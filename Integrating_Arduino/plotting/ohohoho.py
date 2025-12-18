import asyncio
import csv
import os
import struct
import time
from datetime import datetime
from bleak import BleakClient, BleakScanner

LOG_DIR = "log"

DEVICE_NAME  = "GliderIMU"
SERVICE_UUID = "12345678-1234-1234-1234-1234567890AB"
IMU_CHAR_UUID = "12345678-1234-1234-1234-1234567890AC"  # 14 floats
NAV_CHAR_UUID = "12345678-1234-1234-1234-1234567890AD"  # 10 floats

N_IMU = 14
N_NAV = 10

IMU_FIELDS = [
    "t_s",
    "qw","qx","qy","qz",
    "axw_mps2","ayw_mps2","azw_mps2",
    "gx_rads","gy_rads","gz_rads",
    "roll_deg","pitch_deg","yaw_deg",
]
NAV_FIELDS = [
    "vx_mps","vy_mps","vz_mps",
    "px_m","py_m","pz_m",
    "still","zupt","rejected","jerk_mps3",
]

PRINT_HZ = 10.0  # mete 50.0 se quiseres spam total

def prepare_log_file() -> str:
    os.makedirs(LOG_DIR, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(LOG_DIR, f"{ts}_glider_ble_nav.csv")

def unpack_floats(data: bytes, n: int):
    if len(data) != 4*n:
        return None
    return struct.unpack("<" + "f"*n, data)

async def find_device(timeout_s=10.0):
    # Preferir service UUID; fallback para name/local_name
    dev = await BleakScanner.find_device_by_filter(
        lambda d, ad: (
            (ad and ad.service_uuids and SERVICE_UUID.lower() in [s.lower() for s in ad.service_uuids])
            or (ad and ad.local_name == DEVICE_NAME)
            or (d and d.name == DEVICE_NAME)
        ),
        timeout=timeout_s
    )
    return dev

async def connect_with_retry(max_tries=6):
    for k in range(max_tries):
        dev = await find_device(timeout_s=8.0)
        if dev is None:
            wait = min(2 + k, 8)
            print(f"[-] Device not found. Retry in {wait}s...")
            await asyncio.sleep(wait)
            continue

        print(f"[+] Found: name={dev.name!r} address={dev.address}")
        try:
            client = BleakClient(dev, timeout=30.0)
            await client.connect()
            return client
        except Exception as e:
            wait = min(2 + k, 8)
            print(f"[-] Connect failed: {type(e).__name__}: {e}. Retry in {wait}s...")
            try:
                await client.disconnect()
            except Exception:
                pass
            await asyncio.sleep(wait)

    raise RuntimeError("Failed to connect after retries.")

async def main():
    csv_path = prepare_log_file()
    print("[+] CSV:", os.path.abspath(csv_path))

    latest_imu = None
    last_print = 0.0
    print_period = 1.0 / max(1e-6, PRINT_HZ)

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(IMU_FIELDS + NAV_FIELDS)
        f.flush()

        client = await connect_with_retry()
        print("[+] Connected.")

        def on_disconnect(_):
            print("[-] BLE disconnected.")
        try:
            client.set_disconnected_callback(on_disconnect)
        except Exception:
            pass

        def on_imu(_, data: bytearray):
            nonlocal latest_imu
            vals = unpack_floats(bytes(data), N_IMU)
            if vals is None:
                print(f"[-] IMU size {len(data)} != {4*N_IMU}")
                return
            latest_imu = vals

        def on_nav(_, data: bytearray):
            nonlocal latest_imu, last_print
            nav = unpack_floats(bytes(data), N_NAV)
            if nav is None:
                print(f"[-] NAV size {len(data)} != {4*N_NAV}")
                return
            if latest_imu is None:
                return

            w.writerow(list(latest_imu) + list(nav))
            f.flush()

            now = time.time()
            if (now - last_print) >= print_period:
                last_print = now

                (t_s, qw,qx,qy,qz,
                 axw,ayw,azw,
                 gx,gy,gz,
                 roll,pitch,yaw) = latest_imu

                (vx,vy,vz, px,py,pz, still,zupt,rejected, jerk) = nav

                print(
                    f"{t_s:8.3f}s  "
                    f"q=({qw:+.3f},{qx:+.3f},{qy:+.3f},{qz:+.3f})  "
                    f"a_w=({axw:+.3f},{ayw:+.3f},{azw:+.3f})  "
                    f"w=({gx:+.3f},{gy:+.3f},{gz:+.3f})  "
                    f"euler=({roll:+.1f},{pitch:+.1f},{yaw:+.1f})  "
                    f"v=({vx:+.2f},{vy:+.2f},{vz:+.2f})  "
                    f"p=({px:+.2f},{py:+.2f},{pz:+.2f})  "
                    f"still={int(still)} zupt={int(zupt)} rej={int(rejected)} jerk={jerk:.1f}"
                )

        await client.start_notify(IMU_CHAR_UUID, on_imu)
        await client.start_notify(NAV_CHAR_UUID, on_nav)

        print("[+] Notifying... Ctrl+C to stop.")
        try:
            while True:
                await asyncio.sleep(1.0)
        except KeyboardInterrupt:
            print("\n[+] Stopping...")
        finally:
            try:
                await client.stop_notify(IMU_CHAR_UUID)
            except Exception:
                pass
            try:
                await client.stop_notify(NAV_CHAR_UUID)
            except Exception:
                pass
            try:
                await client.disconnect()
            except Exception:
                pass

    print("[+] Done:", csv_path)

if __name__ == "__main__":
    asyncio.run(main())
