import asyncio
import csv
import struct
import time
from datetime import datetime
import os
import glob

from bleak import BleakScanner, BleakClient

# Beep on Windows; on other OS we use a simple fallback
try:
    import winsound
    HAVE_WINSOUND = True
except ImportError:
    HAVE_WINSOUND = False

MAX_LOG_FILES = 10
LOG_DIR = "log"

# Alert moments based on the Arduino clock (ms)
ALERT1_T_IMU_MS = 120000.0   # 2 min
ALERT2_T_IMU_MS = 149000.0   # 2.5 min (end of window)


def prepare_log_file() -> str:
    os.makedirs(LOG_DIR, exist_ok=True)

    pattern = os.path.join(LOG_DIR, "*.csv")
    files = glob.glob(pattern)

    if len(files) >= MAX_LOG_FILES:
        files.sort(key=os.path.getmtime)
        num_to_delete = len(files) - (MAX_LOG_FILES - 1)
        for old_path in files[:num_to_delete]:
            try:
                os.remove(old_path)
                print(f"[LOG] Removed old log: {old_path}")
            except OSError as e:
                print(f"[LOG] Error removing {old_path}: {e}")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{ts}_glider_ble.csv"
    full_path = os.path.join(LOG_DIR, filename)
    return full_path


DEVICE_NAME  = "GliderIMU"
SERVICE_UUID = "12345678-1234-1234-1234-1234567890AB"
CHAR_UUID    = "12345678-1234-1234-1234-1234567890AC"

# Payload: [t_ms, qw, qx, qy, qz, ax_lin_mps2, ay_lin_mps2, az_lin_mps2]
N_FLOATS = 8


async def main():
    print(f"[+] Searching for BLE device '{DEVICE_NAME}'...")

    device = await BleakScanner.find_device_by_filter(
        lambda d, ad: d.name == DEVICE_NAME
    )

    if device is None:
        print("[-] Device not found. Ensure it is powered and advertising.")
        return

    print(f"[+] Found: {device.name} ({device.address})")
    print("[+] Connecting...")

    async with BleakClient(device) as client:
        print("[+] Connected to BLE device.")

        csv_path = prepare_log_file()
        print(f"[+] Writing CSV to: {os.path.abspath(csv_path)}")

        with open(csv_path, mode="w", newline="") as f:
            writer = csv.writer(f)

            writer.writerow([
                "t_imu_ms",
                "pc_timestamp_iso",
                "pc_timestamp_s",
                "qw", "qx", "qy", "qz",
                "ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2",
            ])
            f.flush()

            start_time = time.time()
            alert1_done = False
            alert2_done = False

            def do_beep(freq, dur_ms, label):
                print(f"[ALERT] {label}")
                try:
                    if HAVE_WINSOUND:
                        winsound.Beep(freq, dur_ms)
                    else:
                        print("\a", end="", flush=True)
                except Exception as e:
                    print(f"[ALERT] Beep failed: {e}")

            def notification_handler(sender: int, data: bytearray):
                nonlocal writer, start_time, f, alert1_done, alert2_done

                expected_len = 4 * N_FLOATS
                if len(data) != expected_len:
                    print(
                        f"[-] Unexpected size: {len(data)} bytes "
                        f"(expected {expected_len})"
                    )
                    return

                vals = struct.unpack("<" + "f" * N_FLOATS, data)
                (
                    t_imu_ms,
                    qw, qx, qy, qz,
                    ax_lin_mps2, ay_lin_mps2, az_lin_mps2,
                ) = vals

                # PC time (continua guardado no CSV se quiseres)
                t_pc  = time.time()
                t_rel = t_pc - start_time
                iso   = datetime.fromtimestamp(t_pc).isoformat()

                # Tempo do IMU em segundos
                t_imu_s = t_imu_ms / 1000.0

                # ALERT using IMU clock
                if (not alert1_done) and (t_imu_ms >= ALERT1_T_IMU_MS):
                    do_beep(2000, 400, "IMU reached 120 s (start log window).")
                    alert1_done = True

                if (not alert2_done) and (t_imu_ms >= ALERT2_T_IMU_MS):
                    do_beep(1500, 600, "IMU reached 150 s (end log window).")
                    alert2_done = True

                # CSV: mantém pc_timestamp_s como antes, e t_imu_ms em ms
                writer.writerow([
                    f"{t_imu_ms:.1f}",
                    iso,
                    f"{t_rel:.6f}",
                    f"{qw:.6f}", f"{qx:.6f}", f"{qy:.6f}", f"{qz:.6f}",
                    f"{ax_lin_mps2:.6f}", f"{ay_lin_mps2:.6f}", f"{az_lin_mps2:.6f}",
                ])
                f.flush()

                # Terminal: agora mostra o tempo do IMU em segundos
                print(
                    f"{t_imu_s:8.3f}s (IMU)  "
                    f"q=({qw:+.3f},{qx:+.3f},{qy:+.3f},{qz:+.3f})  "
                    f"a_lin=({ax_lin_mps2:+.3f},{ay_lin_mps2:+.3f},{az_lin_mps2:+.3f})"
                )

            print(f"[+] Subscribing to notifications on characteristic {CHAR_UUID}...")
            await client.start_notify(CHAR_UUID, notification_handler)

            print("[+] Receiving data. Ctrl+C to stop.")
            try:
                while True:
                    await asyncio.sleep(1.0)
            except KeyboardInterrupt:
                print("\n[+] Stopping notifications...")
                await client.stop_notify(CHAR_UUID)

        print(f"[+] Log finished. CSV file: {csv_path}")


if __name__ == "__main__":
    asyncio.run(main())
