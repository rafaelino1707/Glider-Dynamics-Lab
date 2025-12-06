import asyncio
from bleak import BleakScanner, BleakClient
import struct
import csv

# UUIDs reais do serviço e characteristic (visto no inspipi.py)
SERVICE_UUID = "0000180d-0000-1000-8000-00805f9b34fb"
CHAR_UUID    = "00002a37-0000-1000-8000-00805f9b34fb"

# 7 floats: qw, qx, qy, qz, ax, ay, az
PACKET_FLOAT_COUNT = 7
STRUCT_FMT = "<" + "f" * PACKET_FLOAT_COUNT

CSV_FILE = "log/glider_ble_log.csv"


async def main():
    print("Scanning for BLE devices...")
    devices = await BleakScanner.discover(timeout=5.0)

    target = None
    for d in devices:
        if "GliderIMU" in (d.name or ""):
            target = d
            break

    if not target:
        print("Device 'GliderIMU' not found.")
        return

    print(f"Found device: {target.address} ({target.name})")
    print("Connecting...")

    async with BleakClient(target.address) as client:
        print("Connected.")
        print(f"Subscribing to characteristic {CHAR_UUID}")

        f = open(CSV_FILE, "w", newline="")
        writer = csv.writer(f)

        # header
        header = ["qw", "qx", "qy", "qz", "ax", "ay", "az"]
        writer.writerow(header)
        f.flush()

        def notification_handler(sender, data: bytearray):
            # cada pacote deve ter 7 floats -> 28 bytes
            if len(data) != 4 * PACKET_FLOAT_COUNT:
                # print("Unexpected packet len:", len(data))
                return

            vals = struct.unpack(STRUCT_FMT, data)
            writer.writerow(vals)
            f.flush()

        await client.start_notify(CHAR_UUID, notification_handler)

        print("Logging BLE data. Press Ctrl+C to stop.")
        try:
            while True:
                await asyncio.sleep(0.2)
        except KeyboardInterrupt:
            print("Stopping logging...")

        await client.stop_notify(CHAR_UUID)
        f.close()
        print("CSV saved:", CSV_FILE)


if __name__ == "__main__":
    asyncio.run(main())
