import asyncio
import csv
import struct
import time
from datetime import datetime

from bleak import BleakScanner, BleakClient

# Ajusta isto se usares outro nome / UUID
DEVICE_NAME   = "GliderIMU"
SERVICE_UUID  = "12345678-1234-1234-1234-1234567890AB"
CHAR_UUID     = "12345678-1234-1234-1234-1234567890AC"
N_FLOATS      = 7
CSV_FILENAME  = "glider_ble_log2.csv"



async def main():
    print(f"[+] À procura do dispositivo BLE '{DEVICE_NAME}'...")

    device = await BleakScanner.find_device_by_filter(
        lambda d, ad: d.name == DEVICE_NAME
    )

    if device is None:
        print("[-] Dispositivo não encontrado. Garante que está ligado e a anunciar.")
        return

    print(f"[+] Encontrado: {device.name} ({device.address})")
    print("[+] A ligar...")

    async with BleakClient(device) as client:
        print("[+] Ligado ao dispositivo BLE.")

        # Abre o ficheiro CSV
        with open(CSV_FILENAME, mode="w", newline="") as f:
            writer = csv.writer(f)

            # Cabeçalho
            writer.writerow([
                "pc_timestamp_iso",
                "pc_timestamp_s",
                "qw", "qx", "qy", "qz",
                "ax_mps2_raw", "ay_mps2_raw", "az_mps2_raw"
            ])

            start_time = time.time()

            def notification_handler(sender: int, data: bytearray):
                """
                Esta função é chamada sempre que chegam dados novos por BLE.
                'data' é o pacote de bytes enviado pelo micro.
                """
                nonlocal writer

                # Esperamos N_FLOATS floats de 4 bytes, little-endian
                expected_len = 4 * N_FLOATS
                if len(data) != expected_len:
                    print(f"[-] Tamanho inesperado: {len(data)} bytes (esperava {expected_len})")
                    return

                vals = struct.unpack("<" + "f" * N_FLOATS, data)
                qw, qx, qy, qz, ax, ay, az = vals

                t_pc = time.time()
                t_rel = t_pc - start_time
                iso = datetime.fromtimestamp(t_pc).isoformat()

                # Escreve linha no CSV
                writer.writerow([
                    iso,
                    f"{t_rel:.6f}",
                    f"{qw:.6f}", f"{qx:.6f}", f"{qy:.6f}", f"{qz:.6f}",
                    f"{ax:.6f}", f"{ay:.6f}", f"{az:.6f}",
                ])

                # Opcional: também imprime no terminal
                print(
                    f"{t_rel:8.3f}s  "
                    f"q=({qw:+.3f},{qx:+.3f},{qy:+.3f},{qz:+.3f})  "
                    f"a=({ax:+.3f},{ay:+.3f},{az:+.3f})"
                )

            print(f"[+] A subscrever notificações na characteristic {CHAR_UUID}...")
            await client.start_notify(CHAR_UUID, notification_handler)

            print("[+] A receber dados. Ctrl+C para parar.")
            try:
                while True:
                    await asyncio.sleep(1.0)
            except KeyboardInterrupt:
                print("\n[+] A parar notificações...")
                await client.stop_notify(CHAR_UUID)

        print(f"[+] Log terminado. Ficheiro CSV: {CSV_FILENAME}")


if __name__ == "__main__":
    asyncio.run(main())
