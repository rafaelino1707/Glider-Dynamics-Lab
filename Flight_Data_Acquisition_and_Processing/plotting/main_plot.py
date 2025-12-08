import asyncio
import csv
import struct
import time
from datetime import datetime
import os
import glob

from bleak import BleakScanner, BleakClient

MAX_LOG_FILES = 10
LOG_DIR = "log"


def prepare_log_file() -> str:
    """
    Garante que a pasta log/ existe, apaga o ficheiro mais antigo
    se já houver mais de MAX_LOG_FILES, e devolve o caminho para
    um novo ficheiro CSV com timestamp no nome.
    """
    os.makedirs(LOG_DIR, exist_ok=True)

    # lista todos os CSV na pasta log/
    pattern = os.path.join(LOG_DIR, "*.csv")
    files = glob.glob(pattern)

    # se já há demasiados, apaga os mais antigos
    if len(files) >= MAX_LOG_FILES:
        # ordena por data de modificação (mais antigo primeiro)
        files.sort(key=os.path.getmtime)
        num_to_delete = len(files) - (MAX_LOG_FILES - 1)
        for old_path in files[:num_to_delete]:
            try:
                os.remove(old_path)
                print(f"[LOG] Removido log antigo: {old_path}")
            except OSError as e:
                print(f"[LOG] Erro ao remover {old_path}: {e}")

    # cria nome do novo ficheiro com dia e hora
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{ts}_glider_ble.csv"
    full_path = os.path.join(LOG_DIR, filename)
    return full_path


# Dispositivo / UUIDs (iguais ao firmware)
DEVICE_NAME  = "GliderIMU"
SERVICE_UUID = "12345678-1234-1234-1234-1234567890AB"
CHAR_UUID    = "12345678-1234-1234-1234-1234567890AC"

# Payload: [t_ms, qw, qx, qy, qz, ax, ay, az, gx, gy, gz, mx, my, mz]
N_FLOATS = 14


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

        csv_path = prepare_log_file()
        print(f"[+] A escrever CSV em: {os.path.abspath(csv_path)}")

        # abrir CSV NO CAMINHO CERTO
        with open(csv_path, mode="w", newline="") as f:
            writer = csv.writer(f)

            # cabeçalho
            writer.writerow([
                "t_imu_ms",          # timestamp do Arduino (ms, float)
                "pc_timestamp_iso",  # relógio PC
                "pc_timestamp_s",    # tempo relativo PC
                "qw", "qx", "qy", "qz",
                "ax_raw", "ay_raw", "az_raw",       # em g (ou m/s² se mudares no C++)
                "gx_rad_s", "gy_rad_s", "gz_rad_s",
                "mx_uT", "my_uT", "mz_uT",
            ])
            f.flush()

            start_time = time.time()

            def notification_handler(sender: int, data: bytearray):
                nonlocal writer, start_time, f

                expected_len = 4 * N_FLOATS
                if len(data) != expected_len:
                    print(
                        f"[-] Tamanho inesperado: {len(data)} bytes "
                        f"(esperava {expected_len})"
                    )
                    return

                vals = struct.unpack("<" + "f" * N_FLOATS, data)
                (
                    t_imu_ms,
                    qw, qx, qy, qz,
                    ax, ay, az,
                    gx, gy, gz,
                    mx, my, mz,
                ) = vals

                t_pc  = time.time()
                t_rel = t_pc - start_time
                iso   = datetime.fromtimestamp(t_pc).isoformat()

                # escreve SEMPRE cada pacote
                writer.writerow([
                    f"{t_imu_ms:.1f}",
                    iso,
                    f"{t_rel:.6f}",
                    f"{qw:.6f}", f"{qx:.6f}", f"{qy:.6f}", f"{qz:.6f}",
                    f"{ax:.6f}", f"{ay:.6f}", f"{az:.6f}",
                    f"{gx:.6f}", f"{gy:.6f}", f"{gz:.6f}",
                    f"{mx:.6f}", f"{my:.6f}", f"{mz:.6f}",
                ])
                f.flush()  # força escrever no disco já

                # debug no terminal
                print(
                    f"{t_rel:8.3f}s  "
                    f"q=({qw:+.3f},{qx:+.3f},{qy:+.3f},{qz:+.3f})  "
                    f"a=({ax:+.3f},{ay:+.3f},{az:+.3f})  "
                    f"g=({gx:+.3f},{gy:+.3f},{gz:+.3f})  "
                    f"m=({mx:+.3f},{my:+.3f},{mz:+.3f})"
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

        print(f"[+] Log terminado. Ficheiro CSV: {csv_path}")


if __name__ == "__main__":
    asyncio.run(main())
