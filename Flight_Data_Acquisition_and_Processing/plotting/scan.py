import asyncio
from bleak import BleakScanner

async def main():
    print("[+] A fazer scan BLE (5 s)...")
    devices = await BleakScanner.discover(timeout=5.0)

    if not devices:
        print("[-] Nenhum dispositivo BLE encontrado.")
        return

    print("[+] Dispositivos encontrados:")
    for d in devices:
        print(f"  - name={d.name!r:20}  address={d.address}")

if __name__ == "__main__":
    asyncio.run(main())
