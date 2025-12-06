import asyncio
from bleak import BleakScanner, BleakClient

async def main():
    print("Scanning...")
    devices = await BleakScanner.discover(timeout=5.0)

    target = None
    for d in devices:
        if "GliderIMU" in (d.name or ""):
            target = d
            break

    if not target:
        print("GliderIMU not found.")
        return

    print(f"Found: {target.address} ({target.name}), connecting...")

    async with BleakClient(target.address) as client:
        print("Connected.")

        # tentar forçar descoberta de serviços se o método existir
        try:
            get_services = getattr(client, "get_services", None)
            if callable(get_services):
                await get_services()
        except Exception as e:
            print("Warning calling get_services():", e)

        svcs = getattr(client, "services", None)
        if svcs is None:
            print("client.services is None – bleak muito antigo.")
            return

        print("Discovered services and characteristics:")
        for service in svcs:
            print(f"Service: {service.uuid}")
            for char in service.characteristics:
                props = ",".join(char.properties)
                print(f"  Char: {char.uuid}  [{props}]")

if __name__ == "__main__":
    asyncio.run(main())
