import serial
import time
import csv

PORT = "COM7"     # ajusta se for COM8, etc.
BAUD = 115200
TIMEOUT_S = 120   # tempo máximo à espera do dump

OUTFILE = "flight_ram_log.csv"


def main():
    print(f"Opening {PORT} at {BAUD} baud...")
    with serial.Serial(PORT, BAUD, timeout=1) as ser:
        # pequena espera para o boot estabilizar
        time.sleep(2.0)
        ser.reset_input_buffer()

        print("Sending 's' to arm logging sequence...")
        ser.write(b's\n')
        ser.flush()

        print("Waiting for RAM dump markers...")
        in_dump = False
        header = None
        rows = []

        t0 = time.time()

        while True:
            if time.time() - t0 > TIMEOUT_S:
                print("Timeout waiting for dump.")
                break

            line_bytes = ser.readline()
            if not line_bytes:
                continue

            line = line_bytes.decode(errors="ignore").strip()
            # Opcional: ver tudo
            print("RX:", line)

            if line.startswith("=== RAM LOG DUMP BEGIN"):
                in_dump = True
                header = None
                rows = []
                continue

            if line.startswith("=== RAM LOG DUMP END"):
                print("Dump finished.")
                break

            if in_dump:
                if header is None:
                    # primeira linha dentro do dump é o header CSV
                    header = line.split(',')
                else:
                    rows.append(line.split(','))

        if header and rows:
            with open(OUTFILE, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(rows)
            print(f"Saved {len(rows)} samples to {OUTFILE}")
        else:
            print("No dump captured (no BEGIN/END block detected).")


if __name__ == "__main__":
    main()
