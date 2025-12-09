import serial
import csv
import time

PORT = "COM7"         # adjust if using another COM
BAUD = 115200         # same as on the Nano
OUTPUT = "ram_dump.csv"

def main():
    with serial.Serial(PORT, BAUD, timeout=1) as ser, \
         open(OUTPUT, "w", newline="") as f:

        writer = csv.writer(f)
        writer.writerow(["raw_line"])  # header

        # clear any garbage that may be in the buffer
        ser.reset_input_buffer()

        print(f"[+] Port opened {PORT} @ {BAUD}")
        print("[+] Sending command 'D' to request RAM dump...\n")

        # this is equivalent to typing 'D' in Serial Monitor
        ser.write(b"D")      # if your firmware needs ENTER, use b"D\n"
        ser.flush()

        # now we read until it stays silent for a few seconds
        max_quiet_s = 3.0
        last_rx = time.time()

        try:
            while True:
                line = ser.readline()
                now = time.time()

                if line:
                    last_rx = now
                    text = line.decode(errors="ignore").rstrip()
                    print(text)
                    writer.writerow([text])
                else:
                    # if nothing arrives for max_quiet_s, assume dump is finished
                    if now - last_rx > max_quiet_s:
                        print("\n[+] No data for more than "
                              f"{max_quiet_s:.1f}s, assuming end of dump.")
                        break

        except KeyboardInterrupt:
            print("\n[+] Interrupted by user.")

        print(f"[+] Done. File saved to: {OUTPUT}")

if __name__ == "__main__":
    main()
