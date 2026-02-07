import serial
import csv
import time

PORT = "COM11"         # ajusta se for outra porta
BAUD = 115200         # igual ao Arduino
OUTPUT = "ram_dump.csv"

def main():
    with serial.Serial(PORT, BAUD, timeout=1) as ser, \
         open(OUTPUT, "w", newline="") as f:

        writer = csv.writer(f)
        writer.writerow(["raw_line"])  # header simples: uma coluna com a linha completa

        # limpa lixo que possa estar no buffer
        ser.reset_input_buffer()

        print(f"[+] Port opened {PORT} @ {BAUD}")
        print("[+] Sending command 'D' to request RAM dump...\n")

        # equivalente a escrever 'D' no Serial Monitor
        ser.write(b"D")      # se o firmware precisasse de ENTER: ser.write(b"D\n")
        ser.flush()

        # usamos os marcadores do Arduino:
        # "=== RAM LOG DUMP BEGIN ==="
        # "=== RAM LOG DUMP END ==="
        in_dump = False
        max_quiet_s = 3.0
        last_rx = time.time()

        try:
            while True:
                line = ser.readline()
                now = time.time()

                if line:
                    text = line.decode(errors="ignore").rstrip()
                    print(text)

                    # atualiza timestamp sempre que chega qualquer coisa
                    last_rx = now

                    # detetar início do dump
                    if "=== RAM LOG DUMP BEGIN" in text:
                        print("[+] BEGIN marker detected, starting to record dump.")
                        in_dump = True
                        # opcionalmente não guardar esta linha no CSV:
                        continue

                    # detetar fim do dump
                    if "=== RAM LOG DUMP END" in text:
                        print("[+] END marker detected, stopping.")
                        # opcionalmente não guardar esta linha no CSV:
                        break

                    # só gravamos linhas entre BEGIN e END
                    if in_dump:
                        writer.writerow([text])
                        f.flush()
                else:
                    # se estivermos em modo dump e ficar tudo silencioso X s, também paramos
                    if in_dump and (now - last_rx > max_quiet_s):
                        print("\n[+] No data for more than "
                              f"{max_quiet_s:.1f}s while dumping, assuming end.")
                        break

        except KeyboardInterrupt:
            print("\n[+] Interrupted by user.")

        print(f"[+] Done. File saved to: {OUTPUT}")

if __name__ == "__main__":
    main()
