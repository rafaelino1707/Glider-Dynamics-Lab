import serial
import csv
import time

PORT = "COM7"         # ajusta se for outra COM
BAUD = 115200         # igual ao do Nano
OUTPUT = "ram_dump.csv"

def main():
    with serial.Serial(PORT, BAUD, timeout=1) as ser, \
         open(OUTPUT, "w", newline="") as f:

        writer = csv.writer(f)
        writer.writerow(["raw_line"])  # cabeçalho

        # limpa lixo que possa estar no buffer
        ser.reset_input_buffer()

        print(f"[+] Porta aberta {PORT} @ {BAUD}")
        print("[+] A enviar comando 'D' para pedir o dump...\n")

        # aqui é o equivalente a escrever D no Serial Monitor
        ser.write(b"D")      # se no teu código precisa de ENTER, usa b"D\n"
        ser.flush()

        # agora lemos até ficar em silêncio durante alguns segundos
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
                    # se não veio nada durante max_quiet_s, assumimos que acabou o dump
                    if now - last_rx > max_quiet_s:
                        print("\n[+] Sem dados há mais de "
                              f"{max_quiet_s:.1f}s, a assumir fim do dump.")
                        break

        except KeyboardInterrupt:
            print("\n[+] Interrompido pelo utilizador.")

        print(f"[+] Terminado. Ficheiro guardado em: {OUTPUT}")

if __name__ == "__main__":
    main()
