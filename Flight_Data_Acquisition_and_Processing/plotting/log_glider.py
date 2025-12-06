import serial
import time
import sys

# Ajusta aqui se a porta ou baudrate mudarem
PORT = "COM7"
BAUD = 115200
OUTFILE = "log/glider_log.csv"

HEADER = (
    "t_ms,"
    "qw,qx,qy,qz,"
    "ax_mps2,ay_mps2,az_mps2,"
    "gx_deg,gy_deg,gz_deg,"
    "roll_deg,pitch_deg,yaw_deg,"
    "ax_lin_mps2,ay_lin_mps2,az_lin_mps2"
)

def main():
    try:
        with serial.Serial(PORT, BAUD, timeout=1) as ser, open(OUTFILE, "w", newline="") as f:
            print(f"Opened {PORT} at {BAUD} baud")
            print(f"Logging to {OUTFILE}")

            # escreve header uma vez
            f.write(HEADER + "\n")
            f.flush()

            try:
                while True:
                    line = ser.readline().decode(errors="ignore").strip()
                    if not line:
                        continue

                    # opcional: validação simples de formato
                    parts = line.split(',')
                    if len(parts) != 17:
                        # linha estranha, ignora
                        continue

                    f.write(line + "\n")
                    f.flush()

            except KeyboardInterrupt:
                print("\nStopping logging (Ctrl+C pressed).")

    except serial.SerialException as e:
        print(f"Serial error: {e}")
        sys.exit(1)

    print("File closed.")

if __name__ == "__main__":
    main()
