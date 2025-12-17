import pandas as pd
import matplotlib.pyplot as plt

CSV_FILE = "glider_log.csv"

def main():
    df = pd.read_csv(CSV_FILE)

    # tempo em segundos
    t = df["t_ms"].values / 1000.0

    # --- Orientação: roll, pitch, yaw ---
    plt.figure()
    plt.plot(t, df["roll_deg"], label="roll [deg]")
    plt.plot(t, df["pitch_deg"], label="pitch [deg]")
    plt.plot(t, df["yaw_deg"], label="yaw [deg]")
    plt.xlabel("time [s]")
    plt.ylabel("angle [deg]")
    plt.title("Orientation (Euler angles)")
    plt.legend()
    plt.grid(True)

    # --- Aceleração linear (sem gravidade) ---
    plt.figure()
    plt.plot(t, df["ax_lin_mps2"], label="ax_lin [m/s²]")
    plt.plot(t, df["ay_lin_mps2"], label="ay_lin [m/s²]")
    plt.plot(t, df["az_lin_mps2"], label="az_lin [m/s²]")
    plt.xlabel("time [s]")
    plt.ylabel("linear acceleration [m/s²]")
    plt.title("Linear acceleration (gravity removed)")
    plt.legend()
    plt.grid(True)

    # --- Opcional: módulo da aceleração linear ---
    ax = df["ax_lin_mps2"].values
    ay = df["ay_lin_mps2"].values
    az = df["az_lin_mps2"].values
    a_norm = (ax**2 + ay**2 + az**2) ** 0.5

    plt.figure()
    plt.plot(t, a_norm)
    plt.xlabel("time [s]")
    plt.ylabel("|a_lin| [m/s²]")
    plt.title("Linear acceleration magnitude")
    plt.grid(True)

    plt.show()

if __name__ == "__main__":
    main()
