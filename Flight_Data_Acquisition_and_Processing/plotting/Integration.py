import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting

# -----------------------------
# CONFIGURAÇÃO
# -----------------------------
CSV_PATH = "log/20251211_180620_glider_ble.csv"  # mete aqui o teu ficheiro


def main():
    # --- Ler CSV ---
    df = pd.read_csv(CSV_PATH)

    # Verificar colunas mínimas
    required_cols = [
        "t_imu_ms",
        "qw", "qx", "qy", "qz",
        "ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2",
    ]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Coluna obrigatória em falta no CSV: {c}")

    # --- Tempo em segundos (IMU) ---
    t_imu_ms = df["t_imu_ms"].to_numpy(dtype=float)
    t = t_imu_ms / 1000.0  # [s]

    # dt entre samples
    dt = np.diff(t)
    # evitar dt <= 0 (caso de timestamps repetidos)
    dt[dt <= 0] = np.median(dt[dt > 0]) if np.any(dt > 0) else 0.01

    # --- Aceleração (m/s^2) ---
    ax = df["ax_lin_mps2"].to_numpy(dtype=float)
    ay = df["ay_lin_mps2"].to_numpy(dtype=float)
    az = df["az_lin_mps2"].to_numpy(dtype=float)

    # Empacotar em matriz Nx3
    a = np.vstack([ax, ay, az]).T  # shape (N, 3)
    N = a.shape[0]

    # --- Integração simples para velocidade e posição ---
    v = np.zeros((N, 3))  # velocidade
    x = np.zeros((N, 3))  # posição

    # condição inicial: parte com v = 0, x = 0
    for i in range(1, N):
        # integração de aceleração -> velocidade (Euler direto)
        v[i] = v[i-1] + a[i-1] * dt[i-1]
        # integração de velocidade -> posição
        x[i] = x[i-1] + v[i-1] * dt[i-1]

    # --- Distância percorrida (cumulativa) ---
    # velocidade escalar
    speed = np.linalg.norm(v, axis=1)
    dist = np.zeros(N)
    for i in range(1, N):
        dist[i] = dist[i-1] + speed[i-1] * dt[i-1]

    print(f"Distância final (comprimento de trajetória) ≈ {dist[-1]:.2f} m")

    # -----------------------------
    # PLOTS
    # -----------------------------

    # 1) Distância vs tempo
    plt.figure()
    plt.plot(t, dist)
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("Distância percorrida [m]")
    plt.title("Distância percorrida vs tempo")
    plt.grid(True)

    # 2) Trajetória 2D (projeção XY)
    plt.figure()
    plt.plot(x[:, 0], x[:, 1])
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Trajetória 2D (projeção XY)")
    plt.axis("equal")
    plt.grid(True)

    # 3) Trajetória 3D
    fig = plt.figure()
    ax3d = fig.add_subplot(111, projection="3d")
    ax3d.plot(x[:, 0], x[:, 1], x[:, 2])
    ax3d.set_xlabel("x [m]")
    ax3d.set_ylabel("y [m]")
    ax3d.set_zlabel("z [m]")
    ax3d.set_title("Trajetória 3D (IMU integrada)")

    plt.show()


if __name__ == "__main__":
    main()
