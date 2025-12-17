import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------------
# INPUT / OUTPUT
# -----------------------------
CSV_PATH = "log/20251214_172253_glider_ble.csv"  # input
OUT_DIR = "log2"                                 # output folder
OUT_SUFFIX = "_integrated"                       # output filename suffix

# -----------------------------
# ACCEL FIX (bias removal)
# -----------------------------
BIAS_SECONDS = 2.0  # glider MUST be still for the first 2 seconds


def quat_to_euler_zyx(qw, qx, qy, qz):
    """
    Quaternion -> Euler ZYX (yaw-pitch-roll), returns roll, pitch, yaw in radians.
    roll  about X, pitch about Y, yaw about Z.
    """
    # roll
    sinr_cosp = 2.0 * (qw*qx + qy*qz)
    cosr_cosp = 1.0 - 2.0 * (qx*qx + qy*qy)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    # pitch
    sinp = 2.0 * (qw*qy - qz*qx)
    sinp = np.clip(sinp, -1.0, 1.0)
    pitch = np.arcsin(sinp)

    # yaw
    siny_cosp = 2.0 * (qw*qz + qx*qy)
    cosy_cosp = 1.0 - 2.0 * (qy*qy + qz*qz)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return roll, pitch, yaw


def main():
    df = pd.read_csv(CSV_PATH)

    required_cols = [
        "t_imu_ms",
        "qw", "qx", "qy", "qz",
        "ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2",
    ]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    # --- time (IMU) ---
    t_imu_ms = df["t_imu_ms"].to_numpy(dtype=float)
    t = t_imu_ms / 1000.0  # seconds
    N = len(t)
    if N < 2:
        raise ValueError("Not enough samples.")

    # dt (fix non-monotonic/repeated timestamps)
    dt = np.diff(t)
    good = dt > 0
    dt_med = float(np.median(dt[good])) if np.any(good) else 0.01
    dt = np.where(dt > 0, dt, dt_med)
    dt_full = np.concatenate([[dt_med], dt])  # aligned with samples

    # -----------------------------
    # ACCEL BIAS removal (main fix)
    # -----------------------------
    N_bias = int(BIAS_SECONDS / dt_med)
    N_bias = max(1, min(N_bias, N))

    bx = df["ax_lin_mps2"].iloc[:N_bias].mean()
    by = df["ay_lin_mps2"].iloc[:N_bias].mean()
    bz = df["az_lin_mps2"].iloc[:N_bias].mean()

    df["ax_lin_mps2"] -= bx
    df["ay_lin_mps2"] -= by
    df["az_lin_mps2"] -= bz

    # --- accel (m/s^2) ---
    ax = df["ax_lin_mps2"].to_numpy(dtype=float)
    ay = df["ay_lin_mps2"].to_numpy(dtype=float)
    az = df["az_lin_mps2"].to_numpy(dtype=float)
    a = np.vstack([ax, ay, az]).T  # (N,3)

    # --- quaternions ---
    qw = df["qw"].to_numpy(dtype=float)
    qx = df["qx"].to_numpy(dtype=float)
    qy = df["qy"].to_numpy(dtype=float)
    qz = df["qz"].to_numpy(dtype=float)

    # --- orientation (Euler) ---
    roll, pitch, yaw = quat_to_euler_zyx(qw, qx, qy, qz)
    yaw_unwrapped = np.unwrap(yaw)
    rad2deg = 180.0 / np.pi
    roll_deg = roll * rad2deg
    pitch_deg = pitch * rad2deg
    yaw_deg = yaw_unwrapped * rad2deg

    # --- integrate ---
    v = np.zeros((N, 3), dtype=float)
    x = np.zeros((N, 3), dtype=float)

    

    for i in range(1, N):
        v[i] = v[i - 1] + a[i - 1] * dt_full[i]
        x[i] = x[i - 1] + v[i - 1] * dt_full[i]

    # --- métricas em X (AGORA sim, x já existe) ---
    x_net  = float(x[-1, 0] - x[0, 0])                 # deslocamento final (signed)
    x_path = float(np.sum(np.abs(np.diff(x[:, 0]))))   # distância percorrida só em x

    print(f"X net displacement = {x_net:.3f} m")
    print(f"X path distance    = {x_path:.3f} m")


    speed = np.linalg.norm(v, axis=1)
    dist = np.zeros(N, dtype=float)
    for i in range(1, N):
        dist[i] = dist[i - 1] + speed[i - 1] * dt_full[i]

    total_dist = float(dist[-1])
    total_time = float(t[-1] - t[0])
    print(f"Total time ≈ {total_time:.3f} s")
    print(f"Distance (path length) ≈ {total_dist:.3f} m")
    print(f"dt_med ≈ {dt_med:.6f} s  -> freq ≈ {1.0/dt_med:.2f} Hz")
    print(f"Bias removed (m/s^2): bx={bx:.6f}, by={by:.6f}, bz={bz:.6f}")

    # -----------------------------
    # SAVE OUTPUT CSV (log2/)
    # -----------------------------
    os.makedirs(OUT_DIR, exist_ok=True)

    out = pd.DataFrame({
        "t_imu_s": t,
        "t_imu_ms": t_imu_ms,
        "dt_s": dt_full,

        "qw": qw, "qx": qx, "qy": qy, "qz": qz,
        "roll_deg": roll_deg,
        "pitch_deg": pitch_deg,
        "yaw_deg_unwrapped": yaw_deg,

        "ax_lin_mps2": ax,
        "ay_lin_mps2": ay,
        "az_lin_mps2": az,

        "vx_mps": v[:, 0],
        "vy_mps": v[:, 1],
        "vz_mps": v[:, 2],
        "speed_mps": speed,

        "x_m": x[:, 0],
        "y_m": x[:, 1],
        "z_m": x[:, 2],
        "dist_m": dist,
    })

    in_base = os.path.splitext(os.path.basename(CSV_PATH))[0]
    out_path = os.path.join(OUT_DIR, f"{in_base}{OUT_SUFFIX}.csv")

    with open(out_path, "w", newline="") as f:
        f.write(f"# source_file={os.path.basename(CSV_PATH)}\n")
        f.write(f"# total_time_s={total_time:.6f}\n")
        f.write(f"# distance_path_m={total_dist:.6f}\n")
        out.to_csv(f, index=False)

    print(f"Saved: {out_path}")

    # -----------------------------
    # PLOTS
    # -----------------------------
    plt.figure()
    plt.plot(t, dist)
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("Distância percorrida (path) [m]")
    plt.title("Distância percorrida vs tempo")
    plt.grid(True)

    plt.figure()
    plt.plot(t, speed)
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("|v| [m/s]")
    plt.title("Velocidade (magnitude) vs tempo")
    plt.grid(True)

    plt.figure()
    plt.plot(t, v[:, 0], label="vx")
    plt.plot(t, v[:, 1], label="vy")
    plt.plot(t, v[:, 2], label="vz")
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("Velocidade [m/s]")
    plt.title("Componentes da velocidade vs tempo")
    plt.grid(True)
    plt.legend()

    plt.figure()
    plt.plot(t, roll_deg, label="roll")
    plt.plot(t, pitch_deg, label="pitch")
    plt.plot(t, yaw_deg, label="yaw (unwrap)")
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("Ângulo [deg]")
    plt.title("Orientação (Euler) vs tempo")
    plt.grid(True)
    plt.legend()

    plt.figure()
    plt.plot(x[:, 0], x[:, 1])
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Trajetória 2D (XY)")
    plt.axis("equal")
    plt.grid(True)

    plt.figure()
    plt.plot(x[:, 0], x[:, 2])
    plt.xlabel("x [m]")
    plt.ylabel("z [m]")
    plt.title("Trajetória 2D (XZ)")
    plt.axis("equal")
    plt.grid(True)

    plt.figure()
    plt.plot(x[:, 1], x[:, 2])
    plt.xlabel("y [m]")
    plt.ylabel("z [m]")
    plt.title("Trajetória 2D (YZ)")
    plt.axis("equal")
    plt.grid(True)

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
