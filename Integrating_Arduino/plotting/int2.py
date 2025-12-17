import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------------
# INPUT / OUTPUT
# -----------------------------
CSV_PATH = "log/20251214_154520_glider_ble.csv"
OUT_DIR = "log2"
OUT_SUFFIX = "_integrated"

# -----------------------------
# SETTINGS
# -----------------------------
BIAS_SECONDS = 2.0
LPF_CUTOFF_HZ = 7.0

# ZUPT (stationary detection) on LINEAR accel (gravity removed)
ZUPT_ALIN_TOL = 0.35
ZUPT_W_TOL    = 0.15
ZUPT_WIN_S    = 0.25
ZUPT_HOLD_S   = 0.50

# Timestamp glitch protection
DT_MIN = 1e-4
DT_MAX = 0.10

# Despike (Hampel) on accel (kills 1–2 sample glitches, keeps real pulses)
HAMPEL_WIN_S   = 0.15
HAMPEL_NSIGMA  = 4.0

# Jerk gate (blocks isolated insane steps, allows sustained high accel)
JERK_MAX = 80.0
JERK_WIN = 2

# When stationary, force v=0 (aggressive ZUPT)
FORCE_V_ZERO_WHEN_STATIONARY = True


def quat_to_euler_zyx(qw, qx, qy, qz):
    sinr_cosp = 2.0 * (qw*qx + qy*qz)
    cosr_cosp = 1.0 - 2.0 * (qx*qx + qy*qy)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    sinp = 2.0 * (qw*qy - qz*qx)
    sinp = np.clip(sinp, -1.0, 1.0)
    pitch = np.arcsin(sinp)

    siny_cosp = 2.0 * (qw*qz + qx*qy)
    cosy_cosp = 1.0 - 2.0 * (qy*qy + qz*qz)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return roll, pitch, yaw


def lpf_1pole(x, t, fc_hz):
    x = np.asarray(x)
    y = np.zeros_like(x)
    y[0] = x[0]
    for k in range(1, len(x)):
        dt = t[k] - t[k-1]
        if dt <= 0:
            y[k] = y[k-1]
            continue
        alpha = (2*np.pi*fc_hz*dt) / (1 + 2*np.pi*fc_hz*dt)
        y[k] = y[k-1] + alpha*(x[k] - y[k-1])
    return y


def detect_stationary_from_linacc(t, a_lin_mps2, gyro_rads,
                                  alin_tol=ZUPT_ALIN_TOL, w_tol=ZUPT_W_TOL, win_s=ZUPT_WIN_S):
    a_norm = np.linalg.norm(a_lin_mps2, axis=1)
    w_norm = np.linalg.norm(gyro_rads, axis=1)
    raw = (a_norm < alin_tol) & (w_norm < w_tol)

    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    fs = 1.0 / np.mean(dt) if len(dt) else 50.0
    win = max(1, int(win_s * fs))

    stat = np.zeros_like(raw, dtype=bool)
    c = 0
    for i in range(len(raw)):
        c = c + 1 if raw[i] else 0
        stat[i] = (c >= win)
    return stat


def apply_stationary_hold(t, stationary, hold_s):
    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    fs = 1.0 / np.mean(dt) if len(dt) else 50.0
    hold_n = max(1, int(hold_s * fs))

    out = stationary.copy()
    hold = 0
    for i in range(len(out)):
        if stationary[i]:
            hold = hold_n
            out[i] = True
        else:
            if hold > 0:
                out[i] = True
                hold -= 1
    return out


def hampel_1d(x, win, nsigma):
    x = np.asarray(x, dtype=float)
    y = x.copy()
    n = len(x)
    k = 1.4826
    for i in range(n):
        i0 = max(0, i - win)
        i1 = min(n, i + win + 1)
        w = x[i0:i1]
        med = np.median(w)
        mad = k * np.median(np.abs(w - med))
        if mad <= 1e-12:
            continue
        if np.abs(x[i] - med) > nsigma * mad:
            y[i] = med
    return y


def hampel_3d(a, t, win_s, nsigma):
    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    fs = 1.0 / np.mean(dt) if len(dt) else 50.0
    win = max(1, int(win_s * fs))

    out = a.copy()
    for j in range(3):
        out[:, j] = hampel_1d(out[:, j], win=win, nsigma=nsigma)
    return out


def compute_jerk_flags(a_f, dt_full, jerk_max, jerk_win):
    jerk = np.zeros_like(a_f)
    jerk[1:] = (a_f[1:] - a_f[:-1]) / dt_full[1:, None]
    jerk_mag = np.linalg.norm(jerk, axis=1)

    spike = jerk_mag > jerk_max
    isolated = np.zeros_like(spike, dtype=bool)
    run = 0
    N = len(spike)
    for i in range(N):
        if spike[i]:
            run += 1
        else:
            if 1 <= run <= jerk_win:
                isolated[i-run:i] = True
            run = 0
    if 1 <= run <= jerk_win:
        isolated[N-run:N] = True

    return jerk_mag, isolated


def main():
    df = pd.read_csv(CSV_PATH)

    required_cols = ["t_imu_ms", "qw", "qx", "qy", "qz", "ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    # -----------------------------
    # TIME
    # -----------------------------
    t_imu_ms = df["t_imu_ms"].to_numpy(dtype=float)
    t = t_imu_ms / 1000.0
    N = len(t)
    if N < 2:
        raise ValueError("Not enough samples.")

    dt_raw = np.diff(t)
    good = dt_raw > 0
    dt_med = float(np.median(dt_raw[good])) if np.any(good) else 0.01

    t_fix = t.copy()
    for i in range(1, N):
        if (not np.isfinite(t_fix[i])) or (t_fix[i] <= t_fix[i-1]):
            t_fix[i] = t_fix[i-1] + dt_med
    t = t_fix

    dt_full = np.empty(N, dtype=float)
    dt_full[0] = dt_med
    dt_full[1:] = np.diff(t)
    dt_full = np.clip(dt_full, DT_MIN, DT_MAX)

    # -----------------------------
    # QUATS (only for orientation plots)
    # -----------------------------
    qw = df["qw"].to_numpy(dtype=float)
    qx = df["qx"].to_numpy(dtype=float)
    qy = df["qy"].to_numpy(dtype=float)
    qz = df["qz"].to_numpy(dtype=float)

    qnorm = np.sqrt(qw*qw + qx*qx + qy*qy + qz*qz)
    badq = (~np.isfinite(qnorm)) | (qnorm < 0.8) | (qnorm > 1.2)
    qnorm = np.where(qnorm > 0, qnorm, 1.0)
    qw, qx, qy, qz = qw/qnorm, qx/qnorm, qy/qnorm, qz/qnorm
    for i in range(1, N):
        if badq[i]:
            qw[i], qx[i], qy[i], qz[i] = qw[i-1], qx[i-1], qy[i-1], qz[i-1]

    roll, pitch, yaw = quat_to_euler_zyx(qw, qx, qy, qz)
    yaw_unwrapped = np.unwrap(yaw)
    rad2deg = 180.0 / np.pi
    roll_deg = roll * rad2deg
    pitch_deg = pitch * rad2deg
    yaw_deg = yaw_unwrapped * rad2deg

    # -----------------------------
    # ACCEL (linear, already gravity removed)
    # -----------------------------
    a_lin = df[["ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2"]].to_numpy(dtype=float)

    # gyro for ZUPT
    has_gyro = all(c in df.columns for c in ["gx_rads", "gy_rads", "gz_rads"])
    if has_gyro:
        gyro = df[["gx_rads", "gy_rads", "gz_rads"]].to_numpy(dtype=float)
    else:
        gyro = np.zeros((N, 3), dtype=float)
        print("WARNING: gyro columns gx_rads/gy_rads/gz_rads not found -> ZUPT disabled.")

    # bias removal
    N_bias = int(BIAS_SECONDS / dt_med)
    N_bias = max(1, min(N_bias, N))
    b = np.mean(a_lin[:N_bias, :], axis=0)
    a_lin_bias = a_lin - b

    # low-pass + despike
    a_lpf = lpf_1pole(a_lin_bias, t, LPF_CUTOFF_HZ)
    a_f = hampel_3d(a_lpf, t, win_s=HAMPEL_WIN_S, nsigma=HAMPEL_NSIGMA)

    # ZUPT + hold
    if has_gyro:
        stationary0 = detect_stationary_from_linacc(t, a_f, gyro)
        stationary = apply_stationary_hold(t, stationary0, ZUPT_HOLD_S)
    else:
        stationary = np.zeros(N, dtype=bool)

    # jerk flags
    jerk_mag, jerk_isolated = compute_jerk_flags(a_f, dt_full, JERK_MAX, JERK_WIN)

    # -----------------------------
    # INTEGRATION
    # -----------------------------
    v = np.zeros((N, 3), dtype=float)
    x = np.zeros((N, 3), dtype=float)

    for i in range(1, N):
        dt_i = dt_full[i]

        if FORCE_V_ZERO_WHEN_STATIONARY and stationary[i]:
            v[i] = 0.0
            x[i] = x[i-1]
            continue

        # ignore isolated jerk spikes (glitches)
        if jerk_isolated[i-1]:
            dv = np.zeros(3)
        else:
            dv = a_f[i-1] * dt_i

        v[i] = v[i-1] + dv
        x[i] = x[i-1] + v[i] * dt_i

    speed = np.linalg.norm(v, axis=1)
    dist = np.zeros(N, dtype=float)
    for i in range(1, N):
        dist[i] = dist[i-1] + speed[i-1] * dt_full[i]

    # metrics
    x_net = float(x[-1, 0] - x[0, 0])
    x_path = float(np.sum(np.abs(np.diff(x[:, 0]))))
    total_time = float(t[-1] - t[0])
    total_dist = float(dist[-1])

    print(f"Total time ≈ {total_time:.3f} s")
    print(f"Distance (path length) ≈ {total_dist:.3f} m")
    print(f"X net displacement = {x_net:.3f} m")
    print(f"X path distance    = {x_path:.3f} m")
    print(f"dt_med ≈ {dt_med:.6f} s  -> freq ≈ {1.0/dt_med:.2f} Hz")
    print(f"Bias removed (m/s^2): bx={b[0]:.6f}, by={b[1]:.6f}, bz={b[2]:.6f}")
    if has_gyro:
        print(f"Stationary fraction: {100*np.mean(stationary):.1f}%")
    print(f"Hampel: win={HAMPEL_WIN_S}s nsigma={HAMPEL_NSIGMA} | Jerk gate: JERK_MAX={JERK_MAX} win={JERK_WIN}")

    print("mean a_lin:", np.mean(a_lin, axis=0))
    print("mean |a_lin|:", np.mean(np.linalg.norm(a_lin, axis=1)))
    print("median |a_lin|:", np.median(np.linalg.norm(a_lin, axis=1)))
    print("min/max az:", np.min(a_lin[:,2]), np.max(a_lin[:,2]))

    # -----------------------------
    # SAVE OUTPUT CSV
    # -----------------------------
    os.makedirs(OUT_DIR, exist_ok=True)
    in_base = os.path.splitext(os.path.basename(CSV_PATH))[0]
    out_path = os.path.join(OUT_DIR, f"{in_base}{OUT_SUFFIX}.csv")

    out = pd.DataFrame({
        "t_imu_s": t,
        "t_imu_ms": t_imu_ms,
        "dt_s": dt_full,

        "qw": qw, "qx": qx, "qy": qy, "qz": qz,
        "roll_deg": roll_deg,
        "pitch_deg": pitch_deg,
        "yaw_deg_unwrapped": yaw_deg,

        "ax_lin_mps2_raw": a_lin[:, 0],
        "ay_lin_mps2_raw": a_lin[:, 1],
        "az_lin_mps2_raw": a_lin[:, 2],

        "ax_lin_mps2_f": a_f[:, 0],
        "ay_lin_mps2_f": a_f[:, 1],
        "az_lin_mps2_f": a_f[:, 2],

        "jerk_mag": jerk_mag,
        "jerk_isolated": jerk_isolated.astype(int),
        "stationary": stationary.astype(int),

        "vx_mps": v[:, 0],
        "vy_mps": v[:, 1],
        "vz_mps": v[:, 2],
        "speed_mps": speed,

        "x_m": x[:, 0],
        "y_m": x[:, 1],
        "z_m": x[:, 2],
        "dist_m": dist,
    })

    with open(out_path, "w", newline="") as f:
        f.write(f"# source_file={os.path.basename(CSV_PATH)}\n")
        f.write(f"# total_time_s={total_time:.6f}\n")
        f.write(f"# distance_path_m={total_dist:.6f}\n")
        out.to_csv(f, index=False)

    print(f"Saved: {out_path}")

    # -----------------------------
    # PLOTS (mantém os teus + adiciona flags úteis sem substituir)
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
    plt.plot(t, a_lin[:, 2], label="az lin raw")
    plt.plot(t, a_f[:, 2], label="az lin filtered+despiked")
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("a_lin_z [m/s²]")
    plt.title("Linear acceleration Z (raw vs filtered)")
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

    if has_gyro:
        plt.figure()
        plt.plot(t, stationary.astype(int))
        plt.xlabel("Tempo IMU [s]")
        plt.ylabel("stationary (0/1)")
        plt.title("ZUPT stationary detector (with hold)")
        plt.grid(True)

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
    ax3d.set_title("Trajetória 3D (integrada)")

    plt.show()


if __name__ == "__main__":
    main()
