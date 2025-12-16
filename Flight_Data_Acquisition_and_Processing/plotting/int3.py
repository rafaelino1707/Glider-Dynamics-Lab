import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------------
# INPUT / OUTPUT
# -----------------------------
CSV_PATH = "ensaios/20251216_033825_glider_ble.csv"
OUT_DIR = "log2"
OUT_SUFFIX = "_integrated_nozupt_cutTail"

# -----------------------------
# SETTINGS (NO ZUPT)
# -----------------------------
BIAS_SECONDS = 2.0

# Filters on linear accel (gravity removed)
LPF_CUTOFF_HZ = 7.0
HAMPEL_WIN_S  = 0.15
HAMPEL_NSIGMA = 4.0

# Drift-killers WITHOUT ZUPT
HPF_CUTOFF_HZ = 0.10
VEL_LEAK_TAU_S = 6.0

# Jerk gate
JERK_MAX = 80.0
JERK_WIN = 2

# Timestamp glitch protection
DT_MIN = 1e-4
DT_MAX = 0.10

# -----------------------------
# CUT LAST X SECONDS (your fix)
# -----------------------------
CUT_TAIL_SECONDS = 0      # cuts last 4–5 s (where peak happens)
MIN_KEEP_SECONDS = 2.0      # safety: never keep less than this


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


def hpf_1pole(x, t, fc_hz):
    """
    1st-order high-pass filter.
    y[k] = alpha * (y[k-1] + x[k] - x[k-1])
    alpha = 1 / (1 + 2*pi*fc*dt)
    """
    x = np.asarray(x)
    y = np.zeros_like(x)
    y[0] = 0.0
    for k in range(1, len(x)):
        dt = t[k] - t[k-1]
        if dt <= 0:
            y[k] = y[k-1]
            continue
        alpha = 1.0 / (1.0 + 2*np.pi*fc_hz*dt)
        y[k] = alpha * (y[k-1] + (x[k] - x[k-1]))
    return y


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


def build_time_from_t_imu_ms(t_imu_ms):
    t = (t_imu_ms - t_imu_ms[0]) / 1000.0
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

    return t, dt_full, dt_med


def main():
    df_full = pd.read_csv(CSV_PATH)

    required_cols = ["t_imu_ms", "qw", "qx", "qy", "qz",
                     "ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2"]
    for c in required_cols:
        if c not in df_full.columns:
            raise ValueError(f"Missing required column: {c}")

    # -----------------------------
    # FULL TIMEBASE (ALWAYS from t_imu_ms)
    # -----------------------------
    t_imu_ms_full = df_full["t_imu_ms"].to_numpy(dtype=float)
    t_full, dt_full, dt_med = build_time_from_t_imu_ms(t_imu_ms_full)
    N_full = len(t_full)

    total_time_full = float(t_full[-1] - t_full[0])
    print(f"Total time (FULL, via t_imu_ms) ≈ {total_time_full:.3f} s")

    # -----------------------------
    # CUT LAST X SECONDS (remove the peak region)
    # -----------------------------
    cut_t = max(total_time_full - CUT_TAIL_SECONDS, MIN_KEEP_SECONDS)
    keep_mask = t_full <= cut_t
    n_keep = int(np.sum(keep_mask))
    print(f"[CUT] Removing last {CUT_TAIL_SECONDS:.3f} s -> cut at t={cut_t:.3f}s (kept {n_keep}/{N_full} samples)")

    df = df_full.loc[keep_mask].reset_index(drop=True)
    t = t_full[keep_mask]
    dt = dt_full[keep_mask]
    t_imu_ms = t_imu_ms_full[keep_mask]
    N = len(t)

    # -----------------------------
    # QUATS (orientation plots)
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
    # ACCEL (linear, gravity removed)
    # -----------------------------
    a_lin = df[["ax_lin_mps2", "ay_lin_mps2", "az_lin_mps2"]].to_numpy(dtype=float)

    # Initial bias removal (coarse)
    N_bias = int(BIAS_SECONDS / dt_med)
    N_bias = max(1, min(N_bias, N))
    b0 = np.mean(a_lin[:N_bias, :], axis=0)
    a0 = a_lin - b0

    # LPF then Hampel (clean the signal)
    a_lpf = lpf_1pole(a0, t, LPF_CUTOFF_HZ)
    a_clean = hampel_3d(a_lpf, t, win_s=HAMPEL_WIN_S, nsigma=HAMPEL_NSIGMA)

    # High-pass (remove slow drift/bias)
    a_hpf = hpf_1pole(a_clean, t, HPF_CUTOFF_HZ)

    # Jerk flags (on accel used for integration)
    jerk_mag, jerk_isolated = compute_jerk_flags(a_hpf, dt, JERK_MAX, JERK_WIN)

    # -----------------------------
    # INTEGRATION (NO ZUPT)
    # - leaky integrator
    # - trapezoidal position update
    # -----------------------------
    v = np.zeros((N, 3), dtype=float)
    x = np.zeros((N, 3), dtype=float)

    for i in range(1, N):
        dt_i = dt[i]

        # ignore isolated jerk spikes (glitches)
        a_i = np.zeros(3) if jerk_isolated[i-1] else a_hpf[i-1]

        leak = np.exp(-dt_i / VEL_LEAK_TAU_S) if (VEL_LEAK_TAU_S and VEL_LEAK_TAU_S > 1e-6) else 1.0

        v_prev = v[i-1]
        v_i = leak * v_prev + a_i * dt_i
        v[i] = v_i

        x[i] = x[i-1] + 0.5 * (v_prev + v_i) * dt_i

    speed = np.linalg.norm(v, axis=1)

    dist = np.zeros(N, dtype=float)
    for i in range(1, N):
        dist[i] = dist[i-1] + speed[i-1] * dt[i]

    # -----------------------------
    # METRICS
    # -----------------------------
    total_time_used = float(t[-1] - t[0]) if N >= 2 else 0.0
    total_dist = float(dist[-1]) if N else 0.0
    x_net = float(x[-1, 0] - x[0, 0]) if N else 0.0
    x_path = float(np.sum(np.abs(np.diff(x[:, 0])))) if N >= 2 else 0.0

    print(f"Total time (USED after cut) ≈ {total_time_used:.3f} s")
    print(f"Distance (path length) ≈ {total_dist:.3f} m")
    print(f"X net displacement = {x_net:.3f} m")
    print(f"X path distance    = {x_path:.3f} m")
    print(f"dt_med ≈ {dt_med:.6f} s  -> freq ≈ {1.0/dt_med:.2f} Hz")
    print(f"Initial bias removed (m/s^2): bx={b0[0]:.6f}, by={b0[1]:.6f}, bz={b0[2]:.6f}")
    print(f"Cut tail: {CUT_TAIL_SECONDS:.3f} s")

    # -----------------------------
    # SAVE OUTPUT CSV
    # -----------------------------
    os.makedirs(OUT_DIR, exist_ok=True)
    in_base = os.path.splitext(os.path.basename(CSV_PATH))[0]
    out_path = os.path.join(OUT_DIR, f"{in_base}{OUT_SUFFIX}.csv")

    out = pd.DataFrame({
        "t_imu_s": t,
        "t_imu_ms": t_imu_ms,
        "dt_s": dt,

        "qw": qw, "qx": qx, "qy": qy, "qz": qz,
        "roll_deg": roll_deg,
        "pitch_deg": pitch_deg,
        "yaw_deg_unwrapped": yaw_deg,

        "ax_lin_mps2_raw": a_lin[:, 0],
        "ay_lin_mps2_raw": a_lin[:, 1],
        "az_lin_mps2_raw": a_lin[:, 2],

        "ax_lin_mps2_clean": a_clean[:, 0],
        "ay_lin_mps2_clean": a_clean[:, 1],
        "az_lin_mps2_clean": a_clean[:, 2],

        "ax_lin_mps2_used": a_hpf[:, 0],
        "ay_lin_mps2_used": a_hpf[:, 1],
        "az_lin_mps2_used": a_hpf[:, 2],

        "jerk_mag": jerk_mag,
        "jerk_isolated": jerk_isolated.astype(int),

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
        f.write(f"# total_time_full_s={total_time_full:.6f}\n")
        f.write(f"# cut_tail_seconds={CUT_TAIL_SECONDS:.6f}\n")
        f.write(f"# cut_time_s={cut_t:.6f}\n")
        f.write(f"# note=NO_ZUPT; cut last N seconds to remove end-peak region\n")
        out.to_csv(f, index=False)

    print(f"Saved: {out_path}")

    # -----------------------------
    # PLOTS (all the original ones)
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
    plt.plot(t, a_clean[:, 2], label="az lin clean (LPF+Hampel)")
    plt.plot(t, a_hpf[:, 2], label="az lin used (HPF)")
    plt.xlabel("Tempo IMU [s]")
    plt.ylabel("a_lin_z [m/s²]")
    plt.title("Linear acceleration Z (raw vs clean vs used)")
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
    ax3d.set_title("Trajetória 3D (integrada)")

    plt.show()


if __name__ == "__main__":
    main()
