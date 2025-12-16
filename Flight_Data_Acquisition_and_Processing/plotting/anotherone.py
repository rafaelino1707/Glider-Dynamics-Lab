import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------
CSV_PATH = "ensaios/20251216_033100_glider_ble_with_pos.csv"   # <-- mete aqui o path
OUT_DIR  = "plots_out"              # opcional: guarda PNGs

SAVE_PNG = False                    # True para guardar imagens
SHOW_PLOTS = True                   # True para mostrar

# -----------------------------
# LOAD
# -----------------------------
df = pd.read_csv(CSV_PATH)

# Ensure numeric (in case strings)
num_cols = [
    "t_imu_ms","pc_timestamp_s",
    "ax_lin_mps2","ay_lin_mps2","az_lin_mps2",
    "gx_rads","gy_rads","gz_rads",
    "vx_mps","vy_mps","vz_mps",
    "px_m","py_m","pz_m",
    "roll_deg","pitch_deg","yaw_deg"
]
for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["t_imu_ms"]).reset_index(drop=True)

t_s = df["t_imu_ms"].to_numpy() * 1e-3
t_s -= t_s[0]

ax = df["ax_lin_mps2"].to_numpy()
ay = df["ay_lin_mps2"].to_numpy()
az = df["az_lin_mps2"].to_numpy()

vx = df["vx_mps"].to_numpy()
vy = df["vy_mps"].to_numpy()
vz = df["vz_mps"].to_numpy()

px = df["px_m"].to_numpy()
py = df["py_m"].to_numpy()
pz = df["pz_m"].to_numpy()

# -----------------------------
# METRICS
# -----------------------------
a_mag = np.sqrt(ax**2 + ay**2 + az**2)
v_mag = np.sqrt(vx**2 + vy**2 + vz**2)

# distance from position samples (arc length)
dp = np.sqrt(np.diff(px)**2 + np.diff(py)**2 + np.diff(pz)**2)
dist_from_pos = np.nansum(dp)

# distance by integrating speed over time (independent check)
dt = np.diff(t_s)
dist_from_speed = np.nansum(0.5 * (v_mag[1:] + v_mag[:-1]) * dt)

# net displacement
net_disp = np.sqrt((px[-1]-px[0])**2 + (py[-1]-py[0])**2 + (pz[-1]-pz[0])**2)

print("\n--- SUMMARY ---")
print(f"Samples: {len(df)}")
print(f"Duration: {t_s[-1]:.3f} s")
print(f"Total distance (from position): {dist_from_pos:.3f} m")
print(f"Total distance (from speed dt): {dist_from_speed:.3f} m")
print(f"Net displacement: {net_disp:.3f} m")
print(f"Max |a|: {np.nanmax(a_mag):.3f} m/s^2")
print(f"Max |v|: {np.nanmax(v_mag):.3f} m/s")

# -----------------------------
# PLOTS
# -----------------------------
def maybe_save(fig, name):
    if not SAVE_PNG:
        return
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(os.path.join(OUT_DIR, name), dpi=200, bbox_inches="tight")

# 1) Acceleration components + magnitude
fig = plt.figure()
plt.plot(t_s, ax, label="ax")
plt.plot(t_s, ay, label="ay")
plt.plot(t_s, az, label="az")
plt.plot(t_s, a_mag, label="|a|")
plt.xlabel("t [s]")
plt.ylabel("a_lin_world [m/s^2]")
plt.title("Linear acceleration (WORLD)")
plt.legend()
plt.grid(True)
maybe_save(fig, "accel_world.png")

# 2) Velocity components + magnitude
fig = plt.figure()
plt.plot(t_s, vx, label="vx")
plt.plot(t_s, vy, label="vy")
plt.plot(t_s, vz, label="vz")
plt.plot(t_s, v_mag, label="|v|")
plt.xlabel("t [s]")
plt.ylabel("v_world [m/s]")
plt.title("Velocity (WORLD)")
plt.legend()
plt.grid(True)
maybe_save(fig, "vel_world.png")

# 3) Position components
fig = plt.figure()
plt.plot(t_s, px, label="px")
plt.plot(t_s, py, label="py")
plt.plot(t_s, pz, label="pz")
plt.xlabel("t [s]")
plt.ylabel("p_world [m]")
plt.title("Position (WORLD)")
plt.legend()
plt.grid(True)
maybe_save(fig, "pos_world.png")

# 4) 3D trajectory
fig = plt.figure()
ax3 = fig.add_subplot(111, projection="3d")
ax3.plot(px, py, pz)
ax3.set_xlabel("x [m]")
ax3.set_ylabel("y [m]")
ax3.set_zlabel("z [m]")
ax3.set_title("3D trajectory (WORLD)")
maybe_save(fig, "traj_3d.png")

# 5) Total distance over time (cumulative)
cum_dist = np.concatenate([[0.0], np.cumsum(dp)])
fig = plt.figure()
plt.plot(t_s, cum_dist)
plt.xlabel("t [s]")
plt.ylabel("cumulative distance [m]")
plt.title("Cumulative distance (from position)")
plt.grid(True)
maybe_save(fig, "cum_dist.png")

# 6) Euler angles (if present)
if all(c in df.columns for c in ["roll_deg","pitch_deg","yaw_deg"]):
    roll = df["roll_deg"].to_numpy()
    pitch = df["pitch_deg"].to_numpy()
    yaw = df["yaw_deg"].to_numpy()

    fig = plt.figure()
    plt.plot(t_s, roll, label="roll")
    plt.plot(t_s, pitch, label="pitch")
    plt.plot(t_s, yaw, label="yaw")
    plt.xlabel("t [s]")
    plt.ylabel("deg")
    plt.title("Euler angles")
    plt.legend()
    plt.grid(True)
    maybe_save(fig, "euler.png")

# 7) Gyro norm (optional)
if all(c in df.columns for c in ["gx_rads","gy_rads","gz_rads"]):
    gx = df["gx_rads"].to_numpy()
    gy = df["gy_rads"].to_numpy()
    gz = df["gz_rads"].to_numpy()
    w_mag = np.sqrt(gx**2 + gy**2 + gz**2)

    fig = plt.figure()
    plt.plot(t_s, gx, label="gx")
    plt.plot(t_s, gy, label="gy")
    plt.plot(t_s, gz, label="gz")
    plt.plot(t_s, w_mag, label="|w|")
    plt.xlabel("t [s]")
    plt.ylabel("rad/s")
    plt.title("Gyro")
    plt.legend()
    plt.grid(True)
    maybe_save(fig, "gyro.png")

if SHOW_PLOTS:
    plt.show()
