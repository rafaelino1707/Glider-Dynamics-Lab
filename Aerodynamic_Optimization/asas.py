import numpy as np
from math import sin, cos
from scipy.integrate import solve_ivp

# ------------------------------
# Glider parameters (adjust to your report)
# ------------------------------
m   = 0.20      # mass [kg]
S   = 0.0377    # wing area [m^2]
rho = 1.225     # air density [kg/m^3]
g   = 9.81

# Design launch speed
V_stall = 7.1
V_oper  = 8.6
f_rail  = 1.2
Vrail   = max(f_rail * V_stall, V_oper)
print(f"Using V_rail = {Vrail:.2f} m/s")

# ------------------------------
# Aerodynamics at a fixed trim AoA (e.g. alpha ~ 3 deg)
# Get these from your XFOIL polar
# ------------------------------
alpha_trim_deg = 3.0
CL_trim = 1.41   # CL at ~3 deg (example from your polar)
CD_trim = 0.05   # CD at that CL (you must tune this)

def eom(t, y):
    """
    y = [V, gamma, x, z]
    """
    V, gamma, x, z = y

    if V <= 0.5:
        # if it gets too slow, let it just fall
        return [0.0, 0.0, 0.0, -g]

    q = 0.5 * rho * V**2
    L = q * S * CL_trim
    D = q * S * CD_trim

    dVdt    = -D/m - g * np.sin(gamma)
    dgamdt  = (L/(m*V)) - (g/V) * np.cos(gamma)
    dxdt    = V * np.cos(gamma)
    dzdt    = V * np.sin(gamma)

    return [dVdt, dgamdt, dxdt, dzdt]

# Ground-hit event (z=0 downwards)
def hit_ground(t, y):
    return y[3]  # z
hit_ground.terminal  = True
hit_ground.direction = -1

# ------------------------------
# Sweep launch angle gamma0
# ------------------------------
angles_deg = np.arange(5.0, 35.1, 2.5)
results = []

for gam_deg in angles_deg:
    gamma0 = np.deg2rad(gam_deg)

    # initial state at rail exit
    V0 = Vrail
    x0 = 0.0
    z0 = 0.0
    y0 = [V0, gamma0, x0, z0]

    sol = solve_ivp(
        eom,
        t_span=(0.0, 60.0),
        y0=y0,
        events=hit_ground,
        max_step=0.01
    )

    if sol.t_events[0].size > 0:
        T_f = sol.t_events[0][0]
    else:
        T_f = sol.t[-1]

    results.append((gam_deg, T_f))

# Print sweep
for ang, T in results:
    print(f"gamma0 = {ang:5.1f} deg  -> flight time ~ {T:6.2f} s")

best = max(results, key=lambda x: x[1])
print("\nBest angle in this sweep:")
print(f"gamma0 = {best[0]:.1f} deg  -> time ~ {best[1]:.2f} s")
