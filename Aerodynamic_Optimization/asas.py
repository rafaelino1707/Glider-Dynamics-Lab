import numpy as np
from math import sin, cos
from scipy.integrate import solve_ivp

# ------------------------------
# Glider parameters (from report)
# ------------------------------
m   = 0.20        # mass [kg]
S   = 0.0377      # wing area [m^2]
rho = 1.225       # air density [kg/m^3]
g   = 9.81

# Wing incidence (body axis -> wing chord)
i_w_deg = 2.0
i_w     = np.deg2rad(i_w_deg)

# Aerodynamic model based on S-1223 data
# CL(3°) = 1.4106, CL(1°) = 1.1825  -> dCL/dalpha ~ 0.114 per deg
CL_alpha_deg = (1.4106 - 1.1825) / (3.0 - 1.0)  # per degree
CL_alpha     = CL_alpha_deg * np.pi / 180.0     # per rad

# Choose reference AoA and CL at that AoA (around best L/D)
alpha_ref_deg = 3.0
alpha_ref     = np.deg2rad(alpha_ref_deg)
CL_ref        = 1.4106

# Simple linear CL(α) around alpha_ref
def CL_of_alpha(alpha):
    return CL_ref + CL_alpha * (alpha - alpha_ref)

# Drag polar of whole glider
CD0  = 0.03   # tune with your polar / tests
K    = 0.040  # induced drag factor (approx from AR and e)

def CD_of_CL(CL):
    return CD0 + K * CL**2

# Choose a fixed pitch attitude theta0 such that in "typical" glide
# you get alpha ~= alpha_ref when gamma is small.
# For a first approximation, set theta0 = alpha_ref - i_w.
theta0 = alpha_ref - i_w

# ------------------------------
# Equations of motion
# ------------------------------
def eom(t, y):
    V, gamma, x, z = y

    if V <= 0.5:
        # avoid numerical issues; once very slow, just fall nearly vertically
        return [0.0, 0.0, 0.0, -g]

    alpha = theta0 - gamma + i_w
    CL    = CL_of_alpha(alpha)
    CD    = CD_of_CL(CL)

    q  = 0.5 * rho * V**2
    L  = q * S * CL
    D  = q * S * CD

    dVdt    = -D/m - g * np.sin(gamma)
    dgamdt  = (L/(m*V)) - (g/V) * np.cos(gamma)
    dxdt    = V * np.cos(gamma)
    dzdt    = V * np.sin(gamma)

    return [dVdt, dgamdt, dxdt, dzdt]

# Ground impact event (z=0 crossing downward)
def hit_ground(t, y):
    return y[3]  # z
hit_ground.terminal  = True
hit_ground.direction = -1

# ------------------------------
# Choose V_rail from design speeds
# ------------------------------
V_stall   = 7.1   # from report
V_oper    = 8.6   # from report
f_rail    = 1.2

Vrail_min = f_rail * V_stall
Vrail     = max(Vrail_min, V_oper)

print(f"Using V_rail = {Vrail:.2f} m/s")

# ------------------------------
# Sweep launch angle gamma0
# ------------------------------
angles_deg = np.arange(5.0, 35.1, 2.5)
results = []

for gam_deg in angles_deg:
    gamma0 = np.deg2rad(gam_deg)

    # Initial conditions at rail exit
    V0 = Vrail
    x0 = 0.0
    z0 = 0.0
    y0 = [V0, gamma0, x0, z0]

    sol = solve_ivp(
        eom,
        t_span=(0.0, 60.0),  # max 60 s
        y0=y0,
        events=hit_ground,
        max_step=0.01
    )

    if sol.t_events[0].size > 0:
        T_f = sol.t_events[0][0]
    else:
        T_f = sol.t[-1]

    results.append((gam_deg, T_f))

# Print results
for ang, T in results:
    print(f"gamma0 = {ang:5.1f} deg  -> flight time ~ {T:6.2f} s")

best = max(results, key=lambda x: x[1])
print("\nBest angle in this sweep:")
print(f"gamma0 = {best[0]:.1f} deg  -> time ~ {best[1]:.2f} s")
