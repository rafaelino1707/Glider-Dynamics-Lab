import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ------- inputs (edit here) -------
Sw        = 0.03770   # wing area [m^2]
bw        = 0.60      # wing span [m]
Vv_target = 0.05      # vertical tail volume target

Cnb_target   = 0.1  # desired dCn/dβ [1/rad]
Cnb_fus      = -0.05  # fuselage contribution
eta_v        = 0.8   # effectiveness
dSigma_dBeta = 0.20   # sidewash gradient
 
# vertical stabilizer geometry for sweep relations
b_v      = 0.052      # fin height [m] (use your value)
c_r      = 0.05      # root chord [m]
c_t      = 0.0305      # tip chord [m]
LambdaLE = 15.0       # LE sweep [deg]      

# aerodynamic constants
e_oswald = 0.9     # tail efficiency
a0       = 2*np.pi    # 2D lift-curve slope [rad^-1]
ARv_assumed = 3.0     # fin AR used in a_v
M        = 0.0        # Mach (keep 0 for incompressible)

lv = np.linspace(0.15, 0.50, 300)  # tail arm [m]; lv = x_v - x_CG
# -----------------------------------

# ---- sweep conversions (linear chord) ----
def quarter_chord_sweep_deg(LambdaLE_deg, cr, ct, b):
    """tan(Lc/4) = tan(LE) + (ct - cr)/(2 b)"""
    LLE = np.deg2rad(LambdaLE_deg)
    tc4 = np.tan(LLE) + (ct - cr)/(2.0*b)
    return np.rad2deg(np.arctan(tc4))

# ---- lift-curve slope with sweep ----
def a_v_from_AR_sweep(AR, Lambda_c4_deg, e=0.92, a0=2*np.pi, M=0.0):
    mu = np.cos(np.deg2rad(Lambda_c4_deg))
    if M <= 0.0:
        return (a0*mu) / (1.0 + (a0*mu)/(np.pi*e*AR))
    # compressible option if ever needed
    beta = np.sqrt(1.0 - (M*mu)**2)
    return (a0*mu) / (beta + (a0*mu)/(np.pi*e*AR))

#Lambda_c4 = quarter_chord_sweep_deg(LambdaLE, c_r, c_t, b_v)
Lambda_c4 = LambdaLE
a_v = a_v_from_AR_sweep(ARv_assumed, Lambda_c4, e_oswald, a0, M)

# ---- gamma from Vv and from Cn_beta ----
gamma_from_Vv  = Vv_target * (bw / lv)
numerator      = (Cnb_target - Cnb_fus) * bw
denom          = eta_v * a_v * lv * (1.0 - dSigma_dBeta)
gamma_from_Cnb = numerator / denom

# guards
gamma_from_Vv  = np.where(gamma_from_Vv  > 0.0, gamma_from_Vv,  np.nan)
gamma_from_Cnb = np.where(gamma_from_Cnb > 0.0, gamma_from_Cnb, np.nan)

# areas
Sv_from_Vv  = gamma_from_Vv  * Sw
Sv_from_Cnb = gamma_from_Cnb * Sw

# intersection (closest)
mask = np.isfinite(gamma_from_Vv) & np.isfinite(gamma_from_Cnb)
idx  = np.nanargmin(np.abs(gamma_from_Vv[mask] - gamma_from_Cnb[mask]))
lv_i  = lv[mask][idx]
gam_i = 0.5*(gamma_from_Vv[mask][idx] + gamma_from_Cnb[mask][idx])
Sv_i  = gam_i * Sw

# save csv
df = pd.DataFrame({
    'l_v [m]': lv,
    'gamma = Sv/Sw (from Vv)': gamma_from_Vv,
    'gamma = Sv/Sw (from Cn_beta)': gamma_from_Cnb,
    'Sv [m^2] (from Vv)': Sv_from_Vv,
    'Sv [m^2] (from Cn_beta)': Sv_from_Cnb
})
df.to_csv('Log/Vertical_Stabilizer/gamma_vs_lv_Vv_and_Cnb.csv', index=False)
print('Saved: gamma_vs_lv_Vv_and_Cnb.csv')
print(f'Λ_LE = {LambdaLE:.2f} deg, Λ_c/4 = {Lambda_c4:.2f} deg, a_v = {a_v:.3f} rad^-1')
print(f'Intersection ~ l_v={lv_i:.3f} m, gamma={gam_i:.3f}, Sv={Sv_i*1e4:.1f} cm^2')

# plot
plt.figure(figsize=(8,5))
#plt.plot(lv, gamma_from_Vv,  label=rf'$\gamma$ from $V_v$ (target={Vv_target:.2f})', lw=2)
plt.plot(lv, gamma_from_Cnb, label=rf'$\gamma$ from $C_{{n_\beta}}$ (target={Cnb_target:.2f})', lw=2)
#plt.scatter([lv_i], [gam_i], s=25, zorder=5,
           # label=rf'Chosen: $l_v={lv_i:.3f}$ m, $\gamma={gam_i:.3f}$')
plt.xlabel(r'$l_v=(x_v-x_{CG})\ \mathrm{[m]}$')
plt.ylabel(r'$\gamma=S_v/S_w$')
plt.title(r'$\gamma$ vs $(x_h - x_{CG})$'+rf': $C_{{n_\beta}}={Cnb_target:.2f}$, $AR_v={ARv_assumed:.1f}$')
plt.grid(True, ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
