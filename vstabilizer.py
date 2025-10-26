import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ------- inputs (edite aqui) -------
Sw        = 0.03770   # wing area [m^2]
bw        = 0.60      # wing span [m]
Vv_target = 0.05      # vertical tail volume target

Cnb_target   = 0.06   # desired dCn/dβ [1/rad]
Cnb_fus      = -0.05  # fuselage contribution
eta_v        = 0.80   # effectiveness
dSigma_dBeta = 0.20   # sidewash gradient
ARv_assumed  = 3.0    # assumed fin AR for a_v estimate

lv = np.linspace(0.15, 0.50, 300)  # tail arm [m]; lv = x_v - x_CG
# -----------------------------------

def a_v_from_AR(AR):  # lateral lift-curve slope [rad^-1]
    return 2*np.pi * AR / (AR + 2.0)

a_v = a_v_from_AR(ARv_assumed)

# Method 1: from Vv  -> gamma = Sv/Sw
gamma_from_Vv = Vv_target * (bw / lv)

# Method 2: from Cn_beta target (worst-case model)
# Cnβ ≈ η_v a_v (Sv lv)/(Sw bw) (1 - dσ/dβ) + Cnβ_fus
numerator = (Cnb_target - Cnb_fus) * bw
denom     = eta_v * a_v * lv * (1.0 - dSigma_dBeta)
gamma_from_Cnb = numerator / denom

# guards
gamma_from_Vv  = np.where(gamma_from_Vv  > 0.0, gamma_from_Vv,  np.nan)
gamma_from_Cnb = np.where(gamma_from_Cnb > 0.0, gamma_from_Cnb, np.nan)

# areas
Sv_from_Vv  = gamma_from_Vv  * Sw
Sv_from_Cnb = gamma_from_Cnb * Sw

# intersection (closest point where both are finite)
mask = np.isfinite(gamma_from_Vv) & np.isfinite(gamma_from_Cnb)
idx  = np.nanargmin(np.abs(gamma_from_Vv[mask] - gamma_from_Cnb[mask]))
lv_i = lv[mask][idx]
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
df.to_csv('gamma_vs_lv_Vv_and_Cnb.csv', index=False)
print('Saved: gamma_vs_lv_Vv_and_Cnb.csv')
print(f'Intersection ~ l_v={lv_i:.3f} m, gamma={gam_i:.3f}, Sv={Sv_i*1e4:.1f} cm^2')

# plot
plt.figure(figsize=(8,5))
#plt.plot(lv, gamma_from_Vv,  label=rf'$\gamma$ from $V_v$ (target={Vv_target:.2f})', lw=2)
plt.plot(lv, gamma_from_Cnb, color='g', label=rf'$\gamma$ from $C_{{n_\beta}}$ (target={Cnb_target:.2f})', lw=2)
#plt.scatter([lv_i], [gam_i], s=15, color='g', zorder=5, label=rf'Chosen Point: $l_v={lv_i:.3f}$ m, $\gamma={gam_i:.3f}$')

plt.xlabel(r'$l_v=(x_v-x_{CG})\ \mathrm{[m]}$')
plt.ylabel(r'$\gamma=S_v/S_w$')
plt.title(rf'$\text{{Area ratio vs tail arm: }}V_v={Vv_target:.2f},\ C_{{n_\beta}}={Cnb_target:.2f},\ AR_v={ARv_assumed:.1f}$')
plt.grid(True, ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
