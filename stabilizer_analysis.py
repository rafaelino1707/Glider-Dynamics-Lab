import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# -----------------------------
# Given values
# -----------------------------
CLw = 1.29            # Lift coefficient of the wing
CLh = -0.4         # Lift coefficient of the horizontal stabilizer
Sw = 0.02356          # Wing area [m^2]
xac_w = 0.0106        # Aerodynamic center of the wing [m]
xCG = 0.09            # Center of gravity [m]
MAC = 0.0424          # Mean aerodynamic chord [m]
Cm_ac = -0.2494       # Moment coefficient at aerodynamic center

# -----------------------------
# Range for beta (Sh = beta * Sw)
# -----------------------------
beta = np.linspace(0.3, 1, 200)   # From 5% to 50% of the main wing area
Sh = beta * Sw

# -----------------------------
# Compute (x_h - x_CG)
# Equation:
# (x_h - x_CG) = [CLw * Sw * (xac_w - xCG) + Cm_ac * Sw * MAC] / (CLh * Sh)
# -----------------------------
numerator = CLw * Sw * (xac_w - xCG) + Cm_ac * Sw * MAC
xh_minus_xcg = numerator / (CLh * Sh)

data = {
    'beta': beta,
    'S_h [m^2]': Sh,
    '(x_h - x_CG) [m]': xh_minus_xcg
}

df = pd.DataFrame(data)
df.to_csv('xh_vs_Sh.csv', index=False)  # <-- ficheiro CSV guardado no diretório atual
print("File 'xh_vs_Sh.csv' saved successfully.")

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(beta, xh_minus_xcg, color='royalblue', linewidth=2)
plt.title(r'Distance $(x_h - x_{CG})$ as a Function of Stabilizer Area Ratio $\beta$', fontsize=13)
plt.xlabel(r'Area Ratio $\beta = S_h / S_w$', fontsize=12)
plt.ylabel(r'$(x_h - x_{CG})$ [m]', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()


# -----------------------------
# Given values (entrada)
# -----------------------------
CLw = 1.29            # Lift coefficient of the wing
CLh = -0.4            # Lift coefficient assumed for the horizontal stabilizer (para comparar)
Sw = 0.02356          # Wing area [m^2]
xac_w = 0.0106        # Aerodynamic center of the wing [m]
xCG = 0.09            # Center of gravity [m]
MAC = 0.0424          # Mean aerodynamic chord [m]
Cm_ac = -0.2494       # Moment coefficient at aerodynamic center

# geometria / braço estabilizador (usar o valor do projecto; aqui 0.50 m)
l_h = 0.6       # (x_h - x_CG) [m]

# -----------------------------
# Range for beta (Sh = beta * Sw)
# -----------------------------
# beta is the ratio S_h / S_w. Choose a range that makes sense for small models.
beta = np.linspace(0.2, 1.0, 300)   # cobre desde stabilizer muito pequeno até S_h = S_w
Sh = beta * Sw

# -----------------------------
# Numerator of trim equation (independent of Sh)
# numerator = CLw*Sw*(xac_w - xCG) + Cm_ac*Sw*MAC
# -----------------------------
numerator = CLw * Sw * (xac_w - xCG) + Cm_ac * Sw * MAC

# -----------------------------
# Compute (x_h - x_CG) from the equation:
# (x_h - x_CG) = numerator / (CLh * Sh)
# and the CL_h required to trim:
# CLh_required = numerator / (Sh * l_h)
# Also compute V_H = (S_h * l_h) / (S_w * MAC)
# -----------------------------
# protect against division by zero
Sh_safe = np.where(Sh <= 0, np.nan, Sh)

xh_minus_xcg = numerator / (CLh * Sh_safe)           # resultado usando CLh (para plotar a curva que usavas)
CLh_required = numerator / (Sh_safe * l_h)          # que CLh seria necessário para cada Sh dado o l_h
V_H = (Sh_safe * l_h) / (Sw * MAC)                   # coeficiente de volume

# -----------------------------
# Save results to CSV
# -----------------------------
df = pd.DataFrame({
    'beta': beta,
    'S_h [m^2]': Sh_safe,
    '(x_h - x_CG) [m] (with CLh={:.3f})'.format(CLh): xh_minus_xcg,
    'CLh_required_for_trim': CLh_required,
    'V_H': V_H
})
df.to_csv('xh_vs_Sh_and_CLh.csv', index=False)
print("File 'xh_vs_Sh_and_CLh.csv' saved successfully.")

# -----------------------------
# Plots: (1) xh-xCG vs beta  (2) CLh_required vs beta
# -----------------------------
fig, ax = plt.subplots(2, 1, figsize=(8,10), sharex=True)

# plot 1: xh - xCG (using assumed CLh)
ax[0].plot(beta, xh_minus_xcg, linewidth=2)
ax[0].axhline(l_h, color='0.15', linestyle='--', label=r'$(x_h - x_{CG})$'+ f' = {l_h:.2f} m')
ax[0].set_ylabel(r'$(x_h - x_{CG})$ [m]')
ax[0].set_title(r'Distance $(x_h - x_{CG})$ as a Function of Stabilizer Area Ratio $\beta$')
ax[0].grid(True, linestyle='--', alpha=0.5)
ax[0].legend()

# plot 2: CLh_required
ax[1].plot(beta, CLh_required, linewidth=2, label='CLh_required (trim)')
ax[1].axhline(CLh, color='0.15', linestyle='--', label=f'CLh_assumed = {CLh:.2f}')
# mark some realistic CLh bounds
ax[1].set_xlabel(r'Area Ratio $\beta = S_h / S_w$')
ax[1].set_ylabel(r'$C_{L_h,\ required}$')
ax[1].set_title(r'$C_{L_h}$ Required for Trim as Function of Stabilizer Area Ratio $\beta$)')
ax[1].set_ylim(np.nanmin(CLh_required)*1.1, np.nanmax(CLh_required)*0.9)  # ajusta visualização
ax[1].grid(True, linestyle='--', alpha=0.5)
ax[1].legend()

plt.tight_layout()
plt.show()
