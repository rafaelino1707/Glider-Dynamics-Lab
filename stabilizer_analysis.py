import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ==========================
# Dados do teu projeto
# ==========================
Sw   = 0.03770   # área da asa [m^2]
MAC  = 0.0679    # mean aerodynamic chord [m]

# Para a curva de TRIM (opcional, para cruzar com a de VH)
CLw    = 1.29        # CL da asa (condição de referência)
Cm_ac  = -0.2494     # Cm no AC da asa
xac_w  = 0.0170      # AC da asa [m] (desde a referência usada para x)
xCG    = 0.09        # centro de gravidade [m]
CLh_as = -0.40       # CL assumido para o estabilizador (sinal tail-down)

# ==========================
# Escolha do alvo de volume horizontal
# ==========================
VH_target = 0.50     # típico: 0.4–0.7

# ==========================
# Varre l_h = (x_h - x_CG)
# ==========================
l_h = np.linspace(0.15, 0.50, 300)  # braço de cauda [m]

# Curva 1: beta ditada por VH
# VH = (S_h * l_h)/(S_w * MAC)  => beta = S_h/S_w = VH * MAC / l_h
beta_VH = VH_target * MAC / l_h

# Curva 2: beta ditada por TRIM com CLh assumido
# Equilíbrio de momentos em torno do CG (forma linearizada usada no teu código):
# l_h = [CLw*Sw*(xac_w - xCG) + Cm_ac*Sw*MAC] / (CLh * S_h)
# => beta_trim = [CLw*(xac_w - xCG) + Cm_ac*MAC] / (CLh * l_h)
num = CLw*(xac_w - xCG) + Cm_ac*MAC
beta_trim = num / (CLh_as * l_h)

# Restrições simples
beta_VH[beta_VH <= 0] = np.nan
beta_trim[beta_trim <= 0] = np.nan

# ==========================
# Guardar CSV
# ==========================
df = pd.DataFrame({
    'l_h [m]': l_h,
    'beta_from_VH': beta_VH,
    'beta_from_trim(CLh={:.2f})'.format(CLh_as): beta_trim,
    'VH_target': np.full_like(l_h, VH_target),
})
df.to_csv('beta_vs_lh_VH_and_trim.csv', index=False)
print("CSV guardado: beta_vs_lh_VH_and_trim.csv")

# ==========================
# Plot
# ==========================
plt.figure(figsize=(8,5))
plt.plot(l_h, beta_VH,  label=r'$\beta$ from $V_H$ (target={:.2f})'.format(VH_target), lw=2)
plt.plot(l_h, beta_trim,
         label=rf'$\beta$ from trim (assume $C_{{L_h}}$={CLh_as:.2f})',
         lw=2)

plt.xlabel(r'$l_h = (x_h - x_{CG})$ [m]')
plt.ylabel(r'$\beta = S_h/S_w$')
plt.title(r'$\beta$ vs $l_h$: alvo de $V_H$ e requisito de trim')
plt.grid(True, ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()

# ==========================
# Nota prática:
# O ponto de interseção das curvas (se existir) satisfaz simultaneamente:
#   - o V_H escolhido
#   - o trim com CLh assumido
# Caso não haja interseção, ajusta VH_target, CLh_as ou a gama de l_h.
# ==========================
