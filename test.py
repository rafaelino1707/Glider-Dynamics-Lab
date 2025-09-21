# lift_estimates.py
import math
import aerosandbox as asb
import numpy as np
from aerosandbox.aerodynamics.aero_2D.xfoil import XFoil

# constantes
rho = 1.225
mu = 1.81e-5
g = 9.81

# parâmetros do teu planeador (exemplo; adapta)
m = 1.0            # kg
W = m * g
S = 0.25           # m^2 (altera para o teu valor)
b = 1.5            # m (envergadura)
c = S / b          # corda média
AR = b**2 / S

# alpha em graus para testar
alpha_deg = 5.0
alpha_rad = math.radians(alpha_deg)

# 1) thin airfoil (2D) -> CL per rad ~= 2*pi
CL2d_thin = 2 * math.pi * alpha_rad
print("Thin-airfoil CL (2D) at", alpha_deg, "deg:", CL2d_thin)

# 2) lifting-line approx -> wing lift slope
e = 0.9  # eficiência de Oswald aproximada
a_wing = (2 * math.pi * AR) / (2 + AR / e)  # per rad
CL_wing_from_thin = a_wing * alpha_rad
print("Lifting-line approx CL (wing) at", alpha_deg, "deg:", CL_wing_from_thin)

# 3) CL required (para verificar velocidade)
V_from_CL = math.sqrt(2 * W / (rho * S * CL_wing_from_thin))
Re_est = rho * V_from_CL * c / mu
Mach_est = V_from_CL / 340.0
print(f"Estimated V (m/s) = {V_from_CL:.2f}, Re = {Re_est:.2e}, Mach = {Mach_est:.3f}")

# 4) AeroSandbox: NeuralFoil quick estimate (se disponível)
af = asb.Airfoil("s3021-il")  # nome compatível com AeroSandbox
try:
    neural = af.get_aero_from_neuralfoil(alpha=alpha_deg, Re=Re_est)
    print("NeuralFoil estimate:", neural)
except Exception as e:
    print("NeuralFoil failed:", e)

# 5) AeroSandbox XFoil viscous run (requires xfoil binary accessible or download)
xf = XFoil()
xf.Re = Re_est
xf.max_iter = 200
try:
    res = xf.alpha(af, alpha=alpha_deg)   # roda XFoil viscous
    print("XFoil viscous result:", res)
except Exception as e:
    print("XFoil viscous failed (check xfoil available):", e)

# Nota: adapta S, b, m conforme o teu projeto
