# pip install CoolProp
from CoolProp.CoolProp import PropsSI, PhaseSI
import CoolProp.CoolProp as CD
from CoolProp.CoolProp import PropsSI, PhaseSI
import math

fluid = "N2O"
T = 14+273.15          # K
P_bar = 50.0          # bar
P = P_bar * 1e5       # Pa

# --- 1) Densidade a T=290.15 K e P=50 bar ---
rho_TP = PropsSI("D", "T", T, "P", P, fluid)       # kg/m^3
phase_TP = PhaseSI("T", T, "P", P, fluid)          # string da fase prevista


# --- 2) Propriedades de saturação a T=290.15 K ---
# Pressão de saturação
Psat = PropsSI("P", "T", T, "Q", 0, fluid)         # Pa


# Densidades saturadas: líquido (Q=0) e vapor (Q=1)
rho_L_sat = PropsSI("D", "T", T, "Q", 0, fluid)    # kg/m^3
rho_V_sat = PropsSI("D", "T", T, "Q", 1, fluid)    # kg/m^3

mass = rho_L_sat * 0.01147
print(f'MAss={mass}')

print(f"Psat(T={T:.2f} K) = {Psat/1e5:.6g} bar")
print(f"ρ_L_sat(T={T:.2f} K) = {rho_L_sat:.6g} kg/m³")
print(f"ρ_V_sat(T={T:.2f} K) = {rho_V_sat:.6g} kg/m³")

