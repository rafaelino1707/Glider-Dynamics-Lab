import math
import matplotlib.pyplot as plt
import numpy as np

class Glider:
    def __init__(self, mass, wingspan, c_upper, c_lower, oswald_efficiency, cd0):
        self.mass = mass
        self.wingspan = wingspan
        self.c_upper = c_upper
        self.c_lower = c_lower
        self.oswald_efficiency = oswald_efficiency
        self.cd0 = cd0

    # ============================
    # GEOMETRIA
    # ============================

    @property
    def total_chord(self):
        """C = c_upper + c_lower"""
        return self.c_upper + self.c_lower

    @property
    def wing_area(self):
        """S = (π * b / 4) * (c_upper + c_lower)"""
        return (math.pi * self.wingspan / 4) * self.total_chord

    @property
    def mean_aerodynamic_chord(self):
        """MAC = (8 / (3π)) * (c_upper + c_lower)"""
        return (8 / (3 * math.pi)) * self.total_chord

    @property
    def aspect_ratio(self):
        """AR = 4b / [π (c_upper + c_lower)]"""
        return 4 * self.wingspan / (math.pi * self.total_chord)

    # ============================
    # CG GEOMÉTRICO
    # ============================

    @property
    def cg_chordwise(self):
        """x̄ (geométrico) = (4 / (3π)) * (c_upper - c_lower)  (referência: linha média)"""
        return (4 / (3 * math.pi)) * (self.c_upper - self.c_lower)

    @property
    def cg_spanwise_half(self):
        """ȳ_h = (2b) / (3π)"""
        return (2 * self.wingspan) / (3 * math.pi)

    # ============================
    # AERODINÂMICA
    # ============================

    @property
    def induced_drag_factor(self):
        """K = 1 / (π * AR * e)"""
        return 1 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    def minimum_sink_velocity(self, rho=1.225):
        """V_min_sink ≈ sqrt(2/rho) * (K/CD0)^(1/4) * sqrt(W/S)"""
        W = self.mass * 9.81
        term_1 = math.sqrt(2 / rho)
        term_2 = (self.induced_drag_factor / self.cd0) ** 0.25
        term_3 = math.sqrt(W / self.wing_area)
        return term_1 * term_2 * term_3

    def vertical_sink_rate(self, V=None, rho=1.225):
        """Calcula taxa de sink vertical"""
        W = self.mass * 9.81
        S = self.wing_area
        if V is None:
            V = self.minimum_sink_velocity(rho)
        CL = W / (0.5 * rho * V**2 * S)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2
        D = 0.5 * rho * V**2 * S * CD
        sink_rate = V * (D / W)
        glide_ratio = W / D if D != 0 else float("inf")
        return {
            "V_used": V,
            "CL": CL,
            "CD": CD,
            "D_N": D,
            "sink_rate_m_s": sink_rate,
            "glide_ratio": glide_ratio
        }

    def stall_velocity(self, CL_max, rho=1.225):
        """V_stall = sqrt((2W)/(ρ S CL_max))"""
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2 * W) / (rho * S * CL_max))

    # ============================
    # CP e AC (valores em fração da MAC e em metros DESDE o LE da MAC)
    # ============================

    def chordwise_cp(self, Cm, Cl):
        """
        x_CP/MAC ≈ 0.25 - (Cm/Cl)
        Retorna x_CP em metros, MEDIDO DESDE o LE da MAC.
        """
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        """y_CP,h = 2b / (3π)  (válido p/ carga ~elíptica)"""
        return (2 * self.wingspan) / (3 * math.pi)

    @property
    def x_ac(self):
        """AC (25% MAC) em metros, MEDIDO DESDE o LE da MAC"""
        return 0.25 * self.mean_aerodynamic_chord

    # ============================
    # REFERÊNCIA GRÁFICA (origem = linha média)
    # ============================

    @property
    def x_le_mac(self):
        """Bordo de ataque da MAC medido a partir da linha média (referência do gráfico)."""
        # Para a planforma de semi-elipses combinadas: x_LE,MAC = -(8/(3π)) * c_lower
        return -(8.0 / (3.0 * math.pi)) * self.c_lower

    def x_ac_plot(self):
        """AC convertido para a referência do gráfico (somando offset do LE da MAC)."""
        return self.x_le_mac + self.x_ac

    def x_cp_plot(self, Cm, Cl):
        """CP convertido para a referência do gráfico (somando offset do LE da MAC)."""
        return self.x_le_mac + self.chordwise_cp(Cm, Cl)

    # ============================
    # CG FÍSICO (AJUSTÁVEL)
    # ============================

    def x_cg_physical(self, frac_mac):
        """
        CG físico como fração da MAC MEDIDA DESDE o LE da MAC.
        Ex.: frac_mac=0.20 => x = 0.20*MAC (desde o LE da MAC)
        """
        return float(frac_mac) * self.mean_aerodynamic_chord

    def static_margin_cp(self, xcg_phys_from_le, Cm, Cl):
        """
        'Margem' usando CP: (x_CP - x_CG)/MAC, ambos medidos DESDE o LE da MAC.
        >0 => CG à frente do CP (tendência estável).
        """
        xcp_from_le = self.chordwise_cp(Cm, Cl)
        return (xcp_from_le - xcg_phys_from_le) / self.mean_aerodynamic_chord

    # ============================
    # PLOTAGEM (referência coerente)
    # ============================

    def plot_wing_with_cg_cp_ac(self, Cm, Cl, cg_phys_frac_mac=0.20, show_half_cp=True):
        """
        Plota asa + CG geométrico + CP + AC + CG físico (fração da MAC),
        TODOS no mesmo referencial do gráfico (origem = linha média).
        """
        b = self.wingspan
        y = np.linspace(-b/2, b/2, 400)

        x_upper = self.c_upper * np.sqrt(1 - (4 * y**2) / b**2)
        x_lower = -self.c_lower * np.sqrt(1 - (4 * y**2) / b**2)

        # Conversões corretas de referência
        x_cg_geom = self.cg_chordwise                            # já na linha média
        x_cp_plot = self.x_cp_plot(Cm, Cl)                       # CP na ref. do gráfico
        x_ac_plot = self.x_ac_plot()                             # AC na ref. do gráfico
        x_cg_phys_plot = self.x_le_mac + self.x_cg_physical(cg_phys_frac_mac)  # CG físico na ref. do gráfico
        y_center  = 0.0
        y_cp_half = self.spanwise_cp_half

        plt.figure(figsize=(9.5, 4.8))
        plt.plot(y, x_upper, label='Superfície superior')
        plt.plot(y, x_lower, label='Superfície inferior')
        plt.axhline(0, linestyle='--', linewidth=0.8)
        plt.axvline(0, linestyle='--', linewidth=0.8)

        plt.scatter(y_center, x_cg_geom, s=60, label='CG geométrico')
        plt.scatter(y_center, x_cp_plot, s=70, label='CP (asa completa)')
        if show_half_cp:
            plt.scatter(y_cp_half, x_cp_plot, s=55, label='CP (meia-asa)')

        # AC como linha horizontal (mesma referência)
        plt.axhline(x_ac_plot, linestyle=':', linewidth=1.2, label='AC (25% MAC)')

        # CG físico (fração da MAC) – mesma referência
        plt.scatter(y_center, x_cg_phys_plot, s=70, label=f'CG físico ({100*cg_phys_frac_mac:.0f}% MAC)')

        # Margem numérica (independe da origem) usa valores "desde LE da MAC"
        SM = self.static_margin_cp(self.x_cg_physical(cg_phys_frac_mac), Cm, Cl)
        plt.annotate(f"SM≈{SM:.2f} (CP−CG)/MAC", xy=(0, (x_cg_phys_plot + x_cp_plot)/2),
                     xytext=(0.05*b, (x_cg_phys_plot + x_cp_plot)/2))

        plt.title("Asa semi-elíptica: CG geom., CG físico, CP e AC (referência coerente)")
        plt.xlabel("y (spanwise, m)")
        plt.ylabel("x (chordwise, m)")
        plt.legend()
        plt.axis('equal')
        plt.grid(True)
        plt.show()


# ============================
# EXEMPLO DE USO
# ============================

glider = Glider(
    mass=0.8,
    wingspan=0.4,
    c_upper=0.014,
    c_lower=0.026,
    oswald_efficiency=0.75,
    cd0=0.04
)

print("=== GEOMETRIA ===")
print(f"Wing Area S = {glider.wing_area:.5f} m²")
print(f"MAC = {glider.mean_aerodynamic_chord:.5f} m")
print(f"AR = {glider.aspect_ratio:.3f}")
print(f"CG geométrico (x) = {glider.cg_chordwise:.5f} m")
print(f"LE da MAC (x ref. linha média) = {glider.x_le_mac:.5f} m")
print(f"AC (25% MAC) (x, ref. linha média) = {glider.x_ac_plot():.5f} m")

# Dados aerodinâmicos (ajusta aos teus)
Cm_section = -0.05
Cl_section = 1.7

# CP desde LE da MAC e convertido para a ref. do gráfico
x_cp_from_le = glider.chordwise_cp(Cm_section, Cl_section)
x_cp_plot = glider.x_cp_plot(Cm_section, Cl_section)

# CG físico como % da MAC (ex.: 20% MAC)
cg_phys_frac = 0.20
x_cg_phys_from_le = glider.x_cg_physical(cg_phys_frac)
SM = glider.static_margin_cp(x_cg_phys_from_le, Cm_section, Cl_section)

print("\n=== POSIÇÕES AERODINÂMICAS ===")
print(f"x_CP (desde LE da MAC) = {x_cp_from_le:.5f} m  ({0.25 - Cm_section/Cl_section:.3f} × MAC)")
print(f"x_CP (ref. linha média p/ plot) = {x_cp_plot:.5f} m")
print(f"x_CG_físico (desde LE da MAC) = {x_cg_phys_from_le:.5f} m  ({cg_phys_frac:.2f} × MAC)")
print(f"Margem (CP − CG)/MAC ≈ {SM:.3f}  (>0 desejável)")

# Plot geral (tudo na mesma referência)
glider.plot_wing_with_cg_cp_ac(Cm=Cm_section, Cl=Cl_section, cg_phys_frac_mac=cg_phys_frac)
