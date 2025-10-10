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
        """x̄ (geométrico) = (4 / (3π)) * (c_upper - c_lower)"""
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
    # CP e AC
    # ============================

    def chordwise_cp(self, Cm, Cl):
        """
        x_CP/MAC ≈ 0.25 - (Cm/Cl)
        Retorna x_CP em metros (mesma referência usada no gráfico).
        """
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        """y_CP,h = 2b / (3π)  (válido p/ carga ~elíptica)"""
        return (2 * self.wingspan) / (3 * math.pi)

    @property
    def x_ac(self):
        """Aerodynamic Center (25% MAC) em metros no eixo x do gráfico"""
        return 0.25 * self.mean_aerodynamic_chord

    # ============================
    # CG FÍSICO (AJUSTÁVEL)
    # ============================

    def x_cg_physical(self, frac_mac):
        """
        CG físico como fração da MAC medida do bordo de ataque da MAC.
        Ex.: frac_mac=0.20 => x = 0.20*MAC
        """
        return float(frac_mac) * self.mean_aerodynamic_chord

    def static_margin_cp(self, xcg_phys, Cm, Cl):
        """
        'Margem estática' usando CP como proxy: (x_CP - x_CG)/MAC.
        >0 => CG à frente do CP (tendência estável).
        """
        xcp = self.chordwise_cp(Cm, Cl)
        return (xcp - xcg_phys) / self.mean_aerodynamic_chord

    # ============================
    # PLOTAGEM
    # ============================

    def plot_wing_with_cg_cp_ac(self, Cm, Cl, cg_phys_frac_mac=0.20, show_half_cp=True):
        """
        Plota asa + CG geométrico + CP + AC + CG físico (fração da MAC).
        """
        b = self.wingspan
        y = np.linspace(-b/2, b/2, 400)

        x_upper = self.c_upper * np.sqrt(1 - (4 * y**2) / b**2)
        x_lower = -self.c_lower * np.sqrt(1 - (4 * y**2) / b**2)

        # Pontos
        x_cg_geom = self.cg_chordwise
        x_cp = self.chordwise_cp(Cm, Cl)
        x_ac = self.x_ac
        x_cg_phys = self.x_cg_physical(cg_phys_frac_mac)
        y_center = 0.0
        y_cp_half = self.spanwise_cp_half

        # Gráfico
        plt.figure(figsize=(9.5, 4.8))
        plt.plot(y, x_upper, label='Superfície superior')
        plt.plot(y, x_lower, label='Superfície inferior')
        plt.axhline(0, linestyle='--', linewidth=0.8)
        plt.axvline(0, linestyle='--', linewidth=0.8)

        # Marcadores
        plt.scatter(y_center, x_cg_geom, s=60, label='CG geométrico')
        plt.scatter(y_center, x_cp, s=70, label='CP (asa completa)')
        if show_half_cp:
            plt.scatter(y_cp_half, x_cp, s=55, label='CP (meia-asa)')

        # AC (linha vertical pontilhada no 25% MAC)
        plt.axhline(x_ac, linestyle=':', linewidth=1.2, label='AC (25% MAC)')

        # CG físico (ajustável)
        plt.scatter(y_center, x_cg_phys, s=70, label=f'CG físico ({100*cg_phys_frac_mac:.0f}% MAC)')

        # Anotação da margem (CP vs CG físico)
        SM = self.static_margin_cp(x_cg_phys, Cm, Cl)
        plt.annotate(f"SM≈{SM:.2f} (CP−CG)/MAC", xy=(0, (x_cg_phys + x_cp)/2),
                     xytext=(0.05*b, (x_cg_phys + x_cp)/2))

        plt.title("Asa semi-elíptica: CG geom., CG físico, CP e AC")
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
    wingspan=1.0,
    c_upper=0.035,
    c_lower=0.065,
    oswald_efficiency=0.75,
    cd0=0.04
)

print("=== GEOMETRIA ===")
print(f"Wing Area S = {glider.wing_area:.5f} m²")
print(f"MAC = {glider.mean_aerodynamic_chord:.5f} m")
print(f"AR = {glider.aspect_ratio:.3f}")
print(f"CG geométrico (x) = {glider.cg_chordwise:.5f} m")
print(f"AC (25% MAC) (x) = {glider.x_ac:.5f} m")

# Dados aerodinâmicos (ajusta aos teus)
Cm_section = -0.03
Cl_section = 2

# CG físico como % da MAC (ex.: 20% MAC)
cg_phys_frac = 0.20
x_cg_phys = glider.x_cg_physical(cg_phys_frac)
SM = glider.static_margin_cp(x_cg_phys, Cm_section, Cl_section)
print("\n=== POSIÇÕES AERODINÂMICAS ===")
print(f"x_CP = {glider.chordwise_cp(Cm_section, Cl_section):.5f} m")
print(f"x_CG_físico = {x_cg_phys:.5f} m  ({cg_phys_frac:.2f} × MAC)")
print(f"Margem (CP − CG)/MAC ≈ {SM:.3f}  (>0 desejável)")

# Plot geral
glider.plot_wing_with_cg_cp_ac(Cm=Cm_section, Cl=Cl_section, cg_phys_frac_mac=cg_phys_frac)
