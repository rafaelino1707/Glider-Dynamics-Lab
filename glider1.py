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
    # CENTRO DE GRAVIDADE
    # ============================

    @property
    def cg_chordwise(self):
        """x̄ = (4 / (3π)) * (c_upper - c_lower)"""
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
    # CENTRO DE PRESSÃO
    # ============================

    def chordwise_cp(self, Cm, Cl):
        """x_CP/MAC ≈ 0.25 - (Cm/Cl)"""
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        """y_CP,h = 2b / (3π)"""
        return (2 * self.wingspan) / (3 * math.pi)

    # ============================
    # PLOTAGEM
    # ============================

    def plot_wing_with_cg(self):
        """
        Plota o contorno da asa (semi-elíptica combinada)
        e marca o centro de gravidade geométrico (CG)
        """
        b = self.wingspan
        y = np.linspace(-b/2, b/2, 200)

        x_upper = self.c_upper * np.sqrt(1 - (4 * y**2) / b**2)
        x_lower = -self.c_lower * np.sqrt(1 - (4 * y**2) / b**2)

        # CG global
        x_cg = self.cg_chordwise
        y_cg = 0  # centrado no meio da asa

        plt.figure(figsize=(8, 4))
        plt.plot(y, x_upper, 'b', label='Superfície superior')
        plt.plot(y, x_lower, 'r', label='Superfície inferior')
        plt.axhline(0, color='k', linestyle='--', linewidth=0.8)
        plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        plt.scatter(y_cg, x_cg, color='orange', s=20, zorder=5, label='CG')

        plt.title("Asa semi-elíptica combinada com CG")
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
    c_upper=0.05,
    c_lower=0.1,
    oswald_efficiency=0.75,
    cd0=0.04
)

print("=== GEOMETRIA ===")
print(f"Wing Area S = {glider.wing_area:.5f} m²")
print(f"Mean Aerodynamic Chord (MAC) = {glider.mean_aerodynamic_chord:.5f} m")
print(f"Aspect Ratio (AR) = {glider.aspect_ratio:.3f}")
print(f"CG chordwise = {glider.cg_chordwise:.5f} m")
print(f"CG spanwise (half-wing) = {glider.cg_spanwise_half:.5f} m")

# Plota a geometria com o CG
glider.plot_wing_with_cg()
