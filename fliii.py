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

    # Curvas da asa (para plot)
    def chord_upper(self, y):
        b = self.wingspan
        return self.c_upper * np.sqrt(np.maximum(0.0, 1 - (4 * y**2) / b**2))

    def chord_lower(self, y):
        b = self.wingspan
        return -self.c_lower * np.sqrt(np.maximum(0.0, 1 - (4 * y**2) / b**2))

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
    # CP e AC (desde LE da MAC)
    # ============================

    def chordwise_cp(self, Cm, Cl):
        """x_CP/MAC ≈ 0.25 - (Cm/Cl)  → retorna em metros DESDE o LE da MAC"""
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        """y_CP,h = 2b / (3π)  (válido p/ carga ~elíptica)"""
        return (2 * self.wingspan) / (3 * math.pi)

    @property
    def x_ac(self):
        """AC (25% MAC) em metros DESDE o LE da MAC"""
        return 0.25 * self.mean_aerodynamic_chord

    # ============================
    # REFERÊNCIA GRÁFICA (origem = linha média)
    # ============================

    @property
    def x_le_mac(self):
        """Bordo de ataque da MAC medido a partir da linha média (referência do gráfico)."""
        return -(8.0 / (3.0 * math.pi)) * self.c_lower

    @property
    def y_mac(self):
        """Estação spanwise onde c(y)=MAC (para desenhar a MAC)."""
        ratio = (8.0 / (3.0 * math.pi))  # MAC / (c_upper + c_lower)
        return 0.5 * self.wingspan * math.sqrt(max(0.0, 1.0 - ratio**2))

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
        """CG físico como fração da MAC DESDE o LE da MAC (0.20 ⇒ 20% MAC)."""
        return float(frac_mac) * self.mean_aerodynamic_chord

    def static_margin_cp(self, xcg_phys_from_le, Cm, Cl):
        """(CP − CG)/MAC (ambos desde LE da MAC)  >0 ⇒ CG à frente do CP."""
        xcp_from_le = self.chordwise_cp(Cm, Cl)
        return (xcp_from_le - xcg_phys_from_le) / self.mean_aerodynamic_chord

    # ============================
    # PLOTAGEM (MAC como segmento)
    # ============================

    def plot_wing_with_cg_cp_ac(self, Cm, Cl, cg_phys_frac_mac=0.20, show_half_cp=True):
        """
        Plota a asa e, em y = y_MAC, desenha a MAC como segmento (LE→TE).
        Marca AC (25%), CP e CG físico sobre a MAC. CG geométrico também.
        Tudo na mesma referência (origem = linha média).
        """
        b = self.wingspan
        y = np.linspace(-b/2, b/2, 500)

        x_upper = self.chord_upper(y)
        x_lower = self.chord_lower(y)

        # Dados da MAC para desenho
        y_mac = self.y_mac
        x_le_mac_plot = self.x_le_mac                 # ref. linha média
        L_mac = self.mean_aerodynamic_chord
        x_te_mac_plot = x_le_mac_plot + L_mac

        # Posições ao longo da MAC (ref. linha média)
        x_ac_plot  = self.x_ac_plot()
        x_cp_plot  = self.x_cp_plot(Cm, Cl)
        x_cg_phys_plot = self.x_le_mac + self.x_cg_physical(cg_phys_frac_mac)

        # CG geométrico (na origem spanwise y=0)
        x_cg_geom = self.cg_chordwise
        y_center = 0.0

        plt.figure(figsize=(10, 5))
        plt.plot(y, x_upper, label='Superfície superior')
        plt.plot(y, x_lower, label='Superfície inferior')
        plt.axhline(0, linestyle='--', linewidth=0.8)
        plt.axvline(0, linestyle='--', linewidth=0.8)

        # Desenha a MAC (segmento vertical em y=y_mac)
        plt.plot([y_mac, y_mac], [x_le_mac_plot, x_te_mac_plot], linestyle='-', linewidth=2, label='MAC (segmento)')
        # Marcações na MAC
        plt.scatter(y_mac, x_ac_plot, s=60, label='AC (25% MAC)')
        plt.scatter(y_mac, x_cp_plot, s=70, label='CP (asa completa)')
        if show_half_cp:
            plt.scatter(self.spanwise_cp_half, x_cp_plot, s=55, label='CP (meia-asa)')
        plt.scatter(y_mac, x_cg_phys_plot, s=70, label=f'CG físico ({100*cg_phys_frac_mac:.0f}% MAC)')

        # CG geométrico (área)
        plt.scatter(y_center, x_cg_geom, s=60, label='CG geométrico')

        # Margem (numérica, indep. da origem)
        SM = self.static_margin_cp(self.x_cg_physical(cg_phys_frac_mac), Cm, Cl)
        midx = 0.5 * (x_cg_phys_plot + x_cp_plot)
        plt.annotate(f"SM≈{SM:.2f} (CP−CG)/MAC", xy=(y_mac, midx), xytext=(y_mac + 0.05*b, midx))

        plt.title("Asa semi-elíptica: MAC (segmento), AC, CP, CG físico e CG geométrico")
        plt.xlabel("y (spanwise, m)")
        plt.ylabel("x (chordwise, m)")
        plt.legend()
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    

glider = Glider(
    mass=0.8,
    wingspan=0.8,
    c_upper=0.028,
    c_lower=0.062,
    oswald_efficiency=0.75,
    cd0=0.04
)

Glider.plot_wing_with_cg_cp_ac(glider, -0.05, 1.7, 0.2)