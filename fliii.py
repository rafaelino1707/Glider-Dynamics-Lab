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
    # CG GEOMÉTRICO (referência: linha média, x>0 para cima/upper)
    # ============================

    @property
    def cg_chordwise(self):
        """x̄_geom = (4 / (3π)) * (c_upper - c_lower)"""
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
    # CP e AC (DISTÂNCIA DESDE o LE da MAC)
    # ============================

    def chordwise_cp(self, Cm, Cl):
        """
        x_CP/MAC ≈ 0.25 - (Cm/Cl)
        Retorna a distância POSITIVA desde o LE da MAC até ao CP (em metros).
        """
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        """y_CP,h = 2b / (3π)  (válido p/ carga ~elíptica)"""
        return (2 * self.wingspan) / (3 * math.pi)

    @property
    def x_ac(self):
        """Distância desde o LE da MAC até ao AC (25% MAC)"""
        return 0.25 * self.mean_aerodynamic_chord

    # ============================
    # REFERÊNCIA DO GRÁFICO (origem = linha média, upper = +x)
    # ============================

    @property
    def x_le_mac(self):
        """
        Bordo de ataque da MAC medido a partir da linha média.
        Convenção: upper = LE ⇒ x_le_mac > 0.
        """
        return (8.0 / (3.0 * math.pi)) * self.c_upper

    @property
    def y_mac(self):
        """Estação spanwise onde c(y)=MAC (para desenhar a MAC)."""
        ratio = (8.0 / (3.0 * math.pi))  # MAC / (c_upper + c_lower)
        return 0.5 * self.wingspan * math.sqrt(max(0.0, 1.0 - ratio**2))

    def x_ac_plot(self):
        """AC na referência do gráfico (a partir do LE para TRÁS)."""
        return self.x_le_mac - self.x_ac

    def x_cp_plot(self, Cm, Cl):
        """CP na referência do gráfico (a partir do LE para TRÁS)."""
        return self.x_le_mac - self.chordwise_cp(Cm, Cl)

    # ============================
    # CG FÍSICO (AJUSTÁVEL)
    # ============================

    def x_cg_physical(self, frac_mac):
        """Distância desde o LE da MAC até ao CG físico (0.20 ⇒ 20% MAC)."""
        return float(frac_mac) * self.mean_aerodynamic_chord

    def static_margin_cp(self, xcg_phys_from_le, Cm, Cl):
        """(CP − CG)/MAC (ambos desde LE da MAC)  >0 ⇒ CG à frente do CP."""
        xcp_from_le = self.chordwise_cp(Cm, Cl)
        return (xcp_from_le - xcg_phys_from_le) / self.mean_aerodynamic_chord

    # ============================
    # PLOTAGEM (MAC como segmento, percentagens a partir do LE)
    # ============================

    def plot_wing_with_cg_cp_ac(self, Cm, Cl, cg_phys_frac_mac=0.20, show_half_cp=True):
        """
        Plota a asa; desenha a MAC como segmento (em y=y_MAC, só referência),
        e coloca AC, CP (asa completa) e CG físico no centro (y=0).
        O CG geométrico (da área) também aparece em y=0.
        Tudo na referência do gráfico (origem = linha média, upper = +x).
        """
        b = self.wingspan
        y = np.linspace(-b/2, b/2, 500)

        # Contornos
        x_upper = self.chord_upper(y)
        x_lower = self.chord_lower(y)

        # MAC (segmento em y = y_mac, só como "régua" de percentagens)
        y_mac = self.y_mac
        x_le_mac_plot = self.x_le_mac
        L_mac = self.mean_aerodynamic_chord
        x_te_mac_plot = x_le_mac_plot - L_mac

        # Frações na MAC
        f_ac = 0.25
        f_cp = 0.25 - (Cm / Cl)
        f_cg = float(cg_phys_frac_mac)

        # Converter frações para x (LE → para TRÁS)
        x_ac = x_le_mac_plot - f_ac * L_mac
        x_cp = x_le_mac_plot - f_cp * L_mac
        x_cg_phys = x_le_mac_plot - f_cg * L_mac

        # CG geométrico (área)
        x_cg_geom = self.cg_chordwise

        # --- Plot ---
        plt.figure(figsize=(10, 5))
        plt.plot(y, x_upper, label='Superfície superior')
        plt.plot(y, x_lower, label='Superfície inferior')
        plt.axhline(0, linestyle='--', linewidth=0.8)
        plt.axvline(0, linestyle='--', linewidth=0.8)

        # MAC como segmento (em y=y_mac, apenas referência)
        plt.plot([y_mac, y_mac], [x_le_mac_plot, x_te_mac_plot],
                color='g', linewidth=2, label='MAC (segmento, referência)')

        # Marcadores PRINCIPAIS no centro (y=0)
        plt.scatter(0.0, x_ac,      s=20, label='AC (25% MAC)')
        plt.scatter(0.0, x_cp,      s=20, label=f'CP asa ({f_cp*100:.1f}% MAC)')
        plt.scatter(0.0, x_cg_phys, s=20, label=f'CG físico ({f_cg*100:.0f}% MAC)')
        plt.scatter(0.0, x_cg_geom, s=20, label='CG geométrico')

        # Opcional: CP da meia-asa na sua estação spanwise
        if show_half_cp:
            plt.scatter(self.spanwise_cp_half, x_cp, s=20, label='CP (meia-asa)')

        # Margem (numérica, indep. da origem)
        SM = self.static_margin_cp(self.x_cg_physical(cg_phys_frac_mac), Cm, Cl)
        midx = 0.5 * (x_cg_phys + x_cp)
        plt.annotate(f"SM≈{SM:.2f} (CP−CG)/MAC", xy=(0.0, midx), xytext=(0.05*b, midx))

        plt.title("Asa semi-elíptica: MAC (referência), AC, CP e CG em y=0")
        plt.xlabel("y (spanwise, m)")
        plt.ylabel("x (chordwise, m)")
        plt.legend()
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def plot_half_wing_right(self, f_cp=None, Cm=None, Cl=None, cg_phys_frac_mac=None):
        """
        Plota apenas a meia-asa direita (0..b/2) com:
        • CG geométrico da meia-asa (fórmulas integradas)
        • CP da meia-asa no mesmo y (y_CP,h = 2b/(3π))
            - Se f_cp (fração da MAC) for dado, usa-o
            - Caso contrário, se Cm e Cl forem dados, usa f_cp = 0.25 - Cm/Cl
        Opcional: marcar também um CG físico por fração da MAC (cg_phys_frac_mac).
        Convenção: LE no lado upper (x positivo); percentagens contam-se do LE para trás.
        """
        b = self.wingspan
        y = np.linspace(0.0, b/2.0, 400)

        # contornos meia-asa
        x_upper = self.chord_upper(y)
        x_lower = self.chord_lower(y)

        # CG geométrico da meia-asa (mesmo x do total; y integrado da meia-asa)
        x_cg_geom = self.cg_chordwise
        y_cg_half = self.cg_spanwise_half

        # Fração do CP ao longo da MAC
        if f_cp is None:
            if Cm is None or Cl is None:
                raise ValueError("Indica f_cp (fração da MAC) ou então Cm e Cl para calcular o CP.")
            f_cp = 0.25 - (Cm / Cl)

        # Converter frações (LE -> para trás)
        L_mac = self.mean_aerodynamic_chord
        x_le_mac_plot = self.x_le_mac
        x_cp = x_le_mac_plot - float(f_cp) * L_mac
        y_cp_half = y_cg_half  # por simetria / carga ~elíptica

        # (Opcional) CG físico por percentagem da MAC (se quiseres ver no mesmo gráfico)
        x_cg_phys = None
        if cg_phys_frac_mac is not None:
            x_cg_phys = x_le_mac_plot - float(cg_phys_frac_mac) * L_mac

        # --- Plot meia-asa direita ---
        plt.figure(figsize=(8.5, 4.5))
        plt.plot(y, x_upper, label='Upper Surface')
        plt.plot(y, x_lower, label='Lower Surface')
        plt.axhline(0, linestyle='--', linewidth=0.8)
        plt.axvline(0, linestyle='--', linewidth=0.8)

        # CG geométrico da meia-asa
        plt.scatter(y_cg_half, x_cg_geom, s=30, label='Geometric CG (Half-Wing)')

        # CP (meia-asa)
        plt.scatter(y_cp_half, x_cp, s=30, label=f'CP (Half-Wing, {f_cp*100:.1f}% MAC)')

        # (Opcional) CG físico
        if x_cg_phys is not None:
            plt.scatter(y_cg_half, x_cg_phys, s=30, label=f'Physical CG ({cg_phys_frac_mac*100:.0f}% MAC)')

        # MAC desenhada como “régua” (em y=y_MAC, só referência)
        y_mac = self.y_mac
        x_te_mac_plot = x_le_mac_plot - L_mac
        plt.plot([y_mac, y_mac], [x_le_mac_plot, x_te_mac_plot],linestyle='--', color='g', linewidth=2, label='MAC (referência)')

        plt.title("Half Right Wing: Geometric CG and CP")
        plt.xlabel("y (Spanwise, m)")
        plt.ylabel("x (Chordwise, m)")
        plt.legend()
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()



# ============================
# EXEMPLO DE USO
# ============================

if __name__ == "__main__":
    glider = Glider(
        mass=0.8,
        wingspan=0.8,
        c_upper=0.028,   # LE do lado upper (x>0)
        c_lower=0.062,
        oswald_efficiency=0.75,
        cd0=0.04
    )

    # Exemplo: Cm/Cl do teu caso
    Cm_section = -0.05
    Cl_section = 1.7

    # Plot com CG físico em 20% da MAC
    glider.plot_wing_with_cg_cp_ac(Cm=Cm_section, Cl=Cl_section, cg_phys_frac_mac=0.20)

    # CP por Cm/Cl:
    glider.plot_half_wing_right(Cm=-0.05, Cl=1.7, cg_phys_frac_mac=0.20)

    # ou CP dado diretamente por percentagem da MAC (ex.: 30%):
    glider.plot_half_wing_right(f_cp=0.30, cg_phys_frac_mac=0.20)