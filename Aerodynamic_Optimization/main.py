import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

class Glider:
    def __init__(self, mass, wingspan, c_upper, c_lower, oswald_efficiency, cd0):
        self.mass = mass
        self.wingspan = wingspan
        self.c_upper = c_upper
        self.c_lower = c_lower
        self.oswald_efficiency = oswald_efficiency
        self.cd0 = cd0

    @property
    def total_chord(self):
        return self.c_upper + self.c_lower

    @property
    def wing_area(self):
        return (math.pi * self.wingspan / 4) * self.total_chord

    @property
    def mean_aerodynamic_chord(self):
        return (8 / (3 * math.pi)) * self.total_chord

    @property
    def aspect_ratio(self):
        return 4 * self.wingspan / (math.pi * self.total_chord)

    @property
    def cg_chordwise(self):
        return (4 / (3 * math.pi)) * (self.c_upper - self.c_lower)

    @property
    def cg_spanwise_half(self):
        return (2 * self.wingspan) / (3 * math.pi)

    @property
    def induced_drag_factor(self):
        return 1 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    def minimum_sink_velocity(self, rho=1.225):
        W = self.mass * 9.81
        term_1 = math.sqrt(2 / rho)
        term_2 = (self.induced_drag_factor / self.cd0) ** 0.25
        term_3 = math.sqrt(W / self.wing_area)
        return term_1 * term_2 * term_3

    def vertical_sink_rate(self, V=None, rho=1.225):
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
        return {"V_used": V, "CL": CL, "CD": CD, "D_N": D,
                "sink_rate_m_s": sink_rate, "glide_ratio": glide_ratio}

    def stall_velocity(self, CL_max, rho=1.225):
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2 * W) / (rho * S * CL_max))

    def chordwise_cp(self, Cm, Cl):
        xcp_over_mac = 0.25 - (Cm / Cl)
        return xcp_over_mac * self.mean_aerodynamic_chord

    @property
    def spanwise_cp_half(self):
        return (2 * self.wingspan) / (3 * math.pi)

    def print_cp(self, Cm, Cl):
        xcp = self.chordwise_cp(Cm, Cl)
        mac = self.mean_aerodynamic_chord
        ycp = self.spanwise_cp_half
        print(f"x_CP = {xcp:.5f} m  ({xcp/mac:.3f} × MAC)")
        print(f"y_CP (half-span) = {ycp:.5f} m  (~0.212 × b)")

    def _wing_edges(self, b, c_up, c_low, y):
        x_upper = c_up * np.sqrt(1 - (4 * y**2) / b**2)
        x_lower = -c_low * np.sqrt(1 - (4 * y**2) / b**2)
        return x_upper, x_lower

    def plot_glider_model(self,
                          fuse_len=0.50,
                          fuse_w=0.035,
                          tail_span=0.25,
                          tail_scale=0.55,
                          tail_x_offset=0.35):

        b = self.wingspan
        y = np.linspace(-b/2, b/2, 600)
        x_u, x_l = self._wing_edges(b, self.c_upper, self.c_lower, y)

        b_t = tail_span
        y_t = np.linspace(-b_t/2, b_t/2, 400)

        c_up_t = self.c_upper * tail_scale
        c_lo_t = self.c_lower * tail_scale
        x_u_t, x_l_t = self._wing_edges(b_t, c_up_t, c_lo_t, y_t)

        tail_x_offset = -0.285
        tail_y_offset = 0
        x_u_t += tail_x_offset
        x_l_t += tail_x_offset
        y_t += tail_y_offset

        x_np = 0.00321
        y_np = 0
        y_cg = 0
        x_cg = 0.0168
        y_cp_full = 0
        x_cp = -0.0002
        fuse_center_x = -0.25/2

        # --- inversão vertical completa ---
        invert = -1
        x_u *= invert
        x_l *= invert
        x_u_t *= invert
        x_l_t *= invert
        x_cg *= invert
        x_cp *= invert
        x_np *= invert
        fuse_center_x *= invert

        fig, ax = plt.subplots(figsize=(9, 5))

        plt.scatter(y_cg, x_cg, s=15, zorder=5, color='darkred', label=r'$x_{CG} = -0.0168$')
        plt.scatter(y_cp_full, x_cp, s=15, zorder=5, color='darkgreen', label=r'$x_{CP} = 0.0002$')
        plt.scatter(y_np, x_np, s=15, zorder=6, color='darkgoldenrod', label=r'$x_{NP} = -0.00321$')

        # --- MAC numa semi-asa (direita) ---
        MAC = self.mean_aerodynamic_chord
        b = self.wingspan
        ratio = 8 / (3 * math.pi)                     # MAC / (c_upper + c_lower)
        y_mac = (b/2) * math.sqrt(1 - ratio**2)       # posição spanwise da MAC

        # bordos da corda local em y_mac
        x_le_mac = self.c_upper * math.sqrt(1 - (4 * y_mac**2) / b**2)
        x_te_mac = -self.c_lower * math.sqrt(1 - (4 * y_mac**2) / b**2)

        # aplicar a inversão vertical usada no resto da figura
        x_le_mac *= invert
        x_te_mac *= invert

        # desenhar a MAC a tracejado
        ax.plot([y_mac, y_mac], [x_te_mac, x_le_mac],
                linestyle='--', linewidth=1.6, label='MAC (wing)')
        ax.plot(y, x_u, label="LE Wing")
        ax.plot(y, x_l, label="TE Wing")

        body = Ellipse((0.0, fuse_center_x), width=fuse_w, height=fuse_len,
                       fill=False, linewidth=1.4)
        ax.add_patch(body)

        ax.plot(y_t, x_u_t, label="LE Stabilizer")
        ax.plot(y_t, x_l_t, label="TE Stabilizer")

        ax.axhline(0, linestyle='--', linewidth=0.8)
        ax.axvline(0, linestyle='--', linewidth=0.8)

        ax.set_title("Glider Top Visualization")
        ax.set_xlabel("y (spanwise, m)")
        ax.set_ylabel("x (chordwise, m)")
        ax.legend()
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True)
        ax.invert_yaxis()

        plt.show()


# ============================
# EXEMPLO DE USO
# ============================

glider = Glider(
    mass=0.2,
    wingspan=0.6,
    c_upper=0.02,
    c_lower=0.06,
    oswald_efficiency=0.09,
    cd0=0.04
)

print("=== GEOMETRIA ===")
print(f"Wing Area S = {glider.wing_area:.5f} m²")
print(f"Mean Aerodynamic Chord (MAC) = {glider.mean_aerodynamic_chord:.5f} m")
print(f"Aspect Ratio (AR) = {glider.aspect_ratio:.3f}")
print(f"CG chordwise = {glider.cg_chordwise:.5f} m")
print(f"CG spanwise (half-wing) = {glider.cg_spanwise_half:.5f} m")

Cm_section = -0.05
Cl_section = 0.6
print("\n=== CENTRO DE PRESSÃO ===")
glider.print_cp(Cm=Cm_section, Cl=Cl_section)

glider.plot_glider_model(
    fuse_len=0.42,
    fuse_w=0.035,
    tail_span=0.24,
    tail_scale=0.55,
    tail_x_offset=0.35
)
