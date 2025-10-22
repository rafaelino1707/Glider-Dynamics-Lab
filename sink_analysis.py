# sink_analysis.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ----------------------------- Glider class --------------------------------
class Glider:
    def __init__(
        self,
        mass_kg,
        wingspan_m,
        c_upper_m,
        c_lower_m,
        cl_max=1.7187,           # nome padronizado para evitar conflitos
        oswald_efficiency=0.75,
        cd0=0.03,
    ):
        self.mass = mass_kg
        self.wingspan = wingspan_m
        self.c_upper = c_upper_m
        self.c_lower = c_lower_m
        self.cl_max = cl_max
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
    # AERODINÂMICA
    # ============================

    @property
    def induced_drag_factor(self):
        """K = 1 / (π * AR * e)"""
        return 1 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    @property
    def cl_min_sink(self):
        """CL_ms = sqrt(3*CD0/K)"""
        return math.sqrt(3.0 * self.cd0 / self.induced_drag_factor)

    def v_min_sink(self, rho=1.225):
        """V_ms = (2/ρ)^{1/2} * (K/(3*CD0))^{1/4} * (W/S)^{1/2}"""
        W = self.mass * 9.81
        S = self.wing_area
        K = self.induced_drag_factor
        return math.sqrt(2.0 / rho) * (K / (3.0 * self.cd0)) ** 0.25 * math.sqrt(W / S)

    def sink_at_V(self, V, rho=1.225):
        """Taxa de sink a uma velocidade V"""
        W = self.mass * 9.81
        S = self.wing_area
        CL = W / (0.5 * rho * V**2 * S)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2
        D = 0.5 * rho * V**2 * S * CD
        sink = V * (D / W)
        glide = (W / D) if D > 0 else float('nan')
        return {"V": V, "CL": CL, "CD": CD, "D": D, "sink": sink, "glide": glide}

    def sink_min(self, rho=1.225):
        """Taxa de sink mínima usando CL_ms"""
        Vms = self.v_min_sink(rho=rho)
        CL = self.cl_min_sink
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2
        Vy = Vms * (CD / CL)
        return {"Vms": Vms, "CL_ms": CL, "CD_ms": CD, "sink_min": Vy}

    def v_stall(self, rho=1.225):
        """V_stall = sqrt((2W)/(ρ S CLmax))"""
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.cl_max))

    @staticmethod
    def reynolds(V, chord, rho=1.225, mu=1.81e-5):
        """Re = ρ V c / μ"""
        return rho * V * chord / mu

# ----------------------------- Settings ------------------------------------
chord_vals = np.linspace(0.01, 0.15, 41)
span_vals  = np.linspace(0.2, 1.5, 41)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

CL_MAX_DEFAULT = 1.7187
OSWALD_EFF = 0.9
CD0 = 0.03
RHO = 1.225
COLORMAP = 'Greys'

now = datetime.now().strftime('%d_%m_%Y_%H_%M')
out_dir = f'Log/Sink_Analysis/Sink_Analysis_{now}'
os.makedirs(out_dir, exist_ok=True)

# ---------------------------- Compute grid --------------------------------
print("Running sweep over masses:", masses)
rows = []
for mass in masses:
    for chord in chord_vals:
        for span in span_vals:
            c_upper = 0.25 * chord
            c_lower = 0.75 * chord
            g = Glider(
                mass, span, c_upper, c_lower,
                cl_max=CL_MAX_DEFAULT,
                oswald_efficiency=OSWALD_EFF,
                cd0=CD0
            )

            Vms = g.v_min_sink(rho=RHO)
            sink_info = g.sink_at_V(Vms, rho=RHO)
            sink_min = g.sink_min(rho=RHO)
            Vstall = g.v_stall(rho=RHO)

            rows.append({
                "mass_kg": mass,
                "chord_m": chord,
                "span_m": span,
                "wing_area_m2": g.wing_area,
                "AR": g.aspect_ratio,
                "K": g.induced_drag_factor,
                "Vmin_theo_m_s": Vms,
                "sink_rate_m_s": sink_info["sink"],
                "CL_at_Vmin": sink_info["CL"],
                "CD_at_Vmin": sink_info["CD"],
                "glide_at_Vmin": sink_info["glide"],
                "CL_ms": sink_min["CL_ms"],
                "CD_ms": sink_min["CD_ms"],
                "sink_min_m_s": sink_min["sink_min"],
                "Vstall_m_s": Vstall,
                "Re_root_at_Vmin": Glider.reynolds(Vms, chord, rho=RHO)
            })

df = pd.DataFrame(rows)
full_csv = os.path.join(out_dir, 'grid_full.csv')
df.to_csv(full_csv, index=False)
print("Saved full grid CSV to:", full_csv)

# ---------------------------- Plotting ------------------------------------
sink_dir = os.path.join(out_dir, 'sink_heatmaps')
vmin_dir = os.path.join(out_dir, 'vmin_heatmaps')
os.makedirs(sink_dir, exist_ok=True)
os.makedirs(vmin_dir, exist_ok=True)

for mass, group in df.groupby('mass_kg'):
    sub = group.copy()
    pivot_sink = sub.pivot_table(index='chord_m', columns='span_m',
                                 values='sink_min_m_s', aggfunc='mean')
    pivot_vmin = sub.pivot_table(index='chord_m', columns='span_m',
                                 values='Vmin_theo_m_s', aggfunc='mean')
    X = pivot_sink.columns.values
    Y = pivot_sink.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    # --- sink heatmap ---
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(pivot_sink.values, origin='lower', aspect='auto',
                   extent=extent, cmap=COLORMAP)
    ax.set_title(f'sink_min (m/s) — mass={mass} kg')
    ax.set_xlabel('Span (m)')
    ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax)
    c.set_label('sink_min (m/s)')
    try:
        CS = ax.contour(X, Y, pivot_vmin.values, colors='black', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.1f')
    except Exception:
        pass
    fig.savefig(os.path.join(sink_dir, f'heat_sink_min_mass{mass:.1f}.png'), dpi=200)
    plt.close(fig)

    # --- Vmin heatmap ---
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(pivot_vmin.values, origin='lower', aspect='auto',
                   extent=extent, cmap=COLORMAP)
    ax.set_title(f'Vmin_sink (m/s) — mass={mass} kg')
    ax.set_xlabel('Span (m)')
    ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax)
    c.set_label('Vmin (m/s)')
    try:
        CS = ax.contour(X, Y, pivot_sink.values, colors='black', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.2f')
    except Exception:
        pass
    fig.savefig(os.path.join(vmin_dir, f'heat_vmin_mass{mass:.1f}.png'), dpi=200)
    plt.close(fig)

print("\nAll done. Check folder:", out_dir)
