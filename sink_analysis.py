# heatmaps_sink.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ----------------------------- Glider class --------------------------------
class Glider:
    def __init__(self, mass_kg, wingspan_m, c_upper_m, c_lower_m, oswald_efficiency=0.75, cd0=0.04):
        self.mass = mass_kg
        self.wingspan = wingspan_m
        self.c_upper = c_upper_m
        self.c_lower = c_lower_m
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

    def Vmin_theoretical(self, rho=1.225):
        """V_min ≈ sqrt(2/rho) * (K/CD0)^(1/4) * sqrt(W/S)"""
        W = self.mass * 9.81
        S = self.wing_area
        K = self.induced_drag_factor
        return math.sqrt(2.0 / rho) * (K / self.cd0)**0.25 * math.sqrt(W / S)

    def sink_at_V(self, V, rho=1.225):
        """Calcula a taxa de sink a uma velocidade V"""
        W = self.mass * 9.81
        S = self.wing_area
        CL = W / (0.5 * rho * V**2 * S)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2
        D = 0.5 * rho * V**2 * S * CD
        sink = V * (D / W)
        glide = (W / D) if D > 0 else float('nan')
        return {"V": V, "CL": CL, "CD": CD, "D": D, "sink": sink, "glide": glide}

    def Vstall(self, rho=1.225):
        """V_stall = sqrt((2W)/(ρ S CLmax))"""
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.CLmax))

# ----------------------------- Settings ------------------------------------
chord_vals = np.linspace(0.01, 0.15, 41)
span_vals  = np.linspace(0.2, 1.5, 41)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

CL_MAX_DEFAULT = 1.8
OSWALD_EFF = 0.9
CD0 = 0.04
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
            c_upper = 0.25 * chord  # 25% da corda total para c_upper
            c_lower = 0.75 * chord  # 75% da corda total para c_lower
            g = Glider(mass, span, c_upper, c_lower, CLmax=CL_MAX_DEFAULT, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            Vmin = g.Vmin_theoretical(rho=RHO)
            sink_info = g.sink_at_V(Vmin, rho=RHO)
            Vstall = g.Vstall(rho=RHO)

            rows.append({
                "mass_kg": mass,
                "chord_m": chord,
                "span_m": span,
                "wing_area_m2": g.wing_area,
                "Vmin_theo_m_s": Vmin,
                "sink_rate_m_s": sink_info["sink"],
                "CL_at_Vmin": sink_info["CL"],
                "CD_at_Vmin": sink_info["CD"],
                "glide_at_Vmin": sink_info["glide"],
                "Vstall_m_s": Vstall
            })

df = pd.DataFrame(rows)
full_csv = os.path.join(out_dir, 'grid_full.csv')
df.to_csv(full_csv, index=False)
print("Saved full grid CSV to:", full_csv)

# ---------------------------- Plotting (sink only) -------------------------
sink_dir   = os.path.join(out_dir, 'sink_heatmaps')
vmin_dir   = os.path.join(out_dir, 'vmin_heatmaps')
os.makedirs(sink_dir, exist_ok=True)
os.makedirs(vmin_dir, exist_ok=True)

for mass, group in df.groupby('mass_kg'):
    sub = group.copy()
    pivot_sink  = sub.pivot_table(index='chord_m', columns='span_m', values='sink_rate_m_s', aggfunc='mean')
    pivot_vmin  = sub.pivot_table(index='chord_m', columns='span_m', values='Vmin_theo_m_s', aggfunc='mean')
    X = pivot_sink.columns.values
    Y = pivot_sink.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    # -------- 1) sink heatmap with Vmin contours ---------------
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'sink_rate (m/s) at Vmin — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('sink_rate (m/s)')
    try:
        CS = ax.contour(X, Y, pivot_vmin.values, colors='black', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.1f')
    except Exception: pass
    fig.savefig(os.path.join(sink_dir, f'heat_sink_mass{mass:.1f}.png'), dpi=200)
    plt.close(fig)

    # -------- 2) Vmin heatmap with sink contours ---------------
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'Vmin_theoretical (m/s) — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('Vmin (m/s)')
    try:
        CS = ax.contour(X, Y, pivot_sink.values, colors='black', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.2f')
    except Exception: pass
    fig.savefig(os.path.join(vmin_dir, f'heat_vmin_mass{mass:.1f}.png'), dpi=200)
    plt.close(fig)

print("\nAll done. Check folder:", out_dir)
