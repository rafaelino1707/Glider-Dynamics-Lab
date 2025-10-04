# heatmaps_sink.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ----------------------------- Glider class --------------------------------
class Glider:
    def __init__(self, mass_kg, chord_m, span_m, CLmax, oswald_efficiency=0.75, cd0=0.04):
        self.mass = mass_kg
        self.chord = chord_m
        self.span = span_m
        self.CLmax = CLmax
        self.e = oswald_efficiency
        self.cd0 = cd0

    @property
    def wing_area(self):
        return max(self.chord * self.span, 1e-12)

    @property
    def aspect_ratio(self):
        S = self.wing_area
        return (self.span ** 2) / S

    @property
    def induced_drag_factor(self):
        AR = self.aspect_ratio
        return 1.0 / (math.pi * AR * self.e)

    def Vmin_theoretical(self, rho=1.225):
        W = self.mass * 9.81
        S = self.wing_area
        K = self.induced_drag_factor
        return math.sqrt(2.0 / rho) * (K / self.cd0)**0.25 * math.sqrt(W / S)

    def sink_at_V(self, V, rho=1.225):
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
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.CLmax))

# ----------------------------- Settings ------------------------------------
chord_vals = np.linspace(0.02, 0.15, 41)
span_vals  = np.linspace(0.5, 1.5, 41)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

CL_MAX_DEFAULT = 1.2
OSWALD_EFF = 0.75
CD0 = 0.04
RHO = 1.225
COLORMAP = 'Greys'

now = datetime.now().strftime('%d_%m_%Y_%H_%M')
out_dir = f'Sink_Analysis_{now}'
os.makedirs(out_dir, exist_ok=True)

# ---------------------------- Compute grid --------------------------------
print("Running sweep over masses:", masses)
rows = []
for mass in masses:
    for chord in chord_vals:
        for span in span_vals:
            g = Glider(mass, chord, span, CLmax=CL_MAX_DEFAULT, oswald_efficiency=OSWALD_EFF, cd0=CD0)
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
