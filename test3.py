# combined_vmin_sink_heatmaps.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ---------------- Model (mesma formula que usaste) ----------------
class GliderSimple:
    def __init__(self, mass_kg, chord_m, span_m, oswald_efficiency=0.75, cd0=0.04):
        self.mass = mass_kg
        self.chord = chord_m
        self.span = span_m
        self.e = oswald_efficiency
        self.cd0 = cd0

    @property
    def wing_area(self):
        return max(self.chord * self.span, 1e-12)

    @property
    def aspect_ratio(self):
        return (self.span**2) / self.wing_area

    @property
    def induced_drag_factor(self):
        return 1.0 / (math.pi * self.aspect_ratio * self.e)

    def vmin_theoretical(self, rho=1.225):
        W = self.mass * 9.81
        S = self.wing_area
        K = self.induced_drag_factor
        return math.sqrt(2.0 / rho) * (K / self.cd0)**0.25 * math.sqrt(W / S)

    def sink_at_v(self, V, rho=1.225):
        W = self.mass * 9.81
        S = self.wing_area
        CL = W / (0.5 * rho * V**2 * S)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL**2
        D = 0.5 * rho * V**2 * S * CD
        sink = V * (D / W)
        glide = W / D if D > 0 else np.nan
        return {"V": V, "CL": CL, "CD": CD, "D": D, "sink": sink, "glide": glide}

# ---------------- Settings (edita se precisares) ----------------
chord_vals = np.linspace(0.02, 0.10, 41)   # corda 2cm -> 10cm
span_vals = np.linspace(0.45, 0.85, 41)    # span 0.45m -> 0.85m
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)  # 0.1, 0.2, ..., 1.0 kg

OSWALD_EFF = 0.75
CD0 = 0.04
RHO = 1.225

# output dir
now = datetime.now().strftime('%Y%m%d_%H%M%S')
out_dir = f'combined_vmin_sink_heatmaps_{now}'
os.makedirs(out_dir, exist_ok=True)

# quais massas mostrar inline (apenas se corres num notebook)
show_masses_inline = [0.3, 0.5, 1.0]

# ---------------- Compute grid -----------------------------------
rows = []
for mass in masses:
    for chord in chord_vals:
        for span in span_vals:
            g = GliderSimple(mass, chord, span, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            Vmin = g.vmin_theoretical(rho=RHO)
            sink_info = g.sink_at_v(Vmin, rho=RHO)
            rows.append({
                "mass": mass,
                "chord": chord,
                "span": span,
                "Vmin": Vmin,
                "sink": sink_info["sink"],
                "CL": sink_info["CL"],
                "CD": sink_info["CD"],
                "glide": sink_info["glide"]
            })

df = pd.DataFrame(rows)
full_csv = os.path.join(out_dir, 'vmin_sink_grid_full.csv')
df.to_csv(full_csv, index=False)
print("Saved full grid CSV to:", full_csv)

# ---------------- Plot combined figures per mass ----------------------
saved_files = []
for mass, group in df.groupby('mass'):
    sub = group.copy()
    # pivot tables for heatmaps (index = chord, columns = span)
    pivot_sink = sub.pivot_table(index='chord', columns='span', values='sink', aggfunc='mean')
    pivot_vmin = sub.pivot_table(index='chord', columns='span', values='Vmin', aggfunc='mean')

    X = pivot_sink.columns.values
    Y = pivot_sink.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax1, ax2 = axes

    # Left: sink heatmap, overlay Vmin contours (white)
    im1 = ax1.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=extent)
    ax1.set_title(f'sink_rate at Vmin — mass={mass} kg')
    ax1.set_xlabel('Span (m)'); ax1.set_ylabel('Chord (m)')
    c1 = fig.colorbar(im1, ax=ax1); c1.set_label('sink (m/s)')
    try:
        CS = ax1.contour(X, Y, pivot_vmin.values, colors='white', linewidths=0.8)
        ax1.clabel(CS, inline=True, fontsize=8, fmt='%.1f')
    except Exception:
        pass

    # Right: Vmin heatmap, overlay sink contours (black)
    im2 = ax2.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=extent)
    ax2.set_title(f'Vmin_theoretical — mass={mass} kg')
    ax2.set_xlabel('Span (m)'); ax2.set_ylabel('Chord (m)')
    c2 = fig.colorbar(im2, ax=ax2); c2.set_label('Vmin (m/s)')
    try:
        CS2 = ax2.contour(X, Y, pivot_sink.values, colors='black', linewidths=0.8)
        ax2.clabel(CS2, inline=True, fontsize=8, fmt='%.2f')
    except Exception:
        pass

    fig.suptitle(f'Mass = {mass} kg — sink (color) + Vmin (contours)  |  Vmin (color) + sink (contours)', fontsize=12)
    plt.tight_layout(rect=[0,0,1,0.96])

    fname = os.path.join(out_dir, f'combined_mass_{mass}.png')
    fig.savefig(fname, dpi=200)
    saved_files.append(fname)

    # show inline if desired (works in Jupyter)
    plt.close(fig) 


print("Saved combined heatmaps to folder:", out_dir)
print("Example files:", saved_files[:6])

# optional: print top combos per mass (lowest sink)
summary = df.groupby('mass').apply(lambda g: g.nsmallest(6, 'sink')[['chord','span','Vmin','sink']]).reset_index(level=1, drop=True).reset_index()
print("\nTop combos per mass (lowest sink):")
print(summary.head(24))
