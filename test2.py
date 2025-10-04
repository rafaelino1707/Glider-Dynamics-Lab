# heatmaps_multiCL_fixed.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
# grid (ajusta se precisares)
chord_vals = np.linspace(0.02, 0.15, 41)   # corda 2cm->15cm
span_vals  = np.linspace(0.5, 1.5, 41)     # span 0.5->1.5 m
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

# aerodynamics
CL_MAX_DEFAULT = 1.2   # used for baseline Vmin/sink
OSWALD_EFF = 0.75
CD0 = 0.04
RHO = 1.225

# lista dos 4 airfoils (CL_max para cada um) - substitui pelos teus valores reais
CL_LIST = [0.9, 1.2, 1.5, 1.8]

# colors for contour lines (one per CL)
overlay_line_colors = ['#d73027', '#f46d43', '#fdae61', '#fee08b']  # contrastantes (vermelho -> amarelo)
line_width = 1.6

# base cmap and global scale control
BASE_CMAP = 'inferno'   # bom contraste vermelho->amarelo

# output dir
now = datetime.now().strftime('%Y%m%d_%H%M%S')
out_dir = f'heatmaps_multiCL_fixed_{now}'
os.makedirs(out_dir, exist_ok=True)
sink_dir   = os.path.join(out_dir, 'sink_heatmaps')
vmin_dir   = os.path.join(out_dir, 'vmin_heatmaps')
vstall_dir = os.path.join(out_dir, 'vstall_multiCL_lines')
os.makedirs(sink_dir, exist_ok=True)
os.makedirs(vmin_dir, exist_ok=True)
os.makedirs(vstall_dir, exist_ok=True)

# ---------------------------- Compute base grids ---------------------------
print("Computing grid for masses:", masses)
rows = []
for mass in masses:
    for chord in chord_vals:
        for span in span_vals:
            g = Glider(mass, chord, span, CLmax=CL_MAX_DEFAULT, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            Vmin = g.Vmin_theoretical(rho=RHO)
            sink_info = g.sink_at_V(Vmin, rho=RHO)
            rows.append({
                "mass": mass,     # consistent column name 'mass'
                "chord": chord,   # consistent 'chord'
                "span": span,     # consistent 'span'
                "Vmin": Vmin,
                "sink": sink_info["sink"]
            })
df_base = pd.DataFrame(rows)

# precompute Vstall for each CL in CL_LIST (store in dict of dataframes)
vstall_dict = {}
for CLval in CL_LIST:
    rows_vs = []
    for mass in masses:
        for chord in chord_vals:
            for span in span_vals:
                g = Glider(mass, chord, span, CLmax=CLval, oswald_efficiency=OSWALD_EFF, cd0=CD0)
                rows_vs.append({
                    "mass": mass,   # same column name 'mass'
                    "chord": chord,
                    "span": span,
                    "Vstall": g.Vstall(rho=RHO),
                    "CLmax": CLval
                })
    vstall_dict[CLval] = pd.DataFrame(rows_vs)

# compute global vmin/vmax for consistent color scaling
global_vmin_vmin = df_base['Vmin'].min()
global_vmax_vmin = df_base['Vmin'].max()
global_vmin_sink = df_base['sink'].min()
global_vmax_sink = df_base['sink'].max()

print("Global Vmin range:", global_vmin_vmin, "->", global_vmax_vmin)
print("Global sink range:", global_vmin_sink, "->", global_vmax_sink)

# ---------------------------- Plot per mass --------------------------------
for mass, group in df_base.groupby('mass'):
    print(f"Plotting mass {mass} kg ...")
    sub = group.copy()
    # pivot table index=chord, columns=span (names match df_base)
    pivot_sink = sub.pivot_table(index='chord', columns='span', values='sink', aggfunc='mean')
    pivot_vmin = sub.pivot_table(index='chord', columns='span', values='Vmin', aggfunc='mean')

    X = pivot_sink.columns.values
    Y = pivot_sink.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    # ---- sink heatmap (with Vmin contours white) ----
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=extent,
                   cmap=BASE_CMAP, vmin=global_vmin_sink, vmax=global_vmax_sink)
    ax.set_title(f'sink_rate (m/s) @ Vmin — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    cbar = fig.colorbar(im, ax=ax); cbar.set_label('sink (m/s)')

    # Vmin contours (white) - few levels for clarity
    try:
        levels_vmin = np.linspace(np.nanmin(pivot_vmin.values), np.nanmax(pivot_vmin.values), 6)
        CSv = ax.contour(X, Y, pivot_vmin.values, levels=levels_vmin, colors='white', linewidths=0.8)
        ax.clabel(CSv, inline=True, fontsize=8, fmt='%.1f', colors='white')
    except Exception:
        pass

    fname = os.path.join(sink_dir, f'heat_sink_mass{mass:.1f}.png')
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ---- vmin heatmap (with sink contours black) ----
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=extent,
                   cmap=BASE_CMAP, vmin=global_vmin_vmin, vmax=global_vmax_vmin)
    ax.set_title(f'Vmin_theoretical (m/s) — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    cbar = fig.colorbar(im, ax=ax); cbar.set_label('Vmin (m/s)')
    try:
        levels_sink = np.linspace(np.nanmin(pivot_sink.values), np.nanmax(pivot_sink.values), 6)
        CSs = ax.contour(X, Y, pivot_sink.values, levels=levels_sink, colors='black', linewidths=0.8)
        ax.clabel(CSs, inline=True, fontsize=8, fmt='%.2f', colors='black')
    except Exception:
        pass
    fname = os.path.join(vmin_dir, f'heat_vmin_mass{mass:.1f}.png')
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close(fig)

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------- Glider class ----------------
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

    def Vstall(self, rho=1.225):
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.CLmax))

# ---------------- Settings ----------------
chord_vals = np.linspace(0.02, 0.15, 80)
span_vals  = np.linspace(0.5, 1.5, 80)

MASS = 0.2
CL_LIST = [0.9, 1.2, 1.5, 1.8]
RHO = 1.225

# ---------------- Compute ----------------
rows = []
for chord in chord_vals:
    for span in span_vals:
        for CL in CL_LIST:
            g = Glider(MASS, chord, span, CL)
            Vstall = g.Vstall(rho=RHO)
            rows.append({"mass": MASS, "chord": chord, "span": span, "CLmax": CL, "Vstall": Vstall})

df = pd.DataFrame(rows)

# ---------------- Plot ----------------
fig, ax = plt.subplots(figsize=(7,5))

# fundo com Vstall mínimo (em cinza)
pivot_min = df.groupby(["chord","span"])["Vstall"].min().unstack()
X, Y = np.meshgrid(pivot_min.columns.values, pivot_min.index.values)
im = ax.imshow(pivot_min.values, origin='lower', aspect='auto',
               extent=[span_vals.min(), span_vals.max(), chord_vals.min(), chord_vals.max()],
               cmap='Greys', alpha=0.9)

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Vstall mínimo (m/s)")

# cor verde-lima
overlay_line_colors = ["#32CD32"] * len(CL_LIST)

legend_lines = []
for i, CLval in enumerate(CL_LIST):
    vsub = df[df["CLmax"] == CLval]
    pivot_vs = vsub.pivot_table(index='chord', columns='span', values='Vstall', aggfunc='mean')
    Z = pivot_vs.values

    try:
        lvl1 = np.nanpercentile(Z, 30)
        lvl2 = np.nanpercentile(Z, 50)
        lvl3 = np.nanpercentile(Z, 70)
        CS = ax.contour(X, Y, Z, levels=[lvl1, lvl2, lvl3],
                        colors=overlay_line_colors[i], linewidths=1.5, linestyles='-')
        ax.clabel(CS, fmt='%.1f', fontsize=8, inline=True)
    except Exception:
        med = np.nanmedian(Z)
        ax.contour(X, Y, Z, levels=[med], colors=overlay_line_colors[i], linewidths=1.5)

    legend_lines.append(Line2D([0],[0], color=overlay_line_colors[i], lw=2, label=f'CLmax={CLval}'))

ax.set_title(f"Vstall contours (green-lime) — mass={MASS} kg")
ax.set_xlabel("Span (m)")
ax.set_ylabel("Chord (m)")

ax.legend(handles=legend_lines, loc='upper right')

plt.tight_layout()
plt.show()
