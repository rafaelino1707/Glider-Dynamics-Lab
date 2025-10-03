# heatmaps_all_masses.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ----------------------------- Glider class --------------------------------
class Glider:
    """
    Modelo simples com CLmax disponível como self.CLmax.
    Calcula Vmin_teórico (velocidade que minimiza sink segundo o teu modelo),
    sink@V, e Vstall a partir de CLmax.
    """
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
        # K = 1/(pi*AR*e)
        return 1.0 / (math.pi * AR * self.e)

    def Vmin_theoretical(self, rho=1.225):
        """
        Velocidade que teoricamente minimiza a taxa de sink (teoria usada no teu modelo).
        """
        W = self.mass * 9.81
        S = self.wing_area
        K = self.induced_drag_factor
        return math.sqrt(2.0 / rho) * (K / self.cd0)**0.25 * math.sqrt(W / S)

    def sink_at_V(self, V, rho=1.225):
        """
        Calcula sink_rate, CL, CD, D e glide para uma velocidade V.
        """
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
        """
        Velocidade de estol a partir de CLmax.
        Vs = sqrt(2W / (rho * S * CLmax))
        """
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.CLmax))


# ----------------------------- Settings ------------------------------------
# geometry grid
chord_vals = np.linspace(0.02, 0.10, 41)   # corda 2cm -> 10cm
span_vals  = np.linspace(0.45, 0.85, 41)   # span 0.45m -> 0.85m

# masses to sweep (0.1 .. 1.0 kg)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

# aerodynamic params (edit if you want)
CL_MAX_DEFAULT = 1.2     # substitute with your CLmax value (from CL vs AoA curve)
OSWALD_EFF = 0.75
CD0 = 0.04
RHO = 1.225

# colour map (user requested red/yellow style) -> use 'hot' (red->yellow->white)
COLORMAP = 'plasma'

# output directory
now = datetime.now().strftime('%d_%m_%Y_%H_%M')
out_dir = f'Sink_and_Stall_Analysis_{now}'
os.makedirs(out_dir, exist_ok=True)

# ---------------------------- Compute grid --------------------------------
print("Running sweep over masses:", masses)
rows = []
for mass in masses:
    print(f"  mass = {mass} kg")
    for chord in chord_vals:
        for span in span_vals:
            g = Glider(mass, chord, span, CLmax=CL_MAX_DEFAULT, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            S = g.wing_area
            # compute theoretical Vmin and sink at that Vmin
            Vmin = g.Vmin_theoretical(rho=RHO)
            sink_info = g.sink_at_V(Vmin, rho=RHO)
            Vstall = g.Vstall(rho=RHO)

            rows.append({
                "mass_kg": mass,
                "chord_m": chord,
                "span_m": span,
                "wing_area_m2": S,
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

# ---------------------------- Plotting -----------------------------------
# We'll create three separate folders inside out_dir for clarity
sink_dir   = os.path.join(out_dir, 'sink_heatmaps')
vmin_dir   = os.path.join(out_dir, 'vmin_heatmaps')
vstall_dir = os.path.join(out_dir, 'vstall_heatmaps')
os.makedirs(sink_dir, exist_ok=True)
os.makedirs(vmin_dir, exist_ok=True)
os.makedirs(vstall_dir, exist_ok=True)

for mass, group in df.groupby('mass_kg'):
    sub = group.copy()

    # pivot tables (index = chord, columns = span)
    pivot_sink  = sub.pivot_table(index='chord_m', columns='span_m', values='sink_rate_m_s', aggfunc='mean')
    pivot_vmin  = sub.pivot_table(index='chord_m', columns='span_m', values='Vmin_theo_m_s', aggfunc='mean')
    pivot_vstall= sub.pivot_table(index='chord_m', columns='span_m', values='Vstall_m_s', aggfunc='mean')

    # grid coordinates (for contour and extent)
    X = pivot_sink.columns.values
    Y = pivot_sink.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    # -------- 1) sink heatmap with Vmin contours (white) ---------------
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'sink_rate (m/s) at Vmin — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('sink_rate (m/s)')

    # overlay Vmin contours (white)
    try:
        CS = ax.contour(X, Y, pivot_vmin.values, colors='white', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.1f')
    except Exception:
        pass

    fname = os.path.join(sink_dir, f'heat_sink_mass{mass:.1f}.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    # -------- 2) Vmin heatmap with sink contours (black) ---------------
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'Vmin_theoretical (m/s) — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('Vmin (m/s)')

    # overlay sink contours (black)
    try:
        CS = ax.contour(X, Y, pivot_sink.values, colors='white', linewidths=0.8)
        ax.clabel(CS, inline=True, fontsize=8, fmt='%.2f')
    except Exception:
        pass

    fname = os.path.join(vmin_dir, f'heat_vmin_mass{mass:.1f}.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    # -------- 3) Vstall heatmap (from CLmax) ---------------------------
    fig, ax = plt.subplots(figsize=(7,5))
    im = ax.imshow(pivot_vstall.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'Vstall (m/s) — CLmax={CL_MAX_DEFAULT} — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('Vstall (m/s)')

    # optional: overlay sink or Vmin contours if you like (commented)
    # try:
    #     CS = ax.contour(X, Y, pivot_sink.values, colors='white', linewidths=0.6)
    #     ax.clabel(CS, inline=True, fontsize=8)
    # except Exception:
    #     pass

    fname = os.path.join(vstall_dir, f'heat_vstall_mass{mass:.1f}.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    print(f"Saved heatmaps for mass {mass:.1f} kg ->")
    print("  ", fname)

print("\nAll done. Check folder:", out_dir)
