# heatmaps_stall.py
import math, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ----------------------------- Glider class --------------------------------
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
    # AERODINÂMICA
    # ============================

    @property
    def induced_drag_factor(self):
        """K = 1 / (π * AR * e)"""
        return 1 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    def stall_velocity(self, CL_max, rho=1.225):
        """V_stall = sqrt((2W)/(ρ S CL_max))"""
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2 * W) / (rho * S * CL_max))

# ----------------------------- Settings ------------------------------------
chord_vals = np.linspace(0.02, 0.1, 41)  # Gerando o valor total da corda
span_vals = np.linspace(0.5, 1.5, 41)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

CL_MAX_DEFAULT =  1.1836
OSWALD_EFF = 0.9
CD0 = 0.04
RHO = 1.225
COLORMAP = 'Greys'

now = datetime.now().strftime('%d_%m_%Y_%H_%M')
out_dir = f'Log/Stall_Analysis/Stall_Analysis_CL_{CL_MAX_DEFAULT}_{now}'
os.makedirs(out_dir, exist_ok=True)

# ---------------------------- Compute grid --------------------------------
print("Running sweep over masses:", masses)
rows = []
for mass in masses:
    for c in chord_vals:
        c_upper = 0.25 * c  # 25% da corda total para c_upper
        c_lower = 0.75 * c  # 75% da corda total para c_lower
        for span in span_vals:
            g = Glider(mass, span, c_upper, c_lower, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            rows.append({
                "mass_kg": mass,
                "c_upper_m": c_upper,
                "c_lower_m": c_lower,
                "span_m": span,
                "total_chord_m": g.total_chord,  # Agora estamos adicionando a coluna total_chord
                "wing_area_m2": g.wing_area,
                "Vstall_m_s": g.stall_velocity(CL_max=CL_MAX_DEFAULT, rho=RHO)
            })

df = pd.DataFrame(rows)
full_csv = os.path.join(out_dir, 'grid_stall.csv')
df.to_csv(full_csv, index=False)
print("Saved stall grid CSV to:", full_csv)

# ---------------------------- Helper: bilinear interpolation --------------
def bilinear_interpolate(x_vals, y_vals, Z, x, y):
    x = float(x); y = float(y)
    if x <= x_vals[0] or x >= x_vals[-1] or y <= y_vals[0] or y >= y_vals[-1]:
        ix = int(np.abs(x_vals - x).argmin())
        iy = int(np.abs(y_vals - y).argmin())
        return float(Z[iy, ix])
    i = np.searchsorted(x_vals, x) - 1
    j = np.searchsorted(y_vals, y) - 1
    x1, x2 = x_vals[i], x_vals[i+1]
    y1, y2 = y_vals[j], y_vals[j+1]
    z11 = Z[j, i]; z21 = Z[j, i+1]; z12 = Z[j+1, i]; z22 = Z[j+1, i+1]
    tx = (x - x1) / (x2 - x1); ty = (y - y1) / (y2 - y1)
    return float((1-tx)*(1-ty)*z11 + tx*(1-ty)*z21 + (1-tx)*ty*z12 + tx*ty*z22)

# ---------------------------- Plotting (stall only) ------------------------
vstall_dir = os.path.join(out_dir, 'vstall_heatmaps')
os.makedirs(vstall_dir, exist_ok=True)

for mass, group in df.groupby('mass_kg'):
    sub = group.copy()
    pivot_vstall = sub.pivot_table(index='total_chord_m', columns='span_m', values='Vstall_m_s', aggfunc='mean')
    X = pivot_vstall.columns.values
    Y = pivot_vstall.index.values
    extent = [X.min(), X.max(), Y.min(), Y.max()]

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(pivot_vstall.values, origin='lower', aspect='auto', extent=extent, cmap=COLORMAP)
    ax.set_title(f'Vstall (m/s) — CLmax={CL_MAX_DEFAULT} — mass={mass} kg')
    ax.set_xlabel('Span (m)'); ax.set_ylabel('Chord (m)')
    c = fig.colorbar(im, ax=ax); c.set_label('Vstall (m/s)')

    x_ref, y_ref = 0.6, 0.08
    ax.axvline(x=x_ref, color='black', linestyle='--', linewidth=1)
    ax.axhline(y=y_ref, color='black', linestyle='--', linewidth=1)
    ax.plot(x_ref, y_ref, 'ro', markeredgecolor='black', markeredgewidth=1.5,
            markersize=6, label='Stall Velocity \n (b=0.6m, c=0.05m)')
    
    # Interpolando a velocidade de stall no ponto de referência
    vstall_ref = bilinear_interpolate(X, Y, pivot_vstall.values, x_ref, y_ref)
    offset_x = (X.max() - X.min()) * 0.03
    offset_y = (Y.max() - Y.min()) * 0.03
    ax.text(x_ref + offset_x, y_ref + offset_y,
            f"{vstall_ref:.2f} m/s",
            color='black', fontsize=9, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=2))

    ax.legend(loc='best', fontsize=8)

    fig.savefig(os.path.join(vstall_dir, f'heat_vstall_mass{mass:.1f}.png'), dpi=200)
    plt.close(fig)

    print(f"Saved stall heatmap for mass {mass:.1f} kg")

print("\nAll done. Check folder:", out_dir)
