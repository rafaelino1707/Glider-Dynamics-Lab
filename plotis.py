"""
Glider heatmap generator (chord × span × metrics)

Este script varre corda × envergadura para uma lista de massas (0.1 a 1.0 kg passo 0.1 kg)
E gera heatmaps (PNG) para cada massa com as seguintes grandezas:
 - Vmin_teorico (velocidade que minimiza o sink segundo o teu modelo)
 - sink_rate calculado em Vmin_teorico
 - Vstall calculado a partir de CL_max (agora fornecido pela classe como self.CLmax)

Notas importantes:
 - NÃO estou a fazer verificação de estol vs Vmin aqui — apenas calculo e ploto as três quantidades
 - O CL_max usado deve vir da tua curva CL vs AoA (coloca o valor real em CL_MAX_DEFAULT)
 - Ajusta ranges (`chord_vals`, `span_vals`, `masses`) conforme precisares

Uso: salvar como `glider_heatmaps.py` e correr num ambiente com numpy/pandas/matplotlib

"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ---------------------- Classe com CLmax ----------------------------------
class Glider:
    def __init__(self, mass_kg, chord_m, wingspan_m, CLmax, oswald_efficiency=0.75, cd0=0.04):
        self.mass = mass_kg
        self.chord = chord_m
        self.wingspan = wingspan_m
        self.CLmax = CLmax
        self.oswald_efficiency = oswald_efficiency
        self.cd0 = cd0

    @property
    def wing_area(self):
        return max(self.chord * self.wingspan, 1e-12)

    @property
    def aspect_ratio(self):
        return (self.wingspan ** 2) / self.wing_area

    @property
    def induced_drag_factor(self):
        return 1.0 / (math.pi * self.aspect_ratio * self.oswald_efficiency)

    def Vmin_theoretical(self, rho=1.225):
        """Velocidade que minimiza o sink em teoria (não leva em conta CLmax)."""
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
        glide = W / D if D>0 else float('nan')
        return {'V': V, 'CL': CL, 'CD': CD, 'D': D, 'sink': sink, 'glide': glide}

    def Vstall(self, rho=1.225):
        """Velocidade de estol baseada em CLmax: Vs = sqrt(2W/(rho S CLmax))"""
        W = self.mass * 9.81
        S = self.wing_area
        return math.sqrt((2.0 * W) / (rho * S * self.CLmax))

# --------------------------- SETTINGS ------------------------------------
# grid geometry
chord_vals = np.linspace(0.02, 0.10, 41)    # corda 2cm -> 10cm
span_vals = np.linspace(0.45, 0.85, 41)     # span 0.45 -> 0.85 m

# masses to sweep (0.1 -> 1.0 kg step 0.1)
masses = np.round(np.arange(0.1, 1.01, 0.1), 2)

# aerodata - ajusta CL_MAX com o valor da tua curva
CL_MAX_DEFAULT = 1.2   # substitui pelo CLmax real do teu perfil
OSWALD_EFF = 0.75
CD0 = 0.04
RHO = 1.225

# output folder
now = datetime.now().strftime('%Y%m%d_%H%M%S')
out_dir = f'glider_heatmaps_{now}'
os.makedirs(out_dir, exist_ok=True)

# dataframe results
rows = []

print('Running sweep for masses:', masses)
for mass in masses:
    print(f'  mass = {mass} kg')
    for chord in chord_vals:
        for span in span_vals:
            g = Glider(mass, chord, span, CLmax=CL_MAX_DEFAULT, oswald_efficiency=OSWALD_EFF, cd0=CD0)
            S = g.wing_area
            Vmin = g.Vmin_theoretical(rho=RHO)
            sink_info = g.sink_at_V(Vmin, rho=RHO)
            Vstall = g.Vstall(rho=RHO)

            rows.append({
                'mass_kg': mass,
                'chord_m': chord,
                'span_m': span,
                'wing_area_m2': S,
                'Vmin_theo_m_s': Vmin,
                'sink_rate_m_s': sink_info['sink'],
                'CL_at_Vmin': sink_info['CL'],
                'CD_at_Vmin': sink_info['CD'],
                'glide_at_Vmin': sink_info['glide'],
                'Vstall_m_s': Vstall
            })

# create dataframe and save
df = pd.DataFrame(rows)
csv_path = os.path.join(out_dir, 'heatmap_grid_full.csv')
df.to_csv(csv_path, index=False)
print('Saved full grid CSV to', csv_path)

# --------------------------- PLOTTING -----------------------------------
# For each mass, create 3 heatmaps: sink_rate, Vmin, Vstall (chord vs span)
for mass, group in df.groupby('mass_kg'):
    sub = group.copy()
    pivot_sink = sub.pivot_table(index='chord_m', columns='span_m', values='sink_rate_m_s', aggfunc='mean')
    pivot_vmin = sub.pivot_table(index='chord_m', columns='span_m', values='Vmin_theo_m_s', aggfunc='mean')
    pivot_vstall = sub.pivot_table(index='chord_m', columns='span_m', values='Vstall_m_s', aggfunc='mean')

    # common extent for plots
    extent = [pivot_sink.columns.min(), pivot_sink.columns.max(), pivot_sink.index.min(), pivot_sink.index.max()]

    # sink heatmap
    plt.figure(figsize=(7,5))
    plt.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=extent)
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'sink_rate (m/s) at Vmin — mass={mass} kg')
    c = plt.colorbar(); c.set_label('sink_rate (m/s)')
    try:
        X = pivot_sink.columns.values; Y = pivot_sink.index.values; Z = pivot_sink.values
        cs = plt.contour(X, Y, Z, colors='k', linewidths=0.6)
        plt.clabel(cs, inline=True, fontsize=8)
    except Exception:
        pass
    plt.tight_layout()
    fname1 = os.path.join(out_dir, f'heat_sink_mass{mass}.png')
    plt.savefig(fname1, dpi=200); plt.close()

    # Vmin heatmap
    plt.figure(figsize=(7,5))
    plt.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=extent)
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'Vmin_theoretical (m/s) — mass={mass} kg')
    c = plt.colorbar(); c.set_label('Vmin (m/s)')
    try:
        X = pivot_vmin.columns.values; Y = pivot_vmin.index.values; Z = pivot_vmin.values
        cs = plt.contour(X, Y, Z, colors='k', linewidths=0.6)
        plt.clabel(cs, inline=True, fontsize=8)
    except Exception:
        pass
    plt.tight_layout()
    fname2 = os.path.join(out_dir, f'heat_vmin_mass{mass}.png')
    plt.savefig(fname2, dpi=200); plt.close()

    # Vstall heatmap
    plt.figure(figsize=(7,5))
    plt.imshow(pivot_vstall.values, origin='lower', aspect='auto', extent=extent)
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'Vstall (m/s) from CLmax={CL_MAX_DEFAULT} — mass={mass} kg')
    c = plt.colorbar(); c.set_label('Vstall (m/s)')
    try:
        X = pivot_vstall.columns.values; Y = pivot_vstall.index.values; Z = pivot_vstall.values
        cs = plt.contour(X, Y, Z, colors='k', linewidths=0.6)
        plt.clabel(cs, inline=True, fontsize=8)
    except Exception:
        pass
    plt.tight_layout()
    fname3 = os.path.join(out_dir, f'heat_vstall_mass{mass}.png')
    plt.savefig(fname3, dpi=200); plt.close()

    print(f'Wrote heatmaps for mass {mass} kg:')
    print('  ', fname1)
    print('  ', fname2)
    print('  ', fname3)

