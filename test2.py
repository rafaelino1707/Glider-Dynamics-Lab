"""
Glider optimization / grid sweep script (updated)

Este script varre PLA_total de 1.0 kg até 0.1 kg com passo 0.1 kg e gera:
 - heatmaps para sink_rate e Vmin
 - mapas combinados (objetivo normalizado)
 - scatter Vmin vs sink_rate com fronteira de Pareto
 - fatiamentos (slices) de Vmin e sink_rate para cordas selecionadas

Parâmetros no bloco SETTINGS. Use num Jupyter Notebook ou execute como script.
"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

class Glider:
    def __init__(self, mass_kg, chord_m, wingspan_m, oswald_efficiency=0.75, cd0=0.04):
        self.mass = mass_kg
        self.chord = chord_m
        self.wingspan = wingspan_m
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
            V = self.minimum_sink_velocity(rho=rho)
        CL = W / (0.5 * rho * V ** 2 * S)
        K = self.induced_drag_factor
        CD = self.cd0 + K * CL ** 2
        D = 0.5 * rho * V ** 2 * S * CD
        sink_rate = V * (D / W)
        glide_ratio = (W / D) if D > 0 else float('inf')
        return {
            'V_used': V,
            'CL': CL,
            'CD': CD,
            'D_N': D,
            'sink_rate_m_s': sink_rate,
            'glide_ratio': glide_ratio
        }


def wing_mass_from_shell(chord_m, span_m, thickness_m=0.002, rho_pla=1240.0):
    wing_area = chord_m * span_m
    total_shell_area = 2.0 * wing_area
    wing_volume = total_shell_area * thickness_m
    return wing_volume * rho_pla

# ------------------------- SETTINGS ---------------------------------------
chord_vals = np.linspace(0.02, 0.10, 41)   # corda: 2 cm a 10 cm
span_vals = np.linspace(0.45, 0.85, 41)    # span: 0.45 a 0.85 m

# PLA_total de 1.0 kg até 0.1 kg com passo 0.1 kg
PLA_total_list = np.round(np.arange(1.0, 0.0, -0.1), 2)
PLA_total_list = [p for p in PLA_total_list if p >= 0.1]

reserved_mass_kg = 0.36  # reserva fixa (ajusta se quiseres)

rho_PLA = 1240.0
thickness_m = 0.002
oswald_efficiency = 0.75
cd0 = 0.04
rho_air = 1.225

lambda_weight = 0.5  # peso relativo entre Vmin_norm e sink_norm

# Dir de saída
now = datetime.now().strftime('%Y%m%d_%H%M%S')
out_dir = f'glider_opt_outputs_{now}'
os.makedirs(out_dir, exist_ok=True)

# ----------------------- GRID SWEEP --------------------------------------
results = []
for PLA_total in PLA_total_list:
    for chord in chord_vals:
        for span in span_vals:
            wing_mass = wing_mass_from_shell(chord, span, thickness_m=thickness_m, rho_pla=rho_PLA)
            feasible = (wing_mass + reserved_mass_kg) <= PLA_total
            total_mass = wing_mass + reserved_mass_kg

            if feasible:
                g = Glider(total_mass, chord, span, oswald_efficiency=oswald_efficiency, cd0=cd0)
                Vmin = g.minimum_sink_velocity(rho=rho_air)
                sink = g.vertical_sink_rate(V=Vmin, rho=rho_air)
                sink_rate = sink['sink_rate_m_s']
                glide_ratio = sink['glide_ratio']
            else:
                Vmin = np.nan
                sink_rate = np.nan
                glide_ratio = np.nan

            results.append({
                'PLA_total_kg': PLA_total,
                'reserved_mass_kg': reserved_mass_kg,
 'chord_m': chord,
                'span_m': span,
                'wing_area_m2': chord*span,
                'wing_mass_kg': wing_mass,
                'feasible': feasible,
                'total_mass_kg': total_mass,
                'Vmin_m_s': Vmin,
                'sink_rate_m_s': sink_rate,
                'glide_ratio': glide_ratio
            })

df = pd.DataFrame(results)
df_feas = df[df['feasible']].copy()

# Normalizar por grupo PLA_total
df_feas['Vmin_norm'] = df_feas.groupby(['PLA_total_kg'])['Vmin_m_s'].transform(lambda x: (x - x.min())/(x.max()-x.min()) if x.max()!=x.min() else 0.0)
df_feas['sink_norm'] = df_feas.groupby(['PLA_total_kg'])['sink_rate_m_s'].transform(lambda x: (x - x.min())/(x.max()-x.min()) if x.max()!=x.min() else 0.0)
df_feas['combined_obj'] = lambda_weight * df_feas['Vmin_norm'] + (1-lambda_weight) * df_feas['sink_norm']

# Pareto
def pareto_front(df_sub):
    pts = df_sub[['Vmin_m_s','sink_rate_m_s']].dropna().values
    idxs = []
    for i,p in enumerate(pts):
        dominated = False
        for j,q in enumerate(pts):
            if j==i: continue
            if (q[0] <= p[0] and q[1] <= p[1]) and (q[0] < p[0] or q[1] < p[1]):
                dominated = True
                break
        if not dominated:
            idxs.append(i)
    df_nonan = df_sub[['Vmin_m_s','sink_rate_m_s']].dropna().reset_index()
    pareto_indexes = df_nonan.loc[idxs,'index'].values
    return pareto_indexes

# ------------------- Outputs: heatmaps + plots ----------------------------
summary_rows = []
for PLA_total, group in df_feas.groupby('PLA_total_kg'):
    sub = group.copy()
    pivot_sink = sub.pivot_table(index='chord_m', columns='span_m', values='sink_rate_m_s', aggfunc='mean')
    pivot_vmin = sub.pivot_table(index='chord_m', columns='span_m', values='Vmin_m_s', aggfunc='mean')
    pivot_comb = sub.pivot_table(index='chord_m', columns='span_m', values='combined_obj', aggfunc='mean')

    # Heatmaps
    plt.figure(figsize=(7,5))
    plt.imshow(pivot_sink.values, origin='lower', aspect='auto', extent=[pivot_sink.columns.min(), pivot_sink.columns.max(), pivot_sink.index.min(), pivot_sink.index.max()])
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'Sink rate (m/s) — PLA={PLA_total}kg')
    plt.colorbar().set_label('sink rate (m/s)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'heat_sink_PLA{PLA_total}.png'), dpi=200)
    plt.close()

    plt.figure(figsize=(7,5))
    plt.imshow(pivot_vmin.values, origin='lower', aspect='auto', extent=[pivot_vmin.columns.min(), pivot_vmin.columns.max(), pivot_vmin.index.min(), pivot_vmin.index.max()])
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'Vmin (m/s) — PLA={PLA_total}kg')
    plt.colorbar().set_label('Vmin (m/s)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'heat_vmin_PLA{PLA_total}.png'), dpi=200)
    plt.close()

    plt.figure(figsize=(7,5))
    plt.imshow(pivot_comb.values, origin='lower', aspect='auto', extent=[pivot_comb.columns.min(), pivot_comb.columns.max(), pivot_comb.index.min(), pivot_comb.index.max()])
    plt.xlabel('Span (m)'); plt.ylabel('Chord (m)')
    plt.title(f'Combined obj (lambda={lambda_weight}) — PLA={PLA_total}kg')
    plt.colorbar().set_label('combined objective (norm)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'heat_comb_PLA{PLA_total}.png'), dpi=200)
    plt.close()

    # Pareto scatter
    plt.figure(figsize=(6,5))
    plt.scatter(sub['Vmin_m_s'], sub['sink_rate_m_s'], s=10)
    pareto_idx = pareto_front(sub)
    if len(pareto_idx) > 0:
        pareto = sub.loc[pareto_idx]
        plt.scatter(pareto['Vmin_m_s'], pareto['sink_rate_m_s'], s=40)
        pareto_sorted = pareto.sort_values('Vmin_m_s')
        plt.plot(pareto_sorted['Vmin_m_s'], pareto_sorted['sink_rate_m_s'])
    plt.xlabel('Vmin (m/s)'); plt.ylabel('sink_rate (m/s)')
    plt.title(f'Pareto front — PLA={PLA_total}kg')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'pareto_PLA{PLA_total}.png'), dpi=200)
    plt.close()

    # Slices: para algumas cordas de interesse, plot Vmin & sink vs span
    chord_slices = [0.025, 0.05, 0.075, 0.1]
    for chord_val in chord_slices:
        slice_df = sub[np.isclose(sub['chord_m'], chord_val)]
        if slice_df.empty:
            continue
        slice_sorted = slice_df.sort_values('span_m')
        plt.figure()
        plt.plot(slice_sorted['span_m'], slice_sorted['Vmin_m_s'], marker='x')
        plt.xlabel('Span (m)'); plt.ylabel('Vmin (m/s)')
        plt.title(f'Vmin vs span — chord={chord_val} m PLA={PLA_total}kg')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'Vmin_vs_span_chord{int(chord_val*1000)}_PLA{PLA_total}.png'), dpi=200)
        plt.close()

        plt.figure()
        plt.plot(slice_sorted['span_m'], slice_sorted['sink_rate_m_s'], marker='x')
        plt.xlabel('Span (m)'); plt.ylabel('sink_rate (m/s)')
        plt.title(f'sink_rate vs span — chord={chord_val} m PLA={PLA_total}kg')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'sink_vs_span_chord{int(chord_val*1000)}_PLA{PLA_total}.png'), dpi=200)
        plt.close()

    # Top combos
    top_comb = sub.nsmallest(12, 'combined_obj')[['chord_m','span_m','wing_mass_kg','total_mass_kg','Vmin_m_s','sink_rate_m_s','glide_ratio','combined_obj']]
    top_comb['PLA_total_kg'] = PLA_total
    summary_rows.append(top_comb)
    print(f'Outputs OK for PLA_total = {PLA_total} kg')

summary_df = pd.concat(summary_rows, ignore_index=True)
summary_df_sorted = summary_df.sort_values(['combined_obj','PLA_total_kg']).reset_index(drop=True)
summary_csv = os.path.join(out_dir, f'glider_opt_summary_{now}.csv')
summary_df_sorted.to_csv(summary_csv, index=False)

print('\nDone. Files saved in:', out_dir)
print('Summary CSV:', summary_csv)
