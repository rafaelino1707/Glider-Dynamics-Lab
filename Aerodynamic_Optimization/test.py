"""
Glider optimization / grid sweep script

Escolhe ranges de chord/span e várias quantidades de PLA (ou massas reservadas) e
calcula Vmin, sink rate e glide ratio para cada combinação.

Requisitos: Python 3.8+, numpy, pandas, matplotlib

Uso:
  - Abra este ficheiro num Jupyter Notebook ou execute como script:
      python glider_optimization.py
  - Edite os parâmetros no bloco SETTINGS para explorar outras ranges.

Saídas:
  - CSV com todas as combinações e métricas
  - Heatmaps (PNG) para cada valor de PLA_total definido
  - Console com top combos por menor sink rate e menor Vmin

Notas:
  - Modelo simplificado (CD = CD0 + K*CL^2, K = 1/(pi*AR*e)).
  - Wing mass estimada como duas cascas (top+bottom) com espessura thickness_m.
  - Ajusta "reserved_mass_kg" para representar fuselagem + eletrónica.

"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

# -------------------------- CLASSES / FUNÇÕES -------------------------------
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
    """Calcula massa total das duas asas assumindo duas cascas (top+bottom)"""
    wing_area = chord_m * span_m
    total_shell_area = 2.0 * wing_area
    wing_volume = total_shell_area * thickness_m
    return wing_volume * rho_pla


# ----------------------------- SETTINGS ------------------------------------
# Ajusta aqui as ranges e parametros
OUTPUT_DIR = 'glider_results'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Ranges
chord_vals = np.linspace(0.02, 0.10, 41)   # corda: 2 cm a 10 cm
span_vals = np.linspace(0.45, 0.85, 41)    # span: 0.45 m a 0.85 m

# Experimenta vários totais de PLA (kg) — isto simula diferentes quantidades de material
PLA_total_list = [0.8, 1.0, 1.2]  # kg

# Reserva para fuselagem + eletrónica (kg). Podes também passar uma lista se quiseres varrer.
reserved_mass_kg = 0.36  # 360 g

# Materiais e outros parâmetros
rho_PLA = 1240.0  # kg/m^3
thickness_m = 0.002  # 2 mm
oswald_efficiency = 0.75
cd0 = 0.04
rho_air = 1.225

# ---------------------------------------------------------------------------
all_results = []

for PLA_total in PLA_total_list:
    print('\n=== Running sweep for PLA_total =', PLA_total, 'kg ===')

    for chord in chord_vals:
        for span in span_vals:
            wing_mass = wing_mass_from_shell(chord, span, thickness_m=thickness_m, rho_pla=rho_PLA)
            feasible = (wing_mass + reserved_mass_kg) <= PLA_total
            total_mass = wing_mass + reserved_mass_kg

            if feasible:
                g = Glider(total_mass, chord, span, oswald_efficiency=oswald_efficiency, cd0=cd0)
                try:
                    Vmin = g.minimum_sink_velocity(rho=rho_air)
                    sink = g.vertical_sink_rate(V=Vmin, rho=rho_air)
                    sink_rate = sink['sink_rate_m_s']
                    glide_ratio = sink['glide_ratio']
                except Exception:
                    Vmin = np.nan
                    sink_rate = np.nan
                    glide_ratio = np.nan
            else:
                Vmin = np.nan
                sink_rate = np.nan
                glide_ratio = np.nan

            all_results.append({
                'PLA_total_kg': PLA_total,
                'reserved_mass_kg': reserved_mass_kg,
                'chord_m': chord,
                'span_m': span,
                'wing_area_m2': chord * span,
                'wing_mass_kg': wing_mass,
                'feasible': feasible,
                'total_mass_kg': total_mass,
                'Vmin_m_s': Vmin,
                'sink_rate_m_s': sink_rate,
                'glide_ratio': glide_ratio
            })

# Save dataframe
df = pd.DataFrame(all_results)
now = datetime.now().strftime('%Y%m%d_%H%M%S')
csv_path = os.path.join(OUTPUT_DIR, f'glider_sweep_{now}.csv')
df.to_csv(csv_path, index=False)
print('\nSaved CSV with all results to:', csv_path)

# For cada PLA_total, criar heatmap e listar top combos
for PLA_total in sorted(df['PLA_total_kg'].unique()):
    sub = df[df['PLA_total_kg'] == PLA_total].copy()
    feasible = sub[sub['feasible']]

    if feasible.empty:
        print('No feasible combos for PLA_total =', PLA_total)
        continue

    # Pivot para heatmap
    pivot = feasible.pivot_table(index='chord_m', columns='span_m', values='sink_rate_m_s', aggfunc='mean')

    plt.figure(figsize=(8,6))
    plt.imshow(pivot.values, origin='lower', aspect='auto',
               extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()])
    plt.xlabel('Span (m)')
    plt.ylabel('Chord (m)')
    plt.title(f'Sink rate mínimo (m/s) — PLA_total={PLA_total} kg')
    cbar = plt.colorbar()
    cbar.set_label('sink rate (m/s)')
    plt.tight_layout()
    png_path = os.path.join(OUTPUT_DIR, f'heatmap_sinkrate_PLA{PLA_total}_{now}.png')
    plt.savefig(png_path, dpi=200)
    plt.close()
    print('Saved heatmap to', png_path)

    # Top combos
    top_sink = feasible.nsmallest(10, 'sink_rate_m_s')
    top_vmin = feasible.nsmallest(10, 'Vmin_m_s')
    print('\nPLA_total =', PLA_total)
    print('Top 5 (menor sink rate):')
    print(top_sink[['chord_m','span_m','wing_mass_kg','total_mass_kg','Vmin_m_s','sink_rate_m_s','glide_ratio']].head(5).to_string(index=False))

    print('\nTop 5 (menor Vmin):')
    print(top_vmin[['chord_m','span_m','wing_mass_kg','total_mass_kg','Vmin_m_s','sink_rate_m_s','glide_ratio']].head(5).to_string(index=False))

print('\nDone. Checa a pasta', OUTPUT_DIR, 'para CSVs e imagens.')

# ---------------------------- Sugestões de exploração -----------------------
# - Ajusta reserved_mass_kg para refletir melhor a tua eletrónica/bateria.
# - Experimenta outras values para oswald_efficiency e cd0 (melhor perfil/forma reduz CD0).
# - Se quiseres um otimizador (ex: minimizar sink_rate + lambda*Vmin) posso adicionar
#   uma rotina usando differential evolution (scipy) ou um grid search mais fino sobre
#   uma região de interesse.
# ---------------------------------------------------------------------------
