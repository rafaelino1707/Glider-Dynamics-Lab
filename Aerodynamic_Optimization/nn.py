# Regressão linear: V (y) vs R_kOhm (x)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

csv_txt = """Ohmimetro_de_Bancada_kOhm,Circuito_Volts,Erro_Relativo_%
0.991,0.979,1.21
2.176,2.166,0.46
4.67,4.623,1.00
6.65,6.604,0.69
9.87,9.827,0.44
14.7,14.731,0.21
17.8,16.91,5.00
"""

# Ler só kOhm e Volts
df = pd.read_csv(StringIO(csv_txt), usecols=["Ohmimetro_de_Bancada_kOhm","Circuito_Volts"])
x = df["Ohmimetro_de_Bancada_kOhm"].to_numpy()
y = df["Circuito_Volts"].to_numpy()

# Regressão y = m x + b
m, b = np.polyfit(x, y, 1)
y_fit = m * x + b

# R^2
ss_res = np.sum((y - y_fit)**2)
ss_tot = np.sum((y - np.mean(y))**2)
r2 = 1 - ss_res/ss_tot

print(f"y = {m:.6f} x + {b:.6f}")
print(f"R^2 = {r2:.6f}")

# Gráfico
plt.figure()
plt.scatter(x, y, label="Dados")
plt.plot(np.sort(x), m*np.sort(x) + b, label=f"y = {m:.4f}x + {b:.4f}")
plt.xlabel("Resistência (kΩ)")
plt.ylabel("Tensão (V)")
plt.title("Regressão Linear V vs kΩ")
plt.legend()
plt.tight_layout()
plt.show()
