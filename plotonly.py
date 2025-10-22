import matplotlib.pyplot as plt

def read_polar(filename):
    alfas, CLs = [], []
    with open(filename, "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if i >= 12:  # ignora cabeçalho
                values = line.strip().split()
                if len(values) >= 2:
                    alfas.append(float(values[0]))
                    CLs.append(float(values[1]))
    return alfas, CLs

# ======== Lê o ficheiro (ajusta o caminho conforme necessário) ========
filename = "Airfoils/NACA0010/naca0010_Re_5e4.txt"
alfas, CLs = read_polar(filename)

# ======== Plot ========
plt.figure(figsize=(8,5))
plt.plot(alfas, CLs, 'r-', linewidth=2)

plt.title("Lift Curve, Re=50k", fontsize=14)
plt.xlabel("Angle of Attack [deg]")
plt.ylabel("Lift Coefficient (C_M)")
plt.grid(True)
plt.show()
