import matplotlib.pyplot as plt
import numpy as np



def read_polar(filename):
    alfas, CLs, CDs, CDps, CMs, L_D_Quocients = [], [], [], [], [], []
    with open(filename, "r") as file:
        i = 0
        lines = file.readlines()
        for line in lines:
            if i >= 12:
                values = line.strip().split()
                alfas.append(float(values[0]))
                CLs.append(float(values[1]))
                CDs.append(float(values[2]))
                CDps.append(float(values[3]))
                CMs.append(float(values[4]))
                L_D_Quocients.append(float(values[1])/float(values[2]))
            i+=1
    
    return alfas, CLs, CDs, CDps, CMs, L_D_Quocients

def read_polar_inviscid(filename):
    alfas, CLs = [], []
    with open(filename, "r") as file:
        i = 0
        lines = file.readlines()
        for line in lines:
            if i >= 12:
                values = line.strip().split()
                alfas.append(float(values[0]))
                CLs.append(float(values[1]))
            i+=1
    
    return alfas, CLs

a1, b1, c1, d1, e1, f1 = read_polar("Airfoils/selig_s1223/s1223_Re_5e4.txt")
a2, b2, c2, d2, e2, f2 = read_polar("Airfoils/selig_s1223/s1223_Re_1e5.txt")
a3, b3, c3, d3, e3, f3 = read_polar("Airfoils/selig_s1223/s1223_Re_2e5.txt")

aa1, bb1, cc1, dd1, ee1, ff1 = read_polar("Airfoils/selig_s3021/s3021_Re_50000")
aa2, bb2, cc2, dd2, ee2, ff2 = read_polar("Airfoils/selig_s3021/s3021_Re_100000")
aa3, bb3, cc3, dd3, ee3, ff3 = read_polar("Airfoils/selig_s3021/s3021_Re_200000")

a4, b4 = read_polar_inviscid("Airfoils/selig_s1223/s1223_inviscid.txt")

# =======================
# Figura 1: Lift Curves
# =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("Lift Curve Comparison for Different Reynolds", fontsize=16)

# Maximizar a janela automaticamente
try:
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()     # Qt5Agg
    # mng.window.state('zoomed')   # TkAgg (descomenta se precisares)
except Exception as e:
    print("Não consegui maximizar automaticamente:", e)
# Subplot 1
axs[0].plot(a1, b1, "k", label="Selig S-1223")   # linha preta contínua
axs[0].plot(aa1, bb1, "--k", label="Selig S-3021")  # linha preta tracejada
axs[0].set_ylim(0, 2.5)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("Lift Coefficient")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")
axs[0].legend()

# Subplot 2
axs[1].plot(a2, b2, "k", label="Selig S-1223")
axs[1].plot(aa2, bb2, "--k", label="Selig S-3021")
axs[1].set_ylim(0, 2.5)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("Lift Coefficient")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")
axs[1].legend()

# Subplot 3
axs[2].plot(a3, b3, "k", label="Selig S-1223")
axs[2].plot(aa3, bb3, "--k", label="Selig S-3021")
axs[2].set_ylim(0, 2.5)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("Lift Coefficient")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")
axs[2].legend()


# Ajusta margens para centrar bem
fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)

plt.savefig("Lift_Curves")
plt.show()


# =======================
# Figura 2: Drag Curves
# =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("Drag Curve Comparison for Different Reynolds", fontsize=16)

try:
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
except Exception as e:
    print("Não consegui maximizar automaticamente:", e)

# Subplot 1
axs[0].plot(a1, c1, "k")
axs[0].set_ylim(0, 0.15)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("Drag Coefficient")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")

# Subplot 2
axs[1].plot(a2, c2, "k")
axs[1].set_ylim(0, 0.15)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("Drag Coefficient")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")

# Subplot 3
axs[2].plot(a3, c3, "k")
axs[2].set_ylim(0, 0.15)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("Drag Coefficient")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")

# Ajusta margens
fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)

plt.savefig("Drag_Curves")
plt.show()




# =======================
# Figura 2: Drag Curves
# =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("L/D Ratio Comparison for Different Reynolds", fontsize=16)

try:
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
except Exception as e:
    print("Não consegui maximizar automaticamente:", e)

# Subplot 1
axs[0].plot(a1, f1, "k")
axs[0].set_ylim(0, 100)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("L/D Ratio")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")

# Subplot 2
axs[1].plot(a2, f2, "k")
axs[1].set_ylim(0, 100)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("L/D Ratio")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")

# Subplot 3
axs[2].plot(a3, f3, "k")
axs[2].set_ylim(0, 100)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("L/D Ratio")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")

# Ajusta margens
fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)

plt.savefig("L_D_Ratio")
plt.show()

plt.plot(a4, b4, "k")
plt.ylim(0, 4)
plt.ylabel("Lift Coefficient")
plt.xlabel("Angle of Atack [deg]")
plt.grid(True)
plt.title("Lift Curve for Inviscid Model", pad=40)  # <-- sobe o título
plt.show()