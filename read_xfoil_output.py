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

a1, b1, c1, d1, e1, f1 = read_polar("Airfoils/selig_s1223/s1223_Re_5e4.txt")
a2, b2, c2, d2, e2, f2 = read_polar("Airfoils/selig_s1223/s1223_Re_1e5.txt")
a3, b3, c3, d3, e3, f3 = read_polar("Airfoils/selig_s1223/s1223_Re_2e5.txt")
print(b2)
plt.figure(figsize=(15,7))

# Título principal da janela
plt.suptitle("Lift Curve Comparison for Different Reynolds", fontsize=16)

### Plotting the Lift Curve ###
plt.subplot(1, 3, 1)
plt.ylim(0, 2.2)
plt.plot(a1, b1, "k")
plt.xlabel("Angle of Attack [deg]")
plt.ylabel("Lift Coefficient")
plt.grid(True)
plt.title("Re = 5e4")

### Plotting the Drag Curve ###
plt.subplot(1, 3, 2)
plt.ylim(0, 2.2)
plt.plot(a2, b2, "k")
plt.xlabel("Angle of Attack [deg]")
plt.ylabel("Lift Coefficient")
plt.grid(True)
plt.title("Re = 1e5")

### Plotting the Glide Efficiency ###
plt.subplot(1, 3, 3)
plt.ylim(0, 2.2)
plt.plot(a3, b3, "k")
plt.xlabel("Angle of Attack [deg]")
plt.ylabel("Lift Coefficient")
plt.grid(True)
plt.title("Re = 2e5")

plt.tight_layout(rect=[0, 0, 1, 0.95])  # deixa espaço para o título principal
plt.show()


### Plotting the Lift Curve ###
plt.plot(alfas, CLs, "k")
plt.xlabel("Angle of Atack [deg]")
plt.ylabel("Lift Coefficient")
plt.grid(True)
plt.title("Lift Curve")
plt.savefig("Lift_Curve")
plt.show()

### Plotting the Drag Curve ###
plt.plot(alfas, CDs, "k")
plt.xlabel("Angle of Atack [deg]")
plt.ylabel("Drag Coefficient")
plt.grid(True)
plt.title("Drag Curve")
plt.savefig("Drag_Curve")
plt.show()

### Plotting the GLide Efficiency Over Angle of Atack ##
plt.plot(alfas, L_D_Quocients, "k")
plt.xlabel("Angle of Atack [Degrees]")
plt.ylabel("Glide Efficiency [L/D]")
plt.grid(True)
plt.title("Glide Efficiency")
plt.savefig("Glide_Efficiency")

plt.show()