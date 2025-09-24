import matplotlib.pyplot as plt
import numpy as np

alfas, CLs, CDs, CDps, CMs, L_D_Quocients = [], [], [], [], [], []

def read_polar(filename):
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

read_polar("s1223_sim1.txt")
### Plotting the Lift Curve ###
plt.figure(figsize=(11,8))
plt.plot(alfas, CLs, "k")
plt.xlabel("Angle of Atack [deg]")
plt.ylabel("Lift Coefficient")
plt.grid(True)
plt.title("Lift Curve")
plt.savefig("Lift_Curve")

### Plotting the Drag Curve ###
plt.figure(figsize=(11,8))
plt.plot(alfas, CDs, "k")
plt.xlabel("Angle of Atack [deg]")
plt.ylabel("Drag Coefficient")
plt.grid(True)
plt.title("Drag Curve")
plt.savefig("Drag_Curve")

### Plotting the GLide Efficiency Over Angle of Atack ###
plt.figure(figsize=(11,8))
plt.plot(alfas, L_D_Quocients, "k")
plt.xlabel("Angle of Atack [Degrees]")
plt.ylabel("Glide Efficiency [L/D]")
plt.grid(True)
plt.title("Glide Efficiency")
plt.savefig("Glide_Efficiency")

plt.show()