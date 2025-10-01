import matplotlib.pyplot as plt

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
                L_D_Quocients.append(float(values[1]) / float(values[2]))
            i += 1
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
            i += 1
    return alfas, CLs

# ======================= LEITURA DOS DADOS =======================
# S1223
a1, b1, c1, d1, e1, f1 = read_polar("Airfoils/S-1223/s1223_Re_5e4.txt")
a2, b2, c2, d2, e2, f2 = read_polar("Airfoils/S-1223/s1223_Re_1e5.txt")
a3, b3, c3, d3, e3, f3 = read_polar("Airfoils/S-1223/s1223_Re_2e5.txt")

# MH60
aa1, bb1, cc1, dd1, ee1, ff1 = read_polar("Airfoils/mh60/mh60_Re_5e4.txt")
aa2, bb2, cc2, dd2, ee2, ff2 = read_polar("Airfoils/mh60/mh60_Re_1e5.txt")
aa3, bb3, cc3, dd3, ee3, ff3 = read_polar("Airfoils/mh60/mh60_Re_2e5.txt")

# SD7037
aaa1, bbb1, ccc1, ddd1, eee1, fff1 = read_polar("Airfoils/SD7037/sd7037_Re_5e4.txt")
aaa2, bbb2, ccc2, ddd2, eee2, fff2 = read_polar("Airfoils/SD7037/sd7037_Re_1e5.txt")
aaa3, bbb3, ccc3, ddd3, eee3, fff3 = read_polar("Airfoils/SD7037/sd7037_Re_2e5.txt")

# S3021
aaaa1, bbbb1, cccc1, dddd1, eeee1, ffff1 = read_polar("Airfoils/S-3021/s3021_Re_5e4.txt")
aaaa2, bbbb2, cccc2, dddd2, eeee2, ffff2 = read_polar("Airfoils/S-3021/s3021_Re_1e5.txt")
aaaa3, bbbb3, cccc3, dddd3, eeee3, ffff3 = read_polar("Airfoils/S-3021/s3021_Re_2e5.txt")

# Inviscid

# ======================= FUNÇÃO FULL SCREEN =======================
def fullscreen():
    try:
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
    except Exception:
        pass

# ======================= FIGURA 1: Lift Curves =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("Lift Curve Comparison for Different Reynolds", fontsize=16)
fullscreen()

colors = ["r", "g", "b", "m"]
labels = ["Selig S-1223", "MH60", "SD7037", "Selig S-3021"]

# Re = 50k
axs[0].plot(a1, b1, colors[0], label=labels[0])
axs[0].plot(aa1, bb1, colors[1], label=labels[1])
axs[0].plot(aaa1, bbb1, colors[2], label=labels[2])
axs[0].plot(aaaa1, bbbb1, colors[3], label=labels[3])
axs[0].set_ylim(0, 2.5)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("Lift Coefficient")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")
axs[0].legend()

# Re = 100k
axs[1].plot(a2, b2, colors[0], label=labels[0])
axs[1].plot(aa2, bb2, colors[1], label=labels[1])
axs[1].plot(aaa2, bbb2, colors[2], label=labels[2])
axs[1].plot(aaaa2, bbbb2, colors[3], label=labels[3])
axs[1].set_ylim(0, 2.5)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("Lift Coefficient")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")
axs[1].legend()

# Re = 200k
axs[2].plot(a3, b3, colors[0], label=labels[0])
axs[2].plot(aa3, bb3, colors[1], label=labels[1])
axs[2].plot(aaa3, bbb3, colors[2], label=labels[2])
axs[2].plot(aaaa3, bbbb3, colors[3], label=labels[3])
axs[2].set_ylim(0, 2.5)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("Lift Coefficient")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")
axs[2].legend()

fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)
plt.savefig("Lift_Curves")
plt.show()

# ======================= FIGURA 2: CL vs CD =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("CL vs CD Comparison for Different Reynolds", fontsize=16)
fullscreen()

# Re = 50k
axs[0].plot(c1, b1, "r", label=labels[0])
axs[0].plot(cc1, bb1, "g", label=labels[1])
axs[0].plot(ccc1, bbb1, "b", label=labels[2])
axs[0].plot(cccc1, bbbb1, "m", label=labels[3])
axs[0].set_ylim(0, 2.5)
axs[0].set_xlabel("Drag Coefficient")
axs[0].set_ylabel("Lift Coefficient")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")
axs[0].legend()

# Re = 100k
axs[1].plot(c2, b2, "r", label=labels[0])
axs[1].plot(cc2, bb2, "g", label=labels[1])
axs[1].plot(ccc2, bbb2, "b", label=labels[2])
axs[1].plot(cccc2, bbbb2, "m", label=labels[3])
axs[1].set_ylim(0, 2.5)
axs[1].set_xlabel("Drag Coefficient")
axs[1].set_ylabel("Lift Coefficient")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")
axs[1].legend()

# Re = 200k
axs[2].plot(c3, b3, "r", label=labels[0])
axs[2].plot(cc3, bb3, "g", label=labels[1])
axs[2].plot(ccc3, bbb3, "b", label=labels[2])
axs[2].plot(cccc3, bbbb3, "m", label=labels[3])
axs[2].set_ylim(0, 2.5)
axs[2].set_xlabel("Drag Coefficient")
axs[2].set_ylabel("Lift Coefficient")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")
axs[2].legend()

fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)
plt.savefig("CL_vs_CD")
plt.show()

# ======================= FIGURA 3: Drag Curves =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("Drag Curve Comparison for Different Reynolds", fontsize=16)
fullscreen()

# Re = 50k
axs[0].plot(a1, c1, "r", label=labels[0])
axs[0].plot(aa1, cc1, "g", label=labels[1])
axs[0].plot(aaa1, ccc1, "b", label=labels[2])
axs[0].plot(aaaa1, cccc1, "m", label=labels[3])
axs[0].set_ylim(0, 0.15)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("Drag Coefficient")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")
axs[0].legend()

# Re = 100k
axs[1].plot(a2, c2, "r", label=labels[0])
axs[1].plot(aa2, cc2, "g", label=labels[1])
axs[1].plot(aaa2, ccc2, "b", label=labels[2])
axs[1].plot(aaaa2, cccc2, "m", label=labels[3])
axs[1].set_ylim(0, 0.15)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("Drag Coefficient")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")
axs[1].legend()

# Re = 200k
axs[2].plot(a3, c3, "r", label=labels[0])
axs[2].plot(aa3, cc3, "g", label=labels[1])
axs[2].plot(aaa3, ccc3, "b", label=labels[2])
axs[2].plot(aaaa3, cccc3, "m", label=labels[3])
axs[2].set_ylim(0, 0.15)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("Drag Coefficient")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")
axs[2].legend()

fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)
plt.savefig("Drag_Curves")
plt.show()

# ======================= FIGURA 4: L/D Ratio =======================
fig, axs = plt.subplots(1, 3)
fig.suptitle("L/D Ratio Comparison for Different Reynolds", fontsize=16)
fullscreen()

# Re = 50k
axs[0].plot(a1, f1, "r", label=labels[0])
axs[0].plot(aa1, ff1, "g", label=labels[1])
axs[0].plot(aaa1, fff1, "b", label=labels[2])
axs[0].plot(aaaa1, ffff1, "m", label=labels[3])
axs[0].set_ylim(0, 100)
axs[0].set_xlabel("Angle of Attack [deg]")
axs[0].set_ylabel("L/D Ratio")
axs[0].grid(True)
axs[0].set_title("Re = 50 000")
axs[0].legend()

# Re = 100k
axs[1].plot(a2, f2, "r", label=labels[0])
axs[1].plot(aa2, ff2, "g", label=labels[1])
axs[1].plot(aaa2, fff2, "b", label=labels[2])
axs[1].plot(aaaa2, ffff2, "m", label=labels[3])
axs[1].set_ylim(0, 100)
axs[1].set_xlabel("Angle of Attack [deg]")
axs[1].set_ylabel("L/D Ratio")
axs[1].grid(True)
axs[1].set_title("Re = 100 000")
axs[1].legend()

# Re = 200k
axs[2].plot(a3, f3, "r", label=labels[0])
axs[2].plot(aa3, ff3, "g", label=labels[1])
axs[2].plot(aaa3, fff3, "b", label=labels[2])
axs[2].plot(aaaa3, ffff3, "m", label=labels[3])
axs[2].set_ylim(0, 100)
axs[2].set_xlabel("Angle of Attack [deg]")
axs[2].set_ylabel("L/D Ratio")
axs[2].grid(True)
axs[2].set_title("Re = 200 000")
axs[2].legend()

fig.subplots_adjust(left=0.07, right=0.95, wspace=0.3, top=0.88, bottom=0.1)
plt.savefig("L_D_Ratio")
plt.show()
