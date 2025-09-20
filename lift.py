"""
AeroSandbox glider script using a local .dat airfoil (e.g. s3021.dat)
Requirements:
    pip install aerosandbox numpy matplotlib
Place your local airfoil file in the same folder and name it `s3021.dat`
(or change FILENAME below).
"""

import aerosandbox as asb
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --------- CONFIG ----------
FILENAME = "s3021.dat"   # <- change if your dat file has a different name
S_REF = 8.0
SPAN = 10.0
C_REF = 1.0
# --------------------------

def load_dat_coords(path: Path):
    """
    Tries to read Selig/.dat coordinates from a file robustly.
    Returns Nx2 numpy array or raises ValueError.
    """
    raw = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    # Filter lines to those that look like numbers (x y)
    coords = []
    for ln in raw:
        ln = ln.strip()
        if ln == "":
            continue
        # skip obvious header lines (non-numeric start)
        parts = ln.split()
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0])
            y = float(parts[1])
            coords.append([x, y])
        except ValueError:
            # maybe the file has a header line; skip it
            continue
    if len(coords) < 6:
        raise ValueError(f"Too few coordinate lines parsed from {path} ({len(coords)}).")
    return np.array(coords)

# -------------------------
# 1) Try load the local .dat
# -------------------------
af_obj = None
dat_path = Path(FILENAME)

# If load failed, try the internal helper as fallback (optional)
if af_obj is None:
    try:
        from aerosandbox.geometry.airfoil import airfoil_families as af_fam
        # try some common name variants
        for candidate in ["s3021-il", "s3021", "S3021"]:
            try:
                coords = af_fam.get_UIUC_coordinates(candidate)
                if coords is not None and len(coords) > 0:
                    af_obj = asb.Airfoil(coordinates=coords)
                    print(f"[OK] Loaded {candidate} from internal UIUC helper.")
                    break
            except Exception:
                pass
    except Exception:
        pass


# -------------------------
# 2) Define the glider (use the loaded airfoil for the main wing sections)
# -------------------------
glider = asb.Airplane(
    name="Simple Glider",
    xyz_ref=[0, 0, 0],
    s_ref=S_REF,
    b_ref=SPAN,
    c_ref=C_REF,
    wings=[
        asb.Wing(
            name="Main Wing",
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1.2,
                    twist=0.0,
                    airfoil=af_obj,   # <- use loaded airfoil here
                ),
                asb.WingXSec(
                    xyz_le=[3.0, 5.0, 0],
                    chord=0.6,
                    twist=0.0,
                    airfoil=af_obj,   # same airfoil at tip
                ),
            ],
        ),
        asb.Wing(
            name="Horizontal Stabilizer",
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[4.5, 0, -0.8],
                    chord=0.6,
                    twist=0.0,
                    airfoil=asb.Airfoil("naca0012"),
                ),
                asb.WingXSec(
                    xyz_le=[5.5, 1.5, -0.8],
                    chord=0.4,
                    twist=0.0,
                    airfoil=asb.Airfoil("naca0012"),
                ),
            ],
        ),
    ],
)

# Mass (optional)
mass_kg = 2.0
g = 9.80665
weight_N = mass_kg * g

# -------------------------
# 3) Sweep alpha & run lifting line
# -------------------------
alphas_deg = np.linspace(-2, 10, 20)

CL_list = []
CD_list = []
L_over_D_list = []

atm = asb.Atmosphere(altitude=0.0)

for alpha in alphas_deg:
    op_point = asb.OperatingPoint(
        atmosphere=atm,
        velocity=12.0,
        alpha=alpha,
        beta=0.0,
        p=0.0,
        q=0.0,
        r=0.0
    )

    aero = asb.LiftingLine(
        airplane=glider,
        op_point=op_point,
        spanwise_resolution=20
    ).run()

    CL = aero.get("CL", np.nan)

    if "CD" in aero and np.isfinite(aero["CD"]):
        CD = aero["CD"]
    else:
        CD0 = 0.02
        AR = aero.get("AR", (SPAN**2 / S_REF))
        e = 0.9
        k = 1.0 / (np.pi * AR * e)
        CD = CD0 + k * CL**2

    CL_list.append(CL)
    CD_list.append(CD)
    L_over_D_list.append(CL / CD if CD > 0 else np.nan)

# -------------------------
# 4) Plots
# -------------------------
plt.figure(figsize=(11, 6))

plt.subplot(1, 3, 1)
plt.plot(alphas_deg, CL_list, "g")
plt.xlabel("Angle of attack (deg)")
plt.ylabel("C_L")
plt.grid(True)
plt.title("Lift curve")

plt.subplot(1, 3, 2)
plt.plot(CL_list, CD_list, "b")
plt.xlabel("C_L")
plt.ylabel("C_D")
plt.grid(True)
plt.title("Drag polar")

plt.subplot(1, 3, 3)
plt.plot(alphas_deg, L_over_D_list, "r")
plt.xlabel("Angle of attack (deg)")
plt.ylabel("L/D")
plt.grid(True)
plt.title("Glide efficiency")

plt.tight_layout()
plt.show()

best_idx = int(np.nanargmax(L_over_D_list))
print(f"Best L/D ≈ {L_over_D_list[best_idx]:.2f} at α = {alphas_deg[best_idx]:.1f}°")
