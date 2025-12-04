#!/usr/bin/env python3
"""
Iterative minimum-sink estimator for a glider.
- Optionally calls XFOIL (if available in PATH) to re-evaluate CDmin at updated Re.
- If XFOIL is not available, uses a user-provided fallback CD0_2D.
"""
import math
import shutil
import subprocess
import tempfile
import os
import re
from typing import Optional, Tuple, Dict, List

NU_AIR = 1.4607e-5  # m^2/s (kinematic viscosity, ~15°C)


class GliderIter:
    def __init__(self,
                 mass: float,
                 wing_area: Optional[float] = None,
                 chord: Optional[float] = None,
                 wingspan: Optional[float] = None,
                 oswald_eff: float = 0.75,
                 cd0_2d_fallback: float = 0.02817,
                 f3d: float = 1.10,
                 rho: float = 1.225):
        """
        Provide either wing_area OR (chord and wingspan).
        cd0_2d_fallback: initial 2D CD0 from XFOIL (if XFOIL not available, used directly)
        f3d: 3D correction factor (e.g. 1.10)
        """
        self.mass = mass
        self.g = 9.81
        self.cd0_2d_fallback = cd0_2d_fallback
        self.f3d = f3d
        self.oswald_eff = oswald_eff
        self.rho = rho

        if wing_area is None:
            if chord is None or wingspan is None:
                raise ValueError("Provide either wing_area or (chord and wingspan).")
            self.chord = chord
            self.wingspan = wingspan
            self.wing_area = chord * wingspan
        else:
            self.wing_area = wing_area
            self.chord = chord
            self.wingspan = wingspan

        if self.chord is None and self.wingspan is not None:
            self.chord = self.wing_area / self.wingspan
        if self.wingspan is None and self.chord is not None:
            self.wingspan = self.wing_area / self.chord

    @property
    def weight(self):
        return self.mass * self.g

    @property
    def aspect_ratio(self):
        return self.wingspan ** 2 / self.wing_area

    @property
    def induced_drag_factor(self):
        return 1.0 / (math.pi * self.aspect_ratio * self.oswald_eff)

    def _call_xfoil_get_cdmin(self, airfoil_name: str, Re: float,
                              alpha_range: Tuple[float, float] = (-5.0, 15.0),
                              alpha_step: float = 0.25) -> Optional[Tuple[float, float, float]]:
        """
        Versão robusta: se 'airfoil_name' for relativo, tenta encontrá-lo no cwd e na pasta do xfoil.exe,
        depois copia o .dat para o tempdir e roda XFOIL.
        """
        xfoil_path_override = r"C:\Users\rafae\Desktop\Propulsion and AerodynamicSoftwares\XFOIL\xfoil.exe"

        xfoil_path = shutil.which("xfoil")
        if xfoil_path is None:
            if os.path.exists(xfoil_path_override):
                xfoil_path = xfoil_path_override
            else:
                print(f"[DEBUG] xfoil not found in PATH and override not found at: {xfoil_path_override}")
                return None
        else:
            print(f"[DEBUG] xfoil found in PATH: {xfoil_path}")

        print(f"[DEBUG] using xfoil executable: {xfoil_path}")

        # Resolve possible locations for the .dat: as given, cwd, and folder of xfoil.exe
        candidates = []
        if os.path.isabs(airfoil_name):
            candidates.append(airfoil_name)
        else:
            candidates.append(os.path.join(os.getcwd(), airfoil_name))
            # try same folder as xfoil.exe
            xfoil_dir = os.path.dirname(xfoil_path)
            candidates.append(os.path.join(xfoil_dir, airfoil_name))
            # also try the folder you previously mentioned (safe)
            candidates.append(r"C:\Users\rafae\Desktop\Propulsion and AerodynamicSoftwares\XFOIL\{}".format(airfoil_name))

        real_path = None
        for c in candidates:
            if os.path.exists(c):
                real_path = c
                break

        if real_path is None:
            print(f"[DEBUG] airfoil .dat not found in candidates: {candidates}")
            return None

        with tempfile.TemporaryDirectory() as td:
            base_name = os.path.basename(real_path)
            tmp_airfoil = os.path.join(td, base_name)
            try:
                shutil.copyfile(real_path, tmp_airfoil)
            except Exception as e:
                print(f"[DEBUG] failed to copy airfoil file to tempdir: {e}")
                return None

            polar_file = os.path.join(td, "polar.dat")
            cmds = []
            cmds.append(f"LOAD {base_name}\n")
            cmds.append("PPAR\n")
            cmds.append("N 450\n")
            cmds.append(" \n")
            cmds.append("OPER\n")
            cmds.append(f"VISC {Re:.0f}\n")
            cmds.append("M 0.03\n")
            cmds.append("ITER 100\n")
            cmds.append("VPAR\n")
            cmds.append("N 6\n")
            cmds.append(" \n")
            cmds.append("PACC\n")
            cmds.append("polar.dat\n\n")
            cmds.append(f"ASEQ {alpha_range[0]:.3f} {alpha_range[1]:.3f} {alpha_step:.3f}\n")
            cmds.append("PACC\n")
            cmds.append("QUIT\n")
            fullcmd = "".join(cmds)

            try:
                proc = subprocess.run([xfoil_path], input=fullcmd.encode('utf-8'),
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=td, timeout=40)
            except Exception as e:
                print(f"[DEBUG] xfoil run exception: {e}")
                return None

            stdout = proc.stdout.decode('utf-8', errors='replace')
            stderr = proc.stderr.decode('utf-8', errors='replace')
            print(f"[DEBUG] xfoil stdout (truncated):\n{stdout[:1200]}\n---\n[DEBUG] xfoil stderr (truncated):\n{stderr[:1200]}\n---")
            print(f"[DEBUG] polar file exists? {os.path.exists(polar_file)} at {polar_file}")

            if os.path.exists(polar_file):
                try:
                    with open(polar_file, "r") as pf:
                        lines = pf.readlines()
                    data = []
                    for ln in lines:
                        ln = ln.strip()
                        if not ln or ln.lower().startswith("alpha") or ln.lower().startswith("-----"):
                            continue
                        parts = re.split(r'\s+', ln)
                        try:
                            alpha = float(parts[0])
                            cl = float(parts[1])
                            cd = float(parts[2])
                            data.append((alpha, cl, cd))
                        except:
                            continue
                    if not data:
                        print("[DEBUG] no data parsed from polar file")
                        return None
                    cdmin_row = min(data, key=lambda x: x[2])
                    alpha_at_cdmin, cl_at_cdmin, cd_min = cdmin_row
                    print(f"[DEBUG] cdmin found: {cd_min:.5f}, cl: {cl_at_cdmin:.4f}, alpha: {alpha_at_cdmin:.3f}")
                    return cd_min, cl_at_cdmin, alpha_at_cdmin
                except Exception as e:
                    print(f"[DEBUG] exception parsing polar file: {e}")
                    return None
            else:
                return None



    def iterate_vms(self,
                    airfoil_name: Optional[str] = None,
                    initial_re_v: float = 10.0,
                    max_iter: int = 6,
                    tol: float = 1e-3,
                    use_xfoil: bool = False) -> Dict[str, List[float]]:
        """
        Iteratively estimate V_ms. If use_xfoil=True and XFOIL is available, the code will attempt to call XFOIL
        to re-evaluate cd0_2d at the updated Re. Otherwise it uses cd0_2d_fallback.
        Returns a history dict with keys: 'Re', 'cd0_2d', 'cd0_3d', 'CL_ms', 'V_ms'.
        """
        history = {'Re': [], 'cd0_2d': [], 'cd0_3d': [], 'CL_ms': [], 'V_ms': []}

        if self.chord is None:
            raise ValueError("Chord must be defined to calculate Re.")
        Re = initial_re_v * self.chord / NU_AIR

        last_v = None
        for it in range(max_iter):
            cd0_2d = None
            if use_xfoil and airfoil_name is not None:
                res = self._call_xfoil_get_cdmin(airfoil_name=airfoil_name, Re=Re)
                if res is not None:
                    cd0_2d = res[0]

            if cd0_2d is None:
                cd0_2d = self.cd0_2d_fallback

            cd0_3d = self.f3d * cd0_2d
            K = self.induced_drag_factor
            CL_ms = math.sqrt(max(0.0, 3.0 * cd0_3d / K))
            V_ms = math.sqrt((2.0 * self.weight) / (self.rho * self.wing_area * CL_ms))

            history['Re'].append(Re)
            history['cd0_2d'].append(cd0_2d)
            history['cd0_3d'].append(cd0_3d)
            history['CL_ms'].append(CL_ms)
            history['V_ms'].append(V_ms)

            if last_v is not None:
                rel_change = abs(V_ms - last_v) / last_v
                if rel_change < tol:
                    break
            last_v = V_ms
            Re = V_ms * self.chord / NU_AIR

        return history


# ---------------------------
# Example run (edit values as desired)
if __name__ == "__main__":
    gl = GliderIter(
        mass=0.10,          # kg
        chord=0.10,         # m
        wingspan=0.8,       # m
        oswald_eff=0.75,
        cd0_2d_fallback=0.02817,
        f3d=1.10,
        rho=1.225
    )

    # set use_xfoil=True only if xfoil is installed and you provide airfoil .dat name/path
    hist = gl.iterate_vms(use_xfoil=False, initial_re_v=8.0, max_iter=8, tol=1e-4)
    for i in range(len(hist['V_ms'])):
        print(f"Iter {i+1}: Re={hist['Re'][i]:.1f}, cd0_2d={hist['cd0_2d'][i]:.5f}, cd0_3d={hist['cd0_3d'][i]:.5f}, CL_ms={hist['CL_ms'][i]:.4f}, V_ms={hist['V_ms'][i]:.3f} m/s")
    hist = gl.iterate_vms(
    use_xfoil=True,
    airfoil_name="SD7037.dat",
    initial_re_v=8.0,
    max_iter=8,
    tol=1e-4
)