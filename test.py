"""
Robust loader for S3021 airfoil in AeroSandbox.
Tentativas:
  1) asb.Airfoil("s3021-il")
  2) asb.geometry.airfoil.airfoil_families.get_UIUC_coordinates(...)
  3) download from UIUC Selig site (fallback)
Se bem sucedido, imprime shape e escreve coords em CSV.
"""

import aerosandbox as asb
import os
import sys
import numpy as np

# Nome(s) a tentar (vários formatos comuns)
name_variants = [
    "s3021-il",
    "s3021",
    "S3021",
    "S3021-095-84",
    "s3021.dat",
]

def try_create_af_from_name(name):
    try:
        af = asb.Airfoil(name)
        if af is not None and af.coordinates is not None and len(af.coordinates) > 0:
            print(f"[OK] Loaded Airfoil via asb.Airfoil('{name}')")
            return af
        else:
            print(f"[MISS] asb.Airfoil('{name}') returned no coordinates.")
            return None
    except Exception as e:
        print(f"[ERR] asb.Airfoil('{name}') raised: {e}")
        return None

# 1) Try direct names
for n in name_variants:
    af = try_create_af_from_name(n)
    if af is not None:
        break

# 2) If not found, try UIUC helper (if present in this aerosandbox version)
if (('af' not in locals() or af is None)):
    try:
        from aerosandbox.geometry.airfoil import airfoil_families as af_fam
        got = None
        for n in name_variants:
            try:
                coords = af_fam.get_UIUC_coordinates(n)
                if coords is not None and len(coords) > 0:
                    print(f"[OK] Found UIUC coordinates for '{n}' via airfoil_families.get_UIUC_coordinates")
                    af = asb.Airfoil(coordinates=coords)
                    got = True
                    break
            except Exception as e:
                # função pode lançar se não encontrar
                # ignora e tenta a próxima variante
                # print(f"[dbg] get_UIUC_coordinates('{n}') -> {e}")
                pass
        if ('got' not in locals()):
            print("[MISS] airfoil_families.get_UIUC_coordinates didn't find usable coords for tried names.")
    except Exception as e:
        print(f"[WARN] Could not import airfoil_families helper: {e}")

# 3) If still not found, try to locate database directory and list candidates
if (('af' not in locals() or af is None)):
    try:
        db_dir = asb._asb_root / "geometry" / "airfoil" / "airfoil_database"
        print(f"[INFO] Looking into local DB dir: {db_dir}")
        if db_dir.exists():
            files = sorted(os.listdir(db_dir))
            # show matching names with 's3021' substring
            matches = [f for f in files if "s3021" in f.lower()]
            print(f"[INFO] DB files containing 's3021': {matches[:30] or 'NONE'}")
        else:
            print("[INFO] DB dir does not exist in this installation.")
    except Exception as e:
        print(f"[WARN] Error while listing local DB: {e}")

# 4) Fallback: try download from UIUC Selig server (requires internet)
if (('af' not in locals() or af is None)):
    try:
        print("[TRY] Attempting to download from UIUC Selig server (requires internet)...")
        import requests
        # common URL pattern for UIUC Selig coord files:
        possible_urls = [
            "https://m-selig.ae.illinois.edu/ads/coord/s3021.dat",
            "https://m-selig.ae.illinois.edu/ads/coord/s3021-il.dat",
            "https://m-selig.ae.illinois.edu/ads/coord/S3021.dat",
        ]
        coords = None
        for url in possible_urls:
            try:
                r = requests.get(url, timeout=10)
                if r.status_code == 200 and len(r.text) > 50:
                    print(f"[OK] downloaded from {url}")
                    # parse numeric lines (Selig .dat often starts with header + x y lines)
                    lines = [ln.strip() for ln in r.text.splitlines() if ln.strip() != ""]
                    # remove header line if non-numeric
                    parsed = []
                    for ln in lines[1:]:
                        parts = ln.split()
                        if len(parts) >= 2:
                            try:
                                x = float(parts[0])
                                y = float(parts[1])
                                parsed.append([x, y])
                            except:
                                # ignore non-numeric lines
                                pass
                    if len(parsed) > 10:
                        coords = np.array(parsed)
                        break
            except Exception as e:
                # try next URL
                # print(f"[dbg] download fail {url}: {e}")
                pass
        if coords is not None:
            af = asb.Airfoil(coordinates=coords)
            print("[OK] Created Airfoil from downloaded coordinates.")
        else:
            print("[MISS] Could not download/parse coordinates from UIUC URLs tried.")
    except Exception as e:
        print(f"[WARN] Download fallback skipped (requests failed or not installed): {e}")

# Final check
if ('af' in locals() and af is not None and getattr(af, "coordinates", None) is not None):
    print("Final: airfoil loaded. coordinates shape:", af.coordinates.shape)
    # save CSV for inspection
    np.savetxt("s3021_coords.csv", af.coordinates, delimiter=",", header="x,y", comments="")
    print("Coordinates saved to s3021_coords.csv")
else:
    print("FAILED: airfoil not loaded. Try updating AeroSandbox (pip install -U aerosandbox) or place a local s3021.dat and re-run.")
    sys.exit(1)
