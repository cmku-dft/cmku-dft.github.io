#!/usr/bin/env python3
"""
check_dos.py - sanity-check a VASP DOS before using it for carrier concentrations.

Usage (in the VASP run folder):
    python3 check_dos.py                 # reads DOSCAR (+ OUTCAR, EIGENVAL if present)
    python3 check_dos.py DOSCAR-orig
    python3 check_dos.py vasprun.xml
    python3 check_dos.py totdos.dat --nelect 36 [--vbm 0]   # 2-column text (E, DOS)

Checks:
  1. electrons in the valence band  (integrated DOS at the VBM)  == NELECT ?
  2. total states at EMAX           == ISPIN-factor x NBANDS ?   (whole-file scale)
  3. cell volume (DOSCAR line 2 is volume PER ATOM -> x NIONS)
  4. band edges from EIGENVAL vs from the DOS, energy step vs kT
"""
import sys, argparse
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz


def read_outcar(path="OUTCAR"):
    info = {}
    if not Path(path).exists():
        return info
    for line in open(path, errors="ignore"):
        if "NELECT" in line and "total number of electrons" in line:
            info["NELECT"] = float(line.split("=")[1].split()[0])
        elif "NBANDS=" in line:
            info["NBANDS"] = int(line.split("NBANDS=")[1].split()[0])
        elif "ISPIN" in line and "spin polarized" in line:
            info["ISPIN"] = int(line.split("=")[1].split()[0])
        elif "NIONS" in line and "=" in line and "NIONS" not in info:
            try: info["NIONS"] = int(line.split("NIONS =")[1].split()[0])
            except Exception: pass
    return info


def read_eigenval(path="EIGENVAL"):
    """VBM/CBM from occupations (VASP >= 5.4 EIGENVAL has occupations)."""
    if not Path(path).exists():
        return None
    L = Path(path).read_text().split("\n")
    ispin = int(L[0].split()[3])
    nel, nk, nb = map(int, map(float, L[5].split()[:3]))
    occ_e, emp_e, i = [], [], 7
    for _ in range(nk):
        i += 1  # k-point line
        for b in range(nb):
            p = L[i + b].split()
            if ispin == 1 and len(p) >= 3:
                (occ_e if float(p[2]) > 0.5 else emp_e).append(float(p[1]))
            elif ispin == 2 and len(p) >= 5:
                for e, o in ((p[1], p[3]), (p[2], p[4])):
                    (occ_e if float(o) > 0.5 else emp_e).append(float(e))
        i += nb + 1
    if not occ_e or not emp_e:
        return None
    return max(occ_e), min(emp_e)


def read_doscar(path):
    L = Path(path).read_text().split("\n")
    nions = int(L[0].split()[0])
    vol_atom = float(L[1].split()[0])
    nedos, efermi = int(float(L[5].split()[2])), float(L[5].split()[3])
    d = np.array([list(map(float, l.split())) for l in L[6:6 + nedos]])
    if d.shape[1] == 5:
        return d[:, 0], d[:, 1] + d[:, 2], d[:, 3] + d[:, 4], vol_atom * nions, efermi, 2, nions
    return d[:, 0], d[:, 1], d[:, 2], vol_atom * nions, efermi, 1, nions


def read_vasprun(path):
    r = ET.parse(path).getroot()
    ef = float(r.find(".//dos/i[@name='efermi']").text)
    sets = r.findall(".//dos/total/array/set/set")
    a = [np.array([list(map(float, x.text.split())) for x in s.findall("r")]) for s in sets]
    vol = None
    for st in r.findall("structure"):
        v = st.find(".//crystal/i[@name='volume']")
        if v is not None: vol = float(v.text)
    info = {}
    for k in ("NELECT", "NBANDS", "ISPIN"):
        x = r.find(f".//parameters//i[@name='{k}']")
        if x is not None: info[k] = float(x.text)
    nions = len(r.findall(".//atominfo/array[@name='atoms']/set/rc"))
    return a[0][:, 0], sum(x[:, 1] for x in a), sum(x[:, 2] for x in a), vol, ef, len(a), nions, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="?", default="DOSCAR")
    ap.add_argument("--nelect", type=float, help="valence electrons of the DOS cell (else from OUTCAR/vasprun)")
    ap.add_argument("--vbm", type=float, help="VBM energy on the file's energy scale (else EIGENVAL or DOS)")
    a = ap.parse_args()

    f = a.file
    info = read_outcar()
    if f.endswith(".xml"):
        E, D, I, vol, ef, ns, nions, xinfo = read_vasprun(f)
        info.update({k: v for k, v in xinfo.items() if k not in info})
        src = "vasprun.xml"
    elif "DOSCAR" in Path(f).name.upper():
        E, D, I, vol, ef, ns, nions = read_doscar(f)
        src = "DOSCAR"
    else:
        raw = np.loadtxt(f); E, D = raw[:, 0], raw[:, 1]
        I = np.concatenate([[0], np.cumsum(0.5 * (D[1:] + D[:-1]) * np.diff(E))])
        vol, ef, ns, nions, src = None, 0.0, 1, None, "text (E, DOS) - IDOS integrated numerically"
    nel = a.nelect or info.get("NELECT")

    # band edges
    edges = read_eigenval() if src != "text (E, DOS) - IDOS integrated numerically" else None
    tol = 1e-3 * D.max()
    if a.vbm is not None:
        vbm, how = a.vbm, "--vbm"
    elif edges:
        vbm, how = edges[0], "EIGENVAL occupations"
    else:
        occ = (E <= ef + 0.05) & (D > tol)
        vbm, how = E[occ].max(), "DOS threshold (approximate)"
    nval = float(np.interp(vbm, E, I))
    de = float(np.median(np.diff(E)))

    print(f"\nFile             : {f}  [{src}]")
    if vol: print(f"Cell volume      : {vol:.3f} A^3" + (f"  ({nions} atoms; DOSCAR line 2 = {vol/nions:.3f} A^3/atom)" if src == "DOSCAR" else ""))
    print(f"E_fermi (file)   : {ef:.4f} eV")
    print(f"VBM used         : {vbm:.4f} eV  [{how}]")
    if edges: print(f"CBM (EIGENVAL)   : {edges[1]:.4f} eV  -> gap {edges[1]-edges[0]:.3f} eV")
    print(f"Energy step      : {de*1000:.1f} meV  (kT = 25.9 meV at 300 K)" + ("  <-- too coarse, increase NEDOS / narrow EMIN-EMAX" if de > 0.0259 else ""))
    print(f"\nElectrons in VB  : {nval:.4f}")
    ok = True
    if nel:
        r = nval / nel
        print(f"NELECT           : {nel:g}   ratio = {r:.4f}", end="")
        if abs(r - 1) > 0.02:
            ok = False
            print(f"   <-- WRONG by x{1/r:.4g}")
        else:
            print("   OK")
    else:
        print("NELECT           : unknown (give --nelect or run where OUTCAR is)")
    nb = info.get("NBANDS")
    if nb and src != "text (E, DOS) - IDOS integrated numerically":
        expect = 2 * nb  # ISPIN=1: 2 e per band; ISPIN=2: 1 e per band per spin, spins summed
        tot = float(I[-1])
        print(f"States at EMAX   : {tot:.4f}   expected 2 x NBANDS = {expect}", end="")
        print("   OK" if abs(tot / expect - 1) < 0.02 else f"   ratio {tot/expect:.4f}  (if E_max is above all bands)")
    print("\nVERDICT: " + ("DOS normalisation OK" if ok else
          f"DOS scaled wrongly - use NELECT {nel:g} in defconchem (renormalises) or regenerate the DOS") + "\n")


if __name__ == "__main__":
    main()
