#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
defconchem.py  (v11.3)  -  self-consistent defect & carrier concentrations
=========================================================================

Solves charge neutrality  p - n + sum_i q_i c_i = 0  for the Fermi level,
using defect formation energies and the host density of states (DOS).
Optional freeze-in: totals are fixed at a growth temperature, then the
Fermi level, carriers and charge states re-equilibrate at a lower T.

Quick start
-----------
    python3 defconchem.py --template > input.txt   # annotated example input
    python3 defconchem.py input.txt                 # run
    python3 defconchem.py input.txt --check         # only print DOS / site diagnostics

Outputs (per DMU_RANGE "base" name):
    results_<base>.dat        final-temperature results (all species)
    results_<base>_growth.dat growth-temperature results (when freezing)
    plot_<base>.png / .pdf    concentrations and E_F vs chemical potential

Old (fz8/fz10) input files run unchanged. Changes vs. v10 (all physics fixes):
  1. Complexes in equilibrium no longer deplete isolated defects. A complex is a
     defect with E_f = sum(E_f of components) - E_bind (grand canonical).
     Mass conservation is applied only where totals really are fixed
     (frozen defects, or TOTAL_CONC dopants), via a proper solver.
  2. DOS read from VASP DOSCAR or vasprun.xml: both spin channels summed, cell
     volume taken from the DOS file itself, VBM/CBM located automatically,
     scissor shift applied if BANDGAP differs from the DOS gap.
  3. FREEZE_TEMPERATURE without FREEZE_DEFECTS is an error (it silently gave
     results at the wrong temperature before). Use FREEZE_DEFECTS ALL.
  4. Converged solvers with relative tolerances; warnings instead of silent
     clamping; log-space arithmetic (no overflow); no-DMU_RANGE runs a single point.
  5. NELECT keyword: DOS renormalised to NELECT valence electrons (as SC-FERMI does).
     v11.2: DOSCAR volume = line-2 value x NIONS (VASP writes volume per atom).
     v11.3: EFFECTIVE_MASS / CARRIER_MODEL - effective-mass (optionally Kane) bands as an
     alternative to the DOS, with an automatic Nc/Nv comparison when both are given.
  6. Diagnostics: electrons in the valence band (spin/normalisation check),
     effective Nc/Nv and equivalent DOS masses, site densities.
"""

import sys
import math
import argparse
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp

__version__ = "v11.3"
kB = 8.617333262e-5  # eV/K
if not hasattr(np, "trapezoid"):  # numpy < 2.0 (e.g. Pyodide)
    np.trapezoid = np.trapz

TEMPLATE = r"""# ===================== defconchem.py input =====================
# Lines starting with # are comments. Keywords are case-insensitive.
# All formation energies: E_f(q) = Ef0 + q*E_F + n*dmu
#   E_F measured from the VBM; Ef0 = formation energy at E_F = 0 (VBM)
#   at the reference chemical potential (dmu = 0, see DMU_RANGE / STOICH).

# ---------- host & DOS ----------
DOS_FILE      DOSCAR        # DOSCAR, vasprun.xml, or 2-column text (E-VBM, DOS)
BANDGAP       3.29          # eV. If different from the DOS gap, CB is scissor-shifted
# DOS_VOLUME_ANG3 47.57     # only needed for 2-column text DOS (volume of the DOS cell)
# SPIN 2                    # only for text DOS with columns E up down
# NELECT 36                 # valence electrons in the DOS cell: the DOS is renormalised so the VB
#                           #   holds NELECT electrons (like SC-FERMI). Strongly recommended.
# DOS_VBM 4.123             # optional override: absolute VBM energy in the DOS file
# DOS_CBM 7.410             # optional override: absolute CBM energy in the DOS file
#
# Alternative / cross-check: effective-mass bands (no DOS file needed)
# EFFECTIVE_MASS 0.28 0.59  # m_e* m_h* [alpha_e alpha_h]  (DOS masses; alpha = Kane
#                           #   non-parabolicity in 1/eV, default 0 = parabolic)
# CARRIER_MODEL  DOS        # DOS or EMA. If both DOS_FILE and EFFECTIVE_MASS are given,
#                           #   this picks which is used; the other is printed for comparison.

# ---------- defect-site cell ----------
# nsites in DEFECT lines are "sites per cell" of THIS cell (VOLUME_ANG3 or POSCAR):
VOLUME_ANG3   47.57         # e.g. primitive wurtzite ZnO cell (2 O, 2 Zn sites)
# POSCAR_FILE POSCAR

# ---------- temperatures ----------
TEMPERATURE         300     # final (measurement) temperature, K
FREEZE_TEMPERATURE  800     # growth temperature, K (omit for full equilibrium at TEMPERATURE)
FREEZE_DEFECTS      ALL     # ALL, or a list of defect/complex names to freeze

# ---------- defects ----------
# DEFECT <name> <sites per cell>
#   STATE q=<charge> Ef0=<eV> [g=<degeneracy>]
DEFECT ClO 2
  STATE q=1  Ef0=-0.38
  STATE q=0  Ef0=3.43
END
DEFECT VO 2
  STATE q=2  Ef0=-0.10
  STATE q=0  Ef0=4.26
END
DEFECT VZn 2
  STATE q=2  Ef0=2.584
  STATE q=1  Ef0=2.903
  STATE q=0  Ef0=3.800
  STATE q=-1 Ef0=5.154
  STATE q=-2 Ef0=6.964
END

# ---------- complexes (two equivalent ways) ----------
# (a) From binding energy: E_f(complex) = sum nu*E_f(component) - BINDING
#     Several COMPLEX blocks with the same name = charge states of one complex.
COMPLEX ClO-VZn
  COMP ClO_q1 1 VZn_q-2 1     # component charge states and multiplicities
  BINDING 2.09
  G_COMPLEX 1                 # orientational degeneracy per site
# SITES 2                     # sites/cell (default: sites of first component)
END
COMPLEX ClO-VZn
  COMP ClO_q1 1 VZn_q-1 1     # neutral complex
  BINDING 1.37
END
# (b) Directly from its own formation energies (like a DEFECT):
# DEFECT ClO-VZn 2
#   CONTAINS ClO 1 VZn 1      # atom bookkeeping for frozen totals
#   STATE q=0  Ef0=3.40
#   STATE q=-1 Ef0=4.49
# END

# ---------- optional fixed total concentration (e.g. known dopant level) ----------
# TOTAL_CONC ClO 1e19         # cm^-3, total incl. atoms bound in complexes

# ---------- chemical-potential scan ----------
# STOICH <defect> <element> <n> : adds n*dmu to E_f   (n = atoms REMOVED;
#   use effective slopes for linked potentials, e.g. Cl_O vs mu_O = +0.5)
# Complexes inherit the sum of their components' STOICH automatically.
STOICH ClO O  0.5
STOICH VO  O  1
STOICH VZn O -1
# DMU_RANGE <base> <element> <start> <end> <steps>
DMU_RANGE ClZnO O -3.22 0.0 41

# ---------- plotting (optional) ----------
# PLOT_MIN 1e6                # hide species that never exceed this (cm^-3)
# PLOT_SPECIES ClO_q1 VZn_q-2 ClO-VZn_q-1   # explicit list instead
"""


class InputError(Exception):
    pass


def die(msg):
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(1)


def warn(msg):
    print(f"  WARNING: {msg}")


# =====================================================================
# Input parsing
# =====================================================================
def parse_input(path: Path):
    cfg = {"DOS_FILE": None, "BANDGAP": None, "SPIN": 1, "POSCAR_FILE": "POSCAR",
           "VOLUME_ANG3": None, "DOS_VOLUME_ANG3": None, "DOS_VBM": None, "DOS_CBM": None,
           "NSITES_CM3": None, "TEMPERATURE": 300.0, "FREEZE_TEMPERATURE": None,
           "FREEZE_DEFECTS": [], "PLOT_MIN": 1e6, "PLOT_SPECIES": [], "DOS_TOL": None,
           "NELECT": None, "EFFECTIVE_MASS": None, "CARRIER_MODEL": None}
    defects, legacy_cx, dmu_ranges, stoich, total_conc = {}, [], [], {}, {}
    cur, cur_cx = None, None
    if not path.exists():
        die(f"input file '{path}' not found")
    for ln, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        key = parts[0].upper()
        try:
            if cur is not None:
                if key == "END":
                    cur = None
                elif key == "STATE":
                    q, e0, g = None, None, 1.0
                    for tok in parts[1:]:
                        if "=" in tok:
                            k, v = tok.split("=", 1)
                            k = k.lower()
                            if k == "q": q = int(v)
                            elif k == "ef0": e0 = float(v)
                            elif k == "g": g = float(v)
                    if q is None:
                        q, e0 = int(parts[1]), float(parts[2])
                        if len(parts) > 3: g = float(parts[3])
                    if e0 is None:
                        raise InputError("STATE needs Ef0")
                    cur["states"].append({"q": q, "Ef0": e0, "g": g})
                elif key == "CONTAINS":
                    cur["contains"] = {parts[i]: float(parts[i + 1]) for i in range(1, len(parts), 2)}
                else:
                    raise InputError(f"unexpected '{parts[0]}' inside DEFECT block (missing END?)")
                continue
            if cur_cx is not None:
                if key == "END":
                    legacy_cx.append(cur_cx); cur_cx = None
                elif key == "COMP":
                    cur_cx["components"] = [(parts[i], float(parts[i + 1])) for i in range(1, len(parts), 2)]
                elif key == "BINDING": cur_cx["E_bind"] = float(parts[1])
                elif key == "G_COMPLEX": cur_cx["g"] = float(parts[1])
                elif key == "SITES": cur_cx["sites"] = float(parts[1])
                else:
                    raise InputError(f"unexpected '{parts[0]}' inside COMPLEX block (missing END?)")
                continue
            if key == "DEFECT":
                name = parts[1]
                if name in defects:
                    raise InputError(f"defect '{name}' defined twice")
                cur = {"name": name, "sites": float(parts[2]), "states": [], "contains": {}}
                defects[name] = cur
            elif key == "COMPLEX":
                cur_cx = {"name": parts[1], "components": [], "E_bind": 0.0, "g": 1.0, "sites": None}
            elif key == "STOICH":
                stoich.setdefault(parts[1], {})[parts[2]] = float(parts[3])
            elif key == "DMU_RANGE":
                dmu_ranges.append((parts[1], parts[2], float(parts[3]), float(parts[4]), int(parts[5])))
            elif key == "FREEZE_DEFECTS":
                cfg["FREEZE_DEFECTS"].extend(parts[1:])
            elif key == "PLOT_SPECIES":
                cfg["PLOT_SPECIES"].extend(parts[1:])
            elif key == "TOTAL_CONC":
                total_conc[parts[1]] = float(parts[2])
            elif key == "EFFECTIVE_MASS":
                v = [float(x) for x in parts[1:]]
                if len(v) < 2:
                    raise InputError("EFFECTIVE_MASS needs me mh [alpha_e alpha_h]")
                cfg["EFFECTIVE_MASS"] = v + [0.0] * (4 - len(v))
            elif key == "CARRIER_MODEL":
                cfg["CARRIER_MODEL"] = parts[1].upper()
            elif key in cfg:
                cfg[key] = parts[1] if key in ("DOS_FILE", "POSCAR_FILE") else float(parts[1])
            else:
                warn(f"line {ln}: unknown keyword '{parts[0]}' ignored")
        except (IndexError, ValueError) as e:
            die(f"line {ln}: cannot parse '{raw.strip()}' ({e})")
        except InputError as e:
            die(f"line {ln}: {e}")
    if cur is not None or cur_cx is not None:
        die("input ended inside a DEFECT/COMPLEX block (missing END)")
    if cfg["BANDGAP"] is None:
        die("BANDGAP is required")
    if cfg["DOS_FILE"] is None and cfg["EFFECTIVE_MASS"] is None:
        die("give DOS_FILE and/or EFFECTIVE_MASS")
    cm = cfg["CARRIER_MODEL"]
    if cm is None:
        cfg["CARRIER_MODEL"] = "DOS" if cfg["DOS_FILE"] else "EMA"
    elif cm not in ("DOS", "EMA"):
        die("CARRIER_MODEL must be DOS or EMA")
    elif cm == "DOS" and not cfg["DOS_FILE"]:
        die("CARRIER_MODEL DOS needs DOS_FILE")
    elif cm == "EMA" and not cfg["EFFECTIVE_MASS"]:
        die("CARRIER_MODEL EMA needs EFFECTIVE_MASS")

    # --- convert legacy COMPLEX blocks into defects with CONTAINS ---
    for cx in legacy_cx:
        if not cx["components"]:
            die(f"COMPLEX {cx['name']}: no COMP line")
        q_tot, e0_tot, contains = 0, -cx["E_bind"], {}
        for sname, nu in cx["components"]:
            if "_q" not in sname:
                die(f"COMPLEX {cx['name']}: component '{sname}' must look like <defect>_q<charge>, e.g. VZn_q-2")
            par, qs = sname.rsplit("_q", 1)
            if par not in defects:
                die(f"COMPLEX {cx['name']}: unknown defect '{par}'")
            st = [s for s in defects[par]["states"] if s["q"] == int(qs)]
            if not st:
                die(f"COMPLEX {cx['name']}: defect {par} has no charge state {qs}")
            q_tot += nu * int(qs)
            e0_tot += nu * st[0]["Ef0"]
            contains[par] = contains.get(par, 0) + nu
        first = cx["components"][0][0].rsplit("_q", 1)[0]
        d = defects.get(cx["name"])
        if d is None:
            d = {"name": cx["name"], "sites": cx["sites"] or defects[first]["sites"],
                 "states": [], "contains": contains, "legacy": True}
            defects[cx["name"]] = d
        elif d.get("contains") != contains:
            die(f"COMPLEX {cx['name']}: charge states must contain the same defects")
        if any(s["q"] == int(round(q_tot)) for s in d["states"]):
            die(f"COMPLEX {cx['name']}: charge state {q_tot} defined twice")
        q_tot = int(round(q_tot))
        d["states"].append({"q": q_tot, "Ef0": e0_tot, "g": cx["g"]})

    # --- checks & stoichiometry inheritance ---
    for name, d in defects.items():
        if not d["states"]:
            die(f"defect {name} has no STATE lines")
        for p in d["contains"]:
            if p not in defects or defects[p]["contains"]:
                die(f"{name}: CONTAINS '{p}' must be a simple (non-complex) defect")
        if name not in stoich and d["contains"]:
            inh = {}
            for p, nu in d["contains"].items():
                for el, n in stoich.get(p, {}).items():
                    inh[el] = inh.get(el, 0) + nu * n
            stoich[name] = inh
    for p in total_conc:
        if p not in defects or defects[p]["contains"]:
            die(f"TOTAL_CONC: '{p}' must be a simple defect")
    for s in stoich:
        if s not in defects:
            die(f"STOICH refers to unknown defect '{s}'")

    fr = cfg["FREEZE_DEFECTS"]
    if cfg["FREEZE_TEMPERATURE"]:
        if not fr:
            die("FREEZE_TEMPERATURE is set but FREEZE_DEFECTS is empty. "
                "Use 'FREEZE_DEFECTS ALL' or list names (or remove FREEZE_TEMPERATURE).")
        if [x.upper() for x in fr] == ["ALL"]:
            cfg["FREEZE_DEFECTS"] = list(defects)
        for x in cfg["FREEZE_DEFECTS"]:
            if x not in defects:
                die(f"FREEZE_DEFECTS: unknown name '{x}' (names: {', '.join(defects)})")
    elif fr:
        warn("FREEZE_DEFECTS given without FREEZE_TEMPERATURE - ignored (full equilibrium)")
        cfg["FREEZE_DEFECTS"] = []
    if cfg["NSITES_CM3"] is not None:
        warn("NSITES_CM3 is no longer used: complex site densities come from SITES / component sites")
    return cfg, defects, dmu_ranges, stoich, total_conc


# =====================================================================
# DOS reading
# =====================================================================
def read_doscar(path):
    lines = Path(path).read_text().splitlines()
    nions = int(lines[0].split()[0])
    vol = float(lines[1].split()[0]) * nions  # DOSCAR line 2 = volume PER ATOM (A^3)
    emax, emin, nedos, efermi = lines[5].split()[:4]
    nedos, efermi = int(nedos), float(efermi)
    data = np.array([list(map(float, l.split())) for l in lines[6:6 + nedos]])
    if data.shape[1] == 5:  # E up dn intup intdn
        return data[:, 0], data[:, 1] + data[:, 2], data[:, 3] + data[:, 4], vol, efermi, 2
    return data[:, 0], data[:, 1], data[:, 2], vol, efermi, 1


def read_vasprun(path):
    root = ET.parse(path).getroot()
    efermi = float(root.find(".//dos/i[@name='efermi']").text)
    sets = root.findall(".//dos/total/array/set/set")
    arrs = [np.array([list(map(float, r.text.split())) for r in s.findall("r")]) for s in sets]
    E = arrs[0][:, 0]
    dos = sum(a[:, 1] for a in arrs)
    idos = sum(a[:, 2] for a in arrs)
    vol = None
    for st in root.findall("structure"):
        v = st.find(".//crystal/i[@name='volume']")
        if v is not None:
            vol = float(v.text)
    nel = root.find(".//parameters//i[@name='NELECT']")
    if vol is None:
        die("could not find cell volume in vasprun.xml")
    return E, dos, idos, vol, efermi, len(arrs), (float(nel.text) if nel is not None else None)


def poscar_volume(path):
    lines = Path(path).read_text().splitlines()
    s = float(lines[1].split()[0])
    M = np.array([list(map(float, l.split()[:3])) for l in lines[2:5]])
    v = abs(np.linalg.det(M))
    return -s if s < 0 else v * s ** 3


class DOS:
    """Host DOS referenced to VBM = 0, CB scissor-shifted to BANDGAP. Units: states/eV/cm^3."""

    def __init__(self, cfg, site_vol):
        f = cfg["DOS_FILE"]
        if not Path(f).exists():
            die(f"DOS file '{f}' not found")
        name = Path(f).name.lower()
        self.nelect = None
        self.idos = None
        if name.endswith(".xml"):
            E, D, I, vol, ef, ns, self.nelect = read_vasprun(f)
            self.source = f"vasprun.xml (ISPIN={ns}, spins summed)"
            self.idos = I
        elif "doscar" in name:
            E, D, I, vol, ef, ns = read_doscar(f)
            self.source = f"DOSCAR (ISPIN={ns}, spins summed)"
            self.idos = I
        else:  # legacy text: E (already VBM=0), DOS [, DOS_down]
            raw = np.loadtxt(f)
            E = raw[:, 0]
            D = raw[:, 1] + (raw[:, 2] if (int(cfg["SPIN"]) == 2 and raw.shape[1] > 2) else 0)
            ef = 0.0
            vol = cfg["DOS_VOLUME_ANG3"]
            self.source = "text file (assumed: energies relative to VBM, DOS per cell, " + \
                          ("columns up+down summed)" if int(cfg["SPIN"]) == 2 else "total DOS incl. both spins)")
            if vol is None:
                vol = site_vol
                warn("text DOS without DOS_VOLUME_ANG3: assuming the DOS was computed in the "
                     "VOLUME_ANG3/POSCAR cell. Check this - a mismatch scales n and p directly!")
            if cfg["DOS_VBM"] is None:
                cfg["DOS_VBM"] = 0.0
        self.vol = vol
        tol = cfg["DOS_TOL"] or 1e-3 * D.max()
        # --- locate band edges ---
        if cfg["DOS_VBM"] is not None:
            vbm = cfg["DOS_VBM"]
        else:
            occ = (E <= ef + 0.05) & (D > tol)
            if not occ.any():
                die("could not find the VBM in the DOS; set DOS_VBM")
            vbm = E[occ].max()
        if cfg["DOS_CBM"] is not None:
            cbm = cfg["DOS_CBM"]
        else:
            emp = (E > vbm + 0.05) & (D > tol)
            # skip smearing tail: first point above tol after a gap region
            gapstart = np.where((E > vbm) & (D <= tol))[0]
            if len(gapstart) == 0:
                die("no gap found in the DOS; set DOS_VBM / DOS_CBM")
            emp = (E > E[gapstart[0]]) & (D > tol)
            cbm = E[emp].min()
        self.vbm_abs, self.cbm_abs = vbm, cbm
        self.auto_edges = cfg["DOS_VBM"] is None or cfg["DOS_CBM"] is None
        self.eg_dos = cbm - vbm
        self.eg = float(cfg["BANDGAP"])
        self.shift = self.eg - self.eg_dos
        vcm = vol * 1e-24
        iv, ic = E <= vbm, E >= cbm
        self.Ev, self.Dv = E[iv] - vbm, D[iv] / vcm
        self.Ec, self.Dc = E[ic] - vbm + self.shift, D[ic] / vcm
        self.nval = float(np.interp(vbm, E, self.idos)) if self.idos is not None else \
            float(np.trapezoid(D[iv], E[iv]))
        self.de = float(np.median(np.diff(E)))
        self.scale = 1.0
        ne = cfg["NELECT"]
        if ne is not None:
            if self.nelect is None:
                self.nelect = ne
            if abs(self.nval - ne) > 0.02 * ne:
                self.scale = ne / self.nval
                self.Dv *= self.scale
                self.Dc *= self.scale
                warn(f"DOS renormalised x{self.scale:.4g} so that the VB holds NELECT = {ne:g} electrons "
                     f"(file gave {self.nval:.4g}). Check how the DOS file was made!")
                self.nval = ne
        elif self.nval < 1.0:
            warn(f"VB holds only {self.nval:.3g} electrons per cell - the DOS normalisation is almost "
                 "certainly wrong. Add NELECT <valence electrons of the DOS cell> to renormalise.")
        elif "text" in self.source:
            warn("text DOS without NELECT: normalisation cannot be checked. "
                 "Add NELECT <valence electrons of the DOS cell> (e.g. 36 for 4-atom ZnO with Zn 3d)")

    def carriers(self, EF, T):
        kT = kB * T
        n = np.trapezoid(self.Dc * expit(-(self.Ec - EF) / kT), self.Ec)
        p = np.trapezoid(self.Dv * expit((self.Ev - EF) / kT), self.Ev)
        return n, p

    def eff_N(self, T):
        """Effective Nc, Nv (Boltzmann limit) and equivalent DOS masses."""
        kT = kB * T
        Nc = np.trapezoid(self.Dc * np.exp(-(self.Ec - self.eg) / kT), self.Ec)
        Nv = np.trapezoid(self.Dv * np.exp(self.Ev / kT), self.Ev)
        N1 = 2.5094e19 * (T / 300) ** 1.5  # Nc for m* = 1
        return Nc, Nv, (Nc / N1) ** (2 / 3), (Nv / N1) ** (2 / 3)

    def report(self, T_list):
        print(f"  DOS source        : {self.source}")
        print(f"  DOS cell volume   : {self.vol:.3f} A^3")
        if self.scale != 1.0:
            print(f"  DOS rescaled by   : x{self.scale:.4g} (NELECT)")
        print(f"  VBM, CBM in file  : {self.vbm_abs:.3f}, {self.cbm_abs:.3f} eV  (DOS gap {self.eg_dos:.3f} eV)")
        print(f"  BANDGAP used      : {self.eg:.3f} eV  -> CB scissor shift {self.shift:+.3f} eV")
        line = f"  Electrons in VB   : {self.nval:.2f} per DOS cell"
        if self.nelect is not None:
            line += f"  (NELECT = {self.nelect:.2f})"
            if abs(self.nval - self.nelect) > 0.05 * self.nelect:
                line += "  <-- MISMATCH: check DOS normalisation / VBM"
        print(line + "   [should equal the valence electron count of the DOS cell]")
        print(f"  DOS energy step   : {self.de * 1000:.1f} meV")
        for T in T_list:
            Nc, Nv, mc, mv = self.eff_N(T)
            print(f"  T = {T:6.0f} K     : Nc = {Nc:.2e}, Nv = {Nv:.2e} cm^-3  "
                  f"(equivalent DOS masses m_c* = {mc:.2f}, m_v* = {mv:.2f})")
            if self.de > kB * T:
                warn(f"DOS energy step ({self.de*1000:.0f} meV) > kT at {T:.0f} K - carrier integrals may be inaccurate (increase NEDOS)")
        if self.auto_edges and "text" not in self.source:
            print("  (band edges found from a DOS threshold; smearing tails can shift them by ~sigma. For accuracy set\n"
                  "   DOS_VBM / DOS_CBM to the VBM/CBM eigenvalues, e.g. VASPKIT 911 or EIGENVAL, same energy scale as the DOS)")
        if abs(self.shift) > 0.02:
            print("  (CB states shifted rigidly; make sure Ef0 values use the same band edges)")


class EMA:
    """Effective-mass bands (Kane non-parabolic if alpha > 0): E(1 + alpha E) = hbar^2 k^2 / 2m.
    g(E) = (1/2pi^2) (2m/hbar^2)^1.5 sqrt(E(1+aE)) (1+2aE), both spins. Units: states/eV/cm^3.
    m is the DOS mass per band edge: for several valleys/bands use m_dos = (sum_i m_i^1.5)^(2/3)."""
    HB2M = 3.80998  # hbar^2 / 2 m_e  [eV A^2]

    def __init__(self, cfg, Tmax):
        me, mh, ae, ah = cfg["EFFECTIVE_MASS"]
        self.me, self.mh, self.ae, self.ah = me, mh, ae, ah
        self.eg = float(cfg["BANDGAP"])
        emax = max(1.5, 50 * kB * Tmax)
        x = np.linspace(0.0, emax, 20001)
        self.Ec, self.Dc = self.eg + x, self._g(x, me, ae)
        self.Ev, self.Dv = -x[::-1], self._g(x, mh, ah)[::-1]
        self.source = (f"effective-mass model: m_e* = {me:g}, m_h* = {mh:g}"
                       + (f", alpha_e = {ae:g}, alpha_h = {ah:g} 1/eV (Kane)" if (ae or ah) else " (parabolic)"))

    def _g(self, x, m, a):
        return (1 / (2 * np.pi ** 2)) * (m / self.HB2M) ** 1.5 * np.sqrt(x * (1 + a * x)) * (1 + 2 * a * x) * 1e24

    carriers = DOS.carriers
    eff_N = DOS.eff_N

    def report(self, T_list):
        print(f"  Carrier model     : {self.source}")
        print(f"  BANDGAP used      : {self.eg:.3f} eV")
        for T in T_list:
            Nc, Nv, mc, mv = self.eff_N(T)
            print(f"  T = {T:6.0f} K     : Nc = {Nc:.2e}, Nv = {Nv:.2e} cm^-3  "
                  f"(equivalent DOS masses m_c* = {mc:.2f}, m_v* = {mv:.2f})")


# =====================================================================
# Defect model
# =====================================================================
class Model:
    def __init__(self, defects, stoich, site_vol, total_conc):
        self.vcm = site_vol * 1e-24
        self.D = defects
        self.simple = [n for n, d in defects.items() if not d["contains"]]
        self.cx = [n for n, d in defects.items() if d["contains"]]
        self.stoich = stoich
        self.total_conc = total_conc
        self.N = {n: d["sites"] / self.vcm for n, d in defects.items()}

    def species(self):
        return [(n, s["q"]) for n, d in self.D.items() for s in d["states"]]

    def logc_eq(self, EF, T, dmu, elem):
        """log of grand-canonical concentrations, per (defect, q)."""
        kT = kB * T
        out = {}
        for n, d in self.D.items():
            nmu = self.stoich.get(n, {}).get(elem, 0.0)
            lN = math.log(self.N[n])
            out[n] = {s["q"]: lN + math.log(s["g"]) - (s["Ef0"] + s["q"] * EF + nmu * dmu) / kT
                      for s in d["states"]}
        return out

    def concentrations(self, EF, T, dmu, elem, fixed_tot=None, fixed_cx=None):
        """
        fixed_tot : {simple defect: total atoms incl. those bound in complexes}
        fixed_cx  : {complex or simple defect: frozen total of that species (its own)}
        Returns {defect: {q: conc}}.
        """
        fixed_tot = dict(fixed_tot or {})
        fixed_cx = dict(fixed_cx or {})
        lc = self.logc_eq(EF, T, dmu, elem)
        logZ = {n: logsumexp(list(v.values())) for n, v in lc.items()}
        la = {n: 0.0 for n in self.D}  # log activity multipliers

        # simple defects with fixed own total (frozen, not in fixed_tot): scale directly
        for n, tot in fixed_cx.items():
            la[n] = math.log(max(tot, 1e-300)) - logZ[n]
        # defects with fixed atom totals -> solve conservation with free complexes
        cons = [p for p in fixed_tot if p not in fixed_cx]
        free_cx = [c for c in self.cx if c not in fixed_cx]
        if cons:
            rhs = {}
            for p in cons:
                r = fixed_tot[p] - sum(self.D[c]["contains"].get(p, 0) * fixed_cx[c]
                                       for c in self.cx if c in fixed_cx)
                if r <= 0:
                    warn(f"atom balance for {p} non-positive ({r:.3e}); set to tiny")
                    r = 1e-300
                rhs[p] = math.log(r)

            def lcx(c):
                return logZ[c] + sum(nu * la.get(p, 0.0) for p, nu in self.D[c]["contains"].items())

            for it in range(200):
                old = dict(la)
                for p in cons:
                    links = [(c, self.D[c]["contains"][p]) for c in free_cx if p in self.D[c]["contains"]]

                    def f(x):
                        la[p] = x
                        terms = [x + logZ[p]] + [math.log(nu) + lcx(c) for c, nu in links]
                        return logsumexp(terms) - rhs[p]
                    la[p] = brentq(f, -3000, 3000, xtol=1e-12)
                if max(abs(la[p] - old[p]) for p in cons) < 1e-10:
                    break
            else:
                warn("mass-conservation solver did not fully converge")
            for c in free_cx:
                la[c] = sum(nu * la.get(p, 0.0) for p, nu in self.D[c]["contains"].items())
        else:
            for c in free_cx:
                la[c] = sum(nu * la.get(p, 0.0) for p, nu in self.D[c]["contains"].items())

        conc = {n: {q: math.exp(min(v + la[n], 700)) for q, v in lc[n].items()} for n in self.D}
        return conc

    def charge(self, conc):
        return sum(q * c for d in conc.values() for q, c in d.items())

    def atom_totals(self, conc):
        tot = {p: sum(conc[p].values()) for p in self.simple}
        for c in self.cx:
            ctot = sum(conc[c].values())
            for p, nu in self.D[c]["contains"].items():
                tot[p] += nu * ctot
        return tot


def solve_EF(model, dos, T, dmu, elem, fixed_tot=None, fixed_cx=None):
    def Q(EF):
        n, p = dos.carriers(EF, T)
        return p - n + model.charge(model.concentrations(EF, T, dmu, elem, fixed_tot, fixed_cx))

    lo, hi = -1.0, dos.eg + 1.0
    qlo, qhi = Q(lo), Q(hi)
    k = 0
    while qlo * qhi > 0 and k < 4:
        lo, hi = lo - 0.5, hi + 0.5
        qlo, qhi = Q(lo), Q(hi)
        k += 1
    flag = ""
    if qlo * qhi > 0:
        EF = lo if abs(qlo) < abs(qhi) else hi
        flag = "NO_ROOT"
        warn(f"dmu={dmu:+.3f}, T={T:.0f} K: no charge-neutral E_F in [{lo:.1f}, {hi:.1f}] eV - "
             f"reporting edge value {EF:.2f} eV (check inputs)")
    else:
        EF = brentq(Q, lo, hi, xtol=1e-8, rtol=1e-12, maxiter=500)
    conc = model.concentrations(EF, T, dmu, elem, fixed_tot, fixed_cx)
    n, p = dos.carriers(EF, T)
    for nm, d in conc.items():
        if sum(d.values()) > 0.05 * model.N[nm]:
            warn(f"dmu={dmu:+.3f}: [{nm}] > 5% of its sites - dilute-limit formulas are unreliable")
    return {"dmu": dmu, "EF": EF, "n": n, "p": p, "conc": conc,
            "tot": model.atom_totals(conc), "flag": flag}


def run_point(model, dos, cfg, dmu, elem):
    T, Tg = cfg["TEMPERATURE"], cfg["FREEZE_TEMPERATURE"]
    tc = model.total_conc or None
    if not Tg:
        return None, solve_EF(model, dos, T, dmu, elem, fixed_tot=tc)
    g = solve_EF(model, dos, Tg, dmu, elem, fixed_tot=tc)
    frozen = set(cfg["FREEZE_DEFECTS"])
    fixed_tot, fixed_cx = {}, {}
    for p in model.simple:
        if p in frozen:
            fixed_tot[p] = g["tot"][p]          # atoms conserved (incl. in complexes)
        elif p in model.total_conc:
            fixed_tot[p] = model.total_conc[p]
    for c in model.cx:
        if c in frozen:
            fixed_cx[c] = sum(g["conc"][c].values())   # complex population fixed
    # frozen simple defects that appear in no free complex: fix their own total directly
    free_cx = [c for c in model.cx if c not in fixed_cx]
    for p in list(fixed_tot):
        if p in frozen and not any(p in model.D[c]["contains"] for c in free_cx):
            fixed_cx[p] = sum(g["conc"][p].values())
            del fixed_tot[p]
    f = solve_EF(model, dos, T, dmu, elem, fixed_tot=fixed_tot, fixed_cx=fixed_cx)
    return g, f


# =====================================================================
# Output
# =====================================================================
def write_table(fname, res, model, elem):
    sp = model.species()
    hdr = [f"dmu_{elem}", "EF", "n", "p"] + [f"{n}_q{q}" for n, q in sp] + \
          [f"Tot_{p}" for p in model.simple] + [f"Tot_{c}" for c in model.cx]
    with open(fname, "w") as fh:
        fh.write("# E_F in eV above VBM; concentrations in cm^-3. Tot_<simple> includes atoms bound in complexes.\n")
        fh.write("# " + "\t".join(hdr) + "\n")
        for r in res:
            row = [r["dmu"], r["EF"], r["n"], r["p"]] + [r["conc"][n][q] for n, q in sp] + \
                  [r["tot"][p] for p in model.simple] + [sum(r["conc"][c].values()) for c in model.cx]
            fh.write("\t".join(f"{x:.5e}" for x in row) + ("\t# " + r["flag"] if r["flag"] else "") + "\n")


PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def make_plot(base, elem, panels, model, dos, cfg):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        warn("matplotlib not installed - skipping plots")
        return
    sp = model.species()
    pmin = cfg["PLOT_MIN"]
    if cfg["PLOT_SPECIES"]:
        chosen = [s for s in sp if f"{s[0]}_q{s[1]}" in cfg["PLOT_SPECIES"]]
    else:
        chosen = [s for s in sp if any(r["conc"][s[0]][s[1]] > pmin for _, res in panels for r in res)]
    color = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(chosen)}
    styles = ["-", "--", ":", "-."]
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8})
    fig, ax = plt.subplots(2, len(panels), figsize=(3.6 * len(panels) + 1.4, 5.6), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]}, squeeze=False)
    allc = [v for _, res in panels for r in res for v in
            [r["n"], r["p"]] + [r["conc"][s[0]][s[1]] for s in chosen]]
    top = 10 ** math.ceil(math.log10(max(allc + [pmin * 10])))
    for j, (title, res) in enumerate(panels):
        x = [r["dmu"] for r in res]
        a = ax[0, j]
        for i, s in enumerate(chosen):
            y = [r["conc"][s[0]][s[1]] for r in res]
            qq = s[1]
            sup = "0" if qq == 0 else (f"{abs(qq) if abs(qq) > 1 else ''}" + ("+" if qq > 0 else "-"))
            lab = f"{s[0]}$^{{{sup}}}$"
            a.semilogy(x, y, color=color[s], ls=styles[(i // len(PALETTE)) % 4], lw=1.6, label=lab)
        a.semilogy(x, [r["n"] for r in res], color="#222222", lw=1.6, ls="--", label="n")
        a.semilogy(x, [r["p"] for r in res], color="#888888", lw=1.4, ls=":", label="p")
        a.set_ylim(pmin, top)
        a.set_title(title, fontsize=10)
        a.grid(True, which="major", color="#e5e5e5", lw=0.6)
        if j == 0:
            a.set_ylabel(r"Concentration (cm$^{-3}$)")
        b = ax[1, j]
        b.plot(x, [r["EF"] for r in res], color="#222222", lw=1.6)
        b.axhline(dos.eg, color="#888888", lw=0.8, ls="--")
        b.axhline(0, color="#888888", lw=0.8, ls="--")
        b.set_xlabel(rf"$\Delta\mu_{{\mathrm{{{elem}}}}}$ (eV)")
        b.grid(True, color="#e5e5e5", lw=0.6)
        if j == 0:
            b.set_ylabel(r"$E_F$ (eV)")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="center right", frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1 - 1.3 / fig.get_figwidth(), 1))
    for ext in ("png", "pdf"):
        fig.savefig(f"plot_{base}.{ext}", dpi=300)
    plt.close(fig)


# =====================================================================
# Main / API
# =====================================================================
def _row(r, model):
    return {"dmu": r["dmu"], "EF": r["EF"], "n": r["n"], "p": r["p"], "flag": r["flag"],
            "conc": {f"{n}_q{q}": r["conc"][n][q] for n, q in model.species()},
            "tot": {**{p: r["tot"][p] for p in model.simple},
                    **{c: sum(r["conc"][c].values()) for c in model.cx}}}


def run(input_path, check=False, plot=True, write_files=True, progress=None):
    """Run a calculation. Returns a dict with diagnostics and all results (also used by the web app)."""
    cfg, defects, dmu_ranges, stoich, total_conc = parse_input(Path(input_path))

    if cfg["VOLUME_ANG3"]:
        site_vol = cfg["VOLUME_ANG3"]
    elif Path(cfg["POSCAR_FILE"]).exists():
        site_vol = poscar_volume(cfg["POSCAR_FILE"])
    else:
        die("give VOLUME_ANG3 or a POSCAR_FILE for the defect-site cell")

    print(f"\n=== defconchem {__version__} ===")
    Ts = [cfg["TEMPERATURE"]] + ([cfg["FREEZE_TEMPERATURE"]] if cfg["FREEZE_TEMPERATURE"] else [])
    dos_file = DOS(cfg, site_vol) if cfg["DOS_FILE"] else None
    ema = EMA(cfg, max(Ts)) if cfg["EFFECTIVE_MASS"] else None
    dos = dos_file if cfg["CARRIER_MODEL"] == "DOS" else ema
    model = Model(defects, stoich, site_vol, total_conc)
    print(f"  Carriers from     : {'DOS file' if dos is dos_file else 'effective-mass model'}"
          f"  (CARRIER_MODEL {cfg['CARRIER_MODEL']})")
    dos.report(Ts)
    if dos_file is not None and ema is not None:   # comparison of the two band models
        other = ema if dos is dos_file else dos_file
        if other is dos_file:
            print("  --- DOS file (for comparison only) ---")
        else:
            print(f"  --- {ema.source} (for comparison only) ---")
        for T in Ts:
            a1, b1, _, _ = dos.eff_N(T)
            a2, b2, _, _ = other.eff_N(T)
            print(f"  T = {T:6.0f} K     : Nc ratio (used/other) = {a1/a2:.2f}, Nv ratio = {b1/b2:.2f}"
                  + ("   <-- differ by >30%: DOS not converged or masses inconsistent" if abs(a1/a2 - 1) > 0.3 else ""))
    print(f"  Defect-site cell  : {site_vol:.3f} A^3")
    for n, d in defects.items():
        kind = f"complex of {' + '.join(f'{v:g} {p}' for p, v in d['contains'].items())}" if d["contains"] else "defect"
        qs = ", ".join(f"{s['q']:+d}" for s in d["states"])
        print(f"    {n:12s} {kind:28s} sites {model.N[n]:.3e} cm^-3   q = {qs}")
    if cfg["FREEZE_TEMPERATURE"]:
        print(f"  Freeze-in         : grown at {cfg['FREEZE_TEMPERATURE']:.0f} K, evaluated at "
              f"{cfg['TEMPERATURE']:.0f} K; frozen: {', '.join(cfg['FREEZE_DEFECTS'])}")
    else:
        print(f"  Full equilibrium at {cfg['TEMPERATURE']:.0f} K")
    if total_conc:
        print("  Fixed totals      : " + ", ".join(f"{k} = {v:.2e}" for k, v in total_conc.items()))

    out = {"version": __version__, "species": [f"{n}_q{q}" for n, q in model.species()],
           "simple": model.simple, "complexes": model.cx, "bandgap": dos.eg,
           "T": cfg["TEMPERATURE"], "T_growth": cfg["FREEZE_TEMPERATURE"],
           "plot_min": cfg["PLOT_MIN"], "scans": [], "files": []}
    if check:
        return out

    if not dmu_ranges:
        print("  No DMU_RANGE: single point at dmu = 0")
        dmu_ranges = [("run", "X", 0.0, 0.0, 1)]
    total = sum(r[4] for r in dmu_ranges)
    done = 0
    for base, elem, start, end, steps in dmu_ranges:
        print(f"\nScanning dmu_{elem} from {start} to {end} ({steps} points) ...")
        gres, fres = [], []
        for mu in np.linspace(start, end, steps):
            g, f = run_point(model, dos, cfg, float(mu), elem)
            fres.append(f)
            if g:
                gres.append(g)
            done += 1
            if progress:
                progress(done, total)
        files = []
        if write_files:
            write_table(f"results_{base}.dat", fres, model, elem)
            files.append(f"results_{base}.dat")
            if gres:
                write_table(f"results_{base}_growth.dat", gres, model, elem)
                files.append(f"results_{base}_growth.dat")
        T = cfg["TEMPERATURE"]
        print(f"   {'dmu':>7} {'EF(eV)':>8} {'n':>10} {'p':>10}" + ("   [growth: EF, n]" if gres else ""))
        for i, r in enumerate(fres):
            extra = f"   [{gres[i]['EF']:.3f}, {gres[i]['n']:.2e}]" if gres else ""
            print(f"   {r['dmu']:7.3f} {r['EF']:8.3f} {r['n']:10.2e} {r['p']:10.2e}{extra}")
        if plot and steps > 1:
            panels = ([(f"(a) {cfg['FREEZE_TEMPERATURE']:.0f} K, equilibrium", gres),
                       (f"(b) {T:.0f} K, frozen-in", fres)] if gres else
                      [(f"{T:.0f} K, equilibrium", fres)])
            make_plot(base, elem, panels, model, dos, cfg)
            files += [f"plot_{base}.png", f"plot_{base}.pdf"]
        if files:
            print("   Saved: " + ", ".join(files))
        out["files"] += files
        out["scans"].append({"base": base, "elem": elem,
                             "final": [_row(r, model) for r in fres],
                             "growth": [_row(r, model) for r in gres]})
    return out


def main():
    ap = argparse.ArgumentParser(description="Self-consistent defect and carrier concentrations.")
    ap.add_argument("input", nargs="?", help="input file")
    ap.add_argument("--template", action="store_true", help="print an annotated example input and exit")
    ap.add_argument("--check", action="store_true", help="only print DOS and site diagnostics")
    ap.add_argument("--no-plot", action="store_true", help="skip plots")
    a = ap.parse_args()
    if a.template:
        print(TEMPLATE)
        return
    if not a.input:
        ap.print_help()
        return
    run(a.input, check=a.check, plot=not a.no_plot)


if __name__ == "__main__":
    main()
