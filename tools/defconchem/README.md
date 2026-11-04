# defconchem

**Self-consistent Fermi level, defect and carrier concentrations from DFT formation energies — with freeze-in.**

defconchem solves charge neutrality

  p − n + Σᵢ qᵢ cᵢ = 0

for the Fermi level, given defect formation energies, a host density of states (or effective masses) and a temperature. It can freeze defect populations at a growth temperature and re-equilibrate charge states, carriers and E_F at a device temperature, and it handles defect complexes and chemical-potential scans.

It comes in two forms that share **the same Python engine** (`py/defconchem.py`):

* a **web app** (this repository, served with GitHub Pages) that runs entirely in the browser through [Pyodide](https://pyodide.org) — no installation, and **no file ever leaves your computer**;
* a **command-line tool** (`python3 py/defconchem.py input.txt`).

The web app also contains `engine.js`, a JavaScript port of the same engine (same input format and physics; agrees with the Python engine to ~10⁻¹⁰ on the test inputs). It is used automatically when the page is opened directly from disk (`file://`, where browsers block the Python worker), when the Python runtime cannot be downloaded, or with `?engine=js`. The status badge shows which engine is running (`v11.3` = Python, `v11.3-js` = JavaScript).

---

## Using the web app

1. Open the GitHub Pages URL of this repository. The first visit downloads the Python runtime with numpy and scipy (several tens of MB; the browser caches it, so later visits start in a few seconds).
2. Press **Load example (Cl:ZnO)** to see a complete input, or **Open input file** for your own (old `defconchem-fz*` inputs work too).
3. Set the host band structure — **effective masses** or a **DOS file** (DOSCAR, vasprun.xml or two-column text) — the defect-site cell, temperatures and the chemical-potential scan.
4. Describe defects, complexes and stoichiometry in the editor (syntax reference is built in).
5. **Check inputs** prints the diagnostics; **Run** solves charge neutrality, plots concentrations and E_F versus Δμ, and offers the result tables (`.dat`) and the generated input file for download.

If you give both a DOS file and effective masses, the one not selected is reported as an Nc/Nv ratio — a quick convergence check.

## Using the command line

```bash
python3 py/defconchem.py --template > input.txt   # annotated example input
python3 py/defconchem.py input.txt --check        # diagnostics only
python3 py/defconchem.py input.txt                # results_*.dat (+ plot_*.png/.pdf with matplotlib)
python3 py/check_dos.py DOSCAR                    # sanity-check a DOS (electron count, volume, edges)
```

Requirements: Python ≥ 3.8, numpy, scipy (matplotlib optional for plots).

## Input in brief

```
BANDGAP         3.30
EFFECTIVE_MASS  0.28 0.59          # m_e* m_h* [alpha_e alpha_h]   (or DOS_FILE + NELECT)
VOLUME_ANG3     47.67              # cell to which "sites per cell" refer (or POSCAR_FILE)
TEMPERATURE     300
FREEZE_TEMPERATURE 800
FREEZE_DEFECTS  ALL
DMU_RANGE  scan O  -3.19 0.0 41

DEFECT Cl_O 2                      # name, sites per cell
STATE  1  -0.39  1                 # q, Ef0 (eV at E_F = VBM, dmu = 0), degeneracy
STATE  0   3.43  1
END
COMPLEX Cl_O--V_Zn
COMP    Cl_O_q1 1  V_Zn_q-2 1
BINDING 2.09
G_COMPLEX 4
END
STOICH  Cl_O  O  0.5               # E_f += 0.5 * dmu_O   (n = atoms removed)
```

Formation energies follow E_f(q) = Ef0 + q·E_F + n·Δμ, with E_F measured from the VBM. Defect names are labels only; the physics comes from the sites, STOICH, Ef0 and COMP/CONTAINS entries. See `--template` (or *Syntax reference* in the app) for every keyword.

## Method

* **Equilibrium (growth T):** dilute-limit, grand-canonical concentrations c = N_sites · g · exp(−E_f / kT). A complex is a defect with E_f = Σ E_f(components) − E_bind; forming it does **not** deplete the isolated defects (the reservoir supplies atoms).
* **Freeze-in:** total concentrations of frozen defects (including atoms bound in complexes) are fixed at the growth temperature; at the final temperature only charge states, carriers and E_F re-equilibrate (complexes either frozen or re-formed with atom conservation).
* **Fixed totals** (`TOTAL_CONC`) for dopants with known concentration use the same conservation solver.
* **Carriers:** Fermi–Dirac integrals over the host DOS (both spins, scissor-shifted to `BANDGAP`) or over parabolic/Kane effective-mass bands; degenerate doping is handled.
* **Solvers:** bracketed root finding for E_F (Brent), log-space mass conservation; warnings when no neutral E_F exists or a defect exceeds 5 % of its sites.

Validated against an independent grand-canonical/freeze-in solver (agreement ≤ 1 %), and the browser and command-line versions give identical numbers.

## Pitfalls worth knowing (all detected by `--check` / `check_dos.py`)

* **VASP DOS from a run with NSW > 0 is divided by NSW** (the DOS is averaged with `KBLOCK`, default = NSW). Compute the DOS in a static run (`NSW = 0`) or set `KBLOCK = 1`, and always give `NELECT` — defconchem then renormalises the DOS and prints the factor.
* **Line 2 of DOSCAR is the volume per atom**; defconchem multiplies by the number of atoms.
* **Coarse k-meshes** turn a light conduction band into a Γ spike (unphysical equivalent mass); use effective masses from a band-structure fit, or an interpolated DOS. Keep the DOS energy step (NEDOS) below ~kT/3.
* **Band edges** found from a smeared DOS can be off by ~σ — set `DOS_VBM`/`DOS_CBM` from the eigenvalues.
* **Sites per cell** must count the sites in the cell you give (e.g. 2 O sites in 4-atom wurtzite ZnO).

## Limitations

Dilute limit; a single freeze temperature; temperature-independent formation energies and band gap (no vibrational entropy); rigid scissor; concentrations per volume (2D materials need per-area conversion); not for metals.

## Running locally

* **Double-click `index.html`**: works offline with the JavaScript engine (Plotly still loads from the CDN for the charts; the table and result files work without it).
* **Python engine locally**: serve the folder, e.g. `python3 -m http.server 8000`, then open `http://localhost:8000/`.

## Deploying on GitHub Pages

1. Push this folder to a repository (e.g. `cmku/defconchem`).
2. *Settings → Pages → Build and deployment → Deploy from a branch*, branch `main`, folder `/ (root)`.
3. The app is then at `https://<org>.github.io/<repo>/`. The empty `.nojekyll` file keeps GitHub from processing the site.

The page loads Pyodide and Plotly from cdn.jsdelivr.net. For an offline/self-hosted setup, serve a Pyodide distribution (with numpy, scipy, openblas) yourself and open the app with `?pyodide=<URL of that folder>`.

## Repository layout

```
index.html, style.css, app.js   web interface
worker.js                       runs the Python engine in a Web Worker (Pyodide)
engine.js                       JavaScript port of the engine (offline / fallback)
py/defconchem.py                engine + command-line tool
py/check_dos.py                 DOS sanity checker
examples/Cl_ZnO_EMA.input       example (Cl donors in ZnO)
```

## Citation

If you use defconchem, please cite the repository and, for the example data,
S. Charoenphon, A. Sripothongnack, S. Limpijumnong and P. Reunchan, *Thermodynamic and kinetic stability of Cl and F donors in ZnO*, Phys. Rev. Materials **10**, 084602 (2026).

## License

MIT — see `LICENSE`.
