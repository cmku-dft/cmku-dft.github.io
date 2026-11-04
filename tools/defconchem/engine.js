/* defconchem engine — JavaScript port of py/defconchem.py v11.3 (same input format, same physics).
   Exposes DefconEngine.run(inputText, files, {check}) -> {ok, log, err, result}. */
(function (root) {
  "use strict";
  const VERSION = "v11.3-js";
  const kB = 8.617333262e-5;

  class DefconError extends Error {}
  const die = (msg) => { throw new DefconError(msg); };

  /* ---------------- small numerics ---------------- */
  function trapz(y, x) {
    let s = 0;
    for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    return s;
  }
  const expit = (x) => (x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x)));
  function logsumexp(a) {
    let m = -Infinity;
    for (const v of a) if (v > m) m = v;
    if (m === -Infinity) return -Infinity;
    let s = 0;
    for (const v of a) s += Math.exp(v - m);
    return m + Math.log(s);
  }
  function interp(x0, x, y) {
    if (x0 <= x[0]) return y[0];
    if (x0 >= x[x.length - 1]) return y[y.length - 1];
    let lo = 0, hi = x.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x[mid] <= x0) lo = mid; else hi = mid; }
    return y[lo] + (y[hi] - y[lo]) * (x0 - x[lo]) / (x[hi] - x[lo]);
  }
  function median(a) { const s = [...a].sort((p, q) => p - q); const n = s.length; return n % 2 ? s[(n - 1) / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]); }
  // Brent's method (as scipy.optimize.brentq)
  function brentq(f, xa, xb, xtol = 2e-12, rtol = 8.881784197001252e-16, maxiter = 500) {
    let xpre = xa, xcur = xb, xblk = 0, fpre = f(xpre), fcur = f(xcur), fblk = 0, spre = 0, scur = 0;
    if (fpre * fcur > 0) die("brentq: f(a) and f(b) must have different signs");
    if (fpre === 0) return xpre;
    if (fcur === 0) return xcur;
    for (let i = 0; i < maxiter; i++) {
      if (fpre !== 0 && fcur !== 0 && Math.sign(fpre) !== Math.sign(fcur)) { xblk = xpre; fblk = fpre; spre = scur = xcur - xpre; }
      if (Math.abs(fblk) < Math.abs(fcur)) { xpre = xcur; xcur = xblk; xblk = xpre; fpre = fcur; fcur = fblk; fblk = fpre; }
      const delta = (xtol + rtol * Math.abs(xcur)) / 2, sbis = (xblk - xcur) / 2;
      if (fcur === 0 || Math.abs(sbis) < delta) return xcur;
      if (Math.abs(spre) > delta && Math.abs(fcur) < Math.abs(fpre)) {
        let stry;
        if (xpre === xblk) stry = -fcur * (xcur - xpre) / (fcur - fpre);
        else {
          const dpre = (fpre - fcur) / (xpre - xcur), dblk = (fblk - fcur) / (xblk - xcur);
          stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
        }
        if (2 * Math.abs(stry) < Math.min(Math.abs(spre), 3 * Math.abs(sbis) - delta)) { spre = scur; scur = stry; }
        else { spre = sbis; scur = sbis; }
      } else { spre = sbis; scur = sbis; }
      xpre = xcur; fpre = fcur;
      if (Math.abs(scur) > delta) xcur += scur; else xcur += sbis > 0 ? delta : -delta;
      fcur = f(xcur);
    }
    return xcur;
  }

  /* ---------------- python-like formatting ---------------- */
  const pe = (x, d = 2) => { if (!isFinite(x)) return String(x); const s = Number(x).toExponential(d); return s.replace(/e([+-])(\d)$/, "e$10$2"); };
  const pf = (x, d = 3) => Number(x).toFixed(d);
  const sg = (q) => (q >= 0 ? "+" : "") + q;
  const g = (x) => { const s = String(Number(x)); return s; };

  /* ---------------- input parsing ---------------- */
  function parseInput(text, log) {
    const cfg = { DOS_FILE: null, BANDGAP: null, SPIN: 1, POSCAR_FILE: "POSCAR", VOLUME_ANG3: null, DOS_VOLUME_ANG3: null,
      DOS_VBM: null, DOS_CBM: null, NSITES_CM3: null, TEMPERATURE: 300.0, FREEZE_TEMPERATURE: null, FREEZE_DEFECTS: [],
      PLOT_MIN: 1e6, PLOT_SPECIES: [], DOS_TOL: null, NELECT: null, EFFECTIVE_MASS: null, CARRIER_MODEL: null };
    const defects = {}, legacy = [], dmu = [], stoich = {}, totalConc = {};
    let cur = null, cx = null;
    const order = [];
    const num = (s) => { const v = Number(s); if (s === undefined || s === "" || Number.isNaN(v)) throw new DefconError(`cannot parse number '${s}'`); return v; };
    const lines = text.split(/\r?\n/);
    for (let ln = 1; ln <= lines.length; ln++) {
      const raw = lines[ln - 1];
      const line = raw.split("#")[0].trim();
      if (!line) continue;
      const parts = line.split(/\s+/);
      const key = parts[0].toUpperCase();
      try {
        if (cur) {
          if (key === "END") cur = null;
          else if (key === "STATE") {
            let q = null, e0 = null, gg = 1.0;
            for (const tok of parts.slice(1)) {
              if (tok.includes("=")) {
                const [k, v] = tok.split("=", 2); const kk = k.toLowerCase();
                if (kk === "q") q = parseInt(v, 10); else if (kk === "ef0") e0 = num(v); else if (kk === "g") gg = num(v);
              }
            }
            if (q === null) { q = parseInt(parts[1], 10); e0 = num(parts[2]); if (parts.length > 3) gg = num(parts[3]); if (Number.isNaN(q)) throw new DefconError("bad charge"); }
            if (e0 === null) throw new DefconError("STATE needs Ef0");
            cur.states.push({ q, Ef0: e0, g: gg });
          } else if (key === "CONTAINS") {
            cur.contains = {};
            for (let i = 1; i < parts.length; i += 2) cur.contains[parts[i]] = num(parts[i + 1]);
          } else throw new DefconError(`unexpected '${parts[0]}' inside DEFECT block (missing END?)`);
          continue;
        }
        if (cx) {
          if (key === "END") { legacy.push(cx); cx = null; }
          else if (key === "COMP") { cx.components = []; for (let i = 1; i < parts.length; i += 2) cx.components.push([parts[i], num(parts[i + 1])]); }
          else if (key === "BINDING") cx.E_bind = num(parts[1]);
          else if (key === "G_COMPLEX") cx.g = num(parts[1]);
          else if (key === "SITES") cx.sites = num(parts[1]);
          else throw new DefconError(`unexpected '${parts[0]}' inside COMPLEX block (missing END?)`);
          continue;
        }
        if (key === "DEFECT") {
          const name = parts[1];
          if (name in defects) throw new DefconError(`defect '${name}' defined twice`);
          cur = { name, sites: num(parts[2]), states: [], contains: {} };
          defects[name] = cur; order.push(name);
        } else if (key === "COMPLEX") cx = { name: parts[1], components: [], E_bind: 0, g: 1, sites: null };
        else if (key === "STOICH") { (stoich[parts[1]] = stoich[parts[1]] || {})[parts[2]] = num(parts[3]); }
        else if (key === "DMU_RANGE") dmu.push([parts[1], parts[2], num(parts[3]), num(parts[4]), parseInt(parts[5], 10)]);
        else if (key === "FREEZE_DEFECTS") cfg.FREEZE_DEFECTS.push(...parts.slice(1));
        else if (key === "PLOT_SPECIES") cfg.PLOT_SPECIES.push(...parts.slice(1));
        else if (key === "TOTAL_CONC") totalConc[parts[1]] = num(parts[2]);
        else if (key === "EFFECTIVE_MASS") {
          const v = parts.slice(1).map(num);
          if (v.length < 2) throw new DefconError("EFFECTIVE_MASS needs me mh [alpha_e alpha_h]");
          while (v.length < 4) v.push(0);
          cfg.EFFECTIVE_MASS = v;
        } else if (key === "CARRIER_MODEL") cfg.CARRIER_MODEL = parts[1].toUpperCase();
        else if (key in cfg) cfg[key] = (key === "DOS_FILE" || key === "POSCAR_FILE") ? parts[1] : num(parts[1]);
        else log(`  WARNING: line ${ln}: unknown keyword '${parts[0]}' ignored`);
      } catch (e) {
        if (e instanceof DefconError) die(`line ${ln}: ${e.message}`);
        die(`line ${ln}: cannot parse '${raw.trim()}' (${e.message})`);
      }
    }
    if (cur || cx) die("input ended inside a DEFECT/COMPLEX block (missing END)");
    if (cfg.BANDGAP === null) die("BANDGAP is required");
    if (!cfg.DOS_FILE && !cfg.EFFECTIVE_MASS) die("give DOS_FILE and/or EFFECTIVE_MASS");
    const cm = cfg.CARRIER_MODEL;
    if (cm === null) cfg.CARRIER_MODEL = cfg.DOS_FILE ? "DOS" : "EMA";
    else if (cm !== "DOS" && cm !== "EMA") die("CARRIER_MODEL must be DOS or EMA");
    else if (cm === "DOS" && !cfg.DOS_FILE) die("CARRIER_MODEL DOS needs DOS_FILE");
    else if (cm === "EMA" && !cfg.EFFECTIVE_MASS) die("CARRIER_MODEL EMA needs EFFECTIVE_MASS");

    for (const c of legacy) {
      if (!c.components.length) die(`COMPLEX ${c.name}: no COMP line`);
      let qtot = 0, e0 = -c.E_bind; const contains = {};
      for (const [sname, nu] of c.components) {
        const i = sname.lastIndexOf("_q");
        if (i < 0) die(`COMPLEX ${c.name}: component '${sname}' must look like <defect>_q<charge>, e.g. VZn_q-2`);
        const par = sname.slice(0, i), qs = parseInt(sname.slice(i + 2), 10);
        if (!(par in defects)) die(`COMPLEX ${c.name}: unknown defect '${par}'`);
        const st = defects[par].states.filter((s) => s.q === qs);
        if (!st.length) die(`COMPLEX ${c.name}: defect ${par} has no charge state ${qs}`);
        qtot += nu * qs; e0 += nu * st[0].Ef0; contains[par] = (contains[par] || 0) + nu;
      }
      const first = c.components[0][0].slice(0, c.components[0][0].lastIndexOf("_q"));
      let d = defects[c.name];
      if (!d) { d = { name: c.name, sites: c.sites || defects[first].sites, states: [], contains, legacy: true }; defects[c.name] = d; order.push(c.name); }
      else if (JSON.stringify(Object.entries(d.contains).sort()) !== JSON.stringify(Object.entries(contains).sort()))
        die(`COMPLEX ${c.name}: charge states must contain the same defects`);
      qtot = Math.round(qtot);
      if (d.states.some((s) => s.q === qtot)) die(`COMPLEX ${c.name}: charge state ${qtot} defined twice`);
      d.states.push({ q: qtot, Ef0: e0, g: c.g });
    }
    for (const name of order) {
      const d = defects[name];
      if (!d.states.length) die(`defect ${name} has no STATE lines`);
      for (const p in d.contains) if (!(p in defects) || Object.keys(defects[p].contains).length) die(`${name}: CONTAINS '${p}' must be a simple (non-complex) defect`);
      if (!(name in stoich) && Object.keys(d.contains).length) {
        const inh = {};
        for (const [p, nu] of Object.entries(d.contains)) for (const [el, n] of Object.entries(stoich[p] || {})) inh[el] = (inh[el] || 0) + nu * n;
        stoich[name] = inh;
      }
    }
    for (const p in totalConc) if (!(p in defects) || Object.keys(defects[p].contains).length) die(`TOTAL_CONC: '${p}' must be a simple defect`);
    for (const s in stoich) if (!(s in defects)) die(`STOICH refers to unknown defect '${s}'`);
    const fr = cfg.FREEZE_DEFECTS;
    if (cfg.FREEZE_TEMPERATURE) {
      if (!fr.length) die("FREEZE_TEMPERATURE is set but FREEZE_DEFECTS is empty. Use 'FREEZE_DEFECTS ALL' or list names (or remove FREEZE_TEMPERATURE).");
      if (fr.length === 1 && fr[0].toUpperCase() === "ALL") cfg.FREEZE_DEFECTS = [...order];
      for (const x of cfg.FREEZE_DEFECTS) if (!(x in defects)) die(`FREEZE_DEFECTS: unknown name '${x}' (names: ${order.join(", ")})`);
    } else if (fr.length) { log("  WARNING: FREEZE_DEFECTS given without FREEZE_TEMPERATURE - ignored (full equilibrium)"); cfg.FREEZE_DEFECTS = []; }
    if (cfg.NSITES_CM3 !== null) log("  WARNING: NSITES_CM3 is no longer used: complex site densities come from SITES / component sites");
    return { cfg, defects, order, dmu, stoich, totalConc };
  }

  /* ---------------- DOS readers ---------------- */
  const toText = (data) => (typeof data === "string" ? data : new TextDecoder().decode(data));
  const rowsOf = (lines) => lines.map((l) => l.trim().split(/\s+/).map(Number));

  function readDoscar(text) {
    const L = text.split(/\r?\n/);
    const nions = parseInt(L[0].trim().split(/\s+/)[0], 10);
    const vol = Number(L[1].trim().split(/\s+/)[0]) * nions;
    const h = L[5].trim().split(/\s+/);
    const nedos = parseInt(Number(h[2]), 10), ef = Number(h[3]);
    const d = rowsOf(L.slice(6, 6 + nedos));
    if (d[0].length === 5) return { E: d.map((r) => r[0]), D: d.map((r) => r[1] + r[2]), I: d.map((r) => r[3] + r[4]), vol, ef, ns: 2 };
    return { E: d.map((r) => r[0]), D: d.map((r) => r[1]), I: d.map((r) => r[2]), vol, ef, ns: 1 };
  }
  function readVasprun(text) {
    const efm = text.match(/<dos>[\s\S]*?<i\s+name="efermi">\s*([-\d.Ee+]+)\s*<\/i>/);
    const totM = text.match(/<total>([\s\S]*?)<\/total>/);
    if (!efm || !totM) die("could not find the total DOS in vasprun.xml");
    const sets = [...totM[1].matchAll(/<set comment="spin\s*\d+">([\s\S]*?)<\/set>/g)];
    const arrs = sets.map((s) => [...s[1].matchAll(/<r>([^<]*)<\/r>/g)].map((r) => r[1].trim().split(/\s+/).map(Number)));
    const E = arrs[0].map((r) => r[0]);
    const D = E.map((_, i) => arrs.reduce((a, A) => a + A[i][1], 0));
    const I = E.map((_, i) => arrs.reduce((a, A) => a + A[i][2], 0));
    const vols = [...text.matchAll(/<i\s+name="volume">\s*([-\d.Ee+]+)\s*<\/i>/g)];
    if (!vols.length) die("could not find cell volume in vasprun.xml");
    const ne = text.match(/<i[^>]*name="NELECT"[^>]*>\s*([-\d.Ee+]+)\s*<\/i>/);
    return { E, D, I, vol: Number(vols[vols.length - 1][1]), ef: Number(efm[1]), ns: arrs.length, nelect: ne ? Number(ne[1]) : null };
  }
  function poscarVolume(text) {
    const L = text.split(/\r?\n/);
    const s = Number(L[1].trim().split(/\s+/)[0]);
    const M = [2, 3, 4].map((i) => L[i].trim().split(/\s+/).slice(0, 3).map(Number));
    const v = Math.abs(M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]));
    return s < 0 ? -s : v * s ** 3;
  }

  class DOS {
    constructor(cfg, siteVol, files, log) {
      this.log = log;
      const f = cfg.DOS_FILE;
      if (!(f in files)) die(`DOS file '${f}' not found`);
      const name = f.toLowerCase();
      this.nelect = null; this.idos = null;
      let E, D, I = null, vol, ef;
      if (name.endsWith(".xml")) {
        const r = readVasprun(toText(files[f])); ({ E, D, I, vol, ef } = r); this.nelect = r.nelect;
        this.source = `vasprun.xml (ISPIN=${r.ns}, spins summed)`; this.idos = I;
      } else if (name.includes("doscar")) {
        const r = readDoscar(toText(files[f])); ({ E, D, I, vol, ef } = r);
        this.source = `DOSCAR (ISPIN=${r.ns}, spins summed)`; this.idos = I;
      } else {
        const rows = rowsOf(toText(files[f]).split(/\r?\n/).filter((l) => l.trim() && !/^\s*#/.test(l)));
        E = rows.map((r) => r[0]);
        const sp2 = parseInt(cfg.SPIN, 10) === 2 && rows[0].length > 2;
        D = rows.map((r) => r[1] + (sp2 ? r[2] : 0));
        ef = 0;
        vol = cfg.DOS_VOLUME_ANG3;
        this.source = "text file (assumed: energies relative to VBM, DOS per cell, " +
          (parseInt(cfg.SPIN, 10) === 2 ? "columns up+down summed)" : "total DOS incl. both spins)");
        if (vol === null) { vol = siteVol; log("  WARNING: text DOS without DOS_VOLUME_ANG3: assuming the DOS was computed in the VOLUME_ANG3/POSCAR cell. Check this - a mismatch scales n and p directly!"); }
        if (cfg.DOS_VBM === null) cfg.DOS_VBM = 0.0;
      }
      this.vol = vol;
      const Dmax = Math.max(...D);
      const tol = cfg.DOS_TOL || 1e-3 * Dmax;
      let vbm, cbm;
      if (cfg.DOS_VBM !== null) vbm = cfg.DOS_VBM;
      else {
        let m = -Infinity; E.forEach((e, i) => { if (e <= ef + 0.05 && D[i] > tol && e > m) m = e; });
        if (m === -Infinity) die("could not find the VBM in the DOS; set DOS_VBM"); vbm = m;
      }
      if (cfg.DOS_CBM !== null) cbm = cfg.DOS_CBM;
      else {
        const gs = E.findIndex((e, i) => e > vbm && D[i] <= tol);
        if (gs < 0) die("no gap found in the DOS; set DOS_VBM / DOS_CBM");
        let m = Infinity; E.forEach((e, i) => { if (e > E[gs] && D[i] > tol && e < m) m = e; }); cbm = m;
      }
      this.vbm_abs = vbm; this.cbm_abs = cbm;
      this.auto_edges = cfg.DOS_VBM === null || cfg.DOS_CBM === null;
      this.eg_dos = cbm - vbm; this.eg = Number(cfg.BANDGAP); this.shift = this.eg - this.eg_dos;
      const vcm = vol * 1e-24;
      this.Ev = []; this.Dv = []; this.Ec = []; this.Dc = [];
      const Ev0 = [], Dv0 = [];
      E.forEach((e, i) => {
        if (e <= vbm) { this.Ev.push(e - vbm); this.Dv.push(D[i] / vcm); Ev0.push(e); Dv0.push(D[i]); }
        if (e >= cbm) { this.Ec.push(e - vbm + this.shift); this.Dc.push(D[i] / vcm); }
      });
      this.nval = this.idos ? interp(vbm, E, this.idos) : trapz(Dv0, Ev0);
      const dif = []; for (let i = 1; i < E.length; i++) dif.push(E[i] - E[i - 1]);
      this.de = median(dif);
      this.scale = 1.0;
      const ne = cfg.NELECT;
      if (ne !== null) {
        if (this.nelect === null) this.nelect = ne;
        if (Math.abs(this.nval - ne) > 0.02 * ne) {
          this.scale = ne / this.nval;
          this.Dv = this.Dv.map((v) => v * this.scale); this.Dc = this.Dc.map((v) => v * this.scale);
          log(`  WARNING: DOS renormalised x${Number(this.scale.toPrecision(4))} so that the VB holds NELECT = ${g(ne)} electrons (file gave ${Number(this.nval.toPrecision(4))}). Check how the DOS file was made!`);
          this.nval = ne;
        }
      } else if (this.nval < 1.0) {
        log(`  WARNING: VB holds only ${Number(this.nval.toPrecision(3))} electrons per cell - the DOS normalisation is almost certainly wrong. Add NELECT <valence electrons of the DOS cell> to renormalise.`);
      } else if (this.source.startsWith("text")) {
        log("  WARNING: text DOS without NELECT: normalisation cannot be checked. Add NELECT <valence electrons of the DOS cell> (e.g. 36 for 4-atom ZnO with Zn 3d)");
      }
    }
    carriers(EF, T) {
      const kT = kB * T;
      let n = 0, p = 0;
      const Ec = this.Ec, Dc = this.Dc, Ev = this.Ev, Dv = this.Dv;
      let prev = Dc[0] * expit(-(Ec[0] - EF) / kT);
      for (let i = 1; i < Ec.length; i++) { const cur = Dc[i] * expit(-(Ec[i] - EF) / kT); n += 0.5 * (cur + prev) * (Ec[i] - Ec[i - 1]); prev = cur; }
      prev = Dv[0] * expit((Ev[0] - EF) / kT);
      for (let i = 1; i < Ev.length; i++) { const cur = Dv[i] * expit((Ev[i] - EF) / kT); p += 0.5 * (cur + prev) * (Ev[i] - Ev[i - 1]); prev = cur; }
      return [n, p];
    }
    effN(T) {
      const kT = kB * T;
      const Nc = trapz(this.Dc.map((d, i) => d * Math.exp(-(this.Ec[i] - this.eg) / kT)), this.Ec);
      const Nv = trapz(this.Dv.map((d, i) => d * Math.exp(this.Ev[i] / kT)), this.Ev);
      const N1 = 2.5094e19 * (T / 300) ** 1.5;
      return [Nc, Nv, (Nc / N1) ** (2 / 3), (Nv / N1) ** (2 / 3)];
    }
    report(Ts) {
      const log = this.log;
      log(`  DOS source        : ${this.source}`);
      log(`  DOS cell volume   : ${pf(this.vol)} A^3`);
      if (this.scale !== 1.0) log(`  DOS rescaled by   : x${Number(this.scale.toPrecision(4))} (NELECT)`);
      log(`  VBM, CBM in file  : ${pf(this.vbm_abs)}, ${pf(this.cbm_abs)} eV  (DOS gap ${pf(this.eg_dos)} eV)`);
      log(`  BANDGAP used      : ${pf(this.eg)} eV  -> CB scissor shift ${this.shift >= 0 ? "+" : ""}${pf(this.shift)} eV`);
      let line = `  Electrons in VB   : ${pf(this.nval, 2)} per DOS cell`;
      if (this.nelect !== null) {
        line += `  (NELECT = ${pf(this.nelect, 2)})`;
        if (Math.abs(this.nval - this.nelect) > 0.05 * this.nelect) line += "  <-- MISMATCH: check DOS normalisation / VBM";
      }
      log(line + "   [should equal the valence electron count of the DOS cell]");
      log(`  DOS energy step   : ${pf(this.de * 1000, 1)} meV`);
      for (const T of Ts) {
        const [Nc, Nv, mc, mv] = this.effN(T);
        log(`  T = ${pf(T, 0).padStart(6)} K     : Nc = ${pe(Nc)}, Nv = ${pe(Nv)} cm^-3  (equivalent DOS masses m_c* = ${pf(mc, 2)}, m_v* = ${pf(mv, 2)})`);
        if (this.de > kB * T) log(`  WARNING: DOS energy step (${pf(this.de * 1000, 0)} meV) > kT at ${pf(T, 0)} K - carrier integrals may be inaccurate (increase NEDOS)`);
      }
      if (this.auto_edges && !this.source.startsWith("text"))
        log("  (band edges found from a DOS threshold; smearing tails can shift them by ~sigma. For accuracy set\n   DOS_VBM / DOS_CBM to the VBM/CBM eigenvalues, e.g. VASPKIT 911 or EIGENVAL, same energy scale as the DOS)");
      if (Math.abs(this.shift) > 0.02) log("  (CB states shifted rigidly; make sure Ef0 values use the same band edges)");
    }
  }

  class EMA {
    constructor(cfg, Tmax, log) {
      this.log = log;
      const [me, mh, ae, ah] = cfg.EFFECTIVE_MASS;
      this.eg = Number(cfg.BANDGAP);
      const emax = Math.max(1.5, 50 * kB * Tmax), N = 20001;
      const x = Array.from({ length: N }, (_, i) => (emax * i) / (N - 1));
      const gfun = (xx, m, a) => (1 / (2 * Math.PI ** 2)) * (m / 3.80998) ** 1.5 * Math.sqrt(xx * (1 + a * xx)) * (1 + 2 * a * xx) * 1e24;
      this.Ec = x.map((v) => this.eg + v); this.Dc = x.map((v) => gfun(v, me, ae));
      this.Ev = x.map((v) => -v).reverse(); this.Dv = x.map((v) => gfun(v, mh, ah)).reverse();
      this.source = `effective-mass model: m_e* = ${g(me)}, m_h* = ${g(mh)}` + (ae || ah ? `, alpha_e = ${g(ae)}, alpha_h = ${g(ah)} 1/eV (Kane)` : " (parabolic)");
    }
    report(Ts) {
      this.log(`  Carrier model     : ${this.source}`);
      this.log(`  BANDGAP used      : ${pf(this.eg)} eV`);
      for (const T of Ts) {
        const [Nc, Nv, mc, mv] = this.effN(T);
        this.log(`  T = ${pf(T, 0).padStart(6)} K     : Nc = ${pe(Nc)}, Nv = ${pe(Nv)} cm^-3  (equivalent DOS masses m_c* = ${pf(mc, 2)}, m_v* = ${pf(mv, 2)})`);
      }
    }
  }
  EMA.prototype.carriers = DOS.prototype.carriers;
  EMA.prototype.effN = DOS.prototype.effN;

  /* ---------------- defect model ---------------- */
  class Model {
    constructor(P, siteVol) {
      this.vcm = siteVol * 1e-24;
      this.D = P.defects; this.order = P.order;
      this.simple = P.order.filter((n) => !Object.keys(P.defects[n].contains).length);
      this.cx = P.order.filter((n) => Object.keys(P.defects[n].contains).length);
      this.stoich = P.stoich; this.totalConc = P.totalConc;
      this.N = {}; for (const n of P.order) this.N[n] = P.defects[n].sites / this.vcm;
    }
    species() { const s = []; for (const n of this.order) for (const st of this.D[n].states) s.push([n, st.q]); return s; }
    logcEq(EF, T, dmu, elem) {
      const kT = kB * T, out = {};
      for (const n of this.order) {
        const nmu = (this.stoich[n] || {})[elem] || 0, lN = Math.log(this.N[n]);
        out[n] = new Map(this.D[n].states.map((s) => [s.q, lN + Math.log(s.g) - (s.Ef0 + s.q * EF + nmu * dmu) / kT]));
      }
      return out;
    }
    concentrations(EF, T, dmu, elem, fixedTot = {}, fixedCx = {}, log = () => {}) {
      const lc = this.logcEq(EF, T, dmu, elem);
      const logZ = {}; for (const n of this.order) logZ[n] = logsumexp([...lc[n].values()]);
      const la = {}; for (const n of this.order) la[n] = 0;
      for (const [n, tot] of Object.entries(fixedCx)) la[n] = Math.log(Math.max(tot, 1e-300)) - logZ[n];
      const cons = Object.keys(fixedTot).filter((p) => !(p in fixedCx));
      const freeCx = this.cx.filter((c) => !(c in fixedCx));
      const sumLa = (c) => Object.entries(this.D[c].contains).reduce((a, [p, nu]) => a + nu * (la[p] || 0), 0);
      if (cons.length) {
        const rhs = {};
        for (const p of cons) {
          let r = fixedTot[p] - this.cx.filter((c) => c in fixedCx).reduce((a, c) => a + (this.D[c].contains[p] || 0) * fixedCx[c], 0);
          if (r <= 0) { log(`  WARNING: atom balance for ${p} non-positive (${pe(r, 3)}); set to tiny`); r = 1e-300; }
          rhs[p] = Math.log(r);
        }
        const lcx = (c) => logZ[c] + sumLa(c);
        let conv = false;
        for (let it = 0; it < 200; it++) {
          const old = { ...la };
          for (const p of cons) {
            const links = freeCx.filter((c) => p in this.D[c].contains).map((c) => [c, this.D[c].contains[p]]);
            const f = (x) => { la[p] = x; return logsumexp([x + logZ[p], ...links.map(([c, nu]) => Math.log(nu) + lcx(c))]) - rhs[p]; };
            la[p] = brentq(f, -3000, 3000, 1e-12);
          }
          if (Math.max(...cons.map((p) => Math.abs(la[p] - old[p]))) < 1e-10) { conv = true; break; }
        }
        if (!conv) log("  WARNING: mass-conservation solver did not fully converge");
      }
      for (const c of freeCx) la[c] = sumLa(c);
      const conc = {};
      for (const n of this.order) { conc[n] = new Map(); for (const [q, v] of lc[n]) conc[n].set(q, Math.exp(Math.min(v + la[n], 700))); }
      return conc;
    }
    charge(conc) { let s = 0; for (const n of this.order) for (const [q, c] of conc[n]) s += q * c; return s; }
    atomTotals(conc) {
      const tot = {};
      for (const p of this.simple) tot[p] = [...conc[p].values()].reduce((a, b) => a + b, 0);
      for (const c of this.cx) { const ct = [...conc[c].values()].reduce((a, b) => a + b, 0); for (const [p, nu] of Object.entries(this.D[c].contains)) tot[p] += nu * ct; }
      return tot;
    }
  }

  function solveEF(model, dos, T, dmu, elem, fixedTot, fixedCx, log) {
    const Q = (EF) => { const [n, p] = dos.carriers(EF, T); return p - n + model.charge(model.concentrations(EF, T, dmu, elem, fixedTot, fixedCx, log)); };
    let lo = -1.0, hi = dos.eg + 1.0, qlo = Q(lo), qhi = Q(hi), k = 0;
    while (qlo * qhi > 0 && k < 4) { lo -= 0.5; hi += 0.5; qlo = Q(lo); qhi = Q(hi); k++; }
    let EF, flag = "";
    if (qlo * qhi > 0) {
      EF = Math.abs(qlo) < Math.abs(qhi) ? lo : hi; flag = "NO_ROOT";
      log(`  WARNING: dmu=${dmu >= 0 ? "+" : ""}${pf(dmu)}, T=${pf(T, 0)} K: no charge-neutral E_F in [${pf(lo, 1)}, ${pf(hi, 1)}] eV - reporting edge value ${pf(EF, 2)} eV (check inputs)`);
    } else EF = brentq(Q, lo, hi, 1e-8, 1e-12, 500);
    const conc = model.concentrations(EF, T, dmu, elem, fixedTot, fixedCx, log);
    const [n, p] = dos.carriers(EF, T);
    for (const nm of model.order) {
      const s = [...conc[nm].values()].reduce((a, b) => a + b, 0);
      if (s > 0.05 * model.N[nm]) log(`  WARNING: dmu=${dmu >= 0 ? "+" : ""}${pf(dmu)}: [${nm}] > 5% of its sites - dilute-limit formulas are unreliable`);
    }
    return { dmu, EF, n, p, conc, tot: model.atomTotals(conc), flag };
  }

  function runPoint(model, dos, cfg, dmu, elem, log) {
    const T = cfg.TEMPERATURE, Tg = cfg.FREEZE_TEMPERATURE;
    const tc = Object.keys(model.totalConc).length ? { ...model.totalConc } : {};
    if (!Tg) return [null, solveEF(model, dos, T, dmu, elem, tc, {}, log)];
    const gr = solveEF(model, dos, Tg, dmu, elem, tc, {}, log);
    const frozen = new Set(cfg.FREEZE_DEFECTS);
    const fixedTot = {}, fixedCx = {};
    for (const p of model.simple) {
      if (frozen.has(p)) fixedTot[p] = gr.tot[p];
      else if (p in model.totalConc) fixedTot[p] = model.totalConc[p];
    }
    for (const c of model.cx) if (frozen.has(c)) fixedCx[c] = [...gr.conc[c].values()].reduce((a, b) => a + b, 0);
    const freeCx = model.cx.filter((c) => !(c in fixedCx));
    for (const p of Object.keys(fixedTot)) {
      if (frozen.has(p) && !freeCx.some((c) => p in model.D[c].contains)) {
        fixedCx[p] = [...gr.conc[p].values()].reduce((a, b) => a + b, 0);
        delete fixedTot[p];
      }
    }
    return [gr, solveEF(model, dos, T, dmu, elem, fixedTot, fixedCx, log)];
  }

  function row(r, model) {
    const conc = {}; for (const [n, q] of model.species()) conc[`${n}_q${q}`] = r.conc[n].get(q);
    const tot = {}; for (const p of model.simple) tot[p] = r.tot[p];
    for (const c of model.cx) tot[c] = [...r.conc[c].values()].reduce((a, b) => a + b, 0);
    return { dmu: r.dmu, EF: r.EF, n: r.n, p: r.p, flag: r.flag, conc, tot };
  }

  function table(res, model, elem) {
    const sp = model.species();
    const hdr = [`dmu_${elem}`, "EF", "n", "p", ...sp.map(([n, q]) => `${n}_q${q}`), ...model.simple.map((p) => `Tot_${p}`), ...model.cx.map((c) => `Tot_${c}`)];
    const L = ["# E_F in eV above VBM; concentrations in cm^-3. Tot_<simple> includes atoms bound in complexes.", "# " + hdr.join("\t")];
    for (const r of res) {
      const vals = [r.dmu, r.EF, r.n, r.p, ...sp.map(([n, q]) => r.conc[n].get(q)), ...model.simple.map((p) => r.tot[p]),
        ...model.cx.map((c) => [...r.conc[c].values()].reduce((a, b) => a + b, 0))];
      L.push(vals.map((x) => pe(x, 5)).join("\t") + (r.flag ? "\t# " + r.flag : ""));
    }
    return L.join("\n") + "\n";
  }

  function run(inputText, files = {}, opts = {}) {
    const out = [], logf = (s) => out.push(s);
    try {
      const P = parseInput(inputText, logf);
      const cfg = P.cfg;
      let siteVol;
      if (cfg.VOLUME_ANG3) siteVol = cfg.VOLUME_ANG3;
      else if (cfg.POSCAR_FILE in files) siteVol = poscarVolume(toText(files[cfg.POSCAR_FILE]));
      else die("give VOLUME_ANG3 or a POSCAR_FILE for the defect-site cell");
      logf("");
      logf(`=== defconchem ${VERSION} ===`);
      const Ts = [cfg.TEMPERATURE, ...(cfg.FREEZE_TEMPERATURE ? [cfg.FREEZE_TEMPERATURE] : [])];
      const dosFile = cfg.DOS_FILE ? new DOS(cfg, siteVol, files, logf) : null;
      const ema = cfg.EFFECTIVE_MASS ? new EMA(cfg, Math.max(...Ts), logf) : null;
      const dos = cfg.CARRIER_MODEL === "DOS" ? dosFile : ema;
      const model = new Model(P, siteVol);
      logf(`  Carriers from     : ${dos === dosFile ? "DOS file" : "effective-mass model"}  (CARRIER_MODEL ${cfg.CARRIER_MODEL})`);
      dos.report(Ts);
      if (dosFile && ema) {
        const other = dos === dosFile ? ema : dosFile;
        logf(other === dosFile ? "  --- DOS file (for comparison only) ---" : `  --- ${ema.source} (for comparison only) ---`);
        for (const T of Ts) {
          const [a1, b1] = dos.effN(T), [a2, b2] = other.effN(T);
          logf(`  T = ${pf(T, 0).padStart(6)} K     : Nc ratio (used/other) = ${pf(a1 / a2, 2)}, Nv ratio = ${pf(b1 / b2, 2)}` +
            (Math.abs(a1 / a2 - 1) > 0.3 ? "   <-- differ by >30%: DOS not converged or masses inconsistent" : ""));
        }
      }
      logf(`  Defect-site cell  : ${pf(siteVol)} A^3`);
      for (const n of P.order) {
        const d = P.defects[n];
        const kind = Object.keys(d.contains).length ? `complex of ${Object.entries(d.contains).map(([p, v]) => `${g(v)} ${p}`).join(" + ")}` : "defect";
        logf(`    ${n.padEnd(12)} ${kind.padEnd(28)} sites ${pe(model.N[n], 3)} cm^-3   q = ${d.states.map((s) => sg(s.q)).join(", ")}`);
      }
      if (cfg.FREEZE_TEMPERATURE) logf(`  Freeze-in         : grown at ${pf(cfg.FREEZE_TEMPERATURE, 0)} K, evaluated at ${pf(cfg.TEMPERATURE, 0)} K; frozen: ${cfg.FREEZE_DEFECTS.join(", ")}`);
      else logf(`  Full equilibrium at ${pf(cfg.TEMPERATURE, 0)} K`);
      if (Object.keys(P.totalConc).length) logf("  Fixed totals      : " + Object.entries(P.totalConc).map(([k, v]) => `${k} = ${pe(v)}`).join(", "));
      const result = { version: VERSION, species: model.species().map(([n, q]) => `${n}_q${q}`), simple: model.simple, complexes: model.cx,
        bandgap: dos.eg, T: cfg.TEMPERATURE, T_growth: cfg.FREEZE_TEMPERATURE, plot_min: cfg.PLOT_MIN, scans: [], files: {} };
      if (opts.check) return { ok: true, log: out.join("\n"), err: "", result };
      let ranges = P.dmu;
      if (!ranges.length) { logf("  No DMU_RANGE: single point at dmu = 0"); ranges = [["run", "X", 0, 0, 1]]; }
      for (const [base, elem, start, end, steps] of ranges) {
        logf(""); logf(`Scanning dmu_${elem} from ${g(start)} to ${g(end)} (${steps} points) ...`);
        const gres = [], fres = [];
        for (let i = 0; i < steps; i++) {
          const mu = steps === 1 ? start : start + ((end - start) * i) / (steps - 1);
          const [gr, fr] = runPoint(model, dos, cfg, mu, elem, logf);
          fres.push(fr); if (gr) gres.push(gr);
        }
        result.files[`results_${base}.dat`] = table(fres, model, elem);
        if (gres.length) result.files[`results_${base}_growth.dat`] = table(gres, model, elem);
        logf(`   ${"dmu".padStart(7)} ${"EF(eV)".padStart(8)} ${"n".padStart(10)} ${"p".padStart(10)}` + (gres.length ? "   [growth: EF, n]" : ""));
        fres.forEach((r, i) => {
          const extra = gres.length ? `   [${pf(gres[i].EF)}, ${pe(gres[i].n)}]` : "";
          logf(`   ${pf(r.dmu).padStart(7)} ${pf(r.EF).padStart(8)} ${pe(r.n).padStart(10)} ${pe(r.p).padStart(10)}${extra}`);
        });
        result.scans.push({ base, elem, final: fres.map((r) => row(r, model)), growth: gres.map((r) => row(r, model)) });
      }
      return { ok: true, log: out.join("\n"), err: "", result };
    } catch (e) {
      const msg = e instanceof DefconError ? `ERROR: ${e.message}` : `ERROR: ${e && e.message ? e.message : e}`;
      return { ok: false, log: out.join("\n"), err: msg, result: null };
    }
  }

  root.DefconEngine = { run, version: VERSION };
  if (typeof module !== "undefined") module.exports = root.DefconEngine;
})(typeof self !== "undefined" ? self : globalThis);
