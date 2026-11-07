/* defconchem — Grace (xmgrace) project export.
   Builds an .agr file from the plotted results: one column per selected panel
   (e.g. frozen-in / growth), with a log-scale concentration graph above an E_F graph.
   Everything in the file (axis labels, ticks, line styles, legends, view boxes) can be
   edited afterwards in xmgrace. Exposes window.DefconAgr.build(spec) -> string. */
"use strict";
(function () {
  // Colours 16+ are added to Grace's default map; they follow the web-app palette.
  const PAL = [[42, 120, 214], [235, 104, 52], [27, 175, 122], [237, 161, 0],
    [232, 123, 164], [0, 131, 0], [98, 80, 214], [227, 73, 72]];
  const BASE_COLORS = [[255, 255, 255, "white"], [0, 0, 0, "black"], [255, 0, 0, "red"], [0, 255, 0, "green"],
    [0, 0, 255, "blue"], [255, 255, 0, "yellow"], [188, 143, 143, "brown"], [220, 220, 220, "grey"],
    [148, 0, 211, "violet"], [0, 255, 255, "cyan"], [255, 0, 255, "magenta"], [255, 165, 0, "orange"],
    [114, 33, 188, "indigo"], [103, 7, 72, "maroon"], [64, 224, 208, "turquoise"], [0, 139, 0, "green4"]];
  const GREY = 17 + PAL.length;          // colour index of the mid-grey used for band edges / p
  // Grace line styles: 1 solid, 2 dotted, 3 dashed, 4 long dashed, 5 dot-dashed
  const LSTYLE = [1, 3, 2, 5];

  const q = (s) => `"${String(s).replace(/"/g, "'")}"`;
  const f6 = (x) => Number(x).toFixed(6);

  /* Species name such as "V_O_q2" or "Cl_O--V_Zn_q-1" -> Grace string with sub/superscripts. */
  function speciesLabel(sp) {
    const i = sp.lastIndexOf("_q");
    const name = i < 0 ? sp : sp.slice(0, i), qq = i < 0 ? null : parseInt(sp.slice(i + 2), 10);
    const part = (p) => { const j = p.indexOf("_"); return j < 0 ? p : `${p.slice(0, j)}\\s${p.slice(j + 1).replace(/_/g, "")}\\N`; };
    const nm = name.split(/--|-(?=[A-Z])/).map(part).join("-");
    if (qq === null || Number.isNaN(qq)) return nm;
    const sup = qq === 0 ? "0" : `${Math.abs(qq) > 1 ? Math.abs(qq) : ""}${qq > 0 ? "+" : "-"}`;
    return `${nm}\\S${sup}\\N`;
  }

  function niceStep(span, target) {
    if (!(span > 0)) return 1;
    const raw = span / target, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
  }

  function header(lines) {
    lines.push("# Grace project file", "# written by the defconchem web app", "#", "@version 50125",
      "@page size 792, 612", "@page scroll 5%", "@page inout 5%", "@link page off",
      '@map font 0 to "Times-Roman", "Times-Roman"', '@map font 1 to "Times-Italic", "Times-Italic"',
      '@map font 2 to "Times-Bold", "Times-Bold"', '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"',
      '@map font 4 to "Helvetica", "Helvetica"', '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"',
      '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"', '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"',
      '@map font 12 to "Symbol", "Symbol"');
    BASE_COLORS.forEach(([r, g, b, n], i) => lines.push(`@map color ${i} to (${r}, ${g}, ${b}), "${n}"`));
    PAL.forEach(([r, g, b], i) => lines.push(`@map color ${16 + i} to (${r}, ${g}, ${b}), "dcc${i + 1}"`));
    lines.push(`@map color ${16 + PAL.length} to (0, 0, 0), "dccblack"`, `@map color ${GREY} to (128, 128, 128), "dccgrey"`,
      "@reference date 0", "@date wrap off", "@date wrap year 1950", "@default linewidth 1.0", "@default linestyle 1",
      "@default color 1", "@default pattern 1", "@default font 0", "@default char size 1.000000",
      "@default symbol size 1.000000", '@default sformat "%.8g"', "@background color 0", "@page background fill on",
      "@timestamp off");
  }

  function axis(lines, ax, o) {
    const p = `@    ${ax}axis `;
    lines.push(`${p} on`, `${p} type zero false`, `${p} offset 0.000000 , 0.000000`, `${p} bar on`, `${p} bar color 1`,
      `${p} bar linestyle 1`, `${p} bar linewidth 1.5`,
      `${p} label ${q(o.label || "")}`, `${p} label layout para`, `${p} label place auto`,
      `${p} label char size ${f6(o.labelSize || 1.3)}`, `${p} label font 0`, `${p} label color 1`, `${p} label place normal`,
      `${p} tick ${o.ticks === false ? "off" : "on"}`, `${p} tick major ${o.major}`, `${p} tick minor ticks ${o.minor}`,
      `${p} tick default 6`, `${p} tick place rounded true`, `${p} tick in`,
      `${p} tick major size 1.000000`, `${p} tick major color 1`, `${p} tick major linewidth 1.5`, `${p} tick major linestyle 1`,
      `${p} tick major grid off`, `${p} tick minor color 1`, `${p} tick minor linewidth 1.0`, `${p} tick minor linestyle 1`,
      `${p} tick minor grid off`, `${p} tick minor size 0.500000`,
      `${p} ticklabel ${o.ticklabels === false ? "off" : "on"}`, `${p} ticklabel format ${o.format || "general"}`,
      `${p} ticklabel prec ${o.prec ?? 5}`, `${p} ticklabel formula ""`, `${p} ticklabel append ""`, `${p} ticklabel prepend ""`,
      `${p} ticklabel angle 0`, `${p} ticklabel skip 0`, `${p} ticklabel stagger 0`, `${p} ticklabel place normal`,
      `${p} ticklabel offset auto`, `${p} ticklabel offset 0.000000 , 0.010000`, `${p} ticklabel start type auto`,
      `${p} ticklabel start 0.000000`, `${p} ticklabel stop type auto`, `${p} ticklabel stop 0.000000`,
      `${p} ticklabel char size ${f6(o.tickSize || 1.1)}`, `${p} ticklabel font 0`, `${p} ticklabel color 1`,
      `${p} tick place both`, `${p} tick spec type none`);
  }

  function setLines(lines, i, s) {
    const p = `@    s${i}`;
    lines.push(`${p} hidden false`, `${p} type xy`, `${p} symbol 0`, `${p} symbol size 1.000000`, `${p} symbol color ${s.color}`,
      `${p} symbol fill pattern 0`, `${p} line type 1`, `${p} line linestyle ${s.style}`, `${p} line linewidth ${s.width.toFixed(1)}`,
      `${p} line color ${s.color}`, `${p} line pattern 1`, `${p} baseline type 0`, `${p} baseline off`, `${p} dropline off`,
      `${p} fill type 0`, `${p} avalue off`, `${p} errorbar off`, `${p} comment ${q(s.comment || "")}`, `${p} legend  ${q(s.legend || "")}`);
  }

  function graph(lines, g, G) {
    lines.push(`@g${g} on`, `@g${g} hidden false`, `@g${g} type XY`, `@g${g} stacked false`, `@g${g} bar hgap 0.000000`,
      `@g${g} fixedpoint off`, `@with g${g}`,
      `@    world ${G.world.join(", ")}`, `@    view ${G.view.map(f6).join(", ")}`,
      `@    title ${q(G.title || "")}`, "@    title font 0", "@    title size 1.200000", "@    title color 1",
      `@    subtitle ${q(G.subtitle || "")}`, "@    subtitle font 0", "@    subtitle size 1.000000", "@    subtitle color 1",
      "@    xaxes scale Normal", `@    yaxes scale ${G.logy ? "Logarithmic" : "Normal"}`, "@    xaxes invert off", "@    yaxes invert off");
    axis(lines, "x", G.x);
    axis(lines, "y", G.y);
    lines.push("@    altxaxis  off", "@    altyaxis  off",
      `@    legend ${G.legend ? "on" : "off"}`, "@    legend loctype view", `@    legend ${f6(G.legendAt[0])}, ${f6(G.legendAt[1])}`,
      "@    legend box color 1", "@    legend box pattern 1", "@    legend box linewidth 1.0", "@    legend box linestyle 0",
      "@    legend box fill color 0", "@    legend box fill pattern 0", "@    legend font 0", "@    legend char size 0.900000",
      "@    legend color 1", "@    legend length 3", "@    legend vgap 1", "@    legend hgap 1", "@    legend invert false",
      "@    frame type 0", "@    frame linestyle 1", "@    frame linewidth 1.5", "@    frame color 1", "@    frame pattern 1",
      "@    frame background color 0", "@    frame background pattern 0");
    G.sets.forEach((s, i) => setLines(lines, i, s));
  }

  /* spec = {
       xLabel, panels: [{ label, rows: [{ dmu, EF, n, p, conc: {sp: c} }] }],
       species: [...]  (in plotting order), allSpecies: [...] (fixes colours),
       bandgap, plotMin, conc: bool, carriers: bool, ef: bool, legend: bool, bandEdges: bool } */
  function build(spec) {
    const panels = spec.panels;
    const showConc = spec.conc && (spec.species.length > 0 || spec.carriers);
    const showEF = spec.ef;
    if (!panels.length || (!showConc && !showEF)) throw new Error("Nothing selected for the Grace file.");
    const all = spec.allSpecies || spec.species;
    const colorOf = (sp) => 16 + (all.indexOf(sp) % PAL.length);
    const styleOf = (sp) => LSTYLE[Math.floor(all.indexOf(sp) / PAL.length) % LSTYLE.length];

    // common x range and ticks
    const xs = panels.flatMap((p) => p.rows.map((r) => r.dmu));
    let x0 = Math.min(...xs), x1 = Math.max(...xs);
    if (x0 === x1) { x0 -= 0.5; x1 += 0.5; }
    const xmaj = niceStep(x1 - x0, 5);

    // common log-y range for the concentration graphs (same rule as the web chart)
    const pmin = spec.plotMin || 1e6;
    let ymax = 20;
    if (showConc) {
      const ys = panels.flatMap((p) => p.rows.flatMap((r) => [
        ...spec.species.map((sp) => r.conc[sp] || 0), ...(spec.carriers ? [r.n || 0, r.p || 0] : [])])).filter((v) => v > 0);
      if (ys.length) ymax = Math.ceil(Math.log10(Math.max(...ys)));
    }
    const ymin = Math.max(Math.floor(Math.log10(pmin)), ymax - 16);
    const eg = spec.bandgap;
    const efmaj = niceStep(eg + 0.3, 4);

    // page layout (landscape 792 x 612 -> view 1.294 x 1.0): one column per panel
    const ncol = panels.length, L = 0.16, R = 1.24, gap = 0.12;
    const cw = (R - L - gap * (ncol - 1)) / ncol;
    const top = 0.9, bottom = 0.13;
    let concBox = null, efBox = null;
    if (showConc && showEF) { concBox = [0.39, top]; efBox = [bottom, 0.33]; }
    else if (showConc) concBox = [bottom, top];
    else efBox = [bottom, top];

    const xLabel = spec.xLabel;
    const lines = [];
    header(lines);
    const data = [];
    let g = 0;
    panels.forEach((P, c) => {
      const vx0 = L + c * (cw + gap), vx1 = vx0 + cw;
      const x = P.rows.map((r) => r.dmu);
      if (showConc) {
        const sets = [];
        spec.species.forEach((sp) => sets.push({ color: colorOf(sp), style: styleOf(sp), width: 2, legend: speciesLabel(sp),
          comment: sp, y: P.rows.map((r) => r.conc[sp]) }));
        if (spec.carriers) {
          sets.push({ color: 1, style: 3, width: 2, legend: "n", comment: "electrons", y: P.rows.map((r) => r.n) });
          sets.push({ color: GREY, style: 2, width: 2, legend: "p", comment: "holes", y: P.rows.map((r) => r.p) });
        }
        graph(lines, g, {
          world: [x0, `1e${ymin}`, x1, `1e${ymax}`], view: [vx0, concBox[0], vx1, concBox[1]], logy: true,
          title: "", subtitle: P.label,
          x: { label: showEF ? "" : xLabel, major: xmaj, minor: 1, ticklabels: !showEF },
          y: { label: c === 0 ? "Concentration (cm\\S-3\\N)" : "", major: 10, minor: 8, format: "power", prec: 0 },
          legend: spec.legend && c === 0, legendAt: [vx1 - 0.2, concBox[1] - 0.02], sets });
        data.push({ g, sets, x, log: true });
        g++;
      }
      if (showEF) {
        const sets = [{ color: 1, style: 1, width: 2, legend: "", comment: "Fermi level", y: P.rows.map((r) => r.EF) }];
        if (spec.bandEdges) {
          sets.push({ color: GREY, style: 3, width: 1, legend: "", comment: "VBM", y: x.map(() => 0) });
          sets.push({ color: GREY, style: 3, width: 1, legend: "", comment: "CBM", y: x.map(() => eg) });
        }
        graph(lines, g, {
          world: [x0, -0.15, x1, +(eg + 0.15).toFixed(4)], view: [vx0, efBox[0], vx1, efBox[1]], logy: false,
          title: "", subtitle: showConc ? "" : P.label,
          x: { label: xLabel, major: xmaj, minor: 1 },
          y: { label: c === 0 ? "E\\sF\\N (eV)" : "", major: efmaj, minor: 1 },
          legend: false, legendAt: [vx0 + 0.02, efBox[1] - 0.02], sets });
        data.push({ g, sets, x, log: false });
        g++;
      }
    });

    // data blocks
    data.forEach((D) => D.sets.forEach((s, i) => {
      lines.push(`@target G${D.g}.S${i}`, "@type xy");
      D.x.forEach((xv, k) => {
        const yv = s.y[k];
        if (yv === null || yv === undefined || !Number.isFinite(yv) || (D.log && yv <= 0)) return;
        lines.push(`${Number(xv).toPrecision(8)} ${Number(yv).toExponential(6)}`);
      });
      lines.push("&");
    }));
    return lines.join("\n") + "\n";
  }

  window.DefconAgr = { build, speciesLabel };
})();
