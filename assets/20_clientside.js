/* Clientside callbacks for mmcifXplorer.
 * Per-atom styles for the dash-bio Molecule3dViewer are computed in the browser
 * from the per-atom fields shipped in modelData (b, q, c, e, k, s ...), so a
 * colour/representation/highlight change never re-sends the model. */
(function () {
  const SERIES = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"];
  const SEQ = ["#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5", "#5598e7", "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"];
  const MUTED = "#898781", DIM = "#3a3a37";
  const KIND = {protein: SERIES[0], nucleic: SERIES[1], ligand: SERIES[2], ion: SERIES[3], branched: SERIES[4],
                polymer: SERIES[6], water: "#6b6a66"};
  const SS = {H: SERIES[0], E: SERIES[1], T: SERIES[2]};
  const ELEM = {C: "#909090", N: "#3050f8", O: "#ff0d0d", S: "#ffff30", P: "#ff8000", H: "#ffffff", Se: "#ffa100",
                Zn: "#7d80b0", Fe: "#e06633", Mg: "#8aff00", Ca: "#3dff00", Na: "#ab5cf2", K: "#8f40d4",
                Cl: "#1ff01f", Mn: "#9c7ac7", Cu: "#c88033", Co: "#f090a0", Ni: "#50d050", Br: "#a62929",
                I: "#940094", F: "#90e050", Cd: "#ffd98f", Hg: "#b8b8d0"};
  const RES = {  // amino-acid property classes (5 categorical slots)
    ALA: 0, VAL: 0, LEU: 0, ILE: 0, MET: 0, PHE: 0, TRP: 0, PRO: 4, GLY: 4, CYS: 4, MSE: 0,
    SER: 2, THR: 2, ASN: 2, GLN: 2, TYR: 2, HIS: 1, LYS: 1, ARG: 1, ASP: 3, GLU: 3};
  const RES_COLORS = [SERIES[0], SERIES[1], SERIES[2], SERIES[4], SERIES[3]]; // hydrophobic, +, polar, -, special
  function plddt(v) { return v > 90 ? "#0053d6" : v > 70 ? "#65cbf3" : v > 50 ? "#ffdb13" : "#ff7d45"; }

  function seq(v, lo, hi) {
    let t = (v - lo) / (hi - lo); t = Math.max(0, Math.min(1, t));
    return SEQ[Math.round(t * (SEQ.length - 1))];
  }

  function baseColor(a, mode, meta) {
    switch (mode) {
      case "chain": return SERIES[a.c % SERIES.length];
      case "entity": return SERIES[a.e % SERIES.length];
      case "ss": return a.k === "protein" || a.k === "nucleic" || a.k === "polymer" ? (SS[a.s] || MUTED) : KIND[a.k] || MUTED;
      case "bfactor": return meta && meta.plddt ? plddt(a.b) : seq(a.b, meta.b_range[0], meta.b_range[1]);
      case "occupancy": return a.q <= 0 ? "#d03b3b" : (a.q < 0.999 ? "#fab219" : "#5c6f8a");
      case "element": return ELEM[a.elem] || "#ff1493";
      case "residue": {
        const r = RES[a.residue_name];
        return r === undefined ? (KIND[a.k] || MUTED) : RES_COLORS[r];
      }
      case "kind": default: return KIND[a.k] || MUTED;
    }
  }

  function repFor(a, rep) {
    const poly = a.k === "protein" || a.k === "nucleic" || a.k === "polymer";
    if (rep === "stick") {
      if (a.k === "ion") return {visualization_type: "sphere", radius: 0.7};
      if (a.k === "water") return {visualization_type: "sphere", radius: 0.25};
      return {visualization_type: "stick", radius: 0.18};
    }
    if (rep === "sphere") {
      return {visualization_type: "sphere", scale: a.k === "water" ? 0.4 : 1.0};
    }
    // cartoon (+ ligands as sticks, ions/waters as spheres)
    if (poly) return {visualization_type: "cartoon"};
    if (a.k === "ion") return {visualization_type: "sphere", radius: 0.9};
    if (a.k === "water") return {visualization_type: "sphere", radius: 0.28};
    return {visualization_type: "stick", radius: 0.22};
  }

  window.dash_clientside = Object.assign({}, window.dash_clientside, {
    mmcif: {
      styles: function (modelData, rep, mode, hl, meta) {
        if (!modelData || !modelData.atoms) return window.dash_clientside.no_update;
        meta = meta || {b_range: [0, 100]};
        const atoms = modelData.atoms;
        const n = atoms.length;
        let hlSet = null, focusSet = null, stickSet = null;
        if (hl && hl.serials && hl.serials.length) {
          hlSet = new Uint8Array(n); hl.serials.forEach(i => { if (i < n) hlSet[i] = 1; });
          (hl.context || []).forEach(i => { if (i < n) hlSet[i] = 1; });
        }
        if (hl && hl.focus && hl.focus.length) {
          focusSet = new Uint8Array(n); hl.focus.forEach(i => { if (i < n) focusSet[i] = 1; });
        }
        if (hl && hl.sticks && hl.sticks.length) {
          stickSet = new Uint8Array(n); hl.sticks.forEach(i => { if (i < n) stickSet[i] = 1; });
        }
        const out = new Array(n);
        for (let i = 0; i < n; i++) {
          const a = atoms[i];
          let st = repFor(a, rep);
          let color = baseColor(a, mode, meta);
          if (hlSet && !hlSet[i]) color = DIM;
          if (focusSet && focusSet[i]) {
            st = {visualization_type: "sphere", radius: 0.65};
            color = ELEM[a.elem] || "#ffffff";
          } else if (stickSet && stickSet[i] && rep !== "sphere") {
            st = {visualization_type: "stick", radius: 0.2};
            color = ELEM[a.elem] || color;
          }
          out[i] = Object.assign({color: color}, st);
        }
        return out;
      },
      viewerSize: function (_n, _entry) {
        const el = document.getElementById("viewer-wrap");
        const w = el ? Math.max(280, Math.floor(el.getBoundingClientRect().width)) : 600;
        const h = Math.max(360, Math.min(640, Math.floor(window.innerHeight * 0.58)));
        return [w, h];
      }
    }
  });

  // Keep the 3D viewer sized to its column.
  let t = null;
  window.addEventListener("resize", function () {
    clearTimeout(t);
    t = setTimeout(function () {
      if (window.dash_clientside && window.dash_clientside.set_props) {
        window.dash_clientside.set_props("viewer-resize", {data: Date.now()});
      }
    }, 250);
  });
})();
