/*
 * Quantum Wave Lab: browser app. Physics in js/qwave.js (QWave), scenes in js/scenes.js, share links in
 * js/share.js, strings in js/i18n.js. The field shows the central 40 × 40 of a 51.2 × 51.2 periodic box;
 * the rim outside the view is the absorbing layer, so packets start clear of it (docs/porting.md).
 */
(function () {
  "use strict";
  const { Grid, Solver, FluxDetector, gaussianPacket, groundStateSteps, measurePosition, momentumDensity } = QWave;
  const { SCENES, ORDER, becTrap, coverage } = QWaveScenes;
  const Share = QWaveShare;

  // ---------------------------------------------------------------- constants
  const L = 51.2, VIEW = 40;
  const ABSORBER = { width: (L - VIEW) / 2, gammaMax: 6 };   // γ_max ≳ 5k/width for k up to ~6
  const DT = 0.01;                                             // V·Δt ≤ 1 for walls up to V = 100
  const SIGMA_M = 0.8;                                         // resolution of the position measurement
  const K_VIEW = 12;                                           // momentum view: |k_x|, |k_y| ≤ 12 (or k_max)
  const ISO_LEVELS = [1, 2, 4, 8, 16, 32, 64];                 // contour lines of the smooth background V
  const INFERNO = hexToBytes(
    "00000401000501010601010802010a02020c02020e03021004031204031405041706041907051b08051d09061f0a0722" +
    "0b07240c08260d08290e092b10092d110a30120a32140b34150b37160b39180c3c190c3e1b0c411c0c431e0c451f0c48" +
    "210c4a230c4c240c4f260c51280b53290b552b0b572d0b592f0a5b310a5c320a5e340a5f3609613809623909633b0964" +
    "3d09653e0966400a67420a68440a68450a69470b6a490b6a4a0c6b4c0c6b4d0d6c4f0d6c510e6c520e6d540f6d550f6d" +
    "57106e59106e5a116e5c126e5d126e5f136e61136e62146e64156e65156e67166e69166e6a176e6c186e6d186e6f196e" +
    "71196e721a6e741a6e751b6e771c6d781c6d7a1d6d7c1d6d7d1e6d7f1e6c801f6c82206c84206b85216b87216b88226a" +
    "8a226a8c23698d23698f24699025689225689326679526679727669827669a28659b29649d29649f2a63a02a63a22b62" +
    "a32c61a52c60a62d60a82e5fa92e5eab2f5ead305dae305cb0315bb1325ab3325ab43359b63458b73557b93556ba3655" +
    "bc3754bd3853bf3952c03a51c13a50c33b4fc43c4ec63d4dc73e4cc83f4bca404acb4149cc4248ce4347cf4446d04545" +
    "d24644d34743d44842d54a41d74b3fd84c3ed94d3dda4e3cdb503bdd513ade5238df5337e05536e15635e25734e35933" +
    "e45a31e55c30e65d2fe75e2ee8602de9612bea632aeb6429eb6628ec6726ed6925ee6a24ef6c23ef6e21f06f20f1711f" +
    "f1731df2741cf3761bf37819f47918f57b17f57d15f67e14f68013f78212f78410f8850ff8870ef8890cf98b0bf98c0a" +
    "f98e09fa9008fa9207fa9407fb9606fb9706fb9906fb9b06fb9d07fc9f07fca108fca309fca50afca60cfca80dfcaa0f" +
    "fcac11fcae12fcb014fcb216fcb418fbb61afbb81dfbba1ffbbc21fbbe23fac026fac228fac42afac62df9c72ff9c932" +
    "f9cb35f8cd37f8cf3af7d13df7d340f6d543f6d746f5d949f5db4cf4dd4ff4df53f4e156f3e35af3e55df2e661f2e865" +
    "f2ea69f1ec6df1ed71f1ef75f1f179f2f27df2f482f3f586f3f68af4f88ef5f992f6fa96f8fb9af9fc9dfafda1fcffa4");
  const WALL_RGB = [77, 208, 225], DETECTOR_RGB = "245, 182, 66";

  function hexToBytes(hex) {
    const out = new Uint8Array(hex.length / 2);
    for (let i = 0; i < out.length; i++) out[i] = parseInt(hex.substr(2 * i, 2), 16);
    return out;
  }

  // ---------------------------------------------------------------- state
  const $ = (id) => document.getElementById(id);
  const field = $("field"), ctx = field.getContext("2d");
  const chart = $("chart"), cctx = chart.getContext("2d");
  const state = {
    lang: pickLanguage(), N: 256, scene: "doubleSlit", tool: "wall", view: "position", display: "amplitude",
    auto: false, contours: true, running: true, speed: 6, bright: 1, brush: 1, height: 100,
    k0: 4, sigma: 2.5, g: 0, detectorOn: true, detectorX: 0, packet: null, customPacket: false,
    grid: null, solver: null, sceneWalls: null, baseV: null, wallV: null, vDirty: false, detector: null,
    strokes: [], redo: [], stroke: null,
    scale: 1, energy: null, samples: [], frame: 0, msPerStep: null, fps: null, lastFrame: null,
    drag: null, hover: null, issues: [], becCache: {}, becPsi: null, busy: false,
    measurement: null, contours_: null, recorder: null, kIndex: null, hintTimer: null,
  };

  function pickLanguage() {
    try { const s = localStorage.getItem("qwl-lang"); if (s === "en" || s === "ja") return s; } catch (e) { /* no storage */ }
    return (navigator.language || "").startsWith("ja") ? "ja" : "en";
  }
  function t(key, vars = {}) {
    const s = I18N[state.lang][key] ?? I18N.en[key] ?? key;
    return s.replace(/\{(\w+)\}/g, (_, k) => (k in vars ? vars[k] : `{${k}}`));
  }
  const percent = (v) => {
    const p = 100 * v;
    return `${p >= 10 ? p.toFixed(0) : p >= 1 ? p.toFixed(1) : p.toPrecision(2)}%`;
  };

  // ---------------------------------------------------------------- grid and walls
  const viewCells = () => Math.round(VIEW / state.grid.dx);
  const viewOffset = () => (state.N - viewCells()) / 2;

  function setGrid(N) {
    state.N = N;
    state.grid = new Grid([N, N], [L, L]);
    // momentum view: indices with |k| ≤ K_VIEW, in ascending k
    const k = state.grid.kx, kv = Math.min(K_VIEW, state.grid.kmax[0]);
    state.kIndex = [...k.keys()].filter((i) => Math.abs(k[i]) <= kv + 1e-9).sort((a, b) => k[a] - k[b]);
    off = null;
  }

  function totalV() {
    const V = new Float64Array(state.grid.size);
    for (let i = 0; i < V.length; i++) V[i] = state.baseV[i] + state.wallV[i];
    return V;
  }

  /** Walls = the scene's walls with every stroke replayed on top (undo = replay without the last one). */
  function rebuildWalls() {
    state.wallV = Float64Array.from(state.sceneWalls);
    for (const st of state.strokes) applyStroke(st);
    state.vDirty = true;
    state.contours_ = null;
    updateUndoButtons();
  }

  function applyStroke(st) {
    if (st.clear) { state.wallV.fill(0); return; }
    const h = st.erase ? 0 : st.h, pts = st.pts;
    paintDisk(pts[0][0], pts[0][1], st.r, h);
    for (let i = 1; i < pts.length; i++) paintSegment(pts[i - 1], pts[i], st.r, h);
  }

  function paintSegment(a, b, R, h) {
    const steps = Math.max(1, Math.ceil(Math.hypot(b[0] - a[0], b[1] - a[1]) / (R / 3)));
    for (let s = 1; s <= steps; s++) paintDisk(a[0] + ((b[0] - a[0]) * s) / steps, a[1] + ((b[1] - a[1]) * s) / steps, R, h);
  }

  function paintDisk(cx, cy, R, h) {
    const g = state.grid, W = state.wallV;
    const i0 = Math.max(0, Math.floor((cx - R + L / 2) / g.dx) - 1), i1 = Math.min(g.nx - 1, Math.ceil((cx + R + L / 2) / g.dx) + 1);
    const j0 = Math.max(0, Math.floor((cy - R + L / 2) / g.dy) - 1), j1 = Math.min(g.ny - 1, Math.ceil((cy + R + L / 2) / g.dy) + 1);
    for (let i = i0; i <= i1; i++) {
      for (let j = j0; j <= j1; j++) {
        const a = coverage(Math.hypot(g.xs[i] - cx, g.ys[j] - cy), R, g.dx);   // soft 1-cell edge
        if (a > 0) { const idx = i * g.ny + j; W[idx] += (h - W[idx]) * a; }
      }
    }
  }

  function pushStroke(st) {
    // on the 0.05 lattice of share links, so a shared drawing is rasterized exactly like the original
    if (st.pts) st.pts = st.pts.map(([x, y]) => [Math.round(x * 20) / 20, Math.round(y * 20) / 20]);
    state.strokes.push(st);
    state.redo = [];
    rebuildWalls();
  }
  function undo() { if (state.strokes.length) { state.redo.push(state.strokes.pop()); rebuildWalls(); } }
  function redo() { if (state.redo.length) { state.strokes.push(state.redo.pop()); rebuildWalls(); } }
  function updateUndoButtons() { $("undo").disabled = !state.strokes.length; $("redo").disabled = !state.redo.length; }

  // ---------------------------------------------------------------- scenes and launching
  function dtFor(grid, g) {
    // Gross–Pitaevskii: Δt·k_max²/2 < π (Weideman & Herbst), with a 10 % margin
    return g > 0 ? Math.min(DT, (0.9 * 2 * Math.PI) / grid.kmax.reduce((s, k) => s + k * k, 0)) : DT;
  }

  function buildSolver(keepPsi) {
    const old = state.solver;
    let re = null, im = null;
    if (keepPsi && old) { old.norm(); re = old.re; im = old.im; }        // norm() also undoes a pending half-step
    const s = new Solver(state.grid, totalV(), dtFor(state.grid, state.g), { absorber: ABSORBER, nonlinearity: state.g });
    if (re) { s.setPsi(re, im); s.t = old.t; s.nsteps = old.nsteps; }
    state.solver = s;
    state.vDirty = false;
  }

  /** opts (from a share link): strokes, packet, k0, sigma, g, detectorOn, detectorX */
  function loadScene(name, opts = {}) {
    const sc = SCENES[name], g = state.grid;
    state.scene = name;
    state.sceneWalls = sc.walls ? sc.walls(g) : new Float64Array(g.size);
    state.baseV = sc.base ? g.fill(sc.base) : new Float64Array(g.size);
    state.strokes = opts.strokes ?? [];
    state.redo = [];
    rebuildWalls();
    state.detectorOn = opts.detectorOn ?? (sc.detector !== null && sc.detector !== undefined);
    state.detectorX = opts.detectorX ?? sc.detector ?? 0;
    $("detector-on").checked = state.detectorOn;
    document.querySelectorAll("#scenes button").forEach((b) => b.classList.toggle("active", b.dataset.scene === name));
    updateSceneInfo();
    if (sc.bec) { prepareBec(); return; }
    state.customPacket = Boolean(opts.packet);
    state.packet = opts.packet ? { ...opts.packet, aspect: 1 } : { ...sc.packet };
    state.k0 = opts.k0 ?? sc.packet.k0;
    state.sigma = opts.sigma ?? sc.packet.sigma;
    state.g = opts.g ?? 0;
    syncSliders();
    buildSolver(false);
    relaunch();
  }

  function updateSceneInfo() {
    const facts = SCENES[state.scene].facts ?? {};
    const vars = Object.fromEntries(Object.entries(facts).map(([k, v]) => [k, percent(v)]));
    $("scene-info").textContent = t(`info.${state.scene}`, vars);
  }

  function relaunch() {
    if (state.scene === "bec" && state.becPsi && !state.customPacket) { startFrom(state.becPsi[0], state.becPsi[1]); return; }
    const p = state.packet, n = Math.hypot(p.dir[0], p.dir[1]) || 1;
    const k = [(state.k0 * p.dir[0]) / n, (state.k0 * p.dir[1]) / n];
    const [re, im] = gaussianPacket(state.grid, p.center, [state.sigma, state.sigma * (p.aspect ?? 1)], k, Float32Array);
    startFrom(re, im);
  }

  function startFrom(re, im) {
    const s = state.solver;
    if (state.vDirty) { s.setPotential(totalV()); state.vDirty = false; }
    s.setPsi(re, im);
    s.t = 0; s.nsteps = 0;
    state.scale = maxAmplitude(s);
    state.measurement = null;
    state.energy = s.energy();
    state.contours_ = null;
    resetDetector();
    state.samples = [];
    runChecks();
    render();
  }

  function maxAmplitude(s) {
    let m = 0;
    for (let i = 0; i < s.re.length; i++) m = Math.max(m, s.re[i] * s.re[i] + s.im[i] * s.im[i]);
    return Math.sqrt(m) || 1;
  }

  function resetDetector() {
    state.detector = state.detectorOn ? new FluxDetector(state.grid, state.detectorX, { method: "spectral" }) : null;
    if (state.detector) state.detector.record(state.solver);
  }

  // ---------------------------------------------------------------- Bose–Einstein condensate scene
  function prepareBec() {
    const N = state.N;
    state.g = 300; state.customPacket = false; state.packet = { ...SCENES.free.packet };
    syncSliders();
    const release = () => {
      state.baseV = new Float64Array(state.grid.size);        // trap switched off: free expansion
      state.becPsi = state.becCache[N];
      buildSolver(false);
      relaunch();
    };
    if (state.becCache[N]) { release(); return; }
    const g = state.grid;
    // dτ = 0.04 converges in ~330 steps; the state differs from a dτ = 0.005 run by ~6e-5 (L2)
    const it = groundStateSteps(g, g.fill(becTrap), { dtau: 0.04, nonlinearity: 300, tol: 1e-6, chunk: 10 });
    setBusy(true);
    let last = null;
    const work = () => {
      if (state.N !== N || state.scene !== "bec") { setBusy(false); return; }   // the user moved on
      const t0 = performance.now();
      let r;
      do { r = it.next(); if (!r.done) last = r.value; } while (!r.done && !r.value.converged && performance.now() - t0 < 40);
      if (r.done || r.value.converged) {
        state.becCache[N] = [Float32Array.from(last.re), Float32Array.from(last.im)];
        setBusy(false);
        release();
        return;
      }
      $("busy-text").textContent = `${t("preparing")}  μ = ${last.mu.toFixed(4)}`;
      setTimeout(work, 0);
    };
    setTimeout(work, 0);
  }

  function setBusy(on) {
    state.busy = on;
    $("busy").hidden = !on;
    if (on) $("busy-text").textContent = t("preparing");
  }

  // ---------------------------------------------------------------- measurement
  function measure() {
    if (state.busy || !state.solver) return;
    const [x, y] = measurePosition(state.solver, SIGMA_M);
    state.measurement = { x, y, time: performance.now() };
    if (!state.auto) state.scale = maxAmplitude(state.solver);
    state.energy = state.solver.energy();
    state.contours_ = null;
    flashHint(t("measured", { x: x.toFixed(1), y: y.toFixed(1) }));
    render();
  }

  // ---------------------------------------------------------------- pointer tools
  function toPhys(ev) {
    const r = field.getBoundingClientRect();
    return [-VIEW / 2 + ((ev.clientX - r.left) / r.width) * VIEW, VIEW / 2 - ((ev.clientY - r.top) / r.height) * VIEW];
  }

  field.addEventListener("pointerdown", (ev) => {
    if (state.busy || state.view !== "position") return;
    field.setPointerCapture(ev.pointerId);
    const p = toPhys(ev);
    const erase = state.tool === "eraser" || ev.button === 2;
    if (state.tool === "wall" || state.tool === "eraser" || ev.button === 2) {
      if (ev.shiftKey) state.drag = { kind: "line", from: p, to: p, erase };
      else {
        state.stroke = { erase, r: state.brush, h: state.height, pts: [p] };
        state.drag = { kind: "paint" };
        paintDisk(p[0], p[1], state.brush, erase ? 0 : state.height);
        state.vDirty = true;
      }
    } else if (state.tool === "launch") state.drag = { kind: "launch", from: p, to: p };
    else if (state.tool === "detector") { state.detectorX = p[0]; state.detectorOn = true; $("detector-on").checked = true; resetDetector(); }
  });
  field.addEventListener("pointermove", (ev) => {
    const p = toPhys(ev);
    state.hover = p;
    const d = state.drag;
    if (!d) return;
    if (d.kind === "paint") {
      const st = state.stroke, last = st.pts[st.pts.length - 1];
      paintSegment(last, p, st.r, st.erase ? 0 : st.h);
      st.pts.push(p);
      state.vDirty = true;
    } else d.to = p;
  });
  const endDrag = (ev) => {
    const d = state.drag;
    state.drag = null;
    if (!d) return;
    if (d.kind === "paint") {
      const st = state.stroke;
      state.stroke = null;
      st.pts = Share.thin(st.pts, st.r / 4);
      pushStroke(st);                                  // replayed from the (thinned) stroke list
    } else if (d.kind === "line") {
      pushStroke({ erase: d.erase, r: state.brush, h: state.height, pts: [d.from, toPhys(ev)] });
    } else if (d.kind === "launch") {
      const p = toPhys(ev), dir = [p[0] - d.from[0], p[1] - d.from[1]];
      state.customPacket = true;
      state.packet = { center: d.from, dir: Math.hypot(dir[0], dir[1]) > 0.5 ? dir : [1, 0], aspect: 1 };
      relaunch();
      if (!state.running) setRunning(true);
    }
  };
  field.addEventListener("pointerup", endDrag);
  field.addEventListener("pointercancel", () => { if (state.stroke) { state.stroke = null; rebuildWalls(); } state.drag = null; });
  field.addEventListener("pointerleave", () => { state.hover = null; });
  field.addEventListener("contextmenu", (ev) => ev.preventDefault());

  // ---------------------------------------------------------------- rendering
  let img = null, off = null, offCtx = null;

  function ensureBuffer(n) {
    if (!off || off.width !== n) {
      off = document.createElement("canvas"); off.width = off.height = n;
      offCtx = off.getContext("2d");
      img = offCtx.createImageData(n, n);
    }
  }

  function render() {
    if (!state.solver) return;
    if (state.solver.shifted) state.solver.norm();                    // undo a pending half-step
    const S = field.width;
    if (state.view === "momentum") renderMomentum(S); else renderPosition(S);
  }

  function colour(v, a2, re, im, mode, inv, inv2) {
    if (mode === "phase") return hsv(Math.atan2(im, re), Math.min(1, Math.sqrt(a2) * inv));
    const x = mode === "density" ? a2 * inv2 : Math.sqrt(a2) * inv;
    const k = 3 * Math.min(255, Math.max(0, (x * 255) | 0));
    return [INFERNO[k], INFERNO[k + 1], INFERNO[k + 2]];
  }

  function renderPosition(S) {
    const s = state.solver, g = state.grid, n = viewCells(), o = viewOffset(), ny = g.ny;
    ensureBuffer(n);
    const re = s.re, im = s.im, W = state.wallV, data = img.data;
    const inv = state.bright / state.scale, inv2 = inv * inv, mode = state.display;
    let amax = 0;
    for (let r = 0; r < n; r++) {
      const iy = o + n - 1 - r;
      for (let c = 0; c < n; c++) {
        const idx = (o + c) * ny + iy, p = 4 * (r * n + c);
        const a2 = re[idx] * re[idx] + im[idx] * im[idx];
        if (a2 > amax) amax = a2;
        let [R, G, B] = colour(0, a2, re[idx], im[idx], mode, inv, inv2);
        const w = W[idx];
        if (w > 0.5) {                                                   // wall tint, as in the README GIF
          const f = 0.35 + 0.35 * Math.min(w / 100, 1);
          R = R * (1 - f) + WALL_RGB[0] * f; G = G * (1 - f) + WALL_RGB[1] * f; B = B * (1 - f) + WALL_RGB[2] * f;
        }
        data[p] = R; data[p + 1] = G; data[p + 2] = B; data[p + 3] = 255;
      }
    }
    if (state.auto) state.scale = 0.8 * state.scale + 0.2 * (Math.sqrt(amax) || state.scale);
    offCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = "high";
    ctx.drawImage(off, 0, 0, S, S);
    if (state.contours) drawContours(S);
    drawPositionOverlay(S);
  }

  function renderMomentum(S) {
    const s = state.solver, g = state.grid, idx = state.kIndex, m = idx.length, ny = g.ny;
    ensureBuffer(m);
    const phi2 = momentumDensity(s), data = img.data;
    let pmax = 0;
    for (const i of idx) for (const j of idx) pmax = Math.max(pmax, phi2[i * ny + j]);
    const inv2 = state.bright / (pmax || 1), inv = Math.sqrt(inv2), mode = state.display === "density" ? "density" : "amplitude";
    for (let r = 0; r < m; r++) {
      const j = idx[m - 1 - r];
      for (let c = 0; c < m; c++) {
        const a2 = phi2[idx[c] * ny + j], p = 4 * (r * m + c);
        const [R, G, B] = colour(0, a2, 0, 0, mode, inv, inv2);
        data[p] = R; data[p + 1] = G; data[p + 2] = B; data[p + 3] = 255;
      }
    }
    offCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, 0, 0, S, S);
    drawMomentumOverlay(S, g.kx[idx[0]], g.kx[idx[m - 1]], m);
  }

  function hsv(angle, v) {                                             // hue = phase, value = amplitude
    const h = (((angle / (2 * Math.PI)) + 1) % 1) * 6, i = Math.floor(h), f = h - i;
    const q = v * (1 - f), u = v * f;
    const rgb = [[v, u, 0], [q, v, 0], [0, v, u], [0, q, v], [u, 0, v], [v, 0, q]][i % 6];
    return [rgb[0] * 255, rgb[1] * 255, rgb[2] * 255];
  }

  const toPx = (x, S) => ((x + VIEW / 2) / VIEW) * S;
  const toPy = (y, S) => ((VIEW / 2 - y) / VIEW) * S;

  // ---------------------------------------------------------------- contour lines (marching squares)
  /** Path2D of the level set f = level on the view grid (f[r * n + c], row 0 at the top). */
  function isoPath(f, n, level, S, path = new Path2D()) {
    const cell = S / n, px = (c) => (c + 0.5) * cell;
    const cross = (a, b) => (level - a) / (b - a);
    for (let r = 0; r < n - 1; r++) {
      for (let c = 0; c < n - 1; c++) {
        const tl = f[r * n + c], tr = f[r * n + c + 1], bl = f[(r + 1) * n + c], br = f[(r + 1) * n + c + 1];
        const pts = [];
        if ((tl > level) !== (tr > level)) pts.push([px(c + cross(tl, tr)), px(r)]);
        if ((tr > level) !== (br > level)) pts.push([px(c + 1), px(r + cross(tr, br))]);
        if ((bl > level) !== (br > level)) pts.push([px(c + cross(bl, br)), px(r + 1)]);
        if ((tl > level) !== (bl > level)) pts.push([px(c), px(r + cross(tl, bl))]);
        for (let k = 0; k + 1 < pts.length; k += 2) { path.moveTo(...pts[k]); path.lineTo(...pts[k + 1]); }
      }
    }
    return path;
  }

  function viewField(V) {
    const g = state.grid, n = viewCells(), o = viewOffset(), f = new Float64Array(n * n);
    for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) f[r * n + c] = V[(o + c) * g.ny + (o + n - 1 - r)];
    return f;
  }

  function drawContours(S) {
    const E = state.energy;
    let cc = state.contours_;
    if (!cc || cc.S !== S || (E !== null && Math.abs(E - cc.E) > 0.01 * Math.abs(E) + 0.01)) {
      const n = viewCells(), base = viewField(state.baseV), total = viewField(totalV());
      const iso = new Path2D();
      const bmax = base.reduce((a, b) => Math.max(a, b), 0);
      for (const lv of ISO_LEVELS) if (lv < bmax) isoPath(base, n, lv, S, iso);
      const tmax = total.reduce((a, b) => Math.max(a, b), 0);
      const wmax = viewField(state.wallV).reduce((a, b) => Math.max(a, b), 0);
      // V = ⟨E⟩: edge of the classically forbidden region, for smooth traps and barriers comparable to E.
      // Around hard walls (V ≫ E) it would only trace the walls. Not meaningful with g ≠ 0.
      const eLine = E !== null && E > 0 && state.g === 0 && tmax > E && wmax < 4 * E ? isoPath(total, n, E, S) : null;
      cc = state.contours_ = { S, E, iso, eLine };
    }
    ctx.save();
    ctx.strokeStyle = "rgba(210, 220, 230, 0.28)"; ctx.lineWidth = 1.2;
    ctx.stroke(cc.iso);
    if (cc.eLine) { ctx.strokeStyle = "rgba(255, 255, 255, 0.8)"; ctx.lineWidth = 1.6; ctx.setLineDash([2, 5]); ctx.stroke(cc.eLine); }
    ctx.restore();
  }

  // ---------------------------------------------------------------- overlays
  function drawPositionOverlay(S) {
    if (state.detector) {
      const x = toPx(state.detector.position, S);
      ctx.save();
      ctx.strokeStyle = `rgba(${DETECTOR_RGB}, 0.85)`; ctx.lineWidth = 2; ctx.setLineDash([10, 8]);
      ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, S); ctx.stroke();
      ctx.restore();
    }
    const m = state.measurement;
    if (m) {
      const age = (performance.now() - m.time) / 2500;
      if (age < 1) {
        const x = toPx(m.x, S), y = toPy(m.y, S), R = ((2 * SIGMA_M) / VIEW) * S;
        ctx.save();
        ctx.globalAlpha = 1 - age; ctx.strokeStyle = "#fff"; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(x, y, R * (1 + age), 0, 2 * Math.PI); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(x - R * 2.2, y); ctx.lineTo(x - R * 1.3, y); ctx.moveTo(x + R * 1.3, y); ctx.lineTo(x + R * 2.2, y);
        ctx.moveTo(x, y - R * 2.2); ctx.lineTo(x, y - R * 1.3); ctx.moveTo(x, y + R * 1.3); ctx.lineTo(x, y + R * 2.2); ctx.stroke();
        ctx.restore();
      }
    }
    const d = state.drag;
    if (d && d.kind === "launch") {
      const [x0, y0] = [toPx(d.from[0], S), toPy(d.from[1], S)], [x1, y1] = [toPx(d.to[0], S), toPy(d.to[1], S)];
      ctx.save();
      ctx.strokeStyle = "#fff"; ctx.fillStyle = "#fff"; ctx.lineWidth = 3;
      ctx.beginPath(); ctx.arc(x0, y0, (state.sigma / VIEW) * S, 0, 2 * Math.PI); ctx.globalAlpha = 0.25; ctx.fill();
      ctx.globalAlpha = 1; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      const a = Math.atan2(y1 - y0, x1 - x0);
      ctx.beginPath(); ctx.moveTo(x1, y1);
      ctx.lineTo(x1 - 16 * Math.cos(a - 0.4), y1 - 16 * Math.sin(a - 0.4));
      ctx.lineTo(x1 - 16 * Math.cos(a + 0.4), y1 - 16 * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
      ctx.restore();
    }
    if (d && d.kind === "line") {
      ctx.save();
      ctx.strokeStyle = d.erase ? "rgba(255,255,255,0.8)" : `rgb(${WALL_RGB.join(",")})`;
      ctx.lineWidth = (2 * state.brush / VIEW) * S; ctx.lineCap = "round"; ctx.globalAlpha = 0.5;
      ctx.beginPath(); ctx.moveTo(toPx(d.from[0], S), toPy(d.from[1], S)); ctx.lineTo(toPx(d.to[0], S), toPy(d.to[1], S)); ctx.stroke();
      ctx.restore();
    }
    if (state.hover && !d && (state.tool === "wall" || state.tool === "eraser")) {
      ctx.save();
      ctx.strokeStyle = state.tool === "wall" ? `rgb(${WALL_RGB.join(",")})` : "rgba(255,255,255,0.8)";
      ctx.lineWidth = 1.5; ctx.setLineDash(state.tool === "eraser" ? [4, 4] : []);
      ctx.beginPath(); ctx.arc(toPx(state.hover[0], S), toPy(state.hover[1], S), (state.brush / VIEW) * S, 0, 2 * Math.PI); ctx.stroke();
      ctx.restore();
    }
    drawScaleBar(S);
  }

  function drawScaleBar(S) {
    const len = (5 / VIEW) * S, x = 24, y = S - 26;
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.85)"; ctx.fillStyle = "rgba(255,255,255,0.85)"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + len, y); ctx.moveTo(x, y - 6); ctx.lineTo(x, y + 6); ctx.moveTo(x + len, y - 6); ctx.lineTo(x + len, y + 6); ctx.stroke();
    ctx.font = "20px system-ui, sans-serif"; ctx.textBaseline = "bottom";
    ctx.fillText("5", x + len / 2 - 5, y - 6);
    if (state.scene !== "bec" || state.customPacket) ctx.fillText(`λ = 2π/k₀ = ${((2 * Math.PI) / state.k0).toFixed(2)}`, x + len + 16, y + 8);
    ctx.restore();
  }

  function drawMomentumOverlay(S, kFirst, kLast, m) {
    const dk = (kLast - kFirst) / (m - 1), span = m * dk, k0px = (k) => ((k - kFirst + dk / 2) / span) * S;
    const X = k0px, Y = (k) => S - k0px(k);
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.25)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(X(0), 0); ctx.lineTo(X(0), S); ctx.moveTo(0, Y(0)); ctx.lineTo(S, Y(0)); ctx.stroke();
    ctx.fillStyle = "rgba(255,255,255,0.7)"; ctx.font = "18px system-ui, sans-serif";
    for (let k = Math.ceil(kFirst / 5) * 5; k <= kLast; k += 5) {
      if (k === 0) continue;
      ctx.textAlign = "center"; ctx.fillText(String(k), X(k), Y(0) + 24);
      ctx.textAlign = "left"; ctx.fillText(String(k), X(0) + 8, Y(k) + 6);
    }
    const label = (sub, x, y) => {                     // k with a subscript (no reliance on Unicode subscripts)
      ctx.font = "italic 22px system-ui, sans-serif"; ctx.fillText("k", x, y);
      ctx.font = "italic 15px system-ui, sans-serif"; ctx.fillText(sub, x + 12, y + 6);
    };
    ctx.textAlign = "left";
    label("x", S - 30, Y(0) - 12);
    label("y", X(0) - 34, 26);
    ctx.strokeStyle = `rgba(${DETECTOR_RGB}, 0.8)`; ctx.setLineDash([6, 6]); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(X(0), Y(0), (state.k0 / span) * S, 0, 2 * Math.PI); ctx.stroke();
    ctx.restore();
  }

  function drawChart() {
    const dpr = window.devicePixelRatio || 1, W = Math.round(chart.clientWidth * dpr), H = Math.round(chart.clientHeight * dpr);
    if (chart.width !== W || chart.height !== H) { chart.width = W; chart.height = H; }
    cctx.clearRect(0, 0, W, H);
    const padL = 38 * dpr, padR = 8 * dpr, padT = 6 * dpr, padB = 16 * dpr;
    const smp = state.samples, tmax = smp.length ? Math.max(smp[smp.length - 1].t, 1) : 1;
    const X = (tt) => padL + ((W - padL - padR) * tt) / tmax, Y = (v) => H - padB - (H - padT - padB) * Math.min(Math.max(v, 0), 1);
    cctx.font = `${11 * dpr}px system-ui, sans-serif`; cctx.fillStyle = "#8b949e"; cctx.strokeStyle = "#30363d"; cctx.lineWidth = dpr;
    for (const v of [0, 0.5, 1]) {
      cctx.beginPath(); cctx.moveTo(padL, Y(v)); cctx.lineTo(W - padR, Y(v)); cctx.stroke();
      cctx.fillText(`${v * 100}%`, 4 * dpr, Y(v) + 4 * dpr);
    }
    cctx.fillText(`t = ${tmax.toFixed(1)}`, W - padR - 52 * dpr, H - 3 * dpr);
    if (smp.length < 2) return;
    const line = (f, color) => {
      cctx.strokeStyle = color; cctx.lineWidth = 2 * dpr; cctx.beginPath();
      smp.forEach((q, i) => (i ? cctx.lineTo(X(q.t), Y(f(q))) : cctx.moveTo(X(q.t), Y(f(q)))));
      cctx.stroke();
    };
    line((q) => q.norm, "#c9d1d9");
    line((q) => 1 - q.norm, "#f47067");
    if (state.detector) line((q) => q.trans, `rgb(${DETECTOR_RGB})`);
  }

  // ---------------------------------------------------------------- readouts, checks, hints
  function updateReadouts() {
    const s = state.solver;
    $("r-time").textContent = s.t.toFixed(2);
    $("r-energy").textContent = state.energy === null ? "—" : state.energy.toFixed(2);
    $("r-norm").textContent = `${(100 * s.norm()).toFixed(1)}%`;
    $("r-trans").textContent = state.detector ? `${(100 * Math.max(0, state.detector.transmitted)).toFixed(2)}%` : "—";
    const perf = $("r-perf");
    if (state.msPerStep === null) perf.textContent = "—";
    else perf.innerHTML = `${state.msPerStep.toFixed(1)} ms <small>${Math.round(state.fps ?? 0)} fps</small>`;
  }

  function runChecks() {
    // strict (1e-6) on the initial state, 1e-3 while running: hard walls always scatter a trace of
    // probability to high k, which is harmless; one message per kind (the worst axis for k_max)
    const tol = state.solver.t === 0 ? 1e-6 : 1e-3, seen = new Map();
    for (const is of state.solver.diagnose({ tol })) {
      if (!seen.has(is.code) || is.value > seen.get(is.code).value) seen.set(is.code, is);
    }
    state.issues = [...seen.values()];
    const ul = $("checks");
    ul.innerHTML = "";
    const add = (cls, text) => { const li = document.createElement("li"); li.className = cls; li.textContent = text; ul.appendChild(li); };
    if (!state.issues.length) add("ok", t("checksOk"));
    for (const is of state.issues) {
      add("warn", t(`diag.${is.code}`, {
        value: is.code === "kmax" ? percent(is.value) : is.value.toPrecision(3),
        limit: (1 / state.solver.dt).toFixed(0), kmax: is.kmax ? is.kmax.toFixed(1) : "",
      }));
    }
  }

  function defaultHint() { return state.view === "momentum" ? t("hint.momentum") : t(`hint.${state.tool}`); }
  function flashHint(text) {
    const h = $("hint");
    h.textContent = text; h.classList.add("flash");
    clearTimeout(state.hintTimer);
    state.hintTimer = setTimeout(() => { h.textContent = defaultHint(); h.classList.remove("flash"); }, 3000);
  }

  // ---------------------------------------------------------------- main loop
  function frame(now) {
    advance(now);
    requestAnimationFrame(frame);
  }

  function advance(now) {
    if (state.lastFrame !== null) {
      const f = 1000 / Math.max(1, now - state.lastFrame);
      state.fps = state.fps === null ? f : 0.9 * state.fps + 0.1 * f;
    }
    state.lastFrame = now;
    if (!state.solver || state.busy) return;
    const s = state.solver;
    if (state.vDirty) { s.setPotential(totalV()); state.vDirty = false; state.contours_ = null; }
    if (state.running) {
      const t0 = performance.now();
      s.step(state.speed);
      const ms = (performance.now() - t0) / state.speed;
      state.msPerStep = state.msPerStep === null ? ms : 0.9 * state.msPerStep + 0.1 * ms;
      if (state.detector) state.detector.record(s);
      state.frame++;
      if (state.frame % 2 === 0) {
        state.samples.push({ t: s.t, norm: s.norm(), trans: state.detector ? state.detector.transmitted : 0 });
        if (state.samples.length > 1500) state.samples = state.samples.filter((_, i) => i % 2 === 0);
      }
      if (state.frame % 8 === 0) state.energy = s.energy();
      if (state.frame % 60 === 0) runChecks();
    }
    render();
    if (state.frame % 4 === 0 || !state.running) { updateReadouts(); drawChart(); }
  }

  // ---------------------------------------------------------------- export and share
  function download(blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  function saveImage() {
    field.toBlob((b) => b && download(b, `quantum-wave-lab-${state.scene}-t${state.solver.t.toFixed(1)}.png`), "image/png");
  }

  function toggleRecording() {
    if (state.recorder) { state.recorder.stop(); return; }
    const type = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm", "video/mp4"].find((m) => MediaRecorder.isTypeSupported(m));
    const rec = new MediaRecorder(field.captureStream(30), type ? { mimeType: type, videoBitsPerSecond: 8e6 } : undefined);
    const chunks = [];
    rec.ondataavailable = (e) => { if (e.data.size) chunks.push(e.data); };
    rec.onstop = () => {
      const ext = rec.mimeType.includes("mp4") ? "mp4" : "webm";
      download(new Blob(chunks, { type: rec.mimeType }), `quantum-wave-lab-${state.scene}.${ext}`);
      state.recorder = null;
      updateRecordButton();
    };
    rec.start();
    state.recorder = rec;
    updateRecordButton();
  }
  function updateRecordButton() {
    const b = $("record");
    b.classList.toggle("recording", Boolean(state.recorder));
    b.querySelector("span").textContent = state.recorder ? t("stopRecording") : t("record");
  }

  function shareLink() {
    const code = Share.encode({
      scene: state.scene, N: state.N, k0: state.k0, sigma: state.sigma, g: state.g, display: state.display,
      speed: state.speed, detectorOn: state.detectorOn, detectorX: state.detectorX,
      packet: state.customPacket ? state.packet : null, strokes: state.strokes,
    });
    const url = `${location.href.split("#")[0]}#s=${code}`;
    history.replaceState(null, "", url);
    const done = () => flashHint(t("copied"));
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(url).then(done, () => flashHint(url));
    else flashHint(url);
  }

  function loadFromHash() {
    const m = location.hash.match(/^#s=([\w-]+)$/);
    const o = m && Share.decode(m[1]);
    if (!o || !SCENES[o.scene] || ![128, 256, 512].includes(o.N)) return false;
    if (o.N !== state.N || !state.grid) { setGrid(o.N); $("resolution").value = String(o.N); }
    state.display = ["amplitude", "density", "phase"].includes(o.display) ? o.display : "amplitude";
    state.speed = Math.min(20, Math.max(1, Math.round(o.speed)));
    pressed("#display button", "display", state.display);
    loadScene(o.scene, {
      strokes: o.strokes, packet: o.packet, k0: clamp(o.k0, 0.5, 8), sigma: clamp(o.sigma, 0.8, 5),
      g: clamp(o.g, 0, 500), detectorOn: o.detectorOn, detectorX: o.detectorX,
    });
    return true;
  }
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

  // ---------------------------------------------------------------- controls
  function setRunning(on) {
    state.running = on;
    $("play").textContent = on ? t("pause") : t("play");
  }

  function syncSliders() {
    const set = (id, v) => { $(id).value = v; };
    set("k0", state.k0); set("sigma", state.sigma); set("g", state.g); set("speed", state.speed);
    set("brush", state.brush); set("height", state.height); set("bright", state.bright);
    $("o-k0").textContent = state.k0.toFixed(1);
    $("o-sigma").textContent = state.sigma.toFixed(1);
    $("o-g").textContent = String(state.g);
    $("o-speed").textContent = String(state.speed);
    $("o-brush").textContent = state.brush.toFixed(1);
    $("o-height").textContent = String(state.height);
    $("o-bright").textContent = `${state.bright.toFixed(2)}×`;
  }

  function pressed(selector, key, value) {
    document.querySelectorAll(selector).forEach((b) => b.setAttribute("aria-pressed", String(b.dataset[key] === value)));
  }

  function slider(id, key, after) {
    $(id).addEventListener("input", (ev) => { state[key] = parseFloat(ev.target.value); syncSliders(); if (after) after(); });
  }

  function setView(v) {
    state.view = v;
    pressed("#view button", "view", v);
    $("hint").textContent = defaultHint();
    field.style.cursor = v === "momentum" ? "default" : "crosshair";
  }

  function changeResolution(N) {
    setGrid(N);
    const sc = SCENES[state.scene];
    state.sceneWalls = sc.walls ? sc.walls(state.grid) : new Float64Array(state.grid.size);
    state.baseV = sc.base ? state.grid.fill(sc.base) : new Float64Array(state.grid.size);
    rebuildWalls();                                   // the strokes are re-rasterized exactly on the new grid
    if (state.scene === "bec") { loadScene("bec", { strokes: state.strokes }); return; }
    buildSolver(false);
    relaunch();
  }

  function setLanguage(lang) {
    state.lang = lang;
    try { localStorage.setItem("qwl-lang", lang); } catch (e) { /* no storage */ }
    document.documentElement.lang = lang;
    document.querySelectorAll("[data-i18n]").forEach((el) => { el.textContent = t(el.dataset.i18n); });
    pressed("[data-lang]", "lang", lang);
    $("hint").textContent = defaultHint();
    setRunning(state.running);
    updateRecordButton();
    if (state.solver) { runChecks(); updateSceneInfo(); }
  }

  function wireControls() {
    const scenes = $("scenes");
    for (const name of ORDER) {
      const b = document.createElement("button");
      b.type = "button"; b.dataset.scene = name; b.dataset.i18n = `scene.${name}`;
      b.addEventListener("click", () => loadScene(name));
      scenes.appendChild(b);
    }
    slider("speed", "speed");
    slider("brush", "brush");
    slider("height", "height");
    slider("bright", "bright");
    slider("k0", "k0");
    slider("sigma", "sigma");
    slider("g", "g", () => { buildSolver(true); runChecks(); });
    ["k0", "sigma"].forEach((id) => $(id).addEventListener("change", () => relaunch()));
    $("play").addEventListener("click", () => setRunning(!state.running));
    $("relaunch").addEventListener("click", () => relaunch());
    $("measure").addEventListener("click", measure);
    $("undo").addEventListener("click", undo);
    $("redo").addEventListener("click", redo);
    $("save").addEventListener("click", saveImage);
    $("share").addEventListener("click", shareLink);
    if (typeof MediaRecorder === "undefined" || !field.captureStream) $("record").hidden = true;
    else $("record").addEventListener("click", toggleRecording);
    $("clear").addEventListener("click", () => pushStroke({ clear: true }));
    $("auto").addEventListener("change", (ev) => { state.auto = ev.target.checked; if (!state.auto) state.scale = maxAmplitude(state.solver); });
    $("contours").addEventListener("change", (ev) => { state.contours = ev.target.checked; });
    $("detector-on").addEventListener("change", (ev) => { state.detectorOn = ev.target.checked; resetDetector(); });
    $("resolution").addEventListener("change", (ev) => changeResolution(parseInt(ev.target.value, 10)));
    document.querySelectorAll("#tools button").forEach((b) => b.addEventListener("click", () => {
      state.tool = b.dataset.tool;
      pressed("#tools button", "tool", state.tool);
      if (state.view !== "position") setView("position");
      $("hint").textContent = defaultHint();
    }));
    document.querySelectorAll("#view button").forEach((b) => b.addEventListener("click", () => setView(b.dataset.view)));
    document.querySelectorAll("#display button").forEach((b) => b.addEventListener("click", () => {
      state.display = b.dataset.display;
      pressed("#display button", "display", state.display);
    }));
    document.querySelectorAll("[data-lang]").forEach((b) => b.addEventListener("click", () => setLanguage(b.dataset.lang)));
    window.addEventListener("keydown", (ev) => {
      if (ev.target.closest("input, select, textarea")) return;
      const mod = ev.ctrlKey || ev.metaKey;
      if (mod && (ev.key === "z" || ev.key === "Z")) { ev.preventDefault(); if (ev.shiftKey) redo(); else undo(); }
      else if (mod && (ev.key === "y" || ev.key === "Y")) { ev.preventDefault(); redo(); }
      else if (mod) return;
      else if (ev.code === "Space") { ev.preventDefault(); setRunning(!state.running); }
      else if (ev.key === "r" || ev.key === "R") relaunch();
      else if (ev.key === "m" || ev.key === "M") measure();
      else if (ev.key === "k" || ev.key === "K") setView(state.view === "position" ? "momentum" : "position");
    });
    window.addEventListener("hashchange", loadFromHash);
  }

  // ---------------------------------------------------------------- start
  wireControls();
  setGrid(state.N);
  setLanguage(state.lang);
  syncSliders();
  if (!loadFromHash()) loadScene(state.scene);
  requestAnimationFrame(frame);

  // debugging handle (browser console): qwl.state, qwl.tick(n) advances n frames even in a hidden tab
  globalThis.qwl = { state, tick: (n = 1) => { for (let i = 0; i < n; i++) advance(performance.now()); }, measure, undo, redo };
})();
