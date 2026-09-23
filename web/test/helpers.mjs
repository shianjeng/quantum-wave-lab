// Headless replica of the app's setup (js/app.js): same box, absorber, time step, packet and detector.
import "../js/qwave.js";
import "../js/scenes.js";
const { Grid, Solver, FluxDetector, gaussianPacket } = globalThis.QWave;
const { SCENES } = globalThis.QWaveScenes;

export const L = 51.2, ABSORBER = { width: (51.2 - 40) / 2, gammaMax: 6 }, DT = 0.01, STEPS_PER_FRAME = 6;

/** Solver, detector and grid of a scene as the app builds them; `edit(x, y, V)` may change the walls. */
export function setupScene(name, { N = 256, k0 = null, walls = null } = {}) {
  const sc = SCENES[name], g = new Grid([N, N], [L, L]);
  const W = walls ? g.fill(walls) : sc.walls ? sc.walls(g) : new Float64Array(g.size);
  const B = sc.base ? g.fill(sc.base) : new Float64Array(g.size);
  const V = W.map((w, i) => w + B[i]);
  const s = new Solver(g, V, DT, { absorber: ABSORBER });
  const p = sc.packet, n = Math.hypot(...p.dir), kk = k0 ?? p.k0;
  s.setPsi(...gaussianPacket(g, p.center, [p.sigma, p.sigma * (p.aspect ?? 1)], [(kk * p.dir[0]) / n, (kk * p.dir[1]) / n], Float32Array));
  const det = sc.detector === null || sc.detector === undefined ? null : new FluxDetector(g, sc.detector, { method: "spectral" });
  if (det) det.record(s);
  return { g, s, det, W };
}

/** Run until time T with the detector sampled once per frame, as in the app. */
export function runUntil({ s, det }, T) {
  while (s.t < T - 1e-9) {
    s.step(STEPS_PER_FRAME);
    if (det) det.record(s);
  }
  return det ? det.transmitted : null;
}
