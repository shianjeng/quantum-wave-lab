import { test } from "node:test";
import assert from "node:assert/strict";
import "../js/qwave.js";
const { Grid, Solver, gaussianPacket, measurePosition, momentumDensity } = globalThis.QWave;

function seeded(seed) {             // small deterministic PRNG (mulberry32)
  return () => {
    seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

test("momentum density: normalized (Parseval) and peaked at k0", () => {
  const g = new Grid([128, 128], [51.2, 51.2]);
  const s = new Solver(g, new Float64Array(g.size), 0.01, { precision: "float64" });
  s.setPsi(...gaussianPacket(g, [2, -3], [2, 2], [3, -1.5]));
  const phi2 = momentumDensity(s), dk2 = ((2 * Math.PI) / g.lx) ** 2;
  let sum = 0, best = 0;
  for (let i = 0; i < phi2.length; i++) { sum += phi2[i] * dk2; if (phi2[i] > phi2[best]) best = i; }
  assert.ok(Math.abs(sum - 1) < 1e-12, `Σ|φ|² dk² = ${sum}`);
  const kx = g.kx[Math.floor(best / g.ny)], ky = g.ky[best % g.ny];
  assert.ok(Math.abs(kx - 3) < 0.07 && Math.abs(ky + 1.5) < 0.07, `peak at (${kx}, ${ky})`);
});

test("position measurement: normalized, outcomes follow |ψ|², local momentum kept", () => {
  const g = new Grid([256], [51.2]);
  const rand = seeded(1);
  let mean = 0;
  const n = 400;
  for (let k = 0; k < n; k++) {
    const s = new Solver(g, new Float64Array(g.size), 0.01, { precision: "float64" });
    s.setPsi(...gaussianPacket(g, [4], [2], [1.5]));
    const [x0] = measurePosition(s, 0.5, rand);
    mean += x0 / n;
    if (k === 0) {
      assert.ok(Math.abs(s.norm() - 1) < 1e-12);
      // collapsed: narrow around x0, and the phase gradient (momentum 1.5) survives
      let m = 0, p = 0;
      for (let i = 0; i < g.nx; i++) m += g.xs[i] * (s.re[i] ** 2 + s.im[i] ** 2) * g.dx;
      const phi2 = momentumDensity(s), dk = (2 * Math.PI) / g.lx;
      for (let i = 0; i < g.nx; i++) p += g.kx[i] * phi2[i] * dk;
      assert.ok(Math.abs(m - x0) < 0.3, `⟨x⟩ = ${m} vs outcome ${x0}`);
      assert.ok(Math.abs(p - 1.5) < 1e-6, `⟨k⟩ = ${p}`);             // a real window leaves ⟨k⟩ unchanged
    }
  }
  // outcomes: |ψ|² (σ = 2) smeared by σm = 0.5, mean 4; standard error ≈ 2.06/√400 ≈ 0.1
  assert.ok(Math.abs(mean - 4) < 0.35, `mean outcome ${mean}`);
});
