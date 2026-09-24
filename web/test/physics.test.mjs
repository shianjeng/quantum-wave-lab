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

test("flux profile: Σ profile·dy equals the transmitted probability", async () => {
  const { FluxDetector } = globalThis.QWave;
  const g = new Grid([64, 64], [25.6, 25.6]);
  const s = new Solver(g, new Float64Array(g.size), 0.01, { precision: "float64" });
  s.setPsi(...gaussianPacket(g, [-5, 1], [1.5, 1.5], [3, 0.5]));
  const d = new FluxDetector(g, 2, { method: "fd4" });
  d.record(s);
  for (let i = 0; i < 300; i++) { s.step(1); d.record(s); }
  const sum = d.profile.reduce((a, b) => a + b, 0) * g.dy;
  assert.ok(d.transmitted > 0.3, `transmitted ${d.transmitted}`);
  assert.ok(Math.abs(sum - d.transmitted) < 1e-12, `${sum} vs ${d.transmitted}`);
  // the packet moves up (k_y > 0): the profile's centre of mass lies above its start y = 1
  let m = 0; for (let j = 0; j < g.ny; j++) m += g.ys[j] * d.profile[j] * g.dy;
  assert.ok(m / d.transmitted > 1, `centre ${m / d.transmitted}`);
});

test("sampleProfile draws from the positive part of a profile", () => {
  const { sampleProfile } = globalThis.QWave;
  const ys = Float64Array.from({ length: 40 }, (_, j) => -10 + j * 0.5);
  const profile = ys.map((y) => (Math.abs(y - 3) < 1 ? 2 : Math.abs(y + 5) < 1 ? 1 : -0.3));   // two plateaus, 2:1
  const rand = (() => { let s = 7; return () => ((s = (s * 16807) % 2147483647) / 2147483647); })();
  let near3 = 0, nearMinus5 = 0;
  for (let i = 0; i < 3000; i++) {
    const y = sampleProfile(profile, ys, rand);
    if (Math.abs(y - 3) < 1.25) near3++; else if (Math.abs(y + 5) < 1.25) nearMinus5++; else assert.fail(`sample at ${y}`);
  }
  assert.ok(Math.abs(near3 / nearMinus5 - 2) < 0.25, `ratio ${near3 / nearMinus5}`);
  assert.equal(sampleProfile(ys.map(() => -1), ys, rand), null);
  assert.ok(Math.abs(sampleProfile(profile, ys, rand, 0, 10) - 3) < 1.25, "range restriction");
});
