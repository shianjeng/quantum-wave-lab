// The JS solver must reproduce the golden files in ../../reference (written by export_reference.py),
// in float64 and in float32, within the tolerances stored in the manifest.
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { test } from "node:test";
import assert from "node:assert/strict";

import "../js/qwave.js";
const { Grid, Solver, FluxDetector, absorptionRate } = globalThis.QWave;

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "reference");
const manifest = JSON.parse(readFileSync(join(ROOT, "manifest.json"), "utf8"));

function readF64(file) {
  const b = readFileSync(join(ROOT, file));
  return new Float64Array(b.buffer, b.byteOffset, b.byteLength / 8);
}
function readC128(file) {                       // interleaved (re, im) -> [re, im]
  const a = readF64(file), n = a.length / 2, re = new Float64Array(n), im = new Float64Array(n);
  for (let i = 0; i < n; i++) { re[i] = a[2 * i]; im[i] = a[2 * i + 1]; }
  return [re, im];
}
function relL2(re, im, [rr, ri]) {
  let d = 0, s = 0;
  for (let i = 0; i < rr.length; i++) { d += (re[i] - rr[i]) ** 2 + (im[i] - ri[i]) ** 2; s += rr[i] ** 2 + ri[i] ** 2; }
  return Math.sqrt(d / s);
}
const absorberOptions = (a) => (a === null ? null : { width: a.width ?? null, widthFrac: a.width_frac ?? 0.08, gammaMax: a.gamma_max ?? 2.0 });

for (const c of manifest.cases) {
  for (const precision of ["float64", "float32"]) {
    test(`${c.name} (${precision})`, () => {
      const grid = new Grid(c.n, c.length);
      const V = readF64(c.files.V);
      if (c.absorber) {
        const mask = absorptionRate(grid, absorberOptions(c.absorber)).map((g) => Math.exp(-g * c.dt));
        const ref = readF64(c.files.absorber_mask);
        const err = Math.max(...mask.map((m, i) => Math.abs(m - ref[i])));
        assert.ok(err < 1e-14, `absorber mask differs by ${err}`);
      }
      const solver = new Solver(grid, V, c.dt, { absorber: absorberOptions(c.absorber), precision, nonlinearity: c.nonlinearity });
      solver.setPsi(...readC128(c.files.psi0));
      const detectors = c.detectors.map((d) => new FluxDetector(grid, d.position, { span: d.span ?? null, method: d.method }));
      detectors.forEach((d) => d.record(solver));
      const tol = c.tolerance[precision];
      for (const cp of c.checkpoints) {
        solver.step(cp.step - solver.nsteps, (s) => detectors.forEach((d) => d.record(s)));
        assert.ok(Math.abs(solver.t - cp.t) < 1e-9, `t = ${solver.t}, expected ${cp.t}`);
        const err = relL2(solver.re, solver.im, readC128(cp.psi));
        assert.ok(err < tol, `step ${cp.step}: relative L2 error ${err.toExponential(2)} > ${tol.toExponential(2)}`);
        assert.ok(Math.abs(solver.norm() - cp.norm) < tol, `step ${cp.step}: norm ${solver.norm()} vs ${cp.norm}`);
        const e = solver.energy();
        assert.ok(Math.abs(e - cp.energy) < tol * Math.max(1, Math.abs(cp.energy)), `step ${cp.step}: energy ${e} vs ${cp.energy}`);
        (cp.flux_transmitted ?? []).forEach((want, k) => {
          const got = detectors[k].transmitted;
          assert.ok(Math.abs(got - want) < tol, `step ${cp.step}: flux ${k} ${got} vs ${want}`);
        });
      }
    });
  }
}
