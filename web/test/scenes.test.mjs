// Every scene must build on every grid the app offers, start without warnings, and the numbers quoted
// in the scene descriptions (SCENES[name].facts) must be what the app actually measures.
import { test } from "node:test";
import assert from "node:assert/strict";
import { setupScene, runUntil } from "./helpers.mjs";

const { SCENES, ORDER, softBox } = globalThis.QWaveScenes;

for (const N of [128, 256, 512]) {
  test(`all scenes build and start cleanly on ${N}²`, () => {
    for (const name of ORDER) {
      if (SCENES[name].bec) continue;                   // prepared by imaginary time in the page
      const { s } = setupScene(name, { N });
      const issues = s.diagnose();
      // on the coarse 128² grid (k_max = 7.9) a k0 = 4 packet is close to the momentum limit: allowed there
      const serious = issues.filter((is) => !(N === 128 && is.code === "kmax"));
      assert.deepEqual(serious, [], `${name} on ${N}²: ${JSON.stringify(issues)}`);
    }
  });
}

test("slit scenes are mirror-symmetric about y = 0 on every grid", () => {
  for (const N of [128, 256, 512]) {
    for (const name of ["singleSlit", "doubleSlit", "grating"]) {
      const { g, W } = setupScene(name, { N });
      for (let i = 0; i < g.nx; i++)
        for (let j = 1; j < g.ny; j++) assert.equal(W[i * g.ny + j], W[i * g.ny + (g.ny - j)], `${name} ${N}`);
    }
  }
});

const close = (got, want, tol = 0.006) => assert.ok(Math.abs(got - want) < tol, `measured ${got.toFixed(4)}, text says ${want}`);

test("tunneling: transmission quoted in the description", () => {
  close(runUntil(setupScene("tunneling"), 10), SCENES.tunneling.facts.T);
});

test("resonant tunneling: pair, one barrier erased, detuned", () => {
  const f = SCENES.resonant.facts;
  close(runUntil(setupScene("resonant"), 25), f.T2);
  close(runUntil(setupScene("resonant", { walls: (x) => 10 * softBox(x, 0.8, 1.2, 0.1) }), 25), f.T1);
  close(runUntil(setupScene("resonant", { k0: 2.8 }), 25), f.Tdetuned);
});
