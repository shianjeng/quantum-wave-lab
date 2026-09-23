/*
 * Scene definitions for the app (plain script: defines globalThis.QWaveScenes; needs QWave).
 *
 * A scene gives erasable walls (walls(grid)), a smooth background that is not erased (base(x, y)),
 * a packet {center, dir, k0, sigma, aspect}, a detector x position (or null) and an interaction g.
 * Coordinates are physical (the view is -20…20); dimensions are snapped to the grid where needed
 * so openings stay mirror-symmetric (qwave.potentials.slits) at every resolution.
 */
(function (root) {
  "use strict";
  const { centeredCells } = root.QWave;

  const softBox = (u, lo, hi, e) => 0.5 * (Math.tanh((u - lo) / e) - Math.tanh((u - hi) / e));

  /** Fraction of a cell at distance d from the edge of a disk of radius R that lies inside (1-cell ramp). */
  const coverage = (d, R, dx) => Math.min(Math.max((R - d) / dx + 0.5, 0), 1);

  /** Wall of the given thickness at x0 (whole cells) with openings; returns V (Float64Array). */
  function wallWithOpenings(g, x0, thickness, openings, height) {
    const nt = Math.max(1, Math.round(thickness / g.dx));
    return g.fill((x, y) => {
      const inWall = x >= x0 - (nt * g.dx) / 2 - 1e-9 && x < x0 + (nt * g.dx) / 2 - 1e-9;
      return inWall && !openings.some(([c, w]) => centeredCells(y, c, w, g.dy)) ? height : 0;
    });
  }

  /** Odd number of cells closest to `width` (an opening centred on a grid line needs an odd count). */
  const oddCells = (width, d) => { const n = Math.max(1, Math.round(width / d)); return n % 2 ? n : n + 1; };

  const SCENES = {
    tunneling: {
      walls: (g) => g.fill((x) => 10 * softBox(x, -0.5, 0.5, 0.25)),             // ⟨E⟩ ≈ 8 < V0 = 10
      packet: { center: [-8, 0], dir: [1, 0], k0: 4, sigma: 3 }, detector: 3,
      facts: { T: 0.16 },                     // checked by test/scenes.test.mjs
    },
    resonant: {
      // two barriers (V0 = 10, width 0.4, gap 1.6): on resonance (k0 = 3) the pair transmits more than one alone
      walls: (g) => g.fill((x) => 10 * (softBox(x, -1.2, -0.8, 0.1) + softBox(x, 0.8, 1.2, 0.1))),
      packet: { center: [-8, 0], dir: [1, 0], k0: 3, sigma: 3 }, detector: 4,
      facts: { T2: 0.38, T1: 0.25, Tdetuned: 0.19 },   // pair, one barrier erased, k0 = 2.8 (test/scenes.test.mjs)
    },
    singleSlit: {
      walls: (g) => wallWithOpenings(g, -3, 0.4, [[0, oddCells(1.6, g.dy) * g.dy]], 100),
      packet: { center: [-11, 0], dir: [1, 0], k0: 4, sigma: 2.2 }, detector: 1,
    },
    doubleSlit: {
      walls: (g) => {
        // openings of 1.2 and separation 5, snapped so they sit exactly symmetric on this grid
        const nw = Math.max(1, Math.round(1.2 / g.dy));
        let ns = Math.round(5 / g.dy);
        if ((nw + ns) % 2 === 0) ns += 1;
        const w = nw * g.dy, c = (ns * g.dy) / 2;
        return wallWithOpenings(g, -3, 0.4, [[c, w], [-c, w]], 100);
      },
      packet: { center: [-11, 0], dir: [1, 0], k0: 4, sigma: 2.2 }, detector: 1,
    },
    grating: {
      // 7 openings, period ≈ 2.4: sharp diffraction orders at k_y = 2πn/2.4 (see the momentum view)
      walls: (g) => {
        const w = oddCells(0.6, g.dy) * g.dy, p = Math.round(2.4 / g.dy) * g.dy;
        return wallWithOpenings(g, -3, 0.4, [-3, -2, -1, 0, 1, 2, 3].map((m) => [m * p, w]), 100);
      },
      packet: { center: [-11, 0], dir: [1, 0], k0: 4, sigma: 2.2, aspect: 2.2 }, detector: 1,   // σ_y ≈ 4.8: lights all seven
    },
    scattering: {
      walls: (g) => g.fill((x, y) => 100 * coverage(Math.hypot(x, y), 2, g.dx)),
      packet: { center: [-11, 0.8], dir: [1, 0], k0: 4, sigma: 2.5 }, detector: 6,
    },
    harmonic: {
      base: (x, y) => 0.5 * 0.35 ** 2 * (x * x + y * y),                         // coherent state on a circle
      packet: { center: [-8, 0], dir: [0, 1], k0: 2.8, sigma: 1.2 }, detector: null,
    },
    corral: {
      // a ring wall (radius 13) the packet cannot leave: it bounces and fills the corral with interference
      walls: (g) => g.fill((x, y) => { const r = Math.hypot(x, y); return 100 * Math.min(coverage(r, 13.6, g.dx), 1 - coverage(r, 13, g.dx)); }),
      packet: { center: [-5, -2], dir: [1, 0.35], k0: 3.5, sigma: 1.6 }, detector: null,
    },
    bec: { bec: true, detector: null },
    free: { packet: { center: [-10, 0], dir: [1, 0], k0: 4, sigma: 2.5 }, detector: 5 },
  };

  /** The trap the BEC scene is prepared in (released at t = 0). */
  const becTrap = (x, y) => 0.5 * 0.25 ** 2 * (x * x + y * y) + 8 * Math.exp(-(x * x) / 2);

  root.QWaveScenes = { SCENES, ORDER: Object.keys(SCENES), becTrap, softBox, coverage };
})(globalThis);
