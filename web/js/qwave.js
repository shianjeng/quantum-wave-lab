/*
 * qwave.js: browser port of the Python reference solver (qwave 0.3.0).
 *
 *   i ∂ψ/∂t = [ -∇²/2 + V + g|ψ|² ] ψ      (ħ = m = 1), split-step Fourier (Strang), periodic box,
 *   complex absorbing layer, optional Gross–Pitaevskii nonlinearity.
 *
 * Arrays follow reference/README.md: C order with axis 0 (x) slowest, element [ix, iy] at ix*ny + iy.
 * A 1D grid is a 2D grid with ny = 1. ψ is stored as two real arrays (re, im) of the chosen precision
 * (Float32Array like the browser app, Float64Array for the reference tests); the FFT itself works in
 * double precision on one row or column at a time.
 *
 * Plain script (no modules): it defines globalThis.QWave, so index.html can load it with <script>
 * (even from file://) and Node can import it for the tests.
 */
(function (root) {
  "use strict";

  // ------------------------------------------------------------------ FFT (radix 2, in place)
  class FFT {
    constructor(n) {
      if (n < 1 || (n & (n - 1)) !== 0) throw new Error(`FFT size must be a power of 2, got ${n}`);
      this.n = n;
      this.rev = new Uint32Array(n);
      const bits = Math.round(Math.log2(n));
      for (let i = 0; i < n; i++) {
        let r = 0;
        for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b);
        this.rev[i] = r;
      }
      this.cos = new Float64Array(n / 2 || 1);
      this.sin = new Float64Array(n / 2 || 1);
      for (let k = 0; k < n / 2; k++) {
        this.cos[k] = Math.cos((2 * Math.PI * k) / n);
        this.sin[k] = Math.sin((2 * Math.PI * k) / n);
      }
    }

    /** Unnormalized DFT of (re, im) in place; sign -1 forward (e^{-2πi jk/n}), +1 inverse. */
    transform(re, im, sign) {
      const n = this.n, rev = this.rev;
      if (n === 1) return;
      for (let i = 0; i < n; i++) {
        const j = rev[i];
        if (j > i) {
          let t = re[i]; re[i] = re[j]; re[j] = t;
          t = im[i]; im[i] = im[j]; im[j] = t;
        }
      }
      for (let size = 2; size <= n; size <<= 1) {
        const half = size >> 1, step = n / size;
        for (let start = 0; start < n; start += size) {
          for (let j = 0, k = 0; j < half; j++, k += step) {
            const wr = this.cos[k], wi = sign * this.sin[k];
            const a = start + j, b = a + half;
            const tr = re[b] * wr - im[b] * wi;
            const ti = re[b] * wi + im[b] * wr;
            re[b] = re[a] - tr; im[b] = im[a] - ti;
            re[a] += tr; im[a] += ti;
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------ grid
  function fftFreq(n, d) {
    const k = new Float64Array(n);
    for (let j = 0; j < n; j++) k[j] = (2 * Math.PI * (j < Math.ceil(n / 2) ? j : j - n)) / (n * d);
    return k;
  }

  class Grid {
    /** n = [nx] or [nx, ny], length = [lx] or [lx, ly]; coordinates x_j = -L/2 + j·dx. */
    constructor(n, length) {
      if (n.length !== length.length || n.length < 1 || n.length > 2) throw new Error("need 1 or 2 axes");
      this.ndim = n.length;
      this.nx = n[0]; this.ny = n.length === 2 ? n[1] : 1;
      this.lx = length[0]; this.ly = n.length === 2 ? length[1] : 1;
      this.dx = this.lx / this.nx;
      this.dy = this.ndim === 2 ? this.ly / this.ny : 1;
      this.size = this.nx * this.ny;
      this.dV = this.ndim === 2 ? this.dx * this.dy : this.dx;
      this.xs = Float64Array.from({ length: this.nx }, (_, j) => -this.lx / 2 + j * this.dx);
      this.ys = this.ndim === 2 ? Float64Array.from({ length: this.ny }, (_, j) => -this.ly / 2 + j * this.dy)
                                : new Float64Array(1);
      this.kx = fftFreq(this.nx, this.dx);
      this.ky = this.ndim === 2 ? fftFreq(this.ny, this.dy) : new Float64Array(1);
      this.kmax = this.ndim === 2 ? [Math.PI / this.dx, Math.PI / this.dy] : [Math.PI / this.dx];
      this.fftX = new FFT(this.nx);
      this.fftY = new FFT(this.ny);
      this._bufX = [new Float64Array(this.nx), new Float64Array(this.nx)];
    }

    /** Array of f(x, y) over the grid (Float64Array). */
    fill(f) {
      const out = new Float64Array(this.size);
      for (let i = 0; i < this.nx; i++)
        for (let j = 0; j < this.ny; j++) out[i * this.ny + j] = f(this.xs[i], this.ys[j]);
      return out;
    }

    /** In-place 2D DFT of (re, im) (any typed arrays), sign -1 forward / +1 inverse, unnormalized. */
    fft2(re, im, sign) {
      const { nx, ny } = this;
      if (ny > 1) {
        const r = new Float64Array(ny), m = new Float64Array(ny);
        for (let i = 0; i < nx; i++) {
          const o = i * ny;
          for (let j = 0; j < ny; j++) { r[j] = re[o + j]; m[j] = im[o + j]; }
          this.fftY.transform(r, m, sign);
          for (let j = 0; j < ny; j++) { re[o + j] = r[j]; im[o + j] = m[j]; }
        }
      }
      const [r, m] = this._bufX;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) { r[i] = re[i * ny + j]; m[i] = im[i * ny + j]; }
        this.fftX.transform(r, m, sign);
        for (let i = 0; i < nx; i++) { re[i * ny + j] = r[i]; im[i * ny + j] = m[i]; }
      }
    }

    norm(re, im) {
      let s = 0;
      for (let i = 0; i < this.size; i++) s += re[i] * re[i] + im[i] * im[i];
      return s * this.dV;
    }
  }

  // ------------------------------------------------------------------ initial states and potentials
  /** ψ ∝ exp(-(r-r0)²/(4σ²) + i k0·r), normalized; center/sigma/k0 are [x, y] (or [x] in 1D). */
  function gaussianPacket(grid, center, sigma, k0, Float = Float64Array) {
    const re = new Float64Array(grid.size), im = new Float64Array(grid.size);
    const [cx, cy = 0] = center, [sx, sy = 1] = sigma, [kx, ky = 0] = k0;
    for (let i = 0; i < grid.nx; i++) {
      const x = grid.xs[i];
      for (let j = 0; j < grid.ny; j++) {
        const y = grid.ys[j];
        let arg = -((x - cx) ** 2) / (4 * sx * sx), ph = kx * x;
        if (grid.ndim === 2) { arg -= ((y - cy) ** 2) / (4 * sy * sy); ph += ky * y; }
        const a = Math.exp(arg), idx = i * grid.ny + j;
        re[idx] = a * Math.cos(ph); im[idx] = a * Math.sin(ph);
      }
    }
    const s = 1 / Math.sqrt(grid.norm(re, im));
    return [Float.from(re, (v) => v * s), Float.from(im, (v) => v * s)];
  }

  /** γ(r): 0 inside, sin² ramp to gammaMax over `width` (absolute) or widthFrac·L at each edge. */
  function absorptionRate(grid, { width = null, widthFrac = 0.08, gammaMax = 2.0 } = {}) {
    const ramp = (u, L) => {
      const w = width === null ? widthFrac * L : width;
      const e = Math.min(Math.max((w - (L / 2 - Math.abs(u))) / w, 0), 1);
      return gammaMax * Math.sin((Math.PI / 2) * e) ** 2;
    };
    return grid.fill((x, y) => (grid.ndim === 2 ? Math.max(ramp(x, grid.lx), ramp(y, grid.ly)) : ramp(x, grid.lx)));
  }

  /** Openings of `width` symmetric about each centre (qwave.potentials._centered_cells). */
  function centeredCells(u, center, width, du) {
    const n = width / du, m = (2 * center) / du;
    if (Math.abs(n - Math.round(n)) > 1e-9 || Math.abs(m - Math.round(m)) > 1e-9 ||
        (Math.round(n) + Math.round(m)) % 2 === 0) {
      throw new Error(`an opening of width ${width} centred at ${center} (du = ${du}) cannot sit symmetrically ` +
                      `on the grid: (2·centre + width)/du must be an odd integer`);
    }
    return Math.abs(u - center) < width / 2 - 1e-9 * du;
  }

  // ------------------------------------------------------------------ solver
  class Solver {
    /**
     * grid, V (array of grid.size), dt; options: absorber (null or {width|widthFrac, gammaMax}),
     * precision "float32" | "float64", nonlinearity g, renormalizeEvery (null or N).
     */
    constructor(grid, V, dt, { absorber = null, precision = "float32", nonlinearity = 0,
                               renormalizeEvery = null } = {}) {
      this.grid = grid;
      this.Float = precision === "float64" ? Float64Array : Float32Array;
      this.absorberOptions = absorber;
      this.nonlinearity = nonlinearity;
      this.renormalizeEvery = renormalizeEvery;
      if (renormalizeEvery !== null && absorber !== null) throw new Error("renormalizeEvery needs absorber = null");
      this.re = new this.Float(grid.size);
      this.im = new this.Float(grid.size);
      this.t = 0;
      this.nsteps = 0;
      this.shifted = false;          // true while (re, im) hold e^{-iVdt/2}ψ between merged steps
      this.V = Float64Array.from(V);
      this.setDt(dt);
    }

    setPsi(re, im) {
      this.re.set(re); this.im.set(im);
      this.shifted = false;
    }

    setDt(dt) {
      this.dt = dt;
      const { grid } = this, N = grid.size;
      this.kinRe = new this.Float(N); this.kinIm = new this.Float(N);
      for (let i = 0; i < grid.nx; i++) {
        for (let j = 0; j < grid.ny; j++) {
          const k2 = grid.kx[i] ** 2 + (grid.ndim === 2 ? grid.ky[j] ** 2 : 0), ph = -0.5 * k2 * dt;
          this.kinRe[i * grid.ny + j] = Math.cos(ph) / N;      // 1/N of the inverse FFT folded in
          this.kinIm[i * grid.ny + j] = Math.sin(ph) / N;
        }
      }
      this.gamma = this.absorberOptions ? absorptionRate(grid, this.absorberOptions) : null;
      this.mask = this.gamma ? this.gamma.map((g) => Math.exp(-g * dt)) : null;
      this.setPotential(this.V);
    }

    /** Replace V (e.g. after the user painted); keeps ψ. */
    setPotential(V) {
      if (this.shifted) this._unshift();
      this.V = Float64Array.from(V);
      const N = this.grid.size, F = this.Float;
      this.halfRe = new F(N); this.halfIm = new F(N);
      this.halfAbsRe = new F(N); this.halfAbsIm = new F(N);
      this.fullAbsRe = new F(N); this.fullAbsIm = new F(N);
      for (let i = 0; i < N; i++) {
        const ph = -0.5 * this.V[i] * this.dt, a = this.mask ? this.mask[i] : 1;
        const c = Math.cos(ph), s = Math.sin(ph), c2 = Math.cos(2 * ph), s2 = Math.sin(2 * ph);
        this.halfRe[i] = c; this.halfIm[i] = s;
        this.halfAbsRe[i] = a * c; this.halfAbsIm[i] = a * s;
        this.fullAbsRe[i] = a * c2; this.fullAbsIm[i] = a * s2;
      }
    }

    _mul(fr, fi) {
      const { re, im } = this;
      for (let i = 0; i < re.length; i++) {
        const r = re[i], m = im[i];
        re[i] = r * fr[i] - m * fi[i];
        im[i] = r * fi[i] + m * fr[i];
      }
    }

    _kinetic() {
      this.grid.fft2(this.re, this.im, -1);
      this._mul(this.kinRe, this.kinIm);
      this.grid.fft2(this.re, this.im, +1);
    }

    _nonlinearHalf() {
      const { re, im, V, dt } = this, g = this.nonlinearity;
      for (let i = 0; i < re.length; i++) {
        const r = re[i], m = im[i], ph = -0.5 * dt * (V[i] + g * (r * r + m * m));
        const c = Math.cos(ph), s = Math.sin(ph);
        re[i] = r * c - m * s; im[i] = r * s + m * c;
      }
    }

    _unshift() {                       // e^{+iVdt/2} (e^{-iVdt/2} ψ) = ψ
      const { re, im, halfRe, halfIm } = this;
      for (let i = 0; i < re.length; i++) {
        const r = re[i], m = im[i];
        re[i] = r * halfRe[i] + m * halfIm[i];
        im[i] = m * halfRe[i] - r * halfIm[i];
      }
      this.shifted = false;
    }

    /** ψ at one index, correct also between merged steps. */
    valueAt(i) {
      if (!this.shifted) return [this.re[i], this.im[i]];
      const r = this.re[i], m = this.im[i];
      return [r * this.halfRe[i] + m * this.halfIm[i], m * this.halfRe[i] - r * this.halfIm[i]];
    }

    _renormalize() {
      const s = 1 / Math.sqrt(this.grid.norm(this.re, this.im));
      for (let i = 0; i < this.re.length; i++) { this.re[i] *= s; this.im[i] *= s; }
    }

    /**
     * Advance nsteps. observer(solver), if given, runs after every `every`-th step; it reads ψ with
     * valueAt(), which is correct although the potential half-steps are merged in between.
     */
    step(nsteps = 1, observer = null, every = 1) {
      const after = () => {
        this.t += this.dt;
        this.nsteps += 1;
        if (this.renormalizeEvery && this.nsteps % this.renormalizeEvery === 0) this._renormalize();
        if (observer && this.nsteps % every === 0) observer(this);
      };
      if (this.shifted) this._unshift();
      if (this.nonlinearity) {
        for (let n = 0; n < nsteps; n++) {
          this._nonlinearHalf(); this._kinetic(); this._nonlinearHalf();
          if (this.mask) { const { re, im, mask } = this; for (let i = 0; i < re.length; i++) { re[i] *= mask[i]; im[i] *= mask[i]; } }
          after();
        }
        return;
      }
      // merged half-steps: e^{-iVdt/2} K e^{-iVdt/2} · e^{-iVdt/2} K e^{-iVdt/2} ... (see qwave/solver.py)
      this._mul(this.halfRe, this.halfIm);
      for (let n = 0; n < nsteps; n++) {
        this._kinetic();
        if (n < nsteps - 1) {
          this._mul(this.fullAbsRe, this.fullAbsIm);
          this.shifted = true;
        } else {
          this._mul(this.halfAbsRe, this.halfAbsIm);
          this.shifted = false;
        }
        after();
      }
    }

    norm() {
      if (this.shifted) this._unshift();
      return this.grid.norm(this.re, this.im);
    }

    /** E = [⟨T⟩ + ⟨V⟩ + (g/2)∫|ψ|⁴] / norm. */
    energy() {
      if (this.shifted) this._unshift();
      const { grid, re, im, V } = this, N = grid.size;
      const fr = Float64Array.from(re), fi = Float64Array.from(im);
      grid.fft2(fr, fi, -1);
      let ekin = 0, epot = 0, eint = 0;
      for (let i = 0; i < grid.nx; i++) {
        for (let j = 0; j < grid.ny; j++) {
          const idx = i * grid.ny + j, k2 = grid.kx[i] ** 2 + (grid.ndim === 2 ? grid.ky[j] ** 2 : 0);
          ekin += k2 * (fr[idx] ** 2 + fi[idx] ** 2);
        }
      }
      for (let i = 0; i < N; i++) {
        const rho = re[i] * re[i] + im[i] * im[i];
        epot += V[i] * rho; eint += rho * rho;
      }
      const e = (0.5 * ekin * grid.dV) / N + epot * grid.dV + 0.5 * this.nonlinearity * eint * grid.dV;
      return e / grid.norm(re, im);
    }

    /** dt · Σ k_max² / 2: must stay below π when g ≠ 0 (Weideman & Herbst 1986). */
    stabilityNumber() {
      return (this.dt * this.grid.kmax.reduce((s, k) => s + k * k, 0)) / 2;
    }

    /**
     * Same checks as Solver.diagnose in Python. Returns [{code, ...numbers}] (empty = OK), so the
     * app can translate the messages.
     */
    diagnose({ tol = 1e-6, maxLossRate = 1e-4 } = {}) {
      if (this.shifted) this._unshift();
      const { grid, re, im, V } = this, N = grid.size, g = this.nonlinearity, issues = [];
      let rhoMax = 0;
      const rho = new Float64Array(N);
      for (let i = 0; i < N; i++) { rho[i] = re[i] * re[i] + im[i] * im[i]; if (rho[i] > rhoMax) rhoMax = rho[i]; }
      let vmax = 0;
      for (let i = 0; i < N; i++) if (rho[i] > tol * rhoMax) vmax = Math.max(vmax, Math.abs(V[i] + g * rho[i]));
      if (vmax * this.dt > 1 + 1e-9) issues.push({ code: "vdt", value: vmax * this.dt, dtMax: 1 / vmax });
      if (g && this.stabilityNumber() > Math.PI) {
        issues.push({ code: "nonlinear", value: this.stabilityNumber(),
                      dtMax: (2 * Math.PI) / grid.kmax.reduce((s, k) => s + k * k, 0) });
      }
      const fr = Float64Array.from(re), fi = Float64Array.from(im);
      grid.fft2(fr, fi, -1);
      let total = 0;
      const beyond = [0, 0];
      for (let i = 0; i < grid.nx; i++) {
        for (let j = 0; j < grid.ny; j++) {
          const p = fr[i * grid.ny + j] ** 2 + fi[i * grid.ny + j] ** 2;
          total += p;
          if (Math.abs(grid.kx[i]) > (2 / 3) * grid.kmax[0]) beyond[0] += p;
          if (grid.ndim === 2 && Math.abs(grid.ky[j]) > (2 / 3) * grid.kmax[1]) beyond[1] += p;
        }
      }
      for (let a = 0; a < grid.ndim; a++) {
        if (beyond[a] / total > tol) issues.push({ code: "kmax", value: beyond[a] / total, axis: a, kmax: grid.kmax[a] });
      }
      if (this.gamma && this.t === 0) {
        let lost = 0, sum = 0;
        for (let i = 0; i < N; i++) { lost += 2 * this.gamma[i] * rho[i]; sum += rho[i]; }
        if (lost / sum > maxLossRate) issues.push({ code: "absorber", value: lost / sum });
      }
      return issues;
    }
  }

  // ------------------------------------------------------------------ flux detector
  const FD4 = [[-2, 1 / 12], [-1, -8 / 12], [1, 8 / 12], [2, -1 / 12]];

  class FluxDetector {
    /**
     * Probability crossing the line x = position (+x counts positive), optionally only lo <= y < hi.
     * method "spectral" (∂ψ/∂x by an FFT along x of each column) or "fd4" (5-point stencil).
     */
    constructor(grid, position, { span = null, method = "spectral" } = {}) {
      this.grid = grid; this.method = method;
      let best = 0;
      for (let i = 1; i < grid.nx; i++) if (Math.abs(grid.xs[i] - position) < Math.abs(grid.xs[best] - position)) best = i;
      this.index = best;
      this.position = grid.xs[best];
      this.columns = [];
      for (let j = 0; j < grid.ny; j++) {
        const y = grid.ys[j];
        if (span === null || (y >= span[0] && y < span[1])) this.columns.push(j);
      }
      this._r = new Float64Array(grid.nx); this._m = new Float64Array(grid.nx);
      // e^{+i k x_index} for the spectral derivative evaluated at one point
      this._phRe = new Float64Array(grid.nx); this._phIm = new Float64Array(grid.nx);
      for (let k = 0; k < grid.nx; k++) {
        const ph = (2 * Math.PI * k * best) / grid.nx;
        this._phRe[k] = Math.cos(ph); this._phIm[k] = Math.sin(ph);
      }
      this.reset();
    }

    reset() {
      this.lastT = null; this.lastFlux = 0; this.transmitted = 0;       // only the last sample is needed
      // per column: current density j_x(y) now and its time integral ∫ j_x dt (Σ profile·dy = transmitted)
      this.current = new Float64Array(this.grid.ny);
      this.profile = new Float64Array(this.grid.ny);
      this._prevCurrent = null;
    }

    flux(solver) {
      const { grid, index } = this, ny = grid.ny, nx = grid.nx;
      let total = 0;
      for (const j of this.columns) {
        const [r0, m0] = solver.valueAt(index * ny + j);
        let dr = 0, dm = 0;
        if (this.method === "fd4") {
          for (const [o, c] of FD4) {
            const [r, m] = solver.valueAt((((index + o) % nx) + nx) % nx * ny + j);
            dr += c * r; dm += c * m;
          }
          dr /= grid.dx; dm /= grid.dx;
        } else {
          const R = this._r, M = this._m;
          for (let i = 0; i < nx; i++) [R[i], M[i]] = solver.valueAt(i * ny + j);
          grid.fftX.transform(R, M, -1);
          for (let k = 0; k < nx; k++) {       // Σ i k φ_k e^{+ikx} / nx
            const kr = -grid.kx[k] * M[k], km = grid.kx[k] * R[k];
            dr += kr * this._phRe[k] - km * this._phIm[k];
            dm += kr * this._phIm[k] + km * this._phRe[k];
          }
          dr /= nx; dm /= nx;
        }
        const jx = r0 * dm - m0 * dr;          // Im(ψ* ∂ψ)
        this.current[j] = jx;
        total += jx;
      }
      return total * (grid.ndim === 2 ? grid.dy : 1);
    }

    /** Add a sample at solver.t (trapezoid rule over the samples, as in Python). */
    record(solver) {
      const phi = this.flux(solver);
      if (this.lastT !== null) {
        const dt = solver.t - this.lastT;
        this.transmitted += 0.5 * (phi + this.lastFlux) * dt;
        for (const j of this.columns) this.profile[j] += 0.5 * (this.current[j] + this._prevCurrent[j]) * dt;
      }
      this.lastT = solver.t;
      this.lastFlux = phi;
      this._prevCurrent = Float64Array.from(this.current);
      return phi;
    }
  }

  // ------------------------------------------------------------------ imaginary time (ground states)
  /**
   * Ground state of V (+ g|ψ|²) by imaginary-time propagation, run in chunks so a page stays
   * responsive: returns an iterator; each next() does `chunk` steps and yields {mu, converged}.
   * Both nonlinear half-steps use the density at the start of the step (qwave/eigen.py).
   */
  function* groundStateSteps(grid, V, { dtau = 0.005, nonlinearity = 0, tol = 1e-9, chunk = 20, maxSteps = 20000 } = {}) {
    const N = grid.size, g = nonlinearity;
    const re = new Float64Array(N), im = new Float64Array(N);
    for (let i = 0; i < grid.nx; i++)
      for (let j = 0; j < grid.ny; j++)
        re[i * grid.ny + j] = Math.exp(-((grid.xs[i] / (0.25 * grid.lx)) ** 2) - (grid.ndim === 2 ? (grid.ys[j] / (0.25 * grid.ly)) ** 2 : 0));
    const kin = new Float64Array(N);
    for (let i = 0; i < grid.nx; i++)
      for (let j = 0; j < grid.ny; j++)
        kin[i * grid.ny + j] = Math.exp(-0.5 * dtau * (grid.kx[i] ** 2 + (grid.ndim === 2 ? grid.ky[j] ** 2 : 0))) / N;
    const half = new Float64Array(N);
    let muOld = Infinity;
    for (let it = 1; it <= maxSteps; it++) {
      for (let i = 0; i < N; i++) half[i] = Math.exp(-0.5 * dtau * (V[i] + g * (re[i] * re[i] + im[i] * im[i])));
      for (let i = 0; i < N; i++) { re[i] *= half[i]; im[i] *= half[i]; }
      grid.fft2(re, im, -1);
      for (let i = 0; i < N; i++) { re[i] *= kin[i]; im[i] *= kin[i]; }
      grid.fft2(re, im, +1);
      for (let i = 0; i < N; i++) { re[i] *= half[i]; im[i] *= half[i]; }
      const s = 1 / Math.sqrt(grid.norm(re, im));
      for (let i = 0; i < N; i++) { re[i] *= s; im[i] *= s; }
      if (it % chunk === 0) {
        const solver = new Solver(grid, V, dtau, { precision: "float64", nonlinearity: g });
        solver.setPsi(re, im);
        const mu = solver.energy() + 0.5 * g * (() => { let q = 0; for (let i = 0; i < N; i++) q += (re[i] ** 2 + im[i] ** 2) ** 2; return q * grid.dV; })();
        const converged = Math.abs(mu - muOld) < tol;
        muOld = mu;
        yield { mu, converged, re, im, steps: it };
        if (converged) return;
      }
    }
  }

  // ------------------------------------------------------------------ measurement and momentum space
  /**
   * Gaussian position measurement with resolution sigmaM: the outcome r0 is drawn from |ψ|² smeared by
   * a Gaussian of width sigmaM (sample r from |ψ|², add N(0, sigmaM²)); then ψ ← e^{-(r-r0)²/(4σm²)} ψ,
   * renormalized. The window is real, so the local momentum (phase) is kept. `random` returns [0, 1).
   * Returns r0 as [x, y] (y = 0 in 1D).
   */
  function measurePosition(solver, sigmaM, random = Math.random) {
    if (solver.shifted) solver._unshift();
    const { grid, re, im } = solver, N = grid.size;
    let total = 0;
    for (let i = 0; i < N; i++) total += re[i] * re[i] + im[i] * im[i];
    let u = random() * total, idx = N - 1;
    for (let i = 0; i < N; i++) { u -= re[i] * re[i] + im[i] * im[i]; if (u <= 0) { idx = i; break; } }
    const gauss = () => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());
    const x0 = grid.xs[Math.floor(idx / grid.ny)] + sigmaM * gauss();
    const y0 = grid.ndim === 2 ? grid.ys[idx % grid.ny] + sigmaM * gauss() : 0;
    for (let i = 0; i < grid.nx; i++) {
      for (let j = 0; j < grid.ny; j++) {
        const d2 = (grid.xs[i] - x0) ** 2 + (grid.ndim === 2 ? (grid.ys[j] - y0) ** 2 : 0);
        const w = Math.exp(-d2 / (4 * sigmaM * sigmaM)), k = i * grid.ny + j;
        re[k] *= w; im[k] *= w;
      }
    }
    const s = 1 / Math.sqrt(grid.norm(re, im));
    for (let i = 0; i < N; i++) { re[i] *= s; im[i] *= s; }
    return [x0, y0];
  }

  /** |φ(k)|² in FFT order, normalized so that Σ |φ|² (dk)^d = norm of ψ (Parseval). */
  function momentumDensity(solver) {
    if (solver.shifted) solver._unshift();
    const { grid } = solver, N = grid.size;
    const fr = Float64Array.from(solver.re), fi = Float64Array.from(solver.im);
    grid.fft2(fr, fi, -1);
    // Parseval: Σ|F_k|² = N Σ|ψ_j|², so |φ_k|² = |F_k|² dV / (N dk^d) gives Σ|φ|² dk^d = Σ|ψ|² dV
    const dk = ((2 * Math.PI) / grid.lx) * (grid.ndim === 2 ? (2 * Math.PI) / grid.ly : 1);
    const scale = grid.dV / (N * dk);
    const out = new Float64Array(N);
    for (let i = 0; i < N; i++) out[i] = (fr[i] * fr[i] + fi[i] * fi[i]) * scale;
    return out;
  }

  /**
   * Draw one detection position from a detector profile (e.g. FluxDetector.profile): y is chosen with
   * probability ∝ max(profile, 0) among lo <= y < hi and spread uniformly within its cell. null if empty.
   */
  function sampleProfile(profile, ys, random = Math.random, lo = -Infinity, hi = Infinity) {
    const dy = ys.length > 1 ? ys[1] - ys[0] : 1;
    let total = 0;
    for (let j = 0; j < profile.length; j++) if (ys[j] >= lo && ys[j] < hi && profile[j] > 0) total += profile[j];
    if (!(total > 0)) return null;
    let u = random() * total;
    for (let j = 0; j < profile.length; j++) {
      if (!(ys[j] >= lo && ys[j] < hi && profile[j] > 0)) continue;
      u -= profile[j];
      if (u <= 0) return ys[j] + (random() - 0.5) * dy;
    }
    return ys[ys.length - 1];
  }

  root.QWave = { FFT, Grid, Solver, FluxDetector, gaussianPacket, absorptionRate, centeredCells, groundStateSteps,
                 measurePosition, momentumDensity, sampleProfile };
})(globalThis);
