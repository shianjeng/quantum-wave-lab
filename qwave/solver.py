"""
Split-step Fourier (Strang splitting) solver for the time-dependent Schrödinger equation.

    i dψ/dt = [ -1/2 ∇² + V(r, t) + g |ψ|² ] ψ        (units: ħ = m = 1)

g = 0 is the linear Schrödinger equation; g ≠ 0 the Gross–Pitaevskii / nonlinear Schrödinger
equation (g > 0 repulsive, g < 0 attractive; ψ normalized to 1, so g includes the atom number).

One step (exactly unitary without absorber; O(dt²) for smooth V):
    ψ ← e^{-i (V + g|ψ|²) dt/2} ψ
    ψ ← IFFT[ e^{-i k² dt/2} FFT[ψ] ]
    ψ ← e^{-i (V + g|ψ|²) dt/2} ψ
The nonlinear half-steps are exact (they do not change |ψ|), so the scheme stays unitary and O(dt²).
But with g ≠ 0 it is only stable for dt · k_max²/2 < π (k_max² summed over the axes; Weideman &
Herbst 1986): beyond that, grid modes whose kinetic phase per step is a multiple of 2π are pumped by
the nonlinearity and blow up after a while (validate.py). Solver warns, diagnose() reports it.

For a static potential the closing half-step of one step and the opening half-step of the
next are merged into a single full-step factor e^{-i V dt}, so n steps cost n + 1 potential
multiplications instead of 2n. The absorber is diagonal in position space like V, so it is
folded into the same factor.

Works in 1D or 2D (any dimension, via numpy.fft.fftn). This file is the
reference implementation that the C++/WASM port is tested against.
"""
from __future__ import annotations

import warnings

import numpy as np


class Grid:
    """Uniform periodic grid. `n` and `length` are per-axis tuples (or scalars for 1D)."""

    def __init__(self, n, length):
        self.n = tuple(int(v) for v in np.atleast_1d(n))
        self.length = tuple(float(v) for v in np.atleast_1d(length))
        if len(self.n) != len(self.length):
            raise ValueError(f"n has {len(self.n)} axes but length has {len(self.length)}")
        if min(self.n) < 2 or min(self.length) <= 0:
            raise ValueError(f"need n >= 2 and length > 0 on every axis, got n={self.n}, length={self.length}")
        self.ndim = len(self.n)
        self.dx = tuple(L / N for L, N in zip(self.length, self.n))
        self.axes = [(-L / 2 + np.arange(N) * d) for L, N, d in zip(self.length, self.n, self.dx)]
        self.kaxes = [2 * np.pi * np.fft.fftfreq(N, d) for N, d in zip(self.n, self.dx)]
        self.kmax = tuple(np.pi / d for d in self.dx)                 # Nyquist wavenumber per axis
        self.x = np.meshgrid(*self.axes, indexing="ij")          # list of coordinate arrays
        self.k = np.meshgrid(*self.kaxes, indexing="ij")
        self.k2 = sum(ki**2 for ki in self.k)
        self.dV = float(np.prod(self.dx))                         # volume element

    # --- helpers -----------------------------------------------------------
    def norm(self, psi):
        return float(np.sum(np.abs(psi) ** 2) * self.dV)

    def normalize(self, psi):
        return psi / np.sqrt(self.norm(psi))

    # expectation values are divided by the norm, so they stay correct
    # when the absorber has removed part of the wave function
    def expect_x(self, psi, axis=0):
        p = np.abs(psi) ** 2
        return float(np.sum(self.x[axis] * p) / np.sum(p))

    def std_x(self, psi, axis=0):
        p = np.abs(psi) ** 2
        p = p / np.sum(p)
        m = np.sum(self.x[axis] * p)
        return float(np.sqrt(np.sum((self.x[axis] - m) ** 2 * p)))

    def expect_k(self, psi, axis=0):
        p = np.abs(np.fft.fftn(psi)) ** 2
        return float(np.sum(self.k[axis] * p) / np.sum(p))


def gaussian_packet(grid: Grid, center, sigma, k0):
    """ψ ∝ exp(-(r-r0)²/(4σ²) + i k0·r).  σ is the std-dev of |ψ|²."""
    center, sigma, k0 = (np.broadcast_to(np.asarray(a, float), (grid.ndim,)) for a in (center, sigma, k0))
    arg = np.zeros(grid.n, dtype=complex)
    for xi, c, s, k in zip(grid.x, center, sigma, k0):
        arg += -((xi - c) ** 2) / (4 * s**2) + 1j * k * xi
    return grid.normalize(np.exp(arg))


def absorbing_mask(grid: Grid, dt, width_frac=0.08, gamma_max=2.0, width=None):
    """Per-step mask exp(-γ(r) dt) with a smooth sin² ramp γ(r) near the edges
    (equivalent to a complex absorbing potential -iγ). Because it depends on dt,
    the absorption per unit time is independent of the step size.
    Kills waves that would otherwise wrap around the periodic box.

    The layer is `width` long (absolute) or, if width is None, `width_frac` of each box length.
    How well it absorbs depends on the wavenumber k, the layer width and gamma_max: slow waves are
    reflected by a steep ramp, fast waves cross a thin/weak layer and wrap around. Measured residual
    (reflected + wrapped) probability, see `absorber_residual` and validate.py:

        layer width   gamma_max   k = 0.5   k = 1     k = 2     k = 4     k = 8
        3.2           2           3e-1      6e-2      3e-3      4e-2      2e-1    (8 % of a 40-wide box)
        3.2           10          5e-1      2e-1      2e-2      3e-5      4e-4
        32            2           6e-3      2e-6      6e-8      7e-8      4e-7    (8 % of a 400-wide box)

    Rule of thumb: the probability that crosses both edge layers (and wraps around) is
    ≈ exp(-2 · gamma_max · width / k), so pick gamma_max ≳ 5 k / width; much larger values make
    the ramp steep and reflect the slow waves instead."""
    return np.exp(-absorption_rate(grid, width_frac, gamma_max, width) * dt)


def absorption_rate(grid: Grid, width_frac=0.08, gamma_max=2.0, width=None):
    """γ(r): 0 in the interior, rising as sin² to gamma_max at the box edges."""
    gamma = np.zeros(grid.n)
    for xi, L in zip(grid.x, grid.length):
        w = width_frac * L if width is None else width
        d = L / 2 - np.abs(xi)                     # distance to nearest edge
        edge = np.clip((w - d) / w, 0, 1)          # 0 inside, →1 at the boundary
        gamma = np.maximum(gamma, gamma_max * np.sin(np.pi / 2 * edge) ** 2)
    return gamma


def absorber_residual(k, width, gamma_max=2.0, dt=0.01):
    """Probability left after a nearly monochromatic 1D wave packet (wavenumber k, σ = 15) has run
    into an absorbing layer of the given absolute width: reflected + transmitted-and-wrapped part.
    Independent of the box size, so it can be used to choose absorber parameters for any grid."""
    L = 400.0
    g = Grid(4096, L)
    s = Solver(g, np.zeros(g.n), dt, absorber=dict(width=width, gamma_max=gamma_max))
    psi = s.step(gaussian_packet(g, 0.0, 15.0, k), int(round((L / 2 + 60) / k / dt)))
    return g.norm(psi)


class Solver:
    """Propagates ψ on `grid` under the potential `V` with time step `dt`.

    V         array on the grid, or a callable V(t) -> array for a time-dependent potential
              (evaluated at the midpoint t + dt/2 of every step, which keeps the method O(dt²)).
    absorber  None/False (periodic box), True (absorbing layer with default parameters) or a dict
              of keyword arguments for `absorbing_mask`. The mask is built here, for this dt, so it
              can never silently belong to a different step size.
    dtype     np.complex128 (default) or np.complex64 to emulate the single precision of a
              WASM/WebGL port: ψ and all phase factors are stored in that precision between
              operations. (NumPy >= 2 also runs the FFT in single precision; NumPy 1.x upcasts it.)
    renormalize_every
              None (default) or N: rescale ψ to unit norm every N steps. In single precision the norm
              drifts by ~1e-4 per 1000 steps (validate.py), which adds up in a long interactive run.
              Only allowed without absorber, where the norm is supposed to decrease.
    nonlinearity
              g of the Gross–Pitaevskii term g|ψ|² (default 0: linear Schrödinger equation).
    """

    def __init__(self, grid: Grid, V, dt, absorber=None, dtype=np.complex128, renormalize_every=None,
                 nonlinearity=0.0):
        self.grid = grid
        self.nonlinearity = float(nonlinearity)
        self.dtype = np.dtype(dtype)
        if self.dtype.kind != "c":
            raise ValueError(f"dtype must be complex, got {self.dtype}")
        if absorber is True:
            absorber = {}
        elif absorber is False:
            absorber = None
        elif absorber is not None and not isinstance(absorber, dict):
            raise TypeError("absorber must be None, True/False or a dict of absorbing_mask() "
                            "arguments (a precomputed mask would silently tie it to another dt)")
        self._absorber_args = absorber
        if renormalize_every is not None:
            if absorber is not None:
                raise ValueError("renormalize_every cannot be combined with an absorber "
                                 "(the absorbed probability must stay lost)")
            if int(renormalize_every) < 1:
                raise ValueError(f"renormalize_every must be >= 1, got {renormalize_every}")
            renormalize_every = int(renormalize_every)
        self.renormalize_every = renormalize_every
        self.nsteps_done = 0
        self.t = 0.0
        self._V_input = V
        self.set_dt(dt)

    # --- configuration -------------------------------------------------------
    def set_dt(self, dt):
        """Change the time step; rebuilds every dt-dependent factor (kinetic, potential, absorber)."""
        if not dt > 0:
            raise ValueError(f"dt must be positive, got {dt}")
        self.dt = float(dt)
        if self.nonlinearity and self.nonlinear_stability_number() > np.pi:
            warnings.warn(self._instability_message(), RuntimeWarning, stacklevel=2)
        self.kin_phase = np.exp(-0.5j * self.grid.k2 * self.dt).astype(self.dtype)
        self.absorber = (None if self._absorber_args is None
                         else absorbing_mask(self.grid, self.dt, **self._absorber_args))
        self.set_potential(self._V_input)

    def nonlinear_stability_number(self):
        """dt · k_max²/2: the largest kinetic phase per step. Must stay below π when g ≠ 0."""
        return self.dt * sum(k**2 for k in self.grid.kmax) / 2

    def _instability_message(self):
        dt_max = 2 * np.pi / sum(k**2 for k in self.grid.kmax)
        return (f"nonlinear split-step instability: dt·k_max²/2 = {self.nonlinear_stability_number():.3g} > π "
                f"with g = {self.nonlinearity:g}; use dt < {dt_max:.3g} or a coarser grid")

    def set_potential(self, V):
        """Can be called at any time (e.g. when the user draws a new barrier).
        Accepts an array or a callable V(t) -> array."""
        self._V_input = V
        self.time_dependent = callable(V)
        self.V = self._eval_V(self.t)
        self._build_potential_factors(self.V)

    def _eval_V(self, t):
        V = self._V_input(t) if self.time_dependent else self._V_input
        V = np.asarray(V, dtype=float)
        if V.shape != self.grid.n:
            raise ValueError(f"potential has shape {V.shape}, grid is {self.grid.n}")
        return V

    def _build_potential_factors(self, V):
        half = np.exp(-0.5j * V * self.dt)
        full = half * half
        if self.absorber is not None:
            half_abs = half * self.absorber
            full_abs = full * self.absorber
        else:
            half_abs, full_abs = half, full
        self.half_pot_phase = half.astype(self.dtype)
        self._half_abs = half_abs.astype(self.dtype)          # closing half-step (+ absorber)
        self._full_abs = full_abs.astype(self.dtype)          # merged half-steps (+ absorber)

    # --- propagation -------------------------------------------------------
    def _renormalize(self, psi):
        norm = np.sum(np.abs(psi) ** 2, dtype=np.float64) * self.grid.dV
        return psi * psi.real.dtype.type(1 / np.sqrt(norm))

    def _renormalize_due(self):
        return self.renormalize_every is not None and self.nsteps_done % self.renormalize_every == 0

    def _kinetic(self, psi):
        return np.fft.ifftn(self.kin_phase * np.fft.fftn(psi)).astype(self.dtype, copy=False)

    def _potential_half(self, psi):
        """e^{-i (V + g|ψ|²) dt/2} ψ for the nonlinear equation (exact: |ψ| is unchanged)."""
        v = self.V + self.nonlinearity * np.abs(psi) ** 2
        return np.exp(-0.5j * self.dt * v).astype(self.dtype, copy=False) * psi

    def _finish_step(self, psi, callback, every, observed=None):
        self.t += self.dt
        self.nsteps_done += 1
        if self._renormalize_due():
            psi = self._renormalize(psi)
        if callback is not None and self.nsteps_done % every == 0:
            callback(self.t, psi if observed is None else observed(psi))
        return psi

    def step(self, psi, nsteps=1, callback=None, every=1):
        """Advance ψ by `nsteps` steps and return it.

        callback(t, ψ), if given, is called after every `every`-th step (counted over the solver's
        lifetime) with the wave function at that time — for detectors, recording or live plots.
        It must not modify the array it receives.
        """
        psi = np.asarray(psi).astype(self.dtype, copy=False)
        if nsteps <= 0:
            return psi
        if every < 1:
            raise ValueError(f"every must be >= 1, got {every}")
        if self.time_dependent or self.nonlinearity:
            for _ in range(nsteps):
                if self.time_dependent:
                    self.V = self._eval_V(self.t + 0.5 * self.dt)
                    self._build_potential_factors(self.V)
                if self.nonlinearity:
                    psi = self._potential_half(self._kinetic(self._potential_half(psi)))
                    if self.absorber is not None:
                        psi = psi * self.absorber
                else:
                    psi = self._half_abs * self._kinetic(self.half_pot_phase * psi)
                psi = self._finish_step(psi, callback, every)
            if self.time_dependent:
                self.V = self._eval_V(self.t)
            return psi
        # static linear V: e^{-iVdt/2} K e^{-iVdt/2} · e^{-iVdt/2} K e^{-iVdt/2} ... with adjacent halves
        # merged. Between steps the running array is e^{-iVdt/2} ψ(t), so an observer gets it multiplied
        # back by e^{+iVdt/2}; a renormalization is a scalar factor and can be applied at any point.
        psi = self.half_pot_phase * psi
        undo_half = np.conj(self.half_pot_phase) if callback is not None else None
        for i in range(nsteps):
            psi = self._kinetic(psi)
            if i < nsteps - 1:
                psi = self._full_abs * psi
                psi = self._finish_step(psi, callback, every, observed=lambda p: undo_half * p)
            else:
                psi = self._half_abs * psi
                psi = self._finish_step(psi, callback, every)
        return psi

    # --- diagnostics -------------------------------------------------------
    def energy(self, psi):
        """Energy per particle at the current time, E = ⟨T⟩ + ⟨V⟩ + (g/2)∫|ψ|⁴ for normalized ψ
        (the kinetic part evaluated in k-space). Conserved by the exact dynamics."""
        return self._energy_terms(psi, interaction_weight=0.5)

    def chemical_potential(self, psi):
        """μ = ⟨T⟩ + ⟨V⟩ + g∫|ψ|⁴: the eigenvalue of a stationary state, ψ(t) = e^{-iμt} ψ(0).
        Equals energy() for g = 0."""
        return self._energy_terms(psi, interaction_weight=1.0)

    def _energy_terms(self, psi, interaction_weight):
        g = self.grid
        psi = np.asarray(psi, dtype=complex)
        rho = np.abs(psi) ** 2
        phik = np.fft.fftn(psi)
        # Parseval: Σ|ψ|² dV = Σ|φ|² dV / N_total
        ekin = 0.5 * np.sum(g.k2 * np.abs(phik) ** 2) * g.dV / psi.size
        epot = np.sum(self.V * rho) * g.dV
        eint = interaction_weight * self.nonlinearity * np.sum(rho**2) * g.dV
        return float((ekin + epot + eint) / g.norm(psi))

    def diagnose(self, psi, tol=1e-6, max_loss_rate=1e-4):
        """Sanity checks for the current ψ; returns a list of human-readable problems (empty = OK).

        * max |V + g|ψ|²|·dt where ψ actually is: above ~1 the splitting error grows quickly at hard walls
          (make_gif.py: V = 200 with dt = 0.01 gives ~17 % density error, dt = 0.005 gives 0.3 %).
        * momentum content near the Nyquist limit k_max = π/dx: if more than `tol` of |φ(k)|² lies
          beyond 2/3 k_max on any axis, the grid is too coarse for this ψ (aliasing) — or ψ does not
          vanish at the periodic box edge, and the jump there shows up as high-k content.
        * with g ≠ 0: the stability limit dt · k_max²/2 < π of the nonlinear split-step method.
        * at t = 0 only: how fast the absorbing layer is already eating the initial state
          (more than `max_loss_rate` of the probability per unit time).
        Call it on the initial state and, for interactive use, every few frames.
        """
        g = self.grid
        issues = []
        rho = np.abs(psi) ** 2
        present = rho > tol * rho.max()
        v_eff = self.V + self.nonlinearity * rho
        vdt = float(np.max(np.abs(v_eff[present]))) * self.dt if present.any() else 0.0
        if vdt > 1 + 1e-9:
            issues.append(f"max|V + g|ψ|²|·dt = {vdt:.3g} > 1 where ψ is non-negligible: "
                          f"reduce dt below {self.dt / vdt:.3g} or lower the walls / interaction")
        if self.nonlinearity and self.nonlinear_stability_number() > np.pi:
            issues.append(self._instability_message())
        pk = np.abs(np.fft.fftn(psi)) ** 2
        pk /= pk.sum()
        for axis, (k, kmax) in enumerate(zip(g.k, g.kmax)):
            frac = float(pk[np.abs(k) > 2 / 3 * kmax].sum())
            if frac > tol:
                issues.append(f"{frac:.1e} of |φ(k)|² lies beyond 2/3 of k_max = {kmax:.3g} on axis {axis}: "
                              f"either ψ has structure finer than the grid (refine to dx < {g.dx[axis] / 2:.3g} "
                              f"or lower k0) or ψ is cut off at the box edge (move it inwards)")
        if self.absorber is not None and self.t == 0:
            gamma = -np.log(self.absorber) / self.dt
            rate = float(np.sum(2 * gamma * rho) / rho.sum())          # d(norm)/dt at t = 0
            if rate > max_loss_rate:
                issues.append(f"the initial state is already inside the absorbing layer (loses {rate:.1e} "
                              f"of its probability per unit time): start the packet further from the edge")
        return issues
