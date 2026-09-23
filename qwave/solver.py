"""
Split-step Fourier (Strang splitting) solver for the time-dependent Schrödinger equation.

    i dψ/dt = [ -1/2 ∇² + V(r) ] ψ        (units: ħ = m = 1)

One step (exactly unitary without absorber; O(dt²) for smooth V):
    ψ ← e^{-i V dt/2} ψ
    ψ ← IFFT[ e^{-i k² dt/2} FFT[ψ] ]
    ψ ← e^{-i V dt/2} ψ

Works in 1D or 2D (any dimension, via numpy.fft.fftn). This file is the
reference implementation that the C++/WASM port is tested against.
"""
from __future__ import annotations

import numpy as np


class Grid:
    """Uniform periodic grid. `n` and `length` are per-axis tuples (or scalars for 1D)."""

    def __init__(self, n, length):
        self.n = tuple(np.atleast_1d(n).astype(int))
        self.length = tuple(np.atleast_1d(length).astype(float))
        assert len(self.n) == len(self.length)
        self.ndim = len(self.n)
        self.dx = tuple(L / N for L, N in zip(self.length, self.n))
        self.axes = [(-L / 2 + np.arange(N) * d) for L, N, d in zip(self.length, self.n, self.dx)]
        self.kaxes = [2 * np.pi * np.fft.fftfreq(N, d) for N, d in zip(self.n, self.dx)]
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


def gaussian_packet(grid: Grid, center, sigma, k0):
    """ψ ∝ exp(-(r-r0)²/(4σ²) + i k0·r).  σ is the std-dev of |ψ|²."""
    center, sigma, k0 = (np.broadcast_to(np.asarray(a, float), (grid.ndim,)) for a in (center, sigma, k0))
    arg = np.zeros(grid.n, dtype=complex)
    for xi, c, s, k in zip(grid.x, center, sigma, k0):
        arg += -((xi - c) ** 2) / (4 * s**2) + 1j * k * xi
    return grid.normalize(np.exp(arg))


def absorbing_mask(grid: Grid, dt, width_frac=0.08, gamma_max=2.0):
    """Per-step mask exp(-γ(r) dt) with a smooth sin² ramp γ(r) near the edges
    (equivalent to a complex absorbing potential -iγ). Because it depends on dt,
    the absorption per unit time is independent of the step size.
    Kills waves that would otherwise wrap around the periodic box."""
    gamma = np.zeros(grid.n)
    for xi, L in zip(grid.x, grid.length):
        w = width_frac * L
        d = L / 2 - np.abs(xi)                     # distance to nearest edge
        edge = np.clip((w - d) / w, 0, 1)          # 0 inside, →1 at the boundary
        gamma = np.maximum(gamma, gamma_max * np.sin(np.pi / 2 * edge) ** 2)
    return np.exp(-gamma * dt)


class Solver:
    def __init__(self, grid: Grid, V, dt, absorber=None):
        self.grid = grid
        self.dt = dt
        self.set_potential(V)
        self.kin_phase = np.exp(-0.5j * grid.k2 * dt)
        self.absorber = absorber
        self.t = 0.0

    def set_potential(self, V):
        """Can be called at any time (e.g. when the user draws a new barrier)."""
        self.V = np.asarray(V, dtype=float)
        self.half_pot_phase = np.exp(-0.5j * self.V * self.dt)

    def step(self, psi, nsteps=1):
        fftn, ifftn = np.fft.fftn, np.fft.ifftn
        for _ in range(nsteps):
            psi = self.half_pot_phase * psi
            psi = ifftn(self.kin_phase * fftn(psi))
            psi = self.half_pot_phase * psi
            if self.absorber is not None:
                psi = psi * self.absorber
            self.t += self.dt
        return psi

    def energy(self, psi):
        """⟨H⟩ = ⟨T⟩ + ⟨V⟩, kinetic part evaluated in k-space."""
        g = self.grid
        phik = np.fft.fftn(psi)
        # Parseval: Σ|ψ|² dV = Σ|φ|² dV / N_total
        ekin = 0.5 * np.sum(g.k2 * np.abs(phik) ** 2) * g.dV / psi.size
        epot = np.sum(self.V * np.abs(psi) ** 2) * g.dV
        return float((ekin + epot) / g.norm(psi))
