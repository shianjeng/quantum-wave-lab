"""Analytic reference results (ħ = m = 1) used by the validation scripts and tests."""
import numpy as np


def T_rect(E, V0, a):
    """Exact transmission coefficient of a rectangular barrier of height V0 and width a."""
    E = np.asarray(E, float)
    T = np.empty_like(E)
    lo, hi = E < V0, E > V0
    kap = np.sqrt(2 * (V0 - E[lo]))
    T[lo] = 1 / (1 + V0**2 * np.sinh(kap * a) ** 2 / (4 * E[lo] * (V0 - E[lo])))
    q = np.sqrt(2 * (E[hi] - V0))
    T[hi] = 1 / (1 + V0**2 * np.sin(q * a) ** 2 / (4 * E[hi] * (E[hi] - V0)))
    T[~(lo | hi)] = 1 / (1 + V0 * a**2 / 2)
    return T


def T_packet(k0, sigma, V0, a, n=20001):
    """T averaged over the momentum distribution of a Gaussian packet whose |ψ|² has std σ:
    |φ(k)|² ∝ exp(-(k-k0)²/(2σk²)), σk = 1/(2σ)."""
    sk = 1 / (2 * sigma)
    k = np.linspace(k0 - 8 * sk, k0 + 8 * sk, n)
    k = k[k > 0]                                   # negative-k weight is negligible here
    w = np.exp(-((k - k0) ** 2) / (2 * sk**2))
    return float(np.sum(w * T_rect(k**2 / 2, V0, a)) / np.sum(w))


def ho_levels(n, omega=1.0):
    """Lowest n energies of the 1D harmonic oscillator: E_j = ω (j + 1/2)."""
    return omega * (np.arange(n) + 0.5)


def driven_ho_x(t, x0, omega, F0, Omega):
    """⟨x⟩(t) of a harmonic oscillator driven by the force F0 cos(Ωt), starting at rest at x0.
    Exact for any quantum state (Ehrenfest is exact for quadratic potentials); needs Ω ≠ ω."""
    A = F0 / (omega**2 - Omega**2)
    t = np.asarray(t, float)
    return (x0 - A) * np.cos(omega * t) + A * np.cos(Omega * t)


def bloch_period(force, lattice_period):
    """Bloch-oscillation period T_B = 2π / (|F| a) of a particle in a lattice of period a under force F."""
    return 2 * np.pi / (abs(force) * lattice_period)


def fringe_period_ky(slit_sep):
    """Far-field double-slit fringe period in transverse momentum: Δk_y = 2π / d
    (equivalently d sin θ = n λ, since k_y = k sin θ)."""
    return 2 * np.pi / slit_sep
