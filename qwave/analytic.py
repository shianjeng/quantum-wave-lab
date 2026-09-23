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


def T_packet(k0, sigma, V0, a, n=20001, T=None):
    """T averaged over the momentum distribution of a Gaussian packet whose |ψ|² has std σ:
    |φ(k)|² ∝ exp(-(k-k0)²/(2σk²)), σk = 1/(2σ). T(E, V0, a) defaults to T_rect."""
    T = T_rect if T is None else T
    sk = 1 / (2 * sigma)
    k = np.linspace(k0 - 8 * sk, k0 + 8 * sk, n)
    k = k[k > 0]                                   # negative-k weight is negligible here
    w = np.exp(-((k - k0) ** 2) / (2 * sk**2))
    return float(np.sum(w * T(k**2 / 2, V0, a)) / np.sum(w))


def T_sech2(E, V0, w):
    """Exact transmission through the smooth barrier V0 / cosh²(x / w) (Landau & Lifshitz §25),
    valid for 8 V0 w² > 1:  T = sinh²(πkw) / [sinh²(πkw) + cosh²((π/2)√(8V0w² − 1))]."""
    k = np.sqrt(2 * np.asarray(E, float))
    s = np.sinh(np.pi * k * w) ** 2
    return s / (s + np.cosh(0.5 * np.pi * np.sqrt(8 * V0 * w**2 - 1)) ** 2)


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


def bright_soliton(x, t, g, v=0.0, x0=0.0):
    """Exact moving bright soliton of i ψ_t = -ψ_xx/2 + g|ψ|²ψ with g < 0, normalized to 1:
    ψ = √(η/2) sech(η(x - x0 - vt)) exp(i[v(x - x0) + (η² - v²)t/2]),  η = |g|/2."""
    if g >= 0:
        raise ValueError("bright solitons need an attractive nonlinearity g < 0")
    eta = abs(g) / 2
    u = eta * (np.asarray(x, float) - x0 - v * t)
    return np.sqrt(eta / 2) / np.cosh(u) * np.exp(1j * (v * (np.asarray(x, float) - x0) + 0.5 * (eta**2 - v**2) * t))


def bright_soliton_energy(g, v=0.0):
    """Energy per particle of the bright soliton: E = v²/2 - g²/24."""
    return v**2 / 2 - g**2 / 24


def thomas_fermi_mu_1d(g, omega=1.0):
    """Chemical potential of a 1D condensate in V = ω²x²/2 in the Thomas–Fermi limit (kinetic energy
    neglected, n(x) = (μ - V)/g):  μ = (3 g ω / (4√2))^(2/3). Accurate for μ ≫ ω."""
    return (3 * g * omega / (4 * np.sqrt(2))) ** (2 / 3)
