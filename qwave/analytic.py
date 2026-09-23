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
