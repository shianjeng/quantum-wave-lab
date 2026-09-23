"""
Stationary states by imaginary-time propagation (t → -iτ):

    ψ ← e^{-V dτ/2} IFFT[ e^{-k² dτ/2} FFT[ψ] ] e^{-V dτ/2},   then renormalize.

Every component decays as e^{-E_n τ}, so the lowest state survives. Excited states are found one
after another by projecting out the states already found (Gram–Schmidt) after every step.

The fixed point is an eigenvector of the split operator, which differs from the true eigenvector
by O(dτ²); the energy is evaluated as the exact ⟨H⟩ of that vector, so its error is only O(dτ⁴).
"""
from __future__ import annotations

import numpy as np

from .solver import Grid


def _energy(grid: Grid, V, psi):
    phik = np.fft.fftn(psi)
    ekin = 0.5 * np.sum(grid.k2 * np.abs(phik) ** 2) * grid.dV / psi.size
    epot = np.sum(V * np.abs(psi) ** 2) * grid.dV
    return float((ekin + epot) / grid.norm(psi))


def eigenstates(grid: Grid, V, n=1, dtau=0.01, tol=1e-12, max_steps=200_000, check_every=20, seed=0):
    """Lowest `n` eigenstates of H = -∇²/2 + V.

    Returns (energies, states): energies is an array of length n (ascending), states a list of
    normalized arrays. Convergence: the energy changes by less than `tol` over `check_every` steps.
    The rate is set by the gap to the next state (error ∝ e^{-2 ΔE τ}), so nearly degenerate
    states converge slowly; raise max_steps for those. Raises RuntimeError if not converged.
    """
    V = np.asarray(V, dtype=float)
    half_V = np.exp(-0.5 * V * dtau)
    kin = np.exp(-0.5 * grid.k2 * dtau)
    rng = np.random.default_rng(seed)
    # random start inside a broad envelope: overlaps every state, none of them favoured by symmetry
    envelope = np.exp(-sum((xi / (0.25 * L)) ** 2 for xi, L in zip(grid.x, grid.length)))

    energies, states = [], []
    for _ in range(n):
        psi = envelope * (rng.standard_normal(grid.n) + 1j * rng.standard_normal(grid.n))
        E_old = np.inf
        for it in range(1, max_steps + 1):
            psi = half_V * np.fft.ifftn(kin * np.fft.fftn(half_V * psi))
            for phi in states:
                psi -= np.sum(np.conj(phi) * psi) * grid.dV * phi
            psi = grid.normalize(psi)
            if it % check_every == 0:
                E = _energy(grid, V, psi)
                if abs(E - E_old) < tol:
                    break
                E_old = E
        else:
            raise RuntimeError(f"state {len(states)} not converged after {max_steps} steps "
                               f"(last ΔE = {abs(E - E_old):.1e}); increase max_steps or dtau")
        # fix the global phase so real eigenfunctions come out real and positive at their maximum
        psi *= np.exp(-1j * np.angle(psi.flat[np.argmax(np.abs(psi))]))
        states.append(psi)
        energies.append(E)
    return np.array(energies), states


def ground_state(grid: Grid, V, **kwargs):
    """Convenience wrapper: (E0, ψ0)."""
    E, states = eigenstates(grid, V, 1, **kwargs)
    return float(E[0]), states[0]
