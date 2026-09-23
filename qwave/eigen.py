"""
Stationary states by imaginary-time propagation (t → -iτ):

    ψ ← e^{-V dτ/2} IFFT[ e^{-k² dτ/2} FFT[ψ] ] e^{-V dτ/2},   then renormalize.

Every component decays as e^{-E_n τ}, so the lowest state survives. Excited states are found one
after another by projecting out the states already found (Gram–Schmidt) after every step.

With a nonlinearity g (Gross–Pitaevskii), V is replaced by V + g|ψ|² evaluated on the current,
normalized ψ; only the ground state is available then, and the returned eigenvalue is the chemical
potential μ (not the energy per particle). Both nonlinear half-steps use the density at the start
of the step, which keeps the fixed point symmetric and its error O(dτ²); since μ, unlike the energy,
is not variational, that error shows up in μ directly (still < 1e-10 relative for dτ = 0.001–0.004
in validate.py, where the state is stationary in real time to ~1e-7).

The fixed point is an eigenvector of the split operator, which differs from the true eigenvector
by O(dτ²); the energy is evaluated as the exact ⟨H⟩ of that vector, so its error is only O(dτ⁴).
"""
from __future__ import annotations

import numpy as np

from .solver import Grid


def _energy(grid: Grid, V, psi, g=0.0):
    """⟨T⟩ + ⟨V⟩ + g∫|ψ|⁴: the energy for g = 0, the chemical potential μ otherwise."""
    rho = np.abs(psi) ** 2
    phik = np.fft.fftn(psi)
    ekin = 0.5 * np.sum(grid.k2 * np.abs(phik) ** 2) * grid.dV / psi.size
    epot = np.sum((V + g * rho) * rho) * grid.dV
    return float((ekin + epot) / grid.norm(psi))


def eigenstates(grid: Grid, V, n=1, dtau=0.01, tol=1e-12, max_steps=200_000, check_every=20, seed=0,
                nonlinearity=0.0):
    """Lowest `n` eigenstates of H = -∇²/2 + V (+ g|ψ|² with nonlinearity=g, then n must be 1).

    Returns (energies, states): energies is an array of length n (ascending), states a list of
    normalized arrays. Convergence: the energy changes by less than `tol` over `check_every` steps.
    The rate is set by the gap to the next state (error ∝ e^{-2 ΔE τ}), so nearly degenerate
    states converge slowly; raise max_steps for those. Raises RuntimeError if not converged.
    """
    V = np.asarray(V, dtype=float)
    g = float(nonlinearity)
    if g and n != 1:
        raise ValueError("with a nonlinearity only the ground state (n = 1) is defined")
    half_V = np.exp(-0.5 * V * dtau)
    kin = np.exp(-0.5 * grid.k2 * dtau)
    rng = np.random.default_rng(seed)
    # the ground state has no nodes, so it starts from a positive envelope (a random start would also
    # contain e.g. the antisymmetric partner of a double-well ground state, which decays very slowly);
    # the excited states start from random noise inside the envelope, which overlaps every state
    envelope = np.exp(-sum((xi / (0.25 * L)) ** 2 for xi, L in zip(grid.x, grid.length)))

    energies, states = [], []
    for i in range(n):
        if i == 0:
            psi = envelope.astype(complex)
        else:
            psi = envelope * (rng.standard_normal(grid.n) + 1j * rng.standard_normal(grid.n))
        E_old = np.inf
        for it in range(1, max_steps + 1):
            if g:
                # both halves with the density at the start of the step: at the fixed point the step is then
                # the symmetric Strang operator of H[ρ*], so the bias is O(dτ²) (not O(dτ))
                half = np.exp(-0.5 * dtau * (V + g * np.abs(psi) ** 2))
                psi = half * np.fft.ifftn(kin * np.fft.fftn(half * psi))
            else:
                psi = half_V * np.fft.ifftn(kin * np.fft.fftn(half_V * psi))
            for phi in states:
                psi -= np.sum(np.conj(phi) * psi) * grid.dV * phi
            psi = grid.normalize(psi)
            if it % check_every == 0:
                E = _energy(grid, V, psi, g)
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
    """Convenience wrapper: (E0, ψ0), or (μ, ψ0) with nonlinearity=g."""
    E, states = eigenstates(grid, V, 1, **kwargs)
    return float(E[0]), states[0]
