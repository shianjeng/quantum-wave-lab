"""
Observables that can be measured while the simulation runs.

Probability current (ħ = m = 1):   j = Im(ψ* ∇ψ),   ∂|ψ|²/∂t + ∇·j = 0   (without absorber)

The probability that has crossed a line x = x0 up to time t is ∫ Φ dt with the flux
Φ(t) = ∫ j_x(x0, y, t) dy. Unlike summing |ψ|² behind a wall at the end, this

  * gives the transmitted probability live, while the wave is still scattering,
  * is not affected by the absorber eating the transmitted wave later on,
  * can be split into segments (e.g. one detector per slit).
"""
from __future__ import annotations

import numpy as np

from .solver import Grid

# 4th-order central difference: f'(x) ≈ [f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h)] / (12 h)
_FD4 = {-2: 1 / 12, -1: -8 / 12, 1: 8 / 12, 2: -1 / 12}


def gradient(grid: Grid, psi, axis=0):
    """∂ψ/∂x_axis on the whole grid, computed spectrally (exact for band-limited ψ)."""
    return np.fft.ifft(1j * grid.kaxes[axis].reshape([-1 if a == axis else 1 for a in range(grid.ndim)])
                       * np.fft.fft(psi, axis=axis), axis=axis)


def probability_current(grid: Grid, psi, axis=0):
    """j_axis = Im(ψ* ∂ψ/∂x_axis) on the whole grid."""
    return np.imag(np.conj(psi) * gradient(grid, psi, axis))


class FluxDetector:
    """Counts the probability that crosses the plane x_axis = position (in the +x_axis direction; flow
    the other way counts negative).

    In 2D, `span=(lo, hi)` restricts the detector to lo <= y < hi on the other axis (e.g. one slit).
    method="spectral" differentiates with FFTs along `axis` (exact, one extra FFT per sample);
    method="fd4" uses a 4th-order finite difference that only needs 5 grid lines around the detector —
    the cheap variant for the browser port. Its relative error is ≈ (k dx)⁴/30 for wavenumber k
    (5e-5 at k dx = 0.2; validate.py compares both).
    The detector sits on the grid line nearest to `position` (see `self.position`).

    The flux is sampled once per time step, so it needs a well-resolved solution: with sharp potential
    edges (rect_barrier, hard walls) the splitting error puts a little probability into high-k modes
    whose interference with the main wave makes the flux oscillate faster than the time step, and the
    time integral aliases (rect_barrier: ~20 % error at dt = 0.02, still ~6 % at dt = 0.005, depending
    on where the detector sits). Use smooth potentials
    (sech2_barrier, rect_barrier(edge=...)) or a smaller dt; summing |ψ|² behind the barrier is immune.

    Usage:
        det = FluxDetector(grid, position=2.0)
        det.record(solver.t, psi)                              # initial sample
        psi = solver.step(psi, n, callback=det.record)         # accumulates ∫Φ dt (trapezoid rule)
        det.transmitted                                        # probability that crossed so far
    """

    def __init__(self, grid: Grid, position, axis=0, span=None, method="spectral"):
        if method not in ("spectral", "fd4"):
            raise ValueError(f"method must be 'spectral' or 'fd4', got {method!r}")
        if span is not None and grid.ndim != 2:
            raise ValueError("span is only meaningful on a 2D grid")
        self.grid, self.axis, self.method = grid, axis, method
        ax = grid.axes[axis]
        self.index = int(np.argmin(np.abs(ax - position)))
        self.position = float(ax[self.index])
        self.weight = np.prod([d for a, d in enumerate(grid.dx) if a != axis])   # line element dy
        self.mask = None
        if span is not None:
            other = grid.axes[1 - axis]
            self.mask = (other >= span[0]) & (other < span[1])
        self.reset()

    def reset(self):
        self.times, self.fluxes = [], []
        self.transmitted = 0.0

    def _line(self, arr, offset=0):
        idx = (self.index + offset) % self.grid.n[self.axis]
        return np.take(arr, idx, axis=self.axis)

    def flux(self, psi):
        """Φ = ∫ j_axis dy through the detector at this instant."""
        if self.method == "spectral":
            dpsi = self._line(gradient(self.grid, psi, self.axis))
        else:
            dpsi = sum(c * self._line(psi, o) for o, c in _FD4.items()) / self.grid.dx[self.axis]
        j = np.imag(np.conj(self._line(psi)) * dpsi)
        if self.mask is not None:
            j = j[self.mask]
        return float(np.sum(j) * self.weight)

    def record(self, t, psi):
        """Add a sample; integrates Φ over time with the trapezoid rule (O(dt²), like the solver)."""
        phi = self.flux(psi)
        if self.times:
            self.transmitted += 0.5 * (phi + self.fluxes[-1]) * (t - self.times[-1])
        self.times.append(float(t))
        self.fluxes.append(phi)
        return phi
