"""Quantum Wave Lab: split-step Fourier solver for the time-dependent Schrödinger equation."""
from .eigen import eigenstates, ground_state
from .observables import FluxDetector, probability_current
from .solver import Grid, Solver, absorber_residual, absorbing_mask, gaussian_packet

__version__ = "0.3.0"
__all__ = ["FluxDetector", "Grid", "Solver", "absorber_residual", "absorbing_mask", "eigenstates",
           "gaussian_packet", "ground_state", "probability_current"]
