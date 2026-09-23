"""Quantum Wave Lab: split-step Fourier solver for the time-dependent Schrödinger equation."""
from .eigen import eigenstates, ground_state
from .solver import Grid, Solver, absorber_residual, absorbing_mask, gaussian_packet

__version__ = "0.2.0"
__all__ = ["Grid", "Solver", "absorber_residual", "absorbing_mask", "eigenstates", "gaussian_packet",
           "ground_state"]
