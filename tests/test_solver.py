"""Fast unit tests (run by CI). The full, figure-producing check is validate.py."""
import numpy as np

from qwave.potentials import harmonic, rect_barrier_1d
from qwave.solver import Grid, Solver, gaussian_packet


def test_norm_is_conserved():
    g = Grid((128, 128), (30, 30))
    s = Solver(g, harmonic(g), dt=0.01)
    psi = s.step(gaussian_packet(g, (2, 0), 1.0, (0, 3)), 200)
    assert abs(g.norm(psi) - 1) < 1e-12


def test_free_packet_width_matches_analytic():
    g = Grid(2048, 200)
    s = Solver(g, np.zeros(g.n), dt=0.02)
    psi = s.step(gaussian_packet(g, 0, 1.5, 0), 500)          # t = 10
    exact = 1.5 * np.sqrt(1 + (10 / (2 * 1.5**2)) ** 2)
    assert abs(g.std_x(psi) - exact) / exact < 1e-8


def test_ho_ground_state_is_stationary():
    g = Grid(512, 30)
    s = Solver(g, harmonic(g, 1.0), dt=0.005)
    psi0 = gaussian_packet(g, 0, 1 / np.sqrt(2), 0)
    psi = s.step(psi0, 1000)
    overlap = abs(np.sum(np.conj(psi0) * psi) * g.dV)
    assert overlap > 1 - 1e-8
    assert abs(s.energy(psi) - 0.5) < 1e-6


def test_barrier_width_is_exact_on_grid():
    g = Grid(1000, 50)                                          # dx = 0.05
    V = rect_barrier_1d(g, 1.0, 1.0)
    assert np.count_nonzero(V) == 20


def test_absorption_is_independent_of_dt():
    from qwave.solver import absorbing_mask
    g = Grid(1024, 100)
    left = []
    for dt in (0.01, 0.005):
        s = Solver(g, np.zeros(g.n), dt, absorber=absorbing_mask(g, dt))
        psi = s.step(gaussian_packet(g, 30, 2.0, 3.0), int(round(8 / dt)))
        left.append(g.norm(psi))
    assert left[0] < 0.5                                       # most of it was absorbed
    assert abs(left[0] - left[1]) < 1e-3


def test_expectation_values_are_normalized():
    g = Grid(512, 40)
    psi = 0.3 * gaussian_packet(g, 2.0, 1.0, 0)                # norm 0.09, e.g. after absorption
    assert abs(g.expect_x(psi) - 2.0) < 1e-10
    assert abs(g.std_x(psi) - 1.0) < 1e-10


def test_tunneling_matches_analytic():
    from qwave.analytic import T_packet
    g = Grid(2**14, 2**14 * 0.05)
    s = Solver(g, rect_barrier_1d(g, 1.5, 1.0), dt=0.02)
    psi = s.step(gaussian_packet(g, -150, 10.0, 1.5), int(300 / 1.5 / 0.02))
    T = np.sum(np.abs(psi[g.x[0] > 1]) ** 2) * g.dV
    assert abs(T - T_packet(1.5, 10.0, 1.5, 1.0)) < 5e-3


def test_check_on_grid_rejects_off_grid_lengths():
    import pytest
    from qwave.potentials import check_on_grid
    check_on_grid(0.125, 1.0, 1.25)
    with pytest.raises(ValueError):
        check_on_grid(0.15625, 1.0)
