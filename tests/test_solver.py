"""Fast unit tests (run by CI). The full, figure-producing check is validate.py."""
import numpy as np
import pytest

from qwave.analytic import T_packet, driven_ho_x, ho_levels
from qwave.eigen import eigenstates
from qwave.potentials import check_on_grid, double_slit, harmonic, paint_disk, rect_barrier, rect_barrier_1d, slits
from qwave.solver import Grid, Solver, absorber_residual, gaussian_packet


# ---------------------------------------------------------------- grid
def test_grid_rejects_mismatched_axes():
    with pytest.raises(ValueError):
        Grid((64, 64), 10.0)
    with pytest.raises(ValueError):
        Grid(64, -1.0)


def test_expectation_values_are_normalized():
    g = Grid(512, 40)
    psi = 0.3 * gaussian_packet(g, 2.0, 1.0, 1.5)              # norm 0.09, e.g. after absorption
    assert abs(g.expect_x(psi) - 2.0) < 1e-10
    assert abs(g.std_x(psi) - 1.0) < 1e-10
    assert abs(g.expect_k(psi) - 1.5) < 1e-10


# ---------------------------------------------------------------- propagation
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


def test_merged_half_steps_equal_single_steps():
    g = Grid(512, 30)
    V = harmonic(g) + rect_barrier(g, 3.0, 1.0)
    psi0 = gaussian_packet(g, -5, 1.0, 2.0)
    s1, s2 = Solver(g, V, 0.01, absorber=True), Solver(g, V, 0.01, absorber=True)
    many = s1.step(psi0, 300)
    one_by_one = psi0
    for _ in range(300):
        one_by_one = s2.step(one_by_one, 1)
    assert np.max(np.abs(many - one_by_one)) < 1e-12
    assert abs(s1.t - s2.t) < 1e-9


def test_second_order_in_dt_for_smooth_potential():
    g = Grid(256, 30)
    V = harmonic(g)
    psi0 = gaussian_packet(g, 2, 0.7, 1.0)
    ref = Solver(g, V, 1e-4).step(psi0, 10000)                 # t = 1
    err = [np.sqrt(g.norm(Solver(g, V, dt).step(psi0, round(1 / dt)) - ref)) for dt in (0.02, 0.01)]
    assert abs(np.log2(err[0] / err[1]) - 2) < 0.05


def test_time_dependent_potential_driven_oscillator():
    g = Grid(512, 40)
    x = g.x[0]
    drive = dict(omega=1.0, F0=0.5, Omega=1.7)
    s = Solver(g, lambda t: 0.5 * x**2 - drive["F0"] * np.cos(drive["Omega"] * t) * x, dt=0.005)
    psi = s.step(gaussian_packet(g, 1.0, 1 / np.sqrt(2), 0.0), 1000)   # t = 5
    assert abs(g.expect_x(psi) - driven_ho_x(s.t, 1.0, **drive)) < 1e-4
    assert abs(g.norm(psi) - 1) < 1e-12


def test_single_precision_stays_close_to_double():
    g = Grid((128, 128), (30, 30))
    V = paint_disk(np.zeros(g.n), g, (0, 0), 2.0, 100.0)
    psi0 = gaussian_packet(g, (-8, 0.5), 2.0, (3.0, 0))
    p64 = Solver(g, V, 0.005).step(psi0, 500)
    s32 = Solver(g, V, 0.005, dtype=np.complex64)
    p32 = s32.step(psi0, 500)
    assert p32.dtype == np.complex64
    assert np.sqrt(g.norm(p32 - p64)) < 1e-4
    assert abs(g.norm(p32) - 1) < 1e-4


def test_renormalization_removes_single_precision_drift():
    g = Grid((128, 128), (30, 30))
    V = harmonic(g)
    psi0 = gaussian_packet(g, (2, 0), 1.0, (0, 3))
    s = Solver(g, V, 0.01, dtype=np.complex64, renormalize_every=100)
    psi = s.step(psi0, 250)
    psi = s.step(psi, 250)                                     # counting continues across calls
    assert s.nsteps_done == 500
    assert abs(g.norm(psi) - 1) < 1e-6
    # in double precision a renormalization is a no-op: same result as without it
    exact = Solver(g, V, 0.01).step(psi0, 500)
    assert np.max(np.abs(Solver(g, V, 0.01, renormalize_every=100).step(psi0, 500) - exact)) < 1e-12
    with pytest.raises(ValueError):
        Solver(g, V, 0.01, absorber=True, renormalize_every=100)


def test_tunneling_matches_analytic():
    g = Grid(2**14, 2**14 * 0.05)
    s = Solver(g, rect_barrier(g, 1.5, 1.0), dt=0.02)
    psi = s.step(gaussian_packet(g, -150, 10.0, 1.5), int(300 / 1.5 / 0.02))
    T = np.sum(np.abs(psi[g.x[0] > 1]) ** 2) * g.dV
    assert abs(T - T_packet(1.5, 10.0, 1.5, 1.0)) < 5e-3


# ---------------------------------------------------------------- absorber
def test_absorption_is_independent_of_dt():
    g = Grid(1024, 100)
    left = []
    for dt in (0.01, 0.005):
        s = Solver(g, np.zeros(g.n), dt, absorber=True)
        psi = s.step(gaussian_packet(g, 30, 2.0, 3.0), int(round(8 / dt)))
        left.append(g.norm(psi))
    assert left[0] < 0.5                                       # most of it was absorbed
    assert abs(left[0] - left[1]) < 1e-3


def test_set_dt_rebuilds_absorber():
    g = Grid(1024, 100)
    s = Solver(g, np.zeros(g.n), 0.01, absorber=dict(gamma_max=4.0))
    s.set_dt(0.005)
    fresh = Solver(g, np.zeros(g.n), 0.005, absorber=dict(gamma_max=4.0))
    assert np.allclose(s.absorber, fresh.absorber)
    assert np.allclose(s.kin_phase, fresh.kin_phase)


def test_precomputed_absorber_array_is_rejected():
    from qwave.solver import absorbing_mask
    g = Grid(256, 20)
    with pytest.raises(TypeError):
        Solver(g, np.zeros(g.n), 0.01, absorber=absorbing_mask(g, 0.02))


def test_wide_absorber_leaves_little_behind():
    assert absorber_residual(2.0, width=32.0) < 1e-6


# ---------------------------------------------------------------- potentials
def test_barrier_width_is_exact_on_grid():
    g = Grid(1000, 50)                                          # dx = 0.05
    assert np.count_nonzero(rect_barrier(g, 1.0, 1.0)) == 20
    assert np.array_equal(rect_barrier_1d(g, 1.0, 1.0), rect_barrier(g, 1.0, 1.0))


def test_rect_barrier_in_2d_is_a_wall_along_y():
    g = Grid((80, 40), (10.0, 5.0))                             # dx = 0.125
    V = rect_barrier(g, 5.0, 1.0)
    assert np.count_nonzero(V) == 8 * 40
    assert np.all(V == V[:, :1])


def test_smooth_barrier_has_same_area():
    g = Grid(4000, 100)
    sharp, smooth = rect_barrier(g, 2.0, 1.0), rect_barrier(g, 2.0, 1.0, edge=0.2)
    assert abs(np.sum(sharp) - np.sum(smooth)) * g.dx[0] < 1e-9
    assert smooth.max() < 2.0


def test_slits_generalizes_double_slit():
    g = Grid((160, 160), (20, 20))
    assert np.array_equal(slits(g, (2.5, -2.5)), double_slit(g))
    single = slits(g, (0.0,), slit_width=1.25)
    assert np.count_nonzero(single == 0) - np.count_nonzero(slits(g, ()) == 0) == 4 * 10  # 4 x 10 cells open


def test_check_on_grid_rejects_off_grid_lengths():
    check_on_grid(0.125, 1.0, 1.25)
    with pytest.raises(ValueError):
        check_on_grid(0.15625, 1.0)


def test_paint_disk_set_and_add():
    g = Grid((64, 64), (16, 16))
    V = paint_disk(np.zeros(g.n), g, (0, 0), 2.0, 5.0)
    both = paint_disk(V, g, (1, 0), 2.0, 5.0, mode="add")
    over = paint_disk(V, g, (1, 0), 2.0, 5.0)
    assert both.max() == 10.0 and over.max() == 5.0
    assert V.max() == 5.0                                       # input untouched
    erased = paint_disk(V, g, (0, 0), 3.0, 0.0)
    assert not erased.any()
    with pytest.raises(ValueError):
        paint_disk(V, g, (0, 0), 1.0, 1.0, mode="xor")


# ---------------------------------------------------------------- eigenstates
def test_imaginary_time_gives_ho_levels():
    g = Grid(256, 20)
    E, states = eigenstates(g, harmonic(g), 3, dtau=0.01)
    assert np.max(np.abs(E - ho_levels(3))) < 1e-6
    overlaps = [[np.sum(np.conj(a) * b) * g.dV for b in states] for a in states]
    assert np.allclose(overlaps, np.eye(3), atol=1e-8)


# ---------------------------------------------------------------- diagnostics
def test_diagnose_clean_setup():
    g = Grid((256, 256), (40, 40))
    s = Solver(g, double_slit(g, height=50), 0.01, absorber=True)
    assert s.diagnose(gaussian_packet(g, (-10, 0), 1.5, (4, 0))) == []


def test_diagnose_flags_problems():
    g = Grid((160, 160), (20, 20))                              # dx = 0.125, k_max ≈ 25
    wall = paint_disk(np.zeros(g.n), g, (0, 0), 2.0, 200.0)
    # packet sitting on a hard wall with a too-large dt
    msgs = Solver(g, wall, 0.01).diagnose(gaussian_packet(g, (-2, 0), 1.0, (0, 0)))
    assert any("V|·dt" in m for m in msgs)
    # momentum close to the Nyquist limit
    msgs = Solver(g, np.zeros(g.n), 0.005).diagnose(gaussian_packet(g, (-5, 0), 1.0, (20, 0)))
    assert any("k_max" in m for m in msgs)
    # initial packet overlapping the absorbing layer
    msgs = Solver(g, np.zeros(g.n), 0.005, absorber=True).diagnose(gaussian_packet(g, (-8.5, 0), 1.0, (2, 0)))
    assert any("absorbing layer" in m for m in msgs)
