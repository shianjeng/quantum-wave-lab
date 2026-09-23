"""
Numerical validation against analytic results (ħ = m = 1).

  1. Norm conservation (2D, no absorber)           -> should stay 1 to machine precision
  2. Free Gaussian spreading:  σ(t) = σ0 √(1 + (t / 2σ0²)²)
  3. Harmonic oscillator coherent state: ⟨x⟩(t) = x0 cos(ωt), energy constant
  4. Tunneling through a rectangular barrier: transmitted probability vs
     analytic T(E) averaged over the packet's momentum distribution
  5. Double slit: far-field fringe period Δk_y = 2π/d (i.e. d sin θ = nλ)
  6. Bloch oscillations in a tilted lattice: period T_B = 2π/(F a)
  7. Stationary states by imaginary-time propagation: E_n = n + 1/2
  8. Time-step convergence order: O(dt²) for smooth (also time-dependent) V, O(dt) at sharp edges;
     single precision (complex64) vs double
  9. Absorbing boundary: residual (reflected + wrapped) probability vs wavenumber
 10. Flux detector: live transmission through a smooth sech² barrier (exact T(E)), unaffected by
     the absorber; spectral vs 4th-order finite-difference current
 11. Gross–Pitaevskii: moving bright soliton (exact solution), O(dt²), energy
 12. Gross–Pitaevskii ground state (imaginary time): Thomas–Fermi limit, stationary in real time

Each panel is one function; the expensive inner loops (independent runs) go to a process pool.

Run:  python validate.py                 all checks + figures/validation.png; exits non-zero on failure
      python validate.py --quick         all checks, no figure and no illustration-only runs (CI)
      python validate.py --jobs 1        run serially (default: one process per CPU)
      python validate.py --update-docs   also rewrite the tables in docs/validation.md
"""
import argparse
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qwave.analytic import (
    T_packet,
    T_rect,
    T_sech2,
    bloch_period,
    bright_soliton,
    bright_soliton_energy,
    driven_ho_x,
    fringe_period_ky,
    ho_levels,
    thomas_fermi_mu_1d,
)
from qwave.eigen import eigenstates, ground_state
from qwave.observables import FluxDetector
from qwave.potentials import (
    check_on_grid,
    cosine_lattice,
    double_slit,
    harmonic,
    rect_barrier,
    sech2_barrier,
    slits,
    tilt,
)
from qwave.solver import Grid, Solver, absorber_residual, gaussian_packet

results = []          # (name_en, name_ja, error, tol, ok)


def check(name, name_ja, err, tol):
    ok = err < tol
    results.append((name, name_ja, err, tol, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: error = {err:.2e} (tol {tol:.2g})", flush=True)


# ================================================================ workers (run in the process pool)
TUNNEL = dict(V0=1.5, a=1.0, sigma=20.0)


def tunnel_run(k0):
    """Transmitted probability through a rectangular barrier for one packet (dx = 1/16: 16 cells)."""
    V0, a, sigma = TUNNEL["V0"], TUNNEL["a"], TUNNEL["sigma"]
    g = Grid(8192, 512.0)
    s = Solver(g, rect_barrier(g, V0, a), dt=0.02)
    psi = s.step(gaussian_packet(g, -120, sigma, k0), int(240 / k0 / s.dt))   # until fully scattered
    return float(np.sum(np.abs(psi[g.x[0] > a]) ** 2) * g.dV), T_packet(k0, sigma, V0, a)


SLIT = dict(d=5.0, w=1.125, wall_x=-10.0, thick=0.5)


def far_field(centers):
    """P(k_y) of the wave transmitted through openings at the given y-centres."""
    g = Grid((512, 512), (64, 64))                      # dx = 0.125
    V = slits(g, centers, SLIT["wall_x"], SLIT["thick"], SLIT["w"], height=50)
    s = Solver(g, V, dt=0.005, absorber=dict(gamma_max=10))
    psi = s.step(gaussian_packet(g, (-20, 0), (2.0, 8.0), (4.0, 0)), 1000)          # t = 5
    A = np.fft.fftn(np.where(g.x[0] > SLIT["wall_x"] + SLIT["thick"], psi, 0))
    order = np.argsort(g.kaxes[1])
    return g.kaxes[1][order], (np.abs(A) ** 2).sum(axis=0)[order]


CONV_DTS = np.array([0.04, 0.02, 0.01, 0.005])
DRIVE = dict(omega=1.0, F0=0.5, Omega=1.7)
CONV_CASES = ["harmonic (smooth)", "barrier, sharp edges", "barrier, tanh edges (0.25)", "driven oscillator V(x, t)"]


def convergence_errors(label, T=2.0, ref_dt=1e-4):
    """L2 error at t = T for every dt in CONV_DTS against a dt = 1e-4 reference."""
    g = Grid(2048, 100)
    X = g.x[0]
    psi0 = gaussian_packet(g, -10, 2.0, 4.0)
    if label == "harmonic (smooth)":
        V = harmonic(g, 0.3)
    elif label == "barrier, sharp edges":
        V = rect_barrier(g, 10, 1.0)
    elif label == "barrier, tanh edges (0.25)":
        V = rect_barrier(g, 10, 1.0, edge=0.25)
    else:
        def V(t):
            return 0.5 * X**2 - DRIVE["F0"] * np.cos(DRIVE["Omega"] * t) * X
        psi0 = gaussian_packet(g, 1.0, 1 / np.sqrt(2), 0.0)
    ref = Solver(g, V, ref_dt).step(psi0, int(round(T / ref_dt)))
    return np.array([np.sqrt(g.norm(Solver(g, V, dt).step(psi0, int(round(T / dt))) - ref)) for dt in CONV_DTS])


def absorber_point(args):
    k, width, gamma_max = args
    return absorber_residual(k, width, gamma_max)


# ================================================================ panels (main process)
def panel_norm(ax, run, quick):
    g = Grid((320, 320), (40, 40))
    s = Solver(g, double_slit(g, height=50), dt=0.01)
    psi = gaussian_packet(g, (-10, 0), 1.5, (5, 0))
    ts, norms = [], []
    for _ in range(100):
        psi = s.step(psi, 5)
        ts.append(s.t); norms.append(g.norm(psi))
    check("2D norm conservation (double slit)", "2D ノルム保存（二重スリット）", max(abs(np.array(norms) - 1)), 1e-10)
    if ax is not None:
        ax.plot(ts, np.array(norms) - 1)
        ax.set(title="1. Norm conservation (2D, double slit)", xlabel="t", ylabel="‖ψ‖² − 1")


def panel_free(ax, run, quick):
    g = Grid(4096, 400)
    s = Solver(g, np.zeros(g.n), dt=0.02)
    sigma0 = 2.0
    psi = gaussian_packet(g, -50, sigma0, 2.0)
    ts, sig = [0.0], [g.std_x(psi)]
    for _ in range(200):
        psi = s.step(psi, 10)
        ts.append(s.t); sig.append(g.std_x(psi))
    ts = np.array(ts)
    sig_exact = sigma0 * np.sqrt(1 + (ts / (2 * sigma0**2)) ** 2)
    check("Free Gaussian spreading σ(t) (relative)", "自由粒子の波束の広がり σ(t)（相対）",
          max(abs(np.array(sig) - sig_exact) / sig_exact), 1e-6)
    if ax is not None:
        ax.plot(ts, sig_exact, "k-", lw=3, alpha=0.3, label="analytic")
        ax.plot(ts[::10], sig[::10], "o", ms=4, label="numerical")
        ax.set(title="2. Free Gaussian spreading", xlabel="t", ylabel="σ(t)")
        ax.legend()


def panel_harmonic(ax, run, quick):
    g = Grid(1024, 40)
    omega, x0 = 1.0, 3.0
    s = Solver(g, harmonic(g, omega), dt=0.005)
    psi = gaussian_packet(g, x0, 1 / np.sqrt(2 * omega), 0.0)   # ground-state width -> coherent state
    E0 = s.energy(psi)
    ts, xs, Es = [0.0], [g.expect_x(psi)], [E0]
    for _ in range(400):
        psi = s.step(psi, 10)
        ts.append(s.t); xs.append(g.expect_x(psi)); Es.append(s.energy(psi))
    ts = np.array(ts)
    drift = max(abs(np.array(Es) - E0)) / E0
    check("Harmonic-oscillator coherent state ⟨x⟩ = x₀ cos ωt", "調和振動子コヒーレント状態 ⟨x⟩ = x₀ cos ωt",
          max(abs(np.array(xs) - x0 * np.cos(omega * ts))), 1e-4)
    check("Harmonic-oscillator energy conservation (relative)", "調和振動子のエネルギー保存（相対）", drift, 1e-4)
    if ax is not None:
        ax.plot(ts, x0 * np.cos(omega * ts), "k-", lw=3, alpha=0.3, label="x0 cos ωt")
        ax.plot(ts[::8], xs[::8], "o", ms=3, label="⟨x⟩ numerical")
        ax.set(title=f"3. Harmonic oscillator (ΔE/E < {drift:.0e})", xlabel="t")
        ax.legend(loc="upper right")


def panel_tunneling(ax, run, quick):
    k0s = np.linspace(1.0, 2.2, 5 if quick else 13)
    T_num, T_ana = np.array(run(tunnel_run, k0s)).T
    check("Rectangular-barrier transmission T(E) (analytic, averaged over the packet's momentum distribution)",
          "矩形障壁の透過率 T(E)（運動量分布で平均した解析解と比較）", max(abs(T_num - T_ana)), 5e-3)
    if ax is not None:
        V0, a = TUNNEL["V0"], TUNNEL["a"]
        Efine = np.linspace(0.3, 2.6, 400)
        ax.plot(Efine, T_rect(Efine, V0, a), "k-", lw=3, alpha=0.3, label="analytic T(E)")
        ax.plot(k0s**2 / 2, T_num, "o", label="wave-packet simulation")
        ax.axvline(V0, ls=":", c="gray")
        ax.set(title=f"4. Tunneling (V0={V0}, a={a})", xlabel="E = k0²/2", ylabel="T")
        ax.legend()


def panel_double_slit(ax, run, quick):
    # The far field is the momentum distribution of the transmitted wave. With the amplitudes A± of
    # each slit alone, P_double = |A+ + A-|² = P+ + P- + 2 Re(A+* A-), and the interference term
    # oscillates as cos(k_y d): its zeros are spaced exactly π/d, whatever the single-slit envelope.
    d = SLIT["d"]
    check_on_grid(0.125, d / 2, SLIT["w"], SLIT["wall_x"], SLIT["thick"])
    (ky, P2), (_, Pp), (_, Pm) = run(far_field, [(d / 2, -d / 2), (d / 2,), (-d / 2,)])
    interference = (P2 - Pp - Pm) / (2 * np.sqrt(Pp * Pm))
    inside = np.abs(ky) < 4.5
    zeros = [ky[i] - interference[i] * (ky[i + 1] - ky[i]) / (interference[i + 1] - interference[i])
             for i in np.flatnonzero(inside[:-1] & inside[1:]) if interference[i] * interference[i + 1] < 0]
    d_eff = np.pi / np.polyfit(np.arange(len(zeros)), zeros, 1)[0]                 # zero spacing = π/d
    check(f"Double-slit fringe period Δk_y = 2π/d (relative; {len(zeros)} fringe zeros)",
          f"二重スリットの縞の周期 Δk_y = 2π/d（相対、縞のゼロ点 {len(zeros)} 個）", abs(d_eff - d) / d, 5e-3)
    check("Double slit: mirror symmetry of the two single-slit far fields P₊(k_y) = P₋(−k_y) (relative)",
          "二重スリット: 片方ずつの遠方パターンの鏡映対称性 P₊(k_y) = P₋(−k_y)（相対）",
          np.max(np.abs(Pp[1:] - Pm[1:][::-1])) / Pp.max(), 1e-10)
    if ax is not None:
        model = Pp + Pm + 2 * np.sqrt(Pp * Pm) * np.cos(ky * d)
        norm = P2.max()
        ax.plot(ky, model / norm, "k-", lw=3, alpha=0.3, label="P₊ + P₋ + 2√(P₊P₋) cos(k_y d)")
        ax.plot(ky, P2 / norm, ".", ms=4, label="double slit, simulated")
        for n in range(-3, 4):
            ax.axvline(n * fringe_period_ky(d), ls=":", c="gray", lw=0.8)
        ax.set(title=f"5. Double slit far field (d = {d:g}, d_eff = {d_eff:.4f})", xlabel="k_y",
               ylabel="P(k_y) (normalized)", xlim=(-5, 5))
        ax.legend(loc="upper right", fontsize=8)


def panel_bloch(ax, run, quick):
    # Ground state of a lattice + weak trap, then the trap is replaced by a uniform force F.
    # The packet oscillates in space (instead of accelerating) and returns every T_B = 2π/(F a).
    g = Grid(4096, 256)
    lat_a, depth, F = 1.0, 4.0, 0.5
    lattice = cosine_lattice(g, depth, lat_a)
    _, psi = ground_state(g, lattice + harmonic(g, 0.3), dtau=0.005, tol=1e-10)
    s = Solver(g, lattice + tilt(g, F), dt=0.005)
    TB = bloch_period(F, lat_a)
    ts, xs = [0.0], [g.expect_x(psi)]
    while s.t < 3.3 * TB:
        psi = s.step(psi, 10)
        ts.append(s.t); xs.append(g.expect_x(psi))
    ts, xs = np.array(ts), np.array(xs)
    returns = [0.0]                                     # ⟨x⟩ is maximal (back at the start) at t = n T_B
    for i in range(1, len(xs) - 1):
        if xs[i] > xs[i - 1] and xs[i] > xs[i + 1]:
            p, q, r = xs[i - 1], xs[i], xs[i + 1]
            returns.append(ts[i] + 0.5 * (p - r) / (p - 2 * q + r) * (ts[1] - ts[0]))
    TB_num = np.polyfit(np.arange(len(returns)), returns, 1)[0]
    check(f"Bloch-oscillation period T_B = 2π/(F a) (relative; {len(returns) - 1} periods)",
          f"ブロッホ振動の周期 T_B = 2π/(F a)（相対、{len(returns) - 1} 周期）", abs(TB_num - TB) / TB, 1e-3)
    if ax is not None:
        ax.plot(ts, xs)
        for n in range(4):
            ax.axvline(n * TB, ls=":", c="gray")
        ax.set(title=f"6. Bloch oscillation (T_B = {TB:.3f}, numerical {TB_num:.3f})", xlabel="t", ylabel="⟨x⟩")


def panel_eigenstates(ax, run, quick):
    g = Grid(512, 30)
    nstates = 5
    E, states = eigenstates(g, harmonic(g), nstates, dtau=0.01)
    check(f"Imaginary-time eigenstates: harmonic oscillator E₀…E{nstates - 1} = n + 1/2",
          f"虚時間発展による固有状態: 調和振動子 E₀…E{nstates - 1} = n + 1/2",
          max(abs(E - ho_levels(nstates))), 1e-6)
    if ax is not None:
        xg = g.x[0]
        ax.plot(xg, harmonic(g), "k-", lw=1)
        for En, phi in zip(E, states):
            ax.axhline(En, c="gray", lw=0.5)
            ax.plot(xg, En + 0.8 * phi.real)
        ax.set(title="7. Imaginary time: HO eigenstates (offset by Eₙ)", xlabel="x", ylabel="E",
               xlim=(-6, 6), ylim=(0, nstates + 0.5))


def panel_convergence(ax, run, quick):
    errors = dict(zip(CONV_CASES, run(convergence_errors, CONV_CASES)))
    orders = {label: np.polyfit(np.log(CONV_DTS[1:]), np.log(err[1:]), 1)[0] for label, err in errors.items()}
    check("Time-step convergence order, smooth V: |p − 2|", "時間刻みの収束次数（滑らかな V）: |p − 2|",
          abs(orders["harmonic (smooth)"] - 2), 0.05)
    check("Time-step convergence order, time-dependent V: |p − 2|", "時間刻みの収束次数（時間依存 V）: |p − 2|",
          abs(orders["driven oscillator V(x, t)"] - 2), 0.05)
    check("Time-step convergence order, sharp barrier edges: |p − 1| (drops to ≈ first order)",
          "時間刻みの収束次数（障壁の角が鋭い場合）: |p − 1|（ほぼ 1 次に落ちる）",
          abs(orders["barrier, sharp edges"] - 1), 0.25)

    # driven oscillator against the exact ⟨x⟩(t) (Ehrenfest is exact for quadratic V)
    g = Grid(512, 40)
    X = g.x[0]
    s = Solver(g, lambda t: 0.5 * X**2 - DRIVE["F0"] * np.cos(DRIVE["Omega"] * t) * X, dt=0.005)
    psi = gaussian_packet(g, 1.0, 1 / np.sqrt(2), 0.0)
    dev = 0.0
    for _ in range(200):
        psi = s.step(psi, 10)
        dev = max(dev, abs(g.expect_x(psi) - driven_ho_x(s.t, 1.0, **DRIVE)))
    check("Driven oscillator (time-dependent V) ⟨x⟩(t)", "強制振動子（時間依存 V）の ⟨x⟩(t)", dev, 1e-4)

    # single precision, as in a WASM/WebGL port: 2D hard disk, 1000 steps
    # (with NumPy >= 2 the FFT itself runs in single precision too, as it would in the browser)
    g = Grid((320, 320), (40, 40))
    V = np.where(g.x[0] ** 2 + g.x[1] ** 2 <= 4.0, 200.0, 0.0)
    p0 = gaussian_packet(g, (-11, 0.8), 2.5, (4.0, 0))
    p64 = Solver(g, V, 0.005).step(p0, 1000)
    p32 = Solver(g, V, 0.005, dtype=np.complex64).step(p0, 1000)
    check("Single precision (complex64) vs double, 2D, 1000 steps (relative L2)",
          "単精度（complex64）と倍精度の差、2D・1000 ステップ（相対 L2）", np.sqrt(g.norm(p32 - p64)), 2e-4)
    check("Single precision (complex64) norm drift, 1000 steps", "単精度（complex64）のノルムのずれ、1000 ステップ",
          abs(g.norm(p32) - 1), 1e-3)
    s32 = Solver(g, V, 0.005, dtype=np.complex64, renormalize_every=100)
    p32r, drift = p0, 0.0
    for _ in range(100):                                # sample between the renormalizations too
        p32r = s32.step(p32r, 10)
        drift = max(drift, abs(g.norm(p32r) - 1))
    check("Single precision with renormalize_every=100: max norm drift over 1000 steps",
          "単精度＋100 ステップごとの再規格化: 1000 ステップ中のノルムの最大のずれ", drift, 3e-5)
    if ax is not None:
        for label, err in errors.items():
            ax.loglog(CONV_DTS, err, "o-", label=f"{label}: p = {orders[label]:.2f}")
        ax.loglog(CONV_DTS, 0.5 * CONV_DTS**2 / CONV_DTS[0] ** 2 * 1e-3, "k:", lw=1, label="∝ dt²")
        ax.set(title="8. Time-step convergence (L2 error at t = 2)", xlabel="dt", ylabel="‖ψ − ψ_ref‖")
        ax.legend(fontsize=7, loc="lower right")


def panel_absorber(ax, run, quick):
    ks = np.array([0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8])
    layers = [(32.0, 2.0, "width 32, γmax 2"), (3.2, 2.0, "width 3.2, γmax 2"), (3.2, 10.0, "width 3.2, γmax 10")]
    if quick:                                           # only what the check needs
        ks, layers = ks[ks >= 1], layers[:1]
    jobs = [(k, width, gmax) for width, gmax, _ in layers for k in ks]
    values = np.array(run(absorber_point, jobs)).reshape(len(layers), len(ks))
    residual = {label: row for (_, _, label), row in zip(layers, values)}
    wide = residual["width 32, γmax 2"]
    check("Absorbing layer, residual (reflected + wrapped) probability for 1 ≤ k ≤ 8",
          "吸収層の残存確率（反射＋折り返し）、1 ≤ k ≤ 8", wide[(ks >= 1) & (ks <= 8)].max(), 1e-5)
    if ax is not None:
        for label, row in residual.items():
            ax.semilogy(ks, row, "o-", label=label)
        ax.set(title="9. Absorbing layer: residual probability", xlabel="k", ylabel="reflected + wrapped")
        ax.legend(fontsize=8)


def panel_flux(ax, run, quick):
    # Smooth barrier (exact T(E) known) in a box whose absorber eats the transmitted packet before the end:
    # summing |ψ|² behind the barrier then fails, integrating the probability current through a line does not.
    V0, w, sigma, k0 = 1.5, 0.5, 10.0, 1.6
    g = Grid(2**13, 2**13 * 0.05)                        # L = 409.6
    x = g.x[0]
    s = Solver(g, sech2_barrier(g, V0, w), dt=0.02, absorber=dict(width=40.0, gamma_max=1.0))
    psi = gaussian_packet(g, -100, sigma, k0)
    det, det_fd4 = FluxDetector(g, 10.0), FluxDetector(g, 10.0, method="fd4")
    det.record(0, psi); det_fd4.record(0, psi)
    ts, T_live, P_behind = [0.0], [0.0], [0.0]
    for _ in range(125):                                 # t = 250: the transmitted packet is absorbed from t ≈ 160 on
        psi = s.step(psi, 100, callback=lambda t, p: (det.record(t, p), det_fd4.record(t, p)))
        ts.append(s.t); T_live.append(det.transmitted)
        P_behind.append(np.sum(np.abs(psi[x >= det.position]) ** 2) * g.dV)
    T_exact = T_packet(k0, sigma, V0, w, T=T_sech2)
    check("Flux detector: transmission through a sech² barrier vs exact T(E), absorber active",
          "フラックス検出器: sech² 障壁の透過率と厳密解（吸収層あり）", abs(det.transmitted - T_exact), 2e-4)
    check("Flux detector: 4th-order finite-difference vs spectral current",
          "フラックス検出器: 4 次差分と スペクトル法の電流の差", abs(det_fd4.transmitted - det.transmitted), 1e-5)
    if ax is not None:
        ax.plot(ts, T_live, lw=2, label="∫ flux dt (detector at x = 10)")
        ax.plot(ts, P_behind, "--", label="Σ|ψ|² behind the barrier")
        ax.axhline(T_exact, c="k", ls=":", label=f"exact T = {T_exact:.4f}")
        ax.set(title="10. Live transmission (absorber eats the packet)", xlabel="t", ylabel="probability")
        ax.legend(fontsize=8, loc="center right")


def panel_soliton(ax, run, quick):
    gnl, v, T_end = -2.0, 1.0, 20.0
    g = Grid(1024, 102.4)
    x = g.x[0]
    psi0 = bright_soliton(x, 0.0, gnl, v, x0=-20)
    sol_err = []
    for dt in (0.004, 0.002, 0.001):
        s = Solver(g, np.zeros(g.n), dt, nonlinearity=gnl)
        psi = s.step(psi0, int(round(T_end / dt)))
        sol_err.append(np.sqrt(g.norm(psi - bright_soliton(x, T_end, gnl, v, x0=-20))))
    check("Bright soliton (Gross–Pitaevskii, g < 0) vs exact solution at t = 20 (L2)",
          "明るいソリトン（GP 方程式, g < 0）と厳密解の差、t = 20（L2）", sol_err[1], 1e-4)
    check("Bright soliton: convergence order in dt: |p − 2|", "明るいソリトン: 時間刻みの収束次数 |p − 2|",
          abs(np.log2(sol_err[0] / sol_err[1]) - 2), 0.05)
    check("Bright soliton: energy vs v²/2 − g²/24 (relative)", "明るいソリトン: エネルギーと v²/2 − g²/24 の差（相対）",
          abs(s.energy(psi) - bright_soliton_energy(gnl, v)) / abs(bright_soliton_energy(gnl, v)), 1e-8)
    if ax is not None:
        for t_snap in (0.0, 10.0, 20.0):
            ax.plot(x, np.abs(bright_soliton(x, t_snap, gnl, v, x0=-20)) ** 2, "k-", lw=3, alpha=0.3)
        s = Solver(g, np.zeros(g.n), 0.002, nonlinearity=gnl)
        psi = psi0
        for t_snap in (0.0, 10.0, 20.0):
            psi = s.step(psi, int(round((t_snap - s.t) / s.dt)))
            ax.plot(x[::6], np.abs(psi[::6]) ** 2, "o", ms=3, label=f"t = {t_snap:g}")
        ax.set(title=f"11. Bright soliton, g = {gnl:g}, v = {v:g} (grey: exact)", xlabel="x", ylabel="|ψ|²",
               xlim=(-30, 10))
        ax.legend(fontsize=8)


def panel_gp_ground_state(ax, run, quick):
    gnl = 500.0
    g = Grid(256, 40)
    x = g.x[0]
    mu, psi0 = ground_state(g, harmonic(g), dtau=0.002, nonlinearity=gnl)
    mu_tf = thomas_fermi_mu_1d(gnl)
    check(f"Gross–Pitaevskii ground state, μ vs Thomas–Fermi limit (relative; g = {gnl:g}, μ ≈ {mu:.0f}ω)",
          f"GP 方程式の基底状態、μ とトーマス・フェルミ極限の差（相対、g = {gnl:g}, μ ≈ {mu:.0f}ω）",
          abs(mu - mu_tf) / mu_tf, 1e-3)
    s = Solver(g, harmonic(g), 0.005, nonlinearity=gnl)
    psi = s.step(psi0, 2000)                             # t = 10
    overlap = np.sum(np.conj(psi0) * psi) * g.dV         # should be e^{-iμt}
    check("Gross–Pitaevskii ground state is stationary: 1 − |⟨ψ(0)|ψ(10)⟩| and phase error of e^{−iμt}",
          "GP 基底状態の定常性: 1 − |⟨ψ(0)|ψ(10)⟩| と位相 e^{−iμt} のずれ",
          max(1 - abs(overlap), abs(np.angle(overlap * np.exp(1j * mu * s.t)))), 1e-5)
    if ax is not None:
        ax.plot(x, np.maximum(mu_tf - harmonic(g), 0) / gnl, "k-", lw=3, alpha=0.3, label="Thomas–Fermi")
        ax.plot(x, np.abs(psi0) ** 2, label="imaginary time, t = 0")
        ax.plot(x[::4], np.abs(psi[::4]) ** 2, "o", ms=3, label="after real-time t = 10")
        ax.set(title=f"12. GP ground state, g = {gnl:g} (μ = {mu:.3f}, TF {mu_tf:.3f})", xlabel="x",
               ylabel="|ψ|²", xlim=(-15, 15))
        ax.legend(fontsize=8)


PANELS = [panel_norm, panel_free, panel_harmonic, panel_tunneling, panel_double_slit, panel_bloch,
          panel_eigenstates, panel_convergence, panel_absorber, panel_flux, panel_soliton, panel_gp_ground_state]


# ================================================================ docs tables
def sci(v):
    """0.00041 -> '4 × 10⁻⁴'"""
    if v == 0:
        return "0"
    e = int(np.floor(np.log10(abs(v))))
    m = round(v / 10**e)
    if m == 10:
        m, e = 1, e + 1
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m} × 10{sup}"


def update_tables(path):
    text = open(path, encoding="utf-8").read()
    for lang, idx, header in (("en", 0, "| Test | Error |"), ("ja", 1, "| テスト | 誤差 |")):
        rows = [header, "|---|---|"] + [f"| {r[idx].replace('|', chr(92) + '|')} | {sci(r[2])} |" for r in results]
        block = "\n".join(rows)
        pattern = re.compile(rf"(<!-- validation-table:{lang}:start -->\n).*?(<!-- validation-table:{lang}:end -->)",
                             re.S)
        if not pattern.search(text):
            raise SystemExit(f"markers for '{lang}' not found in {path}")
        text = pattern.sub(lambda m, block=block: m.group(1) + block + "\n" + m.group(2), text)
    open(path, "w", encoding="utf-8").write(text)
    print(f"updated tables in {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true", help="all checks, no figure, no illustration-only runs")
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="worker processes (1 = serial)")
    parser.add_argument("--update-docs", action="store_true", help="rewrite the tables in docs/validation.md")
    args = parser.parse_args()
    if args.quick and args.update_docs:
        parser.error("--update-docs needs the full run (drop --quick)")

    t_start = time.time()
    fig, axs = (None, [None] * len(PANELS)) if args.quick else plt.subplots(4, 3, figsize=(15, 16))
    axes = axs if args.quick else axs.flat
    pool = ProcessPoolExecutor(args.jobs) if args.jobs > 1 else None
    try:
        def run(fn, items):
            items = list(items)
            return list(pool.map(fn, items)) if pool else [fn(i) for i in items]

        for panel, ax in zip(PANELS, axes):
            panel(ax, run, args.quick)
    finally:
        if pool:
            pool.shutdown()

    if fig is not None:
        fig.tight_layout()
        os.makedirs("figures", exist_ok=True)
        fig.savefig("figures/validation.png", dpi=110)
        print("saved figures/validation.png")
    if args.update_docs:
        update_tables(os.path.join("docs", "validation.md"))
    failed = sum(not r[4] for r in results)
    print(f"{len(results) - failed}/{len(results)} checks passed in {time.time() - t_start:.0f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
