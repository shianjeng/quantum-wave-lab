"""
Numerical validation against analytic results (ħ = m = 1).

  1. Norm conservation (2D, no absorber)           -> should stay 1 to machine precision
  2. Free Gaussian spreading:  σ(t) = σ0 √(1 + (t / 2σ0²)²)
  3. Harmonic oscillator coherent state: ⟨x⟩(t) = x0 cos(ωt), energy constant
  4. Tunneling through a rectangular barrier: transmitted probability vs
     analytic T(E) averaged over the packet's momentum distribution

Run:  python validate.py        (writes figures/validation.png, exits non-zero on failure)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qwave.analytic import T_packet, T_rect
from qwave.potentials import double_slit, harmonic, rect_barrier_1d
from qwave.solver import Grid, Solver, gaussian_packet

results = []


def check(name, err, tol):
    ok = err < tol
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: error = {err:.2e} (tol {tol:.0e})")


fig, axs = plt.subplots(2, 2, figsize=(11, 8))

# ---------------------------------------------------------------- 1. norm (2D)
g = Grid((256, 256), (40, 40))
s = Solver(g, double_slit(g, height=50), dt=0.01)
psi = gaussian_packet(g, (-10, 0), 1.5, (5, 0))
ts, norms = [], []
for _ in range(100):
    psi = s.step(psi, 5)
    ts.append(s.t); norms.append(g.norm(psi))
check("2D norm conservation", max(abs(np.array(norms) - 1)), 1e-10)
axs[0, 0].plot(ts, np.array(norms) - 1)
axs[0, 0].set(title="1. Norm conservation (2D, double slit)", xlabel="t", ylabel="‖ψ‖² − 1")

# ---------------------------------------------------------------- 2. free spreading
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
check("Free packet width", max(abs(np.array(sig) - sig_exact) / sig_exact), 1e-6)
axs[0, 1].plot(ts, sig_exact, "k-", lw=3, alpha=0.3, label="analytic")
axs[0, 1].plot(ts[::10], sig[::10], "o", ms=4, label="numerical")
axs[0, 1].set(title="2. Free Gaussian spreading", xlabel="t", ylabel="σ(t)")
axs[0, 1].legend()

# ---------------------------------------------------------------- 3. harmonic oscillator
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
check("HO ⟨x⟩(t) = x0 cos ωt", max(abs(np.array(xs) - x0 * np.cos(omega * ts))), 1e-4)
check("HO energy drift (relative)", max(abs(np.array(Es) - E0)) / E0, 1e-4)
print(f"       E_numerical = {E0:.6f}, E_analytic = ω/2 + ω²x0²/2 = {omega/2 + omega**2 * x0**2 / 2:.6f}")
axs[1, 0].plot(ts, x0 * np.cos(omega * ts), "k-", lw=3, alpha=0.3, label="x0 cos ωt")
axs[1, 0].plot(ts[::8], xs[::8], "o", ms=3, label="⟨x⟩ numerical")
axs[1, 0].set(title=f"3. Harmonic oscillator (ΔE/E < {max(abs(np.array(Es)-E0))/E0:.0e})", xlabel="t")
axs[1, 0].legend(loc="upper right")

# ---------------------------------------------------------------- 4. tunneling
V0, a, sigma = 1.5, 1.0, 20.0
g = Grid(2**15, 2**15 * 0.05)                       # dx = 0.05 -> barrier is exactly 20 cells
k0s = np.linspace(1.0, 2.2, 13)
T_num, T_ana = [], []
for k0 in k0s:
    s = Solver(g, rect_barrier_1d(g, V0, a), dt=0.02)
    psi = gaussian_packet(g, -300, sigma, k0)
    psi = s.step(psi, int(600 / k0 / s.dt))           # until the packet has fully scattered
    T_num.append(np.sum(np.abs(psi[g.x[0] > a]) ** 2) * g.dV)
    T_ana.append(T_packet(k0, sigma, V0, a))
check("Tunneling T(E) vs analytic", max(abs(np.array(T_num) - np.array(T_ana))), 5e-3)
Efine = np.linspace(0.3, 2.6, 400)
axs[1, 1].plot(Efine, T_rect(Efine, V0, a), "k-", lw=3, alpha=0.3, label="analytic T(E)")
axs[1, 1].plot(k0s**2 / 2, T_num, "o", label="wave-packet simulation")
axs[1, 1].axvline(V0, ls=":", c="gray")
axs[1, 1].set(title=f"4. Tunneling (V0={V0}, a={a})", xlabel="E = k0²/2", ylabel="T")
axs[1, 1].legend()

fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/validation.png", dpi=130)
print("saved figures/validation.png")
sys.exit(0 if all(results) else 1)
