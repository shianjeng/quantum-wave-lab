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

Run:  python validate.py                  (writes figures/validation.png, exits non-zero on failure)
      python validate.py --update-readme  (also rewrites the error tables in README.md)
"""
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qwave.analytic import T_packet, T_rect, bloch_period, driven_ho_x, fringe_period_ky, ho_levels
from qwave.eigen import eigenstates, ground_state
from qwave.potentials import check_on_grid, cosine_lattice, double_slit, harmonic, rect_barrier, slits, tilt
from qwave.solver import Grid, Solver, absorber_residual, gaussian_packet

results = []          # (name_en, name_ja, error, tol, ok)


def check(name, name_ja, err, tol):
    ok = err < tol
    results.append((name, name_ja, err, tol, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: error = {err:.2e} (tol {tol:.2g})")


fig, axs = plt.subplots(3, 3, figsize=(15, 12))

# ---------------------------------------------------------------- 1. norm (2D)
g = Grid((256, 256), (40, 40))
s = Solver(g, double_slit(g, height=50), dt=0.01)
psi = gaussian_packet(g, (-10, 0), 1.5, (5, 0))
ts, norms = [], []
for _ in range(100):
    psi = s.step(psi, 5)
    ts.append(s.t); norms.append(g.norm(psi))
check("2D norm conservation (double slit)", "2D ノルム保存（二重スリット）", max(abs(np.array(norms) - 1)), 1e-10)
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
check("Free Gaussian spreading σ(t) (relative)", "自由粒子の波束の広がり σ(t)（相対）",
      max(abs(np.array(sig) - sig_exact) / sig_exact), 1e-6)
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
check("Harmonic-oscillator coherent state ⟨x⟩ = x₀ cos ωt", "調和振動子コヒーレント状態 ⟨x⟩ = x₀ cos ωt",
      max(abs(np.array(xs) - x0 * np.cos(omega * ts))), 1e-4)
check("Harmonic-oscillator energy conservation (relative)", "調和振動子のエネルギー保存（相対）",
      max(abs(np.array(Es) - E0)) / E0, 1e-4)
print(f"       E_numerical = {E0:.6f}, E_analytic = ω/2 + ω²x0²/2 = {omega/2 + omega**2 * x0**2 / 2:.6f}")
axs[0, 2].plot(ts, x0 * np.cos(omega * ts), "k-", lw=3, alpha=0.3, label="x0 cos ωt")
axs[0, 2].plot(ts[::8], xs[::8], "o", ms=3, label="⟨x⟩ numerical")
axs[0, 2].set(title=f"3. Harmonic oscillator (ΔE/E < {max(abs(np.array(Es)-E0))/E0:.0e})", xlabel="t")
axs[0, 2].legend(loc="upper right")

# ---------------------------------------------------------------- 4. tunneling
V0, a, sigma = 1.5, 1.0, 20.0
g = Grid(2**15, 2**15 * 0.05)                       # dx = 0.05 -> barrier is exactly 20 cells
k0s = np.linspace(1.0, 2.2, 13)
T_num, T_ana = [], []
for k0 in k0s:
    s = Solver(g, rect_barrier(g, V0, a), dt=0.02)
    psi = gaussian_packet(g, -300, sigma, k0)
    psi = s.step(psi, int(600 / k0 / s.dt))           # until the packet has fully scattered
    T_num.append(np.sum(np.abs(psi[g.x[0] > a]) ** 2) * g.dV)
    T_ana.append(T_packet(k0, sigma, V0, a))
check("Rectangular-barrier transmission T(E) (analytic, averaged over the packet's momentum distribution)",
      "矩形障壁の透過率 T(E)（運動量分布で平均した解析解と比較）",
      max(abs(np.array(T_num) - np.array(T_ana))), 5e-3)
Efine = np.linspace(0.3, 2.6, 400)
axs[1, 0].plot(Efine, T_rect(Efine, V0, a), "k-", lw=3, alpha=0.3, label="analytic T(E)")
axs[1, 0].plot(k0s**2 / 2, T_num, "o", label="wave-packet simulation")
axs[1, 0].axvline(V0, ls=":", c="gray")
axs[1, 0].set(title=f"4. Tunneling (V0={V0}, a={a})", xlabel="E = k0²/2", ylabel="T")
axs[1, 0].legend()

# ---------------------------------------------------------------- 5. double slit
# The far field is the momentum distribution of the transmitted wave. With the amplitudes A± of
# each slit alone, P_double = |A+ + A-|² = P+ + P- + 2 Re(A+* A-), and the interference term
# oscillates as cos(k_y d): its zeros are spaced exactly π/d, whatever the single-slit envelope.
g = Grid((512, 512), (64, 64))                      # dx = 0.125
x, y = g.x
d, w, wall_x, thick = 5.0, 1.25, -10.0, 0.5
check_on_grid(g.dx[0], d / 2, w, wall_x, thick)


def far_field(centers):
    s = Solver(g, slits(g, centers, wall_x, thick, w, height=50), dt=0.005, absorber=dict(gamma_max=10))
    psi = s.step(gaussian_packet(g, (-20, 0), (2.0, 8.0), (4.0, 0)), 1000)          # t = 5
    A = np.fft.fftn(np.where(x > wall_x + thick, psi, 0))
    return (np.abs(A) ** 2).sum(axis=0)                                           # P(k_y)


order = np.argsort(g.kaxes[1])
ky = g.kaxes[1][order]
P2, Pp, Pm = (far_field(c)[order] for c in ((d / 2, -d / 2), (d / 2,), (-d / 2,)))
interference = (P2 - Pp - Pm) / (2 * np.sqrt(Pp * Pm))
inside = np.abs(ky) < 4.5
zeros = [ky[i] - interference[i] * (ky[i + 1] - ky[i]) / (interference[i + 1] - interference[i])
         for i in np.flatnonzero(inside[:-1] & inside[1:]) if interference[i] * interference[i + 1] < 0]
spacing = np.polyfit(np.arange(len(zeros)), zeros, 1)[0]                          # = π/d
d_eff = np.pi / spacing
check(f"Double-slit fringe period Δk_y = 2π/d (relative; {len(zeros)} fringe zeros)",
      f"二重スリットの縞の周期 Δk_y = 2π/d（相対、縞のゼロ点 {len(zeros)} 個）", abs(d_eff - d) / d, 5e-3)
model = Pp + Pm + 2 * np.sqrt(Pp * Pm) * np.cos(ky * d)
norm = P2.max()
axs[1, 1].plot(ky, model / norm, "k-", lw=3, alpha=0.3, label="P₊ + P₋ + 2√(P₊P₋) cos(k_y d)")
axs[1, 1].plot(ky, P2 / norm, ".", ms=4, label="double slit, simulated")
for n in range(-3, 4):
    axs[1, 1].axvline(n * fringe_period_ky(d), ls=":", c="gray", lw=0.8)
axs[1, 1].set(title=f"5. Double slit far field (d = {d}, d_eff = {d_eff:.4f})", xlabel="k_y",
              ylabel="P(k_y) (normalized)", xlim=(-5, 5))
axs[1, 1].legend(loc="upper right", fontsize=8)

# ---------------------------------------------------------------- 6. Bloch oscillations
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
axs[1, 2].plot(ts, xs)
for n in range(4):
    axs[1, 2].axvline(n * TB, ls=":", c="gray")
axs[1, 2].set(title=f"6. Bloch oscillation (T_B = {TB:.3f}, numerical {TB_num:.3f})", xlabel="t",
              ylabel="⟨x⟩")

# ---------------------------------------------------------------- 7. stationary states
g = Grid(512, 30)
nstates = 5
E, states = eigenstates(g, harmonic(g), nstates, dtau=0.01)
check(f"Imaginary-time eigenstates: harmonic oscillator E₀…E{nstates - 1} = n + 1/2",
      f"虚時間発展による固有状態: 調和振動子 E₀…E{nstates - 1} = n + 1/2",
      max(abs(E - ho_levels(nstates))), 1e-6)
xg = g.x[0]
axs[2, 0].plot(xg, harmonic(g), "k-", lw=1)
for En, phi in zip(E, states):
    axs[2, 0].axhline(En, c="gray", lw=0.5)
    axs[2, 0].plot(xg, En + 0.8 * phi.real)
axs[2, 0].set(title="7. Imaginary time: HO eigenstates (offset by Eₙ)", xlabel="x", ylabel="E",
              xlim=(-6, 6), ylim=(0, nstates + 0.5))

# ---------------------------------------------------------------- 8. time-step convergence
dts = np.array([0.04, 0.02, 0.01, 0.005])


def l2_errors(make_solver, psi0, T, ref_dt=1e-4):
    ref = make_solver(ref_dt).step(psi0, int(round(T / ref_dt)))
    return np.array([np.sqrt(g.norm(make_solver(dt).step(psi0, int(round(T / dt))) - ref)) for dt in dts])


def fitted_order(err):
    return np.polyfit(np.log(dts[1:]), np.log(err[1:]), 1)[0]


g = Grid(2048, 100)
psi0 = gaussian_packet(g, -10, 2.0, 4.0)
cases = {
    "harmonic (smooth)": (lambda dt: Solver(g, harmonic(g, 0.3), dt), psi0),
    "barrier, sharp edges": (lambda dt: Solver(g, rect_barrier(g, 10, 1.0), dt), psi0),
    "barrier, tanh edges (0.25)": (lambda dt: Solver(g, rect_barrier(g, 10, 1.0, edge=0.25), dt), psi0),
}
drive = dict(omega=1.0, F0=0.5, Omega=1.7)
X = g.x[0]
cases["driven oscillator V(x, t)"] = (
    lambda dt: Solver(g, lambda t: 0.5 * X**2 - drive["F0"] * np.cos(drive["Omega"] * t) * X, dt),
    gaussian_packet(g, 1.0, 1 / np.sqrt(2), 0.0))
orders = {}
for label, (make, p0) in cases.items():
    err = l2_errors(make, p0, T=2.0)
    orders[label] = fitted_order(err)
    axs[2, 1].loglog(dts, err, "o-", label=f"{label}: p = {orders[label]:.2f}")
axs[2, 1].loglog(dts, 0.5 * dts**2 / dts[0] ** 2 * 1e-3, "k:", lw=1, label="∝ dt²")
check("Time-step convergence order, smooth V: |p − 2|", "時間刻みの収束次数（滑らかな V）: |p − 2|",
      abs(orders["harmonic (smooth)"] - 2), 0.05)
check("Time-step convergence order, time-dependent V: |p − 2|", "時間刻みの収束次数（時間依存 V）: |p − 2|",
      abs(orders["driven oscillator V(x, t)"] - 2), 0.05)
check("Time-step convergence order, sharp barrier edges: |p − 1| (drops to ≈ first order)",
      "時間刻みの収束次数（障壁の角が鋭い場合）: |p − 1|（ほぼ 1 次に落ちる）",
      abs(orders["barrier, sharp edges"] - 1), 0.25)

# driven oscillator against the exact ⟨x⟩(t) (Ehrenfest is exact for quadratic V)
g1 = Grid(512, 40)
X1 = g1.x[0]
s = Solver(g1, lambda t: 0.5 * X1**2 - drive["F0"] * np.cos(drive["Omega"] * t) * X1, dt=0.005)
psi = gaussian_packet(g1, 1.0, 1 / np.sqrt(2), 0.0)
dev = 0.0
for _ in range(200):
    psi = s.step(psi, 10)
    dev = max(dev, abs(g1.expect_x(psi) - driven_ho_x(s.t, 1.0, **drive)))
check("Driven oscillator (time-dependent V) ⟨x⟩(t)", "強制振動子（時間依存 V）の ⟨x⟩(t)", dev, 1e-4)

# single precision, as in a WASM/WebGL port: 2D hard disk, 1000 steps
g2 = Grid((320, 320), (40, 40))
V2 = np.where(g2.x[0] ** 2 + g2.x[1] ** 2 <= 4.0, 200.0, 0.0)
p0 = gaussian_packet(g2, (-11, 0.8), 2.5, (4.0, 0))
p64 = Solver(g2, V2, 0.005).step(p0, 1000)
p32 = Solver(g2, V2, 0.005, dtype=np.complex64).step(p0, 1000)
# (with NumPy >= 2 the FFT itself runs in single precision too, as it would in the browser)
check("Single precision (complex64) vs double, 2D, 1000 steps (relative L2)",
      "単精度（complex64）と倍精度の差、2D・1000 ステップ（相対 L2）", np.sqrt(g2.norm(p32 - p64)), 2e-4)
check("Single precision (complex64) norm drift, 1000 steps", "単精度（complex64）のノルムのずれ、1000 ステップ",
      abs(g2.norm(p32) - 1), 1e-3)
s32 = Solver(g2, V2, 0.005, dtype=np.complex64, renormalize_every=100)
p32r, drift = p0, 0.0
for _ in range(100):                                # sample between the renormalizations too
    p32r = s32.step(p32r, 10)
    drift = max(drift, abs(g2.norm(p32r) - 1))
check("Single precision with renormalize_every=100: max norm drift over 1000 steps",
      "単精度＋100 ステップごとの再規格化: 1000 ステップ中のノルムの最大のずれ", drift, 3e-5)
axs[2, 1].set(title="8. Time-step convergence (L2 error at t = 2)", xlabel="dt", ylabel="‖ψ − ψ_ref‖")
axs[2, 1].legend(fontsize=7, loc="lower right")

# ---------------------------------------------------------------- 9. absorbing boundary
ks = np.array([0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8])
layers = [(32.0, 2.0, "width 32, γmax 2"), (3.2, 2.0, "width 3.2, γmax 2"), (3.2, 10.0, "width 3.2, γmax 10")]
residual = {}
for width, gmax, label in layers:
    residual[label] = np.array([absorber_residual(k, width, gmax) for k in ks])
    axs[2, 2].semilogy(ks, residual[label], "o-", label=label)
wide = residual["width 32, γmax 2"]
check("Absorbing layer, residual (reflected + wrapped) probability for 1 ≤ k ≤ 8",
      "吸収層の残存確率（反射＋折り返し）、1 ≤ k ≤ 8", wide[(ks >= 1) & (ks <= 8)].max(), 1e-5)
axs[2, 2].set(title="9. Absorbing layer: residual probability", xlabel="k", ylabel="reflected + wrapped")
axs[2, 2].legend(fontsize=8)

fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/validation.png", dpi=110)
print("saved figures/validation.png")


# ---------------------------------------------------------------- README tables
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


def update_readme(path="README.md"):
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


if "--update-readme" in sys.argv[1:]:
    update_readme()
sys.exit(0 if all(r[4] for r in results) else 1)
