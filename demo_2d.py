"""
2D demos: double-slit interference and tunneling through a thin wall.
Writes figures/double_slit.png, figures/double_slit.gif, figures/tunneling_2d.png.

Run:  python demo_2d.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from qwave.analytic import T_packet
from qwave.potentials import check_on_grid, double_slit, rect_barrier
from qwave.solver import Grid, Solver, gaussian_packet

os.makedirs("figures", exist_ok=True)
N, L = 480, 60.0                                           # dx = 0.125
g = Grid((N, N), (L, L))
extent = (-L / 2, L / 2, -L / 2, L / 2)
# The default absorber (γmax = 2 over 8 % of the box) lets ~4 % of a k = 4 wave cross it and wrap
# around; γmax = 10 brings that to ~1e-5 (see absorbing_mask / validate.py panel 9).
ABSORBER = dict(gamma_max=10.0)


def render(ax, psi, V, title):
    # show |ψ| (not |ψ|²) so the weak transmitted fringes stay visible next to the reflection
    ax.imshow(np.abs(psi.T), origin="lower", extent=extent, cmap="inferno", vmax=np.max(np.abs(psi)))
    ax.contour(g.axes[0], g.axes[1], V.T, levels=[1], colors="cyan", linewidths=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])


def make_solver(V, psi0, dt, absorber=ABSORBER):
    s = Solver(g, V, dt, absorber=absorber)
    for msg in s.diagnose(psi0):
        print("warning:", msg)
    return s


def run(V, psi0, snapshot_times, dt=0.01, absorber=ABSORBER):
    s = make_solver(V, psi0, dt, absorber)
    psi, frames = psi0, []
    for t in snapshot_times:
        psi = s.step(psi, int(round((t - s.t) / dt)))
        frames.append(psi.copy())
    return frames


# ---------------------------------------------------------------- double slit
check_on_grid(g.dx[0], 5.0, 0.5, 1.125, 5.0)
V = double_slit(g, wall_x=-5, thickness=0.5, slit_width=1.125, slit_sep=5.0, height=200)
psi0 = gaussian_packet(g, (-17, 0), (2.0, 4.0), (4.0, 0))   # clear of the absorbing layer
times = [0, 3.0, 5.0, 8.0]
DT_SLIT = 0.005                                            # V·dt = 1 at the hard walls (dt = 0.01 is too coarse)
frames = run(V, psi0, times, dt=DT_SLIT)
fig, axs = plt.subplots(1, 4, figsize=(14, 3.8))
for ax, psi, t in zip(axs, frames, times):
    render(ax, psi, V, f"t = {t}")
fig.suptitle("Double slit (k0 = 4, λ ≈ 1.57, d = 5)")
fig.tight_layout()
fig.savefig("figures/double_slit.png", dpi=120)

# animation
s = make_solver(V, psi0, DT_SLIT)
state = {"psi": psi0}
fig, ax = plt.subplots(figsize=(5, 5))
im = ax.imshow(np.abs(psi0.T), origin="lower", extent=extent, cmap="inferno")
ax.contour(g.axes[0], g.axes[1], V.T, levels=[1], colors="cyan", linewidths=0.8)
ax.set_xticks([]); ax.set_yticks([])


def update(_):
    state["psi"] = s.step(state["psi"], 16)                # 0.08 time units per frame
    d = np.abs(state["psi"].T)
    im.set_data(d); im.set_clim(0, d.max())
    return (im,)


FuncAnimation(fig, update, frames=90, blit=True).save("figures/double_slit.gif", writer=PillowWriter(fps=24))

# ---------------------------------------------------------------- 2D tunneling
x, y = g.x
V0, a, k0, sig = 9.0, 1.0, 4.0, 3.0
check_on_grid(g.dx[0], a)
V = rect_barrier(g, V0, a)                                # exactly 8 cells = width 1.0
psi0 = gaussian_packet(g, (-12, 0), sig, (k0, 0))
E = Solver(g, V, 0.01).energy(psi0)                       # = k0²/2 + 1/(8σ²)·2 > k0²/2
times = [0, 2.5, 4.0, 6.0]
# nothing reaches the box edge before t = 6, so no absorber: the 2D run stays exactly comparable
# with the 1D runs below (an absorber would also nibble at the packet's far tail)
frames = run(V, psi0, times, absorber=None)
fig, axs = plt.subplots(1, 4, figsize=(14, 3.8))
for ax, psi, t in zip(axs, frames, times):
    right = np.sum(np.abs(psi[x > 0.5]) ** 2) * g.dV
    render(ax, psi, V, f"t = {t}   P(x>wall) = {right:.3f}")
fig.suptitle(f"Tunneling: ⟨E⟩ = {E:.3f} < V0 = {V0}, wall width {a}")
fig.tight_layout()
fig.savefig("figures/tunneling_2d.png", dpi=120)

# cross-check 1 (implementation): the wall depends only on x, so the problem separates and the
# 2D transmitted probability must equal a 1D run on the SAME dx and dt.
# cross-check 2 (accuracy): compare with the analytic T(E) averaged over the packet's k_x spread.
# The 2D demo grid (dx = 0.125) is chosen for speed, so a few-% discretization error is expected;
# this is exactly the accuracy/real-time trade-off the WASM version has to manage.
def transmitted_1d(dx, dt, t_end=6.0):
    n = int(1024 / dx)
    g1 = Grid(n, n * dx)
    s1 = Solver(g1, rect_barrier(g1, V0, a), dt)
    p1 = s1.step(gaussian_packet(g1, -12, sig, k0), int(round(t_end / dt)))
    return np.sum(np.abs(p1[g1.x[0] > a]) ** 2) * g1.dV


T_2d = np.sum(np.abs(frames[-1][x > 0.5]) ** 2) * g.dV    # transmitted probability at the last snapshot
T_exact = T_packet(k0, sig, V0, a)
T_same = transmitted_1d(g.dx[0], 0.01)
T_fine = transmitted_1d(1 / 64, 0.001)
print(f"2D transmitted             = {T_2d:.5f}")
print(f"1D, same dx/dt             = {T_same:.5f}   (|diff| = {abs(T_2d - T_same):.1e}  -> 2D code OK)")
print(f"1D, dx=1/64, dt=0.001      = {T_fine:.5f}")
print(f"analytic (k-averaged)      = {T_exact:.5f}   (demo-grid error {abs(T_2d - T_exact) / T_exact:.1%})")
print("saved figures/double_slit.png, double_slit.gif, tunneling_2d.png")
