"""
Gross–Pitaevskii demos (i ψ_t = [-∇²/2 + V + g|ψ|²] ψ):

  left   two bright solitons (g < 0) collide and pass through each other, shifted but unchanged
  right  a condensate is prepared in a double-well trap by imaginary-time propagation, released,
         and the two halves overlap into interference fringes (as in Andrews et al., Science 1997)

Writes figures/gpe.png.   Run:  python demo_gpe.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qwave.analytic import bright_soliton
from qwave.eigen import ground_state
from qwave.potentials import harmonic
from qwave.solver import Grid, Solver

os.makedirs("figures", exist_ok=True)
fig = plt.figure(figsize=(15, 4.2))
gs = fig.add_gridspec(1, 4, width_ratios=[1.3, 1, 1, 1])

# ---------------------------------------------------------------- soliton collision (1D)
g1, v, x0 = -8.0, 1.5, 10.0
grid = Grid(1024, 80.0)                                      # dx = 0.078: stable for dt < 3.9e-3
x = grid.x[0]
# two solitons of norm 1/2 each: √(1/2)·(soliton normalized for g/2), opposite velocities, phase π/2 apart
psi = np.sqrt(0.5) * (bright_soliton(x, 0, g1 / 2, v, x0=-x0) + 1j * bright_soliton(x, 0, g1 / 2, -v, x0=x0))
s = Solver(grid, np.zeros(grid.n), 0.002, nonlinearity=g1)
for msg in s.diagnose(psi):
    print("warning:", msg)
E0 = s.energy(psi)
rows = [np.abs(psi) ** 2]
while s.t < 14.0 - 1e-9:
    psi = s.step(psi, 50)
    rows.append(np.abs(psi) ** 2)
print(f"soliton collision: norm {grid.norm(psi):.12f}, energy drift {abs(s.energy(psi) - E0):.1e}")
ax = fig.add_subplot(gs[0])
ax.imshow(np.array(rows), origin="lower", aspect="auto", cmap="magma", extent=(x[0], x[-1], 0, s.t),
          vmax=np.max(rows) * 0.8)
ax.set(xlim=(-20, 20), xlabel="x", ylabel="t", title=f"Bright solitons collide (g = {g1:g}, v = ±{v:g})")

# ---------------------------------------------------------------- BEC interference (2D)
g2, omega, barrier = 300.0, 0.25, 8.0
grid = Grid((256, 256), (64.0, 64.0))                        # dx = 0.25: stable for dt < 2.0e-2
x, y = grid.x
trap = harmonic(grid, omega) + barrier * np.exp(-(x**2) / 2)  # harmonic trap split by a Gaussian wall
mu, psi = ground_state(grid, trap, dtau=0.005, nonlinearity=g2, tol=1e-9)
d = 2 * np.sum(np.abs(psi) ** 2 * np.abs(x)) * grid.dV       # distance between the two halves
s = Solver(grid, np.zeros(grid.n), 0.01, nonlinearity=g2, absorber=dict(width=6.0, gamma_max=3.0))
for msg in s.diagnose(psi):
    print("warning:", msg)
print(f"double-well ground state: μ = {mu:.4f}, separation of the halves d ≈ {d:.2f}")
extent = (-32, 32, -32, 32)
for i, t_snap in enumerate((0.0, 3.0, 6.0)):
    psi = s.step(psi, int(round((t_snap - s.t) / s.dt)))
    ax = fig.add_subplot(gs[i + 1])
    ax.imshow(np.abs(psi.T) ** 2, origin="lower", extent=extent, cmap="inferno")
    title = "BEC in a double well"
    if t_snap > 0:
        # fringe spacing along y = 0, against the estimate 2πt/d for two non-interacting point sources;
        # the interaction energy turns into extra expansion, so the measured fringes are ~1.5x wider here
        cut = np.abs(psi[:, grid.n[1] // 2]) ** 2
        peaks = [grid.axes[0][i] for i in range(1, grid.n[0] - 1)
                 if cut[i] > cut[i - 1] and cut[i] > cut[i + 1] and cut[i] > 0.05 * cut.max()]
        spacing = np.median(np.diff(peaks))
        print(f"t = {t_snap:g}: fringe spacing {spacing:.2f}, non-interacting estimate 2πt/d = "
              f"{2 * np.pi * t_snap / d:.2f}")
        title = (f"released, t = {t_snap:g}: fringes {spacing:.2f} apart\n"
                 f"(non-interacting: 2πt/d = {2 * np.pi * t_snap / d:.2f})")
    ax.set(title=title, xlim=(-24, 24), ylim=(-24, 24), xticks=[], yticks=[])
fig.tight_layout()
fig.savefig("figures/gpe.png", dpi=110)
print("saved figures/gpe.png")
