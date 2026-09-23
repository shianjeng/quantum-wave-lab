"""Ready-made potentials. In the web version these are what the user 'draws' with the mouse."""
import numpy as np

from .solver import Grid


def harmonic(grid: Grid, omega=1.0):
    return 0.5 * omega**2 * sum(xi**2 for xi in grid.x)


def rect_barrier(grid: Grid, height, width, center=0.0, edge=0.0):
    """Barrier of the given height on center - width/2 <= x < center + width/2 (x = axis 0).
    Works on 1D and 2D grids; in 2D it is a wall that depends only on x.

    edge = 0: sharp, half-open interval, so (number of cells) × dx == width when width/dx is an
              integer. The time error of the splitting is then only O(dt) (the potential is
              discontinuous), see validate.py.
    edge > 0: smooth tanh edges of that width (same half-height points), which restores O(dt²)
              at the price of a slightly different T(E) from the textbook rectangular barrier.
    """
    x = grid.x[0]
    if edge == 0:
        return np.where(_interval(x, center, width, grid.dx[0]), height, 0.0)
    # half-height points on the same cell faces as the sharp version
    left, right = center - width / 2 - grid.dx[0] / 2, center + width / 2 - grid.dx[0] / 2
    return 0.5 * height * (np.tanh((x - left) / edge) - np.tanh((x - right) / edge))


def sech2_barrier(grid: Grid, height, width, center=0.0):
    """Smooth barrier V = height / cosh²((x - center) / width) along axis 0 (a wall in 2D).
    Its transmission is known exactly (analytic.T_sech2), and being smooth it keeps the splitting
    O(dt²) and the probability current well resolved (observables.FluxDetector)."""
    e = np.exp(-2 * np.abs(grid.x[0] - center) / width)          # sech²u = 4e^{-2|u|}/(1+e^{-2|u|})², no overflow
    return height * 4 * e / (1 + e) ** 2


def rect_barrier_1d(grid: Grid, height, width, center=0.0):
    """Backwards-compatible alias of rect_barrier (which also works in 2D)."""
    return rect_barrier(grid, height, width, center)


def _interval(u, center, width, du):
    """Half-open [center - w/2, center + w/2): contains exactly width/du grid points
    when width/du is an integer and the edges fall on grid points.

    Seen as cells [u_j - du/2, u_j + du/2), the covered region is [center - w/2, center + w/2) shifted
    by -du/2, i.e. walls and barriers built with it sit half a cell below their nominal position
    (a harmless translation for them; openings use _centered_cells instead)."""
    eps = 1e-9 * du
    return (u - center >= -width / 2 - eps) & (u - center < width / 2 - eps)


def _centered_cells(u, center, width, du):
    """The width/du cells placed symmetrically about `center`: exact width, exact position, so a set of
    openings that is mirror-symmetric in the continuum stays mirror-symmetric on the grid.

    That needs the edges center ± width/2 on cell faces, not on grid points, which is the case exactly
    when (2·center + width)/du is an odd integer. Otherwise the opening cannot be both centred and
    width/du cells wide, and a ValueError suggests the nearest widths that work."""
    n, m = width / du, 2 * center / du
    if abs(n - round(n)) > 1e-9 or abs(m - round(m)) > 1e-9 or (round(n) + round(m)) % 2 == 0:
        raise ValueError(
            f"an opening of width {width:g} ({n:g} cells) centred at {center:g} (du = {du:g}) cannot sit "
            f"symmetrically on the grid: (2·center + width)/du must be an odd integer, here it is "
            f"{n + m:g}. Use width {width - du:g} or {width + du:g}, or shift the centre by half a cell")
    return np.abs(u - center) < width / 2 - 1e-9 * du


def double_slit(grid: Grid, wall_x=0.0, thickness=0.5, slit_width=1.125, slit_sep=5.0, height=1e3):
    """Vertical wall at x = wall_x with two openings centered at y = ±slit_sep/2.
    The openings are exactly slit_width wide, slit_sep apart and mirror-symmetric about y = 0, which
    requires (slit_sep + slit_width)/dy to be an odd integer (the defaults fit dy = 0.125); the wall
    thickness should be a multiple of dx (see check_on_grid)."""
    return slits(grid, (slit_sep / 2, -slit_sep / 2), wall_x, thickness, slit_width, height)


def slits(grid: Grid, centers, wall_x=0.0, thickness=0.5, slit_width=1.125, height=1e3):
    """Vertical wall at x = wall_x with openings of width slit_width centered at y = each of `centers`
    (one center: single slit, two: double slit, many: grating). Each opening is the slit_width/dy cells
    symmetric about its centre (see _centered_cells for the condition this puts on the dimensions)."""
    x, y = grid.x
    wall = _interval(x, wall_x, thickness, grid.dx[0])
    hole = np.zeros(grid.n, dtype=bool)
    for c in centers:
        hole |= _centered_cells(y, c, slit_width, grid.dx[1])
    return np.where(wall & ~hole, height, 0.0)


def cosine_lattice(grid: Grid, depth, period, axis=0):
    """Periodic potential V = -depth · cos(2π x / period) along `axis`.
    Choose a period that divides the box length so the lattice is continuous across the boundary."""
    return -depth * np.cos(2 * np.pi * grid.x[axis] / period)


def tilt(grid: Grid, force, axis=0):
    """Linear potential V = force · x (a uniform force -force along `axis`).
    It jumps at the periodic boundary, so keep the wave packet away from the box edges."""
    return force * grid.x[axis]


def check_on_grid(dx, *lengths):
    """Raise if a geometric length is not an integer number of cells —
    otherwise the effective width silently differs from the nominal one."""
    for L in lengths:
        n = L / dx
        if abs(n - round(n)) > 1e-9:
            raise ValueError(f"length {L} is {n:.3f} cells (dx = {dx}); pick a multiple of dx")


def paint_disk(V, grid: Grid, center, radius, height, mode="set"):
    """Mouse-brush primitive: paint a disk of potential. Returns a new array.

    mode="set"  the disk's cells become `height` (overwrites what was there; use height=0 as an eraser)
    mode="add"  `height` is added to the existing potential (overlapping strokes pile up)
    """
    x, y = grid.x
    V = np.array(V, dtype=float)
    disk = (x - center[0]) ** 2 + (y - center[1]) ** 2 <= radius**2
    if mode == "set":
        V[disk] = height
    elif mode == "add":
        V[disk] += height
    else:
        raise ValueError(f"mode must be 'set' or 'add', got {mode!r}")
    return V
