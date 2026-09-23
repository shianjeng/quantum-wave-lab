"""Ready-made potentials. In the web version these are what the user 'draws' with the mouse."""
import numpy as np

from .solver import Grid


def harmonic(grid: Grid, omega=1.0):
    return 0.5 * omega**2 * sum(xi**2 for xi in grid.x)


def rect_barrier_1d(grid: Grid, height, width, center=0.0):
    """Half-open interval so that (number of cells) × dx == width when width/dx is an integer."""
    return np.where(_interval(grid.x[0], center, width, grid.dx[0]), height, 0.0)


def _interval(u, center, width, du):
    """Half-open [center - w/2, center + w/2): contains exactly width/du grid points
    when width/du is an integer and the edges fall on grid points."""
    eps = 1e-9 * du
    return (u - center >= -width / 2 - eps) & (u - center < width / 2 - eps)


def double_slit(grid: Grid, wall_x=0.0, thickness=0.5, slit_width=1.25, slit_sep=5.0, height=1e3):
    """Vertical wall at x = wall_x with two openings centered at y = ±slit_sep/2.
    Choose dimensions that are integer multiples of dx (see check_on_grid)."""
    x, y = grid.x
    wall = _interval(x, wall_x, thickness, grid.dx[0])
    slit1 = _interval(y, slit_sep / 2, slit_width, grid.dx[1])
    slit2 = _interval(y, -slit_sep / 2, slit_width, grid.dx[1])
    return np.where(wall & ~(slit1 | slit2), height, 0.0)


def check_on_grid(dx, *lengths):
    """Raise if a geometric length is not an integer number of cells —
    otherwise the effective width silently differs from the nominal one."""
    for L in lengths:
        n = L / dx
        if abs(n - round(n)) > 1e-9:
            raise ValueError(f"length {L} is {n:.3f} cells (dx = {dx}); pick a multiple of dx")


def paint_disk(V, grid: Grid, center, radius, height):
    """Mouse-brush primitive: add a disk of potential. Returns a new array."""
    x, y = grid.x
    V = V.copy()
    V[(x - center[0]) ** 2 + (y - center[1]) ** 2 <= radius**2] = height
    return V
