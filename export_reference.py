"""
Reference data ("golden files") for the C++ / WASM port.

Every case stores the inputs as arrays (potential, absorber mask, initial ψ), so a port does not
have to reimplement the potentials, plus ψ and a few scalars after 1, 100 and 1000 steps.
The format is described in reference/README.md; tests/test_reference.py checks that this Python
implementation still reproduces the files.

Run:  python export_reference.py          (rewrites reference/)
"""
import json
import os
import shutil

import numpy as np

from qwave.analytic import bright_soliton
from qwave.observables import FluxDetector
from qwave.potentials import double_slit, harmonic, sech2_barrier
from qwave.solver import Grid, Solver, gaussian_packet

OUT = "reference"
CHECKPOINTS = (1, 100, 1000)
FORMAT_VERSION = 1


def cases():
    """name -> (grid, V, psi0, dt, absorber, nonlinearity, detectors, description)"""
    g = Grid(256, 51.2)
    yield ("free_1d", g, np.zeros(g.n), gaussian_packet(g, -10.0, 2.0, 2.0), 0.01, None, 0.0, [],
           "free Gaussian packet, periodic box")

    g = Grid(256, 25.6)
    yield ("harmonic_1d", g, harmonic(g), gaussian_packet(g, 3.0, 1 / np.sqrt(2), 0.0), 0.005, None, 0.0, [],
           "coherent state in V = x²/2")

    g = Grid(512, 102.4)
    # dt = 0.03: after 1000 steps (t = 30) both scattered packets are inside the absorbing layers
    yield ("sech2_barrier_absorber_1d", g, sech2_barrier(g, 1.5, 0.5), gaussian_packet(g, -12.0, 4.0, 1.6), 0.03,
           dict(width=16.0, gamma_max=1.0), 0.0, [dict(position=5.0, method="spectral"),
                                                  dict(position=5.0, method="fd4")],
           "smooth barrier 1.5 sech²(x/0.5), absorbing layer, flux detectors at x = 5")

    g = Grid(512, 51.2)
    yield ("soliton_gp_1d", g, np.zeros(g.n), bright_soliton(g.x[0], 0.0, -2.0, 1.0, x0=-10.0), 0.002, None, -2.0,
           [], "Gross–Pitaevskii bright soliton, g = -2, v = 1")

    g = Grid((128, 128), (32.0, 32.0))
    yield ("double_slit_2d", g, double_slit(g, wall_x=-2.0, thickness=0.5, slit_width=1.25, slit_sep=5.0, height=50.0),
           gaussian_packet(g, (-7.0, 0.0), (1.25, 3.0), (4.0, 0.0)), 0.005, dict(width=2.56, gamma_max=10.0), 0.0,
           [dict(position=0.0, method="spectral", span=(0.0, 16.0)),
            dict(position=0.0, method="spectral", span=(-16.0, 0.0))],
           "double slit with absorbing layer, flux detectors behind the upper / lower half (the on-grid slits "
           "are mirror-symmetric about y = -dy/2, not y = 0, so the two halves differ by a few %)")


def rel_l2(a, b):
    return float(np.sqrt(np.sum(np.abs(a - b) ** 2) / np.sum(np.abs(b) ** 2)))


def run_case(grid, V, psi0, dt, absorber, g, detectors, dtype=np.complex128):
    """Propagate to every checkpoint; returns [(step, t, ψ, {scalars})]."""
    s = Solver(grid, V, dt, absorber=absorber, nonlinearity=g, dtype=dtype)
    dets = [FluxDetector(grid, **d) for d in detectors]
    for d in dets:
        d.record(0.0, psi0)
    out, psi = [], psi0
    for n in CHECKPOINTS:
        psi = s.step(psi, n - s.nsteps_done, callback=lambda t, p: [d.record(t, p) for d in dets])
        scalars = dict(norm=grid.norm(psi), energy=s.energy(psi))
        if dets:
            scalars["flux_transmitted"] = [d.transmitted for d in dets]
        out.append((n, s.t, psi, scalars))
    return s, out


def save(path, arr):
    arr = np.ascontiguousarray(arr)
    if np.iscomplexobj(arr):
        arr = arr.astype("<c16")          # interleaved (re, im) float64 pairs
    else:
        arr = arr.astype("<f8")
    arr.tofile(path)


def main():
    if os.path.isdir(OUT):
        for entry in os.listdir(OUT):
            full = os.path.join(OUT, entry)
            if os.path.isdir(full):
                shutil.rmtree(full)
    manifest = dict(format_version=FORMAT_VERSION, units="hbar = m = 1", cases=[])
    for name, grid, V, psi0, dt, absorber, g, detectors, description in cases():
        os.makedirs(os.path.join(OUT, name), exist_ok=True)
        s, results = run_case(grid, V, psi0, dt, absorber, g, detectors)
        # single-precision run of the same case: sets the tolerance a float32 port has to meet
        _, results32 = run_case(grid, V, psi0, dt, absorber, g, detectors, dtype=np.complex64)
        err32 = max(rel_l2(r32[2], r[2]) for r, r32 in zip(results, results32))

        files = dict(V=f"{name}/V.f64", psi0=f"{name}/psi_0.c128")
        save(os.path.join(OUT, files["V"]), V)
        save(os.path.join(OUT, files["psi0"]), psi0)
        if s.absorber is not None:
            files["absorber_mask"] = f"{name}/absorber_mask.f64"
            save(os.path.join(OUT, files["absorber_mask"]), s.absorber)
        checkpoints = []
        for n, t, psi, scalars in results:
            fname = f"{name}/psi_{n}.c128"
            save(os.path.join(OUT, fname), psi)
            checkpoints.append(dict(step=n, t=t, psi=fname, **scalars))
        manifest["cases"].append(dict(
            name=name, description=description, ndim=grid.ndim, n=list(grid.n), length=list(grid.length),
            dx=list(grid.dx), dt=dt, nonlinearity=g, absorber=absorber, detectors=detectors,
            files=files, checkpoints=checkpoints,
            tolerance=dict(float64=1e-10, float32=float(max(1e-5, 10 * err32))),
            float32_error_measured=err32))
        print(f"{name:28s} norm after {CHECKPOINTS[-1]} steps = {checkpoints[-1]['norm']:.6f}   "
              f"float32 rel. L2 error = {err32:.1e}")
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    size = sum(os.path.getsize(os.path.join(dp, fn)) for dp, _, fns in os.walk(OUT) for fn in fns)
    print(f"wrote {OUT}/ ({size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
