"""The golden files in reference/ (written by export_reference.py) must still be reproduced.
A C++/WASM port runs the same comparison against the same files."""
import json
import os

import numpy as np
import pytest

from export_reference import rel_l2, run_case
from qwave.solver import Grid, Solver

ROOT = os.path.join(os.path.dirname(__file__), "..", "reference")
with open(os.path.join(ROOT, "manifest.json")) as f:
    MANIFEST = json.load(f)


def load(case, key, complex_=False):
    arr = np.fromfile(os.path.join(ROOT, case["files"][key] if key in case["files"] else key),
                      dtype="<c16" if complex_ else "<f8")
    return arr.reshape(case["n"])


@pytest.mark.parametrize("case", MANIFEST["cases"], ids=lambda c: c["name"])
@pytest.mark.parametrize("dtype", [np.complex128, np.complex64], ids=["float64", "float32"])
def test_reference_case(case, dtype):
    grid = Grid(case["n"], case["length"])
    V, psi0 = load(case, "V"), load(case, "psi0", complex_=True)
    if case["absorber"] is not None:
        mask = Solver(grid, V, case["dt"], absorber=case["absorber"]).absorber
        assert np.max(np.abs(mask - load(case, "absorber_mask"))) < 1e-14
    _, results = run_case(grid, V, psi0, case["dt"], case["absorber"], case["nonlinearity"], case["detectors"],
                          dtype=dtype)
    tol = case["tolerance"]["float64" if dtype == np.complex128 else "float32"]
    for (step, t, psi, scalars), ref in zip(results, case["checkpoints"]):
        assert step == ref["step"] and abs(t - ref["t"]) < 1e-9
        assert rel_l2(psi, np.fromfile(os.path.join(ROOT, ref["psi"]), dtype="<c16").reshape(case["n"])) < tol
        assert abs(scalars["norm"] - ref["norm"]) < tol
        assert abs(scalars["energy"] - ref["energy"]) < tol * max(1.0, abs(ref["energy"]))
        for got, want in zip(scalars.get("flux_transmitted", []), ref.get("flux_transmitted", [])):
            assert abs(got - want) < tol
