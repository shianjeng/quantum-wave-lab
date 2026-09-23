# Reference data for ports

Golden files written by `python export_reference.py` (run with NumPy ≥ 2). A C++ / WASM / JS port of
the solver should reproduce them; `tests/test_reference.py` does the same comparison for the Python
implementation. Units ħ = m = 1.

## Files

`manifest.json` lists the cases. Per case:

| key | meaning |
|---|---|
| `n`, `length`, `dx` | grid points, box length and spacing per axis (axis 0 = x, axis 1 = y) |
| `dt`, `nonlinearity` | time step, g of the term g\|ψ\|²ψ (0 = linear) |
| `absorber` | `null` or the `absorbing_mask` parameters; the mask itself is in `files.absorber_mask` |
| `files.V`, `files.psi0` | potential and initial wave function |
| `checkpoints` | after `step` = 1, 100, 1000 steps: time `t`, file `psi`, `norm`, `energy`, `flux_transmitted` |
| `detectors` | flux detectors (see below), one `flux_transmitted` entry each |
| `tolerance` | allowed relative L2 error of ψ (and absolute error of the scalars) for a float64 / float32 port |

Binary layout: raw little-endian, no header, C order (row-major) with **axis 0 slowest**, i.e. element
`[ix, iy]` of a 2D array sits at offset `ix * n[1] + iy`.
`*.f64` are float64 arrays; `*.c128` are complex128 arrays stored as interleaved `(re, im)` float64 pairs.

## Conventions a port must match

- Coordinates: `x_j = -L/2 + j·dx`, `j = 0 … n-1` (periodic box).
- Wavenumbers in FFT order: `k_j = 2π·j/L` for `j < n/2`, `2π·(j - n)/L` otherwise (`numpy.fft.fftfreq`);
  `k² = k_x² + k_y²`. The FFT is unnormalized forward, `1/N` on the inverse.
- One step, linear case:
  `ψ ← P·ψ;  ψ ← IFFT[exp(-i k² dt/2)·FFT[ψ]];  ψ ← A·P·ψ`, with `P = exp(-i V dt/2)` and `A` the
  absorber mask (1 without absorber). The last `P` of a step and the first `P` of the next can be merged
  into `exp(-i V dt)`; the result is identical up to rounding.
- Nonlinear case (g ≠ 0): each `P` is `exp(-i (V + g|ψ|²) dt/2)` evaluated on the ψ it multiplies
  (exact, since it does not change |ψ|); no merging. `A` is applied after the second `P`.
- `norm = Σ|ψ|²·dV`, `dV = Π dx`.
- `energy = [½ Σ k²|FFT[ψ]|² dV / N + Σ V|ψ|² dV + ½ g Σ|ψ|⁴ dV] / norm` (N = total number of points).
- Flux detector at the grid line `x = position` (nearest grid point), optionally restricted to
  `span[0] <= y < span[1]`: `Φ = Σ Im(ψ* ∂ψ/∂x)·dy` along that line. `"spectral"` takes ∂ψ/∂x as
  `IFFT_x[i k_x FFT_x[ψ]]`; `"fd4"` as `[ψ(x-2dx) - 8ψ(x-dx) + 8ψ(x+dx) - ψ(x+2dx)] / (12 dx)`
  (periodic indices). `flux_transmitted` is `∫Φ dt` by the trapezoid rule over samples at `t = 0` and
  after every step.

## Tolerances

`float64` is 1e-10 (different FFT libraries agree to ~1e-13 here). `float32` is 10× the error that a
complex64 run with NumPy's single-precision FFT actually shows (`float32_error_measured`, 2–5 × 10⁻⁵),
so a correct single-precision port passes with margin while a wrong sign or a missing factor does not.
