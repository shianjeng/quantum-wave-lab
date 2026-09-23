# Porting to C++ / WebAssembly

**[English](#english)** | **[日本語](#japanese)** · [← README](../README.md)

<a id="english"></a>
## English

### Reference data

`python export_reference.py` writes [`reference/`](../reference/README.md): five fixed cases (free packet, harmonic oscillator, sech² barrier with absorber and flux detectors, Gross–Pitaevskii soliton, 2D double slit with absorber). Each case stores the inputs as arrays and ψ, norm, energy and flux after 1, 100 and 1000 steps, together with the tolerances a float64 and a float32 port must meet. The binary layout and every convention a port must match (coordinates, FFT ordering, order of the factors, absorber, flux) are described in [reference/README.md](../reference/README.md). `tests/test_reference.py` checks that the Python implementation still reproduces the files. The reference data belong to qwave 0.3.0. The JavaScript port in [`web/`](../web/README.md) already passes them (`cd web && node --test`): float64 to ~10⁻¹³, float32 to ~10⁻⁵.

### Notes from the validation

- **Precision**: single precision is sufficient (relative error ~10⁻⁵ after 1000 steps), but the norm drifts slowly (~10⁻⁴ per 1000 steps). Renormalize every few hundred steps when no absorber is active (`renormalize_every`).
- **Box vs. view**: simulate on a box larger than the visible area, so that packets start clear of the absorbing layer (as `make_gif.py` does: 56 × 56 simulated, 40 × 40 shown).
- **Hard walls**: keep V·Δt ≲ 1 where the wave is, e.g. Δt = 0.005 for V = 200.
- **Absorber**: scale γ_max with the typical k and the layer width, γ_max ≳ 5k / width ([method](method.md#absorbing-boundary)).
- **Live transmission**: the 4th-order finite-difference flux needs only 5 grid lines (error ≈ (kΔx)⁴/30), but it needs smooth potentials. For drawn walls, prefer soft brush edges.
- **Nonlinear mode** (g ≠ 0): keep Δt·k_max²/2 < π, otherwise the solution blows up.
- **Slits**: build openings symmetric about their centres ((2·centre + width)/Δy odd), as `qwave.potentials.slits` does.
- Run the checks of `diagnose` every few frames and show them to the user.

---

<a id="japanese"></a>
## 日本語

### リファレンスデータ

`python export_reference.py` は [`reference/`](../reference/README.md) に 5 つの固定ケースを出力します（自由粒子、調和振動子、吸収層と確率流検出器つきの sech² 障壁、グロス・ピタエフスキーのソリトン、吸収層つきの 2D 二重スリット）。各ケースには入力の配列と、1・100・1000 ステップ後の ψ・ノルム・エネルギー・確率流、そして float64 版と float32 版が満たすべき許容誤差が入っています。バイナリ形式と、移植版が合わせるべき約束事（座標、FFT の並び、因子の順序、吸収層、確率流）は [reference/README.md](../reference/README.md) にまとめています。Python 実装がこれを再現し続けていることは `tests/test_reference.py` で確認しています。このリファレンスデータは qwave 0.3.0 に対応します。[`web/`](../web/README.md#japanese) の JavaScript 版はすでにこれに合格しています（`cd web && node --test`、float64 で ~10⁻¹³、float32 で ~10⁻⁵）。

### 検証から得た注意点

- **精度**: 単精度で十分（1000 ステップ後の相対誤差 ~10⁻⁵）ですが、ノルムがゆっくりずれます（1000 ステップで ~10⁻⁴）。吸収層がないときは数百ステップごとに再規格化します（`renormalize_every`）。
- **計算領域と表示領域**: 計算領域は表示領域より大きく取り、波包の初期位置が吸収層にかからないようにします（`make_gif.py` では 56 × 56 を計算して 40 × 40 を表示）。
- **硬い壁**: 波のある場所で V·Δt ≲ 1 を保ちます（例: V = 200 なら Δt = 0.005）。
- **吸収層**: γ_max は典型的な k と層の幅に合わせ、γ_max ≳ 5k / 幅 とします（[手法](method.md#吸収境界)）。
- **透過率のリアルタイム表示**: 4 次差分の確率流は 5 本の格子線だけで済みます（誤差 ≈ (kΔx)⁴/30）が、滑らかなポテンシャルが前提です。描いた壁には柔らかいブラシの縁を使うのがおすすめです。
- **非線形モード**（g ≠ 0）: Δt·k_max²/2 < π を守らないと解が発散します。
- **スリット**: `qwave.potentials.slits` と同様に、開口を中心について対称に作ります（(2·中心 + 幅)/Δy が奇数）。
- `diagnose` のチェックを数フレームごとに実行し、ユーザーに表示します。
