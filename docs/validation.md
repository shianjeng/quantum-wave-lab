# Validation

**[English](#english)** | **[日本語](#japanese)** · [← README](../README.md)

![validation](../figures/validation.png)

<a id="english"></a>
## English

`python validate.py` compares the solver with analytic results and writes the figure above; CI runs it on every push. `--quick` runs every check without the figure (~20 s on 4 cores), `--jobs N` sets the number of worker processes, and `--update-docs` regenerates the tables below. The tables were produced with NumPy 2, so the single-precision rows include a float32 FFT.

<!-- validation-table:en:start -->
| Test | Error |
|---|---|
| 2D norm conservation (double slit) | 2 × 10⁻¹³ |
| Free Gaussian spreading σ(t) (relative) | 2 × 10⁻¹⁴ |
| Harmonic-oscillator coherent state ⟨x⟩ = x₀ cos ωt | 6 × 10⁻⁵ |
| Harmonic-oscillator energy conservation (relative) | 6 × 10⁻⁶ |
| Rectangular-barrier transmission T(E) (analytic, averaged over the packet's momentum distribution) | 2 × 10⁻³ |
| Double-slit fringe period Δk_y = 2π/d (relative; 14 fringe zeros) | 3 × 10⁻⁴ |
| Double slit: mirror symmetry of the two single-slit far fields P₊(k_y) = P₋(−k_y) (relative) | 3 × 10⁻¹³ |
| Bloch-oscillation period T_B = 2π/(F a) (relative; 3 periods) | 2 × 10⁻⁴ |
| Imaginary-time eigenstates: harmonic oscillator E₀…E4 = n + 1/2 | 4 × 10⁻¹⁰ |
| Time-step convergence order, smooth V: \|p − 2\| | 3 × 10⁻⁴ |
| Time-step convergence order, time-dependent V: \|p − 2\| | 3 × 10⁻⁴ |
| Time-step convergence order, sharp barrier edges: \|p − 1\| (drops to ≈ first order) | 1 × 10⁻¹ |
| Driven oscillator (time-dependent V) ⟨x⟩(t) | 1 × 10⁻⁵ |
| Single precision (complex64) vs double, 2D, 1000 steps (relative L2) | 5 × 10⁻⁵ |
| Single precision (complex64) norm drift, 1000 steps | 8 × 10⁻⁵ |
| Single precision with renormalize_every=100: max norm drift over 1000 steps | 8 × 10⁻⁶ |
| Absorbing layer, residual (reflected + wrapped) probability for 1 ≤ k ≤ 8 | 2 × 10⁻⁶ |
| Flux detector: transmission through a sech² barrier vs exact T(E), absorber active | 4 × 10⁻⁵ |
| Flux detector: 4th-order finite-difference vs spectral current | 9 × 10⁻⁷ |
| Bright soliton (Gross–Pitaevskii, g < 0) vs exact solution at t = 20 (L2) | 7 × 10⁻⁶ |
| Bright soliton: convergence order in dt: \|p − 2\| | 9 × 10⁻⁶ |
| Bright soliton: energy vs v²/2 − g²/24 (relative) | 2 × 10⁻¹² |
| Gross–Pitaevskii ground state, μ vs Thomas–Fermi limit (relative; g = 500, μ ≈ 41ω) | 2 × 10⁻⁴ |
| Gross–Pitaevskii ground state is stationary: 1 − \|⟨ψ(0)\|ψ(10)⟩\| and phase error of e^{−iμt} | 4 × 10⁻⁷ |
<!-- validation-table:en:end -->

In the 2D demo (`demo_2d.py`), the tunneling probability through a wall that depends only on x is compared with
(a) a 1D run on the same Δx and Δt: the difference is 2 × 10⁻⁵, confirming the 2D implementation, and
(b) the analytic result: the difference is 3.7 %, the discretization error of the coarse real-time grid (Δx = 0.125).
The web version will manage this **accuracy vs. real-time trade-off** through the grid size.

---

<a id="japanese"></a>
## 日本語

`python validate.py` で解析解と比較し、上の図を出力します（CI で毎 push 実行）。`--quick` は図を作らずにすべてのチェックだけを行い（4 コアで約 20 秒）、`--jobs N` で並列プロセス数を指定し、`--update-docs` で下の表を再生成します。表は NumPy 2 で作成しているため、単精度の行は FFT も float32 です。

<!-- validation-table:ja:start -->
| テスト | 誤差 |
|---|---|
| 2D ノルム保存（二重スリット） | 2 × 10⁻¹³ |
| 自由粒子の波束の広がり σ(t)（相対） | 2 × 10⁻¹⁴ |
| 調和振動子コヒーレント状態 ⟨x⟩ = x₀ cos ωt | 6 × 10⁻⁵ |
| 調和振動子のエネルギー保存（相対） | 6 × 10⁻⁶ |
| 矩形障壁の透過率 T(E)（運動量分布で平均した解析解と比較） | 2 × 10⁻³ |
| 二重スリットの縞の周期 Δk_y = 2π/d（相対、縞のゼロ点 14 個） | 3 × 10⁻⁴ |
| 二重スリット: 片方ずつの遠方パターンの鏡映対称性 P₊(k_y) = P₋(−k_y)（相対） | 3 × 10⁻¹³ |
| ブロッホ振動の周期 T_B = 2π/(F a)（相対、3 周期） | 2 × 10⁻⁴ |
| 虚時間発展による固有状態: 調和振動子 E₀…E4 = n + 1/2 | 4 × 10⁻¹⁰ |
| 時間刻みの収束次数（滑らかな V）: \|p − 2\| | 3 × 10⁻⁴ |
| 時間刻みの収束次数（時間依存 V）: \|p − 2\| | 3 × 10⁻⁴ |
| 時間刻みの収束次数（障壁の角が鋭い場合）: \|p − 1\|（ほぼ 1 次に落ちる） | 1 × 10⁻¹ |
| 強制振動子（時間依存 V）の ⟨x⟩(t) | 1 × 10⁻⁵ |
| 単精度（complex64）と倍精度の差、2D・1000 ステップ（相対 L2） | 5 × 10⁻⁵ |
| 単精度（complex64）のノルムのずれ、1000 ステップ | 8 × 10⁻⁵ |
| 単精度＋100 ステップごとの再規格化: 1000 ステップ中のノルムの最大のずれ | 8 × 10⁻⁶ |
| 吸収層の残存確率（反射＋折り返し）、1 ≤ k ≤ 8 | 2 × 10⁻⁶ |
| フラックス検出器: sech² 障壁の透過率と厳密解（吸収層あり） | 4 × 10⁻⁵ |
| フラックス検出器: 4 次差分と スペクトル法の電流の差 | 9 × 10⁻⁷ |
| 明るいソリトン（GP 方程式, g < 0）と厳密解の差、t = 20（L2） | 7 × 10⁻⁶ |
| 明るいソリトン: 時間刻みの収束次数 \|p − 2\| | 9 × 10⁻⁶ |
| 明るいソリトン: エネルギーと v²/2 − g²/24 の差（相対） | 2 × 10⁻¹² |
| GP 方程式の基底状態、μ とトーマス・フェルミ極限の差（相対、g = 500, μ ≈ 41ω） | 2 × 10⁻⁴ |
| GP 基底状態の定常性: 1 − \|⟨ψ(0)\|ψ(10)⟩\| と位相 e^{−iμt} のずれ | 4 × 10⁻⁷ |
<!-- validation-table:ja:end -->

2D デモ（`demo_2d.py`）では、x のみに依存する壁でのトンネル確率を
(a) 同じ Δx・Δt の 1D 計算（差 2 × 10⁻⁵：2D 実装の正しさを確認）と、
(b) 解析解（差 3.7%：リアルタイム用の粗い格子 Δx = 0.125 による離散化誤差）の両方と比較しています。
Web 版では、この**精度とリアルタイム性のトレードオフ**を格子サイズで調整します。
