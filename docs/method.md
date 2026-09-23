# Method

**[English](#english)** | **[日本語](#japanese)** · [← README](../README.md)

<a id="english"></a>
## English

### Time step

Units ħ = m = 1. The equation is i ∂ψ/∂t = [−∇²/2 + V(r, t) + g|ψ|²] ψ (g = 0: Schrödinger, g ≠ 0: Gross–Pitaevskii). One step of the split-step Fourier method (Strang splitting):

```
ψ ← e^{-i(V + g|ψ|²) Δt/2} ψ  →  FFT  →  ψ ← e^{-ik² Δt/2} ψ  →  IFFT  →  ψ ← e^{-i(V + g|ψ|²) Δt/2} ψ
```

- Every factor is unitary, so the norm is conserved to machine precision; the nonlinear half-steps are exact because they do not change |ψ|.
- The time error is O(Δt²) for smooth (also time-dependent) potentials, but only ≈ O(Δt) at sharp edges such as a rectangular barrier. `rect_barrier(..., edge=w)` rounds the edges with tanh and restores O(Δt²).
- For a static linear potential the closing half-step of one step and the opening half-step of the next are merged into one factor e^{−iVΔt}. `step(..., callback=f)` still hands the true ψ to `f` after every step.
- A time-dependent potential is passed as a callable `V(t)`, evaluated at the midpoint of each step.
- `Solver(..., dtype=np.complex64)` emulates the float32 arithmetic of the WASM/WebGL port (relative error ~10⁻⁵ after 1000 steps). Its norm drifts slowly; `renormalize_every=N` removes that when no absorber is active.

### Absorbing boundary

Wrap-around at the periodic boundary is suppressed by a complex absorbing potential, applied as the mask e^{−γ(r)Δt} with a smooth sin² ramp; its absorption per unit time does not depend on Δt. How well it works depends on the wavenumber k and the layer width: slow waves are reflected by a steep ramp, fast waves cross a thin or weak layer and wrap around.

| layer width | γ_max | k = 0.5 | k = 1 | k = 2 | k = 4 | k = 8 |
|---|---|---|---|---|---|---|
| 3.2 (8 % of a 40-wide box) | 2 | 3 × 10⁻¹ | 6 × 10⁻² | 3 × 10⁻³ | 4 × 10⁻² | 2 × 10⁻¹ |
| 3.2 | 10 | 5 × 10⁻¹ | 2 × 10⁻¹ | 2 × 10⁻² | 3 × 10⁻⁵ | 4 × 10⁻⁴ |
| 32 (8 % of a 400-wide box) | 2 | 6 × 10⁻³ | 2 × 10⁻⁶ | 6 × 10⁻⁸ | 7 × 10⁻⁸ | 4 × 10⁻⁷ |

(residual = reflected + wrapped probability, `absorber_residual`). Rule of thumb: the probability that crosses both edge layers is ≈ exp(−2 γ_max · width / k), so choose γ_max ≳ 5k / width.

### Geometry on the grid

- Barrier and wall widths are integer multiples of Δx, so the effective width never silently differs from the nominal one (`check_on_grid` raises otherwise). Such a wall is shifted by half a cell from its nominal position, a harmless translation.
- Slit openings are the w/Δy cells placed symmetrically about each centre, so width, separation and mirror symmetry are all exact. This needs the edges on cell faces: (2·centre + width)/Δy must be an odd integer, e.g. width 1.125 and separation 5 on Δy = 0.125. Otherwise `slits`/`double_slit` raise an error that suggests the nearest widths that work.

### Stationary states

`eigenstates(grid, V, n)` propagates in imaginary time (t → −iτ) and projects out the states already found (Gram–Schmidt); the ground state starts from a nodeless envelope. With `nonlinearity=g`, `ground_state` returns the Gross–Pitaevskii ground state and its chemical potential μ. Both nonlinear half-steps use the density at the start of the step, which keeps the fixed point symmetric (bias O(dτ²)).

### Gross–Pitaevskii equation

With g ≠ 0 the split-step Fourier method is only stable for Δt·k_max²/2 < π, where k_max² is summed over the axes (Weideman & Herbst 1986). Beyond that, grid modes whose kinetic phase per step is a multiple of 2π are pumped by the nonlinearity, and after a while the solution blows up (the energy jumped from 8.5 to 411 in one test). `Solver` warns and `diagnose` reports it.

### Live measurements

`FluxDetector` integrates the probability current j = Im(ψ*∇ψ) through a line (trapezoid rule over the per-step samples), optionally over a y-range such as one slit. This gives the transmission while the wave is still scattering, unaffected by the absorber. The derivative is spectral, or a 4th-order finite difference over 5 grid lines (relative error ≈ (kΔx)⁴/30), which is the cheap option for the browser. The flux needs a smooth solution. At sharp potential edges the splitting error puts a little probability into high-k modes, the flux then oscillates faster than the time step, and its time integral aliases (~20 % error for a rectangular barrier at Δt = 0.02). Use smooth potentials such as `sech2_barrier`.

### Sanity checks

`Solver.diagnose(ψ)` reports:
- |V + g|ψ|²|·Δt > 1 where ψ is (hard walls need Δt ≲ 1/V),
- momentum content beyond 2/3 of the Nyquist limit, or ψ cut off at the box edge,
- an initial state that already sits in the absorbing layer,
- the nonlinear stability limit above.

---

<a id="japanese"></a>
## 日本語

### 時間発展

単位系 ħ = m = 1。方程式は i ∂ψ/∂t = [−∇²/2 + V(r, t) + g|ψ|²] ψ（g = 0 でシュレディンガー方程式、g ≠ 0 でグロス・ピタエフスキー方程式）。分割ステップ・フーリエ法（Strang 分解）の 1 ステップは

```
ψ ← e^{-i(V + g|ψ|²) Δt/2} ψ  →  FFT  →  ψ ← e^{-ik² Δt/2} ψ  →  IFFT  →  ψ ← e^{-i(V + g|ψ|²) Δt/2} ψ
```

- 各演算子はユニタリなのでノルムは機械精度で保存されます。非線形の半ステップも |ψ| を変えないので厳密です。
- 時間誤差は滑らかなポテンシャル（時間依存でも）では O(Δt²) ですが、矩形障壁のような鋭い角があるとほぼ O(Δt) に落ちます。`rect_barrier(..., edge=w)` で角を tanh で丸めると O(Δt²) に戻ります。
- 時間に依存しない線形ポテンシャルでは、あるステップの最後の半ステップと次のステップの最初の半ステップを 1 つの因子 e^{−iVΔt} にまとめています。それでも `step(..., callback=f)` の `f` には毎ステップ正しい ψ が渡されます。
- 時間依存ポテンシャルは関数 `V(t)` として渡し、各ステップの中点で評価します。
- `Solver(..., dtype=np.complex64)` で WASM/WebGL 版の float32 演算を再現します（1000 ステップ後の相対誤差 ~10⁻⁵）。ノルムがゆっくりずれるので、吸収層がないときは `renormalize_every=N` で補正します。

### 吸収境界

周期境界での折り返しは、なめらかな sin² 型の複素吸収ポテンシャル（マスク e^{−γ(r)Δt}）で抑えます。単位時間あたりの吸収率は Δt に依存しません。効き具合は波数 k と層の幅で決まり、遅い波は急な立ち上がりで反射され、速い波は薄い・弱い層を通り抜けて反対側から戻ってきます。

| 層の幅 | γ_max | k = 0.5 | k = 1 | k = 2 | k = 4 | k = 8 |
|---|---|---|---|---|---|---|
| 3.2（幅 40 の箱の 8%） | 2 | 3 × 10⁻¹ | 6 × 10⁻² | 3 × 10⁻³ | 4 × 10⁻² | 2 × 10⁻¹ |
| 3.2 | 10 | 5 × 10⁻¹ | 2 × 10⁻¹ | 2 × 10⁻² | 3 × 10⁻⁵ | 4 × 10⁻⁴ |
| 32（幅 400 の箱の 8%） | 2 | 6 × 10⁻³ | 2 × 10⁻⁶ | 6 × 10⁻⁸ | 7 × 10⁻⁸ | 4 × 10⁻⁷ |

（残存確率 = 反射＋折り返し、`absorber_residual`）。目安として、両端の層を通り抜ける確率は ≈ exp(−2 γ_max · 幅 / k) なので、γ_max ≳ 5k / 幅 に取ります。

### 格子上の形状

- 障壁や壁の幅は Δx の整数倍に取り、実効幅が公称値とずれないようにしています（ずれる場合は `check_on_grid` がエラー）。こうした壁は公称位置から半セル平行移動しますが、物理的には影響ありません。
- スリットの開口部は、各中心について対称に並ぶ w/Δy 個のセルなので、幅・間隔・鏡映対称性がすべて厳密です。そのためには開口の端がセルの境界に来る必要があり、(2·中心 + 幅)/Δy が奇数でなければなりません（例: Δy = 0.125 で幅 1.125、間隔 5）。満たさない場合、`slits`/`double_slit` は使える幅を提案してエラーにします。

### 定常状態

`eigenstates(grid, V, n)` は虚時間（t → −iτ）で発展させ、すでに求めた状態を射影で除きながら（グラム・シュミット）順に求めます。基底状態は節のない初期関数から始めます。`nonlinearity=g` を指定すると `ground_state` はグロス・ピタエフスキー方程式の基底状態と化学ポテンシャル μ を返します。非線形の半ステップは 2 回ともステップ開始時の密度を使うので、固定点は対称になります（誤差 O(dτ²)）。

### グロス・ピタエフスキー方程式

g ≠ 0 のとき、分割ステップ・フーリエ法は Δt·k_max²/2 < π（k_max² は各軸の和）でしか安定ではありません（Weideman & Herbst 1986）。これを超えると、1 ステップあたりの運動エネルギーの位相が 2π の倍数になる格子モードが非線形項で増幅され、しばらくすると解が発散します（あるテストではエネルギーが 8.5 から 411 に跳ね上がりました）。`Solver` が警告し、`diagnose` も報告します。

### リアルタイム計測

`FluxDetector` は、ある線を横切る確率流 j = Im(ψ*∇ψ) を時間積分します（各ステップの値を台形公式で積分）。1 本のスリットなど、y の範囲を限定することもできます。散乱の途中でも透過率が分かり、吸収層の影響も受けません。微分はスペクトル法か、5 本の格子線だけを使う 4 次差分（相対誤差 ≈ (kΔx)⁴/30）で、後者はブラウザ向けの軽い方法です。ただし確率流には滑らかな解が必要です。ポテンシャルの角が鋭いと、分割誤差で高波数成分にわずかな確率が乗り、確率流が時間刻みより速く振動して時間積分がエイリアスします（Δt = 0.02 の矩形障壁で約 20% の誤差）。`sech2_barrier` などの滑らかなポテンシャルを使ってください。

### 妥当性チェック

`Solver.diagnose(ψ)` は次のことを報告します。
- ψ のある場所で |V + g|ψ|²|·Δt > 1 になっている（硬い壁では Δt ≲ 1/V が必要）
- ナイキスト限界の 2/3 を超える運動量成分がある、または ψ が箱の端で切れている
- 初期状態がすでに吸収層にかかっている
- 上記の非線形安定条件を満たしていない
