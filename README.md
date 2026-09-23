# Quantum Wave Lab

![CI](https://github.com/shianjeng/quantum-wave-lab/actions/workflows/ci.yml/badge.svg)

**[English](#english)** | **[日本語](#japanese)**

<p align="center">
  <img src="figures/hero.gif" alt="Tunneling, double-slit interference and scattering of a quantum wave packet" width="100%">
</p>

---

<a id="english"></a>
## English

An interactive simulator for the **time-dependent Schrödinger equation**, solved with the **split-step Fourier method (Strang splitting)**. The goal is a browser app where you draw potentials with the mouse and watch tunneling, double-slit interference and scattering in real time.

**Status:** the Python reference implementation and its numerical validation are complete. Next, the solver will be ported to C++ → WebAssembly + WebGL.

<sub>GIF: brightness ∝ |ψ| (fixed scale per panel). Δt = 0.005, Δx = 0.125. Relative density error vs. a 2× finer grid with Δt = 0.001: 1–5 %. Tunneling panel: V₀ = 10 > ⟨E⟩ ≈ 8; over-barrier (E > V₀) components contribute < 1 % of the transmitted probability.</sub>

### Method

Units ħ = m = 1. One time step:

```
ψ ← e^{-iV Δt/2} ψ  →  FFT  →  ψ ← e^{-ik² Δt/2} ψ  →  IFFT  →  ψ ← e^{-iV Δt/2} ψ
```

Every factor is unitary, so the norm is conserved to machine precision; the time error is O(Δt²) for smooth potentials. Wrap-around at the periodic boundary is suppressed by a complex absorbing potential e^{-γ(r)Δt}, whose absorption rate is independent of Δt.
Barrier and slit widths are integer multiples of the grid spacing Δx, so the effective width never silently differs from the nominal one (`check_on_grid`).

### Validation

`python validate.py` compares the solver with analytic results (run by CI on every push).

| Test | Error |
|---|---|
| 2D norm conservation (double slit) | 4 × 10⁻¹⁴ |
| Free Gaussian spreading σ(t) | 2 × 10⁻¹⁴ |
| Harmonic-oscillator coherent state ⟨x⟩ = x₀ cos ωt | 6 × 10⁻⁵ |
| Harmonic-oscillator energy conservation (relative) | 6 × 10⁻⁶ |
| Rectangular-barrier transmission T(E) (analytic, averaged over the packet's momentum distribution) | 9 × 10⁻⁴ |

In the 2D demo (`demo_2d.py`), the tunneling probability through a wall that depends only on x is compared with
(a) a 1D run on the same Δx and Δt — difference 3 × 10⁻⁵, confirming the 2D implementation, and
(b) the analytic result — difference 3.7 %, the discretization error of the coarse real-time grid (Δx = 0.125).
The web version will manage this **accuracy vs. real-time trade-off** through the grid size.

![validation](figures/validation.png)

### Usage

```bash
pip install -r requirements.txt
python -m pytest -q      # fast unit tests
python validate.py       # comparison with analytic results → figures/validation.png
python demo_2d.py        # double slit & 2D tunneling → figures/*.png, *.gif
python make_gif.py       # README GIF → figures/hero.gif
```

### Structure

```
qwave/solver.py      Grid, Solver, gaussian_packet, absorbing_mask
qwave/potentials.py  harmonic, rectangular barrier, double slit, brush painting, grid check
qwave/analytic.py    analytic references (barrier transmission, ...)
validate.py          validation against analytic results (writes figures)
demo_2d.py           2D demos
make_gif.py          README GIF generator
tests/               unit tests for CI
```

### Roadmap

- [x] Python reference implementation and numerical validation
- [ ] Port to C++17 (FFT: pocketfft / KissFFT) with tests that match the Python results
- [ ] Compile to WebAssembly with Emscripten, render with WebGL2
- [ ] Draw potentials with the mouse; presets (tunneling, double slit, scattering)
- [ ] Benchmark: plain JS vs. WASM vs. WASM + SIMD
- [ ] Deploy on GitHub Pages

---

<a id="japanese"></a>
## 日本語

**時間依存シュレディンガー方程式**を**分割ステップ・フーリエ法（Strang 分解）**で解くシミュレータです。マウスでポテンシャルを描き、トンネル効果・二重スリット干渉・散乱をブラウザ上でリアルタイムに観察できる Web アプリを目指しています。

**現状:** Python によるリファレンス実装と数値検証が完成した段階です。今後 C++ → WebAssembly + WebGL へ移植します。

<sub>GIF：明るさ ∝ |ψ|（各パネル固定スケール）。Δt = 0.005, Δx = 0.125。2 倍細かい格子・Δt = 0.001 の計算との密度の相対誤差は 1–5%。トンネル側は V₀ = 10 > ⟨E⟩ ≈ 8 で、透過確率のうち障壁より高いエネルギー成分の寄与は 1% 未満。</sub>

### 手法

単位系 ħ = m = 1。1 ステップは

```
ψ ← e^{-iV Δt/2} ψ  →  FFT  →  ψ ← e^{-ik² Δt/2} ψ  →  IFFT  →  ψ ← e^{-iV Δt/2} ψ
```

各演算子はユニタリなのでノルムは機械精度で保存され、滑らかなポテンシャルでは時間誤差 O(Δt²)。周期境界での折り返しは複素吸収ポテンシャル e^{-γ(r)Δt} で抑制します（吸収率は Δt に依存しません）。
障壁やスリットの幅は格子間隔 Δx の整数倍に取り、実効幅が公称値とずれないようにしています（`check_on_grid`）。

### 数値検証

`python validate.py` で以下を解析解と比較します（CI で毎 push 実行）。

| テスト | 誤差 |
|---|---|
| 2D ノルム保存（二重スリット） | 4 × 10⁻¹⁴ |
| 自由粒子の波束の広がり σ(t) | 2 × 10⁻¹⁴ |
| 調和振動子コヒーレント状態 ⟨x⟩ = x₀ cos ωt | 6 × 10⁻⁵ |
| 調和振動子のエネルギー保存（相対） | 6 × 10⁻⁶ |
| 矩形障壁の透過率 T(E)（運動量分布で平均した解析解と比較） | 9 × 10⁻⁴ |

2D デモ（`demo_2d.py`）では、x のみに依存する壁でのトンネル確率を
(a) 同じ Δx・Δt の 1D 計算（差 3 × 10⁻⁵：2D 実装の正しさを確認）と、
(b) 解析解（差 3.7%：リアルタイム用の粗い格子 Δx = 0.125 による離散化誤差）の両方と比較しています。
Web 版では、この**精度とリアルタイム性のトレードオフ**を格子サイズで調整します。

### 使い方

```bash
pip install -r requirements.txt
python -m pytest -q      # 高速ユニットテスト
python validate.py       # 解析解との比較 → figures/validation.png
python demo_2d.py        # 二重スリット・2D トンネル → figures/*.png, *.gif
python make_gif.py       # README 冒頭の GIF → figures/hero.gif
```

### 構成

```
qwave/solver.py      Grid, Solver, gaussian_packet, absorbing_mask
qwave/potentials.py  調和ポテンシャル, 矩形障壁, 二重スリット, ブラシ描画, 格子整合チェック
qwave/analytic.py    解析解（矩形障壁の透過率など）
validate.py          解析解との比較（図を出力）
demo_2d.py           2D デモ
make_gif.py          README 用 GIF の生成
tests/               CI 用ユニットテスト
```

### ロードマップ

- [x] Python リファレンス実装と数値検証
- [ ] C++17 移植（FFT: pocketfft / KissFFT）、Python と結果を突き合わせるテスト
- [ ] Emscripten で WebAssembly 化、WebGL2 で描画
- [ ] マウスでポテンシャルを描画、プリセット（トンネル・二重スリット・散乱）
- [ ] 性能比較: 純 JS vs WASM vs WASM + SIMD
- [ ] GitHub Pages で公開

---

## License

MIT © 2026 Hank Zhu
