# Quantum Wave Lab: browser version

**[English](#english)** | **[日本語](#japanese)** · [← README](../README.md)

<a id="english"></a>
## English

A static site with no build step: `index.html`, `style.css` and six plain scripts. Open `index.html` directly, or serve the folder:

```bash
python -m http.server --directory web 8000     # then open http://localhost:8000
```

GitHub Pages publishes it with `.github/workflows/pages.yml` (one-time setup: Settings → Pages → Source: *GitHub Actions*).

### What the page does

- **Layout:** the field is the centre of the page. Readouts sit as a HUD on it, the main actions in a dock below, and settings in tabs (Scenes, Draw, Wave, View, Data). The scene gallery shows real simulation thumbnails, computed in the page at 128² in idle time and cached in `localStorage`.
- **Guided tour:** 7 steps on the first visit (the ★ Tour button replays it): wave packet → tunneling → double slit → one particle at a time → momentum space → 3D → your turn.
- **One particle at a time:** detections appear on a screen at x = 15, one dot each. Every dot is drawn from the time-integrated probability current that has reached the screen (`FluxDetector.profile`, `sampleProfile`), and the interference fringes build up as in Tonomura's electron experiment. An orange histogram shows the pattern they converge to; for the double slit its maxima are at 0, ±5.8, ±13.4 (λD/d ≈ 5.65).
- **3D view (WebGL2):** |ψ| as a lit height field with the colours of the 2D view. Walls appear as ridges. Drag to rotate, scroll to zoom; it turns slowly when idle. It renders offscreen and is copied onto the visible canvas, so PNG export and recording work in 3D too.
- **Colour maps:** Inferno, Viridis and Ice. **Fullscreen** (F) for presentations.
- **10 scenes** with a short explanation each: tunneling, resonant tunneling, single slit, double slit, grating, scattering, harmonic trap, quantum corral, BEC interference, empty box.
- **Drawing:** draw and erase walls (Shift + drag for straight walls), with undo/redo. You can also launch a packet by dragging, and move the flux detector.
- **Views:** position (|ψ|, |ψ|², or phase as colour) or momentum (|φ(k)|², where diffraction orders and far-field fringes appear). There is optional auto brightness.
- **Potential lines:** contour lines of smooth potentials, plus the V = ⟨E⟩ line (the edge of the classically forbidden region).
- **Measurements:** time, energy ⟨E⟩, probability in the box, transmission through the detector (live, and unaffected by the absorber), and a chart that also shows what was absorbed.
- **Measure:** a Gaussian position measurement (σ = 0.8). The outcome is drawn from |ψ|², and ψ collapses around it while keeping its local momentum.
- **Sharing and export:** share links that encode the scene, parameters and drawn walls in the URL; save a PNG; record a video (WebM/MP4).
- **Checks, settings and languages:** checks from `diagnose()`; grid 128² / 256² / 512²; interaction g; English and Japanese.

### Files

```
js/qwave.js    solver: FFT, Grid, Solver (float32/float64, absorber, g|ψ|², merged half-steps), FluxDetector,
               diagnose, imaginary-time ground state; a port of the Python package qwave 0.3.0
js/scenes.js   scene definitions (walls, background, packet, detector, verified numbers for the texts)
js/share.js    share links: app state <-> URL-safe string
js/app.js      the page: drawing, launching, rendering, contour lines, particles, tour, gallery, export, checks
js/surface3d.js  the 3D view (WebGL2 height field, orbit camera)
js/i18n.js     UI strings (English / 日本語)
test/          node --test: reference data, scenes, share links, measurement / momentum space
```

### Tests

```bash
cd web && node --test
```

- `reference.test.mjs`: every case in [`reference/`](../reference/README.md) runs in float64 and float32. The float64 results match the Python golden files to ~10⁻¹³, the float32 results to ~10⁻⁵, which is 10× inside the stored tolerance.
- `scenes.test.mjs`: every scene builds on all three grids and starts without warnings, and the slit scenes are exactly mirror-symmetric. The transmissions quoted in the scene texts are re-measured as the app measures them, e.g. resonant tunneling: 38 % with both barriers, 25 % with one, 19 % detuned.
- `share.test.mjs`, `physics.test.mjs`: share links round-trip, and malformed links are rejected. The momentum density obeys Parseval and peaks at k₀. Measurement outcomes follow |ψ|², the post-measurement state is normalized, and ⟨k⟩ is kept.

CI runs them on every push (~40 s).

### Choices (from the Python validation, see [docs/porting.md](../docs/porting.md))

- **Box**: 51.2 × 51.2, periodic. Only the central 40 × 40 is shown; the rim is the absorbing layer (width 5.6, γ_max = 6), so packets start clear of it.
- **Grid**: 128², 256² (default, Δx = 0.2) or 512². One float32 step of 256² takes about 4 ms in Chrome on an M-series Mac, and 6 steps per frame run at ~40 fps.
- **Time step**: Δt = 0.01, so that V·Δt ≤ 1 for walls up to V = 100 (the wall-height slider goes to 200 to show the warning). With g > 0, Δt is lowered to 0.9 × 2π/Σk_max² for stability.
- **Walls**: the brush has a soft one-cell edge. Slit openings are placed symmetrically on the grid at every resolution.
- **Transmission**: a spectral flux detector, sampled once per frame. It matches the region sum to ~0.2 % for the double slit, whereas the 4th-order finite difference would be ~1 % low at kΔx = 0.8.
- **Checks**: `diagnose()` runs at every launch (strict) and every 60 frames. While the wave runs, high-k content only counts if it exceeds 1e-4 of the launched probability: hard walls always scatter a harmless ~1e-5, which relative thresholds would flag as the box empties. With less than 1 % left in the box, the checks and ⟨E⟩ pause.
- **V = ⟨E⟩ line**: drawn only if the drawn walls stay below 4⟨E⟩ (barriers, traps); around hard walls it would only trace the walls.
- **Share links**: stroke points are stored on a 0.05 lattice (finer than the finest grid), and the live drawing uses the same lattice, so a shared drawing is rasterized identically. When the grid resolution changes, the strokes are re-rasterized rather than resampled.
- **BEC scene**: the ground state is prepared by imaginary time with dτ = 0.04 (~330 steps, ~2 s); it differs from a dτ = 0.005 state by ~6 × 10⁻⁵.

In the browser console, `qwl.state` shows the app state and `qwl.tick(n)` advances n frames (also in a background tab).

---

<a id="japanese"></a>
## 日本語

ビルド不要の静的サイトです（`index.html`、`style.css` と 5 つのスクリプト）。`index.html` を直接開くか、フォルダをサーバーで配信します。

```bash
python -m http.server --directory web 8000     # http://localhost:8000 を開く
```

GitHub Pages へは `.github/workflows/pages.yml` で公開します（初回のみ Settings → Pages → Source を *GitHub Actions* に設定）。

### できること

- **レイアウト**: フィールドが主役です。計測値はフィールド上の HUD に、主な操作は下のドックに、設定はタブ（シーン・描画・波束・表示・データ）にまとめています。シーン一覧のサムネイルは、ページ内で 128² の格子を使って空き時間に計算した本物のシミュレーション結果で、`localStorage` にキャッシュします。
- **ガイドツアー**: 初回訪問時に 7 ステップ（★ ツアーでいつでも再生）: 波束 → トンネル効果 → 二重スリット → 粒子を 1 個ずつ → 運動量空間 → 3D → あなたの番。
- **粒子を 1 個ずつ**: x = 15 のスクリーンに、検出が 1 個ずつ点として現れます。点はスクリーンに届いた確率流の時間積分（`FluxDetector.profile`、`sampleProfile`）に従って選ばれ、外村彰らの電子の実験のように干渉縞が積み上がっていきます。オレンジのヒストグラムは点が収束していく先の分布です。二重スリットでは極大が 0、±5.8、±13.4 にあります（λD/d ≈ 5.65）。
- **3D 表示（WebGL2）**: |ψ| を光の当たった地形として、2D と同じ色で表示します。壁は尾根になります。ドラッグで回転、スクロールで拡大縮小し、操作しないとゆっくり回ります。画面外で描画してから表示中のキャンバスにコピーするので、3D でも PNG 保存や録画ができます。
- **カラーマップ**: Inferno・Viridis・Ice。発表用の**全画面表示**（F）。
- **10 のシーン**（それぞれ短い解説つき）: トンネル効果、共鳴トンネル、単スリット、二重スリット、回折格子、散乱、調和ポテンシャル、量子の囲い、BEC の干渉、空の箱。
- **描画**: 壁を描く・消す（Shift + ドラッグで直線）、元に戻す/やり直す。ドラッグで波束を発射し、確率流検出器を動かせます。
- **表示**: 位置（|ψ|、|ψ|²、位相の色表示）と運動量（|φ(k)|²。回折次数や遠方の干渉縞が見える）の切り替え。明るさの自動調整もできます。
- **ポテンシャルの線**: 滑らかなポテンシャルの等高線と、V = ⟨E⟩ の線（古典的に入れない領域の境界）。
- **計測**: 時刻、エネルギー ⟨E⟩、箱の中の確率、検出線の透過率（リアルタイム、吸収層の影響なし）、吸収された分も示すグラフ。
- **「測定」**: ガウス型の位置測定（σ = 0.8）。結果は |ψ|² に従って選ばれ、ψ はその周りに収縮します（局所的な運動量は保たれます）。
- **共有と書き出し**: シーン・パラメータ・描いた壁を URL に含める共有リンク、PNG 保存、動画の録画（WebM/MP4）。
- **チェック・設定・言語**: `diagnose()` のチェック、格子 128² / 256² / 512²、相互作用 g、英語と日本語。

### テスト

```bash
cd web && node --test
```

- `reference.test.mjs`: [`reference/`](../reference/README.md) の全ケースを float64 と float32 で実行。float64 は Python のリファレンスデータと ~10⁻¹³ で一致し、float32 は ~10⁻⁵（許容誤差の 1/10 以下）です。
- `scenes.test.mjs`: 全シーンが 3 種類の格子で作れて、警告なしで始まること、スリットが厳密に鏡映対称であることを確認します。シーンの解説に書いた透過率もアプリと同じ方法で測り直します（例: 共鳴トンネルは障壁 2 枚で 38%、1 枚で 25%、共鳴から外すと 19%）。
- `share.test.mjs`、`physics.test.mjs`: 共有リンクが往復で元に戻ることと、壊れたリンクを拒否することを確認します。運動量密度はパーセバルの等式を満たし、k₀ にピークを持ちます。測定結果は |ψ|² に従い、測定後の状態は規格化され、⟨k⟩ も保たれます。

CI で毎 push 実行します（約 40 秒）。

### 設計上の選択（Python 側の検証に基づく。[docs/porting.md](../docs/porting.md) 参照）

- **計算領域**: 51.2 × 51.2 の周期境界。表示は中央の 40 × 40 だけで、外側は吸収層（幅 5.6、γ_max = 6）なので、波束は吸収層にかからない位置から出発します。
- **格子**: 128²、256²（既定、Δx = 0.2）、512²。256² の float32 の 1 ステップは M シリーズ Mac の Chrome で約 4 ms で、1 フレーム 6 ステップで約 40 fps です。
- **時間刻み**: Δt = 0.01。V = 100 までの壁で V·Δt ≤ 1 になります（壁の高さのスライダーは警告を見せるために 200 まであります）。g > 0 では安定条件のため Δt を 0.9 × 2π/Σk_max² まで下げます。
- **壁**: ブラシの縁は 1 セル分だけ柔らかくしています。スリットの開口はどの格子でも対称に配置します。
- **透過率**: スペクトル法の確率流検出器を 1 フレームに 1 回サンプリングします。二重スリットでは領域の和と ~0.2% で一致します（4 次差分だと kΔx = 0.8 で ~1% 小さくなります）。
- **チェック**: `diagnose()` を発射時（厳しい基準）と 60 フレームごとに実行します。実行中の高波数成分は、発射した確率の 1e-4 を超えたときだけ警告します。硬い壁はいつも ~1e-5 の無害な成分を散乱させ、箱が空になるにつれて相対的なしきい値では警告になってしまうためです。箱に残る確率が 1% 未満になると、チェックと ⟨E⟩ の表示を止めます。
- **V = ⟨E⟩ の線**: 描いた壁が 4⟨E⟩ 未満のとき（障壁やトラップ）だけ表示します。硬い壁の周りでは壁をなぞるだけになるためです。
- **共有リンク**: 線の点は 0.05 刻みの格子（最も細かい格子より細かい）で保存し、描画中も同じ刻みを使うので、共有した絵はまったく同じように格子に載ります。格子の解像度を変えたときは、補間ではなく線を描き直します。
- **BEC シーン**: 基底状態は dτ = 0.04 の虚時間発展で準備します（約 330 ステップ、約 2 秒）。dτ = 0.005 の場合との差は ~6 × 10⁻⁵ です。

ブラウザのコンソールでは、`qwl.state` でアプリの状態を確認でき、`qwl.tick(n)` で n フレーム進められます（バックグラウンドのタブでも動作）。
