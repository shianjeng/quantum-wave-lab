"""
README hero GIF: tunneling / double slit / scattering side by side.

Run:  python make_gif.py        -> figures/hero.gif
Frames are rendered directly with NumPy + Pillow (no matplotlib) so the GIF stays small.

Accuracy notes (checked against a 2x finer grid and dt = 0.001):
  * dt = 0.005 is required: with the hard walls (V = 200) dt = 0.01 gives ~17 % error in the
    scattering density; dt = 0.005 gives 0.3 %.
  * Remaining spatial error at dx = 0.125: 1-5 % in the density (mostly the staircase disk edge).
  * Tunneling: V0 = 10 > E ≈ 8 with σ = 3, so < 1 % of the transmitted probability comes from
    over-barrier (E > V0) momentum components — it really is tunneling.
  * Brightness ∝ |ψ| on a fixed scale per panel (max |ψ| at t = 0), so frames are comparable in time.
    |ψ| rather than |ψ|² keeps weak waves visible; the tunneling panel prints the transmitted probability P(x > wall).
Uses ffmpeg (palettegen) if available for better colors; otherwise falls back to Pillow.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
from matplotlib import colormaps
from PIL import Image, ImageDraw, ImageFont

from qwave.potentials import check_on_grid, double_slit, paint_disk
from qwave.solver import Grid, Solver, absorbing_mask, gaussian_packet

N, L = 320, 40.0                     # dx = 0.125
DT, STEPS_PER_FRAME, N_FRAMES = 0.005, 12, 88     # stop before waves hit the edges
PANEL, PAD, HEADER = 260, 6, 30      # pixel sizes
CMAP = colormaps["inferno"]

g = Grid((N, N), (L, L))
x, y = g.x
check_on_grid(g.dx[0], 1.0, 0.5, 1.25, 5.0, 4.0)


def scenario_tunneling():
    V = np.where((x >= -0.5) & (x < 0.5), 10.0, 0.0)                     # ⟨E⟩ ≈ 8 < V0 = 10
    return "Tunneling", V, gaussian_packet(g, (-10, 0), 3.0, (4.0, 0))


def scenario_double_slit():
    V = double_slit(g, wall_x=-4, thickness=0.5, slit_width=1.25, slit_sep=5.0, height=200)
    return "Double slit", V, gaussian_packet(g, (-13, 0), (2.0, 3.5), (4.0, 0))


def scenario_scattering():
    V = paint_disk(np.zeros(g.n), g, (0, 0), 2.0, 200)                   # hard disk
    return "Scattering", V, gaussian_packet(g, (-11, 0.8), 2.5, (4.0, 0))


def render_panel(psi, V, vmax):
    amp = np.abs(psi).T[::-1]                                            # y up
    rgb = CMAP(np.clip(amp / vmax, 0, 1))[..., :3]
    wall = (V.T[::-1] > 0)[..., None]
    rgb = np.where(wall, 0.35 * rgb + 0.65 * np.array([0.3, 0.85, 0.95]), rgb)
    img = Image.fromarray((rgb * 255).astype(np.uint8)).resize((PANEL, PANEL), Image.BILINEAR)
    return img


def main():
    font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
    font = ImageFont.truetype(font_path, 17) if os.path.exists(font_path) else ImageFont.load_default()
    small = ImageFont.truetype(font_path.replace("-Bold", ""), 13) if os.path.exists(font_path) else font

    sims = []
    for make in (scenario_tunneling, scenario_double_slit, scenario_scattering):
        title, V, psi = make()
        sims.append([title, V, psi, Solver(g, V, DT, absorber=absorbing_mask(g, DT)), np.abs(psi).max()])

    W = 3 * PANEL + 4 * PAD
    H = HEADER + PANEL + PAD + 22
    frames = []
    for f in range(N_FRAMES):
        canvas = Image.new("RGB", (W, H), (13, 17, 23))                   # GitHub dark bg
        draw = ImageDraw.Draw(canvas)
        for i, sim in enumerate(sims):
            title, V, psi, solver, vmax = sim
            if f > 0:
                sim[2] = psi = solver.step(psi, STEPS_PER_FRAME)
            x0 = PAD + i * (PANEL + PAD)
            canvas.paste(render_panel(psi, V, vmax), (x0, HEADER))
            draw.text((x0 + PANEL / 2, HEADER / 2), title, font=font, fill=(230, 237, 243), anchor="mm")
            if title == "Tunneling":
                P = np.sum(np.abs(psi[x > 0.5]) ** 2) * g.dV
                draw.text((x0 + PANEL - 8, HEADER + 8), f"transmitted {P:.3f}",
                          font=small, fill=(230, 237, 243), anchor="ra")
        t = f * STEPS_PER_FRAME * DT
        draw.text((W / 2, H - 12),
                  f"t = {t:4.2f}   ·   brightness ∝ |ψ| (fixed scale)   ·   split-step Fourier, ħ = m = 1",
                  font=small, fill=(139, 148, 158), anchor="mm")
        frames.append(canvas)

    os.makedirs("figures", exist_ok=True)
    out = "figures/hero.gif"
    fps = 20
    if shutil.which("ffmpeg"):
        with tempfile.TemporaryDirectory() as tmp:
            for k, fr in enumerate(frames):
                fr.save(f"{tmp}/f{k:04d}.png")
            src = ["-framerate", str(fps), "-i", f"{tmp}/f%04d.png"]
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", *src,
                            "-vf", "palettegen=max_colors=96:stats_mode=diff", f"{tmp}/pal.png"], check=True)
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", *src, "-i", f"{tmp}/pal.png",
                            "-lavfi", "paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle",
                            "-loop", "0", out], check=True)
    else:
        frames[0].save(out, save_all=True, append_images=frames[1:], duration=1000 // fps,
                       loop=0, optimize=True)
    print(f"saved {out}  ({os.path.getsize(out) / 1e6:.2f} MB, {len(frames)} frames)")


if __name__ == "__main__":
    main()
