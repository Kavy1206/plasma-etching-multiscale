"""Profile-evolution GIFs for the README (the animation the original plan
called for). Uses src/feature_scale/viz.py -- same code path as the app.
Run: python -m scripts.make_gifs"""
import numpy as np
from PIL import Image, ImageDraw

from src.feature_scale.viz import (center_depth, etch_iter, make_synthetic_ions,
                                   open_width_profile, render_rgb)


def make_gif(path, frames_rgb, labels, hold_ms=1800, frame_ms=90):
    pil = []
    for arr, lab in zip(frames_rgb, labels):
        h, w, _ = arr.shape
        canvas = Image.new("RGB", (w, h + 22), "white")
        canvas.paste(Image.fromarray(arr), (0, 22))
        ImageDraw.Draw(canvas).text((4, 5), lab, fill="black")
        pil.append(canvas)
    durations = [frame_ms] * (len(pil) - 1) + [hold_ms]
    pil[0].save(path, save_all=True, append_images=pil[1:], duration=durations,
                loop=0, optimize=True)


def run(label, source, E, th, width, sticking, n_steps, every, path, png):
    frames, labels = [], []
    for step, grid in etch_iter(E, th, source, width, sticking, n_steps, every):
        frames.append(render_rgb(grid, scale=5, pad=25))
        labels.append(f"{label} | t={step:4d} | depth {center_depth(grid)-15:2d}")
    make_gif(path, frames, labels)
    Image.fromarray(frames[-1]).save(png)
    wp = open_width_profile(grid)
    rows = np.arange(15, 15 + 30)
    print(f"{label}: final depth {center_depth(grid)-15} cells below mask; "
          f"open width at mask bottom={wp[15]}, max within etched region={wp[15:center_depth(grid)].max()}")


d = np.load("app/data/argon_sweep.npz")
run("REAL 5 mTorr Ar", "real", d["E_5"], d["ang_5"], 50, 0.3, 3000, 60,
    "figures/profile_evolution_real.gif", "/tmp/last_real.png")
Es, ths = make_synthetic_ions(25.0)
run("SYNTHETIC broad angle", "synthetic", Es, ths, 30, 0.5, 2500, 50,
    "figures/profile_evolution_synthetic.gif", "/tmp/last_synth.png")
