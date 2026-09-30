"""Shared driver + renderer for the feature-scale model, used by BOTH the
Streamlit app (live animation) and scripts/make_gifs.py (README GIFs), so
what people see in the app and on GitHub comes from one code path.

Parameters match the runs reported in PROGRESS.md: "real" = the validated-
pipeline settings used for the ARDE sweep (real simulated 5 mTorr argon
IEDF/IADF); "synthetic" = the labeled broad-angle mechanism demo.
"""
from __future__ import annotations

import numpy as np

from .mc_etch import run_macro_step
from .model import FeatureGrid, YieldModel

SOURCES = {
    "real": dict(flux_per_cell=250 / 90.0, removal=0.012, chem=0.003, grazing=75.0),
    "synthetic": dict(flux_per_cell=4.0, removal=0.012, chem=0.004, grazing=70.0),
}


def make_synthetic_ions(sigma_deg: float = 25.0, n: int = 5000, seed: int = 1):
    """NOT simulated argon: uniform 60-180 eV, |Gaussian| angle spread."""
    rng = np.random.default_rng(seed)
    return rng.uniform(60, 180, n), np.abs(rng.normal(0.0, sigma_deg, n))


def center_depth(grid) -> int:
    c0, c1 = grid.trench_cols
    return int(np.sum(grid.material[:, (c0 + c1) // 2] < 0.5))


def open_width_profile(grid):
    """Open (material<0.5) cells per row, within the trench window +/-20."""
    c0, c1 = grid.trench_cols
    lo, hi = max(c0 - 20, 0), min(c1 + 20, grid.nx)
    return np.sum(grid.material[:, lo:hi] < 0.5, axis=1)


def etch_iter(ion_E, ion_theta, source, width, sticking, n_steps, every,
              nx=160, nz=100, seed=0, yield_model=None):
    """Generator: yields (step, grid) at step 0 and every `every` steps.
    The grid object is mutated in place between yields."""
    cfg = SOURCES[source]
    ym = yield_model or YieldModel(A=0.009, E_th=25.0)
    grid = FeatureGrid(nx=nx, nz=nz, trench_width_cells=width, mask_depth_cells=15)
    rng = np.random.default_rng(seed)
    n_ions = max(int(round(cfg["flux_per_cell"] * width)), 5)
    yield 0, grid
    for step in range(1, n_steps + 1):
        run_macro_step(grid, ym, ion_E, ion_theta, n_ions=n_ions, n_neutrals=n_ions,
                       sticking=sticking, removal_per_ion_hit=cfg["removal"],
                       chem_removal_scale=cfg["chem"], grazing_cutoff_deg=cfg["grazing"],
                       rng=rng)
        if step % every == 0 or step == n_steps:
            yield step, grid


def render_rgb(grid, scale: int = 4, pad: int = 20):
    """Material fraction -> uint8 RGB (solid = dark), cropped around the
    trench and upscaled by nearest-neighbour so cells stay visible."""
    import matplotlib
    cmap = matplotlib.colormaps["gray_r"]
    c0, c1 = grid.trench_cols
    lo, hi = max(c0 - pad, 0), min(c1 + pad, grid.nx)
    rgb = (cmap(grid.material[:, lo:hi])[..., :3] * 255).astype(np.uint8)
    return np.repeat(np.repeat(rgb, scale, axis=0), scale, axis=1)
