import sys, time
sys.path.insert(0, ".")
import numpy as np
from src.feature_scale.model import FeatureGrid, YieldModel
from src.feature_scale.mc_etch import run_macro_step

d = np.load("runs/ar_5mtorr_profiles.npz")
E_gnd, ang_gnd = d["E_gnd"], d["ang_gnd"]
ym = YieldModel(A=0.009, E_th=25.0)

def run(trench_width, n_macro=3000, nx=200, nz=200, sticking=0.3, seed=0,
       removal_per_ion_hit=0.012, chem_removal_scale=0.003, n_ions=250, n_neutrals=250,
       mask_depth_cells=15, report_every=None):
    grid = FeatureGrid(nx=nx, nz=nz, trench_width_cells=trench_width, mask_depth_cells=mask_depth_cells)
    rng = np.random.default_rng(seed)
    history = []
    for step in range(n_macro):
        run_macro_step(grid, ym, E_gnd, ang_gnd, n_ions=n_ions, n_neutrals=n_neutrals,
                       sticking=sticking, removal_per_ion_hit=removal_per_ion_hit,
                       chem_removal_scale=chem_removal_scale, rng=rng)
        if report_every and step % report_every == 0:
            c0, c1 = grid.trench_cols; cc = (c0 + c1) // 2
            history.append((step, np.sum(grid.material[:, cc] < 0.5)))
    return grid, history

if __name__ == "__main__":
    widths = [int(w) for w in sys.argv[1:]] or [20, 35, 55, 90, 130]
    n_macro = 3000
    results = {}
    depths = {}
    for w in widths:
        t0 = time.time()
        grid, _ = run(w, n_macro=n_macro)
        results[w] = grid.material.copy()
        c0, c1 = grid.trench_cols
        cc = (c0 + c1) // 2
        depth = np.sum(grid.material[:, cc] < 0.5)
        depths[w] = depth
        print(f"width={w}: depth={depth} cells, aspect={depth/w:.2f}, elapsed={time.time()-t0:.1f}s")
    np.savez("runs/feature_scale_arde.npz",
            **{f"w{w}": m for w, m in results.items()},
            widths=np.array(widths), depths=np.array([depths[w] for w in widths]))
