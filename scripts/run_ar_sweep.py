"""Argon production sweep: pressure series at fixed voltage, extracting the
IEDF/IADF at the grounded electrode.

SCOPE CUT, stated here and in PROGRESS.md: this is a REDUCED sweep relative
to the original plan (400 RF cycles instead of Turner-scale 1280-15360;
64 cells x 128 particles/cell instead of 128-512) so the full 5-pressure
sweep fits this session's remaining time budget. This trades statistical
smoothness for coverage -- real, independently-simulated argon IEDFs at
each pressure, not helium data and not literature substitutes, but with
wider point-to-point noise than the Case 1 validation run.

Usage: python -m scripts.run_ar_sweep <pressure_mtorr>
"""
import json
import sys
import time

import numpy as np

from src.pic_mcc.config import Config
from src.pic_mcc.constants import K_B
from src.pic_mcc.simulation import Simulation

MTORR_PA = 0.133322368

pressure_mtorr = float(sys.argv[1])
n_gas = pressure_mtorr * MTORR_PA / (K_B * 300.0)

cfg = Config(gas="ar", n_cells=64, particles_per_cell=128,
            n_cycles=400, n_cycles_average=50,
            gamma_see=0.05, T_see=2.0,
            n_gas=n_gas, voltage=300.0, n_plasma_init=8e14,
            seed=20260924 + int(pressure_mtorr))

ckpt = f"runs/ar_{pressure_mtorr:g}mtorr.npz"
sim = Simulation(cfg, ckpt, collect_iedf=True)
t0 = time.time()
res = sim.run(max_wall_seconds=280, progress_every=20000, checkpoint_every=20000)
print(json.dumps({**res, "wall_s": round(time.time() - t0, 1),
                  "pressure_mtorr": pressure_mtorr, "n_gas": n_gas}))

if res["step"] >= res["n_steps"]:
    x, ne, ni, kte = sim.averaged_profiles()
    E_gnd, ang_gnd = sim.iedf_iadf(1)  # grounded/right electrode
    np.savez(f"runs/ar_{pressure_mtorr:g}mtorr_profiles.npz",
            x=x, ne=ne, ni=ni, kte=kte, E_gnd=E_gnd, ang_gnd=ang_gnd,
            pressure_mtorr=pressure_mtorr)
    print(f"n_i(mid)={ni[len(x)//2]:.3e}  kTe(mid)={kte[len(x)//2]:.2f} eV  "
         f"n_ions_collected={E_gnd.size}  mean_E={E_gnd.mean():.1f} eV  DONE")
