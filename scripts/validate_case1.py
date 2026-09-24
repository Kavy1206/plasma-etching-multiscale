"""Compute the Case 1 validation metrics agreed in theory.md sec. 9.5 and
compare against Turner et al. (2013) Table III / the Case 1 reference profile."""
import json

import numpy as np

from src.pic_mcc.config import TURNER_REFERENCE, turner_case
from src.pic_mcc.simulation import Simulation

cfg = turner_case(1)
sim = Simulation(cfg, "runs/case1.npz")
assert sim.step == cfg.n_steps, f"run not finished: {sim.step}/{cfg.n_steps}"

x, n_e, n_i, kTe_ev = sim.averaged_profiles()
mid = len(x) // 2

ref = TURNER_REFERENCE[1]
n_i_mid_ours = n_i[mid]
kTe_mid_ours = kTe_ev[mid]

err_ni = 100 * abs(n_i_mid_ours - ref["n_i_midplane"]) / ref["n_i_midplane"]
err_kTe = 100 * abs(kTe_mid_ours - ref["kTe_midplane"]) / ref["kTe_midplane"]

# Reference profile for the full-gap RMS comparison
refdata = np.loadtxt("data/benchmark/turner2013_case1_reference.csv")
x_ref, ni_ref = refdata[:, 0], refdata[:, 4]
ni_ref_interp = np.interp(x, x_ref, ni_ref)

sheath_frac = 0.15  # outer 15% of the gap on each side, per theory.md sec 3d estimate (~27%);
                    # use a conservative inner definition so "bulk" isn't contaminated
is_sheath = (x < sheath_frac * cfg.length) | (x > (1 - sheath_frac) * cfg.length)
is_bulk = ~is_sheath

def rms_norm(mask):
    diff = n_i[mask] - ni_ref_interp[mask]
    return float(np.sqrt(np.mean(diff ** 2)) / np.max(ni_ref_interp))

rms_bulk = rms_norm(is_bulk)
rms_sheath = rms_norm(is_sheath)
rms_all = rms_norm(np.ones_like(x, dtype=bool))

result = dict(
    n_i_midplane_ours=n_i_mid_ours, n_i_midplane_turner=ref["n_i_midplane"], err_ni_pct=err_ni,
    kTe_midplane_ours=kTe_mid_ours, kTe_midplane_turner=ref["kTe_midplane"], err_kTe_pct=err_kTe,
    rms_norm_bulk_pct=100 * rms_bulk, rms_norm_sheath_pct=100 * rms_sheath, rms_norm_all_pct=100 * rms_all,
    avg_samples=sim.avg.count, n_electrons_final=sim.electrons.n, n_ions_final=sim.ions.n,
)
print(json.dumps(result, indent=2))

np.savez("runs/case1_profiles.npz", x=x, n_e=n_e, n_i=n_i, kTe_ev=kTe_ev,
        x_ref=x_ref, ni_ref=ni_ref)
