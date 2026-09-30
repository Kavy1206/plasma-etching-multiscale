"""Bundle the small, committed data files the Streamlit app needs (runs/ is
gitignored, so the app can't read it on a deployed server).

NOTE: only the 5 mTorr run has CORRECTED IADF angle data (see PROGRESS.md);
the other pressures' saved angle arrays are deliberately NOT bundled."""
import glob
import json
import shutil

import numpy as np

out = {}
for f in glob.glob("runs/ar_*mtorr_profiles.npz"):
    p = float(f.split("ar_")[1].split("mtorr")[0])
    d = np.load(f)
    tag = f"{p:g}"
    for k_out, k_in in (("E", "E_gnd"), ("x", "x"), ("ni", "ni"), ("ne", "ne"), ("kte", "kte")):
        out[f"{k_out}_{tag}"] = d[k_in]
    if p == 5.0:
        out["ang_5"] = d["ang_gnd"]
np.savez_compressed("app/data/argon_sweep.npz", **out)
shutil.copy("runs/case1_profiles.npz", "app/data/case1_profiles.npz")
res = [json.load(open(f)) for f in sorted(glob.glob("runs/lammps/campaign/result_e*.json"))]
json.dump(res, open("app/data/yield_results.json", "w"), indent=1)
print("argon keys:", sorted(out.keys())[:6], "...  yield points:", [r["energy_ev"] for r in res])
