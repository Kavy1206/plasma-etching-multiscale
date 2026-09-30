"""Regenerate app/iedf_lookup.pkl from the real simulated argon runs.
Run this after any change to the argon production sweep in runs/."""
import glob
import pickle

import numpy as np

lut = {}
for f in glob.glob("runs/ar_*mtorr_profiles.npz"):
    p = float(f.split("ar_")[1].split("mtorr")[0])
    d = np.load(f)
    lut[p] = dict(
        mean_E=float(d["E_gnd"].mean()) if d["E_gnd"].size else None,
        n_samples=int(d["E_gnd"].size),
        mid_ni=float(d["ni"][len(d["x"]) // 2]),
        mid_kte=float(d["kte"][len(d["x"]) // 2]),
    )

with open("app/iedf_lookup.pkl", "wb") as fh:
    pickle.dump(lut, fh)
print(f"wrote app/iedf_lookup.pkl with {len(lut)} pressure points")
