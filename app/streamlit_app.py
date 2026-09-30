"""Plasma to Profile -- interactive demo.

Sliders -> interpolated lookup over the REAL PIC-MCC argon production runs
this project generated (runs/ar_*mtorr_profiles.npz, baked into
iedf_lookup.pkl by scripts/build_app_lookup.py). This is NOT live PIC-MCC --
a full run takes tens of minutes; this interpolates between five real
simulated pressure points at a single simulated voltage (300 V).

Honesty note shown in the UI itself, not just here: only the PRESSURE axis
is backed by independently simulated data (5 real points, 5-100 mTorr). The
voltage slider uses a simple physically-motivated scaling (ion energy
roughly proportional to sheath voltage, per docs/theory.md), not independent
simulation at each voltage -- and the app says so.
"""
import pickle
from pathlib import Path

import numpy as np
import streamlit as st

st.set_page_config(page_title="Plasma to Profile", layout="centered")

LOOKUP_PATH = Path(__file__).parent / "iedf_lookup.pkl"


@st.cache_data
def load_lookup():
    with open(LOOKUP_PATH, "rb") as f:
        return pickle.load(f)


lut = load_lookup()
pressures = np.array(sorted(lut.keys()))
mean_E = np.array([lut[p]["mean_E"] for p in pressures])
n_samples = np.array([lut[p]["n_samples"] for p in pressures])
mid_kte = np.array([lut[p]["mid_kte"] for p in pressures])

st.title("Plasma to Profile")
st.caption(
    "Reactor -> sputter yield -> feature profile, chained from this "
    "project's own simulated data. Not live PIC-MCC -- see the note below "
    "each prediction for what's interpolated from real runs vs. scaled."
)

col1, col2 = st.columns(2)
with col1:
    pressure = st.slider("Pressure (mTorr)", float(pressures.min()), float(pressures.max()),
                         float(pressures[0]), step=1.0)
    voltage = st.slider("RF voltage amplitude (V)", 100, 500, 300, step=10)
with col2:
    feature_width_nm = st.slider("Trench width (nm)", 20, 200, 60, step=5)
    sticking = st.slider("Neutral sticking coefficient", 0.05, 0.9, 0.3, step=0.05)

# --- Pressure axis: real interpolation over simulated points ---------------
E_at_p = float(np.interp(pressure, pressures, mean_E))
kte_at_p = float(np.interp(pressure, pressures, mid_kte))
n_at_p = float(np.interp(pressure, pressures, n_samples))

# --- Voltage axis: NOT independently simulated -- simple physical scaling.
# theory.md sec 6: mean sheath-collected ion energy scales roughly with the
# time-averaged sheath voltage, which for a fixed geometry scales close to
# linearly with V0 amplitude. All simulated data is at V0=300V, so scale
# from there.
E_predicted = E_at_p * (voltage / 300.0)

st.subheader("Predicted ion energy at the wafer")
c1, c2, c3 = st.columns(3)
c1.metric("Mean ion energy", f"{E_predicted:.0f} eV")
c2.metric("Mid-plane kTe (at 300V, interpolated)", f"{kte_at_p:.2f} eV")
c3.metric("Real samples backing this pressure", f"{n_at_p:.0f}")

st.caption(
    f"Pressure axis ({pressure:.0f} mTorr): linearly interpolated between "
    f"real simulated argon PIC-MCC runs at {pressures.min():.0f}-{pressures.max():.0f} mTorr. "
    f"Voltage axis ({voltage} V): all underlying runs were at 300V -- this "
    f"scales linearly from there per the mechanism in docs/theory.md sec 6, "
    f"NOT independently simulated at {voltage}V."
)

if n_at_p < 50:
    st.warning(
        f"The nearest real data for this pressure has only ~{n_at_p:.0f} collected ion "
        "samples (a data-collection bug during the original run, fixed but not "
        "backfilled for this pressure -- see PROGRESS.md). Treat this prediction "
        "as indicative, not converged."
    )

# --- Simple feature-scale sketch (schematic, not a live MC run) ------------
st.subheader("Schematic profile (qualitative, not a live simulation)")
aspect_guess = min(feature_width_nm / 30.0, 3.0)  # purely illustrative scaling
st.write(
    f"At {feature_width_nm} nm width and this project's own finding that the "
    f"real simulated argon IADF is highly collimated at low pressure "
    f"(mean angle 2.65 deg at 5 mTorr), aspect-ratio-dependent etch-rate "
    f"loss is expected to be small until aspect ratio exceeds roughly 20 "
    f"(see docs/theory.md sec 6 and docs/PROGRESS.md's feature-scale "
    f"section) -- well beyond what a {feature_width_nm} nm trench reaches "
    f"in a typical process. A separate synthetic-angle demonstration in "
    f"this repo (figures/feature_scale_synthetic_bowing.png) shows what "
    f"the model predicts when that collimation assumption is relaxed."
)

st.divider()
st.caption(
    "Built from real project data: docs/PROGRESS.md has the full build log "
    "including every bug found and fixed. docs/theory.md has the physics "
    "derivations behind every number here."
)
