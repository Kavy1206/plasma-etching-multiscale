"""Plasma to Profile -- interactive walkthrough.

Tabs: pipeline/honesty panel, reactor results (Case 1 validation + argon
IEDF vs pressure), atomic results (MD sputter yield + damage), and a LIVE run
of the actual feature-scale model (the only part computed in the browser
session -- a full PIC-MCC run takes tens of minutes, so reactor and MD
results are precomputed data from this project's own runs, bundled in
app/data/).
"""
import json
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import streamlit as st

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

APP_DIR = Path(__file__).resolve().parent
ROOT = APP_DIR.parent
sys.path.insert(0, str(ROOT))

from src.feature_scale.viz import (center_depth, etch_iter,  # noqa: E402
                                   make_synthetic_ions, render_rgb)

st.set_page_config(page_title="Plasma to Profile", layout="wide")


@st.cache_data
def load_case1():
    d = np.load(APP_DIR / "data" / "case1_profiles.npz")
    return {k: d[k] for k in d.files}


@st.cache_data
def load_argon():
    d = np.load(APP_DIR / "data" / "argon_sweep.npz")
    return {k: d[k] for k in d.files}


@st.cache_data
def load_yield():
    return json.load(open(APP_DIR / "data" / "yield_results.json"))


case1, argon, yld = load_case1(), load_argon(), load_yield()
PRESSURES = [5, 10, 20, 50, 100]

st.title("Plasma to Profile")
st.caption("From a plasma reactor to an etched trench: reactor-scale PIC-MCC -> "
           "atomic-scale MD -> feature-scale Monte Carlo. Built from this project's own simulations.")

tab0, tab1, tab2, tab3 = st.tabs(["Pipeline & what's real", "1 - Reactor (PIC-MCC)",
                                  "2 - Atoms (MD)", "3 - Feature (live etch)"])

# ----------------------------------------------------------------- tab 0
with tab0:
    st.markdown("""
**How the scales connect:** the reactor simulation produces the energy and angle
distribution of ions hitting the wafer (IEDF/IADF) -> MD turns ion energy into a
sputter yield -> the feature-scale model uses both to decide how fast each bit of
a trench wall and floor erodes.
""")
    st.table(pd.DataFrame([
        ["Reactor, helium (Case 1)", "Real 1280-cycle PIC-MCC run", "Validated vs Turner et al. 2013: 1.2% (n_i), 3.8% (kTe), 0.65% (profile RMS)"],
        ["Reactor, argon sweep", "Real PIC-MCC, 5 pressures at 300 V", "Reduced run length; 20 and 100 mTorr have thin ion statistics"],
        ["Sputter yield", "Real LAMMPS MD, 6 energies", "10 sequential shots/energy - only the 500 eV point is solid"],
        ["Feature-scale etch", "LIVE in this page (tab 3)", "Real argon IADF is very collimated; broad-angle option is a labeled synthetic demo"],
    ], columns=["Stage", "What it is", "Caveat"]))
    st.info("Everything in tabs 1-2 is precomputed data from the repo's own runs. "
            "Only tab 3 computes anything here, and it runs the actual model code.")

# ----------------------------------------------------------------- tab 1
with tab1:
    st.subheader("Validation: helium CCP vs Turner et al. (2013), Case 1")
    x, n_i, kte = case1["x"], case1["n_i"], case1["kTe_ev"]
    x_ref, ni_ref = case1["x_ref"], case1["ni_ref"]
    mid = len(x) // 2
    ref_i = np.interp(x, x_ref, ni_ref)
    c1, c2, c3 = st.columns(3)
    c1.metric("Mid-plane n_i error", f"{100 * abs(n_i[mid] - 1.40e14) / 1.40e14:.1f}%")
    c2.metric("Mid-plane kTe error", f"{100 * abs(kte[mid] - 9.36) / 9.36:.1f}%")
    c3.metric("Profile RMS (of peak)", f"{100 * np.sqrt(np.mean((n_i - ref_i) ** 2)) / ref_i.max():.2f}%")
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.plot(x_ref * 1e3, ni_ref / 1e14, "k-", lw=2, label="Turner et al. reference")
    ax.plot(x * 1e3, n_i / 1e14, "o", ms=3.5, color="#d62728", label="this code")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("n_i (1e14 m^-3)"); ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    st.pyplot(fig, clear_figure=True)

    st.subheader("Argon: ion energy at the wafer vs pressure")
    p = st.select_slider("Pressure (mTorr)", options=PRESSURES, value=10, key="pressure")
    E = argon[f"E_{p}"]
    left, right = st.columns(2)
    with left:
        fig, ax = plt.subplots(figsize=(5.5, 3.4))
        for q in PRESSURES:
            ax.hist(argon[f"E_{q}"], bins=30, range=(0, 160), density=True, histtype="step",
                    lw=1, alpha=0.25, color="gray")
        ax.hist(E, bins=30, range=(0, 160), density=True, histtype="stepfilled", alpha=0.6,
                color="#1f77b4")
        ax.set_xlabel("ion energy at grounded electrode (eV)"); ax.set_ylabel("density")
        ax.set_title(f"{p} mTorr (others in gray)")
        st.pyplot(fig, clear_figure=True)
    with right:
        means = [argon[f"E_{q}"].mean() for q in PRESSURES]
        fig, ax = plt.subplots(figsize=(5.5, 3.4))
        ax.plot(PRESSURES, means, "o-", color="#d62728")
        ax.plot([p], [E.mean()], "o", ms=12, mfc="none", mec="k")
        ax.set_xscale("log"); ax.set_xlabel("pressure (mTorr)"); ax.set_ylabel("mean ion energy (eV)")
        ax.grid(alpha=0.3)
        st.pyplot(fig, clear_figure=True)
    m1, m2 = st.columns(2)
    m1.metric("Mean ion energy", f"{E.mean():.0f} eV")
    m2.metric("Ion samples behind this histogram", f"{E.size}")
    if E.size < 500:
        st.warning(f"Only {E.size} ion samples (a checkpoint bug dropped samples on this run; "
                   "fixed but not re-run). Indicative, not converged.")
    st.caption("Higher pressure -> more charge-exchange collisions in the sheath -> "
               "a growing low-energy population and a lower mean energy (docs/theory.md sec. 6).")
    with st.expander("Ion angle distribution (5 mTorr only)"):
        ang = argon["ang_5"]
        fig, ax = plt.subplots(figsize=(5.5, 2.8))
        ax.hist(ang, bins=40, range=(0, 20), color="#2ca02c")
        ax.set_xlabel("angle from surface normal (deg)"); ax.set_ylabel("ions")
        st.pyplot(fig, clear_figure=True)
        st.caption(f"Mean {ang.mean():.2f} deg, median {np.median(ang):.2f} deg. 5 mTorr is the only "
                   "pressure with corrected angle data (see docs/PROGRESS.md).")

# ----------------------------------------------------------------- tab 2
with tab2:
    st.subheader("Ar+ -> Si(100) sputter yield (LAMMPS, Stillinger-Weber + ZBL)")
    Ei = np.array([r["energy_ev"] for r in yld]); Yi = np.array([r["yield_mean"] for r in yld])
    tot = np.array([r["total"] for r in yld]); n = np.array([r["n_shots"] for r in yld])
    order = np.argsort(Ei)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.errorbar(Ei[order], Yi[order], yerr=(np.sqrt(np.maximum(tot, 1)) / n)[order], fmt="o",
                color="#d62728", capsize=3, label="this MD (10 impacts/energy)")
    ax.plot([100, 500, 1000], [0.07, 0.65, 0.93], "k^", ms=9, label="literature anchors")
    ax.set_xlabel("Ar+ energy (eV)"); ax.set_ylabel("yield (atoms/ion)"); ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    st.pyplot(fig, clear_figure=True)
    st.caption("500 eV gives 0.7, matching the literature 0.6-0.7. Lower energies are single-digit "
               "counts and the fitted threshold is not constrained (uncertainty larger than its value).")
    st.image(str(ROOT / "figures" / "lammps_damage_500eV.png"),
             caption="Si(100) after 4 x 500 eV Ar+ impacts: coordination analysis (OVITO) shows the "
                     "damaged, non-4-coordinated region near the surface.")

# ----------------------------------------------------------------- tab 3
REAL = "Real simulated argon, 5 mTorr (collimated, mean 2.65 deg)"
SYNTH = "Synthetic broad-angle (mechanism demo - NOT simulated argon)"
with tab3:
    st.subheader("Live trench etch (the actual feature-scale model)")
    left, right = st.columns([1, 2])
    with left:
        src_label = st.radio("Ion source", [REAL, SYNTH], key="src")
        source = "real" if src_label == REAL else "synthetic"
        width = st.slider("Trench width (cells)", 20, 100, 40, 5, key="width")
        sticking = st.slider("Neutral sticking coefficient", 0.05, 0.9, 0.3, 0.05, key="stick")
        sigma = 25
        if source == "synthetic":
            sigma = st.slider("Angular spread sigma (deg)", 5, 40, 25, key="sigma")
        n_steps = st.slider("Etch time (macro-steps)", 200, 3000, 1500, 100, key="nsteps")
        run = st.button("Run etch", type="primary", key="run")
        st.caption("Cell size is not calibrated to nanometres. Roughly 15-30 ms per step "
                   "on a laptop, so a full run is ~30-90 s.")
    with right:
        live = st.empty()
        if run:
            if source == "real":
                ion_E, ion_th = argon["E_5"], argon["ang_5"]
            else:
                ion_E, ion_th = make_synthetic_ions(float(sigma))
            every = max(n_steps // 40, 1)
            frames, depths = [], []
            bar = st.progress(0.0)
            t0 = time.time()
            for step, grid in etch_iter(ion_E, ion_th, source, width, sticking, n_steps, every):
                frames.append(render_rgb(grid, scale=4, pad=20))
                depths.append((step, center_depth(grid) - 15))
                live.image(frames[-1], caption=f"t = {step} / {n_steps}   depth = {depths[-1][1]} cells")
                bar.progress(step / n_steps)
            bar.empty(); live.empty()
            st.session_state["frames"], st.session_state["depths"] = frames, depths
            st.session_state["meta"] = f"{src_label} | width {width} | {time.time() - t0:.0f} s"
        if "frames" in st.session_state:
            frames, depths = st.session_state["frames"], st.session_state["depths"]
            st.caption(st.session_state["meta"])
            i = st.slider("Scrub through the etch", 0, len(frames) - 1, len(frames) - 1, key="scrub")
            st.image(frames[i], caption=f"t = {depths[i][0]}   depth = {depths[i][1]} cells")
            st.line_chart(pd.DataFrame({"depth below mask (cells)": [d for _, d in depths]},
                                       index=[s for s, _ in depths]))
    st.markdown("""
**What to look for.** With the real argon ions (mean angle ~2.65 deg) the walls stay straight and
depth barely depends on width - ions are too collimated for sidewall shadowing to matter at these
aspect ratios. Switch to the synthetic broad-angle source and widen the spread: the trench bows out
into a bulb. Gray fill inside the etched region is partially eroded material; cells below 0.5
material fraction stop blocking particles, so the model counts them as open.
""")
