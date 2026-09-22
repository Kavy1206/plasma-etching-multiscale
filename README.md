# Plasma to Profile: Multiscale Simulation of Plasma Etching

> **Status:** 🚧 In progress — see [PROGRESS.md](docs/PROGRESS.md) for what's done vs. planned.

A 13.56 MHz capacitively coupled plasma (CCP) reactor model coupled to atomistic
sputter-yield calculations and a feature-scale etch profile simulator — built to
understand, from first principles, how plasma conditions in a reactor turn into
the etch profile on a wafer.

<!-- TODO: once Phase 1 is done, replace this line with the actual validation
     result, e.g. "Validated against Turner et al. (2013) Case 1: electron
     density profile agrees to within X%." Do not claim validation before
     you've actually run the comparison. -->

## What this does

1. **Reactor scale** — a 1D3V particle-in-cell / Monte-Carlo-collision (PIC-MCC)
   solver for an argon CCP discharge, giving the ion energy and angular
   distribution (IEDF/IADF) that hits the wafer.
2. **Atomic scale** — LAMMPS molecular dynamics of Ar⁺ bombarding Si(100),
   giving sputter yield and surface damage as a function of ion energy and angle.
3. **Feature scale** — a 2D Monte Carlo profile simulator that uses (1) and (2)
   to predict trench etch profiles: aspect-ratio-dependent etching, sidewall
   bowing, microtrenching.

## Quick look

<!-- TODO: embed 3–4 key figures here once generated, e.g.:
![IEDF vs pressure](figures/iedf_pressure_sweep.png)
![Sputter yield fit](figures/yield_curve.png)
![Etch profile evolution](figures/profile_animation.gif)
-->

## Try it

<!-- TODO: link once Streamlit app is deployed -->
Live demo: *(coming after Phase 4)*

Run locally:
```bash
pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

## Repo structure

```
src/pic_mcc/        1D3V PIC-MCC reactor solver
src/md_lammps/       LAMMPS input scripts + yield-fitting analysis
src/feature_scale/   2D Monte Carlo feature-profile model
app/                 Streamlit interactive demo
docs/theory.md       Derivations: Bohm criterion, sheath physics, ALE window
docs/STUDY.md        Plain-English explainer, organized as Q&A
docs/PROGRESS.md     What's built, what's validated, what's left
figures/             Generated plots and animations
data/                Cross sections, benchmark reference data
```

## Validation

<!-- TODO: this section is the credibility of the whole repo — fill it in
     honestly as each phase completes. State what you validated against,
     the actual agreement (a number), and what you have NOT validated. -->

## Limitations

<!-- TODO: fill in honestly — e.g. 1D reactor model (no radial effects),
     Ar-only chemistry (no molecular gas chemistry), empirical yield fit
     rather than full reactive MD, 2D feature model (no 3D corner effects) -->

## References

- Turner, M. M. et al., "Simulation benchmarks for low-pressure plasmas:
  Capacitive discharges," *Physics of Plasmas* **20**, 013507 (2013).
- Lieberman, M. A. & Lichtenberg, A. J., *Principles of Plasma Discharges
  and Materials Processing*, 2nd ed., Wiley (2005).
- Kanarik, K. J. et al., on atomic layer etching synergy windows.
