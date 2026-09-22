# Cross-section and benchmark data provenance

## What is here

### `turner_benchmark_he/` — the Turner et al. (2013) prescribed set

Exactly the collision set the benchmark requires. Electron–helium from the
**Biagi v7.1** compilation via LXCat; ion–helium from **A. V. Phelps'** analytic
representation, tabulated as an isotropic-scattering component plus a
backward-scattering component.

| File | Process | Threshold | Frame |
|---|---|---|---|
| `Elastic_He.csv` | e⁻ + He, elastic momentum transfer | — | lab |
| `Excitation1_He.csv` | e⁻ + He → He* | 19.82 eV | lab |
| `Excitation2_He.csv` | e⁻ + He → He** | 20.61 eV | lab |
| `Ionization_He.csv` | e⁻ + He → He⁺ + 2e⁻ | 24.59 eV | lab |
| `Isotropic_He.csv` | He⁺ + He, isotropic component | — | **centre of mass** |
| `Backscattering_He.csv` | He⁺ + He, backward component | — | **centre of mass** |

Format: `energy_eV;sigma_m2`, semicolon-delimited, no header.

### `../benchmark/turner2013_case1_reference.csv`

Reference time-averaged profiles for Case 1, on the $L/128$ grid.
Space-delimited. **Column 0** = $x$ (m), **column 1** = electron density (m⁻³),
**column 4** = ion density (m⁻³). Remaining columns not identified; not needed.

Verification that this is a faithful Case 1 reference: its mid-plane ion density
is $1.4046\times10^{14}$ m⁻³, matching the paper's Table III value of
$0.140\times10^{15}$ m⁻³ to three significant figures.

## Rules the benchmark imposes on how these are used

1. **Linear interpolation** between tabulated points. Not spline, not log-log.
2. Above the maximum tabulated energy, **clamp** to the last tabulated value.
3. **Isotropic scattering in the CM frame for all electron processes.** The Biagi
   momentum-transfer cross sections are only consistent with transport data under
   this assumption. If you later change the scattering model, you must change the
   cross sections too.
4. After an ionising collision, split the residual energy **exactly equally**
   between primary and secondary electron.
5. **Ion tables are functions of centre-of-mass energy.** Convert before lookup.
   For a symmetric-mass collision, $E_{CM} = E_{lab}/2$. Getting this wrong
   changes the ion collision rate by a factor of ~2 and is a common bug.

## Provenance chain — READ BEFORE PUBLISHING

These files were obtained from the public reimplementation
`https://github.com/lase-unb/ccp-benchmark` (which reproduces the original LXCat
and Phelps attribution in `turner_benchmark_he/LICENSE.txt`). **That is a
secondary source.** Before this repo is presented to anyone:

- [ ] Download the electronic supplement attached to the original paper
      (doi:10.1063/1.4775084) and `diff` it against these files. Record the result
      in `docs/PROGRESS.md`.
- [ ] Register at lxcat.net (free) and export Biagi v7.1 He yourself, so the
      first-party download is what lives in the repo.
- [ ] Keep `LICENSE.txt` intact. LXCat's terms require citing both the database
      and the originating contributor.

## Citations required if these data are used

- S. F. Biagi, cross-section compilation version 7.1 (2004), retrieved from
  www.lxcat.net.
- A. V. Phelps, *J. Appl. Phys.* **76**, 747 (1994); analytic compilation at
  jila.colorado.edu/~avp/.
- M. M. Turner et al., *Phys. Plasmas* **20**, 013507 (2013).

## Still missing: argon

The production runs (pressure/voltage sweep, IEDF/IADF for Phase 3) use argon, and
those cross sections are **not** in this directory yet. Needed from LXCat, using
one compilation consistently (Biagi v7.1 *or* Phelps — do not mix):

- e⁻ + Ar: elastic momentum transfer, excitation, ionization ($E_{th}=15.76$ eV)
- Ar⁺ + Ar: elastic/isotropic **and resonant charge exchange**

The charge-exchange channel is not optional — it is the mechanism behind the
low-energy IEDF structure at high pressure (see `docs/theory.md` §6).
