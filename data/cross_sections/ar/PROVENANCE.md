# Argon cross sections — provenance

## Source

LXCat itself was unreachable from both sandboxed environments this project
was built in (registration required; the domain is not in this environment's
network allow-list). Instead of leaving argon blocked, these tables were
generated from the **analytic** cross-section formulas used by **eduPIC**:

> Z. Donko, A. Derzsi, M. Vass, B. Horvath, S. Wilczek, S. Hartmann,
> P. Hartmann, "eduPIC: an educational particle-in-cell/Monte Carlo collision
> code for the study of a low-pressure capacitively coupled radiofrequency
> discharge," *Plasma Sources Sci. Technol.* **30**, 095017 (2021).
> doi:10.1088/1361-6595/ac0b48. Code: github.com/donkozoltan/eduPIC.

eduPIC is a peer-reviewed, purpose-built teaching/reference code for exactly
this kind of argon CCP simulation — not an arbitrary script. Its own citations
for the cross sections:

- **Electron–argon** (elastic, lumped excitation, ionization): analytic fits
  from A. V. Phelps & Z. Lj. Petrović, *Plasma Sources Sci. Technol.* **8**,
  R21 (1999).
- **Ion (Ar⁺)–argon** (isotropic + a "backward" component derived as
  `(qmom - qiso)/2`): analytic fits from A. V. Phelps, *J. Appl. Phys.* **76**,
  747 (1994).

**I read the formulas from the eduPIC source and independently
re-implemented and tabulated them in Python** (`scripts/generate_ar_cross_sections.py`
— the exact code used, kept in the repo so this is reproducible and
auditable) rather than vendoring eduPIC's C++. The formulas themselves —
thresholds, exponents, prefactors — are transcribed exactly; any transcription
error is mine to find, not eduPIC's, which is exactly why the generator script
is kept rather than only the output CSVs.

## What's here

| File | Process | Threshold | Frame |
|---|---|---|---|
| `Elastic_Ar.csv` | e⁻ + Ar, elastic momentum transfer | — | lab |
| `Excitation_Ar.csv` | e⁻ + Ar, lumped effective excitation | 11.5 eV | lab |
| `Ionization_Ar.csv` | e⁻ + Ar → Ar⁺ + 2e⁻ | 15.8 eV | lab |
| `Isotropic_Ar.csv` | Ar⁺ + Ar, isotropic elastic component | — | **centre of mass** |
| `Backscattering_Ar.csv` | Ar⁺ + Ar, backward component | — | **centre of mass** |

Format: `energy_eV;sigma_m2`, semicolon-delimited, 20,001 points, linear grid
(0–1000 eV electrons, 0–2000 eV ions), same interpolation/clamping rules as
the helium set (see `../PROVENANCE.md`).

Note on excitation: this is a **single lumped "effective excitation"**
channel (threshold 11.5 eV), not split into multiple states the way the
Turner helium benchmark's two excitations are. That's a real physical
simplification eduPIC itself makes (argon has a dense manifold of excited
states just above 11.5 eV; lumping them into one effective process with a
fitted cross section is standard practice for CCP/etch modeling, where the
inelastic energy loss channel matters more than resolving individual
states) — carried over here rather than introduced by this project.

Note on the "backward" ion channel: unlike the helium set (where
`Backscattering_He.csv` is a genuinely separate tabulated process), here it
is **derived**, `qback = (qmom - qiso)/2`, from Phelps' momentum-transfer and
isotropic cross sections. This is exactly what eduPIC itself does — it is how
the qmom/qiso decomposition into "isotropic" and "backward" parts is defined
in Phelps (1994) — not an approximation added by this project. One numerical
detail: the raw `(qmom - qiso)/2` goes slightly negative at some energies from
floating-point/fit-extrapolation edge effects; this is floored at zero before
tabulating (a handful of points, negligible cross-section magnitude where it
occurs — check `scripts/generate_ar_cross_sections.py` before relying on this
near the very low-energy end of the table).

## Before this goes on a resume

- [ ] Cross-check a handful of tabulated points against a first-party LXCat
      pull of the Phelps compilation, or against Phelps & Petrović (1999) /
      Phelps (1994) directly, once LXCat access is available.
- [ ] The lumped-excitation simplification should be named explicitly in any
      write-up that discusses excitation-driven effects (it is already named
      in `docs/theory.md` and the README's Limitations section).
