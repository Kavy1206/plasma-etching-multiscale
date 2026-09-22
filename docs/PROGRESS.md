# PROGRESS

## Session 2 — Phase 1, steps 1–4 (grid, push, CIC, Poisson)

### Decisions taken

- **Gas split confirmed.** `Config.gas` switches between `"he"` (Turner
  validation, `gamma_see` forced to 0) and `"ar"` (production sweep, secondaries
  on). `Config.__post_init__` raises if a helium run is given a non-zero
  secondary yield — the benchmark specifies complete absorption with no emitted
  secondaries, so that combination is always a mistake, not a choice.
- **`turner_case(1..4)`** reproduces Table I of the benchmark exactly, including
  the prescribed rounded masses (m_e = 9.109e-31, m_i = 6.67e-27) rather than the
  CODATA values. `TURNER_REFERENCE` holds the Table III results to compare against.

### Built and tested

`src/pic_mcc/{constants,config,grid,poisson,pusher,cross_sections}.py`.
20 tests, all passing (`python -m pytest tests -q`).

| Step | What was verified | Result |
|---|---|---|
| 1 — grid | Case 1 resolution reproduces Table III | lambda_D/dx = 3.722 vs 3.72 |
| 3 — CIC | charge conservation; node/centre weights; flat density at the walls (half-width boundary nodes); gather/scatter partition of unity | exact |
| 4 — Poisson | vacuum case gives linear phi and uniform E; second-order convergence | error matches the analytic truncation error (k dx)^2/12 to 1% at every resolution |
| 0 — data | all 6 helium tables load; tabulated onsets match documented thresholds; inelastic sigma = 0 below threshold; clamping above the last point; nu*dt << 1 | passing |
| 2 + all | cold plasma oscillation | see below |

### Cold plasma oscillation (the integration test)

| Quantity | Value |
|---|---|
| omega_pe (analytic) | 6.586269e8 rad/s |
| omega_leapfrog (expected, incl. integrator dispersion) | 6.589016e8 rad/s |
| omega_measured | 6.587010e8 rad/s |
| error vs omega_pe | +0.011% |
| error vs leapfrog prediction | −0.030% |

The measured frequency lands *between* the two analytic bounds and within 0.03%
of each. The residual below the leapfrog value is the CIC force-interpolation
softening at finite k·dx, which is the expected sign and roughly the expected
size (~(k dx)^2/24 ≈ 0.01% at k dx = 0.049).

### Two things that went wrong and what they were

1. **First attempt measured a frequency 35× too low.** Cause was the test, not the
   solver: a displacement xi = A sin(pi x/L) induces a potential proportional to
   cos(pi x/L), which takes values +C and −C at the two walls. Forcing
   phi(0)=phi(L)=0 then requires a linear homogeneous correction whose gradient is
   a **uniform spurious field of −2C/L** across the domain — measured at −210 V/m
   against a restoring field of ~330 V/m. Switching to the second harmonic
   xi = A sin(2 pi x/L), whose induced potential is equal at both walls, removed
   it entirely (residual offset 2e−15 V/m) and reproduced the analytic cold-slab
   field E = e n0 xi / eps0.
2. **The oscillation then grew 17× by step 4000.** This is the finite-grid
   instability, and it is correct behaviour: a perfectly cold plasma has
   lambda_D = 0, so the accuracy condition lambda_D/dx >~ 2 is maximally violated.
   The envelope is flat to 0.1% for the first 1200 steps, so the frequency is
   measured there (19 periods, by zero-crossing regression rather than an FFT bin)
   and the script asserts the envelope really was flat. This is also a concrete
   demonstration of *why* the benchmark prescribes lambda_D/dx = 3.72.

### One scoped change to the solver

`electric_field()` at the two electrode nodes now uses a second-order one-sided
stencil instead of (phi_1 − phi_0)/dx. The first-order form approximates the field
half a cell inside the domain, which is unacceptable precisely at the electrodes
where the sheath field is largest and the delivered ion energy is set. Measured on
a known analytic field, wall error fell from **2.454% to 0.003%**; bulk unchanged.

### Performance

5.1 ms/step in pure NumPy at Case 1 particle counts (2 species × 65,536). That
extrapolates to **~44 min for Case 1's 512,000 steps** with fields only. With MCC
the estimate is 1–2 h, so **numba is not required for Case 1** — defer it and
revisit if the MCC loop proves to be the bottleneck. (Case 4, at 49.15M steps, is
a different story; it is not a required deliverable.)

### Phase 0 correction carried in

`docs/theory.md` §7 rewritten. The Bohdansky estimate of E_th ≈ 43 eV for Ar⁺→Si
was **wrong by about a factor of two**. Measured values cluster at 15–25 eV
(Wu et al. 2009: 18–20 eV; Wolsky & Zdanuk: 15–20 eV; Kuschel et al. 2025 for
a-Si: 23 eV). Section now carries the table of sources, three reasons the
semi-empirical formulas overshoot for covalent Si, a widened Phase 2 acceptance
band of **15–45 eV**, and — more usefully — absolute-yield anchors
(Y ≈ 0.07 at 100 eV, 0.6–0.7 at 500 eV, 0.93 at 1 keV) as the primary Phase 2
acceptance criterion, since E_th is a fit parameter rather than a measured onset.

### Not yet done (Phase 1 remaining)

- Step 5: null-collision MCC (scaffolding and nu_max are in `cross_sections.py`;
  the scattering kernels are not written).
- Step 6: absorbing walls + secondary emission.
- The driven-electrode boundary condition V(t) = V0 sin(2 pi f t) is implemented in
  the Poisson solver's Dirichlet argument but not yet wired into a time loop.
- Argon cross sections still missing from `data/cross_sections/ar/`.

### Validation metric, fixed in advance

As agreed in theory.md §9.5, before any comparison is run:
- Primary: relative error in mid-plane n_i vs 1.40e14 m^-3.
- Secondary: relative error in mid-plane kT_e vs 9.36 eV.
- Profile: normalised RMS error over the gap, reported separately for bulk and
  sheath regions.
- Target: a few percent at the mid-plane. 20% is a bug, not a tolerance.
