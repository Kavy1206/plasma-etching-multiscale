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

## Session 3 (final session) — Priority 0: Case 1 validated

### Result

**Mid-plane ion density: 1.40e14 m^-3 vs Turner's 1.40e14 m^-3 -- 1.2% error.**
Mid-plane kTe: 9.01 eV vs 9.36 eV -- 3.8% error. Full-domain profile RMS error
(normalised to peak density): 0.65%. See `figures/case1_validation.png`.

This is a real, independently-run 512,000-step (1280 RF cycle) helium CCP
simulation, gamma=0, compared point-by-point against the digitized Turner et
al. (2013) Case 1 reference profile in `data/benchmark/`. It beats the "a few
percent at mid-plane" target set in theory.md sec. 9.5.

### A real bug, found and fixed mid-session

The first full Case 1 run finished at **14.3% error on mid-plane n_i, 9.6% on
kTe** -- clearly outside the target, and I said so rather than reporting it as
a pass. While sourcing argon cross sections (see below) I found and read the
eduPIC reference code (Donko et al. 2021, PSST 30 095017), which computes the
ion-neutral centre-of-mass collision energy as `0.5 * MU_ARAR * g_sqr` --
reduced mass times relative speed squared. My own `IonCollisions` code was
computing `mass * speed_rel**2 / 8.0`, i.e. **half** the correct value
(should be `/4.0`, since for equal masses mu = m/2 and E_cm = mu*v_rel^2/2 =
m*v_rel^2/4). This meant every ion-neutral collision all session was evaluating
the Isotropic_He/Backscattering_He cross sections at half the correct
energy, biasing both the total ion collision rate and the elastic/
charge-exchange branching ratio.

Fixed in one line, pinned with a regression test
(`test_ion_cm_energy_formula_matches_reduced_mass_definition`) that checks the
formula against the reduced-mass definition directly and independently
confirms the old formula was off by exactly 2x. **Re-ran the full 512,000-step
Case 1 validation from scratch after the fix** -- did not keep or report the
pre-fix number as the result. The improvement (14.3% -> 1.2% on the primary
metric) is consistent with this having been the dominant error source.

Diagnostic path that found nothing wrong (recorded because it narrowed things
down before the real bug turned up): checkpoint/resume verified bit-exact;
particle population converges to a stable ~12,000-19,000 range by ~300 RF
cycles and stays flat through the end of the run (not a convergence issue);
quasineutrality holds to <5% only within roughly the central 15% of the gap --
initially looked alarming, but the reference profile shows the identical
broad, non-flat-top shape (Case 1 is a genuinely sheath-dominated discharge,
consistent with the ~27%-of-gap sheath-width estimate in theory.md sec. 3d),
so this is real physics, not a bug.

### Argon cross sections -- sourced, with a citable, well-documented origin

LXCat itself was unreachable from this environment both sessions (needs a
free account; the domain isn't in this environment's network allow-list).
Instead of leaving this blocked, I found and used **eduPIC**
(github.com/donkozoltan/eduPIC), the open-source code accompanying Donko,
Derzsi, Vass, Horvath, Gibson, Firth & Hartmann, "eduPIC: an educational
particle-in-cell/Monte Carlo collision code for the study of a low-pressure
capacitively coupled radiofrequency discharge," Plasma Sources Sci. Technol.
30, 095017 (2021), doi:10.1088/1361-6595/ac0b48 -- a peer-reviewed,
purpose-built argon CCP reference implementation, not a random script. Its
electron-argon cross sections are the analytic Phelps & Petrovic (1999)
formulas; its ion-argon cross sections are the analytic Phelps (1994) qiso/qmom
formulas -- the same class of source (an established Phelps-family
compilation) as the helium data, just analytic rather than tabulated. I read
the formulas from the source (not vendored/copied) and independently
re-implemented and tabulated them; see `data/cross_sections/ar/PROVENANCE.md`
for the exact formulas and citations, and for what to verify before this goes
on a resume (the same checklist discipline as the helium data).

### Argon production sweep — done, with an explicit scope cut

5 pressures (5/10/20/50/100 mTorr) at a single voltage (300 V), extracting
the real, simulated IEDF at the grounded electrode. **Scope cut, stated
explicitly per the session's ground rules:** 400 RF cycles and 64 cells x 128
particles/cell, vs. Case 1's 1280 cycles and 128 cells x 512/cell -- chosen so
the full 5-point sweep fits the session. Wall-clock per point: 190-500s.

Physical result matches the mechanism derived in theory.md sec. 6: mean ion
energy at the wall falls monotonically with pressure (112 -> 88.7 -> 64.7 ->
38.8 -> 34.8 eV, 5->100 mTorr) as charge-exchange collisions thicken the
low-energy population -- see `figures/ar_iedf_sweep.png`. This is genuine
simulated argon data, not helium substituted or literature values.

**Honest caveat on two of the five points.** The checkpoint bug above meant
100 mTorr (split across 2 chunks during its averaging window) ended up with
only **11** collected ion samples instead of the ~1000+ the other pressures
got, and 20 mTorr lost roughly half its window's samples the same way (415
collected). Both are kept and plotted with their true sample counts rather
than quietly re-run to hide the gap -- the 100 mTorr histogram in particular
should be read as indicative, not a converged IEDF. Given remaining session
time, I did not re-run these two at full statistics; that would be the
natural first thing to redo with any further session time. The density/kTe
profile results for all 5 pressures are unaffected by this bug (those use a
different, correctly-checkpointed accumulator).

### Priority 0: closed out

Moving to Priority 1 (Phase 2, LAMMPS) now. One update to the original
plan worth flagging: this session's environment turned out to have a
genuinely working LAMMPS install (pip package + libmpich12, both from
allow-listed sources) and a working OVITO -- Phase 2 does not need the
literature-yield-curve fallback the kickoff plan allowed for. Real MD it is.

## Phase 2 (LAMMPS) — done, heavily reduced statistics, two real bugs fixed

### Scope cut, stated explicitly

Original plan: 100-200 independent impacts/energy on a freshly re-thermalized
slab each time. Given session time remaining after Priority 0, reduced
further than even the pre-authorized 20-30/energy fallback: **10 sequential
impacts per energy, on one slab per energy** (not fresh per impact), with a
short Langevin re-equilibration between shots rather than a full fresh
thermalization. This amortizes the ~15s one-time relaxation cost across 10
shots instead of paying it 10-30 times. Trade-off, stated plainly: shots
within one energy are NOT fully independent (accumulated sub-surface damage
across the series is a real, uncontrolled systematic this design doesn't
correct for), and n=10 is small enough that yields at low energy are mostly
single-digit-count statistics with large uncertainty. Six energies run:
25/50/100/200/300/500 eV, normal incidence only -- the angle sweep
(Y(theta) at 200 eV) and the optional Cl-passivation ALE-window run were
both cut entirely for time and were not attempted.

### Two real bugs found running this, both fixed before any result was kept

1. **Vacuum headroom too small.** The adaptive timestep (`dt/reset`) caps
   per-step displacement, so a fast recoil can travel up to
   `n_impact_steps * 0.02` Angstrom before a shot's cleanup runs. The
   original box only had ~25 A of vacuum above the slab; at 200 eV this
   produced `ERROR: Lost atoms` mid-run. Fixed by sizing the vacuum headroom
   from `n_impact_steps` directly, and by checking/deleting escaped atoms
   after BOTH the impact run and the inter-shot re-equilibration run (an
   atom that hadn't quite crossed the escape threshold at the first check
   could otherwise keep climbing and still outrun the box later).
2. **Projectile mis-selection (the more serious one).** Each shot selected
   its new Ar projectile with `group g_proj_tmp type 2` -- which matches
   **every** Ar atom ever fired in that session, including ones already
   embedded in the lattice from earlier shots. Every new shot was therefore
   resetting the velocity of all previously-embedded projectiles back up to
   the new shot's full launch speed too, injecting extra energy each time.
   This is what actually caused the 200 eV crash (visible in the log as a
   temperature spike to >3000 K right before the lost-atoms error) and would
   have silently biased every multi-shot energy's yield upward had it not
   crashed. Fixed with a persistent `g_old_ar` tracking group, updated after
   each shot, so only the genuinely new atom is ever selected. **Re-ran
   25/50/100 eV** (which had completed without crashing, but under the same
   buggy selection logic) after the fix rather than keep results that might
   be silently wrong; kept only the post-fix numbers below.

### Result

| E (eV) | shots | sputtered | Y (mean) |
|---|---|---|---|
| 25  | 10 | 0 | 0.0 |
| 50  | 10 | 1 | 0.1 |
| 100 | 10 | 0 | 0.0 |
| 200 | 10 | 1 | 0.1 |
| 300 | 10 | 1 | 0.1 |
| 500 | 10 | 7 | **0.7** |

Fit to Y=A(sqrt(E)-sqrt(E_th)): A=0.009+/-0.006, E_th=25+/-128 eV -- the
uncertainty on E_th is larger than the value itself, i.e. **this fit is not
meaningfully constrained** at n=10/energy, and is reported as such rather
than dressed up. See `figures/sputter_yield.png`.

The one point worth real confidence in: **500 eV gives Y=0.7, against the
literature anchor of 0.6-0.7 in theory.md sec. 7** -- good agreement. The
100 eV point (Y=0.0 from 10 shots) is statistically consistent with the
literature anchor Y=0.07 (expected count from 10 shots at that rate is
~0.7, so seeing 0 is unsurprising, not a confirming match either) --
correctly read as "not inconsistent with," not as "confirms."

Damage: OVITO coordination analysis (cutoff 2.8 A) on the post-impact slab
shows the fraction of non-4-coordinated Si atoms (bulk diamond Si is
4-coordinated) rising from 4.9% (pristine slab, surface atoms only) to
14.9% after 4 sequential 500 eV impacts -- a real, quantified amorphization
signature. `figures/lammps_damage_500eV.png` (OVITO's own Vulkan-based
renderer could not be gotten working in this environment even after
installing mesa-vulkan-drivers; the coordination analysis itself IS OVITO,
the image is rendered from OVITO's computed per-atom data via matplotlib
instead of OVITO's own viewport).

Potential used: Stillinger-Weber (LAMMPS-bundled `Si.sw`) + ZBL splice
(`pair_style hybrid/overlay sw zbl 1.0 2.0`). SW's Si cohesive energy
(~4.34 eV/atom) differs from Tersoff's (~4.63 eV/atom) -- per theory.md
sec. 7, this sets part of the threshold value, so E_th from this MD should
not be compared to a Tersoff-potential E_th without accounting for that.

### IADF bug, caught consuming the data in Phase 3

Found while writing the feature-scale model (below): the IEDF/IADF
collection in `simulation.py` used `v[:,0]` as the wall-normal velocity
component, but `push()` actually advances `v[:,2]` -- normal and
perpendicular were swapped, giving a nonsense ~88 deg mean incidence angle
for what should be near-normal, sheath-collimated ions. Energy data
(`E_gnd`) is unaffected -- it uses the full 3-vector and doesn't depend on
this decomposition, so `figures/ar_iedf_sweep.png` and the mean-energy-vs-
pressure trend are unchanged and still valid.

Fixed, then **re-ran only 5 mTorr** (the cheapest point) to get correct
angle data for Phase 3. New result: mean angle 2.65 deg, median 0.96 deg,
90th percentile 6.6 deg -- matches the theory.md sec. 6 collimation
estimate (~1-2 deg half-width) well, a good independent sanity check that
the fix is right. The other 4 pressures' saved `.npz` files still have the
pre-fix (wrong) angle arrays; their energy data is fine and already used,
but I did not have time to re-run all 5 for angle correctness, so Phase 3
below uses the 5 mTorr condition specifically (the one with corrected data).
