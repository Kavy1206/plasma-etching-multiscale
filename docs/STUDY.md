# STUDY

Every answer here is tied to what *this repo's code and results* actually
show, not a generic textbook answer. Numbers are pulled from
[`docs/PROGRESS.md`](PROGRESS.md); if you're asked a follow-up this doesn't
cover, that's where to look next.

---

### "Walk me through this project."

Three coupled scales. A 1D3V particle-in-cell / Monte Carlo collision
(PIC-MCC) code models the plasma reactor itself — electrons and ions in a
13.56 MHz capacitive discharge, self-consistently solving Poisson's equation
every timestep and doing null-collision Monte Carlo for electron-neutral and
ion-neutral collisions. That gives the ion energy and angle distribution
(IEDF/IADF) hitting the wafer. LAMMPS molecular dynamics of single Ar⁺ ions
hitting a Stillinger-Weber Si(100) slab turns that IEDF into a sputter yield
curve. A 2D cell-based Monte Carlo feature-scale model then uses both to
predict the actual trench cross-section as material is etched away —
aspect-ratio-dependent etching, bowing, microtrenching.

The reactor code is validated against a published benchmark (Turner et al.
2013): 1.2% error on mid-plane ion density, 3.8% on electron temperature,
against a from-scratch independent implementation.

### "What does 'validated' actually mean here, precisely?"

Turner et al. (2013) is a fixed reference problem — specific helium CCP
parameters, specific cross sections, specific run length (1280 RF cycles) —
that five independent PIC-MCC codes were run against and shown to agree with
each other to within a few percent at the mid-plane. I ran my code on the
exact same Case 1 parameters (`turner_case(1)` in `config.py`) and compared
my time-averaged ion density and electron temperature profile, point by
point, against their published reference curve
(`data/benchmark/turner2013_case1_reference.csv`). 1.2% / 3.8% / 0.65%
RMS is where I landed. That's a real number from a real 512,000-step run,
not an estimate.

**What it does not mean:** it does not validate the argon production runs
(different gas, different cross sections, no independent reference for
those), the LAMMPS yield curve, or the feature-scale model. Those are
internally consistent and physically reasonable (see below) but not
benchmarked against an external reference the way Case 1 is.

### "Why helium for validation and argon for production? Isn't that inconsistent?"

Because Turner et al.'s benchmark *is* helium — that's what has a published
reference to check against. Argon is what a real etch reactor actually runs
(higher mass, resonant charge exchange, different ionization threshold), so
that's what feeds the sputter-yield and feature-scale stages. This is a
config switch (`Config.gas`), not two different codebases — same PIC-MCC
core, different cross-section tables. `config.py`'s `__post_init__` actively
forbids running helium with secondary emission on, because that combination
would silently invalidate the Case 1 comparison.

### "What's the Bohm criterion and why does it matter here?"

Ions entering the sheath must already be moving at least the Bohm velocity
$u_B=\sqrt{k_BT_e/M}$, or no stable sheath can form — derived in
`docs/theory.md` §2 from requiring the sheath's Poisson equation have a
monotonic (non-oscillatory) solution. It gives a cheap sanity check I
*derived but did not run against the code*: the Bohm flux formula
$J_i=0.61\,n_0\,e\,u_B$ predicts 0.205 A/m² for Case 1's parameters, against
Turner's reported 0.219 A/m² (6.4% off, about what that formula's
flat-presheath assumption costs). My simulation never recorded wall ion flux,
so I can't claim it matched; the Case 1 validation rests on the density and
temperature profiles instead. Comparing the simulated wall flux to 0.219 A/m²
would be a good extra check to add.

### "Explain the null-collision Monte Carlo method."

Rather than compute a real collision probability for every particle every
step (expensive — cross sections aren't cheap to evaluate at every
particle's exact energy), precompute an upper bound
$\nu_{max}=n_{gas}\max_E[\sigma_{tot}(E)v(E)]$. Each step, flag a particle
as a collision *candidate* with probability $1-e^{-\nu_{max}\Delta t}$ (cheap
— one random draw per particle, no cross-section evaluation needed). Only
for the (much smaller) flagged subset, evaluate the *true* local collision
frequency and accept with probability $\nu(E)/\nu_{max}$ — a rejected
candidate is a "null" collision, does nothing. This reproduces the correct
energy-dependent rate while only paying the expensive per-particle
cross-section evaluation for a small fraction of the population each step.
Implemented in `mcc.py`; `test_null_collision_rate_matches_analytic_elastic_only_rate`
checks the resulting rate against the analytic $n\sigma v\Delta t$ prediction
directly.

### "Bug X — walk me through how you found and fixed it." (pick the most relevant one)

The clearest example: the ion-neutral centre-of-mass energy formula. Cross
sections for ion-neutral collisions are tabulated against CM energy,
$E_{cm}=\tfrac12\mu v_{rel}^2$ with reduced mass $\mu=m_1m_2/(m_1+m_2)$. For
equal masses (Ar⁺ on Ar, He⁺ on He), $\mu=m/2$, so $E_{cm}=mv_{rel}^2/4$. My
first implementation had $mv_{rel}^2/8$ — exactly half. It didn't crash or
look obviously wrong; Case 1 still ran and gave a *plausible-looking* result
(14.3% error on mid-plane density) that just happened to be outside my own
pre-declared "a few percent" target. I found the actual bug by chance,
cross-referencing a different reference code (eduPIC) while sourcing argon
cross sections, noticed their formula used a factor of 4 not 8, worked out
by hand which was right, wrote a test pinning the correct formula against
the reduced-mass definition directly, and **re-ran the entire 512,000-step
validation from scratch** rather than patch the old number. Result went from
14.3%/9.6% error to 1.2%/3.8%. The lesson I'd give in an interview: a
plausible-looking wrong answer is more dangerous than a crash, and the fix
was cross-checking the physics against an independent source, not staring
harder at my own code.

(Four other real bugs, similarly documented with what broke and how it was
caught, are in PROGRESS.md — the IEDF/IADF checkpoint data loss, the IADF
normal/perpendicular swap, the LAMMPS vacuum-headroom crash, and the LAMMPS
projectile mis-selection bug that was the real cause of that crash.)

### "Why does the ARDE curve come out flat? Isn't ARDE supposed to happen?"

It is, but it needs enough aspect ratio for sidewall shadowing to matter
geometrically, and that threshold depends on how collimated the incoming
ion angular distribution is. The real simulated IADF at 5 mTorr has a mean
angle of 2.65° from normal (see `docs/theory.md` §6's derivation of why low
pressure means a collimated sheath). By a rough geometric estimate,
shadowing becomes significant around aspect ratio $1/\tan(2.65°)\approx21$ —
the widths I tested only reached aspect ratio 1.9. That estimate is a
hypothesis: I did not run aspect ratios that high, so I haven't shown the
curve turns over there. What I can say is that a flat curve in the tested
range is consistent with that explanation, not a sign the model is broken. I demonstrated the model can produce
ARDE-adjacent effects (specifically bowing) by re-running with a
deliberately broadened synthetic angular distribution — clean, textbook
bowed profile, and I labeled that case as synthetic everywhere, not real
argon data, because presenting it as if it were would be misleading.

### "What would you do differently with more time / compute?"

In priority order: (1) re-run the two thin-statistics argon pressure points
(20, 100 mTorr — 415 and 11 collected ion samples respectively, from a
checkpoint bug fixed mid-session but not backfilled); (2) run the LAMMPS
campaign at full independent-impact statistics (100+ per energy, fresh slab
each time, currently 10 sequential shots per energy) so the yield-curve fit
is actually constrained — right now $E_{th}$'s uncertainty is larger than
its value; (3) extend the feature-scale model to high enough aspect ratio
(~20+) to see the real simulated argon IADF actually produce ARDE, instead
of only demonstrating the mechanism with a synthetic case; (4) source
first-party LXCat cross sections instead of the eduPIC-derived analytic
argon set.

---

## 60-second walkthrough

"I built a three-scale simulation of plasma etching: a kinetic model of the
plasma reactor itself, molecular dynamics of individual ion impacts on
silicon, and a profile model that turns those into an actual etched trench
shape. The reactor code is validated against a published benchmark to
1.2% on ion density — that's a real, independently-run comparison, not an
estimate, and I can show you the bug I found and fixed to get there. The
production argon runs feed a real LAMMPS sputter-yield calculation, which
lands right on the literature value at 500 eV. The feature-scale model then
shows something I think is a genuinely interesting result: at the low
pressure this reactor runs at, the ion angular spread is so narrow — 2.65°
— that aspect-ratio-dependent etching basically doesn't happen at realistic
trench geometries, which is exactly why low pressure is used for anisotropic
etching in the first place. I can also show you a demo case with a broader
angular spread where the model produces a textbook bowed profile, so you can
see the mechanism is really there, just not triggered by this particular
gas condition."

## Honest weak points

Say these before you're asked, not after:

- The LAMMPS yield curve is statistics-limited (10 shots/energy,
  sequential not independent) — the only point I'd defend with real
  confidence is 500 eV (Y=0.7, matches literature); the fitted $E_{th}$ is
  not meaningfully constrained and I'd say so if pushed on it.
- Two of five argon pressure points have thin IEDF statistics from a bug I
  found but didn't have time to fully correct for (the sweep's *density and
  temperature* results are unaffected — only the ion-energy histograms at
  20 and 100 mTorr are thin).
- The feature-scale model's headline result (flat ARDE) is real but is a
  negative result at the tested conditions, not a positive demonstration of
  ARDE, bowing, and microtrenching all together — I have bowing from a
  labeled synthetic case, not from the real simulated argon condition.
- Argon cross sections are analytic literature fits sourced via a reference
  code, not a first-party LXCat pull (LXCat itself was unreachable in the
  build environment).
- No molecular gas chemistry, no 3D effects, 1D reactor model — all stated
  in the README, not discovered by an interviewer first.
