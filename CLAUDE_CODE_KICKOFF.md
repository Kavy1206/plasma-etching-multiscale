# How to use this

1. Unzip the starter repo into a local folder.
2. Do the one-time setup below.
3. Open that folder in **Claude Code**, select the **Opus** model, and paste
   the prompt in the fenced block as your first message.
4. Work through it phase by phase — it's designed to stop and check in with
   you after each phase, not run to completion unsupervised. Don't skip the
   checkpoints even if you're tempted to; that's where you actually learn
   what the code did, which is the entire point of doing this yourself
   before the recruiter conversation.

---

## One-time local setup

```bash
# 1. Python environment
cd plasma-etching-project
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. LAMMPS (needed for Phase 2) — easiest via conda
conda install -c conda-forge lammps

# 3. OVITO (needed for Phase 2 visualization)
pip install ovito
# or use the OVITO Basic desktop app if the Python module gives you trouble

# 4. Git + GitHub auth
git init
gh auth login          # if you have the GitHub CLI
# or set up an SSH key / personal access token the normal way

# 5. Create the empty GitHub repo (or let Claude Code do this in-session
#    with `gh repo create` once gh is authenticated)
gh repo create plasma-etching-multiscale --public --source=. --remote=origin
```

If `gh` isn't installed: `brew install gh` (Mac) or see cli.github.com.
If conda isn't installed: use Miniconda, or ask Claude Code to help you
install LAMMPS a different way for your OS — it's a genuinely fiddly install
and this is a reasonable thing to delegate.

---

## The kickoff prompt — paste this into Claude Code (Opus) as your first message

```
I'm a third-year Materials Science undergrad at IIT Gandhinagar preparing for
an Applied Materials summer internship recruitment process (semiconductor
equipment company — deposition, plasma etch, process engineering). I need a
real, working, well-documented project I can point recruiters and
interviewers to, in this repo, on GitHub.

This is NOT a "write me something that looks impressive" request. It needs
to actually run, actually validate against a published benchmark where I
claim it does, and I need to genuinely understand every phase well enough to
defend it in a technical interview. Prioritize correctness and my
understanding over speed. If something doesn't validate, tell me plainly
instead of smoothing it over or quietly lowering the bar for what counts as
agreement.

THE PROJECT: "Plasma to Profile" — a multiscale simulation of plasma etching.
Three coupled scales:

PHASE 0 — Theory foundation
Before writing any code, fill in docs/theory.md: derive the Bohm criterion,
the Child–Langmuir sheath relation, why an asymmetric CCP develops a negative
DC self-bias, why ICP decouples plasma density from ion energy and CCP
doesn't, and define the ALE synergy window. Source electron–Ar and Ar+–Ar
cross sections (elastic, excitation, ionization for electrons; elastic,
charge exchange for ions) from the LXCat database (Biagi or Phelps set) and
save them to data/cross_sections/.
>>> CHECKPOINT: show me the theory doc and cross-section source before
    writing any simulation code. I want to review it.

PHASE 1 — 1D3V PIC-MCC reactor solver (src/pic_mcc/)
Two parallel plates, 1D in x, argon CCP at 13.56 MHz. Build in this order,
testing after each step: (1) grid with Δx ≲ local Debye length, (2) leapfrog
particle push with ω_pe·Δt < 0.2, (3) cloud-in-cell charge deposition,
(4) Poisson solve via tridiagonal/Thomas algorithm with V(t)=V0·sin(2π·13.56e6·t)
on the driven electrode, (5) null-collision Monte Carlo using the Phase 0
cross sections, (6) absorbing wall boundaries with secondary electron
emission γ≈0.07. Use numba.njit on the hot loops — pure Python will be too
slow for ~1e5 particles/species over 300-500 RF cycles.

VALIDATE against Turner et al., "Simulation benchmarks for low-pressure
plasmas: Capacitive discharges," Physics of Plasmas 20, 013507 (2013), Case
1 — overlay your time-averaged electron density against their published
curve and report the actual percentage agreement in docs/PROGRESS.md. This
number is what the README's validation claim will be built on, so it has to
be real.

Then sweep 5/10/20/50/100 mTorr and 200/400/600 V, extracting the IEDF and
IADF at the grounded electrode. Look for the bimodal IEDF at low pressure
collapsing to a single low-energy peak at high pressure — that's the
signature result. Plot it.
>>> CHECKPOINT: show me the Turner validation overlay and the IEDF sweep
    before moving to Phase 2. If validation doesn't look right, debug it
    with me rather than proceeding — Phase 3 depends on these IEDFs being
    real.
    FALLBACK if PIC-MCC is fighting us after a few real sessions: switch to
    a volume-averaged global model (particle/power balance for n_e, T_e) with
    an analytic collisional-sheath IEDF, and say explicitly in the README
    that it's an analytic model, not a kinetic simulation. Tell me if we hit
    this point — don't silently downgrade scope without flagging it.

PHASE 2 — LAMMPS sputter yield (src/md_lammps/)
Si(100) slab, 4,000-8,000 atoms, Stillinger-Weber or Tersoff potential with
a ZBL splice for high-energy impacts, frozen bottom layers + Langevin
thermostat buffer + free NVE surface region, periodic in x/y, vacuum in z.
Fire single Ar atoms at 25/50/100/200/300/500 eV, 100-200 independent
impacts per energy at randomized lateral positions with re-thermalization
between shots, variable timestep since high-energy impacts need much smaller
dt. Measure sputter yield Y and amorphized-layer depth (via OVITO
coordination/CNA analysis). Fit Y = A(√E − √E_th) and report A, E_th.
Also get Y(θ) at 200 eV for θ = 0/30/45/60°.

OPTIONAL, do it if Phase 1 finished with time to spare: repeat on a
Cl-passivated Si surface and compare E_th against clean Si — that gap is the
ALE process window and is the single most differentiating result in this
project.
>>> CHECKPOINT: show me the fitted yield curve and the OVITO damage
    visualization before moving to Phase 3.
    FALLBACK if compute/time is tight: use a published Y(E), Y(θ) fit for
    Ar+→Si and cite the source explicitly in the README — tell me if we're
    doing this instead of running the MD ourselves.

PHASE 3 — 2D feature-scale profile model (src/feature_scale/)
Cell-based Monte Carlo on a ~200×400 grid (material fraction per cell, 0 to
1) — not a level-set method, this is easier to get topology changes
(bowing, undercut) right without remeshing. Launch pseudo-particles: ions
sampled from the Phase 1 IEDF/IADF, neutrals with cosine angular
distribution and a sticking coefficient. Neutrals re-emit diffusely off
sidewalls when they don't stick; ions reflect specularly at grazing
incidence (this is the mechanism behind microtrenching — make sure that's
actually what's producing it in the sim, not an artifact). Remove material
using the Phase 2 Y(E,θ) plus a chemical term from local radical coverage.
Generate: profile evolution animation, an ARDE curve (etch depth vs aspect
ratio, trench widths 20-200 nm), a bowing case, a microtrenching case, and a
3-panel comparison across at least 3 pressure/bias/passivation "recipes."
>>> CHECKPOINT: show me the ARDE curve and the profile figures before
    packaging.

PHASE 4 — Package (app/, docs/, root)
Streamlit app: sliders for pressure, RF voltage, feature width, passivation
sticking → precomputed lookup table (interpolated, NOT live PIC-MCC in the
browser) → predicted IEDF and profile. Deploy on Streamlit Community Cloud.
Finish README.md (already has TODOs marked — fill them with real numbers,
real figures, and an honest Limitations section: this is 1D reactor + Ar-only
chemistry + empirical yield fit + 2D feature scale, and that should be
stated plainly, not buried). Fully answer every question in docs/STUDY.md —
these need to be specific to what THIS code and THESE results actually show,
not generic textbook answers, because I'm going to be asked about this in an
interview and generic answers will fall apart under a follow-up question.
Also draft the "60-second walkthrough" and "honest weak points" sections at
the bottom of STUDY.md.

THROUGHOUT:
- Commit as you go with real, specific commit messages (not "wip" / "update").
  Don't batch everything into one commit and don't fabricate or backdate
  timestamps — real commits reflecting real work sessions.
- Update docs/PROGRESS.md at the end of every session.
- If a result looks physically wrong (wrong sign, wrong order of magnitude,
  doesn't match the qualitative behavior theory predicts), say so and debug
  it with me rather than adjusting the write-up to fit the wrong result.
- After each phase, give me a plain-English debrief before we move on: what
  you built, what it showed, and what you'd say if I were an interviewer
  asking "walk me through this part."

Start with Phase 0. Show me the theory doc before writing any simulation
code.
```

---

## After Phase 4: pushing and sharing

```bash
git add -A
git commit -m "Phase 4: Streamlit app, final README, complete STUDY.md"
git push -u origin main
```

Add the repo link to your resume's `\plasmarepo` macro and the Streamlit
link to `\fonedemo`-style line once both exist — see the TODO comments at
the top of `resume.tex`.

## One more honest note

If, after real effort, Phase 1 or 2 doesn't fully land by early October —
ship what validated. A repo with Phase 0–2 done properly, clearly documented,
with an honest "Phase 3 in progress" note, is a stronger interview asset than
five phases where you can't answer a follow-up on phases 3–4. The checkpoint
structure above is designed to make that visible to you as it happens rather
than discovering it the night before.
