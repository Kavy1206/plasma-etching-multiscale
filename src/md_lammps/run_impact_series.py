"""Fire N sequential Ar impacts at one energy onto a single Si(100) slab,
via the LAMMPS Python API directly (not a static .in file), so the one-time
box setup + 300 K thermal relaxation cost is paid once and amortized across
all N shots rather than once per shot.

SCOPE CUT, stated explicitly: the original plan called for 100-200
independent impacts per energy, each on a freshly re-thermalized slab
(explicitly allowed to fall to 20-30 per the session's own priority-1
scoping). Session time constraints required going further: 10 sequential
impacts per energy on the SAME slab (brief Langevin re-equilibration between
shots, escaped/sputtered atoms deleted so they don't interfere with the next
shot), 6 energies. This trades slab independence for wall-clock -- accumulated
sub-surface damage across the 10 shots is a real, if modest, systematic effect
this design does not control for, and is reported as a limitation, not hidden.
"""
from __future__ import annotations

import argparse
import json
import os
import random

import lammps
import numpy as np

SI_MASS = 28.0855
AR_MASS = 39.948
EV_TO_J = 1.602176634e-19
AMU_TO_KG = 1.66053906660e-27


def run_series(energy_ev: float, angle_deg: float, n_shots: int, seed: int,
               nx=8, ny=8, nz=10, a0=5.431, frozen_thickness=4.0, thermo_thickness=6.0,
               n_impact_steps=1500, n_rethermalize_steps=400, workdir=".", dump_every=0):
    os.makedirs(workdir, exist_ok=True)
    cwd0 = os.getcwd()
    os.chdir(workdir)
    try:
        rng = random.Random(seed)
        slab_top = nz * a0
        z_frozen = frozen_thickness
        z_thermo = frozen_thickness + thermo_thickness
        # Vacuum headroom: dt/reset caps per-step displacement at 0.02 A, so
        # over n_impact_steps an atom can travel up to n_impact_steps*0.02 A
        # before this shot's cleanup runs. Box must be taller than
        # z_escape + that plus a safety margin, or a fast recoil overruns the
        # box and LAMMPS errors out ("lost atoms") mid-run -- hit exactly this
        # at 200 eV with the original, tighter margin; fixed by sizing the
        # headroom from n_impact_steps instead of a fixed constant.
        max_travel = n_impact_steps * 0.02
        z_escape = slab_top + 15.0
        launch_z = slab_top + 8.0
        nz_box = nz + int((max_travel + 15.0 + 15.0) / a0) + 1

        lmp = lammps.lammps(cmdargs=["-log", f"log.e{energy_ev:g}", "-screen", "none"])
        c = lmp.command
        c("units metal"); c("dimension 3"); c("boundary p p f"); c("atom_style atomic")
        c(f"lattice diamond {a0}")
        c(f"region box block 0 {nx} 0 {ny} 0 {nz_box} units lattice")
        c("create_box 2 box")
        c(f"region r_slab block 0 {nx} 0 {ny} 0 {nz} units lattice")
        c("create_atoms 1 region r_slab")
        c(f"mass 1 {SI_MASS}"); c(f"mass 2 {AR_MASS}")
        c("pair_style hybrid/overlay sw zbl 1.0 2.0")
        c("pair_coeff * * sw Si.sw Si Si")
        c("pair_coeff 1 1 zbl 14.0 14.0")
        c("pair_coeff 1 2 zbl 14.0 18.0")
        c("pair_coeff 2 2 zbl 18.0 18.0")
        c("neighbor 2.0 bin"); c("neigh_modify delay 0 every 1 check yes")

        c(f"region r_frozen block INF INF INF INF INF {z_frozen} units box")
        c(f"region r_thermo block INF INF INF INF {z_frozen} {z_thermo} units box")
        c("group g_frozen region r_frozen")
        c("group g_thermo region r_thermo")
        c("group g_mobile subtract all g_frozen")
        c(f"velocity g_mobile create 300.0 {seed} mom yes rot yes dist gaussian")
        c("fix fzero g_frozen setforce 0.0 0.0 0.0")
        c(f"fix ftherm g_thermo langevin 300.0 300.0 0.1 {seed}")
        c("fix fnve g_mobile nve")
        c("timestep 0.001")
        c("thermo 500")
        c("thermo_style custom step dt temp pe ke")
        c("run 2000")

        speed_mps = (2.0 * energy_ev * EV_TO_J / (AR_MASS * AMU_TO_KG)) ** 0.5
        speed_A_ps = speed_mps * 1.0e-2
        vz0 = -speed_A_ps * np.cos(np.radians(angle_deg))
        vx0 = speed_A_ps * np.sin(np.radians(angle_deg))

        c("fix fdt all dt/reset 1 0.0001 0.002 0.02 units box")
        c(f"region r_escape block INF INF INF INF {z_escape} INF units box")
        c("group g_old_ar type 2")  # empty at this point -- no Ar atoms exist yet
        if dump_every:
            c("compute pe_atom all pe/atom")
            c(f"dump dviz all custom {dump_every} viz.dump id type x y z vx vy vz c_pe_atom")

        sputtered_per_shot = []
        for shot in range(n_shots):
            cx = rng.uniform(0.2, 0.8) * nx * a0
            cy = rng.uniform(0.2, 0.8) * ny * a0
            c(f"create_atoms 2 single {cx} {cy} {launch_z} units box")
            # Isolate ONLY this shot's new atom: all type-2 atoms minus the
            # ones already tracked as "old" (previously-fired, now embedded
            # or otherwise still in the box). Selecting by "type 2" alone was
            # a real bug: it matched every Ar atom ever fired, so each new
            # shot was re-launching earlier embedded projectiles at the new
            # shot's velocity too -- an energy-injection bug that caused a
            # runaway temperature spike and a lost-atoms crash at 200 eV.
            c("group g_proj_tmp type 2")
            c("group g_proj_tmp subtract g_proj_tmp g_old_ar")
            c(f"velocity g_proj_tmp set {vx0} 0.0 {vz0} units box")
            c("fix fproj g_proj_tmp nve")
            c(f"run {n_impact_steps}")
            c("unfix fproj")
            c("group g_proj_tmp delete")

            c("group g_gone region r_escape")
            c("variable n_gone_count equal count(g_gone)")
            n_gone_count = int(lmp.extract_variable("n_gone_count", None, 0))
            sputtered_per_shot.append(n_gone_count)

            if n_gone_count > 0:
                c("delete_atoms group g_gone")
            c("group g_gone delete")
            c("variable n_gone_count delete")

            if shot < n_shots - 1:
                c(f"run {n_rethermalize_steps}")
                # An atom that was fast but hadn't yet crossed z_escape at the
                # check above can keep climbing during re-equilibration; check
                # again so it can't accumulate displacement across shots and
                # eventually outrun even the enlarged box (this is exactly
                # what caused the lost-atoms crash at 200 eV: the first check,
                # once per shot, wasn't tight enough).
                c("group g_gone2 region r_escape")
                c("variable n_gone2 equal count(g_gone2)")
                n_gone2 = int(lmp.extract_variable("n_gone2", None, 0))
                sputtered_per_shot[-1] += n_gone2
                if n_gone2 > 0:
                    c("delete_atoms group g_gone2")
                c("group g_gone2 delete")
                c("variable n_gone2 delete")

            c("group g_old_ar type 2")  # re-sync: whatever type-2 survived is now "old"

        lmp.close()
        return sputtered_per_shot
    finally:
        os.chdir(cwd0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy", type=float, required=True)
    ap.add_argument("--angle", type=float, default=0.0)
    ap.add_argument("--n-shots", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    result = run_series(args.energy, args.angle, args.n_shots, args.seed, workdir=args.workdir)
    out = dict(energy_ev=args.energy, angle_deg=args.angle, n_shots=args.n_shots,
              seed=args.seed, sputtered_per_shot=result, total_sputtered=sum(result),
              yield_mean=sum(result) / len(result))
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out))
