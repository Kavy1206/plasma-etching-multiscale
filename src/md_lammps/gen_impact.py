"""Generate a LAMMPS input script for one Ar -> Si(100) impact simulation.

Slab: Si(100), Stillinger-Weber potential, ZBL splice for the short-range
(high-energy nuclear stopping) part of the Ar-Si and Si-Si interaction where
SW is unphysical. Frozen bottom layers (fixed, absorb no energy but anchor
the lattice), a Langevin thermostat buffer layer beneath the free surface
(removes the impact energy on a realistic timescale instead of letting it
reflect off the frozen boundary and re-heat the surface), and a free NVE
region at the surface where the actual impact physics happens. Periodic in
x/y (bulk-like surface), a vacuum region above +z (the impact side) and
below -z is not needed since -z is the frozen boundary.

Variable timestep: high-energy impacts need small dt for the first ~200 fs
(collision cascade), then can relax to a larger dt.
"""
from __future__ import annotations

import argparse
import os

TEMPLATE = """\
# Ar -> Si(100) single-impact sputtering simulation
# energy={energy} eV, angle={angle} deg from normal, seed={seed}
units           metal
dimension       3
boundary        p p f
atom_style      atomic

# --- build the Si(100) slab -------------------------------------------
lattice         diamond {a0}
region          box block 0 {nx} 0 {ny} 0 {nz_box} units lattice
create_box      2 box
region          r_slab block 0 {nx} 0 {ny} 0 {nz} units lattice
create_atoms    1 region r_slab

mass            1 28.0855      # Si
mass            2 39.948       # Ar (projectile species, added later)

# --- potential: SW for Si-Si, ZBL splice for close-range Ar-Si/Si-Si ---
pair_style      hybrid/overlay sw zbl 1.0 2.0
pair_coeff      * * sw Si.sw Si Si
pair_coeff      1 1 zbl 14.0 14.0    # Si-Si ZBL (splices in at short range only)
pair_coeff      1 2 zbl 14.0 18.0    # Si-Ar ZBL (dominant channel for the impact)
pair_coeff      2 2 zbl 18.0 18.0    # Ar-Ar ZBL (irrelevant, single projectile, but must be defined)

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# --- region definitions for the frozen / thermostat / free layout -----
# Diamond-cubic conventional cell edge = a0; the slab spans z in [0, nz*a0]
# in box units directly from the lattice command above (no runtime lookup
# needed since nx/ny/nz/a0 are fixed at generation time).
variable        z_frozen equal {frozen_thickness}
variable        z_thermo equal {frozen_thickness}+{thermo_thickness}

region          r_frozen block INF INF INF INF INF ${{z_frozen}} units box
region          r_thermo block INF INF INF INF ${{z_frozen}} ${{z_thermo}} units box
region          r_free   block INF INF INF INF ${{z_thermo}} INF units box

group           g_frozen region r_frozen
group           g_thermo region r_thermo
group           g_free   region r_free
group           g_mobile subtract all g_frozen

# --- relax the slab at 300 K before firing -----------------------------
velocity        g_mobile create 300.0 {seed} mom yes rot yes dist gaussian
fix             fzero g_frozen setforce 0.0 0.0 0.0
fix             ftherm g_thermo langevin 300.0 300.0 0.1 {seed}
fix             fnve g_mobile nve

timestep        0.001    # ps  (metal units: ps). 1 fs pre-impact relaxation step.
thermo          200
thermo_style    custom step temp pe ke etotal
run             2000     # 2 ps thermal relaxation

# --- fire the projectile -----------------------------------------------
variable        cx equal (bound(all,xmin)+bound(all,xmax))/2 + {dx}
variable        cy equal (bound(all,ymin)+bound(all,ymax))/2 + {dy}
variable        launch_z equal {zhi_launch}

create_atoms    2 single ${{cx}} ${{cy}} ${{launch_z}} units box
group           g_projectile type 2

variable        speed_mps equal sqrt(2*{energy}*1.602176634e-19/(39.948*1.66053906660e-27))
variable        speed_A_ps equal v_speed_mps*1.0e-2      # m/s -> Angstrom/ps (metal units)
variable        vz0 equal -v_speed_A_ps*cos({angle}*PI/180)
variable        vx0 equal  v_speed_A_ps*sin({angle}*PI/180)

velocity        g_projectile set ${{vx0}} 0.0 ${{vz0}} units box
fix             fproj g_projectile nve

# --- variable timestep: small during the cascade, larger once it settles
fix             fdt all dt/reset 1 0.0001 0.002 0.02 units box

compute         pe_atom all pe/atom
thermo          100
thermo_style    custom step dt temp pe ke etotal
dump            d1 all custom 200 {dumpfile} id type x y z vx vy vz c_pe_atom

run             {n_impact_steps}

# --- report: count atoms that left the slab through +z (sputtered) -----
variable        z_escape equal {slab_top}+15.0
variable        is_escaped atom "z > v_z_escape"
group           g_escaped variable is_escaped
variable        n_sputtered equal count(g_escaped)
print           "SPUTTER_YIELD_RESULT n_sputtered=${{n_sputtered}} energy={energy} angle={angle} seed={seed}"

write_dump      all custom {finalfile} id type x y z vx vy vz c_pe_atom
"""


def write_input(path, dumpfile, finalfile, energy, angle=0.0, seed=12345,
                nx=8, ny=8, nz=10, a0=5.431, frozen_thickness=4.0,
                thermo_thickness=6.0, n_impact_steps=4000, dx=0.0, dy=0.0):
    slab_top = nz * a0
    zhi_launch = slab_top + 8.0
    nz_box = nz + int((25.0 / a0) + 1)  # >= launch height + escape margin + buffer
    text = TEMPLATE.format(energy=energy, angle=angle, seed=seed, a0=a0,
                           nx=nx, ny=ny, nz=nz, nz_box=nz_box, frozen_thickness=frozen_thickness,
                           thermo_thickness=thermo_thickness, dx=dx, dy=dy,
                           zhi_launch=zhi_launch, n_impact_steps=n_impact_steps,
                           dumpfile=dumpfile, finalfile=finalfile, slab_top=slab_top)
    with open(path, "w") as f:
        f.write(text)
    return path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy", type=float, required=True)
    ap.add_argument("--angle", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", type=int, default=4000)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    dump = os.path.basename(args.out).replace(".in", ".dump")
    final = os.path.basename(args.out).replace(".in", ".final.dump")
    write_input(args.out, dump, final, args.energy, args.angle, args.seed,
               n_impact_steps=args.steps)
    print(f"wrote {args.out}")
