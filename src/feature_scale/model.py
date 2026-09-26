"""2D cell-based Monte Carlo feature-scale etch profile model.

Grid: material fraction f in [0,1] per cell (0=empty/etched, 1=solid Si),
not a level-set -- chosen (per the original plan) because topology changes
(bowing, undercut, microtrenching) fall out naturally without remeshing.

Physics, and what's real vs. a stated literature substitute:
- Ion flux: energy E and angle-from-normal theta sampled from the REAL
  argon IEDF/IADF this project simulated in Phase 1
  (runs/ar_<p>mtorr_profiles.npz, E_gnd/ang_gnd arrays).
- Neutral flux: cosine-distributed angle (standard for a thermal, isotropic
  reservoir above the wafer), fixed sticking coefficient. Neutrals that
  don't stick on a sidewall re-emit diffusely (cosine-distributed again),
  consistent with the physical picture of a chemisorption-limited radical
  flux.
- Sputter yield Y(E): the fitted curve from Phase 2's real (if
  statistics-limited) MD, Y=A(sqrt(E)-sqrt(E_th)). Ions with E<E_th deposit
  no physical-sputter removal.
- Angular yield dependence Y(theta)/Y(0): Phase 2's own angle sweep was cut
  for time (see PROGRESS.md) -- this uses the literature Yamamura (1996)
  functional form instead, stated here and in the README, not silently
  substituted: Y(theta)/Y(0) = cos(theta)^-f * exp[f(1-1/cos(theta))],
  f=1.6 (a representative value for this class of ion-target system; not
  independently fit for Ar-Si in this project).
- Chemical (radical-assisted) removal: local rate proportional to the
  ADSORBED radical coverage on that cell (tracked per-cell, Langmuir
  kinetics: adsorption from neutral flux, consumption by ion-assisted
  removal) -- the standard ion-enhanced/chemically-assisted etch picture.
- Ion reflection: specular at grazing incidence (angle-from-surface-normal
  beyond a cutoff), which is the stated mechanism for microtrenching -- an
  ion that reflects off a sloped sidewall lands at the trench base corner
  at a near-yield-maximizing angle and preferentially removes material
  there.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class YieldModel:
    A: float = 0.05           # physical-sputter fit prefactor (atoms/ion/sqrt(eV))
    E_th: float = 20.0         # eV
    yamamura_f: float = 1.6    # angular-dependence shape parameter (literature, see module docstring)
    chem_rate: float = 0.15    # chemical removal rate per unit radical coverage per ion strike

    def physical_yield(self, E, theta_rad):
        Y0 = np.maximum(self.A * np.sqrt(np.maximum(E - self.E_th, 0.0)), 0.0)
        cos_t = np.clip(np.cos(theta_rad), 0.05, 1.0)
        ang_factor = cos_t ** (-self.yamamura_f) * np.exp(self.yamamura_f * (1.0 - 1.0 / cos_t))
        # cap the angular enhancement (Yamamura's own form is only valid well
        # short of grazing incidence; beyond ~80 deg reflection dominates
        # instead, handled separately by the specular-reflection logic)
        return Y0 * np.minimum(ang_factor, 4.0)


class FeatureGrid:
    """nx (lateral) x nz (depth) grid. z=0 is the original mask/surface top;
    increasing z goes DOWN into the wafer. material[iz,ix] in [0,1]."""

    def __init__(self, nx=200, nz=400, dx=1.0, trench_width_cells=60,
                mask_depth_cells=20):
        self.nx, self.nz, self.dx = nx, nz, dx
        self.material = np.ones((nz, nx), dtype=np.float64)
        self.radical = np.zeros((nz, nx), dtype=np.float64)
        # Open a trench in a mask layer at the top: material=0 in the trench
        # opening down to mask_depth_cells, material=1 (mask) elsewhere at top.
        self.mask_depth = mask_depth_cells
        c0 = nx // 2 - trench_width_cells // 2
        c1 = c0 + trench_width_cells
        self.trench_cols = (c0, c1)
        self.material[:mask_depth_cells, :] = 1.0
        self.material[:mask_depth_cells, c0:c1] = 0.0
        # Below the mask: open Si everywhere within the trench footprint down
        # to a "start" surface; solid Si below that.
        self.material[mask_depth_cells:, :] = 1.0
        self.material[mask_depth_cells:, c0:c1] = 0.0
        self.is_mask = np.zeros((nz, nx), dtype=bool)
        self.is_mask[:mask_depth_cells, :] = True
        self.is_mask[:mask_depth_cells, c0:c1] = False

    def surface_normal(self, iz, ix):
        """Estimate the local outward surface normal via the material-fraction
        gradient (standard cell-based approach: normal points from solid
        toward empty, i.e. along -grad(material))."""
        m = self.material
        nzg, nxg = m.shape
        izm, izp = max(iz - 1, 0), min(iz + 1, nzg - 1)
        ixm, ixp = max(ix - 1, 0), min(ix + 1, nxg - 1)
        gz = (m[izp, ix] - m[izm, ix]) / 2.0
        gx = (m[iz, ixp] - m[iz, ixm]) / 2.0
        n = np.array([gx, gz])
        norm = np.linalg.norm(n)
        if norm < 1e-9:
            return np.array([0.0, -1.0])
        return n / norm

    def surface_cells(self):
        """Cells that are partially/fully solid AND have an empty neighbor
        (exposed surface -- the only cells that can be etched this step)."""
        m = self.material
        solid = m > 0.01
        empty_neighbor = np.zeros_like(solid)
        empty_neighbor[:-1, :] |= m[1:, :] < 0.99
        empty_neighbor[1:, :] |= m[:-1, :] < 0.99
        empty_neighbor[:, :-1] |= m[:, 1:] < 0.99
        empty_neighbor[:, 1:] |= m[:, :-1] < 0.99
        return solid & empty_neighbor
