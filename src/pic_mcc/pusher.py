"""Leapfrog particle push, 1D in space, 3V in velocity.

Velocities are staggered half a step behind positions:

    v^{n+1/2} = v^{n-1/2} + (q/m) E(x^n) dt
    x^{n+1}   = x^n + vx^{n+1/2} dt

Only E_x acts (electrostatic, no magnetic field), so vy and vz are untouched
by the push. They exist because collisions scatter into the full 3D velocity
space, and because the transverse velocity at the wall is what sets the ion
angular distribution (IADF) that Phase 3 needs.

Accuracy condition (Turner et al. eq. 1): omega_pe * dt <~ 0.2.
"""
from __future__ import annotations

import numpy as np


def push(x, vx, E_at_particle, qm: float, dt: float):
    """In-place leapfrog advance. qm = q/m for the species. Returns (x, vx)."""
    vx += qm * E_at_particle * dt
    x += vx * dt
    return x, vx


def half_step_back(vx, E_at_particle, qm: float, dt: float):
    """Initialise the leapfrog: retard v by half a step at t=0."""
    return vx - 0.5 * qm * E_at_particle * dt


def kinetic_energy_ev(vx, vy, vz, mass: float) -> np.ndarray:
    from .constants import EV
    return 0.5 * mass * (vx ** 2 + vy ** 2 + vz ** 2) / EV
