"""Step 6: absorbing walls with optional secondary electron emission.

Both electrodes fully absorb every charged particle that reaches them (no
reflection). For gamma_see > 0 (argon production runs only -- config.py's
guard forbids this for helium), each absorbed ION has probability gamma_see
of launching one secondary electron back into the domain, at fixed energy
T_see, with a half-Maxwellian-like inward angular spread (cosine-weighted
polar angle, uniform azimuth -- the standard secondary-emission angular
model, e.g. Birdsall & Langdon).
"""
from __future__ import annotations

import numpy as np

from .constants import E_CHARGE as E


def apply_absorbing_walls(sp, length: float):
    """Remove particles with x < 0 or x > length. Returns (n_left, n_right)
    absorbed counts, and the positions/velocities of the absorbed ones."""
    left = sp.x < 0.0
    right = sp.x > length
    gone = left | right
    absorbed_x = sp.x[gone].copy()
    absorbed_v = sp.v[gone].copy()
    absorbed_wall = np.where(left[gone], 0, 1)  # 0 = left electrode, 1 = right
    sp.remove(~gone)
    return absorbed_wall, absorbed_x, absorbed_v


def emit_secondaries(absorbed_wall: np.ndarray, gamma_see: float, T_see_ev: float,
                     mass: float, length: float, rng: np.random.Generator):
    """For each absorbed ion, emit a secondary electron with probability gamma_see.

    Returns (x_new, v_new) arrays for the emitted electrons (possibly empty).
    """
    n = absorbed_wall.size
    if n == 0 or gamma_see <= 0.0:
        return np.zeros(0), np.zeros((0, 3))
    emit = rng.random(n) < gamma_see
    idx = np.nonzero(emit)[0]
    if idx.size == 0:
        return np.zeros(0), np.zeros((0, 3))

    speed = np.sqrt(2.0 * T_see_ev * E / mass)
    m = idx.size
    # Cosine-weighted polar angle about the inward normal: cos(theta) ~ sqrt(R)
    cos_theta = np.sqrt(rng.random(m))
    sin_theta = np.sqrt(1.0 - cos_theta ** 2)
    phi = 2.0 * np.pi * rng.random(m)

    inward = np.where(absorbed_wall[idx] == 0, 1.0, -1.0)  # left wall -> +x, right wall -> -x
    vx = inward * speed * cos_theta
    vy = speed * sin_theta * np.cos(phi)
    vz = speed * sin_theta * np.sin(phi)
    v_new = np.stack([vx, vy, vz], axis=1)

    x_new = np.where(absorbed_wall[idx] == 0, 1e-6 * length, length * (1.0 - 1e-6))
    return x_new, v_new
