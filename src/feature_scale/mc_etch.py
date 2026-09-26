"""Vectorized Monte Carlo particle tracing + material removal for the 2D
feature-scale model. All particles in one macro-step are traced together
(numpy arrays), stepping by a fixed small increment along each particle's
own direction and checking the grid cell for a hit, so the cost scales with
the number of trace steps (bounded by grid depth) rather than per-particle
Python loops.
"""
from __future__ import annotations

import numpy as np


def trace_batch(material, x0, z0, dirx, dirz, step_frac=0.4, max_steps=None):
    """Trace a batch of straight-line rays through the material grid until
    each either hits a solid cell (fraction > 0.5) or exits the domain.

    Returns (hit_mask, hit_ix, hit_iz, hit_normal) for particles that hit,
    and exits (out-of-domain) are simply excluded (hit_mask False there too
    -- caller doesn't distinguish, which is fine: an exited particle removes
    no material either way).
    """
    nz, nx = material.shape
    n = x0.size
    x = x0.copy().astype(np.float64)
    z = z0.copy().astype(np.float64)
    alive = np.ones(n, dtype=bool)
    hit = np.zeros(n, dtype=bool)
    hit_ix = np.full(n, -1, dtype=np.int64)
    hit_iz = np.full(n, -1, dtype=np.int64)

    if max_steps is None:
        max_steps = int(2 * nz / step_frac)

    for _ in range(max_steps):
        if not np.any(alive):
            break
        x[alive] += dirx[alive] * step_frac
        z[alive] += dirz[alive] * step_frac
        ix = np.clip(np.round(x).astype(np.int64), 0, nx - 1)
        iz = np.round(z).astype(np.int64)

        out_of_domain = alive & ((z < -1) | (z >= nz) | (x < -1) | (x >= nx))
        alive[out_of_domain] = False

        check = alive & (iz >= 0) & (iz < nz)
        if np.any(check):
            idx = np.nonzero(check)[0]
            solid = material[iz[idx], ix[idx]] > 0.5
            newly_hit = idx[solid]
            if newly_hit.size:
                hit[newly_hit] = True
                hit_ix[newly_hit] = ix[newly_hit]
                hit_iz[newly_hit] = iz[newly_hit]
                alive[newly_hit] = False

    return hit, hit_ix, hit_iz


def local_normal(material, ix, iz):
    """Vectorized surface-normal estimate (points from solid to empty)."""
    nz, nx = material.shape
    izm = np.clip(iz - 1, 0, nz - 1); izp = np.clip(iz + 1, 0, nz - 1)
    ixm = np.clip(ix - 1, 0, nx - 1); ixp = np.clip(ix + 1, 0, nx - 1)
    gx = (material[iz, ixp] - material[iz, ixm]) / 2.0
    gz = (material[izp, ix] - material[izm, ix]) / 2.0
    norm = np.sqrt(gx ** 2 + gz ** 2)
    norm = np.where(norm < 1e-9, 1.0, norm)
    return -gx / norm, -gz / norm  # outward normal (solid -> empty)


def run_macro_step(grid, yield_model, ion_E, ion_theta, n_ions, n_neutrals,
                   sticking=0.35, removal_per_ion_hit=0.006, chem_removal_scale=0.004,
                   grazing_cutoff_deg=75.0, rng=None):
    """One macro-step: launch n_ions (sampled from the real IEDF/IADF arrays
    ion_E/ion_theta) and n_neutrals (cosine-distributed) from the trench
    opening, trace them, remove material, update radical coverage. Mutates
    grid.material and grid.radical in place."""
    rng = rng or np.random.default_rng()
    c0, c1 = grid.trench_cols
    nz, nx = grid.material.shape

    # --- ions ---------------------------------------------------------
    idx = rng.integers(0, ion_E.size, size=n_ions)
    E = ion_E[idx]
    theta_deg = ion_theta[idx]
    # random sign for lateral component (IADF is |angle|, symmetric about normal)
    sign = rng.choice([-1.0, 1.0], size=n_ions)
    theta_rad = np.radians(theta_deg) * sign
    x0 = rng.uniform(c0 + 1, c1 - 1, size=n_ions)
    z0 = np.zeros(n_ions)
    dirx = np.sin(theta_rad)
    dirz = np.cos(theta_rad)

    hit, hix, hiz = trace_batch(grid.material, x0, z0, dirx, dirz)
    n_reflected = 0
    if np.any(hit):
        hidx = np.nonzero(hit)[0]
        nx_, nz_ = local_normal(grid.material, hix[hidx], hiz[hidx])
        # angle of incidence from the surface normal
        cos_inc = np.clip(-(dirx[hidx] * nx_ + dirz[hidx] * nz_), -1, 1)
        inc_deg = np.degrees(np.arccos(np.abs(cos_inc)))
        grazing = inc_deg > grazing_cutoff_deg
        n_reflected = int(np.sum(grazing))

        # Non-grazing: remove material (physical sputter + chemical term)
        direct = hidx[~grazing[np.arange(hidx.size)]] if False else hidx[~grazing]
        if direct.size:
            Ei, thi = E[direct], np.abs(theta_rad[direct])
            Y = yield_model.physical_yield(Ei, thi)
            ir, ic = hiz[direct], hix[direct]
            radical_here = grid.radical[ir, ic]
            removal = removal_per_ion_hit * Y + chem_removal_scale * radical_here
            np.subtract.at(grid.material, (ir, ic), removal)
            np.subtract.at(grid.radical, (ir, ic), np.minimum(radical_here, 0.3 * radical_here + 0.05))

        # Grazing: specular reflect and re-trace once (this is the
        # microtrenching mechanism -- ion lands near the sidewall base).
        graze_idx = hidx[grazing]
        if graze_idx.size:
            gx, gz = local_normal(grid.material, hix[graze_idx], hiz[graze_idx])
            dot = dirx[graze_idx] * gx + dirz[graze_idx] * gz
            rdx = dirx[graze_idx] - 2 * dot * gx
            rdz = dirz[graze_idx] - 2 * dot * gz
            rx0 = hix[graze_idx].astype(np.float64) + rdx * 0.6
            rz0 = hiz[graze_idx].astype(np.float64) + rdz * 0.6
            hit2, hix2, hiz2 = trace_batch(grid.material, rx0, rz0, rdx, rdz)
            if np.any(hit2):
                h2 = np.nonzero(hit2)[0]
                Ei2 = E[graze_idx][h2]
                # reflected ions hit near-grazing on the new surface too;
                # approximate their local incidence as near the yield-boosting
                # angle (~70 deg) -- this cell is exactly where microtrenching
                # forms, so we still want it to preferentially remove material.
                Y2 = yield_model.physical_yield(Ei2, np.radians(70.0))
                np.subtract.at(grid.material, (hiz2[h2], hix2[h2]), 0.5 * removal_per_ion_hit * Y2)

    # --- neutrals -------------------------------------------------------
    u1, u2 = rng.random(n_neutrals), rng.random(n_neutrals)
    theta_n = np.arcsin(np.sqrt(u1))  # cosine-weighted polar angle
    sign_n = rng.choice([-1.0, 1.0], size=n_neutrals)
    theta_n = theta_n * sign_n
    xn0 = rng.uniform(c0 + 1, c1 - 1, size=n_neutrals)
    zn0 = np.zeros(n_neutrals)
    dxn = np.sin(theta_n); dzn = np.cos(theta_n)

    hitn, hnix, hniz = trace_batch(grid.material, xn0, zn0, dxn, dzn)
    if np.any(hitn):
        hidx = np.nonzero(hitn)[0]
        stick = rng.random(hidx.size) < sticking
        stick_idx = hidx[stick]
        if stick_idx.size:
            grid.radical[hniz[stick_idx], hnix[stick_idx]] = np.minimum(
                grid.radical[hniz[stick_idx], hnix[stick_idx]] + 0.15, 1.0)
        # non-sticking neutrals: one diffuse re-emission attempt (cosine
        # about the local normal), then given up (absorbed into "lost")
        bounce_idx = hidx[~stick]
        if bounce_idx.size:
            nx_, nz_ = local_normal(grid.material, hnix[bounce_idx], hniz[bounce_idx])
            # sample a cosine direction about the LOCAL normal (small-angle
            # approx: perturb normal by a random cosine-weighted angle)
            u = rng.random(bounce_idx.size)
            dphi = np.arcsin(np.sqrt(u)) * rng.choice([-1.0, 1.0], size=bounce_idx.size)
            # rotate normal by dphi (2D rotation)
            cs, sn = np.cos(dphi), np.sin(dphi)
            rdx = nx_ * cs - nz_ * sn
            rdz = nx_ * sn + nz_ * cs
            rx0 = hnix[bounce_idx].astype(np.float64) + rdx * 0.6
            rz0 = hniz[bounce_idx].astype(np.float64) + rdz * 0.6
            hit2, hix2, hiz2 = trace_batch(grid.material, rx0, rz0, rdx, rdz)
            if np.any(hit2):
                h2 = np.nonzero(hit2)[0]
                stick2 = rng.random(h2.size) < sticking
                s2 = h2[stick2]
                if s2.size:
                    grid.radical[hiz2[s2], hix2[s2]] = np.minimum(
                        grid.radical[hiz2[s2], hix2[s2]] + 0.15, 1.0)

    np.clip(grid.material, 0.0, 1.0, out=grid.material)
    grid.material[grid.is_mask] = 1.0  # mask never erodes (idealized hard mask)
    return dict(n_reflected=n_reflected)
