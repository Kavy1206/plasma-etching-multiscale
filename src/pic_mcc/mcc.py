"""Null-collision Monte Carlo collisions (MCC), step 5.

Electron-neutral: elastic, two excitations, ionization. Target treated as
stationary (standard approximation: v_th,He << v_th,e by ~sqrt(M/m_e), so the
neutral's own motion is negligible next to the electron's). Post-collision
speeds follow the standard formulas (matching Vahedi & Surendra, Comput. Phys.
Commun. 87, 179 (1995), and cross-checked against the algorithm used by the
public Turner-benchmark reimplementation lase-unb/ccp-benchmark, whose
scattering.cpp I read for exact formulas -- I did not vendor or copy that
code; the implementation below is independent, in Python, with its own tests):

  elastic:     E' = E * (1 - (2 m_e/M)(1-cos(chi)))      [energy loss to recoil]
  excitation:  E' = E - E_exc
  ionization:  E_share = (E - E_ion) / 2, split equally between the
               scattered primary and the ejected secondary (both isotropic)

Ion-neutral: isotropic ("elastic") and backscattering ("charge exchange").
Unlike electrons, the neutral's thermal motion is NOT negligible next to a
slow ion, so for ion collisions we sample a target velocity from a Maxwellian
at T_gas for each flagged candidate and use the ion-target RELATIVE velocity
for both cross-section lookup (in the centre-of-mass frame the benchmark's
tables are defined in) and post-collision kinematics, then transform back to
the lab frame. This is the "slow projectile" treatment; electrons use the
"fast projectile" (stationary-target) treatment. Two physically distinct
processes:

  isotropic ("Isotropic_He.csv"):  genuine elastic hard-sphere scattering.
      For EQUAL ion/neutral masses, isotropic-in-CM elastic scattering gives,
      in the lab frame (target initially at the sampled thermal velocity),
      outgoing relative speed v_rel' = v_rel * cos(chi), chi ~ arccos(sqrt(1-R))
      (the sin(chi)cos(chi) solid-angle weighting for a hard-sphere CM
      distribution that is isotropic in the FORWARD hemisphere -- see
      random_chi2 below).
  backscattering ("Backscattering_He.csv"):  in Phelps' He+-He compilation
      this channel represents resonant-charge-transfer-like backscatter: the
      fast ion effectively becomes a slow neutral and a new "ion" appears at
      the target's (thermal) velocity. Implemented by setting the ion's
      post-collision velocity equal to the sampled target velocity.

Null-collision method (Skullerud 1968; Vahedi & Surendra 1995): precompute
nu_max = n_gas * max_E[ sigma_total(E) v(E) ]. Each step, flag a particle as
a collision candidate with probability 1 - exp(-nu_max dt). For flagged
candidates only, evaluate the true local nu(E) and accept with probability
nu(E)/nu_max (a "null" collision does nothing) -- this reproduces the correct
energy-dependent collision rate while using one bulk Bernoulli draw plus a
much smaller local evaluation, rather than evaluating cross sections for
every particle every step.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import E_CHARGE as E
from .constants import K_B, M_E


def random_chi_isotropic(rng: np.random.Generator, n: int) -> np.ndarray:
    """Uniform over the full sphere: cos(chi) ~ Uniform(-1, 1)."""
    return np.arccos(1.0 - 2.0 * rng.random(n))


def random_chi_forward(rng: np.random.Generator, n: int) -> np.ndarray:
    """Hard-sphere elastic CM angle: cos(chi) = sqrt(1 - R), R ~ Uniform(0,1).

    This is the standard distribution for elastic scattering off a hard
    sphere of a Maxwellian-distributed target and is what the equal-mass
    v' = v cos(chi) formula below assumes.
    """
    return np.arccos(np.sqrt(1.0 - rng.random(n)))


def scatter_direction(v: np.ndarray, chi: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Rotate each row-vector in v by polar angle chi about a random azimuth.

    v: (n, 3) array of (unnormalised) direction vectors (only direction is used).
    Returns unit vectors. Standard PIC-MCC rotation (e.g. Birdsall & Langdon
    sec. 14-3): build an orthonormal frame from the incoming direction and
    rotate into it.
    """
    speed = np.linalg.norm(v, axis=1)
    speed = np.where(speed > 0, speed, 1.0)
    u = v / speed[:, None]
    n = v.shape[0]
    phi = 2.0 * np.pi * rng.random(n)

    ux, uy, uz = u[:, 0], u[:, 1], u[:, 2]
    # Guard the pole: if |uz| ~ 1, build the perpendicular frame from x instead of z.
    pole = np.abs(uz) > 0.9999
    perp = np.empty_like(u)
    perp[~pole] = np.stack([uy[~pole], -ux[~pole], np.zeros(np.count_nonzero(~pole))], axis=1)
    perp[pole] = np.stack([np.ones(np.count_nonzero(pole)), np.zeros(np.count_nonzero(pole)),
                           np.zeros(np.count_nonzero(pole))], axis=1)
    perp /= np.linalg.norm(perp, axis=1)[:, None]
    w = np.cross(u, perp)

    cos_chi = np.cos(chi)[:, None]
    sin_chi = np.sin(chi)[:, None]
    cos_phi = np.cos(phi)[:, None]
    sin_phi = np.sin(phi)[:, None]

    return u * cos_chi + perp * (sin_chi * cos_phi) + w * (sin_chi * sin_phi)


@dataclass
class Species:
    """Minimal particle container: positions (m), velocities (m/s, n x 3)."""
    x: np.ndarray
    v: np.ndarray
    mass: float
    charge: float

    @property
    def n(self) -> int:
        return self.x.shape[0]

    def append(self, x_new: np.ndarray, v_new: np.ndarray) -> None:
        self.x = np.concatenate([self.x, x_new])
        self.v = np.concatenate([self.v, v_new])

    def remove(self, keep_mask: np.ndarray) -> None:
        self.x = self.x[keep_mask]
        self.v = self.v[keep_mask]


def kinetic_energy_ev(v: np.ndarray, mass: float) -> np.ndarray:
    return 0.5 * mass * np.sum(v * v, axis=1) / E


class ElectronCollisions:
    """Elastic + 2x excitation + ionization, target at rest."""

    def __init__(self, processes, mass: float, ion_mass: float, n_gas: float,
                 nu_max: float, rng: np.random.Generator):
        self.processes = processes  # [elastic, exc1, exc2, ionization], in this order
        self.mass = mass
        self.ion_mass = ion_mass
        self.n_gas = n_gas
        self.nu_max = nu_max
        self.rng = rng
        self.p_null = None  # set once dt is known

    def set_dt(self, dt: float) -> None:
        self.p_null = 1.0 - np.exp(-self.nu_max * dt)

    def apply(self, sp: Species, ions: Species, t_neutral: float):
        """Mutates sp in place; may append new electrons and ions (ionization)."""
        n = sp.n
        if n == 0:
            return
        flagged = self.rng.random(n) < self.p_null
        idx = np.nonzero(flagged)[0]
        if idx.size == 0:
            return

        v = sp.v[idx]
        E_ev = kinetic_energy_ev(v, self.mass)
        sig = np.stack([p(E_ev) for p in self.processes], axis=1)  # (m, 4)
        sig_total = sig.sum(axis=1)
        speed = np.sqrt(2.0 * E_ev * E / self.mass)
        nu_local = self.n_gas * sig_total * speed
        accept = self.rng.random(idx.size) < (nu_local / self.nu_max)
        idx = idx[accept]
        if idx.size == 0:
            return
        E_ev = E_ev[accept]
        sig = sig[accept]
        sig_total = sig_total[accept]

        # Pick a process per accepted particle, weighted by its own cross section.
        cum = np.cumsum(sig, axis=1) / sig_total[:, None]
        r = self.rng.random(idx.size)
        proc_id = np.sum(r[:, None] > cum, axis=1)  # 0..3

        m = idx.size
        v_dir_scale = np.empty(m)          # new speed magnitude
        chi = random_chi_isotropic(self.rng, m)
        new_ion_x = []
        new_ion_v = []
        new_e_x = []
        new_e_v = []

        for k, proc in enumerate(self.processes):
            sel = proc_id == k
            if not np.any(sel):
                continue
            Ek = E_ev[sel]
            if proc.kind == "elastic":
                delta = (2.0 * self.mass / self.ion_mass) * (1.0 - np.cos(chi[sel]))
                v_dir_scale[sel] = np.sqrt(np.maximum(2.0 * E * Ek * (1.0 - delta), 0.0) / self.mass)
            elif proc.kind == "excitation":
                v_dir_scale[sel] = np.sqrt(np.maximum(2.0 * E * (Ek - proc.threshold_ev), 0.0) / self.mass)
            elif proc.kind == "ionization":
                share = np.maximum(Ek - proc.threshold_ev, 0.0) / 2.0
                v_dir_scale[sel] = np.sqrt(2.0 * E * share / self.mass)
                # Secondary electron: isotropic direction, same speed formula (equal share).
                sec_idx = idx[sel]
                chi2 = random_chi_isotropic(self.rng, sec_idx.size)
                dirs2 = scatter_direction(sp.v[sec_idx], chi2, self.rng)
                v_sec = dirs2 * v_dir_scale[sel][:, None]
                new_e_x.append(sp.x[sec_idx].copy())
                new_e_v.append(v_sec)
                # New ion at (near-)thermal speed, at the ionization site.
                v_th = np.sqrt(K_B * t_neutral / self.ion_mass)
                new_ion_x.append(sp.x[sec_idx].copy())
                new_ion_v.append(self.rng.normal(0.0, v_th, size=(sec_idx.size, 3)))

        dirs = scatter_direction(sp.v[idx], chi, self.rng)
        sp.v[idx] = dirs * v_dir_scale[:, None]

        if new_e_x:
            sp.append(np.concatenate(new_e_x), np.concatenate(new_e_v))
        if new_ion_x:
            ions.append(np.concatenate(new_ion_x), np.concatenate(new_ion_v))


class IonCollisions:
    """Isotropic elastic + backscattering (charge-exchange-like), thermal target."""

    def __init__(self, processes, mass: float, n_gas: float, t_neutral: float,
                 nu_max: float, rng: np.random.Generator):
        self.processes = processes  # [isotropic, backscatter]
        self.mass = mass
        self.n_gas = n_gas
        self.t_neutral = t_neutral
        self.nu_max = nu_max
        self.rng = rng
        self.p_null = None

    def set_dt(self, dt: float) -> None:
        self.p_null = 1.0 - np.exp(-self.nu_max * dt)

    def apply(self, sp: Species):
        n = sp.n
        if n == 0:
            return
        flagged = self.rng.random(n) < self.p_null
        idx = np.nonzero(flagged)[0]
        if idx.size == 0:
            return

        m = idx.size
        v_th = np.sqrt(K_B * self.t_neutral / self.mass)
        v_target = self.rng.normal(0.0, v_th, size=(m, 3))
        v_rel = sp.v[idx] - v_target
        speed_rel = np.linalg.norm(v_rel, axis=1)
        # Ion tables are indexed by CENTRE-OF-MASS energy: E_cm = (1/2) mu v_rel^2
        # with reduced mass mu = m1 m2/(m1+m2). For equal ion/neutral masses,
        # mu = m/2, so E_cm = m v_rel^2 / 4. (Cross-checked directly against the
        # eduPIC code (Donko et al. 2021, PSST 30 095017), which computes the
        # identical quantity for the same physical setup as
        # energy = 0.5 * MU_ARAR * g_sqr -- i.e. (1/2) mu v_rel^2, confirming the
        # factor. A first version of this line used m*v_rel^2/8, i.e. HALF the
        # correct CM energy -- caught during argon cross-section sourcing, see
        # PROGRESS.md, and the full Case 1 validation run was redone after the fix.)
        E_cm = self.mass * speed_rel ** 2 / 4.0 / E

        sig = np.stack([p(E_cm) for p in self.processes], axis=1)
        sig_total = sig.sum(axis=1)
        nu_local = self.n_gas * sig_total * speed_rel
        accept = self.rng.random(m) < (nu_local / self.nu_max)
        idx = idx[accept]
        if idx.size == 0:
            return
        v_rel = v_rel[accept]
        v_target = v_target[accept]
        speed_rel = speed_rel[accept]
        sig = sig[accept]
        sig_total = sig_total[accept]

        cum = np.cumsum(sig, axis=1) / sig_total[:, None]
        r = self.rng.random(idx.size)
        proc_id = np.sum(r[:, None] > cum, axis=1)  # 0=isotropic, 1=backscatter

        elastic_sel = proc_id == 0
        cx_sel = proc_id == 1

        if np.any(elastic_sel):
            n_e = int(np.count_nonzero(elastic_sel))
            chi = random_chi_forward(self.rng, n_e)
            dirs = scatter_direction(v_rel[elastic_sel], chi, self.rng)
            v_rel_new = dirs * (speed_rel[elastic_sel] * np.cos(chi))[:, None]
            sp.v[idx[elastic_sel]] = v_rel_new + v_target[elastic_sel]

        if np.any(cx_sel):
            # Fast ion -> slow neutral; a new "ion" appears at the target velocity.
            sp.v[idx[cx_sel]] = v_target[cx_sel]
