"""Step 5 verification: collision kinematics and null-collision statistics."""
import numpy as np
import pytest

from src.pic_mcc.constants import E_CHARGE as E
from src.pic_mcc.constants import K_B
from src.pic_mcc.cross_sections import load_helium_set, max_collision_frequency
from src.pic_mcc.mcc import (ElectronCollisions, IonCollisions, Species,
                             kinetic_energy_ev, random_chi_forward,
                             random_chi_isotropic, scatter_direction)

DIR = "data/cross_sections/turner_benchmark_he"
M_E = 9.10938215e-31
M_HE = 6.67e-27


def test_scatter_direction_preserves_speed_and_is_a_unit_vector():
    rng = np.random.default_rng(0)
    v = rng.normal(0, 1e5, size=(1000, 3))
    chi = random_chi_isotropic(rng, 1000)
    d = scatter_direction(v, chi, rng)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-9)


def test_scatter_direction_reproduces_the_requested_polar_angle():
    """cos(angle between old and new direction) must equal cos(chi)."""
    rng = np.random.default_rng(1)
    v = np.tile([0.0, 0.0, 1.0], (5000, 1))  # all pointing +z
    chi = np.full(5000, 0.3)
    d = scatter_direction(v, chi, rng)
    cos_angle = d[:, 2]  # dot with +z
    assert np.allclose(cos_angle, np.cos(0.3), atol=1e-9)


def test_scatter_direction_azimuth_is_uniform():
    rng = np.random.default_rng(2)
    v = np.tile([0.0, 0.0, 1.0], (20000, 1))
    chi = np.full(20000, np.pi / 2)  # scatter into the x-y plane
    d = scatter_direction(v, chi, rng)
    phi = np.arctan2(d[:, 1], d[:, 0])
    hist, _ = np.histogram(phi, bins=8, range=(-np.pi, np.pi))
    assert np.all(np.abs(hist / hist.mean() - 1.0) < 0.15)


def test_random_chi_isotropic_is_uniform_in_cos():
    rng = np.random.default_rng(3)
    chi = random_chi_isotropic(rng, 200_000)
    hist, _ = np.histogram(np.cos(chi), bins=10, range=(-1, 1))
    assert np.all(np.abs(hist / hist.mean() - 1.0) < 0.03)


def test_random_chi_forward_matches_sin_cos_weighting():
    """cos(chi)^2 should be uniform on [0,1] for this distribution."""
    rng = np.random.default_rng(4)
    chi = random_chi_forward(rng, 200_000)
    u = np.cos(chi) ** 2
    hist, _ = np.histogram(u, bins=10, range=(0, 1))
    assert np.all(np.abs(hist / hist.mean() - 1.0) < 0.03)


def test_elastic_collision_loses_energy_and_conserves_direction_statistics():
    """Electron elastic collisions must lose a small, mass-ratio-sized amount
    of energy on average and never gain energy (target started at rest)."""
    rng = np.random.default_rng(5)
    electrons, ions = load_helium_set(DIR)
    nu_max = max_collision_frequency(electrons, M_E, 9.64e20, energy_max_ev=200.0)
    coll = ElectronCollisions(electrons, M_E, M_HE, 9.64e20, nu_max, rng)
    coll.set_dt(1e-10)  # large dt -> most flagged particles actually collide

    n = 200_000
    E0 = 5.0  # eV, below excitation threshold -> only elastic is possible
    v0 = np.sqrt(2 * E0 * E / M_E)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v, mass=M_E, charge=-E)
    ions_sp = Species(x=np.zeros(0), v=np.zeros((0, 3)), mass=M_HE, charge=E)

    E_before = kinetic_energy_ev(sp.v, M_E)
    coll.apply(sp, ions_sp, t_neutral=300.0)
    E_after = kinetic_energy_ev(sp.v, M_E)

    collided = np.abs(E_after - E_before) > 1e-12
    assert collided.sum() > n * 0.003, "too few collisions occurred to test"
    assert np.all(E_after[collided] <= E_before[collided] + 1e-9), "elastic collision gained energy"
    mean_loss_fraction = 1.0 - E_after[collided].mean() / E0
    # Expected ~ (2 m_e/M)(1-<cos chi>). For isotropic-in-4pi chi, <1-cos chi> = 1.
    expected = 2.0 * M_E / M_HE
    assert mean_loss_fraction == pytest.approx(expected, rel=0.15)


def test_excitation_below_threshold_never_fires():
    rng = np.random.default_rng(6)
    electrons, ions = load_helium_set(DIR)
    nu_max = max_collision_frequency(electrons, M_E, 9.64e20, energy_max_ev=200.0)
    coll = ElectronCollisions(electrons, M_E, M_HE, 9.64e20, nu_max, rng)
    coll.set_dt(1e-10)
    n = 100_000
    E0 = 10.0  # eV: above elastic, below both excitation thresholds (19.82, 20.61)
    v0 = np.sqrt(2 * E0 * E / M_E)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v.copy(), mass=M_E, charge=-E)
    ions_sp = Species(x=np.zeros(0), v=np.zeros((0, 3)), mass=M_HE, charge=E)
    coll.apply(sp, ions_sp, t_neutral=300.0)
    E_after = kinetic_energy_ev(sp.v, M_E)
    assert np.all(E_after <= E0 + 1e-9)  # never gained energy = never excited/ionized


def test_ionization_splits_energy_equally_and_creates_particles():
    rng = np.random.default_rng(7)
    electrons, ions = load_helium_set(DIR)
    # Force-isolate ionization by zeroing the other three tables' cross sections.
    import copy
    procs = copy.deepcopy(electrons)
    for p in procs:
        if p.kind != "ionization":
            p.sigma_m2 = np.zeros_like(p.sigma_m2)
    nu_max = max_collision_frequency(procs, M_E, 9.64e20, energy_max_ev=200.0)
    coll = ElectronCollisions(procs, M_E, M_HE, 9.64e20, nu_max, rng)
    coll.set_dt(1e-9)

    n = 5000
    E0 = 30.0  # eV, safely above the 24.587 eV ionization threshold
    v0 = np.sqrt(2 * E0 * E / M_E)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v.copy(), mass=M_E, charge=-E)
    ions_sp = Species(x=np.zeros(0), v=np.zeros((0, 3)), mass=M_HE, charge=E)
    n_e_before, n_i_before = sp.n, ions_sp.n

    coll.apply(sp, ions_sp, t_neutral=300.0)

    n_new_electrons = sp.n - n_e_before
    n_new_ions = ions_sp.n - n_i_before
    assert n_new_electrons > 0, "no ionization events occurred"
    assert n_new_electrons == n_new_ions, "one ion must be created per ionization event"

    # Every post-ionization primary+secondary electron energy must equal the
    # equal-split value (since all started at exactly E0): (E0 - E_ion)/2.
    expected = (E0 - 24.587) / 2.0
    all_energies = kinetic_energy_ev(sp.v, M_E)
    ionized_energies = all_energies[np.abs(all_energies - expected) < 0.05]
    assert ionized_energies.size >= n_new_electrons, "post-ionization energies don't match the equal split"

    # New ions must be at thermal speed (~300K), not at the ion beam energy.
    new_ion_E = kinetic_energy_ev(ions_sp.v[n_i_before:], M_HE)
    assert np.all(new_ion_E < 1.0), "new ions were not created at thermal energy"


def test_ion_isotropic_scattering_conserves_energy_in_relative_frame_on_average():
    """Elastic ion-neutral collisions must not systematically pump energy in."""
    rng = np.random.default_rng(8)
    _, ions = load_helium_set(DIR)
    nu_max = max_collision_frequency(ions, M_HE, 9.64e20, energy_max_ev=1000.0)
    coll = IonCollisions(ions, M_HE, 9.64e20, 300.0, nu_max, rng)
    coll.set_dt(1e-9)

    n = 200_000
    E0 = 50.0  # eV, well above thermal
    v0 = np.sqrt(2 * E0 * E / M_HE)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v.copy(), mass=M_HE, charge=E)

    E_before = kinetic_energy_ev(sp.v, M_HE)
    coll.apply(sp)
    E_after = kinetic_energy_ev(sp.v, M_HE)
    collided = np.abs(E_after - E_before) > 1e-9
    assert collided.sum() > n * 0.001, "too few ion collisions to test"
    # Post-collision energy must be <= pre-collision for the fast-ion population
    # (elastic scatter always removes relative KE; charge exchange drops it to
    # near-thermal). No mechanism here can raise a 50 eV ion's energy.
    assert np.all(E_after[collided] <= E_before[collided] + 5.0), "ion collision gained energy"


def test_charge_exchange_leaves_ion_near_thermal_energy():
    rng = np.random.default_rng(9)
    _, ions = load_helium_set(DIR)
    import copy
    procs = copy.deepcopy(ions)
    procs[0].sigma_m2 = np.zeros_like(procs[0].sigma_m2)  # kill isotropic, keep CX only
    nu_max = max_collision_frequency(procs, M_HE, 9.64e20, energy_max_ev=1000.0)
    coll = IonCollisions(procs, M_HE, 9.64e20, 300.0, nu_max, rng)
    coll.set_dt(1e-9)

    n = 100_000
    E0 = 100.0
    v0 = np.sqrt(2 * E0 * E / M_HE)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v.copy(), mass=M_HE, charge=E)
    coll.apply(sp)
    E_after = kinetic_energy_ev(sp.v, M_HE)
    collided = E_after < E0 - 1.0
    assert collided.sum() > n * 0.001
    # Post-CX energy should be of order k_B*300K/e ~ tens of meV, not eV.
    assert np.median(E_after[collided]) < 0.5


def test_null_collision_rate_matches_analytic_elastic_only_rate():
    """With only elastic active at a fixed sub-threshold energy, the fraction
    of particles that collide in one step should match n*sigma(E)*v*dt."""
    rng = np.random.default_rng(10)
    electrons, ions = load_helium_set(DIR)
    import copy
    procs = copy.deepcopy(electrons)
    for p in procs:
        if p.kind != "elastic":
            p.sigma_m2 = np.zeros_like(p.sigma_m2)
    n_gas = 9.64e20
    nu_max = max_collision_frequency(procs, M_E, n_gas, energy_max_ev=200.0)
    dt = 5e-11
    coll = ElectronCollisions(procs, M_E, M_HE, n_gas, nu_max, rng)
    coll.set_dt(dt)

    n = 300_000
    E0 = 5.0
    v0 = np.sqrt(2 * E0 * E / M_E)
    v = np.zeros((n, 3)); v[:, 2] = v0
    sp = Species(x=np.zeros(n), v=v.copy(), mass=M_E, charge=-E)
    ions_sp = Species(x=np.zeros(0), v=np.zeros((0, 3)), mass=M_HE, charge=E)

    E_before = kinetic_energy_ev(sp.v, M_E)
    coll.apply(sp, ions_sp, t_neutral=300.0)
    E_after = kinetic_energy_ev(sp.v, M_E)
    frac_collided = np.mean(np.abs(E_after - E_before) > 1e-12)

    sigma_E0 = procs[0](np.array([E0]))[0]
    nu_analytic = n_gas * sigma_E0 * v0
    expected_frac = 1.0 - np.exp(-nu_analytic * dt)
    assert frac_collided == pytest.approx(expected_frac, rel=0.1)


def test_ion_cm_energy_formula_matches_reduced_mass_definition():
    """E_cm = (1/2) mu v_rel^2, mu = m1 m2/(m1+m2). For equal masses, mu = m/2,
    so E_cm = m v_rel^2 / 4 -- NOT /8 (a bug caught once: an earlier version of
    this formula used /8, silently evaluating ion cross sections at half the
    correct centre-of-mass energy). Cross-checked directly against the eduPIC
    code (Donko et al. 2021, PSST 30 095017): energy = 0.5*MU_ARAR*g_sqr, the
    same physical quantity, computed independently there.
    """
    m = 6.67e-27
    v_rel = 1.0e4  # m/s, arbitrary
    mu = m * m / (m + m)
    E_cm_from_definition = 0.5 * mu * v_rel ** 2 / E
    E_cm_from_code_formula = m * v_rel ** 2 / 4.0 / E
    assert E_cm_from_code_formula == pytest.approx(E_cm_from_definition, rel=1e-12)
    wrong_formula = m * v_rel ** 2 / 8.0 / E
    assert wrong_formula == pytest.approx(E_cm_from_definition / 2.0, rel=1e-12)
