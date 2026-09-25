"""Phase 0 data check: the Turner helium cross sections load correctly."""
import numpy as np
import pytest

from src.pic_mcc.cross_sections import (lab_to_cm_energy, load_helium_set,
                                        max_collision_frequency)

DIR = "data/cross_sections/turner_benchmark_he"


def test_all_six_tables_load():
    electrons, ions = load_helium_set(DIR)
    assert len(electrons) == 4
    assert len(ions) == 2


def test_thresholds_match_the_tables():
    """The tabulated onset energy must agree with the documented threshold."""
    electrons, _ = load_helium_set(DIR)
    expected = {"e-He -> He*": 19.82, "e-He -> He**": 20.61, "e-He ionization": 24.587}
    for xs in electrons:
        if xs.name in expected:
            first_nonzero = xs.energy_ev[np.argmax(xs.sigma_m2 > 0)]
            assert first_nonzero == pytest.approx(expected[xs.name], abs=0.05), xs.name


def test_inelastic_cross_sections_vanish_below_threshold():
    electrons, _ = load_helium_set(DIR)
    for xs in electrons:
        if xs.kind in ("excitation", "ionization"):
            below = np.array([0.1, 1.0, 5.0, 10.0, xs.threshold_ev - 0.5])
            assert np.all(xs(below) == 0.0), xs.name


def test_clamping_above_the_last_tabulated_point():
    electrons, _ = load_helium_set(DIR)
    xs = electrons[0]
    assert xs(np.array([1e9]))[0] == pytest.approx(xs.sigma_m2[-1])


def test_ion_tables_are_flagged_as_centre_of_mass():
    _, ions = load_helium_set(DIR)
    assert all(xs.frame == "cm" for xs in ions)


def test_lab_to_cm_is_a_half_for_equal_masses():
    assert lab_to_cm_energy(np.array([100.0]), 6.67e-27, 6.67e-27)[0] == pytest.approx(50.0)


def test_null_collision_frequency_is_finite_and_positive():
    from src.pic_mcc.config import turner_case
    cfg = turner_case(1)
    electrons, ions = load_helium_set(DIR)
    nu_e = max_collision_frequency(electrons, cfg.m_electron, cfg.n_gas)
    nu_i = max_collision_frequency(ions, cfg.m_ion, cfg.n_gas, energy_max_ev=1000.0)
    assert 0 < nu_e < 1e12 and 0 < nu_i < 1e12
    # Benchmark requires nu*dt << 1 (their Table III: nu_e dt = 0.0158)
    assert nu_e * cfg.dt < 0.2, f"nu_e*dt = {nu_e * cfg.dt:.4f}"
    assert nu_i * cfg.dt < 0.2, f"nu_i*dt = {nu_i * cfg.dt:.4f}"


def test_argon_set_loads_and_thresholds_are_correct():
    from src.pic_mcc.cross_sections import load_gas_set
    electrons, ions = load_gas_set("ar", "data/cross_sections/ar")
    assert len(electrons) == 3 and len(ions) == 2
    names = {p.name: p for p in electrons}
    exc = electrons[[p.kind for p in electrons].index("excitation")]
    ion = electrons[[p.kind for p in electrons].index("ionization")]
    below = np.array([1.0, 5.0, 11.0])
    assert np.all(exc(below) == 0.0)
    assert np.all(ion(np.array([1.0, 5.0, 15.0])) == 0.0)
    assert exc(np.array([20.0]))[0] > 0
    assert ion(np.array([50.0]))[0] > 0


def test_argon_ion_backscatter_is_nonnegative():
    from src.pic_mcc.cross_sections import load_gas_set
    _, ions = load_gas_set("ar", "data/cross_sections/ar")
    back = ions[[p.kind for p in ions].index("backscatter")]
    assert np.all(back(np.linspace(0.01, 2000, 5000)) >= 0.0)


def test_argon_simulation_config_uses_argon_data_dir():
    from src.pic_mcc.config import Config
    cfg = Config(gas="ar", gamma_see=0.05)
    assert "ar" in cfg.cross_section_dir
    assert cfg.m_ion == pytest.approx(39.948 * 1.660538782e-27, rel=1e-6)
