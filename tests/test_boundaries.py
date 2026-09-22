"""Step 6 verification: wall absorption and secondary emission."""
import numpy as np

from src.pic_mcc.boundaries import apply_absorbing_walls, emit_secondaries
from src.pic_mcc.mcc import Species


def test_absorbing_walls_remove_only_out_of_bounds_particles():
    n = 1000
    rng = np.random.default_rng(0)
    x = rng.uniform(-0.5, 1.5, n)  # L = 1.0
    sp = Species(x=x, v=np.zeros((n, 3)), mass=1.0, charge=1.0)
    n_before = sp.n
    wall, ax, av = apply_absorbing_walls(sp, 1.0)
    assert np.all((sp.x >= 0.0) & (sp.x <= 1.0))
    assert sp.n + wall.size == n_before
    assert np.all(ax[wall == 0] < 0.0)
    assert np.all(ax[wall == 1] > 1.0)


def test_gamma_zero_emits_nothing():
    rng = np.random.default_rng(1)
    wall = np.zeros(500, dtype=int)
    x_new, v_new = emit_secondaries(wall, 0.0, 2.0, 1.0, 1.0, rng)
    assert x_new.size == 0 and v_new.shape == (0, 3)


def test_gamma_positive_rate_matches_gamma():
    rng = np.random.default_rng(2)
    n = 200_000
    wall = rng.integers(0, 2, n)
    x_new, v_new = emit_secondaries(wall, 0.08, 2.0, 9.109e-31, 0.067, rng)
    assert x_new.size == pytest_approx_frac(n, 0.08)


def pytest_approx_frac(n, frac, tol=0.01):
    import pytest
    return pytest.approx(n * frac, rel=tol * 10)


def test_secondaries_launch_inward_from_each_wall():
    rng = np.random.default_rng(3)
    n = 5000
    wall = np.array([0] * (n // 2) + [1] * (n // 2))
    L = 0.067
    x_new, v_new = emit_secondaries(wall, 1.0, 2.0, 9.109e-31, L, rng)
    left_emitted = v_new[: n // 2]
    right_emitted = v_new[n // 2:]
    assert np.all(left_emitted[:, 0] > 0), "left-wall secondaries must move in +x"
    assert np.all(right_emitted[:, 0] < 0), "right-wall secondaries must move in -x"
    assert np.all(x_new[: n // 2] < 1e-3)
    assert np.all(x_new[n // 2:] > L - 1e-3)


def test_secondary_speed_matches_the_declared_energy():
    from src.pic_mcc.constants import E_CHARGE as E
    rng = np.random.default_rng(4)
    wall = np.zeros(50000, dtype=int)
    mass = 9.109e-31
    T = 2.0
    x_new, v_new = emit_secondaries(wall, 1.0, T, mass, 1.0, rng)
    speed = np.linalg.norm(v_new, axis=1)
    E_meas = 0.5 * mass * speed ** 2 / E
    assert np.allclose(E_meas, T, rtol=1e-6)
