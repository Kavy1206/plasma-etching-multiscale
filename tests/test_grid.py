"""Steps 1 and 3 verification: grid and cloud-in-cell deposition."""
import numpy as np
import pytest

from src.pic_mcc.grid import Grid1D


def test_debye_resolution_of_benchmark_case_1():
    """Turner Case 1 reports lambda_D/dx = 3.72 at the mid-plane."""
    from src.pic_mcc.config import turner_case
    from src.pic_mcc.constants import EPS0, E_CHARGE
    cfg = turner_case(1)
    n_e, T_e = 1.363e14, 9.36  # Table III mid-plane values
    lam = np.sqrt(EPS0 * T_e / (E_CHARGE * n_e))
    assert lam / cfg.dx == pytest.approx(3.72, rel=0.01)


def test_deposition_conserves_charge():
    g = Grid1D(1.0, 32)
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1.0, 10_000)
    acc = g.deposit(x, weight=1.0)
    assert acc.sum() == pytest.approx(10_000.0, rel=1e-12)


def test_particle_on_node_deposits_entirely_there():
    g = Grid1D(1.0, 10)
    acc = g.deposit(np.array([0.3]), weight=1.0)
    assert acc[3] == pytest.approx(1.0)
    assert acc.sum() == pytest.approx(1.0)


def test_particle_at_cell_centre_splits_evenly():
    g = Grid1D(1.0, 10)
    acc = g.deposit(np.array([0.35]), weight=1.0)
    assert acc[3] == pytest.approx(0.5)
    assert acc[4] == pytest.approx(0.5)


def test_boundary_nodes_use_half_width():
    """A uniform particle distribution must give a UNIFORM density, including
    at the walls. This only works if the boundary nodes are divided by dx/2."""
    g = Grid1D(1.0, 64)
    x = np.linspace(0, 1.0, 640_001)  # uniform, endpoints included
    weight = 1.0 / 640_000            # so total weight = 1
    rho = g.deposit(x, weight) / g.node_width
    assert np.allclose(rho, rho[32], rtol=2e-3), "density is not flat at the walls"


def test_gather_and_scatter_use_the_same_weights():
    """Gathering a field that is 1 everywhere must return exactly 1 for every
    particle -- the partition-of-unity property that kills the self-force."""
    g = Grid1D(1.0, 32)
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1.0, 5_000)
    ones = np.ones(g.n_nodes)
    assert np.allclose(g.gather(ones, x), 1.0)


def test_gather_is_exact_for_a_linear_field():
    g = Grid1D(1.0, 32)
    field = 3.0 * g.x + 7.0
    x = np.array([0.0, 0.123, 0.5, 0.987, 1.0])
    assert np.allclose(g.gather(field, x), 3.0 * x + 7.0)
