"""Step 4 verification: the Poisson solve."""
import numpy as np
import pytest

from src.pic_mcc.constants import EPS0
from src.pic_mcc.poisson import PoissonSolver1D, electric_field


def _manufactured(n_cells, L=6.7e-2, rho0=1e-6):
    """rho(x) = rho0 sin(pi x / L)  =>  phi(x) = rho0 L^2 /(eps0 pi^2) sin(pi x/L)

    with phi(0) = phi(L) = 0. Check by substitution:
      d2phi/dx2 = -(pi/L)^2 phi = -rho0 sin(pi x/L)/eps0 = -rho/eps0.  OK
    """
    dx = L / n_cells
    x = np.linspace(0, L, n_cells + 1)
    rho = rho0 * np.sin(np.pi * x / L)
    phi_exact = rho0 * L ** 2 / (EPS0 * np.pi ** 2) * np.sin(np.pi * x / L)
    return dx, x, rho, phi_exact


def test_dirichlet_boundaries_are_exact():
    solver = PoissonSolver1D(64, 1.0 / 64)
    rho = np.zeros(65)
    phi = solver.solve(rho, phi_left=450.0, phi_right=0.0)
    assert phi[0] == pytest.approx(450.0)
    assert phi[-1] == pytest.approx(0.0)
    # Vacuum between charged plates -> linear potential, uniform field
    assert np.allclose(phi, np.linspace(450.0, 0.0, 65))
    E = electric_field(phi, 1.0 / 64)
    assert np.allclose(E, 450.0, rtol=1e-10)


def test_second_order_convergence():
    errors = {}
    for n in (32, 64, 128, 256):
        dx, x, rho, phi_exact = _manufactured(n)
        phi = PoissonSolver1D(n, dx).solve(rho, 0.0, 0.0)
        errors[n] = np.max(np.abs(phi - phi_exact)) / np.max(np.abs(phi_exact))
    # Halving dx must cut the error by ~4 (second-order stencil)
    for coarse, fine in ((32, 64), (64, 128), (128, 256)):
        ratio = errors[coarse] / errors[fine]
        assert 3.7 < ratio < 4.3, f"order check failed: {coarse}->{fine} ratio {ratio:.2f}"
    # Stronger check: the discrete Laplacian has eigenvalue -(4/dx^2) sin^2(k dx/2)
    # instead of -k^2, so for this manufactured solution the relative error is
    # predicted exactly, to leading order, by (k dx)^2 / 12. Matching that says
    # the only error left is stencil truncation -- no algebra or sign bugs.
    for n, err in errors.items():
        k_dx = np.pi / n
        assert err == pytest.approx(k_dx ** 2 / 12.0, rel=0.01), f"n={n}"


def test_matches_analytic_solution():
    dx, x, rho, phi_exact = _manufactured(128)
    phi = PoissonSolver1D(128, dx).solve(rho, 0.0, 0.0)
    assert np.max(np.abs(phi - phi_exact)) / np.max(np.abs(phi_exact)) < 1e-4
