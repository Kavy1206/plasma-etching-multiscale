"""Step 2 + full-core verification: cold plasma oscillation.

Smaller and shorter than scripts/verify_core.py so it runs in a few seconds,
but the same physics. See that script's docstring for why the second harmonic
and why the measurement window is bounded.
"""
import numpy as np
import pytest

from src.pic_mcc.constants import EPS0, E_CHARGE, M_E
from src.pic_mcc.grid import Grid1D
from src.pic_mcc.poisson import PoissonSolver1D, electric_field
from src.pic_mcc.pusher import half_step_back, push
from scripts.verify_core import _omega_from_zero_crossings

L, N_CELLS, N_DENSITY, AMP = 6.7e-2, 128, 1.363e14, 0.002


def _run(n_steps=1200, ppc=128, omega_dt=0.1):
    grid = Grid1D(L, N_CELLS)
    solver = PoissonSolver1D(N_CELLS, grid.dx)
    omega_pe = np.sqrt(N_DENSITY * E_CHARGE ** 2 / (EPS0 * M_E))
    dt = omega_dt / omega_pe
    n = N_CELLS * ppc
    weight = N_DENSITY * L / n
    x0 = (np.arange(n) + 0.5) * (L / n)
    shape = np.sin(2 * np.pi * x0 / L)
    x = x0 + AMP * L * shape
    vx = np.zeros(n)
    qm, rb = -E_CHARGE / M_E, E_CHARGE * N_DENSITY

    def F(xp):
        rho = rb - E_CHARGE * grid.deposit(xp, weight) / grid.node_width
        return grid.gather(electric_field(solver.solve(rho, 0.0, 0.0), grid.dx), xp)

    vx = half_step_back(vx, F(x), qm, dt)
    amp = np.empty(n_steps)
    for i in range(n_steps):
        amp[i] = np.sum((x - x0) * shape)
        push(x, vx, F(x), qm, dt)
    return amp / (AMP * L * np.sum(shape ** 2)), dt, omega_pe


def test_induced_field_matches_the_analytic_cold_slab_result():
    """E = e n0 xi / eps0 for a displaced cold slab, with no spurious offset."""
    grid = Grid1D(L, N_CELLS)
    solver = PoissonSolver1D(N_CELLS, grid.dx)
    n = N_CELLS * 512
    x0 = (np.arange(n) + 0.5) * (L / n)
    shape = np.sin(2 * np.pi * x0 / L)
    x = x0 + AMP * L * shape
    rho = E_CHARGE * N_DENSITY - E_CHARGE * grid.deposit(x, N_DENSITY * L / n) / grid.node_width
    E = electric_field(solver.solve(rho, 0.0, 0.0), grid.dx)
    E_expect = E_CHARGE * N_DENSITY * AMP * L * np.sin(2 * np.pi * grid.x / L) / EPS0
    scale = np.max(np.abs(E_expect))
    assert abs(np.mean(E - E_expect)) < 1e-9, "spurious uniform field from the boundary condition"
    # Electrode nodes: second-order one-sided stencil, essentially exact.
    assert abs(E[0] - E_expect[0]) / scale < 1e-4
    assert abs(E[-1] - E_expect[-1]) / scale < 1e-4
    # Bulk: dominated by CIC + finite-particle noise, not by the stencil.
    assert np.max(np.abs(E[1:-1] - E_expect[1:-1])) / scale < 1e-2


def test_oscillates_at_the_plasma_frequency():
    amp, dt, omega_pe = _run()
    omega = _omega_from_zero_crossings(amp, dt)
    omega_leapfrog = (2.0 / dt) * np.arcsin(omega_pe * dt / 2.0)
    # Must sit between omega_pe and the leapfrog-dispersion value, close to both.
    assert omega == pytest.approx(omega_pe, rel=1e-3)
    assert omega == pytest.approx(omega_leapfrog, rel=1e-3)


def test_frequency_is_independent_of_timestep():
    """Halving dt must move the answer toward omega_pe by ~4x (second order)."""
    _, dt1, omega_pe = _run(n_steps=600, omega_dt=0.2)
    a1, _, _ = _run(n_steps=600, omega_dt=0.2)
    a2, dt2, _ = _run(n_steps=1200, omega_dt=0.1)
    e1 = abs(_omega_from_zero_crossings(a1, dt1) / omega_pe - 1)
    e2 = abs(_omega_from_zero_crossings(a2, dt2) / omega_pe - 1)
    assert e2 < e1, f"refining dt did not improve the frequency: {e1:.2e} -> {e2:.2e}"
