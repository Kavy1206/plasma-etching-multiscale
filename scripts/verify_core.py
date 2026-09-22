"""Integration check for the Phase 1 electrostatic core (steps 1-4).

The cold plasma oscillation is the standard PIC acceptance test and it is the
right one here because it exercises every piece at once: deposit -> Poisson ->
gather -> push -> deposit. If any of the four has a sign error, a wrong weight,
or a stray factor of dx, the measured frequency will not come out at omega_pe.

Setup: uniform immobile positive background, cold electrons displaced by
    xi(x) = A sin(2 pi x / L)

WHY THE SECOND HARMONIC AND NOT THE FIRST. Cold electrostatic oscillations are
non-propagating (omega = omega_pe for every k), so physically any mode would
do. Not every mode is compatible with the Dirichlet walls this solver uses,
though, and the wrong one produces a spurious force that looks exactly like a
solver bug. A displacement xi = A sin(pi x/L) induces a potential proportional
to cos(pi x/L), which takes the values +C and -C at the two walls; forcing
phi(0)=phi(L)=0 then needs a LINEAR homogeneous correction whose gradient is a
uniform spurious field of -2C/L across the whole domain. It scales with the
mode amplitude, so it drives a growing bulk drift that swamps the oscillation
(measured: -210 V/m against a restoring field of ~330 V/m). With
xi = A sin(2 pi x/L) the induced potential takes the SAME value at both walls,
the correction is a pure constant, and the field reduces exactly to the
analytic cold-slab result E = e n0 xi / eps0 -- verified to 0.05%, residual
uniform offset 2e-15 V/m.

WHY THE MEASUREMENT WINDOW IS LIMITED. A perfectly cold plasma has lambda_D = 0,
so the accuracy condition lambda_D/dx >~ 2 (Turner eq. 2) is maximally violated
and the finite-grid instability is unavoidable. It is present from step zero and
grows exponentially: the mode envelope is flat to ~2% for the first ~1200 steps
and then runs away, reaching an order of magnitude by step 4000. That is the algorithm behaving as
documented, not a bug -- but it means the frequency must be measured inside the
flat window, and the script asserts that the envelope really was flat there.
The real discharge has warm electrons and lambda_D/dx = 3.72, which is exactly
why the benchmark prescribes that resolution.

WHAT THE ANSWER SHOULD BE. Not omega_pe exactly. Leapfrog integration of a
harmonic oscillator has a known dispersion error,
    sin(omega_num dt / 2) = omega_pe dt / 2
so omega_num > omega_pe by +0.042% at omega_pe*dt = 0.1. Checking against THAT
is the sharper test: it says the only error left is the integrator's documented
dispersion.

Run:  python -m scripts.verify_core
"""
from __future__ import annotations

import time

import numpy as np

from src.pic_mcc.constants import EPS0, E_CHARGE, M_E
from src.pic_mcc.grid import Grid1D
from src.pic_mcc.poisson import PoissonSolver1D, electric_field
from src.pic_mcc.pusher import half_step_back, push

L = 6.7e-2
N_CELLS = 128
N_DENSITY = 1.363e14          # m^-3, Turner Case 1 mid-plane electron density
PARTICLES_PER_CELL = 512
AMPLITUDE = 0.002             # fraction of L
OMEGA_DT = 0.1
N_STEPS = 1200                # inside the flat window; see module docstring
N_STEPS_INSTABILITY = 4000    # continued, to measure the growth rate


def plasma_frequency(n: float) -> float:
    return np.sqrt(n * E_CHARGE ** 2 / (EPS0 * M_E))


def _omega_from_zero_crossings(signal: np.ndarray, dt: float) -> float:
    """Successive zero crossings are half a period apart. Linear-regressing
    crossing time against crossing index over many periods is far more precise
    than an FFT bin: it uses the phase information, not just the magnitude."""
    s = signal - signal.mean()
    idx = np.where(np.sign(s[:-1]) != np.sign(s[1:]))[0]
    # sub-step crossing time by linear interpolation
    frac = s[idx] / (s[idx] - s[idx + 1])
    t_cross = (idx + frac) * dt
    k = np.arange(t_cross.size)
    slope = np.polyfit(k, t_cross, 1)[0]     # = T/2
    return np.pi / slope


def run_oscillation():
    grid = Grid1D(L, N_CELLS)
    solver = PoissonSolver1D(N_CELLS, grid.dx)

    omega_pe = plasma_frequency(N_DENSITY)
    dt = OMEGA_DT / omega_pe

    n_particles = N_CELLS * PARTICLES_PER_CELL
    weight = N_DENSITY * L / n_particles

    x0 = (np.arange(n_particles) + 0.5) * (L / n_particles)
    x = x0 + AMPLITUDE * L * np.sin(2 * np.pi * x0 / L)
    vx = np.zeros(n_particles)

    rho_background = E_CHARGE * N_DENSITY
    qm = -E_CHARGE / M_E
    shape = np.sin(2 * np.pi * x0 / L)
    norm = np.sum(shape ** 2)

    def field_at_particles(xp):
        rho = rho_background - E_CHARGE * grid.deposit(xp, weight) / grid.node_width
        phi = solver.solve(rho, 0.0, 0.0)
        return grid.gather(electric_field(phi, grid.dx), xp)

    vx = half_step_back(vx, field_at_particles(x), qm, dt)

    amplitude = np.empty(N_STEPS_INSTABILITY)
    for step in range(N_STEPS_INSTABILITY):
        amplitude[step] = np.sum((x - x0) * shape) / norm
        push(x, vx, field_at_particles(x), qm, dt)
    amplitude /= AMPLITUDE * L

    window = amplitude[:N_STEPS]
    envelope_growth = np.max(np.abs(window)) / np.max(np.abs(amplitude[:200]))
    if envelope_growth > 1.05:
        raise RuntimeError(
            f"envelope grew {envelope_growth:.2f}x inside the measurement window; "
            "the finite-grid instability has reached it. Shorten N_STEPS.")

    omega_measured = _omega_from_zero_crossings(window, dt)
    omega_leapfrog = (2.0 / dt) * np.arcsin(omega_pe * dt / 2.0)
    late_growth = np.max(np.abs(amplitude[3500:])) / np.max(np.abs(amplitude[:200]))

    return dict(omega_pe=omega_pe, omega_measured=omega_measured,
                omega_leapfrog=omega_leapfrog, dt=dt,
                envelope_growth=envelope_growth, late_growth=late_growth,
                n_periods=window.size * dt * omega_pe / (2 * np.pi),
                lambda_D_over_dx=np.sqrt(EPS0 * 9.36 / (E_CHARGE * N_DENSITY)) / grid.dx)


def time_one_step():
    """Cost per step of the electrostatic core at Turner Case 1 particle counts,
    so we know before committing whether numba is actually needed."""
    grid = Grid1D(L, N_CELLS)
    solver = PoissonSolver1D(N_CELLS, grid.dx)
    n = N_CELLS * PARTICLES_PER_CELL
    rng = np.random.default_rng(0)
    xe, xi = rng.uniform(0, L, n), rng.uniform(0, L, n)
    ve, vi = rng.normal(0, 1e6, n), rng.normal(0, 1e3, n)
    dt = 1.0 / (400 * 13.56e6)

    def one_step():
        rho = grid.charge_density([xi, xe], [E_CHARGE, -E_CHARGE], 1e10)
        E = electric_field(solver.solve(rho, 450.0, 0.0), grid.dx)
        push(xe, ve, grid.gather(E, xe), -E_CHARGE / M_E, dt)
        push(xi, vi, grid.gather(E, xi), E_CHARGE / 6.67e-27, dt)
        np.clip(xe, 0.0, L, out=xe)
        np.clip(xi, 0.0, L, out=xi)

    for _ in range(20):
        one_step()
    t0 = time.perf_counter()
    for _ in range(200):
        one_step()
    return (time.perf_counter() - t0) / 200


if __name__ == "__main__":
    r = run_oscillation()
    print("Cold plasma oscillation test")
    print(f"  lambda_D/dx at Case 1 conditions = {r['lambda_D_over_dx']:.3f}   (Turner Table III: 3.72)")
    print(f"  omega_pe * dt                    = {r['omega_pe'] * r['dt']:.4f}")
    print(f"  periods measured over             = {r['n_periods']:.1f}")
    print(f"  envelope growth in window        = {r['envelope_growth']:.3f}x  (must stay < 1.05)")
    print(f"  omega_pe       (analytic)        = {r['omega_pe']:.6e} rad/s")
    print(f"  omega_leapfrog (expected)        = {r['omega_leapfrog']:.6e} rad/s")
    print(f"  omega_measured (simulated)       = {r['omega_measured']:.6e} rad/s")
    print(f"  error vs omega_pe                = {100 * (r['omega_measured'] / r['omega_pe'] - 1):+.4f} %")
    print(f"  error vs leapfrog prediction     = {100 * (r['omega_measured'] / r['omega_leapfrog'] - 1):+.4f} %")
    print(f"  finite-grid instability by step {N_STEPS_INSTABILITY} = {r['late_growth']:.1f}x  (expected; cold plasma)")

    t = time_one_step()
    print("\nPerformance (pure NumPy, 2 species x 65,536 particles)")
    print(f"  cost per step                    = {t * 1e3:.3f} ms")
    print(f"  Turner Case 1 (512,000 steps)    = {t * 512_000 / 60:.0f} min  (fields only, no MCC yet)")
