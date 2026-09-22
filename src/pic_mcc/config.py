"""Run configuration.

One config object describes a whole run. The `gas` field is the switch that
selects between:

  * "he" -- helium, for the Turner et al. (2013) VALIDATION runs. Secondary
            emission MUST be zero here: the benchmark states that fluxes
            reaching the electrodes are completely absorbed with no secondary
            particles emitted.
  * "ar" -- argon, for the PRODUCTION pressure/voltage sweep whose IEDF/IADF
            Phase 3 consumes. Secondary emission is on here.

Never validate an argon run against the helium benchmark. See docs/theory.md
sec. 9.1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .constants import AMU, K_B, M_E

GAS_PROPERTIES = {
    # name: (ion mass kg, ionization threshold eV, default cross-section dir)
    "he": (4.002602 * AMU, 24.587, "data/cross_sections/turner_benchmark_he"),
    "ar": (39.948 * AMU, 15.760, "data/cross_sections/ar"),
}


@dataclass
class Config:
    # --- gas -------------------------------------------------------------
    gas: str = "he"
    m_ion: float | None = None          # kg; None -> from GAS_PROPERTIES
    m_electron: float = M_E
    cross_section_dir: str | None = None

    # --- geometry / operating point --------------------------------------
    length: float = 6.7e-2              # m, electrode gap
    n_gas: float = 9.64e20              # m^-3, fixed neutral density
    T_gas: float = 300.0                # K
    voltage: float = 450.0              # V, amplitude of V(t)=V0 sin(2 pi f t)
    frequency: float = 13.56e6          # Hz

    # --- surfaces ---------------------------------------------------------
    gamma_see: float = 0.0              # secondary electron yield; 0 for Turner
    T_see: float = 2.0                  # eV, emitted secondary energy

    # --- initial conditions ----------------------------------------------
    n_plasma_init: float = 2.56e14      # m^-3
    T_e_init: float = 30000.0           # K
    T_i_init: float = 300.0             # K

    # --- numerics ---------------------------------------------------------
    n_cells: int = 128                  # dx = length / n_cells
    steps_per_cycle: int = 400          # dt = 1 / (steps_per_cycle * frequency)
    particles_per_cell: int = 512
    n_cycles: int = 1280                # total RF cycles to run
    n_cycles_average: int = 32          # averaging window at the end
    seed: int = 20260921

    def __post_init__(self) -> None:
        if self.gas not in GAS_PROPERTIES:
            raise ValueError(f"unknown gas {self.gas!r}; expected one of {list(GAS_PROPERTIES)}")
        m, _, xsdir = GAS_PROPERTIES[self.gas]
        if self.m_ion is None:
            self.m_ion = m
        if self.cross_section_dir is None:
            self.cross_section_dir = xsdir
        if self.gas == "he" and self.gamma_see != 0.0:
            raise ValueError(
                "gamma_see must be 0 for helium benchmark runs: Turner et al. "
                "specify complete absorption with no secondary emission. "
                "Pass gas='ar' for production runs with secondaries."
            )

    # --- derived ----------------------------------------------------------
    @property
    def dx(self) -> float:
        return self.length / self.n_cells

    @property
    def dt(self) -> float:
        return 1.0 / (self.steps_per_cycle * self.frequency)

    @property
    def n_steps(self) -> int:
        return self.n_cycles * self.steps_per_cycle

    @property
    def n_steps_average(self) -> int:
        return self.n_cycles_average * self.steps_per_cycle

    @property
    def pressure_pa(self) -> float:
        return self.n_gas * K_B * self.T_gas

    @property
    def pressure_mtorr(self) -> float:
        return self.pressure_pa / 133.322368e-3

    @property
    def weight(self) -> float:
        """Superparticle weight (physical particles per superparticle, per m^2)."""
        return self.n_plasma_init * self.dx / self.particles_per_cell


def turner_case(case: int) -> Config:
    """Turner et al. (2013) Phys. Plasmas 20, 013507 -- Table I.

    Helium, 13.56 MHz, 6.7 cm gap, no secondary emission.
    """
    params = {
        1: dict(n_gas=9.64e20,  voltage=450.0, n_plasma_init=2.56e14,
                n_cells=128, steps_per_cycle=400,  particles_per_cell=512,  n_cycles=1280),
        2: dict(n_gas=32.1e20,  voltage=200.0, n_plasma_init=5.12e14,
                n_cells=256, steps_per_cycle=800,  particles_per_cell=256,  n_cycles=5120),
        3: dict(n_gas=96.4e20,  voltage=150.0, n_plasma_init=5.12e14,
                n_cells=512, steps_per_cycle=1600, particles_per_cell=128,  n_cycles=5120),
        4: dict(n_gas=321e20,   voltage=120.0, n_plasma_init=3.84e14,
                n_cells=512, steps_per_cycle=3200, particles_per_cell=64,   n_cycles=15360),
    }
    if case not in params:
        raise ValueError("Turner benchmark cases are 1-4")
    return Config(
        gas="he",
        m_ion=6.67e-27,        # prescribed by the benchmark, not the CODATA value
        m_electron=9.109e-31,  # likewise
        gamma_see=0.0,
        **params[case],
    )


# Reference results, Turner et al. Table III. Used by the validation script.
TURNER_REFERENCE = {
    1: dict(n_i_midplane=0.140e15, kTe_midplane=9.36, J_i=0.219,
            omega_pe_dt=0.121, lambda_D_over_dx=3.72),
    2: dict(n_i_midplane=0.828e15, kTe_midplane=4.69, J_i=0.215,
            omega_pe_dt=0.150, lambda_D_over_dx=2.14),
    3: dict(n_i_midplane=1.81e15,  kTe_midplane=3.95, J_i=0.195,
            omega_pe_dt=0.110, lambda_D_over_dx=2.66),
    4: dict(n_i_midplane=2.57e15,  kTe_midplane=3.65, J_i=0.186,
            omega_pe_dt=0.066, lambda_D_over_dx=2.14),
}
