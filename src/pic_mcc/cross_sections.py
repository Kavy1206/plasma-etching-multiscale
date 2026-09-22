"""Cross-section tables and the null-collision setup.

Rules imposed by Turner et al. (2013) sec. II, all of which are enforced here:

  * LINEAR interpolation between tabulated points (not spline, not log-log).
  * Above the last tabulated energy, CLAMP to the last value.
  * Electron tables are functions of LAB energy.
  * Ion tables are functions of CENTRE-OF-MASS energy. For an ion hitting a
    same-mass neutral at rest, E_cm = E_lab / 2. Forgetting this halves or
    doubles the ion collision rate and is the single most common bug in an
    otherwise correct MCC implementation.

File format: "energy_eV;sigma_m2", semicolon-delimited, no header.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class CrossSection:
    name: str
    energy_ev: np.ndarray      # ascending
    sigma_m2: np.ndarray
    kind: str                  # "elastic" | "excitation" | "ionization" | "backscatter"
    threshold_ev: float = 0.0
    frame: str = "lab"         # "lab" for electrons, "cm" for ions

    def __post_init__(self) -> None:
        if self.energy_ev.ndim != 1 or self.energy_ev.shape != self.sigma_m2.shape:
            raise ValueError(f"{self.name}: mismatched table shapes")
        if np.any(np.diff(self.energy_ev) < 0):
            raise ValueError(f"{self.name}: energy column is not monotonically increasing")
        if np.any(self.sigma_m2 < 0):
            raise ValueError(f"{self.name}: negative cross section")

    def __call__(self, energy_ev: np.ndarray) -> np.ndarray:
        """Linear interpolation, clamped at both ends (benchmark requirement)."""
        return np.interp(energy_ev, self.energy_ev, self.sigma_m2,
                         left=self.sigma_m2[0], right=self.sigma_m2[-1])


def load_table(path: Path, name: str, kind: str, frame: str = "lab",
               threshold_ev: float = 0.0) -> CrossSection:
    raw = np.loadtxt(path, delimiter=";")
    return CrossSection(name=name, energy_ev=raw[:, 0], sigma_m2=raw[:, 1],
                        kind=kind, frame=frame, threshold_ev=threshold_ev)


HELIUM_FILES = [
    ("Elastic_He.csv",         "e-He elastic",     "elastic",     "lab", 0.0),
    ("Excitation1_He.csv",     "e-He -> He*",      "excitation",  "lab", 19.82),
    ("Excitation2_He.csv",     "e-He -> He**",     "excitation",  "lab", 20.61),
    ("Ionization_He.csv",      "e-He ionization",  "ionization",  "lab", 24.587),
]
HELIUM_ION_FILES = [
    ("Isotropic_He.csv",       "He+ -He isotropic",    "elastic",      "cm", 0.0),
    ("Backscattering_He.csv",  "He+ -He backscatter",  "backscatter",  "cm", 0.0),
]


def load_helium_set(directory: str | Path):
    """Returns (electron_processes, ion_processes) for the Turner benchmark."""
    d = Path(directory)
    electrons = [load_table(d / f, n, k, fr, th) for f, n, k, fr, th in HELIUM_FILES]
    ions = [load_table(d / f, n, k, fr, th) for f, n, k, fr, th in HELIUM_ION_FILES]
    return electrons, ions


def max_collision_frequency(processes, mass: float, n_gas: float,
                            energy_max_ev: float = 1.0e4, n_samples: int = 20000) -> float:
    """nu_max = n_gas * max_over_E[ sigma_total(E) * v(E) ], for the null-collision method.

    Scanning on a dense grid is deliberate: sigma*v is not monotonic (it peaks
    somewhere above the ionization threshold), so taking the endpoint value
    would under-estimate nu_max and silently lose collisions.
    """
    from .constants import EV
    energy = np.linspace(0.0, energy_max_ev, n_samples)
    sigma_total = np.zeros_like(energy)
    for p in processes:
        sigma_total += p(energy)
    speed = np.sqrt(2.0 * energy * EV / mass)
    return float(n_gas * np.max(sigma_total * speed))


def lab_to_cm_energy(energy_lab_ev: np.ndarray, m_projectile: float,
                     m_target: float) -> np.ndarray:
    """E_cm = E_lab * m_target / (m_projectile + m_target), target at rest."""
    return energy_lab_ev * m_target / (m_projectile + m_target)
