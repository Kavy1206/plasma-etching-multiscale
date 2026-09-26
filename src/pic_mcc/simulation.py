"""Full 1D3V PIC-MCC time loop, wiring together steps 1-6.

Order per step (matches the standard leapfrog PIC cycle, and cross-checked
against the reference implementation's ordering discussed in PROGRESS.md):

  1. deposit electron + ion density to the grid
  2. Poisson solve with V(t) = V0 sin(2 pi f t) on the driven electrode
  3. E = -grad(phi)
  4. gather E to particle positions (at x^n)
  5. leapfrog push: v^{n-1/2} -> v^{n+1/2}, x^n -> x^{n+1}
  6. absorbing walls (+ secondary emission for argon runs)
  7. null-collision MCC for both species

Checkpointing: state is saved to an .npz file periodically so a run can span
multiple process invocations (this benchmark case is ~512,000 steps and does
not reliably finish inside one tool call). Resuming re-loads particle arrays,
RNG state, and the running average accumulators exactly, so the averaged
result is identical to an uninterrupted run.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .boundaries import apply_absorbing_walls, emit_secondaries
from .config import Config
from .constants import E_CHARGE as E
from .constants import K_B
from .cross_sections import load_gas_set, max_collision_frequency
from .grid import Grid1D
from .mcc import ElectronCollisions, IonCollisions, Species, kinetic_energy_ev
from .poisson import PoissonSolver1D, electric_field
from .pusher import half_step_back, push


@dataclass
class Averages:
    """Running sums for time-averaged density and temperature profiles,
    accumulated only during the final n_cycles_average window."""
    n_nodes: int
    sum_ne: np.ndarray = None
    sum_ni: np.ndarray = None
    sum_ne_v2: np.ndarray = None  # for kTe: m<v^2>/3 via grid-binned second moment
    count: int = 0
    iedf_energy: list = field(default_factory=lambda: {0: [], 1: []})  # per wall
    iedf_angle: list = field(default_factory=lambda: {0: [], 1: []})   # deg from normal

    def __post_init__(self):
        if self.sum_ne is None:
            self.sum_ne = np.zeros(self.n_nodes)
            self.sum_ni = np.zeros(self.n_nodes)
            self.sum_ne_v2 = np.zeros(self.n_nodes)


class Simulation:
    def __init__(self, cfg: Config, checkpoint_path: str, collect_iedf: bool = False):
        self.collect_iedf = collect_iedf
        self.cfg = cfg
        self.ckpt_path = Path(checkpoint_path)
        self.grid = Grid1D(cfg.length, cfg.n_cells)
        self.poisson = PoissonSolver1D(cfg.n_cells, self.grid.dx)
        self.electron_procs, self.ion_procs = load_gas_set(cfg.gas, cfg.cross_section_dir)
        self.nu_max_e = max_collision_frequency(self.electron_procs, cfg.m_electron,
                                                 cfg.n_gas, energy_max_ev=1000.0)
        self.nu_max_i = max_collision_frequency(self.ion_procs, cfg.m_ion,
                                                 cfg.n_gas, energy_max_ev=2000.0)
        self.rng = np.random.default_rng(cfg.seed)
        self.e_coll = ElectronCollisions(self.electron_procs, cfg.m_electron, cfg.m_ion,
                                         cfg.n_gas, self.nu_max_e, self.rng)
        self.i_coll = IonCollisions(self.ion_procs, cfg.m_ion, cfg.n_gas, cfg.T_gas,
                                    self.nu_max_i, self.rng)
        self.e_coll.set_dt(cfg.dt)
        self.i_coll.set_dt(cfg.dt)

        self.step = 0
        self.avg = Averages(self.grid.n_nodes)
        self.avg_start_step = cfg.n_steps - cfg.n_steps_average
        self.diag = {"n_e": [], "n_i": [], "step": []}

        if self.ckpt_path.exists():
            self._load_checkpoint()
        else:
            self._init_particles()

    def _init_particles(self):
        cfg = self.cfg
        n0 = cfg.n_cells * cfg.particles_per_cell
        v_th_e = np.sqrt(K_B * cfg.T_e_init / cfg.m_electron)
        v_th_i = np.sqrt(K_B * cfg.T_i_init / cfg.m_ion)
        x0 = self.rng.uniform(0, cfg.length, n0)
        self.electrons = Species(x=x0.copy(),
                                 v=self.rng.normal(0, v_th_e, size=(n0, 3)),
                                 mass=cfg.m_electron, charge=-E)
        self.ions = Species(x=x0.copy(),
                            v=self.rng.normal(0, v_th_i, size=(n0, 3)),
                            mass=cfg.m_ion, charge=E)
        self.weight = cfg.weight
        # half-step back the velocities to start the leapfrog
        E_at = self._field_at(self.electrons.x, self.ions.x, boundary_v=0.0)
        rho, phi, Ex = E_at
        Ee = self.grid.gather(Ex, self.electrons.x)
        Ei = self.grid.gather(Ex, self.ions.x)
        self.electrons.v[:, 2] = half_step_back(self.electrons.v[:, 2], Ee, -E / cfg.m_electron, cfg.dt)
        self.ions.v[:, 2] = half_step_back(self.ions.v[:, 2], Ei, E / cfg.m_ion, cfg.dt)

    def _field_at(self, xe, xi, boundary_v: float):
        rho = self.grid.charge_density([xi, xe], [E, -E], self.weight)
        phi = self.poisson.solve(rho, 0.0, boundary_v)
        Ex = electric_field(phi, self.grid.dx)
        return rho, phi, Ex

    def run(self, max_wall_seconds: float, progress_every: int = 20000,
           checkpoint_every: int = 10000, stop_at_step: int | None = None) -> dict:
        """Advance the simulation until cfg.n_steps is reached, the wall-clock
        budget runs out, or (if given) stop_at_step is reached -- whichever
        comes first. stop_at_step never changes cfg.n_steps or avg_start_step;
        it exists for deterministic chunking/testing independent of wall time.
        Checkpoints either way. Returns a status dict."""
        cfg = self.cfg
        t_start = time.time()
        stopped_reason = "completed"
        target = cfg.n_steps if stop_at_step is None else min(stop_at_step, cfg.n_steps)

        while self.step < target:
            boundary_v = cfg.voltage * np.sin(2 * np.pi * cfg.frequency * cfg.dt * self.step)
            rho, phi, Ex = self._field_at(self.electrons.x, self.ions.x, boundary_v)

            Ee = self.grid.gather(Ex, self.electrons.x)
            Ei = self.grid.gather(Ex, self.ions.x)
            push(self.electrons.x, self.electrons.v[:, 2], Ee, -E / cfg.m_electron, cfg.dt)
            push(self.ions.x, self.ions.v[:, 2], Ei, E / cfg.m_ion, cfg.dt)

            wall_e, xe_abs, ve_abs = apply_absorbing_walls(self.electrons, cfg.length)
            wall_i, xi_abs, vi_abs = apply_absorbing_walls(self.ions, cfg.length)
            if self.step >= self.avg_start_step and wall_i.size and self.collect_iedf:
                E_ion = kinetic_energy_ev(vi_abs, cfg.m_ion)
                # push() advances v[:,2] (see push() calls above) -- that is
                # the wall-normal (1D electrostatic) direction, not v[:,0].
                # An earlier version of this line used v[:,0] as "normal",
                # swapping normal and perpendicular and producing a nonsense
                # ~88 deg mean IADF angle for what should be near-normal
                # sheath-accelerated incidence. Caught while consuming this
                # data in the Phase 3 feature-scale model; see PROGRESS.md.
                v_normal = np.abs(vi_abs[:, 2])
                v_perp = np.sqrt(vi_abs[:, 0] ** 2 + vi_abs[:, 1] ** 2)
                angle_deg = np.degrees(np.arctan2(v_perp, np.maximum(v_normal, 1e-30)))
                for w in (0, 1):
                    sel = wall_i == w
                    if np.any(sel):
                        self.avg.iedf_energy[w].append(E_ion[sel])
                        self.avg.iedf_angle[w].append(angle_deg[sel])
            if cfg.gamma_see > 0.0 and wall_i.size:
                x_new, v_new = emit_secondaries(wall_i, cfg.gamma_see, cfg.T_see,
                                                cfg.m_electron, cfg.length, self.rng)
                if x_new.size:
                    self.electrons.append(x_new, v_new)

            self.e_coll.apply(self.electrons, self.ions, cfg.T_gas)
            self.i_coll.apply(self.ions)

            if self.step >= self.avg_start_step:
                acc_e = self.grid.deposit(self.electrons.x, self.weight) / self.grid.node_width
                acc_i = self.grid.deposit(self.ions.x, self.weight) / self.grid.node_width
                self.avg.sum_ne += acc_e
                self.avg.sum_ni += acc_i
                v2 = self.electrons.v[:, 0] ** 2 + self.electrons.v[:, 1] ** 2 + self.electrons.v[:, 2] ** 2
                acc_ev2 = self.grid.deposit(self.electrons.x, self.weight) / self.grid.node_width  # count only
                # Energy-weighted deposit for kTe: sum(w*v^2) per node / sum(w) per node
                acc_ev2_num = np.zeros(self.grid.n_nodes)
                s = self.electrons.x / self.grid.dx
                j = np.clip(np.floor(s).astype(np.int64), 0, cfg.n_cells - 1)
                wgt = s - j
                np.add.at(acc_ev2_num, j, (1 - wgt) * v2 * self.weight)
                np.add.at(acc_ev2_num, j + 1, wgt * v2 * self.weight)
                self.avg.sum_ne_v2 += acc_ev2_num / self.grid.node_width
                self.avg.count += 1

            self.step += 1
            if self.step % progress_every == 0:
                n_e, n_i = self.electrons.n, self.ions.n
                self.diag["step"].append(self.step)
                self.diag["n_e"].append(n_e)
                self.diag["n_i"].append(n_i)

            if self.step % checkpoint_every == 0:
                self._save_checkpoint()

            if time.time() - t_start > max_wall_seconds:
                stopped_reason = "wall_time_budget"
                break

        self._save_checkpoint()
        elapsed = time.time() - t_start
        return dict(step=self.step, n_steps=cfg.n_steps, reason=stopped_reason,
                   elapsed_s=elapsed, n_electrons=self.electrons.n, n_ions=self.ions.n,
                   steps_this_call=None)

    def _save_checkpoint(self):
        iedf_e0 = np.concatenate(self.avg.iedf_energy[0]) if self.avg.iedf_energy[0] else np.zeros(0)
        iedf_e1 = np.concatenate(self.avg.iedf_energy[1]) if self.avg.iedf_energy[1] else np.zeros(0)
        iedf_a0 = np.concatenate(self.avg.iedf_angle[0]) if self.avg.iedf_angle[0] else np.zeros(0)
        iedf_a1 = np.concatenate(self.avg.iedf_angle[1]) if self.avg.iedf_angle[1] else np.zeros(0)
        np.savez(self.ckpt_path,
                xe=self.electrons.x, ve=self.electrons.v,
                xi=self.ions.x, vi=self.ions.v,
                step=self.step, weight=self.weight,
                sum_ne=self.avg.sum_ne, sum_ni=self.avg.sum_ni,
                sum_ne_v2=self.avg.sum_ne_v2, count=self.avg.count,
                iedf_e0=iedf_e0, iedf_e1=iedf_e1, iedf_a0=iedf_a0, iedf_a1=iedf_a1,
                rng_state=np.void(bytes(str(self.rng.bit_generator.state), "utf-8")),
                diag_step=np.array(self.diag["step"]), diag_ne=np.array(self.diag["n_e"]),
                diag_ni=np.array(self.diag["n_i"]))
        import pickle
        with open(str(self.ckpt_path) + ".rng.pkl", "wb") as f:
            pickle.dump(self.rng.bit_generator.state, f)

    def _load_checkpoint(self):
        d = np.load(self.ckpt_path, allow_pickle=True)
        cfg = self.cfg
        self.electrons = Species(x=d["xe"], v=d["ve"], mass=cfg.m_electron, charge=-E)
        self.ions = Species(x=d["xi"], v=d["vi"], mass=cfg.m_ion, charge=E)
        self.step = int(d["step"])
        self.weight = float(d["weight"])
        self.avg.sum_ne = d["sum_ne"]
        self.avg.sum_ni = d["sum_ni"]
        self.avg.sum_ne_v2 = d["sum_ne_v2"]
        self.avg.count = int(d["count"])
        # Restore as single-element lists so future appends accumulate correctly
        # and iedf_iadf()'s concatenate keeps working either way.
        for w, key_e, key_a in ((0, "iedf_e0", "iedf_a0"), (1, "iedf_e1", "iedf_a1")):
            if key_e in d and d[key_e].size:
                self.avg.iedf_energy[w] = [d[key_e]]
                self.avg.iedf_angle[w] = [d[key_a]]
        self.diag = {"step": list(d["diag_step"]), "n_e": list(d["diag_ne"]), "n_i": list(d["diag_ni"])}
        import pickle
        with open(str(self.ckpt_path) + ".rng.pkl", "rb") as f:
            state = pickle.load(f)
        self.rng.bit_generator.state = state

    def averaged_profiles(self):
        if self.avg.count == 0:
            raise RuntimeError("no averaging samples accumulated yet")
        n_e = self.avg.sum_ne / self.avg.count
        n_i = self.avg.sum_ni / self.avg.count
        mean_v2 = self.avg.sum_ne_v2 / np.maximum(self.avg.sum_ne, 1e-30)
        kTe_ev = (self.cfg.m_electron * mean_v2 / 3.0) / E
        return self.grid.x, n_e, n_i, kTe_ev

    def iedf_iadf(self, wall: int):
        e = self.avg.iedf_energy[wall]
        a = self.avg.iedf_angle[wall]
        if not e:
            return np.zeros(0), np.zeros(0)
        return np.concatenate(e), np.concatenate(a)
