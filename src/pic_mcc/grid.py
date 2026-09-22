"""Grid, cloud-in-cell deposition, and field gather.

CRITICAL INVARIANT: deposition (scatter) and field interpolation (gather) use
the SAME linear weights. If they differ, a particle exerts a force on itself
and the scheme stops conserving momentum. Every PIC bug list starts here.

Node volumes: an interior node "owns" dx of length; the two boundary nodes own
dx/2. Charge density at the boundary nodes is therefore divided by dx/2, not
dx. Getting this wrong biases the sheath fields, which is exactly where this
simulation has to be right.
"""
from __future__ import annotations

import numpy as np


class Grid1D:
    def __init__(self, length: float, n_cells: int):
        self.length = length
        self.n_cells = n_cells
        self.dx = length / n_cells
        self.n_nodes = n_cells + 1
        self.x = np.linspace(0.0, length, self.n_nodes)
        # Node "volume" per unit area, for converting deposited charge -> density
        self.node_width = np.full(self.n_nodes, self.dx)
        self.node_width[0] = self.dx / 2.0
        self.node_width[-1] = self.dx / 2.0

    # -- scatter -----------------------------------------------------------
    def deposit(self, x: np.ndarray, weight: float | np.ndarray) -> np.ndarray:
        """Accumulate particle weight onto nodes with linear (CIC) weights.

        Returns the per-node accumulated weight (NOT a density).
        """
        s = x / self.dx
        j = np.floor(s).astype(np.int64)
        np.clip(j, 0, self.n_cells - 1, out=j)
        w = s - j
        if np.isscalar(weight):
            lo = (1.0 - w) * weight
            hi = w * weight
        else:
            lo = (1.0 - w) * weight
            hi = w * weight
        acc = np.bincount(j, weights=lo, minlength=self.n_nodes)
        acc += np.bincount(j + 1, weights=hi, minlength=self.n_nodes)
        return acc

    def charge_density(self, x_list, q_list, weight: float) -> np.ndarray:
        """Net charge density (C/m^3) at nodes from a list of species."""
        acc = np.zeros(self.n_nodes)
        for x, q in zip(x_list, q_list):
            if x.size:
                acc += q * self.deposit(x, weight)
        return acc / self.node_width

    # -- gather ------------------------------------------------------------
    def gather(self, field: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Interpolate a node field to particle positions with the same weights."""
        s = x / self.dx
        j = np.floor(s).astype(np.int64)
        np.clip(j, 0, self.n_cells - 1, out=j)
        w = s - j
        return (1.0 - w) * field[j] + w * field[j + 1]
