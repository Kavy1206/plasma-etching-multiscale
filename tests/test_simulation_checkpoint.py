"""Verifies checkpoint/resume reproduces an uninterrupted run bit-for-bit.

This matters because the full Turner Case 1 run (512,000 steps) is chunked
across multiple process invocations for wall-clock reasons -- if resume
weren't exact, the averaged validation result would depend on how the run
happened to be chopped up, which would make the reported agreement number
meaningless.
"""
import os

import numpy as np

from src.pic_mcc.config import turner_case
from src.pic_mcc.simulation import Simulation


def _make_cfg():
    cfg = turner_case(1)
    cfg.n_cells = 24
    cfg.particles_per_cell = 32
    cfg.n_cycles = 3
    cfg.n_cycles_average = 1
    cfg.__post_init__()
    return cfg


def test_checkpoint_resume_matches_uninterrupted_run_bitwise(tmp_path):
    sim1 = Simulation(_make_cfg(), str(tmp_path / "a.npz"))
    sim1.run(max_wall_seconds=120)
    x1, ne1, ni1, kte1 = sim1.averaged_profiles()

    ckpt = tmp_path / "b.npz"
    sim2 = Simulation(_make_cfg(), str(ckpt))
    stop_mid_average = sim1.avg_start_step + sim1.avg.count // 3
    sim2.run(max_wall_seconds=120, stop_at_step=stop_mid_average)
    assert sim2.step == stop_mid_average

    sim3 = Simulation(_make_cfg(), str(ckpt))  # reload from checkpoint
    sim3.run(max_wall_seconds=120)
    x2, ne2, ni2, kte2 = sim3.averaged_profiles()

    assert sim1.step == sim3.step
    assert sim1.avg.count == sim3.avg.count
    assert np.array_equal(ne1, ne2)
    assert np.array_equal(ni1, ni2)
    assert np.array_equal(kte1, kte2)
    assert np.array_equal(sim1.electrons.x, sim3.electrons.x)
    assert np.array_equal(sim1.ions.v, sim3.ions.v)
