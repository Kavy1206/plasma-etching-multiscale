import numpy as np

from src.feature_scale.viz import (center_depth, etch_iter, make_synthetic_ions,
                                   open_width_profile, render_rgb)


def test_etch_iter_yields_expected_frames_and_keeps_invariants():
    E, th = make_synthetic_ions(25.0, n=500)
    frames = list(etch_iter(E, th, "synthetic", width=20, sticking=0.4,
                            n_steps=60, every=20, nx=80, nz=60))
    assert [s for s, _ in frames] == [0, 20, 40, 60]
    g = frames[-1][1]
    assert g.material.min() >= 0.0 and g.material.max() <= 1.0
    assert np.all(g.material[g.is_mask] == 1.0)       # mask never erodes
    assert center_depth(g) >= 15                        # at least the open mask window
    assert open_width_profile(g).shape == (60,)


def test_render_rgb_shape_and_dtype():
    E, th = make_synthetic_ions(25.0, n=100)
    _, g = next(iter(etch_iter(E, th, "real", width=20, sticking=0.3,
                               n_steps=1, every=1, nx=80, nz=40)))
    img = render_rgb(g, scale=3, pad=10)
    assert img.dtype == np.uint8 and img.shape == (40 * 3, (20 + 20) * 3, 3)
