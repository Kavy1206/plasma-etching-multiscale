"""Smoke tests for the Streamlit app using streamlit's AppTest harness:
the page renders, the pressure slider works, and the live-etch button runs
the real model for both ion sources without exceptions."""
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app" / "streamlit_app.py")


def _fresh():
    return AppTest.from_file(APP, default_timeout=120).run()


def test_app_renders_without_exceptions():
    at = _fresh()
    assert not at.exception
    assert len(at.tabs) == 4


def test_pressure_slider_switches_histogram_and_warns_on_thin_data():
    at = _fresh()
    at.select_slider(key="pressure").set_value(100).run()
    assert not at.exception
    assert any("ion samples" in w.value for w in at.warning)   # 100 mTorr has only 11


def test_live_etch_runs_for_both_sources():
    at = _fresh()
    at.slider(key="nsteps").set_value(200)
    at.button(key="run").click().run()
    assert not at.exception
    assert len(at.session_state["frames"]) >= 2
    depths_real = [d for _, d in at.session_state["depths"]]
    assert depths_real == sorted(depths_real)                   # depth never decreases

    at.radio(key="src").set_value("Synthetic broad-angle (mechanism demo - NOT simulated argon)").run()
    at.slider(key="nsteps").set_value(200)
    at.button(key="run").click().run()
    assert not at.exception
