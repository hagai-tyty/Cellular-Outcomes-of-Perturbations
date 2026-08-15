"""Unit tests for the paired target audit.

The pure functions are tested on constructed arrays where the right answer is known by
arithmetic, because the whole point of this diagnostic is to decide between three readings and a
branch that is never exercised is not a check (the `verify_1a` lesson).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("dts", ROOT / "experiments" / "diag_target_shift.py")
dts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dts)


# ---- ols ---------------------------------------------------------------------------------- #
def test_ols_recovers_an_exact_line():
    x = np.arange(50, dtype=float)
    slope, intercept, r2 = dts.ols(x, 2.5 * x - 7.0)
    assert slope == pytest.approx(2.5, abs=1e-9)
    assert intercept == pytest.approx(-7.0, abs=1e-8)
    assert r2 == pytest.approx(1.0, abs=1e-12)


def test_ols_on_a_pure_offset_gives_slope_one():
    x = np.linspace(-20, 20, 100)
    slope, intercept, r2 = dts.ols(x, x + 5.0)
    assert slope == pytest.approx(1.0, abs=1e-9)
    assert intercept == pytest.approx(5.0, abs=1e-8)


def test_ols_r2_falls_when_the_change_is_not_linear():
    rng = np.random.default_rng(0)
    x = np.linspace(1, 10, 200)
    _, _, r2 = dts.ols(x, x ** 3 + rng.normal(0, 20, 200))
    assert r2 < 1.0


def test_ols_is_undefined_on_a_constant_x():
    """np.polyfit would emit a rank warning and return garbage rather than refuse."""
    s, i, r = dts.ols(np.full(10, 3.0), np.arange(10, dtype=float))
    assert np.isnan(s) and np.isnan(i) and np.isnan(r)


def test_ols_is_undefined_on_too_few_points():
    s, _, _ = dts.ols(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    assert np.isnan(s)


# ---- verdict ------------------------------------------------------------------------------ #
def test_pure_offset_reads_as_offset():
    assert dts.verdict(slope=1.0, r2=1.0, time_spread=0.0) == "A: OFFSET"


def test_slope_just_inside_the_tolerance_is_offset():
    """Not tested AT the boundary: 1.0 + 0.05 - 1.0 == 0.050000000000000044 in binary floating
    point, so an 'exactly at the edge' assertion tests the representation, not the rule."""
    assert dts.verdict(1.0 + dts.SLOPE_TOL * 0.99, 1.0, 0.0) == "A: OFFSET"
    assert dts.verdict(1.0 - dts.SLOPE_TOL * 0.99, 1.0, 0.0) == "A: OFFSET"


def test_slope_past_the_tolerance_is_scale():
    assert dts.verdict(1.0 + dts.SLOPE_TOL * 1.01, 1.0, 0.0) == "B: SCALE"
    assert dts.verdict(1.0 - dts.SLOPE_TOL * 1.01, 1.0, 0.0) == "B: SCALE"


def test_a_scale_change_reads_as_scale():
    assert dts.verdict(slope=1.8, r2=0.99, time_spread=0.5) == "B: SCALE"


def test_a_poor_linear_fit_reads_as_nonlinear():
    assert dts.verdict(slope=1.0, r2=0.10, time_spread=0.0) == "C: NONLINEAR"


def test_a_time_varying_shift_reads_as_nonlinear_even_when_linear_overall():
    """The load-bearing ordering. A change can be perfectly linear pooled and still be a
    DIFFERENT change at each timepoint -- only the latter forbids a simple re-centring."""
    assert dts.verdict(slope=1.0, r2=1.0, time_spread=dts.TIME_TOL + 1e-9).startswith("C:")


def test_time_spread_at_the_tolerance_edge_is_not_yet_nonlinear():
    assert dts.verdict(1.0, 1.0, dts.TIME_TOL) == "A: OFFSET"


def test_nan_inputs_are_undetermined_not_silently_a_verdict():
    assert dts.verdict(float("nan"), 1.0, 0.0) == "UNDETERMINED"
    assert dts.verdict(1.0, float("nan"), 0.0) == "UNDETERMINED"


# ---- stratum_stats ------------------------------------------------------------------------ #
def test_stratum_stats_on_a_known_offset():
    old = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    s = dts.stratum_stats(old, old + 3.0)
    assert s["n"] == 5
    assert s["mean_shift"] == pytest.approx(3.0)
    assert s["sd_shift"] == pytest.approx(0.0, abs=1e-12)
    assert s["slope"] == pytest.approx(1.0, abs=1e-9)
    assert s["corr"] == pytest.approx(1.0, abs=1e-12)
    assert s["mean_c7"] - s["mean_old"] == pytest.approx(3.0)


def test_stratum_stats_on_a_known_scale_change():
    old = np.linspace(0, 10, 40)
    s = dts.stratum_stats(old, 0.5 * old)
    assert s["slope"] == pytest.approx(0.5, abs=1e-9)
    assert s["sd_c7"] == pytest.approx(0.5 * s["sd_old"], rel=1e-9)
    assert s["sd_shift"] > 0            # a scale change is NOT a constant shift


def test_a_scale_change_is_distinguishable_from_an_offset_by_sd_shift():
    """The two readings must not be confusable on the numbers this reports."""
    old = np.linspace(0, 10, 40)
    off = dts.stratum_stats(old, old + 2.0)
    sca = dts.stratum_stats(old, 1.5 * old)
    assert off["sd_shift"] == pytest.approx(0.0, abs=1e-12)
    assert sca["sd_shift"] > 1.0
    assert dts.verdict(off["slope"], off["r2"], 0.0) == "A: OFFSET"
    assert dts.verdict(sca["slope"], sca["r2"], 0.0) == "B: SCALE"


def test_thresholds_are_stated_not_derived_from_the_data():
    """Pre-registration: these are module constants, so they cannot drift to fit a result."""
    assert dts.SLOPE_TOL == 0.05 and dts.R2_FLOOR == 0.90 and dts.TIME_TOL == 2.0
