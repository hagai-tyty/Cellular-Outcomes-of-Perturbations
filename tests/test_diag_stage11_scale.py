"""Unit tests for Stage 11 -- ΔAge scale calibration.

Three things could have made this stage report a false success: a rescale silently changing the
ordering, `k` leaking from the rows it is scored on, and least-squares shrinkage being mistaken for
accuracy. Each is pinned here.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage11_scale.py"
spec = importlib.util.spec_from_file_location("s11", SRC)
s11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s11)


# ---- the mathematical fact the whole stage rests on ---------------------------------------- #
def test_a_positive_rescale_cannot_change_spearman():
    """Stated in the plan BEFORE the run: rank order is invariant to any positive monotone
    transform, so rho after calibration is arithmetic, not a finding."""
    rng = np.random.default_rng(0)
    p, y = rng.normal(size=60), rng.normal(size=60)
    before = s11.spearman(p, y)
    for k in (0.01, 0.37, 1.0, 5.0, 100.0):
        assert s11.spearman(k * p, y) == pytest.approx(before, abs=1e-12)


def test_a_negative_rescale_flips_it_and_is_therefore_not_a_calibration():
    rng = np.random.default_rng(1)
    p, y = rng.normal(size=60), rng.normal(size=60)
    assert s11.spearman(-p, y) == pytest.approx(-s11.spearman(p, y), abs=1e-12)


# ---- least squares vs variance matching: the shrinkage trade -------------------------------- #
def test_least_squares_shrinks_below_variance_matching_when_correlation_is_imperfect():
    """k_LS = rho * SD(y)/SD(p). With rho < 1 it is strictly smaller than the variance-matching
    k -- which is why it wins on MAE and UNDER-reports magnitude."""
    rng = np.random.default_rng(2)
    truth = rng.normal(0, 10, 400)
    pred = 1.66 * truth + rng.normal(0, 8, 400)      # inflated AND noisy
    k_ls, _ = s11.fit_scale(pred, truth, with_offset=False)
    k_var, _ = s11.fit_scale_variance(pred, truth)
    assert k_ls < k_var


def test_variance_matching_restores_the_spread_exactly():
    rng = np.random.default_rng(3)
    truth = rng.normal(0, 10, 300)
    pred = 1.66 * truth + rng.normal(0, 5, 300)
    k, c = s11.fit_scale_variance(pred, truth)
    assert c == 0.0
    assert np.std(k * pred, ddof=1) == pytest.approx(np.std(truth, ddof=1), rel=1e-9)


def test_least_squares_recovers_the_exact_factor_when_there_is_no_noise():
    truth = np.linspace(-20, 20, 50)
    k, c = s11.fit_scale(2.5 * truth, truth, with_offset=False)
    assert k == pytest.approx(0.4, abs=1e-9) and c == 0.0


def test_the_offset_form_recovers_a_shift():
    truth = np.linspace(-20, 20, 50)
    k, c = s11.fit_scale(2.0 * truth + 9.0, truth, with_offset=True)
    assert k == pytest.approx(0.5, abs=1e-9) and c == pytest.approx(-4.5, abs=1e-8)


def test_scale_fitting_survives_an_all_zero_prediction():
    k, c = s11.fit_scale(np.zeros(10), np.arange(10.0), with_offset=False)
    assert k == 1.0 and c == 0.0
    k2, _ = s11.fit_scale_variance(np.zeros(10), np.arange(10.0))
    assert k2 == 1.0


# ---- the LODO leak that would have guaranteed a false success -------------------------------- #
def test_k_is_never_fitted_on_the_rows_it_is_scored_on():
    """Constructed so each donor needs a DIFFERENT k. If k leaked, every donor would be corrected
    perfectly; held out, each must be corrected by the OTHER donors' factor and stay wrong."""
    truth = np.array([1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0])
    donor = np.array(["A"] * 3 + ["B"] * 3 + ["C"] * 3)
    pred = np.concatenate([truth[:3] * 1.0, truth[3:6] * 10.0, truth[6:] * 100.0])
    cal, ks = s11.lodo_calibrate(pred, truth, donor, with_offset=False)
    assert not np.allclose(cal, truth), "a perfect correction means k leaked"
    assert len({round(v, 6) for v in ks.values()}) > 1


def test_lodo_uses_every_donor_exactly_once():
    rng = np.random.default_rng(4)
    truth = rng.normal(size=30)
    donor = np.repeat(["A", "B", "C"], 10)
    cal, ks = s11.lodo_calibrate(truth * 2.0, truth, donor, with_offset=False)
    assert set(ks) == {"A", "B", "C"} and cal.shape == truth.shape


def test_a_donor_with_too_little_training_data_passes_through_uncalibrated():
    truth = np.array([1.0, 2.0])
    donor = np.array(["A", "B"])
    cal, ks = s11.lodo_calibrate(np.array([5.0, 6.0]), truth, donor, with_offset=False)
    assert ks == {"A": 1.0, "B": 1.0} and np.allclose(cal, [5.0, 6.0])


def test_variance_mode_is_actually_reachable_through_lodo():
    rng = np.random.default_rng(5)
    truth = rng.normal(0, 10, 30)
    donor = np.repeat(["A", "B", "C"], 10)
    # the SAME prediction for both modes -- drawing fresh noise per call would compare two
    # different datasets and the inequality would be luck
    pred = 2.0 * truth + rng.normal(0, 6, 30)
    _, k_ls = s11.lodo_calibrate(pred, truth, donor, False, "ls")
    _, k_var = s11.lodo_calibrate(pred, truth, donor, False, "var")
    assert all(k_var[d] > k_ls[d] for d in k_ls)


# ---- verdict branches ----------------------------------------------------------------------- #
def test_reaching_the_floor_band_reads_as_scale_being_the_problem():
    assert s11.verdict_from(s11.FLOOR) == "SCALE IS THE PROBLEM"
    assert s11.verdict_from(s11.FLOOR * s11.FLOOR_MULT) == "SCALE IS THE PROBLEM"


def test_missing_the_band_reads_as_only_part_of_it():
    assert s11.verdict_from(s11.FLOOR * s11.FLOOR_MULT + 0.01) == "SCALE IS PART OF IT"


def test_a_nan_is_undetermined_not_a_verdict():
    assert s11.verdict_from(float("nan")) == "UNDETERMINED"


def test_the_bars_are_stated_constants():
    assert s11.FLOOR == 7.30 and s11.FLOOR_MULT == 1.5
    assert s11.K_STABILITY_BAR == 2.0


def test_the_stage_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_spearman_guards_degenerate_input():
    assert np.isnan(s11.spearman(np.ones(6), np.arange(6.0)))
    assert np.isnan(s11.spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))
    assert s11.spearman(np.arange(6.0), np.arange(6.0)) == pytest.approx(1.0)
    assert isinstance(pd.Series([1, 2]), pd.Series)
