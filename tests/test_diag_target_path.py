"""Unit tests for the target-path audit / scale-mismatch test.

Only the pure functions. The measurement itself needs built folds and a bundle, so the parts
tested here are the ones that decide a VERDICT -- those must not be discoverable-as-wrong after
the fact.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("dtp", ROOT / "experiments" / "diag_target_path.py")
dtp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dtp)


# ---- slope_of ----------------------------------------------------------------------------- #
def test_slope_recovers_a_known_compression():
    """A model predicting exactly half the truth has slope 0.5 -- the regression form of the
    same claim `compression` makes as an SD ratio."""
    t = np.linspace(-30, 30, 60)
    assert dtp.slope_of(t, 0.5 * t) == pytest.approx(0.5, abs=1e-9)


def test_slope_of_a_perfect_predictor_is_one():
    t = np.linspace(-30, 30, 60)
    assert dtp.slope_of(t, t) == pytest.approx(1.0, abs=1e-9)


def test_slope_of_a_constant_prediction_is_zero():
    """The degenerate case that matters here: a model predicting ~0 regardless of the truth."""
    t = np.linspace(-30, 30, 60)
    assert dtp.slope_of(t, np.zeros_like(t)) == pytest.approx(0.0, abs=1e-9)


def test_slope_is_undefined_when_the_truth_is_constant():
    assert np.isnan(dtp.slope_of(np.full(10, 5.0), np.arange(10, dtype=float)))


def test_slope_is_undefined_on_too_few_points():
    assert np.isnan(dtp.slope_of(np.array([1.0, 2.0]), np.array([1.0, 2.0])))


def test_sd_needs_two_points():
    assert np.isnan(dtp.sd(np.array([1.0])))
    assert dtp.sd(np.array([1.0, 3.0])) == pytest.approx(np.sqrt(2.0))


# ---- hypothesis_verdict ------------------------------------------------------------------- #
def test_uniform_worsening_below_the_ceiling_is_supported():
    old = {"a": 0.9, "b": 0.9, "c": 0.9}
    new = {"a": 0.4, "b": 0.4, "c": 0.4}
    v = dtp.hypothesis_verdict(old, new)
    assert v["verdict"] == "H-SUPPORTED"
    assert v["n_more_compressed"] == 3 and v["n_folds"] == 3


def test_uniform_worsening_that_stays_above_the_ceiling_is_refuted():
    """Direction alone must not carry it -- the effect has to be materially large."""
    old = {"a": 0.99, "b": 0.99, "c": 0.99}
    new = {"a": 0.95, "b": 0.95, "c": 0.95}
    assert dtp.hypothesis_verdict(old, new)["verdict"] == "H-REFUTED"


def test_a_low_median_without_a_majority_is_refuted():
    """Both clauses are required: one fold collapsing must not carry the verdict."""
    old = {"a": 0.5, "b": 0.5, "c": 0.9}
    new = {"a": 0.6, "b": 0.6, "c": 0.1}
    v = dtp.hypothesis_verdict(old, new)
    assert v["n_more_compressed"] == 1
    assert v["verdict"] == "H-REFUTED"


def test_improvement_everywhere_is_refuted():
    old = {"a": 0.4, "b": 0.4, "c": 0.4}
    new = {"a": 0.9, "b": 0.9, "c": 0.9}
    assert dtp.hypothesis_verdict(old, new)["verdict"] == "H-REFUTED"


def test_only_folds_present_in_both_arms_are_compared():
    """N2 has no age-valid test cell under C-7, so it must drop out rather than be counted."""
    old = {"N2": 0.79, "N3": 0.61, "O1": 0.83}
    new = {"N3": 0.69, "O1": 0.44}
    v = dtp.hypothesis_verdict(old, new)
    assert v["folds"] == ["N3", "O1"] and v["n_folds"] == 2


def test_no_shared_folds_is_undetermined_not_a_verdict():
    assert dtp.hypothesis_verdict({"a": 0.5}, {"b": 0.5})["verdict"] == "UNDETERMINED"


def test_the_ceiling_is_a_stated_constant():
    assert dtp.COMPRESSION_CEILING == 0.80


def test_the_recorded_result_reproduces_from_its_own_numbers():
    """Guards the reported verdict against a later silent change to the rule."""
    old = {"N3": 0.609, "O1": 0.826, "O2": 0.880, "Y1": 0.801, "Y2": 0.905}
    new = {"N3": 0.687, "O1": 0.437, "O2": 0.428, "Y1": 0.825, "Y2": 0.534}
    v = dtp.hypothesis_verdict(old, new)
    assert v["verdict"] == "H-SUPPORTED"
    assert v["n_more_compressed"] == 3 and v["n_folds"] == 5
    assert v["median_compression_pre"] == pytest.approx(0.826)
    assert v["median_compression_c7"] == pytest.approx(0.534)
