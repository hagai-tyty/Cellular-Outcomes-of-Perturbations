"""STAGE 1.5 §8.3 E1 — every branch of the within-donor trajectory diagnostic.

Per the `verify_1a` lesson: a branch that never executes is not a check. Nothing here touches data;
`diag_e1_trajectory` keeps all repo-data imports inside `donor_trajectories`.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_e1_trajectory", _ROOT / "experiments" / "diag_e1_trajectory.py")
e1 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = e1
_SPEC.loader.exec_module(e1)


# ------------------------------- donor_trend ------------------------------- #
def test_donor_trend_is_minus_one_for_a_clean_decrease():
    # age falls monotonically as day rises -> perfect negative rank correlation
    assert e1.donor_trend([0, 7, 14, 21, 28], [50, 40, 30, 20, 10]) == pytest.approx(-1.0)


def test_donor_trend_is_plus_one_for_a_clean_increase():
    assert e1.donor_trend([0, 7, 14, 21, 28], [10, 20, 30, 40, 50]) == pytest.approx(1.0)


def test_donor_trend_is_nan_without_enough_points():
    assert np.isnan(e1.donor_trend([0, 7, 14], [3, 2, 1]))          # < MIN_POINTS_PER_DONOR


def test_donor_trend_is_nan_with_no_spread_in_age_or_day():
    assert np.isnan(e1.donor_trend([0, 7, 14, 21], [5, 5, 5, 5]))   # flat age
    assert np.isnan(e1.donor_trend([3, 3, 3, 3], [1, 2, 3, 4]))     # single day


def test_donor_trend_ignores_non_finite_points():
    r = e1.donor_trend([0, 7, np.nan, 21, 28], [50, 40, 99, 20, 10])
    assert r == pytest.approx(-1.0)                                 # the nan pair is dropped


# --------------------------------- e1_power -------------------------------- #
def test_e1_power_rises_with_the_true_effect():
    assert e1.e1_power(-0.2) < e1.e1_power(-0.6) < e1.e1_power(-0.9)


def test_e1_power_is_high_for_a_moderate_consistent_trend():
    assert e1.e1_power(-0.6, donor_sd=0.25, n=6) > 0.95


# -------------------------------- e1_verdict ------------------------------- #
def test_e1_verdict_pass_when_all_donors_fall_and_ci_excludes_zero():
    v = e1.e1_verdict({"N2": -0.8, "N3": -0.7, "O1": -0.9, "O2": -0.85, "Y1": -0.75, "Y2": -0.8})
    assert v["status"] == "PASS" and v["ci95"][1] < 0 and v["n_donors_negative"] == 6


def test_e1_verdict_wrong_direction_when_age_rises_consistently():
    v = e1.e1_verdict(dict(zip("ABCDEF", [0.7, 0.8, 0.75, 0.9, 0.72, 0.85], strict=True)))
    assert v["status"] == "WRONG_DIRECTION" and v["ci95"][0] > 0


def test_e1_verdict_no_trend_when_the_ci_straddles_zero():
    v = e1.e1_verdict(dict(zip("ABCDEF", [-0.6, 0.5, -0.4, 0.6, -0.3, 0.55], strict=True)))
    assert v["status"] == "NO_TREND" and v["ci95"][0] < 0 < v["ci95"][1]


def test_e1_verdict_cannot_verify_below_three_donors():
    assert e1.e1_verdict({"N2": -0.9, "O1": -0.8})["status"] == "CANNOT_VERIFY"


def test_e1_verdict_drops_donors_with_an_undefined_trend():
    v = e1.e1_verdict({"N2": float("nan"), "N3": -0.7, "O1": -0.8, "O2": -0.9})
    assert v["n_donors"] == 3 and "N2" not in v["per_donor"]


# ---------------------------------- bars ----------------------------------- #
def test_bars_registers_e1_and_reports_it_resolvable_for_a_moderate_trend():
    b = {x["id"]: x for x in e1.bars()}
    assert set(b) == {"E1"}
    assert b["E1"]["pass_rate_if_intent_holds"] > 0.95
