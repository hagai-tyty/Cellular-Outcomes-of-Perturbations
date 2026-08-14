"""Unit tests for the ΔAge run on GSE165177 — pure functions only, no GEO files.

`arm_group` decides which samples become the ZERO-POINT and which become the measurement, so a
mis-grouping silently changes every ΔAge in the run rather than raising. The ordering of its
branches is load-bearing: `negative_control_intermediate` contains both "negative_control" and
"intermediate", and `failing_to_transiently_reprogram_intermediate` contains BOTH "fail" and
"transient" — so it must be matched as *failed*, not *transient*.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


A = _load("dage_gse165177", "experiments/dage_gse165177.py")


# ------------------------------------------------------------------------- arm_group ---- #
@pytest.mark.parametrize("arm,want", [("negative_control", "control_fib"),
                                      ("negative_control_intermediate", "control_int")])
def test_the_two_control_states_are_kept_apart(arm, want):
    """Pooling them gave a mixed reference and cost M-E3 its verdict. Regression."""
    assert A.arm_group(arm) == want


@pytest.mark.parametrize("arm,want", [
    ("failed_to_transiently_reprogram", "failed_fib"),
    ("failing_to_transiently_reprogram_intermediate", "failed_int")])
def test_failed_arms_group_as_failed_even_though_they_contain_the_word_transiently(arm, want):
    """The branch-order trap: these contain BOTH 'fail' and 'transient'. Grouping them as
    *transient* would move 33 samples into the treatment arm and invert M-E3."""
    assert A.arm_group(arm) == want


@pytest.mark.parametrize("arm,want", [
    ("transiently_reprogrammed", "transient_fib"),
    ("transient_reprogramming", "transient_fib"),
    ("transient_reprogramming_intermediate", "transient_int")])
def test_returned_fibroblasts_and_intermediates_are_kept_apart(arm, want):
    """THE regression this file exists for. These two populations differ by ~40 yr in our own
    data and ~18 yr in the methylation companion. Pooling them produced a -42 that dissolved to
    -24 once separated, flipping M-E3 from REPRODUCED to NOT REPRODUCED."""
    assert A.arm_group(arm) == want


def test_intermediate_and_fibroblast_never_share_a_group():
    fib = {A.arm_group(a) for a in ("transiently_reprogrammed", "failed_to_transiently_reprogram",
                                    "negative_control")}
    ints = {A.arm_group(a) for a in ("transient_reprogramming_intermediate",
                                     "failing_to_transiently_reprogram_intermediate",
                                     "negative_control_intermediate")}
    assert fib.isdisjoint(ints)


@pytest.mark.parametrize("group,want", [
    ("transient_fib", "control_fib"), ("failed_fib", "control_fib"),
    ("transient_int", "control_int"), ("failed_int", "control_int")])
def test_each_treated_stratum_pairs_to_its_own_matched_control(group, want):
    """Fibroblasts against fibroblast controls, intermediates against intermediate controls.
    A pooled reference is what made the first run's contrast not like-for-like."""
    assert A.control_for(group) == want


def test_day_zero_fibroblasts_are_their_own_group():
    """They are untreated, but they are NOT a contemporaneous control -- they sit at day 0 while
    every control sits at days 10-17. Pooling them into `control` would reintroduce exactly the
    cross-timepoint zero-point that D1 is about."""
    assert A.arm_group("day0_fibroblast") == "day0"


def test_an_unrecognised_arm_is_not_silently_absorbed():
    """It must land in `other` and be excluded, never default into control or treatment."""
    assert A.arm_group("something_new") == "other"


def test_grouping_is_case_insensitive():
    assert A.arm_group("Negative_Control") == "control_fib"
    assert A.arm_group("FAILED_to_reprogram") == "failed_fib"
    assert A.arm_group("Transient_Reprogramming_INTERMEDIATE") == "transient_int"


# -------------------------------------------------------------------------------- ci ---- #
def test_ci_matches_the_t_interval_by_hand():
    m, lo, hi, n = A.ci([1.0, 2.0, 3.0, 4.0])
    assert (m, n) == (2.5, 4)
    se = float(np.std([1.0, 2.0, 3.0, 4.0], ddof=1)) / 2.0
    assert lo == pytest.approx(2.5 - A.T_CRIT[3] * se)
    assert hi == pytest.approx(2.5 + A.T_CRIT[3] * se)


def test_ci_uses_the_t_critical_value_for_n_minus_one():
    """n=3 must use t(.975, df=2)=4.303, not 1.96 — at these sample sizes it is a 2.2x
    difference in the interval width and decides whether a bar fires."""
    _, lo, hi, n = A.ci([10.0, 12.0, 14.0])
    se = float(np.std([10.0, 12.0, 14.0], ddof=1)) / np.sqrt(3)
    assert n == 3
    assert (hi - lo) / 2 == pytest.approx(A.T_CRIT[2] * se)
    assert A.T_CRIT[2] == 4.303


def test_ci_drops_non_finite_values_rather_than_propagating_nan():
    m, lo, hi, n = A.ci([1.0, np.nan, 3.0, np.inf])
    assert n == 2 and m == pytest.approx(2.0)
    assert np.isfinite(lo) and np.isfinite(hi)


def test_ci_on_a_single_value_returns_it_with_no_interval():
    m, lo, hi, n = A.ci([7.0])
    assert (m, n) == (7.0, 1)
    assert np.isnan(lo) and np.isnan(hi)


def test_ci_on_an_empty_input_is_all_nan():
    m, lo, hi, n = A.ci([])
    assert n == 0 and np.isnan(m) and np.isnan(lo) and np.isnan(hi)


def test_ci_on_identical_values_has_zero_width():
    m, lo, hi, n = A.ci([5.0, 5.0, 5.0])
    assert m == 5.0 and lo == pytest.approx(5.0) and hi == pytest.approx(5.0)


def test_ci_falls_back_to_the_normal_quantile_beyond_the_table():
    """T_CRIT is tabulated to df=10; larger n must not KeyError."""
    _, lo, hi, n = A.ci(list(np.arange(30.0)))
    assert n == 30 and np.isfinite(lo) and np.isfinite(hi)


# ----------------------------------------------------------------------- true ages ---- #
def test_the_donor_ages_are_the_ones_geo_declares():
    """Pinned so a future edit cannot quietly move the M-E1 calibration target."""
    assert A.TRUE_AGES == {"O1": 53.0, "O2": 53.0, "O3": 38.0}
    assert float(np.mean(list(A.TRUE_AGES.values()))) == pytest.approx(48.0)
