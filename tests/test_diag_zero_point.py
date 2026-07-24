"""STAGE 1.5 Phase 1 — every branch of the zero-point diagnostics.

Per §5.5 and the `verify_1a` lesson: a branch that never executes is not a check. Every status
each function can emit is driven here, including the ones we HOPE not to see (M1 FAIL, which
escalates past this whole stage) and the ones we EXPECT to see (M3 INDETERMINATE at n=6).

Nothing here touches data; `diag_zero_point` keeps all repo-data imports inside `baseline_ages`.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_zero_point", _ROOT / "experiments" / "diag_zero_point.py")
dzp = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = dzp
_SPEC.loader.exec_module(dzp)

CV = dzp.CLOCK_CV_MAE


# ------------------------------------------------------------------ the bar ---- #
def test_m1_threshold_is_not_merely_greater_than_zero():
    """A 'contrast > 0' bar would be passed half the time by a clock reading pure noise."""
    assert dzp.m1_threshold() > 15.0
    # exactly z_0.95 * cv_mae * sqrt(1/2 + 1/2)
    assert dzp.m1_threshold() == pytest.approx(1.6448536269514722 * CV, rel=1e-9)


def test_m1_threshold_falls_with_more_baselines():
    assert dzp.m1_threshold(n_young=6, n_old=6) < dzp.m1_threshold(n_young=2, n_old=2)


def test_m1_bar_is_resolvable_at_the_real_geometry():
    """§5b: a correct clock must clear its own bar >= 95% of the time. This is why M1 is gated
    and the 29-vs-35 middle contrast is not."""
    assert dzp.m1_power(53.0) >= 0.95            # the 0-vs-53 extreme
    assert dzp.m1_power(6.0) < 0.5               # 29 vs 35 — correctly NOT gated


def test_bars_are_recorded_before_the_run_and_M1_reads_resolvable():
    ids = {b["id"]: b for b in dzp.bars()}
    assert set(ids) == {"M1", "M2", "M3"}
    assert ids["M1"]["verdict"] == "RESOLVABLE"
    assert "CONDITIONAL" in ids["M2"]["verdict"]
    assert "UNRESOLVABLE" in ids["M3"]["verdict"]     # expected, and pre-registered as such


# ----------------------------------------------------------------------- M1 ---- #
def test_m1_passes_when_the_clock_tracks_age():
    pred = {"N2": 5.0, "N3": 2.0, "Y1": 30.0, "Y2": 33.0, "O1": 55.0, "O2": 50.0}
    v = dzp.m1_verdict(pred)
    assert v["status"] == "PASS"
    assert v["contrast_years"] == pytest.approx((55.0 + 50.0) / 2 - (5.0 + 2.0) / 2)
    assert v["true_age_gap"] == 53.0


def test_m1_fails_when_the_clock_reads_nothing():
    """The escalation case: no separation between a 0-year-old and a 53-year-old donor."""
    pred = {"N2": 40.0, "N3": 41.0, "Y1": 39.0, "Y2": 40.5, "O1": 40.2, "O2": 39.8}
    assert dzp.m1_verdict(pred)["status"] == "FAIL"


def test_m1_fails_when_separation_is_real_but_below_the_noise_bar():
    """Ordered correctly yet only ~10 yr apart — indistinguishable from a null clock."""
    pred = {"N2": 30.0, "N3": 30.0, "Y1": 35.0, "Y2": 35.0, "O1": 40.0, "O2": 40.0}
    v = dzp.m1_verdict(pred)
    assert v["contrast_years"] == pytest.approx(10.0)
    assert v["status"] == "FAIL"


def test_m1_cannot_verify_without_overlapping_donors():
    assert dzp.m1_verdict({"ZZ": 10.0})["status"] == "CANNOT_VERIFY"


def test_m1_cannot_verify_without_an_age_contrast():
    assert dzp.m1_verdict({"N2": 10.0, "N3": 12.0})["status"] == "CANNOT_VERIFY"


def test_m1_reports_spearman_but_does_not_gate_on_it():
    """Perfect extreme separation with the MIDDLE pair inverted still PASSES — the middle is
    underpowered by design and must not be able to fail the gate."""
    pred = {"N2": 0.0, "N3": 1.0, "Y1": 40.0, "Y2": 35.0, "O1": 53.0, "O2": 54.0}
    v = dzp.m1_verdict(pred)
    assert v["status"] == "PASS"
    assert np.isfinite(v["spearman_all_donors_REPORTED_NOT_GATED"])


# ----------------------------------------------------------------------- M2 ---- #
def test_m2_not_estimable_with_no_pairs_makes_option_a_impossible():
    """Tightening T4 — a permitted outcome, not a failure to measure."""
    v = dzp.m2_verdict([])
    assert v["status"] == "NOT_ESTIMABLE" and v["n_pairs"] == 0
    assert "option (a) is impossible" in v["reason"]


def test_m2_indeterminate_below_the_pair_floor():
    v = dzp.m2_verdict([1.0, 2.0])
    assert v["status"] == "INDETERMINATE" and v["n_pairs"] == 2


def test_m2_detects_a_consistent_batch_offset():
    v = dzp.m2_verdict([8.0, 7.5, 8.2, 7.9, 8.1])
    assert v["status"] == "BATCH_EFFECT"
    assert v["ci95"][0] > 0


def test_m2_reports_no_effect_when_the_ci_straddles_zero():
    v = dzp.m2_verdict([5.0, -4.0, 3.0, -6.0, 2.0])
    assert v["status"] == "NO_BATCH_EFFECT"
    assert v["ci95"][0] < 0 < v["ci95"][1]


def test_m2_ignores_non_finite_pairs():
    assert dzp.m2_verdict([8.0, np.nan, 7.5, 8.2, np.inf, 7.9, 8.1])["n_pairs"] == 5


# ----------------------------------------------------------------------- M3 ---- #
def test_m3_is_indeterminate_at_the_real_six_donor_geometry():
    """THE pre-registered expectation. The real level shifts give a point estimate near 56%,
    but with 6 donors the chi-square CI spans most of [0,1] — so the honest answer is
    'indeterminate', decided before the run rather than rationalised after it."""
    real = [15.025, -28.346, 0.640, 6.557, -8.134, -20.017]
    v = dzp.m3_verdict(real)
    assert v["n_donors"] == 6
    assert v["observed_sd_years"] == pytest.approx(16.39, abs=0.05)
    assert v["share_of_variance"] == pytest.approx(0.56, abs=0.02)
    assert v["share_ci"][1] - v["share_ci"][0] > 0.5      # too wide to conclude
    assert v["status"] == "INDETERMINATE"


def test_m3_says_baseline_dominates_when_the_spread_is_small_and_well_determined():
    rng = np.random.default_rng(0)
    v = dzp.m3_verdict(list(rng.normal(0, CV, 60)))       # spread == baseline noise, many donors
    assert v["status"] == "BASELINE_DOMINATES"
    assert v["share_of_variance"] >= dzp.M3_MOSTLY


def test_m3_says_baseline_minor_when_the_offset_dwarfs_the_clock_error():
    rng = np.random.default_rng(1)
    v = dzp.m3_verdict(list(rng.normal(0, CV * 4, 60)))
    assert v["status"] == "BASELINE_MINOR"
    assert v["share_of_variance"] < dzp.M3_LITTLE


def test_m3_residual_is_the_variance_left_after_removing_baseline_noise():
    real = [15.025, -28.346, 0.640, 6.557, -8.134, -20.017]
    v = dzp.m3_verdict(real)
    assert v["residual_sd_years"] == pytest.approx(10.9, abs=0.2)


def test_m3_cannot_verify_with_too_few_donors_or_zero_spread():
    assert dzp.m3_verdict([1.0, 2.0])["status"] == "CANNOT_VERIFY"
    assert dzp.m3_verdict([4.0, 4.0, 4.0])["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------------------ decide ---- #
def test_m1_failure_escalates_past_this_stage():
    d = dzp.decide({"status": "FAIL", "reason": "x"}, {"status": "NO_BATCH_EFFECT"},
                   {"status": "BASELINE_MINOR"})
    assert d["action"] == "ESCALATE" and "Stage 4" in d["reason"]


def test_m1_inconclusive_blocks():
    assert dzp.decide({"status": "CANNOT_VERIFY", "reason": "r"}, {}, {})["action"] == "BLOCKED"


def test_clean_baselines_need_only_phase_2():
    d = dzp.decide({"status": "PASS"}, {"status": "NO_BATCH_EFFECT"},
                   {"status": "BASELINE_MINOR"})
    assert d["action"] == "PHASE_2_ONLY"


def test_a_batch_effect_pulls_in_phase_3_led_by_option_a():
    d = dzp.decide({"status": "PASS"}, {"status": "BATCH_EFFECT"}, {"status": "BASELINE_MINOR"})
    assert d["action"] == "PHASE_2_AND_3"
    assert d["phase3_lead"].startswith("(a)")


def test_baseline_dominance_pulls_in_phase_3_led_by_option_b():
    d = dzp.decide({"status": "PASS"}, {"status": "NO_BATCH_EFFECT"},
                   {"status": "BASELINE_DOMINATES"})
    assert d["action"] == "PHASE_2_AND_3"
    assert d["phase3_lead"].startswith("(b)")


def test_phase_3_states_that_it_reopens_stage_1(  ):
    """Tightening T4: changing y_age moves BOTH Stage 1 targets, not just the guards."""
    d = dzp.decide({"status": "PASS"}, {"status": "BATCH_EFFECT"}, {"status": "INDETERMINATE"})
    assert "reopens BOTH" in d["reason"] and "guard record" in d["reason"]


def test_the_expected_real_path_is_phase_2_and_3_with_no_lead():
    """M1 pass + M2 not estimable + M3 indeterminate = the outcome Phase 1 most likely returns."""
    d = dzp.decide({"status": "PASS"}, {"status": "NOT_ESTIMABLE"}, {"status": "INDETERMINATE"})
    assert d["action"] == "PHASE_2_AND_3"
    assert "undetermined by M3" in d["phase3_lead"]


# ------------------------------------------------------------ level shifts ---- #
def test_load_level_shifts_returns_empty_when_no_snapshot(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert dzp.load_level_shifts("scorecard/nope.json") == {}
