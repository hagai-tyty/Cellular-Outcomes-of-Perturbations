"""Every branch of the A1/A3-corrected E1 re-run.

This decides whether the current ΔAge labels are salvageable, so each verdict path is driven
explicitly — including FRAGILE, which exists because E1b (+0.009) and D2 (−0.014) were both
decided by hundredths.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_e1_corrected", _ROOT / "experiments" / "diag_e1_corrected.py")
ec = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = ec
_SPEC.loader.exec_module(ec)


# ------------------------------------------------------------------ paired_ci ---- #
def test_paired_ci_matches_a_hand_worked_case():
    s = ec.paired_ci([-10.0, -20.0, -30.0])
    assert s["n"] == 3 and s["mean"] == pytest.approx(-20.0)
    assert s["n_negative"] == 3
    assert s["ci95"][0] < -20.0 < s["ci95"][1]


def test_paired_ci_drops_non_finite():
    assert ec.paired_ci([1.0, np.nan, 2.0, np.inf])["n"] == 2


def test_paired_ci_handles_too_few():
    assert ec.paired_ci([1.0])["n"] == 1


# -------------------------------------------------------------- window_verdict ---- #
def test_rejuvenation_when_ci_is_entirely_negative():
    v = ec.window_verdict(ec.paired_ci([-30.0, -28.0, -32.0, -29.0]))
    assert v["status"] == "REJUVENATION"
    assert v["ci95"][1] < 0


def test_ageing_when_ci_is_entirely_positive():
    v = ec.window_verdict(ec.paired_ci([30.0, 28.0, 32.0, 29.0]))
    assert v["status"] == "AGEING"


def test_no_effect_when_ci_straddles_zero():
    v = ec.window_verdict(ec.paired_ci([-30.0, 28.0, -32.0, 29.0]))
    assert v["status"] == "NO_EFFECT"


def test_fragile_is_flagged_when_a_bound_hugs_zero():
    """THE §10 lesson: a verdict decided by hundredths must say so in the same breath."""
    # constructed so the upper bound lands just below 0
    v = ec.window_verdict({"n": 6, "mean": -5.0, "ci95": [-9.8, -0.2], "n_negative": 5})
    assert v["status"] == "REJUVENATION_FRAGILE"
    assert "FRAGILE" in v["reason"]


def test_not_fragile_when_the_ci_clears_zero_comfortably():
    v = ec.window_verdict({"n": 6, "mean": -30.0, "ci95": [-50.0, -10.0], "n_negative": 6})
    assert v["status"] == "REJUVENATION"


def test_cannot_verify_with_one_donor():
    assert ec.window_verdict(ec.paired_ci([1.0]))["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------- leave_one_donor_out ---- #
def test_loo_stable_when_all_donors_agree():
    r = ec.leave_one_donor_out({"A": -10.0, "B": -12.0, "C": -11.0, "D": -9.0})
    assert r["status"] == "STABLE"


def test_loo_detects_one_donor_carrying_the_effect():
    """Five near-zero donors and one huge one: dropping it flips the mean's sign."""
    r = ec.leave_one_donor_out({"A": -100.0, "B": 2.0, "C": 3.0, "D": 1.0, "E": 2.0, "F": 3.0})
    assert r["status"] == "SIGN_FLIPS"


def test_loo_needs_three_donors():
    assert ec.leave_one_donor_out({"A": 1.0, "B": 2.0})["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------------------- decide ---- #
def test_rejuvenation_means_labels_are_fine_and_no_control_swap():
    d = ec.decide({"status": "REJUVENATION"}, {}, {"status": "STABLE"})
    assert d["action"] == "LABELS_ARE_FINE"
    assert "do NOT redefine is_control" in d["reason"]


def test_fragile_rejuvenation_still_counts_as_labels_fine():
    d = ec.decide({"status": "REJUVENATION_FRAGILE"}, {}, {"status": "STABLE"})
    assert d["action"] == "LABELS_ARE_FINE"


def test_ageing_means_labels_inadequate_and_indicts_the_day0_control_too():
    d = ec.decide({"status": "AGEING"}, {}, {"status": "STABLE"})
    assert d["action"] == "LABELS_INADEQUATE"
    assert "evidence against the day-0 control" in d["reason"]


def test_no_effect_is_inconclusive_and_reports_loo_stability():
    d = ec.decide({"status": "NO_EFFECT"}, {}, {"status": "SIGN_FLIPS"})
    assert d["action"] == "INCONCLUSIVE"
    assert "n=6 is the binding constraint" in d["reason"]


def test_no_effect_with_stable_loo_says_it_is_a_real_null():
    d = ec.decide({"status": "NO_EFFECT"}, {}, {"status": "STABLE"})
    assert "real null" in d["reason"]


def test_other_windows_are_reported_but_never_selected_on():
    """Guard against fishing: a rejuvenating side-window must be reported, not promoted."""
    d = ec.decide({"status": "AGEING"}, {"13-15": {"status": "REJUVENATION"}}, {"status": "STABLE"})
    assert d["action"] == "LABELS_INADEQUATE"          # NOT flipped by the side window
    assert "13-15" in d["reason"] and "do not select" in d["reason"]


# --------------------------------------------------------------- delta_by_donor ---- #
def test_delta_by_donor_uses_day0_baseline_and_filters_arm_and_window():
    ages = {"A_base": 50.0, "A_r1": 40.0, "A_r2": 30.0, "A_late": 99.0, "A_nr": 80.0,
            "B_base": 60.0, "B_r1": 55.0}
    ct = {"A_base": ec.BASELINE, "A_r1": ec.RESPONDER, "A_r2": ec.RESPONDER,
          "A_late": ec.RESPONDER, "A_nr": ec.NON_RESPONDER,
          "B_base": ec.BASELINE, "B_r1": ec.RESPONDER}
    day = {"A_base": 0.0, "A_r1": 11.0, "A_r2": 13.0, "A_late": 29.0, "A_nr": 11.0,
           "B_base": 0.0, "B_r1": 11.0}
    donor = {k: k.split("_")[0] for k in ages}
    out = ec.delta_by_donor(ages, ct, day, donor, list(ages), ec.RESPONDER, (10.0, 13.0))
    assert out["A"] == pytest.approx((40.0 + 30.0) / 2 - 50.0)   # in-window responders only
    assert out["B"] == pytest.approx(55.0 - 60.0)


def test_delta_by_donor_skips_a_donor_with_no_baseline():
    ages = {"A_r1": 40.0}
    out = ec.delta_by_donor(ages, {"A_r1": ec.RESPONDER}, {"A_r1": 11.0}, {"A_r1": "A"},
                            list(ages), ec.RESPONDER, (10.0, 13.0))
    assert out == {}


def test_the_peak_window_is_pre_committed():
    assert ec.PEAK_WINDOW == (10.0, 13.0)
