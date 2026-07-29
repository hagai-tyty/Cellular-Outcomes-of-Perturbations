"""Every branch of the methylation anchor, plus the Horvath transform against known values.

This diagnostic decides whether ΔAge has a valid anchor, so the maths is checked against published
fixed points rather than trusted: `F(20) = 0` and `F(0) = −log(21)` are properties of Horvath's
transform, not of my implementation, and `anti_trafo` must invert `trafo` exactly.
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_methylation_anchor", _ROOT / "experiments" / "diag_methylation_anchor.py")
ma = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = ma
_SPEC.loader.exec_module(ma)


# ------------------------------------------------- the Horvath transform ---- #
def test_trafo_fixed_points_from_the_published_definition():
    """F(adult_age) = 0 and F(0) = -log(21). These are the clock's own anchors."""
    assert float(ma.trafo(20.0)) == pytest.approx(0.0, abs=1e-12)
    assert float(ma.trafo(0.0)) == pytest.approx(-np.log(21.0), abs=1e-12)


def test_anti_trafo_inverts_trafo_across_the_lifespan():
    ages = np.array([0.0, 1.0, 5.0, 19.9, 20.0, 20.1, 38.0, 53.0, 96.0])
    assert np.allclose(ma.anti_trafo(ma.trafo(ages)), ages, atol=1e-9)


def test_anti_trafo_matches_the_published_branch_form():
    # x >= 0 branch is linear with slope 21 and offset 20
    assert float(ma.anti_trafo(0.0)) == pytest.approx(20.0)
    assert float(ma.anti_trafo(1.0)) == pytest.approx(41.0)
    # x < 0 branch is exponential
    assert float(ma.anti_trafo(-np.log(21.0))) == pytest.approx(0.0, abs=1e-9)


def test_a_53_year_old_transforms_and_back():
    assert float(ma.anti_trafo(ma.trafo(53.0))) == pytest.approx(53.0, abs=1e-9)
    assert float(ma.trafo(53.0)) == pytest.approx((53.0 - 20.0) / 21.0)


# ------------------------------------------------------ linear_predictor ---- #
def test_linear_predictor_sums_only_present_cpgs():
    w = {"cg1": 2.0, "cg2": -1.0, "cg3": 10.0}
    s, n = ma.linear_predictor({"cg1": 0.5, "cg2": 0.25}, w)      # cg3 absent
    assert n == 2
    assert s == pytest.approx(2.0 * 0.5 + (-1.0) * 0.25)


def test_linear_predictor_ignores_non_finite_betas():
    s, n = ma.linear_predictor({"cg1": np.nan, "cg2": 0.5}, {"cg1": 5.0, "cg2": 2.0})
    assert n == 1 and s == pytest.approx(1.0)


# ----------------------------------------------------- implied_intercept ---- #
def test_implied_intercept_is_exact_when_the_clock_is_consistent():
    """Construct linear predictors that DO have a common intercept; it must be recovered with
    zero spread — the property that makes G2 meaningful."""
    true = [53.0, 53.0, 38.0]
    k = -0.4471
    lp = [float(ma.trafo(a)) - k for a in true]
    r = ma.implied_intercept(lp, true)
    assert r["status"] == "OK"
    assert r["mean"] == pytest.approx(k, abs=1e-9)
    assert r["sd"] == pytest.approx(0.0, abs=1e-9)


def test_implied_intercept_exposes_an_inconsistent_clock():
    """If no single intercept fits, the spread is large — which is what disqualifies an anchor."""
    r = ma.implied_intercept([0.0, 0.0, 0.0], [53.0, 20.0, 1.0])
    assert r["status"] == "OK"
    assert r["spread_years"] > 10.0


def test_implied_intercept_needs_two_samples():
    assert ma.implied_intercept([0.1], [53.0])["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------------ g2_verdict ---- #
def test_g2_reproduces_when_predictions_are_close():
    v = ma.g2_verdict([52.0, 54.0, 39.0], [53.0, 53.0, 38.0])
    assert v["status"] == "REPRODUCES" and v["mae_years"] == pytest.approx(1.0)


def test_g2_fails_when_predictions_are_far():
    v = ma.g2_verdict([20.0, 25.0, 80.0], [53.0, 53.0, 38.0])
    assert v["status"] == "FAILS"
    assert "NOT an anchor" in v["reason"]


def test_g2_boundary_is_the_declared_tolerance():
    # MAE is the mean absolute error, so both samples must sit at the tolerance to be AT the bar
    tol = ma.G2_MAE_TOL
    assert ma.g2_verdict([53.0 + tol, 53.0 + tol], [53.0, 53.0])["status"] == "REPRODUCES"
    assert ma.g2_verdict([53.0 + tol + 0.1, 53.0 + tol + 0.1], [53.0, 53.0])["status"] == "FAILS"


def test_g2_cannot_verify_with_one_sample():
    assert ma.g2_verdict([53.0], [53.0])["status"] == "CANNOT_VERIFY"


# -------------------------------------------------------- effect_verdict ---- #
def test_rejuvenation_when_ci_entirely_negative():
    v = ma.effect_verdict(ma.paired_stat([-30.0, -28.0, -32.0, -29.0]), "x")
    assert v["status"] == "REJUVENATION" and v["ci95"][1] < 0


def test_ageing_when_ci_entirely_positive():
    assert ma.effect_verdict(ma.paired_stat([30.0, 28.0, 32.0, 29.0]), "x")["status"] == "AGEING"


def test_no_effect_when_ci_straddles_zero():
    assert ma.effect_verdict(ma.paired_stat([-30.0, 31.0, -29.0, 28.0]), "x")["status"] == "NO_EFFECT"


def test_fragile_flag_when_a_bound_hugs_zero():
    v = ma.effect_verdict({"n": 9, "mean": -5.0, "ci95": [-9.8, -0.2], "n_negative": 8}, "x")
    assert v["status"] == "REJUVENATION_FRAGILE" and "FRAGILE" in v["reason"]


def test_effect_cannot_verify_with_one_pair():
    assert ma.effect_verdict(ma.paired_stat([1.0]), "x")["status"] == "CANNOT_VERIFY"


# --------------------------------------------------------------- decide ---- #
def test_g2_failure_blocks_interpretation_of_everything():
    d = ma.decide({"status": "FAILS"}, {"status": "REJUVENATION"}, {"status": "NO_EFFECT"})
    assert d["action"] == "ANCHOR_INVALID"
    assert "NOT interpreted" in d["reason"]


def test_clean_rejuvenation_with_inert_control_is_the_win_condition():
    d = ma.decide({"status": "REPRODUCES"}, {"status": "REJUVENATION"}, {"status": "NO_EFFECT"})
    assert d["action"] == "ANCHOR_VALID_EFFECT_REAL"


def test_rejuvenation_with_a_moving_control_is_flagged_not_hidden():
    d = ma.decide({"status": "REPRODUCES"}, {"status": "REJUVENATION"}, {"status": "AGEING"})
    assert d["action"] == "ANCHOR_VALID_EFFECT_REAL_CONTROL_MOVES"
    assert "not inert" in d["reason"]


def test_ageing_triggers_a_bug_hunt_not_a_premise_revision():
    d = ma.decide({"status": "REPRODUCES"}, {"status": "AGEING"}, {"status": "NO_EFFECT"})
    assert d["action"] == "CONTRADICTS" and "bug hunt" in d["reason"]


def test_null_is_a_real_finding_not_a_power_failure():
    d = ma.decide({"status": "REPRODUCES"}, {"status": "NO_EFFECT"}, {"status": "NO_EFFECT"})
    assert d["action"] == "NO_EFFECT_AT_THIS_RESOLUTION"
    assert "not a power failure" in d["reason"]


def test_fragile_rejuvenation_still_counts_as_an_effect():
    d = ma.decide({"status": "REPRODUCES"}, {"status": "REJUVENATION_FRAGILE"},
                  {"status": "NO_EFFECT"})
    assert d["action"] == "ANCHOR_VALID_EFFECT_REAL"


# ------------------------------------------------------ pair_by_donor_day ---- #
def test_pairing_matches_on_donor_and_day_only():
    meta = {
        "a": {"donor": "O1", "day": 13.0, "ctype": ma.TR, "age": 53.0},
        "b": {"donor": "O1", "day": 13.0, "ctype": ma.NC, "age": 53.0},
        "c": {"donor": "O2", "day": 15.0, "ctype": ma.TR, "age": 53.0},
        "d": {"donor": "O2", "day": 17.0, "ctype": ma.NC, "age": 53.0},   # day mismatch
    }
    assert ma.pair_by_donor_day(meta, ma.TR, ma.NC) == [("O1", 13.0, ["a"], ["b"])]


def test_pairing_KEEPS_exp1_exp2_replicates_for_averaging():
    """REGRESSION. The first version required a unique sample per (donor, day, arm), which
    silently discarded 6 of 9 M-1 pairs -- GSE165179 runs every condition as exp1 AND exp2.
    Replicates are repeats of one condition and must be averaged, never dropped."""
    meta = {
        "t1": {"donor": "O1", "day": 13.0, "ctype": ma.TR, "age": 53.0},
        "t2": {"donor": "O1", "day": 13.0, "ctype": ma.TR, "age": 53.0},
        "c1": {"donor": "O1", "day": 13.0, "ctype": ma.NC, "age": 53.0},
        "c2": {"donor": "O1", "day": 13.0, "ctype": ma.NC, "age": 53.0},
    }
    got = ma.pair_by_donor_day(meta, ma.TR, ma.NC)
    assert got == [("O1", 13.0, ["t1", "t2"], ["c1", "c2"])]


def test_pairing_drops_a_group_with_no_matching_control():
    meta = {"a": {"donor": "O1", "day": 13.0, "ctype": ma.TR, "age": 53.0}}
    assert ma.pair_by_donor_day(meta, ma.TR, ma.NC) == []


# ---------------------------------------------------------- load_betas ---- #
def test_load_betas_handles_comma_separation_and_detection_pval(tmp_path):
    """The real file's quirks: comma-separated, a Detection Pval column after every sample."""
    p = tmp_path / "m.txt.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write("ID_REF,S1,Detection Pval,S2,Detection Pval\n")
        fh.write("cg111,0.10,0.0,0.20,0.0\n")
        fh.write("cg999,0.90,0.0,0.80,0.0\n")      # not a clock CpG -> dropped
    samples, data = ma.load_betas(p, {"cg111"})
    assert samples == ["S1", "S2"]
    assert data["S1"] == {"cg111": 0.10}
    assert data["S2"] == {"cg111": 0.20}
