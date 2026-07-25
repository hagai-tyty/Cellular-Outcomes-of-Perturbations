"""STAGE 1.5 §9 — every branch of the clock-validity diagnostics.

The stakes: this diagnostic decides whether ΔAge is recoverable (application/domain fix) or the
clock is genuinely broken (replace it). A wrong verdict is the expensive mistake. So every verdict
function is driven down every branch, and the load-bearing ones (coverage math, in-range split,
attribution shares, the decision table) are checked against hand-worked numbers.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_clock_validity", _ROOT / "experiments" / "diag_clock_validity.py")
dcv = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = dcv
_SPEC.loader.exec_module(dcv)


# ----------------------------------------------------------- weighted_coverage ---- #
def test_coverage_counts_weight_not_just_genes():
    """The whole point: a clock can keep MOST genes yet lose most WEIGHT, or vice versa."""
    weights = {"A": 10.0, "B": 0.1, "C": 0.1, "D": 0.1}   # A carries ~97% of |weight|
    # data has 3 of 4 genes but MISSES A -> few genes lost, most weight lost
    cov = dcv.weighted_coverage(weights, ["B", "C", "D", "Z"])
    assert cov["n_overlap"] == 3
    assert cov["frac_genes_present"] == pytest.approx(0.75)
    assert cov["frac_abs_weight_present"] == pytest.approx(0.3 / 10.3, abs=1e-6)  # tiny!


def test_coverage_full_when_all_present():
    weights = {"A": 1.0, "B": -2.0, "C": 3.0}
    cov = dcv.weighted_coverage(weights, ["A", "B", "C", "extra"])
    assert cov["frac_abs_weight_present"] == pytest.approx(1.0)
    assert cov["frac_genes_present"] == pytest.approx(1.0)


def test_coverage_reports_concentration():
    # one gene carries 90% of the weight -> k90 == 1
    weights = {"A": 90.0, **{f"g{i}": 1.0 for i in range(10)}}
    cov = dcv.weighted_coverage(weights, list(weights))
    assert cov["weight_k50"] == 1 and cov["weight_k90"] == 1


# ----------------------------------------------------------- coverage_verdict ---- #
def test_coverage_verdict_ok_degraded_crippled():
    assert dcv.coverage_verdict(0.95)["status"] == "OK"
    assert dcv.coverage_verdict(0.80)["status"] == "DEGRADED"
    assert dcv.coverage_verdict(0.50)["status"] == "CRIPPLED"


def test_coverage_verdict_boundaries():
    assert dcv.coverage_verdict(dcv.COVERAGE_OK)["status"] == "OK"
    assert dcv.coverage_verdict(dcv.COVERAGE_CRIPPLED - 1e-9)["status"] == "CRIPPLED"
    assert dcv.coverage_verdict(dcv.COVERAGE_CRIPPLED)["status"] == "DEGRADED"


# --------------------------------------------------------- intercept_dominance ---- #
def test_intercept_dominance_flags_a_dead_clock():
    v = dcv.intercept_dominance([72.3, 72.5, 72.4, 72.6], intercept=72.4)
    assert v["status"] == "DEAD_NEAR_INTERCEPT"


def test_intercept_dominance_passes_a_moving_clock():
    v = dcv.intercept_dominance([30.0, 45.0, 60.0, 80.0], intercept=72.4)
    assert v["status"] == "MOVES"


def test_intercept_dominance_needs_two_points():
    assert dcv.intercept_dominance([72.4], 72.4)["status"] == "CANNOT_VERIFY"


# --------------------------------------------------------- in_range_age_tracking ---- #
def test_in_range_tracks_when_old_reads_older_than_young():
    """The rescue hypothesis H2: exclude the out-of-range neonatal donors, and the clock tracks."""
    pred = {"N2": 98.7, "N3": 36.4, "Y1": 64.9, "Y2": 57.7, "O1": 79.1, "O2": 79.5}
    v = dcv.in_range_age_tracking(pred, dcv.DONOR_AGE)
    assert v["status"] == "TRACKS_IN_RANGE"
    # median split of in-range ages {29,35,53,53} -> young {Y1,Y2}, old {O1,O2}
    assert v["in_range_contrast_years"] == pytest.approx((79.1 + 79.5) / 2 - (64.9 + 57.7) / 2, abs=1e-6)
    assert set(v["out_of_range_donors"]) == {"N2", "N3"}        # age 0, below [1,96]


def test_in_range_no_tracking_when_order_is_wrong():
    pred = {"Y1": 80.0, "Y2": 78.0, "O1": 40.0, "O2": 42.0}    # young read OLDER than old
    v = dcv.in_range_age_tracking(pred, dcv.DONOR_AGE)
    assert v["status"] == "NO_IN_RANGE_TRACKING"


def test_in_range_cannot_verify_with_too_few_in_range():
    v = dcv.in_range_age_tracking({"N2": 50.0, "N3": 51.0}, dcv.DONOR_AGE)
    assert v["status"] == "CANNOT_VERIFY"


def test_in_range_excludes_out_of_range_from_the_contrast():
    """N2's absurd 98.7 (age 0) must NOT enter the young mean — that was M1's error."""
    pred = {"N2": 98.7, "N3": 36.4, "Y1": 64.9, "Y2": 57.7, "O1": 79.1, "O2": 79.5}
    v = dcv.in_range_age_tracking(pred, dcv.DONOR_AGE)
    # young mean is Y1,Y2 only (in range), not N2/N3
    assert v["in_range_contrast_years"] > 15.0


# --------------------------------------------------------- denominator_sensitivity ---- #
def test_denominator_stable_when_predictions_match():
    v = dcv.denominator_sensitivity([50, 60, 70], [50.1, 60.2, 69.9])
    assert v["status"] == "STABLE"


def test_denominator_sensitive_when_they_diverge():
    v = dcv.denominator_sensitivity([50, 60, 70], [70, 80, 95])
    assert v["status"] == "SENSITIVE"


# --------------------------------------------------------- reproduction_verdict ---- #
def test_reproduction_reproduces_degraded_broken_skipped():
    assert dcv.reproduction_verdict(None, 0)["status"] == "SKIPPED"
    assert dcv.reproduction_verdict(10.0, 133)["status"] == "REPRODUCES"      # < 1.5x CV
    assert dcv.reproduction_verdict(25.0, 133)["status"] == "DEGRADED"        # 1.5x–3x
    assert dcv.reproduction_verdict(60.0, 133)["status"] == "BROKEN"          # > 3x


def test_reproduction_boundary_at_1p5x_cv():
    assert dcv.reproduction_verdict(dcv.CLOCK_CV_MAE * 1.5, 100)["status"] == "REPRODUCES"


# --------------------------------------------------------- attribute_direction ---- #
def test_attribution_flags_out_of_domain_when_oskm_and_cellcycle_dominate():
    contrib = {"POU5F1": 3.0, "MKI67": 2.0, "SOX2": 1.0, "COL1A1": 0.5, "FN1": 0.3}
    v = dcv.attribute_direction(contrib)
    assert v["status"] == "OUT_OF_DOMAIN_CONFOUND"
    assert v["confound_share"] > 0.30


def test_attribution_flags_aging_genes_when_they_drive_it():
    contrib = {"CDKN2A": 3.0, "CDKN1A": 2.0, "SERPINE1": 1.0, "COL1A1": 0.4}
    v = dcv.attribute_direction(contrib)
    assert v["status"] == "AGING_GENES_DRIVE_IT"
    assert v["share_senescence_aging"] >= 0.20


def test_attribution_diffuse_when_no_category_dominates():
    contrib = {f"RANDOM{i}": 1.0 for i in range(20)}
    v = dcv.attribute_direction(contrib)
    assert v["status"] == "DIFFUSE"


def test_attribution_only_counts_positive_contributions():
    """The 'age rises' signal is the positive part; negatives (genes pulling younger) don't dilute."""
    contrib = {"POU5F1": 3.0, "MKI67": 2.0, "COL1A1": -10.0}   # big negative must not swamp shares
    v = dcv.attribute_direction(contrib)
    assert v["status"] == "OUT_OF_DOMAIN_CONFOUND"
    assert v["confound_share"] == pytest.approx(1.0)           # 100% of the POSITIVE rise


def test_attribution_cannot_verify_when_empty():
    assert dcv.attribute_direction({})["status"] == "CANNOT_VERIFY"


# --------------------------------------------------------------------- decide ---- #
def _ok_cov(): return {"status": "OK"}
def _crippled_cov(): return {"status": "CRIPPLED"}
def _tracks(): return {"status": "TRACKS_IN_RANGE"}
def _notrack(): return {"status": "NO_IN_RANGE_TRACKING"}
def _repro(): return {"status": "REPRODUCES"}
def _broken(): return {"status": "BROKEN"}
def _skipped(): return {"status": "SKIPPED"}
def _confound(): return {"status": "OUT_OF_DOMAIN_CONFOUND"}
def _diffuse(): return {"status": "DIFFUSE"}


def test_decide_application_fix_takes_priority():
    d = dcv.decide(_crippled_cov(), _tracks(), _repro(), _confound())
    assert d["action"] == "FIX_APPLICATION"
    d2 = dcv.decide(_ok_cov(), _notrack(), _broken(), _diffuse())
    assert d2["action"] == "FIX_APPLICATION"          # BROKEN reproduction also triggers it


def test_decide_target_recoverable_when_in_range_ok_and_reprogramming_confounded():
    d = dcv.decide(_ok_cov(), _tracks(), _repro(), _confound())
    assert d["action"] == "TARGET_RECOVERABLE_DOMAIN_FIX"


def test_decide_investigate_when_in_range_ok_but_attribution_unclear():
    d = dcv.decide(_ok_cov(), _tracks(), _repro(), _diffuse())
    assert d["action"] == "IN_DOMAIN_OK_INVESTIGATE_REPROGRAMMING"


def test_decide_genuine_limitation_when_clean_but_no_tracking():
    d = dcv.decide(_ok_cov(), _notrack(), _repro(), _diffuse())
    assert d["action"] == "GENUINE_CLOCK_LIMITATION"


def test_decide_inconclusive_when_reproduction_skipped_and_no_tracking():
    d = dcv.decide(_ok_cov(), _notrack(), _skipped(), _diffuse())
    assert d["action"] == "INCONCLUSIVE"


def test_decide_crippled_beats_everything_even_with_tracking():
    """An application defect makes the other axes untrustworthy, so it must win."""
    d = dcv.decide(_crippled_cov(), _tracks(), _repro(), _confound())
    assert d["action"] == "FIX_APPLICATION"


# ------------------------------------------------------------------------ bars ---- #
def test_bars_cover_all_four_hypotheses():
    ids = {b["id"] for b in dcv.bars()}
    assert ids == {"H1-coverage", "H1-reproduction", "H2-in-range", "H3-attribution"}
    for b in dcv.bars():
        assert b["bar"] and b["decides"]
