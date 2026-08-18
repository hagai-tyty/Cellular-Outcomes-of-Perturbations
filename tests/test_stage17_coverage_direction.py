"""Stage 17 -- coverage is target-seeking, and fold agreement must be visible.

Two defects, both found while executing Stage 12 and deliberately left unpatched until that
measurement was finished (changing the instrument mid-measurement invalidates the measurement):

D1  `conformal_coverage` was registered ("higher", ...). Coverage 1.000 is not better than 0.900 --
    it means the intervals are too wide, which is why `conformal_width` sits at 63-81 years. Under
    "higher is better", widening every interval until nothing escaped would have scored ACCEPT.
    4-5 of 6 folds in most committed snapshots are OVER-covering, so this was live, not theoretical.

D2  The table never showed whether the folds AGREED. The scorecard header tells the reader to
    "check the per-fold column", but only `dage_mae_model` ever printed one. `conformal_width` in
    the Stage 12 comparison moved the same way on ALL FIVE folds and still read `noise`.

The most important tests here are the ones that keep D2's fix honest: the tally must never become
a second significance test.
"""
from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import scorecard as sc  # noqa: E402

BASE = ROOT / "scorecard" / "c7_A_keep_hff.json"
TREAT = ROOT / "scorecard" / "c7t_stage12.json"


def _cov(vals: dict, level=0.9):
    return {d: {"conformal_coverage": v, "conformal_level": level} for d, v in vals.items()}


def _real():
    return (json.loads(BASE.read_text(encoding="utf-8"))["folds"],
            json.loads(TREAT.read_text(encoding="utf-8"))["folds"])


# ---- D1: coverage is target-seeking --------------------------------------------------------- #
def test_over_and_under_covering_by_the_same_amount_score_the_same():
    """THE defect. Under ("higher", ...) 1.000 outranked 0.800; both are 0.100 from nominal."""
    f = _cov({"N2": 1.00, "N3": 0.80})
    assert sc._judged(f["N2"], "conformal_coverage", target=True) == pytest.approx(0.10)
    assert sc._judged(f["N3"], "conformal_coverage", target=True) == pytest.approx(0.10)


def test_perfect_coverage_scores_zero_distance():
    assert sc._judged(_cov({"N2": 0.90})["N2"], "conformal_coverage",
                      target=True) == pytest.approx(0.0)


def test_the_target_is_read_per_fold_not_hard_coded():
    """0.90 is data. A fold carrying a different conformal_level must be judged against its own."""
    assert sc._judged({"conformal_coverage": 0.95, "conformal_level": 0.95},
                      "conformal_coverage", target=True) == pytest.approx(0.0)
    assert sc._judged({"conformal_coverage": 0.95, "conformal_level": 0.80},
                      "conformal_coverage", target=True) == pytest.approx(0.15)


def test_conformal_coverage_is_registered_as_target():
    assert sc.METRICS["conformal_coverage"][0] == "target"
    assert sc.TARGET_OF["conformal_coverage"] == "conformal_level"


def test_widening_every_interval_to_full_coverage_no_longer_reads_accept():
    """The exact gaming the old direction permitted: push coverage to 1.0 everywhere and the old
    rule calls it a win. It must now read REGRESSION, because it moves AWAY from nominal."""
    A = _cov({d: 0.90 for d in sc.DONORS})
    B = _cov({d: 1.00 for d in sc.DONORS})
    md_o, (lo_o, hi_o), _ = sc._paired(A, B, "conformal_coverage")
    assert sc._verdict("higher", md_o, lo_o, hi_o) == "ACCEPT (better)", "the defect, reproduced"
    md, (lo, hi), _ = sc._paired(A, B, "conformal_coverage", target=True)
    assert sc._verdict("target", md, lo, hi) == "REGRESSION"


@pytest.mark.parametrize("start", [0.60, 1.00])
def test_moving_toward_nominal_from_either_side_reads_as_better(start):
    A = _cov({d: start for d in sc.DONORS})
    B = _cov({d: 0.90 for d in sc.DONORS})
    md, (lo, hi), _ = sc._paired(A, B, "conformal_coverage", target=True)
    assert sc._verdict("target", md, lo, hi) == "ACCEPT (better)"


def test_a_fold_missing_its_target_drops_out_instead_of_crashing():
    """An old snapshot without `conformal_level` must not take the whole row down."""
    assert sc._fold_ok({"conformal_coverage": 0.9}, "conformal_coverage") is False
    assert sc._fold_ok({"conformal_coverage": 0.9, "conformal_level": 0.9},
                       "conformal_coverage") is True


def test_the_real_snapshots_are_mostly_over_covering():
    """Why D1 mattered in practice rather than in principle."""
    folds = [f for f in json.loads(BASE.read_text(encoding="utf-8"))["folds"].values()
             if isinstance(f, dict) and "_error" not in f and f.get("conformal_coverage")]
    over = sum(1 for f in folds if f["conformal_coverage"] > f["conformal_level"])
    assert over >= 4 and over / len(folds) > 0.7


# ---- D2: the fold-direction tally ----------------------------------------------------------- #
def test_a_unanimous_change_is_flagged_even_when_its_ci_spans_zero():
    """THE reason this exists. `conformal_width` in the Stage 12 comparison: all five folds moved
    the same way, mean -6.97 yr, and the CI still included 0 because the MAGNITUDES vary."""
    A, B = _real()
    t = sc.fold_tally(A, B, "conformal_width", "lower")
    assert (t["better"], t["worse"]) == (5, 0)
    assert t["unanimous"] is True
    md, (lo, hi), _ = sc._paired(A, B, "conformal_width")
    assert lo < 0 < hi, "the CI really does span zero"
    assert sc._verdict("lower", md, lo, hi) == "noise (CI incl. 0)", "and the VERDICT is unchanged"


def test_the_tally_is_symmetric_and_flags_unanimous_regressions_too():
    A = _cov({d: 0.90 for d in sc.DONORS})
    B = _cov({d: 0.60 for d in sc.DONORS})
    t = sc.fold_tally(A, B, "conformal_coverage", "target", target=True)
    assert (t["better"], t["worse"], t["unanimous"]) == (0, 6, True)


def test_ties_count_as_same_and_do_not_make_a_row_unanimous():
    A = _cov({d: 0.90 for d in sc.DONORS})
    t = sc.fold_tally(A, A, "conformal_coverage", "target", target=True)
    assert t["same"] == 6 and t["better"] == 0 and t["worse"] == 0
    assert t["unanimous"] is False, "a row of ties agrees about nothing"


def test_unanimity_needs_at_least_four_folds():
    """Three folds agreeing is not surprising enough to mark."""
    A = _cov({"N2": 0.90, "N3": 0.90, "O1": 0.90})
    B = _cov({"N2": 0.95, "N3": 0.96, "O1": 0.97})
    t = sc.fold_tally(A, B, "conformal_coverage", "target", target=True)
    assert (t["better"], t["worse"]) == (0, 3)
    assert t["unanimous"] is False


def test_neutral_metrics_get_no_tally():
    A = {d: {"ood_rate": 0.3} for d in sc.DONORS}
    B = {d: {"ood_rate": 0.9} for d in sc.DONORS}
    assert sc.fold_tally(A, B, "ood_rate", "neutral")["n"] == 0


def test_the_tally_respects_direction_for_higher_is_better_metrics():
    A = {d: {"fate_prauc": 0.5} for d in sc.DONORS}
    B = {d: {"fate_prauc": 0.9} for d in sc.DONORS}
    assert sc.fold_tally(A, B, "fate_prauc", "higher")["better"] == 6


# ---- the guard: the tally must never become a second significance test ---------------------- #
def test_no_p_value_is_ever_computed_or_printed_for_the_tally():
    """With 5 paired folds the smallest achievable two-sided sign-test p is 0.0625, so a unanimous
    result can never clear 0.05. Printing it would invite exactly the misreading this prevents."""
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "0.0625" in src, "the limitation must be stated in the source"
    assert "p_value" not in src
    assert "binom" not in src
    assert "from math import comb" not in src


def test_the_tally_never_reaches_the_verdict_function():
    """`_verdict` takes only the CI. If the tally could reach it, the decision rule would have
    changed -- and the plan says it does not."""
    assert list(inspect.signature(sc._verdict).parameters) == ["direction", "md", "lo", "hi"]
    v = inspect.getsource(sc._verdict)
    for token in ("tally", "unanimous", "fold_tally"):
        assert token not in v


def test_the_decision_rule_text_is_unchanged():
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "a change is REAL only if the paired 95% CI across folds excludes 0" in src


# ---- what must NOT have changed -------------------------------------------------------------- #
def test_non_target_metrics_are_judged_exactly_as_before():
    A, B = _real()
    for key, (direction, _) in sc.METRICS.items():
        if direction == "target":
            continue
        mag = direction == "abs"
        assert sc._paired(A, B, key, magnitude=mag) == sc._paired(A, B, key, magnitude=mag,
                                                                 target=False)


def test_the_stage_13_invariant_still_holds_on_every_metric():
    """col_B - col_A == mean diff, now including the target metric. It can only hold if the
    columns and the paired statistic take the distance at the same point in the computation."""
    A, B = _real()
    checked = 0
    for key, (direction, _) in sc.METRICS.items():
        mag, tgt = direction == "abs", direction == "target"
        common = sc._common_folds(A, B, key)
        va = sc._agg(A, key, only=common, magnitude=mag, target=tgt)
        vb = sc._agg(B, key, only=common, magnitude=mag, target=tgt)
        md, _, _ = sc._paired(A, B, key, magnitude=mag, target=tgt)
        if va is None or vb is None or md is None:
            continue
        assert vb - va == pytest.approx(md, abs=1e-9), key
        checked += 1
    assert checked == 18


def test_measure_fold_still_stores_raw_coverage_not_a_distance():
    """The stored value must stay raw, or every snapshot becomes unreadable in its own terms and
    the over/under-coverage direction is destroyed permanently."""
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert 'out["conformal_coverage"] = float(((y >= lo[m]) & (y <= hi[m])).mean())' in src


def test_conformal_width_was_deliberately_left_alone():
    """Narrower IS better at equal coverage, and coverage is judged separately. A decision, not
    an oversight."""
    assert sc.METRICS["conformal_width"][0] == "lower"
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "narrower IS genuinely better at equal coverage" in src
