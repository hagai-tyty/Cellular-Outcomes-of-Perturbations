"""Unit tests for the Stage 12 §12.9 rebuild verdict.

The point of this helper is that the scorecard's own `conformal_coverage` row judges the WRONG
quantity for this Change: it is registered `("higher", ...)`, but coverage is target-seeking --
1.000 is not better than 0.900, it means the intervals are too wide. §12.9's target is
|coverage - nominal|, and these tests pin that distinction, because getting it backwards would
turn an over-covering fold into evidence of success.

The verdict branches are tested on synthetic input so all three are reachable regardless of what
the real rebuild produces.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage12_rebuild_verdict.py"
spec = importlib.util.spec_from_file_location("s12v", SRC)
s12v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s12v)

RESULTS = ROOT / "results" / "diag_stage12_rebuild_verdict_results.json"
def _find_plan(name: str) -> Path:
    """Locate a plan by FILENAME, wherever it currently sits under plans/.

    Plans were reorganised into `(older)base plans/` and `(newer)practical plans/` on 2026-08-20,
    which broke every test that read one by a fixed path. Searching by name survives the next
    reshuffle too.
    """
    hits = sorted(ROOT.joinpath("plans").rglob(name))
    if not hits:
        raise AssertionError(f"plan {name!r} not found anywhere under plans/")
    return hits[0]


has_results = pytest.mark.skipif(not RESULTS.exists(), reason="rebuild verdict not yet computed")


def _folds(cov: dict):
    return {d: {"conformal_coverage": c, "conformal_level": 0.9} for d, c in cov.items()}


# ---- the quantity: distance from nominal, not raw coverage ---------------------------------- #
def test_over_covering_and_under_covering_by_the_same_amount_score_the_same():
    """THE distinction. 1.00 and 0.80 are both 0.10 from nominal. A rule preferring 1.00 is
    measuring interval width, not calibration."""
    g = s12v.coverage_gap(_folds({"A": 1.00, "B": 0.80}), ["A", "B"])
    assert g["A"] == pytest.approx(0.10)
    assert g["B"] == pytest.approx(0.10)


def test_perfect_coverage_scores_zero():
    assert s12v.coverage_gap(_folds({"A": 0.90}), ["A"])["A"] == pytest.approx(0.0)


def test_the_scorecard_has_since_been_fixed_and_now_agrees_with_this_helper():
    """This test used to assert `METRICS["conformal_coverage"][0] == "higher"`, with a note that
    if the scorecard were ever fixed the test should be REVISITED rather than left silently
    passing. Stage 17 fixed it, and the note fired. Revisited here.

    The helper is no longer the only thing computing distance-to-target, so the useful assertion
    is now the stronger one: two independently written implementations of the same statistic must
    agree, on real data. If they ever diverge, one of them is wrong."""
    sc = s12v._sc()
    assert sc.METRICS["conformal_coverage"][0] == "target"

    A = json.loads((ROOT / "scorecard" / "c7_A_keep_hff.json").read_text("utf-8"))["folds"]
    B = json.loads((ROOT / "scorecard" / "c7t_stage12.json").read_text("utf-8"))["folds"]

    mine, _, n_mine, common = s12v.paired(s12v.coverage_gap(A, sc.DONORS),
                                          s12v.coverage_gap(B, sc.DONORS))
    theirs, _, n_theirs = sc._paired(A, B, "conformal_coverage", target=True)
    assert n_mine == n_theirs == 5
    assert mine == pytest.approx(theirs, abs=1e-12)

    # and the per-fold distances agree one by one, not merely in the mean
    for d in common:
        assert s12v.coverage_gap(A, [d])[d] == pytest.approx(
            sc._judged(A[d], "conformal_coverage", target=True), abs=1e-12)


def test_errored_and_missing_folds_are_skipped():
    f = _folds({"A": 0.9})
    f["B"] = {"_error": "boom"}
    f["C"] = {"conformal_coverage": None}
    assert set(s12v.coverage_gap(f, ["A", "B", "C", "D"])) == {"A"}


# ---- pairing ---------------------------------------------------------------------------------- #
def test_pairing_uses_only_folds_present_in_both():
    md, _, n, common = s12v.paired({"A": 0.1, "B": 0.2, "C": 0.3}, {"A": 0.0, "B": 0.1})
    assert common == ["A", "B"] and n == 2
    assert md == pytest.approx(-0.1)


def test_pairing_needs_two_folds():
    md, (lo, hi), n, _ = s12v.paired({"A": 0.1}, {"A": 0.0})
    assert md is None and lo is None and hi is None and n == 1


# ---- the three pre-registered branches -------------------------------------------------------- #
def test_a_consistent_move_toward_nominal_reads_as_a_real_improvement():
    a = {d: 0.20 for d in "ABCDEF"}
    b = {d: 0.05 for d in "ABCDEF"}
    md, (lo, hi), _, _ = s12v.paired(a, b)
    assert s12v.verdict_from(md, lo, hi).startswith("COVERAGE MOVED TOWARD NOMINAL")


def test_a_consistent_move_away_from_nominal_demands_investigation():
    a = {d: 0.05 for d in "ABCDEF"}
    b = {d: 0.20 for d in "ABCDEF"}
    md, (lo, hi), _, _ = s12v.paired(a, b)
    v = s12v.verdict_from(md, lo, hi)
    assert v.startswith("COVERAGE MOVED AWAY FROM NOMINAL")
    assert "INVESTIGATE" in v


def test_a_noisy_wash_reads_as_a_publishable_negative_not_a_failure():
    a = {"A": 0.10, "B": 0.20, "C": 0.05, "D": 0.30, "E": 0.02, "F": 0.25}
    b = {"A": 0.30, "B": 0.02, "C": 0.28, "D": 0.04, "E": 0.26, "F": 0.03}
    md, (lo, hi), _, _ = s12v.paired(a, b)
    v = s12v.verdict_from(md, lo, hi)
    assert v.startswith("NO DETECTABLE MOVE")
    assert "publishable negative" in v


def test_an_undetermined_comparison_says_so():
    assert s12v.verdict_from(None, None, None).startswith("UNDETERMINED")


# ---- the guards ------------------------------------------------------------------------------- #
def test_the_registered_guards_are_the_ones_the_plan_names():
    plan = _find_plan("STAGE_12_CELL_ID_UNIQUENESS.md").read_text(encoding="utf-8")
    for g in s12v.GUARDS:
        assert g in plan, f"{g} is checked but not named in the plan"
    assert set(s12v.GUARDS) == {"fate_prauc", "fate_roc", "rank_model_dage"}


def test_the_comparison_targets_the_right_two_snapshots():
    assert s12v.BASELINE_TAG == "c7_A_keep_hff"
    assert s12v.TREATMENT_TAG == "c7t_stage12"
    assert s12v.TREATMENT_TAG != s12v.BASELINE_TAG, "the baseline must not be the treatment"


def test_the_nominal_level_matches_what_the_snapshots_record():
    """0.90 is not a constant to guess -- every fold of the committed baseline records it."""
    snap = json.loads((ROOT / "scorecard" / "c7_A_keep_hff.json").read_text(encoding="utf-8"))
    levels = {f["conformal_level"] for f in snap["folds"].values()
              if isinstance(f, dict) and "_error" not in f}
    assert levels == {s12v.NOMINAL}


def test_a_missing_snapshot_aborts_rather_than_half_reporting():
    r = s12v.run(baseline="c7_A_keep_hff", treatment="does_not_exist")
    assert r["ABORTED"] is True and "missing snapshot" in r["reason"]


# ---- the recorded run -------------------------------------------------------------------------- #
@has_results
def test_the_recorded_verdict_is_the_mechanical_result_of_the_rule(recorded=None):
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert r["verdict"] == s12v.verdict_from(r["mean_diff"], r["ci"][0], r["ci"][1])


@has_results
def test_the_recorded_run_paired_at_least_five_folds():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert r["n_folds"] >= 5, "N2 has no ΔAge under C-7, so 5 is the expected maximum"


# ---- contract ---------------------------------------------------------------------------------- #
def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1
