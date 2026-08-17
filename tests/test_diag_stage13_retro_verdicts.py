"""Unit tests for the Stage 13 retro-verdict pass.

The claim this diagnostic makes is strong — that the old decision rule scored *shuffle controls*
as improvements — so the ways it could be wrong are pinned individually:

* the two defects are DIFFERENT mechanisms and must not be conflated (A1 = the column,
  A3 = the verdict); a test constructs a case where one fires and the other does not;
* the flip count would be inflated ~5x without deduping five value-identical snapshots;
* the correction must not simply flip every verdict.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage13_retro_verdicts.py"
spec = importlib.util.spec_from_file_location("s13", SRC)
s13 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s13)

sc = s13._scorecard()


def _folds(vals: dict, key="level_shift_model"):
    return {d: {key: v} for d, v in vals.items()}


# ---- A1: the column ------------------------------------------------------------------------- #
def test_cancellation_gap_is_infinite_when_the_panel_cancels_exactly():
    """Two donors at +10 and -10 are both badly wrong; |mean| is 0. The old column printed 0."""
    g = s13.cancellation_gap(_folds({"N2": 10.0, "N3": -10.0}), "level_shift_model")
    assert g["abs_of_mean"] == pytest.approx(0.0)
    assert g["mean_of_abs"] == pytest.approx(10.0)
    assert g["understated_by"] == float("inf")


def test_cancellation_gap_is_one_when_every_donor_shares_a_sign():
    g = s13.cancellation_gap(_folds({"N2": 4.0, "N3": 6.0}), "level_shift_model")
    assert g["understated_by"] == pytest.approx(1.0)


def test_cancellation_gap_skips_errored_folds():
    f = _folds({"N2": 4.0, "N3": 6.0})
    f["O1"] = {"_error": "boom", "level_shift_model": 999.0}
    assert s13.cancellation_gap(f, "level_shift_model")["n"] == 2


def test_cancellation_gap_on_real_data_reproduces_the_55x_understatement():
    snaps = s13.load_snapshots()
    g = s13.cancellation_gap(snaps["gc2_A_keep_hff"], "level_shift_ridge")
    assert g["abs_of_mean"] == pytest.approx(0.230, abs=0.001)
    assert g["mean_of_abs"] == pytest.approx(12.723, abs=0.001)
    assert g["understated_by"] == pytest.approx(55.2, abs=0.1)


# ---- A3: the verdict ------------------------------------------------------------------------ #
def test_reverdict_reports_both_rules_on_the_same_folds():
    A = _folds({"N2": -28.0, "N3": -26.0, "O1": -30.0, "O2": -27.0})
    B = _folds({"N2": -22.0, "N3": -20.0, "O1": -24.0, "O2": -21.0})
    r = s13.reverdict(sc, A, B, "level_shift_model")
    assert r["n_folds"] == 4
    assert r["old_diff"] == pytest.approx(+6.0) and r["old_verdict"] == "REGRESSION"
    assert r["new_diff"] == pytest.approx(-6.0) and r["new_verdict"] == "ACCEPT (better)"


def test_reverdict_returns_none_when_the_pair_has_too_few_common_folds():
    assert s13.reverdict(sc, _folds({"N2": 1.0}), _folds({"N2": 2.0}),
                         "level_shift_model") is None


def test_the_two_defects_are_genuinely_different_mechanisms():
    """A1 fires with no comparison at all; A3 fires on a panel where A1's ratio is ~1. If they
    were the same defect this construction would be impossible, and the write-up's separation of
    'the column' from 'the verdict' would be a distinction without a difference."""
    # Every donor shares a sign, so mean(|.|) == |mean| and A1 is INERT (ratio 1.0). A3 still
    # inverts the verdict in both directions, which it could not do if the two were one defect.
    A = _folds({"N2": -20.0, "N3": -22.0, "O1": -21.0, "O2": -19.0})
    B = _folds({"N2": -8.0, "N3": -10.0, "O1": -9.0, "O2": -7.0})
    for panel in (A, B):
        assert s13.cancellation_gap(panel, "level_shift_model")["understated_by"] \
            == pytest.approx(1.0)

    r = s13.reverdict(sc, A, B, "level_shift_model")          # -20.5 -> -8.5: a real improvement
    assert r["old_verdict"] == "REGRESSION", "signed: -20.5 -> -8.5 'increased', so 'worse'"
    assert r["new_verdict"] == "ACCEPT (better)", "magnitude: 20.5 -> 8.5 is a real improvement"

    r2 = s13.reverdict(sc, B, A, "level_shift_model")         # -8.5 -> -20.5: a real worsening
    assert r2["old_verdict"] == "ACCEPT (better)", "signed: -8.5 -> -20.5 'decreased', so 'better'"
    assert r2["new_verdict"] == "REGRESSION", "magnitude: 8.5 -> 20.5 is a real worsening"


# ---- deduping, without which the headline count is ~5x inflated ----------------------------- #
def test_dedupe_groups_snapshots_that_are_identical_on_the_metrics_under_test():
    snaps = {"a": _folds({"N2": 1.0}), "b": _folds({"N2": 1.0}), "c": _folds({"N2": 2.0})}
    groups = s13.dedupe(snaps, ["N2"], ["level_shift_model"])
    assert sorted(len(g) for g in groups) == [1, 2]


def test_the_real_snapshots_collapse_to_five_distinct_ones():
    snaps = s13.load_snapshots()
    groups = s13.dedupe(snaps, sc.DONORS, s13.abs_metrics(sc))
    assert len(snaps) == 9 and len(groups) == 5
    big = max(groups, key=len)
    assert set(big) == {"A_xdonor", "B_fatecal", "B_fatecal_pooled", "baseline",
                        "gc2_A_keep_hff"}


def test_abs_metrics_selects_exactly_the_two_level_shift_rows():
    assert s13.abs_metrics(sc) == ["level_shift_model", "level_shift_ridge"]


# ---- the headline, and the guard against it being an artefact of the fix -------------------- #
def test_the_run_reproduces_the_recorded_counts():
    r = s13.run()
    assert len(r["groups"]) == 5
    assert r["n_comparisons"] == 20
    assert r["n_flips"] == 12
    assert r["shuffle_controls_scored_as_improvements"] == 8


def test_the_correction_does_not_simply_flip_every_verdict():
    """8 of 20 verdicts are UNCHANGED — all of them `noise`. A rule that inverted everything
    would be a different bug, not a fix."""
    r = s13.run()
    unchanged = [c for c in r["comparisons"] if c["old_verdict"] == c["new_verdict"]]
    assert len(unchanged) == 8
    assert {c["new_verdict"] for c in unchanged} == {"noise (CI incl. 0)"}


def test_an_unchanged_verdict_does_not_mean_unchanged_evidence():
    """5 of the 8 unchanged rows had the SIGN of their point estimate flip (e.g. -5.388 ->
    +1.610) and still read `noise` — the same verdict for the opposite reason. Worth pinning so
    'only 12 rows changed' is not read as 'the other 8 were fine'."""
    r = s13.run()
    unchanged = [c for c in r["comparisons"] if c["old_verdict"] == c["new_verdict"]]
    sign_flipped = [c for c in unchanged
                    if np.sign(c["old_diff"]) != np.sign(c["new_diff"])]
    assert len(sign_flipped) == 5


def test_shuffle_controls_are_the_dominant_flip():
    """The direction of the error is the finding: a destroyed-label control cannot IMPROVE the
    per-donor level shift, and under the corrected rule none of them is accepted as better."""
    r = s13.run()
    for c in r["flips"]:
        if "shuffle" in c["b"] and c["old_verdict"] == "ACCEPT (better)":
            assert c["new_verdict"] in ("REGRESSION", "noise (CI incl. 0)")
            assert c["new_diff"] > 0, "magnitude must have INCREASED under the shuffle"


def test_shuffles_really_do_worsen_the_magnitude_in_the_raw_data():
    """Verified from the snapshots directly, not through the corrected code path -- otherwise
    the test would only be proving the fix consistent with itself."""
    snaps = s13.load_snapshots()
    def mag(name):
        v = [f["level_shift_model"] for f in snaps[name].values()
             if isinstance(f, dict) and "_error" not in f]
        return float(np.mean(np.abs(v)))
    assert mag("gc2_A_keep_hff") == pytest.approx(13.120, abs=0.001)
    assert mag("gc2_D_stratshuffle_hff_s0") == pytest.approx(16.245, abs=0.001)
    assert mag("gc2_C_shuffle_hff_s0") == pytest.approx(23.241, abs=0.001)
    assert mag("gc2_A_keep_hff") < mag("gc2_D_stratshuffle_hff_s0") < mag("gc2_C_shuffle_hff_s0")


# ---- read-only contract --------------------------------------------------------------------- #
def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_run_does_not_touch_any_snapshot():
    before = {p.name: p.read_bytes() for p in (ROOT / "scorecard").glob("*.json")}
    s13.run()
    after = {p.name: p.read_bytes() for p in (ROOT / "scorecard").glob("*.json")}
    assert before == after


def test_results_are_json_serialisable():
    json.dumps(s13.run())


# ---- the audit's SCOPE is a closed historical set, not whatever is on disk ------------------- #
def test_the_retro_audit_covers_only_snapshots_the_broken_rule_actually_judged():
    """Originally this globbed `scorecard/*.json`, so its scope grew the moment a new snapshot
    landed -- which happened hours later when the Stage 12 rebuild wrote `c7t_stage12`, and the
    pinned counts failed. That snapshot was produced AFTER Stage 13 shipped, so it was never
    judged by the old rule; counting it would invent comparisons that never happened. The set is
    frozen for that reason, not to keep a test green."""
    assert len(s13.RETRO_SNAPSHOTS) == 9
    assert "c7t_stage12" not in s13.RETRO_SNAPSHOTS
    for name in s13.RETRO_SNAPSHOTS:
        assert (ROOT / "scorecard" / f"{name}.json").exists(), f"{name} is in scope but missing"


def test_snapshots_taken_after_stage_13_exist_but_are_out_of_scope():
    """Guards the claim above: the file really is on disk and really is excluded."""
    on_disk = {p.stem for p in (ROOT / "scorecard").glob("*.json")}
    assert "c7t_stage12" in on_disk, "the Stage 12 rebuild snapshot should exist by now"
    assert set(s13.load_snapshots()) == set(s13.RETRO_SNAPSHOTS)
    assert "c7t_stage12" in set(s13.load_snapshots(names=None)), "names=None still takes all"
