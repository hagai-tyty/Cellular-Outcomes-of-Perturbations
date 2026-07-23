"""Every branch of `dump_pool_diag.verdict`.

Written because of the `verify_1a.py` defect: a decision function whose only exercised path was
the one that said PASS. The asymmetry encoded here is the point -- the pool-only calibrator is
scored IN-SAMPLE, so its advantage is inflated and only a NEGATIVE result is conclusive.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dump_pool_diag import read_fold, verdict  # noqa: E402


# ------------------------------------------------------------------ verdict ---- #
def test_no_data_when_either_side_is_empty():
    assert verdict([], [0.1])[0] == "NO DATA"
    assert verdict([0.1], [])[0] == "NO DATA"
    assert verdict([], [])[0] == "NO DATA"


def test_pool_worse_is_conclusive_against_reverting():
    # pool-only loses even WITH its in-sample advantage -> the union fit is exonerated
    v, why = verdict([0.20, 0.22], [0.26, 0.28])
    assert v == "CONCLUSIVE — do not revert"
    assert "in-sample advantage" in why


def test_exact_tie_is_conclusive_not_worth_testing():
    # a tie still means pool-only failed to win despite the bias; must not read as promising
    assert verdict([0.25], [0.25])[0] == "CONCLUSIVE — do not revert"


def test_small_advantage_is_too_weak_to_spend_a_snapshot_on():
    v, why = verdict([0.250], [0.240])          # gain 0.010 < 0.02
    assert v == "WEAK — not worth a snapshot"
    assert "upper bound" in why


def test_just_under_the_threshold_is_weak():
    assert verdict([0.250], [0.2301])[0] == "WEAK — not worth a snapshot"   # gain 0.0199


def test_large_advantage_is_worth_testing_but_flagged_in_sample():
    v, why = verdict([0.250], [0.180])          # gain 0.070
    assert v == "WORTH TESTING — but in-sample"
    assert "UPPER BOUND" in why
    # it must NOT license shipping straight away
    assert "leave-one-donor-out" in why


def test_threshold_boundary_is_worth_testing():
    assert verdict([0.250], [0.230])[0] == "WORTH TESTING — but in-sample"  # gain exactly 0.02


def test_verdict_uses_means_not_the_first_fold():
    # first fold favours pool-only, but the MEAN does not -- must not be fooled by fold 1
    assert verdict([0.30, 0.10], [0.10, 0.32])[0] == "CONCLUSIVE — do not revert"


# ---------------------------------------------------------------- read_fold ---- #
def test_read_fold_reports_missing_file_instead_of_raising(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    import dump_pool_diag
    monkeypatch.setattr(dump_pool_diag, "RUNS", tmp_path / "runs")
    assert "_error" in dump_pool_diag.read_fold("N2")


def test_read_fold_reports_corrupt_json_instead_of_raising(tmp_path, monkeypatch):
    import dump_pool_diag
    b = tmp_path / "runs" / "cellfate_loocv_N2" / "bundle"
    b.mkdir(parents=True)
    (b / "metrics.json").write_text("{not json", encoding="utf-8")
    monkeypatch.setattr(dump_pool_diag, "RUNS", tmp_path / "runs")
    assert "unreadable" in dump_pool_diag.read_fold("N2")["_error"]


def test_read_fold_returns_none_for_absent_keys_not_keyerror(tmp_path, monkeypatch):
    """A pre-Stage-1 bundle has none of these keys; that must degrade, not crash."""
    import dump_pool_diag
    b = tmp_path / "runs" / "cellfate_loocv_N2" / "bundle"
    b.mkdir(parents=True)
    (b / "metrics.json").write_text(json.dumps({"temperature": 1.0}), encoding="utf-8")
    monkeypatch.setattr(dump_pool_diag, "RUNS", tmp_path / "runs")
    r = dump_pool_diag.read_fold("N2")
    assert "_error" not in r
    assert r["shipped_safe_ece_on_pool"] is None
    assert r["xdonor_only_platt_a"] is None
