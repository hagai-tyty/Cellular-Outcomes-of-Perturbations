"""Stage 16 -- the recalibration, VERIFIED on real artefacts rather than asserted.

The `train_model.py` fix is forward-only: Platt is fitted during training, so every bundle already
on disk kept its soft-fitted coefficients. `local_runners/recalibrate_folds.py` re-runs only that
step against the hard target, producing `_s16` from `_s12` without retraining.

These tests pin what the artefacts actually show. They are the difference between "implemented and
tested" and "empirically validated", and they exist because the first is not evidence for the
second.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RECAL = ROOT / "results" / "stage16_recalibration_results.json"
EV_S12 = ROOT / "results" / "diag_stage16_safety_floor_results_s12.json"
EV_S16 = ROOT / "results" / "diag_stage16_safety_floor_results_s16.json"
SNAP_A = ROOT / "scorecard" / "c7t_stage12.json"
SNAP_B = ROOT / "scorecard" / "c7t_stage16.json"

needs = pytest.mark.skipif(
    not all(p.exists() for p in (RECAL, EV_S12, EV_S16, SNAP_A, SNAP_B)),
    reason="recalibration/verification artefacts not present")

# The fold BUILDS are gitignored (~260 MB each), so they exist on the data machine and never on
# CI. The results/snapshot JSONs above ARE committed, so `needs` alone is true on CI and the two
# tests that open `cellfate_loocv_*/bundle/` would raise FileNotFoundError there. They get their
# own guard; every other test in this file reads committed artefacts and still runs in CI, which
# is where most of the value is.
FOLDS = [ROOT / f"cellfate_loocv_{d}{sfx}" / "bundle" / "temperature.json"
         for d in ("N2", "N3", "O1", "O2", "Y1", "Y2") for sfx in ("_s12", "_s16")]
needs_folds = pytest.mark.skipif(
    not all(p.exists() for p in FOLDS),
    reason="fold builds are gitignored; present only on the data machine")


def _pooled(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))["pooled"]["raw"]


# ---- the recalibration itself ---------------------------------------------------------------- #
@needs
def test_every_fold_was_recalibrated_and_the_slope_roughly_doubled():
    """The soft-fitted calibrator under-sharpens; fitting on the hard class roughly doubles the
    slope, which is what the composed-coefficient prediction said it would do."""
    rows = json.loads(RECAL.read_text(encoding="utf-8"))
    assert {r["donor"] for r in rows} == {"N2", "N3", "O1", "O2", "Y1", "Y2"}
    for r in rows:
        assert 2.5 < r["a_old"] < 2.7, r["donor"]
        assert 4.8 < r["a_new"] < 5.1, r["donor"]
        assert r["a_new"] > 1.7 * r["a_old"]


@needs
@needs_folds
def test_the_original_folds_were_not_mutated():
    """`_s12` and the `c7t_stage12` snapshot taken from it must stay valid: the recalibration
    writes a NEW fold set and only replaces bundle/temperature.json inside it."""
    rows = json.loads(RECAL.read_text(encoding="utf-8"))
    for r in rows:
        old = json.loads((ROOT / f"cellfate_loocv_{r['donor']}_s12" / "bundle" /
                          "temperature.json").read_text(encoding="utf-8"))
        new = json.loads((ROOT / f"cellfate_loocv_{r['donor']}_s16" / "bundle" /
                          "temperature.json").read_text(encoding="utf-8"))
        assert old["platt_a"] == pytest.approx(r["a_old"])
        assert new["platt_a"] == pytest.approx(r["a_new"])
        assert old["platt_a"] != new["platt_a"]


@needs
@needs_folds
def test_each_recalibrated_bundle_records_its_own_provenance():
    for d in ("N2", "O1", "Y2"):
        j = json.loads((ROOT / f"cellfate_loocv_{d}_s16" / "bundle" /
                        "recalibration.json").read_text(encoding="utf-8"))
        assert j["stage"] == 16 and j["target"] == "hard"
        assert j["source_fold"] == f"cellfate_loocv_{d}_s12"


# ---- the held-out safety evaluation ---------------------------------------------------------- #
@needs
def test_sensitivity_more_than_doubles_on_the_real_artefacts():
    a, b = _pooled(EV_S12), _pooled(EV_S16)
    assert a["n"] == b["n"] == 119 and a["n_safe"] == b["n_safe"] == 91
    assert a["sensitivity"] == pytest.approx(0.275, abs=0.005)
    assert b["sensitivity"] == pytest.approx(0.670, abs=0.005)
    assert b["sensitivity"] > 2 * a["sensitivity"]


@needs
def test_the_predicted_specificity_loss_did_not_occur():
    """THE thing that had to be checked. The decision to ship was taken on an expected
    specificity drop 0.929 -> 0.821 and 2 -> 5 false approvals. Neither happened: on the real
    artefacts specificity is UNCHANGED and false approvals are unchanged at 2. The safety posture
    did not loosen."""
    a, b = _pooled(EV_S12), _pooled(EV_S16)
    assert a["specificity"] == pytest.approx(0.929, abs=0.005)
    assert b["specificity"] == pytest.approx(a["specificity"], abs=1e-9)
    assert a["false_approvals"] == b["false_approvals"] == 2


@needs
def test_false_rejections_more_than_halve():
    a, b = _pooled(EV_S12), _pooled(EV_S16)
    assert a["false_rejections"] == 66 and b["false_rejections"] == 30
    assert (a["false_rejections"] - b["false_rejections"]) / a["false_rejections"] >= 0.50


@needs
def test_balanced_accuracy_improves_and_no_fold_gets_worse():
    a = json.loads(EV_S12.read_text(encoding="utf-8"))
    b = json.loads(EV_S16.read_text(encoding="utf-8"))
    assert b["pooled"]["raw"]["balanced_accuracy"] > a["pooled"]["raw"]["balanced_accuracy"]
    for d in a["folds"]:
        assert (b["folds"][d]["raw"]["false_rejections"]
                <= a["folds"][d]["raw"]["false_rejections"]), d


@needs
def test_the_safe_class_now_sits_above_the_bar():
    a = json.loads(EV_S12.read_text(encoding="utf-8"))["pooled"]
    b = json.loads(EV_S16.read_text(encoding="utf-8"))["pooled"]
    assert a["median_S_true_safe"] < a["threshold"] < b["median_S_true_safe"]
    assert b["frac_true_safe_below_threshold"] < 0.4 < a["frac_true_safe_below_threshold"]


# ---- the pre-registered guards, on the scorecard ---------------------------------------------- #
@needs
def test_ranking_metrics_are_bit_identical():
    """Platt is monotone, so PR-AUC and ROC-AUC must not move AT ALL. Any movement would mean the
    implementation reordered cells, not that it helped."""
    A = json.loads(SNAP_A.read_text(encoding="utf-8"))["folds"]
    B = json.loads(SNAP_B.read_text(encoding="utf-8"))["folds"]
    for d, fa in A.items():
        fb = B[d]
        if "_error" in fa or "_error" in fb:
            continue
        for k in ("fate_prauc", "fate_roc"):
            if fa.get(k) is not None:
                assert fb[k] == pytest.approx(fa[k], abs=1e-12), f"{d}/{k}"


@needs
def test_every_delta_age_metric_is_untouched():
    """Fate calibration must not move a single ΔAge number."""
    A = json.loads(SNAP_A.read_text(encoding="utf-8"))["folds"]
    B = json.loads(SNAP_B.read_text(encoding="utf-8"))["folds"]
    for d, fa in A.items():
        fb = B[d]
        if "_error" in fa or "_error" in fb:
            continue
        for k in ("dage_mae_model", "dage_mae_ridge", "level_shift_model", "rank_model_dage",
                  "conformal_coverage", "conformal_width", "ood_rate"):
            if fa.get(k) is not None:
                assert fb[k] == pytest.approx(fa[k], abs=1e-12), f"{d}/{k}"


@needs
def test_res_is_still_exactly_zero():
    """Stage 15's `R_eff` gate is independent of the safety gate, and was stated in advance to be
    unmoved by this. Verified rather than assumed."""
    A = json.loads(SNAP_A.read_text(encoding="utf-8"))["folds"]
    B = json.loads(SNAP_B.read_text(encoding="utf-8"))["folds"]
    for d, f in B.items():
        if "_error" in f:
            continue
        # NOT `== 0.0`: Y1 carries floating-point RESIDUE (8.1e-11 before, 1.2e-12 after) rather
        # than a real score. That dust is the same artefact behind the retracted "Spearman 0.40
        # over RES" headline, so it is asserted as dust explicitly instead of being rounded away.
        assert f["res_max"] < 1e-9, d
        assert f["res_median"] == 0.0, d
        assert f["res_approvals"] == 0, d
        assert A[d]["res_approvals"] == 0, d
    assert B["Y1"]["res_max"] > 0.0, "Y1's residue is real dust, not an exact zero"


@needs
def test_fate_ece_is_the_only_metric_that_moved_and_it_improved():
    """The whole change, in one line: the calibration metric gets better and nothing else moves."""
    A = json.loads(SNAP_A.read_text(encoding="utf-8"))["folds"]
    B = json.loads(SNAP_B.read_text(encoding="utf-8"))["folds"]
    moved = set()
    for d, fa in A.items():
        fb = B[d]
        if "_error" in fa or "_error" in fb:
            continue
        for k, va in fa.items():
            if k.startswith("_") or not isinstance(va, (int, float)):
                continue
            if fb.get(k) is not None and abs(fb[k] - va) > 1e-9:
                moved.add(k)
    # Both CALIBRATION metrics move and nothing else does. `fate_ece_platt` moves per fold
    # because its own calibrator is re-fitted on a changed `S`; its AGGREGATE is unchanged
    # (0.180 -> 0.180, CI includes 0), which is the expected signature of a calibrator that had
    # already been compensating for the defect this fix removes.
    assert moved == {"fate_ece", "fate_ece_platt"},         f"unexpected movement: {moved - {'fate_ece', 'fate_ece_platt'}}"
    ece_a = [f["fate_ece"] for f in A.values() if isinstance(f, dict) and f.get("fate_ece")]
    ece_b = [f["fate_ece"] for f in B.values() if isinstance(f, dict) and f.get("fate_ece")]
    assert sum(ece_b) / len(ece_b) < sum(ece_a) / len(ece_a)
