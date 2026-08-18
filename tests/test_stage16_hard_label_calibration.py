"""Stage 16 -- the fate calibrator must be fitted on the HARD class, not the soft target.

The defect: `train_model.py` fitted Platt on `cal_target[:, SAFE_IDX]`, the SOFT probability
stored in `y_cls`. The calibrator was then excellent at what it was asked to do -- ECE 0.009-0.013
against the soft target on calib -- and wrong for every consumer, all of which read `S` as
P(hard class = safe):

* `res.py`'s safety gate compares it to `tau_safe = 0.85`;
* `scorecard.fate_ece` scores it against `argmax(y_cls)`;
* the served `p_identity_preserved` is read the same way.

Measured cost: ECE against hard labels 0.106-0.113, mean S 0.500 against a hard base rate of
0.540, and 46 of 91 genuinely safe held-out cells carrying a soft target BELOW the 0.76 gate --
so shipped safety sensitivity was 0.297.

These tests pin the fix and, just as importantly, the reasoning: that the alternative repairs
(a second stacked calibrator, or moving `tau_safe`) were rejected on purpose.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cellfate.common.calibration import platt_safe
from cellfate.training.calibrate import fit_platt_binary

ROOT = Path(__file__).resolve().parents[1]
TRAIN_MODEL = ROOT / "src" / "cellfate" / "training" / "train_model.py"


def _ece(p, y, bins=10):
    p, y = np.asarray(p, float), np.asarray(y, float)
    edges, e = np.linspace(0, 1, bins + 1), 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (p >= edges[i]) & (p < hi)
        if m.sum():
            e += (m.sum() / len(p)) * abs(p[m].mean() - y[m].mean())
    return float(e)


# ---- the fix is present, and the defect is gone --------------------------------------------- #
def test_the_calibrator_is_fitted_on_the_hard_class():
    src = TRAIN_MODEL.read_text(encoding="utf-8")
    assert "def _hard_safe(soft: np.ndarray)" in src
    assert "(np.argmax(soft, axis=1) == C.SAFE_IDX).astype(np.float64)" in src
    assert "fate_parts = [(p[:, C.SAFE_IDX], _hard_safe(t)) for p, t in (" in src


def test_the_soft_target_fit_is_gone():
    """The exact expression that caused it."""
    src = TRAIN_MODEL.read_text(encoding="utf-8")
    assert "fate_parts = [(p[:, C.SAFE_IDX], t[:, C.SAFE_IDX]) for p, t in (" not in src


def test_the_xdonor_diagnostic_uses_the_same_target_as_the_shipped_fit():
    """Two calibrators fitted against different targets cannot be compared, and the block exists
    precisely to compare them."""
    src = TRAIN_MODEL.read_text(encoding="utf-8")
    assert "y_xd = _hard_safe(xstats.targets)" in src
    assert "y_xd = xstats.targets[:, C.SAFE_IDX]" not in src


# ---- the mechanism, on constructed data ------------------------------------------------------ #
def test_fitting_on_soft_targets_under_sharpens_relative_to_hard():
    """THE defect, reproduced from scratch: identical probabilities, two targets, two very
    different calibrators -- and the soft-fitted one is systematically flatter."""
    rng = np.random.default_rng(16)
    n = 4000
    hard = (rng.random(n) < 0.54).astype(float)
    # soft labels sit strictly inside (0,1): confident but never certain, as `y_cls` is
    soft = np.where(hard > 0, rng.uniform(0.80, 0.999, n), rng.uniform(0.001, 0.20, n))
    p = np.clip(0.5 + 0.18 * (hard - 0.5) + rng.normal(0, 0.08, n), 0.01, 0.99)

    a_soft, b_soft = fit_platt_binary(p, soft)
    a_hard, b_hard = fit_platt_binary(p, hard)
    assert a_hard > a_soft, "the hard-label fit must sharpen more"
    assert _ece(platt_safe(p, a_hard, b_hard), hard) < _ece(platt_safe(p, a_soft, b_soft), hard)


def test_a_soft_fitted_calibrator_tracks_the_soft_mean_not_the_hard_base_rate():
    """Exactly what was measured on calib: mean S 0.500 against a hard base rate of 0.540."""
    rng = np.random.default_rng(17)
    n = 4000
    hard = (rng.random(n) < 0.60).astype(float)
    soft = np.where(hard > 0, rng.uniform(0.70, 0.95, n), rng.uniform(0.05, 0.30, n))
    p = np.clip(0.5 + 0.2 * (hard - 0.5) + rng.normal(0, 0.05, n), 0.01, 0.99)
    a, b = fit_platt_binary(p, soft)
    cal = platt_safe(p, a, b)
    assert abs(cal.mean() - soft.mean()) < abs(cal.mean() - hard.mean())


def test_the_hard_fit_is_never_worse_than_identity_on_its_own_objective():
    """`fit_platt_binary`'s existing guard must still hold under the new target."""
    rng = np.random.default_rng(18)
    p = rng.uniform(0.05, 0.95, 500)
    y = (rng.random(500) < p).astype(float)
    a, b = fit_platt_binary(p, y)
    assert _ece(platt_safe(p, a, b), y) <= _ece(p, y) + 1e-9


def test_a_single_class_target_still_passes_through_uncalibrated():
    """Hard labels make the single-class case MORE reachable, not less: a fold whose calib is all
    one class now yields a degenerate target rather than a spread of soft values."""
    p = np.linspace(0.05, 0.95, 50)
    assert fit_platt_binary(p, np.ones(50)) == (1.0, 0.0)
    assert fit_platt_binary(p, np.zeros(50)) == (1.0, 0.0)


def test_rank_preservation_survives_the_change():
    """The property that lets this ship without disturbing the fate guards: a > 0 keeps the map
    strictly increasing, so PR-AUC and ROC-AUC are unmoved."""
    rng = np.random.default_rng(19)
    p = rng.uniform(0.05, 0.95, 300)
    y = (rng.random(300) < p).astype(float)
    a, b = fit_platt_binary(p, y)
    assert a > 0
    cal = platt_safe(p, a, b)
    assert np.array_equal(np.argsort(p), np.argsort(cal))


# ---- the rejected alternatives, pinned as decisions ------------------------------------------ #
def test_no_second_calibrator_was_stacked_at_inference():
    """The original 16.8 proposal. Two Platts compose exactly into one, so stacking would have
    worked numerically while hiding that the first was fitted against the wrong target."""
    pred = (ROOT / "src" / "cellfate" / "inference" / "predictor.py").read_text(encoding="utf-8")
    assert pred.count("apply_platt(") == 1, "exactly one calibration step at inference"
    res = (ROOT / "src" / "cellfate" / "inference" / "res.py").read_text(encoding="utf-8")
    assert "platt" not in res.lower(), "the safety gate must not calibrate on the way past"


def test_tau_safe_was_not_moved():
    """Lowering the bar to suit a soft-scale S would be fitting a safety policy to data. The
    empirically-optimal threshold on raw scores was 0.495 against a shipped 0.85 -- and the bar
    stays where it is."""
    cfg = (ROOT / "configs" / "infer" / "default.yaml").read_text(encoding="utf-8")
    assert "tau_safe: 0.85" in cfg
    from cellfate.common.schemas import ResParams
    assert ResParams().tau_safe == 0.85 and ResParams().w == 0.03


def test_the_composition_identity_that_made_the_stacked_fix_look_attractive():
    """sigmoid(a2*logit(sigmoid(a1*logit(p)+b1))+b2) == sigmoid(a1*a2*logit(p) + a2*b1+b2).
    Recorded because it is why the 'add a calibrator' fix would have been indistinguishable from
    'ship different coefficients' -- and therefore why the root cause could have stayed hidden."""
    p = np.linspace(0.02, 0.98, 200)
    a1, b1, a2, b2 = 2.65, 0.42, 1.93, 0.26
    stacked = platt_safe(platt_safe(p, a1, b1), a2, b2)
    composed = platt_safe(p, a1 * a2, a2 * b1 + b2)
    assert np.allclose(stacked, composed, atol=1e-9)


# ---- the recorded evidence -------------------------------------------------------------------- #
def test_the_recorded_stage_16_numbers_are_the_ones_this_fix_targets():
    r = json.loads((ROOT / "results" / "diag_stage16_safety_floor_results.json")
                   .read_text(encoding="utf-8"))
    assert r["pooled"]["raw"]["sensitivity"] == pytest.approx(0.297, abs=0.001)
    assert r["pooled"]["raw"]["false_rejections"] == 64
    assert r["pooled"]["raw"]["n_safe"] == 91


def test_existing_bundles_are_not_silently_assumed_fixed():
    """The change is forward-only: bundles already on disk carry soft-fitted coefficients and
    only a rerun of training re-fits them. Stated in the source so nobody reads a green suite as
    'the shipped folds are now correct'."""
    src = TRAIN_MODEL.read_text(encoding="utf-8")
    assert "STAGE 16 DIAGNOSIS CORRECTED" in src
