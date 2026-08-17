"""Unit tests for the Stage 14 pre-flight.

Two claims are being made and each could be wrong in a way that changes the plan:

* that a rescaled target buys UNITS ONLY on the linear path (exact equivariance), and
* that it does NOT buy units only on the neural path, because `huber_delta` is fixed in years.

The second is the one that inverted the expected answer, so its arithmetic is pinned directly
rather than only through the real-fold run.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage14_calibration_equivariance.py"
spec = importlib.util.spec_from_file_location("s14", SRC)
s14 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s14)

RESULTS = ROOT / "results" / "diag_stage14_calibration_equivariance_results.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="pre-flight has not been run")


@pytest.fixture(scope="module")
def recorded():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


# ---- the constants must match their sources ------------------------------------------------- #
def test_huber_delta_matches_the_shipped_config():
    """The whole E2 argument rests on this number. If the config changes, the finding changes."""
    cfg = (ROOT / "configs" / "train" / "default.yaml").read_text(encoding="utf-8")
    assert "huber_delta: 2.0" in cfg
    assert s14.HUBER_DELTA == 2.0


def test_k_matches_stage_11s_recorded_lodo_factors():
    r = json.loads((ROOT / "results" / "diag_stage11_scale_results.json").read_text("utf-8"))
    raw = r["variants"]["raw"]
    assert np.mean(list(raw["k_per_donor"].values())) == pytest.approx(s14.K_LS, abs=0.001)
    assert np.mean(list(raw["k_var_per_donor"].values())) == pytest.approx(s14.K_VAR, abs=0.001)


def test_the_least_squares_k_is_the_one_that_shrinks():
    """K_LS < K_VAR is not incidental: least squares wins on MAE by UNDER-reporting magnitude.
    For a REPORTING transform that is the wrong trade, which is why the plan prefers K_VAR."""
    assert s14.K_LS < s14.K_VAR
    r = json.loads((ROOT / "results" / "diag_stage11_scale_results.json").read_text("utf-8"))
    raw = r["variants"]["raw"]
    assert raw["sd_ratio_scaled"] < 0.7 < raw["sd_ratio_var_matched"] < 1.05


# ---- E2: the arithmetic that inverted the expected answer ------------------------------------ #
def test_shrinking_a_target_can_only_move_residuals_INSIDE_a_fixed_knee():
    """Monotone, so the DIRECTION of E2 holds regardless of which estimator's residuals are
    used. Only the magnitude depends on that choice."""
    rng = np.random.default_rng(14)
    resid = rng.normal(0, 3.0, 5000)
    prev = None
    for k in (1.0, 0.8, 0.6, 0.4, 0.2):
        h = s14.huber_region(resid, k)
        if prev is not None:
            assert h["frac_inside_after"] >= prev
        prev = h["frac_inside_after"]


def test_huber_region_reports_before_and_after_correctly():
    resid = np.array([1.0, 3.0, 5.0, 9.0])          # 1 of 4 inside a knee of 2.0
    h = s14.huber_region(resid, k=0.2, delta=2.0)
    assert h["frac_inside_before"] == pytest.approx(0.25)
    assert h["frac_inside_after"] == pytest.approx(1.0)   # 0.2, 0.6, 1.0, 1.8 all < 2
    assert h["median_abs_resid_after"] == pytest.approx(h["median_abs_resid_before"] * 0.2)


def test_a_knee_far_beyond_every_residual_makes_the_loss_pure_l1_both_ways():
    """The case I ASSUMED held before measuring: if residuals are all >> delta, rescaling by k
    changes nothing about the loss regime. It is not what the real folds show."""
    resid = np.full(100, 50.0)
    h = s14.huber_region(resid, k=0.5, delta=2.0)
    assert h["frac_inside_before"] == 0.0 and h["frac_inside_after"] == 0.0


# ---- the units-effect guard table ------------------------------------------------------------ #
def test_the_predicted_effect_scales_exactly_the_three_scale_metrics_and_no_others():
    m = {"dage_mae_model": 10.0, "level_shift_model": -8.0, "conformal_width": 40.0,
         "rank_model_dage": 0.9}
    p = s14.predicted_scorecard_effect(m, 0.5)
    assert p["dage_mae_model"] == 5.0
    assert p["level_shift_model"] == -4.0
    assert p["conformal_width"] == 20.0
    assert p["rank_model_dage"] == 0.9, "ranking is rank-invariant and must be EXACTLY unchanged"


def test_the_guard_table_is_derived_from_a_committed_snapshot_not_a_fresh_run():
    src = SRC.read_text(encoding="utf-8")
    assert 'scorecard" / "c7_A_keep_hff.json"' in src


# ---- E1 on the real folds -------------------------------------------------------------------- #
@has_results
def test_ridge_is_exactly_equivariant_on_every_available_fold(recorded):
    """THE demonstration. A linear model's predictions scale exactly, so for the linear path a
    calibrated target buys precisely nothing beyond units -- measured, not argued."""
    assert recorded["folds"], "no folds were measured"
    for d, f in recorded["folds"].items():
        e = f["equivariance_k_ls"]
        assert e["rel_dev"] < 1e-9, f"{d}: ridge is not equivariant"
        assert e["mae_ratio"] == pytest.approx(e["k"], rel=1e-9), f"{d}: MAE did not scale by k"
        assert e["spearman_scaled"] == pytest.approx(e["spearman_unscaled"], abs=1e-12)


@has_results
def test_both_k_values_are_equivariant_not_just_the_one(recorded):
    for f in recorded["folds"].values():
        for key in ("equivariance_k_ls", "equivariance_k_var"):
            assert f[key]["rel_dev"] < 1e-9


@has_results
def test_the_recorded_folds_show_the_loss_regime_actually_changing(recorded):
    """The finding: residuals are COMPARABLE to the 2.0 yr knee, not far beyond it, so rescaling
    pushes the loss from a Huber mix toward near-pure L2. Every fold must show a large jump, or
    the Stage 14 recommendation to rescale huber_delta alongside the target is unfounded."""
    assert recorded["folds"]
    for d, f in recorded["folds"].items():
        h = f["huber"]
        assert h["frac_inside_before"] < 0.70, f"{d}: residuals were already inside the knee"
        assert h["frac_inside_after"] > 0.80, f"{d}: the regime did not move"
        assert h["frac_inside_after"] - h["frac_inside_before"] > 0.25, f"{d}: jump too small"


@has_results
def test_the_held_out_donor_without_age_labels_is_reported_not_silently_dropped(recorded):
    """N2 has no ΔAge under C-7. It must appear in `errors`, not simply be absent."""
    assert "N2" in recorded["errors"]
    assert "N2" not in recorded["folds"]


# ---- contract -------------------------------------------------------------------------------- #
def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_the_script_trains_nothing():
    """It refits a cheap ridge and reads built folds. If it ever imports the trainer, the
    'read-only pre-flight' contract is broken."""
    src = SRC.read_text(encoding="utf-8")
    assert "train_model" not in src and "training.train" not in src
