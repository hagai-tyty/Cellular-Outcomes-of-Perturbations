"""Unit tests for STAGE 1.5.2 gate G-c step 1 — pure functions only, no repo data.

G-c's three-way decision is the thing worth pinning: it routes to KEEP, to a retrain, or to
masking 99.7% of the age labels. Every branch is exercised, plus the day decoding, which is the
one place a silent off-by-one would put the whole trajectory on the wrong axis.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


G = _load("diag_gc_hff_signature", "experiments/diag_gc_hff_signature.py")


# --------------------------------------------------------------- day_from_dose_time ---- #
def test_day_decoding_round_trips_the_real_encoder():
    """Against `encode_dose_time` itself, not a reimplementation of it."""
    from cellfate.data.perturbation import encode_dose_time
    days = [0.0, 2.0, 4.0, 8.0, 14.0, 21.0]
    dt = encode_dose_time([0.0] + [1.0] * 5, [d * 24.0 for d in days])
    assert list(G.day_from_dose_time(dt)) == days


def test_the_day_zero_time_floor_decodes_to_zero_not_to_a_fractional_day():
    """Controls are floored to time_h = 0.01 before the log; that must read as day 0."""
    dt = np.array([[-4.0, np.log(G.TIME_FLOOR)]])
    assert G.day_from_dose_time(dt)[0] == 0.0


# ------------------------------------------------------------------ trajectory_stats ---- #
def _traj(day_means, n=50, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    days, ages = [], []
    for d, m in day_means.items():
        days.append(np.full(n, float(d)))
        ages.append(rng.normal(m, noise, n) if noise else np.full(n, float(m)))
    return np.concatenate(days), np.concatenate(ages)


def test_a_perfect_rejuvenation_trajectory_reads_as_one():
    d, y = _traj({0: 0.0, 2: -6.0, 4: -12.0, 6: -18.0, 8: -24.0}, noise=1.0)
    s = G.trajectory_stats(d, y)
    assert s["rho_timepoint"] == pytest.approx(-1.0)
    assert s["slope_yr_per_day"] == pytest.approx(-3.0, abs=0.3)
    assert s["n_descending_steps"] == s["n_steps"] == 4


def test_the_ipsc_endpoint_is_excluded_not_merely_down_weighted():
    """A cell-type change at day 21 must not enter the trend at all."""
    with_ipsc = {0: 0.0, 2: 1.0, 4: 2.0, 6: 3.0, 21: -400.0}
    d, y = _traj(with_ipsc)
    s = G.trajectory_stats(d, y)
    assert s["n_timepoints"] == 4 and 21.0 not in s["days"]
    assert s["rho_timepoint"] > 0          # the real trend is UP; iPSC would have flipped it


def test_a_flat_trajectory_produces_no_signature():
    d, y = _traj({0: 0.0, 2: 0.0, 4: 0.0, 6: 0.0}, noise=0.01, seed=3)
    s = G.trajectory_stats(d, y)
    assert abs(s["slope_yr_per_day"]) < 0.05


def test_too_few_timepoints_cannot_verify():
    d, y = _traj({0: 0.0})
    assert G.trajectory_stats(d, y)["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------- leave_one_timepoint_out ---- #
def test_leave_one_out_exposes_a_trend_carried_by_a_single_point():
    days = [0.0, 2.0, 4.0, 6.0, 8.0]
    means = [0.0, 0.1, -0.1, 0.0, -50.0]        # everything lives in the last point
    loo = G.leave_one_timepoint_out(days, means)
    assert loo["folds"]["drop_day_8"]["slope"] > -1.0        # collapses without it
    assert loo["folds"]["drop_day_0"]["slope"] < -5.0        # survives without any other


def test_leave_one_out_reports_a_narrow_range_for_a_genuinely_robust_trend():
    days = [0.0, 2.0, 4.0, 6.0, 8.0]
    means = [0.0, -6.0, -12.0, -18.0, -24.0]
    loo = G.leave_one_timepoint_out(days, means)
    lo, hi = loo["slope_range"]
    assert hi - lo < 0.01
    assert loo["rho_range"] == [-1.0, -1.0]


# --------------------------------------------------------------------- gc_verdict ---- #
BAND = (-6.45, -1.61)


def _stats(rho, slope):
    return {"status": "OK", "rho_timepoint": rho, "slope_yr_per_day": slope}


def test_both_criteria_hold_keeps_the_labels():
    v = G.gc_verdict(_stats(-0.90, -3.2), -0.50, *BAND)
    assert v["action"] == "KEEP_HFF_LABELS"


def test_neither_criterion_holds_masks_the_labels():
    v = G.gc_verdict(_stats(-0.20, -0.36), -0.50, *BAND)
    assert v["action"] == "MASK_HFF_IN_PHASE_2"


def test_rho_alone_routes_to_the_retrain_not_to_a_convenient_keep():
    """The measured case, 2026-07-31: rho -0.905 passes, slope -1.526 misses by 0.084."""
    v = G.gc_verdict(_stats(-0.905, -1.526), -0.50, *BAND)
    assert v["action"] == "RUN_STEP_2"
    assert v["rho_ok"] is True and v["slope_ok"] is False


def test_slope_alone_also_routes_to_the_retrain():
    v = G.gc_verdict(_stats(-0.10, -3.2), -0.50, *BAND)
    assert v["action"] == "RUN_STEP_2"


def test_a_slope_too_STEEP_also_leaves_the_band():
    """The band is two-sided: 10x methylation's slope is as unlike the signature as 0.1x."""
    assert G.gc_verdict(_stats(-0.95, -30.0), -0.50, *BAND)["action"] == "RUN_STEP_2"


def test_cannot_verify_propagates():
    assert G.gc_verdict({"status": "CANNOT_VERIFY"}, -0.50, *BAND)["action"] == "CANNOT_VERIFY"


# ------------------------------------------------------------------ the recorded run ---- #
def test_the_recorded_result_is_what_the_rules_produce():
    import json
    p = ROOT / "diag_gc_hff_signature_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    r = json.loads(p.read_text(encoding="utf-8"))
    if "hff_trajectory" not in r:
        pytest.skip("pre-registration only")
    pre = r["preregistration"]
    assert G.gc_verdict(r["hff_trajectory"], pre["rho_bar_used"],
                        *pre["slope_band"])["action"] == r["verdict"]["action"]
    # the finding that refutes §0's cited evidence: rho is not -0.214, it is strongly negative
    assert r["hff_trajectory"]["rho_timepoint"] < -0.80
