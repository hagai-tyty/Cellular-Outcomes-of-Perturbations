"""Unit tests for P3 — pure functions only, no bundle data.

The metric is the whole argument of this experiment. P3 rejects AUC because inside a held-out
timepoint the `day` arm predicts a CONSTANT, which AUC scores at 0.5 by construction and log-loss
scores honestly. These tests pin that property, so the reasoning cannot be quietly undone by a
later metric swap.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


P = _load("p3_progress", "experiments/p3_progress.py")


# ------------------------------------------------------------------------- log_loss ---- #
def test_log_loss_matches_the_definition():
    assert P.log_loss([1, 0], [0.75, 0.25]) == pytest.approx(-np.log(0.75))


def test_a_constant_predictor_naming_the_true_rate_is_scored_well():
    """THE property the metric was chosen for. A held-out timepoint with 30% risk, predicted as a
    flat 0.30 by the day arm, must score at the entropy of the truth -- not be punished for being
    constant, and not be rewarded either."""
    y = np.array([1] * 30 + [0] * 70)
    ideal = P.log_loss(y, np.full(100, 0.30))
    entropy = -(0.3 * np.log(0.3) + 0.7 * np.log(0.7))
    assert ideal == pytest.approx(entropy, abs=1e-9)


def test_a_constant_predictor_beats_a_wrong_per_cell_predictor():
    """Per-cell variation only wins if it is RIGHT. Random per-cell noise must lose to the
    correct constant -- otherwise the metric would reward progress for merely varying."""
    rng = np.random.default_rng(0)
    y = (rng.random(2000) < 0.3).astype(int)
    const = P.log_loss(y, np.full(2000, 0.30))
    noisy = P.log_loss(y, rng.uniform(0.05, 0.95, 2000))
    assert const < noisy


def test_log_loss_is_clipped_so_a_confident_miss_is_finite():
    assert np.isfinite(P.log_loss([1], [0.0]))
    assert np.isfinite(P.log_loss([0], [1.0]))


def test_brier_matches_the_definition():
    assert P.brier([1, 0], [0.75, 0.25]) == pytest.approx((0.0625 + 0.0625) / 2)


# ---------------------------------------------------------------------- fit_predict ---- #
def test_single_class_training_target_falls_back_to_the_base_rate():
    """No decision boundary exists; the honest prediction is the constant base rate, which
    log-loss then scores fairly. It must not raise and must not invent a boundary."""
    out = P.fit_predict(np.arange(10).reshape(-1, 1), np.zeros(10, int), np.zeros((3, 1)))
    assert out.shape == (3,)
    assert np.allclose(out, 0.0)


def test_a_separable_training_set_predicts_in_the_right_direction():
    Xtr = np.array([[0.0], [0.0], [1.0], [1.0]])
    ytr = np.array([0, 0, 1, 1])
    out = P.fit_predict(Xtr, ytr, np.array([[0.0], [1.0]]))
    assert out[0] < out[1]


def test_fit_predict_accepts_multi_column_predictors():
    rng = np.random.default_rng(1)
    Xtr = rng.normal(size=(40, 2))
    ytr = (Xtr[:, 0] > 0).astype(int)
    assert P.fit_predict(Xtr, ytr, rng.normal(size=(5, 2))).shape == (5,)


# ------------------------------------------------------------------------ paired_ci ---- #
def test_paired_ci_uses_the_t_value_for_nine_folds():
    m, lo, hi, n = P.paired_ci([1.0] * 9)
    assert n == 9 and m == pytest.approx(1.0)
    assert lo == pytest.approx(1.0) and hi == pytest.approx(1.0)


def test_paired_ci_width_uses_t_df8_not_the_normal_quantile():
    v = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    m, lo, hi, n = P.paired_ci(v)
    se = float(np.std(v, ddof=1)) / 3.0
    assert n == 9
    assert (hi - lo) / 2 == pytest.approx(P.T8 * se)
    assert P.T8 == pytest.approx(2.306, abs=1e-3)


def test_paired_ci_drops_non_finite_folds():
    m, lo, hi, n = P.paired_ci([1.0, np.nan, 3.0])
    assert n == 2 and m == pytest.approx(2.0)


def test_paired_ci_on_too_few_folds_is_nan():
    assert P.paired_ci([1.0])[3] == 1
    assert np.isnan(P.paired_ci([1.0])[1])
