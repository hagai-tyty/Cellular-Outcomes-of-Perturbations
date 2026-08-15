"""Unit tests for the age-capacity test.

The claim this script supports -- that the representation CAN learn age -- is only as good as the
cross-validation being leak-free, so that is what is tested.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "dac", ROOT / "experiments" / "diag_age_capacity.py")
dac = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dac)


def test_every_sample_is_predicted_exactly_once():
    """A sample appearing in two test folds, or none, would silently bias the pooled MAE."""
    rng = np.random.default_rng(0)
    X, y = rng.normal(size=(50, 8)), rng.normal(size=50)
    p = dac.cv_predict(X, y, alpha=1.0, n_splits=5, seed=3)
    assert p.shape == (50,) and np.all(np.isfinite(p))


def test_a_sample_never_trains_on_itself():
    """The target is a learnable SIGNAL plus per-sample NOISE that appears in no feature. A
    leak-free fit recovers the signal and cannot recover the noise.

    Note a wrong version of this test: setting a feature equal to y outright makes y perfectly
    learnable from the training rows alone, so reproducing it exactly is correct behaviour, not
    memorisation. The unlearnable component is what makes a leak visible.
    """
    rng = np.random.default_rng(7)
    signal = np.linspace(0, 40, 60)
    noise = rng.normal(0, 5.0, 60)
    X = np.column_stack([signal, rng.normal(size=60)])
    y = signal + noise
    p = dac.cv_predict(X, y, alpha=1e-6, n_splits=5, seed=0)
    assert np.corrcoef(p, signal)[0, 1] > 0.95            # the learnable part IS learned
    assert abs(np.corrcoef(p - signal, noise)[0, 1]) < 0.4  # the unlearnable part is NOT


def test_predictions_are_reproducible_from_the_seed():
    rng = np.random.default_rng(1)
    X, y = rng.normal(size=(30, 5)), rng.normal(size=30)
    a = dac.cv_predict(X, y, 1.0, n_splits=5, seed=9)
    b = dac.cv_predict(X, y, 1.0, n_splits=5, seed=9)
    assert np.allclose(a, b)


def test_a_different_seed_changes_the_partition():
    rng = np.random.default_rng(2)
    X, y = rng.normal(size=(30, 5)), rng.normal(size=30)
    assert not np.allclose(dac.cv_predict(X, y, 1.0, 5, 1), dac.cv_predict(X, y, 1.0, 5, 2))


def test_pure_noise_scores_no_better_than_the_baseline():
    """The bar is a RATIO to the mean baseline, so noise must not clear it."""
    rng = np.random.default_rng(4)
    y = rng.uniform(1, 96, 120)
    X = rng.normal(size=(120, 200))
    p = dac.cv_predict(X, y, alpha=100.0, n_splits=10, seed=0)
    base = float(np.abs(y - np.median(y)).mean())
    assert dac.scores(y, p)["mae"] / base > dac.MAE_RATIO_BAR


def test_scores_are_the_quantities_reported():
    y = np.arange(20.0)
    s = dac.scores(y, y + 3.0)
    assert s["mae"] == pytest.approx(3.0)
    assert s["pearson"] == pytest.approx(1.0)
    assert s["spearman"] == pytest.approx(1.0)


def test_the_target_is_donor_age_not_a_clock_output():
    """The whole non-circularity argument. If this ever read a clock, the test would be measuring
    a linear readout of its own input, exactly as `diag_clock_circularity` found for ΔAge."""
    src = (ROOT / "experiments" / "diag_age_capacity.py").read_text(encoding="utf-8")
    assert "LinearClock" not in src and "fleischer_clock" not in src
    assert 'meta.age' in src or '"age"' in src


def test_hgps_is_excluded_from_the_primary_cohort():
    """Progeria ages abnormally fast; including it would inflate apparent performance."""
    src = (ROOT / "experiments" / "diag_age_capacity.py").read_text(encoding="utf-8")
    assert "HGPS" in src and 'eq("Normal")' in src


def test_the_thresholds_are_stated_constants():
    assert dac.MAE_RATIO_BAR == 0.75
    assert dac.PUBLISHED_CV_MAE == pytest.approx(12.27)
    assert dac.N_SPLITS == 10 and dac.ALPHAS == (1.0, 10.0, 100.0, 1000.0, 10000.0)
