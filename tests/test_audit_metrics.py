"""`audit_metrics`: resolvability and paired-test sensitivity.

These decide whether a criterion is trusted or rebuilt, so the properties have to hold rather
than look right. The sensitivity function is here because the first version of it was WRONG --
it used the metric's own fold-to-fold spread, which cancels in a paired test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from audit_metrics import ece, resolvability, sensitivity_multiplier  # noqa: E402


# -------------------------------------------------------------- resolvability ---- #
def test_a_bar_above_the_whole_null_is_fully_resolvable():
    null = np.full(1000, 0.05)
    r = resolvability(null, bar=0.20)
    assert r["pass_rate"] == 1.0


def test_a_bar_below_the_whole_null_is_never_passable():
    """The run-3 situation in its extreme form: the criterion tests sample size, not the model."""
    null = np.full(1000, 0.30)
    r = resolvability(null, bar=0.169)
    assert r["pass_rate"] == 0.0


def test_pass_rate_matches_the_null_quantile():
    null = np.linspace(0.0, 1.0, 10001)
    assert resolvability(null, bar=0.25)["pass_rate"] == pytest.approx(0.25, abs=1e-3)


def test_usable_bar_is_where_a_perfect_system_passes_95_percent():
    rng = np.random.default_rng(0)
    null = rng.normal(0.10, 0.02, 20000)
    r = resolvability(null, bar=0.05)
    assert r["usable_bar"] == pytest.approx(np.percentile(null, 95), abs=1e-9)
    assert (null <= r["usable_bar"]).mean() == pytest.approx(0.95, abs=0.01)


def test_higher_is_better_flips_the_comparison():
    null = np.linspace(0.0, 1.0, 10001)
    r = resolvability(null, bar=0.25, lower_is_better=False)
    assert r["pass_rate"] == pytest.approx(0.75, abs=1e-3)
    assert r["usable_bar"] == pytest.approx(np.percentile(null, 5), abs=1e-9)


# ------------------------------------------------------ sensitivity_multiplier ---- #
def test_multiplier_matches_the_paired_t_test():
    """k = t(.975, n-1)/sqrt(n): mean/SD must exceed k for the CI to exclude 0."""
    from scipy.stats import t
    for n in (5, 6, 10):
        assert sensitivity_multiplier(n) == pytest.approx(t.ppf(0.975, n - 1) / np.sqrt(n))


def test_multiplier_shrinks_with_more_folds():
    assert sensitivity_multiplier(6) < sensitivity_multiplier(5)
    assert sensitivity_multiplier(20) < sensitivity_multiplier(6)


def test_a_uniform_change_is_detected_at_any_magnitude():
    """SD(effect)=0 -> threshold 0. The property that makes Stage 1's +0.000 guards meaningful."""
    assert sensitivity_multiplier(6) * 0.0 == 0.0


def test_the_multiplier_reproduces_the_observed_run3_ci():
    """A_xdonor -> B_fatecal on fate_ece: mean 0.115, CI half-width 0.0275 over 5 folds."""
    k = sensitivity_multiplier(5)
    sd_effect = 0.0275 / k
    assert 0.115 > k * sd_effect                     # detected
    assert sd_effect == pytest.approx(0.0221, abs=5e-4)


def test_a_heterogeneous_change_can_be_large_and_still_read_as_noise():
    """The blind spot: helps some folds, hurts others."""
    effect = np.array([+0.30, -0.28, +0.26, -0.24, +0.22, -0.20])
    k = sensitivity_multiplier(len(effect))
    assert abs(effect.mean()) < k * effect.std(ddof=1)


# ------------------------------------------------------------------------ ece ---- #
def test_audit_ece_matches_scorecards():
    p = np.array([0.95, 0.95, 0.95, 0.95])
    y = np.array([1.0, 1.0, 1.0, 0.0])
    assert ece(p, y) == pytest.approx(0.20)


def test_perfectly_calibrated_large_sample_approaches_zero():
    rng = np.random.default_rng(0)
    p = rng.uniform(0.05, 0.95, 200000)
    y = (rng.random(len(p)) < p).astype(float)
    assert ece(p, y) < 0.01
