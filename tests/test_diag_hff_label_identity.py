"""STAGE 1.5.5 — pure-logic tests for the HFF label-provenance diagnostic. No repo data required."""

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

_SPEC = importlib.util.spec_from_file_location(
    "diag_hff_label_identity", ROOT / "experiments" / "diag_hff_label_identity.py")
hli = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = hli
_SPEC.loader.exec_module(hli)


def test_spearman_is_one_for_a_monotone_transform():
    x = np.arange(1.0, 31.0)
    assert hli.spearman(x, np.exp(x / 7)) == pytest.approx(1.0)


def test_spearman_is_zero_for_independent_data():
    rng = np.random.default_rng(0)
    assert abs(hli.spearman(rng.normal(size=4000), rng.normal(size=4000))) < 0.06


def test_r2_is_one_when_y_is_a_monotone_function_of_a_predictor():
    """r2_on ranks its inputs, so a monotone relationship must come out at 1."""
    x = np.linspace(0, 10, 300)
    assert hli.r2_on(np.exp(x), np.column_stack([x, np.zeros_like(x)])) == pytest.approx(1.0, abs=1e-6)


def test_r2_is_near_zero_for_unrelated_predictors():
    rng = np.random.default_rng(1)
    y = rng.normal(size=3000)
    X = rng.normal(size=(3000, 2))
    assert hli.r2_on(y, X) < 0.01


def test_r2_grows_when_a_real_predictor_is_added():
    rng = np.random.default_rng(2)
    a = rng.normal(size=2000)
    y = a + rng.normal(0, 0.5, 2000)
    noise = rng.normal(size=2000)
    assert hli.r2_on(y, noise.reshape(-1, 1)) < 0.02
    assert hli.r2_on(y, np.column_stack([noise, a])) > 0.5


# ------------------------------------------------- the pre-registered decision rule ---- #
def _pd(rhos):
    return {str(i): {"rho_pluri": r} for i, r in enumerate(rhos)}


def test_verdict_identity_dominated_needs_the_registered_majority():
    assert hli.decide(_pd([0.9] * 6 + [0.0] * 2))["verdict"] == "IDENTITY_DOMINATED"
    assert hli.decide(_pd([0.9] * 5 + [0.0] * 3))["verdict"] == "NOT_DOMINATED"


def test_the_bar_is_on_absolute_rho_so_sign_cannot_rescue_a_label():
    """A label that tracks pluripotency NEGATIVELY is just as much an identity readout."""
    assert hli.decide(_pd([-0.9] * 6 + [0.0] * 2))["verdict"] == "IDENTITY_DOMINATED"


def test_the_bar_is_inclusive_at_exactly_the_threshold():
    assert hli.decide(_pd([hli.RHO_BAR] * 6 + [0.0] * 2))["verdict"] == "IDENTITY_DOMINATED"


def test_nan_timepoints_never_count_toward_the_bar():
    assert hli.decide(_pd([float("nan")] * 8))["n_timepoints_over_bar"] == 0


def test_the_observed_result_is_far_from_the_bar_in_every_timepoint():
    """Regression on the recorded run: the verdict was not a near miss.

    Measured |rho(y_age, pluripotency)| within timepoint ranged 0.036-0.270 against a 0.50 bar, so
    no timepoint came close. Pinned so a future change that inflates these correlations -- e.g. by
    accidentally pooling timepoints, which gives -0.216 -- is caught rather than read as a finding.
    """
    observed = [0.079, -0.079, 0.036, -0.158, -0.130, -0.128, -0.270, -0.194]
    v = hli.decide(_pd(observed))
    assert v["verdict"] == "NOT_DOMINATED"
    assert v["n_timepoints_over_bar"] == 0
    assert max(abs(r) for r in observed) < hli.RHO_BAR
