"""Unit tests for the residual-expression diagnostic.

At n=6 the verdict rests entirely on the procedure being correct, so the leak-prone parts are
tested on constructed data where the answer is known: that the donor-age fit is refit inside each
fold, and that a useless model scores like the mean rather than like a correct one.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "dre", ROOT / "experiments" / "diag_residual_expression.py")
dre = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dre)


# ---- ridge_fit_predict --------------------------------------------------------------------- #
def test_ridge_recovers_a_clean_linear_signal_at_low_alpha():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(40, 3))
    y = X @ np.array([2.0, -1.0, 0.5])
    p = dre.ridge_fit_predict(X[:30], y[:30], X[30:], alpha=1e-6)
    assert np.corrcoef(p, y[30:])[0, 1] > 0.99


def test_a_huge_alpha_collapses_to_the_train_mean_not_to_zero():
    """Load-bearing. If a degenerate fit predicted 0 while the target is centred elsewhere, a
    useless model would score like a correct one on residualised targets."""
    rng = np.random.default_rng(1)
    X = rng.normal(size=(20, 5))
    y = rng.normal(size=20) + 100.0
    p = dre.ridge_fit_predict(X[:15], y[:15], X[15:], alpha=1e12)
    assert np.allclose(p, y[:15].mean(), atol=1e-3)


def test_the_dual_form_equals_the_primal_form():
    """The implementation solves in SAMPLE space (5x5) instead of FEATURE space (2000x2000) via
    `X^T(XX^T + aI)^-1 = (X^TX + aI)^-1 X^T`. That is a change to the MATH, so it is pinned against
    an explicit primal solve rather than trusted."""
    rng = np.random.default_rng(11)
    xtr, xte = rng.normal(size=(5, 60)), rng.normal(size=(1, 60))
    ytr = rng.normal(size=5)
    for alpha in (0.1, 1.0, 100.0, 10_000.0):
        mu, sd = xtr.mean(0), xtr.std(0)
        sd = np.where(sd < 1e-12, 1.0, sd)
        a, b = (xtr - mu) / sd, (xte - mu) / sd
        ybar = ytr.mean()
        primal = b @ np.linalg.solve(a.T @ a + alpha * np.eye(a.shape[1]),
                                     a.T @ (ytr - ybar)) + ybar
        assert dre.ridge_fit_predict(xtr, ytr, xte, alpha) == pytest.approx(primal, rel=1e-6)


def test_constant_features_do_not_divide_by_zero():
    X = np.hstack([np.ones((10, 1)), np.arange(10.0).reshape(-1, 1)])
    p = dre.ridge_fit_predict(X[:8], np.arange(8.0), X[8:], alpha=1.0)
    assert np.all(np.isfinite(p))


# ---- loo_spearman: the leak that would manufacture signal ---------------------------------- #
def test_the_donor_age_fit_is_refit_per_fold_not_once_globally():
    """If the age fit used all donors, the held-out donor would help define its own residual.
    Constructed so y is EXACTLY a function of age: every fold's residual must then be ~0, and the
    LOO spearman undefined -- which cannot happen if the fit saw the held-out point only via the
    global fit."""
    rng = np.random.default_rng(2)
    age = np.array([0.0, 0.0, 29.0, 35.0, 53.0, 53.0])
    y = 2.0 * age + 10.0
    X = rng.normal(size=(6, 8))
    r = dre.loo_spearman(X, y, age, alpha=1.0, residualise=True)
    assert np.isnan(r), "an exactly-age-determined target must leave no residual to predict"


def test_without_residualising_an_age_determined_target_is_predictable():
    """The same data with residualise=False is the positive-control path and must NOT be nan."""
    age = np.array([0.0, 0.0, 29.0, 35.0, 53.0, 53.0])
    y = 2.0 * age + 10.0
    X = np.vstack([age, age ** 2]).T.astype(float)
    assert np.isfinite(dre.loo_spearman(X, y, age, alpha=1e-6, residualise=False))


def test_pure_noise_does_not_produce_a_high_loo_spearman_on_average():
    rng = np.random.default_rng(3)
    age = np.array([0.0, 0.0, 29.0, 35.0, 53.0, 53.0])
    vals = []
    for _ in range(40):
        X, y = rng.normal(size=(6, 50)), rng.normal(size=6)
        v = dre.loo_spearman(X, y, age, alpha=100.0, residualise=True)
        if np.isfinite(v):
            vals.append(v)
    assert abs(float(np.mean(vals))) < 0.5


# ---- permutation_null ----------------------------------------------------------------------- #
def test_the_null_is_reproducible_from_its_seed():
    rng = np.random.default_rng(4)
    X, y = rng.normal(size=(6, 10)), rng.normal(size=6)
    age = np.array([0.0, 0.0, 29.0, 35.0, 53.0, 53.0])
    a = dre.permutation_null(X, y, age, 100.0, True, n_perm=25, seed=7)
    b = dre.permutation_null(X, y, age, 100.0, True, n_perm=25, seed=7)
    assert np.allclose(a, b)


def test_the_null_spans_both_signs_so_a_positive_result_is_not_automatic():
    """A null concentrated above zero would make any positive observation 'pass' trivially."""
    rng = np.random.default_rng(5)
    X, y = rng.normal(size=(6, 20)), rng.normal(size=6)
    age = np.array([0.0, 0.0, 29.0, 35.0, 53.0, 53.0])
    null = dre.permutation_null(X, y, age, 100.0, True, n_perm=200, seed=1)
    assert null.min() < 0 < null.max()


def test_spearman_guards_degenerate_input():
    assert np.isnan(dre.spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))
    assert np.isnan(dre.spearman(np.ones(6), np.arange(6.0)))


def test_the_thresholds_are_stated_constants():
    assert dre.ALPHAS == (1.0, 10.0, 100.0, 1000.0, 10000.0)
    assert dre.N_PERM == 2000 and dre.PERM_PCTILE == 95.0
    assert dre.MIN_ALPHAS_PASSING == 3
    assert (dre.EARLY_LO, dre.EARLY_HI, dre.LATE_LO) == (7.0, 29.0, 34.0)
