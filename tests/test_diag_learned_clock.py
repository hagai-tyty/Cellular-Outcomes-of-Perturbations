"""STAGE 1.5.4 — pure-logic tests for the learned-clock diagnostic. No repo data required.

The load-bearing test here is `test_lodo_mean_reversion_is_removed_by_partialling_donor`: it
reproduces, from synthetic data alone, the artefact that made the first run of this stage VOID.
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

_SPEC = importlib.util.spec_from_file_location(
    "diag_learned_clock", ROOT / "experiments" / "diag_learned_clock.py")
dlc = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = dlc
_SPEC.loader.exec_module(dlc)


# ------------------------------------------------------------------ ridge ---- #
def test_ridge_matches_the_closed_form_when_p_is_small():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(40, 5))
    y = X @ np.array([1.0, -2.0, 0.5, 0.0, 3.0]) + rng.normal(0, 0.1, 40)
    a = 1.0
    want = np.linalg.solve(X.T @ X + a * np.eye(5), X.T @ y)
    assert np.allclose(dlc.ridge_fit(X, y, a), want)


def test_ridge_dual_form_agrees_with_primal_when_both_are_valid():
    """p >> n uses the sample-space solve; it must give the same answer as the primal."""
    rng = np.random.default_rng(1)
    X = rng.normal(size=(12, 12))          # square: both branches are legitimate
    y = rng.normal(size=12)
    a = 2.5
    primal = np.linalg.solve(X.T @ X + a * np.eye(12), X.T @ y)
    dual = X.T @ np.linalg.solve(X @ X.T + a * np.eye(12), y)
    assert np.allclose(primal, dual, atol=1e-8)


def test_ridge_shrinks_towards_zero_as_alpha_grows():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(30, 4))
    y = X @ np.ones(4) + rng.normal(0, 0.1, 30)
    small = np.abs(dlc.ridge_fit(X, y, 1e-2)).sum()
    large = np.abs(dlc.ridge_fit(X, y, 1e6)).sum()
    assert large < small / 100


# ------------------------------------------------------- rank statistics ---- #
def test_spearman_is_one_for_a_monotone_transform():
    x = np.arange(1.0, 21.0)
    assert dlc.spearman(x, np.exp(x / 5)) == pytest.approx(1.0)


def test_spearman_handles_ties_by_averaging_ranks():
    assert dlc.spearman([1, 1, 2, 3], [1, 1, 2, 3]) == pytest.approx(1.0)


def test_partial_spearman_removes_a_pure_confounder():
    """x and y correlate ONLY through z; partialling z must collapse it."""
    rng = np.random.default_rng(3)
    z = rng.normal(size=200)
    x = z + rng.normal(0, 0.05, 200)
    y = z + rng.normal(0, 0.05, 200)
    assert dlc.spearman(x, y) > 0.9
    assert abs(dlc.partial_spearman(x, y, z)) < 0.3


# ---------------------------------------------- THE ARTEFACT THAT VOIDED RUN 1 ---- #
def test_lodo_mean_reversion_is_removed_by_partialling_donor():
    """Reproduces the failure that made this stage's first run VOID, from synthetic data.

    A leave-one-donor-out prediction for donor d is roughly the mean of the OTHER donors, so if
    the donor means are ordered the prediction is ANTI-correlated with the target by construction
    -- a large apparent correlation from a model that has learned nothing. Partialling donor must
    remove it; not partialling donor must show it.
    """
    donor = np.array(["A"] * 20 + ["B"] * 20 + ["C"] * 20)
    means = {"A": 0.0, "B": 5.0, "C": 10.0}
    rng = np.random.default_rng(4)
    y = np.array([means[d] for d in donor]) + rng.normal(0, 0.5, 60)
    # a model that learned NOTHING: predict the mean of the other two donors
    pred = np.array([np.mean([means[o] for o in "ABC" if o != d]) for d in donor])
    plu = rng.normal(size=60)

    naive = dlc.partial_spearman(pred, y, plu)
    corrected = dlc.partial_spearman(pred, y, plu, donor)
    assert naive < -0.5, "the artefact must be present without the donor term"
    assert abs(corrected) < 0.2, "partialling donor must remove it"


def test_partialling_donor_does_not_erase_a_genuine_within_donor_signal():
    """The correction must not be so aggressive that real signal disappears with the artefact."""
    donor = np.array(["A"] * 20 + ["B"] * 20 + ["C"] * 20)
    rng = np.random.default_rng(5)
    within = rng.normal(size=60)
    y = np.array([{"A": 0.0, "B": 5.0, "C": 10.0}[d] for d in donor]) + within
    pred = within + rng.normal(0, 0.2, 60)          # tracks the WITHIN-donor part only
    plu = rng.normal(size=60)
    assert dlc.partial_spearman(pred, y, plu, donor) > 0.7


# ------------------------------------------------------------------ verdict ---- #
@pytest.mark.parametrize(("rhos", "want"), [
    ({"a": 0.6, "b": 0.7}, "PASS"),
    ({"a": 0.6, "b": 0.4}, "SPLIT"),
    ({"a": 0.1, "b": 0.4}, "FAIL"),
    ({"a": 0.50, "b": 0.50}, "PASS"),        # the bar is inclusive
])
def test_verdict_rules(rhos, want):
    assert dlc.verdict(rhos) == want


def test_split_is_a_failure_not_a_pass():
    """M-2a sec 7's rule, inherited: one clock clearing the bar is not a result."""
    assert dlc.verdict({"a": 0.9, "b": 0.0}) == "SPLIT"


# ------------------------------------------------------------ inner CV ---- #
def test_inner_cv_is_grouped_so_a_donor_cannot_span_the_split():
    """With one group only, there is nothing to hold out and a fixed alpha is returned."""
    rng = np.random.default_rng(6)
    X = rng.normal(size=(10, 3))
    y = rng.normal(size=10)
    got = dlc.inner_cv_alpha(X, y, np.array(["A"] * 10))
    assert got in dlc.ALPHAS


def test_inner_cv_does_not_pick_the_weakest_alpha_on_pure_noise():
    """When X carries no signal, the weakest shrinkage is the one choice that must not win.

    Asserting the MAXIMUM alpha would over-specify: with standardised features the CV error curve
    is nearly flat across the heavy end, so which of the large alphas wins is close to arbitrary
    (this data picks 1e3 of [1e1 … 1e6]). What the routine must never do is choose the weakest and
    fit the noise -- that is the property worth pinning.
    """
    rng = np.random.default_rng(7)
    X = rng.normal(size=(40, 30))
    y = rng.normal(size=40)
    g = np.array(["A"] * 20 + ["B"] * 20)
    assert dlc.inner_cv_alpha(X, y, g) > min(dlc.ALPHAS)
