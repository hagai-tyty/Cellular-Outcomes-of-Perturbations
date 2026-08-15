"""Unit tests for the clock-circularity test.

The verdict branches decide whether every recorded ΔAge number describes a prediction or a
tautology, so each one is exercised on constructed input where the answer is known.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "dcc", ROOT / "experiments" / "diag_clock_circularity.py")
dcc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dcc)


# ---- pearson ------------------------------------------------------------------------------ #
def test_pearson_of_an_exact_linear_map_is_one():
    x = np.linspace(-5, 5, 40)
    assert dcc.pearson(x, 3.0 * x + 7.0) == pytest.approx(1.0, abs=1e-12)


def test_pearson_is_sign_aware():
    x = np.linspace(-5, 5, 40)
    assert dcc.pearson(x, -2.0 * x) == pytest.approx(-1.0, abs=1e-12)


def test_pearson_is_undefined_on_a_constant_vector():
    """A constant `clock_panel` would otherwise yield a nan that reads as a low correlation and
    be misreported as NOT CIRCULAR."""
    assert np.isnan(dcc.pearson(np.full(10, 2.0), np.arange(10, dtype=float)))
    assert np.isnan(dcc.pearson(np.arange(10, dtype=float), np.full(10, 2.0)))


def test_pearson_is_undefined_on_too_few_points():
    assert np.isnan(dcc.pearson(np.array([1.0, 2.0]), np.array([1.0, 2.0])))


# ---- clock_weights_for_panel -------------------------------------------------------------- #
def test_weights_align_to_panel_order_not_dict_order():
    """The dot product is against X's COLUMNS, so a mis-ordered weight vector would silently
    score the wrong gene."""
    w, _ = dcc.clock_weights_for_panel(["B", "A"], {"A": 1.0, "B": 2.0})
    assert list(w) == [2.0, 1.0]


def test_genes_the_clock_never_saw_get_zero():
    """Mirrors LinearClock.predict_age (aging.py:55): absent gene -> 0.0, never dropped."""
    w, cov = dcc.clock_weights_for_panel(["A", "UNSEEN"], {"A": 1.0})
    assert list(w) == [1.0, 0.0]
    assert cov["n_panel"] == 2 and cov["n_with_weight"] == 1


def test_retained_mass_is_a_fraction_of_the_whole_clock():
    """Must be measured against the FULL clock, not just the panel, or a panel holding a sliver
    of the clock would report 1.0 and look complete."""
    _, cov = dcc.clock_weights_for_panel(["A"], {"A": 1.0, "B": 3.0})
    assert cov["abs_mass_retained"] == pytest.approx(0.25)


def test_retained_mass_uses_absolute_values():
    """Signed weights would cancel and understate the coverage."""
    _, cov = dcc.clock_weights_for_panel(["A", "B"], {"A": 1.0, "B": -1.0})
    assert cov["abs_mass_retained"] == pytest.approx(1.0)


# ---- verdict ------------------------------------------------------------------------------ #
def test_recoverable_label_and_matching_ridge_is_circular():
    assert dcc.verdict(0.99, 0.99) == "CIRCULAR"


def test_recoverable_label_with_diverging_ridge_is_label_recoverable():
    assert dcc.verdict(0.99, 0.10) == "LABEL-RECOVERABLE"


def test_an_unrecoverable_label_is_not_circular_whatever_ridge_does():
    """T1 gates the verdict: if the panel cannot reconstruct the clock, a high T2 cannot make the
    task a tautology."""
    assert dcc.verdict(0.10, 0.99) == "NOT CIRCULAR"
    assert dcc.verdict(0.10, 0.10) == "NOT CIRCULAR"


def test_the_threshold_is_inclusive_on_both_clauses():
    assert dcc.verdict(dcc.RHO_CIRCULAR, dcc.RHO_CIRCULAR) == "CIRCULAR"


def test_just_below_the_threshold_is_not_circular():
    assert dcc.verdict(dcc.RHO_CIRCULAR - 1e-9, 0.99) == "NOT CIRCULAR"


def test_a_nan_t1_is_undetermined_not_a_verdict():
    assert dcc.verdict(float("nan"), 0.99) == "UNDETERMINED"


def test_a_nan_t2_cannot_produce_circular():
    """An undefined T2 must degrade to LABEL-RECOVERABLE, never silently to CIRCULAR."""
    assert dcc.verdict(0.99, float("nan")) == "LABEL-RECOVERABLE"


def test_the_threshold_is_a_stated_constant():
    assert dcc.RHO_CIRCULAR == 0.95


# ---- the retracted headline ---------------------------------------------------------------- #
def test_the_loocv_runner_no_longer_claims_ranking_generalizes():
    """`run_loocv.py` printed ">>> Ranking generalizes across held-out donors: Spearman 0.40 <<<"
    as a headline for the paper. The RES it ranks maxed at 1.6e-4 in the only folds that produced
    a Spearman at all, and was exactly 0.0 elsewhere -- a correlation over floating-point residue.
    Folds where RES is constant return NaN and vanish from the mean rather than reporting the
    estimator is degenerate, so the surviving number is selected for having rounding noise."""
    src = (ROOT / "local_runners" / "run_loocv.py").read_text(encoding="utf-8")
    assert "Ranking generalizes across held-out donors" not in src.split("# This used to print")[-1] \
        or "no longer allowed" in src, "the retracted generalization headline is being printed again"
    assert "NOT a generalization claim" in src, (
        "the ranking line must carry its caveat: a NaN fold drops out of the mean silently")
    assert "res_max" in src, "the caveat must point at the scale check that exposed this"
