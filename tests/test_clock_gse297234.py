"""Unit tests for the GSE297234 clock test — pure functions, no GEO files.

`dedup_highest` and `per_cell_ages` are the two places a silent error would change every number
in the run without raising: the first decides WHICH row a duplicated gene symbol contributes, and
the second computes ages without ever densifying a 30k x 8k matrix, so an indexing or
normalisation slip would be invisible.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix

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


C = _load("clock_gse297234", "experiments/clock_gse297234.py")


class _Clock:
    """A stand-in with a known answer, so the arithmetic is checkable by hand."""

    def __init__(self, weights, intercept=0.0):
        self.weights, self.intercept = weights, float(intercept)

    def predict_age(self, expr, genes):
        w = np.array([self.weights.get(g, 0.0) for g in genes], float)
        return np.asarray(expr, float) @ w + self.intercept


# ------------------------------------------------------------------- dedup_highest ---- #
def test_duplicate_symbols_keep_the_highest_expressed_row():
    X = csc_matrix(np.array([[1.0, 1.0], [50.0, 50.0], [3.0, 3.0]]))
    Xd, genes = C.dedup_highest(X, ["A", "A", "B"])
    assert genes == ["A", "B"]
    assert Xd.toarray()[0].tolist() == [50.0, 50.0]      # the 50-count row, not the 1-count row


def test_dedup_is_a_no_op_when_symbols_are_unique():
    X = csc_matrix(np.array([[1.0], [2.0], [3.0]]))
    Xd, genes = C.dedup_highest(X, ["A", "B", "C"])
    assert genes == ["A", "B", "C"]
    assert Xd.shape == (3, 1)


def test_dedup_returns_rows_in_a_stable_sorted_order():
    """Row order must track the gene list, or every weight lines up with the wrong gene."""
    X = csc_matrix(np.array([[9.0], [1.0], [5.0], [2.0]]))
    Xd, genes = C.dedup_highest(X, ["B", "A", "B", "A"])
    assert genes == ["B", "A"]                            # original index order 0 then 3? no: 0,3
    assert Xd.toarray().ravel().tolist() == [9.0, 2.0]    # highest B (9) and highest A (2)


def test_dedup_handles_a_gene_that_is_entirely_zero():
    X = csc_matrix(np.array([[0.0, 0.0], [0.0, 0.0]]))
    Xd, genes = C.dedup_highest(X, ["A", "A"])
    assert genes == ["A"] and Xd.shape == (1, 2)


# ------------------------------------------------------------------ per_cell_ages ---- #
def test_per_cell_age_matches_a_hand_computed_value():
    """One gene, weight 1, one cell: age = log1p(1e4 * x / lib) + intercept."""
    X = csc_matrix(np.array([[4.0], [6.0]]))              # lib = 10
    clock = _Clock({"G1": 1.0, "G2": 0.0}, intercept=5.0)
    got = C.per_cell_ages(X, ["G1", "G2"], clock)
    assert got[0] == pytest.approx(np.log1p(1e4 * 4.0 / 10.0) + 5.0)


def test_per_cell_age_is_library_size_normalised():
    """Two cells with the same COMPOSITION but different depth must score identically."""
    X = csc_matrix(np.array([[1.0, 10.0], [3.0, 30.0]]))
    clock = _Clock({"G1": 0.7, "G2": -0.2}, intercept=1.0)
    got = C.per_cell_ages(X, ["G1", "G2"], clock)
    assert got[0] == pytest.approx(got[1])


def test_per_cell_age_ignores_genes_the_clock_does_not_weight():
    X = csc_matrix(np.array([[5.0], [5.0]]))
    a = C.per_cell_ages(X, ["G1", "UNKNOWN"], _Clock({"G1": 1.0}))
    b = C.per_cell_ages(X, ["G1", "G2"], _Clock({"G1": 1.0}))
    assert a[0] == pytest.approx(b[0])


def test_per_cell_age_survives_an_empty_cell_without_dividing_by_zero():
    X = csc_matrix(np.array([[0.0, 4.0], [0.0, 6.0]]))
    got = C.per_cell_ages(X, ["G1", "G2"], _Clock({"G1": 1.0, "G2": 1.0}, intercept=2.0))
    assert np.isfinite(got).all()
    assert got[0] == pytest.approx(2.0)                   # empty cell -> intercept only


def test_per_cell_age_returns_one_value_per_cell():
    X = csc_matrix(np.random.default_rng(0).poisson(2.0, size=(6, 9)).astype(float))
    got = C.per_cell_ages(X, [f"G{i}" for i in range(6)], _Clock({"G0": 1.0}))
    assert got.shape == (9,)


# -------------------------------------------------------------------- pseudobulk ---- #
def test_pseudobulk_sums_across_cells_before_normalising():
    """Not the mean of per-cell values -- summing first is what makes it comparable to bulk."""
    X = csc_matrix(np.array([[1.0, 3.0], [1.0, 1.0]]))    # gene totals 4 and 2, lib 6
    clock = _Clock({"G1": 1.0, "G2": 0.0})
    got = C.pseudobulk_age(X, ["G1", "G2"], clock)
    assert got == pytest.approx(np.log1p(1e4 * 4.0 / 6.0))


def test_pseudobulk_is_invariant_to_cell_order():
    rng = np.random.default_rng(1)
    A = rng.poisson(3.0, size=(5, 7)).astype(float)
    clock = _Clock({f"G{i}": float(i) for i in range(5)})
    g = [f"G{i}" for i in range(5)]
    a = C.pseudobulk_age(csc_matrix(A), g, clock)
    b = C.pseudobulk_age(csc_matrix(A[:, ::-1]), g, clock)
    assert a == pytest.approx(b)


# ------------------------------------------------------------------------ config ---- #
def test_the_two_donors_and_their_ages_are_pinned():
    ages = {s["line"]: s["age"] for s in C.SAMPLES}
    assert ages == {"GM23815": 22.0, "GM00731": 96.0}
    assert abs(ages["GM00731"] - ages["GM23815"]) == 74.0


def test_only_day_zero_files_are_referenced():
    """Any later day would put reprogrammed cells into an age-calibration test."""
    assert all("_D0_" in s["file"] for s in C.SAMPLES)
