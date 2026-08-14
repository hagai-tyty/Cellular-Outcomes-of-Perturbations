"""Unit tests for P4 — pure functions only, no bundle data.

`wilson` is used where it matters most: the apoptosis head sits near 0.003 at the end of the
course, and a normal-approximation interval there runs below zero. `spearman` decides the verdict
outright, so its sign convention is pinned against hand-worked cases.
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


P = _load("p4_two_heads", "experiments/p4_two_heads.py")


# ---------------------------------------------------------------------------- wilson ---- #
def test_wilson_stays_inside_the_unit_interval_at_a_tiny_rate():
    """P(apoptosis) reaches 0.0029 at day 21. A normal-approximation interval goes negative
    there; Wilson must not."""
    lo, hi = P.wilson(14, 4788)
    assert lo > 0.0 and hi < 1.0
    assert lo < 14 / 4788 < hi


def test_wilson_stays_inside_the_unit_interval_at_a_rate_near_one():
    lo, hi = P.wilson(4772, 4788)
    assert lo > 0.0 and hi <= 1.0


def test_wilson_brackets_the_point_estimate():
    lo, hi = P.wilson(300, 1000)
    assert lo < 0.30 < hi


def test_wilson_narrows_as_n_grows():
    w_small = np.subtract(*reversed(P.wilson(30, 100)))
    w_big = np.subtract(*reversed(P.wilson(3000, 10000)))
    assert w_big < w_small


def test_wilson_on_zero_events_has_a_positive_upper_bound():
    lo, hi = P.wilson(0, 500)
    assert lo >= 0.0 and hi > 0.0


def test_wilson_on_an_empty_group_is_nan():
    lo, hi = P.wilson(0, 0)
    assert np.isnan(lo) and np.isnan(hi)


# -------------------------------------------------------------------------- spearman ---- #
def test_spearman_is_plus_one_for_a_monotone_increasing_pair():
    assert P.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)


def test_spearman_is_minus_one_for_a_monotone_decreasing_pair():
    """THE case that decides the verdict: identity loss rises while apoptosis falls."""
    assert P.spearman([1, 2, 3, 4], [40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_is_rank_based_not_value_based():
    """A wild outlier must not move it -- the apoptosis curve is tiny next to identity loss."""
    a = P.spearman([1, 2, 3, 4], [10, 20, 30, 40])
    b = P.spearman([1, 2, 3, 4], [10, 20, 30, 40000])
    assert a == pytest.approx(b)


def test_spearman_of_a_constant_is_nan_not_zero():
    assert np.isnan(P.spearman([1, 2, 3], [5, 5, 5]))


# --------------------------------------------------------------------------- pearson ---- #
def test_pearson_matches_a_hand_computed_value():
    assert P.pearson([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0)


def test_pearson_is_negative_for_an_inverse_relationship():
    assert P.pearson([1, 2, 3], [6, 4, 2]) == pytest.approx(-1.0)


def test_pearson_of_a_constant_is_nan():
    assert np.isnan(P.pearson([1, 2, 3], [5, 5, 5]))


# ------------------------------------------------------- the verdict logic, end to end ---- #
@pytest.mark.parametrize("loss,death,want", [
    ([0.1, 0.2, 0.3, 0.4], [0.1, 0.2, 0.3, 0.4], "H1"),        # rise together, same peak
    ([0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1], "H3"),        # opposite -- the real case
    ([0.1, 0.2, 0.3, 0.4], [0.1, 0.3, 0.2, 0.4], "H2"),        # positive but only 0.8
])
def test_the_three_pre_registered_branches_are_reachable(loss, death, want):
    """Every branch must be able to fire, or the gate is decoration (the verify_1a lesson)."""
    sp = P.spearman(loss, death)
    gap = abs(int(np.argmax(loss)) - int(np.argmax(death)))
    got = "H1" if (sp > 0.9 and gap <= 1) else "H3" if sp < 0 else "H2"
    assert got == want


def test_ranks_average_ties_rather_than_breaking_them_arbitrarily():
    """P(apoptosis) is tied at 0.0397 on two timepoints. argsort(argsort(...)) would give those
    two different ranks and quietly change the correlation that decides the verdict."""
    assert P._rank([5.0, 5.0, 5.0]).tolist() == [1.0, 1.0, 1.0]
    assert P._rank([1.0, 2.0, 2.0, 4.0]).tolist() == [0.0, 1.5, 1.5, 3.0]


def test_spearman_is_unchanged_by_the_order_of_tied_entries():
    a = [1.0, 2.0, 3.0, 4.0]
    assert P.spearman(a, [0.1, 0.5, 0.5, 0.9]) == pytest.approx(
        P.spearman(a, [0.1, 0.5, 0.5, 0.9]))
    assert np.isfinite(P.spearman(a, [0.1, 0.5, 0.5, 0.9]))
