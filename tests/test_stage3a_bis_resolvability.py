"""Unit tests for the Stage 3a-bis §5b resolvability check — pure functions, no repo data.

The load-bearing claim of 3a-bis is a CONTROLLED COMPARISON: regime A and regime B share almost
the same training set (679 vs 643 pairs, both containing every HFF pseudo-replicate) and differ
only in how precisely the HELD-OUT trajectory is measured. That comparison is only valid if the
pieces below behave exactly as claimed, so each is pinned as a property:

  * the pseudo-replicate partition must be balanced WITHIN each timepoint. A partition that split
    cells globally would hand some replicates only early timepoints, and they would no longer be
    trajectories at all -- regime B would then be measuring the partition, not the geometry.
  * `pairs_of` must emit strictly forward pairs. A single backwards pair would leak the future.
  * `hff_curve` must pool by CELL COUNT, not average the replicate means, or the simulated truth
    would not be HFF's measured curve.
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


B = _load("stage3a_bis_resolvability", "experiments/stage3a_bis_resolvability.py")


# ------------------------------------------------- partition_within_timepoints ---- #
def _t(counts):
    """A time vector with `counts[i]` cells at timepoint i."""
    return np.concatenate([np.full(n, float(i)) for i, n in enumerate(counts)])


def test_partition_is_balanced_within_every_timepoint():
    """The whole point: each replicate gets ~1/k of EACH timepoint, not 1/k of the cells."""
    t = _t([100, 100, 100])
    sel = np.ones(len(t), bool)
    rep = B.partition_within_timepoints(t, sel, 10, seed=0)
    for tp in np.unique(t):
        at_tp = np.isclose(t, tp)
        for k in range(10):
            assert (at_tp & (rep == k)).sum() == 10


def test_every_replicate_carries_the_whole_time_course():
    """A replicate missing a timepoint is not a trajectory and would distort regime B."""
    t = _t([31, 47, 23, 55])
    rep = B.partition_within_timepoints(t, np.ones(len(t), bool), 5, seed=1)
    for k in range(5):
        assert len(np.unique(t[rep == k])) == 4


def test_partition_is_a_true_partition_of_the_selected_cells():
    t = _t([20, 30])
    sel = np.ones(len(t), bool)
    rep = B.partition_within_timepoints(t, sel, 4, seed=2)
    assert (rep >= 0).sum() == sel.sum()
    assert sorted(np.bincount(rep[rep >= 0])) == sorted(np.bincount(rep[rep >= 0]))
    assert set(np.unique(rep)) <= set(range(4))


def test_partition_never_touches_unselected_cells():
    """Gill cells must keep -1 so they are not swept into an HFF pseudo-replicate."""
    t = _t([10, 10])
    sel = np.zeros(len(t), bool)
    sel[:10] = True
    rep = B.partition_within_timepoints(t, sel, 3, seed=3)
    assert (rep[10:] == -1).all()
    assert (rep[:10] >= 0).all()


def test_partition_is_deterministic_for_a_seed_and_varies_across_seeds():
    t = _t([40, 40])
    sel = np.ones(len(t), bool)
    a = B.partition_within_timepoints(t, sel, 4, seed=7)
    assert np.array_equal(a, B.partition_within_timepoints(t, sel, 4, seed=7))
    assert not np.array_equal(a, B.partition_within_timepoints(t, sel, 4, seed=8))


def test_partition_handles_a_timepoint_with_fewer_cells_than_replicates():
    """3 cells across 10 replicates: 7 replicates simply miss that timepoint, no crash."""
    t = _t([3, 50])
    rep = B.partition_within_timepoints(t, np.ones(len(t), bool), 10, seed=0)
    assert (rep[:3] >= 0).all()
    assert len({int(v) for v in rep[:3]}) == 3


# ------------------------------------------------------------------- pairs_of ---- #
def _rows(ts, us=None, ns=None):
    us = us if us is not None else [0.0] * len(ts)
    ns = ns if ns is not None else [1] * len(ts)
    return [{"t": float(t), "x": np.array([float(t)]), "u": float(u), "n": int(n)}
            for t, u, n in zip(ts, us, ns, strict=True)]


def test_pairs_are_strictly_forward():
    for p in B.pairs_of(_rows([0.0, 1.0, 2.0, 3.0])):
        assert p["dt"] > 0
        assert p["t_j"] > p["t_j"] - p["dt"]


def test_pair_count_is_n_choose_2():
    for n in (3, 4, 9, 12):
        assert len(B.pairs_of(_rows(list(range(n))))) == n * (n - 1) // 2


def test_pairs_carry_the_endpoint_cell_count_for_the_binomial_draw():
    """`n_j` is what turns Gill's 1-2 cells into the real observation noise."""
    ps = B.pairs_of(_rows([0.0, 1.0], ns=[500, 3]))
    assert len(ps) == 1 and ps[0]["n_j"] == 3


def test_pairs_are_empty_when_every_timepoint_is_identical():
    assert B.pairs_of(_rows([1.0, 1.0, 1.0])) == []


def test_pairs_are_unaffected_by_row_order():
    a = sorted(p["dt"] for p in B.pairs_of(_rows([0.0, 2.0, 5.0])))
    b = sorted(p["dt"] for p in B.pairs_of(_rows([5.0, 0.0, 2.0])))
    assert a == pytest.approx(b)


# ----------------------------------------------------------------- hff_curve ---- #
def test_hff_curve_pools_by_cell_count_not_by_replicate_mean():
    """One replicate with 999 cells at u=1 and one with 1 cell at u=0 pools to ~0.999."""
    fold = {"train": {
        "HFF_r0": _rows([0.0, 1.0], us=[1.0, 1.0], ns=[999, 999]),
        "HFF_r1": _rows([0.0, 1.0], us=[0.0, 0.0], ns=[1, 1]),
    }}
    ts, g = B.hff_curve(fold)
    assert list(ts) == [0.0, 1.0]
    assert g[0] == pytest.approx(999 / 1000)


def test_hff_curve_ignores_non_hff_trajectories():
    """Gill donors must not enter the simulated truth -- it is HFF's curve by construction."""
    fold = {"train": {
        "HFF_r0": _rows([0.0, 1.0], us=[0.2, 0.8], ns=[100, 100]),
        "N3": _rows([0.0, 1.0], us=[1.0, 0.0], ns=[1000, 1000]),
    }}
    _, g = B.hff_curve(fold)
    assert g == pytest.approx([0.2, 0.8])


def test_hff_curve_is_sorted_by_time():
    fold = {"train": {"HFF_r0": _rows([5.0, 0.0, 2.0], us=[0.9, 0.1, 0.5], ns=[10, 10, 10])}}
    ts, g = B.hff_curve(fold)
    assert list(ts) == [0.0, 2.0, 5.0]
    assert g == pytest.approx([0.1, 0.5, 0.9])


# ----------------------------------------------------------- trajectory_rows ---- #
def _traj_inputs(n_per_tp=(4, 4), unsafe=(0.0, 1.0), age_valid=True):
    t = _t(list(n_per_tp))
    X = np.arange(len(t) * 2, dtype=np.float32).reshape(len(t), 2)
    cls = np.zeros(len(t), int)
    y = np.arange(len(t), dtype=float)
    mask = np.full(len(t), age_valid)
    u = np.concatenate([np.full(n, float(v)) for n, v in zip(n_per_tp, unsafe, strict=True)])
    return t, X, cls, y, mask, u


def test_trajectory_rows_reports_the_cell_count_and_unsafe_fraction_per_timepoint():
    t, X, cls, y, mask, u = _traj_inputs((4, 6), (0.0, 1.0))
    rows = B.trajectory_rows(t, X, cls, y, mask, np.ones(len(t), bool), u)
    assert [r["n"] for r in rows] == [4, 6]
    assert [r["u"] for r in rows] == [0.0, 1.0]


def test_trajectory_rows_returns_nan_age_when_every_cell_is_masked():
    """C-7 rule 4 masks a whole donor's dAge; the safety target must still be computed."""
    t, X, cls, y, mask, u = _traj_inputs(age_valid=False)
    rows = B.trajectory_rows(t, X, cls, y, mask, np.ones(len(t), bool), u)
    assert all(np.isnan(r["y"]) for r in rows)
    assert all(r["n_age"] == 0 for r in rows)
    assert all(np.isfinite(r["u"]) for r in rows)


def test_trajectory_rows_averages_expression_in_double_precision():
    t, X, cls, y, mask, u = _traj_inputs()
    rows = B.trajectory_rows(t, X, cls, y, mask, np.ones(len(t), bool), u)
    assert all(r["x"].dtype == np.float64 for r in rows)


def test_trajectory_rows_only_sees_the_selected_cells():
    t, X, cls, y, mask, u = _traj_inputs((4, 4), (0.0, 1.0))
    sel = np.zeros(len(t), bool)
    sel[:4] = True
    rows = B.trajectory_rows(t, X, cls, y, mask, sel, u)
    assert len(rows) == 1 and rows[0]["n"] == 4
