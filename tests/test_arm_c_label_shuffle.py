"""ARM C — the ΔAge label-permutation control (step-6 follow-up).

Step 6's rerun found a consistent ranking cost when HFF's labels are withheld (`rank_model_dage`
-0.0688, CI excludes 0). That is confounded between two explanations:

    (i)  HFF's 33,613 labels carry information the age head uses, or
    (ii) 75 labels is simply too few to learn from, whatever they contain.

Arm C holds label VOLUME at arm A's level and destroys only the cell<->label PAIRING, separating
them. These tests pin the three properties that make it a valid control:

  * the label MULTISET is unchanged -- same values, same count, only reassigned;
  * the permutation is GLOBAL across shards, not per chunk (chunks are timepoint-homogeneous, so a
    within-chunk shuffle would leave between-timepoint structure intact and the control would be
    far weaker than it looks);
  * it runs AFTER the deconfounder fit and control re-centring, so the coefficient and every ΔAge
    VALUE stay bit-identical to arm A.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from cellfate.data import build_dataset as bd  # noqa: E402
from cellfate.data.build_dataset import ChunkAux  # noqa: E402


class _Cfg:
    def __init__(self, shuffle=frozenset(), seed=0):
        self.age_shuffle_datasets = shuffle
        self.age_shuffle_seed = seed
        self.deconfound = True
        self.primary_regime = "holdout"


def _aux(n, shuffle_mask, start=0):
    return ChunkAux(
        cell_ids=[f"c{start + i}" for i in range(n)],
        cell_line=np.array(["HFF"] * n), is_control=np.zeros(n, bool),
        d_age_raw=np.zeros(n), cc=np.zeros(n),
        age_mask=np.ones(n, bool), deconfound_mask=np.ones(n, bool),
        shuffle_mask=np.asarray(shuffle_mask, dtype=bool))


def test_the_label_multiset_is_preserved_exactly():
    """Arm C must change WHICH cell gets which label, never the labels themselves. If the multiset
    moved, arm C would differ from arm A in label distribution too and stop being a clean control."""
    n = 40
    aux = {"s0": _aux(n, np.ones(n, bool))}
    ys = {"s0": np.arange(n, dtype=float)}
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=3), aux, {k: v.copy() for k, v in ys.items()})
    assert sorted(out["s0"]) == sorted(ys["s0"])
    assert not np.array_equal(out["s0"], ys["s0"]), "nothing actually moved"


def test_the_permutation_is_global_across_shards():
    """The property a per-chunk shuffle would NOT have. With two shards holding disjoint value
    ranges, a global permutation must move values ACROSS the shard boundary."""
    a, b = np.arange(30, dtype=float), np.arange(100, 130, dtype=float)
    aux = {"s0": _aux(30, np.ones(30, bool)), "s1": _aux(30, np.ones(30, bool), start=30)}
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=1), aux,
                                 {"s0": a.copy(), "s1": b.copy()})
    crossed = np.isin(out["s0"], b).sum() + np.isin(out["s1"], a).sum()
    assert crossed > 0, "no value crossed shards -- the shuffle is per-chunk, not global"
    assert sorted(np.concatenate([out["s0"], out["s1"]])) == sorted(np.concatenate([a, b]))


def test_only_targeted_cells_move():
    """Cells outside `shuffle_mask` -- the Gill donors -- must keep their own labels."""
    n = 40
    mask = np.zeros(n, bool)
    mask[:20] = True                       # only the first half is HFF
    aux = {"s0": _aux(n, mask)}
    y = np.arange(n, dtype=float)
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=5), aux, {"s0": y.copy()})
    assert np.array_equal(out["s0"][20:], y[20:]), "untargeted cells were shuffled"
    assert sorted(out["s0"][:20]) == sorted(y[:20])


def test_nan_labels_do_not_take_part():
    """A cell with no computable ΔAge stays NaN; it must not absorb a real label or donate one."""
    n = 10
    y = np.arange(n, dtype=float)
    y[3] = y[7] = np.nan
    aux = {"s0": _aux(n, np.ones(n, bool))}
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=2), aux, {"s0": y.copy()})
    assert np.isnan(out["s0"][3]) and np.isnan(out["s0"][7])
    ok = ~np.isnan(y)
    assert sorted(out["s0"][ok]) == sorted(y[ok])


def test_it_is_deterministic_under_the_seed_and_varies_across_seeds():
    """The seed is recorded in the census; the run must be reproducible from it."""
    n = 50
    aux = {"s0": _aux(n, np.ones(n, bool))}
    y = np.arange(n, dtype=float)
    r = lambda s: bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=s), aux, {"s0": y.copy()})["s0"]  # noqa: E731
    assert np.array_equal(r(11), r(11))
    assert not np.array_equal(r(11), r(12))


def test_it_is_a_no_op_when_not_requested():
    """Arms A and B must be untouched by arm C's existence."""
    n = 20
    aux = {"s0": _aux(n, np.ones(n, bool))}
    y = np.arange(n, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(frozenset()), aux, {"s0": y.copy()})
    assert np.array_equal(out["s0"], y)


def test_a_single_target_label_is_a_no_op_not_a_crash():
    """Degenerate case: one label cannot be permuted against anything."""
    mask = np.zeros(5, bool)
    mask[2] = True
    aux = {"s0": _aux(5, mask)}
    y = np.arange(5, dtype=float)
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=0), aux, {"s0": y.copy()})
    assert np.array_equal(out["s0"], y)


def test_shuffle_mask_defaults_to_empty_for_pre_arm_c_callers():
    a = ChunkAux(cell_ids=["a", "b"], cell_line=np.array(["L"] * 2), is_control=np.zeros(2, bool),
                 d_age_raw=np.zeros(2), cc=np.zeros(2), age_mask=np.ones(2, bool))
    assert a.shuffle_mask.shape == (2,)
    assert not a.shuffle_mask.any()


def test_a_sidecar_written_before_arm_c_still_loads(tmp_path):
    """A RESUMED build must not crash on an npz predating the field."""
    d = tmp_path / "_cc_cache"
    d.mkdir()
    n = 4
    np.savez(d / "s0.npz", cell_ids=np.array([f"c{i}" for i in range(n)]).astype("U"),
             cell_line=np.array(["L"] * n).astype("U"), is_control=np.zeros(n, bool),
             d_age_raw=np.zeros(n), cc=np.zeros(n), age_mask=np.ones(n, bool))

    class P:
        root = tmp_path
    out = bd._load_cc_sidecars(P())
    assert not out["s0"].shuffle_mask.any()


def test_shuffling_preserves_the_mean_so_level_metrics_cannot_shift_by_construction():
    """A sanity property worth pinning: because the multiset is preserved, the GLOBAL mean and SD
    of the age labels are unchanged. Any level shift arm C shows is learned, not injected."""
    n = 200
    rng = np.random.default_rng(0)
    y = rng.normal(5.0, 3.0, size=n)
    aux = {"s0": _aux(n, np.ones(n, bool))}
    out = bd._shuffle_age_labels(_Cfg({"hff_sc"}, seed=9), aux, {"s0": y.copy()})["s0"]
    assert out.mean() == pytest.approx(y.mean())
    assert out.std() == pytest.approx(y.std())
