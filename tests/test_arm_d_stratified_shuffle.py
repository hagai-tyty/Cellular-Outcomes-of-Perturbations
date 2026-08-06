"""ARM D — the STRATIFIED ΔAge shuffle.

Arm C permuted HFF's labels globally and destroyed two things at once: the between-timepoint
trajectory (rho(day, ΔAge) = -0.905) AND the within-timepoint cell-level pairing. So it proved the
labels carry exploitable structure without saying which kind.

Arm D permutes **within** each `(cell_line, time_h)` stratum. The property that makes it a valid
separator, and what these tests pin:

    * every stratum's MULTISET of labels is unchanged  -> the between-stratum trajectory survives
      EXACTLY, including each day's mean ΔAge;
    * within a stratum the pairing is destroyed;
    * no label ever crosses a stratum boundary -- that would silently re-introduce arm C.
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
    def __init__(self, strata: bool, seed: int = 0, shuffle=frozenset({"hff_sc"})):
        self.age_shuffle_datasets = shuffle
        self.age_shuffle_seed = seed
        self.age_shuffle_strata = strata
        self.deconfound = True
        self.primary_regime = "holdout"


def _aux(n, strata, start=0):
    return ChunkAux(
        cell_ids=[f"c{start + i}" for i in range(n)],
        cell_line=np.array(["HFF"] * n), is_control=np.zeros(n, bool),
        d_age_raw=np.zeros(n), cc=np.zeros(n),
        age_mask=np.ones(n, bool), deconfound_mask=np.ones(n, bool),
        shuffle_mask=np.ones(n, bool), stratum=np.asarray(strata, dtype=object))


# ======================================================================== #
# THE PROPERTY THAT MAKES ARM D A VALID SEPARATOR                          #
# ======================================================================== #
def test_each_stratum_keeps_its_own_labels():
    """No label crosses a stratum boundary. If one did, arm D would be arm C with extra steps."""
    n = 60
    strata = np.array([f"HFF|{(i % 3) * 48}" for i in range(n)], dtype=object)
    y = np.arange(n, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(strata=True), {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]
    for s in set(strata):
        m = strata == s
        assert sorted(out[m]) == sorted(y[m]), f"labels leaked out of stratum {s}"


def test_the_between_stratum_trajectory_survives_exactly():
    """The whole point: each stratum's MEAN ΔAge is untouched, so rho(day, ΔAge) is preserved.
    Arm C's global shuffle destroys this; arm D must not."""
    n = 90
    days = np.repeat([0, 48, 96], n // 3)
    strata = np.array([f"HFF|{d}" for d in days], dtype=object)
    # a real trajectory: later days are more negative, as HFF's ΔAge actually is
    y = np.concatenate([np.random.default_rng(1).normal(m, 2.0, n // 3)
                        for m in (-2.0, -12.0, -24.0)])
    out = bd._shuffle_age_labels(_Cfg(strata=True), {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]
    for s in sorted(set(strata)):
        m = strata == s
        assert out[m].mean() == pytest.approx(y[m].mean())
        assert out[m].std() == pytest.approx(y[m].std())


def test_the_global_shuffle_does_not_preserve_it():
    """The contrast that makes the previous test meaningful -- arm C breaks the trajectory."""
    n = 90
    days = np.repeat([0, 48, 96], n // 3)
    strata = np.array([f"HFF|{d}" for d in days], dtype=object)
    y = np.concatenate([np.full(n // 3, m, dtype=float) for m in (-2.0, -12.0, -24.0)])
    out = bd._shuffle_age_labels(_Cfg(strata=False), {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]
    moved = [out[strata == s].mean() != pytest.approx(y[strata == s].mean())
             for s in sorted(set(strata))]
    assert any(moved), "the global shuffle left every stratum mean intact -- it is not global"


def test_within_a_stratum_the_pairing_is_actually_destroyed():
    n = 60
    strata = np.array([f"HFF|{(i % 2) * 48}" for i in range(n)], dtype=object)
    y = np.arange(n, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(strata=True, seed=5),
                                 {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]
    assert int((out != y).sum()) > n // 4, "hardly anything moved"


def test_strata_are_respected_across_shards():
    """A stratum spans chunks, so grouping must be global-by-key, not per shard."""
    a = _aux(20, np.array(["HFF|0"] * 10 + ["HFF|48"] * 10, dtype=object))
    b = _aux(20, np.array(["HFF|0"] * 10 + ["HFF|48"] * 10, dtype=object), start=20)
    ya, yb = np.arange(20, dtype=float), np.arange(100, 120, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(strata=True, seed=2), {"s0": a, "s1": b},
                                 {"s0": ya.copy(), "s1": yb.copy()})
    day0 = np.concatenate([out["s0"][:10], out["s1"][:10]])
    truth0 = np.concatenate([ya[:10], yb[:10]])
    assert sorted(day0) == sorted(truth0)
    crossed = np.isin(out["s0"][:10], yb[:10]).sum()
    assert crossed > 0, "no value crossed shards within the same stratum -- grouping is per shard"


def test_a_singleton_stratum_is_left_alone_not_dropped():
    """One cell in a stratum cannot be permuted; it must keep its label, not lose it."""
    strata = np.array(["HFF|0"] * 5 + ["HFF|999"], dtype=object)   # last stratum has n=1
    y = np.arange(6, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(strata=True), {"s0": _aux(6, strata)}, {"s0": y.copy()})["s0"]
    assert out[5] == y[5]
    assert sorted(out[:5]) == sorted(y[:5])


def test_arm_c_is_unchanged_by_arm_ds_existence():
    """Regression: the unstratified path must behave exactly as before."""
    n = 40
    strata = np.array([f"HFF|{(i % 4) * 24}" for i in range(n)], dtype=object)
    y = np.arange(n, dtype=float)
    out = bd._shuffle_age_labels(_Cfg(strata=False, seed=3),
                                 {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]
    assert sorted(out) == sorted(y)
    assert not np.array_equal(out, y)


def test_it_is_deterministic_under_the_seed():
    n = 40
    strata = np.array([f"HFF|{(i % 3) * 24}" for i in range(n)], dtype=object)
    y = np.arange(n, dtype=float)

    def run(seed):
        return bd._shuffle_age_labels(_Cfg(strata=True, seed=seed),
                                      {"s0": _aux(n, strata)}, {"s0": y.copy()})["s0"]

    assert np.array_equal(run(7), run(7))
    assert not np.array_equal(run(7), run(8))


def test_stratum_defaults_to_one_global_group():
    """Pre-arm-D callers get arm C's behaviour, not a crash."""
    a = ChunkAux(cell_ids=["a", "b"], cell_line=np.array(["HFF"] * 2),
                 is_control=np.zeros(2, bool), d_age_raw=np.zeros(2), cc=np.zeros(2),
                 age_mask=np.ones(2, bool))
    assert set(a.stratum) == {"__all__"}


def test_a_sidecar_written_before_arm_d_still_loads(tmp_path):
    d = tmp_path / "_cc_cache"
    d.mkdir()
    n = 4
    np.savez(d / "s0.npz", cell_ids=np.array([f"c{i}" for i in range(n)]).astype("U"),
             cell_line=np.array(["HFF"] * n).astype("U"), is_control=np.zeros(n, bool),
             d_age_raw=np.zeros(n), cc=np.zeros(n), age_mask=np.ones(n, bool))

    class P:
        root = tmp_path
    out = bd._load_cc_sidecars(P())
    assert set(out["s0"].stratum) == {"__all__"}
