"""STAGE 1.5.3 C-I — `y_age` must NOT depend on the training-label policy.

Step 6's first run was confounded because `_deconfound_train_only` fitted the cell-cycle
deconfounder, and re-centred on controls, using ``age_mask`` -- which `AGE_MASKED_DATASETS` changes.
So masking HFF did not merely remove labels, it **redefined the target variable** for every cell,
including the held-out evaluation cells. `dage_mae_ridge` regressed +9.21 yr even though ridge never
touches the trained age head, and the deconfounder slope went -3.93 -> -24.20 on N2 with its
intercept flipping sign. Full write-up: `results/STEP6_REPORT.md` section 3, confound C-I.

The fix separates two questions that had been conflated:

    * ``deconfound_mask`` -- is this cell's ΔAge COMPUTABLE?   -> decides the VALUE of `y_age`
    * ``age_mask``        -- may the age head TRAIN on it?      -> decides which cells the loss uses

These tests pin the invariant that makes step 6 a one-change experiment: **two runs differing only
in `AGE_MASKED_DATASETS` must produce bit-identical `y_age`.** If this ever fails, the arms have
different target variables again and the comparison is void whatever it reports.
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


class _Paths:
    def __init__(self, root: Path) -> None:
        self.root = root

    def shard_file(self, sid: str) -> Path:
        return self.root / f"{sid}.parquet"


class _Cfg:
    deconfound = True
    primary_regime = "holdout"


class _Row:
    def __init__(self, age_mask: bool) -> None:
        self.age_mask = age_mask


def _aux(n: int, *, trainable: np.ndarray, computable: np.ndarray, seed: int = 0) -> ChunkAux:
    rng = np.random.default_rng(seed)
    cc = rng.normal(size=n)
    return ChunkAux(
        cell_ids=[f"c{i}" for i in range(n)],
        cell_line=np.array(["L1"] * (n // 2) + ["L2"] * (n - n // 2)),
        is_control=np.array([i % 4 == 0 for i in range(n)]),
        # a real cell-cycle dependence, so the fitted coefficient is not trivially zero
        d_age_raw=3.0 * cc + rng.normal(scale=0.5, size=n) + 10.0,
        cc=cc,
        age_mask=trainable,
        deconfound_mask=computable,
    )


def _run(monkeypatch, tmp_path, *, trainable, computable, n=40):
    """Run the real `_deconfound_train_only` and capture what it would write."""
    written: dict[str, np.ndarray] = {}
    monkeypatch.setattr(bd.io, "rewrite_shard_yage",
                        lambda path, y: written.__setitem__(Path(path).stem, np.asarray(y)))
    aux = _aux(n, trainable=trainable, computable=computable)
    splits = {"holdout": {f"c{i}": "train" for i in range(n)}}
    rows = [_Row(bool(x)) for x in trainable]
    coef = bd._deconfound_train_only(_Cfg(), _Paths(tmp_path), rows, splits, {"s0": aux})
    return coef, written["s0"]


# ======================================================================== #
# THE INVARIANT                                                            #
# ======================================================================== #
def test_y_age_is_identical_whether_or_not_a_dataset_is_masked(monkeypatch, tmp_path):
    """The one that makes step 6 a one-change experiment.

    Arm A trains on everything; arm B withholds the first 30 of 40 cells (a stand-in for HFF).
    Their ΔAge is still computable, so `deconfound_mask` is unchanged -- and `y_age` must come out
    bit-identical, NaNs in the same places.
    """
    n = 40
    computable = np.ones(n, dtype=bool)
    arm_a = np.ones(n, dtype=bool)
    arm_b = np.ones(n, dtype=bool)
    arm_b[:30] = False                      # withheld from TRAINING only

    coef_a, y_a = _run(monkeypatch, tmp_path, trainable=arm_a, computable=computable, n=n)
    coef_b, y_b = _run(monkeypatch, tmp_path, trainable=arm_b, computable=computable, n=n)

    assert coef_a == pytest.approx(coef_b), (
        f"the deconfounder moved with the training policy: {coef_a} vs {coef_b}")
    assert np.array_equal(np.isnan(y_a), np.isnan(y_b)), "NaN pattern differs between arms"
    ok = ~np.isnan(y_a)
    assert np.array_equal(y_a[ok], y_b[ok]), "y_age is not bit-identical across arms"


def test_the_pre_fix_behaviour_would_have_failed_that_invariant(monkeypatch, tmp_path):
    """Mutation check: reproduce the OLD behaviour by letting `deconfound_mask` follow the
    training mask, and confirm the invariant above genuinely catches it. Without this, the test
    could be passing for a reason unrelated to the fix."""
    n = 40
    arm_a = np.ones(n, dtype=bool)
    arm_b = np.ones(n, dtype=bool)
    arm_b[:30] = False

    _, y_a = _run(monkeypatch, tmp_path, trainable=arm_a, computable=arm_a, n=n)
    _, y_b = _run(monkeypatch, tmp_path, trainable=arm_b, computable=arm_b, n=n)

    both = ~np.isnan(y_a) & ~np.isnan(y_b)
    assert both.any()
    assert not np.allclose(y_a[both], y_b[both]), (
        "the old behaviour did NOT change y_age here, so the invariant test proves nothing")


def test_the_coefficient_is_fitted_on_computable_cells_not_trainable_ones(monkeypatch, tmp_path):
    """C-I's core claim, stated directly: withholding labels must not shrink the fit set."""
    n = 40
    computable = np.ones(n, dtype=bool)
    withheld = np.ones(n, dtype=bool)
    withheld[:38] = False                   # only 2 trainable labels left

    coef_full, _ = _run(monkeypatch, tmp_path, trainable=computable, computable=computable, n=n)
    coef_masked, _ = _run(monkeypatch, tmp_path, trainable=withheld, computable=computable, n=n)
    assert coef_full == pytest.approx(coef_masked)


def test_a_genuinely_uncomputable_cell_still_drops_out(monkeypatch, tmp_path):
    """The fix must not readmit cells excluded for REAL reasons (cancer source, out-of-clock-range).
    Only the `dataset_policy` exclusion is a training-only decision."""
    n = 40
    computable = np.ones(n, dtype=bool)
    computable[:10] = False                 # e.g. cancer_source / donor_out_of_clock_range
    trainable = computable.copy()

    _, y = _run(monkeypatch, tmp_path, trainable=trainable, computable=computable, n=n)
    assert np.isnan(y[:10]).all(), "cells with no computable ΔAge must stay NaN"
    assert not np.isnan(y[10:]).any()


# ======================================================================== #
# backward compatibility                                                   #
# ======================================================================== #
def test_chunkaux_defaults_deconfound_mask_to_age_mask():
    """Every pre-C-I construction site keeps its old meaning rather than crashing."""
    m = np.array([True, False, True])
    a = ChunkAux(cell_ids=["a", "b", "c"], cell_line=np.array(["L"] * 3),
                 is_control=np.zeros(3, bool), d_age_raw=np.zeros(3), cc=np.zeros(3),
                 age_mask=m)
    assert np.array_equal(a.deconfound_mask, m)


def test_a_sidecar_written_before_ci_still_loads(tmp_path):
    """A RESUMED build must not crash on an npz that predates the new field -- the
    `rewrite_shard_yage` lesson, applied before it can bite."""
    d = tmp_path / "_cc_cache"
    d.mkdir()
    n = 5
    np.savez(d / "s0.npz", cell_ids=np.array([f"c{i}" for i in range(n)]).astype("U"),
             cell_line=np.array(["L"] * n).astype("U"),
             is_control=np.zeros(n, bool), d_age_raw=np.zeros(n), cc=np.zeros(n),
             age_mask=np.ones(n, bool))          # NO deconfound_mask

    class P:
        root = tmp_path
    out = bd._load_cc_sidecars(P())
    assert np.array_equal(out["s0"].deconfound_mask, out["s0"].age_mask)
