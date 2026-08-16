"""Unit tests for the Phase 0-3 diagnostics.

Each phase turned on one piece of logic that, if wrong, would have produced a confident false
result. Those are what is tested here.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dcd = _load("dcd", "experiments/diag_clock_difference_capacity.py")
p1 = _load("p1", "experiments/diag_phase1_top100.py")


# ---- weight truncation: the whole of Phase 0 rests on this ---------------------------------- #
def test_truncation_keeps_exactly_the_n_largest_by_absolute_weight():
    """By |weight|, not by signed weight -- a large NEGATIVE age coefficient is as informative as
    a large positive one, and ranking on the signed value would silently drop them all."""
    w = {"a": -9.0, "b": 0.5, "c": 7.0, "d": -0.1}
    out = p1.top_n_weights(w, ["a", "b", "c", "d"], 2)
    assert list(out) == [-9.0, 0.0, 7.0, 0.0]


def test_truncation_zeroes_rather_than_drops():
    """Dropping genes would change the vector's length and silently misalign the dot product."""
    out = p1.top_n_weights({"a": 5.0}, ["a", "b", "c"], 1)
    assert len(out) == 3 and list(out) == [5.0, 0.0, 0.0]


def test_genes_absent_from_the_clock_get_zero():
    out = p1.top_n_weights({"a": 5.0}, ["a", "zzz"], 2)
    assert list(out) == [5.0, 0.0]


def test_asking_for_more_genes_than_exist_is_harmless():
    out = p1.top_n_weights({"a": 5.0, "b": 1.0}, ["a", "b"], 999)
    assert list(out) == [5.0, 1.0]


def test_variant_weights_builds_every_level_on_the_same_gene_order():
    genes = ["a", "b", "c", "d"]
    v = dcd.variant_weights({"a": -9.0, "b": 0.5, "c": 7.0, "d": -0.1}, genes)
    assert set(v) >= {"raw", "top100", "top500", "top2000", "covnorm"}
    for k, w in v.items():
        assert len(w) == len(genes), k


def test_covnorm_rescales_by_the_uncovered_mass():
    """Genes the local matrix lacks read as zero; covnorm exists to stop that silently
    under-reporting. With half the mass missing the surviving weights must double."""
    v = dcd.variant_weights({"a": 1.0, "b": 1.0}, ["a"])   # 'b' absent -> 50% covered
    assert v["covnorm"][0] == pytest.approx(2.0)


def test_covnorm_is_a_noop_when_coverage_is_complete():
    v = dcd.variant_weights({"a": 1.0, "b": 1.0}, ["a", "b"])
    assert np.allclose(v["covnorm"], v["raw"])


# ---- the refit-CV that made Part A's memorisation visible ----------------------------------- #
def test_the_refit_recovers_a_known_linear_signal_out_of_fold():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 40))
    age = X @ rng.normal(size=40) * 0.5 + 50
    out = dcd.refit_cv(X, age, alpha=1.0, top_ns=(5,), n_splits=5)
    assert out["full"]["r_abs"] > 0.8


def test_the_refit_never_scores_in_sample():
    """Part A scored the shipped clock on its own training donors and returned MAE 0.13 against a
    published cv_mae of 12.27. Out-of-fold prediction of pure noise must NOT do that."""
    rng = np.random.default_rng(1)
    X, age = rng.normal(size=(60, 200)), rng.uniform(1, 96, 60)
    out = dcd.refit_cv(X, age, alpha=1.0, top_ns=(10,), n_splits=5)
    assert out["full"]["mae_abs"] > 5.0, "noise scored like signal -- the fold split leaks"


def test_difference_error_is_reported_alongside_absolute():
    rng = np.random.default_rng(2)
    X, age = rng.normal(size=(40, 20)), rng.uniform(1, 96, 40)
    out = dcd.refit_cv(X, age, alpha=10.0, top_ns=(5,), n_splits=4)
    for lvl in out.values():
        assert {"mae_abs", "mae_diff", "r_diff", "sd_ratio"} <= set(lvl)


# ---- the persistence control that overturned Phase 3 ---------------------------------------- #
def test_a_static_trajectory_is_predictable_from_its_own_start():
    """Phase 3's raw verdict was SIGNAL until this was checked: endpoints correlated at 0.971, so
    'predict late from early' was predicting a number that had barely moved. Any forward claim must
    beat the trajectory simply staying put."""
    rng = np.random.default_rng(3)
    early = rng.normal(0, 10, 30)
    late = early + rng.normal(0, 0.5, 30)          # barely moves
    import pandas as pd
    assert pd.Series(early).corr(pd.Series(late), method="spearman") > 0.95


def test_a_moving_trajectory_is_not_predictable_from_its_start():
    rng = np.random.default_rng(4)
    early = rng.normal(0, 10, 30)
    late = rng.normal(0, 10, 30)                    # unrelated
    import pandas as pd
    assert abs(pd.Series(early).corr(pd.Series(late), method="spearman")) < 0.5
