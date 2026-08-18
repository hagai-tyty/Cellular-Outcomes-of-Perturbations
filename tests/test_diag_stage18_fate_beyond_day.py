"""Stage 18 -- separating the fate head's biology from its clock.

`dose_time` is a MODEL INPUT and it encodes the timepoint. On the held-out donors fate is very
nearly a function of that timepoint, so every marginal metric (PR-AUC, ROC-AUC) is inflated by
information the model was handed rather than inferred.

The only honest question is asked WITHIN a timepoint, where `dose_time` is constant and cannot
help. These tests pin that machinery, above all the two ways it could lie: pooling strata across
donors (which would invent pairs that share an hour but not an experiment), and a permutation null
that shuffles globally instead of within strata (which would be far too easy to beat).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage18_fate_beyond_day.py"
spec = importlib.util.spec_from_file_location("s18", SRC)
s18 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s18)

RESULTS = ROOT / "results" / "diag_stage18_fate_beyond_day_results_s16.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="stage 18 has not been run")


@pytest.fixture(scope="module")
def rec():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


# ---- the stratified statistic ---------------------------------------------------------------- #
def test_only_pairs_inside_a_stratum_are_counted():
    """A safe cell at hour 167 and an unsafe cell at hour 1295 form no pair: separating those is
    exactly the clock-reading this test exists to exclude."""
    score = np.array([0.9, 0.1, 0.8, 0.2])
    safe = np.array([True, False, True, False])
    same = np.array(["a", "a", "b", "b"])
    diff = np.array(["a", "b", "c", "d"])
    assert s18.stratified_pairs(score, safe, same)[2] == 2
    assert s18.stratified_pairs(score, safe, diff)[2] == 0


def test_a_perfect_within_stratum_ranking_scores_one():
    score, safe = np.array([0.9, 0.1]), np.array([True, False])
    auc, n = s18.stratified_auc(score, safe, np.array(["a", "a"]))
    assert auc == 1.0 and n == 1


def test_an_inverted_ranking_scores_zero():
    auc, _ = s18.stratified_auc(np.array([0.1, 0.9]), np.array([True, False]),
                                np.array(["a", "a"]))
    assert auc == 0.0


def test_ties_count_as_half():
    auc, _ = s18.stratified_auc(np.array([0.5, 0.5]), np.array([True, False]),
                                np.array(["a", "a"]))
    assert auc == 0.5


def test_a_stratum_with_one_class_contributes_no_pairs():
    """The common case here: 63 of 70 timepoints carry a single class."""
    assert s18.stratified_pairs(np.array([0.9, 0.8]), np.array([True, True]),
                                np.array(["a", "a"]))[2] == 0


def test_no_pairs_yields_nan_rather_than_a_fabricated_score():
    auc, n = s18.stratified_auc(np.array([0.9]), np.array([True]), np.array(["a"]))
    assert n == 0 and np.isnan(auc)


# ---- the null, and the way it could have been too easy --------------------------------------- #
def test_the_permutation_shuffles_within_strata_not_globally():
    """A global shuffle would also destroy the between-timepoint structure, which the model is
    NOT being credited for -- making the null trivially beatable and the p-value meaningless."""
    src = SRC.read_text(encoding="utf-8")
    assert "s[idx] = rng.permutation(s[idx])" in src
    assert "preserving the stratum sizes" in src


def test_a_random_scorer_does_not_beat_the_within_stratum_null():
    rng = np.random.default_rng(0)
    n = 60
    stratum = np.repeat([f"s{i}" for i in range(10)], 6)
    safe = np.tile([True, True, True, False, False, False], 10)
    r = s18.perm_null(rng.normal(size=n), safe, stratum, n_perm=400, seed=1)
    assert 0.02 < r["p"] < 0.98, "a null scorer must not look significant"


def test_a_perfect_scorer_does_beat_it():
    stratum = np.repeat([f"s{i}" for i in range(8)], 4)
    safe = np.tile([True, True, False, False], 8)
    score = np.where(safe, 1.0, 0.0) + np.arange(len(safe)) * 1e-6
    r = s18.perm_null(score, safe, stratum, n_perm=400, seed=2)
    assert r["observed"] == pytest.approx(1.0)
    assert r["p"] < 0.01


def test_the_null_is_deterministic_for_a_seed():
    stratum = np.repeat(["a", "b"], 6)
    safe = np.tile([True, True, True, False, False, False], 2)
    sc = np.linspace(0, 1, 12)
    assert s18.perm_null(sc, safe, stratum, n_perm=200, seed=3) == \
        s18.perm_null(sc, safe, stratum, n_perm=200, seed=3)


# ---- strata must not merge across donors ----------------------------------------------------- #
def test_pooled_strata_are_namespaced_by_donor():
    """The same hour in two donors is not one stratum -- they are different experiments, and
    merging them would manufacture pairs the design never contained."""
    src = SRC.read_text(encoding="utf-8")
    assert 'f"{d}:{v}" for v in stratum' in src
    assert "must not merge across donors" in src


# ---- the recorded result --------------------------------------------------------------------- #
@has_results
def test_fate_is_almost_a_function_of_timepoint(rec):
    """THE structural fact. Only 7 of 70 timepoints carry more than one class, and on three folds
    the number is zero -- so on those folds a lookup table on the hour is unbeatable."""
    p = rec["pooled"]
    assert p["n_strata"] == 70 and p["n_mixed_strata"] == 7
    zero = [d for d, f in rec["folds"].items() if f["n_mixed_strata"] == 0]
    assert set(zero) == {"N2", "O1", "O2"}


@has_results
def test_time_alone_is_perfect_on_two_folds(rec):
    """O1 and O2 reach PR-AUC 1.000 from the timepoint alone. Any model metric on those folds is
    measuring the calendar."""
    for d in ("O1", "O2"):
        assert rec["folds"][d]["baselines"]["time_only"]["prauc"] == pytest.approx(1.0, abs=1e-9)


@has_results
def test_the_within_timepoint_signal_is_real_but_rests_on_twelve_pairs(rec):
    """The finding, with its own power attached. Significant under a within-stratum permutation
    null -- and the entire evidence base is 12 pairs from 7 strata across 3 donors."""
    sm = rec["pooled"]["stratified_model"]
    assert sm["n_pairs"] == 12
    assert sm["observed"] == pytest.approx(0.917, abs=0.005)
    assert sm["p"] < 0.05


@has_results
def test_the_marginal_number_overstates_the_stratified_one(rec):
    """0.931 marginal against 0.917 stratified looks close, but the marginal is computed over 70
    strata and the stratified over 7 -- they are not the same claim on the same evidence."""
    p = rec["pooled"]
    assert p["model_prauc"] == pytest.approx(0.931, abs=0.005)
    assert p["n_mixed_strata"] < 0.15 * p["n_strata"]


@has_results
def test_y1_is_the_fold_carrying_the_evidence(rec):
    """5 of its 11 timepoints are mixed and it supplies 7 of the 12 pairs. It is also the fold
    where the marginal metric is WORST (0.796) -- the two facts are the same fact."""
    y1 = rec["folds"]["Y1"]
    assert y1["n_mixed_strata"] == 5
    assert y1["stratified_pairs"] == 7
    assert y1["model_prauc"] < 0.85


@has_results
def test_n2_contributes_nothing_because_it_has_no_unsafe_cells(rec):
    assert rec["folds"]["N2"]["model_prauc"] is None
    assert rec["folds"]["N2"]["stratified_pairs"] == 0
