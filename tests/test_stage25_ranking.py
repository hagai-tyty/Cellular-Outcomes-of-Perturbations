"""Contracts for Stage 25 — the preregistered ranking test.

This is the one load-bearing new result in Generation 1, and every degree of freedom in it was
frozen before any of its numbers existed. These contracts exist because a preregistration is only
worth something if something checks that the code implements what was registered — the metric, the
population, the weighting, the comparator, the permutation count, and the verdict logic.

The sharpest ones are negative: that an incomplete null is refused rather than silently becoming a
smaller null, and that Stage 25 cannot reach a SUPPORTED verdict without all six criteria.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "stage25"
SRC = ROOT / "experiments" / "run_stage25_ranking.py"
PLAN = ROOT / "plans" / "(newer)practical plans" / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"

A_JSON = OUT / "stage25a_observed.json"
SMOKE = OUT / "stage25_smoke.json"
VERDICT = OUT / "stage25_verdict.json"

ran_a = pytest.mark.skipif(not A_JSON.exists(), reason="25A has not been run")
ran_smoke = pytest.mark.skipif(not SMOKE.exists(), reason="the smoke stage has not been run")
ran_c = pytest.mark.skipif(not VERDICT.exists(), reason="25C has not been run")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def mod():
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_stage25_ranking as S25
    return S25


# ============================================================================================== #
# The metric is the plan's formula
# ============================================================================================== #
def test_within_clone_auc_implements_the_plan_formula(mod):
    """§8.5: pairwise over positive/zero pairs, ties scoring exactly 0.5."""
    f = mod.within_clone_auc
    assert f(np.array([0.9, 0.8, 0.1, 0.2]), np.array([1, 1, 0, 0])) == 1.0
    assert f(np.array([0.1, 0.2, 0.9, 0.8]), np.array([1, 1, 0, 0])) == 0.0
    assert f(np.array([0.5, 0.5]), np.array([1, 0])) == 0.5
    # one positive above one zero and below another -> exactly half the pairs
    assert f(np.array([0.5, 0.9, 0.1]), np.array([1, 0, 0])) == 0.5
    # undefined for a clone with no contrast, which is why such clones are excluded
    assert np.isnan(f(np.array([0.3, 0.4]), np.array([1, 1])))


def test_equal_clone_weighting_not_micro_averaging(mod):
    """§8.5: a clone contributing many pairs must not outweigh one contributing few."""
    import pandas as pd
    df = pd.DataFrame({
        # clone A: 1 positive, 1 zero -> AUC 1.0 from a single pair
        # clone B: 1 positive, 3 zeros -> AUC 0.0 from three pairs
        "clone_id": ["A"] * 2 + ["B"] * 4,
        "treatment": ["Acid", "Cisplatin"] + ["Acid", "Cisplatin", "CoCl2", "Dabrafenib"],
        "y": [1, 0] + [1, 0, 0, 0],
        "s": [0.9, 0.1] + [0.1, 0.9, 0.8, 0.7],
    })
    r, per = mod.rank_score(df, "s")
    assert per["A"] == 1.0 and per["B"] == 0.0
    assert r == 0.5, "equal-clone weighting gives 0.5; micro-averaging would give 0.25"


def test_the_frozen_parameters_are_the_plan_s(mod):
    assert mod.N_PERM == 1000 and mod.SEED_PERM == 23523
    assert mod.N_BOOT == 2000 and mod.SEED_BOOT == 23501
    assert mod.EXPECTED_ELIGIBLE == 892
    assert len(mod.CONDITIONS) == 6


def test_the_module_refuses_an_incomplete_null():
    """A missing draw is an integrity stop, never a smaller null."""
    src = SRC.read_text(encoding="utf-8")
    assert "the null is incomplete" in src
    assert "INTEGRITY STOP, not a smaller null" in src


def test_the_module_verifies_its_inputs_before_reading_a_number():
    src = SRC.read_text(encoding="utf-8")
    assert "the out-of-fold table has changed since Stage 24 hashed it" in src
    assert "the frozen plan has moved since Stage 24 consumed it" in src
    assert "INPUT-INTEGRITY STOP" in src


def test_the_null_is_a_full_refit_not_a_label_shuffle():
    src = SRC.read_text(encoding="utf-8")
    assert "never across it" in src or "never across the boundary" in src.replace("\n", " ")
    assert "Observed-data hyperparameters are NOT reused" in src.replace("\n", " ")
    # the inner selection must actually run inside the null draw
    assert "GroupKFold" in src and "_fit_logistic" in src


# ============================================================================================== #
# 25A / smoke
# ============================================================================================== #
@ran_a
def test_the_observed_statistic_came_from_the_frozen_table():
    a = _json(A_JSON)
    assert a["eligible_clones"] == 892
    assert a["primary_comparator"] == "W4" and a["endpoint"] == "C1"
    assert a["metric"].startswith("equal-clone-weighted within-clone AUROC")
    h = _json(RESULTS / "stage24_handoff_to_stage25.json")
    assert a["oof_table_sha256"] == h["frozen_oof_table"]["sha256"]


@ran_smoke
def test_the_smoke_test_proved_sharding_cannot_change_the_answer():
    s = _json(SMOKE)
    by = {c["check"]: c["pass"] for c in s["checks"]}
    assert by["2 shards run OUT OF ORDER are bit-identical to sequential"] is True
    assert by["draw b is a function of b alone (repeatable)"] is True
    assert by["each shard writes its OWN file"] is True
    assert by["within-clone AUC matches the plan formula"] is True
    assert s["all_passed"] is True
    assert s["measured_seconds_per_draw"] > 0


# ============================================================================================== #
# 25C — the verdict
# ============================================================================================== #
@ran_c
def test_the_verdict_is_one_of_exactly_two_values():
    v = _json(VERDICT)
    assert v["verdict"] in ("STAGE_25_RANKING_SUPPORTED", "STAGE_25_RANKING_NOT_SUPPORTED")
    assert "GEN1_MANDATORY_SHIP" in v["gen1_next"]
    assert "terminal scientific result" in v["terminal"]


@ran_c
def test_supported_requires_all_six_criteria():
    v = _json(VERDICT)
    assert set(v["criteria"]) == {
        "1_delta_RANK_gt_0", "2_bootstrap_lower_endpoint_gt_0", "3_observed_gt_null_p95",
        "4_p_perm_le_0_05", "5_delta_TOP1_ge_0", "6_integrity_leakage_determinism"}
    if v["verdict"] == "STAGE_25_RANKING_SUPPORTED":
        assert all(v["criteria"].values())
        assert v["failing_criteria"] == []
    else:
        assert v["failing_criteria"], "a NOT_SUPPORTED verdict must name what failed"


@ran_c
def test_the_permutation_contract_is_the_frozen_one():
    p = _json(VERDICT)["permutation"]
    assert p["n_perm"] == 1000 and p["base_seed"] == 23523
    assert p["early_stopping"] is False
    assert "full refit" in p["type"]
    assert "hyperparameters never reused" in p["type"]
    # p_perm must equal the finite-sample formula exactly
    assert abs(p["p_perm"] - (1 + p["n_null_ge_observed"]) / 1001) < 1e-12


@ran_c
def test_delta_top1_can_only_withhold_support_never_grant_it():
    d = _json(VERDICT)["delta_TOP1"]
    assert "directional-consistency check, not a significance test" in d["role"]
    assert "withhold support, never grant it" in d["role"]


@ran_c
def test_the_secondary_comparator_cannot_rescue_the_primary():
    s = _json(VERDICT)["secondary"]
    assert "cannot rescue" in s["note"]


@ran_c
def test_the_limitations_survive_into_the_verdict():
    lim = " ".join(_json(VERDICT)["standing_limitations"])
    assert "no independent biological replication" in lim
    assert "3.45x the state contribution" in lim
    assert "Cisplatin is negligible" in lim and "Doxorubicin is negative" in lim
    assert "not death, sensitivity or clinical response" in lim
    assert "gate 18.3 FAILED" in lim


@ran_c
def test_the_descriptives_are_reported_without_rescue_power():
    d = _json(VERDICT)["descriptives"]
    for k in ("by_outer_fold", "by_pretreatment_depth_bin",
              "pairwise_condition_ranking_accuracy_W5", "score_tie_rate_W5",
              "excluded_all_zero_clones", "excluded_all_positive_clones"):
        assert k in d, k
    assert d["excluded_all_zero_clones"] == 472
    assert d["excluded_all_positive_clones"] == 37
    assert len(d["by_outer_fold"]) == 5
