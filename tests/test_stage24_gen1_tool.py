"""Contracts for Stage 24 — Gen-1 Role-B predictor engineering.

Stage 24 is bounded by the frozen Stage-23.5 plan, and its whole value rests on two things being
true: that the freeze it consumed is intact, and that the reproduction it certifies actually
reproduced something. These contracts pin both.

The sharpest one is `test_the_gate_is_self_consistent`. The first 24B run reported a file as
BYTE_IDENTICAL *and* failing a sub-gate, because the gate hardcoded one row count for two endpoints
with different row counts. Byte-identity happened to carry the verdict, so nothing was wrong with
the science -- but the same defect would have forced a spurious INPUT_INTEGRITY_STOP had the
reproduction been merely tolerance-clean. The contract now refuses that combination outright.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "stage24"

A_JSON = OUT / "stage24a_engineering_plan.json"
B_JSON = OUT / "stage24b_reproduction.json"
PROTOCOL = RESULTS / "stage23_5_protocol.json"
HANDOFF = RESULTS / "stage23_5_handoff_to_stage24.json"

ran_a = pytest.mark.skipif(not A_JSON.exists(), reason="24A has not been run")
ran_b = pytest.mark.skipif(not B_JSON.exists(), reason="24B has not been run")
has_freeze = pytest.mark.skipif(not PROTOCOL.exists(), reason="Stage 23.5 is not frozen")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# ============================================================================================== #
# 24A — the freeze it consumed
# ============================================================================================== #
@ran_a
def test_24a_verified_the_freeze_before_any_engineering():
    a = _json(A_JSON)
    assert a["all_checks_pass"] is True
    c = a["checks"]
    for k in ("plan_digest_matches_protocol", "plan_digest_matches_handoff",
              "plan_status_frozen", "stage_24_open", "ranking_protocol_frozen",
              "audit_fully_passed", "compute_budget_accepted", "no_source_artifact_drift"):
        assert c[k] is True, k


@ran_a
def test_24a_confirms_stage_24_may_not_see_the_ranking_metric():
    """Stage 24 generates Stage 25's inputs. It must not be able to look at the answer."""
    a = _json(A_JSON)
    assert a["checks"]["ranking_metric_not_inspected"] is True
    assert a["checks"]["no_ranking_artifact_exists"] is True
    assert a["ranking_artifacts_found"] == []
    assert "inspect the Stage-25 ranking metric" in a["engineering_plan"]["may_not"]


@ran_a
def test_24a_forbids_the_things_that_would_make_stage_24_unbounded():
    may_not = " ".join(_json(A_JSON)["engineering_plan"]["may_not"])
    assert "replace W5 because another architecture scores better" in may_not
    assert "add a dataset" in may_not
    assert "change any Stage-22 or Stage-23 frozen quantity" in may_not


@ran_a
@has_freeze
def test_24a_pins_the_same_digest_the_freeze_recorded():
    a, p, h = _json(A_JSON), _json(PROTOCOL), _json(HANDOFF)
    d = a["plan"]["canonical_lf_sha256"]
    assert d == p["plan_canonical_lf_sha256"] == h["plan_canonical_lf_sha256"]
    assert p["plan_status"] == "FROZEN"


# ============================================================================================== #
# 24B — the reproduction
# ============================================================================================== #
@ran_b
def test_the_reproduction_verdict_is_one_of_the_three_frozen_values():
    b = _json(B_JSON)
    assert b["reproduction_verdict"] in ("BYTE_IDENTICAL", "TOLERANCE_DECLARED",
                                         "INPUT_INTEGRITY_STOP")
    for g in b["gates"].values():
        assert g["verdict"] in ("BYTE_IDENTICAL", "TOLERANCE_DECLARED", "INPUT_INTEGRITY_STOP")


@ran_b
def test_the_gate_is_self_consistent():
    """A byte-identical file that fails a sub-gate means the GATE is broken, not the reproduction.

    This is the contract for the defect found on the first 24B run: EXPECTED_ROWS was a single
    constant, but C1 scores all 8,406 clone x condition rows while C2 is defined only on the 2,256
    detected rows. The module now raises on this combination; this pins it in the artifact too.
    """
    for label, g in _json(B_JSON)["gates"].items():
        assert g["gate_self_consistent"] is True, label
        if g["byte_identical"]:
            assert g["R1_shape_and_key"]["pass"] is True, label
            assert g["R2_every_score"]["pass"] is True, label
            assert g["R3_within_clone_ordering"]["pass"] is True, label


@ran_b
def test_row_counts_are_per_endpoint_not_global():
    """C1 scores every clone x condition row; C2 only the detected ones. 8406 != 2256."""
    g = _json(B_JSON)["gates"]
    c1 = g["C1_W0toW4"]["R1_shape_and_key"]
    c2 = g["C2_W0toW4"]["R1_shape_and_key"]
    assert c1["endpoint"] == "C1" and c1["expected_rows_for_endpoint"] == 8406
    assert c2["endpoint"] == "C2" and c2["expected_rows_for_endpoint"] == 2256
    assert c1["rows_frozen"] == 8406 and c2["rows_frozen"] == 2256
    for r in (c1, c2):
        assert r["frozen_matches_documented_count"] is True
        assert r["repro_matches_frozen_row_count"] is True


@ran_b
def test_r3_checks_within_clone_ordering_for_every_clone_and_both_ranking_models():
    """The Stage-25 ranking test is a function of within-clone orderings and nothing else."""
    r3 = _json(B_JSON)["gates"]["C1_W5"]["R3_within_clone_ordering"]
    assert set(r3["models_checked"]) == {"pred_W4", "pred_W5"}
    for m in ("pred_W4", "pred_W5"):
        assert r3["per_model"][m]["clones_checked"] == 1401
        assert r3["per_model"][m]["clones_with_changed_ordering"] == 0


@ran_b
def test_r2_checks_every_prediction_cell_not_an_aggregate():
    r2 = _json(B_JSON)["gates"]["C1_W5"]["R2_every_score"]
    assert r2["tolerance"] == 1e-12
    assert {"pred_W1", "pred_W4", "pred_W5"} <= set(r2["columns_checked"])
    for col, v in r2["per_column"].items():
        assert v["present_in_repro"] is True, col
        assert v["cells_over_tolerance"] == 0, col
        assert v["cells"] == 8406, col


@ran_b
def test_the_frozen_artifacts_were_not_written_to():
    """The reproduction compares AGAINST the frozen files; it must never write INTO them."""
    b = _json(B_JSON)
    assert b["reproduction_root"].replace("\\", "/").endswith("results/stage24/repro")
    for label, v in b["frozen_artifacts_untouched"].items():
        rel = v["path"].replace("\\", "/")
        blob = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT,
                              capture_output=True).stdout
        if blob:
            assert hashlib.sha256(blob).hexdigest() == v["sha256"], (
                f"{label} differs from git HEAD -- a frozen artifact was modified")


@ran_b
def test_a_tolerance_declared_verdict_would_require_a_named_cause():
    """§7.1 R4. 'Floating point' alone is not a cause."""
    b = _json(B_JSON)
    if b["reproduction_verdict"] == "BYTE_IDENTICAL":
        assert b["R4_cause_named"] is None
    else:
        assert b["R4_cause_named"], "a non-identical reproduction must name its cause"


@ran_b
@ran_a
def test_24b_could_not_have_run_before_24a():
    """The module refuses to reproduce before the freeze is verified."""
    src = (ROOT / "experiments" / "run_stage24_gen1_tool.py").read_text(encoding="utf-8")
    assert 'raise RuntimeError("24A must run before 24B")' in src
    assert 'raise RuntimeError("24A did not pass; 24B may not run")' in src


# ============================================================================================== #
# 24D–24G — the handoff Stage 25 consumes
# ============================================================================================== #
D_JSON = OUT / "stage24d_handoff_table.json"
E_JSON = OUT / "stage24e_determinism_and_leakage.json"
F_JSON = OUT / "stage24f_tool_freeze.json"
G_JSON = RESULTS / "stage24_handoff_to_stage25.json"
D_CSV = OUT / "stage24_oof_for_stage25.csv"

ran_d = pytest.mark.skipif(not D_JSON.exists(), reason="24D has not been run")
ran_e = pytest.mark.skipif(not E_JSON.exists(), reason="24E has not been run")
ran_f = pytest.mark.skipif(not F_JSON.exists(), reason="24F has not been run")
ran_g = pytest.mark.skipif(not G_JSON.exists(), reason="24G has not been run")


@ran_d
def test_the_handoff_table_is_complete_and_the_population_is_892():
    d = _json(D_JSON)
    assert d["all_checks_pass"] is True
    assert d["rows"] == 8406 and d["clones"] == 1401
    assert d["eligible_clones"] == 892 == d["expected_eligible"]
    assert d["excluded"] == {"all_zero": 472, "all_positive": 37}
    for k in ("pred_W1", "pred_W4", "pred_W5", "ranking_eligible", "outer_fold"):
        assert k in d["columns"], k


@ran_d
def test_stage_24_computed_no_ranking_statistic():
    """Stage 24 emits Stage 25's inputs. It is forbidden from looking at the answer."""
    d = _json(D_JSON)
    assert d["ranking_statistic_computed"] is False
    # Scan the DATA, not the prose. `note` and `eligibility_rule` deliberately name the quantities
    # that were NOT computed, so including them would make this contract fail on its own wording.
    data = {k: v for k, v in d.items() if k not in ("note", "eligibility_rule")}
    blob = json.dumps(data).lower()
    for banned in ("delta_rank", "auroc", "auc_i", "top1", "delta_top1"):
        assert banned not in blob, f"a ranking quantity leaked into the handoff table: {banned}"


@ran_e
def test_scoring_is_deterministic_and_reproduces_the_frozen_column():
    e = _json(E_JSON)
    assert e["determinism"]["within_session"] is True
    assert e["determinism"]["across_loads"] is True
    assert e["determinism"]["max_abs_diff_vs_frozen"] < 1e-12
    assert e["clones_sampled"] >= 50


@ran_e
def test_every_fold_component_is_isolated_from_the_clones_it_scores():
    """Verified from the artifact's own recorded training set, not inferred from the OOF."""
    e = _json(E_JSON)
    assert e["leakage"]["every_fold_isolated"] is True
    for f in range(5):
        iso = e["leakage"]["fold_isolation"][f"fold{f}"]
        assert iso["overlap"] == 0, f"fold{f} trained on clones it is used to score"
        assert iso["held_out_clones"] > 0 and iso["train_clones"] > 0
    assert e["leakage"]["deployment_trained_on_all_clones"] is True


@ran_e
def test_the_artifact_carries_no_outcome_and_no_lineage_feature():
    e = _json(E_JSON)
    assert e["leakage"]["no_outcome_length_array_in_artifact"] is True
    assert e["leakage"]["outcome_length_arrays_found"] == []
    assert e["leakage"]["n_expression_features"] == 36601
    assert e["leakage"]["feature_space_is_gene_expression_only"] is True


@ran_f
def test_every_plan_6_5_deliverable_exists():
    f = _json(F_JSON)
    assert f["all_deliverables_present"] is True
    for k in ("python_prediction_api", "command_line_interface", "frozen_model_artifact",
              "machine_readable_schemas", "model_card", "example_dataset", "unit_tests",
              "end_to_end_reproduction"):
        assert k in f["deliverables"], k
    for name in ("MODEL_CARD.md", "io_schema.json", "example_clones.csv"):
        assert name in f["hashes"]


@ran_f
def test_the_model_card_carries_the_limitations_and_the_ranking_caveat():
    card = (OUT / "tool" / "MODEL_CARD.md").read_text(encoding="utf-8")
    flat = " ".join(card.replace("**", "").split())
    assert "not validated on held-out data" in flat
    assert "order is not a validated condition ranking" in flat
    assert "Not a clinical tool" in flat
    assert "3.45x the state contribution" in flat


@ran_g
def test_stage_24_hands_off_ready_and_names_what_stage_25_may_not_do():
    g = _json(G_JSON)
    assert g["stage_24_verdict"] == "STAGE_24_GEN1_TOOL_READY"
    assert g["substage_results"]["24B"] == "BYTE_IDENTICAL"
    assert g["substage_results"]["24C"] == "SERIALIZED_AND_EQUIVALENT"
    assert g["ranking_metric_inspected_by_stage_24"] is False
    assert g["ranking_statistic_computed_by_stage_24"] is False
    assert g["ranking_population"]["eligible_clones"] == 892
    may_not = " ".join(g["stage_25_may_not"])
    for item in ("change the metric", "reduce the permutation count", "rescue a failed C1",
                 "add a dataset", "revise the plan after seeing a result"):
        assert item in may_not, item


@ran_g
def test_the_handoff_table_hash_matches_the_file_on_disk():
    g = _json(G_JSON)
    assert hashlib.sha256(D_CSV.read_bytes()).hexdigest() == g["frozen_oof_table"]["sha256"]
