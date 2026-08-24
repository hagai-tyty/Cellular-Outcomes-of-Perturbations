"""Stage 23.2A — contracts for the resolution protocol and source-design freeze.

23.2A is the substage that decides whether the whole failure decomposition is trustworthy, and it
is unusually easy to get wrong in ways nothing downstream would notice:

1. **Rewriting history.** Stage 23 is closed and read-only. A diagnostic stage that quietly touched
   a Stage-23 artifact, or that let the closed Role-A verdict drift, would invalidate everything it
   then claimed to explain.
2. **A replay that is not the historical pipeline.** The paired 23.2B/C design only works because
   `D00` really is the historical null. If the replay were a fresh sample rather than the frozen
   mappings, "paired" would be a lie and the Monte-Carlo noise it removes would still be there.
3. **Over-claiming the source design.** The two control GSMs are one biological replicate. A stage
   that let `REPLICATE_STRUCTURE_BIOLOGICAL` back in, or that read a lane split into ambiguous
   metadata, would manufacture experimental units that do not exist.
4. **Peeking at reserved evidence.** The reserved ledger names biological replicates 2 and 3. If
   anything beyond declared metadata reached it, the confirmation evidence would already be burned.

These contracts test all four, plus the canonical-JSON provenance design that exists specifically
so 23.2B's code cannot invalidate 23.2A's frozen protocol.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "run_stage23_2_role_a_resolution.py"
spec = importlib.util.spec_from_file_location("s232", SRC)
S232 = importlib.util.module_from_spec(spec)
sys.modules["s232"] = S232
spec.loader.exec_module(S232)

RES = ROOT / "results"
OUT = RES / "stage23_2"
A_RESULTS = OUT / "stage23_2a_results.json"
PROTOCOL = OUT / "stage23_2_protocol.json"
LEDGER = OUT / "stage23_2_reserved_confirmation_candidates.json"
D00 = OUT / "stage23_2_historical_null_d00.json"
SOURCE_DESIGN = OUT / "stage23_2_source_design.json"

ran = pytest.mark.skipif(not A_RESULTS.exists(), reason="23.2A has not been run")


@pytest.fixture(scope="module")
def a():
    return json.loads(A_RESULTS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def design():
    return json.loads(SOURCE_DESIGN.read_text(encoding="utf-8"))


# ---- 1. Stage 23 stays read-only ------------------------------------------------------------- #
@ran
def test_the_closed_stage23_verdicts_are_untouched(a):
    syn = json.loads((RES / "stage23_final_synthesis.json").read_text(encoding="utf-8"))
    assert syn["final_verdicts"]["role_a"] == "ROLE_A_SIGNAL_FAIL"
    assert syn["final_verdicts"]["role_b_additive"] == "ROLE_B_ADDITIVE_PASS"
    assert syn["final_verdicts"]["role_b_interaction"] == "INTERACTION_PASS_MULTI_TREATMENT"
    assert syn["STRUCTURAL_CONTROLS_PASS"] is True
    assert syn["roadmap_gate"]["gate"] == "STAGE_24_BLOCKED_ROLE_A"


@ran
def test_the_historical_artifact_hashes_recorded_by_23_2a_still_match(a):
    """If a Stage-23 artifact changed after 23.2A froze, every diagnostic below it is suspect."""
    for name, digest in a["preflight"]["artifact_hashes"].items():
        assert S232.sha256_file(RES / name) == digest, name


@ran
def test_23_2a_writes_nothing_outside_its_own_directory(a):
    """Every artifact 23.2A claims to have written lives under results/stage23_2/."""
    for name in a["artifacts"]:
        assert (OUT / name).exists(), name
    assert not list(RES.glob("stage23_2_*.json")), "an artifact escaped into results/"


@ran
def test_the_closure_record_still_declares_stage23_closed(a):
    checks = {c["check"]: c["ok"] for c in a["preflight"]["checks"]}
    assert checks["closure record exists and declares Stage 23 formally closed"]
    assert checks["closure record states ROLE_A_SIGNAL_FAIL is permanent"]
    assert checks["legacy STAGE 23R is only an alias, not a competing gate"]


@ran
def test_every_frozen_anchor_was_verified(a):
    assert a["preflight"]["ok"] is True
    assert a["preflight"]["failed"] == []
    assert a["preflight"]["n_checks"] >= 30


# ---- 2. the replay really is the historical pipeline ------------------------------------------ #
@ran
def test_the_committed_d00_array_has_200_values_and_a_matching_digest():
    d = json.loads(D00.read_text(encoding="utf-8"))
    assert d["n_permutations"] == len(d["values"]) == 200
    assert d["base_seed"] == 23323
    assert d["values_sha256"] == hashlib.sha256(np.array(d["values"]).tobytes()).hexdigest()


@ran
def test_the_replayed_array_reproduces_the_committed_stage23_summary_exactly(a):
    """The load-bearing claim: this is the historical null, not a fresh sample of it."""
    rep = a["permutation_recovery"]["replay"]
    assert rep["values_recomputed"] == 200
    assert rep["all_summary_statistics_reproduced"] is True
    assert all(rep["matches_committed_summary"].values()), rep["matches_committed_summary"]

    hist = json.loads((RES / "stage23_permutation_results.json").read_text(encoding="utf-8"))
    t = hist["permutation_tests"]["role_a_delta_AP_state"]
    arr = np.array(json.loads(D00.read_text(encoding="utf-8"))["values"])
    assert round(float(arr.mean()), 12) == round(t["null_mean"], 12)
    assert round(float(np.percentile(arr, 95)), 12) == round(t["null_p95"], 12)
    assert int((arr >= t["observed"]).sum()) == t["n_null_ge_observed"] == 16
    assert round(float((1 + (arr >= t["observed"]).sum()) / 201), 12) == round(t["p_perm"], 12)


@ran
def test_the_mapping_set_is_recorded_by_digest_and_kept_out_of_git(a):
    m = a["permutation_recovery"]
    assert m["n_permutations"] == 200
    assert m["mapping_is_cache_only"] is True
    assert len(m["mapping_set_sha256"]) == 64
    assert m["mapping_rows"] == 200 * 5 * 3147
    assert m["mapping_cache_dir"].startswith("_cc_cache/")
    assert not list(OUT.glob("*mapping*")), "the mapping table must not be committed"


def test_a_regenerated_mapping_obeys_the_frozen_permutation_structure():
    """Recomputed here rather than trusted: no crossing, no stratum change, still a bijection."""
    k = pd.read_csv(RES / "stage22_rewind_clones.csv")
    clones = json.loads((ROOT / "_cc_cache" / "stage23" / "GSE227151_clones.json")
                        .read_text(encoding="utf-8")) if (
        ROOT / "_cc_cache" / "stage23" / "GSE227151_clones.json").exists() else None
    if clones is None:
        pytest.skip("23A clone cache absent")
    k = k.set_index("clone_id").loc[clones]
    strata = S232.S23.rewind_strata(k)
    fold = k["outer_fold"].to_numpy()
    for b in (0, 7, 199):
        rng = np.random.default_rng(S232.S23.SEED_PERMUTATION + b)
        for f in range(S232.S23.N_OUTER):
            pmap = S232.S23.permute_within(strata, fold != f, rng)
            side = fold != f
            assert (side[pmap] == side).all()
            assert (strata[pmap] == strata).all()
            assert sorted(pmap.tolist()) == list(range(len(pmap)))


@ran
def test_the_realized_strata_are_the_five_non_empty_cells(a):
    cells = a["permutation_recovery"]["realized_strata"]["cells"]
    assert cells == {"1|1": 2584, "2|1": 220, "2|2": 196, "3+|1": 37, "3+|2": 110}
    assert "1|2" not in cells, "one pretreatment cell cannot span two lanes"
    assert sum(cells.values()) == 3147


# ---- 3. the source design is not over-claimed ------------------------------------------------- #
@ran
def test_the_within_r1_status_is_one_of_the_three_allowed_values(a, design):
    allowed = {"WITHIN_R1_TECHNICAL_LANES", "WITHIN_R1_SEPARATE_LIBRARIES",
               "WITHIN_R1_STRUCTURE_UNRESOLVED"}
    assert a["within_r1_status"] in allowed
    assert design["source_design"]["within_r1_status"] == a["within_r1_status"]


@ran
def test_the_retired_replicate_vocabulary_never_appears_in_any_artifact():
    """`REPLICATE_STRUCTURE_BIOLOGICAL` was removed as a possible finding in V2."""
    for p in sorted(OUT.glob("*")):
        text = p.read_text(encoding="utf-8", errors="replace")
        for banned in ("REPLICATE_STRUCTURE_BIOLOGICAL", "REPLICATE_STRUCTURE_TECHNICAL_LANES",
                       "REPLICATE_STRUCTURE_UNRESOLVED"):
            assert banned not in text, f"{p.name} contains {banned}"


@ran
def test_the_benchmark_is_recorded_as_one_biological_replicate(design):
    s = design["source_design"]["biological_replicate_count_is_settled"]
    assert s["value"] == 1
    assert s["label"] == "R1"
    assert s["benchmark_evidence"]["biological_replicate_column"] == ["R1"]
    cells = pd.read_csv(RES / "stage22_rewind_cells.csv")
    assert sorted(cells["biological_replicate"].unique()) == ["R1"]


@ran
def test_the_lane_sensitivity_flag_follows_the_status(a):
    """V2 §7.5 may run only under WITHIN_R1_TECHNICAL_LANES."""
    assert a["lane_composition_sensitivity_permitted"] == (
        a["within_r1_status"] == "WITHIN_R1_TECHNICAL_LANES")


@ran
def test_the_sample_numbering_conflict_is_recorded_and_the_benchmark_is_not_re_derived(design):
    """F2: GEO titles and file naming disagree; the benchmark keyed on SampleNum and stays frozen."""
    c = design["source_design"]["sample_numbering_conflict"]
    assert c["conflict_present"] is True
    assert c["geo_title_map"] != c["file_naming_map"]
    assert c["benchmark_agrees_with_resolved_map"] is True
    cells = pd.read_csv(RES / "stage22_rewind_cells.csv")
    live = {str(k): sorted(set(v))[0]
            for k, v in cells.groupby("SampleNum")["gsm"].agg(list).items()}
    assert c["benchmark_SampleNum_to_gsm"] == live, "the frozen mapping was altered"


@ran
def test_non_discriminating_evidence_is_not_used_to_decide_the_status(design):
    """Per-sample GEM loading and 10x indexing are true of both designs, so they must not decide."""
    e = design["source_design"]["declared_evidence"]
    assert e["per_sample_gem_and_indexing_are_non_discriminating"] is True
    assert e["extract_protocol_lines_read"] >= 4, "the protocol spans several series-matrix lines"
    if design["source_design"]["within_r1_status"] == "WITHIN_R1_STRUCTURE_UNRESOLVED":
        assert e["metadata_declares_lane_split"] is False
        assert e["metadata_declares_separate_source_material"] is False
        assert e["characteristics_differ_between_gsms"] is False


# ---- 4. reserved evidence stays untouched ------------------------------------------------------ #
@ran
def test_the_reserved_ledger_holds_declared_metadata_only():
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    allowed = {"accession", "title", "declared_biological_replicate", "library_strategy",
               "library_source", "declared_gating", "role", "locally_downloaded",
               "matching_future_outcome_declared"}
    for e in led["entries"]:
        assert set(e) <= allowed, set(e) - allowed
        assert e["matching_future_outcome_declared"] is None, "outcome status must stay unverified"
    assert led["n_samples_declared"] == 13
    assert "UNVERIFIED" in led["matching_outcome_status"]


@ran
def test_the_ledger_names_the_reserved_replicates_without_evaluating_them():
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    by_acc = {e["accession"]: e for e in led["entries"]}
    assert by_acc["GSM7092515"]["role"] == "USED_BY_STAGE_23"
    assert by_acc["GSM7092516"]["role"] == "USED_BY_STAGE_23"
    assert by_acc["GSM7092517"]["declared_biological_replicate"] == "2"
    assert by_acc["GSM7092519"]["declared_biological_replicate"] == "3"
    assert by_acc["GSM7092520"]["role"] == "RESERVED_DIFFERENT_DESIGN"
    for e in led["entries"]:
        if e["role"].startswith("RESERVED"):
            assert e["locally_downloaded"] is False, f"{e['accession']} was downloaded"


def test_the_module_never_opens_a_reserved_matrix():
    """Source-level: only the two Stage-23 GSMs may be read from disk.

    A reserved accession may be NAMED in prose -- the 23.2F confirmation protocol has to list them
    as forbidden-until-qualified, and 23.2A's ledger describes them. What must never happen is a
    reserved accession appearing in code that reaches the filesystem. So the check is scoped: the
    only function permitted to mention one is the confirmation-protocol writer, which emits
    documentation, and the ledger must still be derived from family.xml rather than hardcoded.
    """
    import ast

    src = SRC.read_text(encoding="utf-8")
    assert 'STAGE23_GSMS = ("GSM7092515", "GSM7092516")' in src
    reserved = ("GSM7092517", "GSM7092518", "GSM7092519", "GSM7092520", "GSM7092521")
    documentation_only = {"_write_confirmation_protocol"}

    tree = ast.parse(src)
    for fn in (n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)):
        if fn.name in documentation_only:
            continue
        body = ast.get_source_segment(src, fn) or ""
        for acc in reserved:
            assert acc not in body, f"{fn.name} names {acc}; data access must not reach a reserved GSM"

    # the ledger is still built from the declared metadata file, not from a hardcoded list
    ledger_fn = src.split("def reserved_confirmation_ledger(")[1].split(chr(10) + "def ")[0]
    assert "FAMILY_XML" in ledger_fn
    for acc in reserved:
        assert acc not in ledger_fn


# ---- 5. the gDNA rule and Bdepth ---------------------------------------------------------------#
@ran
def test_the_gdna_rule_reproduces_the_frozen_label_exactly(a, design):
    g = design["gdna_rule"]
    assert g["reproduces_frozen_positives_exactly"] is True
    assert g["reconstructed_positive_clones"] == 35
    assert g["selected_barcodes"] == 101, "the rank-100 tie yields 101 barcodes"
    assert g["tie_size_at_cutoff"] == 2
    assert g["rank_100_cutoff_counts"] == 2365
    assert g["support_column"] == "counts", "M4: the gDNA table has no nUMI column"
    assert g["sample_num_grouping_is_a_no_op"] is True
    assert g["sample_num_values_in_gdna"] == [3], "gDNA is one pooled library"


@ran
def test_bdepth_is_complete_and_outcome_free(design):
    b = design["bdepth"]
    assert b["clones"] == 3147
    assert b["all_positive_total_umi"] is True
    assert b["detected_matches_normalised_nonzero_pattern"] is True
    tbl = pd.read_csv(OUT / "stage23_2_bdepth.csv")
    assert len(tbl) == 3147
    for banned in ("y_primed", "outcome", "gdna", "counts", "primed"):
        assert not [c for c in tbl.columns if banned in c.lower()], banned


# ---- 6. the canonical-JSON provenance design ---------------------------------------------------#
def test_the_protocol_digest_is_canonical_and_order_independent():
    a = {"b": 1, "a": [1, 2], "c": {"z": 0, "y": 1}}
    b = {"c": {"y": 1, "z": 0}, "a": [1, 2], "b": 1}
    assert S232.canonical_json_sha256(a) == S232.canonical_json_sha256(b)
    assert S232.canonical_json_sha256({"a": 1}) != S232.canonical_json_sha256({"a": 2})


@ran
def test_the_hashed_protocol_payload_excludes_source_and_runtime_provenance(a):
    """V2 §4.1/§4.2 -- the fix for the Stage-23 builder-hash problem, which bit three times."""
    doc = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    payload = doc["protocol"]
    assert doc["canonical_sha256"] == S232.canonical_json_sha256(payload)
    assert doc["canonical_sha256"] == a["stage23_2_protocol_sha256"]
    blob = json.dumps(payload)
    for banned in ("git_commit", "builder", "sha256_of_source", "timestamp", "platform",
                   "D:\\\\", "/mnt/", "runtime_minutes"):
        assert banned not in blob, f"the hashed payload contains {banned}"
    assert "source_provenance" in doc, "provenance is recorded beside the payload, not inside it"
    assert "git_commit" in doc["source_provenance"]


@ran
def test_adding_later_substage_code_cannot_move_the_frozen_protocol_digest(a):
    """The digest depends on the scientific surface only, so 23.2B code cannot invalidate it."""
    doc = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    payload = dict(doc["protocol"])
    first = S232.canonical_json_sha256(payload)
    # a source change is modelled by changing only the provenance block
    doc["source_provenance"]["source_files"]["stage23_2_builder"] = "0" * 64
    assert S232.canonical_json_sha256(payload) == first


@ran
def test_the_frozen_design_matches_the_plan(a):
    doc = json.loads(PROTOCOL.read_text(encoding="utf-8"))["protocol"]
    d = doc["decomposition_design"]
    assert d["no_k_selection_reference"]["fixed_K_arms"] == [10, 20, 50]
    assert "equal" in d["no_k_selection_reference"]["weights"]
    assert d["search_width_ladder"]["conditional"] is False
    assert set(d["cells"]) == {"00", "01", "10", "11"}
    assert doc["label_reliability"]["not_supported_requires_independent_outcome_assay_replication"]
    assert "REMOVED" in doc["label_reliability"]["cross_gsm_gdna_concordance"]
    assert doc["power"]["n_biological_replicates"] == 1
    assert set(doc["power"]["statuses"]) == {"WITHIN_R1_EVENT_COUNT_LIMITATION",
                                             "BIOLOGICAL_REPLICATION_LIMITATION"}


@ran
def test_all_23_2a_gates_pass_and_the_verdict_is_derived(a):
    assert all(a["gates"].values()), [k for k, v in a["gates"].items() if not v]
    assert a["verdict"] == S232.PROTOCOL_FROZEN
    assert a["verdict"] == (S232.PROTOCOL_FROZEN if all(a["gates"].values())
                            else S232.INPUT_BLOCKED)


# ============================================================================================== #
# 23.2B — model-selection null decomposition.
#
# The substage is only meaningful if three things hold, and each can fail silently:
#
#   * the paired basis really is the historical one. `S_j = D00_j - D10_j` is a paired statistic;
#     if `D00` were recomputed rather than read, or the mappings redrawn, the pairing would be
#     fictional and the CI too wide or too narrow for reasons nothing would report.
#   * the no-K-selection reference is the equal-weight mean of the three fixed-K arms, with weights
#     fixed before execution. Promoting whichever arm looks best would be exactly the selection
#     effect the substage claims to measure.
#   * the status is derived from the CI, not chosen. UNRESOLVED must be reachable and reported.
# ============================================================================================== #
B_RESULTS = OUT / "stage23_2_model_selection_decomposition.json"
ran_b = pytest.mark.skipif(not B_RESULTS.exists(), reason="23.2B has not been run")


@pytest.fixture(scope="module")
def b():
    return json.loads(B_RESULTS.read_text(encoding="utf-8"))


@ran_b
def test_cell_00_is_read_from_the_committed_array_and_recomputes_bitwise(b):
    """The paired basis. A recomputation that differed would mean the engine is not the historical
    pipeline, and every paired CI below would be measuring the wrong thing."""
    c = b["cell_00_recomputation"]
    assert c["reproduces_committed_D00_exactly"] is True
    assert c["max_abs_difference"] == 0.0
    committed = np.array(json.loads(D00.read_text(encoding="utf-8"))["values"])
    assert round(b["D00"]["mean"], 12) == round(float(committed.mean()), 12)
    assert b["historical_null_artifact_sha256"] == S232.sha256_file(D00)


@ran_b
def test_the_expression_free_reference_reproduces_the_historical_r1_exactly(b):
    rb = json.loads((RES / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    assert b["reference"]["reproduces_historical_R1_exactly"] is True
    assert b["reference"]["r1_reference_AP"] == rb["pooled_oof_metrics"]["R1"]["AP"]


@ran_b
def test_the_no_k_reference_is_the_equal_weight_mean_of_three_arms(b):
    arms = [b["per_arm"][f"K{k}"]["mean"] for k in (10, 20, 50)]
    assert len(arms) == 3
    assert round(b["D10_no_k_selection"]["mean"], 10) == round(float(np.mean(arms)), 10), (
        "the reference must be the equal-weight mean, not a chosen arm")
    for k in (10, 20, 50):
        assert f"K{k}" in b["per_arm"], f"arm K{k} must be reported as a diagnostic"
    assert b["arm_dispersion"] == pytest.approx(max(arms) - min(arms))


@ran_b
def test_no_arm_was_promoted_to_the_primary_reference(b):
    """If the reference equalled the best-looking arm, the equal weighting would be cosmetic."""
    arms = {k: v["mean"] for k, v in b["per_arm"].items()}
    ref = b["D10_no_k_selection"]["mean"]
    best = max(arms.values())
    assert ref != best or len(set(arms.values())) == 1, "the reference equals the strongest arm"


@ran_b
def test_the_selection_shift_is_the_paired_difference(b):
    s = b["selection_shift"]
    assert round(s["mean"], 10) == round(b["D00"]["mean"] - b["D10_no_k_selection"]["mean"], 10)
    assert s["resamples"] == 10_000
    assert s["seed"] == 23421
    assert s["ci95_low"] < s["mean"] < s["ci95_high"]


@ran_b
def test_the_status_is_derived_from_the_confidence_interval(b):
    s = b["selection_shift"]
    expected = ("SUPPORTED" if s["ci95_low"] > 0
                else "NOT_SUPPORTED" if s["ci95_high"] <= 0 else "UNRESOLVED")
    assert b["MODEL_SELECTION_NULL_INFLATION"] == expected


@ran_b
def test_the_search_width_ladder_ran_unconditionally(b):
    """V2 §6.4 made the ladder unconditional; it must be present whatever the primary status."""
    ladder = b["search_width_ladder"]
    assert set(ladder) == {"4_candidate_fixed_K_mean", "8_candidate", "12_candidate"}
    assert ladder["12_candidate"] == pytest.approx(b["D00"]["mean"])
    assert b["ladder_monotone_increase"] == (
        ladder["4_candidate_fixed_K_mean"] <= ladder["8_candidate"] <= ladder["12_candidate"])


@ran_b
def test_the_fraction_explained_is_only_reported_for_a_positive_null_mean(b):
    frac = b["fraction_null_mean_explained_by_search"]
    if b["D00"]["mean"] > 0:
        assert frac == pytest.approx(b["selection_shift"]["mean"] / b["D00"]["mean"])
    else:
        assert frac is None


@ran_b
def test_the_observed_sensitivity_keeps_the_historical_value_intact(b):
    """V2 §6.5 is effect attribution, not a rescue: the historical observed dAP must be unchanged."""
    o = b["observed_sensitivity"]["delta_AP"]
    assert o["hist12"] == pytest.approx(0.01050162935116511, abs=1e-9)
    assert o["no_k_selection"] == pytest.approx(np.mean([o["K10"], o["K20"], o["K50"]]))
    assert "p_value" not in json.dumps(b["observed_sensitivity"]), \
        "no rescue p-value may be computed here"


@ran_b
def test_23_2b_pins_the_frozen_protocol_and_mapping_set(b):
    proto = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    a = json.loads(A_RESULTS.read_text(encoding="utf-8"))
    assert b["stage23_2_protocol_sha256"] == proto["canonical_sha256"]
    assert b["mapping_set_sha256"] == a["permutation_recovery"]["mapping_set_sha256"]
    assert b["n_permutations"] == 200


# ============================================================================================== #
# 23.2D — outcome-label reliability.
#
# The contract that matters most here is the asymmetry V2 design change 3 introduced. Every one of
# V1's four stability criteria is met by this data, so V1's rule would have returned NOT_SUPPORTED
# and the ledger would have recorded "the label is sound" on the strength of two diagnostics that
# observe only sequencing-count noise and cutoff position. V2 makes NOT_SUPPORTED additionally
# require independent outcome-assay replication, which does not exist in the Rewind materials.
#
# If a later edit ever relaxes that gate, this file should fail loudly.
# ============================================================================================== #
D_RESULTS = OUT / "stage23_2_label_reliability.json"
ran_d = pytest.mark.skipif(not D_RESULTS.exists(), reason="23.2D has not been run")


@pytest.fixture(scope="module")
def d():
    return json.loads(D_RESULTS.read_text(encoding="utf-8"))


@ran_d
def test_reliability_analysis_only_runs_after_the_label_reproduces_exactly(d):
    """V2 §8.1 makes exact reproduction a precondition, not a finding."""
    assert d["source_rule_reproduces_35_positives"] is True
    assert d["cutoff_sensitivity"]["top100"]["positive_clones"] == 35
    assert d["cutoff_sensitivity"]["top100"]["selected_barcodes"] == 101
    assert d["tie_size_at_rank_100"] == 2


@ran_d
def test_not_supported_is_unreachable_without_independent_replication(d):
    """Design change 3, tested literally: stability alone must NOT clear the label."""
    m = d["multinomial_stability"]
    ladder = d["cutoff_sensitivity"]
    all_v1_criteria_met = (m["mean_frozen_positive_retention"] >= 0.90
                           and m["n_positives_below_0_80"] <= 3
                           and ladder["top90"]["jaccard_vs_top100"] >= 0.90
                           and ladder["top110"]["jaccard_vs_top100"] >= 0.90)
    if all_v1_criteria_met and not d["independent_outcome_assay_replication_available"]:
        assert d["OUTCOME_LABEL_LIMITATION"] == "UNRESOLVED", (
            "V1's rule would have said NOT_SUPPORTED here; V2 must not")
    assert d["not_supported_reachable"] == d["independent_outcome_assay_replication_available"]
    if d["OUTCOME_LABEL_LIMITATION"] == "NOT_SUPPORTED":
        assert d["independent_outcome_assay_replication_available"] is True


@ran_d
def test_instability_can_still_reach_supported(d):
    """The asymmetry must be one-directional: fragility is still demonstrable."""
    m = d["multinomial_stability"]
    ladder = d["cutoff_sensitivity"]
    unstable_a = (m["mean_frozen_positive_retention"] < 0.80
                  or m["n_positives_below_0_50"] >= 7)
    unstable_b = min(ladder["top90"]["jaccard_vs_top100"],
                     ladder["top110"]["jaccard_vs_top100"]) < 0.80
    if unstable_a and unstable_b:
        assert d["OUTCOME_LABEL_LIMITATION"] == "SUPPORTED"


@ran_d
def test_the_multinomial_scope_is_declared_and_single_unit(d):
    m = d["multinomial_stability"]
    assert m["selection_units"] == 1, "the gDNA table is one pooled library"
    assert m["resamples"] == 5000
    assert m["seed"] == 23431
    assert m["total_gdna_counts_N"] == 782826
    assert m["distinct_barcodes"] == 1936
    assert "sequencing-count sampling noise only" in m["scope"]


@ran_d
def test_the_cutoff_geometry_covers_ranks_80_to_120(d):
    g = d["cutoff_geometry_ranks_80_to_120"]
    assert [x["rank"] for x in g] == list(range(80, 121))
    at_100 = next(x for x in g if x["rank"] == 100)
    assert at_100["counts"] == d["rank_100_cutoff_counts"] == 2365
    assert at_100["tie_size_at_this_value"] == 2
    assert at_100["ratio_to_rank_100_cutoff"] == 1.0


@ran_d
def test_every_frozen_positive_has_a_recorded_selection_probability(d):
    per = d["multinomial_stability"]["per_frozen_positive"]
    assert len(per) == 35
    assert all(0.0 <= x["P_selected"] <= 1.0 for x in per)
    assert per == sorted(per, key=lambda x: x["P_selected"]), "reported least-stable first"


@ran_d
def test_the_cutoff_ladder_is_complete_and_consistent(d):
    ladder = d["cutoff_sensitivity"]
    assert set(ladder) == {f"top{n}" for n in (80, 90, 100, 110, 120)}
    assert ladder["top100"]["jaccard_vs_top100"] == 1.0
    assert ladder["top100"]["frozen_positives_lost"] == 0
    assert ladder["top100"]["frozen_negatives_gained"] == 0
    for n in (80, 90):
        assert ladder[f"top{n}"]["frozen_negatives_gained"] == 0, "a smaller N cannot add positives"
    for n in (110, 120):
        assert ladder[f"top{n}"]["frozen_positives_lost"] == 0, "a larger N cannot drop positives"


@ran_d
def test_cross_gsm_gdna_concordance_stays_removed(d):
    assert "REMOVED IN V2" in d["cross_gsm_gdna_concordance"]
    blob = json.dumps(d)
    assert "GSM7092515" not in blob and "GSM7092516" not in blob, (
        "per-GSM gDNA quantities are not identifiable and must not appear")


@ran_d
def test_no_predictor_was_fitted_and_no_label_was_shopped(d):
    """V2 §3.7 -- alternate labels may be studied, never selected for predictive performance."""
    assert d["no_predictive_model_fitted"] is True
    assert d["candidate_future_formulations"]["status"] == "EXPLORATORY_PROPOSAL_ONLY"
    if d["OUTCOME_LABEL_LIMITATION"] != "SUPPORTED":
        assert d["candidate_future_formulations"]["candidates"] == []
    # scan the payload for predictive QUANTITIES, after dropping the declarative fields whose
    # own names contain the word ("no_predictive_model_fitted") -- otherwise the check matches
    # the very statement it is verifying
    scan = {k: v for k, v in d.items()
            if k not in ("no_predictive_model_fitted", "candidate_future_formulations",
                         "cross_gsm_gdna_concordance", "status_reason")}
    blob = json.dumps(scan).lower()
    for banned in ("average_precision", "roc_auc", "delta_ap", "auc", "log_loss",
                   "predict_proba"):
        assert banned not in blob, f"23.2D reported a predictive quantity: {banned}"


def test_the_top_n_rule_keeps_ties_at_the_cutoff():
    """The behaviour that turns top-100 into 101 barcodes, checked on a synthetic case."""
    counts = np.array([10, 9, 8, 8, 7, 6])
    mask = S232._select_top_n(counts, 3)
    assert mask.tolist() == [True, True, True, True, False, False]
    assert S232._select_top_n(counts, 6).sum() == 6
    assert S232._select_top_n(counts, 99).sum() == 6


# ============================================================================================== #
# 23.2C — residual depth / nuisance decomposition.
#
# This is the substage that found the dominant mechanism, so the contracts guard the three ways a
# strong-looking result here could be wrong:
#
#   * the 2x2 must be balanced. `D11` has to use the SAME three fixed-K arms as `D10`, or the
#     "factor interaction" term is comparing different constructions and means nothing.
#   * `SUPPORTED` requires outcome-level evidence AND a mechanistic correlation. A depth shift on
#     its own could be an artifact of adding two columns to a nuisance model; the outcome-free
#     retention correlations are what make it mechanistic.
#   * the corrected same-data diagnostic must stay exploratory. It came within a hair of clearing
#     its gate, which is exactly when a stage is tempted to promote it.
# ============================================================================================== #
C_RESULTS = OUT / "stage23_2_depth_decomposition.json"
ran_c = pytest.mark.skipif(not C_RESULTS.exists(), reason="23.2C has not been run")


@pytest.fixture(scope="module")
def c():
    return json.loads(C_RESULTS.read_text(encoding="utf-8"))


@ran_c
def test_the_2x2_is_balanced_across_both_nuisance_blocks(c):
    """D10 and D11 must be built the same way, or the interaction term is meaningless."""
    proto = json.loads(PROTOCOL.read_text(encoding="utf-8"))["protocol"]
    assert proto["decomposition_design"]["no_k_selection_reference"]["fixed_K_arms"] == [10, 20, 50]
    assert set(c["per_arm_Bdepth"]) == {"K10", "K20", "K50"}
    arms = [c["per_arm_Bdepth"][f"K{k}"]["mean"] for k in (10, 20, 50)]
    assert round(c["null_means"]["mu11"], 10) == round(float(np.mean(arms)), 10), (
        "mu11 must be the equal-weight mean of the same three arms used for mu10")
    assert c["arm_dispersion_Bdepth"] == pytest.approx(max(arms) - min(arms))


@ran_c
def test_the_nuisance_block_is_the_frozen_bdepth(c):
    assert c["nuisance_Bdepth"] == [
        "log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
        "log1p(n_detected_GE_features_in_raw_pseudobulk)"]
    joined = " ".join(c["nuisance_Bdepth"]).lower()
    for banned in ("gdna", "y_primed", "primed", "outcome"):
        assert banned not in joined, f"Bdepth contains an outcome-derived term: {banned}"


@ran_c
def test_the_contrast_algebra_is_consistent_with_the_reported_means(c):
    m = c["null_means"]
    assert c["contrasts"]["selection_main_contrast"]["mean"] == pytest.approx(
        m["mu00"] - m["mu10"], abs=1e-9)
    assert c["contrasts"]["depth_main_contrast"]["mean"] == pytest.approx(
        m["mu00"] - m["mu01"], abs=1e-9)
    assert c["contrasts"]["factor_interaction"]["mean"] == pytest.approx(
        m["mu00"] - m["mu10"] - m["mu01"] + m["mu11"], abs=1e-9)
    assert c["depth_shift_full"]["mean"] == pytest.approx(m["mu00"] - m["mu01"], abs=1e-9)


@ran_c
def test_supported_requires_both_outcome_level_and_mechanistic_evidence(c):
    """V2 §7.6 -- a depth shift alone is not enough; a retention correlation must also hold."""
    mech = [v for v in c["technical_retention"].values() if isinstance(v, dict)]
    any_pos = any(v["excludes_zero_positive"] for v in mech)
    all_incl_zero = all(not v["excludes_zero_positive"] for v in mech)
    d = c["depth_shift_full"]
    expected = ("SUPPORTED" if (d["ci95_low"] > 0 and any_pos)
                else "NOT_SUPPORTED" if (d["ci95_high"] <= 0 and all_incl_zero)
                else "UNRESOLVED")
    assert c["RESIDUAL_DEPTH_STRUCTURE"] == expected
    if c["RESIDUAL_DEPTH_STRUCTURE"] == "SUPPORTED":
        assert d["ci95_low"] > 0 and any_pos


@ran_c
def test_the_retention_diagnostic_uses_no_outcome(c):
    assert "no outcome label is used" in c["technical_retention"]["note"]
    for v in c["technical_retention"].values():
        if isinstance(v, dict):
            assert -1.0 <= v["median"] <= 1.0
            assert v["ci95_low"] <= v["mean"] <= v["ci95_high"]


@ran_c
def test_the_corrected_same_data_diagnostic_stays_exploratory(c):
    """It nearly cleared its gate. It must still be unable to emit a confirmatory verdict."""
    d = c["corrected_same_data_diagnostic"]
    expected = ("POSITIVE" if (d["O11"] > 0 and d["O11"] > d["q95_11"] and d["p_diag_11"] <= 0.05)
                else "NEGATIVE")
    assert d["CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC"] == expected
    assert "exploratory" in d["note"]
    assert "ROLE_A_CONFIRMATORY_SUPPORTED" in d["note"]
    blob = json.dumps(c)
    assert "ROLE_A_SIGNAL_PASS" not in blob, "23.2C may not emit a Role-A pass"


@ran_c
def test_the_observed_2x2_preserves_the_historical_cell(c):
    o = c["observed_2x2"]
    assert o["O00_matches_historical_within_tolerance"] is True
    assert o["O00"] == pytest.approx(0.01050162935116511, abs=1e-9)
    assert set(o) >= {"O00", "O10", "O01", "O11"}


@ran_c
def test_lane_composition_sensitivity_respects_the_23_2a_status(c):
    """V2 §7.5 is permitted only under WITHIN_R1_TECHNICAL_LANES."""
    a = json.loads(A_RESULTS.read_text(encoding="utf-8"))
    ls = c["lane_composition_sensitivity"]
    assert ls["within_r1_status"] == a["within_r1_status"]
    assert ls["permitted"] == (a["within_r1_status"] == "WITHIN_R1_TECHNICAL_LANES")
    if not ls["permitted"]:
        assert ls["executed"] is False
        blob = json.dumps(c)
        assert "n_cells_GSM" not in blob, "per-sample counts entered a forbidden diagnostic"


@ran_c
def test_bdepth_is_not_promoted_to_a_production_baseline(c):
    assert "forbids promoting it into a" in c["critical_interpretation"]
    assert "§7.7" in c["critical_interpretation"] or "7.7" in c["critical_interpretation"]


# ============================================================================================== #
# 23.2E — power / identifiability.
#
# The result is a clean monotone ladder, which is exactly when a power study is most dangerous:
# it invites being read as "we just need more clones". The contracts here pin the three boundaries
# V2 draws around that reading.
#
#   * the generating direction must be LABEL-FREE. If `z` had seen `y_primed`, the whole curve
#     would be circular and would overstate detectability.
#   * the two statuses are different KINDS of claim. Event-count detectability is estimated from
#     the curve; biological-replication limitation is a design fact that no simulation can move,
#     and it admits only SUPPORTED / NOT_SUPPORTED.
#   * no bias direction may be claimed for the resampled scales, and the curve may never be
#     described as an independent-sample power curve.
# ============================================================================================== #
E_RESULTS = OUT / "stage23_2_power_identifiability.json"
ran_e = pytest.mark.skipif(not E_RESULTS.exists(), reason="23.2E has not been run")


@pytest.fixture(scope="module")
def e():
    return json.loads(E_RESULTS.read_text(encoding="utf-8"))


@ran_e
def test_all_nine_shards_ran_at_the_frozen_counts(e):
    """V2 §9.4 forbids reducing simulation counts after execution starts."""
    d = e["design"]
    assert d["cohort_scales"] == [1, 2, 4]
    assert sorted(d["target_AUCs"]) == [0.66, 0.70]
    assert d["null_allocations"] == 200
    assert d["alternative_simulations"] == 100
    for s in ("1", "2", "4"):
        v = e["per_scale"][s]
        assert v["n_null"] == 200, f"scale {s} null was shortened"
        assert set(v["power"]) == {"0.66", "0.70"}, f"scale {s} is missing a target AUC"
        for pw in v["power"].values():
            assert pw["n_alt"] == 100


@ran_e
def test_the_frozen_seeds_were_used(e):
    assert e["design"]["seeds"] == {"covariate_resample": 23440, "null": 23441,
                                    "alternative": 23442, "beta_calibration": 23443}


@ran_e
def test_the_cohort_ladder_preserves_the_historical_prevalence(e):
    """Scaling adds events by growing the cohort, not by enriching a fixed one."""
    for s, expect_n, expect_pos in (("1", 3147, 7), ("2", 6294, 14), ("4", 12588, 28)):
        v = e["per_scale"][s]
        assert v["cohort_N"] == expect_n
        assert v["positives_per_fold"] == expect_pos
        prevalence = (v["positives_per_fold"] * 5) / v["cohort_N"]
        assert prevalence == pytest.approx(35 / 3147, rel=0.02), s


@ran_e
def test_beta_calibration_hit_its_targets(e):
    for s in ("1", "2", "4"):
        for target, pw in e["per_scale"][s]["power"].items():
            assert pw["achieved_median_oracle_AUC"] == pytest.approx(float(target), abs=0.02), (
                f"scale {s} target {target} achieved {pw['achieved_median_oracle_AUC']}")


@ran_e
def test_the_null_tightens_and_power_rises_with_cohort_size(e):
    q95 = [e["per_scale"][s]["q95_null"] for s in ("1", "2", "4")]
    assert q95[0] > q95[1] > q95[2], "a larger cohort must tighten the null"
    for target in ("0.66", "0.70"):
        p = [e["per_scale"][s]["power"][target]["estimated_power"] for s in ("1", "2", "4")]
        assert p[0] <= p[1] <= p[2], f"power is not monotone in cohort size at AUC {target}"


@ran_e
def test_the_event_count_status_follows_the_frozen_rule(e):
    p1 = e["per_scale"]["1"]["power"]["0.66"]["estimated_power"]
    reached = any(e["per_scale"][s]["power"]["0.66"]["estimated_power"] >= 0.80
                  for s in ("2", "4"))
    expected = ("NOT_SUPPORTED" if p1 >= 0.80
                else "SUPPORTED" if (p1 < 0.50 and reached) else "UNRESOLVED")
    assert e["WITHIN_R1_EVENT_COUNT_LIMITATION"] == expected


@ran_e
def test_biological_replication_is_a_design_fact_not_an_estimate(e):
    """V2 §9.5.2 -- SUPPORTED while the claim rests on one replicate; UNRESOLVED is not available."""
    assert e["BIOLOGICAL_REPLICATION_LIMITATION"] in {"SUPPORTED", "NOT_SUPPORTED"}
    assert e["BIOLOGICAL_REPLICATION_LIMITATION"] != "UNRESOLVED"
    assert e["BIOLOGICAL_REPLICATION_LIMITATION"] == "SUPPORTED"
    note = e["biological_replication_note"]
    assert "design fact" in note and "n_biological_" in note
    # a high power curve must NOT clear it
    best = max(e["per_scale"][s]["power"][t]["estimated_power"]
               for s in ("1", "2", "4") for t in ("0.66", "0.70"))
    assert best >= 0.80 and e["BIOLOGICAL_REPLICATION_LIMITATION"] == "SUPPORTED", (
        "a power curve reaching 0.80+ must not clear the replication limitation")


@ran_e
def test_no_bias_direction_is_claimed_for_the_resampled_scales(e):
    """The V2 §9.6 wording correction, tested literally."""
    caveat = e["bias_direction_caveat"]
    assert "bias direction is not guaranteed" in caveat
    assert "empirical resamples of biological replicate R1" in caveat
    assert "does not estimate power gained from additional independent biological" in caveat
    blob = json.dumps(e).lower()
    for banned in ("optimistic", "conservative"):
        assert banned not in blob, f"23.2E claimed a bias direction: {banned}"


@ran_e
def test_the_stated_limits_of_the_study_are_recorded(e):
    cannot = " ".join(e["what_this_cannot_do"]).lower()
    for required in ("create biological replicate diversity",
                     "estimate between-replicate variance",
                     "independent-sample power curve"):
        assert required in cannot, required


def test_the_synthetic_direction_never_sees_the_outcome():
    """Source-level: `z` is built from expression residualised on Bdepth, never from y_primed."""
    src = SRC.read_text(encoding="utf-8")
    body = src.split("def synthetic_direction(")[1].split("\ndef ")[0]
    for banned in ("y_primed", "y_prim", "outcome", "_assign_positives", "roc_auc"):
        assert banned not in body, f"the simulation generator touched {banned}"
    assert "BDEPTH_TABLE" in body, "the direction must be residualised on Bdepth"
    assert "PCA" in body


@ran_e
def test_23_2e_pins_the_frozen_protocol(e):
    proto = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    assert e["stage23_2_protocol_sha256"] == proto["canonical_sha256"]
    assert e["design"]["cohort_scales"] == proto["protocol"]["power"]["cohort_scales"]
    assert proto["protocol"]["power"]["n_biological_replicates"] == 1


# ============================================================================================== #
# 23.2F — diagnostic synthesis + confirmatory protocol freeze.
#
# Synthesis is where a stage that found real mechanisms is most tempted to over-deliver. Three
# things must hold and each would fail silently:
#
#   * the corrected hypothesis must be a MECHANICAL consequence of the ledger. Picking whichever
#     correction produced the biggest same-data effect is precisely the failure mode 23.2 exists
#     to prevent.
#   * no same-data result may become a confirmation, however good it looks.
#   * the confirmation protocol must be frozen BEFORE any untouched evidence is inspected, and must
#     forbid the evidence already consumed.
# ============================================================================================== #
F_RESULTS = OUT / "stage23_2_diagnostic_synthesis.json"
HANDOFF_JSON = OUT / "stage23_2_handoff_to_stage24.json"
PLANS = ROOT / "plans" / "(newer)practical plans"
HANDOFF_MD = PLANS / "STAGE_23_2_HANDOFF_TO_STAGE_24.md"
CONFIRMATION_MD = PLANS / "STAGE_23_2_ROLE_A_CONFIRMATION_V1.md"
ran_f = pytest.mark.skipif(not F_RESULTS.exists(), reason="23.2F has not been run")


@pytest.fixture(scope="module")
def f():
    return json.loads(F_RESULTS.read_text(encoding="utf-8"))


def test_23_2f_fits_nothing():
    import ast

    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    names = {"run_23_2f", "_ledger", "_decomposition_table", "_corrected_hypothesis",
             "_minimum_design_requirement", "_write_handoff", "_write_confirmation_protocol"}
    for fn in (n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name in names):
        for call in (c for c in ast.walk(fn) if isinstance(c, ast.Call)):
            if isinstance(call.func, ast.Attribute):
                assert call.func.attr not in {"fit", "fit_transform", "predict", "predict_proba"}
            if isinstance(call.func, ast.Name):
                assert not call.func.id.startswith("_fit_")
                assert call.func.id not in {"_delta_ap_once", "_load_rewind_x", "expression_block"}


@ran_f
def test_the_ledger_is_recomputable_from_the_substage_artifacts(f):
    """Re-derive all six statuses without reading 23.2F's own conclusions."""
    b = json.loads((OUT / "stage23_2_model_selection_decomposition.json").read_text())
    c = json.loads((OUT / "stage23_2_depth_decomposition.json").read_text())
    d = json.loads((OUT / "stage23_2_label_reliability.json").read_text())
    e = json.loads((OUT / "stage23_2_power_identifiability.json").read_text())
    led = f["diagnostic_ledger"]
    assert led["MODEL_SELECTION_NULL_INFLATION"] == b["MODEL_SELECTION_NULL_INFLATION"]
    assert led["RESIDUAL_DEPTH_STRUCTURE"] == c["RESIDUAL_DEPTH_STRUCTURE"]
    assert led["OUTCOME_LABEL_LIMITATION"] == d["OUTCOME_LABEL_LIMITATION"]
    assert led["WITHIN_R1_EVENT_COUNT_LIMITATION"] == e["WITHIN_R1_EVENT_COUNT_LIMITATION"]
    assert led["BIOLOGICAL_REPLICATION_LIMITATION"] == e["BIOLOGICAL_REPLICATION_LIMITATION"]

    corrected = c["corrected_same_data_diagnostic"]["CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC"]
    label = d["OUTCOME_LABEL_LIMITATION"]
    events = e["WITHIN_R1_EVENT_COUNT_LIMITATION"]
    expected = ("SUPPORTED" if (corrected == "POSITIVE" and label != "SUPPORTED")
                else "NOT_SUPPORTED" if (corrected == "NEGATIVE" and label == "NOT_SUPPORTED"
                                         and events == "NOT_SUPPORTED")
                else "UNRESOLVED")
    assert led["ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA"] == expected


@ran_f
def test_no_biological_signal_cannot_be_concluded_from_this_data(f):
    """The intended consequence of the §8.6 and §9.5.2 asymmetries."""
    led = f["diagnostic_ledger"]
    if led["OUTCOME_LABEL_LIMITATION"] != "NOT_SUPPORTED":
        assert led["ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA"] != "NOT_SUPPORTED", (
            "'no robust signal' was concluded while the label limitation is unresolved")
    assert led["_robust_derivation"]["not_supported_reachable"] is False


@ran_f
def test_the_corrected_hypothesis_is_mechanical_not_chosen_by_effect_size(f):
    ch = f["corrected_hypothesis"]
    led = f["diagnostic_ledger"]
    indicated = []
    if led["RESIDUAL_DEPTH_STRUCTURE"] == "SUPPORTED":
        indicated.append("depth_complete_nuisance_control")
    if led["MODEL_SELECTION_NULL_INFLATION"] == "SUPPORTED":
        indicated.append("search_matched_pipeline")
    if led["OUTCOME_LABEL_LIMITATION"] == "SUPPORTED":
        indicated.append("alternative_outcome_representation")
    assert ch["indicated_corrections"] == indicated
    assert ch["status"] == ("ONE_CORRECTION_INDICATED" if len(indicated) == 1
                            else "NO_CORRECTION_INDICATED" if not indicated
                            else "MULTIPLE_CORRECTIONS_EQUALLY_PLAUSIBLE")
    if ch["status"] != "ONE_CORRECTION_INDICATED":
        assert ch["hypothesis"] is None
        assert f["projected_23_2G_exit"] == "ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE"
    # a correction that is not indicated must be explicitly rejected, not silently dropped
    for name, rejected in ch["rejected_because_not_indicated"].items():
        assert rejected == (name not in indicated), name


@ran_f
def test_an_unindicated_correction_is_not_smuggled_into_the_hypothesis(f):
    ch = f["corrected_hypothesis"]
    if ch["hypothesis"] is None:
        pytest.skip("no hypothesis frozen")
    led = f["diagnostic_ledger"]
    if led["MODEL_SELECTION_NULL_INFLATION"] != "SUPPORTED":
        assert "search" not in ch["hypothesis"]["name"]
        assert ch["hypothesis"]["search"].startswith("the frozen historical")
    if led["OUTCOME_LABEL_LIMITATION"] != "SUPPORTED":
        assert "unchanged" in ch["hypothesis"]["outcome"]


@ran_f
def test_no_same_data_result_can_become_a_confirmation(f):
    assert f["same_data_status"] == "MECHANISM_DIAGNOSED_ON_EXISTING_REWIND"
    assert "never ROLE_A_CONFIRMATORY_SUPPORTED" in f["no_same_data_rescue"]
    assert f["projected_23_2G_exit"] != "ROLE_A_CONFIRMATORY_SUPPORTED"
    assert "ROLE_A_CONFIRMATORY_SUPPORTED" not in json.dumps(f["diagnostic_ledger"])


@ran_f
def test_the_historical_verdict_is_restated_as_permanent(f):
    assert "ROLE_A_SIGNAL_FAIL" in f["historical_role_a_verdict"]
    assert "permanent" in f["historical_role_a_verdict"]


@ran_f
def test_the_benchmark_change_firewall_matches_the_correction(f):
    """V2 §10.7 -- a nuisance-block change is not a benchmark semantic; an outcome change is."""
    ch = f["corrected_hypothesis"]
    fw = f["material_benchmark_change_firewall"]
    assert fw["triggered"] == ch["material_benchmark_change"]
    if "alternative_outcome_representation" in ch["indicated_corrections"]:
        assert fw["triggered"] is True
    else:
        assert fw["triggered"] is False
        assert fw["changed_semantics"] == []


@ran_f
def test_the_design_requirement_uses_only_tested_rungs(f):
    dr = f["minimum_design_requirement"]
    e = json.loads((OUT / "stage23_2_power_identifiability.json").read_text())
    tested = {v["positives_per_fold"] * 5 for v in e["per_scale"].values()}
    if dr["minimum_positive_clones"] is not None:
        assert dr["minimum_positive_clones"] in tested, "an untested cohort size was interpolated"
        smallest = dr["smallest_tested_cohort_reaching_0_80"]
        assert smallest["power"] >= 0.80
        below = [v["power"] for k, v in dr["tested_ladder"].items()
                 if v["positive_clones"] < dr["minimum_positive_clones"]]
        assert all(p < 0.80 for p in below), "a smaller tested cohort already reached 0.80"
    assert "forbids extrapolating" in dr["do_not_interpolate"]


@ran_f
def test_the_handoff_is_not_stage_24_ready(f):
    assert f["stage_24_ready"] is False
    h = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    assert h["stage_24_ready"] is False
    assert "23.2G" in h["not_ready_reason"]
    assert "NOT Stage-24-ready" in HANDOFF_MD.read_text(encoding="utf-8")


@ran_f
def test_the_handoff_carries_role_b_forward_unchanged():
    h = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    syn = json.loads((RES / "stage23_final_synthesis.json").read_text(encoding="utf-8"))
    rb = h["role_b"]
    assert rb["frozen_additive_verdict"] == syn["final_verdicts"]["role_b_additive"]
    assert rb["frozen_interaction_verdict"] == syn["final_verdicts"]["role_b_interaction"]
    assert rb["strongest_frozen_simple_baseline"]["must_be_beaten_by_stage_24"] is True
    limits = " ".join(rb["treatment_level_limitations"])
    assert "Doxorubicin" in limits and "Cisplatin" in limits


@ran_f
def test_the_handoff_records_the_evidence_firewall():
    h = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    consumed = " ".join(h["role_a"]["evidence_consumed_by_diagnosis"])
    assert "GSM7092515" in consumed and "GSM7092516" in consumed
    reserved = h["role_a"]["evidence_reserved_for_confirmation"]
    assert set(reserved) >= {"GSM7092517", "GSM7092518", "GSM7092519"}
    assert not (set(reserved) & {"GSM7092515", "GSM7092516"}), "consumed evidence was reserved"
    assert h["role_a"]["replicate_status"]["n_biological_replicates"] == 1


@ran_f
def test_the_handoff_forbids_the_claims_that_would_overreach():
    h = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    forbidden = " ".join(h["role_a"]["forbidden_claims"]).lower()
    for must in ("demonstrated signal", "no biological signal", "confirms anything"):
        assert must in forbidden, must


@ran_f
def test_the_confirmation_protocol_is_frozen_and_excludes_consumed_evidence(f):
    if f["corrected_hypothesis"]["status"] != "ONE_CORRECTION_INDICATED":
        assert not CONFIRMATION_MD.exists(), "a confirmation protocol was written without a "
        return
    text = CONFIRMATION_MD.read_text(encoding="utf-8")
    assert "before any untouched confirmation evidence was inspected" in text.lower()
    assert "GSM7092515" in text and "GSM7092516" in text, "consumed GSMs must be named as forbidden"
    forbidden_block = text.split("Forbidden data")[1]
    assert "GSM7092515" in forbidden_block
    for required in ("Confirmed scientific hypothesis", "Outcome definition", "Nuisance block",
                     "Permutation / null design", "PASS threshold",
                     "Minimum positive-count", "Source-qualification criteria",
                     "Search budget", "Stage-27 firewall"):
        assert required.lower() in text.lower(), f"the protocol omits: {required}"
    assert "observed > null p95" in text and "p_perm <= 0.05" in text


@ran_f
def test_23_2f_pins_every_substage_it_synthesised(f):
    for stage, name in (("23.2A", "stage23_2a_results.json"),
                        ("23.2B", "stage23_2_model_selection_decomposition.json"),
                        ("23.2C", "stage23_2_depth_decomposition.json"),
                        ("23.2D", "stage23_2_label_reliability.json"),
                        ("23.2E", "stage23_2_power_identifiability.json")):
        assert f["source_artifacts"][stage] == S232.sha256_file(OUT / name), stage
    assert f["models_fitted_in_23_2f"] == 0
    assert f["synthesis_is_mechanical"] is True
