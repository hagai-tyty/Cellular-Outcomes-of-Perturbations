"""Stage 22 — the frozen benchmark contracts.

Stage 22's whole job is to make the task un-redefinable after Stage 23 sees performance. These
tests are what makes that stick. They pin the five ways the benchmark could quietly stop being the
task it claims to be:

1. **An ambiguous grouping unit.** 8 Rewind cells carry two clone assignments, and both of their
   clones land in different folds. If they came back, `outer_fold` would be undefined for a cell
   whose expression is one column, and clone information would cross the split.
2. **A different label than the authors'.** The primed set must remain the top-100-with-ties gDNA
   cut — 42 cells / 35 clones — not a read floor, not a presence test.
3. **A benchmark that becomes "among surviving clones, predict abundance."** Every eligible WM989
   clone keeps all six treatment rows, and a zero stays an `observed_zero` against an available
   sample rather than being dropped or relabelled as death.
4. **A shortcut promoted to the headline.** Clone size predicts the WM989 outcome strongly, so the
   feature-eligibility firewall has to name it as nuisance, not as `X`.
5. **A split that drifts.** Folds are frozen here, once, on sorted clone keys — not left to Stage 23
   to reseed until something looks good.

The tests read the committed artifacts; they do not re-run the ~2 minute build.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "build_stage22_prospective_benchmarks.py"
spec = importlib.util.spec_from_file_location("s22", SRC)
S22 = importlib.util.module_from_spec(spec)
sys.modules["s22"] = S22
spec.loader.exec_module(S22)

RES = ROOT / "results"
FILES = {
    "results": RES / "stage22_prospective_benchmark_results.json",
    "rewind_manifest": RES / "stage22_rewind_benchmark_manifest.json",
    "wm989_manifest": RES / "stage22_wm989_benchmark_manifest.json",
    "rewind_cells": RES / "stage22_rewind_cells.csv",
    "rewind_clones": RES / "stage22_rewind_clones.csv",
    "wm989_cells": RES / "stage22_wm989_cell_assignments.csv",
    "wm989_naive": RES / "stage22_wm989_naive_cells.csv",
    "wm989_clones": RES / "stage22_wm989_clones.csv",
    "wm989_ct": RES / "stage22_wm989_clone_treatment.csv",
}
built = pytest.mark.skipif(not all(p.exists() for p in FILES.values()),
                           reason="Stage 22 has not been built")


@pytest.fixture(scope="module")
def art():
    out = {}
    for k, p in FILES.items():
        out[k] = (json.loads(p.read_text(encoding="utf-8")) if p.suffix == ".json"
                  else pd.read_csv(p))
    return out


# ---- constants and helpers, no data required --------------------------------------------------#
def test_the_split_convention_is_frozen_in_code():
    assert S22.N_SPLITS == 5
    assert S22.STAGE22_SPLIT_SEED == 22022


def test_the_author_rule_constants_come_from_stage_21d_not_a_copy():
    """Plan §1.3: one source of truth. A second implementation that merely reproduces the headline
    counts could drift from the source-faithful logic without any test noticing."""
    assert S22.S21D.TOP_N_GDNA == 100
    assert S22.S21D.CELL_LOWER_LIMIT == 100
    assert S22.S21D.COR_THRESHOLD == 0.55
    assert S22.S21D.DIFFERENCE_VAL == 0.2
    assert S22.S21D.POSTERIOR_FLOOR == 0.5
    src = SRC.read_text(encoding="utf-8")
    for fn in ("qc_and_lineage", "barcode_clustering", "barcode_combine",
               "barcoding_posterior_and_assignment", "slice_max_with_ties"):
        assert f"S21D.{fn}" in src, f"{fn} must be imported from 21D, not reimplemented"
    defined = {n.name for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)}
    assert not (defined & {"qc_and_lineage", "barcode_clustering", "barcode_combine",
                           "barcoding_posterior_and_assignment", "slice_max_with_ties"})


def test_folds_are_deterministic_and_independent_of_input_order():
    keys = [f"C{i:04d}" for i in range(97)]
    a = S22.deterministic_folds(keys)
    b = S22.deterministic_folds(list(reversed(keys)))
    assert a == b, "a shuffled input order must not move a single clone"
    assert set(a.values()) == set(range(S22.N_SPLITS))


def test_stratified_folds_keep_both_classes():
    keys = [f"C{i:04d}" for i in range(200)]
    y = {k: int(i < 20) for i, k in enumerate(keys)}
    folds = S22.deterministic_folds(keys, y)
    per = {}
    for k, f in folds.items():
        per.setdefault(f, []).append(y[k])
    assert all(0 < sum(v) < len(v) for v in per.values())


def test_the_builder_takes_roots_as_arguments_rather_than_hard_coding_them():
    """Plan §1.4. The 21D module may keep its D:\\ defaults; the builder must accept explicit
    roots so the benchmark is reproducible off this machine."""
    src = SRC.read_text(encoding="utf-8")
    assert "--rewind-root" in src and "--wm989-root" in src
    for fn in ("build_rewind", "build_wm989"):
        node = next(n for n in ast.walk(ast.parse(src))
                    if isinstance(n, ast.FunctionDef) and n.name == fn)
        assert [a.arg for a in node.args.args] == ["root"], f"{fn} must be root-parameterised"
    assert 'r"D:' not in src and "r'D:" not in src


def test_the_no_modelling_gate_is_computed_from_the_syntax_tree():
    """It is checked on the AST, not on the text -- a text search matches this very sentence, which
    is how the first version of the gate failed on itself."""
    ev = S22.builder_fits_no_model()
    assert ev["no_modelling_imports"] and ev["no_fit_call"]
    assert ev["sklearn_use_is_splitters_only"], "sklearn is used for splitters only"
    assert ev["modelling_imports_found"] == []


# ---- Rewind frozen contract -------------------------------------------------------------------#
@built
def test_rewind_source_audit_reproduces(art):
    s = art["rewind_manifest"]["statistics"]
    assert s["source_author_qc_assignment_records"] == 3921
    assert s["source_unique_cell_uid"] == 3913
    assert s["source_unique_clones"] == 3149
    assert s["bare_cellID_cross_lane_collisions"] == 0


@built
def test_exactly_eight_ambiguous_cells_are_enumerated_and_excluded(art):
    m = art["rewind_manifest"]
    s = m["statistics"]
    assert s["ambiguous_cell_uid_excluded"] == 8
    assert s["source_rows_removed_by_exclusion"] == 16
    ex = m["declared_exclusions"]
    assert ex["n_cell_uid"] == 8 and len(ex["cells"]) == 8
    for c in ex["cells"]:
        assert len(c["clone_ids"]) == 2, "each is a two-clone cell"
        assert c["n_source_rows"] == 2
        assert c["exclusion_reason"] == S22.AMBIGUITY_REASON
        assert c["any_clone_primed"] is False, "none touches a primed clone; positives are intact"


@built
def test_the_retained_rewind_population_is_3905_and_uniquely_keyed(art):
    cells = art["rewind_cells"]
    assert len(cells) == 3905 == 3913 - 8
    assert cells["cell_uid"].is_unique
    assert (cells["cell_uid"] == cells["SampleNum"].astype(str) + ":" + cells["cellID"]).all()
    assert cells.groupby("cell_uid")["clone_id"].nunique().max() == 1


@built
def test_the_post_exclusion_clone_counts_were_recomputed_not_carried_forward(art):
    """Plan §3.5. Two clones consisted only of excluded cells, so 3,149 is NOT the retained count.
    Carrying the Stage-21D figure forward would have been wrong by exactly those two."""
    s = art["rewind_manifest"]["statistics"]
    assert s["retained_benchmark_clones"] == 3147
    assert s["clones_lost_entirely_to_exclusion"] == 2
    assert s["negative_cells"] == 3863 and s["negative_clones"] == 3112
    assert len(art["rewind_clones"]) == 3147


@built
def test_the_author_top_100_with_ties_rule_produced_the_label(art):
    r = art["rewind_manifest"]["author_rule_identifiers"]
    assert r["top_n"] == 100 and r["with_ties"] is True
    assert "slice_max(n=100, with_ties=TRUE)" in r["outcome_rule"]
    s = art["rewind_manifest"]["statistics"]["author_rule"]
    assert s["cutoff_nUMI"] == 2365 and s["barcodes_at_cutoff"] == 2
    assert s["selected_barcodes"] == 101, "the boundary tie is kept, so 101 not 100"


@built
def test_rewind_positive_anchors(art):
    cells, clones = art["rewind_cells"], art["rewind_clones"]
    assert int(cells["y_primed"].sum()) == 42
    assert int(clones["y_primed"].sum()) == 35
    assert cells.loc[cells["y_primed"] == 1, "clone_id"].nunique() == 35
    assert int((cells["y_primed"] == 0).sum()) == 3863


@built
def test_no_rewind_clone_carries_a_contradictory_outcome(art):
    cells = art["rewind_cells"]
    assert cells.groupby("clone_id")["y_primed"].nunique().max() == 1
    clone_y = dict(zip(art["rewind_clones"]["clone_id"], art["rewind_clones"]["y_primed"],
                       strict=True))
    assert (cells["clone_id"].map(clone_y) == cells["y_primed"]).all()


@built
def test_every_rewind_cell_shares_its_clone_group_and_fold(art):
    cells, clones = art["rewind_cells"], art["rewind_clones"]
    assert (cells["outer_group"] == cells["clone_id"]).all()
    assert cells.groupby("clone_id")["outer_fold"].nunique().max() == 1
    fold = dict(zip(clones["clone_id"], clones["outer_fold"], strict=True))
    assert (cells["clone_id"].map(fold) == cells["outer_fold"]).all()


@built
def test_five_rewind_folds_all_carry_positives_and_negatives(art):
    clones = art["rewind_clones"]
    per = clones.groupby("outer_fold")["y_primed"].agg(["size", "sum"])
    assert len(per) == 5
    assert (per["sum"] > 0).all(), "35 positive clones over 5 folds -- this check is mandatory"
    assert (per["sum"] < per["size"]).all()
    assert int(per["sum"].sum()) == 35


@built
def test_rewind_expression_mapping_is_complete_and_consistent(art):
    cells = art["rewind_cells"]
    assert cells["expression_column_index"].notna().all()
    assert (cells["expression_barcode"].str.split("-").str[0] == cells["cellID"]).all()
    assert set(cells["gsm"]) == {"GSM7092515", "GSM7092516"}
    assert (cells.groupby("gsm")["expression_source"].nunique() == 1).all()
    s = art["rewind_manifest"]["statistics"]
    assert s["expression_unmapped_cells"] == 0 and s["expression_mapped_cells"] == 3905


@built
def test_samplenum_is_still_the_containment_mapping_not_the_geo_title(art):
    s = art["rewind_manifest"]["statistics"]
    assert s["samplenum_to_gsm"] == {"1": "GSM7092516", "2": "GSM7092515"}
    c = s["samplenum_containment"]
    assert c["1"]["GSM7092516"] == 1.0 and c["2"]["GSM7092515"] == 1.0


@built
def test_the_rewind_claim_scope_is_frozen_in_the_manifest(art):
    cs = art["rewind_manifest"]["claim_scope"]
    assert cs["biological_replicate"] == "R1"
    assert cs["generalization_scope"] == "within_R1_clone_heldout"
    assert "not independent" in cs["limitation"]
    assert "NOT proven biological" in cs["outcome_semantics_caveat"], \
        "nonprimed is an operational label, never asserted as failure"
    assert "not applicable" in cs["treatment_variation"]
    cells = art["rewind_cells"]
    assert (cells["biological_replicate"] == "R1").all()
    assert (cells["outcome_semantics"] == S22.REWIND_OUTCOME_SEMANTICS).all()


# ---- WM989 frozen contract --------------------------------------------------------------------#
@built
def test_the_author_qc_thresholds_are_recorded_per_condition(art):
    qc = art["wm989_manifest"]["author_rule_identifiers"]["qc"]
    assert qc["Naive1"] == [3000, 75000, 20] and qc["Doxorubicin"] == [1500, 20000, 20]
    per = art["wm989_manifest"]["statistics"]["per_sample_qc"]
    assert len(per) == 9
    for name, row in per.items():
        assert row["cells_post_qc"] < row["cells_raw"], name
    assert art["wm989_manifest"]["statistics"]["post_qc_cells"] == 46891


@built
def test_custom_features_feed_the_lineage_assay(art):
    f = art["wm989_manifest"]["statistics"]["feature_structure"]
    assert f["n_genes"] == 36601 and f["n_lineage_features"] == 153055
    assert f["lineage_assay_source"] == "Custom feature block"


@built
def test_the_clone_call_is_the_author_pipeline_with_no_dominant_fraction_fallback(art):
    r = art["wm989_manifest"]["author_rule_identifiers"]
    assert r["cell_lower_limit"] == 100 and r["cor_threshold"] == 0.55
    assert r["difference_val"] == 0.2 and r["posterior_floor"] == 0.5
    s = art["wm989_manifest"]["statistics"]
    assert "no dominant-fraction" in s["clone_call_source"]
    src = SRC.read_text(encoding="utf-8")
    body = src[src.index('"""', src.index('"""') + 3) + 3:]
    for banned in ("dominant_fraction", "umi_floor", "UMI_FLOOR"):
        assert banned not in body


@built
def test_lineage_assignment_was_run_once_on_the_joint_object(art):
    j = art["wm989_manifest"]["statistics"]["joint_assignment"]
    assert j["run_once_on_all_samples"] is True
    assert j["cells_in_joint_object"] == 46891
    assert j["lineages_in_joint_object"] == 14883
    p = art["wm989_manifest"]["statistics"]["posterior_assignment"]
    assert p["dropped_by_posterior_floor"] == 4007


@built
def test_wm989_cell_uid_is_gsm_plus_barcode_and_the_collisions_are_real(art):
    cells = art["wm989_cells"]
    assert cells["cell_uid"].is_unique
    assert (cells["cell_uid"] == cells["gsm"] + ":" + cells["cell_barcode"]).all()
    s = art["wm989_manifest"]["statistics"]
    # Two populations, both recorded. The plan quotes the raw figure; the uniqueness requirement
    # bites on the post-QC benchmark population, which is the smaller one.
    assert s["reused_bare_barcodes_across_samples_raw"] == 722, "the plan's figure, reproduced"
    assert s["reused_bare_barcodes_across_samples"] == 248, "post-QC benchmark population"
    assert s["reused_bare_barcodes_across_samples"] > 0, \
        "the compound key is required here, not merely defensive"


@built
def test_the_three_naive_lanes_are_one_pretreatment_condition(art):
    cells = art["wm989_cells"]
    assert set(cells.loc[cells["is_naive"], "sample"]) == {"Naive1", "Naive2", "Naive3"}
    assert set(cells.loc[cells["is_naive"], "condition"]) == {"naive"}
    by = art["wm989_manifest"]["statistics"]["by_condition"]
    assert set(by) == {"naive", *S22.S21D.WM989_TREATMENTS}
    assert by["naive"]["post_qc_cells"] == 7226 and by["naive"]["unique_clones"] == 1401
    tot = art["wm989_manifest"]["statistics"]["naive_totals_by_sample"]
    assert tot == {"Naive1": 1433, "Naive2": 1863, "Naive3": 3193}
    assert sum(tot.values()) == art["wm989_manifest"]["statistics"]["naive_pooled_denominator"]


@built
def test_na_lineage_cells_are_documented_and_never_reassigned(art):
    cells = art["wm989_cells"]
    na = cells[~cells["is_assigned"]]
    assert len(na) == 4120
    assert na["assigned_lineage"].isna().all(), "no heuristic backfill"
    ex = art["wm989_manifest"]["declared_exclusions"]
    assert ex["n_cells"] == 4120 and "NOT reassigned" in ex["note"]
    assert not set(na["cell_uid"]) & set(art["wm989_naive"]["cell_uid"])


@built
def test_every_eligible_clone_carries_all_six_treatment_rows(art):
    ct, clones = art["wm989_ct"], art["wm989_clones"]
    assert len(clones) == 1401
    assert len(ct) == 8406 == 1401 * 6
    assert (ct.groupby("clone_id").size() == 6).all()
    assert set(ct["treatment"]) == set(S22.S21D.WM989_TREATMENTS)
    assert set(ct["clone_id"]) == set(clones["clone_id"])


@built
def test_zero_outcomes_are_observed_zeros_not_missing_and_not_relabelled(art):
    ct = art["wm989_ct"]
    zero = ct[ct["n_post_cells"] == 0]
    assert len(zero) == 6150, "73% of rows -- dropping them would change the task"
    assert (zero["outcome_observation_status"] == "observed_zero").all()
    assert not zero["detected_post"].any()
    assert zero["treatment_sample_available"].all()
    assert (ct.loc[ct["n_post_cells"] > 0, "outcome_observation_status"] == "observed_nonzero").all()
    assert ct["n_post_cells"].notna().all(), "no zero was coerced from NA"
    for word in ("death", "failure", "sensitive", "resistant"):
        assert word not in set(ct["outcome_observation_status"])


@built
def test_all_six_treatment_samples_are_available_with_positive_denominators(art):
    per = art["wm989_manifest"]["statistics"]["per_treatment"]
    assert set(per) == set(S22.S21D.WM989_TREATMENTS)
    assert all(v["total_assigned_cells"] > 0 for v in per.values())
    ct = art["wm989_ct"]
    assert ct["treatment_sample_available"].all()
    assert (ct["treatment_total_assigned_cells"] > 0).all()


@built
def test_the_rank_convention_is_min_competition_rank_recomputed_from_the_table(art):
    ct = art["wm989_ct"]
    for t, g in ct.groupby("treatment"):
        expect = g["n_post_cells"].rank(method="min", ascending=False).astype(int)
        assert (g["post_rank"].to_numpy() == expect.to_numpy()).all(), t
        assert np.allclose(g["post_rank_fraction"], g["post_rank"] / 1401)
        tie = g["n_post_cells"].map(g["n_post_cells"].value_counts())
        assert (g["post_tie_size"].to_numpy() == tie.to_numpy()).all(), t
    assert "min" in art["wm989_manifest"]["author_rule_identifiers"]["rank_convention"]


@built
def test_no_binary_resistance_threshold_was_frozen(art):
    ct = art["wm989_ct"]
    assert not [c for c in ct.columns
                if c.startswith("y_") or "resistant" in c or c.startswith("top")]
    assert "n_post_cells" in ct.columns and "post_rank" in ct.columns
    src = SRC.read_text(encoding="utf-8")
    body = src[src.index('"""', src.index('"""') + 3) + 3:]
    assert "num_lin" not in body, "the figure-specific top-5 is discussed in the plan, not coded"


@built
def test_post_fraction_uses_the_full_treatment_denominator(art):
    ct = art["wm989_ct"]
    assert np.allclose(ct["post_fraction"],
                       ct["n_post_cells"] / ct["treatment_total_assigned_cells"], atol=1e-8)
    per = art["wm989_manifest"]["statistics"]["per_treatment"]
    sums = [v["benchmark_fraction_sum"] for v in per.values()]
    assert all(s < 1.0 for s in sums), \
        "the benchmark holds only the eligible subset, so its fractions sum to LESS than one"
    assert "sum to less than one" in art["wm989_manifest"]["claim_scope"]["compositional"]


@built
def test_wm989_clone_folds_are_consistent_across_all_three_tables(art):
    clones, ct, naive = art["wm989_clones"], art["wm989_ct"], art["wm989_naive"]
    fold = dict(zip(clones["clone_id"], clones["outer_fold"], strict=True))
    assert (ct["clone_id"].map(fold) == ct["outer_fold"]).all()
    assert (naive["clone_id"].map(fold) == naive["outer_fold"]).all()
    assert ct.groupby("clone_id")["outer_fold"].nunique().max() == 1, \
        "a clone must not have one treatment arm in train and another in test"
    assert clones["outer_fold"].nunique() == 5
    assert (clones["outer_group"] == clones["clone_id"]).all()


@built
def test_wm989_naive_expression_mapping_is_complete(art):
    naive = art["wm989_naive"]
    assert len(naive) == 6489
    assert naive["expression_column_index"].notna().all()
    assert (naive["expression_barcode"].str.split("-").str[0]
            == naive["cell_uid"].str.split(":").str[1]).all()
    assert set(naive["gsm"]) == {"GSM8562999", "GSM8563000", "GSM8563001"}
    s = art["wm989_manifest"]["statistics"]
    assert s["naive_expression_unmapped_cells"] == 0


@built
def test_the_abundance_confound_is_measured_and_declared(art):
    """Clone size predicts the outcome hard. Stage 23 must clear it before any state claim, so the
    diagnostic is recorded here rather than discovered later."""
    z = art["wm989_manifest"]["statistics"]["zero_rate_by_naive_depth"]
    assert set(z) == {"1", "2", "3-4", "5-9", "10+"}
    rates = [z[k]["zero_rate"] for k in ("1", "2", "3-4", "5-9", "10+")]
    assert rates == sorted(rates, reverse=True), "the zero rate must fall with clone depth"
    assert rates[0] - rates[-1] > 0.3, "and it does so steeply, which is the point"
    assert "abundance-only nuisance baselines" in \
        art["wm989_manifest"]["claim_scope"]["abundance_confound"]


@built
def test_wm989_clone_coverage_survives_the_author_rule(art):
    c = art["wm989_manifest"]["statistics"]["clone_coverage"]
    assert c["in_ge2_treatments"] == 603 and c["in_all_6"] == 37
    assert sum(c["by_n_treatments"].values()) == 1401


# ---- gates, verdicts, provenance ---------------------------------------------------------------#
@built
def test_all_ten_gates_pass_and_none_is_hard_coded(art):
    r = art["results"]
    assert len(r["gates"]) == 10 and all(r["gates"].values())
    assert r["all_gates_pass"] is True
    ev = r["gate_evidence"]
    assert ev["G22-2"]["rewind_target_recomputed_from_outcome_source_alone"] is True
    assert ev["G22-10"]["no_modelling_imports"] and ev["G22-10"]["no_fit_call"]
    src = SRC.read_text(encoding="utf-8")
    assert '"G22-2_no_label_leakage": True' not in src
    assert '"G22-10_no_modelling": True' not in src


@built
def test_the_target_recomputes_from_the_outcome_source_alone(art):
    """G22-2 with teeth: the primed set is rebuilt from the gDNA arm without opening an expression
    file, and must equal what the benchmark wrote."""
    ev = art["results"]["gate_evidence"]["G22-2"]
    assert ev["rewind_target_recomputed_from_outcome_source_alone"]
    assert ev["rewind_written_labels_are_a_subset_of_that_set"]
    assert "no pretreatment expression value enters any target" in ev["wm989_caveat"]


@built
def test_the_feature_firewall_names_clone_size_as_nuisance_not_as_x(art):
    fe = art["results"]["feature_eligibility"]
    assert set(fe) == {"TARGET", "PROVENANCE_ONLY", "BASELINE_NUISANCE", "PRIMARY_X"}
    for col in ("n_naive_cells", "n_pretreatment_cells", "naive_pooled_fraction",
                "treatment_total_assigned_cells"):
        assert col in fe["BASELINE_NUISANCE"], col
    for col in ("outer_fold", "clone_id", "expression_column_index", "gsm"):
        assert col in fe["PROVENANCE_ONLY"], col
    for col in ("y_primed", "n_post_cells", "post_rank"):
        assert col in fe["TARGET"], col
    assert len(fe["PRIMARY_X"]) == 1 and "expression" in fe["PRIMARY_X"][0]
    overlap = set(fe["TARGET"]) & (set(fe["PROVENANCE_ONLY"]) | set(fe["BASELINE_NUISANCE"]))
    assert not overlap, f"a column cannot be both a target and a predictor: {overlap}"


@built
def test_the_verdicts_and_the_stage_23_gate_are_derived(art):
    r = art["results"]
    for d in ("GSE227151", "GSE279162"):
        assert r["datasets"][d]["verdict"] in S22.READY_VERDICTS
    assert r["overall"] == S22.STAGE_23_READY
    role_a = r["datasets"]["GSE227151"]["verdict"]
    assert r["overall"] == (S22.STAGE_23_READY if role_a in S22.READY_VERDICTS
                            else S22.STAGE_23_BLOCKED)
    assert r["role_b_is_non_blocking"] is True
    assert r["role_b_status_preserved"] == r["datasets"]["GSE279162"]["verdict"]


@built
def test_declared_missingness_is_why_neither_is_plain_ready(art):
    """Both datasets carry enumerated exclusions, so the honest verdict is the declared variant."""
    assert art["rewind_manifest"]["verdict"] == S22.BENCHMARK_READY_DECLARED
    assert art["rewind_manifest"]["declared_exclusions"]["n_cell_uid"] == 8
    assert art["wm989_manifest"]["verdict"] == S22.BENCHMARK_READY_DECLARED
    assert art["wm989_manifest"]["declared_exclusions"]["n_cells"] == 4120


@built
def test_provenance_is_complete_and_never_circular(art):
    for key, manifest in (("rewind_manifest", art["rewind_manifest"]),
                          ("wm989_manifest", art["wm989_manifest"])):
        assert manifest["reconstruction_commit"] == "6c2f2d6"
        assert manifest["plan"]["version"] == "V2"
        assert len(manifest["plan"]["sha256"]) == 64
        assert manifest["source_files"] and manifest["author_code_files"]
        for f in manifest["source_files"] + manifest["author_code_files"]:
            assert len(f["sha256"]) == 64 and f["bytes"] > 0
            assert "\\" not in f["name"] and "/" not in f["name"], "identity by basename"
        names = {a["name"] for a in manifest["derived_artifacts"]}
        assert FILES[key].name not in names, "a manifest must not contain its own hash"
    results_names = {a["name"] for a in art["results"]["manifest_hashes"]}
    assert results_names == {FILES["rewind_manifest"].name, FILES["wm989_manifest"].name}
    assert FILES["results"].name not in results_names


@built
def test_the_manifests_carry_no_timestamp_or_repository_head(art):
    """Plan §9: benchmark identity must reproduce byte-for-byte, so mutable state stays out."""
    for m in (art["rewind_manifest"], art["wm989_manifest"], art["results"]):
        blob = json.dumps(m).lower()
        for banned in ("timestamp", "generated_at", "git_head", "\"head\"", "datetime"):
            assert banned not in blob, banned


@built
def test_no_raw_or_author_code_file_was_committed(art):
    assert art["results"]["raw_data_committed"] is False
    assert art["results"]["model_fitted"] is False
    assert art["results"]["src_modified"] is False
    for m in (art["rewind_manifest"], art["wm989_manifest"]):
        assert not (RES / m["source_files"][0]["name"]).exists()
        assert not (RES / m["author_code_files"][0]["name"]).exists()


@built
def test_the_treatment_alias_map_covers_all_three_naming_systems(art):
    a = art["wm989_manifest"]["statistics"]["treatment_aliases"]
    assert set(a) == set(S22.S21D.WM989_TREATMENTS)
    assert a["Dabrafenib"] == {"plan": "dabrafenib", "author_code": "dab"}
    assert all(set(v) == {"plan", "author_code"} for v in a.values())
