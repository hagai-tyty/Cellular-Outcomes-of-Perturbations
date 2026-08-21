"""Stage 21D — reconstruction under the AUTHORS' rules (revision 2).

Revision 1 froze both datasets at `RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE` because the author
code was not on disk. It now is, so these tests pin the transcribed rules and the regressions that
must not creep back:

1. **The rule must be the authors', not a plausible substitute.** Rewind's primed set is the top 100
   gDNA barcodes by summed count, not "present at any count". The two give 42 cells and 102 cells
   respectively, so the difference is not cosmetic.
2. **`slice_max` semantics.** dplyr's default `with_ties = TRUE` returns MORE than n rows when a tie
   straddles the boundary. There IS such a tie here. Silently using `head(n)` is a different rule.
3. **42 is an assertion, never a target.** `n = 100` comes from the author script; nothing was
   tuned. The test asserts the reproduction, and separately asserts the run declares itself untuned.
4. **The superseded exploratory rule must stay superseded.** Revision 1's dominant-fraction / UMI
   floor sweep gave 4,018 naive clones; the author posterior gives 1,401. The floor sweep must not
   return as the production clone call.
5. **Earlier valid findings must not regress**: the transposed `SampleNum` mapping, the clone as the
   outer split unit, and the measured lane sizes rather than the unsupported 6479/6189 figures.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage21d_public_reconstruction.py"
spec = importlib.util.spec_from_file_location("s21d", SRC)
s21d = importlib.util.module_from_spec(spec)
sys.modules["s21d"] = s21d          # @dataclass resolves through sys.modules
spec.loader.exec_module(s21d)

RESULTS = ROOT / "results" / "diag_stage21d_public_reconstruction_results.json"
STAGE_21A = ROOT / "results" / "diag_stage21_data_audit_results.json"
STAGE_21B = ROOT / "results" / "diag_stage21b_source_design_results.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="21D has not been run")

P, AB, UNK = s21d.PRESENT, s21d.ABSENT_PROVEN, s21d.UNKNOWN


@pytest.fixture(scope="module")
def rec():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


# ---- dplyr slice_max semantics, tested on constructed data --------------------------------- #
def test_slice_max_keeps_boundary_ties_like_dplyr():
    """`with_ties = TRUE` returns n+1 rows when two values share the n-th place."""
    df = pd.DataFrame({"bc": list("abcde"), "nUMI": [10, 9, 8, 8, 1]})
    kept, cutoff, at_cut = s21d.slice_max_with_ties(df, "nUMI", 3)
    assert cutoff == 8 and at_cut == 2
    assert sorted(kept["bc"]) == ["a", "b", "c", "d"], "the tied 4th row must be kept, not dropped"


def test_slice_max_returns_everything_when_n_exceeds_the_table():
    df = pd.DataFrame({"bc": ["a", "b"], "nUMI": [3, 1]})
    kept, cutoff, at_cut = s21d.slice_max_with_ties(df, "nUMI", 100)
    assert len(kept) == 2 and cutoff is None and at_cut == 2


def test_the_top_n_is_the_author_constant_not_a_tunable():
    assert s21d.TOP_N_GDNA == 100


# ---- the named stop condition ----------------------------------------------------------------- #
def test_a_missing_gdna_arm_stops_rewind_and_names_that_exact_file(tmp_path):
    r = s21d.audit_rewind(base=tmp_path)
    assert r["verdict"] == s21d.MISSING_FILE
    f = r["findings"]["future_outcome_source"]
    assert s21d.GDNA_FILE in f.value and f.status == AB
    assert "No other" in f.evidence and "substituted" in f.evidence


def test_the_gdna_filename_is_pinned_so_a_rename_cannot_pass_silently():
    assert s21d.GDNA_FILE == "stepThreeStarcodeShavedReads_BC_gDNA.txt"


def test_a_missing_wm989_tree_is_reported_not_guessed(tmp_path):
    r = s21d.audit_gse279162(base=tmp_path)
    assert r["verdict"] == s21d.MISSING_FILE
    assert r["findings"]["required_files"].value["missing"]


def test_without_the_author_code_the_verdict_falls_back_to_pending(monkeypatch, tmp_path):
    """The whole reason revision 1 was PENDING. Removing the author scripts must reproduce that
    state rather than silently reverting to a guessed rule."""
    monkeypatch.setattr(s21d, "REWIND_SCRIPTS", [tmp_path / "not_here.R"])
    monkeypatch.setattr(s21d, "WM989_SCRIPTS", [tmp_path / "not_here.Rmd"])
    for fn in (s21d.audit_rewind, s21d.audit_gse279162):
        r = fn()
        if r["verdict"] == s21d.MISSING_FILE:
            pytest.skip("raw data not on this machine")
        assert r["verdict"] == s21d.RECON_PENDING
        assert r["findings"]["author_code"].status == UNK


def test_the_verdict_constants_are_distinct():
    assert len({s21d.RECON_PASS, s21d.RECON_PENDING, s21d.INVALID_LINKAGE,
                s21d.INSUFFICIENT_UNITS, s21d.MISSING_FILE}) == 5
    assert len({s21d.STAGE22_READY, s21d.STAGE22_PENDING, s21d.STAGE22_BLOCKED}) == 3


# ---- evidence standard -------------------------------------------------------------------------#
def test_a_finding_cannot_carry_an_invalid_status():
    with pytest.raises(ValueError):
        s21d.Finding(value=1, status="probably", evidence="")


@has_results
def test_every_finding_carries_evidence_and_a_legal_status(rec):
    for gse in ("GSE227151", "GSE279162"):
        for k, f in rec[gse]["findings"].items():
            assert f["evidence"].strip(), f"{gse}/{k}"
            assert f["status"] in (P, AB, UNK), f"{gse}/{k}"


# ---- Role A: the author rule ------------------------------------------------------------------ #
@has_results
def test_the_author_code_is_present_and_named(rec):
    f = rec["GSE227151"]["findings"]["author_code"]
    assert f["status"] == P and not f["value"]["missing"]
    names = " ".join(f["value"]["scripts"])
    assert "20220921_R1_primedVersusNonPrimedMarkersAndDistribution.R" in names
    assert "2022.02.14_R1_cellNumberDistributionForPrimedVersusNonPrimed.R" in names


@has_results
def test_our_gdna_file_is_the_authors_dummy_arm_and_the_schema_gap_is_declared(rec):
    """The equivalence the rule depends on. The filter is inert because every row IS the gDNA arm;
    the count column name differs, and that difference is recorded rather than papered over."""
    v = rec["GSE227151"]["findings"]["gdna_arm_equivalence"]["value"]
    assert v["cellID_values"] == ["dummy"]
    assert v["rows_kept_by_filter_dummy"] == v["rows_total"], "filter(cellID=='dummy') is inert"
    assert v["samplenum_values"] == [3]
    assert v["count_column_here"] == ["counts"] and v["count_column_in_author_script"] == ["UMI"]
    assert v["schema_difference"] is True, "the rename must be declared, not silently accepted"


@has_results
def test_slice_max_selected_more_than_100_barcodes_because_of_a_real_tie(rec):
    v = rec["GSE227151"]["findings"]["author_rule_slice_max"]["value"]
    assert v["n"] == 100
    assert v["rows_at_cutoff"] == 2 and v["cutoff_nUMI"] == 2365
    assert v["selected_barcodes"] == 101, "with_ties=TRUE must keep the boundary tie"
    assert v["selected_nUMI_min"] == v["cutoff_nUMI"]


@has_results
def test_the_published_42_is_reproduced_and_declared_untuned(rec):
    f = rec["GSE227151"]["findings"]["published_42_reproduction"]
    assert f["value"]["reconstructed_primed_cells"] == 42
    assert f["value"]["published_primed_cells"] == 42
    assert f["value"]["reproduced"] is True
    assert f["value"]["tuned"] is False
    assert "never targeted" in f["evidence"]


@has_results
def test_the_prospective_label_counts_are_the_author_rules_not_revision_ones(rec):
    v = rec["GSE227151"]["findings"]["prospective_label_counts"]["value"]
    assert v["primed_cells"] == 42 and v["primed_clones"] == 35
    assert v["nonprimed_cells"] == 3879 and v["nonprimed_clones"] == 3114
    assert v["primed_cells"] + v["nonprimed_cells"] == 3921
    assert v["per_sample"]["1"] == {"primed_cells": 24, "primed_clones": 21}
    assert v["per_sample"]["2"] == {"primed_cells": 18, "primed_clones": 17}
    assert v["primed_cells"] != 102, "the superseded presence-at-any-count reading must not return"


@has_results
def test_the_second_scripts_barcode_exclusion_is_inert_and_explained(rec):
    v = rec["GSE227151"]["findings"]["second_script_barcode_exclusion"]["value"]
    assert v["barcode"] == s21d.REWIND_EXCLUDED_BC
    assert v["rows_in_filtered10XCells"] == 0 and v["inert"] is True
    assert v["in_top_selection"] is False


@has_results
def test_the_author_qc_step_is_exactly_that_one_barcode(rec):
    """Revision 1 reported '21.2% of rows dropped' without saying why. It is one barcode."""
    v = rec["GSE227151"]["findings"]["author_qc_step_is_a_single_barcode_exclusion"]["value"]
    assert v["filtered_equals_stepThree_minus_that_barcode"] is True
    assert v["stepThree_rows"] - v["rows_carrying_excluded_barcode"] == v["filtered10XCells_rows"]


# ---- Role A: findings that must NOT regress ----------------------------------------------------#
@has_results
def test_samplenum_is_still_resolved_by_intersection_and_still_transposed(rec):
    f = rec["GSE227151"]["findings"]["samplenum_to_gsm"]
    assert f["value"]["mapping"] == {"1": "GSM7092516", "2": "GSM7092515"}
    cont = f["value"]["containment"]
    assert cont["1"]["GSM7092516"] == 1.0 and cont["1"]["GSM7092515"] < 0.01
    assert cont["2"]["GSM7092515"] == 1.0 and cont["2"]["GSM7092516"] < 0.01
    assert "TRANSPOSED" in f["evidence"]


@has_results
def test_the_measured_lane_sizes_stand_and_the_unsupported_figures_are_not_restored(rec):
    """6479/6189 match nothing in these files. The measured numbers are 7096/6569 raw and
    1878/2035 clone-assigned."""
    e = rec["GSE227151"]["findings"]["expression_cells"]["value"]
    assert e["GSM7092516"]["cells"] == 7096 and e["GSM7092515"]["cells"] == 6569
    a = rec["GSE227151"]["findings"]["cells_with_author_qc_clone"]["value"]
    assert a == {"1": 1878, "2": 2035}
    for bad in (6479, 6189):
        assert bad not in (e["GSM7092516"]["cells"], e["GSM7092515"]["cells"])


@has_results
def test_the_clone_is_still_the_outer_split_unit(rec):
    c = rec["GSE227151"]["findings"]["clone_structure"]["value"]
    assert c["clones_total"] == 3149 and c["shared"] == 311
    assert c["clones_sample1"] + c["clones_sample2"] - c["shared"] == c["clones_total"]
    assert "leak" in rec["GSE227151"]["findings"]["clone_structure"]["evidence"]


@has_results
def test_only_one_biological_replicate_is_locally_reconstructable(rec):
    f = rec["GSE227151"]["findings"]["replicate_structure"]["value"]
    assert f["biological_replicates_local"] == 1
    assert set(f["gsms_local"]) == {"GSM7092515", "GSM7092516"}


# ---- Role B: the Schaff pipeline ---------------------------------------------------------------#
@has_results
def test_the_five_primary_author_files_are_the_provenance_chain(rec):
    v = rec["GSE279162"]["findings"]["author_code"]["value"]
    assert not v["missing"]
    assert v["readme_order"] == ["preprocess_GEX.Rmd", "preprocess_cDNA_BCs.Rmd",
                                 "preprocess_gDNA_BCs.Rmd",
                                 "Find_Markers_Top_Res_lins_in_naive.Rmd"]
    assert "preprocess_gDNA_BCs.Rmd is separate" in rec["GSE279162"]["findings"]["author_code"][
        "evidence"], "gDNA RPM logic must not be confused with the scRNA cell call"


def test_the_qc_filters_match_the_author_script_values():
    """Transcribed from preprocess_GEX.Rmd; a typo here would silently change every count."""
    assert s21d.WM989_QC == {
        "Naive1": (3000, 75000, 20), "Naive2": (3000, 60000, 20), "Naive3": (2500, 50000, 15),
        "Dabrafenib": (2000, 30000, 15), "Trametinib": (1500, 20000, 15),
        "CoCl2": (2500, 50000, 15), "Acid": (1500, 30000, 15),
        "Cisplatin": (1500, 20000, 20), "Doxorubicin": (1500, 20000, 20)}


def test_the_barcode_calling_constants_match_the_author_script():
    assert s21d.CELL_LOWER_LIMIT == 100
    assert s21d.COR_THRESHOLD == 0.55
    assert s21d.DIFFERENCE_VAL == 0.2
    assert s21d.POSTERIOR_FLOOR == 0.5


@has_results
def test_the_author_rna_qc_actually_removed_cells(rec):
    v = rec["GSE279162"]["findings"]["author_rna_qc"]["value"]
    assert set(v) == set(s21d.WM989_QC)
    for name, row in v.items():
        assert row["cells_post_qc"] < row["cells_raw"], name
        assert row["qc"]["nFeature_RNA_gt"] == s21d.WM989_QC[name][0]
    assert sum(r["cells_post_qc"] for r in v.values()) == 46891
    assert sum(r["cells_raw"] for r in v.values()) == 77417


@has_results
def test_the_author_clustering_precondition_held(rec):
    """The author's loop calls stop('Merging happening') if a pair would join two clusters. If that
    fires, our object differs from theirs and the reconstruction is not faithful."""
    v = rec["GSE279162"]["findings"]["barcode_clustering"]["value"]
    assert v["author_stop_merging_happening_fired"] == 0
    assert v["cell_lower_limit"] == 100 and v["cor_threshold"] == 0.55
    assert v["candidate_lineages"] == 604 and v["correlated_pairs"] == 153
    assert v["clusters"] == 90 and v["lineages_merged"] == 204
    assert v["lineages_after_rowsum_filter"] == 14883


@has_results
def test_the_posterior_assignment_reproduces_and_the_floor_bites(rec):
    v = rec["GSE279162"]["findings"]["barcoding_posterior_and_assignment"]["value"]
    assert v["lineages"] == 14769 == 14883 - 204 + 90
    assert v["cells"] == 46891
    assert v["assigned_cells"] == 42771 and v["na_cells"] == 4120
    assert v["dropped_by_posterior_floor"] == 4007
    assert v["assigned_before_posterior_floor"] > v["assigned_cells"]


@has_results
def test_the_naive_lanes_are_pooled_as_the_author_pools_them(rec):
    v = rec["GSE279162"]["findings"]["per_condition_assignment"]["value"]
    assert set(v) == {"naive", *s21d.WM989_TREATMENTS}
    assert "Naive1" not in v and "Naive2" not in v and "Naive3" not in v
    assert v["naive"]["cells_post_qc"] == 7226 and v["naive"]["unique_clones"] == 1401


@has_results
def test_the_superseded_floor_sweep_numbers_do_not_come_back(rec):
    v = rec["GSE279162"]["findings"]["clone_coverage"]["value"]
    assert v["clones_with_naive_observation"] == 1401
    assert v["clones_in_ge2_treatments"] == 603
    assert v["clones_in_all_6"] == 37
    old = rec["supersedes"]["GSE279162"]
    assert v["clones_with_naive_observation"] != old["naive_clones_floor1"]
    assert v["clones_in_ge2_treatments"] != old["clones_in_ge2_treatments"]
    assert "SUPERSEDED" in old["rule_used"]
    assert sum(v["by_n_treatments"].values()) == v["clones_with_naive_observation"]


@has_results
def test_the_interaction_structure_survives_the_stricter_rule(rec):
    """The reason GSE279162 was qualified for Role B at all. It must survive the real rule, not
    only the permissive one."""
    v = rec["GSE279162"]["findings"]["clone_coverage"]["value"]
    assert v["clones_in_ge2_treatments"] >= 500
    assert v["clones_in_all_6"] > 0


@has_results
def test_the_future_outcome_is_abundance_and_num_lin_5_was_not_adopted(rec):
    f = rec["GSE279162"]["findings"]["future_outcome"]
    assert "abundance" in f["value"]["definition"]
    assert "rank" in f["value"]["definition"]
    assert set(f["value"]["per_treatment"]) == set(s21d.WM989_TREATMENTS)
    for t, row in f["value"]["per_treatment"].items():
        assert row["clones_present"] + row["clones_absent"] == 1401, t
    assert "NOT turned into a binary target" in f["evidence"]
    # The figure-specific `num_lin = 5` must be discussed, never implemented: no such binding
    # exists, and the outcome carries no top-k column.
    import ast

    bound = {t.id for node in ast.walk(ast.parse(SRC.read_text(encoding="utf-8")))
             if isinstance(node, ast.Assign) for t in node.targets if isinstance(t, ast.Name)}
    assert not {n for n in bound if "num_lin" in n.lower() or "top_res" in n.lower()}
    assert not [k for k in f["value"] if "top" in k.lower()]


@has_results
def test_the_missing_validation_anchor_is_recorded_as_unknown_not_as_success(rec):
    f = rec["GSE279162"]["findings"]["reproduction_anchor"]
    assert f["status"] == UNK
    assert f["value"]["published_count_available"] is False
    assert f["value"]["is_a_reproduction_failure"] is False
    assert "not a reproduction" in f["evidence"]


@has_results
def test_deviations_from_the_author_object_are_stated(rec):
    f = rec["GSE279162"]["findings"]["known_deviations_from_the_author_object"]
    assert len(f["value"]) >= 3


# ---- verdicts and the audit trail --------------------------------------------------------------#
@has_results
def test_both_datasets_now_pass(rec):
    for gse in ("GSE227151", "GSE279162"):
        assert rec[gse]["verdict"] == s21d.RECON_PASS, gse
        assert rec[gse]["previous_verdict"] == s21d.RECON_PENDING, gse


@has_results
def test_the_overall_gate_is_derived_from_the_two_dataset_verdicts(rec):
    both_pass = all(rec[g]["verdict"] == s21d.RECON_PASS for g in ("GSE227151", "GSE279162"))
    assert rec["overall"] == (s21d.STAGE22_READY if both_pass else s21d.STAGE22_PENDING)
    assert rec["overall"] == "STAGE_22_READY"
    assert rec["previous_overall"] == s21d.STAGE22_PENDING


@has_results
def test_revision_one_is_preserved_rather_than_erased(rec):
    s = rec["supersedes"]
    assert rec["revision"] == 2 and s["revision"] == 1
    assert s["commit"] == "30ca7f0"
    assert s["overall"] == s21d.STAGE22_PENDING
    assert s["GSE227151"]["positive_cells"] == 102 and s["GSE227151"]["positive_clones"] == 82
    assert s["GSE279162"]["naive_clones_floor1"] == 4018
    assert s["why_it_was_reasonable_then"].strip()
    assert len(s["findings_that_still_stand"]) >= 3


# ---- additive contract and engineering constraints ----------------------------------------------#
@has_results
def test_stage_21d_did_not_touch_the_earlier_frozen_results(rec):
    assert STAGE_21A.exists() and STAGE_21B.exists()
    a = json.loads(STAGE_21A.read_text(encoding="utf-8"))
    assert a["verdict"] == "CULTURE_FORWARD_AVAILABLE"
    b = json.loads(STAGE_21B.read_text(encoding="utf-8"))
    assert b["GSE242423"]["verdict"] == "LINEAGE_ABSENT_PROVEN"
    assert b["GSE165176"]["verdict"] == "ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY"
    assert rec["additive_to"] == ["results/diag_stage21_data_audit_results.json",
                                  "results/diag_stage21b_source_design_results.json"]


@has_results
def test_no_model_was_fitted(rec):
    """Checked against the parsed code, not the prose: the docstring may SAY 'no sklearn'."""
    import ast

    assert rec["model_fitted"] is False and rec["src_modified"] is False
    src = SRC.read_text(encoding="utf-8")
    imported = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not (imported & {"sklearn", "torch", "tensorflow", "xgboost", "lightgbm", "statsmodels",
                            "cellfate"}), f"21D must not model or import src/: {sorted(imported)}"
    body = src[src.index('"""', src.index('"""') + 3) + 3:]
    for banned in ("LogisticRegression", "sklearn", "torch", ".fit("):
        assert banned not in body, f"21D must not model: found {banned} in code"


@has_results
def test_neither_the_raw_data_nor_the_author_code_entered_git(rec):
    assert rec["raw_data_committed"] is False and rec["author_code_committed"] is False
    for gse in ("GSE227151", "GSE279162"):
        entries = rec[gse]["manifest"]
        assert entries, gse
        for e in entries:
            assert e["exists"] and e["sha256"] and e["bytes"]
            assert not Path(e["path"]).is_relative_to(ROOT), f"{e['path']} is inside the repo"
    codes = [e for g in ("GSE227151", "GSE279162") for e in rec[g]["manifest"]
             if "author_code" in e["path"]]
    assert len(codes) == 7, "both author-code sets are hashed as provenance"


def test_the_script_writes_only_its_own_results_and_the_two_small_tables():
    src = SRC.read_text(encoding="utf-8")
    assert src.count(".write_text(") == 2 and src.count(".to_csv(WM989_TABLE") == 1
    for name in ("diag_stage21d_public_reconstruction_results.json",
                 "stage21d_rewind_clone_table.tsv",
                 "stage21d_gse279162_clone_table.tsv"):
        assert name in src


@has_results
def test_the_committed_clone_tables_are_small_and_carry_the_author_labels(rec):
    rw = pd.read_csv(ROOT / rec["GSE227151"]["clone_table"], sep="\t")
    assert (ROOT / rec["GSE227151"]["clone_table"]).stat().st_size < 400_000
    assert int(rw["primed"].sum()) == 35, "clone-level positives"
    assert int(rw.loc[rw["primed"] == 1, "cells"].sum()) == 42
    assert int(rw["in_gdna_at_all"].sum()) == 82, "revision 1's set is kept as a column, not a label"

    wm = pd.read_csv(ROOT / rec["GSE279162"]["clone_table"], sep="\t")
    assert (ROOT / rec["GSE279162"]["clone_table"]).stat().st_size < 400_000
    assert list(wm.columns) == ["lineage", "naive", *s21d.WM989_TREATMENTS]
    assert int((wm["naive"] > 0).sum()) == 1401
    assert np.issubdtype(wm[list(s21d.WM989_TREATMENTS)].to_numpy().dtype, np.integer)
