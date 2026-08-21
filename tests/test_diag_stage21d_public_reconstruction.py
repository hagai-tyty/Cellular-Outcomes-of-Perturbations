"""Stage 21D — the reconstruction that turns two qualified accessions into prospective tables.

21D is the stage where a *reconstruction* can quietly become an *invention*. These tests pin the
four ways that could happen here:

1. **Substituting a missing file.** If the future gDNA barcode arm is absent, the Rewind branch has
   to stop and name that exact file. Reaching for another intermediate would silently swap the
   outcome variable for something else.
2. **Trusting a sample label.** The Rewind barcode tables' `SampleNum` is TRANSPOSED relative to
   GEO's "sample N" wording. Believing either label instead of intersecting cellIDs would attach
   every cell's expression to the wrong lane.
3. **Tuning a threshold to a published number.** The paper reports 42 primed cells. No floor may be
   moved to produce it — and in fact no cell-level reading reproduces it, which is recorded rather
   than smoothed over.
4. **Calling a floor-conditional count a fact.** GSE279162's clone counts move with the UMI floor.
   The verdict must stay PENDING while that rule is unresolved.

Plus the tri-state rule carried over from 21A/21B: "not found" is never "absent".
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

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


# ---- the named stop condition ----------------------------------------------------------------- #
def test_a_missing_gdna_arm_stops_rewind_and_names_that_exact_file(tmp_path):
    """Case 1. The future outcome has exactly one source. If it is gone, the branch stops."""
    r = s21d.audit_rewind(base=tmp_path)
    assert r["verdict"] == s21d.MISSING_FILE
    f = r["findings"]["future_outcome_source"]
    assert s21d.GDNA_FILE in f.value
    assert f.status == AB
    assert "No other" in f.evidence and "substituted" in f.evidence


def test_the_gdna_filename_is_pinned_so_a_rename_cannot_pass_silently():
    assert s21d.GDNA_FILE == "stepThreeStarcodeShavedReads_BC_gDNA.txt"


def test_a_missing_wm989_tree_is_reported_not_guessed(tmp_path):
    r = s21d.audit_gse279162(base=tmp_path)
    assert r["verdict"] == s21d.MISSING_FILE
    assert r["findings"]["required_files"].value["missing"], "the absence must be enumerated"


def test_the_verdict_constants_are_distinct():
    assert len({s21d.RECON_PASS, s21d.RECON_PENDING, s21d.INVALID_LINKAGE,
                s21d.INSUFFICIENT_UNITS, s21d.MISSING_FILE}) == 5


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


# ---- Role A: GSE227151 Rewind ------------------------------------------------------------------#
@has_results
def test_samplenum_is_resolved_by_intersection_and_is_transposed_vs_geo(rec):
    """Case 2, and the single most damaging silent error available in this dataset."""
    f = rec["GSE227151"]["findings"]["samplenum_to_gsm"]
    assert f["value"]["mapping"] == {"1": "GSM7092516", "2": "GSM7092515"}
    cont = f["value"]["containment"]
    assert cont["1"]["GSM7092516"] == 1.0 and cont["1"]["GSM7092515"] < 0.01
    assert cont["2"]["GSM7092515"] == 1.0 and cont["2"]["GSM7092516"] < 0.01
    assert "TRANSPOSED" in f["evidence"]


@has_results
def test_the_author_qc_table_is_the_linkage_not_a_rebuilt_one(rec):
    f = rec["GSE227151"]["findings"]["author_qc_table_schema"]
    assert f["value"] == ["cellID", "BC50StarcodeD8", "SampleNum", "nUMI", "fracUMI", "nLineages"]
    assert "NOT rebuilt" in f["evidence"]
    pre = rec["GSE227151"]["findings"]["prefilter_table_rows"]["value"]
    assert pre["filtered10XCells"] < pre["stepThree_BC_10X"], "filtered must be the QC'd subset"


@has_results
def test_the_primed_label_did_not_need_a_threshold(rec):
    """The outcome is invariant across every floor up to the weakest positive clone, so no cutoff
    was chosen. This is what separates Rewind's outcome from GSE279162's."""
    f = rec["GSE227151"]["findings"]["outcome_threshold_sensitivity"]
    assert f["value"]["is_invariant"] is True
    lo, hi = f["value"]["invariant_over"]
    assert lo == 1 and hi >= 500
    for floor, row in f["value"]["sweep"].items():
        if int(floor) <= hi:
            assert row["primed_clones"] == 82 and row["primed_cells"] == 102, floor


@has_results
def test_the_gdna_arm_is_bimodal_and_no_positive_sits_in_the_noise_mode(rec):
    d = rec["GSE227151"]["findings"]["gdna_read_distribution"]["value"]
    assert d["1"] > 1000, "a large 1-read noise mode exists"
    assert d[">=500"] > 300, "and a separate colony mode"
    lo, hi = rec["GSE227151"]["findings"]["outcome_threshold_sensitivity"]["value"]["invariant_over"]
    assert hi > 2, "the weakest positive clone must sit above the 1-2 read noise mode"


@has_results
def test_the_published_42_is_not_reproduced_and_was_not_chased(rec):
    """Case 3. Recording a failed check is the point; adopting the one reading that happens to hit
    42 would be fitting the label to the answer."""
    f = rec["GSE227151"]["findings"]["published_42_check"]
    r = f["value"]["readings"]
    assert r["pooled_cells"] == 102 and r["pooled_clones"] == 82
    assert 42 not in (r["pooled_cells"], r["sample1_cells"], r["sample2_cells"]), \
        "no CELL-level reading reproduces the published 42"
    assert f["value"]["readings_equal_to_42"] == ["sample1_clones"]
    assert "not adopted" in f["evidence"]


@has_results
def test_the_outer_split_unit_must_be_the_clone(rec):
    c = rec["GSE227151"]["findings"]["clone_structure"]["value"]
    assert c["shared"] > 0, "clones really do span both lanes"
    assert c["clones_total"] == 3149
    assert c["clones_sample1"] + c["clones_sample2"] - c["shared"] == c["clones_total"]
    assert "leak" in rec["GSE227151"]["findings"]["clone_structure"]["evidence"]


@has_results
def test_rewind_label_counts_are_rare_event(rec):
    v = rec["GSE227151"]["findings"]["prospective_label_counts"]["value"]
    assert v["positive_clones"] == 82 and v["positive_cells"] == 102
    assert v["negative_clones"] == 3067 and v["negative_cells"] == 3819
    assert v["positive_clones"] + v["negative_clones"] == 3149
    assert v["positive_rate_clones"] < 0.05


@has_results
def test_only_one_biological_replicate_is_locally_reconstructable(rec):
    f = rec["GSE227151"]["findings"]["replicate_structure"]
    assert f["value"]["biological_replicates_local"] == 1
    assert set(f["value"]["gsms_local"]) == {"GSM7092515", "GSM7092516"}
    assert f["value"]["gsms_in_series_not_local"], "the rest of the series is named, not hidden"


# ---- Role B: GSE279162 -------------------------------------------------------------------------#
@has_results
def test_the_clone_id_is_a_feature_row_verified_from_the_matrices(rec):
    f = rec["GSE279162"]["findings"]["feature_structure"]["value"]
    assert f["n_genes"] == 36601 and f["n_lineage_features"] == 153055
    assert f["total"] == 189656 and f["genes_first"] is True
    assert f["by_type"]["Custom"] == 153055


@has_results
def test_per_cell_clone_assignment_is_not_clean_enough_to_call_without_a_rule(rec):
    f = rec["GSE279162"]["findings"]["per_cell_barcode_dominance"]
    assert f["value"]["median_dominant_fraction"] <= 0.6
    assert f["value"]["median_lineage_features_per_cell"] > 1
    assert "REQUIRES an explicit" in f["evidence"]


@has_results
def test_wm989_counts_really_are_floor_conditional(rec):
    """Case 4. The justification for PENDING has to be measured, not asserted."""
    sweep = rec["GSE279162"]["findings"]["outcome_threshold_sensitivity"]["value"]
    naive = [sweep[k]["naive_clones"] for k in sorted(sweep, key=int)]
    assert naive == sorted(naive, reverse=True), "the naive pool must shrink as the floor rises"
    assert naive[0] > 2 * naive[-1], "and it must move materially, not marginally"
    for t in s21d.WM989_TREATMENTS:
        assert sweep["1"][t]["positives"] != sweep["50"][t]["positives"], t


@has_results
def test_the_multi_treatment_structure_role_b_was_qualified_for_actually_exists(rec):
    v = rec["GSE279162"]["findings"]["multi_treatment_coverage"]["value"]
    assert v["clones_in_ge2_treatments"] >= 1000
    assert v["clones_in_all_6"] > 0
    assert sum(v["by_n_treatments"].values()) == v["naive_clones"]


@has_results
def test_the_six_treatments_and_three_naive_lanes_come_from_the_source_design(rec):
    f = rec["GSE279162"]["findings"]["design_from_source"]
    assert f["value"]["treatments"] == list(s21d.WM989_TREATMENTS)
    assert f["value"]["pretreatment_samples"] == list(s21d.WM989_NAIVE)
    assert "Overall-Design" in f["evidence"]


# ---- verdicts ----------------------------------------------------------------------------------#
@has_results
def test_both_verdicts_are_pending_and_each_names_what_is_unresolved(rec):
    for gse in ("GSE227151", "GSE279162"):
        assert rec[gse]["verdict"] == s21d.RECON_PENDING, gse
        f = rec[gse]["findings"]["unresolved_outcome_rule"]
        assert f["status"] == UNK, "an unresolved rule is UNKNOWN, never ABSENT"
        assert len(f["value"]) >= 2, "the open questions are enumerated, not summarised"


@has_results
def test_the_pending_verdict_is_driven_by_the_unknown_not_by_opinion(rec):
    """A PENDING verdict must be structurally tied to an UNKNOWN finding, so it cannot be asserted
    while every finding is resolved, nor withheld while one is not."""
    for gse in ("GSE227151", "GSE279162"):
        unknowns = [k for k, f in rec[gse]["findings"].items() if f["status"] == UNK]
        assert unknowns, f"{gse} is PENDING, so something must be UNKNOWN"


@has_results
def test_stage_22_does_not_open_on_a_pending_reconstruction(rec):
    assert rec["overall"] == s21d.STAGE22_PENDING
    assert s21d.STAGE22_OPEN != s21d.STAGE22_PENDING != s21d.STAGE22_BLOCKED


# ---- additive contract and engineering constraints ----------------------------------------------#
@has_results
def test_stage_21d_did_not_touch_the_earlier_frozen_results(rec):
    assert STAGE_21A.exists() and STAGE_21B.exists()
    a = json.loads(STAGE_21A.read_text(encoding="utf-8"))
    assert a["verdict"] == "CULTURE_FORWARD_AVAILABLE", "21A's frozen verdict is unchanged"
    b = json.loads(STAGE_21B.read_text(encoding="utf-8"))
    assert b["GSE242423"]["verdict"] == "LINEAGE_ABSENT_PROVEN"
    assert b["GSE165176"]["verdict"] == "ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY"
    assert rec["additive_to"] == ["results/diag_stage21_data_audit_results.json",
                                  "results/diag_stage21b_source_design_results.json"]


@has_results
def test_no_model_was_fitted(rec):
    """Checked against the parsed code, not the prose: the module docstring is allowed to SAY
    'no sklearn' without that counting as an import of it."""
    import ast

    assert rec["model_fitted"] is False and rec["src_modified"] is False
    src = SRC.read_text(encoding="utf-8")
    tree = ast.parse(src)
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not (imported & {"sklearn", "torch", "tensorflow", "xgboost", "lightgbm", "statsmodels",
                            "cellfate"}), f"21D must not model or import src/: {sorted(imported)}"

    body = src[src.index('"""', src.index('"""') + 3) + 3:]   # everything after the docstring
    for banned in ("LogisticRegression", "sklearn", "torch", ".fit("):
        assert banned not in body, f"21D must not model: found {banned} in code"


@has_results
def test_the_raw_public_data_stayed_out_of_git(rec):
    """Provenance is committed; the matrices are not."""
    assert rec["raw_data_committed"] is False
    for gse in ("GSE227151", "GSE279162"):
        entries = rec[gse]["manifest"]
        assert entries, gse
        for e in entries:
            assert e["exists"] and e["sha256"] and e["bytes"]
            assert not Path(e["path"]).is_relative_to(ROOT), f"{e['path']} is inside the repo"


def test_the_script_writes_only_its_own_results_and_the_two_small_tables():
    src = SRC.read_text(encoding="utf-8")
    assert src.count(".write_text(") == 3
    for name in ("diag_stage21d_public_reconstruction_results.json",
                 "stage21d_rewind_clone_table.tsv",
                 "stage21d_gse279162_clone_table.tsv"):
        assert name in src


@has_results
def test_the_committed_clone_tables_are_small_and_present(rec):
    for gse, cap in (("GSE227151", 400_000), ("GSE279162", 400_000)):
        p = ROOT / rec[gse]["clone_table"]
        assert p.exists(), p
        assert p.stat().st_size < cap, f"{p} is too large to belong in git"
