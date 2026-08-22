"""Stage 23A — the protocol freeze, and the contracts that keep it frozen.

23A is the substage where a learnability gate can be quietly pre-broken: if `X` is built wrong, or
the protocol leaves a choice open, every later verdict inherits the flaw and no downstream test
will notice. These tests pin the four ways that happens here:

1. **Normalising twice.** Clone pseudobulk must sum RAW counts and apply CP10K + log1p exactly
   once. Summing normalised cells, or a second log1p, silently rescales every distance -- so the
   test inverts the transform and checks the row sums really are 10,000.
2. **A lineage feature reaching X.** WM989 ships 153,055 `Custom` LinNNNN features beside the
   36,601 genes, and they encode clone identity. `X` must be exactly 36,601 columns wide.
3. **A treated cell reaching X_before.** Only the two Rewind control lanes and the three WM989
   naive lanes may be summed.
4. **An open choice.** Anything a later substage could shop for -- K, grids, seeds, the reference
   treatment, the tie-break -- has to be in `stage23_protocol.json` now, not decided later.

Plus the audit contract itself: Stage-22's `overall` string is never trusted on its own.
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
SRC = ROOT / "experiments" / "run_stage23_learnability_gate.py"
spec = importlib.util.spec_from_file_location("s23", SRC)
S23 = importlib.util.module_from_spec(spec)
sys.modules["s23"] = S23
spec.loader.exec_module(S23)

RES = ROOT / "results"
PROTOCOL = RES / "stage23_protocol.json"
REWIND_MAN = RES / "stage23_rewind_clone_expression_manifest.json"
WM989_MAN = RES / "stage23_wm989_clone_expression_manifest.json"
PREP = RES / "stage23_outer_fold_preprocessing.json"
CACHE = ROOT / "_cc_cache" / "stage23"

frozen = pytest.mark.skipif(not PROTOCOL.exists(), reason="23A has not been run")
has_cache = pytest.mark.skipif(
    not (CACHE / "GSE227151_pseudobulk.npz").exists(),
    reason="pseudobulk cache is gitignored and absent (this is the CI condition)")


@pytest.fixture(scope="module")
def proto():
    return json.loads(PROTOCOL.read_text(encoding="utf-8"))


# ---- no data required --------------------------------------------------------------------------#
def test_the_canonical_text_hash_is_lf_crlf_invariant(tmp_path):
    """V2 §1.4. Stage 22 recorded checkout-byte hashes, so its plan digest differs between
    Windows and Linux for identical content. This rule is what stops Stage 23 inheriting that."""
    lf = tmp_path / "lf.txt"
    crlf = tmp_path / "crlf.txt"
    lf.write_bytes(b"alpha\nbeta\ngamma\n")
    crlf.write_bytes(b"alpha\r\nbeta\r\ngamma\r\n")
    assert S23.canonical_text_sha256(lf) == S23.canonical_text_sha256(crlf)
    assert S23.sha256_file(lf) != S23.sha256_file(crlf), "the raw byte hashes really do differ"
    cr = tmp_path / "cr.txt"
    cr.write_bytes(b"alpha\rbeta\rgamma\r")
    assert S23.canonical_text_sha256(cr) == S23.canonical_text_sha256(lf)


def test_the_frozen_constants_match_the_plan():
    assert S23.K_CANDIDATES == (10, 20, 50)
    assert S23.SEED_PROTOCOL == 23023
    assert (S23.SEED_BOOT_REWIND, S23.SEED_BOOT_WM989_C1, S23.SEED_BOOT_WM989_C2) == (
        23123, 23223, 23224)
    assert S23.SEED_PERMUTATION == 23323
    assert S23.N_BOOTSTRAP == 2000 and S23.N_PERMUTATION == 200
    assert S23.LOGISTIC_C == (0.01, 0.1, 1, 10) and S23.RIDGE_ALPHA == (0.1, 1, 10, 100)
    assert S23.N_GENES == 36601


def test_the_treatment_coding_is_case_sensitive_and_matches_stage22():
    """The V1 blocker: V1 froze `acid` as the reference, but the benchmark contains `Acid`."""
    wt = pd.read_csv(RES / "stage22_wm989_clone_treatment.csv")
    assert list(S23.TREATMENT_ORDER) == sorted(wt.treatment.unique())
    assert S23.REFERENCE_TREATMENT == "Acid"
    assert S23.REFERENCE_TREATMENT in set(wt.treatment), "reference must exist verbatim"
    assert len(S23.TREATMENT_ORDER) - 1 == 5, "five non-reference dummies"


def test_the_wm989_nuisance_block_carries_total_captured_depth():
    """V2 §1.2.1. The frozen confound is keyed on TOTAL depth, and log1p is not additive, so the
    three per-lane terms cannot represent it. Omitting it would bias the test toward X."""
    assert "log1p(n_naive_cells)" in S23.WM989_NUISANCE
    assert len(S23.WM989_NUISANCE) == 4
    wk = pd.read_csv(RES / "stage22_wm989_clones.csv")
    for term in S23.WM989_NUISANCE:
        assert term[len("log1p("):-1] in wk.columns


def test_stage23_may_tighten_the_firewall_but_never_loosen_it():
    res = json.loads((RES / "stage22_prospective_benchmark_results.json").read_text("utf-8"))
    fe = res["feature_eligibility"]
    for col, move in S23.STAGE23_TIGHTENED.items():
        stage22 = next(k for k, v in fe.items() if col in v)
        assert stage22 == move["stage22"], f"{col} moved in Stage 22 without notice"
        assert move["stage23"] == "TARGET", "a tightening only ever goes toward TARGET"
    forbidden = set(fe["TARGET"]) | set(fe["PROVENANCE_ONLY"])
    assert not (forbidden & set(S23.REWIND_NUISANCE) & set(S23.WM989_NUISANCE))


# ---- the frozen protocol -----------------------------------------------------------------------#
@frozen
def test_the_protocol_froze_everything_a_later_substage_could_shop_for(proto):
    assert proto["verdict"] == S23.PROTOCOL_FROZEN
    assert proto["model_fitted"] is False
    for key in ("seeds", "outer_cv", "inner_cv", "representation", "grids", "tie_break",
                "treatment_coding", "nuisance_blocks", "inner_selection_metric", "inference",
                "feature_firewall"):
        assert key in proto, key
    assert proto["representation"]["k_candidates"] == [10, 20, 50]
    assert proto["representation"]["pca"]["svd_solver"] == "randomized"
    assert proto["representation"]["pca"]["fit_once_at_max_K_then_reuse_prefixes"] is True
    assert proto["treatment_coding"]["reference"] == "Acid"
    assert proto["treatment_coding"]["standardized"] is False
    assert proto["inference"]["permutation_p"].startswith("(1 + ")


@frozen
def test_the_input_audit_ran_and_did_not_trust_the_overall_string(proto):
    a = proto["input_audit"]
    assert a["verdict"] == S23.AUDITED and not a["failed_checks"]
    for k in ("role_a_ready", "all_gates_pass_true", "every_individual_gate_true",
              "model_fitted_false", "preflight_derived_independently_of_overall_string"):
        assert a["checks"][k]["ok"], k
    assert a["checks"]["stage22_overall_string_for_reference_only"]["detail"] == "STAGE_23_READY"
    assert "without consuming" in a["checks"]["known_stage22_gate_derivation_limitation"]["detail"]
    assert a["checks"]["role_b_status_preserved_separately"]["detail"].startswith("BENCHMARK_READY")


@frozen
def test_the_audit_verified_the_stage22_hash_chain_and_the_raw_files(proto):
    a = proto["input_audit"]["checks"]
    for k in ("six_benchmark_csv_hashes_match_manifests", "two_manifest_hashes_match_results",
              "results_file_omits_its_own_hash", "GSE227151_raw_files_present_and_identical",
              "GSE279162_raw_files_present_and_identical",
              "every_expression_column_index_resolves"):
        assert a[k]["ok"], k
    assert a["every_expression_column_index_resolves"]["detail"] == 0


@frozen
def test_the_wm989_target_dependency_was_independently_rebuilt(proto):
    """V2 §1.3: the literal `True` evidence field from Stage 22 is not accepted as proof."""
    d = proto["input_audit"]["checks"]["wm989_targets_rebuild_from_treated_cells_only"]
    assert d["ok"] and d["detail"]["rows"] == 8406
    assert d["detail"]["max_abs_diff"] == 0
    assert d["detail"]["naive_cells_used"] == 0
    note = proto["input_audit"]["checks"]["wm989_inherited_joint_assignment_dependency"]["detail"]
    assert "0.50%" in note and "is NOT rebuilt" in note


@frozen
def test_the_inherited_crlf_limitation_is_declared_not_hidden(proto):
    d = proto["input_audit"]["checks"]["inherited_crlf_lf_provenance_limitation"]["detail"]
    for who in ("stage22_plan", "stage22_builder"):
        assert d[who]["recorded"] == d[who]["checkout_bytes"], "byte hash still verifies"
        assert d[who]["canonical_lf"] != d[who]["checkout_bytes"], "and LF really differs"
    assert "NOT rewritten" in d["note"]


@frozen
def test_the_frozen_counts_and_folds_were_reverified(proto):
    a = proto["input_audit"]["checks"]
    for k in ("rewind_counts", "wm989_counts", "rewind_positive_clones_per_fold",
              "outer_folds_are_five_each", "rewind_cell_fold_matches_clone_fold",
              "wm989_clone_treatment_fold_matches_clone_fold", "wm989_c2_eligible_clones_929"):
        assert a[k]["ok"], k
    assert a["wm989_c2_eligible_clones_929"]["detail"] == 929


@frozen
def test_the_feature_universe_is_gene_expression_only_and_id_keyed(proto):
    a = proto["input_audit"]["checks"]
    assert a["every_sample_has_36601_gene_expression_features"]["ok"]
    assert a["wm989_custom_block_is_153055_where_present"]["ok"]
    assert a["gene_block_is_the_first_36601_rows"]["ok"]
    assert a["all_samples_share_one_gene_feature_id_list"]["ok"]
    assert a["all_samples_share_one_gene_feature_id_list"]["detail"]["distinct_signatures"] == 1
    for m in (json.loads(REWIND_MAN.read_text("utf-8")), json.loads(WM989_MAN.read_text("utf-8"))):
        assert m["genes"] == 36601
        assert "stable 10x feature id" in m["feature_key"]
        assert "never symbol or row index" in m["feature_key"]


@frozen
def test_one_pseudobulk_row_per_clone_and_only_pretreatment_cells(proto):
    r = json.loads(REWIND_MAN.read_text("utf-8"))
    w = json.loads(WM989_MAN.read_text("utf-8"))
    assert r["clones"] == 3147 and r["pretreatment_cells_summed"] == 3905
    assert w["clones"] == 1401 and w["pretreatment_cells_summed"] == 6489
    # only the two Rewind control lanes and the three WM989 naive lanes were read
    assert all("control" in s for s in r["source_samples"])
    assert all("Naive" in s for s in w["source_samples"])
    assert len(w["source_samples"]) == 3
    for t in ("Dabrafenib", "Trametinib", "CoCl2", "Acid", "Cisplatin", "Doxorubicin"):
        assert not any(t in s for s in w["source_samples"]), f"{t} reached X_before"


@frozen
def test_the_normalization_rule_is_recorded_as_applied_once(proto):
    for m in (json.loads(REWIND_MAN.read_text("utf-8")), json.loads(WM989_MAN.read_text("utf-8"))):
        assert m["normalization"] == (
            "sum raw counts per clone -> CP10K -> log1p, applied exactly once")
        assert m["committed_to_git"] is False
        assert len(m["matrix_content_sha256"]) == 64


@frozen
def test_every_candidate_k_is_feasible_in_every_outer_fold():
    prep = json.loads(PREP.read_text(encoding="utf-8"))["per_dataset"]
    for ds, folds in prep.items():
        assert len(folds) == 5, ds
        for f, v in folds.items():
            assert v["all_K_feasible"], f"{ds} fold {f} cannot support K=50"
            assert v["max_feasible_K"] >= 50
            assert v["retained_genes"] > 0
            assert v["detection_floor"] == max(5, int(np.ceil(0.01 * v["outer_training_clones"])))


@frozen
def test_the_outer_fold_gene_filter_is_only_descriptive():
    note = json.loads(PREP.read_text(encoding="utf-8"))["note"]
    assert "descriptive only" in note and "inner split" in note


# ---- checks that need the gitignored cache ----------------------------------------------------#
@has_cache
def test_the_pseudobulk_really_is_cp10k_log1p_applied_once():
    """Inverting the transform is the only way to catch a double normalisation: a second log1p
    still yields a plausible-looking non-negative matrix."""
    from scipy import sparse
    for name, n_clones in (("GSE227151", 3147), ("GSE279162", 1401)):
        X = sparse.load_npz(CACHE / f"{name}_pseudobulk.npz")
        assert X.shape == (n_clones, 36601)
        back = X.copy()
        back.data = np.expm1(back.data)
        rows = np.asarray(back.sum(axis=1)).ravel()
        assert np.allclose(rows, 1e4, rtol=1e-6), f"{name} row sums are not CP10K"
        assert (X.data >= 0).all() and not np.isnan(X.data).any()


@has_cache
def test_the_cached_matrix_matches_its_committed_content_hash():
    from scipy import sparse
    for name, man in (("GSE227151", REWIND_MAN), ("GSE279162", WM989_MAN)):
        X = sparse.load_npz(CACHE / f"{name}_pseudobulk.npz")
        got = hashlib.sha256(
            X.indptr.tobytes() + X.indices.tobytes() + np.round(X.data, 10).tobytes()).hexdigest()
        assert got == json.loads(man.read_text("utf-8"))["matrix_content_sha256"]


@has_cache
def test_clone_order_matches_the_frozen_benchmark_exactly():
    for name, table in (("GSE227151", "stage22_rewind_clones.csv"),
                        ("GSE279162", "stage22_wm989_clones.csv")):
        clones = json.loads((CACHE / f"{name}_clones.json").read_text(encoding="utf-8"))
        assert clones == sorted(clones), "clone order must be deterministic"
        assert set(clones) == set(pd.read_csv(RES / table)["clone_id"])


# ---- engineering constraints -------------------------------------------------------------------#
def test_23a_fits_no_model_and_reuses_no_production_estimator():
    """V2 §2.6: the existing CellFateNet / fixed baseline registry are not Stage-23 estimators."""
    import ast

    src = SRC.read_text(encoding="utf-8")
    imported = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            imported |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not {m for m in imported if m.split(".")[0] in {"torch", "cellfate"}}
    # 23B onward fits estimators by design; the 23A code path must not.
    tree = ast.parse(src)
    a_functions = {"audit_stage22_inputs", "clone_pseudobulk", "training_fold_gene_filter",
                   "expression_manifest", "build_protocol"}
    for fn in (n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name in a_functions):
        calls = [c for c in ast.walk(fn)
                 if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute)
                 and c.func.attr in {"fit", "fit_transform", "fit_predict"}]
        assert not calls, f"{fn.name} is part of 23A and must fit nothing"
    for banned in ("CellFateNet", "_LinearBase"):
        assert banned not in src


def test_23a_does_not_write_a_large_matrix_into_results():
    for p in RES.glob("stage23_*"):
        assert p.suffix in {".json", ".csv"}, p.name
        assert p.stat().st_size < 2_000_000, f"{p.name} is too large to commit"
    assert not list(RES.glob("*pseudobulk*")), "the matrix must stay in the gitignored cache"


def test_the_pseudobulk_cache_is_gitignored():
    ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    assert "_cc_cache/" in ignore
    assert str(CACHE.relative_to(ROOT)).startswith("_cc_cache")


@frozen
def test_no_stage22_artifact_was_rewritten_by_23a(proto):
    res = json.loads((RES / "stage22_prospective_benchmark_results.json").read_text("utf-8"))
    assert res["overall"] == "STAGE_23_READY" and res["all_gates_pass"] is True
    for a in res["manifest_hashes"]:
        assert S23.sha256_file(RES / a["name"]) == a["sha256"], f"{a['name']} was modified"


# ================================================================================================ #
# 23B — Rewind Role-A learnability
# ================================================================================================ #
REWIND_RESULTS = RES / "stage23_rewind_results.json"
REWIND_OOF = RES / "stage23_rewind_oof_predictions.csv"
ran_23b = pytest.mark.skipif(not REWIND_RESULTS.exists(), reason="23B has not been run")


@pytest.fixture(scope="module")
def rb():
    return json.loads(REWIND_RESULTS.read_text(encoding="utf-8"))


@ran_23b
def test_exactly_the_four_pre_registered_models_are_present(rb):
    assert set(rb["models"]) == {"R0", "R1", "R2", "R3"}
    assert set(rb["pooled_oof_metrics"]) == {"R0", "R1", "R2", "R3"}
    assert rb["models"]["R0"] == "outer-training prevalence"
    assert rb["models"]["R1"] == "nuisance only"


@ran_23b
def test_average_precision_is_primary_and_no_accuracy_gate_exists(rb):
    assert rb["primary_metric"] == "average_precision_score at clone grain"
    blob = json.dumps(rb).lower()
    assert "accuracy" not in blob, "accuracy is reporting-only and must not appear in a gate"
    for m in ("R0", "R1", "R2", "R3"):
        assert set(rb["pooled_oof_metrics"][m]) == {"AP", "ROC_AUC", "log_loss", "brier"}


@ran_23b
def test_one_oof_prediction_per_clone_and_35_positives_exactly_once(rb):
    oof = pd.read_csv(REWIND_OOF)
    clones_tbl = pd.read_csv(RES / "stage22_rewind_clones.csv")
    assert len(oof) == 3147 == rb["clones"]
    assert oof["clone_id"].is_unique
    assert set(oof["clone_id"]) == set(clones_tbl["clone_id"])
    assert int(oof["y_primed"].sum()) == 35 == rb["positives"]
    for m in ("R0", "R1", "R2", "R3"):
        assert oof[f"pred_{m}"].notna().all(), m
        assert ((oof[f"pred_{m}"] >= 0) & (oof[f"pred_{m}"] <= 1)).all(), m


@ran_23b
def test_the_outer_folds_are_exactly_the_frozen_stage22_assignment(rb):
    oof = pd.read_csv(REWIND_OOF)
    frozen_fold = pd.read_csv(RES / "stage22_rewind_clones.csv").set_index("clone_id")["outer_fold"]
    assert (oof.set_index("clone_id")["outer_fold"] == frozen_fold.loc[oof["clone_id"]].values).all()
    assert sorted(oof["outer_fold"].unique()) == [0, 1, 2, 3, 4]
    assert oof.groupby("outer_fold")["y_primed"].sum().tolist() == [7, 7, 7, 7, 7]


@ran_23b
def test_r0_is_a_constant_per_fold_equal_to_the_training_prevalence(rb):
    """R0 must be the outer-TRAINING prevalence, never the test prevalence."""
    oof = pd.read_csv(REWIND_OOF)
    for f, g in oof.groupby("outer_fold"):
        assert g["pred_R0"].nunique() == 1, f
        expected = rb["per_outer_fold"][str(f)]["train_prevalence"]
        assert abs(g["pred_R0"].iloc[0] - expected) < 1e-6, "recorded to 6dp"
        # and it is the TRAINING prevalence, not this fold's test prevalence
        assert abs(g["pred_R0"].iloc[0] - g["y_primed"].mean()) > 0


@ran_23b
def test_hyperparameters_come_only_from_the_frozen_grids(rb):
    for f, sel in rb["selected_hyperparameters_per_outer_fold"].items():
        assert sel["R1"]["K"] is None, "the nuisance-only model has no K"
        for m in ("R1", "R2", "R3"):
            assert sel[m]["C"] in S23.LOGISTIC_C, (f, m)
        for m in ("R2", "R3"):
            assert sel[m]["K"] in S23.K_CANDIDATES, (f, m)


@ran_23b
def test_fold_direction_is_reported_but_is_not_a_pass_requirement(rb):
    """V2 §4.5/§4.6: with seven positive clones per fold this is a high-variance diagnostic.
    V1 made it an AND-condition; the audit showed that was too noisy to gate on."""
    assert rb["fold_direction_is_diagnostic_only"] is True
    assert len(rb["per_fold_average_precision"]) == 5
    for f, a in rb["per_fold_average_precision"].items():
        assert "delta_AP_state" in a and "delta_AP_absolute" in a
        assert abs(a["delta_AP_state"] - (a["R3"] - a["R1"])) < 1e-12, f


@ran_23b
def test_the_bootstrap_is_the_pre_registered_stratified_clone_bootstrap(rb):
    for key in ("delta_AP_state_R3_minus_R1", "delta_AP_absolute_R2_minus_R0"):
        b = rb["inference"][key]
        assert b["replicates"] == 2000 and b["seed"] == 23123
        assert b["ci95_low"] <= b["ci95_high"]
        assert 0.0 <= b["fraction_delta_le_0"] <= 1.0
    s = rb["inference"]["delta_AP_state_R3_minus_R1"]
    assert abs(s["point"] - (rb["pooled_oof_metrics"]["R3"]["AP"]
                             - rb["pooled_oof_metrics"]["R1"]["AP"])) < 1e-12


@ran_23b
def test_the_provisional_verdict_is_derived_not_asserted(rb):
    s = rb["inference"]["delta_AP_state_R3_minus_R1"]
    expected = (S23.ROLE_A_FAIL if s["point"] <= 0
                else S23.ROLE_A_PASS if s["ci95_low"] > 0 else S23.ROLE_A_WEAK)
    assert rb["provisional_verdict"] == expected
    src = SRC.read_text(encoding="utf-8")
    assert '"provisional_verdict": "ROLE_A' not in src, "the verdict must not be hard-coded"


@ran_23b
def test_the_verdict_stays_provisional_until_23e(rb):
    assert "23E" in rb["verdict_is_provisional_until"]
    assert "PERMUTATION" in rb["verdict_is_provisional_until"]


@ran_23b
def test_convergence_warnings_are_surfaced_not_swallowed(rb):
    """V2 §3.7: a convergence warning is a protocol failure to investigate, never a silent drop."""
    assert isinstance(rb["convergence_warnings"], list)
    assert rb["convergence_warnings"] == [], rb["convergence_warnings"]
    src = SRC.read_text(encoding="utf-8")
    assert "ConvergenceWarning" in src and "convergence.append" in src


@ran_23b
def test_no_target_or_provenance_column_reached_the_design_matrix():
    """The design matrix is built from exactly two blocks: PC scores and the nuisance terms.

    Checked on the syntax tree: `_design` is the only place a design matrix is assembled, and it
    concatenates nothing except `pcs` and `nuis`. A target or provenance column could only enter
    by being stuffed into one of those two, and the nuisance block is pinned below.
    """
    import ast

    assert S23.REWIND_NUISANCE == ("log1p(n_pretreatment_cells)", "n_lanes")
    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "_design")
    names = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
    assert names <= {"parts", "pcs", "nuis", "k", "use_x", "use_nuis", "np"}, names
    build = next(n for n in ast.walk(tree)
                 if isinstance(n, ast.FunctionDef) and n.name == "run_23b")
    cols = [ast.unparse(n) for n in ast.walk(build)
            if isinstance(n, ast.Call) and ast.unparse(n).startswith("np.column_stack")]
    assert len(cols) == 1 and "n_pretreatment_cells" in cols[0] and "n_lanes" in cols[0]
    assert "y_primed" not in cols[0] and "outer_fold" not in cols[0]


@ran_23b
def test_the_results_reference_the_protocol_that_was_frozen_in_23a(rb):
    assert rb["protocol_sha256"] == S23.sha256_file(RES / "stage23_protocol.json")
    assert rb["plan"]["version"] == "V2"


# ================================================================================================ #
# 23C — WM989 additive state-signal gate
# ================================================================================================ #
WM989_RESULTS = RES / "stage23_wm989_results.json"
C1_OOF = RES / "stage23_wm989_detection_oof.csv"
C2_OOF = RES / "stage23_wm989_abundance_oof.csv"
ran_23c = pytest.mark.skipif(not WM989_RESULTS.exists(), reason="23C has not been run")


@pytest.fixture(scope="module")
def wb():
    return json.loads(WM989_RESULTS.read_text(encoding="utf-8"))


@ran_23c
def test_exactly_w0_to_w4_are_present_and_no_interaction_yet(wb):
    """23C is the ADDITIVE gate. X x U belongs to 23D and must not leak forward."""
    assert set(wb["models"]) == {"W0", "W1", "W2", "W3", "W4"}
    assert wb["models"] == {"W0": "U", "W1": "B + U", "W2": "X", "W3": "X + U",
                            "W4": "X + B + U"}
    assert wb["interaction_terms_present"] is False
    assert "W5" not in json.dumps(wb)
    assert wb["primary_comparison"] == "W4 vs W1 on both endpoints"


@ran_23c
def test_the_nuisance_block_is_the_four_frozen_terms(wb):
    assert wb["nuisance_block_B"] == list(S23.WM989_NUISANCE)
    assert "log1p(n_naive_cells)" in wb["nuisance_block_B"], "total captured depth is mandatory"


@ran_23c
def test_the_treatment_coding_is_the_frozen_case_sensitive_one(wb):
    assert wb["treatment_coding"]["order"] == list(S23.TREATMENT_ORDER)
    assert wb["treatment_coding"]["reference"] == "Acid"
    d = S23.treatment_dummies(np.array(["Acid", "CoCl2", "Trametinib"]))
    assert d.shape == (3, 5), "five non-reference dummies"
    assert d[0].sum() == 0, "the reference level is all-zero"
    assert d[1].sum() == 1 and d[2].sum() == 1


@ran_23c
def test_c1_covers_every_clone_treatment_row_and_c2_only_the_nonzero_ones(wb):
    c1, c2 = pd.read_csv(C1_OOF), pd.read_csv(C2_OOF)
    assert len(c1) == 8406 == wb["endpoints"]["C1"]["rows"]
    assert c1["clone_id"].nunique() == 1401
    assert (c1.groupby("clone_id").size() == 6).all(), "all six treatments per clone"
    assert len(c2) == 2256 == wb["endpoints"]["C2"]["rows"]
    assert c2["clone_id"].nunique() == 929 == wb["endpoints"]["C2"]["clones"]
    ct = pd.read_csv(RES / "stage22_wm989_clone_treatment.csv")
    nz = ct[ct.n_post_cells > 0]
    assert set(map(tuple, c2[["clone_id", "treatment"]].to_numpy())) == set(
        map(tuple, nz[["clone_id", "treatment"]].to_numpy()))


@ran_23c
def test_c1_target_is_detection_and_is_never_renamed_resistance(wb):
    c1 = pd.read_csv(C1_OOF)
    ct = pd.read_csv(RES / "stage22_wm989_clone_treatment.csv").sort_values(
        ["clone_id", "treatment"]).reset_index(drop=True)
    assert (c1["y"].to_numpy() == (ct["n_post_cells"].to_numpy() > 0).astype(int)).all()
    blob = json.dumps(wb).lower()
    for word in ("resistant", "resistance", "survival", "death"):
        assert word not in blob, f"C1 detection must not be renamed {word!r}"


@ran_23c
def test_c2_target_is_log1p_of_the_nonzero_count(wb):
    c2 = pd.read_csv(C2_OOF)
    ct = pd.read_csv(RES / "stage22_wm989_clone_treatment.csv")
    nz = ct[ct.n_post_cells > 0].sort_values(["clone_id", "treatment"]).reset_index(drop=True)
    assert np.allclose(c2.sort_values(["clone_id", "treatment"])["y"].to_numpy(),
                       np.log1p(nz["n_post_cells"].to_numpy()))


@ran_23c
def test_the_c2_metric_is_clone_balanced_not_row_averaged(wb):
    """A clone seen under five treatments must not count as five independent units."""
    c2 = pd.read_csv(C2_OOF)
    err = np.abs(c2["y"] - c2["pred_W4"])
    row_mean = float(err.mean())
    clone_balanced = float(pd.DataFrame({"c": c2["clone_id"], "e": err})
                           .groupby("c")["e"].mean().mean())
    reported = wb["endpoints"]["C2"]["pooled_oof_metrics"]["W4"]["clone_balanced_MAE"]
    assert abs(reported - clone_balanced) < 1e-9
    assert abs(reported - row_mean) > 1e-9, "the two must actually differ, else the test is vacuous"


@ran_23c
def test_wm989_hyperparameters_come_only_from_the_frozen_grids(wb):
    for ep, grid in (("C1", S23.LOGISTIC_C), ("C2", S23.RIDGE_ALPHA)):
        sel = wb["endpoints"][ep]["selected_hyperparameters_per_outer_fold"]
        assert len(sel) == 5
        for f, models in sel.items():
            for m, v in models.items():
                assert v["hp"] in grid, (ep, f, m, v)
                if m in ("W0", "W1"):
                    assert v["K"] is None, f"{m} uses no expression block"
                else:
                    assert v["K"] in S23.K_CANDIDATES, (ep, f, m)


@ran_23c
def test_the_outer_folds_are_the_frozen_stage22_ones(wb):
    frozen = pd.read_csv(RES / "stage22_wm989_clones.csv").set_index("clone_id")["outer_fold"]
    for path in (C1_OOF, C2_OOF):
        df = pd.read_csv(path)
        assert (df["clone_id"].map(frozen).to_numpy() == df["outer_fold"].to_numpy()).all()
        assert df.groupby("clone_id")["outer_fold"].nunique().max() == 1


@ran_23c
def test_the_bootstrap_resamples_clones_and_carries_their_rows(wb):
    for ep, seed, n in (("C1", 23223, 1401), ("C2", 23224, 929)):
        d = wb["endpoints"][ep]["inference"]["delta_state_W1_minus_W4"]
        assert d["seed"] == seed and d["replicates"] == 2000
        assert d["bootstrap_unit"] == "clone" and d["clones_resampled"] == n
        assert len(d["ci95"]) == 2 and len(d["ci975_two_sided"]) == 2
        assert d["ci975_two_sided"][0] <= d["ci95"][0], "97.5% interval must be the wider one"
        assert d["ci975_two_sided"][1] >= d["ci95"][1]


@ran_23c
def test_the_reported_delta_equals_w1_minus_w4(wb):
    for ep in ("C1", "C2"):
        e = wb["endpoints"][ep]
        d = e["inference"]["delta_state_W1_minus_W4"]
        key = d["metric"]
        expect = e["pooled_oof_metrics"]["W1"][key] - e["pooled_oof_metrics"]["W4"][key]
        assert abs(d["point"] - expect) < 1e-12


@ran_23c
def test_the_verdict_is_derived_from_the_pre_registered_rule(wb):
    ll = wb["endpoints"]["C1"]["inference"]["delta_state_W1_minus_W4"]
    ma = wb["endpoints"]["C2"]["inference"]["delta_state_W1_minus_W4"]
    pass_ll, pass_ma = ll["ci975_two_sided"][0] > 0, ma["ci975_two_sided"][0] > 0
    harm_ll, harm_ma = ll["ci975_two_sided"][1] < 0, ma["ci975_two_sided"][1] < 0
    if (pass_ll and not harm_ma) or (pass_ma and not harm_ll):
        expected = S23.ROLE_B_PASS
    elif ll["point"] <= 0 and ma["point"] <= 0:
        expected = S23.ROLE_B_FAIL
    elif harm_ll or harm_ma:
        expected = S23.ROLE_B_FAIL
    else:
        expected = S23.ROLE_B_WEAK
    assert wb["verdict"] == expected
    src = SRC.read_text(encoding="utf-8")
    assert '"verdict": "ROLE_B' not in src, "the verdict must not be hard-coded"


@ran_23c
def test_the_abundance_baseline_is_load_bearing_and_visible(wb):
    """The whole point of V2 §1.2.1: W1 must be a serious competitor, not a formality."""
    c1 = wb["endpoints"]["C1"]["pooled_oof_metrics"]
    assert c1["W1"]["log_loss"] < c1["W0"]["log_loss"], "B must improve on treatment-only"
    b_gain = c1["W0"]["log_loss"] - c1["W1"]["log_loss"]
    x_gain = c1["W1"]["log_loss"] - c1["W4"]["log_loss"]
    assert b_gain > x_gain, "if X ever exceeds B here, re-read the confound before believing it"


@ran_23c
def test_no_convergence_warning_was_swallowed(wb):
    assert wb["convergence_warnings"] == [], wb["convergence_warnings"]


@ran_23c
def test_the_results_reference_the_23a_protocol(wb):
    assert wb["protocol_sha256"] == S23.sha256_file(RES / "stage23_protocol.json")
    assert wb["plan"]["version"] == "V2"
    assert "23E" in wb["verdict_is_provisional_until"]


# ================================================================================================ #
# 23D — WM989 explicit interaction gate
# ================================================================================================ #
INTERACTION_RESULTS = RES / "stage23_wm989_interaction_results.json"
ran_23d = pytest.mark.skipif(not INTERACTION_RESULTS.exists(), reason="23D has not been run")


@pytest.fixture(scope="module")
def wi():
    return json.loads(INTERACTION_RESULTS.read_text(encoding="utf-8"))


def test_the_interaction_block_is_pc_by_dummy_only_and_zero_for_the_reference():
    """V2 §6.1: only standardized PC score x non-reference dummy. The reference treatment's rows
    must have an all-zero interaction block, so its state contribution stays in the common X
    coefficients rather than being double counted."""
    pcs = np.arange(12.0).reshape(3, 4)
    d = S23.treatment_dummies(np.array(["Acid", "CoCl2", "Trametinib"]))
    blk = S23.interaction_block(pcs, d)
    assert blk.shape == (3, 4 * 5), "exactly 5K columns"
    assert (blk[0] == 0).all(), "the reference row carries no interaction"
    assert blk[1].sum() != 0 and blk[2].sum() != 0
    for row in (1, 2):
        assert int((blk[row] != 0).sum()) == 4, "one treatment block active per row"


@ran_23d
def test_w1_and_w4_are_reused_verbatim_from_23c(wi):
    """The reference must not drift: 23D refits only W5."""
    assert wi["w1_w4_reused_from_23c"] is True
    for name, path in (("C1", RES / "stage23_wm989_interaction_oof.csv"),
                       ("C2", RES / "stage23_wm989_interaction_abundance_oof.csv")):
        new = pd.read_csv(path)
        old = pd.read_csv(RES / ("stage23_wm989_detection_oof.csv" if name == "C1"
                                 else "stage23_wm989_abundance_oof.csv"))
        assert np.allclose(new["pred_W1"], old["pred_W1"]), name
        assert np.allclose(new["pred_W4"], old["pred_W4"]), name
        assert (new["clone_id"] == old["clone_id"]).all()
        assert "pred_W5" in new.columns and new["pred_W5"].notna().all()


@ran_23d
def test_the_design_width_is_exactly_k_plus_b_plus_u_plus_interaction(wi):
    for ep in ("C1", "C2"):
        e = wi["endpoints"][ep]
        for f, sel in e["selected_hyperparameters_per_outer_fold"].items():
            k = sel["K"]
            assert sel["interaction_columns"] == 5 * k
            assert e["per_outer_fold"][f]["design_columns"] == k + 4 + 5 + 5 * k, (ep, f)


@ran_23d
def test_no_gene_level_interaction_was_constructed(wi):
    assert "no gene-level interaction" in wi["interaction_terms"]
    for ep in ("C1", "C2"):
        for f, m in wi["endpoints"][ep]["per_outer_fold"].items():
            assert m["design_columns"] < 400, (ep, f, "a gene-level interaction would be enormous")


@ran_23d
def test_both_required_comparisons_are_present_with_both_intervals(wi):
    for ep in ("C1", "C2"):
        inf = wi["endpoints"][ep]["inference"]
        assert set(inf) == {"interaction_W4_minus_W5", "full_state_W1_minus_W5"}
        for key, d in inf.items():
            assert d["replicates"] == 2000 and d["bootstrap_unit"] == "clone"
            assert d["ci975_two_sided"][0] <= d["ci95"][0]
            assert d["ci975_two_sided"][1] >= d["ci95"][1]
            a, b = d["comparison"].split(" - ")
            expect = wi["endpoints"][ep]["pooled"][a] - wi["endpoints"][ep]["pooled"][b]
            assert abs(d["point"] - expect) < 1e-12, key


@ran_23d
def test_the_full_state_comparison_is_required_not_optional(wi):
    """V2 §6.3: W5 must beat the load-bearing nuisance baseline W1, not merely rearrange error
    relative to W4."""
    for ep in ("C1", "C2"):
        assert "full_state_W1_minus_W5" in wi["endpoints"][ep]["inference"]
    assert "pass_full" in wi["endpoint_families"]["C1"]


@ran_23d
def test_treatment_level_directions_cover_all_six(wi):
    for ep in ("C1", "C2"):
        e = wi["endpoints"][ep]
        assert set(e["by_treatment"]) == set(S23.TREATMENT_ORDER)
        n = sum(v["improved"] for v in e["by_treatment"].values())
        assert n == e["treatments_improved_by_W5_over_W4"]
        for t, v in e["by_treatment"].items():
            assert abs(v["improvement_W4_minus_W5"] - (v["W4"] - v["W5"])) < 1e-12, t
            assert v["improved"] == (v["W4"] - v["W5"] > 0)


@ran_23d
def test_the_verdict_is_derived_from_the_four_pre_registered_criteria(wi):
    fams = wi["endpoint_families"]
    got = S23.INTERACTION_NONE
    passing = None
    for ep, other in (("C1", "C2"), ("C2", "C1")):
        a, b = fams[ep], fams[other]
        if (a["pass_int"] and a["pass_full"] and a["n_treat"] >= 3
                and not b["harm_int"] and not b["harm_full"]):
            got, passing = S23.INTERACTION_PASS, ep
            break
    if got != S23.INTERACTION_PASS:
        if (any(f["point_int"] > 0 or f["point_full"] > 0 for f in fams.values())
                and any(f["n_treat"] >= 1 for f in fams.values())):
            got = S23.INTERACTION_LOCAL
    assert wi["verdict"] == got
    assert wi["passing_endpoint"] == passing
    src = SRC.read_text(encoding="utf-8")
    assert '"verdict": "INTERACTION' not in src, "the verdict must not be hard-coded"


@ran_23d
def test_a_pass_required_beating_w1_on_the_same_endpoint(wi):
    """The specific failure mode V2 §6.5 exists to block: W5 'passing' by rearranging error
    against W4 while still not beating the nuisance baseline."""
    if wi["verdict"] != S23.INTERACTION_PASS:
        pytest.skip("no PASS to check")
    ep = wi["passing_endpoint"]
    f = wi["endpoint_families"][ep]
    assert f["pass_int"] and f["pass_full"], "both bounds are required on the passing endpoint"
    assert f["n_treat"] >= 3, "one favourable treatment is not a broad interaction claim"


@ran_23d
def test_the_outer_folds_are_still_the_frozen_stage22_ones(wi):
    frozen = pd.read_csv(RES / "stage22_wm989_clones.csv").set_index("clone_id")["outer_fold"]
    for path in ("stage23_wm989_interaction_oof.csv",
                 "stage23_wm989_interaction_abundance_oof.csv"):
        df = pd.read_csv(RES / path)
        assert (df["clone_id"].map(frozen).to_numpy() == df["outer_fold"].to_numpy()).all()


@ran_23d
def test_no_convergence_warning_and_the_protocol_is_referenced(wi):
    assert wi["convergence_warnings"] == [], wi["convergence_warnings"]
    assert wi["protocol_sha256"] == S23.sha256_file(RES / "stage23_protocol.json")
    assert "23E" in wi["verdict_is_provisional_until"]


# ============================================================================================== #
# 23E — permutation nulls, structural controls, provenance sentinel, determinism.
#
# 23E is the substage most easily faked, because it is the one that decides whether the three
# earlier PASS verdicts survive. The ways it can be faked are all mechanical:
#
#   * a permutation that leaks -- a profile crossing the outer train/test boundary, or leaving its
#     stratum, quietly rebuilds the very structure the null is supposed to destroy;
#   * caching that is convenient rather than exact -- reusing an outer-training basis is only valid
#     because the permutation preserves the outer-training profile SET, and that has to be checked
#     numerically, not asserted in a comment;
#   * a p-value without the +1 correction, or a "pass" that needs only one of the two gates;
#   * a failure quietly rounded into a pass.
#
# The last one is why `test_role_a_is_recorded_as_a_permutation_failure` exists.
# ============================================================================================== #
PERM_RESULTS = RES / "stage23_permutation_results.json"
DET_RESULTS = RES / "stage23_determinism.json"
ran_23e = pytest.mark.skipif(not PERM_RESULTS.exists(), reason="23E has not been run")


@pytest.fixture(scope="module")
def pe():
    return json.loads(PERM_RESULTS.read_text(encoding="utf-8"))


def test_permutation_p_uses_the_plus_one_correction():
    """V2 §7.3. Without the +1 a null that never reaches the observed value reports p = 0, which
    claims more resolution than 200 draws can carry."""
    null = np.zeros(200)
    got = S23.permutation_p(1.0, null)
    assert got["n_null_ge_observed"] == 0
    assert got["p_perm"] == pytest.approx(1 / 201)
    assert got["p_perm"] > 0, "a finite null can never license p = 0"
    got2 = S23.permutation_p(0.0, null)
    assert got2["n_null_ge_observed"] == 200
    assert got2["p_perm"] == pytest.approx(1.0)


def test_a_pass_needs_both_the_p95_gate_and_the_p_value():
    """Either gate alone is satisfiable by a statistic that is not actually extreme."""
    heavy_tail = np.concatenate([np.zeros(190), np.full(10, 5.0)])
    p95_only = S23.permutation_p(1.0, heavy_tail)
    assert p95_only["exceeds_null_p95"] is True
    assert p95_only["p_perm"] > 0.05
    assert p95_only["passes"] is False, "the p95 gate alone must not promote a claim"

    flat = np.linspace(0.0, 1.0, 200)
    assert S23.permutation_p(0.5, flat)["passes"] is False
    assert S23.permutation_p(2.0, flat)["passes"] is True


def test_the_permutation_never_crosses_the_outer_boundary_or_leaves_its_stratum():
    """V2 §7.1. A profile that moves across the fold boundary, or between strata, reintroduces
    exactly the depth/lane structure the null exists to hold fixed."""
    rng = np.random.default_rng(7)
    strata = np.array(["a", "a", "a", "b", "b", "b", "b", "c"] * 4)
    side = np.array([True, True, False, True, False, True, False, True] * 4)
    for _ in range(50):
        out = S23.permute_within(strata, side, rng)
        assert (side[out] == side).all(), "a profile changed outer train/test side"
        assert (strata[out] == strata).all(), "a profile left its stratum"
        assert sorted(out.tolist()) == list(range(len(strata))), "not a bijection"


def test_the_permutation_actually_moves_something():
    """A 'null' that is the identity map would make every p-value 1.0 by construction."""
    rng = np.random.default_rng(3)
    strata = np.array(["a"] * 40)
    side = np.array([True] * 30 + [False] * 10)
    moved = sum(int((S23.permute_within(strata, side, rng) != np.arange(40)).any())
                for _ in range(20))
    assert moved >= 19, "the permutation is very nearly the identity"


def test_rewind_strata_are_the_frozen_bins():
    """V2 §7.1: n_pretreatment_cells in {1, 2, 3+} crossed with n_lanes."""
    tbl = pd.DataFrame({"n_pretreatment_cells": [1, 2, 3, 9, 1], "n_lanes": [1, 1, 2, 2, 2]})
    assert list(S23.rewind_strata(tbl)) == ["1|1", "2|1", "3+|2", "3+|2", "1|2"]


@frozen
def test_wm989_strata_never_leave_a_cell_under_four_clones():
    """The frozen merge rule exists so a stratum of one cannot shuffle only with itself."""
    wk = pd.read_csv(RES / "stage22_wm989_clones.csv")
    strat = S23.wm989_strata(wk)
    counts = pd.Series(strat).value_counts()
    assert (counts >= 4).all(), counts[counts < 4].to_dict()
    assert len(strat) == len(wk)


def test_the_cached_outer_training_basis_equals_a_fresh_fit():
    """The one claim that makes 23E affordable: the final outer-training transform depends on the
    outer-training profile SET, which the permutation preserves, so it may be cached.

    Checked numerically against `expression_block` rather than argued in a comment.
    """
    from scipy import sparse

    rng = np.random.default_rng(11)
    dense = rng.random((60, 200)) * (rng.random((60, 200)) < 0.4)
    X = sparse.csr_matrix(dense)
    tr = np.arange(0, 45)
    te = np.arange(45, 60)
    ztr, zte, n_keep, _ = S23.expression_block(X, tr, te, 10)
    cache = S23._frozen_pipeline_cache(X, tr, 10)
    assert len(cache["keep"]) == n_keep
    assert np.allclose(S23._apply_cached(X, tr, cache), ztr)
    assert np.allclose(S23._apply_cached(X, te, cache), zte)

    shuffled = tr[rng.permutation(len(tr))]
    reordered = S23._frozen_pipeline_cache(X, shuffled, 10)
    assert list(reordered["keep"]) == list(cache["keep"]), \
        "the basis must depend on the training SET, not on its order"


def test_the_expression_free_models_are_reused_rather_than_refitted():
    """R1/W1 carry no expression, so a permutation of X cannot move them. The null code takes them
    from the frozen observed results instead of refitting -- refitting would burn hours and, worse,
    would let a null run silently disagree with the observed run it is compared against."""
    src = SRC.read_text(encoding="utf-8")
    fam = src.split("def run_23e_family(")[1].split("\ndef ")[0]
    assert 'rew["pooled_oof_metrics"]["R1"]["AP"]' in fam
    assert 'wmc["endpoints"]["C1"]["pooled_oof_metrics"]["W1"]' in fam
    body = src.split("def _wm989_null_once(")[1].split("\ndef ")[0]
    assert '"W4"' in body and '"W5"' in body
    for absent in ('"W0"', '"W2"', '"W3"'):
        assert absent not in body, f"the null refits {absent}, which no permutation can move"


@ran_23e
def test_the_frozen_permutation_count_and_seed_were_used(pe):
    assert pe["n_permutations"] == S23.N_PERMUTATION == 200
    assert pe["permutation_base_seed"] == S23.SEED_PERMUTATION == 23323
    for k, v in pe["permutation_tests"].items():
        assert v["n_permutations"] == 200, f"{k} was tested on a shortened null"


@ran_23e
def test_every_reported_p_value_matches_its_own_null(pe):
    """Recompute p from the cached null arrays rather than trusting the recorded number."""
    nulls: dict = {}
    for fam in ("rewind", "wm989c1", "wm989c2"):
        pth = CACHE / f"stage23e_null_{fam}.json"
        if pth.exists():
            nulls.update(json.loads(pth.read_text(encoding="utf-8"))["nulls"])
    if not nulls:
        pytest.skip("null cache is absent")
    for key, rec in pe["permutation_tests"].items():
        arr = np.array(nulls[key])
        assert len(arr) == 200
        again = S23.permutation_p(rec["observed"], arr)
        assert again["p_perm"] == pytest.approx(rec["p_perm"])
        assert again["passes"] is rec["passes"]
        assert again["null_p95"] == pytest.approx(rec["null_p95"])


@ran_23e
def test_role_a_is_recorded_as_a_permutation_failure(pe):
    """The honest-reporting contract. Role A cleared its bootstrap CI in 23B and then failed here;
    a later edit that quietly flips this to PASS without new evidence must break a test."""
    ra = pe["permutation_tests"]["role_a_delta_AP_state"]
    assert ra["passes"] is False
    assert ra["p_perm"] > 0.05
    assert ra["null_mean"] > 0, ("the null mean is positive -- selection optimism alone produces a "
                                 "gain, which is the whole finding")
    assert pe["claim_permutation_status"]["ROLE_A_PERMUTATION_PASS"] is False


@ran_23e
def test_the_non_candidate_statistic_is_declared_rather_than_silently_dropped(pe):
    """Additive C2 failed its 23C bootstrap, so it is not a PASS candidate. Not testing it is
    legitimate; not saying so would not be."""
    assert "c2_delta_MAE_state" in pe["not_permutation_tested"]
    assert "c2_delta_MAE_state" not in pe["permutation_tests"]


@ran_23e
def test_all_five_structural_controls_ran_and_passed(pe):
    required = {"outer_test_isolation", "feature_firewall", "frozen_fold_identity",
                "canonical_text_hash_lf_crlf", "fresh_clone_determinism"}
    assert required <= set(pe["structural_controls"])
    for name in required:
        assert pe["structural_controls"][name]["ok"] is True, name
    assert pe["STRUCTURAL_CONTROLS_PASS"] is True
    assert pe["STRUCTURAL_CONTROLS_PASS"] == all(
        v["ok"] for v in pe["structural_controls"].values())


@ran_23e
def test_determinism_compared_the_full_artifact_set_against_a_clean_tree():
    if not DET_RESULTS.exists():
        pytest.skip("the determinism check has not been run")
    d = json.loads(DET_RESULTS.read_text(encoding="utf-8"))
    assert d["artifacts_compared"] == len(S23.DETERMINISM_ARTIFACTS)
    assert len(S23.DETERMINISM_ARTIFACTS) >= 13, "the compared set must not shrink"
    assert d["working_tree_clean_for_builder_and_artifacts"] is True, (
        "a dirty builder makes the provenance hashes unreproducible by construction")
    assert d["mismatched"] == {}
    assert d["all_match"] is True
    assert set(d["committed_digests"]) == set(S23.DETERMINISM_ARTIFACTS)


def test_the_determinism_set_covers_every_committed_stage23_artifact():
    """A shrinking artifact list would make determinism trivially true."""
    skip = {"stage23_permutation_results.json", "stage23_determinism.json"}
    on_disk = {p.name for p in RES.glob("stage23_*") if p.name not in skip}
    assert on_disk <= set(S23.DETERMINISM_ARTIFACTS), on_disk - set(S23.DETERMINISM_ARTIFACTS)


def test_the_sentinel_sees_presence_flags_only():
    """V2 §7.4. The sentinel is only meaningful if it is strictly weaker than the real models: no
    expression, no captured counts, no clone identity."""
    body = SRC.read_text(encoding="utf-8").split("def provenance_sentinel(")[1].split("\ndef ")[0]
    assert "> 0).astype(float)" in body, "the sentinel must binarise, not read counts"
    assert "get_dummies" not in body and "OneHot" not in body, \
        "clone identity must never be encoded"
    for banned in ("_load_rewind_x", "_load_wm989_x", "expression_block", "PCA("):
        assert banned not in body, f"the sentinel reached expression via {banned}"


@ran_23e
def test_the_sentinel_does_not_reach_the_models_whose_gain_is_claimed(pe):
    """The alert condition V2 §7.4 actually specifies: if library presence alone reaches R3/W4,
    the claimed gain is library structure rather than biology."""
    s = pe["provenance_sentinel"]
    assert s["rewind"]["sentinel_AP"] < s["rewind"]["R3_AP"]
    assert s["rewind"]["reaches_R3_without_expression"] is False
    assert s["wm989_c1"]["sentinel_log_loss"] > s["wm989_c1"]["W4_log_loss"]
    assert s["wm989_c1"]["reaches_W4_without_expression"] is False
    assert s["rewind"]["alert"] is False and s["wm989_c1"]["alert"] is False


@ran_23e
def test_the_promotion_status_is_mechanical(pe):
    """Every claim's status must follow from its own permutation result -- never from a judgement
    call written into the record."""
    t = pe["permutation_tests"]
    st = pe["claim_permutation_status"]
    assert st["ROLE_A_PERMUTATION_PASS"] == t["role_a_delta_AP_state"]["passes"]
    assert st["ROLE_B_ADDITIVE_PERMUTATION_PASS"] == t["c1_delta_LL_state"]["passes"]
    assert st["ROLE_B_INTERACTION_PERMUTATION_PASS"] == (
        t["c1_delta_LL_interaction"]["passes"] and t["c1_delta_LL_full"]["passes"])
    assert st["C2_INTERACTION_SECONDARY_PERMUTATION_PASS"] == (
        t["c2_delta_MAE_interaction"]["passes"] and t["c2_delta_MAE_full"]["passes"])


@ran_23e
def test_23e_references_the_frozen_protocol_and_plan(pe):
    assert pe["protocol_sha256"] == S23.sha256_file(RES / "stage23_protocol.json")
    assert pe["plan"]["canonical_lf_sha256"] == S23.canonical_text_sha256(S23.PLAN)
    assert pe["stage"] == "23E"


# ============================================================================================== #
# 23F — mechanical synthesis.
#
# 23F is where a stage that produced one failure and three passes gets summarised, and summarising
# is where the failure would go missing. The contracts here are therefore about *arithmetic and
# provenance*, not about performance:
#
#   * 23F must fit nothing -- a synthesis step that trains anything is no longer a synthesis;
#   * every verdict must be recomputable from the frozen artifacts by a test that never reads
#     23F's own conclusions;
#   * the roadmap gate must depend on Role A alone, so no amount of Role-B strength can open
#     Stage 24;
#   * all seven findings the stage must carry forward have to be present AND agree with their
#     source artifact, so a later edit cannot quietly soften one.
# ============================================================================================== #
SYNTHESIS = RES / "stage23_final_synthesis.json"
ran_23f = pytest.mark.skipif(not SYNTHESIS.exists(), reason="23F has not been run")


@pytest.fixture(scope="module")
def sy():
    return json.loads(SYNTHESIS.read_text(encoding="utf-8"))


def test_23f_fits_nothing():
    """A synthesis step that trains anything is not a synthesis. Checked on the AST of the 23F
    code path rather than on its own self-reported `models_fitted_in_23f` field."""
    import ast

    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    f_functions = {"run_23f", "_lower_bound"}
    for fn in (n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name in f_functions):
        for call in (c for c in ast.walk(fn) if isinstance(c, ast.Call)):
            if isinstance(call.func, ast.Attribute):
                assert call.func.attr not in {"fit", "fit_transform", "fit_predict", "predict",
                                              "predict_proba"}, f"{fn.name} fits or predicts"
            if isinstance(call.func, ast.Name):
                assert not call.func.id.startswith("_fit_"), f"{fn.name} calls {call.func.id}"
                assert call.func.id not in {"_load_rewind_x", "_load_wm989_x", "expression_block",
                                            "clone_pseudobulk"}, f"{fn.name} touched expression"


@ran_23f
def test_23f_reports_that_it_fitted_nothing(sy):
    assert sy["models_fitted_in_23f"] == 0
    assert sy["synthesis_is_mechanical"] is True
    assert sy["stage"] == "23F"


@ran_23f
def test_every_verdict_is_recomputable_from_the_frozen_artifacts_alone(sy):
    """Re-derive all four verdicts from 23B/23C/23D/23E without reading 23F's conclusions."""
    rb = json.loads((RES / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    wc = json.loads((RES / "stage23_wm989_results.json").read_text(encoding="utf-8"))
    wi = json.loads((RES / "stage23_wm989_interaction_results.json").read_text(encoding="utf-8"))
    pe = json.loads((RES / "stage23_permutation_results.json").read_text(encoding="utf-8"))
    ctl = pe["STRUCTURAL_CONTROLS_PASS"]
    st = pe["claim_permutation_status"]

    a_ok = (rb["inference"]["delta_AP_state_R3_minus_R1"]["ci95_low"] > 0
            and st["ROLE_A_PERMUTATION_PASS"] and ctl)
    assert sy["final_verdicts"]["role_a"] == (S23.ROLE_A_PASS if a_ok else S23.ROLE_A_FAIL)

    b = wc["endpoints"]["C1"]["inference"]["delta_state_W1_minus_W4"]
    b_ok = b["ci975_two_sided"][0] > 0 and st["ROLE_B_ADDITIVE_PERMUTATION_PASS"] and ctl
    assert sy["final_verdicts"]["role_b_additive"] == (S23.ROLE_B_PASS if b_ok
                                                       else S23.ROLE_B_FAIL)

    e = wi["endpoints"]["C1"]["inference"]
    i_ok = (e["interaction_W4_minus_W5"]["ci975_two_sided"][0] > 0
            and e["full_state_W1_minus_W5"]["ci975_two_sided"][0] > 0
            and st["ROLE_B_INTERACTION_PERMUTATION_PASS"] and ctl)
    broad = wi["endpoints"]["C1"]["treatments_improved_by_W5_over_W4"] >= 3
    assert sy["final_verdicts"]["role_b_interaction"] == (
        S23.INTERACTION_PASS if (i_ok and broad)
        else S23.INTERACTION_LOCAL if i_ok else S23.INTERACTION_NONE)

    e2 = wi["endpoints"]["C2"]["inference"]
    c2_ok = (e2["interaction_W4_minus_W5"]["ci975_two_sided"][0] > 0
             and e2["full_state_W1_minus_W5"]["ci975_two_sided"][0] > 0
             and st["C2_INTERACTION_SECONDARY_PERMUTATION_PASS"] and ctl)
    assert sy["final_verdicts"]["c2_interaction_secondary"] == (
        S23.C2_SECONDARY_CONFIRMED if c2_ok else S23.C2_SECONDARY_NOT_CONFIRMED)


@ran_23f
def test_role_a_failed_and_that_is_what_the_ledger_says(sy):
    """The demotion has to be visible in the ledger, not only in prose. Role A's bootstrap
    criterion PASSED and its permutation gate FAILED -- both facts must survive."""
    a = sy["claims"]["role_a"]
    assert a["required_by_plan"] is True
    assert a["bootstrap_criterion"]["excludes_zero"] is True, "23B's PASS candidacy is not erased"
    assert a["permutation_gate"]["passes"] is False
    assert a["permutation_gate"]["p_perm"] > 0.05
    assert a["final"] == S23.ROLE_A_FAIL


@ran_23f
def test_23f_does_not_agree_with_23b_and_says_so(sy):
    """23B recorded a provisional PASS. 23F must contradict it, or the permutation gate did
    nothing."""
    rb = json.loads((RES / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    assert rb["provisional_verdict"] == S23.ROLE_A_PASS
    assert sy["final_verdicts"]["role_a"] == S23.ROLE_A_FAIL
    assert sy["final_verdicts"]["role_a"] != rb["provisional_verdict"]


@ran_23f
def test_23f_agrees_with_23c_and_23d_where_nothing_overturned_them(sy):
    """The converse contract: a synthesis that silently changed a passing substage's verdict would
    be just as wrong as one that hid a failure."""
    wc = json.loads((RES / "stage23_wm989_results.json").read_text(encoding="utf-8"))
    wi = json.loads((RES / "stage23_wm989_interaction_results.json").read_text(encoding="utf-8"))
    assert sy["final_verdicts"]["role_b_additive"] == wc["verdict"]
    assert sy["final_verdicts"]["role_b_interaction"] == wi["verdict"]


@ran_23f
def test_the_gate_depends_on_role_a_alone(sy):
    """The load-bearing rule: Role B is strong on three separate statistics and still cannot open
    Stage 24."""
    g = sy["roadmap_gate"]
    assert g["role_a_is_mandatory"] is True
    assert g["role_b_may_substitute_for_role_a"] is False
    assert g["gate"] == (S23.STAGE_24_OPEN
                         if sy["final_verdicts"]["role_a"] == S23.ROLE_A_PASS
                         else S23.STAGE_24_BLOCKED_ROLE_A)
    assert g["gate"] == S23.STAGE_24_BLOCKED_ROLE_A
    # Role B really is positive -- the gate is blocked despite that, not because of it.
    assert sy["final_verdicts"]["role_b_additive"] == S23.ROLE_B_PASS
    assert sy["final_verdicts"]["role_b_interaction"] == S23.INTERACTION_PASS
    assert sy["final_verdicts"]["c2_interaction_secondary"] == S23.C2_SECONDARY_CONFIRMED


@ran_23f
def test_the_gate_routes_to_the_failure_resolution_stage(sy):
    """Roadmap V4 sends STAGE_24_BLOCKED_ROLE_A to Stage 23R, not to Stage 24 and not to a rerun
    of Role A inside Stage 23."""
    nxt = sy["roadmap_gate"]["next_stage"]
    assert "23R" in nxt
    assert "STAGE 24" not in nxt.upper().replace("STAGE 23R", "")


@ran_23f
def test_all_seven_preserved_findings_are_present(sy):
    p = sy["preserved_findings"]
    assert len(p) == 7, sorted(p)
    for i, fragment in enumerate(
            ["rewind_absolute_signal", "role_a_permutation_p", "abundance_remains",
             "state_adds_beyond_abundance", "explicit_interaction_adds",
             "doxorubicin", "no_external_generalization"], start=1):
        key = next((k for k in p if k.startswith(f"{i}_")), None)
        assert key is not None, f"finding {i} is missing"
        assert fragment in key, f"finding {i} is not about {fragment}"
        assert p[key], f"finding {i} is empty"


@ran_23f
def test_the_preserved_findings_agree_with_their_source_artifacts(sy):
    """Each carried-forward number is re-read from the artifact it claims to come from."""
    rb = json.loads((RES / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    wc = json.loads((RES / "stage23_wm989_results.json").read_text(encoding="utf-8"))
    wi = json.loads((RES / "stage23_wm989_interaction_results.json").read_text(encoding="utf-8"))
    pe = json.loads((RES / "stage23_permutation_results.json").read_text(encoding="utf-8"))
    p = sy["preserved_findings"]

    f1 = p["1_rewind_absolute_signal_was_weak_before_any_permutation"]
    assert f1["positives"] == rb["positives"] == 35
    assert f1["R3_AP"] == rb["pooled_oof_metrics"]["R3"]["AP"]
    assert f1["prevalence"] == pytest.approx(rb["positives"] / rb["clones"])

    f2 = p["2_role_a_permutation_p"]
    assert f2["p_perm"] == pe["permutation_tests"]["role_a_delta_AP_state"]["p_perm"]
    assert f2["p_perm"] == pytest.approx(0.0846, abs=1e-4)

    f3 = p["3_abundance_remains_the_dominant_wm989_predictor"]
    w0 = wc["endpoints"]["C1"]["pooled_oof_metrics"]["W0"]["log_loss"]
    w1 = wc["endpoints"]["C1"]["pooled_oof_metrics"]["W1"]["log_loss"]
    assert f3["abundance_gain_W0_minus_W1_log_loss"] == pytest.approx(w0 - w1)
    assert f3["ratio"] > 1, "abundance must still be recorded as the dominant predictor"

    f4 = p["4_state_adds_beyond_abundance_on_c1"]
    assert f4["delta_log_loss_state_W1_minus_W4"] == (
        wc["endpoints"]["C1"]["inference"]["delta_state_W1_minus_W4"]["point"])

    f5 = p["5_explicit_interaction_adds_further_signal"]
    assert f5["delta_log_loss_interaction_W4_minus_W5"] == (
        wi["endpoints"]["C1"]["inference"]["interaction_W4_minus_W5"]["point"])

    f6 = p["6_doxorubicin_is_the_consistent_treatment_level_exception"]
    for ep, key in (("C1", "C1_improvement_W4_minus_W5"), ("C2", "C2_improvement_W4_minus_W5")):
        assert f6[key] == wi["endpoints"][ep]["by_treatment"]["Doxorubicin"][
            "improvement_W4_minus_W5"]
        assert f6[key] < 0, f"Doxorubicin must stay recorded as negative on {ep}"
    assert f6["negative_on_both_endpoints"] is True

    f7 = p["7_no_external_generalization_has_been_shown"]
    assert f7["external_biological_replicate_tested"] is False
    assert f7["unseen_treatment_tested"] is False
    assert f7["cross_dataset_transfer_tested"] is False


@ran_23f
def test_23f_pins_the_artifacts_it_synthesised(sy):
    """A synthesis is only reproducible if it names the exact inputs it read."""
    for stage, name in (("23B", "stage23_rewind_results.json"),
                        ("23C", "stage23_wm989_results.json"),
                        ("23D", "stage23_wm989_interaction_results.json"),
                        ("23E", "stage23_permutation_results.json")):
        assert sy["source_artifacts"][stage] == S23.sha256_file(RES / name), stage
    assert sy["protocol_sha256"] == S23.sha256_file(RES / "stage23_protocol.json")
    assert sy["plan"]["canonical_lf_sha256"] == S23.canonical_text_sha256(S23.PLAN)


# ============================================================================================== #
# Formal closure.
#
# The closure declaration is prose, so it is the one part of Stage 23 that could drift away from
# the artifacts without anything breaking. These contracts tie it back down: the verdicts and
# digests written into the record must be the ones the frozen files actually carry, and the
# stage's sections must all be present and in order.
# ============================================================================================== #
RECORD = ROOT / "plans" / "(newer)practical plans" / "RECORDs" / "stage_23_RECORD.md"
closed = pytest.mark.skipif(not RECORD.exists(), reason="the Stage-23 record is absent")


@closed
def test_the_record_carries_every_substage_and_the_closure_in_order():
    import re

    text = RECORD.read_text(encoding="utf-8")
    heads = re.findall(r"^# (23[A-F]|STAGE 23 — FORMAL CLOSURE)", text, re.M)
    assert heads == ["23A", "23B", "23C", "23D", "23E", "23F", "STAGE 23 — FORMAL CLOSURE"], heads


@closed
@ran_23f
def test_the_closure_verdicts_are_the_ones_the_artifacts_carry(sy):
    """A closure section that disagrees with `stage23_final_synthesis.json` would be the single
    most misleading thing in the repository."""
    block = RECORD.read_text(encoding="utf-8").split("# STAGE 23 — FORMAL CLOSURE")[1]
    for name, verdict in sy["final_verdicts"].items():
        assert verdict in block, f"the closure does not state the {name} verdict {verdict}"
    assert sy["roadmap_gate"]["gate"] in block
    assert "STAGE_24_OPEN" not in block, "Stage 24 is blocked; the closure must not say otherwise"
    assert "23R" in block, "the closure must name the stage that follows"


@closed
@ran_23f
def test_the_digests_quoted_in_the_closure_match_the_files(sy):
    block = RECORD.read_text(encoding="utf-8").split("# STAGE 23 — FORMAL CLOSURE")[1]
    assert sy["protocol_sha256"] in block
    for stage, name in (("23B", "stage23_rewind_results.json"),
                        ("23C", "stage23_wm989_results.json"),
                        ("23D", "stage23_wm989_interaction_results.json"),
                        ("23E", "stage23_permutation_results.json")):
        digest = S23.sha256_file(RES / name)
        assert digest == sy["source_artifacts"][stage]
        assert digest in block, f"{stage}'s digest in the closure is not the file's digest"
    assert S23.sha256_file(SYNTHESIS) in block, "23F's own digest is stale in the closure"


@closed
def test_the_closure_states_that_role_a_stays_failed():
    """The rule Stage 23R is bound by. If this sentence is ever softened, a later stage could
    quietly relabel a failed pre-registered test as a pass."""
    block = RECORD.read_text(encoding="utf-8").split("# STAGE 23 — FORMAL CLOSURE")[1]
    assert "permanent" in block.lower()
    assert "may **not** relabel" in block or "may not relabel" in block
    assert "confirmatory evidence" in block
