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
