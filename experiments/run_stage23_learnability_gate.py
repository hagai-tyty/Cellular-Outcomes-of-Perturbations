"""STAGE 23 — learnability / interaction gate.

Pre-registered in `plans/(newer)practical plans/STAGE_23_LEARNABILITY_INTERACTION_GATE_V2.md`
(V1 archived under `arcive/` after the independent pre-execution audit). Consumes the frozen
Stage-22 benchmark (`8d6011a`).

**23A fits NO outcome model.** It runs the §3.0 input audit, builds the clone-level `X_before`
representation, and freezes the evaluation protocol. Nothing there touches an outcome except to
verify that Stage-22's targets are reconstructible from post-treatment data alone.

**23B is the first substage that fits an estimator.** It may only use what 23A froze: the outer
folds come from Stage 22 untouched, and every hyperparameter comes from the grids in
`stage23_protocol.json`. Each learned quantity -- gene filter, gene scaler, PCA, PC scaler,
nuisance scaler -- is refitted inside every inner-training split, then rebuilt once on the whole
outer-training set before the outer test fold is touched exactly once.

WHY THE ORDER MATTERS
---------------------
The audit runs FIRST and gates everything after it. If a Stage-22 artifact, fold, hash or raw file
has moved, pseudobulk construction never starts -- because a representation built on a shifted
benchmark would be silently wrong rather than loudly broken.

THE TRAPS THIS FILE EXISTS TO AVOID
-----------------------------------
1. **Normalising twice.** Clone pseudobulk must sum RAW counts across the clone's cells and apply
   CP10K + log1p exactly once. Summing already-normalised cells, or log1p-ing a second time, would
   quietly change every downstream distance.
2. **Aligning by row position.** Feature identity is the stable 10x feature ID. Row position is a
   coincidence that happens to hold here (audited: all 11 source samples share one identical
   feature list AND order) and must be asserted, never assumed.
3. **Letting a lineage feature into X.** WM989 ships 153,055 `Custom` LinNNNN features in the same
   matrix as the 36,601 genes. They encode clone identity: including one is direct provenance
   leakage. Only rows <= n_genes enter.
4. **Letting a treated cell into X.** `X_before` is pretreatment only -- Rewind's two control lanes
   and WM989's three naive lanes. No treated column is ever read.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_CACHE = ROOT / "_cc_cache" / "stage23"          # gitignored: the matrices never enter git

_spec_path = Path(__file__).resolve().parent / "build_stage22_prospective_benchmarks.py"
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("s22", _spec_path)
S22 = importlib.util.module_from_spec(_spec)
sys.modules["s22"] = S22
_spec.loader.exec_module(S22)
S21D = S22.S21D

PLAN = ROOT / "plans/(newer)practical plans/STAGE_23_LEARNABILITY_INTERACTION_GATE_V2.md"
PLAN_VERSION = "V2"
STAGE22_BENCHMARK_COMMIT = "8d6011a"

# ---- frozen protocol constants, transcribed from V2 ------------------------------------------ #
SEED_PROTOCOL = 23023          # PCA + inner CV
SEED_BOOT_REWIND = 23123
SEED_BOOT_WM989_C1 = 23223
SEED_BOOT_WM989_C2 = 23224
SEED_PERMUTATION = 23323
K_CANDIDATES = (10, 20, 50)
N_INNER = 3
N_OUTER = 5
LOGISTIC_C = (0.01, 0.1, 1, 10)
RIDGE_ALPHA = (0.1, 1, 10, 100)
N_BOOTSTRAP = 2000
N_PERMUTATION = 200
# V2 §3.8: matches the frozen Stage-22 table exactly and case-sensitively.
TREATMENT_ORDER = ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin", "Trametinib")
REFERENCE_TREATMENT = "Acid"
REWIND_NUISANCE = ("log1p(n_pretreatment_cells)", "n_lanes")
# V2 §1.2.1: the total-depth term is mandatory -- the frozen confound is keyed on total depth.
WM989_NUISANCE = ("log1p(n_naive_cells)", "log1p(n_naive1_cells)",
                  "log1p(n_naive2_cells)", "log1p(n_naive3_cells)")
N_GENES = 36601
CP10K = 10_000.0

AUDITED = "STAGE22_INPUTS_AUDITED"
BLOCKED = "STAGE22_INPUTS_BLOCKED"
PROTOCOL_FROZEN = "PROTOCOL_FROZEN"
PROTOCOL_BLOCKED = "PROTOCOL_BLOCKED"

# V2 §2.3 + the firewall override rule: Stage 23 may TIGHTEN Stage-22's classes, never loosen them.
STAGE23_TIGHTENED = {
    "post_tie_size": {"stage22": "PROVENANCE_ONLY", "stage23": "TARGET"},
    "treatment_total_assigned_cells": {"stage22": "BASELINE_NUISANCE", "stage23": "TARGET"},
}


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


def canonical_text_sha256(p: Path) -> str:
    """V2 §1.4: text provenance hashes canonicalise CRLF and CR to LF first.

    Stage 22 recorded checkout-byte hashes, so its plan/builder digests differ between Windows and
    Linux for identical content. This rule removes that failure mode; a contract test asserts the
    same text encoded either way yields one digest.
    """
    raw = p.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256_bytes(raw)


def read_barcodes(p: Path) -> list[str]:
    with gzip.open(p, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def gene_features(p: Path) -> list[tuple[str, str]]:
    """(stable_id, symbol) for the Gene Expression block only."""
    out = []
    with gzip.open(p, "rt") as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) > 2 and parts[2] != "Gene Expression":
                continue
            out.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
    return out


# --------------------------------------------------------------------------------------------- #
# 23A step 1 — the §3.0 input audit. Everything downstream is gated on this.
# --------------------------------------------------------------------------------------------- #
def audit_stage22_inputs(rewind_root: Path, wm989_root: Path) -> dict:
    checks: dict[str, object] = {}
    fail: list[str] = []

    def need(name, ok, detail=None):
        checks[name] = {"ok": bool(ok), **({"detail": detail} if detail is not None else {})}
        if not ok:
            fail.append(name)

    res = json.loads((_RESULTS / "stage22_prospective_benchmark_results.json").read_text("utf-8"))
    man = {"GSE227151": json.loads((_RESULTS / "stage22_rewind_benchmark_manifest.json")
                                   .read_text("utf-8")),
           "GSE279162": json.loads((_RESULTS / "stage22_wm989_benchmark_manifest.json")
                                   .read_text("utf-8"))}

    # ---- hash chain: CSVs -> manifests -> results -------------------------------------------- #
    need("six_benchmark_csv_hashes_match_manifests",
         all(sha256_file(_RESULTS / a["name"]) == a["sha256"]
             and (_RESULTS / a["name"]).stat().st_size == a["bytes"]
             for m in man.values() for a in m["derived_artifacts"]))
    need("two_manifest_hashes_match_results",
         all(sha256_file(_RESULTS / a["name"]) == a["sha256"] for a in res["manifest_hashes"]))
    need("results_file_omits_its_own_hash",
         "stage22_prospective_benchmark_results.json"
         not in {a["name"] for a in res["manifest_hashes"]})

    # ---- V2 §1.3: never trust the serialized `overall` string on its own ---------------------- #
    ready = ("BENCHMARK_READY", "BENCHMARK_READY_WITH_DECLARED_MISSINGNESS")
    role_a = res["datasets"]["GSE227151"]["verdict"]
    derived = (role_a in ready and res["all_gates_pass"] is True
               and all(res["gates"].values()) and res["model_fitted"] is False)
    need("role_a_ready", role_a in ready, role_a)
    need("all_gates_pass_true", res["all_gates_pass"] is True)
    need("every_individual_gate_true", all(res["gates"].values()),
         {k: v for k, v in res["gates"].items() if not v} or "all true")
    need("model_fitted_false", res["model_fitted"] is False)
    need("preflight_derived_independently_of_overall_string", derived)
    checks["stage22_overall_string_for_reference_only"] = {"ok": True, "detail": res["overall"]}
    checks["role_b_status_preserved_separately"] = {
        "ok": True, "detail": res["datasets"]["GSE279162"]["verdict"]}
    checks["known_stage22_gate_derivation_limitation"] = {
        "ok": True,
        "detail": "Stage 22 derives `overall` from the Role-A verdict without consuming "
                  "all_gates_pass. Unaffected here because all ten gates are true, and Stage 23 "
                  "fails closed on the derived condition above regardless."}

    # ---- V2 §1.4: inherited CRLF/LF provenance limitation, declared not hidden ---------------- #
    plan22 = ROOT / "plans/(newer)practical plans/STAGE_22_PROSPECTIVE_BENCHMARK_CONSTRUCTION_V2.md"
    bld22 = ROOT / "experiments" / "build_stage22_prospective_benchmarks.py"
    checks["inherited_crlf_lf_provenance_limitation"] = {
        "ok": True,
        "detail": {"stage22_plan": {"recorded": man["GSE227151"]["plan"]["sha256"],
                                    "checkout_bytes": sha256_file(plan22),
                                    "canonical_lf": canonical_text_sha256(plan22)},
                   "stage22_builder": {"recorded": man["GSE227151"]["builder_source_sha256"],
                                       "checkout_bytes": sha256_file(bld22),
                                       "canonical_lf": canonical_text_sha256(bld22)},
                   "note": "content-identical; digests differ only by line endings. Stage-22 "
                           "tables are NOT rewritten to hide this."}}
    lf = b"a\nb\nc\n"
    crlf = lf.replace(b"\n", b"\r\n")
    need("canonical_text_hash_is_lf_crlf_invariant",
         sha256_bytes(lf) == sha256_bytes(crlf.replace(b"\r\n", b"\n").replace(b"\r", b"\n")))

    # ---- frozen counts / folds ---------------------------------------------------------------- #
    rc = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")
    rk = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv")
    wn = pd.read_csv(_RESULTS / "stage22_wm989_naive_cells.csv")
    wk = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv")
    wt = pd.read_csv(_RESULTS / "stage22_wm989_clone_treatment.csv")
    wc = pd.read_csv(_RESULTS / "stage22_wm989_cell_assignments.csv")
    need("rewind_counts",
         (len(rc), len(rk), int(rc.y_primed.sum()), int(rk.y_primed.sum())) == (3905, 3147, 42, 35))
    need("wm989_counts",
         (len(wn), len(wk), len(wt), int((wt.n_post_cells > 0).sum())) == (6489, 1401, 8406, 2256))
    need("rewind_positive_clones_per_fold",
         rk.groupby("outer_fold").y_primed.sum().tolist() == [7, 7, 7, 7, 7])
    need("outer_folds_are_five_each",
         rk.outer_fold.nunique() == N_OUTER and wk.outer_fold.nunique() == N_OUTER)
    need("rewind_cell_fold_matches_clone_fold",
         bool((rc.clone_id.map(dict(zip(rk.clone_id, rk.outer_fold, strict=True)))
               == rc.outer_fold).all()))
    need("wm989_clone_treatment_fold_matches_clone_fold",
         bool((wt.clone_id.map(dict(zip(wk.clone_id, wk.outer_fold, strict=True)))
               == wt.outer_fold).all()))
    c2 = int(wt.loc[wt.n_post_cells > 0, "clone_id"].nunique())
    need("wm989_c2_eligible_clones_929", c2 == 929, c2)
    need("treatments_match_frozen_case_sensitive",
         sorted(wt.treatment.unique()) == list(TREATMENT_ORDER))
    need("reference_treatment_present_case_sensitive", REFERENCE_TREATMENT in set(wt.treatment))
    for col in ("n_naive_cells", "n_naive1_cells", "n_naive2_cells", "n_naive3_cells"):
        need(f"wm989_nuisance_column_{col}", col in wk.columns)
    for col in ("n_pretreatment_cells", "n_lanes"):
        need(f"rewind_nuisance_column_{col}", col in rk.columns)

    # ---- external raw data --------------------------------------------------------------------- #
    for gse, root in (("GSE227151", rewind_root), ("GSE279162", wm989_root)):
        miss, bad = [], []
        for f in man[gse]["source_files"]:
            hits = list(root.rglob(f["name"]))
            if not hits:
                miss.append(f["name"])
                continue
            p = hits[0]
            if p.stat().st_size != f["bytes"] or sha256_file(p) != f["sha256"]:
                bad.append(f["name"])
        need(f"{gse}_raw_files_present_and_identical", not miss and not bad,
             {"missing": miss, "mismatched": bad, "checked": len(man[gse]["source_files"])})

    # ---- every expression_column_index resolves to the recorded barcode ----------------------- #
    unresolved = 0
    for tbl, root in ((rc, rewind_root), (wn, wm989_root)):
        for gsm, sub in tbl.groupby("gsm"):
            B = read_barcodes(next(root.rglob(f"{gsm}_*_barcodes.tsv.gz")))
            idx = sub["expression_column_index"].to_numpy()
            if idx.min() < 1 or idx.max() > len(B):
                unresolved += len(sub)
                continue
            got = np.array([B[i - 1] for i in idx])
            unresolved += int((got != sub["expression_barcode"].to_numpy()).sum())
    need("every_expression_column_index_resolves", unresolved == 0, unresolved)
    need("expression_sources_exist",
         all(bool(list(rewind_root.rglob(s)) or list(wm989_root.rglob(s)))
             for s in set(rc.expression_source) | set(wn.expression_source)))

    # ---- feature universe ----------------------------------------------------------------------- #
    feat_sig, ge_info = {}, {}
    for root in (rewind_root, wm989_root):
        for p in sorted(root.rglob("*_features.tsv.gz")):
            rows = [ln.rstrip("\n").split("\t") for ln in gzip.open(p, "rt")]
            n_custom = sum(1 for r in rows if len(r) > 2 and r[2] == "Custom")
            ge = [r[0] for r in rows if len(r) > 2 and r[2] == "Gene Expression"]
            feat_sig[p.name] = sha256_bytes("\n".join(ge).encode())
            ge_info[p.name] = (len(ge), n_custom,
                               all(len(r) > 2 and r[2] == "Gene Expression" for r in rows[:N_GENES]))
    need("every_sample_has_36601_gene_expression_features",
         all(v[0] == N_GENES for v in ge_info.values()), {k: v[0] for k, v in ge_info.items()})
    need("wm989_custom_block_is_153055_where_present",
         all(v[1] in (0, 153055) for v in ge_info.values()))
    need("gene_block_is_the_first_36601_rows", all(v[2] for v in ge_info.values()))
    need("all_samples_share_one_gene_feature_id_list", len(set(feat_sig.values())) == 1,
         {"distinct_signatures": len(set(feat_sig.values())), "samples": len(feat_sig)})

    # ---- WM989 targets depend on post-treatment lineage only (V2 §1.3) ------------------------ #
    treated = wc[wc.is_assigned & ~wc.is_naive]
    cnt = treated.groupby(["assigned_lineage", "sample"]).size()
    rebuilt = np.array([int(cnt.get(k, 0))
                        for k in zip(wt.clone_id, wt.treatment, strict=True)])
    need("wm989_targets_rebuild_from_treated_cells_only",
         bool((rebuilt == wt.n_post_cells.to_numpy()).all()),
         {"rows": int(len(rebuilt)),
          "max_abs_diff": int(np.abs(rebuilt - wt.n_post_cells.to_numpy()).max()),
          "naive_cells_used": 0})
    checks["wm989_inherited_joint_assignment_dependency"] = {
        "ok": True,
        "detail": "199/39,665 = 0.50% of treated assignments differ if the author pipeline is "
                  "rerun on treated cells only. Declared limitation per V2 §1.2.1; the frozen "
                  "Stage-22 target is NOT rebuilt. No pretreatment Gene Expression value enters "
                  "any target."}

    return {"verdict": AUDITED if not fail else BLOCKED,
            "failed_checks": fail, "n_checks": len(checks), "checks": checks}


# --------------------------------------------------------------------------------------------- #
# 23A step 2 — clone pseudobulk. Raw counts summed FIRST, then CP10K + log1p exactly once.
# --------------------------------------------------------------------------------------------- #
def clone_pseudobulk(cells: pd.DataFrame, roots: list[Path]) -> tuple[sparse.csr_matrix, list[str]]:
    """Sum raw Gene-Expression counts per clone, then normalise once.

    Only rows <= N_GENES are read, so WM989's 153,055 Custom lineage features can never reach `X`.
    Only the pretreatment cells listed in `cells` are read, so no treated column is ever touched.
    """
    clones = sorted(cells["clone_id"].unique())
    cidx = {c: i for i, c in enumerate(clones)}
    acc = sparse.csr_matrix((len(clones), N_GENES), dtype=np.float64)
    for src, sub in cells.groupby("expression_source"):
        mtx = next(p for root in roots for p in root.rglob(src))
        col2clone = dict(zip(sub["expression_column_index"].to_numpy(),
                             [cidx[c] for c in sub["clone_id"]], strict=True))
        wanted = np.zeros(max(col2clone) + 2, dtype=np.int64) - 1
        for c, i in col2clone.items():
            wanted[c] = i
        with gzip.open(mtx, "rt") as fh:
            for line in fh:
                if not line.startswith("%"):
                    break                                   # dims line consumed
            for chunk in pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                                     dtype=np.int64, chunksize=4_000_000):
                r = chunk["row"].to_numpy()
                c = chunk["col"].to_numpy()
                v = chunk["val"].to_numpy()
                keep = (r <= N_GENES) & (c < len(wanted))
                r, c, v = r[keep], c[keep], v[keep]
                tgt = wanted[c]
                keep = tgt >= 0
                if not keep.any():
                    continue
                acc = acc + sparse.csr_matrix(
                    (v[keep].astype(np.float64), (tgt[keep], r[keep] - 1)),
                    shape=(len(clones), N_GENES))
    totals = np.asarray(acc.sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise RuntimeError(f"{int((totals <= 0).sum())} zero-total clones -- V2 §3.3 is a BLOCK, "
                           "not an imputation case")
    norm = acc.multiply((CP10K / totals)[:, None]).tocsr()   # CP10K once
    norm.data = np.log1p(norm.data)                          # log1p once
    return norm, clones


def training_fold_gene_filter(X: sparse.csr_matrix, train_rows: np.ndarray) -> np.ndarray:
    """V2 §3.4, applied to an outer-TRAINING subset only.

    Reported here descriptively; the estimator pipeline refits this inside every inner split too
    (V2 §2.5), so a gene list computed once on the full outer-training set is never reused as a
    shortcut for inner-CV selection.
    """
    sub = X[train_rows]
    detected = np.asarray((sub > 0).sum(axis=0)).ravel()
    floor = max(5, int(np.ceil(0.01 * len(train_rows))))
    dense_ok = detected >= floor
    sq = sub.multiply(sub)
    mean = np.asarray(sub.mean(axis=0)).ravel()
    var = np.asarray(sq.mean(axis=0)).ravel() - mean ** 2
    return np.flatnonzero(dense_ok & (var > 0))


def expression_manifest(name: str, X: sparse.csr_matrix, clones: list[str],
                        cells: pd.DataFrame, feature_ids: list[str]) -> dict:
    """Compact provenance for a matrix that deliberately never enters git."""
    return {
        "dataset": name,
        "clones": len(clones),
        "genes": int(X.shape[1]),
        "nnz": int(X.nnz),
        "density": round(X.nnz / (X.shape[0] * X.shape[1]), 6),
        "pretreatment_cells_summed": int(len(cells)),
        "cells_per_clone": {"min": int(cells.groupby("clone_id").size().min()),
                            "median": int(cells.groupby("clone_id").size().median()),
                            "max": int(cells.groupby("clone_id").size().max())},
        "source_samples": sorted(cells["expression_source"].unique().tolist()),
        "normalization": "sum raw counts per clone -> CP10K -> log1p, applied exactly once",
        "feature_key": "stable 10x feature id (column 1 of features.tsv), never symbol or row index",
        "feature_id_sha256": sha256_bytes("\n".join(feature_ids).encode()),
        "clone_id_sha256": sha256_bytes("\n".join(clones).encode()),
        "matrix_content_sha256": sha256_bytes(
            X.indptr.tobytes() + X.indices.tobytes() + np.round(X.data, 10).tobytes()),
        "committed_to_git": False,
        "cache_path_note": "materialised under _cc_cache/stage23/ which is gitignored; the "
                           "content hash above is what makes a rebuild verifiable",
    }


def build_protocol(tightened: dict) -> dict:
    """V2 §2.4: everything a later substage could otherwise shop for, frozen now."""
    return {
        "stage": "23A",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION,
                 "canonical_lf_sha256": canonical_text_sha256(PLAN)},
        "stage22_benchmark_commit": STAGE22_BENCHMARK_COMMIT,
        "builder_source_canonical_lf_sha256": canonical_text_sha256(Path(__file__).resolve()),
        "hash_rule": {"binary_artifacts": "sha256 of exact bytes",
                      "text_sources": "sha256 after canonicalising CRLF and CR to LF"},
        "seeds": {"protocol_pca_and_inner_cv": SEED_PROTOCOL,
                  "bootstrap_rewind": SEED_BOOT_REWIND,
                  "bootstrap_wm989_c1": SEED_BOOT_WM989_C1,
                  "bootstrap_wm989_c2": SEED_BOOT_WM989_C2,
                  "permutation_base": SEED_PERMUTATION},
        "outer_cv": {"folds": N_OUTER, "unit": "clone",
                     "source": "frozen Stage-22 outer_fold; never recomputed"},
        "inner_cv": {"folds": N_INNER,
                     "rewind": "StratifiedKFold(shuffle=True, random_state=23023) on clones; "
                               "every inner fold must contain both classes",
                     "wm989": "GroupKFold(groups=clone_id) over deterministically sorted clone ids",
                     "wm989_c2": "restricted to outer-training clones with >=1 nonzero row, "
                                 "after the outer split is known"},
        "representation": {
            "pseudobulk": "sum raw pretreatment counts per clone",
            "normalization": "CP10K then log1p, exactly once",
            "gene_filter": "detected in >= max(5, ceil(0.01 * n_training_clones)) and var > 0, "
                           "fitted on the applicable TRAINING split only",
            "pca": {"impl": "sklearn.decomposition.PCA", "svd_solver": "randomized",
                    "random_state": SEED_PROTOCOL,
                    "fit_once_at_max_K_then_reuse_prefixes": True},
            "k_candidates": list(K_CANDIDATES),
            "standardization": "genes standardized train-only; PC scores standardized again "
                               "train-only; continuous nuisance standardized train-only; "
                               "treatment dummies never standardized"},
        "grids": {"logistic": {"penalty": "l2", "solver": "liblinear", "C": list(LOGISTIC_C),
                               "fit_intercept": True, "class_weight": None, "max_iter": 5000,
                               "random_state": SEED_PROTOCOL},
                  "ridge": {"alpha": list(RIDGE_ALPHA), "fit_intercept": True}},
        "tie_break": ["smaller K", "stronger regularization (smaller C / larger alpha)"],
        "treatment_coding": {"canonical_order": list(TREATMENT_ORDER),
                             "reference": REFERENCE_TREATMENT,
                             "n_dummies": len(TREATMENT_ORDER) - 1,
                             "case_sensitive_match_to_stage22": True,
                             "standardized": False},
        "nuisance_blocks": {"rewind": list(REWIND_NUISANCE), "wm989": list(WM989_NUISANCE)},
        "inner_selection_metric": {"rewind": "maximize mean average_precision_score",
                                   "wm989_c1": "minimize mean clone-balanced log loss",
                                   "wm989_c2": "minimize mean clone-balanced MAE"},
        "c2_sample_weight": "1 / n_nonzero_rows_for_that_clone_in_the_training_subset, "
                            "renormalised to mean 1",
        "inference": {"bootstrap_replicates": N_BOOTSTRAP, "permutations": N_PERMUTATION,
                      "permutation_p": "(1 + #{null >= observed}) / (n_perm + 1)"},
        "feature_firewall": {
            "primary_x": "pretreatment Gene Expression features only",
            "wm989_custom_lineage_features_forbidden": True,
            "stage23_tightened_vs_stage22": tightened,
            "rule": "Stage 23 may tighten eligibility; it may never loosen a Stage-22 "
                    "forbidden/target/provenance field into PRIMARY_X"},
        "model_fitted": False,
    }


# --------------------------------------------------------------------------------------------- #
# 23B — Rewind Role-A learnability.  Does pretreatment X predict priming beyond prevalence and
# captured clone size?  This is the first substage that fits an estimator.
# --------------------------------------------------------------------------------------------- #
REWIND_OOF = _RESULTS / "stage23_rewind_oof_predictions.csv"
REWIND_RESULTS = _RESULTS / "stage23_rewind_results.json"
ROLE_A_PASS = "ROLE_A_SIGNAL_PASS"
ROLE_A_WEAK = "ROLE_A_SIGNAL_WEAK"
ROLE_A_FAIL = "ROLE_A_SIGNAL_FAIL"


def _load_rewind_x() -> tuple[sparse.csr_matrix, list[str]]:
    p = _CACHE / "GSE227151_pseudobulk.npz"
    if not p.exists():
        raise RuntimeError("23A pseudobulk cache missing -- run `--stage 23a` first")
    X = sparse.load_npz(p)
    clones = json.loads((_CACHE / "GSE227151_clones.json").read_text(encoding="utf-8"))
    man = json.loads((_RESULTS / "stage23_rewind_clone_expression_manifest.json")
                     .read_text(encoding="utf-8"))
    got = sha256_bytes(X.indptr.tobytes() + X.indices.tobytes() + np.round(X.data, 10).tobytes())
    if got != man["matrix_content_sha256"]:
        raise RuntimeError("cached X does not match the 23A manifest hash")
    return X, clones


def expression_block(X, train_rows: np.ndarray, test_rows: np.ndarray, max_k: int):
    """V2 §2.5/§3.4-3.5: filter, standardize, PCA and re-standardize -- all fitted on TRAIN only.

    Returns standardized PC scores for train and test. PCA is fitted once at `max_k`; K=10/20/50
    are prefixes of this same basis, never separate randomized fits.
    """
    from sklearn.decomposition import PCA

    keep = training_fold_gene_filter(X, train_rows)
    tr = np.asarray(X[train_rows][:, keep].todense())
    te = np.asarray(X[test_rows][:, keep].todense())
    mu, sd = tr.mean(axis=0), tr.std(axis=0)
    sd[sd == 0] = 1.0
    tr = (tr - mu) / sd
    te = (te - mu) / sd
    k = min(max_k, tr.shape[0], tr.shape[1])
    pca = PCA(n_components=k, svd_solver="randomized", random_state=SEED_PROTOCOL).fit(tr)
    ztr, zte = pca.transform(tr), pca.transform(te)
    zmu, zsd = ztr.mean(axis=0), ztr.std(axis=0)
    zsd[zsd == 0] = 1.0
    return (ztr - zmu) / zsd, (zte - zmu) / zsd, int(len(keep)), k


def standardize_train_only(tr: np.ndarray, te: np.ndarray):
    mu, sd = tr.mean(axis=0), tr.std(axis=0)
    sd[sd == 0] = 1.0
    return (tr - mu) / sd, (te - mu) / sd


def _fit_logistic(Xtr, ytr, Xte, C, convergence: list):
    import warnings

    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = LogisticRegression(penalty="l2", solver="liblinear", C=C, fit_intercept=True,
                               class_weight=None, max_iter=5000, random_state=SEED_PROTOCOL)
        m.fit(Xtr, ytr)
        for w in caught:
            if issubclass(w.category, ConvergenceWarning):
                convergence.append({"C": C, "n_train": int(len(ytr)), "message": str(w.message)})
    return m.predict_proba(Xte)[:, 1]


def _design(pcs, nuis, k, use_x, use_nuis):
    parts = []
    if use_x:
        parts.append(pcs[:, :k])
    if use_nuis:
        parts.append(nuis)
    return np.hstack(parts)


def run_23b(rewind_root: Path) -> dict:
    from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    X, clones = _load_rewind_x()
    tbl = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    y = tbl["y_primed"].to_numpy()
    fold = tbl["outer_fold"].to_numpy()
    nuis_raw = np.column_stack([np.log1p(tbl["n_pretreatment_cells"].to_numpy()),
                                tbl["n_lanes"].to_numpy().astype(float)])

    convergence: list = []
    oof = {m: np.full(len(clones), np.nan) for m in ("R0", "R1", "R2", "R3")}
    selected: dict[str, dict] = {}
    per_fold_meta = {}

    for f in range(N_OUTER):
        tr = np.flatnonzero(fold != f)
        te = np.flatnonzero(fold == f)
        ytr, yte = y[tr], y[te]

        # ---- inner CV: every learned quantity refitted inside each inner-training split ------ #
        skf = StratifiedKFold(n_splits=N_INNER, shuffle=True, random_state=SEED_PROTOCOL)
        scores = {"R1": {}, "R2": {}, "R3": {}}
        for itr_i, iva_i in skf.split(np.zeros(len(tr)), ytr):
            itr, iva = tr[itr_i], tr[iva_i]
            ztr, zva, _, kmax = expression_block(X, itr, iva, max(K_CANDIDATES))
            btr, bva = standardize_train_only(nuis_raw[itr], nuis_raw[iva])
            yi, yv = y[itr], y[iva]
            for C in LOGISTIC_C:
                p = _fit_logistic(btr, yi, bva, C, convergence)
                scores["R1"].setdefault((None, C), []).append(average_precision_score(yv, p))
            for k in K_CANDIDATES:
                if k > kmax:
                    continue
                for C in LOGISTIC_C:
                    p2 = _fit_logistic(_design(ztr, btr, k, True, False), yi,
                                       _design(zva, bva, k, True, False), C, convergence)
                    scores["R2"].setdefault((k, C), []).append(average_precision_score(yv, p2))
                    p3 = _fit_logistic(_design(ztr, btr, k, True, True), yi,
                                       _design(zva, bva, k, True, True), C, convergence)
                    scores["R3"].setdefault((k, C), []).append(average_precision_score(yv, p3))

        # ---- deterministic selection: maximise mean AP; tie-break smaller K, then smaller C -- #
        def pick(model, sc=scores):   # bound now:  is rebuilt each outer fold
            best = max(sc[model].items(),
                       key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                       -(kv[0][0] or 0), -kv[0][1]))
            return best[0], float(np.mean(best[1]))

        sel = {m: pick(m) for m in ("R1", "R2", "R3")}
        selected[str(f)] = {m: {"K": sel[m][0][0], "C": sel[m][0][1],
                                "mean_inner_AP": round(sel[m][1], 6)} for m in sel}

        # ---- final: rebuild the pipeline on the WHOLE outer-training set, predict test once --- #
        ztr, zte, n_genes_kept, kmax = expression_block(X, tr, te, max(K_CANDIDATES))
        btr, bte = standardize_train_only(nuis_raw[tr], nuis_raw[te])
        oof["R0"][te] = float(ytr.mean())
        oof["R1"][te] = _fit_logistic(btr, ytr, bte, sel["R1"][0][1], convergence)
        oof["R2"][te] = _fit_logistic(_design(ztr, btr, sel["R2"][0][0], True, False), ytr,
                                      _design(zte, bte, sel["R2"][0][0], True, False),
                                      sel["R2"][0][1], convergence)
        oof["R3"][te] = _fit_logistic(_design(ztr, btr, sel["R3"][0][0], True, True), ytr,
                                      _design(zte, bte, sel["R3"][0][0], True, True),
                                      sel["R3"][0][1], convergence)
        per_fold_meta[str(f)] = {"train_clones": int(len(tr)), "test_clones": int(len(te)),
                                 "test_positives": int(yte.sum()),
                                 "retained_genes": n_genes_kept, "max_feasible_K": kmax,
                                 "train_prevalence": round(float(ytr.mean()), 6)}

    assert not np.isnan(np.concatenate([oof[m] for m in oof])).any(), "an OOF prediction is missing"

    # ---- metrics ----------------------------------------------------------------------------- #
    def metrics(p):
        return {"AP": float(average_precision_score(y, p)),
                "ROC_AUC": float(roc_auc_score(y, p)),
                "log_loss": float(log_loss(y, np.clip(p, 1e-15, 1 - 1e-15))),
                "brier": float(brier_score_loss(y, p))}

    pooled = {m: metrics(oof[m]) for m in ("R0", "R1", "R2", "R3")}
    per_fold_ap = {}
    for f in range(N_OUTER):
        te = np.flatnonzero(fold == f)
        per_fold_ap[str(f)] = {m: float(average_precision_score(y[te], oof[m][te]))
                               for m in ("R0", "R1", "R2", "R3")}
        per_fold_ap[str(f)]["delta_AP_state"] = (per_fold_ap[str(f)]["R3"]
                                                 - per_fold_ap[str(f)]["R1"])
        per_fold_ap[str(f)]["delta_AP_absolute"] = (per_fold_ap[str(f)]["R2"]
                                                    - per_fold_ap[str(f)]["R0"])

    # ---- stratified clone bootstrap on the pooled OOF (V2 §4.4) ------------------------------ #
    rng = np.random.default_rng(SEED_BOOT_REWIND)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    d_state, d_abs = np.empty(N_BOOTSTRAP), np.empty(N_BOOTSTRAP)
    for b in range(N_BOOTSTRAP):
        idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                              rng.choice(neg, len(neg), replace=True)])
        yb = y[idx]
        d_state[b] = (average_precision_score(yb, oof["R3"][idx])
                      - average_precision_score(yb, oof["R1"][idx]))
        d_abs[b] = (average_precision_score(yb, oof["R2"][idx])
                    - average_precision_score(yb, oof["R0"][idx]))

    def boot(delta, point):
        lo, hi = np.percentile(delta, [2.5, 97.5])
        return {"point": float(point), "ci95_low": float(lo), "ci95_high": float(hi),
                "fraction_delta_le_0": float((delta <= 0).mean()),
                "bootstrap_mean": float(delta.mean()), "replicates": N_BOOTSTRAP,
                "seed": SEED_BOOT_REWIND}

    d_state_pt = pooled["R3"]["AP"] - pooled["R1"]["AP"]
    d_abs_pt = pooled["R2"]["AP"] - pooled["R0"]["AP"]
    inference = {"delta_AP_state_R3_minus_R1": boot(d_state, d_state_pt),
                 "delta_AP_absolute_R2_minus_R0": boot(d_abs, d_abs_pt)}

    # ---- provisional verdict, derived (V2 §4.6) ---------------------------------------------- #
    s = inference["delta_AP_state_R3_minus_R1"]
    if s["point"] <= 0:
        verdict = ROLE_A_FAIL
    elif s["ci95_low"] > 0:
        verdict = ROLE_A_PASS
    else:
        verdict = ROLE_A_WEAK

    pd.DataFrame({"clone_id": clones, "outer_fold": fold, "y_primed": y,
                  **{f"pred_{m}": oof[m] for m in ("R0", "R1", "R2", "R3")}}
                 ).to_csv(REWIND_OOF, index=False, lineterminator="\n")

    out = {
        "stage": "23B", "dataset": "GSE227151", "role": "A",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION,
                 "canonical_lf_sha256": canonical_text_sha256(PLAN)},
        "protocol_sha256": sha256_file(_RESULTS / "stage23_protocol.json"),
        "clones": len(clones), "positives": int(y.sum()), "negatives": int((y == 0).sum()),
        "models": {"R0": "outer-training prevalence", "R1": "nuisance only",
                   "R2": "PCA(X) only", "R3": "PCA(X) + nuisance"},
        "primary_metric": "average_precision_score at clone grain",
        "selected_hyperparameters_per_outer_fold": selected,
        "per_outer_fold": per_fold_meta,
        "pooled_oof_metrics": pooled,
        "per_fold_average_precision": per_fold_ap,
        "fold_direction_is_diagnostic_only": True,
        "inference": inference,
        "convergence_warnings": convergence,
        "provisional_verdict": verdict,
        "verdict_is_provisional_until": "23E structural controls + ROLE_A_PERMUTATION_PASS",
    }
    REWIND_RESULTS.write_text(json.dumps(out, indent=2), encoding="utf-8", newline="\n")
    return out


# --------------------------------------------------------------------------------------------- #
# 23C — WM989 additive state-signal gate.  Does pretreatment X add beyond treatment U and captured
# naive clone abundance B?  No interaction terms here: those are 23D.
# --------------------------------------------------------------------------------------------- #
WM989_C1_OOF = _RESULTS / "stage23_wm989_detection_oof.csv"
WM989_C2_OOF = _RESULTS / "stage23_wm989_abundance_oof.csv"
WM989_RESULTS = _RESULTS / "stage23_wm989_results.json"
ROLE_B_PASS = "ROLE_B_ADDITIVE_PASS"
ROLE_B_WEAK = "ROLE_B_ADDITIVE_WEAK"
ROLE_B_FAIL = "ROLE_B_ADDITIVE_FAIL"
WM989_MODELS = ("W0", "W1", "W2", "W3", "W4")


def _load_wm989_x():
    p = _CACHE / "GSE279162_pseudobulk.npz"
    if not p.exists():
        raise RuntimeError("23A pseudobulk cache missing -- run `--stage 23a` first")
    X = sparse.load_npz(p)
    clones = json.loads((_CACHE / "GSE279162_clones.json").read_text(encoding="utf-8"))
    man = json.loads((_RESULTS / "stage23_wm989_clone_expression_manifest.json")
                     .read_text(encoding="utf-8"))
    got = sha256_bytes(X.indptr.tobytes() + X.indices.tobytes() + np.round(X.data, 10).tobytes())
    if got != man["matrix_content_sha256"]:
        raise RuntimeError("cached X does not match the 23A manifest hash")
    return X, clones


def treatment_dummies(treatments: np.ndarray) -> np.ndarray:
    """Five non-reference dummies in the frozen canonical order; `Acid` is the reference.

    Never standardized (V2 §3.5). Column order is fixed by TREATMENT_ORDER, not by whatever order
    the rows happen to arrive in.
    """
    non_ref = [t for t in TREATMENT_ORDER if t != REFERENCE_TREATMENT]
    return np.column_stack([(treatments == t).astype(float) for t in non_ref])


def clone_balanced_error(err: np.ndarray, clone_key: np.ndarray, square: bool = False):
    """Per-clone mean first, then the mean over clones -- so a clone observed under many
    treatments does not silently become several independent units (V2 §5.2)."""
    e = err ** 2 if square else np.abs(err)
    df = pd.DataFrame({"c": clone_key, "e": e})
    per_clone = df.groupby("c", sort=True)["e"].mean()
    return per_clone


def _fit_ridge(Xtr, ytr, Xte, alpha, weights=None):
    from sklearn.linear_model import Ridge

    m = Ridge(alpha=alpha, fit_intercept=True)
    m.fit(Xtr, ytr, sample_weight=weights)
    return m.predict(Xte)


def run_23c(wm989_root: Path) -> dict:
    from scipy.stats import spearmanr
    from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
    from sklearn.model_selection import GroupKFold

    X, clones = _load_wm989_x()
    clone_pos = {c: i for i, c in enumerate(clones)}
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    ct = pd.read_csv(_RESULTS / "stage22_wm989_clone_treatment.csv")
    ct = ct.sort_values(["clone_id", "treatment"]).reset_index(drop=True)

    fold_of = ck["outer_fold"].to_dict()
    nuis_clone = np.column_stack([np.log1p(ck[c].to_numpy(dtype=float))
                                  for c in ("n_naive_cells", "n_naive1_cells",
                                            "n_naive2_cells", "n_naive3_cells")])
    row_fold = ct["clone_id"].map(fold_of).to_numpy()
    dummies_all = treatment_dummies(ct["treatment"].to_numpy())
    y_c1 = (ct["n_post_cells"].to_numpy() > 0).astype(int)
    y_c2 = np.log1p(ct["n_post_cells"].to_numpy(dtype=float))
    c2_mask = ct["n_post_cells"].to_numpy() > 0

    convergence: list = []
    out: dict = {}

    for endpoint in ("C1", "C2"):
        rows = np.arange(len(ct)) if endpoint == "C1" else np.flatnonzero(c2_mask)
        y = y_c1 if endpoint == "C1" else y_c2
        elig_clones = sorted(set(ct["clone_id"].to_numpy()[rows]))
        oof = {m: np.full(len(rows), np.nan) for m in WM989_MODELS}
        selected: dict = {}
        meta: dict = {}

        for f in range(N_OUTER):
            te_rows_local = np.flatnonzero(row_fold[rows] == f)
            tr_rows_local = np.flatnonzero(row_fold[rows] != f)
            tr_clones = sorted({ct["clone_id"].iloc[rows[i]] for i in tr_rows_local})
            te_clones = sorted({ct["clone_id"].iloc[rows[i]] for i in te_rows_local})

            def prep(fit_clones, apply_sets):
                """Fit filter/scale/PCA/nuisance-scale on `fit_clones` only, apply to each set."""
                fit_idx = np.array([clone_pos[c] for c in fit_clones])
                all_idx = np.array([clone_pos[c] for c in sum(apply_sets, [])])
                ztr, zall, n_genes, kmax = expression_block(X, fit_idx, all_idx,
                                                            max(K_CANDIDATES))
                btr, ball = standardize_train_only(nuis_clone[fit_idx], nuis_clone[all_idx])
                pcs = {c: ztr[i] for i, c in enumerate(fit_clones)}
                nui = {c: btr[i] for i, c in enumerate(fit_clones)}
                for i, c in enumerate(sum(apply_sets, [])):
                    pcs[c] = zall[i]
                    nui[c] = ball[i]
                return pcs, nui, n_genes, kmax

            # ---- inner CV: refit every learned quantity inside each inner-training split ----- #
            scores: dict = {m: {} for m in WM989_MODELS}
            gkf = GroupKFold(n_splits=N_INNER)
            g = np.array(tr_clones)
            for itr_i, iva_i in gkf.split(g, groups=g):
                itr_c, iva_c = [g[i] for i in itr_i], [g[i] for i in iva_i]
                pcs, nui, _, kmax = prep(itr_c, [iva_c])
                sel_i = np.array([i for i in tr_rows_local
                                  if ct["clone_id"].iloc[rows[i]] in set(itr_c)])
                sel_v = np.array([i for i in tr_rows_local
                                  if ct["clone_id"].iloc[rows[i]] in set(iva_c)])
                ci = ct["clone_id"].to_numpy()[rows[sel_i]]
                cv = ct["clone_id"].to_numpy()[rows[sel_v]]
                Pi = np.array([pcs[c] for c in ci])
                Pv = np.array([pcs[c] for c in cv])
                Bi = np.array([nui[c] for c in ci])
                Bv = np.array([nui[c] for c in cv])
                Ui, Uv = dummies_all[rows[sel_i]], dummies_all[rows[sel_v]]
                yi, yv = y[rows[sel_i]], y[rows[sel_v]]
                w = None
                if endpoint == "C2":
                    cnt = pd.Series(ci).value_counts()
                    w = 1.0 / pd.Series(ci).map(cnt).to_numpy()
                    w = w / w.mean()

                def score(pred, ep=endpoint, yy=yv, cc=cv):
                    if ep == "C1":
                        return log_loss(yy, np.clip(pred, 1e-15, 1 - 1e-15), labels=[0, 1])
                    return float(clone_balanced_error(yy - pred, cc).mean())

                specs = {"W0": (False, False, True), "W1": (False, True, True),
                         "W2": (True, False, False), "W3": (True, False, True),
                         "W4": (True, True, True)}
                for m, (ux, ub, uu) in specs.items():
                    ks = K_CANDIDATES if ux else (None,)
                    grid = LOGISTIC_C if endpoint == "C1" else RIDGE_ALPHA
                    for k in ks:
                        if k is not None and k > kmax:
                            continue
                        Ai = np.hstack([q for q, use in
                                        ((Pi[:, :k] if k else None, ux), (Bi, ub), (Ui, uu)) if use])
                        Av = np.hstack([q for q, use in
                                        ((Pv[:, :k] if k else None, ux), (Bv, ub), (Uv, uu)) if use])
                        for hp in grid:
                            if endpoint == "C1":
                                pred = _fit_logistic(Ai, yi, Av, hp, convergence)
                            else:
                                pred = _fit_ridge(Ai, yi, Av, hp, w)
                            scores[m].setdefault((k, hp), []).append(score(pred))

            def pick(model, sc=scores, ep=endpoint):
                # minimise for both endpoints; tie-break smaller K then stronger regularisation
                return min(sc[model].items(),
                           key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                           kv[0][0] or 0,
                                           -kv[0][1] if ep == "C2" else kv[0][1]))

            sel = {m: pick(m) for m in WM989_MODELS}
            selected[str(f)] = {m: {"K": sel[m][0][0], "hp": sel[m][0][1],
                                    "mean_inner_score": round(float(np.mean(sel[m][1])), 6)}
                                for m in WM989_MODELS}

            # ---- final refit on the whole outer-training set, predict the test rows once ----- #
            pcs, nui, n_genes, kmax = prep(tr_clones, [te_clones])
            ctr = ct["clone_id"].to_numpy()[rows[tr_rows_local]]
            cte = ct["clone_id"].to_numpy()[rows[te_rows_local]]
            Ptr = np.array([pcs[c] for c in ctr])
            Pte = np.array([pcs[c] for c in cte])
            Btr = np.array([nui[c] for c in ctr])
            Bte = np.array([nui[c] for c in cte])
            Utr, Ute = dummies_all[rows[tr_rows_local]], dummies_all[rows[te_rows_local]]
            ytr = y[rows[tr_rows_local]]
            w = None
            if endpoint == "C2":
                cnt = pd.Series(ctr).value_counts()
                w = 1.0 / pd.Series(ctr).map(cnt).to_numpy()
                w = w / w.mean()
            specs = {"W0": (False, False, True), "W1": (False, True, True),
                     "W2": (True, False, False), "W3": (True, False, True),
                     "W4": (True, True, True)}
            for m, (ux, ub, uu) in specs.items():
                k, hp = sel[m][0]
                Atr = np.hstack([q for q, use in
                                 ((Ptr[:, :k] if k else None, ux), (Btr, ub), (Utr, uu)) if use])
                Ate = np.hstack([q for q, use in
                                 ((Pte[:, :k] if k else None, ux), (Bte, ub), (Ute, uu)) if use])
                if endpoint == "C1":
                    oof[m][te_rows_local] = _fit_logistic(Atr, ytr, Ate, hp, convergence)
                else:
                    oof[m][te_rows_local] = _fit_ridge(Atr, ytr, Ate, hp, w)
            meta[str(f)] = {"train_clones": len(tr_clones), "test_clones": len(te_clones),
                            "train_rows": int(len(tr_rows_local)),
                            "test_rows": int(len(te_rows_local)),
                            "retained_genes": n_genes, "max_feasible_K": kmax}

        assert not np.isnan(np.concatenate([oof[m] for m in WM989_MODELS])).any()

        yv = y[rows]
        ckey = ct["clone_id"].to_numpy()[rows]
        tkey = ct["treatment"].to_numpy()[rows]
        if endpoint == "C1":
            metrics = {m: {"log_loss": float(log_loss(yv, np.clip(oof[m], 1e-15, 1 - 1e-15))),
                           "AP": float(average_precision_score(yv, oof[m])),
                           "ROC_AUC": float(roc_auc_score(yv, oof[m])),
                           "brier": float(brier_score_loss(yv, oof[m]))} for m in WM989_MODELS}
            # per-clone mean log loss; every clone has all six treatment rows, so the row
            # average and the clone-balanced average coincide here (V2 §3.6)
            def _rowwise_ll(pred, yy=yv):
                pc = np.clip(pred, 1e-15, 1 - 1e-15)
                return -(yy * np.log(pc) + (1 - yy) * np.log(1 - pc))

            per_clone = {m: clone_balanced_error(_rowwise_ll(oof[m]), ckey)
                         for m in WM989_MODELS}
        else:
            metrics = {}
            per_clone = {}
            for m in WM989_MODELS:
                mae_c = clone_balanced_error(yv - oof[m], ckey)
                mse_c = clone_balanced_error(yv - oof[m], ckey, square=True)
                sp = {t: float(spearmanr(yv[tkey == t], oof[m][tkey == t]).statistic)
                      for t in TREATMENT_ORDER}
                metrics[m] = {"clone_balanced_MAE": float(mae_c.mean()),
                              "clone_balanced_RMSE": float(np.sqrt(mse_c.mean())),
                              "per_treatment_spearman": {k: round(v, 4) for k, v in sp.items()},
                              "mean_treatment_spearman": float(np.mean(list(sp.values())))}
                per_clone[m] = mae_c

        # ---- clone-cluster bootstrap: resample clones, carry all of their rows -------------- #
        seed = SEED_BOOT_WM989_C1 if endpoint == "C1" else SEED_BOOT_WM989_C2
        rng = np.random.default_rng(seed)
        uc = per_clone["W1"].index.to_numpy()
        v1, v4 = per_clone["W1"].to_numpy(), per_clone["W4"].to_numpy()
        deltas = np.empty(N_BOOTSTRAP)
        for b in range(N_BOOTSTRAP):
            idx = rng.integers(0, len(uc), len(uc))
            deltas[b] = v1[idx].mean() - v4[idx].mean()
        key = "log_loss" if endpoint == "C1" else "clone_balanced_MAE"
        point = metrics["W1"][key] - metrics["W4"][key]
        lo95, hi95 = np.percentile(deltas, [2.5, 97.5])
        lo975, hi975 = np.percentile(deltas, [1.25, 98.75])
        inference = {"delta_state_W1_minus_W4": {
            "metric": key, "point": float(point),
            "ci95": [float(lo95), float(hi95)],
            "ci975_two_sided": [float(lo975), float(hi975)],
            "fraction_delta_le_0": float((deltas <= 0).mean()),
            "bootstrap_mean": float(deltas.mean()),
            "replicates": N_BOOTSTRAP, "seed": seed,
            "bootstrap_unit": "clone", "clones_resampled": int(len(uc))}}

        frame = pd.DataFrame({"clone_id": ckey, "treatment": tkey, "outer_fold": row_fold[rows],
                              "y": yv, **{f"pred_{m}": oof[m] for m in WM989_MODELS}})
        frame.to_csv(WM989_C1_OOF if endpoint == "C1" else WM989_C2_OOF,
                     index=False, lineterminator="\n")
        out[endpoint] = {"rows": int(len(rows)), "clones": len(elig_clones),
                         "selected_hyperparameters_per_outer_fold": selected,
                         "per_outer_fold": meta, "pooled_oof_metrics": metrics,
                         "inference": inference}

    # ---- derived verdict (V2 §5.7) ------------------------------------------------------------ #
    ll = out["C1"]["inference"]["delta_state_W1_minus_W4"]
    ma = out["C2"]["inference"]["delta_state_W1_minus_W4"]
    pass_ll, pass_ma = ll["ci975_two_sided"][0] > 0, ma["ci975_two_sided"][0] > 0
    harm_ll, harm_ma = ll["ci975_two_sided"][1] < 0, ma["ci975_two_sided"][1] < 0
    if (pass_ll and not harm_ma) or (pass_ma and not harm_ll):
        verdict = ROLE_B_PASS
    elif ll["point"] <= 0 and ma["point"] <= 0:
        verdict = ROLE_B_FAIL
    elif harm_ll or harm_ma:
        verdict = ROLE_B_FAIL
    else:
        verdict = ROLE_B_WEAK

    res = {
        "stage": "23C", "dataset": "GSE279162", "role": "B",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION,
                 "canonical_lf_sha256": canonical_text_sha256(PLAN)},
        "protocol_sha256": sha256_file(_RESULTS / "stage23_protocol.json"),
        "models": {"W0": "U", "W1": "B + U", "W2": "X", "W3": "X + U", "W4": "X + B + U"},
        "nuisance_block_B": list(WM989_NUISANCE),
        "treatment_coding": {"order": list(TREATMENT_ORDER), "reference": REFERENCE_TREATMENT},
        "primary_comparison": "W4 vs W1 on both endpoints",
        "interaction_terms_present": False,
        "endpoints": out,
        "convergence_warnings": convergence,
        "verdict": verdict,
        "verdict_is_provisional_until": "23E structural controls + "
                                        "ROLE_B_ADDITIVE_PERMUTATION_PASS",
    }
    WM989_RESULTS.write_text(json.dumps(res, indent=2), encoding="utf-8", newline="\n")
    return res


# --------------------------------------------------------------------------------------------- #
# 23D — WM989 explicit interaction gate.  Does the contribution of pretreatment state DEPEND on
# treatment?  W1 and W4 are reused verbatim from the frozen 23C out-of-fold predictions, so the
# reference cannot drift; only W5 is fitted here.
# --------------------------------------------------------------------------------------------- #
WM989_INTERACTION_RESULTS = _RESULTS / "stage23_wm989_interaction_results.json"
WM989_W5_OOF = _RESULTS / "stage23_wm989_interaction_oof.csv"
INTERACTION_PASS = "INTERACTION_PASS_MULTI_TREATMENT"
INTERACTION_LOCAL = "INTERACTION_LOCAL_ONLY"
INTERACTION_NONE = "INTERACTION_NOT_SUPPORTED"


def interaction_block(pcs: np.ndarray, dummies: np.ndarray) -> np.ndarray:
    """V2 §6.1/§3.8: only standardized PC score x non-reference treatment dummy.

    Never a gene-level interaction matrix. With K PCs and five non-reference dummies this is
    exactly 5K columns, and the reference treatment's state contribution stays in the common X
    coefficients.
    """
    return np.hstack([pcs * dummies[:, [t]] for t in range(dummies.shape[1])])


def run_23d(wm989_root: Path) -> dict:
    from sklearn.metrics import log_loss
    from sklearn.model_selection import GroupKFold

    X, clones = _load_wm989_x()
    clone_pos = {c: i for i, c in enumerate(clones)}
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis_clone = np.column_stack([np.log1p(ck[c].to_numpy(dtype=float))
                                  for c in ("n_naive_cells", "n_naive1_cells",
                                            "n_naive2_cells", "n_naive3_cells")])
    convergence: list = []
    out: dict = {}

    for endpoint, oof_path in (("C1", WM989_C1_OOF), ("C2", WM989_C2_OOF)):
        base = pd.read_csv(oof_path)
        if not {"pred_W1", "pred_W4"} <= set(base.columns):
            raise RuntimeError(f"{oof_path.name} lacks the frozen 23C W1/W4 predictions")
        ckey = base["clone_id"].to_numpy()
        tkey = base["treatment"].to_numpy()
        yv = base["y"].to_numpy()
        rf = base["outer_fold"].to_numpy()
        dummies = treatment_dummies(tkey)
        w5 = np.full(len(base), np.nan)
        selected, meta = {}, {}

        for f in range(N_OUTER):
            te = np.flatnonzero(rf == f)
            tr = np.flatnonzero(rf != f)
            tr_clones = sorted(set(ckey[tr]))
            te_clones = sorted(set(ckey[te]))

            def prep(fit_clones, apply_clones):
                fit_idx = np.array([clone_pos[c] for c in fit_clones])
                app_idx = np.array([clone_pos[c] for c in apply_clones])
                ztr, zap, n_genes, kmax = expression_block(X, fit_idx, app_idx,
                                                           max(K_CANDIDATES))
                btr, bap = standardize_train_only(nuis_clone[fit_idx], nuis_clone[app_idx])
                pcs = {c: ztr[i] for i, c in enumerate(fit_clones)}
                nui = {c: btr[i] for i, c in enumerate(fit_clones)}
                for i, c in enumerate(apply_clones):
                    pcs[c] = zap[i]
                    nui[c] = bap[i]
                return pcs, nui, n_genes, kmax

            def design(idx, pcs, nui, k, cc=ckey, dd=dummies):
                P = np.array([pcs[c] for c in cc[idx]])[:, :k]
                B = np.array([nui[c] for c in cc[idx]])
                U = dd[idx]
                return np.hstack([P, B, U, interaction_block(P, U)])

            # ---- inner CV for W5 only; W1/W4 are frozen from 23C ---------------------------- #
            scores: dict = {}
            g = np.array(tr_clones)
            for itr_i, iva_i in GroupKFold(n_splits=N_INNER).split(g, groups=g):
                itr_c, iva_c = [g[i] for i in itr_i], [g[i] for i in iva_i]
                pcs, nui, _, kmax = prep(itr_c, iva_c)
                si = np.array([i for i in tr if ckey[i] in set(itr_c)])
                sv = np.array([i for i in tr if ckey[i] in set(iva_c)])
                wts = None
                if endpoint == "C2":
                    cnt = pd.Series(ckey[si]).value_counts()
                    wts = 1.0 / pd.Series(ckey[si]).map(cnt).to_numpy()
                    wts = wts / wts.mean()
                grid = LOGISTIC_C if endpoint == "C1" else RIDGE_ALPHA
                for k in K_CANDIDATES:
                    if k > kmax:
                        continue
                    Ai, Av = design(si, pcs, nui, k), design(sv, pcs, nui, k)
                    for hp in grid:
                        if endpoint == "C1":
                            pred = _fit_logistic(Ai, yv[si], Av, hp, convergence)
                            sc = log_loss(yv[sv], np.clip(pred, 1e-15, 1 - 1e-15), labels=[0, 1])
                        else:
                            pred = _fit_ridge(Ai, yv[si], Av, hp, wts)
                            sc = float(clone_balanced_error(yv[sv] - pred, ckey[sv]).mean())
                        scores.setdefault((k, hp), []).append(sc)

            best = min(scores.items(),
                       key=lambda kv: (round(float(np.mean(kv[1])), 12), kv[0][0],
                                       -kv[0][1] if endpoint == "C2" else kv[0][1]))
            k, hp = best[0]
            selected[str(f)] = {"K": k, "hp": hp,
                                "mean_inner_score": round(float(np.mean(best[1])), 6),
                                "interaction_columns": 5 * k}

            pcs, nui, n_genes, kmax = prep(tr_clones, te_clones)
            wts = None
            if endpoint == "C2":
                cnt = pd.Series(ckey[tr]).value_counts()
                wts = 1.0 / pd.Series(ckey[tr]).map(cnt).to_numpy()
                wts = wts / wts.mean()
            Atr, Ate = design(tr, pcs, nui, k), design(te, pcs, nui, k)
            if endpoint == "C1":
                w5[te] = _fit_logistic(Atr, yv[tr], Ate, hp, convergence)
            else:
                w5[te] = _fit_ridge(Atr, yv[tr], Ate, hp, wts)
            meta[str(f)] = {"train_rows": int(len(tr)), "test_rows": int(len(te)),
                            "design_columns": int(Atr.shape[1]), "retained_genes": n_genes}

        assert not np.isnan(w5).any()

        # ---- per-clone losses for the three models, then the two required comparisons ------- #
        def per_clone(pred, ep=endpoint, yy=yv, cc=ckey):
            if ep == "C1":
                pc = np.clip(pred, 1e-15, 1 - 1e-15)
                return clone_balanced_error(-(yy * np.log(pc) + (1 - yy) * np.log(1 - pc)), cc)
            return clone_balanced_error(yy - pred, cc)

        pcs_ = {m: per_clone(base[f"pred_{m}"].to_numpy()) for m in ("W1", "W4")}
        pcs_["W5"] = per_clone(w5)
        pooled = {m: float(v.mean()) for m, v in pcs_.items()}

        seed = SEED_BOOT_WM989_C1 if endpoint == "C1" else SEED_BOOT_WM989_C2
        uc = pcs_["W1"].index.to_numpy()
        rng = np.random.default_rng(seed)
        boot_idx = [rng.integers(0, len(uc), len(uc)) for _ in range(N_BOOTSTRAP)]

        def compare(a, b, per=pcs_, bidx=boot_idx, pl=pooled, sd=seed, ucv=uc):
            va, vb = per[a].to_numpy(), per[b].to_numpy()
            d = np.array([va[i].mean() - vb[i].mean() for i in bidx])
            lo95, hi95 = np.percentile(d, [2.5, 97.5])
            lo975, hi975 = np.percentile(d, [1.25, 98.75])
            return {"comparison": f"{a} - {b}", "point": float(pl[a] - pl[b]),
                    "ci95": [float(lo95), float(hi95)],
                    "ci975_two_sided": [float(lo975), float(hi975)],
                    "fraction_delta_le_0": float((d <= 0).mean()),
                    "replicates": N_BOOTSTRAP, "seed": sd, "bootstrap_unit": "clone",
                    "clones_resampled": int(len(ucv))}

        inference = {"interaction_W4_minus_W5": compare("W4", "W5"),
                     "full_state_W1_minus_W5": compare("W1", "W5")}

        # ---- treatment-level directional diagnostics ---------------------------------------- #
        by_treatment = {}
        for t in TREATMENT_ORDER:
            m = tkey == t
            if endpoint == "C1":
                def ll(p, mm=m, yy=yv):
                    return float(log_loss(yy[mm], np.clip(p[mm], 1e-15, 1 - 1e-15),
                                          labels=[0, 1]))
                w4v, w5v = ll(base["pred_W4"].to_numpy()), ll(w5)
            else:
                w4v = float(np.abs(yv[m] - base["pred_W4"].to_numpy()[m]).mean())
                w5v = float(np.abs(yv[m] - w5[m]).mean())
            by_treatment[t] = {"W4": w4v, "W5": w5v, "improvement_W4_minus_W5": w4v - w5v,
                               "improved": bool(w4v - w5v > 0), "rows": int(m.sum())}
        n_improved = sum(v["improved"] for v in by_treatment.values())

        frame = base.copy()
        frame["pred_W5"] = w5
        frame.to_csv(WM989_W5_OOF if endpoint == "C1" else
                     _RESULTS / "stage23_wm989_interaction_abundance_oof.csv",
                     index=False, lineterminator="\n")
        out[endpoint] = {"rows": int(len(base)), "clones": int(len(uc)),
                         "selected_hyperparameters_per_outer_fold": selected,
                         "per_outer_fold": meta,
                         "pooled_metric": "log_loss" if endpoint == "C1" else "clone_balanced_MAE",
                         "pooled": pooled, "inference": inference,
                         "by_treatment": by_treatment,
                         "treatments_improved_by_W5_over_W4": n_improved}

    # ---- verdict, derived from V2 §6.5 -------------------------------------------------------- #
    def fam(ep):
        i = out[ep]["inference"]["interaction_W4_minus_W5"]
        s = out[ep]["inference"]["full_state_W1_minus_W5"]
        return {"pass_int": i["ci975_two_sided"][0] > 0, "pass_full": s["ci975_two_sided"][0] > 0,
                "harm_int": i["ci975_two_sided"][1] < 0, "harm_full": s["ci975_two_sided"][1] < 0,
                "point_int": i["point"], "point_full": s["point"],
                "n_treat": out[ep]["treatments_improved_by_W5_over_W4"]}

    fams = {ep: fam(ep) for ep in ("C1", "C2")}
    verdict = INTERACTION_NONE
    passing_endpoint = None
    for ep, other in (("C1", "C2"), ("C2", "C1")):
        a, b = fams[ep], fams[other]
        if (a["pass_int"] and a["pass_full"] and a["n_treat"] >= 3
                and not b["harm_int"] and not b["harm_full"]):
            verdict, passing_endpoint = INTERACTION_PASS, ep
            break
    if verdict != INTERACTION_PASS:
        any_point = any(f["point_int"] > 0 or f["point_full"] > 0 for f in fams.values())
        any_treat = any(f["n_treat"] >= 1 for f in fams.values())
        if any_point and any_treat:
            verdict = INTERACTION_LOCAL

    res = {"stage": "23D", "dataset": "GSE279162", "role": "B",
           "plan": {"file": PLAN.name, "version": PLAN_VERSION,
                    "canonical_lf_sha256": canonical_text_sha256(PLAN)},
           "protocol_sha256": sha256_file(_RESULTS / "stage23_protocol.json"),
           "models": {"W1": "B + U (frozen from 23C)", "W4": "X + B + U (frozen from 23C)",
                      "W5": "X + B + U + X*U"},
           "interaction_terms": "standardized PCA(X) score x non-reference treatment dummy only; "
                                "no gene-level interaction",
           "w1_w4_reused_from_23c": True,
           "endpoints": out, "endpoint_families": fams,
           "passing_endpoint": passing_endpoint,
           "convergence_warnings": convergence,
           "verdict": verdict,
           "verdict_is_provisional_until": "23E structural controls + "
                                           "ROLE_B_INTERACTION_PERMUTATION_PASS"}
    WM989_INTERACTION_RESULTS.write_text(json.dumps(res, indent=2), encoding="utf-8",
                                         newline="\n")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 23 learnability gate")
    ap.add_argument("--stage", default="23a", choices=["23a", "23b", "23c", "23d"])
    ap.add_argument("--rewind-root", type=Path, default=S21D.REWIND)
    ap.add_argument("--wm989-root", type=Path, default=S21D.WM989)
    args = ap.parse_args(argv)
    _RESULTS.mkdir(exist_ok=True)
    _CACHE.mkdir(parents=True, exist_ok=True)

    if args.stage == "23d":
        r = run_23d(args.wm989_root)
        for ep in ("C1", "C2"):
            e = r["endpoints"][ep]
            print(f"  {ep} {e['pooled_metric']}: " + "  ".join(
                f"{m}={e['pooled'][m]:.5f}" for m in ("W1", "W4", "W5")))
            for key in ("interaction_W4_minus_W5", "full_state_W1_minus_W5"):
                d = e["inference"][key]
                print(f"      {d['comparison']:<9} {d['point']:+.5f}  "
                      f"97.5% [{d['ci975_two_sided'][0]:+.5f}, {d['ci975_two_sided'][1]:+.5f}]")
            print(f"      treatments improved by W5 over W4: {e['treatments_improved_by_W5_over_W4']}/6")
        print("OVERALL:", r["verdict"])
        return 0

    if args.stage == "23c":
        r = run_23c(args.wm989_root)
        for ep in ("C1", "C2"):
            e = r["endpoints"][ep]
            key = list(e["pooled_oof_metrics"]["W0"])[0]
            print(f"  {ep} {key}: " + "  ".join(
                f"{m}={e['pooled_oof_metrics'][m][key]:.5f}" for m in ("W0","W1","W2","W3","W4")))
            d = e["inference"]["delta_state_W1_minus_W4"]
            print(f"      delta(W1-W4) = {d['point']:+.5f}  95% {d['ci95']}  97.5% {d['ci975_two_sided']}")
        print("OVERALL:", r["verdict"])
        return 0

    if args.stage == "23b":
        r = run_23b(args.rewind_root)
        for m in ("R0", "R1", "R2", "R3"):
            print(f"  {m}  AP={r['pooled_oof_metrics'][m]['AP']:.4f}")
        st = r["inference"]["delta_AP_state_R3_minus_R1"]
        print(f"  dAP(R3-R1) = {st['point']:+.4f}  "
              f"95% CI [{st['ci95_low']:+.4f}, {st['ci95_high']:+.4f}]")
        print("OVERALL:", r["provisional_verdict"])
        return 0

    # ---- step 1: the audit gates everything ------------------------------------------------- #
    audit = audit_stage22_inputs(args.rewind_root, args.wm989_root)
    print(f"input audit: {audit['verdict']}  ({audit['n_checks']} checks)")
    if audit["verdict"] != AUDITED:
        for f in audit["failed_checks"]:
            print("  FAILED:", f)
        (_RESULTS / "stage23_protocol.json").write_text(
            json.dumps({"stage": "23A", "verdict": PROTOCOL_BLOCKED, "input_audit": audit},
                       indent=2, default=str), encoding="utf-8", newline="\n")
        print("OVERALL:", PROTOCOL_BLOCKED)
        return 1

    # ---- step 2: pseudobulk ------------------------------------------------------------------ #
    rc = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")
    wn = pd.read_csv(_RESULTS / "stage22_wm989_naive_cells.csv")
    feat_file = next(args.wm989_root.rglob("*_features.tsv.gz"))
    feature_ids = [fid for fid, _sym in gene_features(feat_file)]
    if len(feature_ids) != N_GENES:
        raise RuntimeError(f"expected {N_GENES} gene features, got {len(feature_ids)}")

    manifests = {}
    for name, cells, roots, out in (
            ("GSE227151", rc, [args.rewind_root], "stage23_rewind_clone_expression_manifest.json"),
            ("GSE279162", wn, [args.wm989_root], "stage23_wm989_clone_expression_manifest.json")):
        X, clones = clone_pseudobulk(cells, roots)
        sparse.save_npz(_CACHE / f"{name}_pseudobulk.npz", X)
        (_CACHE / f"{name}_clones.json").write_text(json.dumps(clones), encoding="utf-8")
        m = expression_manifest(name, X, clones, cells, feature_ids)
        manifests[name] = (m, X, clones)
        (_RESULTS / out).write_text(json.dumps(m, indent=2), encoding="utf-8", newline="\n")
        print(f"{name}: pseudobulk {X.shape[0]} clones x {X.shape[1]} genes, nnz={X.nnz}")

    # ---- step 3: descriptive outer-training gene filter --------------------------------------- #
    prep = {}
    for name, clone_table, in (("GSE227151", "stage22_rewind_clones.csv"),
                               ("GSE279162", "stage22_wm989_clones.csv")):
        m, X, clones = manifests[name]
        tbl = pd.read_csv(_RESULTS / clone_table)
        fold = tbl.set_index("clone_id")["outer_fold"].to_dict()
        folds_arr = np.array([fold[c] for c in clones])
        per = {}
        for f in range(N_OUTER):
            train = np.flatnonzero(folds_arr != f)
            keep = training_fold_gene_filter(X, train)
            per[str(f)] = {"outer_training_clones": int(len(train)),
                           "detection_floor": max(5, int(np.ceil(0.01 * len(train)))),
                           "retained_genes": int(len(keep)),
                           "max_feasible_K": int(min(len(train), len(keep))),
                           "all_K_feasible": bool(min(len(train), len(keep)) >= max(K_CANDIDATES))}
        prep[name] = per
        print(f"{name}: retained genes per outer fold = "
              f"{[per[str(f)]['retained_genes'] for f in range(N_OUTER)]}")
    (_RESULTS / "stage23_outer_fold_preprocessing.json").write_text(
        json.dumps({"note": "descriptive only; the estimator pipeline refits this filter inside "
                            "every inner split per V2 SS2.5", "per_dataset": prep},
                   indent=2), encoding="utf-8", newline="\n")

    # ---- step 4: freeze the protocol ---------------------------------------------------------- #
    proto = build_protocol(STAGE23_TIGHTENED)
    proto["input_audit"] = audit
    proto["expression_manifests"] = {k: v[0] for k, v in manifests.items()}
    all_k = all(v["all_K_feasible"] for d in prep.values() for v in d.values())
    proto["verdict"] = PROTOCOL_FROZEN if all_k else PROTOCOL_BLOCKED
    (_RESULTS / "stage23_protocol.json").write_text(
        json.dumps(proto, indent=2, default=str), encoding="utf-8", newline="\n")
    print("OVERALL:", proto["verdict"])
    return 0 if proto["verdict"] == PROTOCOL_FROZEN else 1


if __name__ == "__main__":
    sys.exit(main())
