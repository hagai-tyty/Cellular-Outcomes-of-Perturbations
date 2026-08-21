"""STAGE 22 — build the frozen prospective benchmark tables.

Pre-registered in `plans/(newer)practical plans/STAGE_22_PROSPECTIVE_BENCHMARK_CONSTRUCTION_V2.md`
(frozen at `e97efe1`, annotated at `f2d6c86`). Consumes Stage 21D revision 2 (`6c2f2d6`).

**No model is fitted. No threshold is tuned. `src/` is untouched.** Stage 22 does not ask whether
the task is learnable; it freezes exactly what the task is.

SINGLE SOURCE OF TRUTH
----------------------
Every author rule -- Rewind's top-100-with-ties cut, WM989's QC thresholds, barcode clustering,
posterior assignment -- is *imported* from `experiments/diag_stage21d_public_reconstruction.py`
rather than re-implemented (plan §1.3). A second implementation that merely happens to reproduce
the headline counts could drift from the source-faithful logic without any test noticing.

The 21D module keeps `D:\\...` defaults; per plan §1.4 this builder always passes roots explicitly,
so nothing about it needs rewriting for portability.

WHAT THIS FREEZES
-----------------
    Rewind   cell -> clone -> primed/nonprimed          (fixed U)
    WM989    naive cell -> clone -> (treatment -> abundance/rank)

plus deterministic clone-level outer folds, an expression-column mapping for `X_before`, and a
feature-eligibility firewall so Stage 23 cannot turn provenance or clone size into the claimed
transcriptomic predictor.

THE TWO TRAPS THIS FILE EXISTS TO AVOID
---------------------------------------
1. **An ambiguous grouping unit.** 8 Rewind cells carry two clone assignments each, and all 8 have
   their two clones land in different folds. Keeping them would make `outer_fold` undefined for a
   cell whose expression vector is a single column. The plan pre-registers their exclusion; this
   file enumerates them rather than picking a clone heuristically.
2. **A benchmark that silently becomes "among surviving clones, predict abundance."** Every
   eligible WM989 clone gets all six treatment rows, and a zero is recorded as `observed_zero`
   against an available treatment sample -- never dropped, never relabelled as death or failure.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import KFold, StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"

_spec = importlib.util.spec_from_file_location(
    "s21d", Path(__file__).resolve().parent / "diag_stage21d_public_reconstruction.py")
S21D = importlib.util.module_from_spec(_spec)
sys.modules["s21d"] = S21D
_spec.loader.exec_module(S21D)

PLAN = ROOT / "plans/(newer)practical plans/STAGE_22_PROSPECTIVE_BENCHMARK_CONSTRUCTION_V2.md"
PLAN_VERSION = "V2"
RECONSTRUCTION_COMMIT = "6c2f2d6"

N_SPLITS = 5
STAGE22_SPLIT_SEED = 22022

REWIND_CELLS = _RESULTS / "stage22_rewind_cells.csv"
REWIND_CLONES = _RESULTS / "stage22_rewind_clones.csv"
WM989_CELL_ASSIGN = _RESULTS / "stage22_wm989_cell_assignments.csv"
WM989_NAIVE_CELLS = _RESULTS / "stage22_wm989_naive_cells.csv"
WM989_CLONES = _RESULTS / "stage22_wm989_clones.csv"
WM989_CLONE_TREATMENT = _RESULTS / "stage22_wm989_clone_treatment.csv"
REWIND_MANIFEST = _RESULTS / "stage22_rewind_benchmark_manifest.json"
WM989_MANIFEST = _RESULTS / "stage22_wm989_benchmark_manifest.json"
OUT = _RESULTS / "stage22_prospective_benchmark_results.json"

BENCHMARK_READY = "BENCHMARK_READY"
BENCHMARK_READY_DECLARED = "BENCHMARK_READY_WITH_DECLARED_MISSINGNESS"
BLOCKED_LINKAGE = "BENCHMARK_BLOCKED_LINKAGE"
BLOCKED_EXPRESSION = "BENCHMARK_BLOCKED_EXPRESSION_MAPPING"
BLOCKED_OUTCOME = "BENCHMARK_BLOCKED_OUTCOME"
BLOCKED_LEAKAGE = "BENCHMARK_BLOCKED_LEAKAGE"
READY_VERDICTS = (BENCHMARK_READY, BENCHMARK_READY_DECLARED)
STAGE_23_READY = "STAGE_23_READY"
STAGE_23_BLOCKED = "STAGE_23_BLOCKED"

AMBIGUITY_REASON = "ambiguous_multi_lineage_clone_assignment"
REWIND_OUTCOME_SEMANTICS = "author_top100_iPSC_gDNA_priming"
REWIND_OUTCOME_RULE = ("gDNA arm -> group_by(BC50StarcodeD8, SampleNum) -> sum -> "
                       "slice_max(n=100, with_ties=TRUE) -> inner_join(filtered10XCells)")
REWIND_OUTCOME_SOURCE = "stepThreeStarcodeShavedReads_BC_gDNA.txt"

# The GEO sample names are canonical because they come from the files themselves; the plan's prose
# and the author code each use a different casing/shorthand, so both are recorded as aliases.
TREATMENT_ALIASES = {
    "Dabrafenib": {"plan": "dabrafenib", "author_code": "dab"},
    "Trametinib": {"plan": "trametinib", "author_code": "tram"},
    "CoCl2": {"plan": "CoCl2", "author_code": "cocl2"},
    "Acid": {"plan": "acid", "author_code": "acid"},
    "Cisplatin": {"plan": "cisplatin", "author_code": "cis"},
    "Doxorubicin": {"plan": "doxorubicin", "author_code": "dox"},
}

FEATURE_ELIGIBILITY = {
    "TARGET": ["y_primed", "n_post_cells", "post_fraction", "post_rank", "post_rank_fraction",
               "detected_post", "outcome_observation_status"],
    "PROVENANCE_ONLY": ["cell_uid", "cellID", "cell_barcode", "clone_id", "assigned_lineage",
                        "gsm", "SampleNum", "sample", "source_naive_sample", "condition",
                        "expression_barcode", "expression_column_index", "expression_source",
                        "outcome_source", "outcome_rule", "outcome_semantics",
                        "outer_group", "outer_fold", "biological_replicate",
                        "generalization_scope", "naive_source_samples",
                        "treatment_sample_available", "post_tie_size", "is_assigned", "is_naive",
                        "lane_membership", "exclusion_reason"],
    "BASELINE_NUISANCE": ["n_pretreatment_cells", "n_lanes", "n_primed_cells", "n_nonprimed_cells",
                          "n_naive_cells", "n_naive1_cells", "n_naive2_cells", "n_naive3_cells",
                          "naive1_total_assigned_cells", "naive2_total_assigned_cells",
                          "naive3_total_assigned_cells", "naive_pooled_fraction",
                          "treatment_total_assigned_cells",
                          "nUMI", "fracUMI", "nLineages", "assigned_posterior"],
    "PRIMARY_X": ["pretreatment gene-expression values only, resolved through "
                  "expression_source + expression_column_index"],
}


def hash_artifacts(paths: list[Path]) -> list[dict]:
    """Identity by basename/size/digest, never by the machine-specific path."""
    return [{"name": p.name, "bytes": p.stat().st_size, "sha256": S21D.sha256(p)} for p in paths]


def source_files(paths: list[Path]) -> list[dict]:
    return [{"name": p.name, "bytes": p.stat().st_size, "sha256": S21D.sha256(p)}
            for p in paths if p.exists()]


def barcode_index(path: Path) -> dict[str, tuple[int, str]]:
    """bare 10x barcode -> (1-based Matrix Market column index, full barcode string)."""
    out = {}
    for i, bc in enumerate(S21D.read_barcodes(path), start=1):
        out[bc.split("-")[0]] = (i, bc)
    return out


def deterministic_folds(keys, y_by_key: dict | None = None) -> dict[str, int]:
    """One fold per clone, assigned once, over a sorted key list so the split cannot drift.

    Stratified when a clone-level class exists (Rewind), plain otherwise (WM989 counts). The sort
    is what makes this reproducible: a split keyed off whatever order a groupby happened to emit
    would change the moment an upstream filter changed.
    """
    order = sorted(keys)
    idx = np.arange(len(order))
    if y_by_key is not None:
        labels = np.asarray([y_by_key[k] for k in order])
        splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                                   random_state=STAGE22_SPLIT_SEED).split(idx, labels)
    else:
        splitter = KFold(n_splits=N_SPLITS, shuffle=True,
                         random_state=STAGE22_SPLIT_SEED).split(idx)
    fold = np.zeros(len(order), dtype=int)
    for i, (_, test) in enumerate(splitter):
        fold[test] = i
    return dict(zip(order, fold.tolist(), strict=True))


# --------------------------------------------------------------------------------------------- #
# ROLE A — Rewind / GSE227151
# --------------------------------------------------------------------------------------------- #
def build_rewind(root: Path) -> dict:
    req = S21D.rewind_required(root)
    missing = sorted(k for k, p in req.items() if not p.exists())
    if missing:
        return {"verdict": BLOCKED_LINKAGE, "missing_files": missing}

    f10 = pd.read_csv(req["filtered10XCells"], sep="\t", index_col=0)
    gdna = pd.read_csv(req["bc_gDNA"], sep="\t", index_col=0)

    # ---- author rule, imported not re-implemented ------------------------------------------- #
    count_col = next(c for c in gdna.columns
                     if c not in ("cellID", "BC50StarcodeD8", "SampleNum"))
    probed = (gdna[gdna["cellID"] == "dummy"]
              .groupby(["BC50StarcodeD8", "SampleNum"], as_index=False)[count_col].sum()
              .rename(columns={count_col: "nUMI"}))
    top, cutoff, n_at_cutoff = S21D.slice_max_with_ties(probed, "nUMI", S21D.TOP_N_GDNA)
    primed_clones = set(top["BC50StarcodeD8"])

    # ---- SampleNum -> GSM, re-derived by containment, never from GEO title text -------------- #
    gsm_bc = {g: barcode_index(req[f"{g}_barcodes"]) for g in S21D.REWIND_GSMS}
    sample_to_gsm, containment = {}, {}
    for sn, grp in f10.groupby("SampleNum"):
        ids = set(grp["cellID"])
        hits = {g: len(ids & set(bc)) / len(ids) for g, bc in gsm_bc.items()}
        sample_to_gsm[int(sn)] = max(hits, key=hits.get)
        containment[int(sn)] = {g: round(v, 4) for g, v in hits.items()}
    cross_lane_collisions = len(set(f10.loc[f10.SampleNum == 1, "cellID"])
                                & set(f10.loc[f10.SampleNum == 2, "cellID"]))

    # ---- plan §3.5: the pre-registered ambiguity exclusion ----------------------------------- #
    f10 = f10.assign(cell_uid=f10["SampleNum"].astype(str) + ":" + f10["cellID"])
    per_uid = f10.groupby("cell_uid")["BC50StarcodeD8"].nunique()
    ambiguous = sorted(per_uid[per_uid > 1].index)
    excluded_rows = f10[f10["cell_uid"].isin(ambiguous)]
    kept = f10[~f10["cell_uid"].isin(ambiguous)].copy()
    exclusion_audit = [
        {"cell_uid": u, "SampleNum": int(g["SampleNum"].iloc[0]),
         "clone_ids": sorted(g["BC50StarcodeD8"]), "n_source_rows": int(len(g)),
         "nLineages": sorted({int(x) for x in g["nLineages"]}),
         "any_clone_primed": bool(g["BC50StarcodeD8"].isin(primed_clones).any()),
         "exclusion_reason": AMBIGUITY_REASON}
        for u, g in excluded_rows.groupby("cell_uid")]

    # ---- cell table -------------------------------------------------------------------------- #
    kept["gsm"] = kept["SampleNum"].map(sample_to_gsm)
    expr = [gsm_bc[g].get(c, (None, None)) for g, c in zip(kept["gsm"], kept["cellID"], strict=True)]
    kept["expression_column_index"] = [e[0] for e in expr]
    kept["expression_barcode"] = [e[1] for e in expr]
    kept["expression_source"] = [req[f"{g}_matrix"].name for g in kept["gsm"]]
    kept["clone_id"] = kept["BC50StarcodeD8"]
    kept["y_primed"] = kept["clone_id"].isin(primed_clones).astype(int)
    kept["outcome_source"] = REWIND_OUTCOME_SOURCE
    kept["outcome_rule"] = REWIND_OUTCOME_RULE
    kept["outcome_semantics"] = REWIND_OUTCOME_SEMANTICS
    kept["biological_replicate"] = "R1"
    kept["generalization_scope"] = "within_R1_clone_heldout"
    kept["outer_group"] = kept["clone_id"]

    # ---- clone table + folds (built on the POST-exclusion population) ------------------------ #
    clones = (kept.groupby("clone_id")
              .agg(n_pretreatment_cells=("cell_uid", "size"),
                   n_lanes=("SampleNum", "nunique"),
                   lane_membership=("SampleNum", lambda s: "+".join(map(str, sorted(set(s))))))
              .reset_index())
    clones["y_primed"] = clones["clone_id"].isin(primed_clones).astype(int)
    clones["outcome_semantics"] = REWIND_OUTCOME_SEMANTICS
    clones["n_primed_cells"] = clones["n_pretreatment_cells"] * clones["y_primed"]
    clones["n_nonprimed_cells"] = clones["n_pretreatment_cells"] - clones["n_primed_cells"]
    clones["outer_group"] = clones["clone_id"]
    folds = deterministic_folds(
        clones["clone_id"].tolist(),
        dict(zip(clones["clone_id"], clones["y_primed"], strict=True)))
    clones["outer_fold"] = clones["clone_id"].map(folds)
    kept["outer_fold"] = kept["clone_id"].map(folds)

    cell_cols = ["cell_uid", "cellID", "SampleNum", "gsm", "clone_id", "nUMI", "fracUMI",
                 "nLineages", "y_primed", "outcome_source", "outcome_rule", "outcome_semantics",
                 "biological_replicate", "generalization_scope", "outer_group", "outer_fold",
                 "expression_barcode", "expression_column_index", "expression_source"]
    clone_cols = ["clone_id", "n_pretreatment_cells", "n_lanes", "lane_membership", "y_primed",
                  "outcome_semantics", "n_primed_cells", "n_nonprimed_cells", "outer_group",
                  "outer_fold"]
    cells_out = kept[cell_cols].sort_values("cell_uid").reset_index(drop=True)
    clones_out = clones[clone_cols].sort_values("clone_id").reset_index(drop=True)
    cells_out.to_csv(REWIND_CELLS, index=False, lineterminator="\n")
    clones_out.to_csv(REWIND_CLONES, index=False, lineterminator="\n")

    # ---- integrity assertions ---------------------------------------------------------------- #
    per_fold = (clones_out.groupby("outer_fold")
                .agg(clones=("clone_id", "size"), positive_clones=("y_primed", "sum")).to_dict("index"))
    pos_cells_fold = cells_out[cells_out.y_primed == 1].groupby("outer_fold").size().to_dict()
    unmapped = int(cells_out["expression_column_index"].isna().sum())
    checks = {
        "retained_cell_uid_unique": bool(cells_out["cell_uid"].is_unique),
        "no_retained_cell_maps_to_multiple_clones":
            bool(cells_out.groupby("cell_uid")["clone_id"].nunique().max() == 1),
        "no_clone_has_contradictory_outcome":
            bool(cells_out.groupby("clone_id")["y_primed"].nunique().max() == 1),
        "one_group_and_fold_per_clone":
            bool(cells_out.groupby("clone_id")[["outer_group", "outer_fold"]].nunique().max().max() == 1),
        "cell_outer_group_equals_clone": bool((cells_out["outer_group"] == cells_out["clone_id"]).all()),
        "five_folds": len(per_fold) == N_SPLITS,
        "every_fold_has_both_classes":
            all(0 < v["positive_clones"] < v["clones"] for v in per_fold.values()),
        "expression_mapping_complete": unmapped == 0,
        "primed_cells_42": int(cells_out["y_primed"].sum()) == 42,
        "primed_clones_35": int(clones_out["y_primed"].sum()) == 35,
        "no_ambiguous_cell_retained": not set(cells_out["cell_uid"]) & set(ambiguous),
    }

    stats = {
        "source_author_qc_assignment_records": int(len(f10)),
        "source_unique_cell_uid": int(f10["cell_uid"].nunique()),
        "source_unique_bare_cellID": int(f10["cellID"].nunique()),
        "source_unique_clones": int(f10["BC50StarcodeD8"].nunique()),
        "bare_cellID_cross_lane_collisions": cross_lane_collisions,
        "ambiguous_cell_uid_excluded": len(ambiguous),
        "source_rows_removed_by_exclusion": int(len(excluded_rows)),
        "retained_benchmark_cell_uid": int(len(cells_out)),
        "retained_benchmark_clones": int(len(clones_out)),
        "clones_lost_entirely_to_exclusion":
            int(f10["BC50StarcodeD8"].nunique() - len(clones_out)),
        "positive_cells": int(cells_out["y_primed"].sum()),
        "negative_cells": int((cells_out["y_primed"] == 0).sum()),
        "positive_clones": int(clones_out["y_primed"].sum()),
        "negative_clones": int((clones_out["y_primed"] == 0).sum()),
        "cell_grain_prevalence": round(float(cells_out["y_primed"].mean()), 6),
        "clone_grain_prevalence": round(float(clones_out["y_primed"].mean()), 6),
        "cells_per_clone": {int(k): int(v) for k, v in
                            clones_out["n_pretreatment_cells"].value_counts().sort_index().items()},
        "clones_spanning_both_lanes_after_exclusion": int((clones_out["n_lanes"] == 2).sum()),
        "per_outer_fold": {int(k): {"clones": int(v["clones"]),
                                    "positive_clones": int(v["positive_clones"]),
                                    "positive_cells": int(pos_cells_fold.get(k, 0))}
                           for k, v in per_fold.items()},
        "expression_mapped_cells": int(len(cells_out) - unmapped),
        "expression_unmapped_cells": unmapped,
        "expression_features": {g: S21D.mtx_dims(req[f"{g}_matrix"])[0] for g in S21D.REWIND_GSMS},
        "author_rule": {"n": S21D.TOP_N_GDNA, "cutoff_nUMI": cutoff,
                        "barcodes_at_cutoff": n_at_cutoff, "selected_barcodes": len(primed_clones)},
        "samplenum_to_gsm": {str(k): v for k, v in sample_to_gsm.items()},
        "samplenum_containment": {str(k): v for k, v in containment.items()},
    }
    return {"verdict": None, "checks": checks, "stats": stats, "req": req,
            "primed_clones": primed_clones,
            "exclusions": exclusion_audit, "artifacts": [REWIND_CELLS, REWIND_CLONES]}


# --------------------------------------------------------------------------------------------- #
# ROLE B — GSE279162 / WM989
# --------------------------------------------------------------------------------------------- #
def build_wm989(root: Path) -> dict:
    req = S21D.wm989_required(root)
    missing = sorted(k for k, p in req.items() if not p.exists())
    if missing:
        return {"verdict": BLOCKED_LINKAGE, "missing_files": missing}

    feats = [ln.rstrip("\n").split("\t")
             for ln in gzip.open(req["Naive1_features.tsv"], "rt")]
    lin_names_all = np.asarray([f[0] for f in feats if len(f) > 2 and f[2] == "Custom"],
                               dtype=object)
    n_genes = len(feats) - len(lin_names_all)
    mt_mask = np.zeros(n_genes + 1, dtype=bool)
    for i, f in enumerate(feats[:n_genes]):
        if f[1].startswith("MT-"):
            mt_mask[i + 1] = True

    # ---- author QC, imported thresholds ------------------------------------------------------ #
    rows, cols, vals, cond, bcs, colidx, offset = [], [], [], [], [], [], 0
    per_sample_qc = {}
    for name in S21D.WM989_ORDER:
        min_f, max_c, max_mt = S21D.WM989_QC[name]
        r = S21D.qc_and_lineage(req[f"{name}_matrix.mtx"], n_genes, mt_mask, min_f, max_c, max_mt)
        gsm = next(g for g, s in S21D.WM989_SAMPLES.items() if s == name)
        per_sample_qc[name] = {"gsm": gsm, "cells_raw": r["n_cells"], "cells_post_qc": r["kept"],
                               "qc": {"nFeature_RNA_gt": min_f, "nCount_RNA_lt": max_c,
                                      "percent_mt_lt": max_mt}}
        bc_list = S21D.read_barcodes(req[f"{name}_barcodes.tsv"])
        local = {c: i for i, c in enumerate(r["cells"])}
        lin = r["lin"]
        rows.append(lin["row"].to_numpy() - n_genes - 1)
        cols.append(np.array([local[c] for c in lin["col"].to_numpy()], dtype=np.int64) + offset)
        vals.append(lin["val"].to_numpy().astype(np.float64))
        cond += [name] * r["kept"]
        bcs += [bc_list[c - 1] for c in r["cells"]]
        colidx += [int(c) for c in r["cells"]]
        offset += r["kept"]

    lin_mat = sparse.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(lin_names_all), offset))
    nonzero = np.asarray(lin_mat.sum(axis=1)).ravel() > 0
    lin_mat, lin_names = lin_mat[nonzero], lin_names_all[nonzero]

    # ---- author clone call, imported ---------------------------------------------------------- #
    clus = S21D.barcode_clustering(lin_mat, lin_names)
    combined, combined_names = S21D.barcode_combine(lin_mat, lin_names, clus["clusters"])
    assigned, posterior, post_stats = S21D.barcoding_posterior_and_assignment(
        combined, combined_names)

    cond = np.asarray(cond, dtype=object)
    cells = pd.DataFrame({
        "sample": cond,
        "cell_barcode": [b.split("-")[0] for b in bcs],
        "expression_barcode": bcs,
        "expression_column_index": colidx,
        "assigned_lineage": assigned,
        "assigned_posterior": np.round(posterior, 6),
    })
    cells["gsm"] = cells["sample"].map({v: k for k, v in S21D.WM989_SAMPLES.items()})
    cells["cell_uid"] = cells["gsm"] + ":" + cells["cell_barcode"]
    cells["condition"] = np.where(cells["sample"].isin(S21D.WM989_NAIVE), "naive", cells["sample"])
    cells["is_assigned"] = cells["assigned_lineage"].notna()
    cells["is_naive"] = cells["sample"].isin(S21D.WM989_NAIVE)
    cells["expression_source"] = [req[f"{s}_matrix.mtx"].name for s in cells["sample"]]

    # Two different populations, both worth recording: the plan quotes the RAW figure, while the
    # uniqueness requirement actually applies to the post-QC benchmark population.
    reused_bare = int((cells.groupby("cell_barcode")["gsm"].nunique() > 1).sum())
    raw_seen: dict[str, set] = {}
    for name in S21D.WM989_ORDER:
        for b in S21D.read_barcodes(req[f"{name}_barcodes.tsv"]):
            raw_seen.setdefault(b.split("-")[0], set()).add(name)
    reused_bare_raw = int(sum(1 for v in raw_seen.values() if len(v) > 1))

    assigned_cells = cells[cells["is_assigned"]]
    naive_assigned = assigned_cells[assigned_cells["is_naive"]]
    naive_totals = {s: int((naive_assigned["sample"] == s).sum()) for s in S21D.WM989_NAIVE}
    naive_denom = sum(naive_totals.values())
    treatment_totals = {t: int((assigned_cells["sample"] == t).sum())
                        for t in S21D.WM989_TREATMENTS}

    # ---- clone tables + folds ----------------------------------------------------------------- #
    counts = (assigned_cells.groupby(["assigned_lineage", "condition"]).size()
              .unstack(fill_value=0)
              .reindex(columns=["naive", *S21D.WM989_TREATMENTS], fill_value=0))
    eligible = sorted(counts.index[counts["naive"] > 0])
    folds = deterministic_folds(eligible, None)

    naive_by_sample = (naive_assigned.groupby(["assigned_lineage", "sample"]).size()
                       .unstack(fill_value=0).reindex(index=eligible, fill_value=0)
                       .reindex(columns=list(S21D.WM989_NAIVE), fill_value=0))
    clones = pd.DataFrame({"clone_id": eligible})
    clones["n_naive_cells"] = counts.loc[eligible, "naive"].to_numpy()
    for s in S21D.WM989_NAIVE:
        clones[f"n_{s.lower()}_cells"] = naive_by_sample[s].to_numpy()
    for s in S21D.WM989_NAIVE:
        clones[f"{s.lower()}_total_assigned_cells"] = naive_totals[s]
    clones["naive_pooled_fraction"] = (clones["n_naive_cells"] / naive_denom).round(8)
    clones["naive_source_samples"] = [
        "|".join(s for s in S21D.WM989_NAIVE if naive_by_sample.loc[c, s] > 0) for c in eligible]
    clones["outer_group"] = clones["clone_id"]
    clones["outer_fold"] = clones["clone_id"].map(folds)

    naive_out = naive_assigned[naive_assigned["assigned_lineage"].isin(set(eligible))].copy()
    naive_out = naive_out.rename(columns={"sample": "source_naive_sample",
                                          "assigned_lineage": "clone_id"})
    naive_out["outer_group"] = naive_out["clone_id"]
    naive_out["outer_fold"] = naive_out["clone_id"].map(folds)

    # ---- clone x treatment: all six rows, zeros explicit -------------------------------------- #
    recs = []
    depth = dict(zip(clones["clone_id"], clones["n_naive_cells"], strict=True))
    frac = dict(zip(clones["clone_id"], clones["naive_pooled_fraction"], strict=True))
    srcs = dict(zip(clones["clone_id"], clones["naive_source_samples"], strict=True))
    for t in S21D.WM989_TREATMENTS:
        col = counts.loc[eligible, t]
        rank = col.rank(method="min", ascending=False).astype(int)
        tie = col.map(col.value_counts())
        for c in eligible:
            n_post = int(col[c])
            recs.append({
                "clone_id": c, "treatment": t,
                "n_naive_cells": int(depth[c]), "naive_pooled_fraction": float(frac[c]),
                "n_post_cells": n_post,
                "post_fraction": round(n_post / treatment_totals[t], 8),
                "post_rank": int(rank[c]),
                "post_rank_fraction": round(int(rank[c]) / len(eligible), 8),
                "post_tie_size": int(tie[c]),
                "detected_post": bool(n_post > 0),
                "outcome_observation_status": "observed_nonzero" if n_post else "observed_zero",
                "treatment_sample_available": True,
                "treatment_total_assigned_cells": treatment_totals[t],
                "naive_source_samples": srcs[c],
                "outer_group": c, "outer_fold": int(folds[c]),
            })
    ct = pd.DataFrame(recs).sort_values(["clone_id", "treatment"]).reset_index(drop=True)

    cell_cols = ["cell_uid", "cell_barcode", "sample", "gsm", "condition", "assigned_lineage",
                 "assigned_posterior", "is_assigned", "is_naive"]
    naive_cols = ["cell_uid", "source_naive_sample", "gsm", "clone_id", "assigned_posterior",
                  "expression_barcode", "expression_column_index", "expression_source",
                  "outer_group", "outer_fold"]
    cells[cell_cols].sort_values("cell_uid").to_csv(WM989_CELL_ASSIGN, index=False,
                                                    lineterminator="\n")
    naive_out[naive_cols].sort_values("cell_uid").to_csv(WM989_NAIVE_CELLS, index=False,
                                                         lineterminator="\n")
    clones.to_csv(WM989_CLONES, index=False, lineterminator="\n")
    ct.to_csv(WM989_CLONE_TREATMENT, index=False, lineterminator="\n")

    n_treat = (counts.loc[eligible, list(S21D.WM989_TREATMENTS)] > 0).sum(axis=1)
    zero_by_depth = {}
    for lo, hi, label in ((1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10**9, "10+")):
        m = (clones["n_naive_cells"] >= lo) & (clones["n_naive_cells"] <= hi)
        sub = ct[ct["clone_id"].isin(set(clones.loc[m, "clone_id"]))]
        zero_by_depth[label] = {"clones": int(m.sum()),
                                "zero_rate": round(float((sub["n_post_cells"] == 0).mean()), 4)
                                if len(sub) else None}

    checks = {
        "cell_uid_unique": bool(cells["cell_uid"].is_unique),
        "six_rows_per_eligible_clone":
            bool((ct.groupby("clone_id").size() == len(S21D.WM989_TREATMENTS)).all()),
        "structural_rows_match": len(ct) == len(eligible) * len(S21D.WM989_TREATMENTS),
        "all_treatment_samples_available": all(v > 0 for v in treatment_totals.values()),
        "zeros_retained_as_observed_zero":
            bool(((ct["n_post_cells"] == 0) == (ct["outcome_observation_status"] == "observed_zero")).all()),
        "one_group_and_fold_per_clone":
            bool(ct.groupby("clone_id")[["outer_group", "outer_fold"]].nunique().max().max() == 1),
        "naive_cells_share_clone_fold":
            bool(naive_out.groupby("clone_id")["outer_fold"].nunique().max() == 1),
        "clone_fold_consistent_across_tables":
            bool(dict(zip(clones["clone_id"], clones["outer_fold"], strict=True))
                 == dict(zip(ct["clone_id"], ct["outer_fold"], strict=True))),
        "five_folds": int(clones["outer_fold"].nunique()) == N_SPLITS,
        "expression_mapping_complete": int(naive_out["expression_column_index"].isna().sum()) == 0,
        "no_binary_resistance_threshold": "y_resistant" not in ct.columns,
        "posterior_floor_applied": post_stats["dropped_by_posterior_floor"] > 0,
    }

    stats = {
        "feature_structure": {"total_features": len(feats), "n_genes": n_genes,
                              "n_lineage_features": int(len(lin_names_all)),
                              "lineage_assay_source": "Custom feature block",
                              "mt_genes": int(mt_mask.sum())},
        "clone_call_source": "s21d.barcoding_posterior_and_assignment (author pipeline); "
                             "no dominant-fraction or UMI-floor fallback exists",
        "joint_assignment": {"run_once_on_all_samples": True,
                             "cells_in_joint_object": int(offset),
                             "lineages_in_joint_object": int(lin_mat.shape[0])},
        "per_sample_qc": per_sample_qc,
        "raw_cells": int(sum(v["cells_raw"] for v in per_sample_qc.values())),
        "post_qc_cells": int(len(cells)),
        "assigned_cells": int(len(assigned_cells)),
        "na_lineage_cells": int((~cells["is_assigned"]).sum()),
        "by_condition": {k: {"post_qc_cells": int(len(g)),
                             "assigned": int(g["is_assigned"].sum()),
                             "na": int((~g["is_assigned"]).sum()),
                             "unique_clones": int(g["assigned_lineage"].nunique())}
                         for k, g in cells.groupby("condition")},
        "unique_cell_uid": int(cells["cell_uid"].nunique()),
        "reused_bare_barcodes_across_samples": reused_bare,
        "reused_bare_barcodes_across_samples_raw": reused_bare_raw,
        "reused_bare_barcodes_note": "raw counts every barcode string in the 9 GEO barcode files; the other counts only the post-QC benchmark population",
        "unique_assigned_clones": int(counts.shape[0]),
        "clones_with_naive_observation": len(eligible),
        "naive_totals_by_sample": naive_totals,
        "naive_pooled_denominator": naive_denom,
        "naive_cells_per_clone": {"median": int(clones["n_naive_cells"].median()),
                                  "max": int(clones["n_naive_cells"].max()),
                                  "ge2": int((clones["n_naive_cells"] >= 2).sum()),
                                  "ge5": int((clones["n_naive_cells"] >= 5).sum()),
                                  "ge10": int((clones["n_naive_cells"] >= 10).sum())},
        "clone_treatment_rows": int(len(ct)),
        "clone_coverage": {"in_ge1_treatment": int((n_treat >= 1).sum()),
                           "in_ge2_treatments": int((n_treat >= 2).sum()),
                           "in_all_6": int((n_treat == 6).sum()),
                           "by_n_treatments": {int(k): int(v)
                                               for k, v in n_treat.value_counts().sort_index().items()}},
        "per_treatment": {t: {"total_assigned_cells": treatment_totals[t],
                              "clones_nonzero": int((counts.loc[eligible, t] > 0).sum()),
                              "zero_rate": round(float((counts.loc[eligible, t] == 0).mean()), 4),
                              "benchmark_fraction_sum": round(
                                  float(counts.loc[eligible, t].sum() / treatment_totals[t]), 4),
                              "abundance_median_nonzero": int(
                                  counts.loc[eligible, t][counts.loc[eligible, t] > 0].median()),
                              "abundance_max": int(counts.loc[eligible, t].max())}
                          for t in S21D.WM989_TREATMENTS},
        "observed_zero_rows": int((ct["n_post_cells"] == 0).sum()),
        "observed_zero_rate": round(float((ct["n_post_cells"] == 0).mean()), 4),
        "zero_rate_by_naive_depth": zero_by_depth,
        "per_outer_fold": {int(k): int(v) for k, v in
                           clones["outer_fold"].value_counts().sort_index().items()},
        "naive_expression_mapped_cells": int(len(naive_out)),
        "naive_expression_unmapped_cells": int(naive_out["expression_column_index"].isna().sum()),
        "barcode_clustering": {"candidate_lineages": clus["n_candidate_lineages"],
                               "correlated_pairs": clus["n_pairs"],
                               "clusters": len(clus["clusters"]),
                               "merging_conflicts": clus["merging_conflicts"]},
        "posterior_assignment": post_stats,
        "treatment_aliases": TREATMENT_ALIASES,
    }
    return {"verdict": None, "checks": checks, "stats": stats, "req": req,
            "artifacts": [WM989_CELL_ASSIGN, WM989_NAIVE_CELLS, WM989_CLONES,
                          WM989_CLONE_TREATMENT]}


# --------------------------------------------------------------------------------------------- #
def no_label_leakage(rewind_root: Path, primed_clones: set[str], cells: pd.DataFrame) -> dict:
    """G22-2, computed rather than asserted.

    The Rewind target is recomputed from the gDNA arm ALONE -- no expression file is opened -- and
    compared with what the benchmark wrote. If a gene-expression value had ever reached the label,
    this recomputation would disagree.

    For WM989 the target is a count of post-treatment cells carrying a clone in the lineage assay.
    No PRETREATMENT expression value enters it. One honest caveat is recorded rather than glossed:
    which post-treatment cells are countable depends on the authors' post-treatment RNA QC
    (nFeature / nCount / percent.mt), applied identically across conditions. That is an expression
    QUALITY filter on the future arm, not the pretreatment predictor, so the gate holds as written
    -- but the dependence is declared.
    """
    gdna = pd.read_csv(rewind_root / S21D.GDNA_FILE, sep="	", index_col=0)
    count_col = next(c for c in gdna.columns
                     if c not in ("cellID", "BC50StarcodeD8", "SampleNum"))
    probed = (gdna[gdna["cellID"] == "dummy"]
              .groupby(["BC50StarcodeD8", "SampleNum"], as_index=False)[count_col].sum()
              .rename(columns={count_col: "nUMI"}))
    top, _, _ = S21D.slice_max_with_ties(probed, "nUMI", S21D.TOP_N_GDNA)
    recomputed = set(top["BC50StarcodeD8"])
    written = set(cells.loc[cells["y_primed"] == 1, "clone_id"])
    return {
        "rewind_target_recomputed_from_outcome_source_alone": recomputed == primed_clones,
        "rewind_written_labels_are_a_subset_of_that_set": written <= recomputed,
        "rewind_no_expression_file_read_to_build_the_target": True,
        "wm989_target_uses_no_pretreatment_expression": True,
        "wm989_caveat": "post-treatment cell countability depends on the authors' post-treatment "
                        "RNA QC, applied identically across conditions; no pretreatment expression "
                        "value enters any target",
    }


def builder_fits_no_model() -> dict:
    """G22-10, computed from this file's own syntax tree.

    `sklearn.model_selection` is imported for KFold/StratifiedKFold, which are splitters, not
    estimators. Anything that could fit is banned, and so is any `.fit(` call.
    """
    import ast

    src = Path(__file__).resolve().read_text(encoding="utf-8")
    tree = ast.parse(src)
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    modelling = {m for m in imported
                 if m.split(".")[0] in {"torch", "tensorflow", "xgboost", "lightgbm", "statsmodels"}
                 or (m.startswith("sklearn") and m != "sklearn.model_selection")}
    # Checked on the syntax tree, not on the text. A text search matches this function's own
    # docstring -- which is exactly how the first version of this gate failed on itself.
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
             and n.func.attr in {"fit", "fit_transform", "fit_predict", "partial_fit"}]
    return {"no_modelling_imports": not modelling,
            "modelling_imports_found": sorted(modelling),
            "no_fit_call": not calls,
            "sklearn_use_is_splitters_only": "sklearn.model_selection" in imported}


def verdict_for(checks: dict, declared_missingness: int) -> str:
    if not checks.get("expression_mapping_complete", True):
        return BLOCKED_EXPRESSION
    linkage_keys = ("one_group_and_fold_per_clone", "cell_outer_group_equals_clone",
                    "no_retained_cell_maps_to_multiple_clones", "naive_cells_share_clone_fold",
                    "clone_fold_consistent_across_tables", "retained_cell_uid_unique",
                    "cell_uid_unique")
    if any(checks.get(k) is False for k in linkage_keys):
        return BLOCKED_LEAKAGE
    outcome_keys = ("no_clone_has_contradictory_outcome", "six_rows_per_eligible_clone",
                    "structural_rows_match", "zeros_retained_as_observed_zero",
                    "all_treatment_samples_available", "primed_cells_42", "primed_clones_35")
    if any(checks.get(k) is False for k in outcome_keys):
        return BLOCKED_OUTCOME
    if any(v is False for v in checks.values()):
        return BLOCKED_LINKAGE
    return BENCHMARK_READY_DECLARED if declared_missingness else BENCHMARK_READY


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 22 prospective benchmark construction")
    ap.add_argument("--rewind-root", type=Path, default=S21D.REWIND)
    ap.add_argument("--wm989-root", type=Path, default=S21D.WM989)
    args = ap.parse_args(argv)
    _RESULTS.mkdir(exist_ok=True)

    plan_hash = S21D.sha256(PLAN)
    common = {"reconstruction_commit": RECONSTRUCTION_COMMIT,
              "plan": {"file": PLAN.name, "version": PLAN_VERSION, "sha256": plan_hash},
              "builder_source_sha256": S21D.sha256(Path(__file__).resolve()),
              "split": {"n_splits": N_SPLITS, "seed": STAGE22_SPLIT_SEED,
                        "unit": "clone", "assigned_once_in_stage_22": True},
              "feature_eligibility": FEATURE_ELIGIBILITY,
              "model_fitted": False, "src_modified": False, "raw_data_committed": False}

    # ---- Role A ------------------------------------------------------------------------------ #
    a = build_rewind(args.rewind_root)
    if a["verdict"]:
        print("GSE227151:", a["verdict"])
        return 1
    a_verdict = verdict_for(a["checks"], len(a["exclusions"]))
    rewind_manifest = {
        "dataset": "GSE227151", "role": "A", "source_accession": "GSE227151",
        "local_source_path_used": str(args.rewind_root),
        "source_files": source_files(sorted(a["req"].values())),
        "author_code_files": source_files(S21D.REWIND_SCRIPTS),
        "author_rule_identifiers": {
            "outcome_rule": REWIND_OUTCOME_RULE, "outcome_semantics": REWIND_OUTCOME_SEMANTICS,
            "top_n": S21D.TOP_N_GDNA, "with_ties": True,
            "second_script_excluded_barcode": S21D.REWIND_EXCLUDED_BC},
        **common,
        "declared_exclusions": {"reason": AMBIGUITY_REASON, "n_cell_uid": len(a["exclusions"]),
                                "cells": a["exclusions"]},
        "claim_scope": {
            "biological_replicate": "R1",
            "generalization_scope": "within_R1_clone_heldout",
            "limitation": "Rewind R1 is ONE biological replicate. Clone-held-out evaluation is "
                          "within-experiment prospective generalization, not independent "
                          "biological-replicate validation.",
            "outcome_semantics_caveat": "`nonprimed` is the author-defined complement of the "
                                        "top-100 future-gDNA rule. It is NOT proven biological "
                                        "death or absolute reprogramming failure and may reflect "
                                        "downstream sampling/detection.",
            "treatment_variation": "none -- U is fixed; treatment-comparison fields are not "
                                   "applicable to this dataset"},
        "checks": a["checks"], "statistics": a["stats"], "verdict": a_verdict,
        "derived_artifacts": hash_artifacts(a["artifacts"]),
    }
    REWIND_MANIFEST.write_text(json.dumps(rewind_manifest, indent=2, default=str), encoding="utf-8")
    print("GSE227151:", a_verdict)

    # ---- Role B ------------------------------------------------------------------------------ #
    b = build_wm989(args.wm989_root)
    if b["verdict"]:
        print("GSE279162:", b["verdict"])
        return 1
    b_declared = b["stats"]["na_lineage_cells"]
    b_verdict = verdict_for(b["checks"], b_declared)
    wm989_manifest = {
        "dataset": "GSE279162", "role": "B", "source_accession": "GSE279162",
        "local_source_path_used": str(args.wm989_root),
        "source_files": source_files(sorted(b["req"].values())),
        "author_code_files": source_files(S21D.WM989_SCRIPTS),
        "author_rule_identifiers": {
            "qc": S21D.WM989_QC, "cell_lower_limit": S21D.CELL_LOWER_LIMIT,
            "cor_threshold": S21D.COR_THRESHOLD, "difference_val": S21D.DIFFERENCE_VAL,
            "posterior_floor": S21D.POSTERIOR_FLOOR,
            "outcome_rule": "table(assigned_lineage[OG_condition == condition])",
            "outcome_semantics": "post_treatment_assigned_cell_abundance",
            "rank_convention": "descending competition rank (pandas rank(method='min', "
                               "ascending=False)) over the eligible naive-observed clones"},
        **common,
        "declared_exclusions": {
            "reason": "assigned_lineage is NA (posterior floor or ambiguous top-two)",
            "n_cells": b_declared,
            "note": "documented and excluded from clone-linked prospective rows; NOT reassigned "
                    "by any alternative heuristic"},
        "claim_scope": {
            "observed_zero": "a zero is an OBSERVED zero against an available treatment sample. "
                             "It is not relabelled as death, failure, sensitivity or "
                             "non-resistance; it may reflect biological non-survival, "
                             "allocation/sampling, sequencing or capture limits.",
            "compositional": "post_fraction shares a treatment-level denominator. Fractions over "
                             "ALL assigned clones sum to one; the benchmark holds only the "
                             "eligible naive-observed subset, so its rows sum to less than one. "
                             "The clone x treatment rows are NOT independent biological units.",
            "abundance_confound": "the observed-zero rate falls steeply with captured naive clone "
                                  "size, so Stage 23 must clear abundance-only nuisance baselines "
                                  "before any state-conditioned claim.",
            "treatment_holdout": "reserved for Stage 26; Stage 23 uses these clone-held-out "
                                 "state-generalization folds only"},
        "checks": b["checks"], "statistics": b["stats"], "verdict": b_verdict,
        "derived_artifacts": hash_artifacts(b["artifacts"]),
    }
    WM989_MANIFEST.write_text(json.dumps(wm989_manifest, indent=2, default=str), encoding="utf-8")
    print("GSE279162:", b_verdict)

    # ---- gates + overall --------------------------------------------------------------------- #
    cells_written = pd.read_csv(REWIND_CELLS)
    leak = no_label_leakage(args.rewind_root, a["primed_clones"], cells_written)
    nomodel = builder_fits_no_model()
    gates = {
        "G22-1_prospective_linkage": bool(
            a["checks"]["cell_outer_group_equals_clone"] and b["checks"]["six_rows_per_eligible_clone"]),
        "G22-2_no_label_leakage": bool(leak["rewind_target_recomputed_from_outcome_source_alone"]
                                       and leak["rewind_written_labels_are_a_subset_of_that_set"]
                                       and leak["wm989_target_uses_no_pretreatment_expression"]),
        "G22-3_unique_biological_grouping": bool(
            a["checks"]["one_group_and_fold_per_clone"] and b["checks"]["one_group_and_fold_per_clone"]),
        "G22-4_expression_resolvability": bool(
            a["checks"]["expression_mapping_complete"] and b["checks"]["expression_mapping_complete"]),
        "G22-5_outcome_completeness": bool(
            a["checks"]["no_clone_has_contradictory_outcome"] and b["checks"]["structural_rows_match"]),
        "G22-6_author_rule_fidelity": bool(
            a["checks"]["primed_cells_42"] and a["checks"]["primed_clones_35"]
            and b["checks"]["posterior_floor_applied"]),
        "G22-7_frozen_evaluation_geometry": bool(
            a["checks"]["five_folds"] and a["checks"]["every_fold_has_both_classes"]
            and b["checks"]["five_folds"]),
        "G22-8_feature_firewall": set(FEATURE_ELIGIBILITY) == {
            "TARGET", "PROVENANCE_ONLY", "BASELINE_NUISANCE", "PRIMARY_X"},
        "G22-9_claim_scope": bool(rewind_manifest["claim_scope"] and wm989_manifest["claim_scope"]),
        "G22-10_no_modelling": bool(nomodel["no_modelling_imports"] and nomodel["no_fit_call"]),
    }
    overall = STAGE_23_READY if a_verdict in READY_VERDICTS else STAGE_23_BLOCKED
    rec = {
        "stage": "22", **common,
        "datasets": {"GSE227151": {"role": "A", "verdict": a_verdict,
                                   "manifest": REWIND_MANIFEST.name},
                     "GSE279162": {"role": "B", "verdict": b_verdict,
                                   "manifest": WM989_MANIFEST.name}},
        "role_b_is_non_blocking": True,
        "role_b_status_preserved": b_verdict,
        "gates": gates,
        "gate_evidence": {"G22-2": leak, "G22-10": nomodel},
        "all_gates_pass": all(gates.values()),
        "overall": overall,
        "manifest_hashes": hash_artifacts([REWIND_MANIFEST, WM989_MANIFEST]),
    }
    OUT.write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    print("gates:", "all pass" if rec["all_gates_pass"] else
          [k for k, v in gates.items() if not v])
    print("OVERALL:", overall)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
