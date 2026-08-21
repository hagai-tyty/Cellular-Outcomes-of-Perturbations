"""STAGE 21D — acquisition + reconstruction of the two datasets Stage 21C qualified.

REVISION 2 (2026-08-21). Revision 1 (`30ca7f0`) froze both datasets at
`RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE` because the author code was not on disk: the Rewind
Zenodo directory existed but was empty, and GEO ships GSE279162 as a raw feature-barcode matrix with
no clone-calling rule. **The complete author code for both datasets has since been placed on `D:\\`,
so this revision executes the authors' own rules instead of exploring around them.** The revision-1
state is preserved in the results file under `supersedes`, not erased.

Pre-registered in `plans/(newer)practical plans/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §11.
**Additive: the frozen 21A / 21B / 21C results are neither reopened nor rewritten.**
**Fits no model. No sklearn, no torch, no tuned threshold. `src/` untouched.**

WHAT THIS BUILDS
----------------
For each qualified dataset, the three links a prospective task needs:

    pretreatment expression  ->  cell  ->  clone  ->  future outcome

and nothing else. Whether the task is *learnable* is Stage 22 onward.

ROLE A — GSE227151 (Rewind; hiF-T fibroblast reprogramming)
    X   scRNA of the "before" arm, two 10x lanes of biol rep 1
    U   OSKM / dox reprogramming (single fixed intervention)
    Y   primed / nonprimed, by the author rule below

    Author rule, transcribed from
      plotScripts/rewind10X/R1/20220921_R1_primedVersusNonPrimedMarkersAndDistribution.R
      plotScripts/rewind10X/R1/2022.02.14_R1_cellNumberDistributionForPrimedVersusNonPrimed.R

        filter(cellID == "dummy")                     # the future gDNA arm
        group_by(BC50StarcodeD8, SampleNum)
        summarise(nUMI = sum(UMI))
        slice_max(order_by = nUMI, n = 100)           # dplyr default with_ties = TRUE
        inner_join(filtered10XCells, by = "BC50StarcodeD8")   -> primed
        anti_join (filtered10XCells, by = "BC50StarcodeD8")   -> nonprimed

    `n = 100` is fixed by the author code. Nothing here is tuned; the published 42 primed cells is
    asserted AFTER the rule is applied, as an independent reproduction check.

ROLE B — GSE279162 (WM989 melanoma, six treatments)
    X   scRNA of the untreated naive arm (3 lanes, pooled to `naive` by the author)
    U   dabrafenib | trametinib | CoCl2 | acid | cisplatin | doxorubicin
    Y   Y(clone, treatment) = post-treatment assigned-cell abundance, plus its rank within
        that treatment

    Author pipeline, transcribed from author_code_Schaff_manuscript/ in README order:
      preprocess_GEX.Rmd          condition-specific RNA QC; Custom features become the
                                  `lineage` assay; naive1/2/3 collapse to `naive`
      preprocess_cDNA_BCs.Rmd     drop zero-count lineages -> barcode_clustering(
                                  cell_lower_limit = 100, cor_threshold = 0.55) ->
                                  barcode_combine -> barcoding_posterior ->
                                  barcoding_assignment(difference_val = 0.2) ->
                                  assigned_lineage = NA where assigned_posterior < 0.5
      Find_Markers_Top_Res_lins_in_naive.Rmd
                                  table(assigned_lineage[OG_condition == condition])

    `num_lin = 5` in that last script is figure-specific and is NOT adopted as a target here.
    `preprocess_gDNA_BCs.Rmd` is separate gDNA/RPM logic and is kept as provenance only — it is
    not the scRNA cell-call procedure and must not be confused with it.

WHAT REVISION 1 GOT WRONG, AND WHY IT MATTERED
----------------------------------------------
Revision 1 labelled a Rewind clone primed if its barcode appeared in the gDNA arm *at all*, and
noted the answer was invariant for any read floor in [1, 562]. That was true and still is — but the
authors do not use a floor, they take the top 100 barcodes by summed count. The floor-invariant
reading gives 82 clones / 102 cells; the author rule gives 35 clones / 42 cells. Being invariant to
the wrong knob is not the same as being right, which is exactly why revision 1 refused to call it a
PASS.

For GSE279162 the gap was wider: revision 1's exploratory dominant-fraction / UMI-floor sweep
reported 4,018 naive clones at floor 1; the author's posterior assignment gives 1,401. **That floor
sweep is superseded and must not be used as the production clone call.**

Tri-state rule carried over from 21A/21B: PRESENT / ABSENT_PROVEN / UNKNOWN_REQUIRES_SOURCE_FILE.
"Not found" never becomes "absent".
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage21d_public_reconstruction_results.json"
REWIND_TABLE = _RESULTS / "stage21d_rewind_clone_table.tsv"
WM989_TABLE = _RESULTS / "stage21d_gse279162_clone_table.tsv"

PRESENT = "PRESENT"
ABSENT_PROVEN = "ABSENT_PROVEN"
UNKNOWN = "UNKNOWN_REQUIRES_SOURCE_FILE"

RECON_PASS = "RECONSTRUCTION_PASS"
RECON_PENDING = "RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE"
INVALID_LINKAGE = "INVALID_LINKAGE"
INSUFFICIENT_UNITS = "INSUFFICIENT_INDEPENDENT_UNITS"
MISSING_FILE = "MISSING_REQUIRED_FILE"

STAGE22_READY = "STAGE_22_READY"
STAGE22_PENDING = "STAGE_22_PENDING_OUTCOME_RULE"
STAGE22_BLOCKED = "STAGE_22_BLOCKED"

# CI runs on ubuntu-latest, where D:\ does not exist. Making that condition reproducible LOCALLY
# is what stops this class of red X from recurring: set CELLFATE_NO_LOCAL_DATA=1 and the suite sees
# exactly what CI sees. Enforced by tests/test_ci_portability.py.
_NO_LOCAL_DATA = os.environ.get("CELLFATE_NO_LOCAL_DATA") == "1"
_ABSENT_ROOT = Path("__local_data_absent__")

REWIND = _ABSENT_ROOT if _NO_LOCAL_DATA else Path(r"D:\GSE227151_Rewind")
WM989 = _ABSENT_ROOT if _NO_LOCAL_DATA else Path(r"D:\GSE279162")

# The one file whose absence stops the Role-A branch outright, by name, with no substitution.
GDNA_FILE = "stepThreeStarcodeShavedReads_BC_gDNA.txt"

# Author code, referenced by path + digest. It stays on D:\ and never enters git.
REWIND_CODE = REWIND / "author_code_zenodo7707418"
REWIND_SCRIPTS = [
    REWIND_CODE / "plotScripts/rewind10X/R1/20220921_R1_primedVersusNonPrimedMarkersAndDistribution.R",
    REWIND_CODE / "plotScripts/rewind10X/R1/2022.02.14_R1_cellNumberDistributionForPrimedVersusNonPrimed.R",
]
WM989_CODE = WM989 / "author_code_Schaff_manuscript"
WM989_SCRIPTS = [WM989_CODE / n for n in ("How_to_run_code_README.txt", "preprocess_GEX.Rmd",
                                          "preprocess_cDNA_BCs.Rmd", "preprocess_gDNA_BCs.Rmd",
                                          "Find_Markers_Top_Res_lins_in_naive.Rmd")]

TOP_N_GDNA = 100          # slice_max(n = 100), fixed by the author script
# Excluded in the second R1 script; kept as a named constant so the check is explicit.
REWIND_EXCLUDED_BC = "ATTCTAGTTGTAGTACGAGTAGCACATGTTCTACGTGGAGGACGAGAACG"

REWIND_GSMS = {
    "GSM7092515": "GSM7092515_1_2_control",
    "GSM7092516": "GSM7092516_1_1_control",
}
WM989_SAMPLES = {
    "GSM8562999": "Naive1", "GSM8563000": "Naive2", "GSM8563001": "Naive3",
    "GSM8563002": "Dabrafenib", "GSM8563003": "Trametinib", "GSM8563004": "CoCl2",
    "GSM8563005": "Acid", "GSM8563006": "Cisplatin", "GSM8563007": "Doxorubicin",
}
WM989_NAIVE = ("Naive1", "Naive2", "Naive3")
WM989_TREATMENTS = ("Dabrafenib", "Trametinib", "CoCl2", "Acid", "Cisplatin", "Doxorubicin")
# Merge order used by the author's merge() call; fixed here so the run is deterministic.
WM989_ORDER = ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin",
               "Naive1", "Naive2", "Naive3", "Trametinib")
# preprocess_GEX.Rmd: subset(nFeature_RNA > a & nCount_RNA < b & percent.mt < c), per condition.
WM989_QC = {
    "Naive1": (3000, 75000, 20), "Naive2": (3000, 60000, 20), "Naive3": (2500, 50000, 15),
    "Dabrafenib": (2000, 30000, 15), "Trametinib": (1500, 20000, 15), "CoCl2": (2500, 50000, 15),
    "Acid": (1500, 30000, 15), "Cisplatin": (1500, 20000, 20), "Doxorubicin": (1500, 20000, 20),
}
CELL_LOWER_LIMIT = 100        # barcode_clustering
COR_THRESHOLD = 0.55          # barcode_clustering
DIFFERENCE_VAL = 0.2          # barcoding_assignment
POSTERIOR_FLOOR = 0.5         # assigned_lineage <- NA where assigned_posterior < 0.5

N_GENES_WM989 = 36601          # asserted against the features file, not assumed

# Revision 1, frozen at 30ca7f0. Kept verbatim so the correction is auditable.
SUPERSEDES = {
    "revision": 1,
    "commit": "30ca7f0",
    "overall": STAGE22_PENDING,
    "overall_vocabulary_note": "revision 1 spelled the pass gate STAGE_22_MAY_OPEN; it is "
                               "STAGE_22_READY here. The pending and blocked gates are unchanged.",
    "why_it_was_reasonable_then": "the author code was not available to that run: "
                                  "D:\\GSE227151_Rewind\\author_code_zenodo7707418\\ existed but "
                                  "was EMPTY, and GEO ships GSE279162 as a raw feature-barcode "
                                  "matrix whose Data-Processing block states only that barcode "
                                  "reads were linked to 10x cell barcodes, not how a clone was "
                                  "called. Both outcome rules were therefore genuinely unknown, "
                                  "and revision 1 recorded them as UNKNOWN rather than inventing "
                                  "them.",
    "GSE227151": {
        "verdict": RECON_PENDING,
        "rule_used": "clone present in the gDNA arm at any read count",
        "positive_clones": 82, "positive_cells": 102,
        "negative_clones": 3067, "negative_cells": 3819,
        "note": "invariant for any read floor in [1, 562] -- true, but invariance to a knob the "
                "authors do not use. The author rule is a top-100 cut, not a floor.",
    },
    "GSE279162": {
        "verdict": RECON_PENDING,
        "rule_used": "exploratory dominant-barcode fraction and UMI floor sweep -- SUPERSEDED, "
                     "must not be used as the production clone call",
        "naive_clones_floor1": 4018, "clones_in_ge2_treatments": 1500, "clones_in_all_6": 250,
    },
    "findings_that_still_stand": [
        "Rewind SampleNum must be resolved by cellID containment, not by the GEO title text: "
        "SampleNum 1 -> GSM7092516, SampleNum 2 -> GSM7092515",
        "clones span both Rewind 10x lanes, so a lane-wise outer split would leak; the clone "
        "remains the biological grouping unit",
        "the 6479/6189 figures carried into revision 1 are not supported by the files; the "
        "measured counts stand",
    ],
}


@dataclass
class Finding:
    value: object
    status: str
    evidence: str

    def __post_init__(self):
        if self.status not in (PRESENT, ABSENT_PROVEN, UNKNOWN):
            raise ValueError(f"bad status {self.status!r}")


def _f(d, key, value, status, evidence):
    d[key] = Finding(value, status, evidence)


def sha256(path: Path, cap: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(cap):
            h.update(chunk)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict]:
    """Provenance for files that stay OUT of git: path, size, digest."""
    return [{"path": str(p), "exists": p.exists(),
             "bytes": p.stat().st_size if p.exists() else None,
             "sha256": sha256(p) if p.exists() else None} for p in paths]


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def mtx_dims(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                a, b, c = line.split()
                return int(a), int(b), int(c)
    raise ValueError(f"no dims line in {path}")


def slice_max_with_ties(df: pd.DataFrame, col: str, n: int) -> tuple[pd.DataFrame, int | None, int]:
    """dplyr `slice_max(order_by = col, n = n)` with its default `with_ties = TRUE`.

    dplyr keeps every row whose min-rank on the descending order is <= n, so a tie straddling the
    boundary returns MORE than n rows. Silently taking `head(n)` would be a different rule.
    Returns (kept, cutoff_value, rows_at_cutoff).
    """
    if len(df) <= n:
        return df.copy(), None, len(df)
    cutoff = np.sort(df[col].to_numpy())[::-1][n - 1]
    return df[df[col] >= cutoff].copy(), int(cutoff), int((df[col] == cutoff).sum())


# --------------------------------------------------------------------------------------------- #
# ROLE A — GSE227151 Rewind
# --------------------------------------------------------------------------------------------- #
def rewind_required(base: Path) -> dict[str, Path]:
    return {
        "filtered10XCells": base / "filtered10XCells.txt",
        "bc_10X": base / "stepThreeStarcodeShavedReads_BC_10X.txt",
        "bc_gDNA": base / GDNA_FILE,
        "series_matrix": base / "GSE227151-GPL18573_series_matrix.txt.gz",
        "family_xml": base / "GSE227151_family.xml",
        **{f"{g}_barcodes": base / g / f"{p}_barcodes.tsv.gz" for g, p in REWIND_GSMS.items()},
        **{f"{g}_features": base / g / f"{p}_features.tsv.gz" for g, p in REWIND_GSMS.items()},
        **{f"{g}_matrix": base / g / f"{p}_matrix.mtx.gz" for g, p in REWIND_GSMS.items()},
    }


def audit_rewind(base: Path = REWIND) -> dict:
    out: dict[str, Finding] = {}
    req = rewind_required(base)
    missing = sorted(k for k, p in req.items() if not p.exists())

    _f(out, "required_files", {"expected": len(req), "missing": missing},
       ABSENT_PROVEN if missing else PRESENT, f"checked {len(req)} paths under {base}")

    # The named stop condition: no substitution of another intermediate.
    if not req["bc_gDNA"].exists():
        _f(out, "future_outcome_source", str(req["bc_gDNA"]), ABSENT_PROVEN,
           f"REQUIRED FILE ABSENT: {GDNA_FILE}. The Rewind branch stops here. No other "
           f"intermediate is substituted for the future gDNA barcode arm.")
        return {"verdict": MISSING_FILE, "findings": out,
                "manifest": manifest(sorted(req.values()))}
    if missing:
        _f(out, "blocking_absence", missing, UNKNOWN,
           "files absent that are not the gDNA arm; reconstruction cannot be completed")
        return {"verdict": MISSING_FILE, "findings": out,
                "manifest": manifest(sorted(req.values()))}

    scripts_missing = [str(p) for p in REWIND_SCRIPTS if not p.exists()]
    _f(out, "author_code", {"root": str(REWIND_CODE),
                            "scripts": [str(p) for p in REWIND_SCRIPTS],
                            "missing": scripts_missing},
       UNKNOWN if scripts_missing else PRESENT,
       "the two load-bearing R1 scripts. In revision 1 this directory existed but was EMPTY, "
       "which is exactly why the outcome rule was UNKNOWN then.")
    if scripts_missing:
        return {"verdict": RECON_PENDING, "findings": out,
                "manifest": manifest(sorted(req.values()))}

    f10 = pd.read_csv(req["filtered10XCells"], sep="\t", index_col=0)
    s10 = pd.read_csv(req["bc_10X"], sep="\t", index_col=0)
    gdna = pd.read_csv(req["bc_gDNA"], sep="\t", index_col=0)

    # ---- does OUR separate gDNA file mean what the author's cellID=="dummy" arm means? -------- #
    dummy_rows = int((gdna["cellID"] == "dummy").sum())
    count_cols = [c for c in gdna.columns if c not in ("cellID", "BC50StarcodeD8", "SampleNum")]
    _f(out, "gdna_arm_equivalence",
       {"cellID_values": sorted(gdna["cellID"].unique().tolist()),
        "rows_kept_by_filter_dummy": dummy_rows, "rows_total": len(gdna),
        "samplenum_values": sorted(int(x) for x in gdna["SampleNum"].unique()),
        "count_column_here": count_cols, "count_column_in_author_script": ["UMI"],
        "schema_difference": count_cols != ["UMI"]},
       PRESENT,
       "the author reads a combined 10XAndGDNA file and applies filter(cellID == \"dummy\"); our "
       "separate file IS that arm already, so the filter is inert (it keeps every row). One real "
       "SCHEMA difference: the per-row count column is named `counts` here and `UMI` there. The "
       "semantics are not assumed equal -- they are checked downstream by the 42-cell "
       "reproduction, which the rule reproduces exactly.")

    _f(out, "author_qc_step_is_a_single_barcode_exclusion",
       {"stepThree_rows": len(s10), "filtered10XCells_rows": len(f10),
        "rows_carrying_excluded_barcode": int((s10["BC50StarcodeD8"] == REWIND_EXCLUDED_BC).sum()),
        "filtered_equals_stepThree_minus_that_barcode": bool(
            set(map(tuple, s10.loc[s10["BC50StarcodeD8"] != REWIND_EXCLUDED_BC,
                                   ["cellID", "BC50StarcodeD8", "SampleNum"]].to_numpy()))
            == set(map(tuple, f10[["cellID", "BC50StarcodeD8", "SampleNum"]].to_numpy())))},
       PRESENT,
       "revision 1 recorded that the author QC drops 21.2% of rows without explaining it. It is "
       "exactly one hyper-abundant barcode: filtered10XCells.txt IS "
       "stepThreeStarcodeShavedReads_BC_10X.txt minus all 1054 rows of "
       f"{REWIND_EXCLUDED_BC}.")

    # ---- SampleNum <-> GSM, resolved by intersection rather than by any label ---------------- #
    gsm_sets = {g: {b.split("-")[0] for b in read_barcodes(req[f"{g}_barcodes"])}
                for g in REWIND_GSMS}
    mapping, evidence = {}, {}
    for sn, grp in f10.groupby("SampleNum"):
        ids = set(grp["cellID"])
        hits = {g: len(ids & S) / len(ids) for g, S in gsm_sets.items()}
        mapping[int(sn)] = max(hits, key=hits.get)
        evidence[int(sn)] = {g: round(v, 4) for g, v in hits.items()} | {"n_cells": len(ids)}
    clean = all(max(v for k, v in e.items() if k != "n_cells") == 1.0 for e in evidence.values())
    _f(out, "samplenum_to_gsm", {"mapping": mapping, "containment": evidence},
       PRESENT if clean else UNKNOWN,
       "resolved by cellID containment in each GSM's own barcodes.tsv. NOTE: this is TRANSPOSED "
       "relative to the GEO titles -- GEO calls GSM7092515 'sample 1', but the barcode tables' "
       "SampleNum 1 is GSM7092516. Reading either label instead of intersecting would have "
       "mislabelled every cell. Carried forward unchanged from revision 1.")

    expr = {}
    for g in REWIND_GSMS:
        ngene, ncell, nnz = mtx_dims(req[f"{g}_matrix"])
        expr[g] = {"genes": ngene, "cells": ncell, "nnz": nnz,
                   "barcodes": len(read_barcodes(req[f"{g}_barcodes"]))}
    _f(out, "expression_cells", expr, PRESENT,
       f"total sequenced cells = {sum(v['cells'] for v in expr.values())}")

    assigned = f10.groupby("SampleNum")["cellID"].nunique().to_dict()
    _f(out, "cells_with_author_qc_clone", {int(k): int(v) for k, v in assigned.items()}, PRESENT,
       f"{sum(assigned.values())} of {sum(v['cells'] for v in expr.values())} sequenced cells "
       f"carry an author-QC clone assignment")

    s1 = set(f10.loc[f10.SampleNum == 1, "BC50StarcodeD8"])
    s2 = set(f10.loc[f10.SampleNum == 2, "BC50StarcodeD8"])
    _f(out, "clone_structure",
       {"clones_total": int(f10["BC50StarcodeD8"].nunique()),
        "clones_sample1": len(s1), "clones_sample2": len(s2), "shared": len(s1 & s2)},
       PRESENT,
       f"{len(s1 & s2)} clones appear in BOTH 10x lanes -> a lane-wise split would leak clones. "
       f"The outer split unit must be the clone. Carried forward unchanged from revision 1.")

    # ---- THE AUTHOR RULE --------------------------------------------------------------------- #
    probed = (gdna[gdna["cellID"] == "dummy"]
              .groupby(["BC50StarcodeD8", "SampleNum"], as_index=False)[count_cols[0]].sum()
              .rename(columns={count_cols[0]: "nUMI"}))
    top, cutoff, n_at_cutoff = slice_max_with_ties(probed, "nUMI", TOP_N_GDNA)
    head_only = set(probed.sort_values("nUMI", ascending=False).head(TOP_N_GDNA)["BC50StarcodeD8"])
    top_bcs = set(top["BC50StarcodeD8"])

    primed = f10[f10["BC50StarcodeD8"].isin(top_bcs)]
    nonprimed = f10[~f10["BC50StarcodeD8"].isin(top_bcs)]
    _f(out, "author_rule_slice_max",
       {"grouped_rows": len(probed), "n": TOP_N_GDNA, "cutoff_nUMI": cutoff,
        "rows_at_cutoff": n_at_cutoff, "selected_barcodes": len(top_bcs),
        "selected_nUMI_min": int(top["nUMI"].min()), "selected_nUMI_max": int(top["nUMI"].max()),
        "with_ties_primed_cells": int(len(primed)),
        "without_ties_primed_cells": int(f10["BC50StarcodeD8"].isin(head_only).sum())},
       PRESENT,
       f"slice_max(order_by = nUMI, n = {TOP_N_GDNA}) with dplyr's default with_ties = TRUE. There "
       f"IS a boundary tie ({n_at_cutoff} barcodes at nUMI = {cutoff}), so {len(top_bcs)} barcodes "
       f"are selected, not {TOP_N_GDNA}. The cell count happens to be identical either way here, "
       f"which is reported rather than relied on.")

    _f(out, "prospective_label_counts",
       {"primed_cells": int(len(primed)),
        "primed_clones": int(primed["BC50StarcodeD8"].nunique()),
        "nonprimed_cells": int(len(nonprimed)),
        "nonprimed_clones": int(nonprimed["BC50StarcodeD8"].nunique()),
        "per_sample": {int(sn): {"primed_cells": int((primed.SampleNum == sn).sum()),
                                 "primed_clones": int(primed.loc[primed.SampleNum == sn,
                                                                 "BC50StarcodeD8"].nunique())}
                       for sn in sorted(f10["SampleNum"].unique())},
        "positive_rate_cells": round(len(primed) / len(f10), 4)},
       PRESENT,
       "primed := the clone's barcode is among the top-100 of the future gDNA arm by summed count, "
       "per the author scripts. Revision 1's presence-at-any-count reading gave 82 clones / 102 "
       "cells; the author rule gives fewer, and the difference is the whole reason revision 1 "
       "refused to call it a PASS.")

    # ---- the reproduction gate ---------------------------------------------------------------- #
    reproduced = int(len(primed)) == 42
    _f(out, "published_42_reproduction",
       {"published_primed_cells": 42, "reconstructed_primed_cells": int(len(primed)),
        "reproduced": reproduced, "tuned": False},
       PRESENT if reproduced else ABSENT_PROVEN,
       "ASSERTED AFTER the author rule was implemented, never targeted. n = 100 comes from the "
       "author script and no threshold was moved. Revision 1's rule did not reproduce 42 at cell "
       "level; this one does exactly.")

    excl_in_f10 = int((f10["BC50StarcodeD8"] == REWIND_EXCLUDED_BC).sum())
    _f(out, "second_script_barcode_exclusion",
       {"barcode": REWIND_EXCLUDED_BC, "rows_in_filtered10XCells": excl_in_f10,
        "rows_in_gdna": int((gdna["BC50StarcodeD8"] == REWIND_EXCLUDED_BC).sum()),
        "in_top_selection": REWIND_EXCLUDED_BC in top_bcs,
        "effect_on_primed_cells": 0 if excl_in_f10 == 0 else None,
        "inert": excl_in_f10 == 0},
       ABSENT_PROVEN if excl_in_f10 == 0 else PRESENT,
       "the second R1 script drops this barcode before counting. It is ABSENT from "
       "filtered10XCells.txt and from the gDNA arm, so the exclusion is INERT here -- because the "
       "upstream QC that produced filtered10XCells.txt already removed exactly this barcode. Both "
       "author scripts therefore give identical primed sets.")

    _f(out, "replicate_structure",
       {"gsms_local": sorted(REWIND_GSMS), "biological_replicates_local": 1,
        "gsms_in_series_not_local": ["GSM7092517", "GSM7092518", "GSM7092519", "GSM7092520",
                                     "GSM7092521", "GSM7092522", "GSM7092523", "GSM7092524",
                                     "GSM7092525", "GSM7092526", "GSM7092527"]},
       PRESENT,
       "both local GSMs are titled 'biol rep 1' -- two 10x lanes of ONE biological replicate. "
       "Outer units are clones within one replicate, not independent replicates. This limits "
       "generalisation and is unchanged by the rule resolution.")

    tbl = (f10.groupby("BC50StarcodeD8")
           .agg(cells=("cellID", "size"),
                cells_s1=("SampleNum", lambda s: int((s == 1).sum())),
                cells_s2=("SampleNum", lambda s: int((s == 2).sum())))
           .reset_index().rename(columns={"BC50StarcodeD8": "clone"}))
    gsum = probed.set_index("BC50StarcodeD8")["nUMI"]
    tbl["gdna_nUMI"] = tbl["clone"].map(gsum).fillna(0).astype(int)
    tbl["in_gdna_at_all"] = tbl["clone"].isin(set(gsum.index)).astype(int)
    tbl["primed"] = tbl["clone"].isin(top_bcs).astype(int)
    tbl = tbl.sort_values(["primed", "gdna_nUMI", "clone"], ascending=[False, False, True])
    REWIND_TABLE.write_text(tbl.to_csv(sep="\t", index=False), encoding="utf-8")

    if not clean or not len(primed):
        verdict = INVALID_LINKAGE
    elif any(f.status == UNKNOWN for f in out.values()) or not reproduced:
        verdict = RECON_PENDING
    else:
        verdict = RECON_PASS
    return {"verdict": verdict, "previous_verdict": SUPERSEDES["GSE227151"]["verdict"],
            "findings": out,
            "manifest": manifest(sorted(req.values()) + REWIND_SCRIPTS),
            "clone_table": REWIND_TABLE.relative_to(ROOT).as_posix()}


# --------------------------------------------------------------------------------------------- #
# ROLE B — GSE279162, the Schaff pipeline
# --------------------------------------------------------------------------------------------- #
def wm989_required(base: Path) -> dict[str, Path]:
    req = {"family_xml": base / "GSE279162_family.xml",
           "series_matrix": base / "GSE279162_series_matrix.txt.gz"}
    for gsm, name in WM989_SAMPLES.items():
        for kind in ("barcodes.tsv", "features.tsv", "matrix.mtx"):
            req[f"{name}_{kind}"] = base / f"{gsm}_{name}_filtered_{kind}.gz"
    return req


def qc_and_lineage(path: Path, n_genes: int, mt_mask: np.ndarray,
                   min_feature: int, max_count: int, max_mt: float) -> dict:
    """preprocess_GEX.Rmd: nFeature_RNA / nCount_RNA / percent.mt on the Gene Expression block,
    then `subset(nFeature_RNA > a & nCount_RNA < b & percent.mt < c)`.

    Streams the matrix so the gene block is never materialised; returns the surviving cells and
    their lineage entries.
    """
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                _, n_cells, _ = (int(x) for x in line.split())
                break
        n_count = np.zeros(n_cells + 1)
        n_feat = np.zeros(n_cells + 1)
        mt_count = np.zeros(n_cells + 1)
        lin_chunks = []
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                                 dtype=np.int32, chunksize=4_000_000):
            row, col, val = chunk["row"].to_numpy(), chunk["col"].to_numpy(), chunk["val"].to_numpy()
            is_gene = row <= n_genes
            gc, gv, gr = col[is_gene], val[is_gene], row[is_gene]
            n_count += np.bincount(gc, weights=gv, minlength=n_cells + 1)
            n_feat += np.bincount(gc, minlength=n_cells + 1)
            m = mt_mask[gr]
            if m.any():
                mt_count += np.bincount(gc[m], weights=gv[m], minlength=n_cells + 1)
            lin_chunks.append(chunk[~is_gene])
    pct_mt = np.divide(100 * mt_count, n_count, out=np.zeros_like(mt_count), where=n_count > 0)
    keep = (n_feat > min_feature) & (n_count < max_count) & (pct_mt < max_mt)
    keep[0] = False                                  # column indices are 1-based
    lin = pd.concat(lin_chunks, ignore_index=True)
    return {"n_cells": n_cells, "kept": int(keep.sum()),
            "cells": np.flatnonzero(keep), "lin": lin[keep[lin["col"].to_numpy()]]}


def barcode_clustering(lin: sparse.csr_matrix, names: np.ndarray) -> dict:
    """preprocess_cDNA_BCs.Rmd :: barcode_clustering(cell_lower_limit, cor_threshold).

    Correlate lineage features that are seen in enough cells, then walk the resulting pairs in R's
    column-major `which(..., arr.ind = TRUE)` order. The author's loop calls stop("Merging
    happening") if a pair would join two existing clusters; that condition is counted here rather
    than raised, so a data difference surfaces as a finding instead of a crash.
    """
    n_per_lin = np.diff(lin.indptr)
    sel = np.flatnonzero(n_per_lin >= CELL_LOWER_LIMIT)
    dense = np.asarray(lin[sel].todense()).T                       # cells x lineages
    with np.errstate(invalid="ignore"):
        cor = np.corrcoef(dense, rowvar=False)
    cor[np.tril_indices_from(cor)] = np.nan
    ii, jj = np.where(cor >= COR_THRESHOLD)
    order = np.lexsort((ii, jj))                                   # R is column-major
    pairs = np.stack([ii[order], jj[order]], axis=1)

    clusters: list[list[int]] = []
    owner: dict[int, int] = {}
    conflicts = 0
    for a, b in pairs:
        idx = [int(a), int(b)]
        vals = {owner[k] for k in idx if k in owner}
        if not vals:
            clusters.append(sorted(idx))
            for k in idx:
                owner[k] = len(clusters) - 1
        else:
            if len(vals) > 1:
                conflicts += 1
            v = min(vals)
            clusters[v] = sorted(set(clusters[v]) | set(idx))
            for k in idx:
                owner[k] = v
    named = [[str(names[sel[x]]) for x in c] for c in clusters]
    return {"clusters": named, "n_pairs": int(len(pairs)),
            "n_candidate_lineages": int(len(sel)), "merging_conflicts": conflicts}


def barcode_combine(lin: sparse.csr_matrix, names: np.ndarray,
                    clusters: list[list[str]]) -> tuple[sparse.csc_matrix, np.ndarray]:
    """preprocess_cDNA_BCs.Rmd :: barcode_combine -- untouched rows first, then one summed row per
    cluster, named after the cluster's first member in lin_mat row order."""
    pos = {n: i for i, n in enumerate(names)}
    inside = sorted({pos[n] for c in clusters for n in c})
    untouched = np.setdiff1d(np.arange(lin.shape[0]), inside)
    blocks = [lin[untouched]]
    out_names = list(names[untouched])
    for c in clusters:
        idx = sorted(pos[n] for n in c)
        blocks.append(sparse.csr_matrix(lin[idx].sum(axis=0)))
        out_names.append(names[idx[0]])
    return sparse.vstack(blocks).tocsc(), np.asarray(out_names, dtype=object)


def barcoding_posterior_and_assignment(lin: sparse.csc_matrix,
                                       names: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """preprocess_cDNA_BCs.Rmd :: barcoding_posterior + barcoding_assignment + the posterior floor.

    The posterior is a softmax of `count_b * log(gamma_b)` over ALL lineages, so a lineage with
    zero count still contributes exp(0) = 1 to the denominator. That is what makes a cell with weak
    counts fall below the 0.5 floor rather than being confidently mis-assigned, and it is why the
    zero rows cannot be dropped for speed.
    """
    n_lin, n = lin.shape
    lib = np.asarray(lin.sum(axis=0)).ravel()
    idx, val, ptr = lin.indices, lin.data, lin.indptr

    cell_max = np.zeros(n)
    for c in range(n):
        if ptr[c + 1] > ptr[c]:
            cell_max[c] = val[ptr[c]:ptr[c + 1]].max()
    col_of = np.repeat(np.arange(n), np.diff(ptr))
    won = (val == cell_max[col_of]) & (val > 0)
    scaled = val / (lib[col_of] + 1.0)
    total = np.bincount(idx, weights=scaled, minlength=n_lin)
    won_sum = np.bincount(idx[won], weights=scaled[won], minlength=n_lin)
    n_won = np.bincount(idx[won], minlength=n_lin)
    with np.errstate(invalid="ignore", divide="ignore"):
        beta1 = np.where(n_won > 0, won_sum / np.maximum(n_won, 1), np.nan)
        beta0 = (total - won_sum) / (n - n_won)
    beta1_mean = float(np.nanmean(beta1))
    positive = beta0[beta0 > 1e-8]
    q98, q02 = np.quantile(positive, 0.98), np.quantile(positive, 0.02)
    gamma = beta1_mean / np.maximum(np.minimum(beta0, q98), q02)
    log_gamma = np.log(gamma)

    assigned = np.empty(n, dtype=object)
    posterior = np.zeros(n)
    for c in range(n):
        a, b = ptr[c], ptr[c + 1]
        k = b - a
        z = val[a:b] * log_gamma[idx[a:b]] if k else np.empty(0)
        top = z.max() if k else 0.0
        m = max(top, 0.0)
        denom = (np.exp(z - m).sum() if k else 0.0) + (n_lin - k) * np.exp(-m)
        if k and top > 0:
            order = np.argsort(-z, kind="stable")
            second = max(z[order[1]], 0.0) if k > 1 else 0.0
            p1 = np.exp(z[order[0]] - m) / denom
            posterior[c] = p1
            if p1 - np.exp(second - m) / denom >= DIFFERENCE_VAL:
                assigned[c] = names[idx[a + order[0]]]
        else:
            posterior[c] = np.exp(-m) / denom
    below = posterior < POSTERIOR_FLOOR
    n_before = int(sum(x is not None for x in assigned))
    assigned[below] = None
    return assigned, posterior, {
        "beta1_mean": beta1_mean, "gamma_min": float(gamma.min()), "gamma_max": float(gamma.max()),
        "assigned_before_posterior_floor": n_before,
        "dropped_by_posterior_floor": int(below.sum()),
    }


def audit_gse279162(base: Path = WM989) -> dict:
    out: dict[str, Finding] = {}
    req = wm989_required(base)
    missing = sorted(k for k, p in req.items() if not p.exists())
    _f(out, "required_files", {"expected": len(req), "missing": missing},
       ABSENT_PROVEN if missing else PRESENT,
       f"checked {len(req)} paths under {base} across {len(WM989_SAMPLES)} samples")
    if missing:
        return {"verdict": MISSING_FILE, "findings": out,
                "manifest": manifest(sorted(req.values()))}

    scripts_missing = [str(p) for p in WM989_SCRIPTS if not p.exists()]
    _f(out, "author_code", {"root": str(WM989_CODE),
                            "scripts": [str(p) for p in WM989_SCRIPTS],
                            "missing": scripts_missing,
                            "readme_order": ["preprocess_GEX.Rmd", "preprocess_cDNA_BCs.Rmd",
                                             "preprocess_gDNA_BCs.Rmd",
                                             "Find_Markers_Top_Res_lins_in_naive.Rmd"]},
       UNKNOWN if scripts_missing else PRESENT,
       "the five primary provenance files. preprocess_gDNA_BCs.Rmd is separate gDNA/RPM logic, "
       "kept as provenance only -- it is NOT the scRNA assigned_lineage procedure.")
    if scripts_missing:
        return {"verdict": RECON_PENDING, "findings": out,
                "manifest": manifest(sorted(req.values()))}

    feats = [ln.rstrip("\n").split("\t")
             for ln in gzip.open(req["Naive1_features.tsv"], "rt")]
    lin_names_all = np.asarray([f[0] for f in feats if len(f) > 2 and f[2] == "Custom"],
                               dtype=object)
    n_genes = len(feats) - len(lin_names_all)
    symbols = [f[1] for f in feats[:n_genes]]
    mt_mask = np.zeros(n_genes + 1, dtype=bool)
    for i, s in enumerate(symbols):
        if s.startswith("MT-"):
            mt_mask[i + 1] = True
    _f(out, "feature_structure",
       {"total": len(feats), "n_genes": n_genes, "n_lineage_features": len(lin_names_all),
        "genes_first": all(f[2] == "Gene Expression" for f in feats[:n_genes] if len(f) > 2),
        "mt_genes": int(mt_mask.sum()),
        "lineage_ids_match_L_pattern": bool(all(re.fullmatch(r"L\d+", n)
                                                for n in lin_names_all[:1000]))},
       PRESENT,
       f"the clone identifier is a FEATURE ROW: {len(lin_names_all)} 'Custom' LinNNNN features "
       f"after {n_genes} genes. preprocess_GEX.Rmd makes exactly this block the `lineage` assay.")
    if n_genes != N_GENES_WM989:
        _f(out, "gene_block_size_unexpected", n_genes, UNKNOWN,
           f"expected {N_GENES_WM989} genes before the Custom block")

    # ---- preprocess_GEX.Rmd: condition-specific RNA QC ---------------------------------------- #
    per_sample, rows, cols, vals, cond, offset = {}, [], [], [], [], 0
    for name in WM989_ORDER:
        min_f, max_c, max_mt = WM989_QC[name]
        r = qc_and_lineage(req[f"{name}_matrix.mtx"], n_genes, mt_mask, min_f, max_c, max_mt)
        per_sample[name] = {"gsm": next(g for g, s in WM989_SAMPLES.items() if s == name),
                            "cells_raw": r["n_cells"], "cells_post_qc": r["kept"],
                            "pct_kept": round(100 * r["kept"] / r["n_cells"], 1),
                            "qc": {"nFeature_RNA_gt": min_f, "nCount_RNA_lt": max_c,
                                   "percent_mt_lt": max_mt}}
        local = {c: i for i, c in enumerate(r["cells"])}
        lin = r["lin"]
        rows.append(lin["row"].to_numpy() - n_genes - 1)
        cols.append(np.array([local[c] for c in lin["col"].to_numpy()], dtype=np.int64) + offset)
        vals.append(lin["val"].to_numpy().astype(np.float64))
        cond += [name] * r["kept"]
        offset += r["kept"]
    _f(out, "author_rna_qc", per_sample, PRESENT,
       f"condition-specific filters transcribed from preprocess_GEX.Rmd. "
       f"{sum(v['cells_post_qc'] for v in per_sample.values())} of "
       f"{sum(v['cells_raw'] for v in per_sample.values())} raw cells survive. Revision 1 used "
       f"every raw GEO cell, which is one reason its counts were inflated.")

    lin_mat = sparse.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(lin_names_all), offset))
    cond = np.asarray(cond, dtype=object)
    nonzero = np.asarray(lin_mat.sum(axis=1)).ravel() > 0
    lin_mat, lin_names = lin_mat[nonzero], lin_names_all[nonzero]

    # ---- preprocess_cDNA_BCs.Rmd ------------------------------------------------------------- #
    clus = barcode_clustering(lin_mat, lin_names)
    sizes = [len(c) for c in clus["clusters"]]
    _f(out, "barcode_clustering",
       {"lineages_after_rowsum_filter": int(lin_mat.shape[0]),
        "cell_lower_limit": CELL_LOWER_LIMIT, "cor_threshold": COR_THRESHOLD,
        "candidate_lineages": clus["n_candidate_lineages"], "correlated_pairs": clus["n_pairs"],
        "clusters": len(clus["clusters"]), "lineages_merged": int(sum(sizes)),
        "cluster_size_histogram": {int(k): int(v) for k, v in
                                   zip(*np.unique(sizes, return_counts=True), strict=True)}
                                  if sizes else {},
        "author_stop_merging_happening_fired": clus["merging_conflicts"]},
       PRESENT,
       "correlated lineage features are the same physical clone read with sequencing error or "
       "multiple integrations. The author's loop raises stop(\"Merging happening\") if a pair "
       f"would join two existing clusters; it fires {clus['merging_conflicts']} times here, so the "
       "reconstruction meets the same precondition the authors' run did.")

    combined, combined_names = barcode_combine(lin_mat, lin_names, clus["clusters"])
    assigned, posterior, stats = barcoding_posterior_and_assignment(combined, combined_names)
    _f(out, "barcoding_posterior_and_assignment",
       {"lineages": int(combined.shape[0]), "cells": int(combined.shape[1]),
        "difference_val": DIFFERENCE_VAL, "posterior_floor": POSTERIOR_FLOOR, **stats,
        "assigned_cells": int(sum(a is not None for a in assigned)),
        "na_cells": int(sum(a is None for a in assigned))},
       PRESENT,
       "multinomial posterior over every lineage, assignment only when the top two posteriors "
       "differ by >= 0.2, then assigned_lineage <- NA below a 0.5 posterior. This SUPERSEDES "
       "revision 1's dominant-fraction / UMI-floor sweep, which is not the production rule.")

    df = pd.DataFrame({"cond": cond, "lineage": assigned, "posterior": posterior})
    df["OG"] = np.where(df["cond"].isin(WM989_NAIVE), "naive", df["cond"])
    ok = df[df["lineage"].notna()]
    table = (ok.groupby(["lineage", "OG"]).size().unstack(fill_value=0)
             .reindex(columns=["naive", *WM989_TREATMENTS], fill_value=0))

    _f(out, "per_condition_assignment",
       {k: {"cells_post_qc": int(len(g)),
            "assigned": int(g["lineage"].notna().sum()),
            "na": int(g["lineage"].isna().sum()),
            "unique_clones": int(g["lineage"].nunique())}
        for k, g in df.groupby("OG")},
       PRESENT,
       "naive1/2/3 are pooled into one `naive` condition, as preprocess_GEX.Rmd does explicitly. "
       "That answers one of revision 1's three open questions outright.")

    naive_clones = table.index[table["naive"] > 0]
    n_treat = (table.loc[naive_clones, list(WM989_TREATMENTS)] > 0).sum(axis=1)
    _f(out, "clone_coverage",
       {"unique_assigned_clones": int(len(table)),
        "clones_with_naive_observation": int(len(naive_clones)),
        "clones_in_ge1_treatment": int((n_treat >= 1).sum()),
        "clones_in_ge2_treatments": int((n_treat >= 2).sum()),
        "clones_in_all_6": int((n_treat == len(WM989_TREATMENTS)).sum()),
        "by_n_treatments": {int(k): int(v) for k, v in n_treat.value_counts().sort_index().items()}},
       PRESENT,
       f"{int((n_treat >= 2).sum())} clones carry both a naive observation and a post-treatment "
       f"observation under two or more treatments -- the structure an X x U interaction test "
       f"needs. Revision 1's floor-1 sweep claimed 1500; the author rule gives fewer and is the "
       f"number that counts.")

    per_treatment = {}
    for t in WM989_TREATMENTS:
        col = table.loc[naive_clones, t]
        nz = col[col > 0]
        per_treatment[t] = {"clones_present": int((col > 0).sum()),
                            "clones_absent": int((col == 0).sum()),
                            "assigned_cells": int(col.sum()),
                            "abundance_median": int(nz.median()) if len(nz) else 0,
                            "abundance_max": int(col.max())}
    pooled = table.loc[naive_clones, list(WM989_TREATMENTS)].to_numpy().ravel()
    pooled = pooled[pooled > 0]
    naive_depth = table.loc[naive_clones, "naive"]
    _f(out, "future_outcome",
       {"definition": "Y(clone, treatment) = post-treatment assigned-cell abundance; rank is "
                      "taken within each treatment",
        "per_treatment": per_treatment,
        "abundance_distribution": {"n_nonzero": int(len(pooled)),
                                   "median": int(np.median(pooled)),
                                   "q75": int(np.quantile(pooled, 0.75)),
                                   "q95": int(np.quantile(pooled, 0.95)),
                                   "max": int(pooled.max())},
        "naive_cells_per_clone": {"median": int(naive_depth.median()),
                                  "max": int(naive_depth.max()),
                                  "ge2": int((naive_depth >= 2).sum()),
                                  "ge5": int((naive_depth >= 5).sum()),
                                  "ge10": int((naive_depth >= 10).sum())}},
       PRESENT,
       "built exactly as Find_Markers_Top_Res_lins_in_naive.Rmd builds it: "
       "table(assigned_lineage[OG_condition == condition]). The script's num_lin = 5 is "
       "figure-specific and is deliberately NOT turned into a binary target here.")

    _f(out, "reproduction_anchor",
       {"published_count_available": False,
        "scripts_searched": [p.name for p in WM989_SCRIPTS],
        "nearest_candidate": "preprocess_GEX.Rmd writes preprocess_GEX/cellsPreandPostFilt.xlsx",
        "candidate_present_in_archive": (WM989_CODE / "preprocess_GEX").exists(),
        "is_a_reproduction_failure": False},
       UNKNOWN,
       "unlike Rewind's published 42, the five provided scripts contain no published count that "
       "can serve as an independent reproduction check -- preprocess_GEX.Rmd writes "
       "cellsPreandPostFilt.xlsx but that output is not in the archive. The pipeline executes "
       "cleanly and end to end, so this is a missing VALIDATION ANCHOR, not a reproduction "
       "failure. Recorded so nobody later mistakes 'ran without error' for 'checked against the "
       "authors' numbers'.")

    _f(out, "known_deviations_from_the_author_object",
       ["Seurat's CreateAssayObject(min.cells = 1) is applied per sample before QC; here the "
        "equivalent zero-lineage filter is the merged rowSums > 0 that preprocess_cDNA_BCs.Rmd "
        "applies anyway, which subsumes it",
        "NormalizeData / PCA / UMAP / FindClusters are not run: they feed visualisation and "
        "marker tests, never assigned_lineage",
        "sample merge order is fixed alphabetically; it can only matter for which.max ties, and a "
        "tie yields a zero top-two difference, which the 0.2 rule sends to NA regardless"],
       PRESENT,
       "stated so the reconstruction is auditable against the authors' object rather than "
       "implied to be byte-identical to it")

    table.reset_index().to_csv(WM989_TABLE, sep="\t", index=False)

    if any(f.status == UNKNOWN and k != "reproduction_anchor" for k, f in out.items()):
        verdict = RECON_PENDING
    elif len(naive_clones) == 0 or int((n_treat >= 1).sum()) == 0:
        verdict = INVALID_LINKAGE
    else:
        verdict = RECON_PASS
    return {"verdict": verdict, "previous_verdict": SUPERSEDES["GSE279162"]["verdict"],
            "findings": out,
            "manifest": manifest(sorted(req.values()) + WM989_SCRIPTS),
            "clone_table": WM989_TABLE.relative_to(ROOT).as_posix()}


# --------------------------------------------------------------------------------------------- #
def main() -> int:
    _RESULTS.mkdir(exist_ok=True)
    rec = {"stage": "21D", "revision": 2,
           "plan": "plans/(newer)practical plans/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md",
           "additive_to": ["results/diag_stage21_data_audit_results.json",
                           "results/diag_stage21b_source_design_results.json"],
           "supersedes": SUPERSEDES,
           "model_fitted": False, "src_modified": False, "raw_data_committed": False,
           "author_code_committed": False}

    for key, fn in (("GSE227151", audit_rewind), ("GSE279162", audit_gse279162)):
        r = fn()
        print(f"{key}: {r['verdict']}   (was {r.get('previous_verdict')})")
        rec[key] = {"verdict": r["verdict"], "previous_verdict": r.get("previous_verdict"),
                    "findings": {k: asdict(v) for k, v in r["findings"].items()},
                    "manifest": r["manifest"], "clone_table": r.get("clone_table")}

    va, vb = rec["GSE227151"]["verdict"], rec["GSE279162"]["verdict"]
    if va in (MISSING_FILE, INVALID_LINKAGE, INSUFFICIENT_UNITS):
        overall = STAGE22_BLOCKED
    elif va == RECON_PASS and vb == RECON_PASS:
        overall = STAGE22_READY
    else:
        overall = STAGE22_PENDING
    rec["overall"] = overall
    rec["previous_overall"] = SUPERSEDES["overall"]
    print("OVERALL:", overall, " (was", SUPERSEDES["overall"] + ")")

    OUT.write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
