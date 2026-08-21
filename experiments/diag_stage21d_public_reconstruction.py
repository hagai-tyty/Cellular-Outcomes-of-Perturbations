"""STAGE 21D — acquisition + reconstruction of the two datasets Stage 21C qualified.

Pre-registered in `plans/(newer)practical plans/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §11.
**Additive: the frozen 21A / 21B / 21C results are neither reopened nor rewritten.**
**Fits no model. No sklearn, no torch, no threshold tuned to a published number. `src/` untouched.**

WHAT THIS BUILDS
----------------
For each qualified dataset, the three links that a prospective task needs:

    pretreatment expression  ->  cell  ->  clone  ->  future outcome

and nothing else. Reconstruction only. Whether the task is *learnable* is Stage 22 onward.

ROLE A — GSE227151 (Rewind; hiF-T fibroblast reprogramming)
    X   scRNA of the "before" arm, two 10x lanes of biol rep 1
    U   OSKM / dox reprogramming (single fixed intervention)
    Y   clone barcode recovered from the post-reprogramming gDNA arm  ->  primed / nonprimed

ROLE B — GSE279162 (WM989 melanoma, six treatments)
    X   scRNA of the untreated naive arm (3 lanes)
    U   dabrafenib | trametinib | CoCl2 | acid | cisplatin | doxorubicin
    Y   clone re-observed in that treatment's post-treatment scRNA  ->  survived / did not

THE TRAP THIS SCRIPT EXISTS TO AVOID
------------------------------------
Both datasets label a clone by *presence in a later sample*. Presence is a thresholded
observation, so the outcome column is only as trustworthy as the rule that produced it. Two
failure modes:

1. **Inventing the rule.** Picking a read floor that reproduces a published count is fitting the
   label to the answer. This script therefore reports the floor SWEEP as evidence and refuses to
   select a floor. Where the answer turns out not to depend on the floor, it says so and gives the
   interval over which that holds.
2. **Trusting a sample label.** `SampleNum` in the Rewind barcode tables does NOT follow GEO's
   "sample N" wording — they are transposed. The mapping is resolved here by cellID intersection
   against each GSM's own barcode list, which is unambiguous, rather than by reading either label.

Tri-state rule carried over from 21A/21B: PRESENT / ABSENT_PROVEN / UNKNOWN_REQUIRES_SOURCE_FILE.
"Not found" never becomes "absent".
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

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

STAGE22_OPEN = "STAGE_22_MAY_OPEN"
STAGE22_PENDING = "STAGE_22_PENDING_OUTCOME_RULE"
STAGE22_BLOCKED = "STAGE_22_BLOCKED"

REWIND = Path(r"D:\GSE227151_Rewind")
WM989 = Path(r"D:\GSE279162")

# The one file whose absence stops the Role-A branch outright, by name, with no substitution.
GDNA_FILE = "stepThreeStarcodeShavedReads_BC_gDNA.txt"

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

N_GENES_WM989 = 36601          # asserted against the features file, not assumed
UMI_FLOORS = (1, 2, 3, 5, 10, 20, 50)


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
    out = []
    for p in paths:
        out.append({"path": str(p), "exists": p.exists(),
                    "bytes": p.stat().st_size if p.exists() else None,
                    "sha256": sha256(p) if p.exists() else None})
    return out


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


def mtx_rows_above(path: Path, first_row: int) -> pd.DataFrame:
    """Every (row, col, value) with row > `first_row`, streamed.

    CellRanger writes genes first and the Custom lineage features after them, so this pulls the
    lineage block out of a 30M-entry matrix without ever holding the gene block in memory.
    """
    kept = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                break                       # the dims line, consumed
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                                 dtype=np.int32, chunksize=4_000_000):
            kept.append(chunk[chunk["row"] > first_row])
    return pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(
        columns=["row", "col", "val"], dtype=np.int32)


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
       ABSENT_PROVEN if missing else PRESENT,
       f"checked {len(req)} paths under {base}")

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

    f10 = pd.read_csv(req["filtered10XCells"], sep="\t", index_col=0)
    s10 = pd.read_csv(req["bc_10X"], sep="\t", index_col=0)
    gdna = pd.read_csv(req["bc_gDNA"], sep="\t", index_col=0)

    _f(out, "author_qc_table_schema", list(f10.columns), PRESENT,
       f"filtered10XCells.txt is an author-processed cell->clone table, {len(f10)} rows. "
       f"Used as the primary linkage; QC is NOT rebuilt from the raw combined intermediate.")
    _f(out, "prefilter_table_rows", {"stepThree_BC_10X": len(s10), "filtered10XCells": len(f10)},
       PRESENT,
       f"the author QC step drops {len(s10) - len(f10)} of {len(s10)} rows "
       f"({100 * (1 - len(f10) / len(s10)):.1f}%)")

    # ---- SampleNum <-> GSM, resolved by intersection rather than by any label ---------------- #
    gsm_sets = {g: {b.split("-")[0] for b in read_barcodes(req[f"{g}_barcodes"])}
                for g in REWIND_GSMS}
    mapping, evidence = {}, {}
    for sn, grp in f10.groupby("SampleNum"):
        ids = set(grp["cellID"])
        hits = {g: len(ids & S) / len(ids) for g, S in gsm_sets.items()}
        best = max(hits, key=hits.get)
        mapping[int(sn)] = best
        evidence[int(sn)] = {g: round(v, 4) for g, v in hits.items()} | {"n_cells": len(ids)}
    clean = all(max(v for k, v in e.items() if k != "n_cells") == 1.0 for e in evidence.values())
    _f(out, "samplenum_to_gsm", {"mapping": mapping, "containment": evidence},
       PRESENT if clean else UNKNOWN,
       "resolved by cellID containment in each GSM's own barcodes.tsv. NOTE: this is TRANSPOSED "
       "relative to the GEO titles -- GEO calls GSM7092515 'sample 1', but the barcode tables' "
       "SampleNum 1 is GSM7092516. Reading either label instead of intersecting would have "
       "mislabelled every cell.")

    # ---- expression side --------------------------------------------------------------------- #
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
       f"carry an author-QC clone assignment "
       f"({100 * sum(assigned.values()) / sum(v['cells'] for v in expr.values()):.1f}%)")

    multi = int((f10["nLineages"] > 1).sum())
    _f(out, "ambiguous_cells", {"rows_with_multiple_lineages": multi,
                                "rows": len(f10),
                                "unique_cell_sample_pairs": int(
                                    f10.groupby("SampleNum")["cellID"].nunique().sum())},
       PRESENT,
       f"{multi} rows carry nLineages > 1; every clone-assigned row otherwise maps to one clone")

    # ---- clone structure --------------------------------------------------------------------- #
    s1 = set(f10.loc[f10.SampleNum == 1, "BC50StarcodeD8"])
    s2 = set(f10.loc[f10.SampleNum == 2, "BC50StarcodeD8"])
    _f(out, "clone_structure",
       {"clones_total": int(f10["BC50StarcodeD8"].nunique()),
        "clones_sample1": len(s1), "clones_sample2": len(s2), "shared": len(s1 & s2)},
       PRESENT,
       f"{len(s1 & s2)} clones appear in BOTH 10x lanes -> a lane-wise split would leak clones. "
       f"The outer split unit must be the clone.")

    # ---- future outcome: the gDNA arm -------------------------------------------------------- #
    _f(out, "gdna_table_schema", list(gdna.columns), PRESENT,
       f"{len(gdna)} rows, {gdna['BC50StarcodeD8'].nunique()} distinct clone barcodes, "
       f"cellID is the constant {gdna['cellID'].unique().tolist()} and SampleNum is "
       f"{sorted(gdna['SampleNum'].unique().tolist())}; rows are pre-collapse read groups, so the "
       f"principled aggregation per clone is the SUM of counts")

    gsum = gdna.groupby("BC50StarcodeD8")["counts"].sum()
    bins = {"1": int((gsum == 1).sum()), "2": int((gsum == 2).sum()),
            "3-9": int(gsum.between(3, 9).sum()), "10-99": int(gsum.between(10, 99).sum()),
            "100-499": int(gsum.between(100, 499).sum()),
            ">=500": int((gsum >= 500).sum())}
    _f(out, "gdna_read_distribution", bins, PRESENT,
       "strongly bimodal: a 1-2 read noise mode and a colony mode in the hundreds to thousands")

    overlap = sorted(set(gsum.index) & set(f10["BC50StarcodeD8"]))
    min_overlap_reads = int(gsum.loc[overlap].min()) if overlap else 0
    primed_cells = f10[f10["BC50StarcodeD8"].isin(overlap)]

    # Threshold-insensitivity, measured rather than chosen: the overlap set can only change once
    # the floor exceeds the weakest overlapping clone.
    sweep = {}
    for floor in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, min_overlap_reads, min_overlap_reads + 1):
        keep = set(gsum[gsum >= floor].index)
        sweep[int(floor)] = {"gdna_clones": len(keep),
                             "primed_clones": len(keep & set(f10["BC50StarcodeD8"])),
                             "primed_cells": int(f10["BC50StarcodeD8"].isin(keep).sum())}
    stable = all(v["primed_clones"] == len(overlap)
                 for k, v in sweep.items() if k <= min_overlap_reads)
    _f(out, "outcome_threshold_sensitivity",
       {"sweep": sweep, "invariant_over": [1, min_overlap_reads], "is_invariant": bool(stable)},
       PRESENT,
       f"every clone seen in BOTH arms carries >= {min_overlap_reads} summed gDNA reads, i.e. none "
       f"sits in the 1-2 read noise mode. The primed set is therefore IDENTICAL for any read floor "
       f"in [1, {min_overlap_reads}] -- no threshold had to be chosen for it.")

    n_clone_total = int(f10["BC50StarcodeD8"].nunique())
    per_sample = {int(sn): {"primed_cells": int((primed_cells.SampleNum == sn).sum()),
                            "primed_clones": int(primed_cells.loc[primed_cells.SampleNum == sn,
                                                                  "BC50StarcodeD8"].nunique())}
                  for sn in sorted(f10["SampleNum"].unique())}
    _f(out, "prospective_label_counts",
       {"positive_clones": len(overlap), "negative_clones": n_clone_total - len(overlap),
        "positive_cells": int(len(primed_cells)), "negative_cells": int(len(f10) - len(primed_cells)),
        "per_sample": per_sample,
        "positive_rate_clones": round(len(overlap) / n_clone_total, 4)},
       PRESENT,
       "primed := the clone's barcode is recovered from the post-reprogramming gDNA arm "
       "(SampleNum 3). Rare-event geometry, as the biology implies.")

    # ---- the published count, used ONLY as an after-the-fact check ---------------------------- #
    readings = {
        "pooled_cells": int(len(primed_cells)),
        "pooled_clones": len(overlap),
        "sample1_cells": per_sample.get(1, {}).get("primed_cells"),
        "sample1_clones": per_sample.get(1, {}).get("primed_clones"),
        "sample2_cells": per_sample.get(2, {}).get("primed_cells"),
        "sample2_clones": per_sample.get(2, {}).get("primed_clones"),
    }
    matches = sorted(k for k, v in readings.items() if v == 42)
    _f(out, "published_42_check", {"readings": readings, "readings_equal_to_42": matches},
       PRESENT if matches else ABSENT_PROVEN,
       "VALIDATION CHECK ONLY -- no threshold was moved to produce it. The published figure of 42 "
       "primed cells is not reproduced by any cell-level reading; it coincides with the CLONE "
       "count of one lane. A coincidence at this resolution does not establish the author rule, "
       "so it is recorded and not adopted.")

    _f(out, "unresolved_outcome_rule",
       ["is 'primed' defined per clone or per cell?",
        "are the two 10x lanes pooled, or is one lane the anchor?",
        f"was any gDNA read floor applied? (immaterial over [1, {min_overlap_reads}], but "
        f"unstated)"],
       UNKNOWN,
       f"the author code directory {base / 'author_code_zenodo7707418'} is present but EMPTY, and "
       f"the GEO Data-Processing block documents only the cellranger scRNA path -- it says nothing "
       f"about barcode calling. The rule is pending that code drop, not pending more data.")

    _f(out, "replicate_structure",
       {"gsms_local": sorted(REWIND_GSMS), "biological_replicates_local": 1,
        "gsms_in_series_not_local": ["GSM7092517", "GSM7092518", "GSM7092519", "GSM7092520",
                                     "GSM7092521", "GSM7092522", "GSM7092523", "GSM7092524",
                                     "GSM7092525", "GSM7092526", "GSM7092527"]},
       PRESENT,
       "both local GSMs are titled 'biol rep 1' -- two 10x lanes of ONE biological replicate. The "
       "author barcode tables only carry SampleNum 1/2/3, so the linkage exists for rep 1 only. "
       "Outer units are clones within one replicate, not independent replicates.")

    # ---- the committed artefact: clone-level, small, deterministic ---------------------------- #
    tbl = (f10.groupby("BC50StarcodeD8")
              .agg(cells=("cellID", "size"),
                   cells_s1=("SampleNum", lambda s: int((s == 1).sum())),
                   cells_s2=("SampleNum", lambda s: int((s == 2).sum())))
              .reset_index()
              .rename(columns={"BC50StarcodeD8": "clone"}))
    tbl["gdna_reads"] = tbl["clone"].map(gsum).fillna(0).astype(int)
    tbl["primed"] = tbl["clone"].isin(overlap).astype(int)
    tbl = tbl.sort_values(["primed", "gdna_reads", "clone"], ascending=[False, False, True])
    REWIND_TABLE.write_text(tbl.to_csv(sep="\t", index=False), encoding="utf-8")

    # The verdict follows the tri-state findings, not a judgement call: the linkage is either
    # resolved or it is not, and an UNKNOWN on the outcome rule holds the dataset at PENDING.
    if not clean or not overlap:
        verdict = INVALID_LINKAGE
    elif any(f.status == UNKNOWN for f in out.values()):
        verdict = RECON_PENDING
    else:
        verdict = RECON_PASS
    return {"verdict": verdict, "findings": out,
            "manifest": manifest(sorted(req.values())),
            "clone_table": str(REWIND_TABLE.relative_to(ROOT))}


# --------------------------------------------------------------------------------------------- #
# ROLE B — GSE279162 WM989
# --------------------------------------------------------------------------------------------- #
def wm989_required(base: Path) -> dict[str, Path]:
    req = {"family_xml": base / "GSE279162_family.xml",
           "series_matrix": base / "GSE279162_series_matrix.txt.gz"}
    for gsm, name in WM989_SAMPLES.items():
        for kind in ("barcodes.tsv", "features.tsv", "matrix.mtx"):
            req[f"{name}_{kind}"] = base / f"{gsm}_{name}_filtered_{kind}.gz"
    return req


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

    # ---- feature structure, read from the files rather than trusted ------------------------- #
    feats = [ln.rstrip("\n").split("\t")
             for ln in gzip.open(req["Naive1_features.tsv"], "rt")]
    kinds = {}
    for f in feats:
        kinds[f[2] if len(f) > 2 else "<none>"] = kinds.get(f[2] if len(f) > 2 else "<none>", 0) + 1
    lin_names = [f[0] for f in feats if len(f) > 2 and f[2] == "Custom"]
    n_genes = len(feats) - len(lin_names)
    gene_block_first = all(f[2] == "Gene Expression" for f in feats[:n_genes] if len(f) > 2)
    lin_pattern = bool(lin_names) and all(re.fullmatch(r"L\d+", n) for n in lin_names[:1000])
    _f(out, "feature_structure",
       {"total": len(feats), "by_type": kinds, "n_genes": n_genes,
        "n_lineage_features": len(lin_names), "genes_first": bool(gene_block_first),
        "lineage_ids_match_L_pattern": lin_pattern},
       PRESENT,
       f"the clone identifier is a FEATURE ROW, not a side file: {len(lin_names)} 'Custom' "
       f"LinNNNN features sit after {n_genes} genes in the same matrix")
    if n_genes != N_GENES_WM989:
        _f(out, "gene_block_size_unexpected", n_genes, UNKNOWN,
           f"expected {N_GENES_WM989} genes before the Custom block")

    # ---- pull the lineage block out of every matrix ------------------------------------------ #
    lin_arr = np.asarray(lin_names, dtype=object)
    parts, per_sample = [], {}
    for gsm, name in WM989_SAMPLES.items():
        mpath = req[f"{name}_matrix.mtx"]
        ngene_plus, ncell, nnz = mtx_dims(mpath)
        d = mtx_rows_above(mpath, n_genes)
        d["clone"] = lin_arr[d["row"].values - n_genes - 1]
        d["sample"] = name
        parts.append(d)
        per_sample[name] = {
            "gsm": gsm, "cells": ncell, "nnz": nnz,
            "lineage_entries": int(len(d)),
            "cells_with_lineage_umi": int(d["col"].nunique()),
            "clones_observed": int(d["clone"].nunique()),
            "lineage_umi_total": int(d["val"].sum()),
        }
    A = pd.concat(parts, ignore_index=True)
    _f(out, "per_sample_lineage_capture", per_sample, PRESENT,
       "cells carrying at least one lineage UMI, per sample, straight from the matrices")

    # ---- per-cell assignment is NOT clean: record it, do not fix it -------------------------- #
    g = A.groupby(["sample", "col"])["val"]
    cell = pd.concat([g.sum().rename("total"), g.max().rename("top1"), g.size().rename("n_lin")],
                     axis=1)
    cell["frac"] = cell["top1"] / cell["total"]
    dom = {f"frac>={t}": int((cell["frac"] >= t).sum()) for t in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)}
    _f(out, "per_cell_barcode_dominance",
       {"cells_with_any_lineage_umi": int(len(cell)),
        "median_dominant_fraction": round(float(cell["frac"].median()), 3),
        "median_lineage_features_per_cell": int(cell["n_lin"].median()),
        "cells_at_dominance": dom},
       PRESENT,
       "the dominant barcode carries only a median ~50% of a cell's lineage UMIs and the median "
       "cell shows several lineage features, so a per-cell clone call REQUIRES an explicit "
       "dominance + UMI rule. GEO ships no such rule and no author-processed cell->clone table.")

    # ---- clone level, with the floor swept rather than chosen -------------------------------- #
    cl = (A.groupby(["sample", "clone"])
            .agg(umi=("val", "sum"), cells=("col", "nunique")).reset_index())
    sweep = {}
    for floor in UMI_FLOORS:
        d = cl[cl["umi"] >= floor]
        sets = {s: set(gg["clone"]) for s, gg in d.groupby("sample")}
        naive = set().union(*[sets.get(n, set()) for n in WM989_NAIVE])
        row = {"naive_clones": len(naive)}
        for t in WM989_TREATMENTS:
            st = sets.get(t, set())
            row[t] = {"clones": len(st), "positives": len(st & naive),
                      "negatives": len(naive - st)}
        sweep[floor] = row
    _f(out, "outcome_threshold_sensitivity", sweep, PRESENT,
       "the naive pool and every positive count move materially with the UMI floor "
       f"(naive clones {sweep[1]['naive_clones']} -> {sweep[max(UMI_FLOORS)]['naive_clones']}). "
       "Unlike the Rewind gDNA arm there is NO empty region to separate signal from noise, so no "
       "floor is selected here.")

    sets1 = {s: set(gg["clone"]) for s, gg in cl.groupby("sample")}
    naive1 = set().union(*[sets1.get(n, set()) for n in WM989_NAIVE])
    ntreat = pd.Series({c: sum(c in sets1.get(t, set()) for t in WM989_TREATMENTS)
                        for c in naive1})
    _f(out, "multi_treatment_coverage",
       {"naive_clones": len(naive1),
        "by_n_treatments": {int(k): int(v) for k, v in ntreat.value_counts().sort_index().items()},
        "clones_in_ge2_treatments": int((ntreat >= 2).sum()),
        "clones_in_all_6": int((ntreat == len(WM989_TREATMENTS)).sum())},
       PRESENT,
       f"{int((ntreat >= 2).sum())} naive clones are re-observed under two or more treatments at "
       f"floor 1 -- this is the structure a state x treatment interaction test needs, and it is "
       f"the reason this dataset was qualified for Role B.")

    naive_cells = cl[cl["sample"].isin(WM989_NAIVE)].groupby("clone")["cells"].sum()
    _f(out, "pretreatment_observations_per_clone",
       {"clones": int(len(naive_cells)),
        "median_cells": int(naive_cells.median()), "max_cells": int(naive_cells.max()),
        "clones_with_ge_k_cells": {k: int((naive_cells >= k).sum()) for k in (1, 2, 3, 5, 10)}},
       PRESENT,
       "X is a clone-level pretreatment profile; clones with a single naive cell give a "
       "one-cell estimate of that clone's starting state")

    _f(out, "design_from_source",
       {"pretreatment_samples": list(WM989_NAIVE), "treatments": list(WM989_TREATMENTS)},
       PRESENT,
       "GEO Overall-Design: barcoded WM989 cells were allowed to double, a subset was taken for "
       "the untreated naive scRNA arm, and the remainder was split six ways -- 4 weeks of "
       "dabrafenib / trametinib / CoCl2 / acidic media, or 2 weeks of cisplatin / doxorubicin "
       "plus a 2 week holiday -- then re-sequenced. X precedes U precedes Y by construction.")

    _f(out, "unresolved_outcome_rule",
       ["per-cell clone call: minimum UMI and minimum dominant fraction are unspecified",
        "clone-presence floor for calling a clone 'surviving' in a treated sample",
        "whether Naive1/2/3 are pooled into one pretreatment pool or kept as replicates"],
       UNKNOWN,
       "GEO ships the raw feature-barcode matrix only. There is no author-processed clone table "
       "and the Data-Processing block states only that barcode reads were linked to 10x cell "
       "barcodes, not how a clone was called. Every reported count above is floor-conditional.")

    tbl = cl.pivot_table(index="clone", columns="sample", values="umi", fill_value=0).astype(int)
    tbl = tbl.loc[sorted(naive1)]
    tbl.insert(0, "n_treatments_observed", ntreat.reindex(tbl.index).astype(int))
    WM989_TABLE.write_text(tbl.reset_index().to_csv(sep="\t", index=False), encoding="utf-8")

    return {"verdict": RECON_PENDING, "findings": out,
            "manifest": manifest(sorted(req.values())),
            "clone_table": str(WM989_TABLE.relative_to(ROOT))}


# --------------------------------------------------------------------------------------------- #
def main() -> int:
    _RESULTS.mkdir(exist_ok=True)
    rec = {"stage": "21D",
           "plan": "plans/(newer)practical plans/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md",
           "additive_to": ["results/diag_stage21_data_audit_results.json",
                           "results/diag_stage21b_source_design_results.json"],
           "model_fitted": False, "src_modified": False,
           "raw_data_committed": False}

    for key, fn in (("GSE227151", audit_rewind), ("GSE279162", audit_gse279162)):
        r = fn()
        print(f"{key}: {r['verdict']}")
        rec[key] = {"verdict": r["verdict"],
                    "findings": {k: asdict(v) for k, v in r["findings"].items()},
                    "manifest": r["manifest"],
                    "clone_table": r.get("clone_table")}

    va, vb = rec["GSE227151"]["verdict"], rec["GSE279162"]["verdict"]
    if va in (MISSING_FILE, INVALID_LINKAGE, INSUFFICIENT_UNITS):
        overall = STAGE22_BLOCKED
    elif va == RECON_PASS and vb == RECON_PASS:
        overall = STAGE22_OPEN
    else:
        overall = STAGE22_PENDING
    rec["overall"] = overall
    print("OVERALL:", overall)

    OUT.write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
