"""
Harmonization diagnostic: split ΔAge into CONTROL vs REPROGRAMMED per dataset,
and breakdown by phase-dependent timepoint.

USAGE (repo root, venv active):
    python diag_harmonize.py                 # defaults to dataset dir 'cellfate_multi'
    python diag_harmonize.py cellfate_multi
"""
from __future__ import annotations

import sys
from collections import defaultdict

import numpy as np

from cellfate.common import io
from cellfate.common.io import ArtifactPaths

ROOT = sys.argv[1] if len(sys.argv) > 1 else "cellfate_multi"


def dataset_of(cell_line: str) -> str:
    return "hff_sc" if cell_line.upper() == "HFF" else "gill_bulk"


def main() -> None:
    # Force standard output to support UTF-8 characters like 'Δ' in PowerShell
    sys.stdout.reconfigure(encoding='utf-8')
    
    paths = ArtifactPaths.of(ROOT)
    man = io.manifest_rows(io.load_manifest(paths))

    # cell_id -> (cell_line, pert_id, timepoint)
    # Safely extract timepoint from manifest (fallback to '?' if missing)
    meta = {}
    for r in man:
        tp = getattr(r, 'dose_time', getattr(r, 'timepoint', getattr(r, 'time', '?')))
        meta[r.cell_id] = (r.cell_line, r.pert_id, str(tp))

    # cell_id -> ΔAge (age-valid only)
    age = {}
    timepoints_from_shard = {}
    import glob
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        
        # Check for the correct timepoint column name
        has_tp_col = "timepoint" in a
        has_time_col = "time" in a
        has_dose_col = "dose_time" in a
        
        for i in range(len(a["cell_id"])):
            if bool(a["age_mask"][i]):
                cid = a["cell_id"][i]
                age[cid] = float(a["y_age"][i])
                
                # Extract timepoint from shard if available
                if has_dose_col:
                    timepoints_from_shard[cid] = str(a["dose_time"][i])
                elif has_tp_col:
                    timepoints_from_shard[cid] = str(a["timepoint"][i])
                elif has_time_col:
                    timepoints_from_shard[cid] = str(a["time"][i])

    # aggregate by (dataset, cell_line, is_control) and by (dataset, cell_line, timepoint)
    by_line = defaultdict(lambda: {"ctrl": [], "reprog": []})
    by_ds = defaultdict(lambda: {"ctrl": [], "reprog": []})
    by_time = defaultdict(list)

    for cid, dage in age.items():
        cl, pert, manifest_tp = meta.get(cid, ("?", "?", "?"))
        
        # Parquet timepoint takes precedence if available
        tp = timepoints_from_shard.get(cid, manifest_tp)
        
        bucket = "ctrl" if pert == "control" else "reprog"
        ds = dataset_of(cl)
        
        by_line[cl][bucket].append(dage)
        by_ds[ds][bucket].append(dage)
        by_time[(ds, cl, tp)].append(dage)

    def fmt(vals):
        return f"mean={np.mean(vals):+7.2f}  n={len(vals):>5}" if vals else "   (none)      "

    print(f"\n=== ΔAge by cell line: CONTROL vs REPROGRAMMED  ({ROOT}) ===")
    print(f"{'cell line':<10}{'dataset':<11}{'CONTROLS':<26}{'REPROGRAMMED':<26}")
    for cl in sorted(by_line):
        c, r = by_line[cl]["ctrl"], by_line[cl]["reprog"]
        print(f"{cl:<10}{dataset_of(cl):<11}{fmt(c):<26}{fmt(r):<26}")

    # --- NEW: Phase-Dependent Timepoint Diagnostic ---
    print("\n=== BY-TIMEPOINT REJUVENATION TRAJECTORY (PHASE-DEPENDENT) ===")
    print(f"{'dataset':<12}{'cell line':<11}{'timepoint':<12}{'ΔAge':<20}")

    # Custom sort key to sort 'D2', 'D14' chronologically instead of alphabetically, with 'iPSC' at the end
    def tp_sort_key(k):
        ds, cl, tp = k
        try:
            val = float(tp.replace('D', '').replace('Day', '').strip())
        except ValueError:
            val = float('inf') 
        return (ds, cl, val, tp)

    for (ds, cl, tp) in sorted(by_time.keys(), key=tp_sort_key):
        vals = by_time[(ds, cl, tp)]
        print(f"{ds:<12}{cl:<11}{tp:<12}{fmt(vals)}")

    # --- Original A1/A3 Pooled Checks ---
    print("\n=== ΔAge pooled by DATASET (the acceptance test) ===")
    print(f"{'dataset':<12}{'CONTROLS (A1: want ~0)':<30}{'REPROGRAMMED':<26}")
    for ds in sorted(by_ds):
        c, r = by_ds[ds]["ctrl"], by_ds[ds]["reprog"]
        print(f"{ds:<12}{fmt(c):<30}{fmt(r):<26}")

    hff_c = by_ds["hff_sc"]["ctrl"]
    gill_c = by_ds["gill_bulk"]["ctrl"]
    if hff_c and gill_c:
        gap = abs(np.mean(hff_c) - np.mean(gill_c))
        print(f"\n[A1] HFF control mean   = {np.mean(hff_c):+.2f}")
        print(f"[A1] Gill control mean  = {np.mean(gill_c):+.2f}")
        print(f"[A3] control scale gap  = {gap:.2f}  -> {'PASS (<~3)' if gap < 3 else 'FAIL (gap remains)'}")
        print("\nRead: if BOTH control means are ~0 and the gap is small, harmonization worked and")
        print("the sanity-table +14..+32 is reprogramming signal (sign is a separate question).")
        print("If the Gill CONTROL mean is itself large, a residual batch gap remains.")


if __name__ == "__main__":
    main()