"""
Harmonization diagnostic: split ΔAge into CONTROL vs REPROGRAMMED per dataset,
and show the ΔAge trajectory BY TIMEPOINT (to tell phase-dependent biology from noise).

The timepoint is recovered from each cell's dose_time = [log10(dose), log(time_h)]
in the shard: day = exp(dose_time[1]) / 24.

USAGE (repo root, venv active):
    python diag_harmonize.py                 # defaults to 'cellfate_multi'
    python diag_harmonize.py cellfate_multi
"""
from __future__ import annotations

import glob
import sys
from collections import defaultdict

import numpy as np

from cellfate.common import io
from cellfate.common.io import ArtifactPaths

try:
    sys.stdout.reconfigure(encoding="utf-8")   # so 'Δ' prints on Windows
except Exception:
    pass

ROOT = sys.argv[1] if len(sys.argv) > 1 else "cellfate_multi"


def dataset_of(cell_line: str) -> str:
    return "hff_sc" if cell_line.upper() == "HFF" else "gill_bulk"


def day_of(dose_time_row) -> float:
    """Decode timepoint (days) from dose_time = [log10 dose, log time_h]."""
    t_h = float(np.exp(np.asarray(dose_time_row)[1]))
    return round(t_h / 24.0, 1)


def fmt(vals):
    return f"mean={np.mean(vals):+7.2f}  n={len(vals):>5}" if vals else "   (none)"


def main() -> None:
    paths = ArtifactPaths.of(ROOT)
    meta = {r.cell_id: (r.cell_line, r.pert_id) for r in io.manifest_rows(io.load_manifest(paths))}

    by_line = defaultdict(lambda: {"ctrl": [], "reprog": []})
    by_ds = defaultdict(lambda: {"ctrl": [], "reprog": []})
    by_time = defaultdict(list)   # (dataset, cell_line, day) -> [ΔAge]

    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i in range(len(a["cell_id"])):
            if not bool(a["age_mask"][i]):
                continue
            cid = a["cell_id"][i]
            dage = float(a["y_age"][i])
            cl, pert = meta.get(cid, ("?", "?"))
            ds = dataset_of(cl)
            bucket = "ctrl" if pert == "control" else "reprog"
            by_line[cl][bucket].append(dage)
            by_ds[ds][bucket].append(dage)
            by_time[(ds, cl, day_of(a["dose_time"][i]))].append(dage)

    print(f"\n=== ΔAge by cell line: CONTROL vs REPROGRAMMED  ({ROOT}) ===")
    print(f"{'cell line':<10}{'dataset':<11}{'CONTROLS':<26}{'REPROGRAMMED':<26}")
    for cl in sorted(by_line):
        c, r = by_line[cl]["ctrl"], by_line[cl]["reprog"]
        print(f"{cl:<10}{dataset_of(cl):<11}{fmt(c):<26}{fmt(r):<26}")

    print("\n=== ΔAge BY TIMEPOINT (does the sign track reprogramming phase?) ===")
    print(f"{'dataset':<11}{'cell line':<10}{'day':>6}   ΔAge")
    for key in sorted(by_time, key=lambda k: (k[0], k[1], k[2])):
        ds, cl, day = key
        print(f"{ds:<11}{cl:<10}{day:>6}   {fmt(by_time[key])}")

    print("\n=== pooled by DATASET ===")
    print(f"{'dataset':<12}{'CONTROLS (want ~0)':<28}{'REPROGRAMMED':<26}")
    for ds in sorted(by_ds):
        c, r = by_ds[ds]["ctrl"], by_ds[ds]["reprog"]
        print(f"{ds:<12}{fmt(c):<28}{fmt(r):<26}")

    hff_c, gill_c = by_ds["hff_sc"]["ctrl"], by_ds["gill_bulk"]["ctrl"]
    if hff_c and gill_c:
        gap = abs(np.mean(hff_c) - np.mean(gill_c))
        print(f"\n[controls] HFF {np.mean(hff_c):+.2f} | Gill {np.mean(gill_c):+.2f} | "
              f"gap {gap:.2f} -> {'PASS' if gap < 3 else 'FAIL'} (controls ~0 by construction)")
        print("Read the BY-TIMEPOINT table: a smooth sign that tracks day = phase-dependent biology;")
        print("signs jumping around randomly across days = noise / clock artifact.")


if __name__ == "__main__":
    main()
