"""STAGE 1.5.6 step 1b — WHERE do HFF's ΔAge years come from? Decompose the pipeline.

    python experiments/diag_pipeline_decompose.py "D:\\GSE242423\\GSE242423"

READ-ONLY. Writes `results/diag_pipeline_decompose_results.json`. `src/` untouched, no labels move.

WHY
---
Applying the clock DIRECTLY to HFF counts gives day-14 ΔAge = **-10.62 yr**. The pipeline's own
`y_age`, read from built shards by G-c step 1, gives **-24.02 yr**. A **13.4-year gap**, and nobody
has audited it.

That matters more than the sparse-clock work it interrupted: 1.5.6 found a **-14.10 yr bias inside
the clock**, and this gap is the same size and sits DOWNSTREAM of it. Fixing one while the other is
unmeasured would be arithmetic theatre.

THE CHAIN, FROM `build_dataset.py`
----------------------------------
Harmonization is OFF (no config sets `harmonize: true`), so for HFF the chain is exactly:

  S1  age        = clock.predict_age(normalize_counts(counts))
  S2  d_age_raw  = age - control_baseline(age, cell_line, is_control)     <- day-0 is the control
  S3  d          = d_age_raw - (a*cc + b)          `deconfound_age`, coef fit by OLS on TRAIN cells
  S4  y_age      = recenter_on_control_arrays(d)   re-subtract the control mean after deconfounding

**S1+S2 is what the direct run measured (-10.62).** So the 13.4 yr must be created by S3+S4 — the
cell-cycle deconfounder and the re-centring that follows it. This script measures each step.

ONE APPROXIMATION, STATED
-------------------------
The real `(a, b)` is fit on the primary regime's TRAIN cells across **all** datasets, and the train
split is not reconstructible without the built shards. Here it is fit on **HFF's own cells**. HFF is
**99.8 %** of the age-labelled corpus, so the global fit is dominated by these cells anyway — but
this is an approximation and the number it produces is therefore indicative, not exact. The SHAPE of
the contribution (which step moves ΔAge, and in which direction) does not depend on it.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

IPSC_DAY = 21.0
SHARD_DAY14 = -24.02          # results/diag_gc_hff_signature_results.json


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def by_day(day: np.ndarray, v: np.ndarray) -> dict:
    return {int(d): float(v[day == d].mean()) for d in sorted(np.unique(day))}


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    from cellfate.data.aging import (
        LinearClock,
        _control_baseline,
        recenter_on_control_arrays,
    )
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.proliferation import cell_cycle_score, deconfound_age, fit_deconfounder
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")

    hff = Path(sys.argv[1])
    genes_file = next(hff.glob("*genes.tsv.gz"))
    samples = []
    for mtx in sorted(hff.glob("*.matrix.mtx.gz")):
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            samples.append({"matrix": str(mtx), "barcodes": str(bc),
                            "label": mtx.name.split(".")[0].split("_")[-1]})
    src = GSE242423SingleCellSource(samples, str(genes_file), min_genes=hli.MIN_GENES,
                                    max_cells_per_sample=hli.MAX_CELLS,
                                    cells_per_run=hli.CELLS_PER_RUN)
    qc = QCConfig(max_mito_frac=0.20, min_genes=hli.MIN_GENES)
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")

    ages, ccs, days = [], [], []
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        norm = normalize_counts(raw.counts)
        ages.append(np.asarray(clock.predict_age(norm, raw.genes), float))
        ccs.append(np.asarray(cell_cycle_score(norm, raw.genes), float))
        days.append(raw.obs["time_h"].to_numpy(dtype=float) / 24.0)

    age = np.concatenate(ages)
    cc = np.concatenate(ccs)
    day = np.concatenate(days)
    keep = day != IPSC_DAY
    age, cc, day = age[keep], cc[keep], day[keep]
    line = np.array(["HFF"] * len(age))
    is_ctrl = day == 0.0
    print(f"\n[shape before statistic] {len(age)} cells, {len(np.unique(day))} timepoints")

    # ---- S1 -> S4, each step applied on top of the last ------------------------------------- #
    s1 = age
    s2 = age - _control_baseline(age, line, is_ctrl)
    coef = fit_deconfounder(s2, cc)
    s3 = deconfound_age(s2, cc, coef)
    s4 = recenter_on_control_arrays(s3, line, is_ctrl)

    steps = {"S1_clock_absolute_age": s1, "S2_control_relative": s2,
             "S3_deconfounded": s3, "S4_recentred_y_age": s4}
    out = {"script": "diag_pipeline_decompose", "utc": datetime.now(UTC).isoformat(),
           "n_cells": int(len(age)), "deconfounder_coef": {"a": coef[0], "b": coef[1]},
           "cc_by_day": by_day(day, cc), "shard_day14_reference": SHARD_DAY14, "steps": {}}

    print(f"\n  cell-cycle deconfounder fitted on HFF: dAge ~ {coef[0]:+.4f}*cc {coef[1]:+.4f}")
    print(f"\n  {'step':>26} {'day-0':>8} {'day-6':>8} {'day-14':>8}   what it adds at day 14")
    prev14 = None
    for name, v in steps.items():
        d = by_day(day, v)
        out["steps"][name] = d
        delta = "" if prev14 is None else f"{d[14] - prev14:+8.2f} yr"
        print(f"  {name:>26} {d[0]:8.2f} {d[6]:8.2f} {d[14]:8.2f}   {delta}")
        prev14 = d[14]

    d2, d4 = by_day(day, s2)[14], by_day(day, s4)[14]
    out["verdict"] = {
        "day14_after_S2_control_relative": d2,
        "day14_after_S4_pipeline_y_age": d4,
        "processing_contribution": float(d4 - d2),
        "gap_to_recorded_shards": float(d4 - SHARD_DAY14),
        "reproduces_shards": bool(abs(d4 - SHARD_DAY14) < 5.0),
    }
    print(f"\n  S2 (clock + control baseline) day-14 = {d2:+.2f}")
    print(f"  S4 (full pipeline y_age)      day-14 = {d4:+.2f}")
    print(f"  => PROCESSING CONTRIBUTES {d4 - d2:+.2f} yr at day 14")
    print(f"  recorded shards say {SHARD_DAY14:+.2f}; this reproduces it: "
          f"{out['verdict']['reproduces_shards']} (gap {d4 - SHARD_DAY14:+.2f})")
    (_RESULTS / "diag_pipeline_decompose_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print("  wrote results/diag_pipeline_decompose_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
