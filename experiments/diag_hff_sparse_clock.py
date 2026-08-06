"""STAGE 1.5.6 step 1 — THE FALSIFIER. Does the sparse clock move HFF the predicted amount?

    python experiments/diag_hff_sparse_clock.py "D:\\GSE242423\\GSE242423"

READ-ONLY. Writes `results/diag_hff_sparse_clock_results.json`. `src/` untouched, no labels move.

THE PREDICTION, MADE BEFORE THIS RAN
------------------------------------
On Gill, restricting the Fleischer clock to its ~100 largest-|weight| genes removed a **-14.10 yr
systematic bias** in ΔAge (MAE 16.61 -> 5.36, validated leave-one-donor-out).

That bias is a property of the CLOCK's dense weights, not of any dataset -- so it is applied to
HFF's 33,613 labels too. G-c step 1 measured HFF's trajectory reaching **-24.0 yr at day 14** under
the full clock. If the same artefact is at work:

    PREDICTED   day-14 ΔAge moves from -24.0 to roughly -10 yr
    PREDICTED   the trajectory SHAPE survives -- rho(day, ΔAge) stays near -0.9

**If the shape collapses, the hypothesis is wrong and Stage 1.5.6 stops here.** That is the point of
running this first: it is the cheapest thing that can falsify the whole idea, and it needs no
retrain and no new data.

WHY THE SHAPE MATTERS MORE THAN THE MAGNITUDE
---------------------------------------------
A bias correction should shift a trajectory, not destroy it. If sparsifying merely deleted signal,
rho would fall towards zero and the day-ordering would break. If it removed an offset, rho survives
and the curve slides up. Those two outcomes look nothing alike, which is what makes this a test
rather than a demonstration.
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

KS = (50, 100, 150, None)          # None = the full clock, the incumbent
IPSC_DAY = 21.0
PRED_DAY14 = -10.0                 # the pre-registered prediction
PRED_RHO = -0.90


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    b0 = float(clock.get("intercept", 0.0))

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

    days, mats, gene_ref = [], [], None
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        gene_ref = raw.genes
        mats.append(normalize_counts(raw.counts))
        days.append(raw.obs["time_h"].to_numpy(dtype=float) / 24.0)
    day = np.concatenate(days)
    keep = day != IPSC_DAY
    day = day[keep]
    print(f"\n[shape before statistic] {int(keep.sum())} cells, {len(np.unique(day))} timepoints")

    wv = np.array([W.get(g, 0.0) for g in gene_ref])
    ordr = np.argsort(-np.abs(wv))
    out = {"script": "diag_hff_sparse_clock", "utc": datetime.now(UTC).isoformat(),
           "n_cells": int(keep.sum()),
           "prediction": {"day14_dage": PRED_DAY14, "rho": PRED_RHO,
                          "source": "Gill bias -14.10 yr removed at k~100; G-c step 1 measured -24.0"},
           "k": {}}

    for k in KS:
        sub = wv.copy()
        if k is not None:
            sub = np.zeros_like(wv)
            sub[ordr[:k]] = wv[ordr[:k]]
        age = np.concatenate([m @ sub for m in mats])[keep] + b0
        base = age[day == 0.0].mean()
        dage = age - base
        uk = sorted(np.unique(day))
        means = [float(dage[day == d].mean()) for d in uk]
        rho = hli.spearman(np.array(uk, float), np.array(means))
        out["k"][str(k)] = {"days": [float(d) for d in uk], "mean_dage": means,
                            "rho_day": rho, "day14": means[-1]}
        lab = "all 33,155" if k is None else f"top{k}"
        print(f"\n  {lab:>12}  rho(day, dAge) = {rho:+.3f}   day-14 dAge = {means[-1]:+7.2f} yr")
        print("               " + "  ".join(f"d{int(d)}:{m:+6.1f}" for d, m in zip(uk, means, strict=True)))

    full = out["k"]["None"]
    best = out["k"]["100"]
    shape_ok = abs(best["rho_day"]) >= 0.70
    moved = best["day14"] - full["day14"]
    # The prediction was "-24.0 -> ~-10". It is only meaningful if THIS run reproduces the -24.0
    # baseline it was predicated on. G-c step 1 measured -24.0 from BUILT SHARDS (pipeline `y_age`,
    # which adds harmonization, cell-cycle deconfounding and control re-centring); this script
    # applies the clock DIRECTLY to counts and does none of that. Checking the baseline first is
    # what stops "landed near the predicted number" being read as a confirmation when the thing it
    # was supposed to have moved FROM never held.
    RECORDED_FULL_DAY14 = -24.02        # results/diag_gc_hff_signature_results.json
    baseline_reproduced = abs(full["day14"] - RECORDED_FULL_DAY14) < 5.0
    out["verdict"] = {
        "shape_survived": bool(shape_ok),
        "day14_full": full["day14"], "day14_top100": best["day14"], "shift": float(moved),
        "recorded_full_day14": RECORDED_FULL_DAY14,
        "baseline_reproduced": bool(baseline_reproduced),
        "baseline_gap": float(full["day14"] - RECORDED_FULL_DAY14),
        "status": ("BASELINE_NOT_REPRODUCED" if not baseline_reproduced
                   else ("SHAPE_LOST" if not shape_ok
                         else ("CONFIRMED" if abs(best["day14"] - PRED_DAY14) < 6.0
                               else "SHIFTED_BUT_OFF_PREDICTION"))),
    }
    if not baseline_reproduced:
        print(f"\n  !! this run's FULL-clock day-14 is {full['day14']:+.2f}, not the recorded "
              f"{RECORDED_FULL_DAY14:+.2f} ({full['day14'] - RECORDED_FULL_DAY14:+.2f} yr gap).")
        print("     The prediction assumed the recorded value, so it cannot be confirmed here.")
        print("     The gap is the PIPELINE's processing, not the clock -- see plan sec 5.1.")
    print(f"\n  PREDICTED day-14 ~ {PRED_DAY14:+.1f} with shape preserved (|rho| >= 0.70)")
    print(f"  ACTUAL    day-14 = {best['day14']:+.2f}   rho = {best['rho_day']:+.3f}")
    print(f"  => {out['verdict']['status']}")
    (_RESULTS / "diag_hff_sparse_clock_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print("  wrote results/diag_hff_sparse_clock_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
