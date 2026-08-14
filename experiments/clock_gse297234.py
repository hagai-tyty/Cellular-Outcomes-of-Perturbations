"""CLOCK ON GSE297234 — can our clock order two donors 74 years apart?

    python experiments/clock_gse297234.py

Graded against `plans/STAGE_1_5_8_CLOCK_ON_GSE297234_PREREG.md`, committed BEFORE this file
existed. READ-ONLY with respect to the repo: no build, no retrain, `src/` untouched.

WHY
---
M-E1 found the clock reads fresh 38-53 yr fibroblasts as 72-82 -- a +30 yr floor. Our own data
cannot diagnose it: 3 donors spanning 15 years gives P(correct ordering) = 0.755, so nothing there
is testable. `GSE297234` (Lu et al. 2025) has two UNTREATED adult donors 74 years apart --
GM23815 (22) and GM00731 (96) -- same lab, same protocol, 10x. At that gap the same error model
gives P = 0.9997.

THE DISCRIMINATING FACT
-----------------------
Fleischer's training lines are CORIELL repository stock (AG09599, AG04054, ... NIA/NIGMS).
GM23815 and GM00731 are ALSO Coriell. GSE165177's cells were purchased from LONZA. So if the +30
floor is driven by cell source/supplier, this dataset should show a SMALLER bias; if the bias is
unchanged, source is not the driver and Gill's ComBat-harmonisation route is the only one left.

CLOCKING
--------
The clock is bulk-trained, so the PRIMARY estimate is PSEUDOBULK -- sum raw counts across cells,
then the project's own path: normalize_counts(target_sum=1e4) -> predict_age. `normalize_counts`
applies CP10k AND log1p itself (`normalize.py:29`) -- wrapping it in a second `np.log1p` is a bug,
and it is the one that corrupted the first version of this run. A per-cell distribution is reported
as a clearly-labelled SECONDARY quantity, because per-cell dropout is a domain shift the clock was
never fitted for.
"""
from __future__ import annotations

import contextlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "clock_gse297234_results.json"

DATA = Path(r"D:\GSE297234")
SAMPLES = [
    {"gsm": "GSM8986590", "line": "GM23815", "age": 22.0,
     "file": "GSM8986590_GM23815_D0_filtered_feature_bc_matrix.h5"},
    {"gsm": "GSM8986586", "line": "GM00731", "age": 96.0,
     "file": "GSM8986586_GM00731_D0_filtered_feature_bc_matrix.h5"},
]
CLOCK = REPO / "configs" / "clocks" / "fleischer_clock.json"
CV_MAE = 12.26879346460328
GSE165177_DAY0_BIAS = 30.0          # measured in M-E1's bias decomposition
BOOTSTRAP = 200
SEED = 0


def load_10x_h5(path: Path):
    """Cell Ranger filtered_feature_bc_matrix.h5 -> (csc genes x cells, gene symbols)."""
    import h5py
    from scipy.sparse import csc_matrix

    with h5py.File(path, "r") as f:
        g = f["matrix"]
        data = np.asarray(g["data"], dtype=np.float64)
        indices = np.asarray(g["indices"])
        indptr = np.asarray(g["indptr"])
        shape = tuple(np.asarray(g["shape"]))          # (n_genes, n_cells)
        names = [s.decode() if isinstance(s, bytes) else str(s)
                 for s in np.asarray(g["features"]["name"])]
    return csc_matrix((data, indices, indptr), shape=shape), names


def dedup_highest(X, genes: list[str]):
    """Duplicate symbols: keep the highest-expressed row. The project's own rule."""
    tot = np.asarray(X.sum(axis=1)).ravel()
    best: dict[str, int] = {}
    for i, g in enumerate(genes):
        if g not in best or tot[i] > tot[best[g]]:
            best[g] = i
    keep = np.array(sorted(best.values()))
    return X[keep], [genes[i] for i in keep]


def pseudobulk_age(X, genes, clock) -> float:
    """Sum counts across cells, then the project's normalisation path."""
    from cellfate.data.normalize import normalize_counts

    counts = np.asarray(X.sum(axis=1)).ravel()[None, :]
    expr = normalize_counts(counts, target_sum=1e4)          # applies CP10k AND log1p itself
    return float(clock.predict_age(expr, genes)[0])


def per_cell_ages(X, genes, clock) -> np.ndarray:
    """age per cell, without densifying: log1p(CP10k) keeps the sparsity pattern."""
    from scipy.sparse import diags
    lib = np.asarray(X.sum(axis=0)).ravel()
    lib[lib == 0] = 1.0
    Y = (X @ diags(1e4 / lib)).tocsr()
    Y.data = np.log1p(Y.data)
    w = np.array([clock.weights.get(g, 0.0) for g in genes], dtype=np.float64)
    return np.asarray(Y.T @ w).ravel() + clock.intercept


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    install_pretty_console()

    print("\n" + "=" * 92)
    print("CLOCK ON GSE297234 — two untreated donors 74 years apart (22 vs 96)")
    print("=" * 92)
    print("Graded against plans/STAGE_1_5_8_CLOCK_ON_GSE297234_PREREG.md, committed BEFORE")
    print("this script existed. Day 0 only — no reprogramming, no ΔAge, no safety target.")

    clock = LinearClock.from_json(CLOCK)
    rng = np.random.default_rng(SEED)
    out: dict = {"script": "clock_gse297234",
                 "prereg": "plans/STAGE_1_5_8_CLOCK_ON_GSE297234_PREREG.md",
                 "cv_mae": CV_MAE, "samples": []}

    rows, per_cell_rows = [], []
    for s in SAMPLES:
        X, genes = load_10x_h5(DATA / s["file"])
        n_genes_raw, n_cells = X.shape
        X, genes = dedup_highest(X, genes)
        age = pseudobulk_age(X, genes, clock)

        gset = set(genes)
        present = sum(1 for g in clock.weights if g in gset)
        w_tot = sum(abs(v) for v in clock.weights.values())
        w_pres = sum(abs(v) for g, v in clock.weights.items() if g in gset)

        # bootstrap over CELLS -- pseudobulk stability ONLY, NOT a donor-level interval
        boot = []
        for _ in range(BOOTSTRAP):
            idx = rng.integers(0, X.shape[1], X.shape[1])
            counts = np.asarray(X[:, idx].sum(axis=1)).ravel()[None, :]
            boot.append(float(clock.predict_age(
                normalize_counts(counts, target_sum=1e4), genes)[0]))
        blo, bhi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))

        pc = per_cell_ages(X, genes, clock)
        rows.append([s["line"], f"{s['age']:.0f}", str(n_cells), f"{age:.1f}",
                     f"[{blo:.1f},{bhi:.1f}]", f"{age - s['age']:+.1f}",
                     f"{100 * present / len(clock.weights):.1f}%",
                     f"{100 * w_pres / w_tot:.1f}%"])
        per_cell_rows.append([s["line"], f"{s['age']:.0f}", str(len(pc)), f"{pc.mean():.1f}",
                              f"{np.median(pc):.1f}", f"{pc.std():.1f}",
                              f"[{np.percentile(pc, 2.5):.1f},{np.percentile(pc, 97.5):.1f}]"])
        out["samples"].append({**s, "n_cells": int(n_cells), "n_genes_raw": int(n_genes_raw),
                               "n_genes_dedup": len(genes), "pseudobulk_age": age,
                               "boot_ci": [blo, bhi], "bias": age - s["age"],
                               "clock_genes_present_frac": present / len(clock.weights),
                               "clock_weight_mass_frac": w_pres / w_tot,
                               "per_cell_mean": float(pc.mean()),
                               "per_cell_median": float(np.median(pc)),
                               "per_cell_sd": float(pc.std())})
        del X

    print("\n  PSEUDOBULK — the primary estimate (the clock is bulk-trained)")
    print(render_table(["line", "true age", "cells", "predicted", "cell-boot 95%", "bias",
                        "clock genes", "weight mass"], rows,
                       aligns=["l", "r", "r", "r", "r", "r", "r", "r"]))
    print("   the bootstrap resamples CELLS -- it measures pseudobulk stability ONLY and is NOT")
    print("   a donor-level interval. With n=2 donors no donor-level interval exists.")

    print("\n  PER-CELL — secondary, and a domain the clock was never fitted for")
    print(render_table(["line", "true age", "cells", "mean", "median", "SD", "95% spread"],
                       per_cell_rows, aligns=["l", "r", "r", "r", "r", "r", "r"]))

    young, old = out["samples"][0], out["samples"][1]
    yp, op = young["pseudobulk_age"], old["pseudobulk_age"]
    gap = old["age"] - young["age"]

    # ---- N1 ordering -------------------------------------------------------------------------
    n1 = op > yp
    print("\n" + "-" * 92)
    print("PRE-REGISTERED OUTCOMES")
    print("-" * 92)
    print(f"N1 ORDERING at a {gap:.0f}-year gap: pred({old['age']:.0f}) = {op:.1f} vs "
          f"pred({young['age']:.0f}) = {yp:.1f}  ->  "
          f"{'CORRECT' if n1 else 'WRONG -- ESCALATE, this is not a noise story'}")

    # ---- N2 calibration ----------------------------------------------------------------------
    n2 = abs(yp - young["age"]) <= CV_MAE and abs(op - old["age"]) <= CV_MAE
    print(f"N2 CALIBRATION: |bias| = {abs(yp - young['age']):.1f} and "
          f"{abs(op - old['age']):.1f} against cv_mae {CV_MAE:.2f}  ->  "
          f"{'PASS-CALIBRATION' if n2 else 'FAIL-CALIBRATION'}")

    # ---- N3 the discriminating test ----------------------------------------------------------
    mean_bias = float(np.mean([young["bias"], old["bias"]]))
    if mean_bias <= 15.0:
        n3 = ("SOURCE MATTERS — Coriell-to-Coriell transfers far better than Coriell-to-Lonza; "
              "supplier becomes a first-class variable in the acquisition spec")
    elif abs(mean_bias - GSE165177_DAY0_BIAS) <= 10.0:
        n3 = ("SOURCE IS NOT THE DRIVER — a general cross-study/platform failure; "
              "ComBat-style harmonisation is the route, and supplier drops out of the spec")
    elif mean_bias >= 45.0:
        n3 = ("SINGLE-CELL ADDS ITS OWN PENALTY on top — pseudobulk-from-scRNA is a separate "
              "domain shift and must be handled as one")
    else:
        n3 = (f"BETWEEN the pre-registered bands (mean bias {mean_bias:+.1f} vs GSE165177's "
              f"+{GSE165177_DAY0_BIAS:.0f}) — report the number, claim neither branch")
    print(f"N3 SOURCE TEST: mean bias here {mean_bias:+.1f} yr vs GSE165177 day-0 "
          f"+{GSE165177_DAY0_BIAS:.0f} yr\n   -> {n3}")

    # ---- N4 slope ----------------------------------------------------------------------------
    slope = (op - yp) / gap
    n4 = ("tracks age at the right RATE — a pure intercept offset, fixable additively"
          if 0.7 <= slope <= 1.3 else
          "COMPRESSED dynamic range — a scale error as well as an offset; an intercept fix is "
          "not enough" if 0 < slope < 0.7 else
          "NO age signal on this data")
    print(f"N4 SLOPE: ({op:.1f} − {yp:.1f}) / {gap:.0f} = {slope:.3f}  ->  {n4}")
    print("   two points determine a slope exactly; this has NO error bar.")

    # ---- N5 coverage -------------------------------------------------------------------------
    print(f"N5 COVERAGE: genes {100 * young['clock_genes_present_frac']:.1f}% / weight mass "
          f"{100 * young['clock_weight_mass_frac']:.1f}%   vs GSE165177's 57.1% / 89.2%")

    # ---- expectations, graded ----------------------------------------------------------------
    exp = {
        "P-N1 clock orders 22 < 96 correctly": bool(n1),
        "P-N2 calibration FAILS with a positive bias": bool(not n2 and mean_bias > 0),
        "P-N3 bias smaller than +30 (Coriell hypothesis)": bool(mean_bias < GSE165177_DAY0_BIAS),
        "P-N4 slope well below 1": bool(slope < 0.7),
    }
    print("\n  PRE-REGISTERED EXPECTATIONS, GRADED")
    print(render_table(["expectation", "held?"],
                       [[k, "YES" if v else "NO"] for k, v in exp.items()], aligns=["l", "l"]))

    out["outcomes"] = {"N1_ordering_correct": bool(n1), "N2_calibration_pass": bool(n2),
                       "N3_verdict": n3, "N3_mean_bias": mean_bias,
                       "N4_slope": float(slope), "N4_verdict": n4}
    out["expectations"] = exp
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   LIMITS: 2 donors, so N1 is an OBSERVATION not a test and no donor-level interval")
    print("   exists; pseudobulk-from-scRNA is itself a domain shift; day 0 only.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
