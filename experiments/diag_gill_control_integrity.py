"""
STAGE 1.5.6 step 3b — ROOT CAUSE: is a Gill control sample DEGENERATE in the raw GEO matrix?

    python experiments/diag_gill_control_integrity.py "D:/Gill"

READ-ONLY. Writes `results/diag_gill_control_integrity_results.json`. `src/` untouched, no build
touched, no label moved. Reads the raw GEO matrix only.

WHY THIS EXISTS
---------------
§5.5 eliminated T2 (the variance floor) and showed T1 (the mask) is not a carrier. §5.6 then showed
T3 (per-gene `sigma_gill`) survives, with the fold ordering reproduced exactly
(Spearman +1.000). Both stopped at *statistics over* `sigma_gill` without ever asking the prior
question:

> **Are the six Gill control samples that `sigma_gill` is estimated from actually sound?**

They are not. One is degenerate, and it is the one whose removal defines the anomalous fold.

WHAT IS CHECKED
---------------
For every sample column in `GSE165176_Log2_RPM_*`, on the RAW Log2 RPM values before any transform:

  * `min`, `median`, `mean`, `max` and the log2 DYNAMIC RANGE (`max - min`)
  * the `mean - min` gap: for real RNA-seq the mean sits well above the floor, because a minority
    of genes carry most of the signal. A column whose MEAN equals its FLOOR is nearly constant.
  * the implied linear library size after the pipeline's own inversion `2**x - 1`
    (`sources.py`: *"The matrix is Log2 RPM (already normalised). We invert it to linear RPM"*)

and the same for the normalised profiles, plus the rank correlation of each control to the others.

WHY IT MATTERS WHEREVER IT LANDS
---------------------------------
The day-0 `_Fib_` sample of each donor is `is_control` (`sources.py:417`), so it is BOTH:

  1. that donor's entire ΔAge **zero-point** -- Stage 1.5 audit §5.2's `n = 1` finding, and
  2. one of the **five or six** control samples the harmonizer estimates `sigma_gill` from,
     which sets the gain `sigma_gill / sigma_hff` applied to **HFF's** labels
     (`fit_harmonizer`, `build_dataset.py`).

So a bad control does not stay in its own donor. It reaches HFF -- 99.7 % of the age-labelled
corpus -- through the harmonizer, in **every fold that does not hold that donor out**.

READ:
  - no degenerate column -> the fold instability is ordinary estimator variance; §5.3's
    reconstruction proceeds as written.
  - a degenerate column, and it is a CONTROL -> the instability has a named data defect behind it,
    and the fold that EXCLUDES it is the one whose harmonizer is clean, not the one that is odd.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
REPO = Path(__file__).resolve().parents[1]
OUT = _RESULTS / "diag_gill_control_integrity_results.json"

# A column is FLAT if its mean sits essentially on its own floor. Real RNA-seq never does this:
# the mean is pulled well above the minimum by the highly-expressed minority.
FLAT_MEAN_MINUS_MIN = 0.05      # log2 units
NARROW_RANGE = 4.0              # log2 units of max-min; real columns span 12-17


def main() -> int:
    import pandas as pd
    from scipy.stats import spearmanr

    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    gill = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("D:/Gill")
    hits = glob.glob(str(gill / "*Log2_RPM*.txt.gz"))
    if not hits:
        print(f"FATAL: no *Log2_RPM*.txt.gz under {gill}")
        return 2
    src = hits[0]
    df = pd.read_csv(src, sep="\t", low_memory=False)
    cols = [c for c in df.columns[1:] if "_Sendai_" in c]

    print("\n" + "=" * 78)
    print("GILL CONTROL INTEGRITY — raw Log2 RPM, before any pipeline transform")
    print("=" * 78)
    print(f"source: {Path(src).name}\ncolumns scanned: {len(cols)}")

    stats = {}
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        lin = np.power(2.0, v) - 1.0
        lin[lin < 0] = 0.0
        stats[c] = {
            "min": float(v.min()), "median": float(np.median(v)),
            "mean": float(v.mean()), "max": float(v.max()),
            "log2_range": float(v.max() - v.min()),
            "mean_minus_min": float(v.mean() - v.min()),
            "implied_library_size": float(lin.sum()),
            "is_control": "_Fib_" in c,
        }

    degenerate = [c for c, s in stats.items()
                  if s["mean_minus_min"] < FLAT_MEAN_MINUS_MIN or s["log2_range"] < NARROW_RANGE]
    ctrl = [c for c in stats if stats[c]["is_control"]]

    print("\n  the six CONTROL columns (day-0 `_Fib_`, `is_control` per sources.py:417)")
    print(render_table(
        ["column", "min", "median", "mean", "max", "range", "mean-min", "library"],
        [[c, f"{stats[c]['min']:.3f}", f"{stats[c]['median']:.3f}", f"{stats[c]['mean']:.3f}",
          f"{stats[c]['max']:.3f}", f"{stats[c]['log2_range']:.2f}",
          f"{stats[c]['mean_minus_min']:.3f}", f"{stats[c]['implied_library_size']:.3g}"]
         for c in sorted(ctrl)],
        aligns=["l"] + ["r"] * 7))

    print(f"\n  DEGENERATE columns over all {len(stats)} samples "
          f"(mean-min < {FLAT_MEAN_MINUS_MIN} or range < {NARROW_RANGE} log2):")
    if degenerate:
        for c in degenerate:
            s = stats[c]
            print(f"    {c}   mean-min {s['mean_minus_min']:.4f}   range {s['log2_range']:.3f}"
                  f"   {'** IS A CONTROL **' if s['is_control'] else '(treatment sample)'}")
    else:
        print("    none")

    # normalised-profile agreement between the controls
    order = sorted(ctrl)
    M = np.stack([np.power(2.0, pd.to_numeric(df[c], errors="coerce")
                           .to_numpy(float)) - 1.0 for c in order])
    M[~np.isfinite(M)] = 0.0
    M[M < 0] = 0.0
    M = M / np.maximum(M.sum(axis=1, keepdims=True), 1e-12) * 1e4
    M = np.log1p(M)
    S = spearmanr(M.T).correlation
    print("\n  rank agreement of each control profile with the other five (log1p-CP10k)")
    print(render_table(
        ["donor", "mean rho to others"],
        [[c.split("_")[0], f"{S[i, [j for j in range(len(order)) if j != i]].mean():.4f}"]
         for i, c in enumerate(order)], aligns=["l", "r"]))

    payload = {"script": "diag_gill_control_integrity", "source": Path(src).name,
               "n_columns": len(stats), "thresholds":
                   {"flat_mean_minus_min": FLAT_MEAN_MINUS_MIN, "narrow_range": NARROW_RANGE},
               "degenerate_columns": degenerate,
               "degenerate_includes_a_control": any(stats[c]["is_control"] for c in degenerate),
               "per_column": stats,
               "control_rank_agreement": {c.split("_")[0]: float(
                   S[i, [j for j in range(len(order)) if j != i]].mean())
                   for i, c in enumerate(order)}}
    OUT.write_text(json.dumps(payload, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    if payload["degenerate_includes_a_control"]:
        print("   VERDICT: a CONTROL sample is degenerate.")
        print("   That sample is both its donor's entire ΔAge zero-point AND one of the five or")
        print("   six controls sigma_gill is estimated from, so it reaches HFF's labels through")
        print("   the harmonizer in every fold that does NOT hold its donor out.")
        print("   The fold that EXCLUDES it is the one with the clean harmonizer.")
    else:
        print("   VERDICT: no control is degenerate; the fold instability is estimator variance.")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
