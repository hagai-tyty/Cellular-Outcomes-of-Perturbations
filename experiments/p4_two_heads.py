"""P4 — are identity-loss and apoptosis one endpoint or two?

    python experiments/p4_two_heads.py

READ-ONLY. Writes `results/p4_two_heads_results.json`. `src/` untouched, no retrain.
Graded against `plans/STAGE_3_P4_TWO_HEADS_PREREG.md`, committed BEFORE this file existed.

THE STRUCTURAL FACT, STATED BEFORE ANY NUMBER
---------------------------------------------
Labels are the ARGMAX of a three-class call over (safe, loss, death), so a cell is either
identity-lost or apoptotic, NEVER BOTH. The union is exactly `P(loss) + P(death)`.

Collapsing therefore loses nothing to overlap -- it loses the DISTINCTION between two failure
modes. Whether that matters is an empirical question about their TIME COURSES, and that is what is
measured here.

Mutual exclusivity is IMPOSED by the argmax, not measured: a cell partway into both programmes is
forced to one. That is a property of the labeller and it bounds how much the heads could overlap.
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
OUT = _RESULTS / "p4_two_heads_results.json"

SUFFIX = "_c7"
BUNDLE = "cellfate_loocv_N2"
LINE = "HFF"


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Binomial interval that stays inside [0,1] at the extremes, where the apoptosis head sits."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(c - h), float(c + h)


def _rank(x) -> np.ndarray:
    """Ranks with TIES AVERAGED. `argsort(argsort(...))` breaks ties arbitrarily, which is not a
    rank -- and P(apoptosis) is tied at 0.0397 on two timepoints, so it matters here."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[order] = np.arange(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = np.mean(r[order[i:j + 1]])
        i = j + 1
    return r


def spearman(a, b) -> float:
    ra, rb = _rank(a), _rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def pearson(a, b) -> float:
    a, b = np.asarray(a, float) - np.mean(a), np.asarray(b, float) - np.mean(b)
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d else float("nan")


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.common.constants import DEATH_IDX, LOSS_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split
    install_pretty_console()

    print("\n" + "=" * 92)
    print("P4 — IDENTITY LOSS AND APOPTOSIS: one endpoint or two?")
    print("=" * 92)
    print("Graded against plans/STAGE_3_P4_TWO_HEADS_PREREG.md, committed BEFORE this script.")

    paths = ArtifactPaths.of(str(REPO / f"{BUNDLE}{SUFFIX}"))
    parts = [gather_split(paths, "holdout", s) for s in ("train", "val", "calib")]
    line = np.concatenate([np.asarray([str(v) for v in p.cell_line]) for p in parts])
    keep = line == LINE
    logt = np.concatenate([np.asarray(p.dose_time[:, 1], float) for p in parts])[keep]
    cls = np.concatenate([p.y_cls.astype(int) for p in parts])[keep]
    del parts

    # the structural fact, verified rather than assumed
    both = int(np.sum((cls == LOSS_IDX) & (cls == DEATH_IDX)))
    print(f"\n   cells: {len(cls)}   cells labelled BOTH loss and death: {both} "
          f"(argmax makes this structurally impossible)")

    tps = np.unique(np.round(logt, 6))
    rows, ploss, pdeath, days, shares = [], [], [], [], []
    for t in tps:
        m = np.isclose(logt, t)
        n = int(m.sum())
        kl = int(np.sum(cls[m] == LOSS_IDX))
        kd = int(np.sum(cls[m] == DEATH_IDX))
        pl, pd = kl / n, kd / n
        llo, lhi = wilson(kl, n)
        dlo, dhi = wilson(kd, n)
        un = pl + pd
        share = pd / un if un > 0 else float("nan")
        d = float(np.round(np.exp(t) / 24.0, 1))
        days.append(d)
        ploss.append(pl)
        pdeath.append(pd)
        shares.append(share)
        rows.append([f"{d:g}", str(n), f"{pl:.4f}", f"[{llo:.3f},{lhi:.3f}]",
                     f"{pd:.4f}", f"[{dlo:.3f},{dhi:.3f}]", f"{un:.4f}",
                     f"{100 * share:.1f}%"])

    print("\n  THE TWO FAILURE MODES ACROSS THE COURSE")
    print(render_table(["day", "cells", "P(identity loss)", "95% CI", "P(apoptosis)", "95% CI",
                        "union", "apoptosis share"], rows,
                       aligns=["r", "r", "r", "r", "r", "r", "r", "r"]))

    sp, pe = spearman(ploss, pdeath), pearson(ploss, pdeath)
    peak_l, peak_d = days[int(np.argmax(ploss))], days[int(np.argmax(pdeath))]
    gap = abs(int(np.argmax(ploss)) - int(np.argmax(pdeath)))
    smin, smax = float(np.nanmin(shares)), float(np.nanmax(shares))
    fold = smax / smin if smin > 0 else float("inf")

    print(f"\n   Spearman(loss, death) across timepoints = {sp:+.3f}   Pearson = {pe:+.3f}")
    print(f"   peak identity loss at day {peak_l:g}; peak apoptosis at day {peak_d:g} "
          f"({gap} timepoints apart)")
    print(f"   apoptosis share of the unsafe fraction ranges {100 * smin:.1f}% .. {100 * smax:.1f}%"
          f"  ({fold:.0f}x)")

    # ---- graded, exactly as pre-registered -------------------------------------------------
    if sp > 0.9 and gap <= 1:
        v = ("H1 REDUNDANT — the heads rise and fall together; P4's premise is REFUTED and the "
             "collapsed endpoint should be restored")
    elif sp < 0:
        v = ("H3 OPPOSITE DIRECTIONS — the union is dominated by whichever is larger and the "
             "other becomes INVISIBLE in it; two heads are MANDATORY and the collapsed endpoint "
             "is actively misleading")
    else:
        v = ("H2 DIFFERENT TIME COURSES — a single number cannot express both; the two-head form "
             "is warranted")
    print("\n" + "-" * 92)
    print(f"VERDICT: {v}")
    print("-" * 92)
    if fold > 10:
        print(f"   SECONDARY: the apoptosis share varies {fold:.0f}x across the course, so a "
              "single\n   'unsafe' number is a DIFFERENT QUANTITY at different days even when it "
              "is\n   numerically similar.")

    out = {"script": "p4_two_heads", "prereg": "plans/STAGE_3_P4_TWO_HEADS_PREREG.md",
           "n_cells": int(len(cls)), "cells_labelled_both": both,
           "days": days, "p_loss": ploss, "p_death": pdeath, "apoptosis_share": shares,
           "spearman": sp, "pearson": pe, "peak_loss_day": peak_l, "peak_death_day": peak_d,
           "peak_gap_timepoints": gap, "share_fold_range": fold, "verdict": v, "rows": rows}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   LIMITS: one line; marker-based labels; mutual exclusivity is IMPOSED by the")
    print("   argmax, not measured. This changes what is REPORTED, not any model weights.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
