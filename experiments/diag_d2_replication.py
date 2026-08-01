"""STAGE 1.5 §10 D2 — does the reprogramming "age rises" effect REPLICATE on independent data?

    python experiments/diag_d2_replication.py                    # defaults to D:\\GSE242423
    python experiments/diag_d2_replication.py "D:\\GSE242423"

READ-ONLY. Writes `diag_d2_replication_results.json`. Nothing is rebuilt or refitted; `src/` is not
touched.

WHY (STAGE_1_5 §9-§10)
----------------------
§9 settled that the clock is NOT broken: it reproduces age on its own domain (MAE 0.8 yr, rho 0.99)
and tracks in-range adult fibroblast age on Gill (+18 yr for a 21 yr gap). The one live finding is
E1b: during reprogramming (day <= 15) predicted age RISES (rho +0.205, CI [+0.009, +0.401]) -- weak,
but the wrong direction. H3 could not say why (`DIFFUSE`).

Before explaining that effect we check it is REAL. GSE242423 shares nothing with Gill except the
clock: different lab (Kundaje), different modality (single-cell, not bulk), different donor (HFF),
different protocol. If the rise reappears there it is a property of the clock on reprogramming
biology; if it does not, Gill's n=6 bulk design is implicated and the escalation weakens sharply.

DESIGN, PRE-COMMITTED (see the lab notebook, §10 D2)
  * 9 timepoints: D0 D2 D4 D6 D8 D10 D12 D14 iPSC.
  * PRIMARY WINDOW = D0-D14 (8 points), matching E1b's day <= 15. iPSC is EXCLUDED from the primary
    (a completed cell-type change -- the confound E1 already excluded) and reported as sensitivity.
  * Pseudobulk per timepoint through the production path (normalize_counts -> frozen clock), capped
    at 2000 cells/timepoint, seed 0. Pseudobulk (not per-cell) because the clock is bulk-fit and
    Gill is bulk -- like-for-like.
  * 5 disjoint pseudo-replicate pools per timepoint are reported as SPREAD only. They are NOT
    treated as independent samples for significance (that would be pseudo-replication).
  * Statistic: Spearman(predicted_age, day) over the 8 timepoint-level pseudobulks.

POWER, STATED UP FRONT. n = 8 timepoints, so |rho| ~ 0.74 is needed for p < 0.05. This test is
UNDERPOWERED for significance and is judged on DIRECTION, which is the honest bar for a replication
check and is fixed before the run:

    replicates  : rho > 0        (same direction as Gill's E1b +0.205)
    contradicts : rho <= -0.2
    ambiguous   : -0.2 < rho <= 0

WHAT THIS CANNOT DO. It anchors direction, not magnitude, on one donor at n=8. A positive rho shows
the effect is not a Gill artefact; it does NOT establish that ΔAge is invalid.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Pure logic — data-free, fully unit-tested; nothing below imports repo data.  #
# --------------------------------------------------------------------------- #
GILL_E1B_RHO = 0.205          # the effect being replicated (Stage 1.5 §8.6)
MAX_CELLS_PER_TIMEPOINT = 2000
N_PSEUDO_REPLICATES = 5
CONTRADICTS_AT = -0.2         # pre-committed boundary for "opposite direction"
IPSC_LABEL = "IPSC"


def day_of_label(label: str) -> float | None:
    """`D0`->0, `D14`->14, `iPSC`->None (no day; it is an endpoint, not a timepoint). Pure."""
    s = str(label).strip().upper()
    if s == IPSC_LABEL:
        return None
    if s.startswith("D") and s[1:].isdigit():
        return float(s[1:])
    return None


def spearman(x: list[float], y: list[float]) -> float:
    """Spearman rank correlation. NaN when undefined. Pure."""
    a, b = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).correlation)


def d2_verdict(rho: float, n_points: int, gill_rho: float = GILL_E1B_RHO) -> dict:
    """Direction-based replication verdict. Pure; no I/O.

    Judged on DIRECTION, not significance: n=8 timepoints cannot resolve a modest effect, and
    pretending otherwise would be the error this project keeps catching.
    """
    if not np.isfinite(rho) or n_points < 3:
        return {"status": "CANNOT_VERIFY", "n_points": int(n_points),
                "reason": f"only {n_points} usable timepoints; a trend needs >=3"}
    if rho > 0:
        s = "REPLICATES"
        r = (f"predicted age RISES with reprogramming day (rho {rho:+.3f} over {n_points} "
             f"timepoints), the same direction as Gill's E1b ({gill_rho:+.3f}) on fully independent "
             "data -> the effect is not a Gill design artefact")
    elif rho <= CONTRADICTS_AT:
        s = "CONTRADICTS"
        r = (f"predicted age FALLS with reprogramming day (rho {rho:+.3f}), the OPPOSITE of Gill's "
             f"E1b ({gill_rho:+.3f}) -> the Gill result is implicated as design-specific and the "
             "E1b escalation largely dissolves")
    else:
        s = "AMBIGUOUS"
        r = (f"rho {rho:+.3f} over {n_points} timepoints -- no clear replication either way at this "
             "sample size; underpowered, do not over-read in either direction")
    return {"status": s, "rho": float(rho), "n_points": int(n_points),
            "gill_e1b_rho": float(gill_rho), "reason": r}


def bars() -> list[dict]:
    """Pre-registered (ground rule §5b), with the power limitation stated rather than hidden."""
    return [{
        "id": "D2",
        "bar": f"Spearman(predicted_age, day) over D0-D14 is positive (replicates Gill's "
               f"E1b {GILL_E1B_RHO:+.3f}); <= {CONTRADICTS_AT} contradicts; between is ambiguous",
        "null": "the clock reads no consistent age trend across reprogramming (rho ~ 0)",
        "power_note": ("n=8 timepoints needs |rho| ~ 0.74 for p<0.05, so this is UNDERPOWERED for "
                       "significance and is judged on DIRECTION only -- fixed before the run"),
    }]


# --------------------------------------------------------------------------- #
# Real-data wiring (imports repo machinery only when actually run)            #
# --------------------------------------------------------------------------- #
def timepoint_ages(gse_dir: str) -> tuple[dict, dict]:
    """Pseudobulk predicted age per timepoint, via the production normalisation + frozen clock."""
    root = Path(__file__).resolve().parents[1]
    for p in (root, root / "local_runners", root / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    from run_multi_local import discover_gse  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GSE242423SingleCellSource

    samples, genes_file = discover_gse(gse_dir)
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")

    # One source per timepoint so each is loaded and pseudobulked independently. cells_per_run=None
    # keeps a timepoint whole; max_cells_per_sample caps memory and is pre-committed.
    out: dict[str, dict] = {}
    meta: dict = {"n_timepoints": len(samples), "max_cells": MAX_CELLS_PER_TIMEPOINT,
                  "n_pseudo_replicates": N_PSEUDO_REPLICATES, "errors": {}}
    for s in samples:
        label = str(s["label"])
        try:
            src = GSE242423SingleCellSource([s], genes_file, cell_line="HFF", min_genes=500,
                                            max_cells_per_sample=MAX_CELLS_PER_TIMEPOINT,
                                            cells_per_run=None, seed=0)
            chunk = src.plan()[0]
            raw = src.fetch(chunk)
            counts = np.asarray(raw.counts, dtype=np.float64)
            n = counts.shape[0]
            # timepoint-level pseudobulk = the primary unit
            pooled = counts.sum(axis=0, keepdims=True)
            age = float(clock.predict_age(normalize_counts(pooled), raw.genes)[0])
            # disjoint pseudo-replicates -> SPREAD only, never used for significance
            rng = np.random.default_rng(0)
            idx = rng.permutation(n)
            reps = [float(clock.predict_age(
                normalize_counts(counts[part].sum(axis=0, keepdims=True)), raw.genes)[0])
                for part in np.array_split(idx, N_PSEUDO_REPLICATES) if len(part)]
            out[label] = {"day": day_of_label(label), "n_cells": int(n),
                          "pseudobulk_age": age,
                          "replicate_ages": reps,
                          "replicate_sd": float(np.std(reps, ddof=1)) if len(reps) > 1 else None}
        except Exception as exc:  # noqa: BLE001 — recorded per timepoint, never aborts the scan
            meta["errors"][label] = repr(exc)[:160]
    return out, meta


def main() -> int:
    gse_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
    print("STAGE 1.5 §10 D2 — independent replication of the reprogramming age-rise (read-only)\n")
    print("  PRE-REGISTERED BAR (ground rule §5b):")
    for b in bars():
        print(f"    {b['id']}: {b['bar']}")
        print(f"        vs null: {b['null']}")
        print(f"        POWER  : {b['power_note']}")

    if not Path(gse_dir).exists():
        print(f"\n   !! {gse_dir} not found. Pass it: "
              'python experiments/diag_d2_replication.py "D:\\GSE242423"')
        return 1

    print(f"\n  loading {gse_dir} (streaming 10x MTX per timepoint; this takes a few minutes)...")
    tps, meta = timepoint_ages(gse_dir)
    if not tps:
        print("   !! no timepoints could be loaded.")
        return 1

    print(f"\n  {'timepoint':<11}{'day':>6}{'cells':>8}{'pseudobulk age':>16}{'rep SD':>9}")
    print("  " + "-" * 50)
    for lbl in sorted(tps, key=lambda k: (tps[k]["day"] is None, tps[k]["day"] or 0)):
        v = tps[lbl]
        d = "iPSC" if v["day"] is None else f"{v['day']:.0f}"
        sd = "n/a" if v["replicate_sd"] is None else f"{v['replicate_sd']:.2f}"
        print(f"  {lbl:<11}{d:>6}{v['n_cells']:>8}{v['pseudobulk_age']:>16.1f}{sd:>9}")
    for lbl, e in meta["errors"].items():
        print(f"  [!] {lbl}: {e}")

    # PRIMARY: D0-D14 only (matches E1b's day<=15 window); iPSC excluded as a cell-type change.
    prim = [(v["day"], v["pseudobulk_age"]) for v in tps.values() if v["day"] is not None]
    prim.sort()
    rho = spearman([d for d, _ in prim], [a for _, a in prim])
    verdict = d2_verdict(rho, len(prim))

    # SENSITIVITY: include iPSC, mapped past the last day, to show what the endpoint does.
    if any(v["day"] is None for v in tps.values()):
        last = max((d for d, _ in prim), default=14.0)
        withi = prim + [(last + 7.0, v["pseudobulk_age"])
                        for v in tps.values() if v["day"] is None]
        rho_i = spearman([d for d, _ in withi], [a for _, a in withi])
    else:
        rho_i = float("nan")

    print(f"\n  D2 PRIMARY (D0-D14, iPSC excluded) : {verdict['status']}")
    print(f"      {verdict['reason']}")
    print(f"  D2 sensitivity (+iPSC)             : rho {rho_i:+.3f}"
          "   [reported only; iPSC is a cell-type change, not aging]")

    out = {"script": "diag_d2_replication", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "gse_dir": gse_dir, "bars": bars(), "timepoints": tps, "meta": meta,
           "verdict_primary_D0_D14": verdict, "rho_with_ipsc_sensitivity": float(rho_i)}
    _RESULTS / "diag_d2_replication_results.json".write_text(json.dumps(out, indent=2, default=str),
                                                        encoding="utf-8")
    print("\n  wrote diag_d2_replication_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
