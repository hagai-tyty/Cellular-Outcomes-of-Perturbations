"""LIKE-FOR-LIKE: is RNA-ΔAge's disagreement with methylation inside methylation's OWN?

    python experiments/diag_instrument_floor.py

READ-ONLY. Reads `results/dage_ledger.csv` only. Writes
`results/diag_instrument_floor_results.json`. `src/` untouched, no label moves, no data re-read.

WHY THIS EXISTS
---------------
`DELTA_AGE_WHERE_IT_STANDS.md` §2 put two numbers side by side and then refused to compare them:

  * §1's RNA-vs-methylation **MAE 5.36 yr**, on 68 conditions, transient arm, vs multi-tissue
  * `dage_meth_concordance`'s inter-clock **RMS 9.07 yr**, on 9 control groups, sb vs mt

**Different statistic, different sample set, different pairing.** §0's ERROR 2 is exactly that
mistake — comparing a condition-level ΔAge error to a donor-level absolute-age error — and it was
withdrawn. This script does the comparison properly instead of asserting it.

THE FIX: one table, one statistic, one set of rows
--------------------------------------------------
`results/dage_ledger.csv` already carries, **per condition**, both methylation truths and every RNA
variant. So every pairing below is computed on **exactly the same rows**, with **exactly the same
statistic**, paired condition-by-condition:

    d(A,B) = A_dage - B_dage      per condition, then MAE and RMS over conditions

  reference-vs-reference   mt  vs  sb      <- methylation's own self-disagreement = THE FLOOR
  RNA-vs-reference         rna vs  mt
                           rna vs  sb

Rows are restricted to those where **both** methylation truths exist, so the floor and the RNA
comparisons cannot be computed on different samples. Controls are excluded: a control's ΔAge is 0
against itself by construction and would deflate every number identically.

THE BAR, PRE-REGISTERED BEFORE THE NUMBERS
-------------------------------------------
**PASS** for a given RNA variant and reference clock R:

    MAE(rna, R)  <=  MAE(mt, sb)

i.e. the RNA readout disagrees with that methylation clock **no more than the two methylation clocks
disagree with each other**. A 95 % bootstrap CI (10 000 resamples over conditions, paired) is
reported on the difference `MAE(rna,R) - MAE(mt,sb)`; **PASS requires the point estimate to clear
and the CI to be reported, not to be significant** — with n this size, demanding significance would
be a bar no correct system could clear (§5b).

WHAT A PASS WOULD AND WOULD NOT MEAN
-------------------------------------
It would mean: **the remaining disagreement is at the scale of the references' own disagreement, so
it cannot be attributed to the RNA readout.** It would NOT mean the RNA clock is "correct" — two
instruments can agree and both be biased, and §1.5.2's factor-loading arithmetic already showed
these three instruments are *not* jointly consistent with one common age factor.

It also does NOT rescue same-timepoint ΔAge *prediction*, which is circular at ρ 0.96
(`diag_clock_circularity`). This is about the **measurement**, not about a model.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

LEDGER = ROOT / "results" / "dage_ledger.csv"
N_BOOT = 10000
SEED = 0
# The RNA variants worth the comparison: the shipped clock, and §1's sparse candidate.
VARIANTS = ("raw", "top100", "top500", "top2000")


def _f(row: dict, key: str) -> float:
    v = row.get(key, "")
    if v is None or str(v).strip() == "":
        return float("nan")
    try:
        return float(v)
    except ValueError:
        return float("nan")


def mae(d: np.ndarray) -> float:
    return float(np.mean(np.abs(d)))


def rms(d: np.ndarray) -> float:
    return float(np.sqrt(np.mean(d ** 2)))


def boot_ci(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95 % CI on MAE(a) - MAE(b), resampling CONDITIONS (paired), not residuals."""
    n = len(a)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    diffs = np.abs(a)[idx].mean(axis=1) - np.abs(b)[idx].mean(axis=1)
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if not LEDGER.exists():
        print(f"FATAL: {LEDGER} not found.")
        return 1

    rows = list(csv.DictReader(LEDGER.open(encoding="utf-8")))
    print(f"\n[ledger] {len(rows)} rows")

    # ---- the row set: non-control, BOTH methylation truths present -------------------- #
    keep = []
    for r in rows:
        if str(r.get("is_control", "")).strip().lower() in ("true", "1", "yes"):
            continue
        mt, sb = _f(r, "TRUTH_meth_dage_mt"), _f(r, "TRUTH_meth_dage_sb")
        if np.isfinite(mt) and np.isfinite(sb):
            keep.append(r)
    if not keep:
        print("FATAL: no non-control rows carry BOTH methylation truths.")
        return 1

    arms = sorted({str(r.get("arm_set", "?")) for r in keep})
    print(f"[rows]   {len(keep)} non-control conditions with BOTH truths | arms: {arms}")

    rng = np.random.default_rng(SEED)
    out: dict = {"script": "diag_instrument_floor", "utc": datetime.now(UTC).isoformat(),
                 "n_boot": N_BOOT, "seed": SEED, "ledger_rows": len(rows), "by_arm": {}}

    for arm in [*arms, "__POOLED__"]:
        sub = keep if arm == "__POOLED__" else [r for r in keep if str(r.get("arm_set")) == arm]
        if len(sub) < 4:
            print(f"\n[{arm}] only {len(sub)} conditions — skipped, too few to compare")
            out["by_arm"][arm] = {"n": len(sub), "skipped": "fewer than 4 conditions"}
            continue

        mt = np.array([_f(r, "TRUTH_meth_dage_mt") for r in sub])
        sb = np.array([_f(r, "TRUTH_meth_dage_sb") for r in sub])
        floor = mt - sb                                   # the reference self-disagreement
        rec: dict = {"n": len(sub),
                     "floor_mt_vs_sb": {"mae": mae(floor), "rms": rms(floor),
                                        "mean_signed": float(np.mean(floor))},
                     "meth_dage_mean": {"mt": float(mt.mean()), "sb": float(sb.mean())},
                     "meth_dage_sd": {"mt": float(mt.std(ddof=1)), "sb": float(sb.std(ddof=1))},
                     "variants": {}}

        print(f"\n{'='*78}\n[{arm}]  n = {len(sub)} conditions")
        print(f"  THE FLOOR — methylation vs methylation (mt − sb): "
              f"MAE {mae(floor):6.2f}   RMS {rms(floor):6.2f}   mean {np.mean(floor):+6.2f}")
        print(f"\n  {'variant':>9} {'vs':>3} {'MAE':>7} {'RMS':>7} {'ΔMAE vs floor':>14} "
              f"{'95% CI':>20}  {'':>6}")

        for v in VARIANTS:
            col = f"ACTUAL_rna_dage_{v}"
            if col not in sub[0]:
                continue
            rna = np.array([_f(r, col) for r in sub])
            if not np.isfinite(rna).all():
                print(f"  {v:>9}  -- non-finite values, skipped")
                continue
            vrec = {}
            for name, ref in (("mt", mt), ("sb", sb)):
                d = rna - ref
                lo, hi = boot_ci(d, floor, rng)
                dm = mae(d) - mae(floor)
                ok = mae(d) <= mae(floor)
                vrec[name] = {"mae": mae(d), "rms": rms(d), "mean_signed": float(np.mean(d)),
                              "delta_mae_vs_floor": dm, "ci95": [lo, hi], "pass": bool(ok)}
                print(f"  {v:>9} {name:>3} {mae(d):7.2f} {rms(d):7.2f} {dm:+14.2f} "
                      f"  [{lo:+7.2f},{hi:+7.2f}]  {'PASS' if ok else 'fail':>6}")
            vrec["pass_both_clocks"] = bool(vrec["mt"]["pass"] and vrec["sb"]["pass"])
            rec["variants"][v] = vrec
        out["by_arm"][arm] = rec

    # ---- the headline ----------------------------------------------------------------- #
    pooled = out["by_arm"].get("__POOLED__", {})
    winners = [v for v, r in pooled.get("variants", {}).items() if r.get("pass_both_clocks")]
    any_one = {v: [k for k in ("mt", "sb") if r[k]["pass"]]
               for v, r in pooled.get("variants", {}).items() if isinstance(r, dict)}
    out["verdict"] = {
        "variants_inside_floor_on_BOTH_clocks": winners,
        "variants_inside_floor_per_clock": any_one,
        "interpretation": ("PASS means the RNA readout disagrees with that methylation clock no "
                           "more than the two methylation clocks disagree with each other, on the "
                           "SAME conditions with the SAME statistic. It does NOT mean the RNA "
                           "clock is correct, and it does not rescue same-timepoint dAge "
                           "PREDICTION, which is circular at rho 0.96."),
    }
    print(f"\n{'='*78}")
    print(f"  inside the floor on BOTH reference clocks (pooled): "
          f"{winners if winners else 'NONE'}")
    print(f"  inside the floor on at least one:                   {any_one}")

    (_RESULTS / "diag_instrument_floor_results.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\n  wrote results/diag_instrument_floor_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
