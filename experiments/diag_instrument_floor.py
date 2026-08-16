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
**PASS** for a given RNA variant and reference clock R requires BOTH:

    (i)   MAE(rna, R)  <=  MAE(mt, sb)         inside the floor
    (ii)  MAE(rna, R)  <   MAE(0, R)           beats a CONSTANT-ZERO predictor

**(ii) is not optional and the test is close to worthless without it.** A predictor that collapses
toward zero scores a low MAE while carrying no information, and §0 caught exactly that on the Sendai
arm: `top100`'s apparent gain there was the prediction sliding onto the clock's intercept
b0 = 72.43, its SD falling to 6.50 against truth's 14.75. So a constant-zero predictor is scored
against both references on the SAME rows, and a variant inside the floor that does NOT beat it is
reported as `shrinkage!` rather than as a pass.

Reported per variant, ungraded but diagnostic: **SD ratio** against the truth (a shrunk predictor
has SD << 1) and **Spearman rho** — because §2's lesson was that a correlation cannot see a bias,
and its inverse is equally true: an MAE cannot see a lost ordering.

A 95 % bootstrap CI (10 000 resamples over conditions, paired) is reported on
`MAE(rna,R) - MAE(mt,sb)`; **PASS requires the point estimate to clear, NOT the CI to exclude
zero** — at n = 44 demanding significance would be a bar no correct system could clear (§5b).

**All nine ledger variants are scored, not a subset.** Choosing which to report after seeing the
numbers is how a family-wise result becomes a cherry-picked one.

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
# ALL nine ledger variants, not a chosen subset -- picking which to score after seeing the
# numbers is how a family-wise result becomes a cherry-picked one.
VARIANTS = ("raw", "top100", "top500", "top2000", "covnorm", "ranknorm",
            "resid_cc", "resid_pluri", "resid_both")


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
        # ---- THE SHRINKAGE CONTROL, without which the floor test is gameable ---------- #
        # A predictor that collapses toward zero scores a LOW MAE while carrying no
        # information. §0 caught exactly this on the Sendai arm, where top100's "improvement"
        # was the prediction sliding onto the clock's intercept. So a constant-zero predictor
        # is scored against the floor on the SAME rows, and any variant that fails to beat it
        # is disqualified no matter how well it does against the floor.
        zero = np.zeros(len(sub))
        z_mt, z_sb = zero - mt, zero - sb
        rec: dict = {"n": len(sub),
                     "floor_mt_vs_sb": {"mae": mae(floor), "rms": rms(floor),
                                        "mean_signed": float(np.mean(floor))},
                     "meth_dage_mean": {"mt": float(mt.mean()), "sb": float(sb.mean())},
                     "meth_dage_sd": {"mt": float(mt.std(ddof=1)), "sb": float(sb.std(ddof=1))},
                     "constant_zero_control": {"mae_vs_mt": mae(z_mt), "mae_vs_sb": mae(z_sb)},
                     "variants": {}}

        print(f"\n{'='*100}\n[{arm}]  n = {len(sub)} conditions")
        print(f"  THE FLOOR — methylation vs methylation (mt − sb): "
              f"MAE {mae(floor):6.2f}   RMS {rms(floor):6.2f}   mean {np.mean(floor):+6.2f}")
        print(f"  truth SD: mt {mt.std(ddof=1):5.2f}  sb {sb.std(ddof=1):5.2f}    "
              f"CONSTANT-ZERO control: MAE {mae(z_mt):5.2f} (mt)  {mae(z_sb):5.2f} (sb)")
        print(f"\n  {'variant':>11} {'vs':>3} {'MAE':>7} {'Δ floor':>8} {'95% CI':>18} "
              f"{'SDratio':>8} {'rho':>6} {'beats0':>7}  {'verdict':>9}")

        for v in VARIANTS:
            col = f"ACTUAL_rna_dage_{v}"
            if col not in sub[0]:
                print(f"  {v:>11}  -- column absent from the ledger, skipped")
                continue
            rna = np.array([_f(r, col) for r in sub])
            if not np.isfinite(rna).all():
                print(f"  {v:>11}  -- non-finite values, skipped")
                continue
            vrec: dict = {"sd": float(rna.std(ddof=1)), "mean": float(rna.mean())}
            for name, ref, zc in (("mt", mt, z_mt), ("sb", sb, z_sb)):
                d = rna - ref
                lo, hi = boot_ci(d, floor, rng)
                dm = mae(d) - mae(floor)
                inside = mae(d) <= mae(floor)
                beats0 = mae(d) < mae(zc)
                sd_ratio = float(rna.std(ddof=1) / ref.std(ddof=1))
                # Spearman without scipy: Pearson on ranks.
                ra = np.argsort(np.argsort(rna)).astype(float)
                rb = np.argsort(np.argsort(ref)).astype(float)
                rho = float(np.corrcoef(ra, rb)[0, 1])
                ok = bool(inside and beats0)
                vrec[name] = {"mae": mae(d), "rms": rms(d), "mean_signed": float(np.mean(d)),
                              "delta_mae_vs_floor": dm, "ci95": [lo, hi],
                              "inside_floor": bool(inside), "beats_constant_zero": bool(beats0),
                              "sd_ratio_vs_truth": sd_ratio, "spearman": rho, "pass": ok}
                verdict = "PASS" if ok else ("shrinkage!" if inside and not beats0 else "fail")
                print(f"  {v:>11} {name:>3} {mae(d):7.2f} {dm:+8.2f} "
                      f"[{lo:+6.2f},{hi:+6.2f}] {sd_ratio:8.2f} {rho:+6.2f} "
                      f"{'yes' if beats0 else 'NO':>7}  {verdict:>9}")
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
