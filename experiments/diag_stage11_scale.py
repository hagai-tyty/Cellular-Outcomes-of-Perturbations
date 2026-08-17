"""STAGE 11 -- is the ΔAge scale error fixable by calibration?  (read-only)

    python experiments/diag_stage11_scale.py

Pre-registered in `plans/STAGE_11_DAGE_SCALE_CALIBRATION.md`, written before this ran.
`src/` is NOT touched by this stage under any outcome.

WHY
---
Stage 10 established that `raw`'s problem is NOT pluripotency -- it is SCALE. Raw ΔAge has the right
ORDERING (Spearman 0.770 vs methylation) and a 66% magnitude inflation (SD ratio 1.66). Nobody had
tried simply correcting the scale.

A MATHEMATICAL FACT, STATED BEFORE THE RUN
------------------------------------------
A pure linear rescale y -> k*y CANNOT change Spearman: rank order is invariant to any positive
monotone transform. So rescaled-raw's rho WILL be exactly 0.770. That is arithmetic, not a result,
and is not reported as a finding. It is what makes this stage clean: any MAE improvement is
attributable to SCALE ALONE, and the residual gap to top100's 0.810 is the part scale cannot
explain.

THE LEAK THAT MUST NOT HAPPEN
-----------------------------
`k` is fitted LEAVE-ONE-DONOR-OUT. Fitting it on all 44 conditions and scoring the same conditions
would guarantee an improvement and measure nothing.

PRE-REGISTERED READING (plan 11.3)
  SCALE IS THE PROBLEM   LODO-calibrated MAE <= 1.5 x the 7.30 floor (= 10.95)
  SCALE IS PART OF IT    improves but stays above 10.95; report how much of the gap closes
  NOT SCALE              no improvement
Plus, independently: `k` STABILITY across folds. Varying by more than 2x means calibration is not
transferable, and that caveat outranks any MAE gain.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))

LEDGER = ROOT / "results" / "dage_ledger.csv"
FLOOR = 7.30            # methylation vs methylation MAE, 44 conditions
FLOOR_MULT = 1.5        # 11.3: <= FLOOR_MULT * FLOOR  ->  SCALE IS THE PROBLEM
K_STABILITY_BAR = 2.0   # max/min of k across folds; above this, not transferable
VARIANTS = ("raw", "top100", "top500", "top2000", "resid_pluri")


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def fit_scale(pred: np.ndarray, truth: np.ndarray, with_offset: bool) -> tuple[float, float]:
    """k (and c) minimising squared error of k*pred + c against truth.

    Without an offset this is the least-squares scale through the origin, which is the right form
    for a magnitude-inflation defect: ΔAge is already control-relative, so zero must stay zero.

    CAUTION, and this is why `fit_scale_variance` exists beside it: least squares is a SHRINKAGE
    estimator. k_LS = rho * SD(truth)/SD(pred), so with imperfect correlation it deliberately
    shrinks BELOW the variance-matching value to minimise MSE. That lowers MAE while making the
    calibrated ΔAge systematically UNDER-report magnitude -- the same trade that made `ranknorm`
    and `resid_pluri` look good on MAE. A number used to claim "these cells got N years younger"
    must not be silently shrunk.
    """
    if with_offset:
        k, c = np.polyfit(pred, truth, 1)
        return float(k), float(c)
    denom = float(pred @ pred)
    return (float(pred @ truth / denom) if denom > 1e-12 else 1.0), 0.0


def fit_scale_variance(pred: np.ndarray, truth: np.ndarray) -> tuple[float, float]:
    """k matching SPREAD rather than minimising error: k = SD(truth)/SD(pred).

    Preserves magnitude (SD ratio -> 1.0) at the cost of a higher MAE than least squares. This is
    the calibration to use if the number is reported as an amount of rejuvenation."""
    sp = float(np.std(pred, ddof=1))
    return (float(np.std(truth, ddof=1) / sp) if sp > 1e-12 else 1.0), 0.0


def lodo_calibrate(pred: np.ndarray, truth: np.ndarray, donor: np.ndarray,
                   with_offset: bool, mode: str = "ls") -> tuple[np.ndarray, dict]:
    """Leave-one-DONOR-out calibration. k is never fitted on the rows it is scored on.

    mode="ls"  least squares (minimises MAE, SHRINKS magnitude)
    mode="var" variance matching (preserves magnitude, higher MAE)
    """
    out = np.empty_like(pred)
    ks: dict[str, float] = {}
    for d in sorted(set(donor)):
        te = donor == d
        tr = ~te
        if tr.sum() < 2:
            out[te], ks[str(d)] = pred[te], 1.0
            continue
        k, c = (fit_scale_variance(pred[tr], truth[tr]) if mode == "var"
                else fit_scale(pred[tr], truth[tr], with_offset))
        out[te] = k * pred[te] + c
        ks[str(d)] = k
    return out, ks


def verdict_from(mae: float) -> str:
    if not np.isfinite(mae):
        return "UNDETERMINED"
    if mae <= FLOOR_MULT * FLOOR:
        return "SCALE IS THE PROBLEM"
    return "SCALE IS PART OF IT"


def main() -> None:
    led = pd.read_csv(LEDGER)
    t = led[(~led.is_control.astype(bool)) & led.TRUTH_meth_dage_mt.notna()].copy()
    mt = t.TRUTH_meth_dage_mt.to_numpy(float)
    donor = t.donor.to_numpy()

    print("=" * 104)
    print("STAGE 11 -- can the ΔAge scale error be calibrated away?")
    print(f"  floor (methylation vs methylation) MAE {FLOOR}   bar for 'SCALE IS THE PROBLEM': "
          f"<= {FLOOR_MULT * FLOOR:.2f}")
    print(f"  n={len(t)} conditions, donors {sorted(set(donor))} -- k fitted LEAVE-ONE-DONOR-OUT")
    print("  NOTE: a pure rescale cannot change Spearman. Any rho below is unchanged BY ARITHMETIC.")
    print("=" * 104)
    print(f"  {'':<12}{'UNCALIBRATED':>16} | {'LEAST SQUARES':>14} | {'VAR-MATCHED':>14} |")
    print(f"  {'variant':<12}{'MAE':>9}{'SD':>7} |{'MAE':>8}{'SD':>7} |{'MAE':>8}{'SD':>7} |"
          f"{'rho':>7}   k_ls per donor")
    res: dict = {"n": int(len(t)), "floor": FLOOR, "variants": {}}

    for v in VARIANTS:
        col = f"ACTUAL_rna_dage_{v}"
        if col not in t.columns:
            continue
        p = t[col].to_numpy(float)
        ok = np.isfinite(p) & np.isfinite(mt)
        p, y, d = p[ok], mt[ok], donor[ok]
        mae0 = float(np.abs(p - y).mean())
        sd0 = float(np.std(p, ddof=1) / np.std(y, ddof=1))
        rho = spearman(p, y)

        cal_k, ks = lodo_calibrate(p, y, d, with_offset=False)
        cal_kc, _ = lodo_calibrate(p, y, d, with_offset=True)
        cal_v, ksv = lodo_calibrate(p, y, d, with_offset=False, mode="var")
        mae_k = float(np.abs(cal_k - y).mean())
        mae_kc = float(np.abs(cal_kc - y).mean())
        mae_v = float(np.abs(cal_v - y).mean())
        sd_k = float(np.std(cal_k, ddof=1) / np.std(y, ddof=1))
        sd_v = float(np.std(cal_v, ddof=1) / np.std(y, ddof=1))
        kvals = list(ks.values())
        kspread = max(kvals) / min(kvals) if min(kvals) > 1e-9 else float("inf")

        print(f"  {v:<12}{mae0:>9.2f}{sd0:>7.2f} |{mae_k:>8.2f}{sd_k:>7.2f} |{mae_v:>8.2f}"
              f"{sd_v:>7.2f} |{rho:>7.3f}   " + " ".join(f"{k}={val:.2f}" for k, val in ks.items()))
        res["variants"][v] = {"mae_raw": mae0, "mae_scaled": mae_k, "mae_scaled_offset": mae_kc,
                              "mae_var_matched": mae_v, "sd_ratio_raw": sd0,
                              "sd_ratio_scaled": sd_k, "sd_ratio_var_matched": sd_v,
                              "spearman": rho, "k_per_donor": ks, "k_var_per_donor": ksv,
                              "k_spread": kspread, "verdict": verdict_from(mae_k)}

    r = res["variants"]["raw"]
    print(f"\n  RAW: {r['mae_raw']:.2f} -> {r['mae_scaled']:.2f} after a ONE-PARAMETER LODO rescale")
    gap_before = r["mae_raw"] - res["variants"]["top100"]["mae_raw"]
    gap_after = r["mae_scaled"] - res["variants"]["top100"]["mae_raw"]
    closed = (1 - gap_after / gap_before) * 100 if abs(gap_before) > 1e-9 else float("nan")
    print(f"  gap to top100: {gap_before:+.2f} -> {gap_after:+.2f} yr   "
          f"({closed:.0f}% of the gap closed by scale alone)")
    print(f"  k stability across donors: spread {r['k_spread']:.2f}x "
          f"(bar {K_STABILITY_BAR}x)  -> {'STABLE' if r['k_spread'] <= K_STABILITY_BAR else 'NOT TRANSFERABLE'}")
    print(f"  -> {r['verdict']}")
    print(f"\n  ordering, unchanged by any rescale: raw {r['spearman']:.3f} vs top100 "
          f"{res['variants']['top100']['spearman']:.3f}  -- the part scale CANNOT explain")
    res["gap_closed_pct"] = closed
    res["k_stable"] = bool(r["k_spread"] <= K_STABILITY_BAR)
    res["verdict"] = r["verdict"]

    _RESULTS.mkdir(exist_ok=True)
    p_out = _RESULTS / "diag_stage11_scale_results.json"
    p_out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p_out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
