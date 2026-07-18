"""
Test 7.4.1 (ΔAge lab notebook) — PATCH of Test 7.4's broken half.

WHY THIS EXISTS. Test 7.4 tried to test whether recalibrating fate helps RES as an *actionable*
score, but defined the decision cutoff as a CALIB-SET QUANTILE. A quantile is
MONOTONE-EQUIVARIANT: Platt scaling shifts the scores and the quantile shifts identically, so
exactly the same cells are flagged — by construction. (Signature: C-B precision = +0.000 with a
zero-width CI.) That half of 7.4 measured nothing, the same structural blindness as 7.3's
rank-invariant Spearman.

THE FIX. RES is documented as living in [0, 1), so a FIXED ABSOLUTE cutoff is meaningful and is
NOT monotone-equivariant: Platt changes the score VALUES, so more/fewer cells cross a fixed
numeric bar. Two genuinely value-based decision tests:

  (1) FIXED CUTOFF SWEEP — flag cells with RES >= c for c in {0.05, 0.10, 0.20, 0.30, 0.50};
      report how many are flagged and what fraction are truly SAFE-AND-REJUVENATING.
  (2) THE APPROVED GATE — compute_res_batch also returns a status, and REJECTED_UNSAFE fires on
      the ABSOLUTE test S < tau_safe - 3*w. Recalibration changes S, so it moves this gate.
      This is the product's real accept/reject decision, not a proxy.

BUILT-IN SENSITIVITY CHECK (the lesson from 7.4): the script explicitly reports whether the
flagged SETS differ between B and C. If they never differ, the test is blind again and says so
loudly instead of reporting a fake "tied".

Same isolation as 7.3/7.4: ΔAge held identical (ridge) for both arms, so the ONLY difference is
fate calibration.
  B = RES(ridge ΔAge + RAW model fate)
  C = RES(ridge ΔAge + Platt-RECALIBRATED fate)
(A = plain ΔAge sort is NOT included here: it lives on a different scale (years), so absolute
RES cutoffs are not comparable to it. A-vs-C was already settled in 7.3: C-A = -0.298, CI
excludes 0. This test is purely "does recalibration help RES on a value-based decision?")

READ:
  - C flags a better set than B (higher precision at the same cutoffs, or a better APPROVED set)
    -> recalibration DOES improve RES as an actionable score.
  - sets differ but precision is tied -> recalibration moves the decision without improving it.
  - sets NEVER differ -> the test is blind; report that, do not claim "tied".

USAGE (repo root, venv active; needs cellfate_loocv_* bundles + calib/test splits):
    python test7_4_1_res_fixed_threshold.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import LOSS_IDX, SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
CUTOFFS = [0.05, 0.10, 0.20, 0.30, 0.50]   # FIXED absolute cutoffs on the native RES scale
APPROVED = "APPROVED"


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def ridge_pred(tr, target_split) -> np.ndarray:
    """Ridge ΔAge on [state, fingerprint, dose_time] — best ΔAge predictor (Test 6)."""
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(target_split.X), np.asarray(target_split.fp, float),
                     sdt.transform(target_split.dose_time)])
    m = tr.mask
    return Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m]).predict(fte)


def _platt(p_cal, is_pos, p_te) -> np.ndarray:
    """Platt scaling fit on the calib split, applied to test probs. Identity if calib is
    single-class."""
    p_cal = np.asarray(p_cal, float).reshape(-1, 1)
    is_pos = np.asarray(is_pos, int)
    if not (0 < is_pos.sum() < len(is_pos)):
        return np.asarray(p_te, float)
    lr = LogisticRegression(max_iter=1000).fit(p_cal, is_pos)
    return lr.predict_proba(np.asarray(p_te, float).reshape(-1, 1))[:, 1]


def paired_ci(diffs):
    diffs = [d for d in diffs if np.isfinite(d)]
    n = len(diffs)
    if n < 2:
        return float("nan"), (float("nan"), float("nan")), n
    md = sum(diffs) / n
    sd = math.sqrt(sum((x - md) ** 2 for x in diffs) / (n - 1))
    se = sd / math.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n


def prec_rec(flag, good):
    """Precision/recall of retrieving truly safe-and-rejuvenating cells."""
    if flag.sum() == 0:
        return float("nan"), 0.0
    return float(good[flag].mean()), float(good[flag].sum() / max(good.sum(), 1))


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        cal = gather_split(paths, REGIME, "calib")
    except Exception:  # noqa: BLE001
        return None
    m = te.mask
    if m.sum() < 3:
        return None

    pred = Predictor(root)
    est = ModelEstimator(pred)

    te_rows = est.rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in te_rows])
    P_loss = np.array([r["P_loss"] for r in te_rows])
    sig = np.array([r["sigma_age"] for r in te_rows])
    ind = np.array([r["in_dist"] for r in te_rows])

    cal_rows = est.rows(cal.X, cal.fp, cal.dose_time)
    S_cal = np.array([r["S"] for r in cal_rows])
    Pl_cal = np.array([r["P_loss"] for r in cal_rows])
    cls_cal = cal.y_cls.astype(int)
    S_re = _platt(S_cal, cls_cal == SAFE_IDX, S)
    Pl_re = _platt(Pl_cal, cls_cal == LOSS_IDX, P_loss)

    mu = ridge_pred(tr, te)                       # SAME ΔAge in both arms

    B, statB = compute_res_batch(S, P_loss, mu, sig, ind, pred.res_params)
    C, statC = compute_res_batch(S_re, Pl_re, mu, sig, ind, pred.res_params)

    good = ((te.y_cls[m].astype(int) == SAFE_IDX) & (te.y_age[m] < 0)).astype(float)
    if good.sum() == 0 or good.sum() == len(good):
        return {"degenerate": True}

    Bm, Cm = B[m], C[m]
    out = {"degenerate": False, "n": int(m.sum()), "n_good": int(good.sum()),
           "B_range": (float(Bm.min()), float(Bm.max())),
           "C_range": (float(Cm.min()), float(Cm.max()))}

    for c in CUTOFFS:
        fB, fC = Bm >= c, Cm >= c
        out[f"setdiff@{c}"] = int(np.sum(fB != fC))     # SENSITIVITY: do the sets differ at all?
        out[f"nB@{c}"], out[f"nC@{c}"] = int(fB.sum()), int(fC.sum())
        pB, rB = prec_rec(fB, good)
        pC, rC = prec_rec(fC, good)
        out[f"pB@{c}"], out[f"rB@{c}"] = pB, rB
        out[f"pC@{c}"], out[f"rC@{c}"] = pC, rC

    aB, aC = (statB[m] == APPROVED), (statC[m] == APPROVED)
    out["appr_setdiff"] = int(np.sum(aB != aC))
    out["nB_appr"], out["nC_appr"] = int(aB.sum()), int(aC.sum())
    out["pB_appr"], out["rB_appr"] = prec_rec(aB, good)
    out["pC_appr"], out["rC_appr"] = prec_rec(aC, good)
    return out


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7.4.1 — does recalibration help RES at a FIXED ABSOLUTE cutoff? (patches 7.4)")
    print("7.4's cutoff was a calib QUANTILE (monotone-equivariant => blind). RES lives in [0,1),")
    print("so fixed numeric cutoffs CAN move under Platt. B = RES(raw fate), C = RES(recal fate);")
    print("ΔAge identical (ridge) in both arms. 'good' = truly SAFE and rejuvenating (ΔAge<0).")

    per = {d: r for d in DONORS if (r := one_fold(d)) is not None and not r.get("degenerate")}
    if not per:
        print("\n   No usable folds (need calib+test with safe-rejuv variation).")
        return
    skipped = [d for d in DONORS if d not in per]
    if skipped:
        print(f"   (folds without safe-rejuv variation, skipped: {', '.join(skipped)})")

    print("\n  SENSITIVITY CHECK — do B and C flag DIFFERENT cells? (0 everywhere => test is blind)")
    rows = [[d] + [str(per[d][f"setdiff@{c}"]) for c in CUTOFFS] + [str(per[d]["appr_setdiff"])]
            for d in per]
    print(render_table(["fold"] + [f"RES>={c}" for c in CUTOFFS] + ["APPROVED"],
                       rows, aligns=["l"] + ["r"] * (len(CUTOFFS) + 1)))
    total_diff = sum(per[d][f"setdiff@{c}"] for d in per for c in CUTOFFS) + \
        sum(per[d]["appr_setdiff"] for d in per)
    if total_diff == 0:
        print("   !! B and C flag IDENTICAL cells everywhere -> this test is BLIND too.")
        print("      Do NOT read the numbers below as 'tied' — recalibration did not move any")
        print("      decision boundary at these cutoffs (likely the scores barely shifted).")
    else:
        print(f"   OK: {total_diff} cell-level decision flips total -> the metric CAN see recal.")

    for c in CUTOFFS:
        rows = [[d, f"{per[d][f'nB@{c}']}/{per[d]['n']}", f"{per[d][f'nC@{c}']}/{per[d]['n']}",
                 f"{per[d][f'pB@{c}']:.2f}" if np.isfinite(per[d][f'pB@{c}']) else "n/a",
                 f"{per[d][f'pC@{c}']:.2f}" if np.isfinite(per[d][f'pC@{c}']) else "n/a",
                 f"{per[d][f'rB@{c}']:.2f}", f"{per[d][f'rC@{c}']:.2f}"] for d in per]
        print(f"\n  FIXED CUTOFF RES >= {c}")
        print(render_table(["fold", "n flagged B", "n flagged C", "prec B", "prec C",
                            "rec B", "rec C"], rows,
                           aligns=["l", "r", "r", "r", "r", "r", "r"]))
        diffs = [per[d][f'pC@{c}'] - per[d][f'pB@{c}'] for d in per
                 if np.isfinite(per[d][f'pC@{c}']) and np.isfinite(per[d][f'pB@{c}'])]
        md, (lo, hi), n = paired_ci(diffs)
        v = ("recal IMPROVES precision" if lo > 0 else "recal WORSENS" if hi < 0 else "tied")
        print(f"   paired C-B precision: mean={md:+.3f} CI=[{lo:+.3f},{hi:+.3f}] (n={n}) -> {v}")

    print("\n  THE APPROVED GATE (the product's real accept/reject decision)")
    rows = [[d, f"{per[d]['nB_appr']}/{per[d]['n']}", f"{per[d]['nC_appr']}/{per[d]['n']}",
             f"{per[d]['pB_appr']:.2f}" if np.isfinite(per[d]['pB_appr']) else "n/a",
             f"{per[d]['pC_appr']:.2f}" if np.isfinite(per[d]['pC_appr']) else "n/a",
             f"{per[d]['rB_appr']:.2f}", f"{per[d]['rC_appr']:.2f}"] for d in per]
    print(render_table(["fold", "n appr B", "n appr C", "prec B", "prec C", "rec B", "rec C"],
                       rows, aligns=["l", "r", "r", "r", "r", "r", "r"]))
    diffs = [per[d]['pC_appr'] - per[d]['pB_appr'] for d in per
             if np.isfinite(per[d]['pC_appr']) and np.isfinite(per[d]['pB_appr'])]
    md, (lo, hi), n = paired_ci(diffs)
    v = ("recal IMPROVES approved-set precision" if lo > 0 else "recal WORSENS" if hi < 0
         else "tied")
    print(f"   paired C-B precision: mean={md:+.3f} CI=[{lo:+.3f},{hi:+.3f}] (n={n}) -> {v}")

    print("\n  RES score ranges per fold (context: are the scores even near the cutoffs?)")
    rows = [[d, f"[{per[d]['B_range'][0]:.3f}, {per[d]['B_range'][1]:.3f}]",
             f"[{per[d]['C_range'][0]:.3f}, {per[d]['C_range'][1]:.3f}]",
             f"{per[d]['n_good']}/{per[d]['n']}"] for d in per]
    print(render_table(["fold", "B RES range", "C RES range", "good cells"],
                       rows, aligns=["l", "r", "r", "r"]))

    print("\n   WHAT THIS MEANS:")
    print("     - Read the SENSITIVITY CHECK first. If all zeros, the test is blind and the")
    print("       precision numbers carry no information (7.4's mistake — do not repeat it).")
    print("     - If sets differ and C's precision is higher -> recalibration DOES improve RES")
    print("       as an actionable score, even though it cannot change rank (7.3).")
    print("     - If sets differ but precision is tied -> recal moves the decision boundary")
    print("       without improving the decision -> RES's problem is the formula, not its inputs.")
    print("     - NOTE: 7.3 already showed recalibrated RES still ranks far BELOW a plain ΔAge")
    print("       sort (C-A = -0.298, CI excludes 0). This test can only change 'how much RES")
    print("       underperforms', not whether it does.")


if __name__ == "__main__":
    main()
