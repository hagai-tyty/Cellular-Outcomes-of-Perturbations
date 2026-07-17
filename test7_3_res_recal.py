"""
Test 7.3 (ΔAge lab notebook) — does RECALIBRATING FATE rescue RES for ranking?

Tests 7 / 7.1 / 7.2 showed the RES score ranks rejuvenation WORSE than a plain ΔAge sort, and
Test 7.2 isolated the cause to the RES *formula* — its fate terms (S, P_loss), which Test 8.2 showed
are well-ranked but MISCALIBRATED out-of-donor (ECE ~0.28; Platt ~halves it).

Salvage question: if we Platt-recalibrate the fate probabilities on the calib split BEFORE feeding
them to RES, does RES stop scrambling the ranking?

Clean isolation (same design as 7.2): hold ΔAge fixed = ridge ΔAge for every ranking, so the ONLY
thing that changes from B to C is the fate calibration.
  A = rank by ridge ΔAge directly           (current champion, ~0.95 vs true ΔAge)
  B = RES(ridge ΔAge + RAW model fate)       (reproduces Test 7.2's losing RES — a built-in check)
  C = RES(ridge ΔAge + RECALIBRATED fate)    (this test)

NOTE ON RECALIBRATION SCOPE: RES uses both the safe prob S and the loss prob P_loss, so this
recalibrates BOTH independently (Platt per class on the calib split). Test 8.2 validated Platt only
on the SAFE probability; recalibrating P_loss too is the faithful "recalibrated fate for RES" choice,
but means S_re + P_loss_re + P_death no longer sum to 1 (RES doesn't require it). A safe-prob-only
variant is a one-line change (drop the Ploss_re line, pass P_loss) if you want to match 8.2 exactly.

Scored vs (i) true ΔAge and (ii) safe-rejuvenation (gated), per fold + aggregate Spearman + paired
95% CIs for (C-B) [did recal help?] and (C-A) [does recal make RES earn its keep?].

READ:
  - C > B beyond noise      -> recalibration HELPS RES (miscalibration was part of the damage).
  - C >= A (tie or better)  -> recalibrated RES finally matches/earns its keep -> RES salvageable.
  - C < A still             -> even calibrated, RES ranks worse than a plain ΔAge sort -> for ranking
                               use ΔAge directly; RES needs redesign, not just recalibration.

USAGE (repo root, venv active; needs cellfate_loocv_* bundles + their calib/test splits):
    python test7_3_res_recal.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import DEATH_IDX, LOSS_IDX, SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def ridge_pred(tr, te) -> np.ndarray:
    """Ridge ΔAge on [state, fingerprint, dose_time] — the best ΔAge predictor (Test 6)."""
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(te.X), np.asarray(te.fp, float), sdt.transform(te.dose_time)])
    m = tr.mask
    return Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m]).predict(fte)


def _platt(p_cal, is_pos, p_te) -> np.ndarray:
    """Platt scaling: logistic(calib prob -> calib outcome), applied to test probs.
    Identity fallback if the calib split lacks both classes."""
    p_cal = np.asarray(p_cal, float).reshape(-1, 1)
    is_pos = np.asarray(is_pos, int)
    if not (0 < is_pos.sum() < len(is_pos)):
        return np.asarray(p_te, float)
    lr = LogisticRegression(max_iter=1000).fit(p_cal, is_pos)
    return lr.predict_proba(np.asarray(p_te, float).reshape(-1, 1))[:, 1]


def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


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


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
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
    sigma_age = np.array([r["sigma_age"] for r in te_rows])
    in_dist = np.array([r["in_dist"] for r in te_rows])

    # Platt-recalibrate the fate probabilities on the calib split (Test 8.2's fix).
    # If the calib split is missing, fall back to raw (then C == B for this fold).
    recal = False
    try:
        cal = gather_split(paths, REGIME, "calib")
        cal_rows = est.rows(cal.X, cal.fp, cal.dose_time)
        S_cal = np.array([r["S"] for r in cal_rows])
        Ploss_cal = np.array([r["P_loss"] for r in cal_rows])
        cls_cal = cal.y_cls.astype(int)
        S_re = _platt(S_cal, cls_cal == SAFE_IDX, S)
        Ploss_re = _platt(Ploss_cal, cls_cal == LOSS_IDX, P_loss)
        recal = True
    except Exception:  # noqa: BLE001
        S_re, Ploss_re = S, P_loss

    mu = ridge_pred(tr, te)                            # SAME ΔAge fed to A, B and C

    A = -mu                                            # rank by ridge ΔAge directly
    B, _ = compute_res_batch(S, P_loss, mu, sigma_age, in_dist, pred.res_params)        # raw fate
    C, _ = compute_res_batch(S_re, Ploss_re, mu, sigma_age, in_dist, pred.res_params)   # recal fate

    y_age = te.y_age[m]
    y_cls = te.y_cls[m].astype(int)
    unsafe = (y_cls == LOSS_IDX) | (y_cls == DEATH_IDX)
    gated = (-y_age).astype(float)
    if unsafe.any() and (~unsafe).any():
        gated[unsafe] = (-y_age)[~unsafe].min() - 1.0

    return {
        "recal": recal,
        "A_true": spearman(A[m], -y_age), "B_true": spearman(B[m], -y_age),
        "C_true": spearman(C[m], -y_age),
        "A_safe": spearman(A[m], gated), "B_safe": spearman(B[m], gated),
        "C_safe": spearman(C[m], gated),
    }


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7.3 — does RECALIBRATED fate rescue RES for ranking? (ΔAge held = ridge)")
    print("A = ridge ΔAge sort.  B = RES(raw fate).  C = RES(Platt-recalibrated fate).")
    per_fold = {d: r for d in DONORS if (r := one_fold(d)) is not None}
    if not per_fold:
        print("\n   No folds found (need cellfate_loocv_* with calib+test splits).")
        return
    n_recal = sum(v["recal"] for v in per_fold.values())
    if n_recal < len(per_fold):
        print(f"   (note: {len(per_fold) - n_recal}/{len(per_fold)} folds had no usable calib split -> "
              f"C=B there)")

    for tgt, label in [("true", "vs TRUE ΔAge"), ("safe", "vs SAFE-REJUVENATION (gated)")]:
        rows = [[d, f"{per_fold[d][f'A_{tgt}']:+.3f}", f"{per_fold[d][f'B_{tgt}']:+.3f}",
                 f"{per_fold[d][f'C_{tgt}']:+.3f}",
                 f"{per_fold[d][f'C_{tgt}'] - per_fold[d][f'B_{tgt}']:+.3f}",
                 f"{per_fold[d][f'C_{tgt}'] - per_fold[d][f'A_{tgt}']:+.3f}"] for d in per_fold]
        print(f"\n  {label}")
        print(render_table(["fold", "A ridge ΔAge", "B RES(raw)", "C RES(recal)", "C - B", "C - A"],
                           rows, aligns=["l", "r", "r", "r", "r", "r"]))
        aggA = np.nanmean([per_fold[d][f'A_{tgt}'] for d in per_fold])
        aggB = np.nanmean([per_fold[d][f'B_{tgt}'] for d in per_fold])
        aggC = np.nanmean([per_fold[d][f'C_{tgt}'] for d in per_fold])
        print(f"   aggregate: A={aggA:.3f}   B={aggB:.3f}   C={aggC:.3f}")
        for name, hi, lo in [("C - B (did recal help?)", "C", "B"),
                             ("C - A (does RES now earn its keep?)", "C", "A")]:
            diffs = [per_fold[d][f'{hi}_{tgt}'] - per_fold[d][f'{lo}_{tgt}'] for d in per_fold]
            md, (clo, chi), nn = paired_ci(diffs)
            sign = ("HIGHER" if clo > 0 else "LOWER" if chi < 0 else "tied (noise)")
            print(f"     {name}: mean={md:+.3f} 95% CI=[{clo:+.3f},{chi:+.3f}] (n={nn}) -> {sign}")

    print("\n   WHAT THIS MEANS:")
    print("     - C > B  -> recalibrating fate HELPS RES (miscalibration was part of the damage).")
    print("     - C >= A -> recalibrated RES finally matches/earns its keep -> RES is salvageable.")
    print("     - C < A  -> even calibrated, RES ranks worse than a plain ΔAge sort -> for ranking")
    print("                 use ΔAge directly; RES needs redesign, not just recalibration.")


if __name__ == "__main__":
    main()
