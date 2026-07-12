"""
Test 7.2 (ΔAge lab notebook) — clean isolation of the RES FORMULA (user's design).

Test 7 confounded two things: model_RES used the MODEL's ΔAge, ridge_dAge used RIDGE's ΔAge.
So we couldn't tell if RES lost because of the FORMULA or because the model's ΔAge is worse.

Fix: hold ΔAge constant. Feed the SAME ΔAge (ridge's — the best predictor, Test 6) into both:
  A = rank by that ΔAge directly.
  B = rank by RES(that same ridge ΔAge as mu, + model fate S/P_loss, + model uncertainty, + OOD).
The ONLY difference is the RES transform. If B loses to A, it is unambiguously the RES formula
degrading the ranking — not a ΔAge-quality difference.

Scored against (i) true ΔAge and (ii) safe-rejuvenation (gated: unsafe cells ranked last),
per fold + aggregate Spearman + paired (B - A) 95% CI.

READ:
  - B < A  -> the RES FORMULA hurts ranking even with a good ΔAge -> for ranking use ΔAge
             directly; RES needs redesign (or serves a different purpose).
  - B ~= A -> RES is redundant for ranking; ΔAge suffices.
  - B > A  -> RES formula is sound; Test 7's loss was the model's worse ΔAge -> pair RES with
             the better ΔAge. (Flips the earlier verdict.)

USAGE (repo root, venv active):
    python test7_2_res_isolation.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import DEATH_IDX, LOSS_IDX
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
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(te.X), np.asarray(te.fp, float), sdt.transform(te.dose_time)])
    m = tr.mask
    return Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m]).predict(fte)


def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    paths = ArtifactPaths.of(root)
    tr = gather_split(paths, REGIME, "train")
    te = gather_split(paths, REGIME, "test")
    m = te.mask
    if m.sum() < 3:
        return None

    pred = Predictor(root)
    rows = ModelEstimator(pred).rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    sigma_age = np.array([r["sigma_age"] for r in rows])
    in_dist = np.array([r["in_dist"] for r in rows])

    r_pred = ridge_pred(tr, te)                     # the SAME ΔAge for A and B

    # A = rank by ridge ΔAge directly (higher score = more rejuvenation)
    A = -r_pred
    # B = RES with mu_age = ridge ΔAge (fate/uncertainty/OOD from the model)
    B, _ = compute_res_batch(S, P_loss, r_pred, sigma_age, in_dist, pred.res_params)

    y_age = te.y_age[m]
    y_cls = te.y_cls[m].astype(int)
    unsafe = (y_cls == LOSS_IDX) | (y_cls == DEATH_IDX)
    gated = (-y_age).astype(float)
    if unsafe.any() and (~unsafe).any():
        gated[unsafe] = (-y_age)[~unsafe].min() - 1.0

    return {
        "A_true": spearman(A[m], -y_age),           # vs true ΔAge
        "B_true": spearman(B[m], -y_age),
        "A_safe": spearman(A[m], gated),            # vs safe-rejuvenation
        "B_safe": spearman(B[m], gated),
    }


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


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7.2 — RES formula vs plain ΔAge, SAME ΔAge (ridge) fed to both")
    print("A = rank by ridge ΔAge.  B = RES(ridge ΔAge + model fate/uncertainty/OOD).")
    per_fold = {d: r for d in DONORS if (r := one_fold(d)) is not None}
    if not per_fold:
        print("\n   No folds found (need cellfate_loocv_* or runs/cellfate_loocv_*).")
        return

    for tgt, label in [("true", "vs TRUE ΔAge"), ("safe", "vs SAFE-REJUVENATION (gated)")]:
        rows = [[d, f"{per_fold[d][f'A_{tgt}']:+.3f}", f"{per_fold[d][f'B_{tgt}']:+.3f}",
                 f"{per_fold[d][f'B_{tgt}'] - per_fold[d][f'A_{tgt}']:+.3f}"] for d in per_fold]
        print(f"\n  {label}")
        print(render_table(["fold", "A = ridge ΔAge", "B = RES(ridge ΔAge)", "B − A"],
                           rows, aligns=["l", "r", "r", "r"]))
        aggA = np.nanmean([per_fold[d][f'A_{tgt}'] for d in per_fold])
        aggB = np.nanmean([per_fold[d][f'B_{tgt}'] for d in per_fold])
        diffs = [per_fold[d][f'B_{tgt}'] - per_fold[d][f'A_{tgt}'] for d in per_fold]
        md, (lo, hi), nn = paired_ci(diffs)
        verdict = ("RES(B) BEATS plain ΔAge(A)" if lo > 0
                   else "RES(B) WORSE than plain ΔAge(A)" if hi < 0 else "tied (noise)")
        print(f"   aggregate: A={aggA:.3f}   B={aggB:.3f}")
        print(f"   paired (B − A): mean={md:+.3f}  95% CI=[{lo:+.3f},{hi:+.3f}] (n={nn}) "
              f"-> {verdict}")

    print("\n   WHAT THIS MEANS (same ΔAge in both, so any gap is the RES FORMULA):")
    print("     - B < A -> the RES transform DEGRADES ranking even fed a good ΔAge -> for")
    print("       ranking, use ΔAge directly; RES needs redesign.")
    print("     - B ~= A -> RES is redundant for ranking.")
    print("     - B > A -> RES formula is sound; Test 7's loss was the model's worse ΔAge.")


if __name__ == "__main__":
    main()
