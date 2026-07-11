"""
Test 7 (ΔAge lab notebook) — does the model's RES RANKING beat ranking-by-ridge-ΔAge?

We proved ridge matches the model on ΔAge magnitude (Test 6). The model's headline is
RANKING (Spearman 0.69) via its RES score (fate + ΔAge + uncertainty + OOD). Does that RES
ranking actually beat simply sorting perturbations by ridge's predicted ΔAge? If not, RES
isn't earning its keep for ranking.

Ranking quality = Spearman(score, -true_ΔAge): a good score ranks the most-rejuvenating
(most-negative ΔAge) cells highest. Per held-out donor we compute THREE rankings on the SAME
cells:
  1. model_RES     — the model's full RES score            (the product's ranking)
  2. model_dAge     — sort by the model's own ΔAge only     (does RES beat model-ΔAge-sort?)
  3. ridge_dAge     — sort by ridge's predicted ΔAge        (the baseline in question)

Report per-fold Spearman, aggregate mean +/- std, and the paired (model_RES - ridge_dAge)
difference across folds with a 95% CI.

READ:
  - model_RES beats ridge_dAge beyond noise -> ranking IS a real contribution over a linear
    baseline. The model earns its keep on ranking.
  - tied -> ranking is driven by ΔAge, which ridge predicts equally well -> ranking is NOT a
    unique contribution over linear (honest, important).
  - model_RES loses -> RES hurts ranking vs a simple ΔAge sort -> reconsider RES for ranking.

USAGE (repo root, venv active):
    python test7_ranking.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.evaluation.metrics import ranking_metrics
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
    """Fit ridge (same feature view as pipeline) on train age-masked cells, predict test."""
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(te.X), np.asarray(te.fp, float), sdt.transform(te.dose_time)])
    m = tr.mask
    reg = Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m])
    return reg.predict(fte)


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    paths = ArtifactPaths.of(root)
    tr = gather_split(paths, REGIME, "train")
    te = gather_split(paths, REGIME, "test")
    m = te.mask
    if m.sum() < 3:
        return None

    pred = Predictor(root)
    est = ModelEstimator(pred)
    rows = est.rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    mu_age = np.array([r["mu_age"] for r in rows])
    sigma_age = np.array([r["sigma_age"] for r in rows])
    in_dist = np.array([r["in_dist"] for r in rows])
    res, _ = compute_res_batch(S, P_loss, mu_age, sigma_age, in_dist, pred.res_params)

    r_pred = ridge_pred(tr, te)
    true = te.y_age

    # ranking_metrics(score, measured): quality = -measured; Spearman(score, quality)
    sp_res = ranking_metrics(res[m], true[m])["spearman"]              # model RES
    sp_mage = ranking_metrics(-mu_age[m], true[m])["spearman"]         # model ΔAge sort
    sp_ridge = ranking_metrics(-r_pred[m], true[m])["spearman"]        # ridge ΔAge sort
    return {"model_RES": sp_res, "model_dAge": sp_mage, "ridge_dAge": sp_ridge}


def paired_ci(diffs):
    n = len(diffs)
    if n < 2:
        return float("nan"), (float("nan"), float("nan"))
    md = sum(diffs) / n
    sd = math.sqrt(sum((x - md) ** 2 for x in diffs) / (n - 1))
    se = sd / math.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se)


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7 — does the model's RES ranking beat ranking-by-ridge-ΔAge?")
    print("Spearman vs true ΔAge (higher = better ordering). same held-out cells for all three.")

    names = ["model_RES", "model_dAge", "ridge_dAge"]
    per_fold = {}
    for d in DONORS:
        r = one_fold(d)
        if r is not None:
            per_fold[d] = r

    if not per_fold:
        print("\n   No folds found (need cellfate_loocv_* or runs/cellfate_loocv_*).")
        return

    rows = []
    for d in DONORS:
        if d in per_fold:
            rows.append([d] + [f"{per_fold[d][n]:+.3f}" for n in names])
        else:
            rows.append([d] + ["n/a"] * 3)
    print("\n" + render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))

    agg = {n: np.nanmean([per_fold[d][n] for d in per_fold]) for n in names}
    std = {n: np.nanstd([per_fold[d][n] for d in per_fold]) for n in names}
    print("\n   aggregate Spearman:  " +
          "   ".join(f"{n}={agg[n]:.3f}±{std[n]:.3f}" for n in names))

    diffs = [per_fold[d]["model_RES"] - per_fold[d]["ridge_dAge"] for d in per_fold]
    md, (lo, hi) = paired_ci(diffs)
    print(f"\n   paired (model_RES − ridge_dAge): mean={md:+.3f}  95% CI=[{lo:+.3f},{hi:+.3f}]")
    wins = sum(1 for x in diffs if x > 0)
    print(f"   model_RES ranks better on {wins}/{len(diffs)} folds")

    print("\n   WHAT THIS MEANS:")
    if lo > 0:
        print("     -> model_RES BEATS ranking-by-ridge-ΔAge (CI above 0). Ranking IS a real")
        print("        contribution over a linear baseline — the model earns its keep on ranking.")
    elif hi < 0:
        print("     -> model_RES is WORSE than ranking-by-ridge-ΔAge (CI below 0). RES hurts")
        print("        ranking vs a simple ΔAge sort — reconsider RES for ranking.")
    else:
        print("     -> TIED (CI includes 0): ranking is driven by ΔAge, which ridge predicts")
        print("        equally well. The model's ranking is NOT a unique contribution over a")
        print("        linear ΔAge sort. Honest, important finding for the writeup.")
    print("   (compare model_RES vs model_dAge to see if RES's fate/uncertainty adds anything")
    print("    over sorting by the model's own ΔAge.)")


if __name__ == "__main__":
    main()
