"""
Test 7.1 (ΔAge lab notebook) — is RES actually WORSE, or was Test 7 scoring it against the
wrong target?

Test 7 scored rankings against TRUE ΔAge only, and RES lost badly (0.69 vs 0.95). But RES is
designed to rank by SAFE REJUVENATION (rejuvenating AND safe), deliberately down-weighting
cells that rejuvenate but are unsafe. So Test 7 may have penalized RES for doing its job.

Here we score the SAME three rankings (model_RES, model_dAge, ridge_dAge) against a
SAFE-REJUVENATION ground truth built from TRUE fate + TRUE ΔAge on the held-out cells:
  - a good outcome = truly SAFE and strongly rejuvenating (negative ΔAge);
  - truly UNSAFE cells (loss/death) are bad regardless of how much they "rejuvenate".
Two reasonable definitions (to avoid gaming one):
  gated     : quality = -ΔAge for safe cells; unsafe cells ranked at the bottom.
  penalized : quality = -ΔAge - 3*std(ΔAge) for unsafe cells (soft penalty).
Plus precision@k: of the top-k by score, how many are truly safe-AND-rejuvenating.

We also report the FATE COMPOSITION of each held-out donor's age-valid cells — if there are
no unsafe cells, safety can't matter and the test degenerates to Test 7 (noted).

READ:
  - RES beats the ΔAge sorts against safe-rejuvenation -> RES does its job; Test 7 used the
    wrong target; report the right metric. (Requires the model's fate predictions to align
    with true safety.)
  - RES still loses even against safe-rejuvenation -> RES genuinely doesn't help ranking
    (its fate weighting doesn't align with true safety) -> rank by ΔAge; reconsider RES.

USAGE (repo root, venv active):
    python test7_1_safe_ranking.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import DEATH_IDX, LOSS_IDX, SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
K_TOP = 5


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


def safe_rejuv_targets(y_age, y_cls):
    """Two 'safe rejuvenation' quality targets. Higher = better outcome."""
    unsafe = (y_cls == LOSS_IDX) | (y_cls == DEATH_IDX)
    base = -y_age                                   # rejuvenation (more negative ΔAge = better)
    # gated: unsafe cells pushed below every safe cell
    gated = base.copy().astype(float)
    if unsafe.any() and (~unsafe).any():
        gated[unsafe] = base[~unsafe].min() - 1.0
    # penalized: soft penalty proportional to ΔAge spread
    pen = base - (3.0 * (np.std(y_age) + 1e-9)) * unsafe
    return gated, pen, unsafe


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
    mu_age = np.array([r["mu_age"] for r in rows])
    sigma_age = np.array([r["sigma_age"] for r in rows])
    in_dist = np.array([r["in_dist"] for r in rows])
    res, _ = compute_res_batch(S, P_loss, mu_age, sigma_age, in_dist, pred.res_params)
    r_pred = ridge_pred(tr, te)

    y_age = te.y_age[m]
    y_cls = te.y_cls[m].astype(int)
    scores = {"model_RES": res[m], "model_dAge": -mu_age[m], "ridge_dAge": -r_pred[m]}
    gated, pen, unsafe = safe_rejuv_targets(y_age, y_cls)

    out = {"n": int(m.sum()),
           "n_safe": int((y_cls == SAFE_IDX).sum()),
           "n_unsafe": int(unsafe.sum())}
    for name, sc in scores.items():
        out[f"{name}|gated"] = spearman(sc, gated)
        out[f"{name}|pen"] = spearman(sc, pen)
        # precision@k: of top-k by score, fraction truly safe AND rejuvenating (ΔAge<0)
        good = (y_cls == SAFE_IDX) & (y_age < 0)
        k = min(K_TOP, len(sc))
        topk = np.argsort(-sc)[:k]
        out[f"{name}|p@{K_TOP}"] = float(good[topk].mean()) if k else float("nan")
    return out


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

    print("\nTEST 7.1 — rank against SAFE REJUVENATION (RES's actual objective), not raw ΔAge")
    per_fold = {d: r for d in DONORS if (r := one_fold(d)) is not None}
    if not per_fold:
        print("\n   No folds found (need cellfate_loocv_* or runs/cellfate_loocv_*).")
        return

    # fate composition
    comp = [[d, f"{per_fold[d]['n']}", f"{per_fold[d]['n_safe']}", f"{per_fold[d]['n_unsafe']}"]
            for d in per_fold]
    print("\n" + render_table(["fold", "age-valid cells", "safe", "unsafe (loss/death)"],
                              comp, aligns=["l", "r", "r", "r"]))

    names = ["model_RES", "model_dAge", "ridge_dAge"]
    for tgt, label in [("gated", "GATED (unsafe ranked last)"),
                       ("pen", "PENALIZED (soft unsafe penalty)")]:
        rows = [[d] + [f"{per_fold[d][f'{n}|{tgt}']:+.3f}" for n in names] for d in per_fold]
        print("\n  safe-rejuv target = " + label)
        print(render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
        agg = {n: np.nanmean([per_fold[d][f'{n}|{tgt}'] for d in per_fold]) for n in names}
        print("   aggregate: " + "   ".join(f"{n}={agg[n]:.3f}" for n in names))
        diffs = [per_fold[d][f'model_RES|{tgt}'] - per_fold[d][f'ridge_dAge|{tgt}']
                 for d in per_fold]
        md, (lo, hi), nn = paired_ci(diffs)
        verdict = ("RES BEATS ΔAge-sort" if lo > 0 else "RES WORSE than ΔAge-sort"
                   if hi < 0 else "tied (noise)")
        print(f"   paired (model_RES − ridge_dAge): mean={md:+.3f} 95% CI=[{lo:+.3f},{hi:+.3f}]"
              f" (n={nn}) -> {verdict}")

    # precision@k
    rows = [[d] + [f"{per_fold[d][f'{n}|p@{K_TOP}']:.2f}" for n in names] for d in per_fold]
    print(f"\n  precision@{K_TOP}: of top-{K_TOP} by score, fraction truly SAFE & rejuvenating")
    print(render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
    aggp = {n: np.nanmean([per_fold[d][f'{n}|p@{K_TOP}'] for d in per_fold]) for n in names}
    print("   aggregate: " + "   ".join(f"{n}={aggp[n]:.2f}" for n in names))

    print("\n   WHAT THIS MEANS:")
    print("     - If any unsafe-cell count is ~0, safety can't matter -> this reduces to Test 7")
    print("       (RES loses because it only adds noise to a clean ΔAge signal).")
    print("     - If RES beats the ΔAge sorts here (but lost in Test 7) -> RES optimizes safe")
    print("       rejuvenation, Test 7 used the wrong target -> report THIS metric.")
    print("     - If RES loses here too -> its fate weighting doesn't align with TRUE safety")
    print("       (model fate predictions off, or RES mis-designed) -> rank by ΔAge; fix/redesign RES.")


if __name__ == "__main__":
    main()
