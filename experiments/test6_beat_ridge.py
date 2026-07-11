"""
Test 6 (ΔAge lab notebook) — can ANY model beat ridge on the REAL ΔAge?

The structural argument (linear clock -> linear ΔAge -> ridge optimal) has a gap: the
model predicts ΔAge from x_pert WITHOUT being given x_ctrl, so it must infer the per-donor
control offset from x_pert — which need not be linear. So a more powerful model COULD beat
ridge on real ΔAge. We test it empirically, on the real Gill data, same leave-one-donor-out
protocol, SAME features as the pipeline's ridge (only the estimator changes):

  ridge (linear)  vs  gradient-boosted trees  vs  random forest  vs  kernel ridge (RBF)

READ:
  - no nonlinear model beats ridge beyond noise -> ridge is at the ceiling; real ΔAge is
    linearly predictable from this input; the neural-net tie is CORRECT -> closed empirically.
  - a nonlinear model beats ridge beyond noise -> real exploitable structure exists -> the
    neural net only tying ridge is a REAL underperformance to diagnose (Test 6.1).

Deterministic, CPU, off-the-shelf. Reads the same LOOCV folds the other scripts use.

USAGE (repo root, venv active):
    python test6_beat_ridge.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.data import gather_split

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
SUBSAMPLE_TRAIN = 20000   # cap train rows for the heavy models (speed); ridge uses all


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def features(Xtr, fptr, dttr, Xte, fpte, dtte):
    """Replicate the pipeline ridge's feature view: standardized X + fp + standardized
    dose_time. Scalers fit on train only."""
    sx = StandardScaler().fit(Xtr)
    sdt = StandardScaler().fit(dttr)
    ftr = np.hstack([sx.transform(Xtr), np.asarray(fptr, float), sdt.transform(dttr)])
    fte = np.hstack([sx.transform(Xte), np.asarray(fpte, float), sdt.transform(dtte)])
    return ftr, fte


def models():
    return {
        "ridge": Ridge(alpha=1.0),
        "boosted_trees": HistGradientBoostingRegressor(max_depth=4, learning_rate=0.05,
                                                       max_iter=300),
        "random_forest": RandomForestRegressor(n_estimators=80, n_jobs=-1, max_depth=12),
        "kernel_rbf": make_pipeline(Nystroem(gamma=None, n_components=300, random_state=0),
                                    Ridge(alpha=1.0)),
    }


def one_fold(donor: str):
    paths = ArtifactPaths.of(resolve_root(f"cellfate_loocv_{donor}"))
    tr = gather_split(paths, REGIME, "train")
    te = gather_split(paths, REGIME, "test")
    mtr, mte = tr.mask, te.mask
    if mtr.sum() < 10 or mte.sum() < 1:
        return None
    ftr, fte = features(tr.X, tr.fp, tr.dose_time, te.X, te.fp, te.dose_time)
    Xtr, ytr = ftr[mtr], tr.y_age[mtr]
    Xte, yte = fte[mte], te.y_age[mte]

    # subsample train for the heavy models only (ridge is cheap, gets all rows)
    rng = np.random.default_rng(0)
    if len(Xtr) > SUBSAMPLE_TRAIN:
        sub = rng.choice(len(Xtr), SUBSAMPLE_TRAIN, replace=False)
    else:
        sub = np.arange(len(Xtr))

    out = {}
    for name, mdl in models().items():
        Xfit, yfit = (Xtr, ytr) if name == "ridge" else (Xtr[sub], ytr[sub])
        mdl.fit(Xfit, yfit)
        out[name] = float(np.abs(mdl.predict(Xte) - yte).mean())
    return out


def paired_ci(diffs):
    n = len(diffs)
    if n < 2:
        return float("nan"), (float("nan"), float("nan"))
    mean_d = sum(diffs) / n
    sd = math.sqrt(sum((x - mean_d) ** 2 for x in diffs) / (n - 1))
    se = sd / math.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return mean_d, (mean_d - t * se, mean_d + t * se)


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 6 — can any nonlinear model beat ridge on REAL ΔAge? (leave-one-donor-out)")
    print("same features as the pipeline ridge; only the estimator changes. lower MAE better.")

    names = list(models().keys())
    per_fold = {}
    for d in DONORS:
        res = one_fold(d)
        if res is not None:
            per_fold[d] = res

    if not per_fold:
        print("\n   No folds found. Ensure cellfate_loocv_* (or runs/cellfate_loocv_*) exist.")
        return

    # per-fold table
    rows = []
    for d in DONORS:
        if d in per_fold:
            r = per_fold[d]
            rows.append([d] + [f"{r[n]:.2f}" for n in names])
        else:
            rows.append([d] + ["n/a"] * len(names))
    print("\n" + render_table(["fold"] + names, rows,
                              aligns=["l"] + ["r"] * len(names)))

    # aggregate + paired-vs-ridge verdict
    agg = {n: np.mean([per_fold[d][n] for d in per_fold]) for n in names}
    print("\n   aggregate MAE:  " + "   ".join(f"{n}={agg[n]:.2f}" for n in names))

    print("\n   IS ANY NONLINEAR MODEL BETTER THAN RIDGE (beyond noise)?")
    beat = False
    for n in names:
        if n == "ridge":
            continue
        diffs = [per_fold[d][n] - per_fold[d]["ridge"] for d in per_fold]
        mean_d, (lo, hi) = paired_ci(diffs)
        verdict = ("BEATS ridge" if hi < 0 else "worse than ridge" if lo > 0 else "tied (noise)")
        if hi < 0:
            beat = True
        print(f"     {n:<14} mean(model−ridge)={mean_d:+.2f}  95% CI=[{lo:+.2f},{hi:+.2f}]  -> {verdict}")

    print("\n   WHAT THIS MEANS:")
    if not beat:
        print("     -> NO nonlinear model beats ridge beyond noise. Ridge is at the ceiling:")
        print("        real ΔAge IS linearly predictable from this input. The neural net tying")
        print("        ridge is CORRECT, not an underperformance. ΔAge question closed EMPIRICALLY.")
    else:
        print("     -> A nonlinear model BEATS ridge -> real exploitable structure exists.")
        print("        The neural net only tying ridge is then a REAL underperformance. REOPEN:")
        print("        Test 6.1 — why can't the net capture what a tree/kernel can? (optimization")
        print("        / multi-task tradeoff / architecture).")


if __name__ == "__main__":
    main()
