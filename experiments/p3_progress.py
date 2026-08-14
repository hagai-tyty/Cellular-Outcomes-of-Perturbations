"""P3 — does MOLECULAR PROGRESS predict risk better than CALENDAR DAY? (HFF only)

    python experiments/p3_progress.py

READ-ONLY. Writes `results/p3_progress_results.json`. `src/` untouched, no retrain.
Graded against `plans/STAGE_3_P3_PROGRESS_PREREG.md`, committed BEFORE this file existed.

WHY THE METRIC IS LOG-LOSS AND NOT AUC
--------------------------------------
All cells at day 8 share a day. They do NOT share a molecular progress -- reprogramming is
asynchronous. Inside a held-out timepoint the `day` arm therefore predicts a CONSTANT, so AUC would
hand progress an automatic win (AUC 0.5 by construction for a constant), while a per-timepoint
metric would hand day one. Log-loss is a proper scoring rule: a constant predictor that names the
held-out timepoint's true risk exactly scores well, and progress only wins if per-cell variation is
real AND predictable.

LEAKAGE
-------
For each of the 9 leave-one-timepoint-out folds, the scaler, the progress coordinate and the risk
model are ALL fit on training timepoints only; held-out cells are projected into the frozen
coordinate. Nothing derived from the held-out day touches any of it.

MODEL CLASS IS HELD CONSTANT ACROSS THE GATED ARMS
--------------------------------------------------
Every gated arm is the same logistic regression, differing only in which coordinate it is given.
That isolates the INFORMATION in the coordinate rather than the flexibility of the fit. Because
that could handicap `day` -- 8 discrete training values fitted by a monotone function -- an
empirical-interpolation day predictor is reported as a robustness check, and it is the STRONGEST
fair day arm: the observed risk at each training day, interpolated to the held-out day.
"""
from __future__ import annotations

import contextlib
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "p3_progress_results.json"

SUFFIX = "_c7"
BUNDLE = "cellfate_loocv_N2"          # HFF lives in train/val/calib of every fold; one read is enough
LINE = "HFF"
CONTEXT_MAX = 10000                   # cells subsampled for the context arm only
SEED = 0
EPS = 1e-6
T8 = 2.3060041350333704               # t(0.975, df=8) -- nine folds


def log_loss(y, p) -> float:
    p = np.clip(np.asarray(p, float), EPS, 1 - EPS)
    y = np.asarray(y, float)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


def brier(y, p) -> float:
    return float(np.mean((np.asarray(p, float) - np.asarray(y, float)) ** 2))


def fit_predict(Xtr, ytr, Xte, seed: int = SEED):
    """One logistic regression. A single-class training target has no boundary: fall back to the
    base rate, which is the correct constant prediction and is scored honestly by log-loss."""
    ytr = np.asarray(ytr, int)
    if ytr.min() == ytr.max():
        return np.full(len(Xte), float(ytr.mean()))
    m = LogisticRegression(max_iter=1000, random_state=seed)
    m.fit(np.asarray(Xtr, float).reshape(len(Xtr), -1), ytr)
    return m.predict_proba(np.asarray(Xte, float).reshape(len(Xte), -1))[:, 1]


def paired_ci(diffs):
    v = np.asarray([d for d in diffs if np.isfinite(d)], float)
    if len(v) < 2:
        return float("nan"), float("nan"), float("nan"), len(v)
    m = float(v.mean())
    se = float(v.std(ddof=1)) / np.sqrt(len(v))
    t = T8 if len(v) == 9 else 1.96
    return m, m - t * se, m + t * se, len(v)


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.common.constants import DEATH_IDX, LOSS_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split
    install_pretty_console()

    print("\n" + "=" * 92)
    print("P3 — MOLECULAR PROGRESS vs CALENDAR DAY (HFF only, leave-one-timepoint-out)")
    print("=" * 92)
    print("Graded against plans/STAGE_3_P3_PROGRESS_PREREG.md, committed BEFORE this script.")

    paths = ArtifactPaths.of(str(REPO / f"{BUNDLE}{SUFFIX}"))
    parts = [gather_split(paths, "holdout", s) for s in ("train", "val", "calib")]
    line = np.concatenate([np.asarray([str(v) for v in p.cell_line]) for p in parts])
    keep = line == LINE
    X = np.concatenate([np.asarray(p.X, np.float32) for p in parts])[keep]
    logt = np.concatenate([np.asarray(p.dose_time[:, 1], float) for p in parts])[keep]
    cls = np.concatenate([p.y_cls.astype(int) for p in parts])[keep]
    del parts

    day = np.round(np.exp(logt) / 24.0, 1)
    tps = np.unique(np.round(logt, 6))
    heads = {"identity_loss": (cls == LOSS_IDX).astype(int),
             "apoptosis": (cls == DEATH_IDX).astype(int)}
    print(f"\n   {X.shape[0]} {LINE} cells x {X.shape[1]} genes; {len(tps)} timepoints")
    print(render_table(["day", "cells", "P(identity loss)", "P(apoptosis)"],
                       [[f"{day[np.isclose(logt, t)][0]:g}", str(int(np.isclose(logt, t).sum())),
                         f"{heads['identity_loss'][np.isclose(logt, t)].mean():.3f}",
                         f"{heads['apoptosis'][np.isclose(logt, t)].mean():.3f}"] for t in tps],
                       aligns=["r", "r", "r", "r"]))

    out: dict = {"script": "p3_progress", "prereg": "plans/STAGE_3_P3_PROGRESS_PREREG.md",
                 "n_cells": int(X.shape[0]), "n_timepoints": int(len(tps)), "heads": {}}
    rng = np.random.default_rng(SEED)

    for hname, y in heads.items():
        print("\n" + "-" * 92)
        print(f"HEAD: {hname}   (run separately -- no composite endpoint, per P4)")
        print("-" * 92)
        rows, per_arm = [], {a: [] for a in ("day", "progress", "day+progress",
                                             "context", "day_interp")}
        for t in tps:
            te = np.isclose(logt, t)
            tr = ~te
            d_out = float(day[te][0])
            if len(np.unique(y[tr])) < 2:
                rows.append([f"{d_out:g}", str(int(te.sum())), "-", "-", "-", "-", "-", "skip"])
                continue

            # --- everything below is fit on TRAINING cells only -------------------------------
            sc = StandardScaler().fit(X[tr])
            Ztr, Zte = sc.transform(X[tr]), sc.transform(X[te])
            prog_model = Ridge(alpha=1.0).fit(Ztr, logt[tr])       # expression -> day
            ptr = prog_model.predict(Ztr).reshape(-1, 1)
            pte = prog_model.predict(Zte).reshape(-1, 1)           # frozen projection
            dtr = logt[tr].reshape(-1, 1)
            dte = logt[te].reshape(-1, 1)

            preds = {
                "day": fit_predict(dtr, y[tr], dte),
                "progress": fit_predict(ptr, y[tr], pte),
                "day+progress": fit_predict(np.hstack([dtr, ptr]), y[tr], np.hstack([dte, pte])),
            }
            # context arm -- NOT gated, subsampled for cost
            idx = np.flatnonzero(tr)
            sub = rng.choice(idx, min(CONTEXT_MAX, len(idx)), replace=False)
            preds["context"] = fit_predict(sc.transform(X[sub]), y[sub], Zte)
            # strongest fair day arm: empirical risk per training day, interpolated
            tdays = np.unique(np.round(logt[tr], 6))
            emp = np.array([y[np.isclose(logt, u)].mean() for u in tdays])
            preds["day_interp"] = np.full(int(te.sum()), float(np.interp(t, tdays, emp)))

            ll = {k: log_loss(y[te], v) for k, v in preds.items()}
            for k, v in ll.items():
                per_arm[k].append(v)
            rows.append([f"{d_out:g}", str(int(te.sum())), f"{y[te].mean():.3f}",
                         f"{ll['day']:.4f}", f"{ll['progress']:.4f}", f"{ll['day+progress']:.4f}",
                         f"{ll['day_interp']:.4f}", f"{ll['context']:.4f}"])

        print(render_table(["day", "cells", "true risk", "day", "progress", "day+prog",
                            "day interp", "context*"], rows,
                           aligns=["r", "r", "r", "r", "r", "r", "r", "r"]))
        print("   log-loss, lower is better. *context is NOT part of the gate.")

        res = {"per_arm_mean": {k: float(np.mean(v)) for k, v in per_arm.items() if v},
               "per_fold": {k: v for k, v in per_arm.items()}}
        # --- the gate: paired (progress - day) across folds --------------------------------
        dif = [p - d for p, d in zip(per_arm["progress"], per_arm["day"], strict=True)]
        m, lo, hi, n = paired_ci(dif)
        verdict = ("G1 PROGRESS BEATS DAY" if np.isfinite(hi) and hi < 0 else
                   "G3 PROGRESS LOSES" if np.isfinite(lo) and lo > 0 else "G2 TIES")
        dif2 = [p - d for p, d in zip(per_arm["progress"], per_arm["day_interp"], strict=True)]
        m2, lo2, hi2, _ = paired_ci(dif2)
        print(f"\n   means: day={np.mean(per_arm['day']):.4f}  "
              f"progress={np.mean(per_arm['progress']):.4f}  "
              f"day+prog={np.mean(per_arm['day+progress']):.4f}  "
              f"day_interp={np.mean(per_arm['day_interp']):.4f}  "
              f"context={np.mean(per_arm['context']):.4f}")
        print(f"   paired (progress − day):        {m:+.4f}  95% CI [{lo:+.4f},{hi:+.4f}] (n={n})")
        print(f"   paired (progress − day_interp): {m2:+.4f}  95% CI [{lo2:+.4f},{hi2:+.4f}]"
              "   <- robustness, day given its strongest form")
        print(f"   -> {verdict}")
        res["gate"] = {"paired_mean": m, "ci": [lo, hi], "n": n, "verdict": verdict,
                       "vs_day_interp": {"paired_mean": m2, "ci": [lo2, hi2]}}
        out["heads"][hname] = res

    v = {h: out["heads"][h]["gate"]["verdict"] for h in out["heads"]}
    eligible = any(x.startswith("G1") for x in v.values())
    print("\n" + "=" * 92)
    print(render_table(["head", "verdict"], [[k, x] for k, x in v.items()], aligns=["l", "l"]))
    print(f"   P5 ELIGIBLE: {'YES' if eligible else 'NO'}  "
          f"({'G1 on at least one head' if eligible else 'both heads tie or lose -> STOP'})")
    print("   The gate was fixed in the pre-registration and is not adjusted to this result.")
    out["p5_eligible"] = bool(eligible)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
