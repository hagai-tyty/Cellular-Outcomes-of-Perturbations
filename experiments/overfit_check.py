"""
Overfitting check for ΔAge: compare the model's ΔAge error on data it TRAINED on
vs the held-out donor. This tells us which of three explanations we're in:

  train-MAE << test-MAE   -> OVERFITTING (model memorizes train, fails to generalize)
                             fix = regularization / simpler model  (no input change)
  train-MAE ~= test-MAE   -> NOT overfitting; the model genuinely can't fit better.
                             Points to: signal is linear (nothing to fix) OR training
                             capacity. Next test: the synthetic nonlinear-signal rung.

By default this runs across ALL SIX leave-one-donor-out folds (the cellfate_loocv_*
folders already on disk) and reports the per-fold + aggregate picture -- because a
single fold's overfitting number could be a fluke, and we already know the folds
behave very differently (N2 rejuvenates, N3 doesn't). If the model overfits on some
donors but not others, that's a finding one fold would hide.

USAGE (repo root, venv active):
    python overfit_check.py                 # all 6 LOOCV folds (default)
    python overfit_check.py cellfate_multi  # a single named run instead
"""
from __future__ import annotations

import sys

import numpy as np

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor

REGIME = "holdout"
DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]


def mae(est: ModelEstimator, sd, subset=None) -> tuple[float, int]:
    """Model ΔAge MAE over age-valid cells (optionally a boolean subset mask)."""
    m = sd.mask.copy()
    if subset is not None:
        m &= subset
    if not m.any():
        return float("nan"), 0
    rows = est.rows(sd.X, sd.fp, sd.dose_time)
    age = np.array([r["mu_age"] for r in rows], dtype=np.float64)
    return float(np.abs(age[m] - sd.y_age[m]).mean()), int(m.sum())


def one_run(root: str) -> dict | None:
    """Return {train_gill, train_hff, test} MAEs for one finished run, or None."""
    try:
        paths = ArtifactPaths.of(root)
        est = ModelEstimator(Predictor(root))
        train = gather_split(paths, REGIME, "train")
        test = gather_split(paths, REGIME, "test")
    except Exception as exc:  # noqa: BLE001
        print(f"   ! could not read {root}: {exc}")
        return None
    is_hff = np.array([str(c).upper() == "HFF" for c in train.cell_line])
    return {
        "train_hff": mae(est, train, subset=is_hff),
        "train_gill": mae(est, train, subset=~is_hff),
        "test": mae(est, test),
    }


def resolve_root(name: str) -> str:
    """Find a run folder whether it's at repo root or moved into runs/.
    Lets you reorganize outputs into runs/ without breaking this script."""
    from pathlib import Path
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name  # fall through (will error informatively downstream)


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    # single-run mode if a folder is named explicitly
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        roots = [(sys.argv[1], resolve_root(sys.argv[1]))]
    else:
        roots = [(d, resolve_root(f"cellfate_loocv_{d}")) for d in DONORS]

    print("\nΔAge OVERFITTING CHECK — model error on TRAINING Gill donors vs the HELD-OUT donor")
    print("(train-Gill << test  => memorizing; train-Gill ~= test  => genuinely can't fit better)")

    rows, gaps = [], []
    for label, root in roots:
        r = one_run(root)
        if r is None:
            rows.append([label, "  n/a", "  n/a", "  n/a", "fold not found"])
            continue
        tg, th, te = r["train_gill"][0], r["train_hff"][0], r["test"][0]
        verdict = ("memorizing" if np.isfinite(tg) and np.isfinite(te) and tg < 0.4 * te
                   else "cannot-fit" if np.isfinite(tg) and np.isfinite(te) and tg > 0.75 * te
                   else "partial")
        rows.append([label, f"{th:.2f}", f"{tg:.2f}", f"{te:.2f}", verdict])
        if np.isfinite(tg) and np.isfinite(te):
            gaps.append((tg, te))

    print("\n" + render_table(
        ["held-out donor", "train-HFF", "train-Gill", "TEST (held-out)", "read"],
        rows, aligns=["l", "r", "r", "r", "l"]))
    print("   MAE in years. train-HFF is context (single-cell, huge n). The Gill-vs-test")
    print("   comparison is the overfitting signal for the bulk ΔAge question.")

    # ---- aggregate verdict ------------------------------------------------- #
    if gaps:
        tgs = np.array([g[0] for g in gaps])
        tes = np.array([g[1] for g in gaps])
        print("\n   AGGREGATE across folds:")
        print(f"     mean train-Gill MAE = {tgs.mean():.2f}    mean held-out MAE = {tes.mean():.2f}")
        ratio = tes.mean() / tgs.mean() if tgs.mean() > 0 else float("inf")
        n_memo = int(np.sum(tgs < 0.4 * tes))
        n_cant = int(np.sum(tgs > 0.75 * tes))
        print(f"     folds reading 'memorizing': {n_memo}/{len(gaps)}   "
              f"'cannot-fit': {n_cant}/{len(gaps)}")
        print("\n   WHAT THIS MEANS:")
        if n_memo >= n_cant and n_memo >= len(gaps) / 2:
            print(f"     -> OVERFITTING dominates ({ratio:.1f}x worse on held-out). The model memorizes")
            print("        the training donors. Fix is in the MODEL (regularization / simpler / more")
            print("        ensembling) — NOT the input. This is a real, fixable model problem.")
        elif n_cant >= len(gaps) / 2:
            print("     -> NOT overfitting on most folds: the model fits training donors about as")
            print("        poorly as the held-out one. It genuinely can't fit ΔAge better on this")
            print("        input. Next: is the signal linear (nothing to fix) or is training the")
            print("        limit? -> run the synthetic nonlinear-signal rung to decide.")
        else:
            print("     -> MIXED: overfits on some donors, not others (matches the known donor")
            print("        heterogeneity). Both a model-regularization tweak AND the signal-linearity")
            print("        question are worth testing.")


if __name__ == "__main__":
    main()

