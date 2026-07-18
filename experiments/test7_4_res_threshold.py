"""
Test 7.4 (ΔAge lab notebook) — does recalibration help RES on a NON-rank-invariant metric?

Test 7.3 measured ranking with Spearman, which is RANK-INVARIANT. Platt recalibration is
monotonic (rescales probabilities, preserves order), so RES's rank correlation was unchanged
BY CONSTRUCTION (C-B = +0.000 on 5/6 folds). That does not test whether recalibration helps —
it tests the one axis recalibration cannot move.

Recalibration's real benefit is on VALUE-based (not rank-based) metrics. NOTE: precision@k /
top-k is ALSO rank-based (top-k by score = order), so it is ALSO rank-invariant and would ALSO
show C=B — it does NOT test the point either. The metrics that CAN move under a monotonic
recalibration are:
  (1) SCORE CALIBRATION — bin cells by score VALUE; does score value predict the true
      safe-rejuvenation rate? Recal changes the values, so this can change.
  (2) FIXED-ABSOLUTE-THRESHOLD decision — choose a cutoff on the calib set, apply that VALUE
      to test; recal shifts the distribution so which test cells pass changes -> precision/
      recall of retrieving truly safe-AND-rejuvenating cells changes.

Same A / B / C as 7.3 (ridge ΔAge / RES-raw / RES-recal), same held-out cells, same ΔAge fed to
all. Core comparison is B vs C (only fate calibration differs). A shown for context.

READ:
  - C better than B on calibration / fixed-threshold  -> recal DOES help RES as an actionable
    score; 7.3's null was a metric artifact -> report RES on calibrated-score terms.
  - C ~= B even here  -> recal doesn't help RES on any axis -> 7.3's conclusion is robust.

USAGE (repo root, venv active; needs cellfate_loocv_* bundles + calib/test splits):
    python test7_4_res_threshold.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
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
THRESH_Q = 0.70   # calib-set quantile that defines the "act on it" cutoff (top 30%)


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


def _platt(p_cal, is_pos, p_te) -> np.ndarray:
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


def norm01(x_cal, x_te):
    """Min-max normalize test scores using the CALIB range (value-based, so recal-sensitive)."""
    lo, hi = np.min(x_cal), np.max(x_cal)
    if hi <= lo:
        return np.zeros_like(x_te, dtype=float), 0.5
    z = (np.asarray(x_te, float) - lo) / (hi - lo)
    thr = (np.quantile(x_cal, THRESH_Q) - lo) / (hi - lo)
    return np.clip(z, 0, 1), thr


def score_ece(score01, good, bins=5):
    """Calibration of a [0,1] score vs binary 'good' outcome (value-based)."""
    score01 = np.asarray(score01, float)
    good = np.asarray(good, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (score01 >= edges[i]) & (score01 < hi)
        if m.sum():
            e += (m.sum() / len(score01)) * abs(score01[m].mean() - good[m].mean())
    return float(e)


def prec_recall_at(score01, thr, good):
    """Precision/recall of retrieving 'good' cells when flagging score01 >= thr (fixed value)."""
    flag = score01 >= thr
    if flag.sum() == 0:
        return float("nan"), 0.0
    prec = float(good[flag].mean())
    rec = float(good[flag].sum() / max(good.sum(), 1))
    return prec, rec


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        cal = gather_split(paths, REGIME, "calib")
    except Exception:  # noqa: BLE001
        return None
    if te.mask.sum() < 3 or cal.mask.sum() < 3:
        return None
    pred = Predictor(root)
    est = ModelEstimator(pred)

    def res_scores(sd, S_over=None, Pl_over=None):
        rows = est.rows(sd.X, sd.fp, sd.dose_time)
        S = np.array([r["S"] for r in rows]) if S_over is None else S_over
        Pl = np.array([r["P_loss"] for r in rows]) if Pl_over is None else Pl_over
        sig = np.array([r["sigma_age"] for r in rows])
        ind = np.array([r["in_dist"] for r in rows])
        return S, Pl, sig, ind

    # calib fate + recal maps
    cal_rows = est.rows(cal.X, cal.fp, cal.dose_time)
    S_cal = np.array([r["S"] for r in cal_rows])
    Pl_cal = np.array([r["P_loss"] for r in cal_rows])
    cls_cal = cal.y_cls.astype(int)
    # test fate
    S_te, Pl_te, sig_te, ind_te = res_scores(te)
    S_te_re = _platt(S_cal, cls_cal == SAFE_IDX, S_te)
    Pl_te_re = _platt(Pl_cal, cls_cal == LOSS_IDX, Pl_te)
    S_cal_re = _platt(S_cal, cls_cal == SAFE_IDX, S_cal)
    Pl_cal_re = _platt(Pl_cal, cls_cal == LOSS_IDX, Pl_cal)

    mu_te = ridge_pred(tr, te)
    mu_cal = ridge_pred(tr, cal)
    sig_cal = np.array([r["sigma_age"] for r in cal_rows])
    ind_cal = np.array([r["in_dist"] for r in cal_rows])

    # A/B/C scores on test AND calib (need calib to set the value-threshold + normalization)
    A_te, A_cal = -mu_te, -mu_cal
    B_te, _ = compute_res_batch(S_te, Pl_te, mu_te, sig_te, ind_te, pred.res_params)
    B_cal, _ = compute_res_batch(S_cal, Pl_cal, mu_cal, sig_cal, ind_cal, pred.res_params)
    C_te, _ = compute_res_batch(S_te_re, Pl_te_re, mu_te, sig_te, ind_te, pred.res_params)
    C_cal, _ = compute_res_batch(S_cal_re, Pl_cal_re, mu_cal, sig_cal, ind_cal, pred.res_params)

    m = te.mask
    good = ((te.y_cls[m].astype(int) == SAFE_IDX) & (te.y_age[m] < 0)).astype(float)
    if good.sum() == 0 or good.sum() == len(good):
        return {"degenerate": True}

    out = {"degenerate": False}
    for name, s_te, s_cal in [("A", A_te, A_cal), ("B", B_te, B_cal), ("C", C_te, C_cal)]:
        z_te, thr = norm01(s_cal[cal.mask], s_te[m])
        out[f"{name}_ece"] = score_ece(z_te, good)
        p, r = prec_recall_at(z_te, thr, good)
        out[f"{name}_prec"] = p
        out[f"{name}_rec"] = r
    return out


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 7.4 — does recalibration help RES on VALUE-based metrics? (7.3 used rank-only)")
    print("B=RES raw, C=RES recal, A=ridge ΔAge (context). score calibration (lower better) +")
    print(f"precision/recall of true safe-rejuv cells at the calib top-{int((1-THRESH_Q)*100)}% cutoff.")
    per = {d: r for d in DONORS if (r := one_fold(d)) is not None and not r.get("degenerate")}
    if not per:
        print("\n   No usable folds (need calib+test with safe-rejuv variation).")
        return

    # calibration table
    rows = [[d, f"{per[d]['A_ece']:.3f}", f"{per[d]['B_ece']:.3f}", f"{per[d]['C_ece']:.3f}",
             f"{per[d]['C_ece'] - per[d]['B_ece']:+.3f}"] for d in per]
    print("\n  SCORE CALIBRATION error (lower better)")
    print(render_table(["fold", "A ridge", "B RES raw", "C RES recal", "C−B"],
                       rows, aligns=["l", "r", "r", "r", "r"]))
    dB = [per[d]['C_ece'] - per[d]['B_ece'] for d in per]
    md, (lo, hi), n = paired_ci(dB)
    v = ("recal IMPROVES calibration" if hi < 0 else "recal WORSENS" if lo > 0 else "tied")
    print(f"   paired C−B (calibration): mean={md:+.3f} CI=[{lo:+.3f},{hi:+.3f}] (n={n}) -> {v}")

    # precision/recall table
    rows = [[d, f"{per[d]['B_prec']:.2f}/{per[d]['B_rec']:.2f}",
             f"{per[d]['C_prec']:.2f}/{per[d]['C_rec']:.2f}"] for d in per]
    print("\n  PRECISION/RECALL at fixed calib cutoff (retrieving true safe-rejuv cells)")
    print(render_table(["fold", "B RES raw", "C RES recal"], rows, aligns=["l", "r", "r"]))
    dP = [per[d]['C_prec'] - per[d]['B_prec'] for d in per
          if np.isfinite(per[d]['C_prec']) and np.isfinite(per[d]['B_prec'])]
    md, (lo, hi), n = paired_ci(dP)
    v = ("recal IMPROVES precision" if lo > 0 else "recal WORSENS" if hi < 0 else "tied")
    print(f"   paired C−B (precision): mean={md:+.3f} CI=[{lo:+.3f},{hi:+.3f}] (n={n}) -> {v}")

    print("\n   WHAT THIS MEANS:")
    print("     - C better than B here (calibration lower / precision higher) -> recalibration")
    print("       DOES help RES as an actionable score; 7.3's null was a rank-metric artifact.")
    print("     - C ~= B even on these value-based metrics -> recal genuinely doesn't help RES")
    print("       -> 7.3's conclusion (drop RES for ranking) is robust across metric types.")


if __name__ == "__main__":
    main()
