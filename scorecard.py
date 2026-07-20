"""
CellFate-Rx SCORECARD — freeze a full metric snapshot, then diff snapshots across changes.

Purpose: every future change to the model must be judged against the SAME battery of numbers,
per fold, so we can see exactly what improved and what regressed. No more "it feels better".

    python scorecard.py snapshot --tag baseline        # run everything, save scorecard/baseline.json
    python scorecard.py compare baseline after_recal   # diff two snapshots, per metric per fold
    python scorecard.py list                           # show saved snapshots

WHAT IT MEASURES (per leave-one-donor-out fold, plus aggregate):

  ΔAge          dage_mae_model, dage_mae_ridge      (lower better)
                level_shift_model, level_shift_ridge = med(pred) - med(true)   (|.| lower better)
                  ^ Test 7.4.3's core finding: per-donor level shift, +-12.7 yr, cancels on average
  RANKING       rank_res, rank_model_dage, rank_ridge_dage  (Spearman vs true ΔAge, higher better)
  FATE          fate_prauc, fate_roc (higher better), fate_ece, fate_ece_platt (lower better)
  RES           res_approvals, res_approvals_oracle  (composition matters -- see notes)
                res_median, res_max                  (Test 7.4.2: raw RES collapses to ~0)
  UNCERTAINTY   conformal_coverage vs conformal_level, interval_width   <- NEVER VALIDATED BEFORE
  OOD           ood_rate = fraction of held-out cells flagged out-of-distribution  <- LIKEWISE

Uncertainty and OOD have never been isolated in any test; they are included here so the baseline
captures them before any code changes.

NOTE ON res_approvals: more is NOT better. Test 7.4.3 showed the model approves MORE cells than
the oracle (14 vs 11) and N3 approves 7 where truth says 0. The meaningful quantity is approvals
RELATIVE to oracle, so compare shows both and flags over-approval.

USAGE (repo root, venv active). Snapshots are written to scorecard/<tag>.json.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
SNAP_DIR = Path("scorecard")
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
APPROVED = "APPROVED"

# metric -> ("lower"|"higher"|"neutral", pretty label)
METRICS = {
    "dage_mae_model":      ("lower",  "ΔAge MAE (model)"),
    "dage_mae_ridge":      ("lower",  "ΔAge MAE (ridge)"),
    "level_shift_model":   ("abs",    "level shift (model)"),
    "level_shift_ridge":   ("abs",    "level shift (ridge)"),
    "rank_res":            ("higher", "rank: RES"),
    "rank_model_dage":     ("higher", "rank: model ΔAge"),
    "rank_ridge_dage":     ("higher", "rank: ridge ΔAge"),
    "fate_prauc":          ("higher", "fate PR-AUC"),
    "fate_roc":            ("higher", "fate ROC-AUC"),
    "fate_ece":            ("lower",  "fate ECE"),
    "fate_ece_platt":      ("lower",  "fate ECE (Platt)"),
    "res_median":          ("neutral", "RES median"),
    "res_max":             ("neutral", "RES max"),
    "res_approvals":       ("neutral", "RES approvals"),
    "res_approvals_oracle": ("neutral", "RES approvals (oracle)"),
    "conformal_coverage":  ("higher", "conformal coverage"),
    "conformal_width":     ("lower",  "interval width"),
    "ood_rate":            ("neutral", "OOD flag rate"),
    "n_cells":             ("neutral", "held-out cells"),
}


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def _ridge(tr, targets):
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    reg = Ridge(alpha=1.0).fit(ftr[tr.mask], tr.y_age[tr.mask])
    out = []
    for s in targets:
        f = np.hstack([sx.transform(s.X), np.asarray(s.fp, float), sdt.transform(s.dose_time)])
        out.append(reg.predict(f))
    return out


def _platt(p_cal, is_pos, p_te):
    p_cal = np.asarray(p_cal, float).reshape(-1, 1)
    is_pos = np.asarray(is_pos, int)
    if not (0 < is_pos.sum() < len(is_pos)):
        return np.asarray(p_te, float)
    lr = LogisticRegression(max_iter=1000).fit(p_cal, is_pos)
    return lr.predict_proba(np.asarray(p_te, float).reshape(-1, 1))[:, 1]


def _ece(p, y, bins=10):
    p, y = np.asarray(p, float), np.asarray(y, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (p >= edges[i]) & (p < hi)
        if m.sum():
            e += (m.sum() / len(p)) * abs(p[m].mean() - y[m].mean())
    return float(e)


def _sp(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(spearmanr(a, b).correlation)


def measure_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        cal = gather_split(paths, REGIME, "calib")
    except Exception as exc:  # noqa: BLE001
        return {"_error": repr(exc)[:120]}
    m = te.mask
    if m.sum() < 3:
        return {"_error": "too few age-valid cells"}

    pred = Predictor(root)
    est = ModelEstimator(pred)
    p = pred.res_params

    rows = est.rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    mu = np.array([r["mu_age"] for r in rows])
    sig = np.array([r["sigma_age"] for r in rows])
    ind = np.array([r["in_dist"] for r in rows])

    crows = est.rows(cal.X, cal.fp, cal.dose_time)
    S_cal = np.array([r["S"] for r in crows])
    cls_cal = cal.y_cls.astype(int)

    (r_te,) = _ridge(tr, [te])
    y = te.y_age[m]
    out = {"n_cells": int(m.sum())}

    # ---- ΔAge ----
    out["dage_mae_model"] = float(np.abs(mu[m] - y).mean())
    out["dage_mae_ridge"] = float(np.abs(r_te[m] - y).mean())
    out["level_shift_model"] = float(np.median(mu[m]) - np.median(y))
    out["level_shift_ridge"] = float(np.median(r_te[m]) - np.median(y))

    # ---- ranking (vs true ΔAge; higher score = more rejuvenation) ----
    res, stat = compute_res_batch(S, P_loss, mu, sig, ind, p)
    out["rank_res"] = _sp(res[m], -y)
    out["rank_model_dage"] = _sp(-mu[m], -y)
    out["rank_ridge_dage"] = _sp(-r_te[m], -y)

    # ---- fate ----
    st = (te.y_cls.astype(int) == SAFE_IDX).astype(int)
    if 0 < st.sum() < len(st):
        out["fate_prauc"] = float(average_precision_score(st, S))
        out["fate_roc"] = float(roc_auc_score(st, S))
        out["fate_ece"] = _ece(S, st)
        out["fate_ece_platt"] = _ece(_platt(S_cal, cls_cal == SAFE_IDX, S), st)
    else:
        for k in ("fate_prauc", "fate_roc", "fate_ece", "fate_ece_platt"):
            out[k] = None

    # ---- RES ----
    out["res_median"] = float(np.median(res[m]))
    out["res_max"] = float(np.max(res[m]))
    out["res_approvals"] = int((np.asarray(stat)[m] == APPROVED).sum())
    _, stat_o = compute_res_batch(S, P_loss, np.where(m, te.y_age, mu), sig, ind, p)
    out["res_approvals_oracle"] = int((np.asarray(stat_o)[m] == APPROVED).sum())

    # ---- uncertainty (conformal) — NEVER VALIDATED BEFORE ----
    lo, hi = mu - pred.q, mu + pred.q
    out["conformal_coverage"] = float(((y >= lo[m]) & (y <= hi[m])).mean())
    out["conformal_level"] = float(pred.conformal_level)
    out["conformal_width"] = float(2.0 * pred.q)

    # ---- OOD — NEVER VALIDATED BEFORE ----
    out["ood_rate"] = float((~ind[m]).mean())
    return out


def cmd_snapshot(tag: str):
    SNAP_DIR.mkdir(exist_ok=True)
    print(f"\nSCORECARD snapshot '{tag}' — measuring {len(DONORS)} folds...")
    folds = {}
    for d in DONORS:
        r = measure_fold(d)
        folds[d] = r
        status = r.get("_error", "ok")
        print(f"   {d}: {status}")
    snap = {
        "tag": tag,
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "folds": folds,
    }
    path = SNAP_DIR / f"{tag}.json"
    path.write_text(json.dumps(snap, indent=2), encoding="utf-8")
    print(f"\n   saved -> {path}")
    _print_snapshot(snap)


def _agg(folds, key):
    vals = [f[key] for f in folds.values()
            if isinstance(f, dict) and f.get(key) is not None and "_error" not in f]
    return float(np.mean(vals)) if vals else None


def _print_snapshot(snap):
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    folds = snap["folds"]
    ok = [d for d in DONORS if d in folds and "_error" not in folds[d]]
    print(f"\n  SNAPSHOT '{snap['tag']}' — {len(ok)} folds")
    rows = []
    for key, (_, label) in METRICS.items():
        cells = []
        for d in ok:
            v = folds[d].get(key)
            cells.append("n/a" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v)))
        a = _agg(folds, key)
        rows.append([label] + cells + ["n/a" if a is None else f"{a:.3f}"])
    print(render_table(["metric"] + ok + ["mean"], rows,
                       aligns=["l"] + ["r"] * (len(ok) + 1)))


def _paired(A_folds, B_folds, key):
    """Per-fold paired differences (B - A) and their 95% CI. This is the accept/reject
    statistic: a change is REAL only if the CI excludes zero."""
    diffs = []
    for d in DONORS:
        fa, fb = A_folds.get(d), B_folds.get(d)
        if not isinstance(fa, dict) or not isinstance(fb, dict):
            continue
        if "_error" in fa or "_error" in fb:
            continue
        va, vb = fa.get(key), fb.get(key)
        if va is None or vb is None:
            continue
        diffs.append(float(vb) - float(va))
    n = len(diffs)
    if n < 2:
        return None, (None, None), n
    md = float(np.mean(diffs))
    sd = float(np.std(diffs, ddof=1))
    se = sd / np.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n


def _verdict(direction, md, lo, hi):
    """Pre-committed decision rule: accept only if the paired CI excludes 0 in the
    improving direction. Everything else is noise or a regression."""
    if md is None:
        return "n/a"
    if direction == "neutral":
        return "(context)"
    better_is_down = direction in ("lower", "abs")
    if lo > 0:                       # significantly increased
        return "REGRESSION" if better_is_down else "ACCEPT (better)"
    if hi < 0:                       # significantly decreased
        return "ACCEPT (better)" if better_is_down else "REGRESSION"
    return "noise (CI incl. 0)"


def cmd_compare(tag_a: str, tag_b: str):
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    pa, pb = SNAP_DIR / f"{tag_a}.json", SNAP_DIR / f"{tag_b}.json"
    for x in (pa, pb):
        if not x.exists():
            print(f"   missing snapshot: {x}")
            return
    A = json.loads(pa.read_text())
    B = json.loads(pb.read_text())
    print(f"\nSCORECARD compare:  {tag_a}  ->  {tag_b}")
    print("  DECISION RULE: a change is REAL only if the paired 95% CI across folds excludes 0.")
    print("  'noise (CI incl. 0)' means the change is not distinguishable from fold variation.")

    rows = []
    for key, (direction, label) in METRICS.items():
        va, vb = _agg(A["folds"], key), _agg(B["folds"], key)
        if va is None or vb is None:
            rows.append([label, "n/a", "n/a", "", "", "n/a"])
            continue
        md, (lo, hi), n = _paired(A["folds"], B["folds"], key)
        if direction == "abs":       # judge |level shift|, not signed
            va, vb = abs(va), abs(vb)
        ci = "" if md is None else f"[{lo:+.3f},{hi:+.3f}]"
        rows.append([label, f"{va:.3f}", f"{vb:.3f}",
                     "" if md is None else f"{md:+.3f}", ci,
                     _verdict(direction, md, lo, hi)])
    print("\n  AGGREGATE + PAIRED TEST (n folds)")
    print(render_table(["metric", tag_a, tag_b, "mean diff", "95% CI", "verdict"], rows,
                       aligns=["l", "r", "r", "r", "r", "l"]))

    print("\n  PER-FOLD ΔAge MAE (model) — where did the change land?")
    ok = [d for d in DONORS if d in A["folds"] and d in B["folds"]
          and "_error" not in A["folds"][d] and "_error" not in B["folds"][d]]
    rows = []
    for d in ok:
        a, b = A["folds"][d].get("dage_mae_model"), B["folds"][d].get("dage_mae_model")
        if a is None or b is None:
            continue
        rows.append([d, f"{a:.2f}", f"{b:.2f}", f"{b - a:+.2f}",
                     "+ better" if b < a else ("- worse" if b > a else "same")])
    print(render_table(["fold", tag_a, tag_b, "delta", "verdict"], rows,
                       aligns=["l", "r", "r", "r", "l"]))

    ra = _agg(A["folds"], "res_approvals")
    ro = _agg(A["folds"], "res_approvals_oracle")
    rb = _agg(B["folds"], "res_approvals")
    rob = _agg(B["folds"], "res_approvals_oracle")
    if None not in (ra, ro, rb, rob):
        print(f"\n  RES over-approval (approvals - oracle):  {tag_a}: {ra - ro:+.2f}   "
              f"{tag_b}: {rb - rob:+.2f}   (closer to 0 is better)")
    print("\n  NOTE: 'RES approvals' alone is NOT a quality metric — Test 7.4.3 showed the model")
    print("  approves MORE than the oracle. Judge it by the over-approval gap above.")
    print("\n  ACCEPT the change only if the TARGET metric says ACCEPT and no guard metric says")
    print("  REGRESSION. See MASTER_PLAN.md §7b for the pre-registered criteria per change.")


def main():
    ap = argparse.ArgumentParser(description="CellFate-Rx metric scorecard")
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("snapshot")
    s.add_argument("--tag", required=True)
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    sub.add_parser("list")
    args = ap.parse_args()
    if args.cmd == "snapshot":
        cmd_snapshot(args.tag)
    elif args.cmd == "compare":
        cmd_compare(args.a, args.b)
    elif args.cmd == "list":
        SNAP_DIR.mkdir(exist_ok=True)
        snaps = sorted(SNAP_DIR.glob("*.json"))
        print("\n  saved snapshots:" if snaps else "\n  no snapshots yet")
        for x in snaps:
            j = json.loads(x.read_text())
            print(f"   {j['tag']:<20} {j.get('utc', '')}")
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
