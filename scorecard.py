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
import os
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

# metric -> ("lower"|"higher"|"abs"|"target"|"neutral", pretty label)
#
# STAGE 17. "target" means the metric should APPROACH a per-fold target rather than climb or
# fall. `conformal_coverage` was registered ("higher", ...), but coverage 1.000 is not better
# than 0.900 -- it means the intervals are too wide, which is why `conformal_width` sits at
# 63-81 years. Under "higher is better" a change that simply widened every interval until
# nothing escaped would have scored ACCEPT. Across the committed snapshots 4-5 of 6 folds are
# OVER-covering, so this was not hypothetical.
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
    "conformal_coverage":  ("target", "conformal coverage"),
    # NOT changed to "target": narrower IS genuinely better at equal coverage, and coverage is
    # judged separately above. Noted in plans/STAGE_17_... §17.5 rather than altered.
    "conformal_width":     ("lower",  "interval width"),
    "ood_rate":            ("neutral", "OOD flag rate"),
    "n_cells":             ("neutral", "held-out cells"),
}


# STAGE 17: a "target" metric names the per-fold key holding its target. Read from the fold,
# never hard-coded -- 0.90 is data, and a future run may set a different conformal level.
TARGET_OF = {"conformal_coverage": "conformal_level"}


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
    # STAGE 1.5.3 step 6: `CELLFATE_FOLD_SUFFIX` lets the two arms keep SEPARATE fold builds
    # (`cellfate_loocv_N2_armA` / `_armB`) instead of the second overwriting the first. The first
    # run lost arm A's `scalers.json`, so its deconfounder coefficient had to be reported from a
    # proxy build. Unset -> the historical name, so every existing snapshot still resolves.
    root = resolve_root(f"cellfate_loocv_{donor}{os.environ.get('CELLFATE_FOLD_SUFFIX', '')}")
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

    # Loading the bundle belongs in the SAME error contract as loading the splits above: a
    # missing, incomplete or schema-mismatched bundle is a per-fold condition, not a reason to
    # abandon the snapshot. Outside the try, one un-retrained fold aborted the whole run and
    # discarded the folds that had already succeeded -- which is how a partial retrain (or a
    # single crashed fold) costs you every other measurement.
    try:
        pred = Predictor(root)
    except Exception as exc:  # noqa: BLE001
        return {"_error": repr(exc)[:120]}
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
        # RAW rows, kept so ECE can also be POOLED across folds at aggregation time.
        # `fate_ece` above is a per-fold ECE over ~21 cells in 10 bins, which is biased upward
        # hard enough that a PERFECTLY calibrated model scores 0.183 and clears the 0.169 bar
        # only 26.9% of the time -- it measures the sample size, not the model
        # (audit_metrics.py). Pooling is the more correct LOOCV estimate: every cell is still
        # predicted by a model that never saw it, and the pass rate for a correct system rises
        # to 99.6%. Underscored: raw data, not a metric, so METRICS-driven tables ignore it.
        out["_fate_S"] = [float(v) for v in S]
        out["_fate_y"] = [int(v) for v in st]
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


def _fold_ok(f, key) -> bool:
    if not (isinstance(f, dict) and "_error" not in f and f.get(key) is not None):
        return False
    # STAGE 17: a "target" metric is only usable if its per-fold target is present too. An old
    # snapshot missing `conformal_level` must drop out of that row rather than crash it.
    tk = TARGET_OF.get(key)
    return tk is None or f.get(tk) is not None


def _judged(fold, key, magnitude=False, target=False) -> float:
    """The scalar a row is judged on.

    STAGE 17 adds `target`: distance from the fold's own target, so 1.000 and 0.800 both score
    0.100 against a 0.900 nominal. Under the old "higher is better" the first outranked the
    second, and widening every interval until nothing escaped would have read as ACCEPT.
    """
    v = float(fold[key])
    if target:
        return abs(v - float(fold[TARGET_OF[key]]))
    return abs(v) if magnitude else v


def _agg(folds, key, only=None, magnitude=False, target=False):
    """Mean of `key` across folds.

    STAGE 13. Two parameters, each repairing a defect that inverted a decision:

    `only` restricts to a fold set. Without it the two columns of a comparison each average
    over whatever folds ARE VALID IN THEIR OWN SNAPSHOT, so they can describe different donors
    while sitting side by side. In the C-7 comparison N2 errors out in one snapshot but not the
    other, and 13 of 18 rows printed a 6-fold mean next to a 5-fold mean.

    `magnitude` averages |value| instead of the caller taking |mean|. For a metric whose sign
    VARIES BY DONOR -- level shift is a per-donor bias -- the signed mean measures how far the
    donor panel CANCELS, not how large the error is. It printed 0.230 for a shift whose true
    mean magnitude is 12.72 yr: the founding measurement of Stage 2, rendered as zero.
    """
    names = list(folds) if only is None else only
    vals = [_judged(folds[d], key, magnitude, target)
            for d in names if d in folds and _fold_ok(folds[d], key)]
    return float(np.mean(vals)) if vals else None


def pooled_fate_ece(folds, trials: int = 4000, seed: int = 0):
    """ECE over ALL held-out cells at once, with the floor a perfect model would score.

    Returns None when the snapshot predates this (no `_fate_S`), so old snapshots still load.

    `floor` is the median ECE for a PERFECTLY calibrated model with this exact probability
    vector (`y ~ Bernoulli(p)`, so every bit of it is estimator bias). `excess = ece - floor` is
    the only quantity comparable ACROSS calibrators: raw ECE moves when a calibrator merely
    sharpens, because sharper probabilities sit in extreme bins where Bernoulli variance -- and
    so the floor -- is smaller. Measured on run 3: 75% of one apparent improvement was that.
    """
    S, Y = [], []
    for d in DONORS:
        f = folds.get(d)
        if not isinstance(f, dict) or "_error" in f:
            continue
        s, y = f.get("_fate_S"), f.get("_fate_y")
        if s and y:
            S.extend(s); Y.extend(y)
    if len(S) < 10:
        return None
    S, Y = np.asarray(S, float), np.asarray(Y, float)
    obs = _ece(S, Y)
    rng = np.random.default_rng(seed)
    sims = np.array([_ece(S, (rng.random(len(S)) < S).astype(float)) for _ in range(trials)])
    return {"n": int(len(S)), "ece": float(obs), "floor": float(np.median(sims)),
            "excess": float(obs - np.median(sims)),
            "pctile": float((sims < obs).mean())}


def _print_pooled(snap, label="POOLED"):
    p = pooled_fate_ece(snap["folds"])
    if p is None:
        print(f"\n  {label} fate ECE: n/a (snapshot predates pooled scoring)")
        return None
    print(f"\n  {label} fate ECE over all {p['n']} held-out cells:"
          f"  ECE {p['ece']:.3f}   floor {p['floor']:.3f}   EXCESS {p['excess']:+.3f}"
          f"   (pctile of null {p['pctile']:.1%})")
    return p


def _print_snapshot(snap):
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    folds = snap["folds"]
    ok = [d for d in DONORS if d in folds and "_error" not in folds[d]]
    print(f"\n  SNAPSHOT '{snap['tag']}' — {len(ok)} folds")
    rows = []
    for key, (direction, label) in METRICS.items():
        # STAGE 13: per-fold cells stay SIGNED -- the direction of a donor's shift is real
        # information. Only the aggregate becomes a magnitude, and says so in its label.
        mag = direction == "abs"
        tgt = direction == "target"
        cells = []
        for d in ok:
            v = folds[d].get(key)
            cells.append("n/a" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v)))
        a = _agg(folds, key, magnitude=mag, target=tgt)
        shown = f"|{label}|" if mag else (f"{label} (dist. to target)" if tgt else label)
        rows.append([shown] + cells + ["n/a" if a is None else f"{a:.3f}"])
    print(render_table(["metric"] + ok + ["mean"], rows,
                       aligns=["l"] + ["r"] * (len(ok) + 1)))
    print("  per-fold cells are SIGNED; the mean of a |metric| row is mean(|per-fold|), which for")
    print("  a per-donor bias is the size of the error rather than how far the panel cancels.")
    _print_pooled(snap)


def _common_folds(A_folds, B_folds, key):
    """Folds where BOTH snapshots carry a usable value for `key` — the only set on which a
    paired statistic, or a before/after pair of column means, is defined."""
    return [d for d in DONORS
            if _fold_ok(A_folds.get(d), key) and _fold_ok(B_folds.get(d), key)]


def _paired(A_folds, B_folds, key, magnitude=False, target=False):
    """Per-fold paired differences (B - A) and their 95% CI. This is the accept/reject
    statistic: a change is REAL only if the CI excludes zero.

    STAGE 13: `magnitude` takes |value| PER FOLD before differencing. For an "abs" metric this
    is the difference between judging a real quantity and judging noise -- the old signed form
    read `-28 -> -22` (a 6 yr improvement in magnitude) as a `+6` INCREASE, and then
    `_verdict`'s better_is_down turned that into `REGRESSION`.
    """
    diffs = [_judged(B_folds[d], key, magnitude, target)
             - _judged(A_folds[d], key, magnitude, target)
             for d in _common_folds(A_folds, B_folds, key)]
    n = len(diffs)
    if n < 2:
        return None, (None, None), n
    md = float(np.mean(diffs))
    sd = float(np.std(diffs, ddof=1))
    se = sd / np.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n


def fold_tally(A_folds, B_folds, key, direction, magnitude=False, target=False) -> dict:
    """STAGE 17 -- how many folds moved which way. DESCRIPTIVE ONLY; never a verdict.

    `scorecard.py`'s own header tells the reader to "check the per-fold column before trusting an
    aggregate verdict" -- but only `dage_mae_model` ever printed one, so for 17 of 18 metrics the
    check was impossible. In the Stage 12 comparison `conformal_width` moved the same way on
    ALL FIVE folds (-6.97 yr mean) and still read `noise`, because the paired t is driven by the
    consistency of MAGNITUDE, not of direction. Nothing in the table showed that.

    Deliberately reports COUNTS and no p-value. With 5 paired folds the smallest achievable
    two-sided sign-test p is 2/2**5 = 0.0625, so a unanimous result can NEVER clear 0.05 at this
    n; printing a p-value would invite exactly the "second bite at the apple" this is meant to
    prevent. The accept/reject rule remains the paired 95% CI, unchanged.
    """
    if direction == "neutral":
        return {"better": 0, "worse": 0, "same": 0, "n": 0, "unanimous": False}
    down_is_better = direction in ("lower", "abs", "target")
    better = worse = same = 0
    for d in _common_folds(A_folds, B_folds, key):
        diff = (_judged(B_folds[d], key, magnitude, target)
                - _judged(A_folds[d], key, magnitude, target))
        if diff == 0:
            same += 1
        elif (diff < 0) == down_is_better:
            better += 1
        else:
            worse += 1
    n = better + worse + same
    # Unanimity is only meaningful once there are enough folds to be surprised by; and a row of
    # ties is not agreement about anything.
    unanimous = n >= 4 and (better == n or worse == n)
    return {"better": better, "worse": worse, "same": same, "n": n, "unanimous": unanimous}


def _verdict(direction, md, lo, hi):
    """Pre-committed decision rule: accept only if the paired CI excludes 0 in the
    improving direction. Everything else is noise or a regression."""
    if md is None:
        return "n/a"
    if direction == "neutral":
        return "(context)"
    better_is_down = direction in ("lower", "abs", "target")
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
    print("  NOTE: the paired CI is built on DIFFERENCES, so its sensitivity is set by how")
    print("  CONSISTENT a change is across folds (min detectable mean ~ 1.05 x SD of the effect,")
    print("  6 folds). A uniform change is caught at any size; one that helps some folds and")
    print("  hurts others can be large in the mean and still read as noise -- check the")
    print("  per-fold column before trusting an aggregate verdict.")

    # STAGE 13. Both columns are now averaged over the SAME folds that the paired test uses, so
    # `col_B - col_A == mean diff` exactly, on every row. That identity was false for 13 of 18
    # rows before this change (a 6-fold mean printed beside a 5-fold one), and it is the
    # invariant the tests assert -- it can only hold if the column means and the paired
    # statistic take the magnitude at the same point in the computation.
    rows = []
    for key, (direction, label) in METRICS.items():
        mag = direction == "abs"
        tgt = direction == "target"
        common = _common_folds(A["folds"], B["folds"], key)
        va = _agg(A["folds"], key, only=common, magnitude=mag, target=tgt)
        vb = _agg(B["folds"], key, only=common, magnitude=mag, target=tgt)
        if va is None or vb is None:
            rows.append([label, "n/a", "n/a", "", "", "", "", "n/a"])
            continue
        md, (lo, hi), n = _paired(A["folds"], B["folds"], key, magnitude=mag, target=tgt)
        ci = "" if md is None else f"[{lo:+.3f},{hi:+.3f}]"
        ft = fold_tally(A["folds"], B["folds"], key, direction, magnitude=mag, target=tgt)
        bw = ("" if direction == "neutral"
              else f"{ft['better']}/{ft['worse']}" + ("*" if ft["unanimous"] else ""))
        shown = f"|{label}|" if mag else (f"{label} vs target" if tgt else label)
        rows.append([shown, f"{va:.3f}", f"{vb:.3f}",
                     "" if md is None else f"{md:+.3f}", ci, str(n), bw,
                     _verdict(direction, md, lo, hi)])
        if tgt:
            # Signed, so over- vs under-coverage stays visible. The distance alone cannot say
            # WHICH SIDE of nominal a fold sits on, and that is the actionable half.
            sa = _agg(A["folds"], key, only=common)
            sb = _agg(B["folds"], key, only=common)
            lv = next((A["folds"][d][TARGET_OF[key]] for d in common), None)
            rows.append([f"   ^ raw {label} (target {lv})", f"{sa:.3f}", f"{sb:.3f}",
                         "", "", str(len(common)), "", "(context, never judged)"])
        if mag:
            # The signed mean answers a DIFFERENT question -- "is there a global offset, or do
            # donors cancel?" -- and is worth keeping. It is context, never a verdict: judging
            # it is the defect this stage removes.
            sa = _agg(A["folds"], key, only=common)
            sb = _agg(B["folds"], key, only=common)
            rows.append([f"   ^ signed mean ({label})", f"{sa:+.3f}", f"{sb:+.3f}",
                         "", "", str(len(common)), "", "(context, never judged)"])
    print("\n  AGGREGATE + PAIRED TEST")
    print(render_table(["metric", tag_a, tag_b, "mean diff", "95% CI", "folds", "b/w", "verdict"],
                       rows, aligns=["l", "r", "r", "r", "r", "r", "r", "l"]))
    print("  'b/w' = folds BETTER / WORSE; '*' marks a UNANIMOUS direction across >=4 folds.")
    print("  CONTEXT, never a verdict -- the rule is still the paired 95% CI. A unanimous")
    print("  run of 5 folds is a sign-test p of 0.0625 and can NEVER clear 0.05 at this n,")
    print("  so no p-value is printed: it would invite a second bite at the apple. It exists")
    print("  because the note below asks you to check per-fold agreement, and until now only")
    print("  ONE metric printed a per-fold table.")
    dropped = sorted({d for key in METRICS
                      for d in DONORS
                      if _fold_ok(A["folds"].get(d), key) != _fold_ok(B["folds"].get(d), key)})
    if dropped:
        print(f"\n  NOTE: {', '.join(dropped)} carry a usable value in one snapshot but not the")
        print("  other for at least one metric, and are excluded from that metric's row in BOTH")
        print("  columns. Per-metric fold counts are in the 'folds' column.")
    print("\n  Rows shown as |metric| are judged on the PER-FOLD MAGNITUDE. Level shift is a")
    print("  per-donor bias whose sign varies by donor, so a signed mean measures how far the")
    print("  panel cancels, not how large the error is -- it read 0.230 for a 12.72 yr shift.")

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

    # STAGE 13: the same fold set for all four terms, or the over-approval GAP is a difference
    # between two different donor panels.
    res_folds = [d for d in DONORS
                 if all(_fold_ok(F.get(d), k)
                        for F in (A["folds"], B["folds"])
                        for k in ("res_approvals", "res_approvals_oracle"))]
    ra = _agg(A["folds"], "res_approvals", only=res_folds)
    ro = _agg(A["folds"], "res_approvals_oracle", only=res_folds)
    rb = _agg(B["folds"], "res_approvals", only=res_folds)
    rob = _agg(B["folds"], "res_approvals_oracle", only=res_folds)
    if None not in (ra, ro, rb, rob):
        print(f"\n  RES over-approval (approvals - oracle):  {tag_a}: {ra - ro:+.2f}   "
              f"{tag_b}: {rb - rob:+.2f}   (closer to 0 is better)")
    print("\n  NOTE: 'RES approvals' alone is NOT a quality metric — Test 7.4.3 showed the model")
    print("  approves MORE than the oracle. Judge it by the over-approval gap above.")
    # POOLED fate ECE, with the floor. The per-fold `fate_ece` row above is measured on ~21
    # cells in 10 bins, where a PERFECTLY calibrated model scores 0.183 and clears the 0.169 bar
    # only 26.9% of the time; pooled, the same model clears it 99.6% of the time
    # (audit_metrics.py). Judge the calibration target here, not on the per-fold row.
    print("\n  POOLED fate ECE (the resolvable form of the calibration target):")
    pa_, pb_ = _print_pooled(A, f"  {tag_a}"), _print_pooled(B, f"  {tag_b}")
    if pa_ and pb_:
        print(f"\n    raw ECE   {pa_['ece']:.3f} -> {pb_['ece']:.3f}   "
              f"({pb_['ece'] - pa_['ece']:+.3f})")
        print(f"    EXCESS    {pa_['excess']:+.3f} -> {pb_['excess']:+.3f}   "
              f"({pb_['excess'] - pa_['excess']:+.3f})   <- compare calibrators on THIS")
        print("    (raw ECE also moves when a calibrator merely sharpens, because sharper")
        print("     probabilities lower the floor; excess subtracts that off.)")

    print("\n  ACCEPT the change only if the TARGET metric says ACCEPT and no guard metric says")
    print("  REGRESSION. See MASTER_PLAN.md §7b for the pre-registered criteria per change.")


def main():
    # cp1255 (Hebrew) and other legacy Windows codepages cannot encode the characters this
    # script prints. Without this the script does ALL its work, prints its verdict, and then
    # dies on a print -- BEFORE writing its JSON, so the result is lost. Same guard as
    # plan_tests/verify_stage1_5.py:270; stderr is included too, because a traceback whose
    # source line contains one of those characters would fail the same way.
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
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
