"""STAGE 16 — the safety floor rejects cells that are demonstrably safe. Why, and is it fixable?

Stage 15 found that >half of held-out cells are `REJECTED_UNSAFE` (`S < tau_safe - 3w = 0.76`)
while `fate_prauc` is 0.965-0.992. A head that RANKS that well while its probabilities sit that
low is the signature of miscalibration, not of danger.

The gate that could have killed this (measured before the plan): 75-100% of held-out cells are
TRUE SAFE. N2 is 19/19 safe with 15 rejected. So the rejections are not correct, and the question
is which of three remaining mechanisms produces them:

  H1  miscalibration  -- probabilities compressed toward the middle
  H3  threshold       -- 0.76 is the wrong bar for a well-ranked head
  H4  prior shift     -- trained on HFF (D0->iPSC, mixed), scored on a donor that is 75-100% safe

H1, H3 and H4 are different diseases with different cures. The decisive separation is fitting the
calibrator TWO ways: on `calib` (deployable) and on `test` (ORACLE, not deployable, an upper
bound). The gap between them IS the measurement of H4 -- if the oracle repairs it and the
deployable one does not, calibration works but does not transfer, and Platt-on-calib is the wrong
fix.

Pooled across folds, because n~20 per fold is too few for a rate and every cell is predicted by a
model that never saw it -- the same argument `scorecard.pooled_fate_ece` already relies on.
Per-fold numbers are reported beside the pooled ones so a single-fold artefact cannot hide.

Read-only. Runs inference on ~20 held-out cells per fold. Trains nothing, writes one results file.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
# STAGE 16 VERIFICATION: the fold set is selectable so the SAME evaluation can be run against the
# soft-fitted bundles (`_s12`) and the hard-refitted ones (`_s16`) without editing the script --
# two runs of one instrument rather than two instruments. Default unchanged.
SUFFIX = os.environ.get("CELLFATE_FOLD_SUFFIX", "_c7t")

OUT = _RESULTS / f"diag_stage16_safety_floor_results{'' if SUFFIX == '_c7t' else SUFFIX}.json"

# Pre-registered in plans/STAGE_16_SAFETY_FLOOR_MISCALIBRATION.md §16.4, before the run.
FALSE_REJECTION_DROP_BAR = 0.50


def _sc():
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def platt(p_fit, y_fit, p_score):
    """Logistic recalibration. Returns `p_score` unchanged when the fit set has one class -- an
    unidentifiable boundary must pass through, not silently produce a constant."""
    from sklearn.linear_model import LogisticRegression
    p_fit = np.asarray(p_fit, float).reshape(-1, 1)
    y_fit = np.asarray(y_fit, int)
    if not (0 < y_fit.sum() < len(y_fit)):
        return np.asarray(p_score, float), False
    lr = LogisticRegression(max_iter=1000).fit(p_fit, y_fit)
    return lr.predict_proba(np.asarray(p_score, float).reshape(-1, 1))[:, 1], True


def split_safe_fraction(paths, regime: str, split: str, safe_idx: int) -> float:
    """T5's train prior, without paying for the expression matrix.

    `gather_split` reads every shard's full 2000-dim X; calling it for `train` costs ~34k x 2000
    floats per fold and dominates the runtime of this script. The class prior needs only `y_cls`,
    so read that column alone. The split lookup is the same one `gather_split` performs
    (`cell_id -> split`), including the pre-Stage-12 collision behaviour, so the fraction
    describes the split the model was actually trained on.
    """
    import pandas as pd

    from cellfate.evaluation.data import load_splits
    assign = load_splits(paths, regime)
    n = safe = 0
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        df = pd.read_parquet(shard, columns=["cell_id", "y_cls"])
        keep = np.fromiter((assign.get(c) == split for c in df["cell_id"]), bool, len(df))
        if not keep.any():
            continue
        soft = np.stack(df["y_cls"].to_numpy()[keep])
        hard = np.argmax(soft, axis=1)
        n += len(hard)
        safe += int((hard == safe_idx).sum())
    return safe / n if n else float("nan")


def confusion(S: np.ndarray, truth: np.ndarray, thr: float) -> dict:
    """The cost at a threshold. BOTH directions, because a calibrator that simply pushes every
    probability up would 'fix' false rejections by destroying the false-approval rate."""
    S, truth = np.asarray(S, float), np.asarray(truth, bool)
    rejected = S < thr
    fr = int((rejected & truth).sum())          # truly safe, rejected -- the cost of interest
    fa = int((~rejected & ~truth).sum())        # truly unsafe, approved -- the trap
    tp = int((~rejected & truth).sum())
    tn = int((rejected & ~truth).sum())
    n_safe, n_unsafe = int(truth.sum()), int((~truth).sum())
    sens = tp / n_safe if n_safe else float("nan")
    spec = tn / n_unsafe if n_unsafe else float("nan")
    bal = np.nanmean([sens, spec])
    return {"false_rejections": fr, "false_approvals": fa,
            "true_approvals": tp, "true_rejections": tn,
            "n_safe": n_safe, "n_unsafe": n_unsafe, "n": int(len(S)),
            "false_rejection_rate": fr / n_safe if n_safe else float("nan"),
            "sensitivity": sens, "specificity": spec,
            "balanced_accuracy": float(bal)}


def best_threshold(S: np.ndarray, truth: np.ndarray) -> dict:
    """T4, ORACLE: the threshold maximising balanced accuracy on the data it is scored on. Bounds
    H3 -- how much of the problem is bar placement rather than calibration."""
    S, truth = np.asarray(S, float), np.asarray(truth, bool)
    if len(np.unique(truth)) < 2:
        # Balanced accuracy is undefined with one class present; returning a "best" threshold
        # here would be fitting a bar to a set that cannot disagree with it.
        return {"oracle_threshold": None, "oracle_balanced_accuracy": None}
    cands = np.unique(np.concatenate([S, [0.0, 1.0]]))
    best, best_t = -1.0, None
    for t in cands:
        b = confusion(S, truth, t)["balanced_accuracy"]
        if b > best:
            best, best_t = b, float(t)
    return {"oracle_threshold": best_t, "oracle_balanced_accuracy": float(best)}


def run(donors=None) -> dict:
    from cellfate.common.constants import SAFE_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.baselines import ModelEstimator
    from cellfate.evaluation.data import gather_split
    from cellfate.inference import Predictor
    sc = _sc()
    per_fold, errs = {}, {}
    pool: dict[str, list] = {"S": [], "S_cal_deployable": [], "S_cal_oracle": [], "truth": []}
    thr = None
    priors = {}
    for d in donors or DONORS:
        root = sc.resolve_root(f"cellfate_loocv_{d}{SUFFIX}")
        try:
            paths = ArtifactPaths.of(root)
            cal = gather_split(paths, "holdout", "calib")
            te = gather_split(paths, "holdout", "test")
            pred = Predictor(root)
        except Exception as exc:                              # noqa: BLE001
            errs[d] = repr(exc)[:120]
            continue
        p = pred.res_params
        thr = float(p.tau_safe - 3.0 * p.w)
        est = ModelEstimator(pred)
        S = np.array([r["S"] for r in est.rows(te.X, te.fp, te.dose_time)])
        S_cal = np.array([r["S"] for r in est.rows(cal.X, cal.fp, cal.dose_time)])
        truth = te.y_cls.astype(int) == SAFE_IDX
        truth_cal = cal.y_cls.astype(int) == SAFE_IDX

        # (a) DEPLOYABLE: fitted on calib only -- the fit set and the score set are disjoint.
        s_dep, fitted_dep = platt(S_cal, truth_cal, S)
        # (b) ORACLE: fitted on the very cells it scores. NOT deployable. Upper bound only.
        s_orc, fitted_orc = platt(S, truth, S)

        # T5: the priors. If train/calib differ from test, H4 has a mechanism, not just a symptom.
        priors[d] = {
            "train_safe_frac": split_safe_fraction(paths, "holdout", "train", SAFE_IDX),
            "calib_safe_frac": float(truth_cal.mean()),
            "test_safe_frac": float(truth.mean()),
        }
        per_fold[d] = {
            "threshold": thr,
            "raw": confusion(S, truth, thr),
            "deployable": confusion(s_dep, truth, thr),
            "oracle": confusion(s_orc, truth, thr),
            "platt_fitted_deployable": bool(fitted_dep),
            "platt_fitted_oracle": bool(fitted_orc),
            "median_S_true_safe": float(np.median(S[truth])) if truth.any() else None,
            "median_S_true_unsafe": float(np.median(S[~truth])) if (~truth).any() else None,
            "priors": priors[d],
        }
        pool["S"].extend(S.tolist())
        pool["S_cal_deployable"].extend(np.asarray(s_dep, float).tolist())
        pool["S_cal_oracle"].extend(np.asarray(s_orc, float).tolist())
        pool["truth"].extend(truth.tolist())

    P = {k: np.asarray(v) for k, v in pool.items()}
    pooled = {}
    if len(P["S"]):
        pooled = {
            "raw": confusion(P["S"], P["truth"], thr),
            "deployable": confusion(P["S_cal_deployable"], P["truth"], thr),
            "oracle": confusion(P["S_cal_oracle"], P["truth"], thr),
            "threshold": thr,
            "oracle_threshold_raw": best_threshold(P["S"], P["truth"]),
            "frac_true_safe_below_threshold": float((P["S"][P["truth"]] < thr).mean()),
            "median_S_true_safe": float(np.median(P["S"][P["truth"]])),
            "median_S_true_unsafe": float(np.median(P["S"][~P["truth"]])),
        }
        fr0 = pooled["raw"]["false_rejections"]
        for arm in ("deployable", "oracle"):
            fr = pooled[arm]["false_rejections"]
            pooled[f"{arm}_fr_drop"] = (fr0 - fr) / fr0 if fr0 else float("nan")
            pooled[f"{arm}_clears_bar"] = bool(fr0 and (fr0 - fr) / fr0 >= FALSE_REJECTION_DROP_BAR)

    return {"folds": per_fold, "pooled": pooled, "errors": errs,
            "bar_false_rejection_drop": FALSE_REJECTION_DROP_BAR,
            "prior_gap_correlation_POST_HOC": prior_gap_correlation(per_fold),
            "verdict": verdict_from(pooled) if pooled else "UNDETERMINED"}


def prior_gap_correlation(folds: dict) -> dict:
    """POST-HOC, EXPLORATORY. Not pre-registered — noticed in the per-fold table after the run.

    If H4's mechanism is real, a fold whose test prior sits FURTHER from the training prior
    should suffer MORE false rejections. Reported with its power stated, because n=6 cannot
    establish it: the two-sided critical Spearman at n=6, alpha=0.05, is 0.886.
    """
    import pandas as pd
    gaps, rates, names = [], [], []
    for d, f in folds.items():
        raw, p = f["raw"], f["priors"]
        if not raw["n_safe"]:
            continue
        names.append(d)
        gaps.append(p["test_safe_frac"] - p["train_safe_frac"])
        rates.append(raw["false_rejection_rate"])
    if len(gaps) < 3:
        return {"n": len(gaps), "spearman": None}
    rho = float(pd.Series(gaps).corr(pd.Series(rates), method="spearman"))
    return {"n": len(gaps), "folds": names, "prior_gaps": gaps,
            "false_rejection_rates": rates, "spearman": rho,
            "critical_rho_n6_alpha05": 0.886,
            "significant": bool(abs(rho) >= 0.886) if len(gaps) == 6 else None,
            "note": "POST-HOC, not pre-registered; n=6 is underpowered"}


def verdict_from(pooled: dict) -> str:
    """The pre-registered reading (plan §16.4), applied mechanically."""
    dep = pooled.get("deployable_clears_bar")
    orc = pooled.get("oracle_clears_bar")
    if dep and orc:
        return "H1 - PLAIN MISCALIBRATION (a deployable calibrator fixes it)"
    if orc and not dep:
        return "H4 - PRIOR/COHORT SHIFT (calibration works but does not transfer from calib)"
    return "H3 OR STRUCTURAL (calibration is not the lever)"


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                     # noqa: BLE001
            pass
    r = run()
    pl = r["pooled"]
    print("\nSTAGE 16 — does the safety floor reject cells that are actually safe, and why?")
    if r["errors"]:
        print(f"  folds unavailable: {r['errors']}")
    if not pl:
        print("  nothing measured")
        return 1

    print(f"\n  threshold = tau_safe - 3w = {pl['threshold']:.2f}")
    print(f"  POOLED n={pl['raw']['n']}  ({pl['raw']['n_safe']} truly safe, "
          f"{pl['raw']['n_unsafe']} truly unsafe)")
    print(f"  median S: true-safe {pl['median_S_true_safe']:.3f}   "
          f"true-unsafe {pl['median_S_true_unsafe']:.3f}   "
          f"-> the head separates; the safe class just sits low")
    print(f"  {pl['frac_true_safe_below_threshold']:.1%} of TRULY SAFE cells fall below the bar")

    print("\n  T1/T3 — the cost, and whether calibration repairs it (pooled)")
    print(f"     {'arm':<24}{'false rej':>10}{'drop':>8}{'false appr':>12}{'sens':>8}"
          f"{'spec':>8}{'bal acc':>9}")
    for arm, label in (("raw", "raw (as shipped)"),
                       ("deployable", "Platt on calib [DEPLOY]"),
                       ("oracle", "Platt on test [ORACLE]")):
        c = pl[arm]
        drop = "" if arm == "raw" else f"{pl[f'{arm}_fr_drop']:>7.1%}"
        print(f"     {label:<24}{c['false_rejections']:>10}{drop:>8}"
              f"{c['false_approvals']:>12}{c['sensitivity']:>8.3f}"
              f"{c['specificity']:>8.3f}{c['balanced_accuracy']:>9.3f}")
    print(f"     bar: false rejections must drop by >= {r['bar_false_rejection_drop']:.0%}"
          f"   deployable {'CLEARS' if pl['deployable_clears_bar'] else 'MISSES'}"
          f"   oracle {'CLEARS' if pl['oracle_clears_bar'] else 'MISSES'}")

    ot = pl["oracle_threshold_raw"]
    print(f"\n  T4 — ORACLE best threshold on the RAW scores: {ot['oracle_threshold']}"
          f"  (balanced acc {ot['oracle_balanced_accuracy']})   shipped bar {pl['threshold']:.2f}")

    print("\n  T5 — class priors (H4's mechanism, if it has one)")
    print(f"     {'fold':<6}{'train safe':>12}{'calib safe':>12}{'test safe':>11}")
    for d, f in r["folds"].items():
        p = f["priors"]
        print(f"     {d:<6}{p['train_safe_frac']:>11.1%}{p['calib_safe_frac']:>12.1%}"
              f"{p['test_safe_frac']:>11.1%}")

    print("\n  PER-FOLD false rejections (raw -> deployable -> oracle), so the pool hides nothing")
    print(f"     {'fold':<6}{'n safe':>8}{'raw':>6}{'deploy':>8}{'oracle':>8}")
    for d, f in r["folds"].items():
        print(f"     {d:<6}{f['raw']['n_safe']:>8}{f['raw']['false_rejections']:>6}"
              f"{f['deployable']['false_rejections']:>8}{f['oracle']['false_rejections']:>8}")

    pg = r["prior_gap_correlation_POST_HOC"]
    if pg.get("spearman") is not None:
        print("\n  POST-HOC (not pre-registered): does the false-rejection rate track how far a")
        print("  fold's test prior sits from the TRAINING prior? -- H4's mechanism, if real")
        print(f"     Spearman {pg['spearman']:+.3f} over n={pg['n']} folds; two-sided critical "
              f"rho at n=6 is {pg['critical_rho_n6_alpha05']} -> "
              f"{'SIGNIFICANT' if pg.get('significant') else 'NOT significant, suggestive only'}")

    print(f"\n  VERDICT: {r['verdict']}")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
