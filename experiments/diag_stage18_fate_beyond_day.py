"""STAGE 18 -- is the fate head predicting biology, or reading a clock?

`fate_prauc` is 0.958 and `fate_roc` 0.954 on held-out donors. But `dose_time` -- a direct MODEL
INPUT -- encodes the timepoint, and the record already warned that "day > 31 -> unsafe" gets most
of the way there. This measures how much of that 0.958 survives once the timepoint is known.

THE STRUCTURE THAT MAKES THIS URGENT. On the held-out cells, fate is almost a function of
timepoint: for N2, O1 and O2, ZERO timepoints carry more than one class; N3 and Y2 have one each;
only Y1 has several. A lookup table on `time_h` therefore scores near-perfectly on four of six
donors without using a single gene.

THE DECISIVE TEST IS STRATIFIED, NOT MARGINAL. Any marginal metric (PR-AUC, ROC-AUC) is inflated
by the timepoint the model was handed. The question "does expression add anything?" is only asked
by comparing cells that SHARE a timepoint -- there, the day input is constant and cannot help. So
the headline number here is a WITHIN-TIMEPOINT concordance: over all pairs (safe, unsafe) drawn
from the same timepoint, how often does the model rank the safe one higher?

That set is small (a handful of strata) and the answer is reported with a permutation null and an
exact pair count, because a concordance over a dozen pairs is not a result unless it is stated
with its own power.

Read-only. Runs inference on ~20 held-out cells per fold. Trains only cheap sklearn baselines.
"""
from __future__ import annotations

import itertools
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
SUFFIX = os.environ.get("CELLFATE_FOLD_SUFFIX", "_s16")
OUT = _RESULTS / f"diag_stage18_fate_beyond_day_results{SUFFIX}.json"
N_PERM = 20000
SEED = 18


def _sc():
    import importlib.util
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def stratified_pairs(score, safe, stratum):
    """Every (safe, unsafe) pair drawn from the SAME timepoint.

    This is the whole point: inside a stratum the `dose_time` input is constant, so anything the
    model gets right here it got from expression. Returns (n_concordant, n_tied, n_pairs).
    """
    score, safe, stratum = np.asarray(score, float), np.asarray(safe, bool), np.asarray(stratum)
    conc = tied = total = 0
    for u in np.unique(stratum):
        m = stratum == u
        pos, neg = score[m & safe], score[m & ~safe]
        for a, b in itertools.product(pos, neg):
            total += 1
            if a > b:
                conc += 1
            elif a == b:
                tied += 1
    return conc, tied, total


def stratified_auc(score, safe, stratum):
    c, t, n = stratified_pairs(score, safe, stratum)
    return ((c + 0.5 * t) / n) if n else float("nan"), n


def perm_null(score, safe, stratum, n_perm=N_PERM, seed=SEED):
    """Shuffle the scores WITHIN each stratum, preserving the stratum sizes and the label layout.

    A global shuffle would also destroy the between-timepoint structure and make the null far too
    easy to beat; permuting within strata tests exactly the claim being made.
    """
    rng = np.random.default_rng(seed)
    score = np.asarray(score, float).copy()
    obs, n_pairs = stratified_auc(score, safe, stratum)
    if not n_pairs:
        return {"observed": None, "p": None, "n_pairs": 0}
    strata = [np.flatnonzero(stratum == u) for u in np.unique(stratum)]
    ge = 0
    for _ in range(n_perm):
        s = score.copy()
        for idx in strata:
            s[idx] = rng.permutation(s[idx])
        if stratified_auc(s, safe, stratum)[0] >= obs:
            ge += 1
    return {"observed": float(obs), "p": float((ge + 1) / (n_perm + 1)), "n_pairs": int(n_pairs)}


def baselines(tr, te, safe_tr, safe_te):
    """Marginal PR-AUC/ROC for three nested feature sets, fitted on train, scored on the held-out
    donor. `time` is the timepoint column of `dose_time`."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    out = {}
    feats = {
        "time_only": (tr.dose_time[:, 1:2], te.dose_time[:, 1:2]),
        "expr_only": (tr.X, te.X),
        "expr_plus_time": (np.hstack([tr.X, tr.dose_time[:, 1:2]]),
                           np.hstack([te.X, te.dose_time[:, 1:2]])),
    }
    for name, (ftr, fte) in feats.items():
        if len(set(safe_tr.tolist())) < 2 or len(set(safe_te.tolist())) < 2:
            out[name] = {"prauc": None, "roc": None}
            continue
        sx = StandardScaler().fit(ftr)
        lr = LogisticRegression(max_iter=2000).fit(sx.transform(ftr), safe_tr)
        p = lr.predict_proba(sx.transform(fte))[:, 1]
        out[name] = {"prauc": float(average_precision_score(safe_te, p)),
                     "roc": float(roc_auc_score(safe_te, p)),
                     "_score": [float(v) for v in p]}
    return out


def run(donors=None) -> dict:
    from sklearn.metrics import average_precision_score, roc_auc_score

    from cellfate.common.constants import SAFE_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.baselines import ModelEstimator
    from cellfate.evaluation.data import gather_split
    from cellfate.inference import Predictor
    sc = _sc()
    per_fold, errs = {}, {}
    pool = {"S": [], "safe": [], "stratum": [], "time_only": []}
    for d in donors or DONORS:
        root = sc.resolve_root(f"cellfate_loocv_{d}{SUFFIX}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, "holdout", "train")
            te = gather_split(paths, "holdout", "test")
            pred = Predictor(root)
        except Exception as exc:                              # noqa: BLE001
            errs[d] = repr(exc)[:120]
            continue
        S = np.array([r["S"] for r in ModelEstimator(pred).rows(te.X, te.fp, te.dose_time)])
        safe_te = (te.y_cls.astype(int) == SAFE_IDX)
        safe_tr = (tr.y_cls.astype(int) == SAFE_IDX)
        stratum = np.round(te.dose_time[:, 1].astype(float), 4)

        n_strata = len(np.unique(stratum))
        mixed = sum(1 for u in np.unique(stratum)
                    if len(set(safe_te[stratum == u].tolist())) > 1)
        bl = baselines(tr, te, safe_tr, safe_te)
        both = len(set(safe_te.tolist())) > 1
        s_auc, s_pairs = stratified_auc(S, safe_te, stratum)
        per_fold[d] = {
            "n": int(len(S)), "n_safe": int(safe_te.sum()),
            "n_strata": int(n_strata), "n_mixed_strata": int(mixed),
            "model_prauc": float(average_precision_score(safe_te, S)) if both else None,
            "model_roc": float(roc_auc_score(safe_te, S)) if both else None,
            "baselines": {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
                          for k, v in bl.items()},
            "stratified_auc_model": None if not s_pairs else float(s_auc),
            "stratified_pairs": int(s_pairs),
        }
        pool["S"].extend(S.tolist())
        pool["safe"].extend(safe_te.tolist())
        # strata must not merge across donors: the same hour in two donors is not one stratum
        pool["stratum"].extend([f"{d}:{v}" for v in stratum])
        pool["time_only"].extend(bl["time_only"].get("_score", [np.nan] * len(S)))

    S = np.asarray(pool["S"], float)
    safe = np.asarray(pool["safe"], bool)
    stratum = np.asarray(pool["stratum"])
    t_only = np.asarray(pool["time_only"], float)

    pooled = {
        "n": int(len(S)), "n_safe": int(safe.sum()),
        "model_prauc": float(average_precision_score(safe, S)),
        "model_roc": float(roc_auc_score(safe, S)),
        "time_only_prauc": (float(average_precision_score(safe, t_only))
                            if np.isfinite(t_only).all() else None),
        "time_only_roc": (float(roc_auc_score(safe, t_only))
                          if np.isfinite(t_only).all() else None),
        "n_strata": int(len(np.unique(stratum))),
        "n_mixed_strata": int(sum(1 for u in np.unique(stratum)
                                  if len(set(safe[stratum == u].tolist())) > 1)),
    }
    pooled["stratified_model"] = perm_null(S, safe, stratum)
    return {"suffix": SUFFIX, "folds": per_fold, "pooled": pooled, "errors": errs}


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                     # noqa: BLE001
            pass
    r = run()
    p = r["pooled"]
    print(f"\nSTAGE 18 -- is the fate head predicting biology, or reading a clock?  ({r['suffix']})")
    if r["errors"]:
        print(f"  folds unavailable: {r['errors']}")

    print("\n  MARGINAL metrics -- inflated by the timepoint the model was handed")
    print(f"     {'fold':<6}{'n':>4}{'strata':>8}{'mixed':>7}{'model PR':>10}{'time-only':>11}"
          f"{'expr-only':>11}{'expr+time':>11}")
    for d, f in r["folds"].items():
        b = f["baselines"]

        def g(k, _b=b):
            v = _b[k]["prauc"]
            return f"{v:.3f}" if v is not None else "n/a"
        mp = f"{f['model_prauc']:.3f}" if f["model_prauc"] is not None else "n/a"
        print(f"     {d:<6}{f['n']:>4}{f['n_strata']:>8}{f['n_mixed_strata']:>7}{mp:>10}"
              f"{g('time_only'):>11}{g('expr_only'):>11}{g('expr_plus_time'):>11}")
    print(f"\n     POOLED  model PR-AUC {p['model_prauc']:.3f}  ROC {p['model_roc']:.3f}")
    if p["time_only_prauc"] is not None:
        print(f"             time-only PR-AUC {p['time_only_prauc']:.3f}  "
              f"ROC {p['time_only_roc']:.3f}")

    print("\n  THE DECISIVE TEST -- WITHIN-TIMEPOINT concordance (dose_time is constant inside a")
    print("  stratum, so anything right here came from EXPRESSION)")
    print(f"     {p['n_mixed_strata']} of {p['n_strata']} timepoints carry more than one class")
    sm = p["stratified_model"]
    if sm["n_pairs"]:
        print(f"     stratified AUC {sm['observed']:.3f} over {sm['n_pairs']} (safe, unsafe) pairs")
        print(f"     permutation p = {sm['p']:.4f}  (scores shuffled WITHIN strata, "
              f"{N_PERM} draws)")
    else:
        print("     NO within-timepoint pairs exist: fate is a function of timepoint on this data,")
        print("     and the marginal numbers above cannot separate biology from the calendar.")

    print("\n     per fold:")
    for d, f in r["folds"].items():
        s = f["stratified_auc_model"]
        print(f"       {d:<5} pairs {f['stratified_pairs']:>3}   "
              f"AUC {'n/a' if s is None else f'{s:.3f}'}")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
