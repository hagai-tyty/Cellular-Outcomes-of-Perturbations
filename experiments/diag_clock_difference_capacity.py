"""PHASE 0 -- which clock variant is best at DIFFERENCES?  Ground truth, n=133.  (read-only)

    python experiments/diag_clock_difference_capacity.py

Pre-registered in `plans/THREE_TESTS_PREREG.md` Phase 0, before the run.

THE MECHANISM UNDER TEST
------------------------
`top100` is not a different instrument -- it is the 100 largest-|weight| genes of the SAME
Fleischer clock. A dense ridge over 33,155 genes fit from 133 samples carries many small weights
fit to noise. In ABSOLUTE age those largely cancel: the fit was optimised for that target on that
cohort. In a DIFFERENCE between two samples they do NOT cancel -- the two noise terms are
independent and compound. ΔAge is a difference. Hence:

    PREDICTION: dense wins on absolute age; sparse wins on differences.

If that holds, "use top100 everywhere" is the wrong rule. The right one is *dense for absolute,
sparse for differences*, and each downstream phase should use whichever it actually needs.

WHY THIS COHORT
---------------
GSE113957 has 133 donors with declared ages, hence **8,778 pairs with a KNOWN age difference** --
ground truth for exactly the quantity ΔAge is, at a scale nothing else here offers, and with no
methylation clock in the loop.

It is also the cohort the clock was FITTED on. That makes absolute performance optimistic for
`raw` **by construction**, which makes a `raw` loss on differences stronger evidence, not weaker.

TWO QUESTIONS, REPORTED SEPARATELY, BECAUSE THEY ARE NOT THE SAME
----------------------------------------------------------------
AS-USED    the variant's own weights, no refitting. This is how the ledger and the instrument-floor
           comparison use them. Differences need no intercept, so this is well defined for every
           variant; absolute age is NOT, because truncating weights invalidates the clock's
           intercept. Reported for differences only.
CALIBRATED k-fold CV fit of `age ~ a*score + b` per variant. This separates INFORMATION CONTENT
           from calibration -- a variant could carry the signal and merely be on the wrong scale.

PRE-REGISTERED READING
  MECHANISM CONFIRMED  the ranking by difference-error differs from the ranking by absolute error,
                       AND a sparse variant beats `raw` on difference-error.
  MECHANISM REFUTED    the two rankings agree.

Pair differences are NOT independent (each donor appears in 132 pairs), so pair statistics are
descriptive and no p-value is computed from them.

PART A IS INVALID FOR THIS PURPOSE, AND PART B IS WHY THIS FILE HAS TWO PARTS
----------------------------------------------------------------------------
Part A ran first and returned `raw` at **MAE_abs 0.13 yr, r = 1.000** on 133 donors. The clock's own
published cross-validated MAE is **12.27**. A 94x gap is not a good clock -- it is MEMORISATION:
33,155 weights fitted on these exact 133 samples, then scored on them.

The pre-registration flagged this cohort as "optimistic for raw by construction". That was right in
direction and wrong by two orders of magnitude. In-sample, `raw` cannot lose, and truncating a
memorised fit destroys the memorisation rather than removing noise -- so Part A's ranking says
nothing about the mechanism and its verdict is NOT reported as a finding.

Part A is kept, not deleted, because the memorisation itself is worth recording: it quantifies how
overfit the shipped dense clock is, which is the same property the mechanism hypothesis blames for
its poor DIFFERENCE behaviour out of sample.

PART B -- the honest test. Refit a dense clock INSIDE each CV fold, truncate THAT fit, and score
out of fold. Nothing is in-sample, so truncation is tested on its merits. Truncated weights get
their scale and offset refitted on the TRAINING fold only, so absolute error is comparable across
truncation levels without leaking the held-out ages.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dac = _load("dac", "experiments/diag_age_capacity.py")

CLOCK_PATH = ROOT / "configs" / "clocks" / "fleischer_clock.json"
N_SPLITS = 10
SEED = 0
TOP_NS = (100, 500, 2000)


def variant_weights(weights: dict[str, float], genes: list[str]) -> dict[str, np.ndarray]:
    """One weight vector per variant, aligned to `genes`. Absent genes get 0.0 -- the same rule
    `LinearClock.predict_age` uses (aging.py:55)."""
    w_full = np.array([float(weights.get(g, 0.0)) for g in genes], dtype=np.float64)
    out = {"raw": w_full}
    order = np.argsort(-np.abs(w_full))
    for n in TOP_NS:
        w = np.zeros_like(w_full)
        keep = order[:n]
        w[keep] = w_full[keep]
        out[f"top{n}"] = w
    # covnorm: rescale by the share of total |weight| mass this gene space actually covers, so the
    # uncovered fraction stops silently reading as zero.
    total = sum(abs(v) for v in weights.values())
    covered = float(np.abs(w_full).sum())
    out["covnorm"] = w_full * (total / covered) if covered > 0 else w_full
    return out


def rank_transform(expr: np.ndarray) -> np.ndarray:
    """Rank each SAMPLE across genes, scaled to [0,1]. Immune to library size and to the
    log/linear choice -- and, per the earlier shrinkage control, prone to collapsing."""
    out = np.empty_like(expr)
    for i in range(expr.shape[0]):
        out[i] = pd.Series(expr[i]).rank().to_numpy() / expr.shape[1]
    return out


def cv_linear(score: np.ndarray, age: np.ndarray, n_splits: int = N_SPLITS,
              seed: int = SEED) -> tuple[np.ndarray, float]:
    """Held-out predictions from `age ~ a*score + b`, plus the mean fitted slope. One free scale
    and offset per variant, so INFORMATION content is compared rather than calibration."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(age))
    pred = np.empty(len(age))
    slopes = []
    for te in np.array_split(order, n_splits):
        tr = np.setdiff1d(order, te)
        a, b = np.polyfit(score[tr], age[tr], 1)
        pred[te] = a * score[te] + b
        slopes.append(a)
    return pred, float(np.mean(slopes))


def pair_stats(pred_diff: np.ndarray, true_diff: np.ndarray) -> dict:
    return {"mae_diff": float(np.abs(pred_diff - true_diff).mean()),
            "r_diff": float(np.corrcoef(pred_diff, true_diff)[0, 1])}


def fit_dense_ridge(xtr: np.ndarray, ytr: np.ndarray, alpha: float) -> tuple[np.ndarray, float]:
    """Dense clock refit from scratch, dual form (n x n, not 39k x 39k). Returns (weights, ybar)."""
    ybar = float(ytr.mean())
    n = xtr.shape[0]
    dual = np.linalg.solve(xtr @ xtr.T + alpha * np.eye(n), ytr - ybar)
    return xtr.T @ dual, ybar


def refit_cv(expr: np.ndarray, age: np.ndarray, alpha: float, top_ns: tuple[int, ...],
             n_splits: int = N_SPLITS, seed: int = SEED) -> dict:
    """PART B. Refit inside each fold, truncate THAT fit, score out of fold.

    The truncated weights get a scale+offset refitted on the TRAINING fold only -- truncation
    changes the score's scale, and without that rescaling the comparison would measure calibration
    rather than information. Held-out ages are never touched.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(age))
    levels = ("full", *[f"top{n}" for n in top_ns])
    pred = {k: np.empty(len(age)) for k in levels}
    for te in np.array_split(order, n_splits):
        tr = np.setdiff1d(order, te)
        w, ybar = fit_dense_ridge(expr[tr], age[tr], alpha)
        rank = np.argsort(-np.abs(w))
        for lvl in levels:
            if lvl == "full":
                wl = w
            else:
                wl = np.zeros_like(w)
                k = rank[: int(lvl[3:])]
                wl[k] = w[k]
            s_tr, s_te = expr[tr] @ wl, expr[te] @ wl
            if np.std(s_tr) < 1e-12:
                pred[lvl][te] = ybar
                continue
            a, b = np.polyfit(s_tr, age[tr], 1)      # train-fold calibration only
            pred[lvl][te] = a * s_te + b
    iu = np.triu_indices(len(age), k=1)
    true_diff = (age[:, None] - age[None, :])[iu]
    out = {}
    for lvl, p in pred.items():
        d = (p[:, None] - p[None, :])[iu]
        out[lvl] = {"mae_abs": float(np.abs(p - age).mean()),
                    "r_abs": float(np.corrcoef(p, age)[0, 1]),
                    "mae_diff": float(np.abs(d - true_diff).mean()),
                    "r_diff": float(np.corrcoef(d, true_diff)[0, 1]),
                    "sd_ratio": float(np.std(d, ddof=1) / np.std(true_diff, ddof=1))}
    return out


def main() -> None:
    meta = dac.load_meta()
    expr, genes = dac.load_expression(meta.gsm.tolist())
    keep = meta.disease.eq("Normal").to_numpy() & meta.age.notna().to_numpy()
    expr, age = expr[keep], meta.age.to_numpy(float)[keep]
    n = len(age)

    clock = json.loads(CLOCK_PATH.read_text(encoding="utf-8"))
    wv = variant_weights(clock["weights"], genes)
    ranked = rank_transform(expr)

    iu = np.triu_indices(n, k=1)
    true_diff = (age[:, None] - age[None, :])[iu]

    print("=" * 108)
    print(f"PHASE 0 -- clock variants on ABSOLUTE age vs DIFFERENCES.  n={n} donors, "
          f"{len(true_diff)} pairs, ground-truth ages")
    print("  prediction under test: dense wins ABSOLUTE, sparse wins DIFFERENCES")
    print("=" * 108)
    print(f"  {'variant':<10}{'n_genes':>9}{'AS-USED':>26}{'CALIBRATED (CV)':>34}")
    print(f"  {'':<10}{'':>9}{'MAE_diff':>12}{'r_diff':>8}{'sd_ratio':>6}"
          f"{'MAE_abs':>10}{'r_abs':>8}{'MAE_diff':>10}{'ratio':>8}")
    res: dict = {"n_donors": int(n), "n_pairs": int(len(true_diff)), "variants": {}}

    for name in ["raw", "covnorm", "top2000", "top500", "top100", "ranknorm"]:
        if name == "ranknorm":
            score = ranked @ wv["raw"]
        else:
            score = expr @ wv[name]
        nz = int((wv["raw"] != 0).sum()) if name in ("raw", "covnorm", "ranknorm") \
            else int((wv[name] != 0).sum())

        # AS-USED: the variant's own weights, no refit. Differences need no intercept.
        d_used = (score[:, None] - score[None, :])[iu]
        used = pair_stats(d_used, true_diff)
        sd_ratio = float(np.std(d_used, ddof=1) / np.std(true_diff, ddof=1))

        # CALIBRATED: one free scale+offset, cross-validated.
        pred, slope = cv_linear(score, age)
        mae_abs = float(np.abs(pred - age).mean())
        r_abs = float(np.corrcoef(pred, age)[0, 1])
        d_cal = (pred[:, None] - pred[None, :])[iu]
        cal = pair_stats(d_cal, true_diff)
        # <1 means errors CANCEL in a difference, >1 means they COMPOUND
        ratio = cal["mae_diff"] / (np.sqrt(2.0) * mae_abs)

        print(f"  {name:<10}{nz:>9}{used['mae_diff']:>12.2f}{used['r_diff']:>8.3f}"
              f"{sd_ratio:>6.2f}{mae_abs:>10.2f}{r_abs:>8.3f}{cal['mae_diff']:>10.2f}"
              f"{ratio:>8.2f}")
        res["variants"][name] = {"n_genes": nz, "as_used": {**used, "sd_ratio": sd_ratio},
                                 "calibrated": {"mae_abs": mae_abs, "r_abs": r_abs,
                                                "slope": slope, **cal, "compound_ratio": ratio}}

    order_abs = sorted(res["variants"], key=lambda k: res["variants"][k]["calibrated"]["mae_abs"])
    order_dif = sorted(res["variants"], key=lambda k: res["variants"][k]["as_used"]["mae_diff"])
    print(f"\n  ranking by CALIBRATED absolute error : {order_abs}")
    print(f"  ranking by AS-USED difference error  : {order_dif}")
    sparse_beats_raw = any(
        res["variants"][v]["as_used"]["mae_diff"] < res["variants"]["raw"]["as_used"]["mae_diff"]
        for v in ("top100", "top500", "top2000"))
    verdict = ("MECHANISM CONFIRMED" if (order_abs != order_dif and sparse_beats_raw)
               else "MECHANISM REFUTED")
    print(f"\n  rankings differ: {order_abs != order_dif}   a sparse variant beats raw on "
          f"differences: {sparse_beats_raw}")
    print(f"  -> {verdict}")
    res["ranking_absolute"], res["ranking_difference"] = order_abs, order_dif
    res["partA_verdict_INVALID"] = verdict
    res["partA_note"] = ("raw scores in-sample here: the clock was fitted on these 133 donors, so "
                         "its 0.13 yr MAE against a published cv_mae of 12.27 is memorisation. "
                         "Part A cannot test the mechanism; Part B does.")
    print("\n  !! PART A IS NOT A FINDING. `raw` MAE_abs "
          f"{res['variants']['raw']['calibrated']['mae_abs']:.2f} yr against a published cv_mae of "
          "12.27 is\n     MEMORISATION -- the clock was fitted on these exact donors. Truncating a "
          "memorised fit\n     destroys the memorisation, not noise. Part B refits inside folds.")

    # ---- PART B: refit inside folds, truncate the refit, score out of fold ------------------ #
    print("\n" + "=" * 108)
    print("PART B -- HONEST TEST: dense clock REFIT inside each CV fold, truncations scored "
          "out of fold")
    print("=" * 108)
    res["partB"] = {}
    for alpha in (1.0, 10.0, 100.0, 1000.0):
        r = refit_cv(expr, age, alpha, TOP_NS)
        res["partB"][str(alpha)] = r
        print(f"\n  alpha={alpha:.0f}")
        print(f"    {'level':<9}{'MAE_abs':>10}{'r_abs':>8}{'MAE_diff':>11}{'r_diff':>9}"
              f"{'sd_ratio':>10}")
        for lvl, v in r.items():
            print(f"    {lvl:<9}{v['mae_abs']:>10.2f}{v['r_abs']:>8.3f}{v['mae_diff']:>11.2f}"
                  f"{v['r_diff']:>9.3f}{v['sd_ratio']:>10.2f}")
        best_abs = min(r, key=lambda k: r[k]["mae_abs"])
        best_dif = min(r, key=lambda k: r[k]["mae_diff"])
        print(f"    best absolute: {best_abs}   best difference: {best_dif}")

    # Verdict on the mechanism, from Part B only.
    wins = {a: (min(v, key=lambda k: v[k]["mae_abs"]), min(v, key=lambda k: v[k]["mae_diff"]))
            for a, v in res["partB"].items()}
    sparse_wins_diff = sum(1 for _, (_, d) in wins.items() if d != "full")
    differs = sum(1 for _, (a_, d) in wins.items() if a_ != d)
    mech = ("MECHANISM SUPPORTED" if sparse_wins_diff > len(wins) / 2
            else "MECHANISM NOT SUPPORTED")
    print(f"\n  sparse wins DIFFERENCES in {sparse_wins_diff}/{len(wins)} alphas; "
          f"absolute and difference winners differ in {differs}/{len(wins)}")
    print(f"  -> {mech}")
    res["partB_verdict"] = mech

    # ---- PART C: the SHIPPED clock, full vs truncated, OUT of cohort ------------------------ #
    # Part B refit in-cohort, so truncation had no cohort-specific overfitting to remove and could
    # only discard signal. The instrument-floor win for top100 used the SHIPPED clock on a
    # DIFFERENT cohort, where its 33,155 weights carry GSE113957-specific structure. That predicts
    # truncation helps out of cohort and not in it. This is the direct test, and it is what decides
    # which clock the later phases use.
    print("\n" + "=" * 108)
    print("PART C -- SHIPPED clock on OUT-OF-COHORT day-0 fibroblasts with known ages")
    print("=" * 108)
    dat = _load("dat", "experiments/diag_age_transfer.py")
    res["partC"] = {}
    for cname, (ex, se) in (("GSE165176 Gill", dat.GILL), ("GSE165177", dat.TRANS)):
        d = dat.load_bulk_day0(ex, se, 12)
        if d.empty:
            continue
        ages = d["__age__"].to_numpy(float)
        te = d.drop(columns="__age__")
        cols = [g for g in te.columns if g in set(genes)]
        gi = {g: i for i, g in enumerate(genes)}
        X = te[cols].to_numpy(dtype=np.float64)
        print(f"\n  [{cname}]  n={len(ages)}  ages {sorted(ages.astype(int))}  genes {len(cols)}")
        print(f"    {'level':<9}{'rho(score,age)':>16}{'sd(score)':>12}")
        res["partC"][cname] = {"n": len(ages), "ages": ages.tolist(), "levels": {}}
        for lvl in ["raw", "top2000", "top500", "top100"]:
            w = np.array([wv[lvl][gi[g]] for g in cols], dtype=np.float64)
            s = X @ w
            rho = (float(pd.Series(s).corr(pd.Series(ages), method="spearman"))
                   if np.std(s) > 0 and np.std(ages) > 0 else float("nan"))
            print(f"    {lvl:<9}{rho:>16.3f}{np.std(s, ddof=1):>12.2f}")
            res["partC"][cname]["levels"][lvl] = {"spearman": rho,
                                                  "sd_score": float(np.std(s, ddof=1))}

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_clock_difference_capacity_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
