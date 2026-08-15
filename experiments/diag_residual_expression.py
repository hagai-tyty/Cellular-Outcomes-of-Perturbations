"""Does EARLY EXPRESSION predict the late residual, after donor age is removed?  (read-only)

    python experiments/diag_residual_expression.py

THE LAST LIVE OPTION
--------------------
`diag_early_late_forward` found the early->late signal is donor chronological age: the partial
correlation of early CD13 clock age with the late plateau, given donor age, is -0.064. But that
used a one-number SUMMARY of the early window (its clock age), which throws away 2,000 dimensions.
This asks whether the early EXPRESSION carries anything about the late outcome that donor age does
not already give.

If yes -- even weakly -- there is a real forward signal and the project has a target worth pursuing.
If no, the conclusion is that at 6 donors spanning ages 0-53 the outcome is donor age, and the
constraint is DATA, not architecture.

DESIGN, and why each piece is necessary at n=6
----------------------------------------------
* LEAVE-ONE-DONOR-OUT. The donor is the independent unit; there are 6.
* THE DONOR-AGE FIT IS REFIT INSIDE EACH FOLD. Residualising on all 6 first and then doing LOO
  would let the held-out donor help define its own residual -- a leak that manufactures signal.
* NO ANALYTIC p-VALUE. At n=6, with 2,000 features and ridge, no closed-form null is credible.
  Significance comes from a PERMUTATION NULL: the same full procedure, with the target shuffled
  across donors, repeated N_PERM times.
* A GRID OF ALPHAS, ALL REPORTED. Selecting alpha by nested CV on 6 points is theatre. Every
  alpha's result is printed, and the verdict requires the effect to survive across the grid rather
  than at one hand-picked value.
* A POSITIVE CONTROL. The same procedure against the RAW late age (not residualised) must show
  signal, because donor age predicts it at +0.931 and expression tracks donor age. If the control
  is null too, the procedure is broken and the residual result means nothing.

PRE-REGISTERED READING (constants below, fixed before running)
  SIGNAL   observed LOO Spearman exceeds the PERM_PCTILE-th percentile of its own permutation null
           for at least MIN_ALPHAS_PASSING of the alphas.
  NULL     otherwise.

  A NULL result here does NOT establish that no forward signal exists. With 6 donors, one covariate
  and 2,000 features, the power to detect anything short of a very strong effect is close to nil.
  It bounds what THIS dataset can answer -- which is the useful output either way.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))

from cellfate.data.aging import LinearClock  # noqa: E402
from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402
from cellfate.data.sources import GillReprogrammingSource  # noqa: E402

GILL_EXPR = r"D:\Gill\GSE165176_Log2_RPM_Sendai_reprogramming (1).txt.gz"
GILL_SERIES = r"D:\Gill\GSE165176_series_matrix.txt.gz"
CLOCK_PATH = ROOT / "configs" / "clocks" / "fleischer_clock.json"
PANEL_PATH = ROOT / "cellfate_loocv_N2_c7t" / "panel.json"

EARLY_LO, EARLY_HI, LATE_LO = 7.0, 29.0, 34.0
ALPHAS = (1.0, 10.0, 100.0, 1_000.0, 10_000.0)
N_PERM = 2000
PERM_PCTILE = 95.0
MIN_ALPHAS_PASSING = 3
SEED = 0


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def ridge_fit_predict(xtr: np.ndarray, ytr: np.ndarray, xte: np.ndarray,
                      alpha: float) -> np.ndarray:
    """Ridge in closed form on TRAIN-only standardisation.

    DUAL FORM. With 5 training donors and 2,000 features the primal solve is 2000x2000; the
    identity `X^T(XX^T + aI)^-1 = (X^TX + aI)^-1 X^T` gives the same answer from a 5x5 solve. At
    2,000 permutations x 5 alphas x 2 runs the primal version needs ~120,000 large solves and does
    not finish; this is the same result in seconds.

    Centring y on the TRAIN mean matters: a degenerate fit then predicts the train mean rather
    than zero, so a useless model scores like the mean baseline instead of accidentally scoring
    like a correct one on a residualised (near-zero-mean) target.
    """
    mu, sd = xtr.mean(0), xtr.std(0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    a, b = (xtr - mu) / sd, (xte - mu) / sd
    ybar = ytr.mean()
    n = a.shape[0]
    dual = np.linalg.solve(a @ a.T + alpha * np.eye(n), ytr - ybar)
    return b @ (a.T @ dual) + ybar


def loo_spearman(X: np.ndarray, y: np.ndarray, age: np.ndarray, alpha: float,
                 residualise: bool) -> float:
    """Leave-one-donor-out. The donor-age fit is refit on TRAIN donors only inside every fold."""
    n = len(y)
    pred, actual = np.zeros(n), np.zeros(n)
    for i in range(n):
        tr = np.arange(n) != i
        ytr, yte = y[tr].copy(), y[i]
        if residualise:
            # linear fit late ~ donor_age on TRAIN donors only
            A = np.vstack([age[tr], np.ones(tr.sum())]).T
            coef, *_ = np.linalg.lstsq(A, ytr, rcond=None)
            ytr = ytr - (coef[0] * age[tr] + coef[1])
            yte = yte - (coef[0] * age[i] + coef[1])
        pred[i] = float(ridge_fit_predict(X[tr], ytr, X[i:i + 1], alpha)[0])
        actual[i] = yte
    # A residual that is numerically zero relative to the quantity it came from is NOT a target.
    # `spearman`'s std==0 guard tests EXACT zero, so residuals of ~1e-14 would sail through and
    # produce a correlation over floating-point dust -- which is precisely the defect that made
    # "Ranking generalizes: Spearman 0.40" out of RES values maxing at 1.6e-4. Caught by
    # `test_the_donor_age_fit_is_refit_per_fold_not_once_globally`, which constructs exactly that.
    scale = float(np.std(y)) or 1.0
    if float(np.std(actual)) < 1e-9 * scale:
        return float("nan")
    return spearman(actual, pred)


def permutation_null(X, y, age, alpha, residualise, n_perm=N_PERM, seed=SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(n_perm)
    for k in range(n_perm):
        p = rng.permutation(len(y))
        out[k] = loo_spearman(X, y[p], age, alpha, residualise)
    return out[np.isfinite(out)]


def load_donor_table() -> tuple[pd.DataFrame, np.ndarray]:
    src = GillReprogrammingSource(GILL_EXPR, GILL_SERIES)
    src.bulk_integrity_gate = True
    src._load()
    clock = LinearClock.from_json(str(CLOCK_PATH))
    rpm = src._rpm
    # normalize_counts applies CP10k AND log1p itself (normalize.py:29) -- not double-wrapped.
    expr = normalize_counts(np.clip(rpm.to_numpy(dtype=np.float64).T, 0.0, None), target_sum=1e4)
    ages = clock.predict_age(expr, list(rpm.index))

    panel = [g for g in GenePanel.load(str(PANEL_PATH)).genes if g in set(rpm.index)]
    idx = [list(rpm.index).index(g) for g in panel]

    meta = []
    for j, c in enumerate(rpm.columns):
        m = src._meta[c]
        meta.append({"j": j, "donor": m["donor"], "day": float(m["day"]),
                     "donor_age": float(m["age"]), "age": float(ages[j])})
    M = pd.DataFrame(meta)

    donors, rows, y, da = [], [], [], []
    for d, g in M.groupby("donor"):
        e = g[(g.day >= EARLY_LO) & (g.day <= EARLY_HI)]
        late = g[g.day >= LATE_LO]
        if e.empty or late.empty:
            continue
        donors.append(d)
        rows.append(expr[e.j.to_numpy(), :][:, idx].mean(0))
        y.append(late.age.mean())
        da.append(g.donor_age.iloc[0])
    return (pd.DataFrame({"donor": donors, "late": y, "donor_age": da}),
            np.vstack(rows))


def main() -> None:
    T, X = load_donor_table()
    y, age = T.late.to_numpy(float), T.donor_age.to_numpy(float)
    print("=" * 96)
    print("EARLY EXPRESSION -> LATE RESIDUAL (donor age removed)   leave-one-donor-out")
    print(f"  donors: {list(T.donor)}   features: {X.shape[1]} panel genes")
    print(f"  pre-registered: SIGNAL iff observed > {PERM_PCTILE:.0f}th pct of its own permutation "
          f"null\n                  for >= {MIN_ALPHAS_PASSING} of {len(ALPHAS)} alphas "
          f"({N_PERM} permutations each)")
    print("=" * 96)
    res: dict = {"donors": list(T.donor), "n_features": int(X.shape[1]), "runs": {}}

    for label, residualise in (("POSITIVE CONTROL: raw late age", False),
                               ("TEST: late residual | donor age", True)):
        print(f"\n[{label}]")
        print(f"  {'alpha':>10}{'LOO spearman':>15}{'null p95':>11}{'null max':>11}"
              f"{'pctile':>9}   pass")
        n_pass, rows = 0, {}
        for a in ALPHAS:
            obs = loo_spearman(X, y, age, a, residualise)
            null = permutation_null(X, y, age, a, residualise)
            p95 = float(np.percentile(null, PERM_PCTILE))
            pct = float((null < obs).mean() * 100) if np.isfinite(obs) else float("nan")
            ok = bool(np.isfinite(obs) and obs > p95)
            n_pass += ok
            rows[str(a)] = {"observed": obs, "null_p95": p95, "null_max": float(null.max()),
                            "percentile": pct, "pass": ok}
            print(f"  {a:>10.0f}{obs:>15.3f}{p95:>11.3f}{float(null.max()):>11.3f}"
                  f"{pct:>9.1f}   {'YES' if ok else 'no'}")
        verdict = "SIGNAL" if n_pass >= MIN_ALPHAS_PASSING else "NULL"
        print(f"  -> {n_pass} of {len(ALPHAS)} alphas pass  ->  {verdict}")
        res["runs"][label] = {"alphas": rows, "n_pass": n_pass, "verdict": verdict}

    _RESULTS.mkdir(exist_ok=True)
    out = _RESULTS / "diag_residual_expression_results.json"
    out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
