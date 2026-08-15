"""Is the early-expression -> late-residual result FRAGILE?  (read-only, post-hoc)

    python experiments/diag_residual_robustness.py

WHAT THIS CAN AND CANNOT DO
---------------------------
`diag_residual_expression` found SIGNAL: LOO Spearman 0.943, 5 of 5 alphas above their permutation
nulls. This is the SAME six donors, so nothing here can CONFIRM that -- confirmation needs donors
not used to find it. What it can do is show whether the effect survives reasonable changes to
choices that were made arbitrarily. An effect that evaporates when the early window moves by a few
days, or when a different gene set is used, is noise dressed as a result.

The procedure is IMPORTED from `diag_residual_expression`, not reimplemented, so every variant runs
the identical leave-one-donor-out + per-fold age refit + permutation null.

ONE AXIS AT A TIME
------------------
Baseline is the published configuration: early window 7-29 d, all markers, 1,903 panel genes.
Each variant changes exactly ONE thing. A full grid would be 48 cells of multiple-testing soup and
would let a reader pick whichever supports their prior.

  early window   (7,15)  (7,21)  (11,29)
  marker         CD13-only   SSEA4-only
  feature set    clock genes   top-500 variable   random-500

PRE-REGISTERED READING (constants below, fixed before running)
  ROBUST    the baseline shows SIGNAL and at least MIN_ROBUST of the RUNS (baseline + variants)
            show SIGNAL.
  FRAGILE   otherwise -- which would mean the headline result should not be relied on.

Some variants are EXPECTED to fail legitimately: SSEA4-only was the weakest predictor in the
correlation analysis (+0.42 vs +0.83 for CD13), and a 7-15 d window has a third of the samples.
That is why the bar is a majority rather than unanimity, and why every variant is printed.

NOTE ON `top-500 variable`: it selects features using X from all six donors, including the held-out
one. That is unsupervised (it never touches y), but it is not blind, and it is reported as such.
`random-500` is not a control on validity -- the permutation null already covers that -- it asks
whether the effect is specific to a gene set or just to donors resembling each other.
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

_spec = importlib.util.spec_from_file_location(
    "dre", ROOT / "experiments" / "diag_residual_expression.py")
dre = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dre)

from cellfate.data.aging import LinearClock  # noqa: E402
from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402
from cellfate.data.sources import GillReprogrammingSource  # noqa: E402

MIN_ROBUST = 6          # of 9 runs (baseline + 8 variants)
RANDOM_GENES_SEED = 12345
N_RANDOM = 500
N_HVG = 500


def build(markers: tuple[str, ...], early: tuple[float, float], feature_set: str):
    """Per-donor early-mean expression over `feature_set`, plus late age and donor age."""
    src = GillReprogrammingSource(dre.GILL_EXPR, dre.GILL_SERIES)
    src.bulk_integrity_gate = True
    src._load()
    clock = LinearClock.from_json(str(dre.CLOCK_PATH))
    rpm = src._rpm
    expr = normalize_counts(np.clip(rpm.to_numpy(dtype=np.float64).T, 0.0, None), target_sum=1e4)
    ages = clock.predict_age(expr, list(rpm.index))
    idx_of = {g: i for i, g in enumerate(rpm.index)}

    if feature_set == "panel":
        genes = [g for g in GenePanel.load(str(dre.PANEL_PATH)).genes if g in idx_of]
    elif feature_set == "clock":
        genes = [g for g in clock.weights if g in idx_of]
    elif feature_set == "random":
        rng = np.random.default_rng(RANDOM_GENES_SEED)
        all_g = list(rpm.index)
        genes = [all_g[i] for i in rng.choice(len(all_g), size=N_RANDOM, replace=False)]
    elif feature_set == "hvg":
        genes = None            # resolved after the donor means are formed
    else:
        raise ValueError(feature_set)

    rows = []
    for j, c in enumerate(rpm.columns):
        m = src._meta[c]
        mk = "CD13" if "_CD13_" in c else ("SSEA4" if "_SSEA4_" in c else "Fib")
        rows.append({"j": j, "donor": m["donor"], "day": float(m["day"]), "marker": mk,
                     "donor_age": float(m["age"]), "age": float(ages[j])})
    M = pd.DataFrame(rows)

    donors, X, y, da = [], [], [], []
    for d, g in M.groupby("donor"):
        e = g[(g.day >= early[0]) & (g.day <= early[1]) & (g.marker.isin(markers))]
        late = g[g.day >= dre.LATE_LO]
        if e.empty or late.empty:
            continue
        donors.append(d)
        X.append(expr[e.j.to_numpy(), :].mean(0))
        y.append(late.age.mean())
        da.append(g.donor_age.iloc[0])
    Xf = np.vstack(X)
    if feature_set == "hvg":
        cols = np.argsort(-Xf.std(0))[:N_HVG]
    else:
        cols = np.array([idx_of[g] for g in genes])
    return donors, Xf[:, cols], np.asarray(y, float), np.asarray(da, float)


def run_one(name: str, markers, early, feature_set: str) -> dict:
    donors, X, y, age = build(markers, early, feature_set)
    if len(donors) < 5:
        return {"name": name, "n_donors": len(donors), "verdict": "SKIPPED (too few donors)"}
    n_pass, alphas = 0, {}
    for a in dre.ALPHAS:
        obs = dre.loo_spearman(X, y, age, a, residualise=True)
        null = dre.permutation_null(X, y, age, a, residualise=True)
        p95 = float(np.percentile(null, dre.PERM_PCTILE))
        ok = bool(np.isfinite(obs) and obs > p95)
        n_pass += ok
        alphas[str(a)] = {"observed": obs, "null_p95": p95, "pass": ok}
    verdict = "SIGNAL" if n_pass >= dre.MIN_ALPHAS_PASSING else "null"
    med = float(np.nanmedian([v["observed"] for v in alphas.values()]))
    return {"name": name, "n_donors": len(donors), "n_features": int(X.shape[1]),
            "median_observed": med, "n_pass": n_pass, "verdict": verdict, "alphas": alphas}


RUNS = [
    ("BASELINE  7-29d, all markers, panel", ("CD13", "SSEA4", "Fib"), (7.0, 29.0), "panel"),
    ("window    7-15d",                     ("CD13", "SSEA4", "Fib"), (7.0, 15.0), "panel"),
    ("window    7-21d",                     ("CD13", "SSEA4", "Fib"), (7.0, 21.0), "panel"),
    ("window    11-29d",                    ("CD13", "SSEA4", "Fib"), (11.0, 29.0), "panel"),
    ("marker    CD13 only",                 ("CD13",),                (7.0, 29.0), "panel"),
    ("marker    SSEA4 only",                ("SSEA4",),               (7.0, 29.0), "panel"),
    ("features  clock genes",               ("CD13", "SSEA4", "Fib"), (7.0, 29.0), "clock"),
    ("features  top-500 variable",          ("CD13", "SSEA4", "Fib"), (7.0, 29.0), "hvg"),
    ("features  random-500",                ("CD13", "SSEA4", "Fib"), (7.0, 29.0), "random"),
]


def main() -> None:
    print("=" * 100)
    print("ROBUSTNESS OF THE EARLY-EXPRESSION -> LATE-RESIDUAL RESULT   (post-hoc, same 6 donors)")
    print(f"  pre-registered: ROBUST iff the baseline shows SIGNAL and >= {MIN_ROBUST} of "
          f"{len(RUNS)} runs do")
    print("  cannot confirm the result -- only test whether it survives arbitrary choices")
    print("=" * 100)
    print(f"  {'run':<38}{'feats':>7}{'median rho':>12}{'alphas':>8}   verdict")
    out, n_signal, baseline_ok = [], 0, False
    for name, mk, win, fs in RUNS:
        r = run_one(name, mk, win, fs)
        out.append(r)
        if r["verdict"] == "SIGNAL":
            n_signal += 1
            if name.startswith("BASELINE"):
                baseline_ok = True
        print(f"  {name:<38}{r.get('n_features', 0):>7}{r.get('median_observed', float('nan')):>12.3f}"
              f"{r.get('n_pass', 0):>4}/{len(dre.ALPHAS)}   {r['verdict']}")
    verdict = "ROBUST" if (baseline_ok and n_signal >= MIN_ROBUST) else "FRAGILE"
    print(f"\n  {n_signal} of {len(RUNS)} runs show SIGNAL; baseline "
          f"{'holds' if baseline_ok else 'FAILS'}  ->  {verdict}")
    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_residual_robustness_results.json"
    p.write_text(json.dumps({"runs": out, "n_signal": n_signal, "verdict": verdict},
                            indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
