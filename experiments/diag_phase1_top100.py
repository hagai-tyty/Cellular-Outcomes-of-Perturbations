"""PHASE 1 -- redo the two ΔAge forward tests on the clock Phase 0 licensed.  (read-only)

    python experiments/diag_phase1_top100.py

Pre-registered in `plans/THREE_TESTS_PREREG.md` Phase 1.

WHICH CLOCK, AND WHY -- decided by Phase 0, not by preference
-------------------------------------------------------------
Phase 0 tested and REFUTED the reason I originally gave for switching to top100 (that sparsity
helps DIFFERENCES). Truncation is worse than dense on in-cohort differences (Part B: top100 MAE_diff
27.73 vs full 17.80) and worse on out-of-cohort absolute age (Part C: rho 0.616 vs 0.872).

What top100 DOES fix is contamination by the perturbation. On 44 reprogramming conditions scored
against methylation: raw 22.69 (SD ratio 1.66), resid_cc 23.65 (cell cycle is not the culprit),
resid_pluri 13.00 (SD 1.14), top100 7.15 (SD 0.98). The dense clock reads the reprogramming
programme -- specifically pluripotency -- on top of age, and over-reports magnitude by 66%.

Both windows here are REPROGRAMMING samples, so top100 is the licensed readout. Dense would be
licensed for resting day-0 samples, and is not used.

WHAT IS RE-RUN, UNCHANGED IN EVERY OTHER RESPECT
------------------------------------------------
  A  early -> late, and the partial correlation given donor age   (was: partial -0.064)
  B  late residual from early EXPRESSION, LOO + permutation null  (was: SIGNAL, then FRAGILE 3/9)

Scale note: top100's truncated weights have no valid intercept, so its scores are not years. Every
statistic here (Pearson, Spearman, partial correlation, LOO Spearman) is invariant to a linear
rescaling, so this is sound -- but no MAE in years is reported for top100, because that number
would be meaningless.

PRE-REGISTERED READING (from the plan, fixed before this ran)
  CHANGED    a verdict flips: the partial clears |r| > 0.878 (df=3), or robustness reaches 6 of 9.
  UNCHANGED  otherwise -- and the earlier conclusions then stand on a BETTER instrument, which
             strengthens them rather than merely repeating them.
  The partial's magnitude is reported either way. A move from -0.064 to +0.3 that still misses the
  bar is informative and must not be reported as "unchanged".
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


delf = _load("delf", "experiments/diag_early_late_forward.py")
dre = _load("dre", "experiments/diag_residual_expression.py")
drr = _load("drr", "experiments/diag_residual_robustness.py")

from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402
from cellfate.data.sources import GillReprogrammingSource  # noqa: E402

TOP_N = 100
EARLY_LO, EARLY_HI, LATE_LO = delf.EARLY_LO, delf.EARLY_HI, delf.LATE_LO


def top_n_weights(weights: dict, genes: list[str], n: int) -> np.ndarray:
    w = np.array([float(weights.get(g, 0.0)) for g in genes], dtype=np.float64)
    out = np.zeros_like(w)
    keep = np.argsort(-np.abs(w))[:n]
    out[keep] = w[keep]
    return out


def build() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    src = GillReprogrammingSource(delf.GILL_EXPR, delf.GILL_SERIES)
    src.bulk_integrity_gate = True
    src._load()
    clock = json.loads(Path(delf.CLOCK_PATH).read_text(encoding="utf-8"))
    rpm = src._rpm
    genes = list(rpm.index)
    expr = normalize_counts(np.clip(rpm.to_numpy(dtype=np.float64).T, 0.0, None), target_sum=1e4)
    score = expr @ top_n_weights(clock["weights"], genes, TOP_N)
    rows = []
    for j, c in enumerate(rpm.columns):
        m = src._meta[c]
        mk = "CD13" if "_CD13_" in c else ("SSEA4" if "_SSEA4_" in c else "Fib")
        rows.append({"donor": m["donor"], "day": float(m["day"]), "marker": mk,
                     "donor_age": float(m["age"]), "score": float(score[j]), "j": j})
    return pd.DataFrame(rows), expr, genes


def main() -> None:
    df, expr, genes = build()
    print("=" * 100)
    print(f"PHASE 1 -- both ΔAge forward tests redone on top{TOP_N} (Phase 0 licensed it for "
          "PERTURBED samples)")
    print("=" * 100)
    res: dict = {"clock": f"top{TOP_N}", "n_gill_samples": int(len(df))}

    # ---- A: early -> late, and the partial given donor age -------------------------------- #
    per = {}
    for d, g in df.groupby("donor"):
        e = g[(g.day >= EARLY_LO) & (g.day <= EARLY_HI)]
        late = g[g.day >= LATE_LO]
        if e.empty or late.empty:
            continue
        per[d] = {"donor_age": g.donor_age.iloc[0], "late": late.score.mean(),
                  "early_all": e.score.mean(),
                  "early_cd13": e[e.marker == "CD13"].score.mean(),
                  "early_ssea4": e[e.marker == "SSEA4"].score.mean()}
    P = pd.DataFrame(per).T.astype(float)
    print("\n[A] early -> late  (top100 scores; linear scale, so correlations are the statistic)")
    print(P.round(3).to_string())
    print(f"\n  {'relation':<28}{'pearson':>10}{'spearman':>10}")
    res["A"] = {"per_donor": P.to_dict(), "corr": {}}
    for a, b in [("early_cd13", "late"), ("early_all", "late"), ("early_ssea4", "late"),
                 ("donor_age", "late"), ("donor_age", "early_cd13")]:
        r = delf.pearson(P[a], P[b])
        rs = float(pd.Series(P[a]).corr(pd.Series(P[b]), method="spearman"))
        print(f"  {a + ' ~ ' + b:<28}{r:>10.3f}{rs:>10.3f}")
        res["A"]["corr"][f"{a}~{b}"] = {"pearson": r, "spearman": rs}
    p1 = delf.partial_corr(P["early_cd13"], P["late"], P["donor_age"])
    p2 = delf.partial_corr(P["early_all"], P["late"], P["donor_age"])
    print(f"\n  PARTIAL early_cd13 ~ late | donor_age = {p1:+.3f}   "
          f"(bar |r| > {delf.PARTIAL_CRIT_DF3}, df=3)   was -0.064 on the dense clock")
    print(f"  PARTIAL early_all  ~ late | donor_age = {p2:+.3f}")
    flipped_a = bool(np.isfinite(p1) and abs(p1) > delf.PARTIAL_CRIT_DF3)
    print(f"  -> {'CHANGED' if flipped_a else 'UNCHANGED'}")
    res["A"]["partial_early_cd13"] = p1
    res["A"]["partial_early_all"] = p2
    res["A"]["verdict"] = "CHANGED" if flipped_a else "UNCHANGED"

    # ---- B: late residual from early EXPRESSION, plus the 9-variant robustness sweep ------- #
    print(f"\n[B] late residual | donor age, from early EXPRESSION  (LOO + "
          f"{dre.N_PERM}-draw permutation null)")
    panel = [g for g in GenePanel.load(str(dre.PANEL_PATH)).genes if g in set(genes)]
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi[g] for g in panel])
    donors = sorted(per)
    X = np.vstack([expr[df[(df.donor == d) & (df.day >= EARLY_LO) & (df.day <= EARLY_HI)].j
                        .to_numpy(), :][:, cols].mean(0) for d in donors])
    y = np.array([per[d]["late"] for d in donors], float)
    age = np.array([per[d]["donor_age"] for d in donors], float)
    print(f"  donors {donors}   features {X.shape[1]}")
    print(f"  {'alpha':>10}{'LOO spearman':>15}{'null p95':>11}{'pctile':>9}   pass")
    n_pass, per_a = 0, {}
    for a in dre.ALPHAS:
        obs = dre.loo_spearman(X, y, age, a, residualise=True)
        null = dre.permutation_null(X, y, age, a, residualise=True)
        p95 = float(np.percentile(null, dre.PERM_PCTILE))
        pct = float((null < obs).mean() * 100) if np.isfinite(obs) else float("nan")
        ok = bool(np.isfinite(obs) and obs > p95)
        n_pass += ok
        per_a[str(a)] = {"observed": obs, "null_p95": p95, "percentile": pct, "pass": ok}
        print(f"  {a:>10.0f}{obs:>15.3f}{p95:>11.3f}{pct:>9.1f}   {'YES' if ok else 'no'}")
    base_verdict = "SIGNAL" if n_pass >= dre.MIN_ALPHAS_PASSING else "NULL"
    print(f"  -> baseline {base_verdict} ({n_pass}/{len(dre.ALPHAS)} alphas)")
    res["B"] = {"baseline": {"alphas": per_a, "n_pass": n_pass, "verdict": base_verdict}}

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_phase1_top100_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
