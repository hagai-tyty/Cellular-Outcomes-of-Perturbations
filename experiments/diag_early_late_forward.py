"""Is the early->late signal REPROGRAMMING, or just donor chronological age?  (read-only)

    python experiments/diag_early_late_forward.py

BACKGROUND
----------
`diag_clock_circularity` showed the same-timepoint ΔAge regression is CIRCULAR: the label is a
linear readout of the model's own input. An early->late formulation escapes that, because the
target is measured from a DIFFERENT sample taken later, which the input cannot contain.

A first look found early CD13 clock age predicting the late plateau at Pearson +0.830 (n=6) on RAW
ages -- no control subtraction, so the shared-zero-point artifact cannot explain it.

But BOTH ends track donor chronological age (late plateau vs donor age is nearly monotone). If the
whole correlation is mediated by donor age, this predicts the DONOR, not the reprogramming --
the same defect as the circularity, one level up. That is the question here.

WHAT IS MEASURED
  r(early, late)                 the raw forward correlation
  r(early, late | donor_age)     partial correlation, donor age removed
  r(donor_age, late | early)     the reverse: does donor age add anything beyond early?

  With n=6 and one covariate the partial has df = 3, so an inconclusive answer is a LIKELY and
  legitimate outcome. It would still be informative: it bounds how many donors would be needed.

PRE-REGISTERED TRANSITION-TIMING RULE  (written BEFORE it was computed)
----------------------------------------------------------------------
Stated as a rule rather than chosen after inspection, because the impression that "older donors
transition later" was formed by EYE, after seeing the trajectories. That makes it a post-hoc
hypothesis. This gives it a fair test; it cannot confirm it. Confirmation needs donors not used to
form it.

  For each donor, over SSEA4 samples only (the only marker present in the late window), ordered by
  day:
     mid    = ( mean(age over days <= EARLY_HI) + mean(age over days >= LATE_LO) ) / 2
     T_day  = the FIRST measured day whose age is below `mid` AND from which EVERY later measured
              day is also below `mid`.
  "Sustained" is load-bearing: the early trajectory swings 40-50 yr between adjacent single
  samples (e.g. O1 d9=121 -> d11=77), so a first-crossing rule without it would report noise.
  Undefined (nan) if the final measured day is not below `mid`.

  HYPOTHESIS UNDER TEST: Spearman(donor chronological age, T_day) > 0.
  At n=6 the critical |Spearman| for p<0.05 is 0.886. Anything below that is NOT a result.

NOT CLAIMED: that a surviving partial correlation would make this a working tool. n=6 donors is
the binding constraint on everything here.
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
from cellfate.data.normalize import normalize_counts  # noqa: E402
from cellfate.data.sources import GillReprogrammingSource  # noqa: E402

GILL_EXPR = r"D:\Gill\GSE165176_Log2_RPM_Sendai_reprogramming (1).txt.gz"
GILL_SERIES = r"D:\Gill\GSE165176_series_matrix.txt.gz"
CLOCK_PATH = ROOT / "configs" / "clocks" / "fleischer_clock.json"

EARLY_LO, EARLY_HI = 7.0, 29.0
LATE_LO = 34.0
SPEARMAN_CRIT_N6 = 0.886      # |rho| needed for p<0.05, two-tailed, n=6
PEARSON_CRIT_N6 = 0.811       # |r| needed for p<0.05, two-tailed, n=6
PARTIAL_CRIT_DF3 = 0.878      # |r| needed for p<0.05, two-tailed, df=3 (n=6, one covariate)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def partial_corr(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """corr(x, y) with z removed from both. nan when any input is degenerate or the
    denominator vanishes (z explaining one variable completely)."""
    rxy, rxz, ryz = pearson(x, y), pearson(x, z), pearson(y, z)
    if not all(np.isfinite([rxy, rxz, ryz])):
        return float("nan")
    den = np.sqrt(max(0.0, 1 - rxz ** 2) * max(0.0, 1 - ryz ** 2))
    if den < 1e-12:
        return float("nan")
    return float((rxy - rxz * ryz) / den)


def transition_day(days, ages, early_hi: float = EARLY_HI,
                   late_lo: float = LATE_LO) -> float:
    """The PRE-REGISTERED rule. First SUSTAINED crossing of the donor's own early/late midpoint."""
    d, a = np.asarray(days, float), np.asarray(ages, float)
    o = np.argsort(d)
    d, a = d[o], a[o]
    early, late = a[d <= early_hi], a[d >= late_lo]
    if len(early) == 0 or len(late) == 0:
        return float("nan")
    mid = (float(early.mean()) + float(late.mean())) / 2.0
    below = a < mid
    for i in range(len(d)):
        if below[i:].all():
            return float(d[i])
    return float("nan")


def load_ages() -> pd.DataFrame:
    src = GillReprogrammingSource(GILL_EXPR, GILL_SERIES)
    src.bulk_integrity_gate = True
    src._load()
    clock = LinearClock.from_json(str(CLOCK_PATH))
    rpm = src._rpm
    # normalize_counts applies CP10k AND log1p (normalize.py:29), matching the clock's declared
    # `log1p_cp10k`. NOT wrapped in a second log1p -- that bug has cost this project a run before.
    expr = normalize_counts(np.clip(rpm.to_numpy(dtype=np.float64).T, 0.0, None), target_sum=1e4)
    ages = clock.predict_age(expr, list(rpm.index))
    rows = []
    for j, c in enumerate(rpm.columns):
        m = src._meta[c]
        mk = "CD13" if "_CD13_" in c else ("SSEA4" if "_SSEA4_" in c else "Fib")
        rows.append({"donor": m["donor"], "day": float(m["day"]), "marker": mk,
                     "donor_age": float(m["age"]) if m.get("age") is not None else float("nan"),
                     "age": float(ages[j])})
    return pd.DataFrame(rows)


def main() -> None:
    df = load_ages()
    per = {}
    for d, g in df.groupby("donor"):
        e = g[(g.day >= EARLY_LO) & (g.day <= EARLY_HI)]
        late = g[g.day >= LATE_LO]
        s = g[g.marker == "SSEA4"].groupby("day").age.mean()
        per[d] = {
            "donor_age": float(g.donor_age.dropna().iloc[0]) if g.donor_age.notna().any() else np.nan,
            "early_cd13": e[e.marker == "CD13"].age.mean(),
            "early_ssea4": e[e.marker == "SSEA4"].age.mean(),
            "early_all": e.age.mean(),
            "late": late.age.mean(),
            "T_day": transition_day(s.index.to_numpy(), s.to_numpy()),
        }
    P = pd.DataFrame(per).T.astype(float)
    print("=" * 96)
    print("EARLY -> LATE FORWARD SIGNAL  (raw clock ages, NO control subtraction)")
    print("=" * 96)
    print(P.round(2).to_string())

    print(f"\nCORRELATIONS  (n={len(P)};  |r| must exceed {PEARSON_CRIT_N6} for p<0.05)")
    pairs = [("early_cd13", "late"), ("early_all", "late"), ("early_ssea4", "late"),
             ("donor_age", "late"), ("donor_age", "early_cd13")]
    res: dict = {"per_donor": P.to_dict(), "corr": {}, "partial": {}}
    for a, b in pairs:
        r = pearson(P[a], P[b])
        rs = float(pd.Series(P[a]).corr(pd.Series(P[b]), method="spearman"))
        flag = "*" if abs(r) > PEARSON_CRIT_N6 else " "
        print(f"  {a:<12} ~ {b:<12} pearson {r:+.3f}{flag}  spearman {rs:+.3f}")
        res["corr"][f"{a}~{b}"] = {"pearson": r, "spearman": rs}

    print(f"\nPARTIAL CORRELATIONS  (df=3;  |r| must exceed {PARTIAL_CRIT_DF3} for p<0.05)")
    p1 = partial_corr(P["early_cd13"], P["late"], P["donor_age"])
    p2 = partial_corr(P["donor_age"], P["late"], P["early_cd13"])
    print(f"  early_cd13 ~ late | donor_age  = {p1:+.3f}"
          f"   {'clears' if abs(p1) > PARTIAL_CRIT_DF3 else 'DOES NOT clear'} the threshold")
    print(f"  donor_age  ~ late | early_cd13 = {p2:+.3f}"
          f"   {'clears' if abs(p2) > PARTIAL_CRIT_DF3 else 'DOES NOT clear'} the threshold")
    res["partial"] = {"early_given_age": p1, "age_given_early": p2}

    print(f"\nPRE-REGISTERED TIMING TEST  (|spearman| must exceed {SPEARMAN_CRIT_N6} for p<0.05)")
    print("  T_day per donor:", {k: (None if np.isnan(v) else v) for k, v in P["T_day"].items()})
    v = P[["donor_age", "T_day"]].dropna()
    if len(v) >= 3:
        rs = float(pd.Series(v.donor_age).corr(pd.Series(v.T_day), method="spearman"))
        ok = abs(rs) > SPEARMAN_CRIT_N6
        print(f"  spearman(donor_age, T_day) = {rs:+.3f}  (n={len(v)})  ->  "
              f"{'PASSES' if ok else 'DOES NOT PASS'} the pre-registered bar")
        print("  NOTE: post-hoc hypothesis. Even a pass is hypothesis-GENERATING; confirming it "
              "needs donors\n        that were not used to form it.")
        res["timing"] = {"spearman": rs, "n": len(v), "passes": bool(ok)}

    _RESULTS.mkdir(exist_ok=True)
    out = _RESULTS / "diag_early_late_forward_results.json"
    out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
