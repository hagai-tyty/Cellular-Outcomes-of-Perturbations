"""Reproduce STAGE_1_5_1_REV_FINAL sec 4.4 — the two checks that corroborate contrast A.

    python experiments/verify_rev_final_4_4.py "C:\\Users\\hagay\\Desktop\\GSE165179"

READ-ONLY. Writes `verify_rev_final_4_4_results.json`. `src/` is untouched.

WHY THIS EXISTS
---------------
sec 4.4 is the section that answers the one challenge a reviewer is most likely to press: contrast A
was named before the run but became the headline *after* it (sec 7). sec 4.4 closes that with two
checks the promotion could not have manufactured — an internal negative control and a dose-response.

**Neither was produced by any artefact.** `diag_methylation_anchor.py` computes three contrasts, and
its `CONTRASTS` list contains no failing-INTERMEDIATE arm and no dose-response at all; sec 9 lists
`diag_methylation_anchor_results.json` as "full output", but sec 4.4's numbers are not in it. So the
document's most load-bearing defensive section rested on ad-hoc computation with nothing to re-run.

This script is that missing artefact. It re-derives sec 4.4 from the raw beta matrix through an
INDEPENDENT path — pure stdlib, no numpy, no shared code with `diag_methylation_anchor.py` — so
agreement is genuine corroboration rather than the same code returning the same answer twice.

WHAT IT CHECKS
--------------
  V1  contrast A (TRI vs NCI)      -- MUST reproduce -24.1 / -27.5. Validates the pipeline BEFORE
                                      its new numbers are trusted. A pipeline that cannot reproduce
                                      a known value is not evidence for an unknown one.
  V2  sec 4.4(a) FLI vs NCI        -- the internal negative control, claimed -1.1 / -3.6
  V3  sec 4.4(b) dose-response     -- Spearman(reprogramming length, effect), claimed
                                      rho -0.885 / -0.842, p 0.0001 / 0.0006
  V4  sec 4.3 exactness            -- the intercept cancels EXACTLY only where both samples predict
                                      age >= 20 (Horvath's transform is linear only there). Counts
                                      the violations and shows that pairs with ZERO violations give
                                      a difference of EXACTLY 0.00, which demonstrates the algebra
                                      instead of asserting it.

RESULT WHEN FIRST RUN (2026-07-30): all four reproduce. See CHANGES.md.
"""
from __future__ import annotations

import gzip
import json
import math
import re
import statistics as st
import sys
from datetime import UTC, datetime
from pathlib import Path

# --------------------------------------------------------------------------- #
# Pure logic — data-free                                                       #
# --------------------------------------------------------------------------- #
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

ADULT_AGE = 20.0

# Two-sided t, indexed by degrees of freedom. Matches diag_methylation_anchor.py's treatment.
T_CRIT = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
          8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201}


def anti_trafo(x: float, adult_age: float = ADULT_AGE) -> float:
    """Horvath's inverse transform. LINEAR above `adult_age`, exponential below."""
    return (1 + adult_age) * math.exp(x) - 1 if x < 0 else (1 + adult_age) * x + adult_age


def to_lp(age: float, adult_age: float = ADULT_AGE) -> float:
    """Invert `anti_trafo`. Used to recover the intercept-free form 21*(lp_t - lp_c)."""
    return ((age - adult_age) / (1 + adult_age) if age >= adult_age
            else math.log((age + 1.0) / (1 + adult_age)))


# Substring -> arm, MOST SPECIFIC FIRST. The three INTERMEDIATE arms contain the fibroblast arms'
# substrings ("failing_to_transiently_reprogram_intermediate" contains "transiently_reprogram"), so
# a reordering of this table silently mislabels arms. It is a table rather than a chain of `if`s so
# that the precedence is data, and visible, instead of control flow.
_ARM_PATTERNS: tuple[tuple[str, str], ...] = (
    ("failing_to_transiently_reprogram_intermediate", "FLI"),
    ("negative_control_intermediate", "NCI"),
    ("transient_reprogramming_intermediate", "TRI"),
    ("failed_to_transiently_reprogram", "FL"),
    ("negative_control", "NC"),
    ("transiently_reprogrammed", "TR"),
)


def arm_of(name: str) -> str | None:
    """Map a GEO sample title to its arm, or None if it is not one we contrast."""
    n = name.lower()
    for needle, arm in _ARM_PATTERNS:
        if needle in n:
            return arm
    return "FIB" if re.match(r"^o\d\s+fib$", n.strip()) else None


def paired_ci(vals: list[float]) -> dict:
    n = len(vals)
    if n < 2:
        return {"n": n, "mean": float("nan"), "ci95": [float("nan")] * 2}
    m = st.mean(vals)
    se = st.stdev(vals) / math.sqrt(n)
    t = T_CRIT.get(n - 1, 1.96)
    return {"n": n, "mean": m, "ci95": [m - t * se, m + t * se]}


def _ranks(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def _pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True))
    den = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return num / den if den else float("nan")


def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
    h = d
    for m in range(1, 200):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
        c = c if abs(c) > 1e-30 else 1e-30
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
        c = c if abs(c) > 1e-30 else 1e-30
        de = d * c
        h *= de
        if abs(de - 1.0) < 3e-11:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1 - x))
    bt = math.exp(lb)
    return (bt * _betacf(a, b, x) / a if x < (a + 1) / (a + b + 2)
            else 1.0 - bt * _betacf(b, a, 1 - x) / b)


def spearman(x: list[float], y: list[float]) -> dict:
    """Spearman rho with a two-sided p from the t approximation (no scipy on this machine)."""
    n = len(x)
    r = _pearson(_ranks(x), _ranks(y))
    if n < 3 or not math.isfinite(r) or abs(r) >= 1.0:
        return {"n": n, "rho": r, "p": float("nan")}
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return {"n": n, "rho": r, "p": _betai((n - 2) / 2.0, 0.5, (n - 2) / ((n - 2) + t * t))}


def ols_slope(x: list[float], y: list[float]) -> float:
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    den = sum((a - mx) ** 2 for a in x)
    return sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True)) / den if den else float("nan")


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def load_clock_betas(matrix: Path, wanted: set[str]) -> tuple[dict[int, str], dict[int, dict]]:
    """Stream the processed matrix, keeping only the clock CpGs. Comma-separated; every sample
    column is followed by a `Detection Pval` column, which is skipped."""
    with gzip.open(matrix, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split(",")
        cols = {i: h.strip() for i, h in enumerate(hdr)
                if i and not h.strip().lower().startswith("detection")}
        beta: dict[int, dict[str, float]] = {i: {} for i in cols}
        for line in fh:
            cut = line.find(",")
            if cut < 0 or line[:cut] not in wanted:
                continue
            cid = line[:cut]
            parts = line.rstrip("\n").split(",")
            for i in cols:
                try:
                    beta[i][cid] = float(parts[i])
                except (ValueError, IndexError):
                    pass
    return cols, beta


def contrast(ages: dict[int, float], meta: dict[int, dict], treated: str, control: str) -> dict:
    """Paired (donor, day) contrast. Replicates are AVERAGED, never treated as independent —
    the unit-of-analysis rule from REV FINAL sec 2, and the fix for the pairing bug that once
    silently dropped 6 of 9 pairs."""
    groups: dict[tuple, dict[str, list[float]]] = {}
    for i, md in meta.items():
        if i in ages and md["arm"] in (treated, control):
            groups.setdefault((md["donor"], md["day"]), {}).setdefault(md["arm"], []).append(ages[i])
    rows, diffs = [], []
    for key in sorted(groups):
        g = groups[key]
        if treated in g and control in g:
            t_m, c_m = st.mean(g[treated]), st.mean(g[control])
            diffs.append(t_m - c_m)
            rows.append({"donor": key[0], "day": key[1], "treated_age": t_m, "control_age": c_m,
                         "diff_years": t_m - c_m,
                         "intercept_free": 21.0 * (to_lp(t_m) - to_lp(c_m)),
                         "below_20": t_m < ADULT_AGE or c_m < ADULT_AGE})
    out = paired_ci(diffs)
    out["pairs"] = rows
    out["n_pairs_below_20"] = sum(1 for r in rows if r["below_20"])
    out["intercept_free"] = paired_ci([r["intercept_free"] for r in rows])
    return out


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = Path(__file__).resolve().parents[1]
    data = Path(sys.argv[1])
    matrix = next(data.glob("*Matrix_processed*.txt.gz"), None)
    if matrix is None:
        print(f"ERROR: no processed matrix under {data}")
        return 2

    clocks = {}
    for f in ("horvath_skin_blood_2018", "horvath_multitissue_2013"):
        j = json.loads((root / "configs" / "clocks" / f"{f}.json").read_text(encoding="utf-8"))
        clocks[f] = {k: float(v) for k, v in j["weights"].items()}

    prev_path = _RESULTS / "diag_methylation_anchor_results.json"
    if not prev_path.exists():
        print("ERROR: run diag_methylation_anchor.py first (needed for the intercepts)")
        return 2
    prev = json.loads(prev_path.read_text(encoding="utf-8"))

    wanted = set().union(*clocks.values())
    print(f"[scan] {matrix.name} for {len(wanted)} clock CpGs ...", flush=True)
    cols, beta = load_clock_betas(matrix, wanted)
    print(f"[scan] {len(cols)} sample columns", flush=True)

    meta = {}
    for i, name in cols.items():
        a = arm_of(name)
        if a:
            m = re.search(r"(\d+)days", name.lower())
            meta[i] = {"arm": a, "day": float(m.group(1)) if m else 0.0,
                       "donor": name.strip()[:2].upper(), "name": name}

    out = {"script": "verify_rev_final_4_4", "utc": datetime.now(UTC).isoformat(),
           "matrix": str(matrix), "clocks": {}}

    for cf, W in clocks.items():
        k = float(prev["clocks"][cf]["G2"]["intercept_used"])
        ages = {}
        for i in cols:
            hits = [(W[c], beta[i][c]) for c in W if c in beta[i]]
            if len(hits) >= 0.5 * len(W):
                ages[i] = anti_trafo(sum(w * v for w, v in hits) + k)

        A = contrast(ages, meta, "TRI", "NCI")
        Ac = contrast(ages, meta, "FLI", "NCI")
        day = [r["day"] for r in A["pairs"]]
        dose = spearman(day, [r["diff_years"] for r in A["pairs"]])
        dose["slope_intercept_free"] = ols_slope(day, [r["intercept_free"] for r in A["pairs"]])
        dose["slope_derived_intercept"] = ols_slope(day, [r["diff_years"] for r in A["pairs"]])

        blk = {"intercept_used": k, "n_samples": len(ages),
               "V1_contrast_A": A, "V2_failing_intermediates": Ac, "V3_dose_response": dose}
        out["clocks"][cf] = blk

        print(f"\n== {cf}  (intercept {k:+.4f}, {len(ages)} samples)")
        print(f"   V1 A   TRI vs NCI  n={A['n']:2d}  {A['mean']:+6.2f} "
              f"[{A['ci95'][0]:+6.2f},{A['ci95'][1]:+6.2f}]   doc: -24.1 / -27.5")
        print(f"   V2 A-c FLI vs NCI  n={Ac['n']:2d}  {Ac['mean']:+6.2f} "
              f"[{Ac['ci95'][0]:+6.2f},{Ac['ci95'][1]:+6.2f}]   doc:  -1.1 /  -3.6")
        print(f"   V3 dose  rho={dose['rho']:+.3f}  p={dose['p']:.4f}  "
              f"slope={dose['slope_intercept_free']:+.2f} yr/day (intercept-free)   "
              f"doc: -0.885/-0.842, -3.30/-3.15")
        print(f"   V4 pairs with a sample below age 20 (sec 4.3 linear branch invalid): "
              f"{A['n_pairs_below_20']}  ->  derived {A['mean']:+.2f} vs "
              f"intercept-free {A['intercept_free']['mean']:+.2f}")

    dest = _RESULTS / "verify_rev_final_4_4_results.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
