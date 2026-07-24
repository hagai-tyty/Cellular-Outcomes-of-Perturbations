"""STAGE 1.5 PHASE 1 — is the ±12.7 yr per-donor offset biology, batch, or one noisy baseline?

    python experiments/diag_zero_point.py                 # defaults to D:\\Gill
    python experiments/diag_zero_point.py "D:\\Gill"

READ-ONLY. Loads the Gill series + expression matrix, runs the FROZEN clock on the six day-0
baselines, and writes `diag_zero_point_results.json`. Nothing is rebuilt, refitted or written back;
`src/` is not touched. Phase 1 exists to decide whether the expensive Phase 3 is needed at all.

WHY (STAGE_1_5_HARMONIZATION_AUDIT.md §5.2-§5.3, tightened in §6.2)
------------------------------------------------------------------
Group E ruled out the self-centring fallback. It surfaced something else: every Gill donor's
zero-point is **one unreplicated day-0 sample**, and the frozen clock's own cross-validated error is
`cv_mae = 12.27 yr` (`configs/clocks/fleischer_clock.json`) against a per-donor offset of
**±12.7 yr**. Any error in that single baseline propagates 1:1 into every ΔAge for that donor, i.e.
lands as exactly the per-donor additive offset Stage 2 is premised on. The two are currently
indistinguishable. Three measurements separate them:

  M1  does the clock read chronological age on THIS data?   (GEO donor ages: 0,0,29,35,53,53)
  M2  is there an Exp1/Exp2 batch effect?                    (all six baselines are Exp2)
  M3  how much of the offset VARIANCE could one noisy baseline explain?

BARS ARE PRE-REGISTERED AND RESOLVABILITY-CHECKED BEFORE RUNNING (ground rule §5b, tightening T1).
Each bar below is stated with the null it is tested against AND the rate at which a system that
meets the intent would pass it. A measurement whose bar a correct system cannot pass is a
description, not a test -- `bars()` reports that for all three, and `tests/test_bars_resolvable.py`
holds the entries.

M1 is the one that can escalate past this whole stage: if the clock cannot separate a 0-year-old
donor from a 53-year-old one, ΔAge's target is unvalidated and that reaches into Stage 4.
"""

from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Pure logic — data-free, fully unit-tested; nothing below imports repo data.  #
# --------------------------------------------------------------------------- #

# GEO GSE165176 declared donor chronological ages. NOT currently parsed by the pipeline
# (`_parse_series` reads only `days of reprogramming` and `cell type`) -- finding D3.
DONOR_AGE: dict[str, float] = {"N2": 0.0, "N3": 0.0, "Y1": 29.0, "Y2": 35.0, "O1": 53.0, "O2": 53.0}

# The frozen clock's own cross-validated error (configs/clocks/fleischer_clock.json,
# Fleischer 2018 / GSE113957, 133 samples). Every bar below is scaled by it.
CLOCK_CV_MAE = 12.26879346460328

MIN_PAIRS_FOR_M2 = 3          # below this a paired CI is not worth stating
M3_MOSTLY = 0.50              # share of offset variance that counts as "baseline noise dominates"
M3_LITTLE = 0.25


def m1_threshold(cv_mae: float = CLOCK_CV_MAE, n_young: int = 2, n_old: int = 2,
                 alpha: float = 0.05) -> float:
    """Smallest age contrast that beats "the clock reads nothing" at one-sided `alpha`.

    A bar of merely "contrast > 0" would be passed half the time by a clock that reads pure noise,
    so it tests nothing. Under the null the contrast is `N(0, cv_mae*sqrt(1/n_young + 1/n_old))`;
    the bar is `z_alpha` of those standard errors.
    """
    se = cv_mae * np.sqrt(1.0 / n_young + 1.0 / n_old)
    z = 1.6448536269514722          # Phi^-1(0.95); scipy not needed for a fixed alpha
    return float(z * se)


def m1_power(true_contrast: float, cv_mae: float = CLOCK_CV_MAE,
             n_young: int = 2, n_old: int = 2, alpha: float = 0.05) -> float:
    """P(pass) for a clock that IS correct — the §5b resolvability figure for M1."""
    se = cv_mae * np.sqrt(1.0 / n_young + 1.0 / n_old)
    from math import erf, sqrt
    z = (true_contrast - m1_threshold(cv_mae, n_young, n_old, alpha)) / se
    return float(0.5 * (1.0 + erf(z / sqrt(2.0))))


def m1_verdict(pred: dict[str, float], chrono: dict[str, float] | None = None,
               cv_mae: float = CLOCK_CV_MAE) -> dict:
    """Does the clock separate the age extremes? `pred` maps donor -> predicted age.

    Only the EXTREME contrast is gated. The middle donors (29 vs 35) differ by half the clock's
    error and are deliberately underpowered — reported, never claimed.
    """
    chrono = chrono or DONOR_AGE
    common = [d for d in pred if d in chrono]
    if not common:
        return {"status": "CANNOT_VERIFY", "reason": "no donors overlap the age table"}
    ages = np.array([chrono[d] for d in common], float)
    vals = np.array([pred[d] for d in common], float)
    young, old = vals[ages == ages.min()], vals[ages == ages.max()]
    if len(young) == 0 or len(old) == 0 or ages.min() == ages.max():
        return {"status": "CANNOT_VERIFY", "reason": "no age contrast among the donors present"}

    contrast = float(old.mean() - young.mean())
    thr = m1_threshold(cv_mae, len(young), len(old))
    true_gap = float(ages.max() - ages.min())
    # Spearman over ALL donors: reported, NOT gated (n=6, and 29 vs 35 is half the clock error)
    order_pred = np.argsort(np.argsort(vals))
    order_true = np.argsort(np.argsort(ages))
    rho = float(np.corrcoef(order_pred, order_true)[0, 1]) if len(common) > 2 else float("nan")

    return {
        "status": "PASS" if contrast >= thr else "FAIL",
        "contrast_years": contrast,
        "threshold_years": float(thr),
        "true_age_gap": true_gap,
        "n_young": int(len(young)), "n_old": int(len(old)),
        "power_if_clock_correct": m1_power(true_gap, cv_mae, len(young), len(old)),
        "spearman_all_donors_REPORTED_NOT_GATED": rho,
        "reason": (f"extreme contrast {contrast:.1f} yr vs bar {thr:.1f} yr "
                   f"(true gap {true_gap:.0f} yr, clock cv_mae {cv_mae:.2f})"),
    }


def parse_title(title: str) -> dict | None:
    """Split a Gill sample title into `(donor, day, marker, exp)`. Pure.

    Titles look like `N2_d11_CD13_Sendai_Exp1` (treatment) and `N2_Fib_Sendai_Exp2` (the day-0
    baseline). The pipeline's `obs` discards batch and marker identity (finding D1), so the
    series-matrix titles are the ONLY place they survive — which is why M2 reads them here rather
    than from `obs`.

    Returns `None` for anything not fully specified — notably the `_Fib_` baselines, which carry no
    day or marker and so cannot form a matched pair.
    """
    p = title.split("_")
    if len(p) < 4:
        return None
    exp = p[-1]
    if not exp.upper().startswith("EXP"):
        return None
    donor, day_tok, marker = p[0], p[1], p[2]
    if not (day_tok.startswith("d") and day_tok[1:].isdigit()):
        return None                      # excludes the '_Fib_' day-0 baselines by construction
    return {"donor": donor, "day": float(day_tok[1:]), "marker": marker, "exp": exp}


def group_matched_pairs(parsed: dict[str, dict]) -> list[tuple[str, str]]:
    """`(donor, day, marker)` groups holding BOTH an Exp1 and an Exp2 sample. Pure.

    Returns `(exp1_title, exp2_title)` so the caller can difference them in that order. Groups with
    only one batch, or with duplicates within a batch, are skipped — a pair must be unambiguous.
    """
    buckets: dict[tuple[str, float, str], dict[str, list[str]]] = {}
    for title, v in parsed.items():
        key = (v["donor"], v["day"], v["marker"])
        buckets.setdefault(key, {}).setdefault(v["exp"], []).append(title)
    pairs: list[tuple[str, str]] = []
    for _key, by_exp in sorted(buckets.items()):
        e1, e2 = by_exp.get("Exp1", []), by_exp.get("Exp2", [])
        if len(e1) == 1 and len(e2) == 1:
            pairs.append((e1[0], e2[0]))
    return pairs


def m2_verdict(diffs: list[float], min_pairs: int = MIN_PAIRS_FOR_M2) -> dict:
    """Exp1 − Exp2 offset from matched `(donor, day, marker)` pairs.

    Tightening T4: the pair COUNT is reported first and gates everything. If no matched pairs
    exist the Exp1↔Exp2 offset is unidentifiable, and fix option (a) is off the menu regardless of
    anything else — that is a legitimate, plan-permitted outcome, not a failure to measure.
    """
    d = np.asarray([x for x in diffs if np.isfinite(x)], float)
    n = len(d)
    if n == 0:
        return {"status": "NOT_ESTIMABLE", "n_pairs": 0,
                "reason": "no matched (donor, day, marker) pairs span Exp1 and Exp2; the batch "
                          "offset is unidentifiable and fix option (a) is impossible"}
    if n < min_pairs:
        return {"status": "INDETERMINATE", "n_pairs": n,
                "reason": f"only {n} matched pair(s) (< {min_pairs}); a paired CI is not worth stating"}
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    se = sd / np.sqrt(n)
    t = {2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571}.get(n - 1, 1.96)
    lo, hi = mean - t * se, mean + t * se
    return {
        "status": "BATCH_EFFECT" if (lo > 0 or hi < 0) else "NO_BATCH_EFFECT",
        "n_pairs": n, "mean_offset_years": mean, "ci95": [float(lo), float(hi)],
        "reason": (f"paired Exp1-Exp2 offset {mean:+.2f} yr, 95% CI [{lo:+.2f}, {hi:+.2f}] "
                   f"over {n} pairs"),
    }


def m3_verdict(offsets: list[float], cv_mae: float = CLOCK_CV_MAE) -> dict:
    """What share of the per-donor offset VARIANCE could one unreplicated baseline explain?

    A single day-0 sample carries the clock's own error, so it contributes `cv_mae**2` of variance
    to the per-donor offset by construction. Compare that with the observed spread.

    The chi-square CI is load-bearing, not decoration: with 6 donors the variance is known only to
    within roughly a factor of six, so the point estimate alone is not a finding (tightening T2).
    """
    o = np.asarray([x for x in offsets if np.isfinite(x)], float)
    n = len(o)
    if n < 3:
        return {"status": "CANNOT_VERIFY", "n_donors": n,
                "reason": f"{n} donors; a variance needs at least 3"}
    var_obs = float(o.var(ddof=1))
    if var_obs <= 0:
        return {"status": "CANNOT_VERIFY", "n_donors": n, "reason": "zero observed variance"}
    var_base = float(cv_mae ** 2)
    share = float(min(var_base / var_obs, 1.0))
    # chi-square bounds on the observed variance -> bounds on the share
    chi = {3: (0.216, 9.348), 4: (0.484, 11.143), 5: (0.831, 12.833), 6: (1.237, 14.449),
           7: (1.690, 16.013), 8: (2.180, 17.535)}.get(n - 1, (0.831, 12.833))
    var_hi = var_obs * (n - 1) / chi[0]
    var_lo = var_obs * (n - 1) / chi[1]
    share_lo = float(min(var_base / var_hi, 1.0))
    share_hi = float(min(var_base / var_lo, 1.0))
    wide = (share_hi - share_lo) > 0.5

    if wide:
        status = "INDETERMINATE"
    elif share >= M3_MOSTLY:
        status = "BASELINE_DOMINATES"
    elif share < M3_LITTLE:
        status = "BASELINE_MINOR"
    else:
        status = "INDETERMINATE"
    resid = float(np.sqrt(max(var_obs - var_base, 0.0)))
    return {
        "status": status, "n_donors": n,
        "observed_sd_years": float(np.sqrt(var_obs)),
        "baseline_sd_years": float(cv_mae),
        "share_of_variance": share, "share_ci": [share_lo, share_hi],
        "residual_sd_years": resid,
        "reason": (f"observed offset SD {np.sqrt(var_obs):.1f} yr vs {cv_mae:.1f} yr from a single "
                   f"baseline -> share {share:.0%} (95% CI {share_lo:.0%}-{share_hi:.0%}); "
                   f"{resid:.1f} yr SD remains for biology+batch+model"),
    }


def decide(m1: dict, m2: dict, m3: dict) -> dict:
    """Fold M1/M2/M3 into the one action Phase 1 licenses (§5.4 table, extended by T2)."""
    if m1.get("status") == "FAIL":
        return {"action": "ESCALATE",
                "reason": "the clock does not separate the age extremes on this data, so ΔAge's "
                          "target is unvalidated. This reaches past Stage 1.5 into Stage 4 and "
                          "Stage 2's premise is void as stated. Do not proceed to Phase 2/3."}
    if m1.get("status") != "PASS":
        return {"action": "BLOCKED", "reason": f"M1 inconclusive: {m1.get('reason', '')}"}

    batch = m2.get("status") == "BATCH_EFFECT"
    lead = {"BASELINE_DOMINATES": "(b) shrinkage baseline — the n=1 variance is the larger problem",
            "BASELINE_MINOR": "(a) remove the Exp1/Exp2 offset — the batch term leads",
            }.get(m3.get("status"), "undetermined by M3; let M2 lead")
    if not batch and m3.get("status") == "BASELINE_MINOR":
        return {"action": "PHASE_2_ONLY",
                "reason": "baselines are informative and unconfounded; instrument them and let "
                          "Stage 2 proceed as planned"}
    return {"action": "PHASE_2_AND_3", "phase3_lead": lead,
            "reason": ("a real defect is quantified: "
                       + ("Exp1/Exp2 batch confound; " if batch else "")
                       + f"M3 says {m3.get('status', '?')}. Exactly ONE Phase 3 option ships, as "
                         "its own pre-registered Change (it changes y_age, so it reopens BOTH "
                         "Stage 1 targets and restarts the +0.000 guard record — tightening T4)")}


def bars() -> list[dict]:
    """The pre-registered bars, with the rate a system meeting the intent would pass each.

    Ground rule §5b: this is computed and recorded BEFORE the run, not after. M3's entry is the
    interesting one — it is expected to be UNRESOLVABLE at n=6, which is itself worth knowing in
    advance rather than discovering in the results.
    """
    thr = m1_threshold()
    power = m1_power(53.0)
    return [
        {"id": "M1", "bar": f"extreme age contrast >= {thr:.1f} yr",
         "null": "a clock that reads nothing (contrast ~ N(0, 12.27))",
         "pass_rate_if_intent_holds": power,
         "verdict": "RESOLVABLE" if power >= 0.95 else "UNRESOLVABLE"},
        {"id": "M2", "bar": "paired 95% CI on the Exp1-Exp2 offset excludes 0",
         "null": "no batch effect",
         "pass_rate_if_intent_holds": None,
         "verdict": "CONDITIONAL — depends on the matched-pair count, reported at runtime (T4)"},
        {"id": "M3", "bar": f"share of offset variance >= {M3_MOSTLY:.0%} or < {M3_LITTLE:.0%}",
         "null": "n/a — an estimate with a chi-square CI, not a test",
         "pass_rate_if_intent_holds": None,
         "verdict": "EXPECTED UNRESOLVABLE at n=6 — the variance is known only to ~a factor of "
                    "six, so INDETERMINATE is the pre-registered likely outcome (T2)"},
    ]


# --------------------------------------------------------------------------- #
# Real-data wiring (imports repo machinery only when actually run)            #
# --------------------------------------------------------------------------- #
def baseline_ages(gill_dir: str) -> tuple[dict[str, float], dict]:
    """Predicted age of each donor's single day-0 control, via the production path.

    Uses the SAME normalisation the build uses (`normalize_counts` -> log1p-CP10k), which is the
    space the frozen clock was fitted in, so these ages are comparable to the clock's own CV MAE.
    """
    root = Path(__file__).resolve().parents[1]
    for p in (root, root / "local_runners", root / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    from run_multi_local import discover_gill  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GillReprogrammingSource

    expr, series = discover_gill(gill_dir)
    src = GillReprogrammingSource(expr_tsv=expr, series_matrix=series)
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")

    ages: dict[str, float] = {}
    meta: dict = {"donors": {}, "errors": {}}
    for chunk in src.plan():
        donor = chunk["cell_line"]
        try:
            raw = src.fetch(chunk)
            is_ctrl = np.asarray(raw.obs["is_control"].to_numpy(), dtype=bool)
            if not is_ctrl.any():
                meta["errors"][donor] = "no day-0 control"
                continue
            norm = normalize_counts(raw.counts[is_ctrl])
            a = clock.predict_age(norm, raw.genes)
            ages[donor] = float(np.mean(a))
            meta["donors"][donor] = {"n_baseline": int(is_ctrl.sum()),
                                     "predicted_age": float(np.mean(a)),
                                     "chronological_age": DONOR_AGE.get(donor)}
        except Exception as exc:  # noqa: BLE001 — recorded per donor, never aborts the scan
            meta["errors"][donor] = repr(exc)[:160]
    return ages, meta


def matched_exp_offsets(gill_dir: str) -> tuple[list[float], dict]:
    """Exp1 − Exp2 predicted-age difference for every matched `(donor, day, marker)` pair.

    Batch identity is parsed from the SERIES-MATRIX TITLES, because the pipeline's `obs` discards it
    (finding D1). Same normalisation (`normalize_counts` → log1p-CP10k) and same frozen clock as
    `baseline_ages`, so the differences are in the clock's own units.
    """
    root = Path(__file__).resolve().parents[1]
    for p in (root, root / "local_runners", root / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    from run_multi_local import discover_gill  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GillReprogrammingSource

    expr, series = discover_gill(gill_dir)
    src = GillReprogrammingSource(expr_tsv=expr, series_matrix=series)
    src._load()
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")

    parsed = {t: v for t, v in ((t, parse_title(t)) for t in src._rpm.columns) if v}
    pairs = group_matched_pairs(parsed)
    meta: dict = {"n_titles": int(len(src._rpm.columns)), "n_parsed": len(parsed),
                  "n_pairs": len(pairs)}
    if not pairs:
        return [], meta

    needed = sorted({t for pr in pairs for t in pr})
    mat = src._rpm[needed].to_numpy(dtype=np.float64).T          # samples x genes (linear RPM)
    ages = dict(zip(needed, clock.predict_age(normalize_counts(mat), list(src._genes)),
                    strict=True))
    diffs = [float(ages[a] - ages[b]) for a, b in pairs]
    meta["pairs"] = [{"exp1": a, "exp2": b, "diff_years": float(ages[a] - ages[b])}
                     for a, b in pairs]
    return diffs, meta


def main() -> int:
    gill_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\Gill"
    print("STAGE 1.5 PHASE 1 — zero-point diagnostics (read-only)\n")

    print("  PRE-REGISTERED BARS (ground rule §5b — recorded before the numbers):")
    for b in bars():
        rate = ("n/a" if b["pass_rate_if_intent_holds"] is None
                else f"{b['pass_rate_if_intent_holds']:.1%}")
        print(f"    {b['id']}: {b['bar']}")
        print(f"        vs null: {b['null']}")
        print(f"        a correct system passes: {rate}   -> {b['verdict']}")
    print()

    try:
        ages, meta = baseline_ages(gill_dir)
    except Exception as exc:  # noqa: BLE001
        print(f"  [!] could not load Gill data from {gill_dir}: {exc!r}")
        print("      (pass the directory as the first argument)")
        return 1

    print(f"  {'donor':<7}{'chrono':>8}{'predicted':>11}{'n_baseline':>12}")
    print("  " + "-" * 38)
    for d in sorted(ages, key=lambda x: DONOR_AGE.get(x, -1)):
        m = meta["donors"][d]
        print(f"  {d:<7}{str(m['chronological_age']):>8}{m['predicted_age']:>11.1f}"
              f"{m['n_baseline']:>12}")
    for d, e in meta["errors"].items():
        print(f"  [!] {d}: {e}")

    m1 = m1_verdict(ages)
    # M2: batch identity is absent from `obs` (finding D1) but present in the series-matrix titles,
    # so it is parsed from there and MEASURED. An earlier revision handed in an empty list, which
    # made M2 report "no matched pairs ... option (a) is impossible" -- a false claim, since
    # N2_d11_CD13_Sendai_Exp1 and _Exp2 both exist.
    try:
        diffs, m2_meta = matched_exp_offsets(gill_dir)
    except Exception as exc:  # noqa: BLE001 — a failed M2 must not hide the M1 verdict
        diffs, m2_meta = [], {"error": repr(exc)[:200]}
    m2 = m2_verdict(diffs)
    m2["evidence"] = m2_meta
    # M3 uses the per-donor level shifts already measured by the scorecard.
    shifts = load_level_shifts()
    m3 = m3_verdict(list(shifts.values())) if shifts else {"status": "CANNOT_VERIFY",
                                                          "reason": "no snapshot with level shifts"}
    decision = decide(m1, m2, m3)

    for name, r in (("M1", m1), ("M2", m2), ("M3", m3)):
        print(f"\n  {name}: {r.get('status')}\n      {r.get('reason', '')}")
    print(f"\n  ACTION: {decision['action']}\n      {decision['reason']}")

    out = {"script": "diag_zero_point", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "gill_dir": gill_dir, "bars": bars(), "baselines": meta,
           "M1": m1, "M2": m2, "M3": m3, "decision": decision}
    p = Path("diag_zero_point_results.json")
    p.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"\n  wrote {p}")
    return 0


def load_level_shifts(snapshot: str = "scorecard/B_fatecal_pooled.json") -> dict[str, float]:
    """Per-donor `level_shift_model` from a scorecard snapshot — the ±12.7 yr quantity itself."""
    p = Path(snapshot)
    if not p.exists():
        for alt in ("scorecard/B_fatecal.json", "scorecard/baseline.json"):
            if Path(alt).exists():
                p = Path(alt)
                break
        else:
            return {}
    folds = json.loads(p.read_text(encoding="utf-8")).get("folds", {})
    return {d: float(f["level_shift_model"]) for d, f in folds.items()
            if isinstance(f, dict) and f.get("level_shift_model") is not None}


if __name__ == "__main__":
    raise SystemExit(main())
