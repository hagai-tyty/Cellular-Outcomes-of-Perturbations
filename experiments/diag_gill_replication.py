"""STAGE 1.5.1 VERIFICATION — can we reproduce Gill's ~30 yr rejuvenation with OUR clock?

    python experiments/diag_gill_replication.py "D:\\Gill"

READ-ONLY. Writes `diag_gill_replication_results.json`. `src/` untouched.

WHY. Gill et al. 2022 (eLife 71624) report ~30 yr transcriptomic rejuvenation at days 10-13 of MPTR
in these exact samples. Our E1/E1b/D2 found no trend (or the wrong sign). Before blaming the clock,
two flaws in OUR analysis have to be ruled out -- both found by reading the metadata, not the code:

  F1  CELL-TYPE POOLING. GSE165176 labels each sample `Reprogramming fibroblast` (65) vs
      `Failing to reprogram fibroblast` (47) vs `iPSC` (6) vs `Dermal fibroblast` (6). E1/E1b/D2
      excluded ONLY iPSC, so 47 samples that by definition did NOT reprogram were averaged in with
      the ones that did. Non-responders cannot rejuvenate; pooling them dilutes any real effect.

  F2  WRONG STATISTIC FOR THE SHAPE. Gill's effect is a DIP peaking at day 10-13, after which cells
      re-differentiate (OSKM is withdrawn ~day 13). E1/E1b used a monotonic Spearman over days,
      which is near-blind to a non-monotonic dip: a fall-then-rise can rank-correlate to ~0.

This script tests the contrast Gill actually reported -- age at the PEAK window vs the donor's own
day 0 -- separately per cell type. `Failing to reprogram` is a built-in NEGATIVE CONTROL: if the
clock is reading real rejuvenation, responders should drop and non-responders should not.

No bar is being gated here; this is a replication check of a published result on our pipeline.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "local_runners", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

PEAK_DAYS = (10.0, 13.0)      # Gill's reported optimum; our data has days 11 and 13 in-window
GILL_REPORTED = -30.0
REPROG = "Reprogramming fibroblast"
FAILING = "Failing to reprogram fibroblast"
BASELINE = "Dermal fibroblast"


def contrast(day0: list[float], peak: list[float]) -> dict:
    """Mean predicted-age change from a donor's own day-0 baseline to the peak window. Pure."""
    if not day0 or not peak:
        return {"n_day0": len(day0), "n_peak": len(peak), "delta": float("nan")}
    return {"n_day0": len(day0), "n_peak": len(peak),
            "age_day0": float(np.mean(day0)), "age_peak": float(np.mean(peak)),
            "delta": float(np.mean(peak) - np.mean(day0))}


def summarise(deltas: dict) -> dict:
    """Aggregate per-donor deltas: mean, paired 95% CI (t, small n). Pure."""
    v = np.array([d for d in deltas.values() if np.isfinite(d)], float)
    n = len(v)
    if n < 2:
        return {"n": int(n), "mean": float(v[0]) if n else float("nan"), "ci95": [None, None]}
    m, sd = float(v.mean()), float(v.std(ddof=1))
    t = {2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571}.get(n - 1, 1.96)
    se = sd / np.sqrt(n)
    return {"n": int(n), "mean": m, "sd": sd, "ci95": [m - t * se, m + t * se],
            "n_negative": int((v < 0).sum())}


def load_samples(gill_dir: str) -> list[dict]:
    """Every Gill sample with donor, day, cell type and the clock's predicted age."""
    from run_multi_local import discover_gill  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GillReprogrammingSource

    expr, series = discover_gill(gill_dir)
    src = GillReprogrammingSource(expr_tsv=expr, series_matrix=series)
    src._load()
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    cols = list(src._rpm.columns)
    ages = clock.predict_age(
        normalize_counts(src._rpm[cols].to_numpy(dtype=np.float64).T), list(src._genes))
    out = []
    for col, age in zip(cols, ages, strict=True):
        m = src._meta.get(col)
        if m:
            out.append({"title": col, "donor": m["donor"], "day": float(m["day"]),
                        "ctype": str(m.get("ctype", "")).strip(), "age": float(age)})
    return out


def main() -> int:
    gill = sys.argv[1] if len(sys.argv) > 1 else r"D:\Gill"
    print("STAGE 1.5.1 — can OUR clock reproduce Gill's ~30 yr rejuvenation? (read-only)\n")
    if not Path(gill).exists():
        print(f"  !! {gill} not found")
        return 1

    s = load_samples(gill)
    from collections import Counter
    print(f"  {len(s)} samples | cell types: {dict(Counter(x['ctype'] for x in s))}")
    print(f"  peak window = days {PEAK_DAYS[0]:.0f}-{PEAK_DAYS[1]:.0f} "
          f"(Gill's reported optimum); baseline = each donor's own day 0\n")

    results = {}
    for label, keep in (("REPROGRAMMING (responders)", REPROG),
                        ("FAILING to reprogram (negative control)", FAILING),
                        ("POOLED, as E1/E1b/D2 did", None)):
        per_donor: dict[str, float] = {}
        detail = {}
        for donor in sorted({x["donor"] for x in s}):
            # baseline: the donor's day-0 dermal fibroblast, regardless of arm
            d0 = [x["age"] for x in s if x["donor"] == donor and x["day"] == 0.0]
            pk = [x["age"] for x in s if x["donor"] == donor
                  and PEAK_DAYS[0] <= x["day"] <= PEAK_DAYS[1]
                  and x["ctype"] != "iPSC"
                  and (keep is None or x["ctype"] == keep)]
            c = contrast(d0, pk)
            detail[donor] = c
            if np.isfinite(c["delta"]):
                per_donor[donor] = c["delta"]
        agg = summarise(per_donor)
        results[label] = {"per_donor": detail, "aggregate": agg}
        print(f"  {label}")
        for d, c in detail.items():
            if np.isfinite(c["delta"]):
                print(f"     {d}: day0 {c['age_day0']:6.1f} -> peak {c['age_peak']:6.1f}"
                      f"   delta {c['delta']:+7.1f} yr   (n_peak={c['n_peak']})")
        if agg.get("ci95") and agg["ci95"][0] is not None:
            print(f"     MEAN {agg['mean']:+.1f} yr  95% CI [{agg['ci95'][0]:+.1f}, "
                  f"{agg['ci95'][1]:+.1f}]  ({agg.get('n_negative')}/{agg['n']} negative)")
        print()

    # ---- the proper statistic: PAIRED responder-minus-nonresponder, per donor ----
    print("  PAIRED per-donor: responder MINUS non-responder (cancels the identity artefact)")
    rd = results["REPROGRAMMING (responders)"]["per_donor"]
    fd = results["FAILING to reprogram (negative control)"]["per_donor"]
    paired = {d: rd[d]["delta"] - fd[d]["delta"] for d in rd
              if np.isfinite(rd[d].get("delta", np.nan)) and np.isfinite(fd[d].get("delta", np.nan))}
    for d, v in paired.items():
        print(f"     {d}: {v:+7.1f} yr")
    pag = summarise(paired)
    print(f"     MEAN {pag['mean']:+.1f} yr  95% CI [{pag['ci95'][0]:+.1f}, {pag['ci95'][1]:+.1f}]"
          f"   ({pag['n_negative']}/{pag['n']} negative)")
    sig = pag["ci95"][1] < 0
    print(f"     -> CI {'EXCLUDES' if sig else 'includes'} zero"
          f"{' => a significant rejuvenation signal at n=6' if sig else ''}\n")
    results["PAIRED responder-minus-nonresponder"] = {"per_donor": paired, "aggregate": pag}

    # ---- time course: does the separation peak where Gill says it does? ----
    print("  TIME COURSE of the responder-minus-nonresponder gap (Gill: optimum day 10-13)")
    tc = {}
    days = sorted({x["day"] for x in s if x["day"] > 0 and x["ctype"] in (REPROG, FAILING)})
    for day in days:
        per = {}
        for donor in sorted({x["donor"] for x in s}):
            a = [x["age"] for x in s if x["donor"] == donor and x["day"] == day
                 and x["ctype"] == REPROG]
            b = [x["age"] for x in s if x["donor"] == donor and x["day"] == day
                 and x["ctype"] == FAILING]
            if a and b:
                per[donor] = float(np.mean(a) - np.mean(b))
        if len(per) >= 3:
            g = summarise(per)
            tc[day] = g
            star = "  <-- Gill's optimum" if PEAK_DAYS[0] <= day <= PEAK_DAYS[1] else ""
            print(f"     day {day:5.0f}: gap {g['mean']:+7.1f} yr  "
                  f"95% CI [{g['ci95'][0]:+6.1f}, {g['ci95'][1]:+6.1f}]  n={g['n']}{star}")
    results["time_course_gap"] = tc
    print()

    r = results["REPROGRAMMING (responders)"]["aggregate"].get("mean", float("nan"))
    f = results["FAILING to reprogram (negative control)"]["aggregate"].get("mean", float("nan"))
    p = results["POOLED, as E1/E1b/D2 did"]["aggregate"].get("mean", float("nan"))
    print(f"  Gill reported ~{GILL_REPORTED:.0f} yr at this window.")
    print(f"  ours: responders {r:+.1f} | non-responders {f:+.1f} | pooled {p:+.1f}")
    if np.isfinite(r) and np.isfinite(f):
        print(f"  responder-minus-control separation = {r - f:+.1f} yr "
              "(the quantity a working clock should make negative)")

    (_RESULTS / "diag_gill_replication_results.json").write_text(json.dumps(
        {"script": "diag_gill_replication", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
         "peak_days": PEAK_DAYS, "gill_reported": GILL_REPORTED, "results": results},
        indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_gill_replication_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
