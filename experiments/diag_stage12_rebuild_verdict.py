"""STAGE 12 §12.9 — apply the pre-registered reading to the rebuild, mechanically.

`scorecard.py compare` gives the paired CIs, but §12.9's target is **coverage moving TOWARD
nominal**, and that is not what the scorecard's `conformal_coverage` row measures.

    METRICS["conformal_coverage"] = ("higher", ...)

Higher is better only UP TO the nominal level. Coverage 1.000 is not a triumph -- it means the
intervals are too wide, which is why `conformal_width` sits at 63-81 years. A coverage of 1.00 and
a coverage of 0.80 are both 0.10 away from a nominal 0.90, and a rule that prefers the former is
measuring width, not calibration.

So the target statistic here is **|coverage - nominal| per fold**, paired across folds -- exactly
the shape Stage 13 established for `level_shift`, applied to a metric the scorecard still judges
directionally. This script computes it, alongside the guards, and prints the pre-registered
verdict without a human choosing which branch to read.

Read-only: two committed snapshots in, one results file out. Runs in milliseconds.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage12_rebuild_verdict_results.json"

BASELINE_TAG = "c7_A_keep_hff"     # the _c7t folds
TREATMENT_TAG = "c7t_stage12"      # the _s12 folds
NOMINAL = 0.90                     # conformal_level, as recorded in every fold of the snapshots

# §12.9 guards. `rank_model_dage` may be noise or better, never a regression.
GUARDS = ("fate_prauc", "fate_roc", "rank_model_dage")


def _sc():
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def coverage_gap(folds: dict, donors, nominal: float = NOMINAL) -> dict:
    """|coverage - nominal| per fold. The quantity §12.9 actually cares about."""
    out = {}
    for d in donors:
        f = folds.get(d)
        if isinstance(f, dict) and "_error" not in f and f.get("conformal_coverage") is not None:
            out[d] = abs(float(f["conformal_coverage"]) - nominal)
    return out


def paired(a: dict, b: dict) -> tuple:
    """Paired mean difference (b - a) and 95% CI over the folds present in both."""
    common = sorted(set(a) & set(b))
    diffs = [b[d] - a[d] for d in common]
    n = len(diffs)
    if n < 2:
        return None, (None, None), n, common
    md = float(np.mean(diffs))
    sd = float(np.std(diffs, ddof=1))
    se = sd / np.sqrt(n)
    t = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n, common


def verdict_from(md, lo, hi) -> str:
    """§12.9's three pre-registered outcomes, applied to the coverage GAP (lower = closer)."""
    if md is None:
        return "UNDETERMINED (too few paired folds)"
    if hi < 0:
        return ("COVERAGE MOVED TOWARD NOMINAL -- the split defect was materially degrading "
                "calibration; Stage 12 is a real improvement, not only a correctness tidy-up")
    if lo > 0:
        return ("COVERAGE MOVED AWAY FROM NOMINAL -- INVESTIGATE BEFORE ACCEPTING. A more "
                "representative calib set cannot make calibration worse by itself")
    return ("NO DETECTABLE MOVE -- the composition shift was too small to matter at n=6. "
            "Stage 12 remains correct-but-inert, and that is a publishable negative")


def run(baseline: str = BASELINE_TAG, treatment: str = TREATMENT_TAG) -> dict:
    sc = _sc()
    pa = ROOT / "scorecard" / f"{baseline}.json"
    pb = ROOT / "scorecard" / f"{treatment}.json"
    for p in (pa, pb):
        if not p.exists():
            return {"ABORTED": True, "reason": f"missing snapshot: {p.name}"}
    A = json.loads(pa.read_text(encoding="utf-8"))["folds"]
    B = json.loads(pb.read_text(encoding="utf-8"))["folds"]

    ga, gb = coverage_gap(A, sc.DONORS), coverage_gap(B, sc.DONORS)
    md, (lo, hi), n, common = paired(ga, gb)

    guards = {}
    for key in GUARDS:
        direction = sc.METRICS[key][0]
        gm, (glo, ghi), gn = sc._paired(A, B, key, magnitude=(direction == "abs"))
        guards[key] = {"mean_diff": gm, "ci": [glo, ghi], "n": gn,
                       "verdict": sc._verdict(direction, gm, glo, ghi)}

    raw_cov = {d: {"baseline": A[d].get("conformal_coverage") if d in A else None,
                   "treatment": B[d].get("conformal_coverage") if d in B else None}
               for d in common}
    return {
        "ABORTED": False, "baseline": baseline, "treatment": treatment, "nominal": NOMINAL,
        "coverage_gap_baseline": ga, "coverage_gap_treatment": gb,
        "paired_folds": common, "n_folds": n,
        "mean_diff": md, "ci": [lo, hi],
        "verdict": verdict_from(md, lo, hi),
        "raw_coverage": raw_cov,
        "guards": guards,
        "guard_breach": [k for k, v in guards.items() if v["verdict"] == "REGRESSION"],
    }


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                     # noqa: BLE001
            pass
    r = run()
    if r.get("ABORTED"):
        print(f"\n  ABORTED — {r['reason']}")
        return 1
    print(f"\nSTAGE 12 §12.9 — {r['baseline']}  ->  {r['treatment']}")
    print(f"  TARGET: |conformal_coverage - {r['nominal']}| per fold, paired. Lower is closer to")
    print("  nominal. The scorecard's own row judges coverage as 'higher is better', which is")
    print("  wrong for a target-seeking metric -- 1.000 means the intervals are too wide.")
    print(f"\n     {'fold':<6}{'cov base':>10}{'cov new':>10}{'gap base':>10}{'gap new':>10}"
          f"{'change':>10}")
    for d in r["paired_folds"]:
        c = r["raw_coverage"][d]
        ga, gb = r["coverage_gap_baseline"][d], r["coverage_gap_treatment"][d]
        print(f"     {d:<6}{c['baseline']:>10.3f}{c['treatment']:>10.3f}{ga:>10.3f}{gb:>10.3f}"
              f"{gb - ga:>+10.3f}")
    print(f"\n  paired mean change in gap: {r['mean_diff']:+.4f}  "
          f"95% CI [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}]  (n={r['n_folds']} folds)")
    print(f"\n  VERDICT: {r['verdict']}")

    print("\n  GUARDS (§12.9)")
    for k, v in r["guards"].items():
        ci = "" if v["mean_diff"] is None else f"[{v['ci'][0]:+.3f},{v['ci'][1]:+.3f}]"
        md = "n/a" if v["mean_diff"] is None else f"{v['mean_diff']:+.3f}"
        print(f"     {k:<20}{md:>9}  {ci:>18}  n={v['n']}  {v['verdict']}")
    print(f"\n  GUARD BREACHES: {r['guard_breach'] or 'none'}")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
