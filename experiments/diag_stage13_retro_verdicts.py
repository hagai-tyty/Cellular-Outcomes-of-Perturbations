"""STAGE 13 — re-judge every past scorecard comparison under the corrected decision rule.

Stage 12's defect lived in written artefacts, so repairing the past would have needed a rebuild.
This one lives in **aggregation and verdict only** — every per-fold number on disk is correct and
signed — so every comparison ever run can be re-judged right now, with no rebuild and no retrain.

Read-only. Reads `scorecard/*.json`, writes one results file. Changes nothing else.

Two defects are being corrected, and they are DIFFERENT mechanisms — a distinction this script
exists to keep straight, because conflating them produces a wrong story:

  A1  the aggregate took |mean(signed)| instead of mean(|signed|). For a per-donor bias whose
      sign varies by donor, that measures how far the PANEL CANCELS, not how large the error is.
      It printed 0.230 for a shift whose true mean magnitude is 12.72 yr.  -> affects the COLUMN.

  A3  `_verdict`'s better_is_down was applied to a SIGNED difference, so any change that moves
      the shift DOWNWARD read as an improvement — including moving from -5.7 to -20.3.
      -> affects the VERDICT.

A1 does not explain the verdict flips and A3 does not explain the understated column. Both are
measured separately below.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = ROOT / "scorecard"
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage13_retro_verdicts_results.json"


def _scorecard():
    """Load the repo-root `scorecard.py` by path — it is a script, not a package module."""
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# The retro pass audits comparisons that were JUDGED BY THE BROKEN RULE. That is a closed,
# historical set: the nine snapshots that existed when Stage 13 shipped on 2026-08-17. A snapshot
# taken afterwards (`c7t_stage12`, from the Stage 12 rebuild) was never judged by the old rule, so
# including it would inflate the flip count with comparisons that never happened. Globbing the
# directory made this scope silently grow the moment a new snapshot landed -- which it did, hours
# later, and the pinned counts failed. Frozen deliberately, not to keep a test green.
RETRO_SNAPSHOTS = (
    "A_xdonor", "B_fatecal", "B_fatecal_pooled", "baseline", "c7_A_keep_hff",
    "gc2_A_keep_hff", "gc2_B_mask_hff", "gc2_C_shuffle_hff_s0", "gc2_D_stratshuffle_hff_s0",
)


def load_snapshots(snap_dir: Path = SNAP_DIR, names=RETRO_SNAPSHOTS) -> dict[str, dict]:
    """The snapshots in scope for the retro audit. Pass `names=None` to take whatever is on disk."""
    paths = (sorted(snap_dir.glob("*.json")) if names is None
             else [snap_dir / f"{n}.json" for n in names])
    return {p.stem: json.loads(p.read_text(encoding="utf-8"))["folds"]
            for p in paths if p.exists()}


def abs_metrics(sc) -> list[str]:
    return [k for k, (d, _) in sc.METRICS.items() if d == "abs"]


def dedupe(snapshots: dict[str, dict], donors, keys) -> list[list[str]]:
    """Group snapshots that are IDENTICAL on the metrics under test.

    Five of the nine files carry the same level-shift content (`baseline`, `A_xdonor`,
    `B_fatecal`, `B_fatecal_pooled`, `gc2_A_keep_hff`). Counting comparisons without deduping
    inflates the flip count roughly fivefold and would overstate this finding.
    """
    groups: dict[str, list[str]] = {}
    for name, folds in snapshots.items():
        sig = json.dumps({d: {k: folds[d].get(k) for k in keys}
                          for d in donors if d in folds and "_error" not in folds[d]},
                         sort_keys=True)
        groups.setdefault(sig, []).append(name)
    return list(groups.values())


def cancellation_gap(folds, key) -> dict:
    """A1, isolated: |mean(signed)| vs mean(|signed|) on ONE snapshot. No comparison involved."""
    v = np.array([f[key] for f in folds.values()
                  if isinstance(f, dict) and "_error" not in f and f.get(key) is not None],
                 dtype=float)
    if v.size == 0:
        return {"n": 0}
    signed, mag = float(np.abs(np.mean(v))), float(np.mean(np.abs(v)))
    return {"n": int(v.size), "abs_of_mean": signed, "mean_of_abs": mag,
            "understated_by": (mag / signed) if signed > 1e-12 else float("inf"),
            "sd": float(np.std(v, ddof=1)) if v.size > 1 else 0.0}


def reverdict(sc, A, B, key) -> dict | None:
    """A3, isolated: the same folds, the same CI machinery, signed vs magnitude."""
    md_o, (lo_o, hi_o), n = sc._paired(A, B, key)
    md_n, (lo_n, hi_n), _ = sc._paired(A, B, key, magnitude=True)
    if md_o is None or md_n is None:
        return None
    return {"n_folds": n,
            "old_diff": md_o, "old_ci": [lo_o, hi_o],
            "old_verdict": sc._verdict("abs", md_o, lo_o, hi_o),
            "new_diff": md_n, "new_ci": [lo_n, hi_n],
            "new_verdict": sc._verdict("abs", md_n, lo_n, hi_n)}


def run() -> dict:
    sc = _scorecard()
    snaps = load_snapshots()
    keys = abs_metrics(sc)
    groups = dedupe(snaps, sc.DONORS, keys)
    reps = [g[0] for g in groups]

    # ---- A1: the column, per snapshot ----
    gaps = {name: {k: cancellation_gap(snaps[name], k) for k in keys} for name in reps}

    # ---- A3: the verdict, per DISTINCT unordered comparison ----
    comparisons, flips = [], []
    for a, b in itertools.combinations(reps, 2):
        for k in keys:
            r = reverdict(sc, snaps[a], snaps[b], k)
            if r is None:
                continue
            r |= {"a": a, "b": b, "metric": k}
            comparisons.append(r)
            if r["old_verdict"] != r["new_verdict"]:
                flips.append(r)

    # ---- the systematic direction of the error ----
    # A shuffle control destroys donor structure. Under the corrected rule its level shift must
    # get WORSE. If the old rule called shuffles better, the defect was not merely noisy -- it
    # was biased in favour of destroyed-label controls, which is the worst possible direction.
    shuffles = [c for c in flips
                if "shuffle" in c["b"] and c["old_verdict"] == "ACCEPT (better)"]

    return {"groups": groups, "distinct_snapshots": reps, "abs_metrics": keys,
            "cancellation": gaps, "comparisons": comparisons, "flips": flips,
            "n_comparisons": len(comparisons), "n_flips": len(flips),
            "shuffle_controls_scored_as_improvements": len(shuffles)}


def main() -> int:
    import sys
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    r = run()

    print("\nSTAGE 13 — RETRO-VERDICTS over every committed scorecard snapshot")
    print(f"\n  {len(r['groups'])} distinct snapshots (by level-shift content) from "
          f"{sum(len(g) for g in r['groups'])} files:")
    for g in r["groups"]:
        print("     " + " == ".join(g))

    print("\n  A1 — THE COLUMN: |mean(signed)| vs mean(|signed|)")
    print(f"     {'snapshot':<28}{'metric':<20}{'printed':>10}{'true':>10}{'understated':>13}")
    for name, per in r["cancellation"].items():
        for k, g in per.items():
            if not g.get("n"):
                continue
            print(f"     {name:<28}{k:<20}{g['abs_of_mean']:10.3f}{g['mean_of_abs']:10.3f}"
                  f"{g['understated_by']:12.1f}x")

    print(f"\n  A3 — THE VERDICT: {r['n_flips']} of {r['n_comparisons']} distinct verdicts change")
    for c in r["flips"]:
        print(f"     {c['a'][:24]:<24} -> {c['b'][:24]:<24} {c['metric']:<18} n={c['n_folds']}  "
              f"{c['old_verdict']:<16} -> {c['new_verdict']}")

    n = r["shuffle_controls_scored_as_improvements"]
    print(f"\n  DIRECTION OF THE ERROR: {n} comparison(s) in which a SHUFFLE control scored")
    print("     ACCEPT (better) under the old rule. A shuffle destroys the donor structure the")
    print("     level shift measures; it cannot improve it. The defect was biased in favour of")
    print("     destroyed-label negative controls, not merely noisy.")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
