"""STAGE 1.5.2 §17 — the re-audit: is the meth<->meth CEILING the same in every arm?

    python experiments/diag_m2a_per_arm_ceiling.py

READ-ONLY. Reads `diag_m2a_calibratability_results.json` (no raw data, no re-measurement) and
writes `diag_m2a_per_arm_ceiling_results.json`. `src/` untouched.

WHY THIS EXISTS
---------------
§12-R measured one pooled number: the two Horvath clocks agree with each other at rho_partial
**+0.568** over all 68 conditions. §11 separately reported RNA<->methylation **per arm**, and read
the pattern as:

    "The RNA clock tracks methylation age in cells that are NOT reprogramming, and stops tracking --
     or inverts -- in exactly the cells that are."

**Those two results were never put side by side, and they must be.** A low RNA<->meth correlation
inside an arm means one of two very different things:

  * the methylation reference is SHARP there and the RNA clock disagrees with it  -> the RNA clock
    is failing, and that is the strongest evidence the stage can produce;
  * the methylation reference is BLUNT there (the two clocks do not even agree with each other)
    -> nothing can be concluded about the RNA clock in that arm at all.

§11 did not distinguish them. This does, by computing the **per-arm** meth<->meth agreement as the
denominator for the per-arm RNA<->meth figure already reported.

WHAT IS AND IS NOT COMPARABLE
-----------------------------
These per-arm correlations are **unpartialled** Spearman within an arm, whereas §12-R's +0.568 is
partialled for pluripotency across all arms. They are not the same statistic and are not compared
numerically to it. Conditioning on arm already removes most of the between-arm reprogramming axis
that the partialling existed to handle, which is why the within-arm form is the right one here.

Every arm's n is 9-12, which the §6 freeze established is **UNRESOLVABLE** for a rho bar. So nothing
here is a criterion, and no verdict moves. It is a **reading correction** to §11, which is exactly
what a re-audit is for.
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
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

# The arms in which cells are actively reprogramming. Taken from the arm vocabulary, not inferred
# from any measurement.
REPROGRAMMING_ARMS = {"transiently_reprogrammed", "transient_reprogramming_intermediate"}
MIN_N = 4                 # below this a Spearman is not computable in any useful sense
SHARP_CEILING = 0.70      # what M-2a's own null assumed a real agreement looks like


def per_arm_ceiling(rows_a: list[dict], rows_b: list[dict]) -> dict:
    """Per-arm meth<->meth agreement, and RNA<->meth as a fraction of it. Pure.

    `rows_a` / `rows_b` are the per-condition rows M-2a already wrote for the two clocks, keyed by
    (donor, arm, day). Nothing is re-measured; this is a different cut of the same numbers.
    """
    from scipy.stats import spearmanr
    A = {(r["donor"], r["arm"], r["day"]): r for r in rows_a}
    B = {(r["donor"], r["arm"], r["day"]): r for r in rows_b}
    keys = sorted(set(A) & set(B))
    out = {}
    for arm in sorted({k[1] for k in keys}):
        ks = [k for k in keys if k[1] == arm]
        if len(ks) < MIN_N:
            continue
        ma = [A[k]["age_meth"] for k in ks]
        mb = [B[k]["age_meth"] for k in ks]
        ra = [A[k]["age_rna"] for k in ks]
        ceil = float(spearmanr(ma, mb).correlation)
        r_a = float(spearmanr(ra, ma).correlation)
        r_b = float(spearmanr(ra, mb).correlation)
        out[arm] = {
            "n": len(ks), "meth_vs_meth": ceil, "rna_vs_clock_a": r_a, "rna_vs_clock_b": r_b,
            "rna_mean": (r_a + r_b) / 2.0,
            "frac_of_ceiling": ((r_a + r_b) / 2.0 / ceil) if ceil > 0 else float("nan"),
            "reprogramming": arm in REPROGRAMMING_ARMS,
            # An arm is only INTERPRETABLE about the RNA clock if the reference agrees with itself
            # there. Below the level M-2a's own null assumed, it does not.
            "interpretable": ceil >= SHARP_CEILING,
        }
    return out


def reading_correction(arms: dict) -> dict:
    """What §11's per-arm table may and may not be read as. Pure."""
    interp = {a: v for a, v in arms.items() if v["interpretable"]}
    if not interp:
        return {"status": "NO_ARM_IS_INTERPRETABLE",
                "reason": "in no arm do the two methylation clocks agree with each other well "
                          "enough to arbitrate the RNA clock. §11's per-arm table cannot be read "
                          "as evidence about the RNA clock in either direction."}
    sharpest = max(arms.items(), key=lambda kv: kv[1]["meth_vs_meth"])
    bluntest = min(arms.items(), key=lambda kv: kv[1]["meth_vs_meth"])
    fails_where_sharp = [a for a, v in interp.items() if v["rna_mean"] < 0.2]
    return {
        "status": "PARTIAL",
        "n_interpretable": len(interp), "n_arms": len(arms),
        "sharpest_arm": sharpest[0], "sharpest_ceiling": sharpest[1]["meth_vs_meth"],
        "sharpest_rna": sharpest[1]["rna_mean"],
        "bluntest_arm": bluntest[0], "bluntest_ceiling": bluntest[1]["meth_vs_meth"],
        "rna_fails_in_interpretable_arms": fails_where_sharp,
        "reason": (f"only {len(interp)} of {len(arms)} arms have a methylation reference sharp "
                   f"enough to arbitrate anything (meth<->meth >= {SHARP_CEILING}). The rest are "
                   "uninterpretable about the RNA clock, in either direction."),
    }


def main() -> int:
    src = ROOT / "results" / "diag_m2a_calibratability_results.json"
    if not src.exists():
        print(f"[!] {src} not found — run diag_m2a_calibratability.py first.")
        return 1
    m = json.loads(src.read_text(encoding="utf-8"))
    names = list(m["clocks"])
    arms = per_arm_ceiling(m["clocks"][names[0]]["rows"], m["clocks"][names[1]]["rows"])

    print("STAGE 1.5.2 §17 — per-arm methylation CEILING (re-audit, read-only)\n")
    print("  A low RNA<->meth correlation means nothing until you know whether the methylation")
    print("  reference is sharp in that arm. §11 reported the numerator without the denominator.\n")
    print(f"  {'arm':<46}{'n':>4}{'meth<->meth':>13}{'RNA (mean)':>12}{'% of ceil':>11}  ")
    print("  " + "-" * 90)
    for a, v in sorted(arms.items(), key=lambda kv: -kv[1]["meth_vs_meth"]):
        tag = "REPROG" if v["reprogramming"] else "  -   "
        interp = "" if v["interpretable"] else "   <- reference too blunt to interpret"
        pct = f"{v['frac_of_ceiling']:.0%}" if np.isfinite(v["frac_of_ceiling"]) else "n/a"
        print(f"  {a:<46}{v['n']:>4}{v['meth_vs_meth']:>13.3f}{v['rna_mean']:>12.3f}"
              f"{pct:>11}  {tag}{interp}")

    dec = reading_correction(arms)
    print(f"\n  ==> {dec['status']}: {dec['reason']}")
    if dec["status"] == "PARTIAL":
        print(f"      sharpest arm: {dec['sharpest_arm']} "
              f"(meth<->meth {dec['sharpest_ceiling']:+.3f}), RNA there "
              f"{dec['sharpest_rna']:+.3f}")
        print(f"      bluntest arm: {dec['bluntest_arm']} "
              f"(meth<->meth {dec['bluntest_ceiling']:+.3f}) — uninterpretable")

    out = {"script": "diag_m2a_per_arm_ceiling",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "source": src.name, "sharp_ceiling": SHARP_CEILING,
           "arms": arms, "reading_correction": dec}
    (_RESULTS / "diag_m2a_per_arm_ceiling_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_m2a_per_arm_ceiling_results.json")
    print("  NOTE: every arm's n is 9-12, which §6 froze as UNRESOLVABLE for a rho bar.")
    print("  Nothing here is a criterion and no verdict moves. It is a READING correction to §11.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
