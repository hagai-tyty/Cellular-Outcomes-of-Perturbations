"""Register the bar for STAGE 1.5.3 step 6 / Stage 1.5.2 G-c step 2 -- the arm comparison.

    python plan_tests/register_gc_step2_bar.py

READ-ONLY. Writes `results/register_gc_step2_bar_results.json`. `src/` untouched.

WHY THIS EXISTS
---------------
Step 6 decides whether HFF's 33,613 age labels -- **99.7 % of every age label in the project** --
are kept or discarded. It compares two arms on `dage_mae_model`, paired across 6 donor folds, and
the outcome table reads the result as "A better, CI excludes 0 => keep them".

**That criterion had no registered bar.** Every bar registered for this stage (B1/B2 for C-5,
A1/A2/A3 for C-5c) grades a MECHANISM; none grades the comparison the stage actually turns on.
Ground rule sec 5b is explicit: *"a bar with no such test is not considered pre-registered."*

That gap is the exact failure that bit Stage 1 twice on `fate_ece` -- a bar honestly set in advance
that a correct system still fails, because the estimator is noisy at the geometry it is graded on.
At 6 folds the paired CI's minimum detectable effect is `t(.975,5)/sqrt(6) = 1.049` standard
deviations of the PER-FOLD DIFFERENCE, and that SD has never been measured. Baseline
`dage_mae_model` already ranges 5.39 -> 29.69 across the six folds.

THE ASYMMETRY THAT SHAPES THIS BAR
----------------------------------
The two directions are not equally dangerous.

  * **"B better"** (masking helps) is self-limiting: it says the labels were harmful, and the
    action -- drop them -- is what G-c step 1's artefact finding already suggests.
  * **"CI includes 0"** is the trap. Read as *"HFF's labels contribute nothing, so discard them"*,
    it licenses throwing away 99.7 % of the age labels **on a null that may simply be underpowered.**
    An underpowered null is not evidence of absence.

So the protection this bar provides is mostly about how a NULL may be read, and that is
pre-registered below rather than left to whoever reads the output.

DELTA_STAR -- DERIVED, NOT PICKED
---------------------------------
The smallest effect worth acting on is taken from the project's OWN existing bar for a meaningful
change in this same metric: `STAGE_2_LEVEL_CORRECTION.md` sec 12 registers a **>= 25 % drop in
`dage_mae_model`** as its TARGET. Against the recorded baseline mean of 14.29 yr that is
**3.57 yr**. Using the established threshold for the same metric avoids inventing a number here.
"""
from __future__ import annotations

import json
import math
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# Spelled out rather than `ROOT / "results"`: tests/test_results_paths.py checks this form by
# regex and cannot follow an indirection through another name.
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

from audit_metrics import MIN_PASS_RATE, sensitivity_multiplier  # noqa: E402

N_FOLDS = 6
N_SIM = 40_000
SEED = 0

# The six recorded baseline folds. Used only to characterise how heterogeneous this metric already
# is -- NOT as the difference's SD, which is a different (paired) quantity and is unmeasured.
BASELINE_FOLDS = {"N2": 21.7936, "N3": 29.6950, "O1": 5.3876,
                  "O2": 7.5350, "Y1": 7.2791, "Y2": 14.0567}

# Candidate per-fold SDs of the PAIRED DIFFERENCE. 0.5-2 is the regime where the arms move together
# (most fold variance cancels in the pairing); 13.7 is the pessimistic bound where they are
# effectively independent, sqrt(2) * SD(baseline folds).
CANDIDATE_SDS = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 13.7)


def delta_star() -> float:
    """Smallest effect worth acting on: Stage 2 sec 12's registered >=25 % drop, applied to the
    recorded baseline mean. Derived from an existing bar rather than chosen for this test."""
    return 0.25 * float(np.mean(list(BASELINE_FOLDS.values())))


def paired_ci_excludes_zero(diffs: np.ndarray) -> np.ndarray:
    """Vectorised: for each simulated run (a row of `N_FOLDS` paired differences), does the
    two-sided 95 % t CI on the mean exclude zero? This is exactly the rule the outcome table uses.
    """
    from scipy.stats import t as student_t

    n = diffs.shape[1]
    mean = diffs.mean(axis=1)
    se = diffs.std(axis=1, ddof=1) / math.sqrt(n)
    tcrit = float(student_t.ppf(0.975, n - 1))
    return np.abs(mean) > tcrit * se


def power_at(delta: float, sd: float, rng: np.random.Generator) -> dict:
    """Fraction of runs where the CI excludes 0 AND points the right way, when the TRUE effect is
    `delta`. With delta = 0 this is the false-positive rate and must sit near 0.05."""
    d = rng.normal(delta, sd, size=(N_SIM, N_FOLDS))
    excl = paired_ci_excludes_zero(d)
    correct = excl & (np.sign(d.mean(axis=1)) == np.sign(delta)) if delta else excl
    return {"delta": delta, "sd": sd,
            "pass_rate": float(correct.mean()),
            "mde": sensitivity_multiplier(N_FOLDS) * sd}


def max_resolvable_sd(delta: float, rng: np.random.Generator, *,
                      lo: float = 0.1, hi: float = 8.0, iters: int = 24) -> float:
    """Largest SD(per-fold difference) at which DELTA* still clears MIN_PASS_RATE. Pure-ish.

    Solved by bisection rather than read off `CANDIDATE_SDS`. The grid is for display and jumps
    1.0 -> 2.0, so its largest passing gridpoint is 1.0 while the true crossover is near 1.9 --
    taking the gridpoint understates the usable SD by about 2x and would wrongly label a
    perfectly well-powered run INCONCLUSIVE. Power falls monotonically in SD at fixed delta,
    which is what makes bisection valid here.
    """
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if power_at(delta, mid, rng)["pass_rate"] >= MIN_PASS_RATE:
            lo = mid
        else:
            hi = mid
    return lo


def main() -> int:
    rng = np.random.default_rng(SEED)
    ds = delta_star()
    k = sensitivity_multiplier(N_FOLDS)
    base_sd = float(np.std(list(BASELINE_FOLDS.values()), ddof=1))

    print("\nSTEP 6 / G-c step 2 -- registering the arm-comparison bar")
    print(f"  geometry            : {N_FOLDS} donor folds, PAIRED")
    print(f"  MDE multiplier      : t(.975,{N_FOLDS-1})/sqrt({N_FOLDS}) = {k:.4f}")
    print(f"  baseline fold spread: SD {base_sd:.2f} yr (range "
          f"{min(BASELINE_FOLDS.values()):.2f}-{max(BASELINE_FOLDS.values()):.2f})")
    print(f"  DELTA* (Stage 2 sec12's 25% of the {np.mean(list(BASELINE_FOLDS.values())):.2f} yr "
          f"baseline) = {ds:.2f} yr\n")

    rows = []
    print(f"  {'SD(diff)':>9} {'MDE':>7} {'P(detect DELTA*)':>17}  verdict")
    for sd in CANDIDATE_SDS:
        r = power_at(ds, sd, rng)
        r["verdict"] = "RESOLVABLE" if r["pass_rate"] >= MIN_PASS_RATE else "UNRESOLVABLE"
        rows.append(r)
        print(f"  {sd:9.1f} {r['mde']:7.2f} {r['pass_rate']:17.4f}  {r['verdict']}")

    fp = power_at(0.0, 3.0, rng)
    resolvable = [r["sd"] for r in rows if r["verdict"] == "RESOLVABLE"]
    grid_sd = max(resolvable) if resolvable else 0.0
    max_sd = max_resolvable_sd(ds, rng)

    print(f"\n  false-positive check (true delta = 0): {fp['pass_rate']:.4f}  "
          f"(must be ~0.05; a CI that over-rejects would manufacture a verdict)")
    print(f"\n  => DELTA* is detectable at >= {MIN_PASS_RATE:.0%} whenever "
          f"SD(per-fold difference) <= {max_sd:.2f} yr.")
    print(f"     (The grid above only shows {grid_sd:.1f}, its largest PASSING gridpoint -- it "
          f"jumps {grid_sd:.1f} -> 2.0 and\n      never samples between, so reading the crossover "
          "off it understates the usable SD by ~2x.\n      The figure above is solved for, not "
          "read off the grid. Corrected 2026-08-02.)")
    print("     The run MUST report its own observed SD and MDE. If |observed effect| <= MDE the")
    print("     result is INCONCLUSIVE -- NOT 'the labels make no difference'.\n")

    out = {
        "script": "register_gc_step2_bar",
        "utc": datetime.now(UTC).isoformat(),
        "geometry": {"n_folds": N_FOLDS, "paired": True, "n_sim": N_SIM, "seed": SEED},
        "mde_multiplier": k,
        "baseline_folds": BASELINE_FOLDS,
        "baseline_fold_sd": base_sd,
        "delta_star_years": ds,
        "max_resolvable_sd_years": max_sd,
        "max_resolvable_sd_gridpoint": grid_sd,
        "delta_star_derivation":
            "Stage 2 sec 12 registers >=25% drop in dage_mae_model as its TARGET; applied to the "
            "recorded baseline mean of 14.29 yr. Derived from an existing bar, not chosen here.",
        "checks": {f"SD={r['sd']}": r for r in rows},
        "false_positive_rate_at_delta_0": fp["pass_rate"],
        "max_resolvable_sd": max_sd,
        "decision_rule": {
            "conclusive_A_better": "mean(A-B) < 0 favouring A AND CI excludes 0 AND |effect| > MDE",
            "conclusive_B_better": "mean(A-B) > 0 favouring B AND CI excludes 0 AND |effect| > MDE",
            "INCONCLUSIVE":
                "CI includes 0 while MDE > DELTA*. The comparison could not have detected an "
                "effect worth acting on, so it licenses NOTHING -- in particular it does NOT "
                "license discarding HFF's labels. Report as underpowered and say so.",
            "null_is_only_evidence_of_absence_when": "MDE <= DELTA* (i.e. observed SD <= "
                f"{max_sd:.1f} yr at {N_FOLDS} folds)",
        },
    }
    (_RESULTS / "register_gc_step2_bar_results.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("  wrote results/register_gc_step2_bar_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
