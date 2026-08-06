"""Register ARM D's bar — the STRATIFIED shuffle. Written BEFORE the run.

    python plan_tests/register_arm_d_bar.py

READ-ONLY. Writes `results/register_arm_d_bar_results.json`. `src/` untouched.

THE QUESTION
------------
Arm C permuted HFF's ΔAge labels **globally** and ranking collapsed
`rank_model_dage` 0.9476 -> 0.5765. That proved the labels carry *consistent exploitable structure*
but could not say what kind, because a global permutation destroys **two** things at once:

    (a) the BETWEEN-timepoint trajectory   -- rho(day, ΔAge) = -0.905, which the identity artefact
                                              would also produce, and
    (b) the WITHIN-timepoint cell-level pairing -- which only real per-cell signal would carry.

Stage 1.5.5 then removed the two mundane candidates for (b): within a timepoint, identity explains
2-16 % of ΔAge variance and technical covariates ~0-9 %, so 83-97 % is neither. Clock noise and real
signal both remain live.

**Arm D separates (a) from (b).** It permutes WITHIN each `(cell_line, time_h)` stratum, so the
between-stratum trajectory survives intact and only within-stratum pairing is destroyed.

    D ranks like A  ->  the exploitable structure is the BETWEEN-timepoint trajectory. Cell-level
                        pairing carries little. Consistent with a day-level systematic artefact.
    D ranks like C  ->  the structure is WITHIN-timepoint and cell-level. Consistent with real
                        per-cell signal; a day-level artefact cannot produce it.
    D between       ->  both contribute; report the split, claim neither pure account.
    D OUTSIDE [A,C] ->  registered below. Arm C's table did NOT do this and the result landed
                        540 % of the way from A to B, so no branch fitted.

TWO REGISTRATION LESSONS FROM ARM C, APPLIED HERE
--------------------------------------------------
1. **An outcome for landing outside the bracket.** Arm C's table offered only "like A" or "like B"
   and both assumed C in [A, B]. It was not, so the reading had to be labelled beyond-registration.
   §OUTSIDE below fixes that in advance.
2. **"No change from arm A" is an EQUIVALENCE claim.** A paired CI containing zero is absence of
   evidence. It needs a margin fixed before the run and a TOST -- and, as arm C showed, the margin's
   resolvability depends on an SD that is not known until the run, so that dependency is reported
   here rather than discovered afterwards.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

from audit_metrics import MIN_PASS_RATE, bar_verdict  # noqa: E402

N_FOLDS = 6
N_SIM = 40_000
SEED = 0

# Measured, per fold, in the step-6 runs. rank_model_dage.
ARM_A = [0.9104, 0.9091, 0.9896, 0.9701, 0.9596, 0.9468]
ARM_C = [0.4909, 0.5273, 0.8429, 0.5455, 0.4561, 0.5961]
# arm C - arm A, per fold: the full bracket this experiment must partition.
AC_DIFFS = [c - a for a, c in zip(ARM_C, ARM_A, strict=True)]


def bracket() -> tuple[float, float]:
    """The A->C ranking gap and its per-fold SD. Pure."""
    d = np.asarray(AC_DIFFS, dtype=float)
    return float(d.mean()), float(d.std(ddof=1))


def delta_eq() -> float:
    """Equivalence margin for 'D is like A'.

    Half the A->C bracket, so the margin partitions [A, C] at its midpoint -- the same rule used for
    arm C's Δ_eq (half the A-B gap), applied to this experiment's own measured bracket. Derived,
    not chosen, and NOT widened to whatever would pass.
    """
    return abs(bracket()[0]) / 2.0


def tost_equivalent(diffs: np.ndarray, margin: float, alpha: float = 0.05) -> np.ndarray:
    """Vectorised TOST: is the (1-2a) CI wholly inside +/- margin? Pure."""
    n = diffs.shape[1]
    mean = diffs.mean(axis=1)
    se = diffs.std(axis=1, ddof=1) / np.sqrt(n)
    half = float(stats.t.ppf(1 - alpha, n - 1)) * se
    return (mean - half > -margin) & (mean + half < margin)


def resolvability(margin: float, sd: float, true_effect: float = 0.0,
                  rng: np.random.Generator | None = None) -> float:
    """P(the bar declares equivalence) when the truth IS equivalence. Pure."""
    rng = rng or np.random.default_rng(SEED)
    return float(tost_equivalent(rng.normal(true_effect, sd, size=(N_SIM, N_FOLDS)), margin).mean())


def detectability(effect: float, sd: float, rng: np.random.Generator | None = None) -> float:
    """P(a paired 95 % CI excludes 0) for a true `effect`. The DIFFERENCE branch. Pure."""
    rng = rng or np.random.default_rng(SEED + 1)
    d = rng.normal(effect, sd, size=(N_SIM, N_FOLDS))
    m, s = d.mean(axis=1), d.std(axis=1, ddof=1)
    half = float(stats.t.ppf(0.975, N_FOLDS - 1)) * s / np.sqrt(N_FOLDS)
    return float((np.abs(m) > half).mean())


def main() -> int:
    rng = np.random.default_rng(SEED)
    gap, sd_ac = bracket()
    margin = delta_eq()
    mean_a, mean_c = float(np.mean(ARM_A)), float(np.mean(ARM_C))

    print("\nARM D — registering the STRATIFIED-shuffle bar, before the run\n")
    print(f"  arm A mean rank_model_dage : {mean_a:.4f}")
    print(f"  arm C mean rank_model_dage : {mean_c:.4f}")
    print(f"  bracket A -> C             : {gap:+.4f}  (per-fold SD {sd_ac:.4f}, n={N_FOLDS})")
    print(f"  DELTA_EQ = |bracket| / 2   : {margin:.4f}\n")

    # ---- the DIFFERENCE branch: can we detect D sitting at C? ------------------------- #
    det = detectability(gap, sd_ac, rng=rng)
    print(f"  'D is like C' is a DIFFERENCE test. P(detect a shift of the full bracket) = {det:.1%}")
    print(f"     -> {'WELL POWERED' if det >= MIN_PASS_RATE else 'UNDERPOWERED'}"
          f" at the A->C SD of {sd_ac:.4f}\n")

    # ---- the EQUIVALENCE branch: can we ever conclude 'D is like A'? ------------------ #
    print("  'D is like A' is an EQUIVALENCE claim -- TOST, 90 % CI inside +/-DELTA_EQ.")
    print("  Its resolvability depends on SD(D-A), unknown until the run, so it is swept:")
    print(f"     {'SD(D-A)':>9}{'P(declare equivalence | D truly = A)':>40}   verdict")
    rows = []
    for sd in (0.01, 0.02, 0.03, 0.05, 0.08, sd_ac):
        rate = resolvability(margin, sd, rng=rng)
        v = bar_verdict(np.array([rate]), MIN_PASS_RATE, lower_is_better=False)
        ok = rate >= MIN_PASS_RATE
        tag = "  <- the A->C SD (pessimistic)" if abs(sd - sd_ac) < 1e-9 else ""
        rows.append({"sd": sd, "pass_rate": rate, "resolvable": ok, "verdict": v["verdict"]})
        print(f"     {sd:>9.4f}{rate:>40.1%}   {'RESOLVABLE' if ok else 'UNRESOLVABLE'}{tag}")
    # The LARGEST SD that still resolves -- not the first swept one. Reporting the smallest
    # passing gridpoint is the exact bug caught in `register_gc_step2_bar.py`, where
    # `max(passing gridpoint)` understated the usable SD by 2x. Solved by bisection so the figure
    # is not an artefact of which SDs happen to be in the sweep.
    lo, hi = 0.001, 1.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if resolvability(margin, mid, rng=rng) >= MIN_PASS_RATE:
            lo = mid
        else:
            hi = mid
    sd_needed = lo
    print(f"\n     -> the equivalence branch resolves for ANY SD(D-A) up to ~{sd_needed:.3f}"
          f" (solved, not read off the sweep).")
    print("        The achieved SD will be reported alongside the verdict either way.")

    # ---- false-equivalence: a D that really sits at C must not read as 'like A' ------- #
    fp = resolvability(margin, sd_ac, true_effect=gap, rng=rng)
    print(f"\n  false-equivalence check: a D truly AT arm C is declared 'like A' {fp:.1%} "
          "of the time (must be ~0)")

    # ---- the registered outcome table, INCLUDING outside the bracket ------------------ #
    print("\n  PRE-REGISTERED OUTCOMES (position = (A - D) / (A - C), 0 % = at A, 100 % = at C):")
    table = [
        ("D equivalent to A by TOST",
         "the exploitable structure is the BETWEEN-timepoint trajectory; cell-level pairing "
         "carries little. Consistent with a day-level systematic artefact."),
        ("D - C not detectable AND D - A detectable",
         "the structure is WITHIN-timepoint and cell-level. A day-level artefact cannot produce "
         "it; real per-cell signal remains the live explanation."),
        ("both differences detectable (D strictly between)",
         "both components contribute. Report the split as a proportion; claim neither pure "
         "account. This is NOT a null."),
        ("neither difference detectable",
         "INCONCLUSIVE. Underpowered, licenses nothing -- the same rule that governed step 6."),
        ("D OUTSIDE [A, C] -- better than A, or worse than C",
         "REGISTERED BECAUSE ARM C LANDED OUTSIDE ITS OWN TABLE. Not a graded branch: report the "
         "position, state that no pre-registered reading applies, and treat any interpretation as "
         "beyond-registration. A D worse than C would mean stratification destroys MORE than the "
         "global shuffle, which no current account predicts and would need its own investigation."),
    ]
    for i, (cond, reading) in enumerate(table, 1):
        print(f"   {i}. IF {cond}")
        print(f"      -> {reading}")

    out = {"script": "register_arm_d_bar",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "arm_a_mean": mean_a, "arm_c_mean": mean_c,
           "bracket_a_to_c": gap, "bracket_sd": sd_ac,
           "delta_eq": margin, "delta_eq_derivation": "half the measured A->C rank bracket",
           "test_equivalence": "TOST, 90 % CI inside +/- delta_eq",
           "test_difference": "paired 95 % CI excludes 0",
           "difference_branch_power_at_ac_sd": det,
           "equivalence_resolvability_by_sd": rows,
           "sd_at_which_equivalence_resolves": sd_needed,
           "false_equivalence_rate_if_D_is_really_C": fp,
           "outcomes": [{"condition": c, "reading": r} for c, r in table],
           "registers_outside_bracket_outcome": True}
    (_RESULTS / "register_arm_d_bar_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote register_arm_d_bar_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
