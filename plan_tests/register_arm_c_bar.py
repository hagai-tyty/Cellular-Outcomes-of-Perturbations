"""Register ARM C's equivalence bar — BEFORE the run. Stage 1.5.3, step 6 follow-up.

    python plan_tests/register_arm_c_bar.py

READ-ONLY. Writes `results/register_arm_c_bar_results.json`. `src/` untouched.

THE QUESTION ARM C ANSWERS
--------------------------
Step 6's rerun found a consistent ranking cost when HFF's labels are withheld:
`rank_model_dage` **-0.0688** (95% CI [-0.1004, -0.0371]), with ridge moving nearly the same amount.
That is real, but it is CONFOUNDED between two explanations:

    (i)  HFF's 33,613 labels carry information the age head uses, or
    (ii) 75 labels is simply too few to learn from, whatever they contain.

Arm C separates them. It is arm A in every respect -- same cells, same 33,688 labels, same
`age_window_k = 4`, same seeds -- with only the cell<->label PAIRING destroyed inside HFF. Label
volume is held constant; information content is destroyed.

    C ranks like A  ->  the gain was volume / trunk regularisation; HFF's labels are uninformative
    C ranks like B  ->  HFF's labels carry real signal despite the artefact

WHY THIS NEEDS AN EQUIVALENCE MARGIN, NOT A CI CONTAINING ZERO
--------------------------------------------------------------
"C ranks like A" is an EQUIVALENCE claim. A paired CI that happens to contain zero is *absence of
evidence*, not evidence of absence -- exactly the trap that made step 6's null uninterpretable.
Equivalence needs a margin fixed in advance and a TOST: declare equivalence only if the **90% CI
lies entirely inside [-DELTA_EQ, +DELTA_EQ]**.

DELTA_EQ IS DERIVED, NOT PICKED
-------------------------------
The interval this experiment must partition is the measured A-B ranking gap, **0.0688**. Half of it
places the decision boundary at the midpoint: if |C-A| < DELTA_EQ then C sits in A's half, and if
|C-B| < DELTA_EQ it sits in B's half.

    DELTA_EQ = |A-B ranking gap| / 2 = 0.0344

Using the project's own measured gap rather than a round number keeps this consistent with how
DELTA* was derived from Stage 2 sec 12 for the primary metric.

WHY THE MARGIN IS NOT WIDENED TO WHATEVER WOULD PASS
----------------------------------------------------
The smallest margin resolvable at the pessimistic SD is **0.0523 = 0.76x the A-B gap**. Adopting it
would declare arm C "equivalent to A" even if C sat three quarters of the way to B -- a bar almost
nothing could fail, which answers nothing. The principled midpoint is kept and its resolvability
reported as CONDITIONAL instead. A margin reverse-engineered to pass is not a bar.

RESOLVABILITY IS THE POINT OF THIS SCRIPT
-----------------------------------------
Ground rule sec 5b: a bar with no resolvability test is not pre-registered. So: **if arm C truly
equals arm A, how often would this bar actually say so at 6 folds?** If the answer is not >= 95%,
the bar is UNRESOLVABLE and arm C cannot deliver the conclusion it is being run for -- which must be
known BEFORE spending ~5 h, not after.
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

# Measured in the step-6 rerun: per-fold rank_model_dage, arm B - arm A.
AB_RANK_DIFFS = [-0.0649, -0.0338, -0.0662, -0.1091, -0.0982, -0.0403]


def ab_gap() -> tuple[float, float]:
    """The measured A-B ranking gap and its per-fold SD. Pure."""
    a = np.asarray(AB_RANK_DIFFS, dtype=float)
    return float(a.mean()), float(a.std(ddof=1))


def delta_eq() -> float:
    """Equivalence margin: half the measured A-B gap, so it partitions [A, B] at the midpoint."""
    return abs(ab_gap()[0]) / 2.0


def tost_equivalent(diffs: np.ndarray, margin: float, alpha: float = 0.05) -> np.ndarray:
    """Vectorised TOST: is the (1-2*alpha) CI entirely inside +/- margin? Pure.

    Rows are simulated runs of `N_FOLDS` paired differences. The 90% CI (alpha=0.05 each side) is
    the standard equivalence interval; using the 95% CI here would be conservative in the wrong
    direction for an equivalence claim.
    """
    n = diffs.shape[1]
    mean = diffs.mean(axis=1)
    se = diffs.std(axis=1, ddof=1) / np.sqrt(n)
    half = float(stats.t.ppf(1 - alpha, n - 1)) * se
    return (mean - half > -margin) & (mean + half < margin)


def resolvability(margin: float, sd: float, true_effect: float = 0.0,
                  rng: np.random.Generator | None = None) -> float:
    """P(the bar declares equivalence) when the truth IS equivalence. Pure."""
    rng = rng or np.random.default_rng(SEED)
    d = rng.normal(true_effect, sd, size=(N_SIM, N_FOLDS))
    return float(tost_equivalent(d, margin).mean())


def smallest_resolvable_margin(sd: float, rng: np.random.Generator | None = None,
                               lo: float = 0.001, hi: float = 0.5, iters: int = 30) -> float:
    """Smallest margin at which a truly-equivalent arm C clears MIN_PASS_RATE. Pure, bisection."""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if resolvability(mid, sd, rng=rng) >= MIN_PASS_RATE:
            hi = mid
        else:
            lo = mid
    return hi


def main() -> int:
    rng = np.random.default_rng(SEED)
    gap, sd_ab = ab_gap()
    margin = delta_eq()

    print("\nARM C — registering the EQUIVALENCE bar, before the run\n")
    print(f"  measured A-B ranking gap : {gap:+.4f}  (per-fold SD {sd_ab:.4f}, n={N_FOLDS})")
    print(f"  DELTA_EQ = |gap| / 2     : {margin:.4f}")
    print("  test                     : TOST, 90% CI must lie inside +/-DELTA_EQ")
    t90 = float(stats.t.ppf(0.95, N_FOLDS - 1))
    print(f"  resolution floor         : 90% CI half-width = {t90 * sd_ab / np.sqrt(N_FOLDS):.4f}"
          f"  at SD {sd_ab:.4f}\n")

    # The C-A SD is unknown. It is plausibly SMALLER than the A-B SD (both arms carry the full
    # 33,688 labels and differ only in pairing), so sweep rather than assume one value.
    print(f"  {'SD(C-A)':>9}{'P(declare equivalence | C truly = A)':>40}   verdict")
    rows = []
    for sd in (0.005, 0.010, 0.015, 0.020, sd_ab, 0.040):
        rate = resolvability(margin, sd, rng=rng)
        v = bar_verdict(np.array([rate]), MIN_PASS_RATE, lower_is_better=False)
        ok = rate >= MIN_PASS_RATE
        tag = "  <- the A-B SD (pessimistic)" if abs(sd - sd_ab) < 1e-9 else ""
        rows.append({"sd": sd, "pass_rate": rate, "resolvable": ok, "verdict": v["verdict"]})
        print(f"  {sd:>9.4f}{rate:>40.1%}   {'RESOLVABLE' if ok else 'UNRESOLVABLE'}{tag}")

    sd_needed = None
    for r in rows:
        if r["resolvable"]:
            sd_needed = r["sd"]
            break
    smallest = smallest_resolvable_margin(sd_ab, rng=rng)

    print(f"\n  At the pessimistic SD {sd_ab:.4f}, DELTA_EQ = {margin:.4f} is "
          f"{'RESOLVABLE' if resolvability(margin, sd_ab, rng=rng) >= MIN_PASS_RATE else 'NOT RESOLVABLE'}.")
    print(f"  Smallest margin that WOULD be resolvable at that SD: {smallest:.4f} "
          f"({smallest / abs(gap):.2f}x the A-B gap).")
    if sd_needed is not None:
        print(f"  DELTA_EQ = {margin:.4f} becomes resolvable once SD(C-A) <= ~{sd_needed:.4f}.")

    # A false-positive check, mirroring the DELTA* bar: a C that truly sits at B must NOT be
    # declared equivalent to A.
    fp = resolvability(margin, sd_ab, true_effect=gap, rng=rng)
    print(f"\n  false-equivalence check: C truly AT arm B is declared 'like A' "
          f"{fp:.1%} of the time (must be small)")

    out = {"script": "register_arm_c_bar",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "ab_rank_gap": gap, "ab_rank_sd": sd_ab, "n_folds": N_FOLDS,
           "delta_eq": margin, "delta_eq_derivation": "half the measured A-B rank_model_dage gap",
           "test": "TOST, 90% CI inside +/- delta_eq",
           "resolution_floor_at_ab_sd": t90 * sd_ab / np.sqrt(N_FOLDS),
           "resolvability_by_sd": rows,
           "smallest_resolvable_margin_at_ab_sd": smallest,
           "sd_at_which_delta_eq_resolves": sd_needed,
           "false_equivalence_rate_if_C_is_really_B": fp,
           "resolvable_at_pessimistic_sd": bool(
               resolvability(margin, sd_ab, rng=rng) >= MIN_PASS_RATE)}
    (_RESULTS / "register_arm_c_bar_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote register_arm_c_bar_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
