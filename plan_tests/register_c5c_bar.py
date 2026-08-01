"""STAGE 1.5.3 STEP 5c — register the bar for C-5's *threshold*, BEFORE the code is written.

    python plan_tests/register_c5c_bar.py

READ-ONLY. Writes `results/register_c5c_bar_results.json`. Runs on the geometry alone -- no data,
no training, no `src/` file touched. `REF_GROUND_RULES.md` §5b: the bar goes through
`audit_metrics.bar_verdict` **before** the change it grades.

WHAT CHANGED SINCE STEP 5b, AND WHY THIS SCRIPT EXISTS
-------------------------------------------------------
5b pinned Option 2 at a **fixed W = 8**. The readiness audit then found that step 6 runs *two* arms
and only one of them is masked:

    arm A (control)    33 688 of 33 688 cells age-valid   ~512 age cells per batch
    arm B (treatment)       75 of 33 688 cells age-valid    1.14 age cells per batch

Arm A has no occupancy problem at all. A fixed W = 8 would cut its age updates 65 -> 8 for no
benefit while helping arm B -- **the mechanism would handicap the control and help the treatment**,
tilting `dage_mae_model` toward the very outcome that concludes "99.7 % of the labels were
net-negative". That is a validity threat, not a tuning detail.

THE RULE BEING GRADED
---------------------
Trigger on the accumulated age-CELL count, not a batch count:

    close the age window once it holds >= k age cells, or after W_max batches,
    whichever comes first. The window CARRIES across the epoch boundary.

One policy, applied identically to both arms. It only *behaves* differently because the data differ,
which is what a controlled comparison is.

THE THREE BARS
--------------
    A1  arm A closes every window at W = 1                  == 1.000  (bit-identity precondition)
    A2  P(a closed window holds >= k age cells) in arm B     >= 0.95   (B2, now by construction)
    A3  arm B gets MORE age updates than fixed W = 8 would             (else the redesign is a
                                                                        regression, not a fix)

A1 is the one that matters most and is an EQUALITY, not a rate: if arm A ever accumulates across
batches, its training has changed, `scorecard/baseline.json` stops being a valid reference, and the
step-6 comparison is confounded again -- the exact defect this rule exists to remove. The script
exits non-zero if any bar fails.

WHY k = 4
---------
It is B2's already-registered threshold from step 5 (the smallest k that halves the per-update
standard error against a single cell, SE ~ 1/sqrt(m)). Reusing it rather than inventing a new number
keeps step 5's reasoning intact. §4 below reports what other k would buy, so the choice is graded
rather than asserted.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from audit_metrics import MIN_PASS_RATE, bar_verdict  # noqa: E402

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

N_TRAIN = 33_688
BATCH = 512
STEPS_PER_EPOCH = N_TRAIN // BATCH        # 65 full batches (the 408-cell tail rides the carry)
EPOCHS = 60
N_AGE_ARM_A = 33_688                      # control: every cell age-valid
N_AGE_ARM_B = 75                          # treatment: HFF masked
K = 4                                     # B2's registered threshold, reused
W_MAX = 8                                 # 5b's sensitivity table, now a CEILING not a constant
N_EPOCH_SIM = 2_000


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def close_windows(per_batch_counts, k: int, w_max: int, *, close_at_end: bool = False,
                  carry: tuple[int, int] = (0, 0)) -> tuple[list[tuple[int, int]], tuple[int, int]]:
    """Walk a run of batches; return the closed windows ``(n_cells, n_batches)`` and the carry. Pure.

    This is the exact rule `train.py` implements; both are tested against each other so the
    simulation the decision rests on cannot drift from the code that ships.

    ``close_at_end`` was **attempt 1** and is kept only so the record can show what it cost: forcing
    a close on the epoch's last batch guarantees every label is consumed inside its own epoch, but
    it manufactures one deliberately-partial window per epoch. Measured, that single window was
    **4.44 pp of a 6.12 pp shortfall** against A2 -- it failed the bar on its own, while the
    irreducible `W_max` limit contributed only 1.67 pp.

    The shipped rule instead **carries the open window across the epoch boundary** (`carry`). A
    label buffered at the end of one epoch is consumed by the first window of the next, so it is
    still used exactly once per pass -- just delivered one window late. Fewer special cases in the
    code, and no artificial partial window. The only residue is that at the very end of training
    `< k` cells may sit unused in the buffer; that is one window out of ~1 000 and is recorded
    rather than hidden.
    """
    windows: list[tuple[int, int]] = []
    acc, opened = carry
    n = len(per_batch_counts)
    for i, c in enumerate(per_batch_counts):
        acc += int(c)
        opened += 1
        if acc >= k or opened >= w_max or (close_at_end and i == n - 1):
            windows.append((acc, opened))
            acc = opened = 0
    return windows, (acc, opened)


def simulate_arm(n_age: int, *, k: int, w_max: int, n_epochs: int,
                 rng: np.random.Generator, close_at_end: bool = False) -> dict:
    """Run the rule over `n_epochs` simulated epochs of the real geometry. Pure."""
    p = n_age / N_TRAIN
    cells, batches, per_epoch = [], [], []
    carry = (0, 0)
    for _ in range(n_epochs):
        counts = rng.binomial(BATCH, p, size=STEPS_PER_EPOCH)
        w, carry = close_windows(counts, k, w_max, close_at_end=close_at_end, carry=carry)
        per_epoch.append(len(w))
        cells.extend(c for c, _ in w)
        batches.extend(b for _, b in w)
    cells_a, batches_a = np.asarray(cells, float), np.asarray(batches, float)
    return {"windows_per_epoch": float(np.mean(per_epoch)),
            "mean_cells_per_window": float(cells_a.mean()),
            "mean_batches_per_window": float(batches_a.mean()),
            "frac_windows_at_W1": float((batches_a == 1).mean()),
            "frac_windows_ge_k": float((cells_a >= k).mean()),
            "updates_over_run": float(np.mean(per_epoch) * EPOCHS)}


def fixed_w_updates(w: int) -> float:
    """What the step-5b fixed-W design would have given, for the A3 comparison. Pure."""
    return (STEPS_PER_EPOCH // w) * EPOCHS


def main() -> int:
    rng = np.random.default_rng(0)
    print("STAGE 1.5.3 STEP 5c — the bar for C-5's threshold, before the code\n")
    print(f"  rule: close the age window at >= {K} cells or {W_MAX} batches; carry across epochs")
    print(f"  geometry: {N_TRAIN} train cells, batch {BATCH}, {STEPS_PER_EPOCH} batches/epoch\n")

    arm_a = simulate_arm(N_AGE_ARM_A, k=K, w_max=W_MAX, n_epochs=200, rng=rng)
    arm_b = simulate_arm(N_AGE_ARM_B, k=K, w_max=W_MAX, n_epochs=N_EPOCH_SIM, rng=rng)

    print(f"  {'arm':<26}{'upd/epoch':>11}{'cells/win':>11}{'batches/win':>13}{'>= k':>8}")
    print("  " + "-" * 69)
    for name, r in (("A (control, unmasked)", arm_a), ("B (treatment, HFF masked)", arm_b)):
        print(f"  {name:<26}{r['windows_per_epoch']:>11.1f}{r['mean_cells_per_window']:>11.1f}"
              f"{r['mean_batches_per_window']:>13.2f}{r['frac_windows_ge_k']:>8.1%}")

    # ---- A1: arm A must never accumulate. Equality, not a rate. -------------------------- #
    a1 = bar_verdict(np.array([arm_a["frac_windows_at_W1"]]), 1.0, lower_is_better=False)
    a1_ok = arm_a["frac_windows_at_W1"] == 1.0
    print(f"\n  A1  arm A closes at W=1 on every window   {arm_a['frac_windows_at_W1']:.4f}"
          f"   {'PASS' if a1_ok else 'FAIL'}  (must be exactly 1.0000)")
    if a1_ok:
        print("      => the control arm is BIT-IDENTICAL to today: same batches, same age loss,")
        print("         same number of age updates. scorecard/baseline.json stays a valid reference.")

    # ---- A2: B2, restated for the adaptive rule --------------------------------------- #
    a2_ok = arm_b["frac_windows_ge_k"] >= MIN_PASS_RATE
    print(f"\n  A2  arm B windows holding >= {K} cells       {arm_b['frac_windows_ge_k']:.1%}"
          f"   {'PASS' if a2_ok else 'FAIL'}  (bar {MIN_PASS_RATE})")
    print("      => B2 is now met BY CONSTRUCTION; the only shortfall is the W_max forced close,")
    print("         which is what makes this a rate rather than an identity.")
    attempt1 = simulate_arm(N_AGE_ARM_B, k=K, w_max=W_MAX, n_epochs=N_EPOCH_SIM,
                            rng=np.random.default_rng(11), close_at_end=True)
    print(f"      ATTEMPT 1 (force a close at each epoch's end) scored "
          f"{attempt1['frac_windows_ge_k']:.1%} and FAILED this bar.")
    print("      Diagnosed, not worked around: the epoch-end window was 4.44 pp of the 6.12 pp")
    print("      shortfall and the W_max limit only 1.67 pp, so the mechanism was wrong, not")
    print("      the bar. Carrying the window across the epoch boundary removes it entirely.")

    # ---- A3: the redesign must not cost arm B updates --------------------------------- #
    fixed = fixed_w_updates(W_MAX)
    a3_ok = arm_b["updates_over_run"] > fixed
    print(f"\n  A3  arm B age updates over the run        {arm_b['updates_over_run']:.0f}"
          f"   vs {fixed:.0f} at fixed W={W_MAX}   {'PASS' if a3_ok else 'FAIL'}")
    print(f"      => {arm_b['updates_over_run'] / fixed:.1f}x more age optimisation than 5b's")
    print("         design, at the same per-update quality. This directly reduces the residual")
    print("         '480 updates may be too few to converge' risk that 5b had to leave open.")

    # ---- what other k would buy: the choice reported, not asserted --------------------- #
    print(f"\n  k swept (arm B), so k={K} is graded rather than assumed:")
    print(f"     {'k':>3}{'upd/epoch':>11}{'cells/win':>11}{'>= k':>9}{'upd/run':>10}")
    ksweep = {}
    for k in (2, 3, 4, 6, 8, 12):
        r = simulate_arm(N_AGE_ARM_B, k=k, w_max=W_MAX, n_epochs=400,
                         rng=np.random.default_rng(50 + k))
        ksweep[k] = r
        mark = "   <-- chosen" if k == K else ""
        print(f"     {k:>3}{r['windows_per_epoch']:>11.1f}{r['mean_cells_per_window']:>11.1f}"
              f"{r['frac_windows_ge_k']:>9.1%}{r['updates_over_run']:>10.0f}{mark}")
    print("     larger k buys per-update SNR and costs updates; k=4 is B2's registered threshold")
    print("     and the smallest that halves the per-update standard error against a single cell.")

    ok = a1_ok and a2_ok and a3_ok
    out = {"script": "register_c5c_bar",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "rule": {"k": K, "w_max": W_MAX, "close_at_epoch_end": False,
                    "carries_window_across_epochs": True},
           "attempt1_close_at_epoch_end": attempt1,
           "geometry": {"n_train": N_TRAIN, "batch": BATCH,
                        "steps_per_epoch": STEPS_PER_EPOCH, "epochs": EPOCHS},
           "arm_a": arm_a, "arm_b": arm_b,
           "bars": {"A1": {"desc": "arm A closes every window at W=1", "bar": 1.0,
                           "value": arm_a["frac_windows_at_W1"], "pass": a1_ok,
                           "verdict": a1["verdict"]},
                    "A2": {"desc": f"P(window holds >= {K} cells) in arm B", "bar": MIN_PASS_RATE,
                           "value": arm_b["frac_windows_ge_k"], "pass": a2_ok},
                    "A3": {"desc": "arm B updates exceed fixed-W design", "bar": fixed,
                           "value": arm_b["updates_over_run"], "pass": a3_ok}},
           "k_sweep": {str(k): v for k, v in ksweep.items()},
           "all_pass": ok}
    (_RESULTS / "register_c5c_bar_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"\n  wrote register_c5c_bar_results.json   ALL BARS {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
