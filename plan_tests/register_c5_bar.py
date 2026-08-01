"""STAGE 1.5.3 STEP 5 — register C-5's bar, and let it choose between the three options.

    python plan_tests/register_c5_bar.py

READ-ONLY. Writes `results/register_c5_bar_results.json`. Touches no `src/` file and moves no
label. Runs entirely on the geometry — no training, no data files.

WHY A BAR AT ALL, AND WHY THIS ONE
-----------------------------------
`STAGE_1_5_3_EXECUTE.md` §5 step 5: *"C-5 design + its bar through `audit_metrics.bar_verdict`;
bar RESOLVABLE **before any retrain**."* That constraint decides what the bar can measure.
`dage_mae_model` is unavailable — it needs the retrain that step 6 performs — so the bar has to
grade the **mechanism** rather than the outcome.

The mechanism's job is stated exactly in C-5: with HFF masked, 75 age labels sit among 33 688
training cells, uniform shuffling puts **1.14** of them in a 512-cell batch, and `losses.py:55-57`
returns a hard zero for a batch with none. So the measurable question is:

    per optimiser update, does the age task receive a gradient, and is that gradient computed
    over more than one or two cells?

TWO BARS, because "non-zero" is not the same as "usable"
--------------------------------------------------------
    B1  P(update contributes ANY age gradient)        >= 0.95
    B2  P(update's age gradient uses >= 4 cells)      >= 0.95

**B1 alone would be too easy.** C-5's diagnosis is not only the 32 % of empty batches, it is also
that the surviving ones carry *"a Huber loss over one or two cells"* and that `MultiTaskLoss` learns
`s_age` from them. A mechanism can clear B1 and still feed the optimiser per-cell noise, so B2 is
what separates "a gradient exists" from "a gradient means something".

`MIN_PASS_RATE = 0.95` is not invented here — it is `audit_metrics`' own project-wide standard, and
both bars are graded through `bar_verdict` like every other registered bar.

WHY k = 4 IN B2, RATHER THAN A ROUNDER NUMBER
----------------------------------------------
A perfectly even spread of 75 labels over the 66 updates in an epoch gives **1.14** per update, so
any k >= 2 is a bar only an **oversampling or accumulating** mechanism can meet — which is the whole
point of asking. k = 4 halves the per-update standard error relative to a single cell (SE scales as
1/sqrt(m), so 1 -> 4 cells is a 2x reduction) and is the smallest value that does so. Stating the
reasoning matters more than the value: the bar is chosen to *discriminate between the options*, and
§3 below reports what each one scores rather than only whether it passed.

RESOLVABILITY, AND THE DIRECTION THAT ACTUALLY MATTERS
------------------------------------------------------
The system that "meets the intent exactly" is today's **dense** regime — before masking, every cell
carries an age label, so every update is full. That passes both bars at 100 %, which makes them
resolvable but says little.

**The useful check is discrimination:** a bar that every candidate passes, or that none can, decides
nothing. §3 therefore reports all three options against both bars, and the run FAILS if the bars do
not separate them.
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

# --- the real geometry, every number cited from the plan's evidence table ------------- #
N_TRAIN = 33_688          # E11: training-split cells
N_AGE_MASKED = 75         # E11: age-valid labels once HFF is masked
BATCH = 512               # E14: configs/train/default.yaml
N_SIM = 20_000
RNG = np.random.default_rng(0)

MIN_CELLS = 4             # B2's k -- see the module docstring
B1_BAR = MIN_PASS_RATE    # 0.95
B2_BAR = MIN_PASS_RATE


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def age_cells_per_update(n_age: int, n_train: int, batch: int, *, weight: float = 1.0,
                         accumulate: int = 1, n_sim: int = N_SIM,
                         rng: np.random.Generator | None = None) -> np.ndarray:
    """How many age-valid cells reach one optimiser update. Pure.

    Covers all three candidate mechanisms through two knobs, so they are compared on one
    footing rather than three ad-hoc simulations:

      * ``weight``     -- Option 1. A `WeightedRandomSampler` weight of `w` on age-valid cells
                          makes each draw hit one with probability `n_age*w / (n_age*w + rest)`.
                          `weight=1` is today's uniform shuffling.
      * ``accumulate`` -- Option 2. The age loss is summed over `W` batches before it is
                          stepped, so ONE update sees `batch*W` draws.
      * Option 3 changes neither: it pins `s_age` and leaves the sampling alone, so it is
        `weight=1, accumulate=1` -- identical to the status quo by construction, which is
        exactly the criticism C-5 makes of it.
    """
    rng = rng or RNG
    rest = n_train - n_age
    p = (n_age * weight) / (n_age * weight + rest)
    return rng.binomial(batch * accumulate, p, size=n_sim).astype(float)


def required_weight(n_age: int, n_train: int, batch: int, target_cells: float) -> float:
    """The Option-1 sampler weight that puts `target_cells` age cells in an average batch. Pure.

    Solving `batch * n_age*w / (n_age*w + rest) = target` for `w`.
    """
    rest = n_train - n_age
    if target_cells >= batch:
        return float("inf")
    return float(target_cells * rest / (n_age * (batch - target_cells)))


def oversampling_cost(n_age: int, n_train: int, weight: float) -> dict:
    """What Option 1 costs the FATE task. Pure.

    The fate head consumes every cell in the batch, so raising the age cells' draw probability
    changes the fate task's training distribution too. C-5 flags this and refuses to wave it
    through; this quantifies it before the run rather than after.
    """
    rest = n_train - n_age
    p_base = n_age / n_train
    p_new = (n_age * weight) / (n_age * weight + rest)
    return {"weight": weight,
            "age_share_before": p_base, "age_share_after": p_new,
            "fold_oversampled": p_new / p_base,
            "extra_share_of_batch": p_new - p_base}


def discriminates(scores: dict[str, dict]) -> bool:
    """Does the bar separate the candidates? A bar everything passes decides nothing. Pure."""
    verdicts = {k: (v["B1"]["verdict"], v["B2"]["verdict"]) for k, v in scores.items()}
    return len(set(verdicts.values())) > 1


def main() -> int:
    print("STAGE 1.5.3 STEP 5 — C-5's bar, registered before any retrain\n")
    print(f"  geometry: {N_AGE_MASKED} age labels among {N_TRAIN} cells, batch {BATCH}")
    print(f"            -> {N_AGE_MASKED * BATCH / N_TRAIN:.2f} age cells in an average batch")
    print(f"  B1  P(update contributes ANY age gradient)   >= {B1_BAR}")
    print(f"  B2  P(age gradient uses >= {MIN_CELLS} cells)          >= {B2_BAR}\n")

    w_needed = required_weight(N_AGE_MASKED, N_TRAIN, BATCH, target_cells=8.0)
    candidates = {
        "status quo (uniform shuffling)":      dict(weight=1.0, accumulate=1),
        "Option 3 (pin s_age only)":           dict(weight=1.0, accumulate=1),
        "Option 2 (accumulate, W=8)":          dict(weight=1.0, accumulate=8),
        f"Option 1 (sampler, w={w_needed:.1f})": dict(weight=w_needed, accumulate=1),
    }

    print(f"  {'candidate':<34}{'mean cells':>11}{'B1':>9}{'B2':>9}   verdict")
    print("  " + "-" * 78)
    scores: dict[str, dict] = {}
    for name, kw in candidates.items():
        sim = age_cells_per_update(N_AGE_MASKED, N_TRAIN, BATCH, **kw)
        b1 = bar_verdict((sim >= 1).astype(float), 1.0, lower_is_better=False)
        b2 = bar_verdict((sim >= MIN_CELLS).astype(float), 1.0, lower_is_better=False)
        scores[name] = {"config": kw, "mean_cells": float(sim.mean()), "B1": b1, "B2": b2}
        both = "PASS" if b1["verdict"] == b2["verdict"] == "RESOLVABLE" else "FAIL"
        print(f"  {name:<34}{sim.mean():>11.2f}{b1['pass_rate']:>9.1%}{b2['pass_rate']:>9.1%}   {both}")

    # the intent-satisfying reference: the DENSE regime, before any masking
    dense = age_cells_per_update(N_TRAIN, N_TRAIN, BATCH)
    d1 = bar_verdict((dense >= 1).astype(float), 1.0, lower_is_better=False)
    d2 = bar_verdict((dense >= MIN_CELLS).astype(float), 1.0, lower_is_better=False)
    print("\n  reference -- the DENSE regime a correct system enjoys today:")
    print(f"     B1 {d1['pass_rate']:.1%}   B2 {d2['pass_rate']:.1%}   "
          f"-> {'RESOLVABLE' if d1['verdict'] == d2['verdict'] == 'RESOLVABLE' else 'NOT RESOLVABLE'}")

    sep = discriminates(scores)
    print(f"\n  does the bar DISCRIMINATE between the options? {'yes' if sep else 'NO'}")
    if not sep:
        print("  ==> ABORT: a bar every candidate passes (or none can) decides nothing.")

    cost = oversampling_cost(N_AGE_MASKED, N_TRAIN, w_needed)
    print("\n  Option 1's cost to the FATE task, quantified before the run:")
    print(f"     age cells are {cost['fold_oversampled']:.1f}x oversampled "
          f"({cost['age_share_before']:.4%} -> {cost['age_share_after']:.4%} of each batch)")
    print(f"     so the fate head's training mix shifts by "
          f"{cost['extra_share_of_batch']:.2%} of the batch.")
    print("     GUARD: fate_prauc, fate_roc and fate_ece must read 'noise' in the paired")
    print("     6-fold comparison at step 6. A move there is a finding, not a trade-off.")

    out = {"script": "register_c5_bar",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "geometry": {"n_train": N_TRAIN, "n_age_masked": N_AGE_MASKED, "batch": BATCH,
                        "mean_cells_uniform": N_AGE_MASKED * BATCH / N_TRAIN},
           "bars": {"B1": {"desc": "P(update contributes any age gradient)", "bar": B1_BAR},
                    "B2": {"desc": f"P(age gradient uses >= {MIN_CELLS} cells)", "bar": B2_BAR,
                           "min_cells": MIN_CELLS}},
           "reference_dense": {"B1": d1, "B2": d2},
           "candidates": scores, "discriminates": sep,
           "required_weight_for_8_cells": w_needed,
           "option1_fate_cost": cost}
    (_RESULTS / "register_c5_bar_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote register_c5_bar_results.json")
    return 0 if sep else 1


if __name__ == "__main__":
    raise SystemExit(main())
