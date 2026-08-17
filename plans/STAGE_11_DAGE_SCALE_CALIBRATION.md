# Stage 11 — Is the ΔAge scale error fixable by calibration?

**Status:** PLAN. `src/` is not touched by this stage under any outcome; the most a positive result
buys is the right to PROPOSE a Change, separately pre-registered (same rule as Stage 10 §10.5).

**Opened because** Stage 10 established that `raw`'s problem is **not** pluripotency — it is
**scale**. Raw ΔAge has the right ordering (Spearman 0.770 vs methylation) and a 66 % magnitude
inflation (SD ratio 1.66). Nobody has tried simply *correcting the scale*.

---

## 11.0 What is already established, and what follows from it

| variant | MAE | SD ratio | ρ vs methylation |
|---|---|---|---|
| raw | 22.69 | **1.66** | **0.770** |
| resid_pluri | 13.00 | 1.14 | 0.354 ← Stage 10: signal deleted |
| ranknorm | 10.15 | 0.30 | 0.138 ← collapsed |
| top100 | **7.15** | **0.98** | **0.810** |
| *methylation-vs-methylation floor* | *7.30* | *1.00* | *0.613* |

The gap between `raw` and `top100` on MAE is 15.5 yr. Two things could produce it: **scale** (raw is
1.66× too big) or **ordering** (top100 orders better, 0.810 vs 0.770). This stage separates them.

---

## 11.1 A mathematical fact, stated BEFORE the run

**A pure linear rescale `y → k·y` cannot change Spearman.** Rank order is invariant to any positive
monotone transform. Therefore:

- Rescaled-raw's ρ **will be exactly 0.770**. That is arithmetic, not a result, and it must not be
  reported as a finding.
- Whatever MAE improvement calibration produces is therefore attributable to **scale alone**.
- The residual difference from top100's 0.810 is the part scale **cannot** explain.

This is why the stage is worth running: it cleanly partitions the 15.5 yr gap.

---

## 11.2 The calibration, and the leak that must not happen

`ΔAge_cal = k · ΔAge_raw`, with `k` fitted to minimise error against methylation ΔAge.

**`k` must be fitted LEAVE-ONE-DONOR-OUT.** Fitting `k` on all 44 conditions and then scoring those
same conditions is circular — it would guarantee an improvement and measure nothing. Fold on
**donor** (3 in the transient arm), fit `k` on two, apply to the third.

Also fitted and reported: an **offset** `ΔAge_cal = k·raw + c`, because a scale error and a mean
shift are different defects and the data should say which is present. Both variants reported; the
offset one is NOT preferred merely for fitting better.

---

## 11.3 PRE-REGISTERED READING

Let `FLOOR = 7.30` (methylation-vs-methylation MAE) and `FLOOR_MULT = 1.5`.

| outcome | verdict |
|---|---|
| LODO-calibrated raw MAE ≤ `FLOOR_MULT × FLOOR` (= 10.95) | **SCALE IS THE PROBLEM** — a one-parameter fix recovers most of the gap; the dense clock was never broken, only mis-scaled |
| MAE improves but stays above 10.95 | **SCALE IS PART OF IT** — report how much of the 15.5 yr gap closes |
| MAE does not improve | **NOT SCALE** — the inflation is not a simple multiplicative error |

Separately, and independent of the above:

- **`k` STABILITY.** Report `k` per fold. If it varies by more than 2× across three donors,
  calibration is not transferable and that caveat outranks any MAE gain.
- **top100 STILL WINS ON ORDERING** is expected (0.810 > 0.770) and is *not* evidence about scale.
  It is reported so the two effects stay separated.

---

## 11.4 What this stage cannot settle

- **n = 3 donors** for the LODO calibration. A `k` fitted on two donors is a weak estimate.
- **Transient arm only** (44 conditions carrying both methylation truths). No Sendai condition
  carries both, so there is no independent cohort to confirm `k` on.
- Whether `k` transfers to a *new* dataset is untestable here and must not be claimed.
- A successful rescale does **not** rescue same-timepoint ΔAge PREDICTION, which is circular
  (ρ 0.96–0.99) regardless of how the target is scaled.

---

## 11.5 Verification

| item | how |
|---|---|
| diagnostic | `experiments/diag_stage11_scale.py`, read-only, writes `results/diag_stage11_scale_results.json` |
| tests | `tests/test_diag_stage11_scale.py` — the rank-invariance fact, the LODO no-leak property, and every verdict branch, on constructed input |
| record | `CHANGES.md` with the verdict and the `k` stability caveat |
| suite | full `pytest` green + `ruff` clean before commit |
