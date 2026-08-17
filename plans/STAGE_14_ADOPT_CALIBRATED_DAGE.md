# Stage 14 — adopting the calibrated ΔAge

**Status:** PRE-REGISTRATION. `src/` is **not** touched by this document. It states what the
Change would be, what it would and would not buy, and the guards it must pass — so the decision is
made on measured grounds rather than on the appeal of a 63 % drop in ΔAge MAE.

**A pre-flight was run before writing this plan** (`experiments/diag_stage14_calibration_
equivariance.py`, 14 tests). It changed the plan's central recommendation. §14.3 records what I
expected and what the data said instead.

---

## 14.1 What Stage 11 licenses, and what it does not

Stage 11 measured, LODO over 3 donors, against the 7.30 yr methylation-vs-methylation floor:

| | MAE | SD ratio | ρ |
|---|---|---|---|
| raw ΔAge | 22.69 | 1.66 | 0.770 |
| **raw × k, least squares** (k̄ = 0.364) | **6.78** | **0.597** | 0.770 |
| **raw × k, variance-matched** (k̄ = 0.599) | 10.68 | **0.990** | 0.770 |
| *methylation floor* | *7.30* | *1.00* | *0.613* |

Verdict **SCALE IS THE PROBLEM**; `k` stable across donors (spread 1.19×, bar was 2×).

**ρ is 0.770 in all three rows by arithmetic** — a positive rescale cannot change rank order.
Calibration therefore cannot improve ordering, cannot improve model skill, and cannot rescue the
same-timepoint prediction circularity. **It corrects units.** That is worth doing, because every
absolute claim the project makes — "OSKM rejuvenates by X years", and every conformal interval in
years — is currently inflated by roughly 2.7×.

**What Stage 11 explicitly did not license** (§11.4): that `k` transfers. It was fitted on donors
**O1/O2/O3 of the transient arm only** — the only rows in `results/dage_ledger.csv` carrying
`TRUTH_meth_dage_mt` (68 of 90 rows; the Sendai cohort carries none). Applying it to the
Sendai + HFF training pipeline is an out-of-cohort transfer, and out-of-cohort transfer is the
thing this project has already watched fail catastrophically once.

## 14.2 Which `k` — and why the better-fitting one is the wrong one

`k_LS = ρ · SD(truth)/SD(pred)`. With ρ < 1 it is **strictly smaller** than the variance-matching
factor: it wins on MAE by **under-reporting magnitude** (SD ratio 0.597 — it shrinks the spread by
40 %). That is the shrinkage trap this project already hit once, where MAE improved while the
ordering degraded.

For a **reporting** transform the objective is an unbiased magnitude, not a minimal MAE.
**`k_var = 0.599` is the correct choice**, and its worse MAE (10.68 vs 6.78) is the honest number,
not a defeat. Adopting `k_LS` because 6.78 beats the 7.30 floor would be selecting the estimator
that flatters the headline.

## 14.3 THE PRE-FLIGHT — what I expected, and what the data said

**Expected:** `huber_delta = 2.0` yr sits far below a ΔAge SD of 13–23, so essentially every
residual would be outside the knee, the loss would be effectively L1, and rescaling the target
would be a near-pure units change.

**Measured** (ridge in-sample training residuals, five folds):

| fold | med \|resid\| | → after ×k | % inside the 2.0 yr knee | → after |
|---|---|---|---|---|
| N3 | 2.27 | 0.83 | 44.9 % | **87.2 %** |
| O1 | 2.39 | 0.87 | 42.9 % | **85.0 %** |
| O2 | 2.37 | 0.86 | 43.2 % | **85.4 %** |
| Y1 | 1.36 | 0.50 | 66.6 % | **97.3 %** |
| Y2 | 2.40 | 0.87 | 42.9 % | **85.1 %** |

**The expectation was wrong.** Training residuals are **comparable to** the knee, not far beyond
it. The loss today is a genuine Huber mix (~43–67 % quadratic); after rescaling it becomes
**85–97 % quadratic — effectively plain MSE.**

**Consequence:** rescaling the training target is **not** a units change for the neural model. It
silently converts a robust loss into an outlier-sensitive one, at the same time as changing units,
and the two effects would be inseparable in the result. That is a two-change experiment wearing the
costume of a one-change experiment.

*Limitation, stated:* the fractions come from ridge in-sample residuals standing in for the
network's training residuals; the exact percentages depend on the network's fit. The **direction**
is monotone and therefore robust (a test pins this) — only the magnitude is indicative.

**Established in the same run, exactly:** ridge is equivariant to 1e-12. `MAE(k·y) = k·MAE(y)`,
`Δρ = 0.00e+00`, on every fold. So for the linear path, calibration buys units and nothing else —
measured, not argued.

## 14.4 The three options, and the recommendation

### (a) Rescale the training target `y_age → k·y_age`. **NOT RECOMMENDED.**

Requires a rebuild + retrain + full re-score, restarts the guard record, and — per §14.3 — changes
the loss regime as a side effect. If it were done anyway, `huber_delta` **must** be rescaled by the
same `k` (2.0 → 1.20 at `k_var`), or the units change and the loss change are confounded. Even
done correctly it buys nothing over (c), because the model learns the same function up to scale.

### (b) Refit the clock. **REJECTED, and not open for discussion.**

`fleischer_clock.json` is a frozen external artefact with published provenance. Reweighting it to
improve our numbers is fitting the test. Standing rule since the Stage 1.5 fix plan.

### (c) Calibrate at the REPORTING boundary. **RECOMMENDED.**

Apply `k` where ΔAge is reported — `mu_age`, the conformal half-width `q`, and any ΔAge quoted in
years — leaving the target, the training, and every stored artefact untouched.

**Why this is the right shape:**

- **No rebuild, no retrain, no re-score.** The guard record is not restarted.
- It achieves the entire actual goal: absolute ΔAge claims and interval widths become honest.
- It cannot change model skill, and therefore cannot be mistaken for having done so.
- It is trivially reversible — one factor, one place.
- Because ρ is invariant, every ranking-based result stands unchanged and needs no re-audit.

## 14.5 The mandatory guard, computed in advance

A pure rescale must produce **exactly** these, off `scorecard/c7_A_keep_hff.json` at `k_LS`
(recomputed at whichever `k` ships):

| fold | ΔAge MAE | level shift | conformal width |
|---|---|---|---|
| N3 | 7.74 | −8.02 | 25.69 |
| O1 | 4.63 | +0.53 | 25.06 |
| O2 | 5.01 | +3.57 | 25.07 |
| Y1 | 3.72 | +0.25 | 29.32 |
| Y2 | 7.47 | −5.12 | 22.92 |

**`rank_model_dage` must be EXACTLY unchanged.** Any deviation from this table means the change
did something other than convert units — investigate before keeping. This is the same shape of
guard as Stage 2's "a level shift is rank-invariant by construction".

## 14.6 How the result must be reported, whatever it is

**Calibrated ΔAge is reported ALONGSIDE raw, never silently in place of it**, and always carrying
its provenance: *k estimated on 3 donors of the transient arm against methylation truth; transfer
to the Sendai and HFF cohorts is UNTESTED.*

This is not excessive caution — it is the only honest form available. The cohort carrying
methylation truth and the cohort the model trains on are disjoint, and nothing in this project has
established that a scale factor crosses between them.

**A 63 % drop in ΔAge MAE must never be reported as an improvement.** It is a change of units. The
`compare` table will show it as a large ACCEPT, and that verdict is meaningless here.

## 14.7 What would license the transfer

Recorded so the open question has an owner rather than being left as a caveat:

1. **A Sendai-arm condition carrying methylation truth.** None exists in the current data; this is
   a data-acquisition requirement, not an analysis one (cf.
   `plans/DATA_REQUIREMENT_SECOND_TIMECOURSE.md`).
2. **Failing that**: estimate `k` independently within each variant and check the factors agree
   across cohorts that *are* shared. Weaker, and it does not establish transfer — it only fails to
   refute it.

## 14.8 Verification, when the Change is executed

| item | how |
|---|---|
| the factor | `k_var` from `results/diag_stage11_scale_results.json`, not a fresh fit |
| units only | the §14.5 table reproduced exactly; `rank_model_dage` bit-identical |
| nothing retrained | no fold rebuilt, no bundle rewritten, `scalers.json` untouched |
| reversibility | removing the factor restores every number exactly |
| provenance travels | the transfer caveat appears wherever a calibrated number is printed |
| record | `CHANGES.md`, stating plainly that the MAE drop is units, not skill |
