# Stage 17 — the scorecard rewards over-coverage, and hides fold agreement

**Status:** PLAN, then EXECUTE. Changes `scorecard.py` only (metric direction, aggregation,
printing). **No model change, no rebuild, no snapshot rewritten.**

Found while executing Stage 12 §12.9. I deliberately did **not** patch either defect mid-Change —
changing the instrument during a measurement invalidates the measurement. This is that Change.

---

## 17.1 Defect D1 — `conformal_coverage` is judged "higher is better"

`METRICS["conformal_coverage"] = ("higher", ...)`. But coverage is **target-seeking**: it should
approach `conformal_level` (0.90, stored per fold), not climb. Coverage 1.000 is not a triumph —
it means the intervals are too wide, which is exactly why `conformal_width` sits at **63–81 years**.

**The harm is not hypothetical.** Across the distinct snapshots, most folds are *over*-covering:

| snapshot | mean coverage | mean \|cov − 0.90\| | folds over-covering |
|---|---|---|---|
| A_xdonor | 0.889 | 0.178 | **5/6** |
| gc2_B_mask_hff | 0.849 | 0.170 | **5/6** |
| c7_A_keep_hff | 0.923 | 0.098 | **4/5** |
| c7t_stage12 | 0.914 | 0.107 | **4/5** |
| gc2_C_shuffle_hff_s0 | 0.904 | 0.096 | 4/6 |
| gc2_D_stratshuffle | 0.896 | 0.120 | 4/6 |
| baseline | 0.401 | 0.499 | 0/6 |

Mean coverage looks respectable (0.85–0.92) while the mean *distance from nominal* is 0.10–0.18 —
the folds scatter widely and mostly **above** the target. Under "higher is better", **a change that
simply widened every interval until nothing escaped would score ACCEPT.**

**Measured effect:** 4 of 21 distinct coverage verdicts change under the corrected rule.

## 17.2 Defect D2 — the table never shows whether the folds AGREE

`scorecard.py`'s own header says:

> *"A uniform change is caught at any size; one that helps some folds and hurts others can be large
> in the mean and still read as noise — check the per-fold column before trusting an aggregate
> verdict."*

**But there is no per-fold column.** Only `dage_mae_model` gets a per-fold table; the other 17
metrics print a mean and a CI and nothing else. The instrument asks the reader to perform a check
it does not make available.

**This already cost a correct reading.** In the Stage 12 comparison:

| metric | mean diff | 95% CI | verdict | folds better/worse |
|---|---|---|---|---|
| `conformal_width` | **−6.97 yr** | [−19.93, +5.98] | noise | **5 / 0 — unanimous** |

Every fold moved the same way and the CI still spanned zero, because the paired t is driven by the
*consistency of magnitude*, not of direction. I recorded that result as "directional, not
significant, not claimed", which is defensible but understated: **5/5 unanimous (sign test
p = 0.062) is a near miss, not a shrug.** The Stage 12 record is corrected accordingly (§17.6).

## 17.3 The fix

**F1 — a `"target"` direction.** `conformal_coverage` is judged on `|coverage − conformal_level|`
per fold, aggregated as `mean(|·|)` and paired on the per-fold gaps, lower-is-better. The target is
read **per fold** from `conformal_level`, never hard-coded — it is data, and a future run may
change it.

**F2 — the signed context row is kept.** Exactly as Stage 13 did for `level_shift`: an `"abs"`-style
metric also prints its signed mean, marked `(context, never judged)`. Here it answers a genuinely
different question — *are we over- or under-covering?* — which the magnitude alone destroys.

**F3 — a fold-direction tally on every row.** Each row gains `better/worse/same` counts across the
paired folds, so a unanimous-but-heterogeneous change is visible without hunting.

## 17.4 The guard that keeps F3 honest

A direction tally is a hair's breadth from a second significance test, and a second test on the same
data is how a project talks itself into findings the first test rejected.

**The tally is DESCRIPTIVE and never produces a verdict.** The accept/reject rule stays exactly what
it was: the paired 95 % CI. The tally sits beside it as context, in the same class as the
`(context)` rows.

Two further reasons it cannot be used as a rescue:

- **The sign test cannot reach significance at these n.** With 5 paired folds the smallest
  achievable two-sided p is 2/2⁵ = **0.0625**; with 6 it is 2/2⁶ = 0.031. A unanimous 5/5 is
  *by construction* p = 0.0625 and can never clear 0.05. Printing the p-value would invite
  exactly the misreading it is meant to prevent, so **the tally prints counts, not a p-value.**
- **It is symmetric.** It flags unanimous *regressions* as loudly as unanimous improvements.

## 17.5 What this stage does NOT do

- **No `conformal_width` direction change.** Narrower is genuinely better *at equal coverage*, and
  coverage is separately judged. Left as `("lower", ...)`; noted, not changed.
- **No model change, no rebuild, no re-snapshot.** Aggregation and printing only.
- **No stored snapshot modified.** All ten remain readable and are re-judged in place.
- **No past record rewritten.** Corrections are appended beside the originals.
- **The decision rule is unchanged.** Still the paired 95 % CI. F3 adds context, not authority.

## 17.6 Correction owed to the Stage 12 record

The 2026-08-18 entry recorded `conformal_width` −6.97 yr as *"directional, NOT significant, not
claimed"*. That stands as written, and is **appended to**: the change was **unanimous across all
five folds (5/0)**, which the table could not show at the time. It remains not significant by the
pre-registered rule, and Stage 12's verdict — the pre-registered null — is **unchanged**. What
changes is only how strong the unclaimed directional hint was.

## 17.7 Verification

| item | how |
|---|---|
| D1 — target direction | test: coverage 1.00 and 0.80 score equally against a 0.90 target |
| D1 — per-fold target | test: the target is read from `conformal_level`, not a literal |
| D1 — the real flips | test reproducing the 4 flips on the committed snapshots |
| D1 — over-coverage no longer rewarded | test: widening every interval to coverage 1.0 must NOT read ACCEPT |
| D2 — the tally | test: a unanimous change reports n/0 while its CI still spans zero |
| the guard | test: the tally never appears in a verdict string, and no p-value is printed |
| symmetry | test: a unanimous regression is flagged as loudly as a unanimous improvement |
| nothing else moved | full `pytest` green, `ruff` clean; non-target metrics' verdicts unchanged |
| record | `CHANGES.md` + the §17.6 correction, originals left intact |
