# Stage 18 — is the fate head predicting biology, or reading a clock?

**Status:** 🟩 **RETROSPECTIVE RECORD, written 2026-08-20. EXECUTED 2026-08-18, commit `b60ee0f`.**

> ⚠️ **This is not a plan and must not be read as one.** Stage 18 ran without a plan file: its
> design and its decisive statistic were fixed in the docstring of
> `experiments/diag_stage18_fate_beyond_day.py` before the numbers, and the result went straight to
> `CHANGES.md`. This document exists so `plans/` is not silently missing a stage.

---

## The question

`fate_prauc` ≈ 0.96 is the strongest number the project has. But **`dose_time` is a model input**
and it encodes the timepoint, and the record already warned that "day > 31 → unsafe" gets most of
the way there. How much of that 0.96 survives once the clock is taken away?

## The structural fact found first

On the held-out donors, **fate is very nearly a function of timepoint**:

| fold | timepoints | carrying >1 class | time-only PR-AUC | model PR-AUC |
|---|---|---|---|---|
| N2 | 11 | **0** | n/a (19/19 safe) | n/a |
| N3 | 12 | 1 | 0.892 | 1.000 |
| **O1** | 12 | **0** | **1.000** | 0.994 |
| **O2** | 12 | **0** | **1.000** | 1.000 |
| Y1 | 11 | **5** | 0.660 | 0.796 |
| Y2 | 12 | 1 | 0.988 | 1.000 |

**Only 7 of 70 timepoints carry more than one class.** On O1 and O2 the timepoint *alone* reaches
PR-AUC **1.000** — a lookup table on the hour, using no genes, is unbeatable there.

## The decisive statistic, fixed before the run

Any marginal metric is inflated by an input the model was handed. The question is only asked
**within a timepoint**, where `dose_time` is constant and cannot help. Over every (safe, unsafe)
pair drawn from the same timepoint **of the same donor**:

> **stratified AUC 0.917 over 12 pairs, permutation p = 0.0091**

The null shuffles scores **within strata**, preserving stratum sizes. A global shuffle would also
destroy the between-timepoint structure the model is not being credited for, and would be trivially
beaten. Per fold: N3 1.000 (4 pairs), Y1 0.857 (7), Y2 1.000 (1). N2/O1/O2 contribute **zero**.

## The verdict, both halves

**The marginal 0.93–0.96 is very largely the clock** — measured, not suspected.

**But there IS signal underneath.** 11 of 12 same-timepoint pairs ranked correctly is significant
against a proper null. The fate head is not only reading a calendar.

**And the entire evidence base for that is 12 pairs, from 7 strata, across 3 donors.** A real
result on a very thin base, and the honest phrasing is exactly that.

## The fold carrying the evidence is the one that looks worst

Y1 supplies **7 of the 12 pairs** and has **5 of the 7 mixed timepoints** — and it has the *lowest*
marginal PR-AUC (0.796 against 1.000 elsewhere). The same fact seen twice: Y1 is the only donor
whose fate does not track its clock, so it is simultaneously the hardest fold for a clock-reader
and the only fold that can test one.

**Consequence for the LOOCV design:** aggregate `fate_prauc` is dominated by folds that cannot
discriminate the hypothesis. A future evaluation should weight, or at least report, folds by how
much within-timepoint contrast they contain.

## What it does NOT say

It does not say the fate head is worthless — it says its headline metric is inflated by an input it
was handed. It does not settle whether the within-timepoint signal generalises: **12 pairs cannot.**
And it does not touch ΔAge, which was already established as circular for a different reason
(`diag_clock_circularity`, ρ 0.96–0.99).

## The open requirement this creates

**More mixed-timepoint data is the only thing that grows those 12 pairs.** This is a
data-acquisition requirement, not an analysis one — see
`plans/DATA_REQUIREMENT_SECOND_TIMECOURSE.md`.

## Artefacts

`experiments/diag_stage18_fate_beyond_day.py` ·
`results/diag_stage18_fate_beyond_day_results_s16.json` ·
`tests/test_diag_stage18_fate_beyond_day.py` (17 tests) · `CHANGES.md` 2026-08-18.
