# STEP 6 / G-c step 2 — full report

**CellFate-Rx · Stage 1.5.3 · runs of 2026-08-02 and 2026-08-03**

> **The question.** HFF contributes **33 613** of the project's **33 688** ΔAge training labels —
> 99.7 %. Stage 1.5.2's G-c step 1 found those labels carry an artefact signature. So: *do they help
> the age head, or does it learn artefact from them while the 75 usable Gill labels are drowned out?*
>
> **The answer, after two full 6-fold LOOCV runs.** **Inconclusive, and not for want of trying.** The
> first run was invalid — the treatment silently redefined the target variable. That was found,
> fixed, and verified. The second run is a clean one-change comparison and returns **+0.661 yr
> (95 % CI [−4.384, +5.707])** — an effect of 4.6 % of baseline that the design cannot resolve.
> **Nothing licenses discarding the labels.**
>
> The binding constraint is no longer a bug. It is that **two of the six donors sit outside the
> clock's validated age range and carry 4.3× the variance of the other four.** No re-run fixes that.

---

## Contents

1. [What had to be fixed before it could run at all](#1-what-had-to-be-fixed-before-it-could-run-at-all)
2. [Run 1 — the confounded result](#2-run-1--the-confounded-result)
3. [The three confounds](#3-the-three-confounds)
4. [C-I — the fix, and how it was verified](#4-c-i--the-fix-and-how-it-was-verified)
5. [Run 2 — the de-confounded result](#5-run-2--the-de-confounded-result)
6. [Power, honestly](#6-power-honestly)
7. [Guards and secondary metrics](#7-guards-and-secondary-metrics)
8. [What is and is not licensed](#8-what-is-and-is-not-licensed)
9. [Open anomalies](#9-open-anomalies)
10. [Options forward](#10-options-forward)
11. [Reproduction](#11-reproduction)

---

## 1. What had to be fixed before it could run at all

A pre-flight was run instead of starting the compute. It found three blockers, the first of which
would not have failed — it would have returned a plausible answer.

| | blocker | why it mattered |
|---|---|---|
| **B-1** | `retrain_stage1.py` **reuses existing shards** and has no build step, but `age_mask` is written at *build* time | Both arms would have trained on identical data. Measured directly: with the constant set either way, N2's shards read **127 815 / 127 815** age-valid. The paired CI would have included 0, and the pre-registered table reads that as *"the labels contribute nothing → mask them"* — **discarding 99.7 % of the age labels on a run where the treatment was never applied.** |
| **B-2** | `age_window_k` was not plumbed into any driver's `TrainConfig` | C-5 Option 2 would have run at `k = 1` = OFF regardless of what the plan said. |
| **B-3** | arm B is not a mask flip — `_deconfound_train_only` moves `y_age` itself | Foreshadowed C-I below. |

**Fixes:** the arm switch set before the build (`aging.py:304` reads it during `delta_age`);
`age_window_k` / `age_window_max_batches` passed into `TrainConfig`; and a **B-1 guard placed before
training** that asserts the arm took effect and writes `step6_arm_census.json`.

> A correction on the record: the pre-flight also claimed *"no raw GEO input on this machine."* That
> was wrong — it searched inside the repo (`find . -maxdepth 2`) when the pipeline's defaults are
> `D:\GSE242423` and `D:\Gill`, outside it. The data was there all along.

### The contrast proof

Both arms were then built at smoke scale with identical geometry, differing only in the mask:

| arm | age-valid / train | | `k` |
|---|---:|---:|---:|
| A (control) | **5 718 / 5 718** | 100.00 % | 4 |
| B (treatment) | **78 / 5 718** | 1.36 % | 4 |

Same train-cell count; 98.64 % of labels removed. Both guard branches executed and passed.

---

## 2. Run 1 — the confounded result

Both arms rebuilt from raw GEO and retrained, 6 folds each. The B-1 guard passed on all 12 folds:
arm A **33 688 / 33 688** age-valid, arm B **75 / 33 688** — 75 exactly, the plan's E11 prediction,
reached independently by a fresh build.

**Arm A reproduced `scorecard/baseline.json` to 3 decimals on all six folds**, confirming bar A1 in
production: `age_window_k = 4` costs the control arm nothing.

| fold | donor age | arm A | arm B | Δ |
|---|---:|---:|---:|---:|
| N2 | 0 | 21.79 | 43.17 | **+21.37** |
| N3 | 0 | 29.69 | 22.48 | −7.22 |
| O1 | 53 | 5.39 | 9.02 | +3.63 |
| O2 | 53 | 7.54 | 9.37 | +1.84 |
| Y1 | 29 | 7.28 | 12.63 | +5.35 |
| Y2 | 35 | 14.06 | 12.91 | −1.15 |

**effect +3.971 yr · SD 9.599 · MDE 10.074 · CI [−6.102, +14.045] · power for Δ\* 11.3 %**

Inconclusive — and dropping the N2 outlier did not rescue it (n=5: effect +0.491, SD 4.933,
MDE 6.125, still inconclusive). But the deeper problem was not power.

---

## 3. The three confounds

### C-I — the ΔAge target itself moved *(a bug; fixed)*

`_deconfound_train_only` fitted the cell-cycle deconfounder, **and re-centred on controls**, using
`age_mask` — which `AGE_MASKED_DATASETS` changes. Masking HFF dropped the fit from **33 613 single
cells to 75 bulk samples**, and the refitted transform was applied to *every* cell's `y_age`,
including the held-out evaluation targets.

| fold | arm A (slope, intercept) | arm B (slope, intercept) |
|---|---|---|
| N2 | −3.93, −3.42 | **−24.20, +10.12** |
| O1 | −9.27, −10.09 | **−24.80, +6.62** |
| Y2 | −9.66, −9.92 | **−24.88, +3.06** |

Slope 2.5–6× steeper, intercept sign-flipped in all three. **Arm B was not "arm A minus HFF's
labels" — it was a different target variable.** Two changes at once: a one-change-rule violation in
the *design*, not the execution.

The tell was that `dage_mae_ridge` regressed **+9.21 yr** even though ridge never touches the trained
age head.

**This was predicted in writing before the run** — the plan's C-5 "second consequence."

### C-II — the label pool shifts into the clock's extrapolation zone *(not a bug; unfixed)*

`fleischer_clock.json` declares `age_range = [1.0, 96.0]`. **N2 and N3 are donor age 0 — outside it.**
Their ΔAge labels are extrapolations.

| | share of age labels from out-of-range donors |
|---|---|
| arm A | **0.09 %** (30 of 33 688) |
| arm B | **40 %** (30 of 75 — N2 14, N3 16) |

Masking HFF does not merely *reduce* the labels; it **up-weights out-of-clock-range donors ~400×**.
This is a genuine property of the treatment, not an implementation error — which is why it survives
into run 2 and becomes the binding constraint.

### C-III — ridge is not a control *(a reasoning error of mine, corrected)*

I initially framed ridge as isolating the target change. It does not: `scorecard.py:95` fits it on
`tr.y_age[tr.mask]`, the **same masked labels**, so it suffers both changes too. What it does show
is that under identical damage the neural model degraded *less* than the linear baseline (mean excess
−5.24 yr) — suggestive, but the difference-in-differences CI [−12.58, +2.10] includes 0.

---

## 4. C-I — the fix, and how it was verified

The fix separates two questions that had been conflated:

| | |
|---|---|
| `deconfound_mask` | *is this cell's ΔAge **computable**?* → decides the **value** of `y_age` |
| `age_mask` | *may the age head **train** on it?* → decides which cells the loss uses |

**Three call sites had to change, not one:**

1. the deconfounder coefficient fit
2. pass 2's **control re-centring**
3. `if cfg.deconfound and age_mask.any()` — an all-HFF chunk has no trainable label under arm B, so
   it was dropped from the sidecar cache **entirely** and never reached the deconfounder at all

`deconfound_mask` is recomputed by calling `age_label_policy` with `masked_datasets=frozenset()`,
**not** by inverting `age_mask_reason`: the policy records only the *first* reason that fires, so a
cell both out-of-range and dataset-masked reads `"dataset_policy"` and inverting would wrongly
readmit it.

### Verification — three independent checks

| check | result |
|---|---|
| deconfounder coef, both arms, all 6 folds | **identical** (N2 `−3.9289, −3.4207` in both; was `−3.93,−3.42` vs `−24.20,+10.12`) |
| `y_age` across arms, **row-exact** over 7 062 rows | **`max|Δ| = 0.000e+00`**, 0 NaN-pattern mismatches, 6 938 rows differing in `age_mask` and nothing else |
| arm A rerun vs pre-C-I arm A, all 6 folds | **`max|Δ| = 0.000e+00`** — C-I is a no-op on the control, as predicted |

> The pre-flight's own gate had a defect worth recording: the first version keyed a dict on
> `cell_id` and compared **1 024 of 7 062 rows** without saying so — `cell_id` is not unique across
> the dataset. It reported "1,024 cells" and passed, reading like full coverage. A 15 % sample is
> not a gate for a 10 h run; it now compares positionally per shard.

**Downstream audit.** C-I changes shard content — withheld cells now carry a finite `y_age` instead
of NaN. `schemas.py:142` enforces `age_mask=False ⟹ y_age None/NaN` but never fires
(`assemble.py:49` passes `None` when masked; there is no shard→`Sample` read path). All twelve
`y_age` consumers gate on the mask. Nothing reads a withheld label.

---

## 5. Run 2 — the de-confounded result

Both arms rebuilt and retrained, 6 folds each, `age_window_k = 4`, **arm-suffixed roots** so both
survive (`CELLFATE_FOLD_SUFFIX`, honoured by `scorecard.py` too).

| fold | arm A | arm B | Δ | *(run 1 Δ)* |
|---|---:|---:|---:|---:|
| N2 | 21.79 | 29.55 | **+7.76** | *+21.37* |
| N3 | 29.69 | 22.53 | −7.16 | *−7.22* |
| O1 | 5.39 | 7.61 | +2.22 | *+3.63* |
| O2 | 7.54 | 7.36 | −0.17 | *+1.84* |
| Y1 | 7.28 | 8.59 | +1.31 | *+5.35* |
| Y2 | 14.06 | 14.07 | +0.02 | *−1.15* |

### Primary — all 6 folds

| quantity | run 2 | run 1 |
|---|---|---|
| **effect** | **+0.661 yr** | +3.971 |
| **observed SD** | **4.808 yr** | 9.599 |
| **MDE** (1.049 × SD) | **5.045 yr** | 10.074 |
| **95 % CI** | **[−4.384, +5.707]** | [−6.102, +14.045] |
| power for Δ\* = 3.572 | 31.5 % | 11.3 % |

**The SD halved and the effect shrank six-fold, to 4.6 % of the 14.29 yr baseline.** Most of run 1's
+3.97 was the confound, not the labels.

**Verdict: CI includes 0 and MDE 5.045 > Δ\* 3.572 → INCONCLUSIVE.**

### Secondary — the pre-registered 4 in-range folds

Registered **before arm B ran**, on C-II grounds. Folds O1, O2, Y1, Y2.

| | |
|---|---|
| per-fold Δ | `+2.22, −0.17, +1.31, +0.02` |
| effect | **+0.843 yr** |
| **observed SD** | **1.130 yr** |
| MDE (1.591 × SD) | 1.799 yr |
| 95 % CI | **[−0.956, +2.642]** |

### 🔬 C-II is now measured, not hypothesised

**Dropping N2 and N3 collapses the SD 4.808 → 1.130 — a factor of 4.3, on 2 of 6 folds.** The two
donors outside the clock's validated range carry almost all the fold-to-fold variance.

---

## 6. Power, honestly

At face value the secondary's **MDE 1.799 ≤ Δ\* 3.572**, which under the registered outcome table
reads *"the labels are genuinely not contributing → mask them."* **That reading is not taken here.**

With n = 4 the SD is itself a noisy estimate. The χ² 95 % interval on the true σ:

| | n | SD observed | 95 % CI on σ | MDE at σ_high | vs Δ\* = 3.572 |
|---|---:|---:|---:|---:|---|
| run 1, all folds | 6 | 9.599 | [5.992, 23.543] | 24.706 | ≫ |
| **run 2, all folds** | 6 | 4.808 | [3.001, 11.792] | **12.375** | ≫ |
| **run 2, in-range** | 4 | 1.130 | [0.640, **4.215**] | **6.707** | **>** |

**Neither analysis is robustly powered once σ is admitted to be an estimate.** The secondary was
pre-registered as "underpowered by construction"; it turned out *better* powered than the primary,
and still not enough.

This is the single most important line in the report: **the temptation was to take the secondary's
face-value MDE and declare the labels expendable. The uncertainty on σ does not support it.**

---

## 7. Guards and secondary metrics

| metric | arm A | arm B | mean diff | 95 % CI | verdict | *(run 1)* |
|---|---:|---:|---:|---|---|---|
| `rank_model_dage` | 0.948 | 0.879 | −0.069 | [−0.100, −0.037] | **REGRESSION** | *−0.186* |
| `rank_ridge_dage` | 0.955 | 0.891 | −0.064 | [−0.110, −0.018] | **REGRESSION** | *−0.146* |
| `dage_mae_ridge` | 14.05 | 17.54 | +3.485 | [−0.130, +7.099] | noise | *REGRESSION +9.21* |
| `interval_width` | 65.89 | 71.99 | +6.104 | [−8.477, +20.684] | noise | *REGRESSION* |
| `conformal_coverage` | 0.889 | 0.849 | −0.040 | [−0.485, +0.406] | noise | noise |
| `fate_prauc` | 0.992 | 0.978 | −0.013 | [−0.043, +0.017] | noise ✅ | noise ✅ |
| `fate_roc` | 0.983 | 0.966 | −0.017 | [−0.051, +0.018] | noise ✅ | noise ✅ |
| `fate_ece` | 0.249 | 0.320 | +0.071 | [−0.019, +0.161] | noise ✅ | noise ✅ |
| `fate_ece` (Platt) | 0.140 | 0.236 | **+0.096** | **[+0.011, +0.182]** | **REGRESSION** ⚠️ | *REGRESSION* |
| `ood_flag_rate` | 0.273 | 0.516 | **+0.243** | **[+0.068, +0.419]** | context | *+0.162* |

**The three registered fate guards hold.** The plan requires `fate_prauc`, `fate_roc` and `fate_ece`
not to move, and they do not.

**Ranking is the one consistent cost.** `rank_model_dage` −0.069 with the CI excluding 0 — and ridge
degrades by nearly the same amount, so this is two learners both ranking worse on 75 labels than on
33 688. Consistent, and no longer confounded by a moving target.

---

## 8. What is and is not licensed

### ❌ Not licensed

- **Discarding HFF's 33 613 labels.** Neither the primary nor the secondary is robustly powered for
  Δ\*. A null that could not have detected the effect is not evidence of absence.
- **Any of the three branches of the pre-registered outcome table**, as written — the first run
  because it was confounded, the second because MDE > Δ\* once σ uncertainty is admitted.

### 📌 The bounded estimate — the reportable result of step 6

Reported as an **estimate with limits**, not a verdict.

> **Masking HFF's 33 613 ΔAge labels changes `dage_mae_model` by**
> ### +0.661 yr · 95 % CI [−4.384, +5.707]
> **observed SD 4.808 · MDE 5.045 · Δ\* 3.572 · n = 6 paired donor folds**
>
> The point estimate is **4.6 % of the 14.29 yr baseline**. The interval is wide enough to contain
> both a meaningful improvement and a meaningful regression, so **the sign is not established**.
>
> **MDE 5.045 > Δ\* 3.572: the design could not have detected the smallest effect worth acting on.**
> This therefore **licenses nothing about discarding HFF's labels** — in either direction. It is a
> bound on the effect's size, not a measurement of its value.
>
> **The confounded first run's +3.971 yr was mostly the confound, not the labels.** After C-I the
> estimate fell six-fold and the SD halved (9.599 → 4.808). Any use of the earlier figure should be
> replaced by this one.
>
> *Carry forward as:* `dage_mae_model` effect of masking HFF = **+0.66 yr [−4.38, +5.71]**,
> undetermined at this geometry.

### ✅ Supported

- **The effect on `dage_mae_model` is small.** Point estimate +0.661 yr, 4.6 % of baseline. Run 1's
  +3.97 was mostly the confound.
- **C-I is fixed and verified**, three independent ways.
- **C-5 Option 2 at `k = 4` costs the control arm nothing** — arm A bit-identical across two
  independent full rebuilds three weeks apart.
- **C-II is real and quantified**: 2 of 6 donors carry 4.3× the variance.
- **Ranking degrades detectably** when the labels are removed (−0.069, CI excludes 0).

---

## 9. Open anomalies

> ### 🔵 Resolved 2026-08-03 — items 1 and 2 were **not** anomalies. One mechanism explains four.
>
> I flagged `fate_ece` (Platt) as a standing anomaly on the grounds that *"the fate head consumes no
> ΔAge."* That is true of the **labels** and false of the **network**. `models/network.py:60-62`:
>
> ```python
> self.trunk = nn.Sequential(_mlp_block(d_cell + d_u, latent_dim, p_drop),
>                            _mlp_block(latent_dim, latent_dim, p_drop))
> ...
> z = self.trunk(torch.cat([self.cell(x), self.pert(u, dose_time)], dim=1))
> return self.cls_head(z), self.age_head(z), z
> ```
>
> **One shared trunk feeds both heads**, and `MultiTaskLoss` sums their losses. Masking 99.7 % of the
> age labels changes the age loss → changes trunk gradients → changes the representation `z` that the
> fate head reads. So **`rank_model_dage`, `rank_ridge_dage`, `fate_ece` (Platt) and `ood_flag_rate`
> moving together is one mechanism, not four separate puzzles**: a degraded shared representation.
>
> This matters beyond tidiness. It gives the ranking result a **causal account** rather than a
> correlation, which is what makes option 3 worth considering at all — and it explains why the three
> registered fate guards (`prauc`, `roc`, raw `ece`) can hold while the *calibrated* ECE moves: the
> calibration path is more sensitive to representation drift than the ranking metrics are.

1. ~~**`fate_ece` (Platt) regressed in both runs** — unexplained.~~ **Explained above**: shared-trunk
   coupling. It remains worth watching as a *magnitude* — +0.096 [+0.011, +0.182] — but it is no
   longer a puzzle about how ΔAge could possibly reach the fate head.
2. ~~**`ood_flag_rate` nearly doubled** — not established.~~ **Same mechanism.** 0.273 → 0.516 is
   what a drifted shared representation predicts.
3. **`scorecard.py`'s `level shift` row prints the mean without its sign.** In run 1 it read `5.713`
   when the value was **−5.713**; signed means moved −5.713 → +2.267 (looks better) while magnitudes
   moved 13.12 → 18.66 (worse). **A reader trusting that row draws the opposite conclusion.**
   Not yet fixed.
4. **Snapshot tags collide across runs.** Run 2 reused `gc2_A_keep_hff` and overwrote run 1's
   snapshot; it had to be recovered from git. The fold *roots* are now suffixed, but the tags are not.

---

## 10. Options forward

**Another 10 h of this design will not help.** The binding limit is six paired folds with two outside
the clock's validated range, and no re-run touches that.

| | option | cost | what it buys |
|---|---|---|---|
| **1** | **Accept and report.** State that step 6 cannot resolve Δ\* at this geometry; carry +0.661 [−4.384, +5.707] forward as a bounded estimate, not a verdict. | none | honesty, and no further spend |
| **2** | **Fix C-II at source** — **more donors, and/or a clock validated at age 0.** See the correction below: a 4-donor re-run is *not* among the options. | new data | removes the dominant variance source |
| **3** | **Change the estimand to ranking.** `rank_model_dage` shows a consistent, detectable effect where MAE does not. If ranking is what Stage 2 needs, register it as primary. | a **new** pre-registration | an estimand this geometry can actually resolve |

> ### 🔵 Correction to option 2, 2026-08-03
>
> The earlier wording offered *"a 4-donor study, accepting the power that implies."* **That option
> does not exist, and this run already measured why.** The pre-registered in-range secondary **is**
> the 4-donor design (O1, O2, Y1, Y2), and it was not powered either: σ's 95 % upper bound is
> **4.213**, giving **MDE 6.704 > Δ\* 3.572**. Dropping the neonatal donors shrinks the *point*
> estimate of the SD but loses two folds, and the two effects cancel.
>
> **Option 2 therefore reduces entirely to: more donors, and/or a clock validated at age 0.**
> Nobody should read a 4-donor re-run as viable — it has been run.

Option 3 is the most promising on the evidence, and is explicitly **not** a re-read of this run — the
ranking result here was a secondary metric, and promoting it after seeing it would be exactly the
post-hoc selection the pre-registration exists to prevent.

---

## 11. Reproduction

```bash
python plan_tests/step6_preflight.py          # gates the run; exits non-zero on any failure
./run_step6_arm.sh A                          # control,   snapshot gc2_A_keep_hff
./run_step6_arm.sh B                          # treatment, snapshot gc2_B_mask_hff
python scorecard.py compare gc2_A_keep_hff gc2_B_mask_hff
```

| | |
|---|---|
| cost | ~5 h per arm — 6 folds × (~25 min build + ~35 min train), 6 ensembles/fold, GTX 1080 |
| artefacts | `scorecard/gc2_{A,B}_*.json`, `cellfate_loocv_<donor>_arm{A,B}/`, `results/step6_preflight_results.json` |
| bars | `plan_tests/register_gc_step2_bar.py`, rows in `tests/test_bars_resolvable.py` |
| invariants | `tests/test_ci_deconfounder_arm_invariance.py` (6 tests, incl. a mutation check) |
| raw record | `CHANGES.md`, and `results/STEP6_REPORT.md` for the two runs as they were written |

---
---

# ARM C — the label-permutation control (2026-08-03)

**Purpose.** The A−B ranking gap (−0.0688) was confounded between *(i)* HFF's labels carry
information, and *(ii)* 75 labels is simply too few, whatever they contain. Arm C holds label
**volume** at arm A's level and destroys only the cell↔label **pairing**.

**Validity, pre-flight, all six checks passed on the real build:** 33 688 trainable labels (= arm A
exactly), HFF label multiset identical, **42 481 of 42 481 labels moved**, non-HFF `y_age`
bit-identical to arm A (`max|Δ| = 0.0`), `age_mask` identical everywhere. Shuffle seed **0**,
recorded in every fold's census. The permutation runs after the deconfounder fit, so arm C's
coefficient equals arm A's.

## The result

| `rank_model_dage` | N2 | N3 | O1 | O2 | Y1 | Y2 | **mean** |
|---|---:|---:|---:|---:|---:|---:|---:|
| **A** true labels | 0.9104 | 0.9091 | 0.9896 | 0.9701 | 0.9596 | 0.9468 | **0.9476** |
| **B** no HFF labels | 0.8455 | 0.8753 | 0.9234 | 0.8610 | 0.8614 | 0.9065 | **0.8788** |
| **C** shuffled labels | 0.4909 | 0.5273 | 0.8429 | 0.5455 | 0.4561 | 0.5961 | **0.5765** |

| comparison | effect | SD | 95 % CI | TOST vs ±0.0344 |
|---|---:|---:|---|---|
| **C − A** | **−0.3711** | 0.1213 | [−0.4985, −0.2438] excludes 0 | not equivalent |
| **C − B** | **−0.3024** | 0.1139 | [−0.4219, −0.1829] excludes 0 | not equivalent |
| B − A *(reference)* | −0.0688 | 0.0302 | [−0.1004, −0.0371] | — |

**C sits 540 % of the way from A to B — far outside the interval the pre-registration assumed.**

## ⚠️ Neither pre-registered branch fired

The registered table offered *"C ranks like A"* → volume effect, or *"C ranks like B"* → labels
informative. **Both assumed C ∈ [A, B]. It is not.** Shuffled labels are **dramatically worse than no
labels at all**. What follows is therefore an **inference beyond the registered outcomes** and is
labelled as such — it is not a pre-registered conclusion.

The equivalence bar was never usable either: achieved **SD(C−A) = 0.1213** against the ≤ 0.020 the
bar needed. Moot, given effect sizes 10× the margin, but recorded — the equivalence branch was
unavailable, exactly the conditional the registration warned about.

## It is not a broken run — the damage is age-specific

| | A | C |
|---|---:|---:|
| `fate_prauc` | 0.9915 | **0.9898** |
| `fate_roc` | 0.983 | 0.979 |
| `conformal_coverage` | 0.889 | 0.904 |
| `dage_mae_model` | 14.29 | **28.93** |
| `rank_ridge_dage` | 0.955 | **0.742** |
| `level_shift_model` | −5.71 | **−20.28** |

Training did not diverge: the **fate head is essentially untouched** (0.9915 → 0.9898) while every
age metric collapses. That is precisely the signature of corrupting age labels and nothing else — and
it also **bounds the shared-trunk coupling**: real, but modest. Ridge, refit on the same shuffled
labels, collapses alongside the model.

## What this establishes

**Explanation (ii) is eliminated.** If HFF's 33 613 labels were uninformative filler whose only
contribution was volume and trunk regularisation, permuting them would have left ranking near arm A.
Instead ranking fell by **−0.371, 5.4× the entire A−B gap**. The labels carry **structure the model
actively exploits**; the A−B gap is not a volume artefact.

### 🔴 But this does NOT establish that the labels are *correct*

Arm C destroys **any** consistent structure — real biological signal and **systematic artefact
alike**. A consistent artefact is exactly as shuffleable as a true signal, and G-c step 1 already
found HFF's labels carry an artefact signature.

So the honest statement is: **HFF's labels contain consistent, exploitable structure — but arm C
cannot say whether that structure is biology or a reproducible artefact.** It narrows the original
G-c question without closing it.

## Consequences

1. **Discarding HFF's labels is now harder to justify, not easier.** They are demonstrably load-
   bearing for the model, whatever their provenance.
2. **Option 3 (promote ranking to primary) has a mechanism and a well-powered channel** — the arm C
   effect is 12× the ranking MDE. It remains a **new** pre-registration, not a re-read.
3. **A follow-up worth registering:** shuffle *within* donor/timepoint strata rather than globally.
   Global shuffling destroys both biological signal and artefact; a stratified shuffle that preserves
   the artefact's structure while destroying the within-stratum pairing would separate them. That is
   the experiment that could close G-c.
