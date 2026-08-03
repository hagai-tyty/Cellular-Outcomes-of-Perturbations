# STEP 6 / G-c step 2 — full report

**Run date** 2026-08-02. Both arms rebuilt from raw GEO and retrained, 6 LOOCV folds each,
`age_window_k = 4` in both. Snapshots: `scorecard/gc2_A_keep_hff.json`, `scorecard/gc2_B_mask_hff.json`.

> **The question:** do HFF's 33,613 ΔAge labels help the age head, or does it learn artefact from
> them while 75 usable labels are drowned out?
>
> **The answer this run supports:** none. The primary metric is **INCONCLUSIVE**, and the design
> turns out to be **confounded three ways**. The run is informative — but about the experiment, not
> about the hypothesis.

---

## 1. The pre-flight held

| | arm A (control) | arm B (treatment) |
|---|---|---|
| `AGE_MASKED_DATASETS` | `frozenset()` | `{"hff_sc"}` |
| age-valid / train cells (every fold) | **33 688 / 33 688** | **75 / 33 688** |
| `age_window_k` | 4 | 4 |

**75 exactly** — the plan's E11 prediction, reached independently by a fresh build. The B-1 guard
passed on all 12 folds. Arm A reproduced `scorecard/baseline.json` to 3 decimals on all six folds,
confirming **bar A1 in production**: `age_window_k = 4` costs the control arm nothing.

---

## 2. The primary metric — effect, SD and MDE together

`dage_mae_model`, paired B − A across the six donor folds. Lower is better, so **positive = arm B worse**.

| fold | donor age | arm A | arm B | Δ |
|---|---:|---:|---:|---:|
| N2 | 0 | 21.79 | 43.17 | **+21.37** |
| N3 | 0 | 29.69 | 22.48 | −7.22 |
| O1 | 53 | 5.39 | 9.02 | +3.63 |
| O2 | 53 | 7.54 | 9.37 | +1.84 |
| Y1 | 29 | 7.28 | 12.63 | +5.35 |
| Y2 | 35 | 14.06 | 12.91 | −1.15 |

| quantity | value |
|---|---|
| **observed effect** | **+3.971 yr** |
| **observed SD of the per-fold difference** | **9.599 yr** |
| **MDE** = `t(.975,5)/√6 × SD` = 1.0494 × SD | **10.074 yr** |
| 95 % paired CI | **[−6.102, +14.045]** — includes 0 |
| Δ\* (smallest effect worth acting on) | 3.572 yr |
| **power to detect Δ\* at the observed SD** | **11.3 %** |
| SD needed for ≥95 % power (registered) | **≤ 1.91 yr** |

### The reading, per the rule registered *before* any number existed

**|effect| 3.971 ≤ MDE 10.074 → INCONCLUSIVE.**

`plan_tests/register_gc_step2_bar.py`: *"A null with MDE > Δ\* is INCONCLUSIVE and licenses
nothing."* This is **not** evidence that HFF's labels contribute nothing, and it does **not** license
discarding 99.7 % of the project's age labels.

The observed SD is **5× the 1.91 yr** the design needed. The arms do not track each other — which is
the precise risk the bar was registered to expose, and it could only be measured by running it.

### It is not just the outlier

N2 (+21.37 against a −7…+5 spread) dominates the SD. Dropping it does not rescue the null:

| | n | effect | SD | MDE | verdict |
|---|---:|---:|---:|---:|---|
| all folds | 6 | +3.971 | 9.599 | 10.074 | INCONCLUSIVE |
| drop N2 | 5 | +0.491 | 4.933 | 6.125 | **still INCONCLUSIVE** |

---

## 3. Three confounds — why the outcome table cannot be applied

### C-I. The ΔAge target itself moved (the dominant one)

`_deconfound_train_only` (`build_dataset.py:448`) fits `ΔAge ~ a·cc + b` on **age-valid TRAIN cells**
and re-applies it to *every* shard. Masking HFF drops that fit from **33 613 single cells to 75 bulk
samples**. The coefficient does not merely shift — it changes character:

| fold | arm A (slope, intercept) | arm B (slope, intercept) |
|---|---|---|
| N2 | −3.93, −3.42 | **−24.20, +10.12** |
| O1 | −9.27, −10.09 | **−24.80, +6.62** |
| Y2 | −9.66, −9.92 | **−24.88, +3.06** |

Slope 2.5–6× steeper; intercept sign-flipped in all three. **Arm B is not "arm A minus HFF's
labels" — it is a different target variable**, including for the held-out evaluation cells.

This is C-5's *second consequence*, which the plan predicted in writing and which this run confirms.

### C-II. The label pool's composition shifts toward the clock's extrapolation zone

`fleischer_clock.json` declares `age_range = [1.0, 96.0]`. **N2 and N3 are donor age 0 — outside it.**
Their ΔAge labels are extrapolations.

| | share of age labels from out-of-range donors |
|---|---|
| arm A | **0.09 %** (30 of 33 688) |
| arm B | **40 %** (30 of 75 — N2 14, N3 16) |

Masking HFF does not just *reduce* the labels; it **up-weights out-of-clock-range neonatal donors
from negligible to 40 % of the entire training signal.** Nothing registered this, and it is a
second change riding along with the treatment.

*(Stated as a compositional fact, not a causal claim. That it explains N2's blow-up — MAE 21.79 →
43.17, conformal coverage **1.00 → 0.095** — is a plausible hypothesis this run cannot test.)*

### C-III. The ridge baseline is not a control

`scorecard.py:95` fits ridge on `tr.y_age[tr.mask]` — the **same masked labels**. So ridge suffers
both changes too and cannot decompose target-vs-labels. What it does show:

| fold | Δ model | Δ ridge | excess (model − ridge) |
|---|---:|---:|---:|
| N2 | +21.37 | +23.82 | −2.45 |
| N3 | −7.22 | +11.09 | −18.30 |
| O1 | +3.63 | +3.08 | +0.55 |
| O2 | +1.84 | +2.93 | −1.09 |
| Y1 | +5.35 | +13.21 | −7.86 |
| Y2 | −1.15 | +1.15 | −2.30 |
| **mean** | **+3.97** | **+9.21** | **−5.24** (CI [−12.58, +2.10], includes 0) |

**Ridge degraded on all six folds; the model on four.** Under identical damage the neural model
degrades *less* than the linear baseline. Suggestive of robustness — but the difference-in-differences
CI includes 0, so it is **not** an established result.

---

## 4. Secondary and guard metrics

| metric | arm A | arm B | verdict |
|---|---:|---:|---|
| `rank_model_dage` | 0.948 | 0.761 | **REGRESSION** [−0.358, −0.015] |
| `rank_ridge_dage` | 0.955 | 0.808 | **REGRESSION** [−0.222, −0.070] |
| `dage_mae_ridge` | 14.05 | 23.27 | **REGRESSION** [+0.12, +18.30] |
| `interval_width` | 65.9 | 91.9 | **REGRESSION** [+3.59, +48.53] |
| `fate_prauc` | 0.992 | 0.981 | noise ✅ |
| `fate_roc` | 0.983 | 0.961 | noise ✅ |
| `fate_ece` | 0.249 | 0.326 | noise ✅ |
| `fate_ece` (Platt) | 0.140 | 0.288 | **REGRESSION** [+0.008, +0.288] ⚠️ |
| `conformal_coverage` | 0.889 | 0.833 | noise |

**Fate guards: three of four hold.** The plan requires `fate_prauc`, `fate_roc`, `fate_ece` not to
move, and they do not. But **`fate_ece` (Platt) regressed**, and the plan says a move there is *"a
finding to explain, not a trade-off."* Consistent with the same root cause — `y_age` moving changes
which cells the calibration path sees. Recorded, not waved through.

⚠️ **A reporting trap worth naming:** `scorecard.py`'s `level shift` row prints the mean **without its
sign** (A reads `5.713`, actually **−5.713**). The signed means move −5.713 → +2.267, which looks like
an improvement; the **magnitudes move 13.12 → 18.66, i.e. worse**. Read that row per-fold, never
aggregate.

---

## 5. Conclusion

1. **The primary result is INCONCLUSIVE**, by the rule fixed in advance. Nothing about keeping or
   discarding HFF's labels is licensed.
2. **The experiment as designed cannot answer its question.** The treatment is entangled with a refit
   of the ΔAge target (C-I) and with a 400× up-weighting of out-of-clock-range donors (C-II). That is
   a violation of the one-change rule **in the design**, discovered by running it.
3. **The machinery is sound.** Guards fired, both arms landed on their predicted label counts, arm A
   reproduced baseline exactly, C-5 Option 2 cost the control arm nothing. What failed is the
   comparison's *validity*, not its execution.
4. **The plan's own prediction was right.** C-5's second consequence was written down before the run
   and is exactly what dominates the result.

## 6. What would make it answerable

| | change | cost |
|---|---|---|
| **1** | **Decouple the deconfounder from the training-label policy.** It needs cells where ΔAge is *computable*, not cells whose labels are *used*. `age_mask_reason` (C-6) already distinguishes `dataset_policy` from `cancer_source` / `donor_out_of_clock_range` — fit the deconfounder on everything not excluded for `dataset_policy`, and both arms share one transform by construction. Removes C-I. | small code change + rerun arm B |
| **2** | **Decide C-2 before, not after.** With HFF masked, 40 % of the signal is out-of-clock-range. C-2 stops being a separate later experiment and becomes part of this one's specification. | a decision |
| **3** | **Report per-fold, and pre-register the outlier rule.** N2 carries the SD. Whatever is done about it must be fixed in advance. | a decision |
| **4** | **Accept the power ceiling honestly.** Even perfectly de-confounded, 6 paired folds resolve Δ\* = 3.57 only if SD ≤ 1.91. Observed is 9.60. If fixing C-I does not collapse the SD, the design is under-powered whatever else is done, and that should be said rather than run repeatedly. | — |

---
---

# STEP 6 RERUN (post-C-I) — full report

**Run date** 2026-08-03. Both arms rebuilt from raw GEO and retrained, 6 folds each,
`age_window_k = 4` in both, **arm-suffixed roots** so both survive. C-I fixed: `y_age` no longer
depends on the training-label policy.

## 1. The fix worked — verified three ways, not asserted

| check | first run | rerun |
|---|---|---|
| deconfounder coef, N2 | A `−3.93, −3.42` vs B **`−24.20, +10.12`** | A `−3.9289, −3.4207` vs B **`−3.9289, −3.4207`** |
| …all six folds | wildly different | **identical in all 6** |
| `y_age` across arms (pre-flight, row-exact) | — | **`max|Δ| = 0.000e+00` over 7 062 rows** |
| arm A vs pre-C-I arm A | — | **`max|Δ| = 0.000e+00` on all 6 folds** |

C-I was a no-op on the control (as predicted — `deconfound_mask == age_mask` there) and removed the
target-variable drift in the treatment. **The comparison is now genuinely one-change.**

## 2. Primary — all 6 folds

| fold | arm A | arm B | Δ | *(first run Δ)* |
|---|---:|---:|---:|---:|
| N2 | 21.79 | 29.55 | **+7.76** | *+21.37* |
| N3 | 29.69 | 22.53 | −7.16 | *−7.22* |
| O1 | 5.39 | 7.61 | +2.22 | *+3.63* |
| O2 | 7.54 | 7.36 | −0.17 | *+1.84* |
| Y1 | 7.28 | 8.59 | +1.31 | *+5.35* |
| Y2 | 14.06 | 14.07 | +0.02 | *−1.15* |

| quantity | rerun | first run |
|---|---|---|
| **effect** | **+0.661 yr** | +3.971 |
| **observed SD** | **4.808 yr** | 9.599 |
| **MDE** (1.049 × SD) | **5.045 yr** | 10.074 |
| 95 % CI | **[−4.384, +5.707]** | [−6.102, +14.045] |
| power for Δ\* = 3.57 | 31.5 % | 11.3 % |

**The SD halved** — C-I accounted for about half the variance. The effect shrank six-fold, to
**+0.661 yr, 4.6 % of the 14.29 yr baseline.**

**Verdict: CI includes 0 and MDE 5.045 > Δ\* 3.572 → INCONCLUSIVE. Licenses nothing.**

## 3. Secondary — the pre-registered 4 in-range folds

Registered *before* arm B ran, on the C-II grounds that N2/N3 are donor age 0, outside the clock's
`age_range = [1.0, 96.0]`.

| | |
|---|---|
| per-fold Δ (O1, O2, Y1, Y2) | `+2.22, −0.17, +1.31, +0.02` |
| effect | **+0.843 yr** |
| **observed SD** | **1.130 yr** |
| MDE (1.591 × SD) | 1.799 yr |
| 95 % CI | **[−0.956, +2.642]** |

### 🔬 C-II is confirmed, and it is the dominant remaining variance source

**Dropping N2 and N3 collapses the SD 4.808 → 1.130 — a factor of 4.3, on 2 of 6 folds.** The two
donors outside the clock's validated range carry almost all the fold-to-fold variance. That was a
hypothesis when the first run ended; it is now measured.

### ⚠️ But the "we were powered" reading does NOT survive scrutiny

At face value MDE 1.799 ≤ Δ\* 3.572, which under the registered table would read *"the labels are
genuinely not contributing → mask them."* **That reading is not safe**, because with n = 4 the SD is
itself a noisy estimate. The χ² 95 % interval on the true σ:

| | SD observed | 95 % CI on σ | MDE at σ_high | vs Δ\* |
|---|---:|---:|---:|---|
| primary (n=6) | 4.808 | [3.001, 11.792] | 12.375 | ≫ Δ\* |
| secondary (n=4) | 1.130 | [0.640, **4.213**] | **6.704** | **> Δ\*** |

**Neither analysis is robustly powered once σ is admitted to be an estimate.** So the secondary is
INCONCLUSIVE too — for a different and more interesting reason than expected. It was pre-registered
as "underpowered by construction"; it turned out *better* powered than the primary, and still not
enough.

## 4. Secondary and guard metrics

| metric | arm A | arm B | verdict | *(first run)* |
|---|---:|---:|---|---|
| `rank_model_dage` | 0.948 | 0.879 | **REGRESSION** [−0.100, −0.037] | *−0.186* |
| `rank_ridge_dage` | 0.955 | 0.891 | **REGRESSION** [−0.110, −0.018] | *−0.146* |
| `dage_mae_ridge` | 14.05 | 17.54 | noise [−0.130, +7.099] | *REGRESSION +9.21* |
| `interval_width` | 65.9 | 72.0 | noise [−8.48, +20.68] | *REGRESSION* |
| `fate_prauc` | 0.992 | 0.978 | noise ✅ | noise ✅ |
| `fate_roc` | 0.983 | 0.966 | noise ✅ | noise ✅ |
| `fate_ece` | 0.249 | 0.320 | noise ✅ | noise ✅ |
| `fate_ece` (Platt) | 0.140 | 0.236 | **REGRESSION** [+0.011, +0.182] ⚠️ | *REGRESSION* |
| `ood_flag_rate` | 0.273 | 0.516 | **+0.243** [+0.068, +0.419] | *+0.162* |

**Ranking still degrades detectably** (−0.069 Spearman, CI excludes 0) — and ridge degrades by
nearly the same amount (−0.064). Since `scorecard.py:95` refits ridge on the *same* masked labels,
that is two learners both ranking worse on 75 labels than on 33 688. Consistent, and no longer
confounded by a moving target.

**`fate_ece` (Platt) regressed again**, +0.096 [+0.011, +0.182], surviving C-I. It is now a
**standing anomaly**: the fate head consumes no ΔAge, so a persistent move needs explaining rather
than repeating. **`ood_flag_rate` doubled** — arm B's model finds the held-out donor much more
out-of-distribution.

## 5. Conclusion of the rerun

1. **C-I is fixed and verified.** Identical deconfounders, bit-identical `y_age`, control untouched.
   The first run's dominant confound is gone.
2. **The primary is still INCONCLUSIVE** — but for an honest reason now (limited power), not because
   the arms measured different things.
3. **The effect is small.** +0.661 yr, 4.6 % of baseline, six-fold smaller than the confounded
   estimate. Masking 99.7 % of the age labels costs little on ΔAge MAE.
4. **C-II is real and is now the binding constraint.** Two out-of-range donors carry 4.3× the
   variance of the other four. **No further re-running of this design will fix that.**
5. **Ranking is the one place a consistent cost shows up** (−0.069, CI excludes 0), in both learners.

## 6. What this does and does not license

**Does not license:** discarding HFF's labels. Neither analysis is robustly powered for Δ\*.

**Does support, weakly:** the effect on ΔAge MAE is small — the confounded run's +3.97 was mostly the
confound, not the labels.

**The honest next step is not another 10 h of the same design.** The binding limit is 6 paired folds
with 2 of them out of the clock's validated range. Options, in order of what they cost:

| | |
|---|---|
| **Accept and report** | State that step 6 cannot resolve Δ\* at this geometry, and carry both arms' numbers forward as a bounded estimate rather than a verdict. Costs nothing. |
| **Fix C-II at the source** | Get a clock validated at age 0, or exclude the neonatal donors from the *design* — not just the analysis — which means a 4-donor study and accepting the power that implies. |
| **Change the estimand** | The ranking metric shows a consistent, detectable effect where MAE does not. If ranking is what Stage 2 actually needs, register it as primary — but that is a new pre-registration, not a re-read of this one. |
