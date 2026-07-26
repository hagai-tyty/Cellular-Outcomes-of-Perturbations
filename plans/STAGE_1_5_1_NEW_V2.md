# STAGE 1.5.1 — NEW **V2** (independent adjudication + corrected plan)

**Supersedes** `STAGE_1_5_1_NEW_CHANGES.md` as the execution guide. Both earlier documents —
`STAGE_1_5_1_CLOCK_PRECISION.md` (V1) and `STAGE_1_5_1_NEW_CHANGES.md` (the review) — are left
**byte-unmodified**; everything here is additive, per the project's annotate-never-rewrite rule.

**Status:** adjudication **EXECUTED** 2026-07-26 on the real GSE113957 + Gill. Every claim below was
re-derived from raw data by `experiments/stage_1_5_1_tests.py`
(→ `stage_1_5_1_tests_results.json`); **nothing was accepted on either document's report.** The
execution steps in §6 are PLANNED, not run. `git diff --stat src/` is empty.

---

## 0. Bottom line — what actually changed

| | Finding | Consequence |
|---|---|---|
| ✅ | The review's **R4 refutation is correct** — I reproduced all three of its tests independently | V1's Step 1 leak-audit is genuinely unnecessary |
| ✅ | The review's **C3 elimination is correct** | C3 drops to a reported control |
| ❌ | The review's **P3 ("🔴 CORRECTNESS", ten samples unexplained)** is **wrong in severity and its suggested remedy is dangerous** | §2 — the ten are **progeria patients** |
| 🔴 | **V1's primary bar (`cv_mae ≤ 4.0`) is unreachable**: the published result on this exact dataset is **mean MAE 7.7 yr** | §4 — **the single most important correction in this document** |
| 🔴 | **V1's lead candidate C1 (ElasticNet) is worse than the control** — measured 13.01 vs 12.27 | §5 — the plan's main hope is already falsified |
| 🆕 | The published method (**LDA ensemble on ~4.8k filtered genes**) is absent from both documents, and the paper **explicitly tested and rejected ridge** | §5 — this becomes the lead candidate |
| ✅ | The bar's `√2·cv_mae` label-noise assumption — never checked by either document — is **empirically correct** (measured 17.90 vs assumed 17.35) | §3 — V1's arithmetic is vindicated |

**The precision diagnosis (SNR ≈ 1) survives everything.** What changes is *which bar to hold it to*,
*which method to try*, and *what outcome to expect*.

---

## 1. Adjudication of the review's three claims

Every number below is mine, from my own run — not the review's.

### P1 / R4 — "no scaler leak exists" → ✅ **UPHELD**

| Test | Method | My result |
|---|---|---|
| T1 | `normalize_counts(X[:5])` vs `normalize_counts(X)[:5]`; then rescale rows 10+ ×7 and re-check rows 0–4 | `0.000e+00` **both** ways → strictly per-row |
| T2 | token scan of `clock_fit.py` for `StandardScaler`/`scaler`/`fit_transform`/… | **NONE** |
| T3 | donor identity across both series matrices | **143 values, 143 unique, 0 repeated** → no group leakage either |
| T4b | exact protocol replay on the correct sample set | **cv_mae 12.27, pearson 0.837, alpha 0.2721** — matches the artefact to 4 s.f. |

**R4 is refuted, and V1's Step 1 leak-audit would indeed find nothing.** The review's reasoning
(normalisation ≠ standardisation) is correct.

**But keep the guard it motivated.** The review says this and is right: gene selection by
age-correlation *is* a cross-sample statistic and **will** leak outside the fold. V1's §3 guard stays,
and §5's new candidates make it load-bearing. *(Nuance neither document notes: the published feature
filter — expression range and abundance — is **age-blind**, so it does not leak and may sit outside
the fold. Supervised selection may not. The guard must distinguish the two.)*

### P2 / C3 — "slope recalibration eliminated" → ✅ **UPHELD**

Out-of-fold recalibration, my run: **12.67 → 13.16 (+3.9% worse)**. The review measured
12.67 → 12.78 (+1%). Magnitudes differ; the conclusion is identical and directionally robust.
Correcting the slope of a model that is memorising rescales errors without removing them.

### P3 — "ten samples unexplained; 🔴 CORRECTNESS" → ❌ **OBSERVATION RIGHT, SEVERITY WRONG, REMEDY DANGEROUS**

The **arithmetic is right** (143 in GEO, 133 in the artefact) and worth documenting. Everything built
on it is not. Neither document read the `disease` field:

```
T10:  disease across all 143 samples  ->  {'Normal': 133, 'HGPS': 10}
```

**The ten are Hutchinson–Gilford Progeria (HGPS) patients.** Corroborated externally: Fleischer et
al. 2018 trained on **133 healthy** individuals and used **10 progeria patients as a separate
validation set**. And confirmed decisively by measurement:

```
T4b (Normal-only, n=133):  cv_mae 12.27  |  pearson 0.837  |  alpha 0.2721
shipped artefact        :  cv_mae 12.27  |  pearson 0.837  |  alpha 0.2721
```

**An exact reproduction.** So:

1. The exclusion is **not** unrecorded — it is in the GEO metadata, and it is **scientifically
   required**: HGPS fibroblasts show accelerated ageing at young chronological ages, so training on
   them teaches the clock that old-looking transcriptomes are young.
2. **P3's suggested remedy is actively harmful.** The review floats that the ten may have been
   "dropped by an accident of parsing" and are "**7.5% more training data** for a p≫n problem."
   Adding them would corrupt the clock.
3. **Both replications reported so far were contaminated.** The review's 12.67 *and* my own first
   run's 12.67 both included the 10 HGPS. The review attributed the 12.67-vs-12.27 gap to "143 vs
   133 and a newer NCBI annotation" — the sample count was the cause, but the annotation was not,
   and the correct set reproduces exactly.

**P3 downgrades from 🔴 CORRECTNESS to 🟡 documentation.** Step 0 is still worth doing — but it is
"filter `disease == Normal`", not "investigate an unknown exclusion".

---

## 2. The label-noise assumption — checked for the first time

Both documents derive the bar from `ΔAge noise ≈ √2 · cv_mae`. That is an **assumption**, and a
questionable one: `cv_mae` is a *between-donor* error, while ΔAge is a *within-donor* difference, and
any per-donor systematic offset **cancels** in a difference (the same cancellation proved in Stage
1.5 §2 Group A). If most of `cv_mae` were per-donor bias, the true label noise would be far smaller
and **the bar would be far too strict**.

Measured directly (T7), with no reference to `cv_mae`: Gill's Exp1/Exp2 pairs are the **same donor,
same day, same marker**, so their true ages are identical and the spread of their predicted-age
differences *is* the ΔAge label noise.

```
12 matched pairs ->  SD of differences        = 17.90 yr   [MEASURED]
                     sqrt(2) * cv_mae         = 17.35 yr   [ASSUMED by both documents]
                     ratio                    = 1.03
```

**The assumption is correct to 3%.** Single-measurement within-donor SD is **12.66 yr** ≈ `cv_mae`
12.27 — i.e. the clock's error is *not* mostly per-donor bias; it is per-measurement noise that does
**not** cancel. V1's bar arithmetic is vindicated, now on evidence rather than assumption.

*(Scope: this is the noise on **Gill bulk** labels — the age-valid arm. D2's pseudobulk replicate SD
was ~1 yr because averaging 400 cells suppresses sampling noise; that is technical reproducibility,
not accuracy, and does not transfer to Gill's one-sample-per-donor-day design.)*

---

## 3. Root causes — what the evidence now says

| | V1's claim | Verdict on my evidence |
|---|---|---|
| **R1** dense ridge at p/n≈250 overfits | ✅ **confirmed, dramatically** — in-sample MAE **0.05** vs CV **12.67** on n=143 (**252×**). Near-total memorisation |
| **R2** shrinkage compression | ✅ **confirmed and quantified** — slope **0.717**, and T9's decile bias is textbook regression-to-the-mean: young read **old** (+7.6, +11.7), old read **young** (−14.1, −8.1) |
| **R3** fragile to gene-set mismatch | ⏳ still plausible, untested here |
| **R4** CV optimism | ❌ **refuted** (§1) |
| **R5** age-0 out of range | ✅ **confirmed structural** — T6: **0 samples below age 1**. No refit on this dataset can fix it |

### ⚠️ A hypothesis of my own that I must withdraw

My first run showed `alpha = 0.1` sitting exactly on the grid's lower edge, and I was ready to
report "the alpha grid is mis-specified" as a new finding. It was an artefact of including the 10
HGPS samples. On the correct set, `alpha = 0.2721`, **not** at the edge, and widening the grid
changes almost nothing (T11: 12.27 → 12.25 → 12.20). **The alpha grid is fine; my hypothesis was
wrong.** Recorded because a review document that only reports the other side's errors is not a
review.

---

## 4. 🔴 THE BAR IS UNREACHABLE — the most important correction

V1's primary bar: **PASS = `cv_mae ≤ 4.0 yr`**, derived (correctly, per §2) from what the tool needs.
Neither document asked what is *achievable*. The published result on **this exact dataset**:

> **Fleischer et al. 2018** (Genome Biology), *same 133 samples*: **median absolute error 4.0 yr**,
> **mean absolute error 7.7 yr**, R² 0.81, using an **ensemble of LDA classifiers** over ~4,852
> filtered genes — having **explicitly tested and rejected ridge regression**, linear regression and
> SVR.

| | mean MAE | median |
|---|---|---|
| our shipped clock | **12.27** | 9.47 |
| **published, same data** | **7.7** | **4.0** |
| V1's PASS bar | **≤ 4.0** | — |

**V1's PASS bar demands beating the published state of the art by ~1.9× on its own data.** A perfect
reimplementation of the best known method would score 7.7 and be labelled **FAIL**. Wider context:
field-leading transcriptomic clocks report ~4–6 yr *median* (BiT age 5.55 yr on human cortex;
ATAC-clock 5.27 yr median) — so a **mean** ≤ 4.0 is at or beyond the frontier.

By this project's own ground rule §5b — *a bar a correct system cannot pass is a description, not a
test* — **the V1 bar fails resolvability and must be corrected before running.**

### The fix: split the bar in two (do NOT simply lower it)

V1 conflated two different questions. Separating them keeps the science honest without
goalpost-moving:

| Bar | Question | Threshold | Achievable? |
|---|---|---|---|
| **E — Engineering** | is our clock as good as the method allows? | **mean ≤ 7.7 AND median ≤ 4.0** (match published) | **Yes** — a known method reaches it |
| **S — Sufficiency** | is the clock good enough for per-cell rejuvenation claims? | **mean ≤ 4.0** (from §2's verified arithmetic) | **Probably not on this dataset** |

Both get reported. **The likely outcome is E pass / S fail** — a real 37% error reduction that still
does not make per-cell quantification measurable. That routes to V1 §6's fallback, which is
therefore the **expected** path rather than a remote contingency.

**Also fixed:** V1 compares a mean to a literature median. **Report both, always**; the distributions
are skewed (ours: mean 12.27 vs median 9.47; published: 7.7 vs 4.0).

---

## 5. 🔴 The lead candidate is already falsified — and the right one is missing

V1's Step 2 leads with **C1 (ElasticNet/Lasso)**, "sparse selection, targets R1 + R3". Tested, with
all selection inside each fold, on the correct Normal-only set:

| candidate | mean MAE | median | pearson | slope | non-zero genes |
|---|---|---|---|---|---|
| **C4 dense ridge (control)** | **12.27** | 9.47 | 0.837 | 0.717 | 33,155 |
| **C1 ElasticNet (l1_ratio 0.5)** | **13.01** ❌ | 10.92 | 0.823 | 0.684 | 1,166 |
| C3 slope recalibration | 13.16 ❌ | — | — | — | — |

**Sparse linear regression is *worse* than the dense control.** So the fix is *not* "sparse instead
of dense" — the linear-regression family itself is the limit, which is exactly what Fleischer found
when they rejected ridge, linear regression and SVR in favour of LDA.

### Revised candidate list

| | Candidate | Rationale |
|---|---|---|
| **C5 (NEW, lead)** | **Ensemble of LDA classifiers over staggered age bins** + the paper's age-blind gene filter (5-fold expression range, >5 FPKM ⇒ ~4.8k genes) | the published method; the only one demonstrated to reach 7.7/4.0 on this data. **Absent from both documents** |
| **C2 (keep)** | ridge on a fold-internal filtered gene set | untested, cheap, isolates R3 from R1 |
| C1 (demote) | ElasticNet | **tested, worse than control** — report as a result, not a hope |
| C3 (demote) | slope recalibration | tested, worse — reported control |
| C4 (keep) | dense ridge | control |

**Why C5 should work where C1 failed:** age is being predicted from a compressed, noisy signal;
discretising it into overlapping bins and ensembling classifiers is more robust to that than fitting
one continuous linear map, and it directly attacks R2 (the classifier ensemble has no single
regression slope to compress). This is mechanism, not hope — and it is the mechanism the source
paper reports.

---

## 6. Corrected execution plan

### Step 0 (rewritten) — pin the sample set: **filter `disease == Normal`**
Not "discover which ten". Record the rule (`disease == Normal` ⇒ n=133) and the 10 HGPS GSM ids in
the clock metadata so Stage 5 can state the training set. **Reserve the 10 HGPS as a held-out
validation set** — reproducing the paper's accelerated-ageing result on them is a strong, free
external check of any new clock. *(Neither document proposed this.)*

### Step 1 (mostly complete) — nothing left but confirmation
R4 ✅ refuted, slope ✅ 0.717, decile profile ✅ done (worst/best ratio 1.94, monotonic bias). The
leak audit V1 budgeted is **cancelled**.

### Step 2 (rewritten) — evaluate **C5** and C2 against C4, reporting C1/C3 as measured failures
One leak-free harness; **age-blind** filters may sit outside the fold, **supervised** selection may
not. Report mean **and** median for every candidate.

### Step 3 (rewritten) — gate on **two** bars
Report E and S separately. **E fail ⇒ stop and diagnose** (we cannot match a published method on its
own data). **E pass / S fail ⇒ ship the better clock AND invoke V1 §6's fallback** — that is a
success and a scope finding simultaneously, not a defeat.

### Step 4 — rebuild + revalidate
Only after Step 3, and **only if the new clock beats the shipped one on E**. Everything in V1 §5
(V1–V7) stands unchanged.

---

## 7. Honest expectation, recorded before Step 2

| outcome | my estimate | note |
|---|---|---|
| **E pass** (mean ≤7.7, median ≤4.0) | **~55%** | reimplementing a published method on its own data — the main risks are the FPKM-vs-counts filter mismatch and bin-width tuning |
| **S pass** (mean ≤4.0) | **~10%** | would require beating the published mean by ~1.9×. The review's 25–35% PASS estimate was made without the literature comparison and is, I think, too optimistic |
| **Neither** | **~35%** | C5 fails to transfer; the fallback triggers on weaker evidence |

**Most likely single outcome: a materially better clock that is still not sufficient for per-cell
rejuvenation claims.** Recorded now so that result is read as the pre-registered outcome it is.

---

## 8. What stands unchanged from V1

Not re-litigated: the **SNR ≈ 1 diagnosis** (§0), the **revalidation suite** (§5 V1–V7), the
**fallback** (§6), the **discipline** (§7), the guard-streak warning (`y_age` moves ⇒ the `+0.000`
record ends by construction), and `FRAGILE` reporting within 0.5 yr of any boundary.

---

## 9. Reproduction

```bash
python experiments/stage_1_5_1_tests.py "D:\GSE113957" "D:\Gill"   # -> stage_1_5_1_tests_results.json
```

| Test | Question | Result |
|---|---|---|
| T1/T2 | is there a cross-sample scaler? | no — `0.000e+00`, no tokens |
| T3 | group leakage? | 143/143 unique donors |
| T4 | replay CV on all 143 | 12.67 / 0.841 / slope 0.717 / 252× in-sample gap |
| **T4b** | replay on Normal-only 133 | **12.27 / 0.837 / alpha 0.2721 — exact match** |
| T5 | does C3 help? | no, +3.9% |
| T6 | sample count | 143 parsed vs 133 in artefact; 0 below age 1 |
| **T7** | ΔAge label noise, measured | **17.90 yr** vs 17.35 assumed (ratio 1.03) |
| **T8** | is ≤4.0 reachable by sparse? | **13.01 — FAIL, worse than control** |
| T9 | error by decile | ratio 1.94; bias +7.6 → −14.1 |
| **T10** | who are the 10? | **Normal 133 / HGPS 10** |
| T11 | alpha at grid edge? | no (0.2721); widening gains 0.07 yr |

**Sources for the external comparison:**
[Fleischer et al. 2018, *Genome Biology*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6300908/) ·
[BiT age (Meyer & Schumacher 2021)](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13320) ·
[ATAC-clock (2023)](https://link.springer.com/article/10.1007/s11357-023-00986-0)
