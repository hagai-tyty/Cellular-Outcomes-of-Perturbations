# STAGE 1.5.1 — NEW **V2** (independent adjudication + corrected plan)

**Supersedes** `STAGE_1_5_1_NEW_CHANGES.md` as the execution guide. Both earlier documents —
`STAGE_1_5_1_CLOCK_PRECISION.md` (V1) and `STAGE_1_5_1_NEW_CHANGES.md` (the review) — are left
**byte-unmodified**, per the project's annotate-never-rewrite rule.

**Status:** adjudication **EXECUTED** 2026-07-26 against real GSE113957 + Gill. Sixteen tests in
`experiments/stage_1_5_1_tests.py` → `stage_1_5_1_tests_results.json`. **Nothing was accepted on
either document's report — including my own first draft of this file, whose errors are recorded in
§2.** Execution steps (§7) are PLANNED. `git diff --stat src/` is empty.

---

## 0. Bottom line

| | Finding | Consequence |
|---|---|---|
| ✅ | The review's **R4 refutation** and **C3 elimination** are both correct — reproduced independently | V1's leak-audit is cancelled; C3 demoted |
| ❌ | The review's only 🔴 CORRECTNESS claim (**P3**) is wrong in severity, and its suggested remedy would **corrupt the clock** | §1.3 — the ten are **progeria patients** |
| 🔴 | **V1's PASS bar (mean ≤ 4.0) is unreachable** — published SOTA on this exact data is **mean 7.7** | §4 |
| 🔴 | **V1's lead candidate C1 (ElasticNet) is worse than the control** — 12.93 vs 12.27, and *not* a convergence artefact | §6 |
| 🔴 | **Matching published SOTA would still leave SNR ≈ 1.0.** The clock fix, executed perfectly, does not on its own make the effect resolvable | §5 |
| 🔴 | **`MASTER_PLAN` §5b-ter already prescribes the cheaper, stronger lever** (condition-level scoring) and it is **untested** | §5 |
| 🆕 | The published **LDA-ensemble family — now tested by me** — beats ridge (11.28 / 7.50 vs 12.27 / 9.47) but does not reach 7.7 / 4.0 | §6 |
| ✅ | The `√2·cv_mae` label-noise assumption both documents relied on is **empirically sound**, though imprecise (CI [12.7, 30.4]) | §3 |

**The SNR ≈ 1 diagnosis survives.** What changes is that clock precision is **one of two necessary
levers, and the more expensive one** — a material amendment to V1's framing of it as *the*
root-cause fix.

---

## 1. Adjudication of the review's three claims

Every number is from my own run.

### 1.1 P1 / R4 "no scaler leak" → ✅ **UPHELD**

| Test | Method | Result |
|---|---|---|
| T1 | `normalize_counts(X[:5])` vs `…(X)[:5]`, then rescale rows 10+ ×7 and re-check rows 0–4 | **0.000e+00 both ways** → strictly per-row |
| T2 | token scan of `clock_fit.py` | **no** scaler tokens |
| T3 | donor identity, both series matrices | **143 values, 143 unique, 0 repeated** → no group leakage |
| T4b | exact protocol replay, correct sample set | **12.27 / 0.837 / alpha 0.2721** — matches the artefact to 4 s.f. |

**R4 refuted; V1's Step 1 leak-audit is cancelled.** The review's reasoning (normalisation ≠
standardisation) is right.

**Keep the guard it motivated** — the review says this and is correct. One refinement neither
document makes: the published gene filter (dynamic range + abundance) is **age-blind**, so it does
not leak and may sit outside the fold; **supervised** selection may not. The guard must distinguish
the two, or it will forbid a legitimate step.

### 1.2 P2 / C3 "slope recalibration eliminated" → ✅ **UPHELD**
Out-of-fold recalibration: **12.67 → 13.16 (+3.9% worse)**; review measured +1%. Same conclusion.

### 1.3 P3 "ten samples unexplained, 🔴 CORRECTNESS" → ❌ **SEVERITY WRONG, REMEDY DANGEROUS**

Arithmetic right (143 in GEO, 133 in the artefact); everything built on it wrong. Neither document
read the `disease` field:

```
T10:  disease across all 143  ->  {'Normal': 133, 'HGPS': 10}
```

**The ten are Hutchinson–Gilford Progeria patients.** Fleischer et al. trained on **133 healthy** and
used **10 progeria as a separate validation set**. Confirmed decisively:

```
T4b Normal-only n=133 :  cv_mae 12.27 | pearson 0.837 | alpha 0.2721
shipped artefact      :  cv_mae 12.27 | pearson 0.837 | alpha 0.2721      <- exact
```

1. The exclusion is **recorded** (GEO metadata) and **scientifically required** — HGPS fibroblasts
   look old at young chronological ages; training on them teaches the clock backwards.
2. **The review's remedy is harmful.** It floats that the ten may be "**7.5% more training data**".
3. **Both prior replications were contaminated** — the review's 12.67 *and my own first run's 12.67*
   included the 10 HGPS. The review blamed "a newer NCBI annotation"; the cause was solely the
   sample set, and the correct set reproduces exactly.

**P3 → 🟡 documentation.** Step 0 becomes "filter `disease == Normal`", not "investigate".

---

## 2. Adjudication of my own first V2 draft — three errors

A review that only audits the other side is not a review. My draft asserted things I had not tested.

| # | My error | What testing showed |
|---|---|---|
| **M1** | I criticised V1 for leading with an **untested** candidate (C1) — then proposed **C5 (LDA ensemble) as the new lead on the paper's authority alone.** Same error. | Now tested (T12): C5 **does** beat ridge — 11.28 / 7.50 vs 12.27 / 9.47 — but reaches **nowhere near** the published 7.7 / 4.0. It is promising, **not** proven |
| **M2** | I hypothesised the alpha grid was mis-specified (alpha pinned at the 0.1 edge) | **Artefact of the HGPS contamination.** On the correct set alpha = 0.2721, not at the edge; widening the grid gains **0.07 yr** (T11). Withdrawn |
| **M3** | I claimed the paper's **leave-one-out** vs our **5-fold** was "not apples-to-apples" and might explain part of the gap | **It explains −0.03 yr** (T16: 12.27 vs 12.30). The gap is method, not protocol. Withdrawn |
| **M4** | I quoted the measured label noise as "17.90, ratio 1.03" — implying precision it does not have | n=12 pairs ⇒ 95% CI **[12.68, 30.40]** (T14). The assumption sits inside it, but the measurement is coarse |

---

## 3. The label-noise assumption — checked, and its uncertainty stated

Both documents derived the bar from `ΔAge noise ≈ √2·cv_mae` without testing it. That deserved a
check: `cv_mae` is *between-donor*, while ΔAge is a *within-donor difference*, so any per-donor
systematic offset **cancels** (as proved in Stage 1.5 §2 Group A). If most of `cv_mae` were per-donor
bias, the bar would be far too strict.

Measured with no reference to `cv_mae` — Gill's Exp1/Exp2 pairs share donor, day and marker, so their
true ages are identical:

```
12 pairs -> SD of differences = 17.90 yr   [MEASURED]   95% CI [12.68, 30.40]
            sqrt(2)*cv_mae    = 17.35 yr   [ASSUMED]    -> inside the CI
```

**The assumption is sound.** Single-measurement within-donor SD is **12.66 yr ≈ cv_mae 12.27** — the
clock's error is *not* mostly per-donor bias, so it does **not** cancel in ΔAge. V1's arithmetic is
vindicated on evidence. *(Caveat: the CI is wide; treat 17.9 as "≈ the clock's own error", not as a
precise constant.)*

---

## 4. 🔴 The bar is unreachable

V1's **PASS = mean MAE ≤ 4.0**, derived (correctly, §3) from what the tool needs. Neither document
asked what is *achievable*:

> **Fleischer et al. 2018**, *these exact 133 samples*: **median 4.0 yr, mean 7.7 yr**, R² 0.81, via
> an **ensemble of LDA classifiers** over ~4,852 filtered genes — having **explicitly tested and
> rejected ridge regression, linear regression and SVR**. (Confirmed by a second source.)

| | mean | median |
|---|---|---|
| our shipped clock | 12.27 | 9.47 |
| **published, same data** | **7.7** | **4.0** |
| V1's PASS bar | **≤ 4.0** | — |

**V1's bar demands beating published SOTA by ~1.9× on its own data**; field-leading transcriptomic
clocks report 4–6 yr *median*. By ground rule §5b — *a bar a correct system cannot pass is a
description, not a test* — **it fails resolvability and must be corrected before running.**

**Fix — split the bar (do not simply lower it):**

| Bar | Question | Threshold | Achievable? |
|---|---|---|---|
| **E — Engineering** | is our clock as good as the method allows? | mean ≤ 7.7 **and** median ≤ 4.0 | plausibly — a published method reaches it |
| **S — Sufficiency** | is it good enough for the tool? | mean ≤ 4.0 | **almost certainly not on this dataset** |

**Also fixed:** V1 compares our *mean* to a literature *median*. **Report both, always** — the
distributions are skewed (ours 12.27/9.47; published 7.7/4.0).

---

## 5. 🔴 The deeper problem: the clock fix cannot succeed alone

Neither document asked the obvious follow-up — *if we hit the bar, is the problem solved?* Measured
(T15), using the **measured** single-sample SD of 12.66 yr against the recorded **−11.35 yr** effect
(`MASTER_PLAN` §5b-ter). SNR ≥ 2.0 is needed for a 2σ-resolvable effect:

| clock | k=1 | k=3 | k=10 |
|---|---|---|---|
| **current** | 0.61 | 1.06 | **1.94** |
| **published SOTA (E bar)** | **0.98** | 1.70 | 3.10 |
| sufficiency target (S bar) | 1.88 | 3.26 | 5.96 |

Three readings, all consequential:

1. **Achieving the E bar leaves SNR ≈ 1.0.** Matching the published clock — the best realistic
   outcome of this entire stage — **does not make the effect resolvable.**
2. **Even the S bar reaches only 1.88 at k=1**, short of its own 2.0 target. V1's arithmetic was
   marginally optimistic.
3. **The current clock at k=10 (SNR 1.94) ≈ a perfect clock at k=1 (1.88).** **Replication is as
   powerful a lever as a 3× better clock — and far cheaper.**

### The project already knew this, and it is untested

`MASTER_PLAN` §5b-ter contains the same arithmetic and its conclusion:

> *"Uncertainty on a mean shrinks by √n … n=21 cells → SE 3.7–4.6 yr → comfortably detectable"*
> *"**RES should score CONDITIONS (populations of cells), not individual cells.**"*

So the cheaper lever is **already specified in the project's own master plan and has never been
run.** V1 proposed the expensive lever without referencing it.

### But the two levers fix *different* defects — both are needed

This is the reason Stage 1.5.1 still earns its place:

| lever | fixes | cannot fix |
|---|---|---|
| **Aggregation** (§5b-ter) | **variance** — SE shrinks by √n | **bias**: the slope is **0.717**, so ΔAge magnitudes are compressed ~28% low. Averaging a biased estimator gives a *precise wrong answer* |
| **Clock precision** (this stage) | **bias/compression**, and variance | it cannot reach SNR 2 alone (row 2 above) |

**Amendment to V1's framing:** Stage 1.5.1 is **necessary but not sufficient**, and it is the more
expensive of the two levers. It should not be sold as *the* root-cause fix.

⚠️ **Scope limit on aggregation, so it is not oversold:** averaging helps where many cells share a
condition (HFF single-cell; a patient sample at deployment). It does **not** help Gill's *training
labels*, which are ~1 bulk sample per donor-timepoint and cannot be averaged across timepoints
without mixing biology.

---

## 6. Candidates — now tested rather than argued

All on the correct Normal-only set, selection inside each fold.

| candidate | mean | median | pearson | slope | verdict |
|---|---|---|---|---|---|
| **C5 LDA ensemble** (approximated) | **11.28** | **7.50** | 0.811 | 0.709 | **best tested** — beats control, far short of published |
| C5b LDA, single binning | 11.82 | **6.00** | 0.786 | 0.795 | best *median*; ensemble wins on mean |
| C4 dense ridge (control) | 12.27 | 9.47 | 0.837 | 0.717 | control |
| C1 ElasticNet (best of l1 0.1/0.5/0.9) | 12.93 | 10.92 | 0.823 | 0.684 | ❌ **worse than control** |
| C3 slope recalibration | 13.16 | — | — | — | ❌ worse |

**C1's failure is robust, not a convergence artefact** (T13): loose (`max_iter` 3000, `tol` 1e-3) and
tight (60000, 1e-6) give **identical** MAEs with **zero** convergence warnings, across three
`l1_ratio` values. **Sparse linear regression does not beat dense** — the *linear family* is the
limit, exactly what Fleischer found when they rejected ridge, linear regression and SVR.

**On C5, honestly:** my implementation approximates the paper (age-blind filter → 6,625 genes →
in-fold PCA-50 → 20 staggered LDA classifiers). It improves the mean 8% and the median 21% over
ridge, confirming the family is the right direction — but it lands at 11.28, not 7.7. The remaining
gap is the paper's exact recipe (FPKM units, ~4.8k genes, no PCA). **C5 is the best lead available
and it is not yet close to the target.**

---

## 7. Corrected execution plan

### Step 0 — pin the sample set: filter `disease == Normal` (n=133)
Record the rule and the 10 HGPS GSM ids in the clock metadata so Stage 5 can state its training set.
**Reserve the 10 HGPS as a held-out check** — reproducing the paper's accelerated-ageing result on
them is a free, strong external validation. *(Neither earlier document proposed this.)*

### Step 0b (NEW, and it should come first) — test the cheaper lever
Before spending this stage's effort, run the **condition-level aggregation** `MASTER_PLAN` §5b-ter
already prescribes and quantify SNR at realistic n. If it delivers, Stage 1.5.1's *urgency* drops
even though its *necessity* (the compression bias) remains.

### Step 1 — complete
R4 ✅ refuted, slope ✅ 0.717, decile profile ✅ (ratio 1.94; bias +7.6 young → −14.1 old). The
budgeted leak-audit is **cancelled**.

### Step 2 — pursue **C5 exactly**, report C1/C3/C4 as measured
The family is validated (11.28); the gap to 7.7 is recipe fidelity. Reproduce the paper's units,
gene filter and bin width before concluding anything about achievability.

### Step 3 — gate on **both** bars
**E fail ⇒ stop and diagnose.** **E pass / S fail ⇒ ship the better clock AND invoke the fallback** —
success and scope finding simultaneously, which §5 says is the *expected* outcome.

### Step 4 — rebuild + revalidate
Only after Step 3, and only if the new clock beats the shipped one on E. V1 §5's V1–V7 stand.

---

## 8. Honest expectation, recorded before Step 2

| outcome | estimate | basis |
|---|---|---|
| **E pass** (mean ≤7.7, median ≤4.0) | **~40%** | C5 tested at 11.28; closing to 7.7 needs recipe fidelity that may not survive the counts-vs-FPKM difference. Lower than my draft's 55% because I have now measured the gap instead of assuming the paper transfers |
| **S pass** (mean ≤4.0) | **~5%** | requires beating published SOTA by ~1.9×. The review's 25–35% was made without the literature comparison |
| **Neither** | **~55%** | |

**Most likely: a materially better clock that still does not make per-cell quantification
measurable** — which §5 shows was never achievable by this lever alone.

---

## 9. Unchanged from V1

The SNR ≈ 1 diagnosis; the revalidation suite (§5 V1–V7); the fallback (§6); the discipline (§7);
the guard-streak warning (`y_age` moves ⇒ the `+0.000` record ends **by construction**); and
`FRAGILE` reporting within 0.5 yr of any boundary.

---

## 10. Reproduction

```bash
python experiments/stage_1_5_1_tests.py "D:\GSE113957" "D:\Gill"
```

| Test | Question | Result |
|---|---|---|
| T1/T2 | cross-sample scaler? | none — `0.000e+00`, no tokens |
| T3 | group leakage? | 143/143 unique donors |
| T4 | replay CV on all 143 | 12.67 / 0.841 / slope 0.717 / **252× in-sample gap** |
| **T4b** | replay on Normal-only 133 | **12.27 / 0.837 / 0.2721 — exact match** |
| T5 | does C3 help? | no, +3.9% |
| T6 | sample count | 143 parsed vs 133 shipped; **0 below age 1** (R5) |
| **T7** | ΔAge label noise, measured | **17.90 yr** (assumed 17.35) |
| T8 | sparse feasibility | 13.01 — FAIL |
| T9 | error by decile | ratio 1.94; bias **+7.6 → −14.1** |
| **T10** | who are the 10? | **Normal 133 / HGPS 10** |
| T11 | alpha at grid edge? | no (0.2721); widening gains 0.07 |
| **T12** | does LDA-ensemble beat regression? | **11.28 / 7.50** vs control 12.27 / 9.47 |
| **T13** | is C1's failure a convergence artefact? | **no** — loose ≡ tight, 0 warnings |
| **T14** | how precise is the label noise? | 95% CI **[12.68, 30.40]** |
| **T15** | does hitting the bar fix SNR? | **no** — E bar ⇒ SNR **0.98** at k=1 |
| **T16** | LOO vs 5-fold? | **−0.03 yr** — protocol is not the gap |

**Sources:** [Fleischer et al. 2018, *Genome Biology*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6300908/) ·
[BiT age (Meyer & Schumacher 2021)](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13320) ·
[ATAC-clock (2023)](https://link.springer.com/article/10.1007/s11357-023-00986-0)
