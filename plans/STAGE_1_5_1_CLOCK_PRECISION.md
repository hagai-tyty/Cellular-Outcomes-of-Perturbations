# STAGE 1.5.1 — Clock precision: make the instrument sharper than the effect

**Implements:** the root-cause fix identified at the end of Stage 1.5 §10 (option **B**).
**Depends on:** Stage 1.5 closed (§9 clock validity, §10 D2 replication).
**Blocking for:** Stage 2, and every quantitative rejuvenation claim in Stage 5.
**Scope:** refit + revalidate the aging clock. Touches `src/cellfate/data/clock_fit.py` and the
clock artefact only. **No model, training, calibration or inference code changes.**

**Status:** PLANNED — nothing executed. This document is the pre-registration.

---

## 0. Why this stage exists — the one number

Stage 1.5 ran five independent measurements. Their results only cohere one way:

| measurement | result | §  |
|---|---|---|
| M1 — absolute age on Gill | FAIL: 11.8 yr contrast for a 53 yr true gap | §7 |
| E1 — within-donor trajectory (full) | NO_TREND: mean rho −0.064 | §8.5 |
| E1b — trajectory, reprogramming phase | WRONG_DIRECTION: +0.205, cleared bar by **0.009** | §8.6 |
| §9 — clock validity | clock **reproduces** on its own domain (MAE 0.8 yr, rho 0.99); **tracks** in-range adult age (+18 yr for a 21 yr gap) | §9 |
| D2 — independent replication | CONTRADICTS: −0.214, cleared bar by **0.014** | §10 |

The clock is **not broken** (§9 proves it), yet nothing about reprogramming replicates and the two
trajectory verdicts point in **opposite directions**, each decided by hundredths. Put the magnitudes
side by side and the reason is unambiguous:

| quantity | size |
|---|---|
| clock's own cross-validated error (`cv_mae`, in the artefact's metadata) | **±12.27 yr** |
| per-donor offset Stage 2 exists to correct | **±12.7 yr** |
| D2 spread across the **entire** D0→D14 trajectory (77.1 → 90.2) | **13.1 yr** |
| rejuvenation effect the tool must grade | **~11 yr** |

> **Every effect this project measures is the same size as the measuring instrument's error.**
> Signal-to-noise ≈ 1. At SNR 1 nothing replicates, signs flip between datasets, and verdicts hinge
> on the third decimal — which is precisely the pattern §7–§10 produced.

This is not a modelling problem, a calibration problem, or a target-definition problem. **It is a
precision problem in the instrument**, and it is upstream of all of them: ΔAge is the label the model
trains on, so the label's noise is a hard floor on everything downstream — `sigma_age`, conformal
width, RES, and every quantitative claim in the manuscript.

### Why not the alternatives

- **Option A (proceed to Stage 2)** would spend wet-lab money (k≈3 reference cells/donor) correcting
  a ±12.7 yr offset that is **one clock-error in size** — i.e. possibly measuring nothing.
- **Option C (ship fate-only)** discards a target that §9 showed is *sound in-domain*, and invites
  the reviewer question we cannot currently answer: *"why is your aging clock 3× less accurate than
  the published result on the same dataset?"*

---

## 1. Root-cause hypotheses — why the clock is imprecise

Read from the shipped artefact and `src/cellfate/data/clock_fit.py`:

```
model      : RidgeCV(alphas=logspace(-1, 4, 24))     # DENSE — no feature selection
n_samples  : 133
n_genes    : 33,155        -> every one carries a non-zero weight
cv_mae     : 12.27 yr
cv_pearson : 0.837
intercept  : 72.43
age_range  : [1, 96]
```

**R1 — dense ridge at p/n ≈ 250 is the wrong estimator.** 33,155 features from 133 samples, with
*no* selection: the signal is spread across every gene including pure noise. Ridge controls this only
by shrinking everything toward zero, which drags predictions toward the intercept.

**R2 — the shrinkage signature is visible in our own data.** `cv_pearson = 0.837` (the clock *does*
track age) against `cv_mae = 12.27` (but the magnitudes are wrong) is the classic
compression-to-the-mean pattern. It shows up directly: **every** Gill donor and the HFF line read
*high* and cluster near the 72.4 intercept — N2(0)→98.7, Y2(35)→57.7, O1(53)→79.1, HFF(neonatal)→84.5.

**R3 — a dense clock is fragile to gene-set mismatch.** §9 measured only **57%** of clock genes
present in Gill (89% of weight mass). With weight spread thin over 33k genes, every absent gene
silently contributes 0 and pushes the prediction toward the intercept. A **sparse** clock whose
selected genes are present is structurally more robust to this.

**R4 — possible CV optimism.** `cross_val_predict` is run on `Xn`, which is standardised **before**
the split. If the scaler sees all samples, `cv_mae = 12.27` is itself an *optimistic* estimate and
the true error is worse. To be audited, not assumed.

**R5 — out-of-range is a DATA limit, not a fit limit.** GSE113957 spans ages **1–96**. N2/N3 (age 0)
and HFF (neonatal) sit below it. **No refit on this dataset can fix that** — it must be solved with
data or declared as a scope condition. Stated here so it is not later mistaken for a fixable defect.

---

## 2. What we will do

Ordered so that **cheap measurement gates expensive computation**. Nothing downstream is rebuilt
until the new clock has demonstrably cleared its bar.

### Step 1 — Audit the current fit *(hours, no rebuild)*
Establish an honest baseline before changing anything.
- Re-measure `cv_mae` under **leak-free nested CV** (scaler fitted inside each fold) → tests **R4**.
  If the honest number is worse than 12.27, the SNR problem is *larger* than stated above.
- Report the error profile **by age decile** — quantify the compression of **R2**.
- Report predicted-vs-true slope. A slope well below 1 confirms shrinkage compression.

### Step 2 — Refit candidates *(hours, no rebuild)*
Same data, same normalisation, same production path. Evaluated under one identical leak-free CV
harness so the comparison is apples-to-apples.
- **C1 ElasticNet / Lasso** — sparse selection, directly targeting **R1** and **R3**.
- **C2 Ridge on a pre-filtered gene set** (e.g. top-*k* by age correlation, selected **inside** each
  CV fold to avoid leakage).
- **C3 Ridge + slope recalibration** — the cheapest possible fix for **R2** alone: rescale the
  predictions so the predicted-vs-true slope is 1. Tests how much of the error is pure compression.
- **C4 The current dense ridge** — carried through unchanged as the control.

**All candidates are specified here, before any is run.** Whichever wins, it wins on the
pre-registered bar in §3 — not on a bar chosen afterwards.

### Step 3 — Gate *(decision point)*
Compare candidates on §3's bar. **If none clears the minimum, STOP** and go to §6's fallback. Do not
proceed to the rebuild on a clock that has not earned it.

### Step 4 — Rebuild and revalidate *(expensive: ~4 h GPU + rescore)*
Only if Step 3 passes. Rebuild the dataset with the new clock, retrain the six LOOCV folds, re-run
the full §5 revalidation suite.

---

## 3. THE BARS — set now, before any fit

### Primary bar — clock accuracy

Derived from the science, not chosen for reachability. The tool must resolve a **~11 yr** rejuvenation
effect. ΔAge is a *difference* of two clock readings, so its noise is ≈ `√2 × cv_mae`. For an 11 yr
effect to sit at ≥2σ of the label noise we need `√2 · cv_mae ≤ 5.5`, i.e. **`cv_mae ≤ 3.9`**.

| verdict | `cv_mae` (leak-free CV) | meaning |
|---|---|---|
| **PASS** | **≤ 4.0 yr** | ΔAge label noise ≈5.7 yr; an 11 yr effect is ~2σ. The tool becomes measurable |
| **MARGINAL** | 4.0 – 6.0 yr | label noise ≈8.5 yr; effect is ~1.3σ. Usable for *ranking* only; per-cell quantification still not supported |
| **FAIL** | > 6.0 yr | SNR remains ≲2. The transcriptomic-clock route is not viable on n=133 → §6 fallback |

**Guard against gaming:** `cv_mae` must be measured with the **scaler fitted inside each fold** and
with **any feature selection inside each fold**. A number produced any other way does not count.

### Secondary bars — all must hold

| # | Bar | Why |
|---|---|---|
| S1 | predicted-vs-true **slope ∈ [0.85, 1.15]** on held-out folds | kills the R2 compression rather than hiding it in the MAE |
| S2 | `cv_pearson` **≥ 0.90** | rank tracking must improve too, not just be recentred |
| S3 | **≥ 90% of the new clock's |weight| present in Gill AND in GSE242423** | a sparse clock is worthless if its genes are absent from our data (R3) |
| S4 | error must **not be concentrated in one age decile** | a clock that is precise only in mid-range does not solve our problem |
| S5 | `git diff` touches **only** `clock_fit.py` + the clock artefact | no silent changes elsewhere while a target moves |

### Bar-resolvability check (ground rule §5b)

Before running: confirm via `audit_metrics.bar_verdict` that a clock which genuinely achieves 4.0 yr
would be *reported* as passing at n=133. Per §10's lesson, also report **`FRAGILE`** if `cv_mae`
lands within **0.5 yr** of a boundary — E1b and D2 were both decided by hundredths, and this plan
will not repeat that.

---

## 4. Why this will improve things — the mechanism, not a hope

| finding today | after a 4 yr clock | why |
|---|---|---|
| M1 fails: 11.8 yr contrast for a 53 yr gap | contrast should approach the true gap | the bar was `≥20.2 yr`, derived from `cv_mae`; a 3× smaller error shrinks the bar *and* sharpens the estimate |
| E1/E1b/D2 disagree, each by hundredths | a real trajectory signal becomes resolvable, or its absence becomes conclusive | at SNR 3 a null means *no effect*, not *no power* |
| per-donor offset ±12.7 yr, CI includes zero | offset separates from noise — or is shown to be noise | Stage 2's premise finally becomes testable |
| `sigma_age` ~19 yr ⇒ RES ≡ 0 (per-cell approval arithmetically impossible) | label-noise floor drops ~3× | RES can become non-zero without any change to the RES code |
| conformal width 65.9 yr | narrows with the label noise | the width is honest today *because* the label is noisy |

**Note the sign of the risk.** If the refit succeeds, several Stage 1 conclusions become *stale in a
good way* — they were correct given a noisy label. If it fails, we will have established that
transcriptomic age on n=133 fibroblasts cannot support per-cell rejuvenation claims, which is itself
a publishable scope finding and directly informs Stage 5.

---

## 5. What must be re-cleared after the change

Changing the clock changes `y_age`, so **everything downstream must be re-earned.** The harness
already exists — Stage 1.5 built it — which is what makes this affordable.

| # | Test | Bar |
|---|---|---|
| V1 | `diag_clock_validity.py` — own-domain reproduction | `REPRODUCES`, and **on held-out samples**, not in-sample |
| V2 | `diag_clock_validity.py` — H2 in-range tracking on Gill | `TRACKS_IN_RANGE`, contrast ≥ the (now smaller) bar |
| V3 | `diag_zero_point.py` — M1 absolute age | should now **PASS**; if it still fails, R5 (out-of-range) is the cause and must be shown to be |
| V4 | `diag_e1_trajectory.py` — E1 / E1b | a **conclusive** verdict either way; "no trend" at SNR 3 is a real result |
| V5 | `diag_d2_replication.py` — independent replication | E1/E1b and D2 must now **agree**. Disagreement at high SNR would mean a real biological difference between datasets, not noise |
| V6 | full `pytest` + `ruff check src/ tests/ scripts/` | green |
| V7 | `scorecard.py snapshot` + `compare baseline <tag>` | recorded, with the guard caveat below |

### ⚠️ The guard record restarts — stated in advance

`y_age` changes, so `dage_mae_model`, `rank_model_dage`, `fate_prauc`, `fate_roc` **will move**. The
four-run `+0.000` streak ends **by construction, not by defect**. Interpreting those movements as
regressions would be wrong. Stage 1's PARTIAL verdict is measured against the *old* label and does
not automatically carry over; whether Stage 1 must be re-scored is a decision for §7, taken after V1–V5.

---

## 6. If it fails — the pre-registered fallback

**`cv_mae > 6.0` on every candidate** means transcriptomic age on 133 fibroblast samples cannot reach
the precision this tool needs. That is a finding, not a dead end:

1. **More/better training data for the clock** — the honest fix for a p≫n problem (links to Stage 6).
   Note this also addresses **R5**: a training set including neonatal samples is the *only* way to
   make N2/N3/HFF in-range.
2. **Restrict claims to what the precision supports** — large effects only (adult age range,
   fibroblast→iPSC), explicitly *not* per-cell reprogramming increments. This is option C arrived at
   *with evidence*, which is a far stronger position than adopting it now on suspicion.
3. **Ranking-only framing** — if S2 (`cv_pearson`) improves while `cv_mae` does not, the clock orders
   correctly without calibrated magnitudes; the tool ships as a ranker with quantification disclaimed.

**Not on the menu:** tuning the bar to whatever the refit achieves. The §3 bars are fixed now
precisely so that outcome is unavailable later.

---

## 7. Discipline

- **Pre-register before running.** Predictions for Step 1 and Step 2 go in the lab notebook *before*
  either executes, per the §7–§10 pattern.
- **One change at a time.** The clock is the only thing that moves. No calibration, model or RES
  changes ride along — that rule is what let Stage 1 attribute its results at all.
- **Annotate, never rewrite.** Stage 1.5's conclusions stay as written; anything this stage
  supersedes gets a dated note beside it.
- **Report `FRAGILE`.** Any verdict landing within 0.5 yr of a boundary is labelled fragile in the
  same breath as the verdict (§10's lesson).
- **Record failures.** If the refit does not beat dense ridge, that goes in the notebook with the
  same prominence a success would get.

---

## 8. Deliverables

| Artefact | Purpose |
|---|---|
| `experiments/diag_clock_refit.py` (new, read-only) | the leak-free CV harness: audits the current fit and evaluates C1–C4 on one identical protocol. Writes `diag_clock_refit_results.json` |
| `tests/test_diag_clock_refit.py` (new) | every branch of the verdict logic, per the `verify_1a` lesson |
| `configs/clocks/<new>_clock.json` | the new clock — **only if** §3's bar is cleared |
| notebook + `CHANGES.md` entries | pre-registration, then results, then the §7 decision |

**Nothing in Steps 1–3 changes `src/`.** `git diff --stat src/` must be empty until Step 4, and even
then only `clock_fit.py` may move.

---

# 9. REVIEW (2026-07-26) — verified against the tree. Diagnosis endorsed; **R4 is factually wrong**.

Checked rather than accepted. §0–§8 are left as written; this section annotates.

## 9.1 Verified as stated

| Claim | Check | Result |
|---|---|---|
| dense ridge, 33,155 genes, no selection | read the artefact | ✅ all 33,155 weights non-zero |
| `cv_mae` 12.27, `cv_pearson` 0.837, range [1,96] | read the artefact | ✅ exact |
| primary bar arithmetic (`√2·cv_mae ≤ 5.5` ⇒ `cv_mae ≤ 3.9`) | recomputed | ✅ correct |
| **R1** dense ridge at p/n≈250 overfits | **new evidence, see 9.3** | ✅ **strongly corroborated** |

## 9.2 ❌ R4 is FALSE — there is no scaler leak to audit

R4 claims *"`cross_val_predict` is run on `Xn`, which is standardised **before** the split"*, and
Step 1 budgets time to re-measure `cv_mae` "leak-free". **`clock_fit.py` contains no cross-sample
scaler at all** (`grep StandardScaler|scaler|fit_transform` → nothing). The only transform is
`normalize_counts`, which is **per-row** (`lib = counts.sum(axis=1, keepdims=True)`) — library-size
CP10k+log1p uses one sample's own total and no cross-sample statistic, so it cannot leak. And
`cross_val_predict(RidgeCV(...), ...)` refits `RidgeCV` inside each fold, so alpha selection is
already nested.

**Consequence:** `cv_mae = 12.27` is already an honest leak-free estimate. Step 1's "re-measure
leak-free" will reproduce 12.27 and find nothing. **This does not weaken the plan** — it removes the
one hope that the SNR problem was overstated. The instrument really is ±12.27 yr. Step 1 should be
kept for the *other* two items (error-by-decile and predicted-vs-true slope, which test R2), and
the guard in §3 ("scaler fitted inside each fold") remains correct for the **new** candidates C1/C2,
where feature selection genuinely must be inside the fold or it *will* leak.

## 9.3 New evidence that sharpens R1/R2 — the in-sample/CV gap is 16×

§9's reproduction check (run on the real GSE113957) measured the clock **in-sample**:

| | |
|---|---|
| in-sample MAE (143 samples the clock was fit on) | **0.77 yr** |
| cross-validated MAE (same data, held-out folds) | **12.27 yr** |
| ratio | **≈16×** |

A model that is 16× better on its own training data than on held-out data is not "slightly
over-regularised" — it is **memorising**. This is the most direct evidence for R1 in the record, it
was not available when §1 was written, and it means the ridge penalty is effectively not binding at
33k features / 133 samples. It also explains R2's compression: with alpha too small, the fit chases
training noise, and held-out predictions collapse toward the intercept.

**Implication for Step 2:** C1 (sparse) and C2 (pre-filtered) are attacking the right defect, and
C3 (slope recalibration alone) is unlikely to be sufficient — recalibrating the slope of a
memorising model rescales its errors without removing them. C3 should be kept as the *cheap
control* it is, not as a likely winner.

## 9.4 Why the OOD-detector idea (proposed before §10) is **withdrawn**

Recorded because it was proposed in this thread and should not be quietly dropped. After §9 I
suggested: *treat reprogramming intermediates as out-of-domain and gate ΔAge on the model's existing
OOD detector.* **Three independent reasons it fails:**

1. **It measures the wrong distribution.** The detector is a Gaussian over the *model's* latent `z`,
   fitted on `train_ds` (`train_model.py:291`). The model's training set **contains** the Gill
   reprogramming intermediates — so those cells are *in*-distribution for the detector by
   construction. It would flag nothing.
2. **The detector is already known to be uninformative.** Measured OOD AUC ≈ **0.47** (chance), and
   `train_model.py:288-290` records the pre-existing decision that this is a property of the
   representation, with "disable the gate" as the anticipated outcome.
3. **§10 removed the premise.** D2 found the reprogramming trajectory sign **flips** between datasets
   (Gill +0.205, GSE242423 −0.214, each clearing its bar by hundredths). A stable "reprogramming is
   out-of-domain, so age reads high" story predicts the *same* sign in both. It is not a domain
   property — it is noise at SNR≈1, which is exactly this stage's thesis.

Even had all three been fine, gating would flag the reprogramming cells — i.e. the entire use case —
making it option C (retreat) wearing a domain-condition costume. **The precision fix is the correct
target; the OOD route is closed.**

## 9.5 R4 TESTED, NOT ARGUED (2026-07-26) — and Step 1 is effectively complete

§9.2 refuted R4 by reading. It has now been refuted by **running**, on the real GSE113957, and two
of Step 1's three items are answered as a by-product. Total compute: ~3 seconds.

### Test 1 — does any cross-sample statistic exist? (the stated mechanism)

```
normalize_counts(X[:5])  vs  normalize_counts(X)[:5]      ->  max |difference| = 0.0
```

A sample's normalised value is **identical** whether or not the other 138 samples are present.
`normalize_counts` divides each row by its **own** library size (`lib = counts.sum(axis=1,
keepdims=True)`), so no cross-sample quantity is ever computed. A train/test split cannot leak
through a function that never looks across samples. `grep StandardScaler|scaler|fit_transform`
over `clock_fit.py` → **no match**. **The scaler R4 describes does not exist.**

### Test 2 — is there a DIFFERENT leak R4 might have been sensing? (group leakage)

The plausible alternative: if GSE113957 contained several samples per donor, `KFold(shuffle=True)`
would split replicates across folds and the CV *would* be optimistic — a real leak, unrelated to
scalers. Measured:

```
143 samples   ->   143 unique `cell id`s   ->   0 donors appearing twice
```

Every sample is a distinct donor. **No group leakage either.**

### Test 3 — reproduce the CV directly (the decisive one)

Ran the exact protocol (`KFold(5, shuffle, seed=0)`, `RidgeCV` refit inside each fold, per-row
normalisation) on the real data:

| | replicated | shipped artefact |
|---|---|---|
| `cv_mae` | **12.67 yr** | 12.27 yr |
| `cv_pearson` | **0.841** | 0.837 |

It reproduces. (The small gap is expected: 143 samples here vs the artefact's 133, and a slightly
different gene set from the newer NCBI annotation.) **`cv_mae = 12.27` is an honest, leak-free
number.** R4 is wrong in its mechanism *and* in its conclusion.

**Why the error is an understandable one:** scaler-before-split is *the* classic CV leak, and the
code does read `Xn = normalize_counts(X)` immediately above `cross_val_predict`. But
**normalisation ≠ standardisation**: per-sample library-size normalisation uses one row's own total,
while per-gene standardisation uses a column statistic across samples. They look alike in code and
leak completely differently. The instinct was right; the identification was not. R4 was correctly
hedged ("possible", "to be audited, not assumed") — the cost of leaving it in would have been an
audit that could not find anything.

### What the same run DID find — R2 confirmed and quantified

| quantity | value | reading |
|---|---|---|
| predicted-vs-true **slope** on held-out folds | **0.717** | **R2 confirmed.** Well below 1.0 → real shrinkage compression, exactly the S1 bar's target (`slope ∈ [0.85, 1.15]`) |
| `alpha` chosen by RidgeCV | **0.272** | near the **bottom** of `logspace(-1, 4)` → the penalty is barely binding → **R1 confirmed** |
| in-sample MAE vs CV MAE | **0.77 vs 12.67 (16×)** | memorisation, not mild over-regularisation |

### C3 pre-emptively eliminated

§9.3 predicted slope recalibration alone would underperform. Tested (recalibration fitted
out-of-fold so it cannot cheat):

| candidate | `cv_mae` | verdict vs §3 bar |
|---|---|---|
| C4 dense ridge (control) | 12.67 yr | FAIL |
| **C3 + slope recalibration** | **12.78 yr (−1%, i.e. no help)** | **FAIL** |

Correcting the slope of a model that is memorising rescales its errors without removing them.
**C3 is answered — it is not a candidate.** Step 2 should run C1 (sparse) and C2 (pre-filtered)
against C4 only.

### Net effect on the plan

| Step 1 item | Status |
|---|---|
| re-measure `cv_mae` leak-free (R4) | ✅ **done — no leak exists; 12.27 is honest** |
| predicted-vs-true slope (R2) | ✅ **done — 0.717, compression confirmed** |
| error profile by age decile | ⏳ still worth running (S4) |

**The plan's thesis is untouched and strengthened.** Removing R4 removes the one route by which the
SNR problem could have been *overstated*: the ruler really is ±12.3 yr, honestly measured. The
remaining work is Step 2 with C1/C2 vs C4.
