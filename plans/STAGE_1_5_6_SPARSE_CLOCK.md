# STAGE 1.5.6 — The clock's density IS the defect. Sparsify it.

**Status:** ⚠️ **MEASURED 2026-08-04 and CONFINED. The Gill-side result is validated leave-one-donor-out; it does NOT transfer to HFF — on HFF it inverts (§4.6). Not applied to any label, and it must not be until step 4.**

**Scope of what has run:** 3 new read-only scripts, **0 lines changed in `src/`**, **no label moved.**
**Scope of what this proposes:** one config change, gated, with a pre-registered bar.


> ## 🆕 2026-08-07 — **STEP 3b ADDED, and it GATES step 4**
>
> *Additive. The status line above and every section below are unmodified.*
>
> Reproducing July's HFF result on arm A found that **HFF's day-14 ΔAge swings 16.67 yr across the
> six LOOCV folds** — N2 reads −7.35 where the median is −22.51, a **3.1× compression** — even
> though HFF is never the held-out line and is 99.7% of the age-labelled corpus. §4.3–§4.5's own
> closed form predicts this: `σ_ref` is refit per fold from **five single Gill control cells**.
>
> **Step 4 compares sparse vs dense clock across those same six folds.** Run as written, its paired
> CI would carry a 16.67 yr nuisance term. **Step 3b measures whether the clock-weighted gain
> accounts for that spread, and step 4 does not start until it answers.** Free, read-only.
>
> New: **§4.7** (the finding) and **§5.1** (step 3b's full pre-registration, bar, and decision
> branches). Steps 2 and 3 are unaffected and remain free to run now.

---

## 0. ⚠️ SELF-AUDIT 2026-08-04 — three errors in this document, and one thing that got stronger

*Run before anything here is acted on. Read this before §1.*

### ✅ The core result SURVIVES the hardest test I could put to it

The obvious way a MAE gain can be fake is **shrinkage** — a predictor that collapses toward zero
scores well on MAE while carrying no information. It does not happen here:

| | mean | SD | ratio to truth |
|---|---:|---:|---:|
| TRUTH (methylation ΔAge) | −6.57 | 11.41 | — |
| raw (full clock) | −20.67 | 22.99 | **2.02× over-dispersed** |
| **top100** | −8.19 | **11.82** | **1.04× — matched** |

And the decisive control:

> **A constant-zero predictor scores MAE 8.45. `top100` scores 5.36 — it beats the shrinkage floor.
> The full clock scores 16.61 — it is WORSE THAN PREDICTING NOTHING AT ALL.**

That last line is stronger than anything §1 originally claimed, and it is the honest headline: the
dense clock is not merely imprecise on this data, it is **actively worse than a constant**.

### ❌ ERROR 1 — §1.3's "the Sendai arm agrees, independently" is WRONG

Sendai is scored on **absolute** age, where the intercept does **not** cancel. Sparsifying shrinks
the weighted sum, so every prediction slides toward the clock's intercept **b0 = 72.43**:

| | mean prediction | distance from b0 | SD |
|---|---:|---:|---:|
| raw | 98.86 | 26.42 | 24.83 |
| top100 | 67.71 | **4.73** | **6.50** |
| truth | 28.29 | 44.14 | 14.75 |

`top100` did not move toward truth — it **collapsed onto the intercept**, and its SD fell to 6.50
against truth's 14.75. The MAE "improvement" is that artefact. **This is not independent
corroboration and §1.3 must not be read as such.** The transient arm is unaffected: it scores ΔAge,
where the intercept cancels exactly.

### ❌ ERROR 2 — the "below methylation's own ±7 yr error" comparison is not like-for-like

The **±7 yr** figure (1.5.2 §12-R) is a **donor-level ABSOLUTE-age** error — two donors of true age
53 reading 44.0 and 58.5. **MAE 5.36 is a condition-level ΔAge error.** Different quantities.
It is suggestive of the right order of magnitude; it is **not** a demonstration that the RNA clock
now matches methylation's precision.

### ❌ ERROR 3 — on skin & blood, `top100` is worse than predicting zero

§3 said ordering "stays poor". The sharper and more damning statement:

| skin & blood, transient | SD | ratio | MAE |
|---|---:|---:|---:|
| truth | 11.46 | — | — |
| raw | 22.99 | 2.01 | 17.84 |
| top100 | 11.82 | 1.03 | **8.79** |
| **constant zero** | — | — | **6.83** |

**`top100` (8.79) loses to a constant-zero predictor (6.83) on skin & blood.** The spread is right;
the *ordering* is wrong, and getting the spread right without the ordering is worse than abstaining.
So the sb/mt split is not "one clock is weaker" — **on one clock the sparse ΔAge is harmful.**

### What survives, precisely

| claim | status |
|---|---|
| dense clock carries a −14.10 yr bias, removed at k ≈ 100 | ✅ **holds** (transient, ΔAge, intercept cancels) |
| MAE 16.61 → 5.36 on multi-tissue | ✅ **holds**, and beats the 8.45 zero-floor |
| spread is preserved, not shrunk (ratio 1.04) | ✅ **holds** |
| leave-one-donor-out generalisation | ✅ **holds** (6.70 / 6.84 / 5.45 vs 16.6) |
| the full clock is worse than a constant | ✅ **new, and stronger than the original claim** |
| Sendai corroborates | ❌ **withdrawn** — intercept artefact |
| beats methylation's ±7 yr | ❌ **withdrawn** — not like-for-like |
| skin & blood merely "poor" | ❌ **understated** — it loses to a constant |

**Net: the finding is real and confined to multi-tissue ΔAge on the transient arm.** It is one clock,
one arm, one estimand — not the two-arm agreement §1.3 claimed.

---

## 1. The finding

The Fleischer clock is a dense RidgeCV over **33,155 genes fitted from 133 samples**. Restricting it
to its **~100 largest-|weight| genes** and changing nothing else:

| | MAE vs methylation | bias | ρ | sign agreement |
|---|---:|---:|---:|---:|
| **full clock (33,155 genes)** | **16.61 yr** | **−14.10** | +0.703 | 0.62 |
| **top-100 genes** | **5.36 yr** | **−1.61** | **+0.835** | **0.94** |

*68 conditions, Gill transient arm, Horvath multi-tissue as truth, ΔAge vs ΔAge, replicates averaged.*

**MAE 5.36 yr is below the reference instrument's own donor-level error of ±7 yr** (1.5.2 §12-R: two
donors of identical age 53 read 44.0 and 58.5). The RNA clock now agrees with methylation about as
closely as methylation agrees with itself.

### 1.1 The mechanism, visible in the sweep

| k | 20 | 50 | **100** | 150 | 300 | 1000 | all 33,155 |
|---|---:|---:|---:|---:|---:|---:|---:|
| MAE | 6.19 | 5.99 | **5.36** | 5.42 | 9.99 | 12.33 | **16.61** |
| bias | +5.06 | +3.22 | **−1.61** | −1.99 | −7.69 | −10.15 | **−14.10** |

**The bias crosses zero exactly where MAE bottoms out.** That is not a tuning coincidence — it is the
mechanism. Thousands of near-zero weights each contribute a little drift; summed over 33,155 genes
they become a **−14 yr systematic offset**. Dropping them removes the offset rather than merely
shrinking noise.

### 1.2 It generalises — leave-one-donor-out

| held-out donor | k chosen on the other two | held-out MAE | full clock |
|---|---:|---:|---:|
| O1 | 50 | **6.70** | 16.59 |
| O2 | 100 | **6.84** | 16.48 |
| O3 | 100 | **5.45** | 16.75 |

**k is stable at 50–100 and the ~2.5–3× improvement survives donor-level cross-validation.** The full
clock is worse than every sparse variant for **every donor individually**, so this is not selection.

### 1.3 The Sendai arm agrees, independently

22 conditions, absolute age (GSE165178 has no untreated control, so ΔAge cannot be formed there).
MAE minimum at **k = 50** — 28.52 yr against the full clock's 65.63. Different protocol, partly
different donors, same conclusion about where the optimum lies.

---

## 2. Why nine months of correlation tests never saw this

Every earlier test in this arc scored `ρ_partial`. **Spearman is shift- and scale-free, so a uniform
−14 yr offset is completely invisible to it.** M-2a, 1.5.4 and the first variant sweep all measured
ordering and all reported "weak but present". None of them could see the bias, because the statistic
was blind to it by construction.

**MAE and bias against a real instrument are what exposed it.** The lesson generalises: a
correlation is the wrong summary for a quantity whose *units* are the claim.

---

## 3. What is NOT fixed

**Horvath skin & blood does not come right at any k.** Its MAE improves (17.84 → 6.69 at k=50) but
ordering stays poor throughout (ρ ≤ 0.43, sign agreement 0.41–0.68, sometimes below chance). So the
sb/mt asymmetry recorded in 1.5.2 and 1.5.4 **survives** — sparsification fixes the *offset*, not the
disagreement between the two reference clocks.

**This does not clear M-2a's SPLIT rule.** A variant passing on one clock and not the other is still
a SPLIT. What has changed is that we now know *what kind* of failure it is: not a weak signal, a
biased one on one axis and a mis-ordered one on the other.

---

## 4. 🔑 Why this plausibly fixes HFF too

**The −14 yr bias is a property of the clock's dense weights, not of any dataset.** It is applied to
every cell the pipeline labels — including HFF's **33,613 labels, 99.8 % of the total**. Sparsifying
the clock changes all of them, not just the 90 with paired methylation.

**A concrete, falsifiable prediction:** G-c step 1 measured HFF's trajectory reaching **−24.0 yr at
day 14** under the full clock. If the −14 yr offset is the same artefact, the sparse clock should
read HFF's day-14 at roughly **−10 yr**, with the same ρ ≈ −0.9 trajectory shape.

**That is checkable today, on data already on disk, with no retrain.** It is step 1 below.

### 4.1 ⚠️ RUN 2026-08-04 — the prediction is VOID, and the reason is a bigger finding

| k | ρ(day, ΔAge) | day-14 ΔAge |
|---|---:|---:|
| top50 | −0.905 | −8.36 |
| **top100** | −0.881 | **−10.72** |
| top150 | −0.857 | −11.22 |
| **all 33,155** | −0.857 | **−10.62** |

**The shape survives sparsification** (ρ −0.857 → −0.881) — that was the falsification condition and
it passed. **But the predicted magnitude shift did not occur, because the full clock in this run
already reads −10.62, not the −24.0 the prediction was built on.**

**My script initially reported CONFIRMED. That was wrong** — it checked the result against the
prediction without checking that the baseline the prediction was predicated on still held. Verdict
logic corrected to `BASELINE_NOT_REPRODUCED`.

### 4.2 ✅ STEP 1b RAN — and it REFUTED my own attribution

**I wrote that "about half of HFF's ΔAge magnitude is contributed by pipeline processing". That is
false.** Decomposing the chain measured each step:

| step | day-0 | day-6 | day-14 | contributes at day 14 |
|---|---:|---:|---:|---:|
| S1 clock, absolute age | 78.65 | 78.30 | 68.04 | — |
| S2 control-relative (ΔAge) | −0.00 | −0.35 | **−10.62** | — |
| S3 cell-cycle deconfounded | −0.32 | +0.04 | −8.83 | **+1.79** |
| S4 re-centred = `y_age` | 0.00 | +0.36 | **−8.51** | +0.32 |

**Deconfounding and re-centring contribute +2.11 yr in total, and they move ΔAge TOWARD zero.** They
are not the source of the gap. The recorded shard value of **−24.02** remains **15.51 yr** away.

### 4.3 🔑 The actual source — and Stage 1.5 derived it in closed form nine days ago

**Harmonization was ON in the real build**, and I had ruled it out on bad evidence: no
`configs/data/*.yaml` sets `harmonize: true`, but the build is driven by the runner —
`local_runners/run_multi_local.py:161`, `harmonize=HARMONIZE, harmonize_ref_dataset="gill_bulk"`.

`STAGE_1_5_HARMONIZATION_AUDIT.md` §2 **Group B** already proved what that does:

```
ΔAge = Σ_g (x_pert,g − x_ctrl,g) · sigma_ref,g / (sigma_d,g + EPS) · w_g
```

> **`sigma_d` does not cancel. It survives as a per-dataset multiplicative GAIN, and HFF carries
> `sigma_gill / sigma_hff`.** That is why "batch-immune by construction" was recorded as an
> overstatement: ΔAge is immune to *additive* batch effects, not to *scale* ones.

**Implied gain: 24.02 / 10.62 ≈ 2.26.** A ~2.3× scale factor on HFF is exactly the shape Group B
predicts, and it is consistent with the direct route's SD being about half the pipeline's.

**Not yet confirmed — this is attribution by elimination plus a matching closed form**, not a
measurement of the gain itself. The confirming test is one number: compute
`Σ_g |w_g| · sigma_gill,g / sigma_hff,g` over the clock's genes and check it lands near 2.26.

### 4.4 What this does to the sparse-clock plan

**The two effects are not the same size after all, and they are not independent.**

| | size at day 14 | status |
|---|---|---|
| clock density bias (Gill, ΔAge) | **−14.10 yr** | ✅ measured, LODO-validated |
| deconfound + re-centre | **+2.11 yr** | ✅ measured — small, and the *wrong direction* to explain anything |
| **harmonization gain** | **≈ ×2.26** | ⚠️ attributed, closed form known, **not yet measured** |

A **gain** and a **bias** compose differently: sparsifying the clock changes the weighted sum, and
harmonization then *multiplies* it. **So the sparse clock must be evaluated with harmonization ON**,
which the Gill-side work in §1 did not do — Gill is the reference dataset, so its own gain is ≈ 1 and
the effect is invisible there. **That is a real limitation of §1's result on HFF specifically**, and
it is why step 4's rebuild cannot be skipped.


### 4.5 ✅ STEP 1c — the gain is MEASURED, and it is not a scale factor

| | |
|---|---|
| HFF day-14 ΔAge, clock applied directly | **−9.96 yr** |
| HFF day-14 ΔAge, with the harmonization gain applied | **−21.43 yr** |
| **measured gain** | **2.152** — predicted **2.26** ✅ |
| recorded shard `y_age` | −24.02 |

The attribution holds: harmonization closes **~11.5 of the ~14 yr**. The residual ≈2.6 yr is within
what differing cell subsets (50,241 cells here vs 37,693 in the recorded build) and a deconfounder
refit on harmonized data would move.

#### 🔴 The part that is worse than a gain

**The median σ ratio is 0.608 and the mean is 0.560** — for most genes `sigma_gill < sigma_hff`,
which would *shrink* ΔAge. The effective gain is nonetheless **2.152**.

> **So harmonization is not rescaling ΔAge. It is REWEIGHTING it** — the ratio is applied per gene,
> and the clock's heavy-weight genes happen to sit where `sigma_gill / sigma_hff` is large. A
> majority of genes are damped while a minority are amplified, and the amplified ones carry the
> clock.

That is a stronger statement than Group B's closed form implies on its face. Group B said `sigma_d`
survives as a per-dataset gain; measured, **it survives as a per-GENE reweighting whose net effect on
this clock is ×2.15**, with the median gene pulling the other way.

#### Why this compounds with §1's sparsification rather than adding to it

Sparsifying to the top-100 weights **changes which genes carry ΔAge** — and the gain is *per gene*.
So the two do not compose additively: **a sparse clock has a different harmonization gain from the
dense one**, and neither §1's number (measured with Gill as reference, where the gain is ≈1) nor the
2.152 above transfers to the combination.

**Consequence for the plan, and it is not optional:** step 4's rebuild must measure the sparse
clock's gain on HFF directly. §1's Gill-side result cannot be extrapolated to HFF at all — not
because it is wrong, but because Gill is the reference dataset and is the one place where this effect
is invisible by construction.

### 4.6 🔴 STEP 1d — the sparse clock and harmonization interact ADVERSARIALLY on HFF

| k | direct | harmonized | **gain** | median σ-ratio on kept genes |
|---:|---:|---:|---:|---:|
| 50 | −8.22 | −28.17 | **3.429** | 0.838 |
| **100** | −10.73 | **−29.70** | **2.769** | 0.836 |
| 150 | −11.68 | −30.97 | 2.651 | 0.780 |
| 300 | −15.45 | −32.99 | 2.135 | 0.690 |
| 1000 | −12.78 | −24.83 | 1.944 | 0.608 |
| **all 33,155** | −9.96 | **−21.43** | **2.152** | 0.608 |

**The sparse clock's gain is HIGHER than the dense clock's — 2.769 against 2.152 — and the gain rises
monotonically as k falls.** The reason is visible in the last column: the clock's largest-|weight|
genes sit precisely where `sigma_gill / sigma_hff` is largest (**0.836** among the top 100 against
**0.608** over all genes). **Sparsification concentrates the clock onto exactly the genes
harmonization amplifies most.**

#### What that does to §1's conclusion, on HFF specifically

Under the pipeline as it actually runs, `top100` gives HFF day-14 **−29.70 yr** — **further from
plausible than the dense clock's −21.43.**

> **§1's finding does not transfer to HFF. It inverts.** Sparsifying removes a −14 yr bias on Gill
> and *increases* the magnitude on HFF, because Gill is the harmonization reference (gain ≈ 1) and
> HFF is not. The two effects do not compose — they compound in opposite directions on the two
> datasets.

#### The plan consequence, and it is a real constraint

**A single clock cannot be adopted globally on this evidence.** The same change improves Gill and
degrades HFF, and HFF is 99.8 % of the age labels. The options are now:

| | |
|---|---|
| adopt `top100` **and** disable harmonization | harmonization exists to align bulk with single-cell; removing it is its own pre-registered change with its own guards, not a side effect |
| adopt `top100` **and** re-fit the harmonizer on the sparse gene space | the σ ratio is computed over the *admissible* gene set; restricting that set changes the variance floor and therefore every ratio |
| keep the dense clock for HFF | asymmetric labels across datasets — a new confound of exactly the kind C-I was fixed to remove |
| do neither yet | ✅ **the honest position until step 4 measures the combination end-to-end** |

**This is why step 4's rebuild was never optional, and it is now the gate rather than a formality.**

---

## 4.7 🆕 2026-08-07 — the gain is NOT STABLE ACROSS FOLDS, and that changes what step 4 measures

*Additive. Nothing in §4.1–§4.6 is modified. Measured while reproducing July's HFF result on arm A
(`experiments/repro_hff_signature_armA.py`; notebook section "MATCHED-DATA HFF REPRODUCTION").*

Same script, same build family, **only the held-out Gill donor differs**:

| fold | N2 | N3 | O1 | O2 | Y1 | Y2 |
|---|---|---|---|---|---|---|
| HFF day-14 ΔAge | **−7.35** | −22.12 | −24.02 | −22.89 | −22.05 | −23.87 |
| slope yr/day | **−0.49** | −1.26 | −1.53 | −1.44 | −1.40 | −1.47 |
| harmonized genes | 5026 | 5258 | 5328 | 5304 | 5402 | 5305 |

**day-14 spread 16.67 yr. N2 is a 3.1× compression off a median of −22.51.**

### Why §4.3–§4.5's own closed form predicts this

Stage 1.5's audit replaced the false "batch-immune by construction" claim with

```
ΔAge = Σ_g δ_g · σ_ref,g / (σ_d,g + EPS) · w_g
```

and §4.5 confirmed the gain is a **per-gene reweighting**, not a scale. Two facts make that
fold-dependent:

* `Harmonizer.fit` takes **training control cells only** — *"already excludes the held-out donor"*
  ([harmonize.py:59-60](src/cellfate/data/harmonize.py)) — and `ref_dataset = "gill_bulk"`.
* Every Gill donor's zero-point is **one unreplicated control** (Stage 1.5 audit §5.2).

So `σ_ref` is estimated from **five single control samples**, and *which* five changes every fold.
**HFF's entire label scale is multiplied by a factor estimated from five control cells.** HFF is
42481 of 42605 age-labelled cells in the arm-A build (99.7%; the 99.8% quoted elsewhere in this
document is the earlier 33,688-cell corpus).

A mean-of-ratios does **not** explain it — N2 (σ ratio 0.425) and Y1 (0.410) are similar and only
N2 collapses. That is expected from §4.5: the gain is per-gene and §4.6 showed it concentrates on
the clock's heavy genes. **The right statistic is the clock-weighted gain, which nobody has
computed per fold.** That is step 3b.

### The consequence for step 4, and it is the reason 3b comes first

Step 4 compares **sparse vs dense clock** across a 6-fold LOOCV. If HFF's labels already swing 3×
between those folds for reasons unrelated to the clock, step 4's paired CI is built on differences
contaminated by a 16.67 yr nuisance term. **The retrain would be spent on a confounded comparison.**

It also sharpens §4.6's option table: *"re-fit the harmonizer on the sparse gene space"* was listed
as one of four ways forward. 3b measures whether the harmonizer's **fit protocol** — not its gene
space — is the live problem.

---

## 5. The plan

| # | step | gate | cost |
|---|---|---|---|
| **1** | ✅ **DONE.** Shape survives (ρ −0.881), but the baseline was not reproduced — see §4.1/§4.2 | — | done |
| **1b** | ✅ **DONE.** Deconfound + re-centre = **+2.11 yr** and move ΔAge *toward* zero — **not** the source. Gap attributed to the **harmonization gain** (§4.3) | — | done |
| **1c** | ✅ **DONE — CONFIRMED.** Measured gain **2.152** (predicted 2.26). And it is a per-GENE reweighting, not a scale: median σ ratio 0.608 while the net effect is ×2.15 (§4.5) | — | done |
| **1d** | ✅ **DONE — and it INVERTS §1 on HFF.** top-100 gain **2.769** vs dense **2.152**; harmonized day-14 **−29.70** vs dense **−21.43**. The clock's heavy genes sit where the σ ratio is largest (0.836 vs 0.608). See §4.6 | — | done |
| **2** | **Pre-register the bar** for adopting a sparse clock: MAE ≤ 8 yr **and** sign agreement ≥ 0.80 vs methylation, on **both** arms, k fixed at 100 in advance | `bar_verdict` row in `tests/test_bars_resolvable.py` | free |
| **3** | **Write `configs/clocks/fleischer_clock_top100.json`** — the same coefficients, 33,055 zeroed. Provenance in `meta`, original untouched | ships as a **new file**; nothing switches automatically | free |
| **3b** | 🆕 **GATE ON STEP 4 — what carries the 16.67 yr fold spread in HFF's labels (§4.7)?** RECONSTRUCT `d_f` from each fold's own harmonizer inputs and decompose the residual across three named terms — **T1 mask, T2 variance floor, T3 σ_gill** — which also answers §4.6 **option 2** (same lever). Read-only, **0 lines in `src/`** | `bar_verdict` row in `tests/test_bars_resolvable.py`; pre-registration in **§5.3** (§5.1 superseded by A1, left visible). **Runs only if 3c does not settle it** — §5.8 | free |
| **3b-audit** | 🆕 **Independent audit of 3b — 3 defects stand; A5 withdrawn and A2 downgraded after review (§5.2).** `G_f` is the statistic §4.5 disproved and its exact form is a tautology; the question is settled by elimination so the instrument should be a per-fold RECONSTRUCTION; the ATTRIBUTED branch's remedy reintroduces donor leakage; the gene-set/variance-floor mechanism in §4.7's own table is untested — **and it is the same lever as §4.6's option 2** | — | free, done |
| **3c** | 🆕 **RUNS BEFORE 3b — is the DEGENERATE CONTROL the carrier (§5.7)?** Leave-one-CONTROL-out on the O1 fold's harmonizer, all five, recomputing HFF's day-14 each time, with the four healthy controls as the built-in negative control. Read-only, **0 lines in `src/`** | `bar_verdict` row in `tests/test_bars_resolvable.py`; pre-registration in **§5.8** | free |
| **4** | **One rebuild + LOOCV under the sparse clock**, full scorecard, snapshot and rollback | every Stage 1 guard reported before/after | one retrain |
| **5** | Only then decide on the label change | — | — |

**Step 1 is the falsifier and it costs nothing. Do it first.**

---

## 5.1 🆕 STEP 3b — PRE-REGISTRATION (written 2026-08-07, before the measurement)

> ## ⛔ **2026-08-08 — SUPERSEDED BY §5.3. Do not run the metric below.**
>
> *Left visible on purpose. §5.2's A1 showed this section's primary metric is broken, and the
> broken bar stays in the record rather than being edited away — the same convention §5.2 applies
> to its own withdrawn A5.*
>
> **What is wrong (A1, verified independently):** `G_f` here is a `|w|`-weighted **mean of
> ratios**, which drops the per-gene deltas `δ_g`. §4.5 established that exact failure — median σ
> ratio **0.608** against a net gain of **×2.152**. And substituting the *exact* gain does not
> rescue it: with `G_f = Σδ·r·w / Σδ·w`, then `d_f/G_f ≡ Σδ·w`, **fold-invariant by algebra**, so
> `R = 0` for any data whatsoever. The proxy measures the proxy; the exact form measures nothing.
>
> **What else is wrong:** the ATTRIBUTED branch below prescribes *"fit-once-and-freeze"*, which
> **reintroduces leakage** (A3) — every Gill donor has exactly one control, so no subset is
> in-train for all six folds, and fitting once means the held-out donor's control is in the fit.
> That is precisely what [harmonize.py:59-60](src/cellfate/data/harmonize.py) exists to prevent.
>
> **What was missed:** §4.7 attributed the spread to `σ_ref` alone and never mentioned the
> **variance floor** (A4) — see §5.3's T2 term.
>
> §5.3 replaces the metric with a **reconstruction**. Everything below is the record of the first
> attempt, not an instruction.


**Owner:** this stage. **Cost:** free, read-only. **Scope:** 1 new script, 1 test row,
**0 lines changed in `src/`**, **no label moves under any outcome.**
**Blocking for:** step 4 only. It does not gate steps 2 or 3, which are independent and free.

### The question

> Does the **clock-weighted harmonization gain** account for the 16.67 yr fold-to-fold spread in
> HFF's day-14 ΔAge (§4.7)?

### The statistic

Per fold `f`, over the genes admissible in that fold, with `w_g` the clock's coefficients:

```
G_f  =  Σ_g |w_g| · σ_gill,g^(f) / (σ_hff,g^(f) + EPS)   /   Σ_g |w_g|
```

`G_f` is the closed form of §4.3 restricted to what the clock actually reads — deliberately **not**
the mean σ ratio, which §4.7 already showed does not discriminate.

**Primary metric — residual spread ratio.** Let `d_f` be HFF's day-14 mean ΔAge in fold `f`, and
`spread(x) = (max x − min x) / |median x|`:

```
R  =  spread( d_f / G_f )  /  spread( d_f )
```

`R → 0` the gain explains the instability; `R → 1` it explains none of it.

### The bar

| | |
|---|---|
| **B1 ATTRIBUTED** | `R ≤ 0.35` — the gain removes ≥ 65% of the relative spread |
| **direction** | lower is better |
| **resolvability** | `bar_verdict` **must be run before the measurement**, on a null simulated at this geometry: 6 folds, gains drawn to match the observed `G_f` dispersion, `d_f = G_f × constant + N(0, s)` with `s` the observed per-fold SEM (0.19–0.24 yr). RESOLVABLE requires a correct system to clear `R ≤ 0.35` at `MIN_PASS_RATE = 0.95`. **If the simulation returns UNRESOLVABLE, the bar moves to `usable_bar` BEFORE the run and the move is recorded** — `REF_GROUND_RULES.md` §5b |

**Secondary, reported not graded:** Spearman ρ(`G_f`, `d_f`) across the 6 folds, and the correlation
between `G_f` and the **per-fold arm A − arm B differences already in
`scorecard/gc2_A_keep_hff.json` / `gc2_B_mask_hff.json`.** Reported to bound how much of step 6's
SD 4.808 this could carry. **Not graded, and it cannot settle that question** — settling it needs a
rebuild, which step 3b deliberately does not do.

### Decision branches, fixed in advance

| outcome | what it means | what happens to step 4 |
|---|---|---|
| **ATTRIBUTED** (`R ≤ bar`) | the fit protocol — `σ_ref` from five single controls, refit per fold — is the instability | **Step 4 is BLOCKED.** Fixing the protocol (fit-once-and-freeze vs refit-per-fold) becomes its own pre-registered Change with its own bar and snapshot, per the one-change rule. Only then is a sparse-vs-dense comparison interpretable |
| **NOT ATTRIBUTED** (`R > bar`) | the gain is not the mechanism | **Step 4 proceeds as written**, and §4.7's spread is carried into its interpretation as a stated nuisance term. A new owner is needed for the instability; it is not this stage's |
| **UNRESOLVABLE** even at `usable_bar` | 6 folds cannot decide it | Step 4 proceeds, §4.7 is recorded as an unexplained confound, and the item is handed to `STAGE_6_NEW_DATA_REV.md` §3 beside D2 — the same n = 1 control problem, which is already donor-blocked |

### What step 3b does NOT license

* **It is not a fix.** It measures attribution. Changing `Harmonizer.fit`'s protocol is a separate
  Change with its own bar — this step only decides whether that Change is needed.
* **It does not re-open step 6.** Step 6 (= G-c step 2, snapshots `gc2_A`/`gc2_B`) ran and returned
  INCONCLUSIVE. Re-running it requires its own pre-registration and is **not** authorised here.
* **It does not touch the clock.** `k = 100` stays frozen per §6.
* **It says nothing about whether HFF's labels are CORRECT age** — still open after arms C and D.

### Falsifiability check (the `verify_1a` lesson)

The script must **self-test that it can fail**: fed synthetic folds whose `d_f` is constant while
`G_f` varies, it must return NOT ATTRIBUTED. A branch that never executes is not a check.

---

## 5.2 🆕 2026-08-07 — AUDIT OF STEP 3b (second machine), and two items our own work left unrecorded

*Additive. §5.1 is unmodified — this records what an independent check of it found. Every claim
below was verified against the code or the raw data, not against the document.*

### ⚠️ 2026-08-08 — CORRECTION to this section, after review by the second machine

*Additive. A1–A5 below are left **exactly as written**; two of them are wrong or overstated and are
marked here rather than deleted, so the record stays visible.*

| finding | status after review |
|---|---|
| A1, A3, A4 | ✅ **stand** |
| **A5** | ❌ **WITHDRAWN — my misreading, not an inconsistency in §4.7** |
| **A2** | ⚠️ **DOWNGRADED — the instrument it recommends stands; "settled by elimination" does not** |

#### A5 — withdrawn

§4.7 reads *"the 99.8% quoted elsewhere in this document is the earlier **33,688-cell corpus**"*.
**33,688 is the CORPUS size, not HFF's label count.** §4's 33,613 is HFF's count, and

```
33,613 / 33,688 = 99.777 %  ≈  99.8 %
```

Consistent, and checked. There was no contradiction — I read a sentence about a corpus as a second
claim about HFF's count. The sentence invites the misreading and is worth tightening, but the
document was right and A5 was wrong.

#### A2 — downgraded, and the reason is a repeat of this arc's signature error

A2 argued the deconfounder cannot carry the fold spread *"because step 1b bounded S3+S4 at
+2.11 yr"*. **Step 1b ran on the UNHARMONIZED path.** Its own docstring opens:

> *"Harmonization is OFF (no config sets `harmonize: true`), so for HFF the chain is exactly:"*
> — `experiments/diag_pipeline_decompose.py`

That line was written **before §4.3 established the real build has harmonization ON**
(`run_multi_local.py:45`, `HARMONIZE = True`). So +2.11 yr is the deconfounder's contribution to a
ΔAge roughly **2.15× smaller** than the one it actually sees, and S3's inputs in the real build are
harmonized values. **Using it to bound fold-to-fold variance in a harmonized build is an
extrapolation, not a proof** — the same species of error as ruling out harmonization from the
absence of `harmonize: true` in the YAML (§4.3), committed against the very measurement that
corrected it.

**What survives is the part that matters, and the correction STRENGTHENS it.** The elimination
framing falls. The **instrument** does not — and it is now better motivated, not worse:

> If the deconfounder's fold-to-fold behaviour is **unknown** rather than bounded, then a statistic
> that presumes harmonization is the mechanism is exactly the wrong tool. A **reconstruction**
> separates them by construction: recompute `d_f` from each fold's own harmonizer and check it
> reproduces the recorded `d_f`. Where the harmonization-only reconstruction **tracks** `d_f`,
> harmonization is attributed; **the residue where it fails to track is where the deconfounder — or
> anything else — lives.**

It is a **measurement**, not an assumption that the answer is already known. A2 should be read as
recommending the reconstruction, and **not** as having pre-decided its outcome.

---

### ✅ What checks out

* **`σ_ref` really is estimated from five single control samples.** Independently confirmed from the
  raw matrix: `GSE165176_Log2_RPM_Sendai_reprogramming` carries **124 sample columns**, and exactly
  **6** are day-0 dermal fibroblasts — `N2/N3/O1/O2/Y1/Y2_Fib_Sendai_Exp2`, one per donor, matching
  `sources.py:417`'s definition of `is_control`. Leave-one-donor-out leaves five.
* **The spread is arithmetically right.** −7.35 against −24.02 is **16.67 yr**.
* **`_clock_range` is NOT the fail-open it looks like.** It reads `getattr(clock, "age_range", None)`
  while the JSON stores the range under `meta`. Checked directly: `LinearClock.from_json` lifts
  `meta.age_range` into a real attribute (`aging.py:80`) and returns **(1.0, 96.0)**. The wiring is
  correct.

### 🔴 A1 — `G_f` is the statistic §4.5 disproved, and its exact form is a tautology

§5.1 defines `G_f = Σ_g |w_g|·σ_gill,g/σ_hff,g ÷ Σ_g |w_g|`. §4.7 is right to reject the plain
mean-of-ratios — but this is the same family with different weights. **§4.5's finding was that the
median ratio is 0.608 while the net effect is ×2.152**, and the reason is `δ_g`, the per-gene
perturbation deltas, which `G_f` drops entirely. It also drops the *sign* of `w_g`.

And substituting the exact gain does not repair it. The exact gain is

```
G_f = Σ_g δ_g·ratio_g·w_g  /  Σ_g δ_g·w_g
```

so `d_f / G_f ≡ Σ_g δ_g·w_g`, which is **fold-invariant by algebra** — `R = 0` for any data
whatsoever. The statistic is caught between a proxy that measures the proxy and an exact form that
measures nothing.

### ⚠️ A2 — ~~the question is already settled by elimination~~ — **DOWNGRADED 2026-08-08. The reconstruction instrument stands; the elimination claim does not (see the correction box above)**

Across folds the clock is frozen and HFF's cells are fixed. The only inputs to HFF's ΔAge that can
change are:

| input | can it move across folds? |
|---|---|
| `σ_gill` (5 of 6 donors) | **yes** |
| the admissible gene set (5026–5402) | **yes** |
| the deconfounder coefficients `(a, b)` | yes — **but step 1b bounded the whole S3+S4 contribution at +2.11 yr** |
| `w_g` (clock), HFF `δ_g`, `μ` terms | no — frozen, or cancel in a difference |

**By elimination, harmonization is the only remaining candidate that can move 16.67 yr, and no
statistic is needed to reach that.** What is actually unknown is not *whether* but *by what path*.

> **The decisive test is a reconstruction, not a correlation:** recompute `d_f` from each fold's own
> harmonizer (`σ_gill^(f)`, `σ_hff`, gene set, floor) and check it reproduces the recorded `d_f`.
> Deterministic, read-only, and it cannot return "six folds cannot decide it". Where the
> reconstruction *fails* to track is the informative residue.

### 🔴 A3 — the ATTRIBUTED branch's remedy would reintroduce donor leakage

§5.1's ATTRIBUTED branch names the fix as *"fit-once-and-freeze vs refit-per-fold"*. **Fit-once
means the harmonizer has seen the held-out donor's control sample before that donor's cells are
transformed as test data.** That is precisely what `harmonize.py:59-60` (*"already excludes the
held-out donor"*) and Stage 1.5's Group C leak-safety test exist to prevent. There is no subset of
Gill controls that is in-train for every fold, because every donor is held out in one of them.

**So the branch as written points at a fix that cannot be adopted.** Candidates that do not leak,
none of them named in §5.1:

| | |
|---|---|
| shrink `σ_ref` toward a pooled/global estimate | keeps the per-fold refit, damps a 5-sample estimator |
| change `ref_dataset` away from `gill_bulk` | HFF has thousands of control cells; Gill has six samples. But the clock was fit on bulk, which is why the reference is bulk — this is a real trade, not a free swap |
| more control replicates per Gill donor | correct, and **donor-blocked** — this is Stage 6's D2 |
| accept and report it as a nuisance term | the NOT-ATTRIBUTED branch, applied honestly |

### 🔴 A4 — the gene-space mechanism is in §4.7's own table and goes untested

§4.7 records **harmonized genes 5026–5402 per fold**, then attributes the instability to `σ_ref`
alone. But `harmonize.py:112` computes

```
floor = float(np.median(sigma))      # median over the ADMISSIBLE genes
sigma = np.maximum(sigma, floor)
```

`median(sigma)` is a **set-level** statistic. Change the admissible set and the floor moves for
**every gene at once**, including genes whose own σ did not change. **N2 has the fewest admissible
genes (5026) and is the fold that collapses.** That correlation is sitting in §4.7's table untested,
and it is a different mechanism from "σ_ref estimated from five samples" — it would survive even if
σ_ref were perfectly stable.

### ❌ A5 — factual inconsistency — **WITHDRAWN 2026-08-08, I MISREAD §4.7 (see the correction box above)**

§4 states HFF carries **33,613** labels; §4.7 refers to *"the earlier **33,688**-cell corpus"*. One
of the two is wrong.

---

### 📌 Two items our own 1b–1d work left unrecorded

**U1 — `diag_harmonization_gain.py` (the script behind 1c and 1d) has a fail-open control selector.**

```python
g_day = [float(m.group(1)) if (m := re.search(r"_d(\d+)_", s)) else 0.0 for s in samples]
g_ctrl = g_norm[g_day == 0.0]        # anything UNPARSEABLE silently becomes a control
```

Checked: **118 of 124** sample names carry `_dNN_`, and **not one carries `_d0_`**. The 6 that do not
parse are exactly the fibroblast baselines. **So 1c/1d used the right six rows — but by accident, not
by design.** The recorded numbers stand and the script is deliberately **left unmodified** so those
results stay reproducible; the defect is recorded here instead, and
`experiments/diag_harmonizer_refit_sparse.py` selects on `sources.py:417`'s definition and refuses to
guess.

**U2 — §4.6's option 2 is mechanically narrower than the option table states, and it is the SAME
question as step 3b.**

§4.6 offered *"re-fit the harmonizer on the sparse gene space"* on the reasoning that *"restricting
that set changes the variance floor and therefore every ratio"*. Reading `Harmonizer.fit` pins what
re-fitting can and cannot do:

| quantity | does restricting the gene set move it? |
|---|---|
| `mu_g` (`harmonize.py:110`) | **no** — per-gene |
| `sigma_g` (`harmonize.py:111`) | **no** — per-gene |
| `floor = median(sigma)` (`harmonize.py:112`) | **yes** — set-level |
| the admissible mask (`harmonize.py:88, 91`) | **yes** — drops genes entirely |

> **Option 2 is, mechanically, a change to the variance floor and the admissible mask, and nothing
> else.** If the top-100 clock genes already sit above both floors, re-fitting is a **no-op** and
> option 2 cannot work at all.

**And that is the identical lever as A4.** Option 2 varies the gene set deliberately; 3b's folds vary
it incidentally (5026–5402). **They are one measurement, and running them as two spends two efforts
on one question.**

`experiments/diag_harmonizer_refit_sparse.py` is written to answer it — three floor regimes
(pipeline / clock-space / sparse-refit) with the reconciliation to 1c and 1d built in — and is
**NOT YET RUN.**

---

## 5.3 🆕 STEP 3b — PRE-REGISTRATION v2, THE RECONSTRUCTION (written 2026-08-08, supersedes §5.1)

**Owner:** this stage. **Cost:** free, read-only. **Scope:** 1 script, 1 test row, **0 lines
changed in `src/`**, **no label moves under any outcome.** **Blocking for:** step 4 only.
**Supersedes:** §5.1, which stays visible and must not be run.
**Merges in:** §4.6 **option 2** — see "one lever, not two" below.

### Why a reconstruction and not a statistic

§5.2's A1 killed the scalar-gain metric, and A2's downgrade removed the ground for *assuming*
harmonization is the mechanism: step 1b's +2.11 yr bound was measured on the **unharmonized** path
(`diag_pipeline_decompose.py`'s docstring opens *"Harmonization is OFF"*, written before §4.3
established `HARMONIZE = True` at `run_multi_local.py:45`), so it cannot bound the deconfounder's
fold-to-fold behaviour in the real build.

Those two together decide the instrument. If the deconfounder's fold behaviour is **unknown**
rather than bounded, any statistic that presumes harmonization is the mechanism is the wrong tool.
**A reconstruction separates them by construction:** recompute `d_f` from each fold's own
harmonizer inputs; where it tracks the recorded `d_f`, harmonization is attributed, and **the
residue where it fails to track is where the deconfounder — or anything else — lives.**

### The structural fact that makes the decomposition sharp

**HFF's control cells are fold-invariant.** Every LOOCV fold holds out a *Gill* donor; HFF's
day-0 controls are identical in all six. So `σ_raw^(hff)` is **the same number, per gene, in every
fold**. HFF's side of the ratio can therefore move through exactly two channels — which genes are
in the admissible set, and the floor applied to them — and through nothing else. Any HFF-side
variation **is** the mask or the floor, not σ estimation. That is not an assumption; it follows
from `Harmonizer.fit` taking controls per dataset.

### What is reconstructed

Per fold `f`, from `harmonize.py`'s own chain, with `w_g` the frozen clock coefficients:

```
d̂_f  =  Σ_{g ∈ G^(f)}  δ_g · ( σ_gill,g^(f) / σ_hff,g^(f) ) · w_g

  G^(f)      admissible set: intersection of per-dataset {mean control expr ≥ 0.1}   (:88, :91)
  σ_raw      per-gene std over that dataset's control observations                    (:111)
  floor^(ds,f) = median( σ_raw^(ds,f) ) over G^(f)          ← SET-LEVEL                (:112)
  σ         = max( σ_raw, floor )                            ← clamps ~half the genes  (:113)
  δ_g        mean_day14(x_hff,g) − mean_control(x_hff,g)     ← fold-INVARIANT
```

`d̂_f` reconstructs **S1+S2 only** (harmonized, control-relative). The recorded `d_f` is post-S3+S4.
The difference is therefore the deconfounder-and-re-centring contribution plus reconstruction error
— which is exactly the residue we want named, not a nuisance to be explained away.

### One lever, not two — §4.6 option 2 is folded in here

`mu_g` (:110) and `σ_g` (:111) are **per-gene**: restricting the gene set moves neither. The gene
set enters in exactly two places — the **variance floor** (a set-level median) and the
**admissible mask**. So §4.6's option 2, *"re-fit the harmonizer on the sparse gene space"*, is
**mechanically a change to the variance floor and nothing else** — the same lever the folds vary
*incidentally* (5026–5402 genes) and option 2 varies *deliberately*. Carrying the floor as its own
term answers both in one measurement.

### The three terms, carried explicitly

| term | what varies across folds | mechanism | whose claim |
|---|---|---|---|
| **T1 mask** | `G^(f)` — which genes enter the sum at all | admissible intersection moves with the Gill controls | A4 / §4.6 option 2 |
| **T2 floor** | `floor^(gill,f)`, `floor^(hff,f)` | `median(σ)` is set-level, so it moves for **every** gene at once — including genes whose own σ did not move | **A4** |
| **T3 σ_gill** | `σ_raw^(gill,f)` | estimated from **five single control samples**, and which five changes per fold | §4.7's original claim |

**Ablation ladder:** hold two terms at the `O1` reference and vary one; then all three. Report the
`d̂` spread each induces. **T1 and T2 are not orthogonal** — the mask determines the floor — so the
individual effects need not sum to the total, and the write-up must say so rather than presenting a
clean variance decomposition it does not have. The single-term runs are well-defined counterfactuals
(apply `O1`'s floor to fold `f`'s gene set), and they are labelled as counterfactuals.

### Gate, metric and bar

**G0 — FIDELITY GATE, must pass before `F` is read.** The reconstruction at fold `O1` must agree
with `experiments/diag_harmonizer_refit_sparse.py`'s **regime A (pipeline floor)** to **≤ 0.5 yr**.
Two independent implementations of the same quantity. *Deliberately not* a check against `d_O1`
itself — that would bake in an assumed S3+S4 magnitude, which A2's downgrade says is unknown.
**G0 fails ⇒ the reconstruction is wrong and nothing downstream is read.**

**Primary metric — residual spread ratio.** Over the six folds, `spread(x) = max(x) − min(x)` in
years:

```
F  =  spread( d_f − d̂_f )  /  spread( d_f )                 spread(d_f) = 16.671 yr
```

Not algebraically degenerate: `d̂` is computed from raw harmonizer inputs, never from `d`.

| | bar |
|---|---|
| **ATTRIBUTED** | `F ≤ 0.25` — the reconstruction removes ≥ 75 % of the spread (residual ≤ 4.2 yr) |
| **direction** | lower is better |
| **resolvability** | `bar_verdict` **run before the measurement**, on a correct system simulated at this geometry: 6 folds, `d̂` reproducing `d` up to per-fold noise at the observed SEM (0.19–0.24 yr). If UNRESOLVABLE, the bar moves to `usable_bar` **before** the run and the move is recorded — `REF_GROUND_RULES.md` §5b |

**Secondary, reported not graded:** the T1/T2/T3 ladder; Spearman(`d̂_f`, `d_f`); and the
correlation of `d̂_f` with the per-fold arm A − arm B differences already in
`scorecard/gc2_A_keep_hff.json` / `gc2_B_mask_hff.json`, to **bound** how much of step 6's SD 4.808
this could carry. **It cannot settle that** — settling it needs a rebuild, which 3b does not do.

### Decision branches, fixed in advance

| outcome | reading | step 4 |
|---|---|---|
| **ATTRIBUTED** `F ≤ 0.25` | harmonization carries the spread; the ladder names which term | **BLOCKED.** The fix targets the named term and ships as its own Change with its own bar |
| **PARTIAL** `0.25 < F ≤ 0.60` | harmonization is a major but not sole contributor | **BLOCKED.** The residue gets an owner before step 4 is interpretable |
| **NOT ATTRIBUTED** `F > 0.60` | harmonization is **not** the mechanism | Step 4 proceeds with the spread carried as a **stated nuisance term**. The residue is now the target: the S1→S2→S3→S4 decomposition, for which `experiments/diag_pipeline_decompose.py` already has the machinery and needs only repointing at built shards per fold |

### The remedies, corrected for A3's leakage constraint

**"Fit-once-and-freeze" is struck** — it leaks. Non-leaking candidates, in the order the ladder
would motivate them:

| if the ladder names | candidate remedy | leakage-safe? |
|---|---|---|
| **T2 floor** | a floor that is not a function of the current fold's gene set — a fixed reference set, or a fixed quantile of a pooled estimate | ✅ **yes** — the floor never touches the held-out donor's values |
| **T1 mask** | freeze `G` on a fold-independent admissibility rule | ✅ yes, same reason |
| **T3 σ_gill** | shrink `σ_ref` toward a pooled estimate; or move `ref_dataset` off `gill_bulk` (a real trade — the clock was fitted on bulk); or more control replicates per donor | shrinkage/ref-change ✅; replicates ⛔ **donor-blocked**, `STAGE_6_NEW_DATA_REV.md` §3 / D2 |
| any | accept it as a stated nuisance term and size step 4 around it | ✅ |

**This asymmetry is itself a reason to run the ladder.** A4's mechanisms (T1, T2) have
**leakage-free fixes**; §4.7's original mechanism (T3) largely does not. Which term carries the
variance therefore decides whether this is cheap to fix or donor-blocked.

### Reuse, and one defect inherited with it

**Reuse `experiments/diag_harmonizer_refit_sparse.py`** (304 lines, in the repo, **NOT YET RUN**).
It already selects Gill's controls by `sources.py:417`'s own definition with an asserted count, and
already computes the three floor regimes (A pipeline / B clock-space / C sparse-refit) with
reconciliation back to 1c and 1d. 3b extends it **per fold**; it does not replace it.

> **U1 — carried forward, not fixed.** `diag_harmonization_gain.py` selected Gill's controls with a
> regex whose `else 0.0` default makes any **unparseable** name a control. It happened to be right:
> 118 of 124 names carry `_dNN_`, none carries `_d0_`, and the 6 that fail to parse are exactly the
> day-0 dermal fibroblasts. **1c/1d got the right 6 rows by accident.** That script is left
> unmodified so 1c/1d stay reproducible; the defect is recorded, and 3b must not inherit the
> pattern.
>
> **Consequence for §4.5/§4.6's numbers:** 1c/1d floored at the median over the **clock genes**,
> while the pipeline floors over the **full admissible space** `genes_G` and applies an expression
> floor the clock-gene set never had. So **2.152 and 2.769 are near-pipeline gains, not the
> pipeline's own.** Regime A is what reconciles them. Nothing in §4.5/§4.6 is withdrawn on this —
> it is a stated precision limit on those two numbers, and 3b reports the reconciliation.

### Falsifiability self-test (mandatory, the `verify_1a` lesson)

Fed synthetic folds whose harmonizer inputs vary while `d_f` is held **constant**, the script must
return **NOT ATTRIBUTED** (`F ≈ 1`). And fed folds where `d_f` is generated *from* the reconstruction,
it must return **ATTRIBUTED**. Both branches must execute in the test suite. A branch that never
runs is not a check.

### What §5.3 does NOT license

* **It is not a fix.** It names a term. Changing `Harmonizer.fit` is a separate Change with its own
  bar and snapshot, per the one-change rule.
* **It does not re-open step 6.** Step 6 (= G-c step 2, snapshots `gc2_A`/`gc2_B`) ran and returned
  INCONCLUSIVE. Re-running it needs its own pre-registration and is **not** authorised here.
* **It does not touch the clock.** `k = 100` stays frozen per §6.
* **It says nothing about whether HFF's labels are CORRECT age** — still open after arms C and D.

---

## 5.4 🆕 2026-08-08 — AUDIT OF §5.3, and a MEASUREMENT that reorders its three terms

*Additive. §5.3 is unmodified. Everything below was checked against the code or against an artifact
already on disk — nothing here required a run.*

### ✅ §5.3's structural fact is right, and the guarantee is STRONGER than its citation

§5.3 attributes HFF's fold-invariance to *"`Harmonizer.fit` taking controls per dataset"*. The real
guarantee is upstream and harder:

```python
is_ctrl  = obs["is_control"].to_numpy().astype(bool)
not_test = ~obs["cell_line"].isin(heldout).to_numpy()
keep     = is_ctrl & not_test                      # build_dataset.py:352-354
```

> *"'Training control' = a control cell whose cell_line is NOT held out ... **decidable from
> cell_line alone, before the full split**."* — `fit_harmonizer`'s own docstring

**The train/val/calib split never reaches the harmonizer at all.** So HFF's control set is not
merely "the same donors" — it is **bit-identical** in all six folds, and no cell-level split seed
can perturb it.

**A sharpening §5.3 does not draw.** `admissible[ds]` is computed per dataset from that dataset's
own pooled controls (`harmonize.py:87-88`). HFF's controls are fold-invariant, so **HFF's
admissible set is too** — and `genes_G` is their intersection. Therefore:

> **`G^(f)` moves through Gill's side and nothing else. T1, T2 and T3 all trace back to the same
> five Gill control samples** — they are three channels out of one estimate, not three independent
> sources.

### ✅ `F` is well-conditioned — checked, no issue

§5.1's `spread` divided by `|median|`, which would blow up exactly in the success regime
(`d − d̂ → 0`). §5.3 redefines `spread(x) = max(x) − min(x)` **in years**, so `F → 0` cleanly. The
defect does not carry over.

### 🔴 G0 compares two different quantities

`diag_harmonizer_refit_sparse.py` selects Gill's controls as **every** `_Fib_` sample and asserts
the count — **all six donors, no fold exclusion.** The O1 fold fits `σ_gill` on **five**. So its
regime A is an all-six-donor fit, **not fold O1**, and a ≤ 0.5 yr gate between them can fail for a
reason that has nothing to do with implementation fidelity. Extending that script per fold fixes the
quantity but forfeits the independence, since both sides then run the same code path.

### 🟢 The genuinely independent G0 reference is already on disk

`runs/cellfate_multi/harmonization.json` **is the O1 fold's shipped harmonizer.** Confirmed two
independent ways:

| check | result |
|---|---|
| gene count | **5328** — matches §4.7's O1 column exactly |
| held-out donor | 21 split-map entries marked `test`; **Y1 has only 19**, so the fold is not Y1 |

It carries the pipeline's own per-gene `mu` and **post-floor `sigma`, per dataset**. Reconstructing
against *that* validates T1, T2 and T3's inputs **gene by gene against ground truth**, instead of
one scalar to 0.5 yr — and it is independent in the way G0 wants, because it is the pipeline's
output rather than a second script.

### 🔑 THE MEASUREMENT — T2 is not a second-order term. It IS the transform's centre.

Read directly from that artifact, no run required:

| | |
|---|---|
| floor, `gill_bulk` | **0.15821** |
| floor, `hff_sc` | **0.42388** |
| **floor ratio** `floor_gill / floor_hff` | **0.3732** |
| genes clamped at the floor | **2664 / 5328 = 50.0 %** in *each* dataset — mechanical, since `floor = median(σ)` |
| clamped in **BOTH** datasets | **1848 = 34.7 %**, and their ratio is **exactly 0.3732**, min = max |
| **median ratio over all 5328 genes** | **0.3732** — the median gene's ratio **is** the floor constant |
| median ratio, clamped in neither | 0.5335 (n = 1848) |

> **More than a third of the harmonizer's genes carry one identical ratio, and that ratio is the
> median of the entire distribution.** The variance floor does not nudge the transform at the
> margins — it sets its central tendency.

### Why that reorders §5.3's ladder before the ladder is run

If `floor_gill / floor_hff` shifts between folds, **1848 genes' ratios move in lockstep** — a
coherent, non-averaging perturbation. Per-gene `σ_gill` noise (T3) is estimated from five samples
and is large, but it is *independent across genes* and therefore largely averages out over
thousands of them in a weighted sum.

> **This predicts T2 > T3 as the carrier of the fold spread — and unlike the ladder, it is checkable
> from each fold's two floor scalars alone, before any reconstruction is written.**

**It also raises the stakes on §5.3's own asymmetry argument.** §5.3 notes T1/T2 have leakage-free
fixes and T3 largely does not. If T2 additionally turns out to be the dominant term, then the
instability is both **the cheapest to fix and the one not donor-blocked** — which would move this
off Stage 6's critical path entirely.

**Stated as a prediction with a stated mechanism, not a result.** Nothing here measures the
per-fold floors; only the O1 fold's artifact exists on disk.

### ⚪ Checked and clean — recording so nobody re-checks it

`dataset_summary.json` reports `split_sizes = {train 852, calib 115, val 117, test 21}` against
`n_samples = 42605`, which reads like 41,500 unassigned cells. It is not: those are **split-map
entries**, not cells — 1105 of them, `HFF 981` plus ~21 per Gill donor (`Y1` 19). No defect.

---

## 5.5 🆕 2026-08-08 — §5.4's PRECHECK RUN. Its prediction is NOT supported, and §5.3 is corrected

*Additive. §5.3 and §5.4 are unmodified. This section RAN a measurement; §5.4 proposed it.*
**Artefacts:** `experiments/diag_fold_floor_precheck.py`,
`results/diag_fold_floor_precheck_results.json`. Read-only, `src/` untouched, no label moved.

### First, §5.4's three code findings — all verified independently, all correct

| finding | verified how | verdict |
|---|---|---|
| `fit_harmonizer` is the real guarantee, not `Harmonizer.fit`'s docstring | read `build_dataset.py`: `not_test = ~obs["cell_line"].isin(heldout)`, `keep = is_ctrl & not_test` — decidable from `cell_line` alone, before the split | ✅ **stronger than §5.3's citation.** HFF's control set is **bit-identical**, not merely the same donors |
| HFF's admissible set is fold-invariant too, so `G^(f)` moves through Gill alone | `admissible[ds]` is computed inside the per-dataset loop from that dataset's own pooled controls (`harmonize.py:87-88`) | ✅ **correct, and §5.3 missed it.** T1/T2/T3 are three channels out of **one** five-sample estimate |
| **G0 compares two different quantities** | `diag_harmonizer_refit_sparse.py` selects controls as `no _dNN_ and "_Fib_" in name` — **all six donors, no fold exclusion** | ✅ **a real defect in §5.3**, corrected below |
| §5.4's floor measurement | recomputed from `runs/cellfate_multi/harmonization.json` on this machine | ✅ **reproduces exactly** — floor_gill 0.15821, floor_hff 0.42388, R 0.3732, 1848/5328 = 34.7 % clamped in both with ratio min = max to 12 dp, median ratio over all genes = 0.3732 |

`runs/cellfate_multi` is confirmed as the O1 fold independently of §5.4's split-map argument: its
**5328 genes are unique to O1** among the six folds (N2 5026, N3 5258, O1 5328, O2 5304, Y1 5402,
Y2 5305).

### §5.4's prediction — and the two-scalar test it proposed

> *"if `floor_gill/floor_hff` shifts between folds, 1848 genes' ratios move IN LOCKSTEP … That
> predicts **T2 > T3** … checkable from each fold's two floor scalars alone, before any
> reconstruction is written. I'd do that first — it is two numbers per fold."*

Run:

| fold | genes | floor_gill | floor_hff | **R = ratio** | both % | \|w\|-cover | day-14 |
|---|---|---|---|---|---|---|---|
| **N2** | 5026 | 0.16328 | 0.43154 | **0.3784** | 33.3 % | 0.4681 | **−7.352** |
| N3 | 5258 | 0.15611 | 0.42492 | 0.3674 | 34.8 % | 0.4778 | −22.121 |
| O1 | 5328 | 0.15821 | 0.42388 | 0.3732 | 34.7 % | 0.4825 | −24.023 |
| O2 | 5304 | 0.15603 | 0.42443 | 0.3676 | 34.6 % | 0.4834 | −22.891 |
| **Y1** | 5402 | 0.10334 | 0.42234 | **0.2447** | 34.7 % | 0.4878 | −22.049 |
| Y2 | 5305 | 0.15796 | 0.42408 | 0.3725 | 35.3 % | 0.4827 | −23.869 |

**They are anti-aligned.** The fold that collapses (**N2**) has a floor ratio **1.4 % off O1's**.
The fold with the genuinely anomalous floor ratio (**Y1**, 34 % low) has **completely normal
labels**. Spearman(R, day-14) = **−0.14**.

### The maximum-leverage bound — why two scalars can falsify a mechanism

`d = Σ_g δ_g r_g w_g`. Split into the set `B` clamped in both datasets, whose ratio is exactly `R_f`:

```
d_f  =  R_f · Σ_{g∈B} δ_g w_g   +   Σ_{g∉B} δ_g r_g w_g
```

`d` is affine in `R_f`. Zeroing the second term gives T2 its **largest possible** leverage, and then
`d_f = d_O1 · R_f/R_O1`. Under that ceiling:

| term | Spearman | F (spread surviving) | explains at MAX leverage | worst miss | state |
|---|---|---|---|---|---|
| **T2 variance floor** | −0.14 | **1.398** | **−39.8 %** (worse than not correcting) | 17.00 yr | ⛔ **ELIMINATED** |
| **T1 mask** (clock-weight coverage) | −0.20 | **0.957** | **+4.3 %** | 15.96 yr | ⚠️ **NOT A CARRIER** |

For N2 the max-leverage T2 prediction is **−24.35** against an actual **−7.35** — a 17.00 yr miss on
a 16.671 yr spread. **T2 cannot produce N2's collapse at any leverage fraction.** T1 is not formally
eliminated but removes 4.3 % of the spread at its ceiling; N2 in fact carries *more* top-100 clock
genes than O1 (49 vs 48).

**Verdict: `NO_SCALAR_TERM_IS_A_CARRIER`. §5.4's T2 > T3 prediction is not supported.** The
precheck was exactly the right instrument and it cost two numbers per fold — the prediction it was
built to test simply did not survive it.

### What this does NOT establish — the scope limit is the point

`R_f` and `C_f` are **scalars**. This eliminates the **lockstep-constant** channel and the
**total-coverage** channel. Floor effects acting through *which* genes clamp, and mask effects
acting through gene **identity** rather than summed weight, are **not scalar and are not tested
here**. They live inside the reconstruction. **This narrows the ladder; it does not close it.**

### Where it leaves the ladder — and an inversion worth stating

Within harmonization only **T3** (per-gene `σ_gill`) survives as a scalar-inaccessible channel. But
**T3 is the term §5.4's own averaging argument says should largely cancel** — independent per-gene
noise summed over ~5000 weighted genes. If T1 and T2 are out by measurement and T3 is expected to
average out, the residue points **outside harmonization**, at the deconfounder — which **A2's
downgrade left unbounded** precisely because step 1b measured it on the unharmonized path.

That is an argument, not a result: dropping 1 of 5 samples can move a per-gene σ a long way, and
those moves are **not** independent of gene identity — they are driven by *which* donor is dropped.
So T3 must be measured, not reasoned away. **The reconstruction is now more necessary, not less: it
is the only instrument that separates the surviving harmonization channel from the deconfounder.**

### Corrections to §5.3, carried here rather than by editing it

| §5.3 item | correction |
|---|---|
| **G0** | **Replaced.** Do NOT gate against `diag_harmonizer_refit_sparse.py` regime A — it is an all-six-donor fit and fold O1 fits on five, so a ≤ 0.5 yr gate can fail for reasons unrelated to fidelity. **Use `runs/cellfate_multi/harmonization.json`** — the O1 fold's *shipped* harmonizer, carrying the pipeline's own per-gene `mu` and post-floor `sigma` per dataset. It validates the reconstruction's inputs **gene by gene against ground truth** instead of one scalar, and it is genuinely independent because it was produced by the pipeline, not by the diagnostic |
| **the structural-fact citation** | cite `fit_harmonizer` (`build_dataset.py`), not `Harmonizer.fit`'s docstring; and add that **HFF's admissible set is fold-invariant too**, so `G^(f)` moves through Gill's side alone |
| **the three terms** | they are **not independent sources** — three channels out of one five-sample Gill estimate |
| **the ladder order** | T2 **eliminated**, T1 **not a carrier**, both by this precheck. The reconstruction should spend its effort on **T3 and the residue**, and report T1/T2 only to confirm the precheck at gene level |
| **the branch table** | unchanged in structure. NOT ATTRIBUTED is now the *expected* branch on current evidence, which makes `diag_pipeline_decompose.py`'s S1→S4 machinery — repointed at built shards per fold — the likely next step rather than a fallback |

### Housekeeping, verified and recorded so nobody re-checks it

* `dataset_summary.json`'s `split_sizes {852, 115, 117, 21}` against `n_samples 42605` is **not**
  41,500 unassigned cells — those are split-map **entries** (1105 total: HFF 981, ~21 per Gill
  donor, Y1 19). Confirmed: every fold's `splits/holdout.json` has `map` of length 1105.
  **No defect.**
* `loocv_results/{folds,summary}.json` are **tracked and committed on main.** The copies modified
  on this machine are local run output, not a missing commit.

---

## 5.6 🆕 2026-08-08 — T3 SURVIVES the same test that killed T2, and my averaging argument is DEAD

*Additive. §5.3 and §5.5 unmodified. `experiments/diag_t3_sigma_gill_leverage.py`, read-only, no
HFF stream — Gill's six control samples and the shipped O1 harmonizer are enough.*

### First: my §5.4 prediction was WRONG, and §5.5's precheck is what killed it

§5.4 predicted **T2 > T3** from the floor's coherence. §5.5 measured the floor ratios and they are
**anti-aligned** with the labels — N2, the fold that collapses, sits **1.4 %** off O1; Y1, whose
ratio is genuinely anomalous (34 % low), has entirely normal labels. The max-leverage bound then
eliminated T2 outright: **−24.35 predicted against −7.35 actual, a 17.00 yr miss on a 16.671 yr
spread.**

**The prediction is withdrawn.** It cost two numbers per fold and it died on the first fold it
touched — which is what a precheck is for.

### 🟢 G0 PASSES BIT-EXACTLY, and settles the fold identity by construction

Before any counterfactual: recompute `σ_gill` from the raw Gill matrix with **O1 held out**, floor
it, and compare gene by gene to `runs/cellfate_multi/harmonization.json`.

| | |
|---|---|
| unclamped genes compared | **2664** |
| **median relative error** | **0.0000** |
| p90 relative error | **0.0000** |
| shipped genes found in the Gill matrix | **5328 / 5328** |

**Bit-exact.** Three things follow at once, none of them assumed:

1. **The shipped artifact IS the O1 fold** — proven by reconstruction, not inferred from gene count.
2. **`sources.py:417`'s control definition is the pipeline's** — six `_Fib_` samples, one per donor.
3. **The Gill side of the reconstruction is validated end to end**, so a T3 counterfactual built on
   it is standing on checked ground. This is the G0 §5.3 asked for, and it is stronger than a 0.5 yr
   scalar agreement.

### The T3-only counterfactual — vary `σ_raw^(gill)`, hold mask and floor at O1

Per fold, `rho_g = ratio_g^(f) / ratio_g^(O1)`, and `d̂_f/d̂_O1` is a `δ_g·w_g`-weighted average of
`rho_g`:

| fold | observed day-14 | `d_f/d_O1` | `rho` range on clock genes | `|w|`-mean | observed inside? |
|---|---:|---:|---|---:|:---:|
| **N2** | **−7.35** | **0.306** | [0.108, 1.766] | **0.870** | ✅ |
| Y1 | −22.05 | 0.918 | [0.542, 1.926] | 0.915 | ✅ |
| N3 | −22.12 | 0.921 | [0.144, 2.021] | 0.987 | ✅ |
| O2 | −22.89 | 0.953 | [0.426, 1.999] | 0.995 | ✅ |
| Y2 | −23.87 | 0.994 | [0.427, 2.007] | 0.999 | ✅ |
| O1 | −24.02 | 1.000 | [1.000, 1.000] | 1.000 | ✅ |

**T3 is NOT eliminated.** Unlike T2, every fold's observation lies inside T3's containment interval,
and N2's interval reaches down to **0.108** — far more than the 0.306 it needs.

### 🔑 And the ordering is EXACT

```
Spearman( |w|-weighted rho , observed d_f/d_O1 )  over 6 folds  =  +1.000
```

**T3's leverage rank-orders all six folds exactly as the labels are ordered.** The magnitude is
under-predicted (0.870 against 0.306) — which is precisely what §5.2 A1 says a `|w|`-weighted mean
must do, since it drops `δ_g`. **The ordering statistic is scale-free, so A1's objection does not
reach it.** What the two together say is: right channel, wrong estimator.

### Why my averaging argument was wrong — the mechanism, not just the verdict

§5.4 argued T3 "largely averages out" over ~5000 weighted genes because per-gene `σ` noise is
independent. **It is not independent.** Dropping one of five controls removes **one entire donor
expression profile**, so every gene's `σ` moves in response to **a single shared latent** — which
donor was dropped. The perturbation is therefore **coherent across genes in exactly the way I
attributed to the floor**, and the coherence is donor-specific rather than gene-specific.

That is visible in the table: dropping **N2** moves the clock-weighted `σ_gill` more than dropping
any other donor, and N2 is the fold that collapses. It also explains the **N2 / N3 asymmetry** that
donor age alone could not (§`00_START_HERE.md`, both are donor age 0): in the N3 fold, **N2 is still
in the control set** and still carries the spread; only in the N2 fold is it removed.

### Where this leaves §5.5's inversion

§5.5 argued that with T1 and T2 out and T3 expected to average out, the residue points **outside**
harmonization at the unbounded deconfounder. **The second half of that premise is now measured
false.** T3 does not average out, it survives its max-leverage test, and it orders the folds
perfectly.

> **The reconstruction is still necessary — but its most likely outcome has changed.** The live
> hypothesis is no longer "harmonization is exonerated"; it is **T3, the five-sample `σ_ref`
> estimate**, which is the one term whose leakage-free remedies are shrinkage or a reference change,
> and whose replication remedy is donor-blocked (§5.3's asymmetry table).

**Stated at its true strength: this NARROWS the field to T3 and shows the ordering is right. It does
not establish magnitude** — that needs `δ`, hence the HFF stream, hence the reconstruction.

---

## 5.7 🚨 2026-08-08 — ROOT CAUSE: a Gill CONTROL sample is DEGENERATE in the raw GEO matrix

*Additive. §5.3, §5.5 and §5.6 are unmodified.* **Artefacts:**
`experiments/diag_gill_control_integrity.py`, `results/diag_gill_control_integrity_results.json`.
Read-only, raw GEO only, `src/` untouched.

§5.5 eliminated T2 and T1. §5.6 showed T3 survives with the fold ordering exact. Both were
statistics **over** `sigma_gill`. Neither asked the prior question: **are the six control samples
`sigma_gill` is estimated from actually sound?** They are not.

### The measurement — raw `Log2 RPM`, before any pipeline transform

| control column | min | median | mean | max | **log2 range** | mean−min | implied library |
|---|---|---|---|---|---|---|---|
| **N2_Fib_Sendai_Exp2** | 11.490 | 11.490 | **11.490** | 13.227 | **1.74** | **0.0008** | **1.03e+08** |
| N3_Fib_Sendai_Exp2 | 0.760 | 0.760 | 2.323 | 15.221 | 14.46 | 1.563 | 1.52e+06 |
| O1_Fib_Sendai_Exp2 | 2.218 | 2.218 | 3.182 | 15.313 | 13.10 | 0.964 | 1.66e+06 |
| O2_Fib_Sendai_Exp2 | 1.052 | 1.052 | 2.468 | 14.340 | 13.29 | 1.416 | 1.51e+06 |
| Y1_Fib_Sendai_Exp2 | 2.166 | 2.166 | 2.440 | 16.600 | 14.43 | 0.275 | 1.51e+06 |
| Y2_Fib_Sendai_Exp2 | 0.705 | 0.705 | 2.267 | 14.371 | 13.67 | 1.562 | 1.48e+06 |

**N2's day-0 control is nearly constant.** Its mean sits **0.0008 log2 above its own floor** and its
entire dynamic range is **1.74 log2 units** where every other control spans **13–15**. Real RNA-seq
cannot look like this — the mean is always pulled well above the floor by the highly-expressed
minority. After the pipeline's own inversion `2**x − 1` it implies a library **68× larger** than
every other control, and its `log1p`-CP10k profile has SD **0.011** against ~**0.58**.

Rank agreement with the other five controls: **N2 0.096**, against N3 0.679, O1 0.677, O2 0.685,
Y1 0.503, Y2 0.684.

**Six of the 124 Gill sample columns are degenerate**, and exactly one of them is a control:

    N2_Fib_Sendai_Exp2        ** IS A CONTROL **   range 1.737
    N2_d21_CD13_Sendai_Exp2   treatment            range 7.261
    N3_d21_SSEA4_Sendai_Exp2  treatment            range 2.474
    O2_d40_SSEA4_Sendai_Exp2  treatment            range 9.135
    O2_d9_SSEA4_Sendai_Exp1   treatment            range 2.152
    Y1_d7_CD13_Sendai_Exp1    treatment            range 0.152

### Why one bad control does not stay in its own donor

The day-0 `_Fib_` sample is `is_control` (`sources.py:417`), which makes it **two things at once**:

1. **that donor's entire ΔAge zero-point** — Stage 1.5 audit §5.2's `n = 1` finding, now with a
   concrete instance: N2's 21 ΔAge labels are measured against a constant vector, **in every
   fold**, because a donor's baseline does not depend on which fold is running;
2. **one of the five or six controls `sigma_gill` is estimated from** (`fit_harmonizer`), and
   `sigma_gill/sigma_hff` is the gain applied to **HFF's** labels — 99.7 % of the age-labelled
   corpus — **in every fold that does not hold N2 out.**

A near-constant column sits far from the other donors at most genes, so including it **inflates**
`sigma_gill`, inflating the gain and therefore `|ΔAge|` on HFF. Removing it — which happens in
**exactly one fold, N2's** — deflates it. Direction and identity both match §4.7: five folds read
−22 to −24 and the N2 fold reads **−7.35**.

> ### The uncomfortable consequence, stated plainly
>
> **The N2 fold is the only one whose harmonizer excludes the degenerate control.** On this
> evidence the anomalous-looking number may be the *clean* one, and the five agreeing folds —
> including **O1, which is July's −24.02 reference and the source of every recorded HFF ΔAge** —
> may be the contaminated ones. Agreement among five folds is not corroboration when all five
> share the same contaminant.

### What is ESTABLISHED, and what is NOT

**Established, measured directly from the raw GEO file, no inference:** the six columns above are
degenerate; `N2_Fib_Sendai_Exp2` is one of them and is a control; its rank agreement with the other
five controls is 0.096; all 20 of N2's other samples are normal.

**NOT established:** that this *fully accounts* for the 16.67 yr spread. Leave-one-donor-out on
`sigma_gill` moves the `|w|`-weighted mean by −11.8 % when N2 is dropped but **−20.1 % when Y1 is
dropped**, and Y1's labels are normal — so a scalar σ argument still does not deliver the
magnitude, exactly as §5.6 found. **The reconstruction is still required**, and its question is now
sharper: does removing this specific column reproduce N2's `d/d_O1 = 0.306`?

**Also not established:** whether the defect is in GEO's deposited matrix or in how it is read.
Either way the pipeline consumes it as a control.

### Gate gap

`apply_qc` runs on every fetched chunk and this column survived it. A control that is constant to
within 1.7 log2 units, with a library 68× the cohort, should not be able to reach the harmonizer
silently. Stage 1.5's G-a made `n = 1` **visible**; nothing checks whether that `n = 1` is *sound*.

### Scope — this reaches past Stage 1.5.6

| affected | how |
|---|---|
| N2's 21 Gill ΔAge labels | zero-point is a constant vector, in every fold |
| HFF's labels in 5 of 6 folds | through `sigma_gill` in the harmonizer |
| every step-6 arm (`gc2_A/B/C/D`) | all six folds aggregated, five contaminated |
| July's −24.02 HFF reference | the O1 fold, which includes N2's control |
| 5 further Gill treatment samples | degenerate in their own right, listed above |

**Nothing is withdrawn on the strength of this section.** It names a defect and its reach; the
reconstruction quantifies it.

---

## 5.8 🆕 STEP 3c — PRE-REGISTRATION: is the degenerate control the carrier? (2026-08-08)

**Owner:** this stage. **Cost:** free, read-only. **Scope:** 1 script, 1 test row, **0 lines
changed in `src/`**, **no label moves under any outcome.**
**Runs BEFORE step 3b** and may make it unnecessary. **Blocking for:** step 4, same as 3b.

### Why a new step rather than a rung of 3b's ladder

§5.7 established by direct measurement that `N2_Fib_Sendai_Exp2` — N2's day-0 control — is nearly a
constant vector in the raw GEO matrix, and that it enters `σ_gill` in **every fold that does not
hold N2 out**. That is a **named, specific** candidate. 3b's ladder was designed when the candidate
was a diffuse property of a five-sample estimate; testing a named contaminant is cheaper and
strictly more decisive, so it goes first.

**What §5.7 did NOT establish, and what this step exists to settle:**

> Does removing that column actually reproduce the spread? A scalar σ argument does not deliver the
> magnitude — leave-one-donor-out moves the `|w|`-weighted `σ_gill` by **−11.8 %** dropping N2 but
> **−20.1 %** dropping Y1, and **Y1's labels are normal**. Direction and identity match; magnitude
> is unmeasured.

### The test

On the **O1 fold** — July's reference, `d_O1 = −24.023`, and one of the five folds that *include*
the contaminant — refit the harmonizer **five times**, each time dropping **one** of its five
control samples (N2, N3, O2, Y1, Y2), everything else held exactly: same admissible rule, same
floor rule, same clock, same HFF stream. Recompute HFF's day-14 ΔAge each time.

`MIN_REPLICATES = 3` (`harmonize.py:27`), so four controls remain legal in every arm.

**This requires the HFF stream** — `δ_g` is needed for a magnitude, and §5.6 established that the
mixed-sign weights (2648 +, 2589 −) make every `δ`-free shortcut unbounded. That is the cost, and
it is the same machinery 3b needs, used once instead of eighteen times.

### The built-in negative control — the reason this can fail

Dropping **any** control changes `σ_gill`. The claim is not that N2's removal moves the number; it
is that **N2's removal is an OUTLIER among the five**. The four healthy drops are the negative
control, and they are measured in the same run, by the same code, on the same fold.

### Metrics and bar

Let `d_O1^(−k)` be HFF's day-14 ΔAge with control `k` dropped, and `Δ_k = d_O1^(−k) − d_O1`.

| | |
|---|---|
| **B1 — outlier (PRIMARY)** | `\|Δ_N2\|` is the **largest** of the five **and** `≥ 2 ×` the second largest |
| **B2 — magnitude (PRIMARY)** | gap closed `A = Δ_N2 / (d_N2 − d_O1) ≥ 0.70`, where `d_N2 − d_O1 = +16.671 yr` |
| **B3 — direction (gate)** | `Δ_N2 > 0` — removing the contaminant must move O1's ΔAge **toward zero**, the direction §5.7 predicts. `Δ_N2 < 0` falsifies the mechanism outright |
| **resolvability** | `bar_verdict` run **before** the measurement on a correct system simulated at this geometry: five drops, `Δ` distributed as the healthy drops plus a contaminant term. If UNRESOLVABLE, the bar moves to `usable_bar` **before** the run and the move is recorded — `REF_GROUND_RULES.md` §5b |

**Stated in advance so it cannot be read as a shortfall later:** `O1 minus N2's control` is **not**
the N2 fold. The N2 fold holds out N2 and therefore *includes O1's* control, and fits on a different
admissible set. Exact reproduction of `−7.352` is **not** predicted and is not the bar; B2 asks what
fraction of the gap the contaminant carries.

### Decision branches, fixed in advance

| outcome | reading | consequence |
|---|---|---|
| **B1 ∧ B2 ∧ B3** — ATTRIBUTED | the degenerate control is the carrier | **Step 4 stays BLOCKED** until the control is handled. **Step 3b becomes unnecessary** — its ladder was searching for a mechanism now named. The remedy is a *data* fix, and see the asymmetry note below |
| **B1 ∧ B3, B2 fails** — PARTIAL | N2's removal is special but carries < 70 % of the gap | Step 4 blocked; **3b runs** on the residue, with the contaminant now a known term rather than a hypothesis |
| **B1 fails** — GENERIC | every control drop moves it comparably | The contaminant is **not** the carrier; §5.7 stands as a data defect with its own owner, and **3b runs as written** |
| **B3 fails** — FALSIFIED | removing the contaminant moves ΔAge the wrong way | §5.7's mechanism is wrong. Record it, do not rescue it; **3b runs as written** |

### The asymmetry this could resolve, and why it matters more than the ladder

§5.3's remedy table ranked the fixes by leakage-safety: T1/T2 leakage-free, T3 largely not, and
replication **donor-blocked** behind `STAGE_6_NEW_DATA_REV.md` §3 / D2. A **degenerate input sample**
sits outside that table entirely:

* excluding or repairing one bad column **never touches the held-out donor**, so it is
  **leakage-free**;
* it needs **no new donors**, so it is **not donor-blocked**;
* it is a data-ingest fix, not a change to the estimator.

**If 3c attributes, the instability is the cheapest thing on this page to fix and comes off Stage
6's critical path entirely.** That is why it runs before 3b and before step 4.

### Carried with 3c, but NOT part of its measurement

| item | what | owner |
|---|---|---|
| **3c.2** | Is the defect GEO's deposited matrix or our read of it? Check the GSM's own record and the series matrix against the supplementary file. Cheap, read-only, and it decides whether the fix is *exclude* or *re-read* | this stage, reported beside 3c |
| **3c.3** | `apply_qc` passed a control that is constant to 1.74 log2 with a library 68× the cohort. G-a made `n = 1` **visible**; nothing checks it is **sound**. A guard is a `src/` change and therefore **its own Change with its own bar and snapshot**, per the one-change rule | **not this stage** — named so it is not lost |
| **3c.4** | Five further degenerate Gill **treatment** samples (`N2_d21_CD13`, `N3_d21_SSEA4`, `O2_d40_SSEA4`, `O2_d9_SSEA4`, `Y1_d7_CD13`) | **not this stage** — they do not enter `σ_gill`, but they do enter their own ΔAge and any prior Gill analysis |

### Falsifiability self-test (mandatory)

Fed a synthetic control set in which **no** sample is degenerate, the script must return **GENERIC**.
Fed one where a known column is replaced by a constant vector, it must return **ATTRIBUTED**. Both
branches must execute in the test suite.

### What step 3c does NOT license

* **It is not a fix.** It measures attribution. Excluding or repairing the column changes `y_age`
  and is its own Change, with a rebuild and a restarted guard record.
* **It does not re-open step 6.** `gc2_A/B/C/D` ran and returned INCONCLUSIVE; re-running needs its
  own pre-registration and is **not** authorised here.
* **It does not withdraw July's −24.02**, or any recorded arm result. It measures how much of the
  fold spread one named column carries.
* **It says nothing about whether HFF's labels are CORRECT age** — still open after arms C and D.

---

## 5.9 🆕 2026-08-08 — §5.7 CONFIRMED, and the defect is LARGER: 12 columns, and TWO controls

> ### ⚠️ **RENUMBERED 5.8 → 5.9, and PARTIALLY CORRECTED by §5.10.**
> *Two machines pushed a `## 5.8` concurrently (78fd8a9 and 09bc61f) and neither saw the other's.
> §5.8 is step 3c's pre-registration; this section is renumbered to 5.9. **Only the header number
> changed — every word below is as written.** The step-3b and step-3c table rows point at §5.8,
> which now resolves unambiguously.*
>
> **Two claims below do not reproduce and are corrected in §5.10 — read it before acting on this
> section:**
>
> * **"3.4× margin, no false positives"** — the next SOUND column is `N3_d11_SSEA4` at
>   **0.2967**, not 0.964, so the margin is **1.08×** and there is **no separation**. 0.964 is the
>   next sound *control*, a different population from the 124 the screen ran on.
> * **`Y1_Fib` as a SECOND defective control** — **not established.** Its library (1.51e6) and log2
>   range (14.43) are both **normal**; only `mean−min` is low, on a continuum.
>
> Consequently the **§5.6 supersession below is held at PARTIAL** (it requires `Y1_Fib` to be
> defective) and **the proposed gate statistic is changed** — §5.10 moves it to the library tell,
> which is the half of this section's own argument that does hold.

*Additive. §5.7 unmodified. Census over the raw `GSE165176_Log2_RPM_Sendai_reprogramming` matrix,
before any transform. No run, no HFF stream.*

### ✅ First, a correction I owe on §5.6

§5.7 is right that **the containment interval cannot falsify.** `d̂_f/d̂_O1 = Σ c_g·ρ_g / Σ c_g` with
`c_g = δ_g·r_g·w_g`, and the clock's weights are near-balanced in sign (2648 +, 2589 −), so `c` is
mixed-sign and the average is **not** bounded by `[min ρ, max ρ]`. T2's bound was valid because `d`
is affine in **one scalar**; T3's is not.

`diag_t3_sigma_gill_leverage.py`'s docstring says exactly this — and then its output table prints an
**"inside?"** column reading `YES`, and a verdict reading `T3_STILL_LIVE`. **A stated limit that the
output contradicts is a defect, not a caveat.** "T3 not eliminated" is true only in the sense that
nothing eliminated it, and the summary implied a test had been passed. **The Spearman +1.000 is the
real evidence in §5.6; the interval is not evidence at all.**

### ✅ §5.7's root cause reproduces here exactly

Raw Log2 RPM, the six controls, sorted by `mean − min`:

| control | `mean − min` | range | % of genes at the column min | linear RPM sum |
|---|---:|---:|---:|---:|
| **N2_Fib_Sendai_Exp2** | **0.0008** | **1.74** | **99.7 %** | **1.03e+08** |
| **Y1_Fib_Sendai_Exp2** | **0.2745** | 14.43 | **88.5 %** | 1.51e+06 |
| O1_Fib_Sendai_Exp2 | 0.9639 | 13.10 | 68.6 % | 1.66e+06 |
| O2_Fib_Sendai_Exp2 | 1.4155 | 13.29 | 62.8 % | 1.51e+06 |
| Y2_Fib_Sendai_Exp2 | 1.5625 | 13.67 | 60.6 % | 1.48e+06 |
| N3_Fib_Sendai_Exp2 | 1.5628 | 14.46 | 60.7 % | 1.52e+06 |

**N2_Fib confirmed to the digit.** A library **68×** the cohort is the mechanical tell: these are
*reads per million*, so linear RPM must sum to ≈1e6, and every sound column does (1.48–1.66e6).

**One screen that does NOT work, recorded so it is not tried again:** "fraction of genes at the
column min" flags **all 124 columns** at >50 %, because ~60 % zero-inflation is normal for this
assay and every zero-count gene lands on the log2 floor. `mean − min` is the discriminating
statistic; `%at-min` alone is not.

### 🔴 The defect is larger than §5.7 states — 12 columns, and TWO controls

Screening all 124 at `mean − min < ⅕ × cohort median` (cohort median **1.4196**):

| sample | `mean − min` | range | role |
|---|---:|---:|---|
| **Y1_d7_CD13_Sendai_Exp1** | **0.0000** | **0.15** | ⚠️ **entirely constant**, library 2.15e+09 |
| N3_d21_SSEA4_Sendai_Exp2 | 0.0004 | 2.47 | |
| O2_d9_SSEA4_Sendai_Exp1 | 0.0005 | 2.15 | |
| **N2_Fib_Sendai_Exp2** | **0.0008** | 1.74 | 🔴 **CONTROL** — §5.7's finding |
| N2_d21_CD13_Sendai_Exp2 | 0.0142 | 7.26 | |
| O2_d40_SSEA4_Sendai_Exp2 | 0.0254 | 9.13 | |
| O2_d34_SSEA4_Sendai_Exp2 | 0.0799 | 9.69 | |
| O1_d34_SSEA4_Sendai_Exp2 | 0.0902 | 9.79 | |
| N2_d11_CD13_Sendai_Exp2 | 0.0965 | 9.02 | |
| Y2_d34_SSEA4_Sendai_Exp2 | 0.0995 | 10.05 | |
| O1_d11_CD13_Sendai_Exp2 | 0.1034 | 9.00 | |
| **Y1_Fib_Sendai_Exp2** | **0.2745** | 14.43 | 🔴 **CONTROL** — a second one |

**Clean separation, no overlap:** every flagged column is below **0.284**; the next sound column is
**0.964**. A 3.4× gap.

### 🔑 The second control explains what §5.5 left dangling

§5.5 flagged Y1's floor ratio as *"genuinely anomalous — 34 % low"* with entirely normal labels, and
left it unexplained. **Y1_Fib is the second defective control.** Removing it (the Y1 fold) deflates
`σ_gill` for the same reason removing N2_Fib does.

And §5.6's own table already carried the signature without my seeing it:

| fold | `|w|`-weighted ρ | removes a defective control? |
|---|---:|---|
| **N2** | **0.870** ← lowest | ✅ N2_Fib (`mean−min` 0.0008) |
| **Y1** | **0.915** ← 2nd lowest | ✅ Y1_Fib (`mean−min` 0.2745) |
| N3 | 0.987 | — |
| O2 | 0.995 | — |
| Y2 | 0.999 | — |

> **The two folds whose clock-weighted `σ_gill` drops most are exactly the two folds that remove a
> defective control, ranked in the same order as the severity of the defect.** That is a mechanism,
> not a coincidence — and it means §5.6's Spearman +1.000 was reading the contamination, not
> donor biology. **§5.6's "N2 is the atypical donor" reading is superseded: it was never the donor,
> and it was never donor age 0. It is the sample.**

### 🔴 Ten degenerate NON-control columns are in the build as Gill training labels

The other ten are perturbation samples, and `gill_bulk` is a **training source** — its 124 columns
carry age labels into the corpus. `Y1_d7_CD13_Sendai_Exp1` is **entirely constant** (range 0.15
across ~20k genes). **`apply_qc` passes all of them**: its gates are `min_genes` and
`max_mito_frac`, designed for single cells, and a constant bulk column clears both trivially.

**This does not reach §1's headline.** §1's MAE 16.61 → 5.36 is the **transient** arm
(GSE165177/165178/165179), a different matrix. It does reach anything scored on the Sendai arm.

### The gate this argues for — one statistic, already validated on this cohort

> **`mean − min` in log2 space, per bulk sample, with a floor at ⅕ of the cohort median.** On this
> matrix it flags 12/124 with a 3.4× margin and no false positives. Equivalently: **linear RPM must
> sum to ≈1e6** — sound columns land at 1.48–1.66e6, the four worst at 1.0e8–2.1e9.

Stage 1.5.2's gate G-a made `n = 1` **visible**; nothing checks whether that `n = 1` is **sound**.
That is the hole, and it is one assertion wide.

**Not established here:** that removing these columns reproduces `d/d_O1 = 0.306`. That is step 3c.

---

## 5.10 🆕 2026-08-08 — VERIFICATION of §5.9: the separation claim does not reproduce

*Additive. §5.9's body is unmodified — only its header number changed, and the reason is recorded
in its own banner. Everything below was recomputed here from the raw matrix.*

### ✅ What reproduces exactly

| §5.9 claim | verified |
|---|---|
| 12 of 124 columns flagged at `mean−min < ⅕ × cohort median` (median 1.4196, threshold 0.2839) | ✅ **12/124**, same twelve |
| the `%at-min` screen is useless | ✅ **124/124 columns exceed 50 %** — a genuinely useful negative result, and worth the record so neither machine retries it |
| `N2_Fib` to the digit | ✅ `mean−min` 0.0008, range 1.74, 99.7 % at the column min, library 1.03e+08 |
| ten degenerate NON-control columns are live Gill training labels | ✅ `gill_bulk` is a training source; `Y1_d7_CD13_Sendai_Exp1` is entirely constant (range 0.15, library 2.15e+09) |
| `apply_qc` passes them — `min_genes` / `max_mito_frac` are single-cell gates a constant bulk column clears trivially | ✅ |
| §5.6's containment interval cannot falsify | ✅ agreed, and §5.9's framing of it as a **defect** rather than a caveat is the right call |

### 🔴 What does NOT reproduce — the separation

§5.9: *"Every flagged column is below 0.284; the next sound one is 0.964. A 3.4× gap, no overlap."*

Sorted by `mean−min` over all 124 columns, across the threshold:

| rank | column | `mean−min` | library |
|---|---|---|---|
| 11 | **Y1_Fib_Sendai_Exp2** | 0.2745 | 1.51e+06 |
| — | *threshold* | *0.2839* | |
| 12 | **N3_d11_SSEA4_Sendai_Exp2** | **0.2967** | 8.47e+05 |
| 13 | Y1_d11_CD13_Sendai_Exp1 | 0.3069 | 2.29e+06 |

**The next sound column is 0.2967, not 0.964. The margin is 1.08×, not 3.4×, and there is no gap** —
the threshold falls between two adjacent values 8 % apart. **0.964 is `O1_Fib`, the next sound
CONTROL.** Comparing `Y1_Fib` against the five other controls gives 3.4×; comparing it against the
124-column population the screen actually ran on gives 1.08×. The screen and the separation claim
were computed on different populations.

### 🔴 `Y1_Fib` as a second defective control — NOT ESTABLISHED

Three tells, and only one is even suggestive:

| tell | `N2_Fib` | `Y1_Fib` | sound population |
|---|---|---|---|
| linear library | **1.03e+08** | **1.51e+06** | 8.47e+05 – 2.29e+06 |
| log2 dynamic range | **1.74** | **14.43** | 13.1 – 14.5 |
| `mean−min` | **0.0008** | 0.2745 | continuous from 0.2967 |

**`Y1_Fib`'s library and dynamic range are both entirely normal.** Only `mean−min` is low, and it is
the lowest of six controls — but `n = 6`, and against 124 columns it is unremarkable. Its `max` is
the **highest of any control** (16.600), which reads more like a library dominated by a few
transcripts than a constant column. That may still be a quality concern; it is **not** the defect
`N2_Fib` has, and §5.9's "same direction, milder" does not follow from these numbers.

### ✅ What survives: FIVE columns, ONE control

Flagged by **both** mechanical tells — `mean−min < 0.015` **and** library ≥ 1.69e+07:

| column | `mean−min` | range | library | |
|---|---|---|---|---|
| Y1_d7_CD13_Sendai_Exp1 | 0.0000 | 0.15 | 2.15e+09 | |
| N3_d21_SSEA4_Sendai_Exp2 | 0.0004 | 2.47 | 2.38e+08 | |
| O2_d9_SSEA4_Sendai_Exp1 | 0.0005 | 2.15 | 1.63e+08 | |
| **N2_Fib_Sendai_Exp2** | 0.0008 | 1.74 | 1.03e+08 | **CONTROL** |
| N2_d21_CD13_Sendai_Exp2 | 0.0142 | 7.26 | 1.69e+07 | |

The next column by library is 3.88e+06 — a **4.4× gap**. **That** is the separation §5.9 was
reaching for, and it is real; it just is not where the `mean−min` threshold put it. **Exactly one
of the five is a control, and it is `N2_Fib`** — §5.7 unchanged.

### Consequence for §5.9's §5.6 supersession — held at PARTIAL

§5.9 argues the two lowest `|w|`-weighted ρ folds (N2 0.870, Y1 0.915) are *"exactly the two folds
that remove a defective control, ranked by severity"*. That requires `Y1_Fib` to be defective. It
also rests on Y1 (0.918) versus N3 (0.921) in observed `d/d_O1` — a **0.003 margin**, the same
fragile pair flagged when §5.6's Spearman was checked (drop the O1 anchor and N2 and the ordering
holds at exact `p = 0.042`, hinging on this pair).

**What stands:** §5.6's "N2 is the atypical donor / donor age 0" reading **is** superseded —
`N2_Fib` is a sample defect and that much is established. **What does not:** that Y1 is the second
instance. Recorded as **N2_Fib explains N2; Y1 remains unexplained.**

### The gate — moved to the library tell

§5.9 proposes `mean−min` floored at ⅕ of the cohort median, *"validated on this cohort — 12/124,
3.4× margin, no false positives"*. **That validation does not hold**, and a threshold cutting a
continuous distribution 8 % from its neighbour will flag or miss arbitrarily on a new cohort.

**Use the library instead** — which is §5.9's own stronger half: *"the library is the mechanical
tell — these are reads per million, so a sound column sums to ~1e6."*

> **Gate: assert each bulk sample's linear RPM sum lies within a stated band of 1e6.**
> Mechanically justified by the matrix's own units rather than by a quantile of this cohort, and it
> separates the five degenerate columns from the rest by **4.4×**.

One correction to §5.9's supporting figure: *"sound columns land at 1.48–1.66e6"* is the range of
the six **controls**. Over all 124 columns the sound range is **8.47e+05 – 2.29e+06**, so the band
must be set from the full population, not the controls.

### Step 3c is UNCHANGED, and deliberately

§5.9's covering note asks that 3c test *"removing the TWO defective controls"* jointly. **Declined,
for two reasons:**

1. **3c already answers it.** Its design (§5.8) is leave-one-**CONTROL**-out over **all five** O1-fold
   controls, individually. `Y1_Fib` is one of the five. If it matters, 3c returns it as a **second
   outlier**, B1's *"≥ 2× the second largest"* clause fails, and the run routes to **PARTIAL** — the
   correct branch, reached by measurement.
2. **Joint removal would destroy the separability** that makes the question answerable, and would
   bake in a premise this section shows is unestablished.

**No change to §5.8, and none needed.** The individual design is what lets 3c settle §5.9's claim
rather than assume it.

### Also carried forward

§5.9's reconciliation request stands and 3c settles it: the leave-one-out figures in §5.7
(−11.8 % dropping N2, −20.1 % dropping Y1) are over **all genes, unweighted**, while §5.6's ρ is
**`|w|`-weighted on clock genes** and orders the two the other way. They are different quantities;
3c computes actual ΔAge and therefore adjudicates between them.

---

## 6. What this does not license

* **It is not yet a label change.** Nothing in `src/` or `configs/` has moved.
* **It does not resolve the sb/mt split** (§3).
* **It does not make ΔAge "correct"** — it makes it agree with methylation to within methylation's
  own error on one arm, one clock. That is a large improvement on a 16.6 yr disagreement, and it is
  not the same as validation.
* **k = 100 is chosen from a 12-point sweep.** LODO says the choice generalises, but k should be
  **frozen at 100 before step 4**, not re-tuned on the retrain.

---

## 7. Artefacts

| file | role |
|---|---|
| `experiments/diag_dage_variants.py` | step 1 of the search — 9 variants vs the FACS outcome (all null; ceiling effect at AUC 0.9221) |
| `experiments/diag_dage_variants_meth.py` | the same 9 vs methylation, both arms |
| `experiments/diag_dage_ledger.py` | the per-condition ledger: truth / expected / actual / error |
| `experiments/diag_dage_ksweep.py` | the k-sweep that located the optimum |
| `results/DAGE_LEDGER.md` | the readable record, 68 per-condition rows |
| 🆕 `experiments/repro_hff_signature_armA.py` | the 2026-08-07 measurement that found §4.7's fold instability (matched-fold reproduction is EXACT; R2 fold-stability FAILS) |
| 🆕 `results/repro_hff_signature_armA_results.json` | its per-fold output — the `d_f` values step 3b consumes |
| 🆕 `experiments/diag_fold_floor_precheck.py` | §5.5's two-scalar falsifier: T2 ELIMINATED, T1 not a carrier |
| 🆕 `results/diag_fold_floor_precheck_results.json` | its output |
| `results/dage_ledger.csv` | full table, 90 rows × 60 columns |
| `results/diag_dage_ksweep_results.json` | every k, both clocks, both arms |
