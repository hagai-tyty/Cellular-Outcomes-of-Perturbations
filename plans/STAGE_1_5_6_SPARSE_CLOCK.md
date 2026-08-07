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
| **3b** | 🆕 **GATE ON STEP 4 — is the harmonization gain stable across folds?** Compute the clock-weighted gain per fold and test whether it accounts for the 16.67 yr day-14 spread (§4.7). Read-only, **0 lines in `src/`** | `bar_verdict` row in `tests/test_bars_resolvable.py`; pre-registration in §5.1 | free |
| **3b-audit** | 🆕 **Independent audit of 3b — 4 defects + 1 inconsistency (§5.2).** `G_f` is the statistic §4.5 disproved and its exact form is a tautology; the question is settled by elimination so the instrument should be a per-fold RECONSTRUCTION; the ATTRIBUTED branch's remedy reintroduces donor leakage; the gene-set/variance-floor mechanism in §4.7's own table is untested — **and it is the same lever as §4.6's option 2** | — | free, done |
| **4** | **One rebuild + LOOCV under the sparse clock**, full scorecard, snapshot and rollback | every Stage 1 guard reported before/after | one retrain |
| **5** | Only then decide on the label change | — | — |

**Step 1 is the falsifier and it costs nothing. Do it first.**

---

## 5.1 🆕 STEP 3b — PRE-REGISTRATION (written 2026-08-07, before the measurement)

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

### 🔴 A2 — the question is already settled by elimination; a 6-point spread ratio is the wrong instrument

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

### ⚠️ A5 — factual inconsistency

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
| `results/dage_ledger.csv` | full table, 90 rows × 60 columns |
| `results/diag_dage_ksweep_results.json` | every k, both clocks, both arms |
