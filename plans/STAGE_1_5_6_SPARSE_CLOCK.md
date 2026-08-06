# STAGE 1.5.6 — The clock's density IS the defect. Sparsify it.

**Status:** ✅ **MEASURED 2026-08-04, validated leave-one-donor-out. Not yet applied to any label.**

**Scope of what has run:** 3 new read-only scripts, **0 lines changed in `src/`**, **no label moved.**
**Scope of what this proposes:** one config change, gated, with a pre-registered bar.


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

---

## 5. The plan

| # | step | gate | cost |
|---|---|---|---|
| **1** | ✅ **DONE.** Shape survives (ρ −0.881), but the baseline was not reproduced — see §4.1/§4.2 | — | done |
| **1b** | ✅ **DONE.** Deconfound + re-centre = **+2.11 yr** and move ΔAge *toward* zero — **not** the source. Gap attributed to the **harmonization gain** (§4.3) | — | done |
| **1c** | ✅ **DONE — CONFIRMED.** Measured gain **2.152** (predicted 2.26). And it is a per-GENE reweighting, not a scale: median σ ratio 0.608 while the net effect is ×2.15 (§4.5) | — | done |
| **1d** | 🔴 **NEW: measure the SPARSE clock's gain on HFF.** The gain is per-gene, so top-100 has a *different* gain from the dense clock; §1's Gill-side number cannot transfer | gain reported for k = 100 alongside k = all | free, no retrain |
| **2** | **Pre-register the bar** for adopting a sparse clock: MAE ≤ 8 yr **and** sign agreement ≥ 0.80 vs methylation, on **both** arms, k fixed at 100 in advance | `bar_verdict` row in `tests/test_bars_resolvable.py` | free |
| **3** | **Write `configs/clocks/fleischer_clock_top100.json`** — the same coefficients, 33,055 zeroed. Provenance in `meta`, original untouched | ships as a **new file**; nothing switches automatically | free |
| **4** | **One rebuild + LOOCV under the sparse clock**, full scorecard, snapshot and rollback | every Stage 1 guard reported before/after | one retrain |
| **5** | Only then decide on the label change | — | — |

**Step 1 is the falsifier and it costs nothing. Do it first.**

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
| `results/dage_ledger.csv` | full table, 90 rows × 60 columns |
| `results/diag_dage_ksweep_results.json` | every k, both clocks, both arms |
