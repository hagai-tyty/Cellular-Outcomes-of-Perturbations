# WHERE ΔAge ACTUALLY STANDS — and why "10 % error" is the wrong target

**Status:** 🟢 **AUDIT-3, 2026-08-16.** Read-only. `src/` untouched, nothing withdrawn.
**Conclusion: the work is sound, the model did not "get worse", and the accuracy goal is
mis-specified in a way that is now measured rather than argued.**

---

## 1. The three findings that decide everything, all his

### 🔑 A. ΔAge prediction is CIRCULAR — and this is the most valuable result of the whole arc

`diag_clock_circularity`: of the model's **2 000-gene panel**, **1 956 carry Fleischer clock
weights**, retaining **21.4 %** of the clock's absolute weight mass. The clock's own weights on that
panel reconstruct the label at **ρ 0.96–0.97**, and ridge reproduces that readout at **ρ 0.96–0.99**
on all five evaluable folds.

> **The label is a linear functional of the input.** "Predicting ΔAge" from expression is reading
> back a linear combination of the vector handed in. **The 5.84 yr MAE measures how well the panel
> preserves the clock's direction — not how well anything predicts age.**

**Consequence for a 10 % target: it is achievable and meaningless.** Put more clock genes in the
panel and the error falls toward zero, having learned nothing. **Any same-timepoint ΔAge accuracy
number is a panel-fidelity measurement.**

### 🔴 B. "Use ΔAge to predict later ΔAge" is already dead — he measured it

`diag_early_late_forward`, 6 donors:

| | Pearson | Spearman |
|---|---:|---:|
| early ΔAge → late ΔAge | 0.741 | 0.829 |
| **donor age → late ΔAge** | **0.931** | **0.971** |
| **early ΔAge → late, PARTIALLING OUT donor age** | **−0.064** | — |
| donor age → late, partialling out early | **0.759** | — |

**Early ΔAge predicts late ΔAge only because both are driven by donor chronological age.** Control
for age and the signal is **−0.064** — nothing. **Donor age is known at t = 0 and needs no model.**

The follow-up (`diag_residual_robustness`) tested whether *early expression* predicts the late
residual after donor age: **3 of 9 configurations, "effectively 1 of 9", verdict FRAGILE.**

### C. Regime E fired P0 — and my GSE165177 recommendation was WRONG for a structural reason

I argued GSE165177's 7–8 samples per timepoint would fix `p_unsafe`'s saturation. **It does not, and
the reason is structural, not statistical:**

> *"`p_unsafe` is a fraction of CELLS; a bulk sample is already a population average, so a per-sample
> hard label collapses the fraction to 0/1 before it can be estimated."*

Measured: `unsafe_sd_by_donor` = O1 **0.10**, O2 **0.00**, O3 **0.00**. `P0_void = True`.

**Sample replication cannot recover a per-cell fraction from bulk.** My argument confused replication
at the sample level with resolution at the cell level. **He is right; I was wrong.** What GSE165177
*does* still fix — replicated contemporaneous controls, so no `n = 1` zero-point — stands, and is
why `dage_gse165177` reproduced Gill 2022.

---

## 2. 🔵 The ceiling, and it is not the model

`dage_meth_concordance`:

| | |
|---|---|
| **inter-clock RMS — Horvath skin&blood vs Horvath multi-tissue, same samples** | **9.07 yr** |
| per-clock SD across control groups | 3.36 / 3.41 yr |
| typical HFF day-14 ΔAge after C-7 | ≈ **−6.5 yr** |

> **The two reference instruments disagree with each other by more than the effect being measured.**

A 10 % error target on a −6.5 yr ΔAge is **±0.65 yr**. The gold standards differ by **9.07 yr**.
**You cannot validate a quantity to a precision finer than the disagreement between the instruments
that define it.** That is not a modelling limitation and no architecture, dataset or loss function
moves it.

### And §1's number should be read against that ceiling, not against zero

§1 measured RNA-ΔAge against Horvath multi-tissue at **MAE 5.36 yr** (68 conditions, transient arm,
sparse clock). **5.36 < 9.07.** On the face of it the RNA readout agrees with one methylation clock
*better than the two methylation clocks agree with each other*.

**⚠️ Not yet established like-for-like** — 5.36 is an MAE on 68 conditions against multi-tissue;
9.07 is an RMS on 9 control groups between two clocks. **Different statistic, different sample set,
different pairing.** §0's ERROR 2 was exactly this mistake and it must not be repeated. **But it is
one cheap analysis away from being established, and it is the single highest-value number left.**

---

## 3. So: did the model get worse?

**No — and "worse" is measuring the wrong thing.**

`diag_target_path` reports `median_compression` **0.826 → 0.534** after C-7, verdict `H-SUPPORTED`.
The prediction became *more compressed relative to the label*. That reads as "worse" and is not:

* **C-7 removed a contaminant that inflated HFF's labels ~3×** (fold spread 16.67 → 3.69 yr). The
  labels got **smaller and more honest**. A model that was partly fitting an artefact now has less
  artefact to fit, so its apparent explanatory power drops. **That is the fix working.**
* And his own target-path audit says it directly: *"N3's 'improvement' is not the model improving."*

**Degrading against a contaminated baseline is the expected sign of a real correction.**

---

## 4. What to do — the honest options, none of which drop ΔAge

**ΔAge stays.** The standing constraint holds and nothing here argues otherwise. What changes is
**which claim about it is defensible.**

| | |
|---|---|
| **1. STOP optimising same-timepoint ΔAge accuracy** | it is circular at ρ 0.96 (§1A). Every hour spent lowering that number buys panel fidelity, not biology. **This is the single most important change of direction** |
| **2. Establish the instrument-floor comparison LIKE-FOR-LIKE** | RNA-ΔAge vs each methylation clock, and clock-vs-clock, **same conditions, same statistic, same pairing**. If RNA sits inside the inter-clock envelope, the honest headline is *"agreement with methylation at the limit of methylation's own reproducibility"* — a strong claim, **not a scoped-down one**, and possibly already true |
| **3. Drop "early ΔAge → late ΔAge"** | measured at partial **−0.064** after donor age. Donor age does all the work and is free |
| **4. Keep what demonstrably works** | `fate_roc` **0.983**, within-donor ranking Spearman **0.925–0.983**. Untouched by every ΔAge problem in this document |

### On the 10 % target specifically

**Stated against the RNA clock's own output it is trivial (circular). Stated against methylation it
is below the references' mutual disagreement.** Neither is a bar a correct system can clear, which
is exactly what `REF_GROUND_RULES.md` §5b exists to catch — applied, this once, to a project goal
rather than to a test.

**The reachable version of the same ambition:** *ΔAge agreement with methylation at or inside the
inter-clock RMS of 9.07 yr, demonstrated like-for-like on both reference clocks.* §1's 5.36 yr
suggests that may already be met on one clock. **Item 2 settles it, and it is cheap.**
