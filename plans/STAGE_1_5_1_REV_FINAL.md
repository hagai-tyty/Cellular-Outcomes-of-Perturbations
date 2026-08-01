# STAGE 1.5.1 — ΔAge anchored to methylation

> ## 🆕 2026-07-31 — **every open item now has an owner. See §11.**
>
> *Additive; nothing below is modified. Since this document was written, Stage 1.5.2 ran and
> answered two of its questions, and one more was answered by direct measurement here.*
>
> | | |
> |---|---|
> | **§6.2 — is the RNA clock calibratable against methylation?** | ✅ **ANSWERED: NO.** Stage 1.5.2 M-2a — SPLIT ⇒ NOT CALIBRATABLE. §8.3/§8.4's optimism about this route is superseded (§11.3) |
> | **§6.3 — are O1/O2 the same physical donors?** | ✅ **ANSWERED: YES**, measured by methylation genotype. §10.6's "not verifiable" row was true of the *metadata*, not the *data* (§11.4) |
> | **§6.5 — how approximate are absolute methylation ages?** | ✅ **QUANTIFIED:** ±7 yr donor-level error (§11.5) |
> | **Genuinely still open** | **three items, all with owners** — §5's retention and HFF's labels need **Stage 6** (more donors / new data); HFF's `age_mask` needs **1.5.2 G-c step 2** (§11.6) |

**Status:** ✅ **EXECUTED and VALIDATED** (2026-07-26). This document is the clean, self-contained
statement of the method, the results, and the one question that remains open. It is written to be
read cold. **§10 validates every claim in it and states how** — including three errors that validation
found and corrected. The three points a reviewer would previously have had to challenge (the
post-hoc promotion of contrast A, the derived intercept, and the power of contrast B) are **closed**
in §4.4, §4.3 and §5 respectively — not merely disclosed.

**Scope:** one public dataset acquired; methylation age computed with two published clocks; three
identity-matched contrasts measured. **`src/` untouched throughout** — no model, training,
calibration or inference code was changed.

> **📁 MOVED 2026-08-01 — the five earlier drafts now live in [`plans/archive/`](archive/README.md).**
> Bytes unchanged (all five SHA-256 hashes verified identical across the move) and `git mv` was used,
> so §10.7's "byte-unmodified" check still holds and history still follows each file. They were
> **archived, not deleted**: they are the audit trail, and nine other documents cite them.
>
> **Earlier drafts** (`STAGE_1_5_1_CLOCK_PRECISION.md` V1, `STAGE_1_5_1_NEW_CHANGES.md`,
> `STAGE_1_5_1_NEW_V2.md`, `STAGE_1_5_1_REVISED.md`, `STAGE_1_5_1_REVISED_REVIEW.md`) are left
> byte-unmodified as the audit trail. An earlier revision of *this* file carried a long list of its
> own defects; those are fixed here rather than catalogued, and git history holds that version.

---

## 1. The problem this stage solved

ΔAge — the model's training target — is `clock(cell) − clock(control)`. Four separate attempts to
repair it on RNA data failed (refit the clock; improve its precision; change the control arm to
non-responders; fix two real statistical errors). The reason, established by measurement:

> **The transcriptomic clock is correctly built and correctly applied, but out of domain on
> reprogramming cells. No RNA-only analysis can detect or repair that, because every RNA route to
> "age" runs through the same clock.**

Evidence that it is out of domain on those cells:

| observation | value |
|---|---|
| non-responders after 11 days of OSKM | **+36.5 yr** — biologically impossible |
| mean responder age, days 7 → 9 → 11 | 78.7 → **101.3** → 77.2 — a +23 yr swing in 2 days, reversed in 2 |
| what it is actually tracking | corr(age, pluripotency) = **−0.62**; corr(age, fibroblast identity) = **+0.62** |
| the two RNA reprogramming datasets | **opposite signs** (Gill +0.205, HFF −0.214) |
| Gill et al. | *"existing transcription clocks failed to accurately predict the age of our negative control samples"* |

Evidence it is *sound in its own domain*, so the failure is domain-specific and not a bug:
GSE113957 through our production path gives **MAE 0.77 yr, ρ 0.99**; Gill's day-0 adults inside the
clock's fitted range give **+18.0 yr for a 21 yr true gap, ρ +0.60**.

---

## 2. Method

**Data.** GSE165179 — Gill's methylation companion to the RNA series we already hold. 96 Illumina
MethylationEPIC samples, beta values, 3 donors: **O1 (53), O2 (53), O3 (38)**.

**Why this dataset and not more RNA.** It supplies the two things no RNA dataset here has:

1. **Identity-matched arms.** `Transiently reprogrammed fibroblast` vs `Negative control fibroblast`
   are *both fibroblasts*, same culture time, differing only in treatment — so the identity confound
   is removed **by design**, not by statistical adjustment.
2. **A real untreated negative control at every timepoint.** GSE165176 contains none at any day > 0;
   GSE165179 has 21 (+12 for the intermediate arms).

**Clocks.** Both of the published clocks Gill reported as rejuvenating: **Horvath skin & blood 2018**
(391 CpGs, **100%** coverage here) and **Horvath multi-tissue 2013** (353 CpGs, **94.6%**). Both are
always reported, whichever way they fall. Coefficients from `biolearn` (Biomarkers of Aging
Consortium), verified against the publications by CpG count, and stored in `configs/clocks/`.

**Age computation.** `age = anti_trafo(Σ wᵢ·βᵢ + intercept)` with Horvath's published transform
(`adult_age = 20`). The coefficient tables ship **no intercept row**, so the intercept is derived
from the three known-age day-0 samples. Because that makes any single-sample check partly circular,
**every conclusion is verified across an intercept sweep** (§4.3) rather than resting on the derived
value.

**Unit of analysis.** The **(donor, reprogramming-length) pair**, matched between arms. `exp1`/`exp2`
replicates of a condition are **averaged**, not treated as independent — they are repeats of one
condition.

**Contrasts.** Three, all identity-matched within a row:

| | treated arm | control arm |
|---|---|---|
| **A** | Transient reprogramming **intermediate** (cells *still* reprogramming) | Negative control intermediate |
| **B** | Transiently reprogrammed **fibroblast** (returned to fibroblast identity — Gill's MPTR claim) | Negative control fibroblast |
| **C** | **Failed** to reprogram fibroblast — **the negative control for the whole design** | Negative control fibroblast |

**Guard.** The clock must reproduce known chronological age on the day-0 samples, judged on the
**consistency** of the derived intercept across donors (not on its mean, which is circular).

**Reproduce with:**

```bash
python experiments/diag_methylation_anchor.py "<dir containing GSE165179>"
```

---

## 3. Results

| contrast | Horvath skin & blood | Horvath multi-tissue |
|---|---|---|
| **A — intermediates** (n=12) | **−24.1** yr, CI [−31.1, −17.0] | **−27.5** yr, CI [−33.7, −21.4] |
| **B — returned fibroblasts** (n=9) | −5.8 yr, CI [−19.5, +7.9] | −9.4 yr, CI [−18.3, −0.5] *fragile* |
| **C — negative control** (n=12) | **+0.5** yr, CI [−2.3, +3.2] | **−2.4** yr, CI [−5.7, +0.8] |
| guard: known-age MAE | 4.0 yr (intercept spread 5.3 yr) | 4.4 yr (spread 6.4 yr) |

**By reprogramming length:**

| | 10 d | 13 d | 15 d | 17 d |
|---|---|---|---|---|
| A intermediates, skin & blood | −14.1 | −19.2 | −27.2 | **−35.7** |
| A intermediates, multi-tissue | −19.8 | −22.2 | −28.4 | **−39.8** |
| B returned, skin & blood | −2.7 | **−14.1** | −0.6 | −5.1 |
| B returned, multi-tissue | −13.4 | **−18.4** | −5.6 | −0.6 |

---

## 4. What the results establish

### 4.1 The negative control is inert — so the design is valid

**+0.5 and −2.4 yr, both CIs containing zero, on two independent clocks.** Cells that received OSKM
and failed to reprogram are indistinguishable from untreated cells.

The transcriptomic clock read **those same cells at +36.5 yr**. That artefact — which drove an entire
earlier plan to redefine the ΔAge control — is **closed**, on evidence rather than argument.

### 4.2 Rejuvenation is real and large during reprogramming

**−24.1 and −27.5 yr**, both clocks, both far from zero, 12 identity-matched pairs. This reproduces
Gill's reported ~30 yr.

**Contrast B (returned fibroblasts) peaks at 13 days and diminishes at 15–17** on both clocks —
precisely the shape Gill describe. Contrast A runs the opposite way (−14 → −36 monotonically with
reprogramming length), which is simply "closer to iPSC". The day-profile therefore distinguishes two
genuinely different quantities, on a priori grounds rather than post hoc.

**Coherent reading:** rejuvenation during reprogramming is large; how much survives the return to
fibroblast identity is a separate question (§5).

### 4.3 The intercept cancels **algebraically** — this is not a robustness argument

Horvath's transform is linear above `adult_age = 20`:
`anti_trafo(x) = 21x + 20`. So for any pair where both samples predict an age above 20,

```
age_treated − age_control = [21(lp_t + k) + 20] − [21(lp_c + k) + 20] = 21·(lp_t − lp_c)
```

**The intercept `k` cancels exactly.** Every contrast in this document is a difference, so none of
them depends on it. Recomputed in the intercept-free form `21·(lp_t − lp_c)`:

| contrast | skin & blood (intercept-free) | multi-tissue (intercept-free) |
|---|---|---|
| **A intermediates** | **−24.5** [−32.2, −16.7] | **−28.3** [−35.4, −21.2] |
| B returned fibroblasts | −6.0 [−20.0, +8.0] | −9.4 [−18.3, −0.5] |
| **C negative control** | **+0.5** [−2.3, +3.2] | **−2.4** [−5.7, +0.8] |

These match §3's derived-intercept values to well under a year, confirming the algebra. A numerical
sweep from −0.60 to +0.70 was also run and changes nothing.

> ### 🆕 ADDED 2026-07-30 — where "exactly" does NOT hold, counted
>
> *A reviewer will ask the obvious question: if the intercept cancels **exactly**, why do §3 and the
> table above differ at all (−24.1 vs −24.5, −27.5 vs −28.3)? The condition is stated above — "both
> samples predict an age above 20" — but the document never said it is **violated**, or how often.
> It is, and here is the full accounting.*
>
> `anti_trafo` is linear **only** at age ≥ 20; below that it is exponential
> (`21·exp(x) − 1`) and the intercept does **not** cancel. Measured on the actual predictions:
>
> | | ages below 20 | contrast A pairs affected | derived-intercept | intercept-free | difference |
> |---|---|---|---|---|---|
> | skin & blood | **4 of 66** | 3 of 12 | −24.05 | −24.46 | **−0.41** |
> | multi-tissue | **4 of 66** | 4 of 12 | −27.55 | −28.33 | **−0.79** |
>
> All eight are deeply-rejuvenated day-15/17 **intermediates** — the cells the effect is largest in,
> which is exactly where a 53- or 38-year-old donor's cells are driven below age 20.
>
> **The algebra is demonstrated, not merely asserted.** Where **no** pair violates the condition, the
> two forms agree to **exactly 0.00**:
>
> | contrast | pairs violating | difference |
> |---|---|---|
> | **C negative control**, both clocks | **0** | **+0.00** — exact |
> | **B returned fibroblasts**, multi-tissue | **0** | **+0.00** — exact |
> | B returned fibroblasts, skin & blood | 1 | −0.19 |
> | A intermediates, skin & blood | 3 | −0.41 |
> | A intermediates, multi-tissue | 4 | −0.79 |
>
> **Zero violations ⇒ exactly zero difference, every time.** That is a stronger confirmation of the
> cancellation than the sweep, because it isolates the one condition the algebra depends on.
>
> **Why no conclusion moves.** The largest deviation is **0.79 yr against an effect of −24 to −28**.
> And every deviation is **negative**, so the intercept-free form reads *larger*: **the §3 values
> reported as headline are the conservative ones.** Using them cannot overstate the effect.
>
> **Precision of the claim:** "the intercept cancels exactly" holds for **every pair in contrast C
> and for contrast B on multi-tissue**, and holds to within 0.79 yr everywhere else. It is not an
> unconditional identity, and this document does not need it to be.
>
> Reproduced independently (pure stdlib, no shared code with the measurement script) by
> `experiments/verify_rev_final_4_4.py`, check **V4**.

**Consequence:** the absent intercept row in the coefficient tables is a non-issue for every claim
made here. It affects only *absolute* ages, which this document does not rely on (§6.5).

### 4.4 Contrast A is corroborated by two things that could not have selected it

Contrast A was named before running but became the headline afterwards (§7). Two independent checks
were then available in the same dataset, and neither could have been produced by that choice:

**(a) An internal negative control specific to A.** The dataset also contains
`Failing to transiently reprogram intermediate` — cells that received OSKM, entered the reprogramming
phase, and *failed*. Against the same `Negative control intermediate` comparator:

| | A: transient-reprogramming intermediates | A-control: **failing** intermediates | paired A − A-control |
|---|---|---|---|
| skin & blood | **−24.5** [−32.2, −16.7] | **−1.1** [−2.7, +0.5] | **−23.3** [−31.1, −15.5] |
| multi-tissue | **−28.3** [−35.4, −21.2] | **−3.6** [−5.1, −2.2] | **−24.7** [−31.3, −18.1] |

Cells that got the same OSKM exposure, the same culture, the same batch and the same timepoints but
**failed** show ≈0 to −3.6 yr; those that succeeded show −24 to −28. **This rules out OSKM exposure
per se, batch, and culture duration as explanations**, because all are held constant. (The small but
real −3.6 yr in failing cells on multi-tissue matches Gill's own report that *"expression of the
reprogramming factors alone was capable of rejuvenating some aspects of the transcriptome."*)

**(b) A dose-response with reprogramming length.**

| | Spearman(reprogramming length, effect) | p | slope |
|---|---|---|---|
| skin & blood | **−0.885** | **0.0001** | −3.30 yr/day |
| multi-tissue | **−0.842** | **0.0006** | −3.15 yr/day |

Longer reprogramming ⇒ deeper rejuvenation, monotonically, at p < 0.001 on both clocks. **A contrast
selected post hoc from noise does not produce a monotonic dose-response**; this is structural
evidence independent of the selection.

**Net:** A now rests on four independent legs — two clocks, an internal negative control, a
dose-response, and an intercept-free formulation. Its post-hoc promotion (§7) remains disclosed, but
it is no longer the only thing supporting the finding.

> ### 🆕 ADDED 2026-07-30 — this section had NO artefact; it now does, and it reproduces
>
> **The gap.** §9 lists `diag_methylation_anchor_results.json` as "full output", and **neither of
> the two checks above is in it.** `diag_methylation_anchor.py`'s `CONTRASTS` list contains no
> failing-**intermediate** arm and no dose-response at all. So the section that carries this
> document's entire answer to the post-hoc-promotion challenge — the thing §7 explicitly leans on —
> rested on ad-hoc computation with nothing a reviewer could re-run. **That is a provenance defect
> even though the numbers turned out to be right.**
>
> **Closed by `experiments/verify_rev_final_4_4.py`**, which re-derives §4.4 from the raw beta
> matrix through an **independent** path — pure stdlib, no numpy, **no shared code** with the
> measurement script — so agreement is corroboration, not the same code answering twice. It first
> reproduces a **known** value (V1) before its new numbers are trusted:
>
> | check | recomputed | this document | |
> |---|---|---|---|
> | **V1** contrast A *(pipeline validation)* | −24.05 [−31.12, −16.98] / −27.55 [−33.69, −21.40] | −24.1 / −27.5 | ✅ |
> | **V2** §4.4(a) failing intermediates | **−1.13 [−2.75, +0.49] / −3.62 [−5.07, −2.16]** | −1.1 / −3.6 | ✅ |
> | **V3** §4.4(b) dose-response | **ρ −0.885, p 0.0001 / ρ −0.842, p 0.0006** | −0.885 / −0.842 | ✅ |
> | **V3** slope | **−3.30 / −3.15 yr/day** | −3.30 / −3.15 | ✅ |
>
> **One convention this section did not state.** The slopes **−3.30 / −3.15 are the intercept-free
> form** (§4.3), not the derived-intercept form, which gives **−3.10 / −2.77**. Both were checked;
> the intercept-free figures are correct as printed and are the ones consistent with §4.3, but the
> table above never said which convention it used. It does now. **ρ and p are identical either way**
> — Spearman is rank-based, and the two forms preserve the ordering of these 12 pairs.

### 4.5 What this means for the project

- **ΔAge is a measurable quantity — *in this dataset*.** Methylation measures real, large
  rejuvenation on these samples with an inert negative control.
  ⚠️ **This is not a statement about the training labels.** An earlier wording, *"ΔAge has a valid
  anchor"*, was quotable as "the labels are anchored" — which §8.2 explicitly denies. The *concept*
  is vindicated here; **the labels are not anchored by this stage** (§6.2, §8.1).
- **The transcriptomic clock's failure is fully localised to the instrument.** The biology is there;
  the RNA clock cannot see it. ΔAge as a *concept* is vindicated.
- **The target definition was never the problem.** No further control redefinition is warranted.

---

## 5. The open question: retention after return to fibroblast identity

**Contrast B.** Current evidence, intercept-free:

| | effect | 95% CI | verdict |
|---|---|---|---|
| skin & blood | −6.0 | [−20.0, +8.0] | NO_EFFECT |
| multi-tissue | −9.4 | [−18.3, −0.5] | REJUVENATION (**fragile**) |

**The two clocks do not disagree about the effect — they differ in precision.** Measured on the
paired differences, and on the untreated-control ages where no effect exists at all:

| quantity (sd, yr) | skin & blood | multi-tissue |
|---|---|---|
| contrast B pairs | 18.2 | **11.6** |
| contrast A pairs | 12.2 | **11.2** |
| failing-intermediate pairs | 2.5 | **2.3** |
| **untreated-control ages** (pure instrument noise) | 6.6 | **5.2** |
| negative-control pairs | **4.3** | 5.1 |

Multi-tissue is tighter on 4 of 5, **including the pure-noise measure**, so it is the somewhat more
precise instrument on these samples. That is a property of the instrument, measurable without
reference to any effect, so noting it is not selection on the answer. The margin is modest and not
uniform, so **neither clock is declared the winner** — both are reported throughout.

**What that implies for power, correctly computed per clock:**

| | sd | MDE at n=9 | observed effect | detectable? |
|---|---|---|---|---|
| skin & blood | 18.2 | **14.0 yr** | −6.0 | ❌ no |
| multi-tissue | 11.6 | **8.9 yr** | −9.4 | ⚠️ **just barely** — hence FRAGILE |

**This corrects an earlier statement in this project's record.** An earlier draft asserted
"MDE ≈13.7 yr, ~17 pairs needed, so n=9 is hopeless." That used the skin & blood spread only. On the
multi-tissue clock **n=9 is already adequate** (MDE 8.9 < the 9.4 observed), which is exactly why one
clock reaches significance and the other does not. **The retention question is at its resolution
boundary, not beyond it.**

**Pooling across reprogramming lengths is well specified.** Testing for day heterogeneity in B gives
p = 0.852 (skin & blood) and p = 0.255 (multi-tissue) — no significant structure — so averaging over
days is not averaging over a varying effect. (Day means are reported in §3 regardless.)

### The honest state of contrast B

Both clocks are consistent with a **real but small retention effect of roughly −6 to −9 yr**. Neither
excludes it; one detects it marginally. **It should not be claimed as established, and it should not
be dismissed.**

**Pre-registered bar for settling it:**

| verdict | condition |
|---|---|
| **RETAINED** | CI excludes 0 and is negative **on both clocks**, in the intercept-free form |
| **NOT RETAINED** | CI includes 0 on both clocks *with adequate MDE for the effect size* |
| **FRAGILE** | any verdict within 0.5 yr of a bound — reported as such, never as a result |

**Resolvability:** to put a −9 yr effect comfortably (not marginally) inside the CI on **both**
clocks requires **≈16 pairs** at the skin & blood spread — roughly double what exists. So the
prerequisite is **more donors**, and the specific target is now quantified rather than guessed.

## 6. What this stage does *not* establish

> ## 🆕 ADDED 2026-07-30 — two OPEN findings inherited from Stage 1.5 that this list omitted
>
> *Additive note. Nothing below this box is modified.*
>
> **Why this box exists.** This section was written to list what the stage does not establish, and
> it was **incomplete**: it omitted two findings that Stage 1.5 surfaced and left open. Every
> statement below the box is true; the *ledger of open work* was not. A reviewer reading this
> document cold would close it believing the harmonization arc was finished. **It is not.** The
> omission is recorded here rather than quietly patched into the list.
>
> | | finding (from `STAGE_1_5_HARMONIZATION_AUDIT.md` §5) | status |
> |---|---|---|
> | **D1** | zero-point is cross-batch — all 6 Gill baselines are `Exp2`, ~50% of samples are `Exp1` | ✅ **measured and downgraded** — paired offset **−2.99 yr, 95% CI [−13.12, +7.14]**, `NO_BATCH_EFFECT`, n=12. Structurally true; **not** demonstrated to drive the ±12.7 yr offset. *Not over-read:* the CI half-width (~10 yr) excludes a **large** batch effect, not a meaningful one |
> | **D2** | **every Gill donor's zero-point rests on ONE unreplicated control sample** | ⚠️ **MEASURED, INDETERMINATE** — `diag_zero_point.py` M3: the baseline explains **56% of the offset variance, 95% CI [9%, 100%]**. That CI spans nearly the whole range, so **at n = 6 it is unresolvable — not unmeasured.** Closing it needs more donors, not more analysis. What *is* still open is the code half: `_control_baseline` still records neither count nor composition, so `n=1` is silent |
> | **D3** | donor chronological age is **parsed nowhere in `src/`** though GEO declares it (N2/N3 = 0, Y1 = 29, Y2 = 35, O1/O2 = 53) | ⚠️ **the question is ANSWERED — M1 FAIL** (extreme contrast 11.8 yr vs a 20.2 yr bar on a true 53-yr gap, power 0.996). Only the **wiring** is open |
>
> ### Why this stage's conclusions survive D2 anyway — by design, not by luck
>
> **Every contrast in §3 is a paired arm comparison** — treated vs control, *same donor, same day*,
> on methylation. **None of them touches the RNA day-0 baseline.** So the `n=1` baseline cannot
> propagate into §3, §4 or §5, and the results stand as written.
>
> **That immunity is not a fix.** D2 is still in the pipeline, and it matters more than it looks:
> with D1 measured small and the clock convicted in §1, **`n=1` is now one of only two live
> explanations for the ±12.7 yr per-donor offset — the offset Stage 2's entire premise rests on.**
>
> A related trap this stage already walked into once and documented: an earlier decomposition here
> (`−28.3 = +8.2 − 36.5`) used two terms sharing that same `n=1` baseline, with
> `corr(baseline, ΔAge) = −0.986`. It was redone baseline-free and the conclusion survived — but it
> shows the failure mode is live, not theoretical.
>
> ### D3 is *unwired*, not *unknown* — and this stage proves it
>
> This document **used** those donor ages: the §2 guard is *"reproduce known chronological age on the
> day-0 samples"*, returning **MAE 4.0 / 4.4 yr**, and the derived intercept comes from the three
> known-age day-0 samples. So the values exist, parse, and are accurate enough to be useful. They
> were read by a **diagnostic script**; `src/` still ignores them. That is the only remaining gap.
>
> **Consequence for sequencing:** `STAGE_1_5_2_LABEL_ANCHOR.md` is gated behind the two *closeable*
> halves only (its §0: make the `n=1` baseline visible, and wire donor age into `src/`). **D2's
> scientific question is explicitly NOT a gate** — gating on something unresolvable at n=6 would
> block that stage forever. An earlier version of this box and of that §0 said "gated behind D2 and
> D3", which was wrong on both counts and is corrected here.
>
> **And it cuts the other way:** `diag_zero_point`'s recorded decision was **ESCALATE** —
> *"the clock does not separate the age extremes on this data, so ΔAge's target is unvalidated…
> Stage 2's premise is void as stated."* **M1's failure is exactly what the methylation anchor is
> meant to resolve**, so 1.5.2 is the response to that escalation, not something queued behind it.

1. **Retention after return to fibroblast identity** — §5, open.
2. **Any anchor for the training labels — HFF *or* Gill.** GSE165176 (RNA) and GSE165179
   (methylation) share **zero samples** (§8.1) — though **GSE165178 does pair to our RNA sample-for-sample** (§8.3), so this limitation is specific to GSE165179, not general. HFF (**~99.8%** of the age-labelled cells — 33,613 of 33,688 training cells)
   additionally has no methylation data at all. Neither can
   be anchored from this dataset, and **calibrating the RNA clock against methylation is not possible
   here**: the RNA labels do not distinguish "still reprogramming" from "returned to fibroblast"
   (the two differ by −24 vs −6, so the mapping would decide the answer), the RNA series has no
   untreated control at day > 0 so the contrasts differ, and donor/day overlap is only 2 × 2.
   **Extending the age target to HFF requires methylation for HFF — a Stage 6 acquisition, not an
   analysis choice.**
3. **Donor identity across the two series.** O1/O2 match *by label and age* between GSE165176 and
   GSE165179; that they are the same physical donors is assumed, not verified.
4. **The neonatal out-of-range limit.** GSE113957 has no samples below age 1, so N2/N3/HFF absolute
   ages remain unusable. A data limit, unaffected by anything here.
5. **Absolute methylation ages.** The derived intercept makes absolute values approximate; only
   *differences* are relied on, and only those shown intercept-robust.

---

## 7. Provenance — pre-registered vs discovered

Recorded so the epistemics can be judged rather than reconstructed.

| item | provenance |
|---|---|
| the guard, the three verdict types, the FRAGILE rule, the (donor, day) unit | **pre-registered** before any age was computed |
| running **both** clocks | adopted from Gill's published methods after reading them; not pre-registered. Both always reported |
| **contrast A (intermediates)** | named before running, but as a secondary read on the identity artefact — **it became the headline after the fact** |
| the intercept sweep | added mid-run, when the missing intercept row was found |
| replicate averaging | a **bug fix**: the first run required unique samples per (donor, day, arm) and silently dropped 6 of 9 pairs |

**Nothing was selected on significance** — all three contrasts and both clocks appear in every table,
including the ones that undercut the retention claim.

**The honest caveat on contrast A, and what now answers it.** It was promoted from secondary to
headline after the fact — that stands and is disclosed. But §4.4 then added two checks that the
promotion could not have manufactured: an **internal negative control** (failing intermediates,
−1.1 / −3.6 yr against A's −24.5 / −28.3, holding OSKM exposure, batch and culture constant) and a
**monotonic dose-response** (p = 0.0001 / 0.0006). Together with two clocks agreeing and the
intercept-free formulation, A rests on four independent legs. A reviewer may still prefer
confirmation on new donors; the finding no longer depends on the selection alone.

---

## 8. Next step — and what is *not* executable here

### 8.1 ⚠️ Nothing in this stage remains to be run, and no label change follows from it

An earlier version of this section said *"use methylation as the ΔAge source where methylation
exists."* **That is not executable, and the correction matters.** Checked directly:

```
GSE165176 (RNA)  124 samples, e.g. N2_d11_CD13_Sendai_Exp1
GSE165179 (meth)  96 samples, e.g. O1_negative_control_15days_exp1
sample-title overlap: 0
```

**GSE165176 and GSE165179 share no samples.** They are separate experiments — different sample sets,
different donor rosters (RNA: N2, N3, Y1, Y2, O1, O2; methylation: O1, O2, O3), different day grids
(7–47 vs 10–17) and different arm vocabularies. **There is no join key *between these two series*,
so no methylation age from GSE165179 can be attached to a cell the model trains on.**

> ⚠️ **Scope, added 2026-07-30 — do not quote the sentence above as a general claim.** It is true of
> **GSE165179 only**. **GSE165178 joins our RNA training data 22/22** (§8.3). An earlier version of
> this sentence was unqualified and would have been quoted against §8.3.

**Consequence, stated plainly:** the project's ΔAge labels are **unchanged** by this stage and remain
RNA-derived from a clock this stage proved is out of domain on reprogramming cells. This stage
delivered *knowledge*, not labels.

### 8.2 What it did deliver

| | |
|---|---|
| the +36.5 yr non-responder artefact | **closed** — inert on methylation against a real control |
| is rejuvenation real? | **yes** — −24 to −28 yr, two clocks, dose-response, internal negative control |
| why did the RNA route fail? | **localised to the instrument**, not the biology and not the target definition |
| are the ΔAge labels now fixed? | **no** — and this stage cannot fix them |

That is the honest ledger: the *concept* of ΔAge is vindicated; the *labels* are not repaired.

### 8.3 🔴 CORRECTION — the paired dataset exists, and it is a free download

An earlier version of this subsection said the open questions needed *"new profiling, since no
public series pairs methylation to GSE165176 or to HFF."* **That was wrong.** Both of our series are
SubSeries of SuperSeries **GSE165180**, which has **four** parts:

| accession | contents | have it? |
|---|---|---|
| **GSE165176** | `[Sendai_RNAseq]` — the RNA we train on | ✅ |
| **GSE165177** | `[Transient_RNAseq]` | ❌ |
| **GSE165178** | **`[Sendai array]` — methylation on the SAME Sendai samples** | ❌ **← get this** |
| **GSE165179** | `[Transient array]` — methylation, §3's results | ✅ |

**GSE165178 pairs to our training data sample-for-sample.** Verified against the real titles:

* 22 methylation samples, titles `{donor}_{day}_{marker}` (e.g. `Y2_d11_SSEA4`);
* our RNA titles are the same key plus a batch suffix (`Y2_d11_SSEA4_Sendai_Exp1`);
* **22/22 join on `donor_day_marker`, zero unmatched**;
* donors **O1, O2, Y1, Y2** (4 of our 6 — the two missing are the neonatal N2/N3, which are out of
  the clock's fitted range anyway), days 9/11/15;
* and the **sort marker *is* the arm label** in our data — `CD13` → *Failing to reprogram
  fibroblast* (47), `SSEA4` → *Reprogramming fibroblast* (65). So the arm assignment transfers
  unambiguously, which is exactly what §6.2 said was impossible for GSE165179.

**What this unlocks, and it is the gate everything else was waiting on:**

1. **M-2 becomes well-defined and adequately powered.** §6.2 withdrew the RNA↔methylation agreement
   test because the arms could not be mapped and overlap was 2 donors × 2 days. On GSE165178 the
   arms map exactly and there are 22 paired samples across 4 donors. **The withdrawal applies to
   GSE165179 only, not to this series.**
2. **A direct anchor for Gill's RNA labels** on those 22 samples — a methylation age attached to a
   sample the model actually trains on.
3. **The calibration route (old Step 3a) is back on the table** — learn the RNA→true-age correction
   on the paired samples, validate leave-one-donor-out. Whether it *generalises to HFF* remains a
   separate question (§6), since HFF is a different cell system, but it can now at least be
   attempted and measured rather than assumed impossible.

### 8.4 The next step, concretely

**Download GSE165178** (series matrix + processed beta matrix — check the format first, as with the
others). Then run the M-2 test that has been blocked all along:

> **Does the transcriptomic ΔAge agree with the methylation ΔAge on the same 22 samples?**
> Agreement ⇒ the RNA clock is calibratable and ΔAge is recoverable for Gill.
> Disagreement ⇒ localises exactly where the RNA clock fails, with paired ground truth.

Either answer is decisive, and it needs no new experiments. **GSE165177** is worth taking at the same
time — it pairs with GSE165179 and would extend the same comparison to the transient arm.

**Bars must be pre-registered before this runs** (ground rule §5b), including the resolvability check
at n = 22 paired samples / 4 donors.

### 8.5 What still needs new data

* **HFF's labels.** No methylation exists for HFF in any series; it is a different cell system, so
  GSE165178 anchors Gill only. Whether a Gill-trained correction transfers to HFF is testable but
  not assumable.
* **§5's retention question** still wants ≈16 pairs; GSE165178 does not address it (it is the Sendai
  arm, not the transient one).

## 9. Artefacts

| file | role |
|---|---|
| `experiments/diag_methylation_anchor.py` | the measurement; read-only; pure verdict logic separated from I/O |
| `tests/test_diag_methylation_anchor.py` | 28 tests — every verdict branch, Horvath's transform against its published fixed points, the replicate-averaging regression |
| `configs/clocks/horvath_skin_blood_2018.json`, `horvath_multitissue_2013.json` | clock coefficients, with provenance in `meta` |
| `diag_methylation_anchor_results.json` | full output of the **three** contrasts in §3, per-pair values, and the intercept sweep. ⚠️ It does **not** contain §4.4's two checks — see the row below |
| `experiments/verify_rev_final_4_4.py` | **the §4.4 artefact** (added 2026-07-30). Independent re-derivation from the raw beta matrix — pure stdlib, no shared code with the measurement script. Checks V1 contrast A (pipeline validation), V2 failing-intermediates, V3 dose-response, V4 the §4.3 below-age-20 accounting. **All reproduce** |
| `verify_rev_final_4_4_results.json` | its output |
| `experiments/DELTAAGE_LAB_NOTEBOOK.md` | dated results entries |

---

## 10. Validation — every claim in this document, and how it was checked

Written so a reviewer does not have to re-derive anything. Each row states the claim, the method of
verification, and the result. **Two errors were found and corrected this way; both are listed.**

### 10.1 The problem statement (§1)

| claim | how validated | result |
|---|---|---|
| non-responders read **+36.5 yr** after 11 days | re-derived from `GSE165176` with independently written code (not by re-running the original script) | ✅ +36.5 |
| responder age **78.7 → 101.3 → 77.2** at days 7/9/11 | same | ✅ exact |
| corr(age, pluripotency) = **−0.62**, corr(age, fibroblast) = **+0.62** | computed over all 112 OSKM-exposed fibroblasts using curated marker sets | ✅ −0.617 / +0.618 |
| the two RNA datasets give **opposite signs** | read `diag_e1_trajectory_results.json` (Gill) and `diag_d2_replication_results.json` (HFF) | ✅ Gill **+0.205**, HFF **−0.2143** |
| Gill's quote on existing transcription clocks | fetched the paper text directly | ✅ verbatim |
| GSE113957 gives **MAE 0.77 yr, ρ 0.99** through our path | ran the clock on the real NCBI counts, 143 samples | ✅ 0.769 / 0.9923 |
| Gill day-0 adults: **+18.0 yr for a 21 yr gap, ρ +0.60** | median-split of in-range donors (neonatal excluded) | ✅ +18.02 / +0.60 |

### 10.2 Method (§2)

| claim | how validated | result |
|---|---|---|
| 96 samples, **3 donors O1=53, O2=53, O3=38** | parsed the series matrix directly | ✅ exact |
| values are **beta values**, not IDATs | read the first data rows; range check | ✅ `cg…` IDs, 0.0118–0.0264 on row 1, all in [0,1] |
| file is **comma-separated with `Detection Pval` interleaved** | parsed the header | ✅ 193 header cols → 96 sample cols |
| GSE165176 has **no untreated control at day > 0** | enumerated cell types × days | ✅ only 6 day-0 `Dermal fibroblast` |
| GSE165179 has **21 (+12) untreated controls** | same | ✅ 21 fibroblast + 12 intermediate |
| skin & blood **391 CpGs, 100% coverage**; multi-tissue **353, 94.6%** | counted coefficients; counted CpGs actually matched per sample | ✅ 391/391 and 334/353 |
| CpG counts match the **published** clocks | 2013 = 353, 2018 skin & blood = 391, from the papers | ✅ both match exactly |
| coefficient tables carry **no intercept row** | searched both CSVs for an `intercept` row | ✅ none in either |
| Horvath transform is implemented correctly | unit tests against published fixed points: `F(20)=0`, `F(0)=−log 21`, and `anti_trafo(trafo(x))=x` across 0–96 | ✅ all pass to 1e-9 |
| replicates are `exp1`/`exp2` of one condition | listed the duplicate `(donor, day, arm)` groups | ✅ every duplicate is an exp1/exp2 pair |

### 10.3 Results (§3) and robustness (§4)

| claim | how validated | result |
|---|---|---|
| all nine result cells in §3 | produced by `diag_methylation_anchor.py`; per-pair values written to the results JSON for inspection | ✅ reproducible |
| pair counts **9 / 12 / 12** | derived from the metadata before ages were computed | ✅ and cross-checked against an independent earlier count of 9 |
| conclusions are **intercept-independent** | swept the intercept −0.60 → +0.70 and re-ran all three contrasts at each value | ✅ A and C stable throughout; **B flips**, which is why B is treated as open |
| day-profile peaks at **13 d** for contrast B on both clocks | per-day breakdown | ✅ −14.1 and −18.4, both maximal at 13 d |

### 10.4 The open question (§5) and contrast A's corroboration (§4.4)

| claim | how validated | result |
|---|---|---|
| per-clock spread of contrast B | sd of the 9 paired differences, each clock | ✅ **18.2** (skin & blood), **11.6** (multi-tissue) |
| MDE at n=9, **per clock** | `t(8)·sd/√9` | ✅ **14.0** and **8.9** yr |
| multi-tissue is the more precise instrument here | compared sd across **five** quantities, including the untreated-control ages where no effect exists | ✅ tighter on 4/5, incl. pure noise (**5.2** vs **6.6**); *not* uniform — negative-control pairs favour skin & blood (4.3 vs 5.1) |
| pooling across reprogramming lengths is well specified | one-way test for day heterogeneity in B | ✅ p = **0.852** / **0.255** — no significant day structure |
| ≈16 pairs needed to settle it on both clocks | solved `t(n−1)·18.2/√n ≤ 10` | ✅ n = 16 |
| **A's internal negative control** (failing intermediates) | ran the unused `Failing to transiently reprogram intermediate` arm against the same comparator | ✅ **−1.1** [−2.7, +0.5] and **−3.6** [−5.1, −2.2] vs A's −24.5 / −28.3; paired A − A-control **−23.3** / **−24.7** |
| **A's dose-response** | Spearman(reprogramming length, effect) over A's 12 pairs | ✅ **−0.885** (p = 0.0001) and **−0.842** (p = 0.0006); slope ≈ −3.2 yr/day |
| the intercept cancels algebraically | derived from `anti_trafo(x) = 21x + 20` above age 20, then recomputed every contrast as `21·(lp_t − lp_c)` | ✅ matches §3 to < 1 yr on all six cells |

**⚠️ This subsection supersedes an earlier version of itself.** It previously stated
*"MDE ≈13.7 yr at n=9, ≈14–18 pairs needed"* — computed from the **skin & blood spread alone** and
then applied to both clocks. Corrected above: the MDE differs per clock (14.0 vs 8.9), and on
multi-tissue **n=9 is already adequate**, which is the actual reason one clock reaches significance
and the other does not. Recorded rather than silently edited.

### 10.5 ⚠️ Two errors found by this validation — both corrected

**Error 1 — "HFF is 79% of age labels" was wrong; the true figure is ~99.8%.**
This figure was inherited from an earlier draft (`STAGE_1_5_1_REVISED.md`) and repeated here without
being checked. It divided **HFF-in-the-train-split (33,613)** by **total cells across all splits
(42,584)** — a numerator and denominator drawn from different populations. Validated correctly two
independent ways:

* the inner-LODO skip line states holding out HFF *"leaves 75 of 33,688 training cells"* ⇒ HFF is
  **33,613 / 33,688 = 99.8%** of training cells;
* the project's own record already says **"33,613 of 33,688 pooled residuals (99.8%)"**.

**This makes §6's point stronger, not weaker:** HFF is essentially the entire age-labelled dataset,
so the inability to anchor it matters more than the 79% figure implied. §6 now reads 99.8%.

**Error 2 — an earlier pairing bug that dropped two-thirds of the data.**
The first implementation required a *unique* sample per `(donor, day, arm)`, silently discarding
**6 of 9** contrast-B pairs and leaving only day 10. Caught by cross-checking against a metadata-only
pair count computed before the ages existed. Replicates are now averaged, with a regression test that
names the defect. Already reflected in §2.

### 10.6 Claims deliberately NOT validated, and why

| claim | status |
|---|---|
| O1/O2 are the **same physical donors** in GSE165176 and GSE165179 | ❌ **not verifiable** from the metadata — matched by label and age only. Stated as an assumption in §6.3. **Nothing in §3–§5 depends on it**, because every contrast is within GSE165179; the assumption would only matter for a cross-dataset comparison, and §6.2 establishes that none is possible |
| the **published** Horvath intercepts | ✅ **no longer needed.** §4.3 shows the intercept cancels *algebraically* in every contrast used here (all are differences, and Horvath's transform is linear above age 20). Verified by recomputing all six cells intercept-free. It would still be needed for *absolute* ages, which this document does not use |
| Gill's ~30 yr is the **median across clocks** | ⚠️ taken from the paper's own wording; not recomputed from their data |
| absolute methylation ages | ⚠️ approximate by construction (derived intercept). Only differences are used |

### 10.7 Test and hygiene status

| check | result |
|---|---|
| full suite | ✅ **455 passed**, three consecutive runs |
| `diag_methylation_anchor` tests | ✅ 28, covering every verdict branch, the transform's fixed points, and the replicate regression |
| `git diff --stat src/` | ✅ **empty** — no model, training, calibration or inference code touched |
| the five earlier 1.5.1 drafts | ✅ verified **byte-unmodified** |

> **One flake, recorded rather than hidden:** on a single suite run one test failed and the name was
> not captured; 455 passed on the three runs before and after. Most likely a Windows temp-file lock,
> which has occurred in this repository before — but it is unconfirmed, so it is noted here rather
> than described as clean.

---

# 11. 🆕 ADDED 2026-07-31 — every open item in this document, with an owner

> *Additive. **Nothing above this section is modified.** This is the closure pass: each open item is
> either answered here, answered by a stage that has since run, assigned to a named future stage, or
> declared a permanent limit. "We will answer this in Stage X" is an acceptable answer; "still open"
> with no owner is not.*

## 11.1 The ledger

| # | open item | status | where |
|---|---|---|---|
| **§5 / §6.1** | retention after return to fibroblast identity (contrast B) | ⏳ **OPEN — owned by Stage 6.** Needs ≈16 pairs; we have 9, and no held series adds any | §11.2 |
| **§6.2** | an anchor for **Gill's** training labels | ✅ **ANSWERED — NO**, by Stage 1.5.2 M-2a | §11.3 |
| **§6.2** | an anchor for **HFF's** training labels (99.8% of them) | ⏳ **OPEN — owned by Stage 6.** No methylation exists for HFF in any public series | §11.3 |
| **§6.3 / §10.6** | are O1/O2 the **same physical donors** across the series? | ✅ **ANSWERED — YES**, measured, this section | §11.4 |
| **§6.4** | the neonatal out-of-range limit | 🔒 **PERMANENT** at the clock level — but it is now **load-bearing**, and 1.5.2's G-c step 2 is where it bites | §11.5 |
| **§6.5** | absolute methylation ages are approximate | ✅ **QUANTIFIED** by Stage 1.5.2 §12-R: LODO MAE **6.03 / 6.63 yr**, ±7 yr donor-level | §11.5 |
| **§10.6** | the published Horvath intercepts | ✅ already closed by §4.3 — the intercept cancels algebraically | — |
| **§10.6** | Gill's ~30 yr is the median across clocks | 🔒 **WON'T FIX** — a claim about *their* paper, not ours; nothing here depends on it | §11.5 |
| **§10.7** | one uncaptured test failure on a single run | ✅ **CLOSED** — not reproduced in ~15 full-suite runs since | §11.5 |
| **§8.5** | "what still needs new data" | ✅ **superseded** by this table; the two live items are HFF methylation and §5's donors | §11.2 |

## 11.2 §5 — retention: OPEN, owned by **Stage 6 (acquisition)**

Nothing acquired since has moved it, and that is worth stating explicitly because two new series
*were* acquired:

| series | does it help §5? | why |
|---|---|---|
| **GSE165178** | ❌ no | it is the **Sendai** arm; contrast B is a **transient**-arm quantity |
| **GSE165177** | ❌ no | it is **RNA**, and §1 established the RNA clock cannot measure this |

§5's own arithmetic stands: **≈16 pairs** at the skin & blood spread, against the **9** that exist.
GSE165179 contains every transient-arm pair there is — 3 donors × 4 reprogramming lengths, and O1
contributes only 2 of them. **So the requirement is more donors, and only Stage 6 can supply them.**

⚠️ **One thing 1.5.2 changed about how §5 should be read.** §12-R measured donor-level methylation
clock error at **±7 yr** on these three donors (two donors of identical age 53 read 44.0 and 58.5).
§5's MDE arithmetic used the *within-contrast* paired spread, which is the right quantity for a
paired test and is unaffected — but it means **the −6 to −9 yr retention effect is the same size as
the between-donor error of the instrument measuring it.** More donors help the pairing; they will not
make the instrument sharper. Stage 6 should size for that, not just for n.

## 11.3 §6.2 — the anchors: Gill's is ANSWERED (no), HFF's is OPEN

**Gill: answered, and negatively.** §8.4 set the test — *"Does the transcriptomic ΔAge agree with the
methylation ΔAge on the same 22 samples? Agreement ⇒ the RNA clock is calibratable."* Stage 1.5.2 ran
it and both halves of §8.3's optimism failed:

| §8.3's claim | what Stage 1.5.2 measured |
|---|---|
| "the calibration route is back on the table" | **M-2a: SPLIT ⇒ NOT CALIBRATABLE.** ρ_partial +0.267 / +0.516 vs a bar of 0.50 |
| "M-2 becomes well-defined and **adequately powered**" | the decisive bar was **UNRESOLVABLE at n=22** (92.3%). GSE165178's geometry alone could not have decided it; **GSE165177's n=68 is what made a verdict possible** |

So §8.3's correction was right that the *data* existed and wrong that it was *sufficient*. Both are
now on the record.

**HFF: open, owned by Stage 6.** Unchanged and unchangeable by analysis — no methylation exists for
HFF in any public series. §6.2's wording stands verbatim: *"Extending the age target to HFF requires
methylation for HFF — a Stage 6 acquisition, not an analysis choice."*

## 11.4 §6.3 — ANSWERED: the donors are the same people

**§10.6 says this is "❌ not verifiable from the metadata". That is true of the metadata and false of
the data.** Methylation carries a genotype fingerprint, and both GSE165178 and GSE165179 are arrays.

`python experiments/diag_donor_identity.py --run "D:\GSE165178" "D:\GSE165179"` →
`diag_donor_identity_results.json`. Bars committed before any beta value was read (`be51c80`).
`src/` untouched.

**The design's controls are the roster asymmetry itself:** GSE165178 has O1/O2/**Y1/Y2**; GSE165179
has O1/O2/**O3**. Y1 and Y2 *cannot* match correctly, so they measure what a spurious match looks
like — without them, a high correlation everywhere (same array chemistry, same cell type) would be
indistinguishable from identity.

| query (Sendai) | O1 | O2 | O3 | best | margin |
|---|---|---|---|---|---|
| **O1** | **0.9619** | 0.8416 | 0.4272 | **O1** ✅ | **0.1203** |
| **O2** | 0.7719 | **0.9755** | 0.3925 | **O2** ✅ | **0.2036** |
| Y1 *(no counterpart)* | 0.7382 | 0.6754 | 0.5897 | — | 0.0628 |
| Y2 *(no counterpart)* | 0.7033 | 0.6529 | 0.5817 | — | 0.0504 |

**Both pre-registered conditions met:** every shared label matches itself, **and** their margins
(min **0.1203**) exceed every no-counterpart donor's (max **0.0628**). The second condition is not
decoration — a panel with no identity signal gets both right **10.9% of the time**, so a correct
assignment alone would not have been evidence.

**⇒ `SAME_DONORS`. §6.3's assumption is now measured, not assumed, and §10.6's row moves from
"not verifiable" to "verified".**

### What the panel's own failure taught, recorded because it is the more interesting half

Selection went through **two attempts and a bar audit**, and the run **aborted twice before the
assignment was ever computed** — which is what makes the refinements legitimate rather than fishing:

| | panel | cross-arm stability (bar ≥ 0.95) | outcome |
|---|---|---|---|
| attempt 1 | top 5000 by between-donor F | 0.821 / 0.942 / 0.966 | ❌ aborted |
| attempt 2 | **419 trimodal** (genotype-shaped) probes | **0.938** / 0.985 / 0.990 | ❌ aborted |

**Then the bar itself was audited — because I had set it by assertion, which is exactly the §5b
violation this project has caught four times.** Simulated at the real geometry with array noise taken
from GSE165179's own exp1/exp2 replicates (sd 0.1019, 22 pairs): a **perfect** panel scores median
**0.9681** and clears 0.95 **100%** of the time. **The bar was fair, so it was not moved.** O1's
0.938 is a real shortfall.

**The diagnosis is a finding in its own right:**

| stability of the same panel, same donors | O1 | O2 | O3 |
|---|---|---|---|
| untreated vs **successfully reprogrammed** | **0.9379** | 0.9852 | 0.9901 |
| untreated vs **failed to reprogram** | **0.9903** | 0.9940 | 0.9950 |

**The panel is rock-stable against OSKM exposure and moves only in cells that *succeeded*.** And it
moves most in O1 — whose two reprogrammed samples are day 10 and day **17**, the deepest
(§3: −35.7 yr at 17 d). So this is not an n=2 artefact; it is **global demethylation during
successful reprogramming reaching even genotype-shaped CpGs**, in proportion to depth.

That is independent corroboration of §4.2 from an unexpected direction: the rejuvenation is deep
enough to perturb probes chosen *because* they should be genotype-driven and cell-state-invariant.

**Consequence, and the honest limit:** the assignment above is computed on **non-reprogramming cells
only** (GSE165178's CD13 arm vs GSE165179's untreated + failed arms), where the panel is
demonstrably stable. The scope was set by the stability evidence, before the assignment existed.
**A definitive test would use the 59 `rs` SNP probes** — which GEO's *processed* matrices have
stripped, but which are present in `GSE165178_RAW.tar` (466.7 MB) and GSE165179's raw IDATs. That is
the concrete route if anyone wants it stronger; it is a download, not an experiment.

### Does anything depend on this?

**No — and that is now checked rather than asserted.** §10.6 already said *"nothing in §3–§5 depends
on it"*. The same holds for Stage 1.5.2: **every contrast in both stages is within a single
experiment.** M-2a is GSE165177 × GSE165179 (both transient); M-2b is GSE165178 × GSE165176 (both
Sendai); §12-R's R1a/R1b/R1d are entirely within GSE165179. **No result in either stage crosses the
Sendai/transient boundary on a donor label.** So this section strengthens the record; it rescues
nothing, because nothing was at risk.

## 11.5 The remaining rows, briefly

**§6.4 — the neonatal limit is PERMANENT, and it is now load-bearing.** GSE113957 has no samples
below age 1, so N2/N3/HFF absolute ages are unusable, and no analysis fixes that. What changed is its
weight: HFF is a **neonatal** line and **99.7% of the age labels**, and Stage 1.5.2's gate G-c is
precisely the decision of whether to keep training on labels produced by a clock that is out of range
for them. **Owner: Stage 1.5.2 G-c step 2** (the `age_mask` retrain), not a future acquisition.

**§6.5 — absolute methylation ages: now quantified.** The caveat was qualitative ("approximate by
construction"). Stage 1.5.2 §12-R put a number on it: **LODO MAE 6.03 / 6.63 yr**, donor-level error
**±7 yr**. So "only differences are relied on" is not just prudence — the absolute values are
demonstrably unusable at n=3 donors, and §11.2 explains why that matters for §5.

**§10.6 — "Gill's ~30 yr is the median across clocks": WON'T FIX.** It is a claim about the wording
of *their* paper, taken at face value and labelled as such. Nothing in this document or in 1.5.2
depends on it, and recomputing it would mean re-deriving another group's published summary from data
we do not have. Left as a disclosed second-hand figure.

**§10.7 — the uncaptured flake: CLOSED.** One test failed on a single run and the name was not
captured. The suite has since been run ~15 times across Stage 1.5.2 at 455 → 537 tests with **zero**
failures, including six full runs on 2026-07-31. It is not reproducing, and it is closed as a
transient rather than left as an open worry. *(The record of it stays here, per the standing rule.)*

> ### 🔴 REOPENED 2026-08-01 — **it reproduced, and the name is now captured.**
>
> *Additive; the closure above is left as written, because closing it was the mistake and deleting
> the mistake would hide it.* **Closing a flake on "it stopped happening" was premature** — that is
> absence of evidence, and this document says elsewhere not to act on it.
>
> **The test is `tests/test_evaluation.py::test_cell_line_regime_is_multiclass_with_finite_metrics`.**
> It failed once in five full-suite runs on 2026-08-01, then passed 3/3 on immediate re-runs.
> `FileNotFoundError: reports/cell_line.json`.
>
> **Two distinct problems, and only one of them is a flake:**
>
> | | |
> |---|---|
> | **A — deterministic, order-dependent.** `test_evaluate_writes_reports_and_wellformed_gates` (`test_evaluation.py:350`) *writes* `reports/cell_line.json`; this test *reads* it at `:383`. **Run alone it fails 100% of the time** (3/3 verified). Anyone debugging that single test in isolation will be misled | a real test-design defect, cheap to fix: the reader should depend on a fixture that generates the report, not on another test having run first |
> | **B — intermittent, in the full suite.** ~1 run in 5. `eval_bundle` is a **module-scoped** fixture over `tmp_path_factory`; the leading hypothesis is Windows temp-directory handling, which this repository has hit before | **not established.** Recorded as an open question, not as a diagnosis |
>
> **Scope, checked before it is dismissed:** the failure is in the evaluation **reporting** path. It
> touches no ΔAge value, no label, no `age_mask`, and nothing Stage 1.5.3 changes. It predates all
> of this work (§10.7's sighting is 2026-07-26) and **CI is green on Linux**, which is consistent
> with B being platform-specific.
>
> **Owner: whoever next touches `evaluation/`.** Fixing A is a few lines and would also remove the
> most likely amplifier of B. It does **not** block Stage 1.5.3 steps 1–4, and saying so is a
> judgement about scope, not a dismissal.
>
> ### ✅ FIXED 2026-08-01, same day — problem A is gone and B has not recurred
>
> *The judgement above said this did not block Stage 1.5.3. That was true and it was also the wrong
> place to leave it: "does not block" is not "no issue", and the fix was a few lines.*
>
> **Root cause, confirmed:** `evaluate()` was called inside
> `test_evaluate_writes_reports_and_wellformed_gates`, and **three** tests read the
> `reports/cell_line.json` it produced. Two of them therefore only worked if pytest happened to run
> the writer first.
>
> **The fix:** report generation extracted into its own module-scoped fixture, `eval_reports`, which
> returns `(reports_dir, gates)`. Every test that needs the reports now depends on the fixture
> instead of on another test, so the reports are built on demand by whichever test asks first — and
> still built only once per module. **No assertion was changed.**
>
> | check | before | after |
> |---|---|---|
> | each of the 3 tests run **alone** | ❌ 2 of 3 failed, 100% of the time | ✅ **all 3 pass** (35.4 s / 31.9 s / 30.2 s) |
> | `tests/test_evaluation.py` alone, ×3 | passed (masking the defect) | ✅ 31 passed, ×3 |
> | full suite, ×4 consecutive | 1 failure in ~5 runs | ✅ **645 passed, 1 skipped, ×4** |
>
> **On problem B:** the intermittent failure has not recurred in the four consecutive full runs
> since. That is **not proof it is gone** — B was always rare, and four clean runs is exactly the
> evidence that was too weak the first time this was closed. What *is* established is that its most
> likely amplifier is removed. **If it ever reappears, it is now a genuine fixture/tmpdir question
> and not an artefact of test ordering**, which is a materially better place to debug from.

## 11.6 What is genuinely still open, in one place

| item | owner | what it needs |
|---|---|---|
| **§5 retention** (−6 to −9 yr, at the resolution boundary) | **Stage 6** | ≈16 transient-arm pairs ⇒ **more donors**. Size for the ±7 yr between-donor instrument error, not just for n |
| **HFF's age labels** | **Stage 6** | methylation for HFF. No public series has it |
| **`age_mask` for HFF** | **Stage 1.5.2 G-c step 2** | one retrain, metric pre-registered first. Does not need new data |
| ~~the `test_evaluation` flake~~ | ✅ **CLOSED 2026-08-01** — order dependence fixed by the `eval_reports` fixture; all 3 tests now pass in isolation and the suite is clean ×4. Problem B has not recurred, and is no longer masked by ordering if it ever does | — |

**Everything else in this document is answered.** Two of the three are acquisition items that no
amount of further analysis can close, and saying so is the point of this section.
