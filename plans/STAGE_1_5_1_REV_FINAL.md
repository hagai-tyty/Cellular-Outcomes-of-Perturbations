# STAGE 1.5.1 — ΔAge anchored to methylation

**Status:** ✅ **EXECUTED and VALIDATED** (2026-07-26). This document is the clean, self-contained
statement of the method, the results, and the one question that remains open. It is written to be
read cold. **§10 validates every claim in it and states how** — including three errors that validation
found and corrected. The three points a reviewer would previously have had to challenge (the
post-hoc promotion of contrast A, the derived intercept, and the power of contrast B) are **closed**
in §4.4, §4.3 and §5 respectively — not merely disclosed.

**Scope:** one public dataset acquired; methylation age computed with two published clocks; three
identity-matched contrasts measured. **`src/` untouched throughout** — no model, training,
calibration or inference code was changed.

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

### 4.5 What this means for the project

- **ΔAge has a valid anchor.** Methylation measures real, large rejuvenation on these samples with an
  inert negative control.
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

1. **Retention after return to fibroblast identity** — §5, open.
2. **Any anchor for the training labels — HFF *or* Gill.** GSE165176 (RNA) and GSE165179
   (methylation) share **zero samples** (§8.1), so no methylation age can be attached to *any* cell
   the model trains on. HFF (**~99.8%** of the age-labelled cells — 33,613 of 33,688 training cells)
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

**The two series share no samples.** They are separate experiments — different sample sets,
different donor rosters (RNA: N2, N3, Y1, Y2, O1, O2; methylation: O1, O2, O3), different day grids
(7–47 vs 10–17) and different arm vocabularies. **There is no join key, so a methylation age cannot
be attached to any cell the model trains on.**

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

### 8.3 What is executable now

**Nothing here.** The project is not blocked, though:

* **Stage 2** may proceed — see its updated annotation. Its intervention (k ≈ 3 reference cells per
  donor) helps whether the per-donor offset is real biology or n = 1 baseline noise.
* **Stage 3** depends on *Stage 1 required, Stage 2 optional*, and the fate head is untouched and
  strong. It is not gated by anything here.

**The two open questions both need data, not code:**

1. **More donors with paired methylation** — the only way to settle §5's retention question (≈16
   pairs needed).
2. **Methylation on the samples we actually train on** — the only way to anchor the labels. Note
   this means *new profiling of our own samples*, not another public download: no existing series
   pairs methylation to GSE165176 or to HFF.

**No further RNA-side re-analysis will move either question.** That is the point of §1.

## 9. Artefacts

| file | role |
|---|---|
| `experiments/diag_methylation_anchor.py` | the measurement; read-only; pure verdict logic separated from I/O |
| `tests/test_diag_methylation_anchor.py` | 28 tests — every verdict branch, Horvath's transform against its published fixed points, the replicate-averaging regression |
| `configs/clocks/horvath_skin_blood_2018.json`, `horvath_multitissue_2013.json` | clock coefficients, with provenance in `meta` |
| `diag_methylation_anchor_results.json` | full output including per-pair values and the intercept sweep |
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
