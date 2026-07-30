# STAGE 1.5.1 — ΔAge anchored to methylation

**Status:** ✅ **EXECUTED** (2026-07-26). This document is the clean, self-contained statement of the
method, the results, and the one question that remains open. It is written to be read cold.

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

### 4.3 The conclusions do not depend on the derived intercept

Swept from −0.60 to +0.70 — wider than any plausible published value:

| | intermediates (A) | negative control (C) | returned fibroblasts (B) |
|---|---|---|---|
| skin & blood | REJUVENATION throughout | inert throughout | NO_EFFECT throughout |
| multi-tissue | REJUVENATION throughout | inert throughout | **flips** NO_EFFECT ↔ fragile |

**A and C are intercept-robust. B is not** — which is why §5 treats B as open.

The reason is structural: most predictions sit above age 20, where Horvath's transform is linear, so
a constant cancels in a difference.

### 4.4 What this means for the project

- **ΔAge has a valid anchor.** Methylation measures real, large rejuvenation on these samples with an
  inert negative control.
- **The transcriptomic clock's failure is fully localised to the instrument.** The biology is there;
  the RNA clock cannot see it. ΔAge as a *concept* is vindicated.
- **The target definition was never the problem.** No further control redefinition is warranted.

---

## 5. The one open question, with its bar set now

**How much rejuvenation survives the return to fibroblast identity?** (Contrast B.)

Current evidence: right sign on both clocks, Gill's exact day-shape, largest at 13 days
(−14.1 / −18.4) — but not significant on skin & blood, fragile on multi-tissue, and
intercept-sensitive. **This must not be claimed either way.**

**Why it cannot be settled at n=9.** The observed pair-to-pair spread is sd ≈ 17.8 yr, so at 9 pairs
the minimum detectable effect is **≈ 13.7 yr** — almost exactly the size of the effect being looked
for. The measurement is at its own resolution limit.

**Pre-registered bar for the next attempt:**

| verdict | condition |
|---|---|
| **RETAINED** | 95% CI excludes 0 and is negative, **on both clocks**, and stable across the intercept sweep |
| **NOT RETAINED** | CI includes 0 on both clocks at adequate power |
| **FRAGILE** | any verdict within 0.5 yr of a bound, or one that flips under the sweep — reported as such, never as a result |

**Resolvability (ground rule §5b): to reach a 10 yr MDE at sd 17.8 requires ≈ 14–18 pairs — roughly
double what exists.** So the honest prerequisite is **more donors**, not another analysis of these
nine. Stating that now prevents a re-analysis being mistaken for new evidence.

---

## 6. What this stage does *not* establish

1. **Retention after return to fibroblast identity** — §5, open.
2. **Any anchor for HFF.** HFF (79% of the project's age labels) has no methylation data. It cannot
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

**The honest caveat on contrast A:** because it was promoted from secondary to headline, a reviewer
may reasonably require it be re-registered as primary and re-run on new donors before it supports a
manuscript claim. The arguments that it is nonetheless sound: it was named pre-run, two independent
clocks agree, the negative control is inert, it is intercept-robust, and the day-profile separates it
from contrast B on a priori grounds.

---

## 8. Next step

**Use methylation as the ΔAge source where methylation exists** (Gill's three donors). **Treat HFF's
age labels as unanchored** — the fate head is unaffected and remains strong.

Two things unlock the rest, both data acquisitions rather than analyses:

1. **More donors with paired methylation** — the only way to settle §5's retention question.
2. **Methylation for HFF** — the only way to anchor the 79% of age labels it holds.

Neither is a code change, and no further RNA-side re-analysis will move either question.

---

## 9. Artefacts

| file | role |
|---|---|
| `experiments/diag_methylation_anchor.py` | the measurement; read-only; pure verdict logic separated from I/O |
| `tests/test_diag_methylation_anchor.py` | 28 tests — every verdict branch, Horvath's transform against its published fixed points, the replicate-averaging regression |
| `configs/clocks/horvath_skin_blood_2018.json`, `horvath_multitissue_2013.json` | clock coefficients, with provenance in `meta` |
| `diag_methylation_anchor_results.json` | full output including per-pair values and the intercept sweep |
| `experiments/DELTAAGE_LAB_NOTEBOOK.md` | dated results entries |
