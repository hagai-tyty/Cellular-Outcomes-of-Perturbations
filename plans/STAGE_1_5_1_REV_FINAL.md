# STAGE 1.5.1 — ΔAge anchored to methylation

**Status:** ✅ **EXECUTED and VALIDATED** (2026-07-26). This document is the clean, self-contained
statement of the method, the results, and the one question that remains open. It is written to be
read cold. **§10 validates every claim in it and states how** — including two errors that validation
found and corrected.

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
2. **Any anchor for HFF.** HFF (**~99.8%** of the age-labelled cells — 33,613 of 33,688 training cells) has no methylation data. It cannot
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
2. **Methylation for HFF** — the only way to anchor the ~99.8% of age labels it holds.

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

### 10.4 The open question (§5)

| claim | how validated | result |
|---|---|---|
| pair-to-pair spread **sd ≈ 17.8 yr** | recomputed from contrast B's own CI: half-width 13.7 at t(8)=2.306 ⇒ SE 5.94 ⇒ sd = SE·√9 | ✅ 17.8 |
| **MDE ≈ 13.7 yr at n=9** | `t(8)·sd/√9` = 2.306 · 17.8 / 3 | ✅ 13.7 |
| **≈14–18 pairs** needed for a 10 yr MDE | solve `t(n−1)·17.8/√n = 10` | ✅ n ≈ 17 (16 gives 10.3; 18 gives 9.6) |

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
| O1/O2 are the **same physical donors** in GSE165176 and GSE165179 | ❌ **not verifiable** from the metadata — matched by label and age only. Stated as an assumption in §6.3. Nothing in §3–§5 depends on it, because no cross-dataset comparison is made |
| the **published** Horvath intercepts | ❌ not obtained. Mitigated by the sweep (§4.3), which shows the two robust conclusions hold for any value in a range wider than any published one |
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
