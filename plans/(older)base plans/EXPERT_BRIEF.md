# BRIEF FOR AN EXTERNAL EXPERT — one conversation, two separate asks

> ## ⛔ **RETIRED — 2026-08-14. Not sent, and every ΔAge figure below is superseded TWICE.**
>
> The expert consultation was dropped; `plans/EXPERT_MESSAGE.md` (the sendable version) is deleted.
> **This file is kept only for its provenance index**, which records which artefact produced which
> claim, and which is still accurate.
>
> **Do not quote any ΔAge number from this document.** Every one predates two corrections:
>
> 1. **Pseudoreplication** — the intervals were computed over 12 donor-day cells that come from
>    only 3 donors.
> 2. **A double-`log1p`** — the expression was log-compressed twice, shrinking every magnitude by
>    roughly 2.4×.
>
> | quoted below | current |
> |---|---|
> | ΔAge(transient) −17.88 [−21.13, −14.64] | **−42.45 [−67.39, −17.51]**, three donor-level estimates |
> | transient − failed −9.58 [−12.77, −6.39] | **−24.19 [−49.53, +1.14]** — no longer excludes zero; reported as "direction consistent in 3/3 donors" |
> | "our −18 sits below Gill's ~30" | **reversed** — at −42.45 we are *above* it |
>
> Current records: `plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md` §6–§9,
> `plans/STAGE_1_5_8_CLOCK_ON_GSE297234_PREREG.md` §7, and `CHANGES.md`.
>
> Also superseded: the *safety*-half framing here predates the finding that a molecular-progress
> coordinate, not calendar day, is the coordinate worth testing — see
> `plans/WORK_ORDER_2026_08_14.md` P3.

**Written 2026-08-12.** Supersedes the questions in `DATA_REQUIREMENT_SECOND_TIMECOURSE.md`, which
were written before the two runs below and would have wasted the conversation.
**Self-contained** — assumes no knowledge of this project.

---

## 0. What we are building, and where it stands

A model that takes the transcriptome of a cell culture partway through OSKM reprogramming and
predicts, for each candidate **withdrawal day**, (a) the transcriptomic age change, ΔAge, and
(b) the probability the culture has passed into an unsafe state (identity loss or apoptosis).

**Two halves, and they have just come apart.** One works and one is blocked, for a reason we can
now state mechanically rather than as a hunch:

| half | status |
|---|---|
| **ΔAge / rejuvenation** | ✅ **works.** We reproduce the published rejuvenation effect on data that was in no training config |
| **forward safety (`p_unsafe`)** | ❌ **blocked**, and not by sample size — the quantity is not expressible in the data type we have |

Everything below is measured, with the artefact named. Nothing here is an estimate.

---

## 1. Data in hand

| dataset | modality | design | depth |
|---|---|---|---|
| **GSE242423** | **10x scRNA-seq** | 1 line (HFF, neonatal foreskin), OSKM, 9 timepoints D0→D14→iPSC | 42,481 cells, ~4,700/timepoint |
| **GSE165176** (Gill 2022) | bulk RNA-seq | 6 donors, transient/Sendai, 11–12 timepoints | **~1.7 samples per (donor, timepoint)**; **one** day-0 control per donor |
| **GSE165177** (Gill 2022) | bulk RNA-seq | **3 donors aged 53, 53, 38**, 4 timepoints (d10–17) + day 0 | **4–6 treated samples** per (donor, day) and **33 contemporaneous negative controls**, 2–3 per donor **per timepoint** |
| GSE113957 (Fleischer 2018) | bulk RNA-seq | 133 fibroblast samples, ages 1–96 | the aging clock's training data; cross-sectional |

The clock is a linear model over 33,155 genes, `cv_mae = 12.27 yr`, fitted range `[1, 96]`.

---

## 2. What we measured (the two findings that shape the questions)

### 2.1 `p_unsafe` is not expressible in bulk — at any replication

`p_unsafe` is defined as **the fraction of cells** in an unsafe state. We label cells from
pluripotency / somatic-identity / apoptosis marker programs, z-scored against that line's own
controls.

On GSE165177 the target came back **1.000 at 11 of 12 (donor, day) cells** — no time variation at
all. The tell: the **untreated day-0 fibroblasts also labelled "identity lost"** (P = 0.97, 0.88,
0.73), which is impossible biologically. The reason is structural:

> **A bulk sample is already a population average.** A hard label per sample collapses the
> within-sample fraction to 0 or 1 *before it can be counted*, so the "fraction of cells" silently
> becomes a "fraction of samples".

This single mechanism explains both the GSE165177 saturation **and** why the 6-donor GSE165176 set
had 63 of 70 values pinned at 0 or 1 — which we had previously blamed on its 1.7 samples per
timepoint. GSE165177 has **4–6× more replication and is *more* saturated, not less**. More bulk
cannot fix it.

### 2.2 ΔAge works — and reproduces the published effect

Same dataset, same normalisation, same clock. ΔAge is measured against **contemporaneous,
replicated controls of the same donor at the same day** — something GSE165176 cannot provide.

| | mean | 95 % CI | n |
|---|---|---|---|
| ΔAge(transiently reprogrammed) vs its own contemporaneous control | **−17.88 yr** | [−21.13, −14.64] | 12 cells |
| paired: transient − failed-to-reprogram | **−9.58 yr** | [−12.77, −6.39] | 12 cells |

11 of 12 (donor, day) cells negative. **Gill 2022's central rejuvenation claim, recovered by our
own clock, on data that was in no training config.**

**But the clock's absolute age is badly wrong here, and it decomposes cleanly:**

| samples | n | predicted | true | bias |
|---|---|---|---|---|
| **day-0 fibroblasts** (never in reprogramming media) | 3 | 78.0 | 48.0 | **+30.0 yr** |
| **negative controls** (cultured alongside 10–17 d) | 33 | 95.6 | 48.0 | **+47.6 yr** |

So: a **~+30 yr cross-study floor even on fresh cells**, plus a **further ~+18 yr that tracks time
in culture**. Not a gene-coverage artefact — 57 % of clock genes are present but that is **89 % of
the clock's total absolute weight**.

ΔAge survives all of this because it is a **difference**: the bias appears in the sample and in its
control and cancels. Absolute age needs the clock to be *accurate*; ΔAge needs it only to be
*consistent*.

---

## 3. ❓ THE QUESTIONS — ordered by what a good answer is worth

### Track A — the blocked half (forward safety)

**A1.** Do you know of a human OSKM reprogramming **single-cell** RNA-seq time course covering
**more than one donor or cell line**? Multi-line within a single study is worth far more to us than
several single-line studies, because a cross-study comparison re-confounds line with batch.

**A2.** Besides GSE242423, is there *any* published human somatic→iPSC **scRNA-seq** time course
with **≥ 6 sampled days and a day-0 control**? We searched GEO by title and could not verify one,
but our search would miss anything filed under "iPSC generation", bundled in a SuperSeries, or
deposited outside GEO (HCA, ArrayExpress/BioStudies, GSA/CNGBdb, Synapse).

**A3.** ⭐ **Is it biologically plausible that the *timing* of identity loss during OSKM
reprogramming transfers across donor lines at all?** If the answer is "no — it is strongly line-,
passage- and batch-specific", that is decisive and **saves the entire acquisition**. We would stop
and ship the retrospective model. *A well-founded "no" here is worth more to us than a dataset.*

**A4.** If nothing suitable exists, what would it take to generate the minimum viable version —
1–2 donor lines distinct from HFF, 6–8 timepoints from day 0 to ~day 14 plus a terminal sample,
~3–5k cells/timepoint, standard 10x 3′ (roughly 8–16 lanes)? Rough cost and turnaround.

### Track B — the working half (ΔAge), and the anomaly we most want explained

**B1.** ⭐ **Is a +30 yr offset expected when a transcriptomic age clock trained on one lab's
fibroblast RNA-seq is applied to another's?** Ours reads fresh 53-, 53- and 38-year-old fibroblasts
as 72–82. Is cross-study absolute transfer simply not a thing for expression clocks — as opposed to
methylation clocks — and is control-relative use the only defensible mode?

**B2.** ⭐ **Do fibroblasts cultured 10–17 days in reprogramming media genuinely "age"
transcriptomically by ~18 years, or is that a confluence/media/passage artefact?** Our negative
controls drift by that much relative to day 0. It cancels in our contemporaneous-control design, so
it does not threaten the result — but if it is real biology it is interesting in its own right, and
if it is an artefact we would like to know what drives it.

**B3.** Our clock reads the 38-year-old donor as **older** than one 53-year-old. With `cv_mae`
12.27 yr against a 15 yr age gap we treat this as unresolvable and claim nothing. **Is 3 donors
across 2 distinct ages simply hopeless for validating age discrimination — and what donor-age
spread would be needed?**

**B4.** We reproduce ~**−18 yr** ΔAge for transient reprogramming (~−9.6 yr against
failed-to-reprogram controls) at days 10–17. **Is that magnitude consistent with the literature, or
suspiciously large?** We would rather hear "too big, check X" now than later.

**B5.** Given that ΔAge works on bulk with contemporaneous replicated controls: **is a second bulk
dataset of that design a sensible buy for the ΔAge half**, and what would we look for — more
donors, wider age range, or more timepoints?

---

## 4. What we do NOT need (so nothing is over-specified)

- **No fate/safety annotation and no age annotation.** Both are computed by us from expression.
- **No cell-type annotation, clustering, or embedding.** Raw counts plus a day label suffice.
- **No matched ATAC, methylation or protein.**
- **No more bulk for the safety half** — §2.1 shows it cannot help at any depth.

## 5. Deal-breakers for Track A

Bulk RNA-seq of any depth; mouse/MEF (including the Waddington-OT time course — our clock is
human); differentiation *from* iPSC (the most common false hit when searching GEO for "iPSC");
direct lineage conversion (fibroblast→neuron); a single endpoint with no intermediate timepoints;
microarray (e.g. GSE241435 has exactly the right design but the wrong technology).

---

## 6. Provenance

| claim | artefact |
|---|---|
| `p_unsafe` saturation, day-0 fibroblasts labelled "loss" | `experiments/stage3a_regime_e.py`, `results/stage3a_regime_e_results.json` |
| bulk cannot carry a per-cell fraction; GSE165176's 63/70 | same, and `plans/STAGE_3A_REGIME_E_PREREG.md` §6 |
| ΔAge −17.88 [−21.13, −14.64]; −9.58 [−12.77, −6.39] | `experiments/dage_gse165177.py`, `results/dage_gse165177_results.json` |
| bias +30.0 (day 0) / +47.6 (cultured controls) | same, "BIAS DECOMPOSED" |
| 89.2 % of clock weight mass present | same |
| clock `cv_mae` 12.27, range [1, 96], 33,155 genes | `configs/clocks/fleischer_clock.json` |
| GSE165177 passes our integrity gate 0/93 rejected | `experiments/dage_gse165177.py` M-E0 |
