# Message to send — copy from the line below

---

Hi [NAME],

I'm working on a model that reads the transcriptome of a fibroblast culture partway through OSKM
reprogramming and predicts, for each candidate withdrawal day, (a) the transcriptomic age change
and (b) whether the culture has passed into an unsafe state (identity loss / apoptosis). I've hit
two things I can't resolve from the data or the literature on my own, and I'd rather ask once and
ask well. Everything below is measured — happy to send the numbers, code or plots for any of it.

**What I have:** GSE242423 (10x scRNA-seq, HFF neonatal foreskin fibroblasts, OSKM, 9 timepoints
D0→D14→iPSC, ~42k cells, ~4,700 per timepoint), plus the two Gill 2022 bulk RNA-seq sets
(GSE165176: 6 donors, 11–12 timepoints, ~1.7 samples per donor-timepoint, one day-0 control per
donor; GSE165177: 3 donors aged 53/53/38, days 0/10/13/15/17, 4–6 treated samples per donor-day and
33 *contemporaneous* negative controls, 2–3 per donor per timepoint). The age clock is a linear
expression clock trained on Fleischer 2018 (GSE113957, 133 fibroblast samples, ages 1–96, CV MAE
12.3 yr).

---

## Finding 1 — the safety half is blocked, and I think it's structural

I define "unsafe" as the **fraction of cells** past an identity-loss / apoptosis threshold. On
bulk data that fraction can't be measured: each sample is already a population average, so a
per-sample call collapses it to 0 or 1 before it can be counted, and "fraction of cells" silently
becomes "fraction of samples".

The tell: on GSE165177 the target came back saturated at 1.0 in 11 of 12 donor-day groups, and the
**untreated day-0 fibroblasts were also called "identity lost"** — which can't be right. More
replication doesn't help; GSE165177 has 4–6× the replication of GSE165176 and is *more* saturated,
not less. So I think this half needs single-cell data, and I only have one single-cell line.

## Finding 2 — the ΔAge half works, but the clock's absolute numbers are badly off

Measuring ΔAge against contemporaneous, replicated same-day controls, I reproduce the published
rejuvenation effect: **−17.9 yr (95% CI −21.1 to −14.6)** for transiently reprogrammed vs its own
control, and **−9.6 yr (CI −12.8 to −6.4)** against failed-to-reprogram cells, with 11 of 12
donor-day groups negative.

But absolute ages are wrong, and they split in an interesting way:

| samples | n | clock reads | true age | bias |
|---|---|---|---|---|
| day-0 fibroblasts (never in reprogramming media) | 3 | 78.0 | 48.0 | **+30 yr** |
| negative controls (cultured alongside, 10–17 d) | 33 | 95.6 | 48.0 | **+48 yr** |

So a ~**+30 yr floor even on fresh cells**, plus a further ~**+18 yr that tracks time in culture**.
It's not missing genes — 57% of clock genes are present but that's 89% of the clock's total
absolute weight. ΔAge survives it because the bias cancels between a sample and its same-day
control.

---

## Questions

Answer any subset — the starred three are what I'd most like your view on.

**On the safety half:**

1. Do you know of a human OSKM reprogramming **single-cell** RNA-seq time course covering **more
   than one donor or line**? Multi-line within one study is worth far more to me than several
   single-line studies, since comparing across studies re-confounds line with batch.
2. Failing that, is there *any* second human somatic→iPSC scRNA-seq time course with ≥6 sampled
   days and a day-0 control? I searched GEO by title and couldn't verify one, but I'd miss things
   filed as "iPSC generation", bundled in a SuperSeries, or deposited outside GEO.
3. ⭐ **Is it even biologically plausible that the *timing* of identity loss during OSKM
   reprogramming transfers across donor lines?** If it's strongly line-, passage- and
   batch-specific, that's decisive and I'd stop — a well-founded "no" here is more useful to me
   than a dataset.
4. If nothing exists: roughly what would it cost to generate a minimum version — 1–2 lines other
   than HFF, 6–8 timepoints day 0 to ~day 14 plus a terminal sample, ~3–5k cells/timepoint,
   standard 10x 3′ (~8–16 lanes)?

**On the ΔAge half:**

5. ⭐ **Is a +30 yr offset simply expected when an expression-based age clock trained on one lab's
   fibroblast RNA-seq is applied to another's?** Mine reads fresh 53-, 53- and 38-year-old
   fibroblasts as 72–82. Is cross-study *absolute* transfer just not a thing for expression clocks
   (unlike methylation clocks), making control-relative use the only defensible mode?
6. ⭐ **Do fibroblasts cultured 10–17 days in reprogramming media genuinely age ~18 transcriptomic
   years, or is that a confluence / media / passage artefact?** It cancels in my design, so it
   doesn't threaten the result — but I'd like to know which it is.
7. My clock reads the 38-year-old donor as *older* than one of the 53-year-olds. With a 12.3 yr CV
   MAE against a 15 yr age gap I treat this as unresolvable and claim nothing. Is 3 donors across
   2 distinct ages simply hopeless for validating age discrimination, and what donor-age spread
   would actually be needed?
8. Is a ~−18 yr ΔAge for transient reprogramming at days 10–17 consistent with the literature, or
   suspiciously large? I'd rather hear "too big, check X" now than later.
9. Given ΔAge does work on bulk with contemporaneous replicated controls — is a second *bulk*
   dataset of that design a sensible buy, and what would you prioritise: more donors, wider age
   range, or more timepoints?

---

## If you do know of a dataset — what it needs to have

**Must have (all of them):**

- Human, not mouse — the clock is human and an ortholog mapping would reintroduce the confound
- **Single-cell** RNA-seq (droplet / 10x or equivalent) — bulk cannot express the target, per
  Finding 1
- Whole-transcriptome, not a targeted panel
- **A different donor / line from HFF** (neonatal foreskin fibroblast) — a second HFF experiment
  tests nothing
- **≥6 sampled timepoints** spanning day 0 → terminal (3 is the bare minimum for any forward pairs)
- **A day-0 untransduced control of the same line** — non-negotiable; it's both the age zero-point
  and the reference the safety labels are scored against
- ≥~300 cells per timepoint after QC
- Raw or filtered count matrices with the day recoverable (not only a processed embedding or
  Seurat object)

**Strongly prefer:** ≥3 independent lines in a single study (this is the single highest-value
item — it's the difference between one anecdote and a real confidence interval); adult donors of
known chronological age; a similar OSKM delivery method to GSE242423's; a comparable timepoint grid.

**Deal-breakers:** bulk RNA-seq at any depth; mouse/MEF (including the Waddington-OT time course);
differentiation *from* iPSC (by far the most common false hit when searching GEO for "iPSC");
direct lineage conversion (fibroblast→neuron); a single endpoint with no intermediate timepoints;
microarray (GSE241435 has exactly the right design but the wrong technology).

**I do NOT need** any fate/safety annotation, age annotation, cell-type labels, clustering, or
matched ATAC/methylation — all of that I compute myself. Raw counts plus a day label are enough.

Thanks very much for reading this far — even a partial answer, especially on 3, 5 or 6, would save
me a lot of wasted effort.

[YOUR NAME]
