# Message to send — PLAIN TEXT, no tables, no markdown

The recipient reads this in email or chat, not a markdown renderer, so it carries no
tables and nothing that degrades when pasted. Copy everything below the line.

Internal record with full provenance: `plans/EXPERT_BRIEF.md`.

---

Hi [NAME],

I'm working on a model that reads the transcriptome of a fibroblast culture partway
through OSKM reprogramming and predicts, for each candidate withdrawal day, (a) the
transcriptomic age change and (b) whether the culture has passed into an unsafe state
(identity loss / apoptosis). I've hit two things I can't resolve from the data or the
literature on my own, and I'd rather ask once and ask well. Everything below is
measured - happy to send numbers, code or plots for any of it.

WHAT I HAVE

- GSE242423: 10x scRNA-seq, HFF neonatal foreskin fibroblasts, OSKM, 9 timepoints
  D0 to D14 plus iPSC, ~42k cells, ~4,700 per timepoint.
- GSE165176 (Gill 2022): bulk RNA-seq, 6 donors, 11-12 timepoints, but only ~1.7
  samples per donor-timepoint and a single day-0 control per donor.
- GSE165177 (Gill 2022): bulk RNA-seq, 3 donors aged 53/53/38, days 0/10/13/15/17,
  4-6 treated samples per donor-day, and 33 CONTEMPORANEOUS negative controls at 2-3
  per donor per timepoint.
- The age clock is a linear expression clock trained on Fleischer 2018 (GSE113957,
  133 fibroblast samples, ages 1-96, cross-validated MAE 12.3 yr).

FINDING 1 - THE SAFETY HALF IS BLOCKED, AND I THINK IT'S STRUCTURAL

I define "unsafe" as the FRACTION OF CELLS past an identity-loss / apoptosis
threshold. On bulk data that fraction can't be measured: each sample is already a
population average, so a per-sample call collapses it to 0 or 1 before it can be
counted, and "fraction of cells" silently becomes "fraction of samples".

The tell: on GSE165177 the target came back saturated at 1.0 in 11 of 12 donor-day
groups, and the untreated day-0 fibroblasts were also called "identity lost", which
can't be right. More replication doesn't help - GSE165177 has 4-6x the replication of
GSE165176 and is MORE saturated, not less. So I think this half needs single-cell
data, and I only have one single-cell line.

FINDING 2 - THE dAGE HALF WORKS, BUT THE CLOCK'S ABSOLUTE NUMBERS ARE BADLY OFF

Measuring dAge against contemporaneous, replicated same-day controls, I reproduce the
published rejuvenation effect: -17.9 yr (95% CI -21.1 to -14.6) for transiently
reprogrammed cells vs their own control, and -9.6 yr (CI -12.8 to -6.4) against
failed-to-reprogram cells, with 11 of 12 donor-day groups negative.

But absolute ages are wrong, and they split in an interesting way. The day-0
fibroblasts, which never saw reprogramming media (n=3), read 78.0 yr against a true
mean of 48.0 - a bias of +30 yr. The negative controls, cultured alongside the
experiment for 10-17 days (n=33), read 95.6 against the same true 48.0 - a bias of
+48 yr.

So there's a ~+30 yr floor even on fresh cells, plus a further ~+18 yr that tracks
time in culture. It isn't missing genes: 57% of the clock's genes are present in this
data, but those carry 89% of the clock's total absolute weight. dAge survives all of
it because the bias cancels between a sample and its same-day control.

QUESTIONS - answer any subset; the three marked (*) are what I'd most value

On the safety half:

1. Do you know of a human OSKM reprogramming SINGLE-CELL RNA-seq time course covering
   MORE THAN ONE donor or line? Multi-line within one study is worth far more to me
   than several single-line studies, since comparing across studies re-confounds line
   with batch.

2. Failing that, is there any second human somatic-to-iPSC scRNA-seq time course with
   at least 6 sampled days and a day-0 control? I searched GEO by title and couldn't
   verify one, but I'd miss anything filed as "iPSC generation", bundled inside a
   SuperSeries, or deposited outside GEO.

3. (*) Is it even biologically plausible that the TIMING of identity loss during OSKM
   reprogramming transfers across donor lines? If it's strongly line-, passage- and
   batch-specific, that's decisive and I'd stop. A well-founded "no" here is more
   useful to me than a dataset.

4. If nothing exists: roughly what would it cost to generate a minimum version - 1-2
   lines other than HFF, 6-8 timepoints from day 0 to ~day 14 plus a terminal sample,
   ~3-5k cells per timepoint, standard 10x 3-prime (roughly 8-16 lanes)?

On the dAge half:

5. (*) Is a +30 yr offset simply expected when an expression-based age clock trained
   on one lab's fibroblast RNA-seq is applied to another's? Mine reads fresh 53-, 53-
   and 38-year-old fibroblasts as 72-82. Is cross-study ABSOLUTE transfer just not a
   thing for expression clocks, unlike methylation clocks, making control-relative use
   the only defensible mode?

6. (*) Do fibroblasts cultured 10-17 days in reprogramming media genuinely age ~18
   transcriptomic years, or is that a confluence / media / passage artefact? It
   cancels in my design so it doesn't threaten the result, but I'd like to know which
   it is.

7. My clock reads the 38-year-old donor as OLDER than one of the 53-year-olds. With a
   12.3 yr CV MAE against a 15 yr age gap I treat this as unresolvable and claim
   nothing. Is 3 donors across 2 distinct ages simply hopeless for validating age
   discrimination, and what donor-age spread would actually be needed?

8. Is a ~-18 yr dAge for transient reprogramming at days 10-17 consistent with the
   literature, or suspiciously large? I'd rather hear "too big, check X" now than
   later.

9. Given that dAge does work on bulk with contemporaneous replicated controls, is a
   second BULK dataset of that design a sensible buy, and what would you prioritise:
   more donors, wider age range, or more timepoints?

IF YOU DO KNOW OF A DATASET - WHAT IT NEEDS TO HAVE

Must have, all of them:
- Human, not mouse. The clock is human; an ortholog mapping reintroduces the confound.
- SINGLE-CELL RNA-seq (droplet / 10x or equivalent). Bulk cannot express the target,
  per Finding 1.
- Whole-transcriptome, not a targeted panel.
- A DIFFERENT donor or line from HFF (neonatal foreskin fibroblast). A second HFF
  experiment tests nothing.
- At least 6 sampled timepoints spanning day 0 to a terminal point. Three is the bare
  minimum for any forward pairs at all.
- A day-0 untransduced control of the same line. Non-negotiable - it is both the age
  zero-point and the reference the safety labels are scored against.
- At least ~300 cells per timepoint after QC.
- Raw or filtered count matrices with the day recoverable, not only a processed
  embedding or Seurat object.

Strongly prefer:
- Three or more independent lines in a single study. This is the highest-value item -
  it's the difference between one anecdote and a real confidence interval.
- Adult donors of known chronological age.
- A similar OSKM delivery method to GSE242423's.
- A comparable timepoint grid.

Deal-breakers:
- Bulk RNA-seq at any depth.
- Mouse or MEF, including the Waddington-OT time course.
- Differentiation FROM iPSC. This is by far the most common false hit when searching
  GEO for "iPSC".
- Direct lineage conversion, e.g. fibroblast to neuron.
- A single endpoint with no intermediate timepoints.
- Microarray. GSE241435 has exactly the right design but the wrong technology.

I do NOT need any fate or safety annotation, age annotation, cell-type labels,
clustering, or matched ATAC/methylation - I compute all of that myself. Raw counts
plus a day label are enough.

Thanks very much for reading this far. Even a partial answer, especially on 3, 5 or 6,
would save me a lot of wasted effort.

[YOUR NAME]
