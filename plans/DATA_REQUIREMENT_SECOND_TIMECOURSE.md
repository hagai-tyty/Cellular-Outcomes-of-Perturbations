# DATA REQUIREMENT — a second single-cell reprogramming time course

**Status:** open data-acquisition requirement. Written 2026-08-12, after `STAGE_3_TOOL.md` "3a-bis"
measured that the current corpus cannot answer Stage 3a's question at any effect size.
**Audience:** an external domain expert. Self-contained — assumes no knowledge of this project.

---

## 1. What is being built, in three sentences

A model that takes a transcriptome of a cell culture undergoing OSKM reprogramming partway through
the protocol, and predicts, for each candidate **withdrawal day**, (a) the expected transcriptomic
age change (ΔAge, from a published aging clock) and (b) the probability the culture has passed into
an unsafe state — loss of somatic identity or apoptosis.

The product question is **forward and cross-culture**: *given a culture I have not seen before, in
its current state, when should I withdraw the factors?*

That requires a model trained on some cultures to work on a **new** one.

---

## 2. What we already have (please do not suggest these)

| Dataset | What it is | Modality | Trajectories | Timepoints | Depth |
|---|---|---|---|---|---|
| **GSE242423** (Kundaje lab) | human fibroblast OSKM reprogramming to iPSC | **10x scRNA-seq** | **1 line** (HFF, neonatal foreskin) | **9** — D0, D2, D4, D6, D8, D10, D12, D14, iPSC | **42,481 cells**, ~4,700/timepoint |
| **GSE165176** + GSE165177/8/9 (Gill 2022) | transient / Sendai reprogramming | **bulk RNA-seq** | 6 donors (N2, N3, O1, O2, Y1, Y2) | 11–12 each | 18–21 **samples** total per donor → **~1.7 per timepoint** |
| **GSE113957** (Fleischer 2018) | 133 fibroblast samples, donor ages 1–96 | bulk RNA-seq | n/a — cross-sectional | n/a | the aging clock's training data, not a time course |

Also already held and **not** useful for this: `GSE242421`, `GSE242424`, `GSE242419` — the ATAC and
joint-assay siblings of GSE242423, same study, same cells.

---

## 3. The measured blocker (why this request exists)

We ran the pre-registered resolvability check required before grading the gate: simulate a forward
signal that is **real by construction**, using HFF's own measured unsafe-fraction curve
(0.0835 at D0 rising to 0.9996 at D21), observe it at the **real cell counts** via
`Binomial(n_j, p)/n_j`, and ask how often a correct system is detected. Pass rate over 2000 trials,
`MIN_PASS_RATE = 0.95`:

| held-out trajectory | scored at **2** cells/timepoint | scored at **472** cells/timepoint |
|---|---|---|
| **HFF pseudo-replicate** (same line, same modality) | **0.965** ✅ | **1.000** ✅ |
| **Gill donor** (different line, bulk) | **0.000** ❌ | **0.000** ❌ |

Reading across the bottom row: giving the bulk donors 472 cells per timepoint — a counterfactual
they can never have — **changes nothing**. So the blocker is **not** sampling depth. It is that the
held-out trajectory is a different cell line measured on a different platform.

**And we cannot say which of those two, because in this corpus they are perfectly confounded:**
every Gill donor is bulk, and the only single-cell dataset is a single line. That is the gap.

The dense within-line direction works — held out, the forward signal is recovered at pass rate
1.000 — but holding out a *pseudo-replicate of the same culture* proves the curve is recoverable,
never that it transfers to a new culture. **Transfer is the untested thing, and it is the thing the
product depends on.**

---

## 4. THE ASK — hard requirements

A published (or obtainable) **single-cell RNA-seq time course of human somatic cells being
reprogrammed toward iPSC**, on **at least one biological source that is not HFF**.

| # | Requirement | Why it is hard, not soft |
|---|---|---|
| **H1** | **Human.** *Homo sapiens*, not mouse. | The aging clock is human (33,155 human gene weights). A mouse time course cannot be scored at all; an ortholog mapping would introduce exactly the kind of confound this request exists to remove. |
| **H2** | **Single-cell RNA-seq** (droplet, 10x or equivalent). | Bulk is the modality already confounded. A bulk dataset adds nothing to what Gill already provides. |
| **H3** | **Whole-transcriptome**, not a targeted panel. | The clock needs genome-wide coverage, and the fate labels are computed from pluripotency / somatic-identity / apoptosis marker programs. |
| **H4** | **An independent biological source** — a different donor and/or cell line from HFF (neonatal human foreskin fibroblast). | This is the entire point. A second HFF experiment does not test transfer. |
| **H5** | **≥ 6 sampled timepoints** spanning day 0 → terminal (iPSC or late). | Forward pairs are (tᵢ → tⱼ) within a trajectory. 3 is the absolute minimum for any pairs; 6 gives 15 pairs; GSE242423's 9 gives 36. Below 6 the trajectory contributes too little. |
| **H6** | **A day-0 (untransduced / uninfected) control sample of the same line.** | Non-negotiable, and it is needed **twice**: it is the zero-point for control-relative ΔAge, *and* the reference distribution the fate labels are z-scored against. A time course starting at day 2 is unusable. |
| **H7** | **≥ ~300 cells per timepoint after QC** (min 500 genes/cell). | From the table in §3: depth is not the blocker, but it does set the smallest effect we could see. At ~470/timepoint a real effect of half HFF's amplitude is caught every time; at 2/timepoint it is missed 57% of the time. |
| **H8** | **Raw or filtered count matrices**, per timepoint, with the day recoverable. | We need counts, not a processed embedding, Seurat object, or normalised expression matrix only. |

---

## 5. STRONGLY PREFERRED (materially increases the value)

| # | Preference | Why |
|---|---|---|
| **P1** | **≥ 3 independent lines/donors in one study.** | With 1 extra line we can test transfer once (2 folds). With 3+ we get leave-one-line-out with a real confidence interval — which is what actually grades the gate rather than giving a single anecdote. **This is the single highest-value item on the page.** |
| **P2** | **Adult donors with known chronological age.** | HFF is neonatal foreskin — donor age 0, which sits **outside the clock's fitted range of [1, 96] years**. Every ΔAge label on 99.7% of our cells is therefore an extrapolation. Adult donors of stated age would fix a *second*, independent weakness. |
| **P3** | **OSKM delivery similar to GSE242423's.** | If the new line also differs in reprogramming method (e.g. Sendai vs lentiviral vs mRNA), we have broken one confound and introduced another. A different method is still usable, but the method must be **stated** so it can be recorded as a caveat. |
| **P4** | **Timepoints on a similar grid** (roughly every 2 days through ~day 14, plus a terminal iPSC sample). | Forward pairs compare Δt across trajectories; wildly different sampling grids reduce the overlap. |
| **P5** | **Multiple donors × the same protocol in one study**, rather than two separate studies. | Cross-study batch effect would be confounded with the line difference all over again. |
| **P6** | 10x Cell Ranger output as `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` with **HGNC symbol in column 2**. | Matches the existing loader exactly; anything else needs a new reader. Not a blocker, just cheaper. |

---

## 6. DEAL-BREAKERS

- **Bulk RNA-seq** of any depth or donor count. Does not address the modality confound.
- **Mouse / MEF** reprogramming, including the well-known Waddington-OT / Schiebinger time course.
- **Differentiation *from* iPSC** (iPSC → neuron, cardiomyocyte, organoid, etc.). This is the
  opposite direction and by far the most common false hit when searching GEO for "iPSC".
- **Direct lineage conversion** (fibroblast → neuron, fibroblast → cardiomyocyte). Not
  reprogramming to pluripotency; different biology and no iPSC endpoint.
- **Trans-differentiation, metabolic reprogramming, immune-cell reprogramming.** The word
  "reprogramming" in GEO titles is dominated by these; ~717 human hits, almost none relevant.
- **A single endpoint** (day 0 vs iPSC only). No intermediate timepoints means no forward pairs.
- **Microarray.** For example `GSE241435` has exactly the right design — fibroblast reprogramming
  sampled at days 0, 3, 5, 7, 14, 21 — but is microarray, so it cannot serve.

---

## 7. WHAT WE DO **NOT** NEED (so nothing is over-specified)

- **No fate/safety annotation.** Labels are computed by us from expression: pluripotency,
  somatic-identity and apoptosis marker programs, z-scored **against that line's own day-0
  controls**. A raw count matrix plus a day label is sufficient.
- **No ΔAge or clock annotation.** Also computed by us.
- **No cell-type annotation, clustering, or embedding.**
- **No matched ATAC, methylation or protein.** Nice, unused.
- **No reprogramming-efficiency or colony-count metadata.**
- **No matched bulk.** We have plenty.

---

## 8. QUESTIONS TO PUT TO THE EXPERT

Ordered so the most valuable answer comes first.

1. **Do you know of a human OSKM reprogramming scRNA-seq time course covering more than one
   donor or cell line?** Multi-line in a single study is worth far more to us than several
   single-line studies (§P1, §P5).
2. **Besides GSE242423, is there *any* published human somatic → iPSC scRNA-seq time course with
   ≥6 sampled days and a day-0 control?** We searched GEO by title and could not verify one — but
   our search would miss anything filed under "iPSC generation", bundled inside a SuperSeries, or
   deposited outside GEO.
3. **Is there anything with adult donors of known chronological age?** Our clock is fitted on ages
   1–96 and our only dense dataset is neonatal, so every age label we have is an extrapolation
   (§P2).
4. **Are there resources outside GEO we should be searching?** Human Cell Atlas, ArrayExpress /
   BioStudies, GSA / CNGBdb, Synapse, Zenodo, or consortium data behind managed access.
5. **Is there an unpublished or in-house dataset that could be shared or collaborated on?**
6. **If nothing suitable exists — what would it take to generate the minimum viable version?**
   Our floor is: 1–2 donor lines distinct from HFF, 6–8 timepoints from day 0 to ~day 14 plus a
   terminal sample, ~3–5k cells per timepoint, standard 10x 3' — i.e. roughly 8–16 10x lanes.
   Rough cost and turnaround would let us decide between acquiring and abandoning the forward tool.
7. **Sanity check on our reasoning:** is it biologically plausible that the *timing* of identity
   loss during OSKM reprogramming transfers across donor lines at all? If the answer is "no, it is
   strongly line- and batch-specific", that is itself decisive — it would mean the forward tool is
   not merely unmeasured but not worth building, and we would stop.

**Question 7 matters as much as question 1.** A well-founded "this does not transfer" saves the
acquisition entirely.

---

## 9. IF THE ANSWER IS "IT DOES NOT EXIST"

Recorded here so the fallback is not reinvented later:

- **(a) Ship the retrospective scoring model without the forward tool.** The fate/safety head is
  unaffected by all of this and scores well — ROC 0.983, PR-AUC 0.992. *Read those two numbers
  with the same caution as everything else on this page:* they are the **mean over the 5 graded
  Gill folds** of `scorecard/baseline.json`, each fold ~20 bulk samples, with two folds reading
  exactly 1.000 and one (Y1) reading 0.932. They also **predate the C-7 label fix** and have not
  been re-measured since. This is the current default, not a validated claim.
- **(b) Report a ΔAge trajectory readout with no safety recommendation** — the recommendation is
  what needs forward prediction; the readout does not.
- **(c) Generate the data** (§8 Q6).
- **(d) Abandon forward prediction.** Legitimate, and cheaper than any of the above if Q7 comes
  back negative.

None of these is chosen yet. What is settled is that the forward gate **cannot be graded on the
data currently in hand**, and that this is a property of the corpus rather than of the estimator.

---

## 10. PROVENANCE OF EVERY NUMBER ABOVE

| Claim | Source |
|---|---|
| 42,481 HFF cells, 9 timepoints, ~4,700/timepoint | `results/stage3a_bis_resolvability_results.json`; probe over the C-7 fold bundles |
| HFF unsafe fraction 0.0835 → 0.9996, per-timepoint SE 0.006 | same, "3a-bis" §"What 3a was actually reading" |
| Gill 1.7 samples per timepoint | same, geometry table |
| The 2×2 pass-rate table | `experiments/stage3a_bis_resolvability.py`, regimes A/B/C/D, 2000 trials each |
| Clock fitted range [1, 96] yr, 33,155 gene weights | `configs/clocks/fleischer_clock.json` → `meta.age_range` |
| HFF donor age 0, asserted not parsed | `src/cellfate/data/sources.py`, `GSE242423SingleCellSource.DONOR_AGE_PROVENANCE` |
| Fate labels are control-relative marker-program z-scores | `src/cellfate/data/labels.py::fate_labels` |
| min 500 genes/cell QC gate | `GSE242423SingleCellSource(min_genes=500)` |
| Fate head ROC 0.983 / PR-AUC 0.992 | mean over the 5 graded folds of `scorecard/baseline.json` (N3 .981/.997, O1 1.0/1.0, O2 1.0/1.0, Y1 .932/.961, Y2 1.0/1.0) — verified from the file, not taken from the prose that quotes it; **pre-C-7** |
