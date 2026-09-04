# CellFate-Rx Generation 1 — CLAIM LOCK

```text
  GEN1_CLAIMS_LOCKED

  evidence lock digest   213593c3b71e7db7064a5bb704288d3d51bddf725d5002aeefaa5b936c0e53b8
  allowed claims         3
  forbidden claims       9
  adversarial sentences  15 of 15 caught
```

Every previous stage produced a number. This one fixes **sentences** — because that is where a
project of this kind actually fails: not in the statistics, but in the abstract, where "six
observed experimental conditions in one melanoma line" becomes "treatments in cancer".

## The ceiling may be lowered. It may not be raised.

---

## Allowed claims

### Primary
> Within the existing multi-condition WM989 lineage system, pretreatment Gene Expression contains treatment-specific information about future clonal detection beyond treatment identity and captured pretreatment clone abundance under frozen clone-held-out evaluation.

### Ranking — selected by STAGE_25_RANKING_SUPPORTED
> A frozen state x treatment model improves clone-specific ordering of the six observed experimental conditions over a non-interactive additive model.

```text
  delta_RANK  +0.051605
  CI95        [+0.037197, +0.065571]
  null        0 of 1000 draws reached the observed value
  p           p < 0.001 (0 of 1,000); never a point estimate
  clones      892 eligible
```

### Supporting — Role A
> A separately reconstructed reprogramming system showed positive but underpowered evidence that pretreatment transcriptional state carries prospective information about a later lineage outcome.

**This claim must travel with** the word SUPPORTING. Rewind does not confirm anything, and its own confirmation gate failed.
Its own evidence: confirmation gate 18.3 FAILED at 0.64; audited to ~0.45.

---

## Mandatory qualifiers

A claim quoted without its qualifier is a claim this lock did not grant.

```text
  system       one BRAF-V600E melanoma cell line, WM989 (GSE279162), 1,401 lineage-traced clones
  vocabulary   six observed experimental conditions: Acid, Cisplatin, CoCl2, Dabrafenib, Doxorubicin, Trametinib
  outcome      an observed post-treatment clone DETECTION proxy -- not death, sensitivity, resistance, clinical response or patient benefit
  evaluation   clone-held-out, five outer folds, frozen before any result existed
  replication  NONE. Independent biological replication is Generation 2 and was not performed.
```

---

## Forbidden claims — the nine, unchanged

```text
  1  NEVER  unseen-treatment generalization
  2  NEVER  cross-cell-line or cross-patient generalization
  3  NEVER  clinical treatment recommendation
  4  NEVER  causal treatment-effect estimation
  5  NEVER  calibrated probability unless calibration is separately frozen and passed
  6  NEVER  independent biological replication of Role B
  7  NEVER  uniform benefit across all six conditions
  8  NEVER  confirmed Role-A prediction
  9  NEVER  single-cell input equivalence when the model was trained on clone pseudobulk
```

Parsed from §3.5 of the frozen ship plan, not from a copy. This lock may add to the list. It may
not subtract from it or reword an entry.

---

## Where the boundary is

A writer told only "no" learns nothing. Each forbidden sentence below was fired at the same scanner
Stage 26 used on the shipped tool, and every one was caught. Each is paired with the nearest
phrasing that is permitted — and every permitted phrasing was itself scanned clean, so this table
never teaches one refused sentence to be swapped for another.

```text
  FORBIDDEN   The model generalises to new treatments.
  PERMITTED   The model is supported for the six observed experimental conditions; anything else returns UNSUPPORTED_TREATMENT.

  FORBIDDEN   Our approach works for any drug.
  PERMITTED   Our approach was evaluated on six observed experimental conditions in one cell line.

  FORBIDDEN   We predict outcomes in cancer cells.
  PERMITTED   We predict a clone-detection outcome in one BRAF-V600E melanoma cell line, WM989.

  FORBIDDEN   Validated on other cell lines.
  PERMITTED   Evaluated clone-held-out within WM989; no other cell line was tested.

  FORBIDDEN   The model predicts patient response.
  PERMITTED   The model predicts an observed post-treatment clone-detection proxy, not clinical response.

  FORBIDDEN   The tool identifies the best treatment for each clone.
  PERMITTED   The tool orders six observed experimental conditions for each clone; this is experimental-condition selection, not treatment recommendation.

  FORBIDDEN   This supports clinical decision-making.
  PERMITTED   This is a research tool for the WM989 experiment and is not a clinical tool.

  FORBIDDEN   A therapeutic recommendation tool for melanoma.
  PERMITTED   A research predictor for six observed conditions in one melanoma cell line.

  FORBIDDEN   We estimate the causal effect of each treatment.
  PERMITTED   We estimate a predictive association between pretreatment state and a later observed outcome, under clone-held-out evaluation.

  FORBIDDEN   The model outputs a calibrated probability of death.
  PERMITTED   The model outputs an uncalibrated future_detection_score; no calibration was frozen or passed.

  FORBIDDEN   Independently replicated in an external cohort.
  PERMITTED   No independent biological replication was performed; clone-held-out folds and two endpoint families are not replication.

  FORBIDDEN   The interaction helps uniformly across all six conditions.
  PERMITTED   Four of six conditions carry meaningful interaction; Cisplatin is negligible on C1 and Doxorubicin is negative on both endpoints.

  FORBIDDEN   The finding was confirmed in a second system.
  PERMITTED   A separately reconstructed reprogramming system gave positive but underpowered supporting evidence; its confirmation gate FAILED.

  FORBIDDEN   Rewind confirms the result.
  PERMITTED   Rewind supports the result without confirming it; gate 18.3 FAILED at 0.64, audited ~0.45.

  FORBIDDEN   The tool scores a single cell.
  PERMITTED   The tool scores a clone-level pseudobulk profile; a single cell is not an equivalent input.
```

---

## A permitted abstract

Assembled only from locked claims and their mandatory qualifiers, and scanned by the same
instrument. A demonstration of the ceiling, not a mandated abstract.

> In one BRAF-V600E melanoma cell line (WM989, 1,401 lineage-traced clones), pretreatment gene expression carries treatment-specific information about a later observed clonal detection outcome, beyond treatment identity and captured pretreatment clone abundance, under clone-held-out evaluation frozen before any result existed. Under a preregistered test, a frozen state-by-treatment interaction model improves clone-specific ordering of the six observed experimental conditions over a non-interactive additive model (delta +0.0516 in equal-clone-weighted within-clone AUROC, 95% CI [+0.0372, +0.0656]; no draw of 1,000 full-refit permutations reached the observed value, p < 0.001). The outcome is an observed post-treatment clone-detection proxy and is not death, sensitivity, resistance or clinical response. The six conditions are the entire supported vocabulary; the model makes no claim about unseen treatments, other cell lines, or patients, and emits no calibrated probability. Independent biological replication was not performed and remains future work.

---

## What locking a claim does not do

It grants nothing. It fixes the ceiling of what may be said about evidence that was locked
separately. No claim-lock outcome reopens an earlier stage, changes a recorded number, or
authorizes new data, a new condition or a new model.
