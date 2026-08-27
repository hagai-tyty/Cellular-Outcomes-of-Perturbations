# CellFate-Rx Generation 1 — SCOPE LIMIT

`KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`, recorded by Stage 26 under
`STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md`.

This is the authoritative scope document. The evidence lock, the claim lock and the manuscript are
written against it. Where any other document disagrees with this one, this one governs.

## The system

```text
  WM989 (GSE279162), 1,401 lineage-traced clones, one BRAF-V600E melanoma cell line
  six observed experimental conditions
  endpoint C1, post-treatment clone DETECTION, an observed proxy
  evaluation: clone-held-out, five outer folds, frozen before any result
```

## The vocabulary — closed

```text
  Acid   Cisplatin   CoCl2   Dabrafenib   Doxorubicin   Trametinib
```

Anything else returns `UNSUPPORTED_TREATMENT` and no score. 56 adversarial strings were tried against the
shipped tool -- case variants, whitespace, dose formats, unicode confusables, controls, and sixteen
real oncology drugs including `Vemurafenib`, the drug for this exact mutation, and `Carboplatin`, a
platinum agent one substitution from a condition that IS supported.
**56 of 56 were refused.**

The vocabulary is closed by geometry, not only by a list:

```text
  309 design columns = 50 PCs + 4 nuisance + 5 dummies + 250 interaction terms
```

A seventh condition cannot be added without changing that number, which cannot happen without
refitting.

## What Generation 1 MAY claim

> Within the existing multi-condition WM989 lineage system, pretreatment Gene Expression contains
> treatment-specific information about future clonal detection beyond treatment identity and
> captured pretreatment clone abundance, under frozen clone-held-out evaluation.

And, because Stage 25 recorded `STAGE_25_RANKING_SUPPORTED`:

> A frozen state x treatment model improves clone-specific ordering of the six observed
> experimental conditions over a non-interactive additive model.

```text
  delta_RANK   +0.051605
  CI95         [+0.037197, +0.065571]
  null         0 of 1000 full-refit permutation draws reached the observed value
```

Rewind (GSE227151) may support only:

> A separately reconstructed reprogramming system showed positive but underpowered evidence that
> pretreatment transcriptional state carries prospective information about a later lineage outcome.

## What Generation 1 MAY NOT claim

Each line is written as its own prohibition rather than as an item under a heading, so that no line
can be quoted out of this document and read as a claim.

```text
 1  NEVER  unseen-treatment generalization
 2  NEVER  cross-cell-line or cross-patient generalization
 3  NEVER  clinical treatment recommendation
 4  NEVER  causal treatment-effect estimation
 5  NEVER  calibrated probability
 6  NEVER  independent biological replication of Role B
 7  NEVER  uniform benefit across all six conditions
 8  NEVER  confirmed Role-A prediction
 9  NEVER  single-cell input equivalence
```

Every shipped surface was scanned for all nine. Not one appears except inside a negation.

## What the tool is not applicable to

The nuisance block `B` counts a clone's cells in WM989's three specific naive libraries
(Naive1/2/3). Those libraries are the structure of one experiment, not a property of melanoma.
Data from another lab, cell line or library design cannot produce a valid `B` and cannot be scored.
`clone_input_from_cells` removes a chore for someone working with WM989-structured data; it does
not make the model transferable. That is a Generation-2 modelling change, not packaging.

## Outcome semantics

An observed zero means *no assigned post-treatment cell was observed for that clone-condition row*.
It is not proven death, not sensitivity, not resistance, not clinical response, and not patient
benefit. The tool therefore says `future_detection_score` and `low-persistence condition`.

## Standing limitations, carried on every response

  1. NOT APPLICABLE TO ANOTHER EXPERIMENT. The nuisance block counts a clone's cells in WM989's three specific naive libraries (Naive1/2/3). Those libraries are the structure of one experiment, not a property of melanoma, so data from another lab, cell line or library design cannot produce a valid B and cannot be scored.
  2. known conditions only; UNSUPPORTED_TREATMENT for anything outside the six
  3. requires the complete frozen nuisance block B; it may not be imputed
  4. trained on clone-level pseudobulk, so a single cell is not an equivalent input
  5. captured pretreatment clone abundance remains ~3.45x the state contribution
  6. no calibrated probability; the score is not a calibrated risk
  7. no independent biological replication of the Role-B finding
  8. ranking is NOT validated until Stage 25 records RANKING_SUPPORTED

## Ranking status

`ranking_status = SUPPORTED` only when the Stage-25 verdict file is supplied. Without it the tool
reports `NOT_SUPPORTED` and withholds `validated_condition_order`. **The six scores are identical
either way** -- the verdict unlocks a claim, not a computation.

Validated order for the shipped example clone, lowest predicted detection first:

```text
  Trametinib > Dabrafenib > CoCl2 > Cisplatin > Doxorubicin > Acid
```

This is experimental-condition selection within one benchmark, not a treatment recommendation.

## What Stage 26 does not do

It grants no claim. It records that the existing claim is enforced in the code that ships. No
Stage-26 outcome reopens an earlier stage, changes a recorded number, or authorizes new data, a new
condition or a new model.
