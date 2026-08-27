# CellFate-Rx Gen-1 — Model Card

**Model** `gen1-w5-c1-v1` — W5 = `X + B + U + X*U`
**Feature contract** `wm989-ge-36601-v1`
**Endpoint** C1, post-treatment clone detection
**System** WM989 (GSE279162), 1401 lineage-traced clones, six conditions

## What it does
For one starting clone it returns a `future_detection_score` for each of the six observed
conditions: the model's propensity that the clone is still detected after that condition.

## What it does NOT do
```text
  NOT APPLICABLE TO ANOTHER EXPERIMENT. The nuisance block counts a clone's cells in WM989's three specific naive libraries (Naive1/2/3). Those libraries are the structure of one experiment, not a property of melanoma, so data from another lab, cell line or library design cannot produce a valid B and cannot be scored.
  known conditions only; UNSUPPORTED_TREATMENT for anything outside the six
  requires the complete frozen nuisance block B; it may not be imputed
  trained on clone-level pseudobulk, so a single cell is not an equivalent input
  captured pretreatment clone abundance remains ~3.45x the state contribution
  no calibrated probability; the score is not a calibrated risk
  no independent biological replication of the Role-B finding
  ranking is NOT validated until Stage 25 records RANKING_SUPPORTED
```

## Validation
```text
  frozen out-of-fold reproduction   max |diff| vs the Stage-23 frozen column   2.498e-16
  determinism, same session         True
  determinism, across loads         True
  fold isolation                    every fold component verified disjoint from its test clones
  eligible ranking clones           892 of 1401
```

The **deployment** component is packaging, not validation: it is the same specification and the same
selection rule fitted once on all clones, and it is **not** validated on held-out data. Its
performance is estimated by the frozen out-of-fold result, which came from the fold components.

## Ranking
`ranking_status` is `NOT_SUPPORTED` until Stage 25 records `STAGE_25_RANKING_SUPPORTED` under its
pre-registered test. Until then the six scores are returned but their **order is not a validated
condition ranking** and `validated_condition_order` is withheld.

## Intended use
Research use on **the WM989 experiment itself** -- reproducing its frozen out-of-fold predictions,
or scoring a clone from that experiment. The nuisance block counts cells per WM989 naive library
(Naive1/2/3), so **data from another experiment cannot produce a valid B and cannot be scored**.
Making the model transferable would require a dataset-independent nuisance definition, which is a
Generation-2 modelling change, not packaging.

Supported for the six observed experimental conditions. **Not** a clinical tool, **not** a treatment recommendation, **not** a calibrated
probability, and **not** applicable to unseen treatments, other cell lines, or patients.


<!-- STAGE-26 SCOPE LOCK -- appended, nothing above this line altered -->

## Stage 26 — scope lock (STAGE_25_RANKING_SUPPORTED)

Stage 25 has since run. Its preregistered ranking test recorded **`STAGE_25_RANKING_SUPPORTED`**:
`delta_RANK = +0.051605`, CI95
[+0.037197, +0.065571],
0 of 1000 full-refit
permutation draws reached the observed value. So `validated_condition_order` **is** exposed --
but only when the Stage-25 verdict file is supplied to `Gen1Predictor.load`, and the six scores
are identical whether or not it is. The verdict unlocks a claim, not a computation.

The **Ranking** section above was written before that run and speaks in the future tense. It is
left byte-for-byte as it was frozen; this section is the current state.

### The vocabulary is closed, and it was attacked

Stage 26 fired adversarial condition strings at this tool -- case variants, whitespace, dose
formats, unicode confusables, controls, and real oncology drugs including `Vemurafenib`, the drug
for this exact BRAF-V600E mutation, and `Carboplatin`, one substitution from a supported
condition. Every one was refused with `UNSUPPORTED_TREATMENT` and no score.

This matters more than it looks. `Acid` is the reference level and is encoded as five zero
dummies, so the dummy encoder on its own maps **any** unknown string to the Acid row and would
return the Acid score under another name. Stage 26 confirmed that hazard is real and confirmed the
vocabulary filter in `predict()` blocks it.

### Scope

Authoritative scope document: `results/stage26/GEN1_SCOPE_LIMIT.md`. Nine claims are forbidden in
Generation 1, each written as its own prohibition so that no line can be quoted out of context and
read as a claim:

```text
  NEVER  unseen-treatment generalization
  NEVER  cross-cell-line or cross-patient generalization
  NEVER  clinical treatment recommendation
  NEVER  causal treatment-effect estimation
  NEVER  calibrated probability
  NEVER  independent biological replication of Role B
  NEVER  uniform benefit across all six conditions
  NEVER  confirmed Role-A prediction
  NEVER  single-cell input equivalence
```

Every shipped surface was scanned for all nine. Not one appears except inside a negation.
