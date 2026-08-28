# Pretreatment transcriptional state carries condition-specific information about future clonal detection in a lineage-traced melanoma line

**CellFate-Rx, Generation 1.**

```text
  evidence lock   455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9
  claim lock      0b3c7f038a41c3b58b4c47cc5769d1e5e7be27e4a01d6c5f790a8e7da296ae5f
```

Both digests are verifiable from the repository. See **Availability**.

---

## Abstract

Whether a cell's molecular state before a perturbation predicts what happens to it afterwards is
usually asked retrospectively, after outcome and state have been measured in the same cells. We ask
it prospectively at clone level, in one BRAF-V600E melanoma cell line (WM989, GSE279162) in which
1,401 barcoded clones were split and exposed to six observed experimental conditions: Acid,
Cisplatin, CoCl2, Dabrafenib, Doxorubicin and Trametinib.

Within this system, pretreatment gene expression contains condition-specific information about
future clonal detection beyond condition identity and captured pretreatment clone abundance, under
clone-held-out evaluation frozen before any result existed. Under a test preregistered in full — the
metric, population, weighting, comparator, null and verdict rule all fixed before the numbers
existed — a frozen state-by-condition interaction model improves clone-specific ordering of the six
conditions over a non-interactive additive model: +0.051605 in equal-clone-weighted within-clone
AUROC, 95% CI [+0.037197, +0.065571], with 0 of 1000 full-refit permutation draws reaching the
observed value (p < 0.001).

The outcome is an observed post-treatment clone-detection proxy and is **not death**, sensitivity,
resistance or clinical response. The six conditions are the entire supported vocabulary. We make no
claim about unseen conditions, other cell lines, or patients, and the model emits no calibrated
probability. Independent biological replication has not been performed and is Generation 2 work,
not a gate on this result.

---

## Introduction

The question is easy to state and hard to evaluate honestly: before a perturbation is applied, does
a cell's transcriptional state say anything about which perturbation it will survive?

Most evidence for state-dependent outcomes is retrospective. Outcome and state are read from the
same cells, or state is inferred from cells that already experienced the perturbation, and the
apparent predictive signal can be a consequence of the outcome rather than a cause. Barcoded
lineage tracing offers a way around this: a clone is split, part of it is profiled before anything
is done to it, and the remaining parts are exposed to different conditions. The pretreatment
profile is then genuinely prior to the outcome.

Two things make this hard to evaluate rather than merely to observe. The first is that clone
abundance dominates. A clone that was large before treatment is more likely to be detected after
it, for reasons that have nothing to do with state, so any comparison that does not hold abundance
fixed will find a signal that is really a headcount. The second is that the interesting claim is
not "state predicts outcome" but "state predicts *which condition* this clone survives" — a
clone-specific ordering claim, which a model can only satisfy through an explicit interaction
between state and condition.

This work evaluates exactly that, once, under a protocol frozen in advance.

---

## Data

```text
ROLE B, primary       GSE279162 (WM989)      benchmark, model, tool and ranking analysis
ROLE A, supporting    GSE227151 (Rewind)     historical supporting evidence only
```

Role B is one BRAF-V600E melanoma cell line with 1,401 lineage-traced clones split across six
observed experimental conditions. It carries the whole primary claim.

Role A is a separately reconstructed reprogramming system. It contributes one supporting sentence
and nothing else; its own confirmation gate failed. It is not a replication of Role B, and it does
not provide the same multi-condition task or the same outcome.

No additional dataset was searched, downloaded, qualified or used. Raw sequencing data is not
vendored; accessions are given above.

---

## Methods

### Benchmark construction

Clone assignments, outcome construction, the five outer folds, feature rules, condition aliases and
exclusions were all fixed before any model was fitted, and were not revised afterwards.

Two endpoint families were built. **C1** is post-treatment clone detection: an observed zero means
no assigned post-treatment cell was seen for that clone-condition row. **C2** is a clone-balanced
abundance endpoint. Everything reported below is C1.

Expression is clone-level pseudobulk: raw pretreatment counts are summed over a clone's cells, then
CP10K-normalised and `log1p`-transformed exactly once. Applying the transform twice, or summing
already-normalised cells, produces a different feature space and a model the benchmark never
evaluated.

### Models

```text
  W1   B + U                nuisance and condition identity only
  W4   X + B + U            plus an additive expression term
  W5   X + B + U + X*U      plus an explicit state-by-condition interaction
```

`X` is the clone expression profile reduced to 50 principal components on a train-only basis. `B`
is the captured-abundance nuisance block: `log1p` cell counts, total and per pretreatment library.
`U` is condition identity, encoded as five non-reference indicators with Acid as reference. The
design has 309 columns: 50 components, 4 nuisance terms, 5 indicators and 250 interaction terms.

`B` is not optional. Abundance is the confounder that makes a naive version of this question
trivial, so it is in every model including the null, and the tool refuses to score without it.

### Evaluation

Five outer folds, held out **by clone**, so no model ever scores a clone it trained on. Gene
filtering, the PCA basis and every scaler are refitted inside each training fold. Hyperparameters
are selected by an inner grouped split within the training folds only.

### The preregistered ranking test

The primary question is whether W5 orders the six conditions *within a clone* better than W4. The
metric is equal-clone-weighted within-clone AUROC: for each clone, the mean over positive/zero
condition pairs with ties scoring exactly 0.5; then a plain mean over clones, so a clone
contributing many pairs cannot outweigh one contributing few.

The evaluable population is clones with at least one detected and one undetected condition; a clone
with no contrast has an undefined within-clone AUROC. This yields 892 of 1,401 clones — 472 were
never detected under any condition and 37 were always detected.

The comparator is W4, not W1. An additive expression term cannot by construction create
clone-specific ordering, so W4 isolates the interaction as the only thing that could.

The null is a **full refit**: profiles are permuted within stratum on each side of the outer-fold
boundary and never across it, and the entire pipeline — filtering, PCA, scaler fitting and
hyperparameter selection — is re-run inside every draw. Observed-data hyperparameters are never
reused. 1000 draws, no early stopping, with a completeness assertion that refuses an incomplete
null rather than silently reporting a smaller one.

Uncertainty on the observed statistic is a 2,000-replicate clone bootstrap, conditional on the
fitted models.

---

## Results

### The interaction improves clone-specific ordering

```text
  R(W1)   0.692654      nuisance + condition
  R(W4)   0.692176      + additive X
  R(W5)   0.743781      + explicit X x U

  delta_RANK   +0.051605     CI95 [+0.037197, +0.065571]
```

`R(W4)` sits *below* `R(W1)` by 0.0005. The additive expression term contributes nothing to
ordering, which is precisely why W4 was preregistered as the comparator. The entire ordering gain is
the interaction.

### The separation, not the p-value, is the result

```text
  null p95                    0.008672
  observed / null p95         6.0x
  observed vs null mean       11.8 null standard deviations above it
  largest of 1,000 draws      0.013722      -- the observed value exceeds EVERY null draw
  draws reaching observed     0 of 1000 full-refit permutations
```

`p < 0.001` is the floor of a 1,000-draw permutation test and is reported as such, never as a point
estimate. The number that carries weight is the separation: nothing the null produced came close.

### It holds in every stratum it was broken down by

```text
  BY OUTER FOLD                     BY PRETREATMENT DEPTH
    fold 0   +0.0435                  1 cell     +0.0535
    fold 1   +0.0548                  2 cells    +0.0528
    fold 2   +0.0658                  3-4        +0.0314
    fold 3   +0.0506                  5-9        +0.0462
    fold 4   +0.0435                  10+        +0.0779
```

Positive in all five folds and all five depth strata. These were preregistered as descriptive and
could not have rescued a failed primary gate; they were not asked to.

### Choosing the lowest-scoring condition

```text
  delta_TOP1   +0.115471    CI95 [+0.082960, +0.145740]
```

Selecting each clone's lowest predicted detection score finds a genuine zero for 82.8% of evaluable
clones under W5 against 71.3% under W4. This was preregistered as a directional-consistency check,
not a significance test: it could withhold support, never grant it. It did not withhold.

---

## The tool

A frozen predictor ships with this work. For one starting clone it returns a
`future_detection_score` for each of the six observed conditions, reproducing the frozen
out-of-fold predictions to within 5e-16.

What it refuses is as much of the specification as what it returns:

```text
  an unknown condition        -> UNSUPPORTED_TREATMENT, and no score
  a missing nuisance block    -> MISSING_REQUIRED_NUISANCE, never imputed
  a wrong feature schema      -> UNSUPPORTED_FEATURE_SCHEMA
  a validated ordering        -> withheld unless the preregistered verdict file is supplied
```

The refusal was tested adversarially rather than asserted. 56 of 56 hostile condition strings were
refused: case and whitespace variants, dose formats, unicode look-alikes, controls, and sixteen real
oncology drugs including Vemurafenib — the drug for this exact mutation — and Carboplatin, one
substitution from a condition that *is* supported.

That test found something worth stating. Acid is the reference level and is encoded as five zero
indicators, so the indicator encoder alone maps *any* unrecognised string to the Acid row and would
return the Acid score under another name. The vocabulary filter is the only thing preventing this,
and it was verified to hold rather than assumed to.

---

## Limitations

Carried verbatim from the preregistered verdict.

```text
  1  No independent biological replication of the Role-B finding. Clone-held-out folds and
     two endpoint families are not replication.
  2  Captured pretreatment clone abundance remains ~3.45x the whole state contribution.
     The ordering is abundance first, then condition-specific state.
  3  Four of six conditions carry meaningful interaction. Cisplatin is negligible on C1
     and Doxorubicin is negative on both endpoint families.
  4  C1 is an observed detection proxy -- not death, sensitivity or clinical response.
  5  The bootstrap interval is conditional on the fitted models; only the null refits.
  6  Role A remains positive-but-underpowered supporting evidence. Its confirmation gate
     18.3 FAILED at 0.64 against a 0.80 threshold, and a later audit of the instrument
     put the true power at 0.45 -- lower still. Its effect size must not be quoted as an
     estimate, because an underpowered design that reaches significance inflates it.
```

Limitation 2 is the one most likely to be misread. This work does not show that state dominates
outcome. It shows that state adds something specific, on top of an abundance term that is several
times larger.

Limitation 6 is the honest position on Role A. Its own gate failed, we audited our own power
calculation and found it had been too generous, and we report the worse number.

---

## What this does not show

Separate from the limitations above, these are claims this work may not make, in any form.

```text
  1  NEVER  unseen-condition generalization
  2  NEVER  cross-cell-line or cross-patient generalization
  3  NEVER  clinical treatment recommendation
  4  NEVER  causal treatment-effect estimation
  5  NEVER  a calibrated probability
  6  NEVER  independent biological replication of Role B
  7  NEVER  uniform benefit across the six conditions
  8  NEVER  a confirmed Role-A result
  9  NEVER  single-cell input equivalence, the model being trained on clone pseudobulk
```

The result is bounded to one cell line, six observed experimental conditions, and one observed
detection proxy. It is not evidence about therapy, and the six conditions include non-clinical
stress contexts that no one would administer to anything.

---

## Availability

### Verify before reading anything else

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
```

The first re-hashes every locked artifact and refuses if one has moved. The second does the same
for the claim set. Both were shown to refuse a one-bit change before either was issued.

```text
  evidence lock digest   455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9
  claim lock digest      0b3c7f038a41c3b58b4c47cc5769d1e5e7be27e4a01d6c5f790a8e7da296ae5f
```

### What is in the repository

Benchmark tables, frozen out-of-fold predictions, the serialized model metadata, the prediction API
and CLI, the model card and schema, every stage protocol, every stage record, and every executor
and contract file. Full inventory and per-file hashes: `results/evidence_lock/`.

### What is not

```text
  stage24_w5_artifact.npz   44 MB, gitignored. A fresh clone does NOT contain it. Its hash
                            is locked and it rebuilds in about half a minute:
                              python experiments/run_stage24_gen1_tool.py --stage 24c
  raw sequencing data       GSE279162, GSE227151. Accessions are locked; bytes are not
                            vendored.
```

Naming a gap is not closing it. Both remain open.

Full reproduction instructions, environment and runtimes: `results/manuscript/REPRODUCIBILITY.md`.

---

## Generation 2

What would actually test this, none of which was a gate on the present result:

```text
  independent biological replication in a different lineage-traced system
  transfer to conditions the model has never seen, which requires a condition
    representation the present design does not have
  a dataset-independent nuisance definition, since the current abundance block counts
    cells in this experiment's three specific pretreatment libraries and therefore
    cannot be computed anywhere else
  calibration, frozen and tested separately
  out-of-distribution behaviour
```

The first is the one that matters. Everything here rests on a single lineage-traced system, and one
system is one system however carefully it is evaluated.
