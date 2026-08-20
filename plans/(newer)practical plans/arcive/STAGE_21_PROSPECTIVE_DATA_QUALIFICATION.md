# STAGE 21 — Prospective Data Qualification

**Status:** PRE-REGISTERED PLAN — freeze before execution.  
**Purpose:** Establish whether CellFate-Rx has access to a scientifically valid, independently grounded, statistically resolvable prospective dataset before any new prospective model is trained.  
**Current model:** `_s16` remains frozen.  
**This stage does not modify `src/`.**  
**Primary output:** one qualified prospective dataset, or a documented conclusion that no currently available dataset is adequate.

---

## 1. Why Stage 21 exists

The original CellFate-Rx goal is prospective:

> Given an initial molecular state `X` and a proposed perturbation `U`, predict a later cellular outcome `Y`.

The current project has not yet demonstrated that claim.

Two structural facts matter:

1. The current fate target `y_cls` is derived from transcriptomic signature scores. It is therefore a molecular surrogate, not an independently observed future biological outcome.
2. The datasets already used by the project are not obviously capable of forming a strong prospective training/evaluation task:
   - `GSE242423` is one HFF timecourse with many cells but no demonstrated lineage linkage across time.
   - `GSE165177` is bulk and donor-limited.

Stage 21 therefore moves **data qualification onto the critical path**.

The stage asks:

> **What prospective dataset can actually support the paper-making experiment?**

Only after that question is answered do we train a prospective baseline or change CellFate-Rx.

---

# 2. Core distinction: forecasting is not automatically independent fate prediction

The stage must preserve three different tasks.

### Same-state reconstruction

```text
X_t -> f(X_t)
```

This is not prospective. If `f` is directly derived from the same expression vector, a high score can be circular or reconstructive.

### Future molecular-state forecasting

```text
X_t -> f(X_(t+k))
```

This is a legitimate temporal forecast.

However, if `f(X_(t+k))` is itself only an expression-derived fate signature, the result supports:

> **prediction of a later transcriptomic fate surrogate**

not:

> **independently validated biological fate prediction**

### Desired prospective biological-outcome task

```text
X_t + U -> Y_(t+k)
```

where `Y_(t+k)` is measured independently of the RNA vector used as input.

Examples of potentially qualifying outcomes include:

- clone survives / disappears;
- colony forms / does not form;
- later experimentally measured reprogramming success;
- viability or survival assay;
- imaging phenotype;
- orthogonal sorting / surface phenotype;
- lineage-resolved terminal fate.

This third geometry is the target for the main prospective paper.

---

# 3. Stage 21 structure

Stage 21 has three parts.

```text
STAGE 21A
Audit existing local data
        ↓

STAGE 21B
Qualify public prospective datasets
        ↓

STAGE 21C
Acquire + reconstruct + reproduce one qualified dataset
        ↓

PASS
Qualified prospective dataset ready for Stage 22
```

Local-data modelling is **not** the default next step.

If Stage 21A confirms that the local data cannot support the final prospective claim, Stage 21 proceeds immediately to public-data qualification.

---

# 4. Stage 21A — Local data audit

## 4.1 Goal

Turn the current qualitative understanding into a reproducible artifact.

Do not fit a model.

The audit must determine, for each local dataset:

- what the early biological unit is;
- whether early RNA is measured before the outcome;
- whether a real early→later linkage exists;
- whether treatment / dose / time varies;
- what the later outcome actually is;
- whether that outcome is independent of RNA;
- how many independent biological units exist;
- whether a statistically meaningful prospective evaluation is even possible.

---

## 4.2 Local datasets

Audit at minimum:

```text
GSE242423
GSE165176
GSE165177
GSE165178
GSE165179
GSE113957
GSE297234
```

The audit must not assume that thousands of cells imply thousands of independent trajectories.

---

## 4.3 Required classification — linkage

Every dataset receives one strongest linkage level.

### LINK 3 — `LINEAGE_LINKED`

Early transcriptome and later outcome are connected by a true lineage / clone / sister-cell identifier.

### LINK 2 — `CULTURE_LINKED`

Early culture or biological replicate is explicitly followed to a later outcome.

### LINK 1 — `TRAJECTORY_ONLY`

Early and late populations come from the same broad trajectory, but individual lineage/culture linkage is absent.

### LINK 0 — `NO_FORWARD_LINK`

No legitimate early→later relation can be constructed.

---

## 4.4 Required classification — outcome truth

Every candidate target also receives one truth level.

### TRUTH 2 — `ORTHOGONAL_OUTCOME`

The later outcome is measured independently of the transcriptomic vector used as the prediction target.

### TRUTH 1 — `EXPRESSION_DERIVED_SURROGATE`

The later target is a deterministic or engineered function of later RNA.

Example:

```text
future_y_cls = fate_labels(X_future)
```

This may support a molecular forecast, but not independently validated biological fate.

### TRUTH 0 — `NO_USABLE_OUTCOME`

No defensible later target exists.

---

## 4.5 Required audit table

For each local dataset report:

```text
dataset
modality
donors
cell_lines
biological_replicates
timepoints
treatments
dose_information

early_rna_available
lineage_id_available
clone_id_available
culture_followup_available

future_target_candidate
future_target_source
future_target_is_expression_derived
future_target_is_orthogonal

linkage_level
truth_level

outer_split_unit
n_independent_outer_units
mixed_outcome_units
minimum_attainable_p
resolvable_for_pass

licensed_claim
verdict
```

---

# 5. Resolvability gate

Stage 21 must not spend serious modelling effort on an experiment that cannot change the scientific decision.

The effective `n` is the number of **independent biological outer units**, not the number of cells.

Examples:

```text
independent donors
independent cultures
independent clones
independent trajectories
```

For a simple two-sided sign-flip comparison across `n` independent units:

```text
p_min = 2 / 2^n
```

So:

```text
n = 3 -> p_min = 0.25
n = 4 -> p_min = 0.125
n = 5 -> p_min = 0.0625
n = 6 -> p_min = 0.03125
```

If the planned primary test cannot theoretically cross its frozen evidence bar, mark:

```text
UNRESOLVABLE_FOR_PASS
```

and do not proceed to a full local modelling exercise.

---

# 6. Expected local-data outcome — recorded before the audit

These are predictions, not verdicts.

## `GSE242423`

Expected:

```text
LINK 1
TRUTH 1 at best for current y_cls
one HFF trajectory
no true lineage linkage
no independent replicate cultures
```

Likely conclusion:

```text
TRAJECTORY_PILOT_ONLY
```

## `GSE165177`

Expected:

```text
bulk RNA
3 donors
culture/sample-level information only
current fate proxy unsuitable as per-cell truth
```

Likely conclusion:

```text
UNRESOLVABLE_FOR_PASS
```

## `GSE165176`

Experimental arm metadata such as reprogramming / failing-to-reprogram labels may exist.

Stage 21A must explicitly determine whether those labels can be attached to an **earlier** molecular state without future-arm leakage.

Do not assume that an arm label automatically creates a prospective target.

---

# 7. Stage 21A verdicts

## `LOCAL_ORTHOGONAL_TASK_READY`

Requirements:

```text
LINK >= 2
TRUTH = 2
RESOLVABLE_FOR_PASS
```

If this unexpectedly occurs, the dataset may move directly into Stage 21C reconstruction.

## `LOCAL_SURROGATE_ONLY`

A valid future molecular surrogate exists, but:

```text
TRUTH = 1
```

This can support a feasibility forecast, not the main biological-fate claim.

## `LOCAL_TASK_UNRESOLVABLE`

A forward task exists but the independent-unit count cannot support the frozen statistical decision.

## `NO_VALID_LOCAL_FORWARD_TASK`

No defensible local prospective task exists.

For the last three outcomes, continue immediately to Stage 21B.

---

# 8. Stage 21B — Public prospective-dataset qualification

## 8.1 Goal

Find a public dataset that actually supplies the missing geometry:

```text
early molecular state X
+
known perturbation U
+
valid lineage / culture linkage
+
later independently measured biological outcome Y
+
enough independent units for evaluation
```

This is now the main second half of Stage 21, not a fallback.

---

# 9. First candidate set

Audit these first:

```text
Rewind / reprogramming lineage data
ReSisTrace
CellTag / CellTag-multi
other lineage-resolved reprogramming datasets discovered during qualification
```

The exact accession and files must be verified during Stage 21B.

No dataset is accepted because of its paper abstract alone.

---

# 10. Candidate priority

## Priority 1 — Reprogramming + orthogonal outcome

Best case:

```text
pre-treatment fibroblast transcriptome
+
OKSM / reprogramming condition
->
later experimentally measured reprogramming success
```

with clone or lineage linkage.

This most directly connects to CellFate-Rx's biological use case.

## Priority 2 — Multiple perturbations + orthogonal response

Best case:

```text
pre-treatment transcriptome
+
treatment identity
->
future survival / resistance / fate
```

with clone/sister-cell linkage.

This is especially valuable for testing the original `X + U` thesis.

## Priority 3 — Lineage-linked but transcriptomic endpoint

A dataset with excellent lineage linkage but a later outcome still defined mainly from transcriptomic state is useful, but it remains a surrogate benchmark rather than the final orthogonal-fate dataset.

---

# 11. Mandatory public-dataset qualification table

For every candidate:

```text
dataset_name
accession
species
cell_type
biological_system

pre_outcome_rna_available
pre_outcome_time
treatment_identity_available
dose_available

lineage_or_clone_link_available
linkage_method

future_outcome_available
future_outcome_measurement_method
outcome_is_orthogonal_to_rna

n_independent_experiments
n_independent_clones_or_lineages
n_positive_outcomes
n_negative_outcomes

multiple_treatments
multiple_doses
multiple_donors

processed_files_available
raw_files_available
barcode_mapping_available
metadata_join_possible

outer_split_unit
minimum_attainable_p
resolvable_for_pass

paper_relevance
qualification_verdict
```

---

# 12. Stage 21B qualification requirements

A dataset is **QUALIFIED_FOR_PRIMARY** only if all required fields pass.

### Q1 — Temporal validity

RNA is measured before the later outcome.

### Q2 — Biological linkage

The early molecular state can be linked to the later outcome by clone, lineage, sister-cell, or tracked culture.

### Q3 — Orthogonal future outcome

The later endpoint is not merely `f(X_future)`.

### Q4 — Prediction-time validity

No future information is already encoded in the early metadata in a way that makes the task trivial.

### Q5 — Independent units

There are enough independent biological units for the frozen evaluation to be resolvable.

### Q6 — Reconstructability

The public files contain the identifiers needed to actually build the `X + U -> Y` table.

### Q7 — Outcome variation

Both outcomes occur in useful numbers.

A dataset that fails one required criterion is not silently repaired by relaxing the criterion.

---

# 13. Qualification classes

## `QUALIFIED_PRIMARY`

Suitable for the main prospective experiment.

Expected shape:

```text
LINK 3 or strong LINK 2
TRUTH 2
resolvable
joinable
```

## `QUALIFIED_REPLICATION`

Useful as a second system but not ideal as the first dataset.

## `QUALIFIED_SURROGATE`

Strong prospective lineage geometry but outcome is still transcriptomic / molecular.

Useful for method development, not sufficient for the final biological-fate headline.

## `NOT_QUALIFIED`

Fails required geometry, truth, or reconstruction.

---

# 14. Do not train during Stage 21B

Stage 21B is a data qualification stage.

Forbidden:

```text
CellFateNet training
logistic-regression performance fishing
gene-selection tuning
trying multiple outcomes until one predicts well
```

The dataset is accepted or rejected based on experimental geometry and resolvability, not on whether a model score looks attractive.

---

# 15. Stage 21C — Acquire and reconstruct one dataset

Once one candidate is `QUALIFIED_PRIMARY`, acquire it.

Stage 21C must:

1. download the exact required public files;
2. checksum / record provenance;
3. reconstruct the clone / lineage mapping;
4. reconstruct treatment metadata;
5. reconstruct the future outcome;
6. produce one frozen analysis table;
7. reproduce one central published descriptive result from the source study;
8. verify no leakage in the join;
9. freeze the final task definition for Stage 22.

---

# 16. Required Stage-21C analysis table

The final table passed to Stage 22 should look conceptually like:

```text
unit_id
clone_or_lineage_id
biological_replicate
donor_or_cell_line

X_time
X_1
X_2
...
X_p

treatment_id
dose
exposure_time

Y_time
future_outcome
future_outcome_measurement

split_group
```

The exact schema may vary by dataset, but the biological meaning cannot.

---

# 17. Reproduce before modelling

Before CellFate-Rx is allowed to use the dataset, reproduce at least one source-study result that verifies the reconstruction.

Examples:

```text
reported number of successful clones
reported resistant / sensitive split
reported reprogramming-success enrichment
reported barcode/lineage counts
```

The purpose is not to reproduce the whole paper.

It is to prove that our metadata/barcode/outcome join corresponds to the experiment the original authors actually ran.

---

# 18. Leakage audit

The reconstructed dataset must hard-fail if:

```text
future outcome leaks into X
future cluster label appears in metadata
post-treatment expression enters the early input
same clone appears in train and test when clone is the split unit
same biological replicate appears in both sides when replicate is the split unit
treatment label is reconstructed from the outcome
```

Every guard must have a synthetic regression test.

---

# 19. Stage 21 PASS condition

Stage 21 passes only when:

```text
1. one public/local dataset is QUALIFIED_PRIMARY;
2. its required files are acquired;
3. the early→later biological linkage is reconstructed;
4. the future endpoint is independently measured;
5. the task is statistically resolvable;
6. one source-study sanity result is reproduced;
7. the final Stage-22 table is frozen;
8. leakage guards pass;
9. src/ remains unchanged.
```

Console:

```text
STAGE 21 — PROSPECTIVE DATA QUALIFICATION
=========================================

LOCAL DATA:
  final prospective claim supported: NO / YES

PUBLIC QUALIFICATION:
  candidates audited: ...
  qualified primary: ...

SELECTED DATASET:
  ...
  linkage: ...
  outcome truth: ORTHOGONAL
  independent units: ...
  resolvable: YES
  reconstruction: PASS
  source-study sanity check: PASS

VERDICT: PASS
NEXT: STAGE 22 — PROSPECTIVE BASELINE EXPERIMENT
```

---

# 20. Stage 21 FAIL / BLOCKED condition

If no candidate qualifies:

```text
VERDICT: DATA_BLOCKED
```

The result must state exactly which requirement fails most often:

```text
no pre-outcome RNA
no lineage link
no independent outcome
too few independent units
missing barcode mapping
no treatment variation
outcome not reconstructable
```

At that point the project knows what collaboration or new experiment is required.

No model development is allowed to substitute for missing experimental geometry.

---

# 21. Files owned by Stage 21

Suggested:

```text
plans/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION.md

experiments/diag_stage21_local_audit.py
tests/test_diag_stage21_local_audit.py
results/diag_stage21_local_audit.json

experiments/diag_stage21_public_qualification.py
tests/test_diag_stage21_public_qualification.py
results/diag_stage21_public_qualification.json

experiments/repro_stage21_selected_dataset.py
tests/test_repro_stage21_selected_dataset.py
results/repro_stage21_selected_dataset.json

data/stage21/<selected_dataset>/...
```

Do not create model-training code in this stage.

---

# 22. What Stage 21 does NOT prove

A PASS means:

> **The project now possesses a legitimate and reconstructable dataset on which prospective fate prediction can be tested.**

It does not prove:

- that RNA predicts the outcome;
- that treatment adds information;
- that CellFate-Rx beats simple baselines;
- held-out-treatment generalization;
- rejuvenation;
- working RES;
- future ΔAge accuracy;
- clinical utility.

Those questions begin at Stage 22.

---

# 23. Why this stage is now the critical path

The current project's limiting factor is no longer primarily architecture.

The missing object is the right experiment:

```text
X_before + U -> independently measured Y_after
```

The system already contains:

- transcriptome encoding;
- perturbation encoding;
- calibration machinery;
- OOD machinery;
- deterministic evaluation;
- extensive audit infrastructure.

The next scientific bottleneck is therefore whether we can attach those tools to a valid prospective ground truth.

That is what Stage 21 must settle.

---

# 24. One-line definition

> **Stage 21 qualifies and reconstructs the prospective dataset that the original CellFate-Rx claim actually requires; local data are audited first, then public lineage-resolved data become the main path, and no prospective model is trained until one independently grounded future-outcome dataset is ready.**
