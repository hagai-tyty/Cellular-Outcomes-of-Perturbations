# STAGE 21 — Prospective Data Qualification v2

**Status:** PRE-REGISTERED PLAN — freeze before execution.  
**Purpose:** Acquire and qualify the data geometry required for a genuinely new prospective CellFate-Rx paper.  
**Current model:** `_s16` remains frozen.  
**No model training in this stage. No `src/` changes.**

---

## 0. The strategic change

Stage 21 is no longer a local-data feasibility exercise whose fallback is acquisition.

**Acquisition / public-data qualification is now the critical path.**

The paper we want is not:

> "CellFate-Rx can classify another transcriptomic dataset."

The paper we want is built around a new capability:

```text
CURRENT MOLECULAR STATE X
          +
PROPOSED PERTURBATION U
          ↓
INDEPENDENTLY OBSERVED FUTURE OUTCOME Y
```

and, where the data allow it:

```text
same / comparable starting state
        ↓
candidate treatment U1 -> predicted outcome
candidate treatment U2 -> predicted outcome
candidate treatment U3 -> predicted outcome
        ↓
rank which treatment is best for that state
```

Stage 21 therefore exists to secure the experimental ground truth needed to build and test that system.

---

# 1. Why the current data are not enough

The current fate target `y_cls` is produced from transcriptomic signature scores.

That means a same-timepoint prediction:

```text
X_t -> fate_labels(X_t)
```

is not an independently grounded biological-outcome experiment.

A forward version:

```text
X_t -> fate_labels(X_(t+k))
```

is a legitimate temporal molecular forecast, but still only predicts a **future transcriptomic surrogate**.

For the main paper we require:

```text
X_t + U -> Y_(t+k)
```

where `Y_(t+k)` is obtained independently of the RNA vector used as the input/target.

Examples include:

```text
clone survives / disappears
colony forms / does not form
experimentally measured reprogramming success
future viability
independent imaging phenotype
orthogonal cell-surface / sorting phenotype
lineage-resolved terminal outcome
```

---

# 2. Stage 21 has three phases

```text
STAGE 21A
LOCAL DATA AUDIT
        ↓

STAGE 21B
PUBLIC DATA QUALIFICATION
        ↓

STAGE 21C
DOWNLOAD + RECONSTRUCT + VERIFY
        ↓

PASS
READY FOR STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

Stage 21 does **not** fit predictive models.

---

# 3. Stage 21A — Local data audit

Run the existing-data audit once to convert the current understanding into a reproducible result.

Audit:

```text
GSE242423
GSE165176
GSE165177
GSE165178
GSE165179
GSE113957
GSE297234
```

For each dataset record:

```text
pre-outcome RNA?
real early -> late linkage?
lineage / clone id?
tracked culture?
orthogonal future outcome?
treatment variation?
dose variation?
independent biological units?
statistically resolvable?
```

Do not fit a prospective model here unless the audit unexpectedly discovers a strong orthogonal task.

---

# 4. Linkage classes

## LINK 3 — `LINEAGE_LINKED`

Early molecular state and future outcome are joined by a real clone, lineage, sister-cell or equivalent biological identifier.

## LINK 2 — `CULTURE_LINKED`

A tracked culture / replicate is measured early and later.

## LINK 1 — `TRAJECTORY_ONLY`

Early and late populations are part of the same trajectory but cannot be linked as independent biological units.

## LINK 0 — `NO_FORWARD_LINK`

No defensible early→later mapping.

---

# 5. Outcome-truth classes

## TRUTH 2 — `ORTHOGONAL_OUTCOME`

Future outcome is measured independently of the RNA vector.

This is required for the main biological-fate claim.

## TRUTH 1 — `EXPRESSION_DERIVED_SURROGATE`

Future outcome is computed from future transcriptomics.

This is useful for method development but is not independent biological fate.

## TRUTH 0 — `NO_USABLE_OUTCOME`

No defensible future target.

---

# 6. Resolvability

Cell count is not the effective sample size.

The effective `n` is the number of independent:

```text
donors
cultures
clones / lineages
experiments
```

For a simple two-sided sign-flip comparison:

```text
p_min = 2 / 2^n
```

Examples:

```text
n=3 -> 0.25
n=4 -> 0.125
n=5 -> 0.0625
n=6 -> 0.03125
```

The audit must report whether the proposed outer-unit test is capable of crossing the frozen statistical bar **before model work begins**.

---

# 7. Expected local result — recorded in advance

Expected, not assumed:

```text
GSE242423:
    one HFF trajectory
    no lineage across time
    no independent culture replicates
    LINK 1
    current fate endpoint TRUTH 1

GSE165177:
    bulk
    3 donors
    donor-level n too small for decisive outer-unit inference
    current p_unsafe unsuitable as per-cell truth

GSE165176:
    experimentally meaningful arm/sort metadata may exist
    but a legitimate EARLY -> LATER mapping is not yet established
```

If the audit proves any of this wrong, the audit wins.

---

# 8. Stage 21B — Public prospective-data qualification

This is the main work of Stage 21.

Search for and qualify public datasets with:

```text
X_before
+
U
+
valid biological linkage
+
independently measured Y_after
```

The dataset must be accepted based on experimental geometry, not on whether a model later performs well.

---

# 9. We need TWO complementary data capabilities

A strong paper should not depend on only one fixed-treatment lineage dataset.

Stage 21 therefore searches for two complementary capabilities.

---

## DATA ROLE A — Reprogramming anchor

Purpose:

> Show that the prospective framework applies to the biological problem CellFate-Rx was originally built for.

Preferred geometry:

```text
pre-reprogramming fibroblast state
+
OKSM / reprogramming protocol
->
later independently measured reprogramming success / failure
```

Desired:

```text
lineage or clone linkage
orthogonal future endpoint
multiple independent clones
multiple biological replicates if available
```

Candidate family to audit first:

```text
Rewind-style reprogramming datasets
CellTag / CellTag-multi reprogramming datasets
other lineage-traced reprogramming studies
```

No candidate is accepted from reputation or abstract alone.

---

## DATA ROLE B — Multi-perturbation interaction dataset

Purpose:

> Prove that CellFate-Rx can do something more interesting than predict whether a starting state is generally "good" or whether a treatment is generally "strong."

Desired geometry:

```text
early molecular state X
+
multiple possible perturbations U1, U2, ...
->
later orthogonal response Y
```

Best case:

the same clone / comparable starting population is split across several treatments.

That lets the eventual model answer:

> **Which treatment is best for THIS starting molecular state?**

Candidate families to audit first:

```text
ReSisTrace-like sister-cell treatment-response datasets
multi-treatment clonal tracing datasets
other barcoded drug-response studies
```

---

# 10. Optional DATA ROLE C — Locked replication candidate

If possible, identify a third dataset during Stage 21 and reserve it.

Do not use it to choose model architecture or hyperparameters.

Its purpose is a later locked external replication.

This is strongly preferred but not required to begin Stage 22.

---

# 11. Mandatory qualification table

For every public candidate:

```text
dataset_name
accession
paper

species
cell_type
biological_system

pre_outcome_rna_available
pre_outcome_time

treatments
n_treatments
dose_available
exposure_available

lineage_link_available
clone_link_available
sister_cell_link_available
linkage_method

future_outcome
future_outcome_measurement
future_outcome_orthogonal_to_rna

n_independent_experiments
n_independent_clones
n_positive
n_negative

same_clone_across_multiple_treatments
multi_treatment_interaction_test_possible

processed_files_available
barcode_files_available
outcome_files_available
join_key_reconstructable

outer_split_unit
minimum_attainable_p
resolvable

role_A_reprogramming_anchor
role_B_multi_perturbation
role_C_replication

verdict
```

---

# 12. Hard qualification criteria

A primary dataset must pass:

### Q1 — temporal validity

`X` is measured before `Y`.

### Q2 — biological linkage

The early state is linked to the later outcome through a real biological unit.

### Q3 — orthogonal outcome

The later endpoint is not simply a function of the RNA target.

### Q4 — prediction-time validity

No future label or future-derived field is already present in the predictors.

### Q5 — reconstructability

The public files contain the join keys needed to build the table ourselves.

### Q6 — outcome variation

Useful positive and negative outcomes exist.

### Q7 — independent-unit resolvability

The intended evaluation can, in principle, reach its pre-registered evidence bar.

---

# 13. Additional criterion for the multi-perturbation role

A dataset counts as the main **Role B** dataset only if the perturbation variable contains scientifically meaningful variation.

Preferred:

```text
same clone / sister population exposed to multiple U
```

Acceptable:

```text
multiple treatments with enough shared biological structure to test X + U
```

Not enough:

```text
one treatment + elapsed time only
```

The eventual paper needs evidence that `U` is actually load-bearing.

---

# 14. Qualification verdicts

## `QUALIFIED_ROLE_A`

Valid reprogramming prospective dataset.

## `QUALIFIED_ROLE_B`

Valid multi-treatment prospective dataset.

## `QUALIFIED_BOTH`

One dataset supports both roles.

## `QUALIFIED_SURROGATE`

Excellent forward lineage geometry, but outcome is still expression-derived.

Useful for development, not the independent-fate headline.

## `QUALIFIED_REPLICATION`

Suitable as a locked second/third system.

## `NOT_QUALIFIED`

Fails temporal, linkage, outcome, reconstruction or resolvability requirements.

---

# 15. Stage 21B PASS requirement

The preferred PASS is:

```text
>= 1 QUALIFIED_ROLE_A
AND
>= 1 QUALIFIED_ROLE_B
```

They may be the same dataset if it genuinely supports both.

Why require both?

Because otherwise the final paper risks becoming one of two weak stories:

```text
"we predicted reprogramming success under one fixed treatment"
```

or:

```text
"we predicted generic drug response but lost the rejuvenation/reprogramming connection"
```

The two roles together preserve both the biological motivation and the treatment-conditioned novelty.

---

# 16. Partial qualification branch

If only Role A or only Role B is found:

```text
VERDICT = PARTIAL_DATA
```

Do not start expensive architecture work.

Continue qualification/acquisition for the missing role.

A simple reconstruction may proceed, but Stage 22 does not become the full paper benchmark until the missing role is resolved.

---

# 17. Stage 21C — download and reconstruct

For each dataset needed to satisfy PASS:

1. download the exact public files;
2. record accession, URLs/identifiers and checksums;
3. reconstruct lineage / clone mapping;
4. reconstruct treatment metadata;
5. reconstruct future outcome;
6. verify early-vs-late temporal ordering;
7. build one frozen analysis table;
8. reproduce a central descriptive result from the source study;
9. run leakage guards;
10. freeze the dataset version.

---

# 18. Standard table schema

Each reconstructed dataset should map into:

```text
dataset_id
unit_id
clone_or_lineage_id
replicate_id
donor_or_cell_line

X_time
X_gene_1
...
X_gene_p

treatment_id
dose
exposure_time

Y_time
future_outcome
future_outcome_measurement

outer_group
```

Not every dataset needs every field, but every missing field must be explicit.

---

# 19. Source-study reproduction gate

Before using a reconstructed dataset for our own model, reproduce at least one published descriptive fact.

Examples:

```text
reported clone count
reported outcome class count
reported responder/non-responder enrichment
reported lineage success fraction
reported barcode mapping statistic
```

Purpose:

> prove our reconstruction corresponds to the biological experiment the authors actually performed.

No prospective model is fitted until this passes.

---

# 20. Leakage guards

Hard fail if:

```text
future RNA enters X
future outcome appears in metadata
future cluster appears in predictors
same clone leaks across the outer split
same replicate leaks across the outer split
treatment is reconstructed from the outcome
post-treatment feature selection influences X
```

Synthetic regression tests are mandatory.

---

# 21. Stage 21 final PASS

Stage 21 passes when:

```text
1. Role A is qualified and reconstructed.
2. Role B is qualified and reconstructed.
3. Both use true pre-outcome molecular state.
4. Both have valid biological linkage.
5. Their primary future outcomes are orthogonal to RNA.
6. Both are statistically resolvable for their intended benchmark role.
7. Source-study reconstruction checks pass.
8. Standardized frozen tables exist.
9. Leakage guards pass.
10. src/ is unchanged.
```

Preferred extra:

```text
a third locked replication candidate is identified.
```

---

# 22. Stage 21 output

Example:

```text
STAGE 21 — PROSPECTIVE DATA QUALIFICATION
=========================================

LOCAL DATA
final paper geometry available locally: NO

ROLE A — REPROGRAMMING
dataset: ...
linkage: ...
outcome: ...
orthogonal: YES
resolvable: YES
reconstruction: PASS

ROLE B — MULTI-PERTURBATION
dataset: ...
treatments: ...
same-clone multi-U: YES/NO
outcome: ...
orthogonal: YES
resolvable: YES
reconstruction: PASS

ROLE C — LOCKED REPLICATION
candidate: ...

VERDICT: PASS
NEXT: STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

---

# 23. What Stage 21 proves

A PASS proves only:

> **The project now possesses the datasets required to test a genuinely prospective, treatment-conditioned future-outcome model.**

It does not prove the model works.

---

# 24. What Stage 21 explicitly does NOT do

No:

```text
CellFateNet retraining
deep architecture search
prospective performance claim
RES work
ΔAge work
treatment ranking result
held-out treatment claim
publication metric selection
```

Those begin later.

---

# 25. Next stage

```text
STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

Stage 22 turns the qualified datasets into a standardized benchmark with frozen tasks, splits, metrics and leakage controls.

Only after that benchmark exists do we ask whether the prospective problem is learnable.

---

# 26. One-sentence definition

> **Stage 21 secures the two complementary prospective datasets needed for a non-trivial CellFate-Rx paper: a lineage-grounded reprogramming anchor and a multi-perturbation dataset capable of testing state-dependent treatment response, both with independently measured future outcomes and reconstructable public ground truth.**
