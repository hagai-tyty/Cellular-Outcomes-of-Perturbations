# STAGE 21 — Prospective Fate Ground-Truth & Feasibility Audit on Existing Data

**Status:** PRE-REGISTERED PLAN — freeze in git before execution.  
**Purpose:** Determine, using the datasets already on disk, whether the project can pose a scientifically valid and statistically decision-capable **early-state → later-outcome** task before any prospective CellFate-Rx development begins.  
**Mandatory phase:** Stage 21A — metadata / ground-truth / resolvability audit.  
**Conditional phase:** Stage 21B — cheap forward modelling only if Stage 21A proves the task can change a decision.  
**Fallback:** Stage 21C — prospective-data acquisition handoff.  
**Current model:** `_s16` remains frozen and untouched.  
**Out of scope:** ΔAge optimisation, RES rescue, architecture tuning, retraining `_s16`, target fishing, or manuscript claims.

---

# 0. Why this revision exists

The first Stage-21 draft had the correct evaluation logic but mixed together two different questions:

1. **Can early expression forecast a molecular state observed later?**
2. **Can early expression predict an independently verified biological fate observed later?**

Those are not equivalent.

The project's current fate target, `y_cls`, is generated from expression itself: the pipeline computes somatic-identity / pluripotency / apoptosis-related expression programs and converts them into the safe / identity-loss / death label. Therefore a label at time `t+k` can be a valid **future molecular-state surrogate**, but it is not automatically an independently verified future fate.

This is **not the same literal circularity** as same-timepoint ΔAge:

```text
same-timepoint:
X_t -> f(X_t)
```

is circular / reconstructive.

Whereas:

```text
forward surrogate:
X_t -> f(X_{t+k})
```

is a real forecasting problem because the future transcriptome is not contained in the early transcriptome.

However, the latter still licenses only:

> **early expression forecasts a later expression-derived fate surrogate**

unless the later outcome is anchored by an orthogonal experimental observation.

Stage 21 is revised so that those two levels can never be silently conflated.

---

# 1. Current facts that Stage 21 must respect

These are treated as known constraints to be mechanically verified by Stage 21A, not as assumptions to be rediscovered after modelling.

## 1.1 Current fate labels are expression-derived

The current CellFate-Rx `y_cls` is not an orthogonal biological assay.

It is derived from molecular programs in the same transcriptome being labelled.

Therefore:

```text
current y_cls != independent biological fate ground truth
```

A future `y_cls(t+k)` may be used as a **future molecular-state surrogate**, but must never be described as an independently observed fate.

---

## 1.2 `GSE242423` is one HFF trajectory

Known geometry to verify from the local files:

```text
D0
D2
D4
D6
D8
D10
D12
D14
iPSC
```

The review of the source-file listing found:

- one HFF line;
- one GSM/sample per timepoint;
- no independent culture replicate per day;
- no varying treatment arm — the trajectory is OSKM exposure over time;
- ordinary 10x cell barcodes, not lineage barcodes;
- no demonstrated cross-time cell/clone linkage.

Therefore, unless the audit discovers contrary metadata:

```text
STRICT_LINEAGE      -> expected impossible
CULTURE_REPLICATE   -> expected impossible
TRAJECTORY_FORWARD  -> possible
```

The many single cells **do not create many independent future trajectories**.

---

## 1.3 `GSE165177` is bulk and donor-limited

Known geometry:

- bulk RNA;
- donors O1/O2/O3;
- transient-reprogramming experiment;
- replicated samples and contemporaneous controls;
- only 3 donor-level independent units for donor-held-out inference.

The project already established that a bulk sample cannot provide a per-cell `p_unsafe`.

Stage 21 must not reinterpret bulk sample replication as independent cell-fate observations.

---

## 1.4 Experimental arm annotations may exist and must be audited separately

The project metadata contain experimental marker/arm information such as:

```text
CD13  -> failing-to-reprogram arm
SSEA4 -> reprogramming arm
```

These annotations are **not the same object as `y_cls`**.

Stage 21A must explicitly determine:

- how the arm label was experimentally obtained;
- whether it is orthogonal to the RNA expression vector used as input;
- whether it is available only contemporaneously or can define a later outcome;
- whether an early biological unit can legitimately be linked to that later arm/outcome.

The presence of an experimental label does **not** by itself create a prospective pair.

---

# 2. Stage 21's actual question

The mandatory Stage-21 question is:

> **Do the datasets already held locally contain both (a) a legitimate early→later mapping and (b) enough independent biological units to support a decision-capable prospective test?**

Only if the answer is yes do we ask the modelling question:

> **Does early transcriptomic state improve prediction of the later outcome beyond time, treatment, donor and other allowed metadata?**

---

# 3. Two axes must be audited independently

The original plan used one Level-0-to-Level-3 hierarchy. That is retained for **linkage strength**, but outcome truth is now a separate axis.

This prevents a strong lineage link to a weak surrogate from being mistaken for independent fate ground truth.

---

## 3.1 AXIS A — early→later linkage

### LINK 3 — `LINEAGE_LINKED`

```text
early molecular state
+
known intervention
+
lineage / clone / sister-cell identifier
->
later outcome from the same biological lineage
```

Strongest prospective geometry.

---

### LINK 2 — `CULTURE_LINKED`

```text
early culture / replicate molecular state
+
known intervention
->
later outcome from that same tracked culture / replicate
```

Prospective at culture level.

Does not license a single-cell claim.

---

### LINK 1 — `TRAJECTORY_ONLY`

```text
early population state
->
later population state
```

One continuous trajectory with no direct lineage/replicate link.

Valid for descriptive/pilot forecasting only.

Does not create independent cell-level outcomes.

---

### LINK 0 — `NO_FORWARD_LINK`

Examples:

- cross-sectional only;
- unrelated cells arbitrarily paired across days;
- no later measurement;
- only D0;
- outcome precedes input;
- linkage created solely from matching donor + day without a tracked biological unit.

No forward model may run.

---

## 3.2 AXIS B — outcome ground truth

### TRUTH 2 — `ORTHOGONAL_OUTCOME`

Later outcome comes from information not deterministically computed from the RNA vector being predicted.

Examples:

```text
lineage survival / extinction
colony formation
imaging phenotype
experimentally measured viability
cell-surface sorting / FACS category
functional assay
independent lineage fate
other orthogonal endpoint
```

The audit must document the exact source.

This is the target class required for a future paper claim about **biological fate prediction**.

---

### TRUTH 1 — `EXPRESSION_DERIVED_SURROGATE`

Later outcome is computed from the later transcriptome, for example:

```text
y_cls(t+k) = fate_labels(X_(t+k))
```

This is legitimate for:

> future molecular-state forecasting

but not sufficient for:

> independently validated biological fate prediction.

---

### TRUTH 0 — `NO_USABLE_OUTCOME`

No later target exists or it is contaminated by information unavailable at prediction time.

No forward model may run.

---

# 4. What combinations mean

| Linkage | Truth | Meaning | Can support Paper-1 fate headline? |
|---|---|---|---|
| LINK 3 | TRUTH 2 | lineage-resolved independent future fate | **potentially yes** |
| LINK 2 | TRUTH 2 | culture-level independent future outcome | **potentially yes, culture-level wording** |
| LINK 1 | TRUTH 2 | one trajectory with independent later assay | pilot only unless independent trajectories exist |
| LINK 3 | TRUTH 1 | future molecular-state surrogate on lineage | useful surrogate forecast, **not independent fate** |
| LINK 2 | TRUTH 1 | future culture molecular surrogate | useful surrogate forecast, **not independent fate** |
| LINK 1 | TRUTH 1 | future expression-derived state on one trajectory | descriptive pilot only |
| any | TRUTH 0 | no usable target | stop |
| LINK 0 | any | no valid forward geometry | stop |

---

# 5. Pre-registered hypotheses

Stage 21 distinguishes the **ground-truth hypothesis** from the **molecular-signal hypothesis**.

## H-GT — usable prospective ground truth exists locally

At least one existing dataset provides:

```text
LINK >= 2
AND
TRUTH = 2
AND
statistically resolvable independent units
```

If this fails, the local data are insufficient for the paper-making prospective claim.

---

## H-RNA — incremental early molecular signal

Conditional on a licensed task:

```text
Model(metadata + early RNA) > Model(metadata only)
```

on held-out independent biological units.

This is tested only if the task is resolvable.

---

## H0-RNA

Once allowed metadata are known, adding early RNA gives no reproducible improvement, or any gain is reproduced by the restricted molecular-shuffle null.

---

# 6. Stage 21 remains diagnostic

The following are forbidden:

- modifying `src/`;
- retraining `_s16`;
- changing `fate_labels`;
- redefining the future outcome after seeing model results;
- choosing among several horizons because one works;
- trying several safe/unsafe definitions until one passes;
- architecture tuning;
- embeddings;
- CellFateNet retraining;
- RES changes;
- ΔAge changes;
- random-cell splits as evidence of generalization;
- pseudo-pairing cells across timepoints;
- treating repeated cells from one culture as independent trajectories;
- treating later expression-derived `y_cls` as an independent fate assay;
- lowering the evidence bar after Stage 21A reveals small `n`.

Stage 21 may add only diagnostic / plan / result files.

`src/` must remain byte-unmodified.

---

# 7. Files owned by Stage 21

Mandatory:

```text
plans/STAGE_21_PROSPECTIVE_FATE_EXISTING_DATA.md
experiments/diag_stage21_data_audit.py
tests/test_diag_stage21_data_audit.py
results/diag_stage21_data_audit_results.json
```

Conditional, only if Stage 21A says modelling is decision-capable:

```text
experiments/diag_stage21_forward_fate.py
tests/test_diag_stage21_forward_fate.py
results/diag_stage21_forward_fate_results.json
```

If acquisition is required, Stage 21 opens an acquisition handoff under:

```text
plans/STAGE_21C_PROSPECTIVE_DATA_REQUIREMENT.md
```

**Do not call acquisition "Stage 22".**

This resolves the numbering conflict.

`Stage 22` remains reserved for:

```text
STAGE 22 — Prospective CellFate-Rx
```

and begins only after suitable prospective ground truth has been qualified.

---

# 8. Datasets to audit

## 8.1 `GSE242423` — first audit target

Audit:

```text
n_GSM
n_timepoints
n_cell_lines
n_independent_cultures
n_biological_replicates
lineage_id_available
clone_id_available
cross_time_barcode_link_available
treatment_arms
dose_variation
future_outcome_sources
orthogonal_outcome_available
expression_surrogate_available
```

Known expectation:

```text
one HFF line
nine sequential timepoints
one sample/GSM per timepoint
no lineage
no independent culture replicate
one OSKM trajectory
```

This expectation must be checked mechanically rather than merely copied into the result.

---

## 8.2 `GSE165177` — second audit target

Audit:

```text
n_donors
n_samples
n_replicates
n_timepoints
arm labels
how arm labels were experimentally obtained
whether an early culture is tracked to a later outcome
whether donor is the only independent split unit
whether any orthogonal outcome exists
whether the target is merely fate_labels(expression)
```

Known expectation:

```text
bulk
3 donors
not valid for per-cell p_unsafe
```

---

## 8.3 `GSE165176` — metadata-only prospective audit

The first draft focused on GSE242423 and GSE165177.

This revision adds a **metadata audit only** for GSE165176 because the existing records contain CD13/SSEA4 experimental arm labels.

Question:

> Is there any legitimate early sample / culture whose later experimental arm outcome can be linked forward without already knowing that arm at prediction time?

If no, record that explicitly.

Do not train on it merely because an arm label exists.

---

## 8.4 Excluded from forward modelling

### `GSE165178`

Methylation anchor/reference.

May help verify how experimental arms are named, but is not by itself a forward single-cell fate dataset.

### `GSE165179`

Methylation reference/twin.

Not a forward fate dataset.

### `GSE113957`

Clock-training cohort.

No reprogramming future-fate task.

### `GSE297234`

D0-only.

No later outcome.

---

# 9. STAGE 21A — mandatory data / ground-truth audit

No model fitting is permitted before Stage 21A produces a verdict.

---

## 9.1 Required audit table

Each candidate dataset must produce:

```text
dataset
modality
species
cell_line_count
donor_count
biological_replicate_count
culture_count
n_observations
n_timepoints
n_treatment_arms
dose_available

lineage_id_available
clone_id_available
cross_time_link_available
culture_followup_available

current_y_cls_source
future_target_candidate
future_target_source
future_target_is_expression_derived
future_target_is_orthogonal

linkage_level
truth_level

split_unit
n_independent_outer_units
class_balance_by_outer_unit
mixed_outcome_outer_units

minimum_attainable_p
pass_statistically_resolvable

leakage_risk
licensed_claim
verdict
```

---

# 10. The missing piece in the first plan: RESOLVABILITY

Stage 21A must determine whether Stage 21B can possibly change a scientific decision **before any model is fitted**.

Cell count is not the effective `n`.

The effective `n` is the number of independent outer biological units relevant to the claim:

```text
donors
independent cultures
independent clones / lineages
independent experimental trajectories
```

Repeated cells or timepoints within one trajectory do not automatically increase this `n`.

---

## 10.1 Exact minimum-p check

For the primary paired outer-unit comparison:

```text
Δ_i = loss(M+X)_i - loss(M)_i
```

Stage 21A must report the smallest attainable two-sided unit-level sign-flip/permutation p-value if every independent unit favors `M+X`.

For a simple sign-flip test:

```text
p_min = 2 / 2^n
```

subject to the exact permutation design actually used.

Examples:

```text
n = 1 -> p_min = 1.000
n = 2 -> p_min = 0.500
n = 3 -> p_min = 0.250
n = 4 -> p_min = 0.125
n = 5 -> p_min = 0.0625
n = 6 -> p_min = 0.03125
```

Therefore an `n=3` donor comparison cannot produce a conventional two-sided `p < 0.05` even under perfect directional agreement.

This must be known **before** Stage 21B.

---

## 10.2 Resolvability verdict

### `RESOLVABLE_FOR_PASS`

The planned independent-unit test can, in principle, cross the frozen statistical bar and the geometry contains outcome variation.

Stage 21B may run.

---

### `UNRESOLVABLE_FOR_PASS`

Even a perfect directional result cannot cross the frozen bar, or the independent units do not contain enough outcome variation.

Default action:

```text
SKIP STAGE 21B
```

Do not spend a week generating an estimate that cannot alter the project's decision.

A descriptive pilot may be run only under a separately recorded decision, and it can never be upgraded to PASS.

---

# 11. Stage 21A verdicts

Stage 21A returns exactly one top-level verdict.

## `LOCAL_ORTHOGONAL_TASK_READY`

Requirements:

```text
LINK >= 2
TRUTH = 2
RESOLVABLE_FOR_PASS
```

Action:

```text
Run Stage 21B.
```

A strong Stage-21B result may license Stage 22.

---

## `LOCAL_SURROGATE_TASK_READY`

Requirements:

```text
LINK >= 2
TRUTH = 1
RESOLVABLE_FOR_PASS
```

Action:

Stage 21B may run if the purpose is to determine whether early transcriptomics carries forward molecular-state information.

Even a strong result does **not** establish independently validated fate.

It does not by itself license the Paper-1 biological-fate headline.

---

## `TRAJECTORY_PILOT_ONLY`

Typical geometry:

```text
LINK = 1
```

or only one continuous biological trajectory.

Action:

Default:

```text
do not run expensive Stage 21B modelling
```

Record the geometry and proceed to Stage 21C acquisition unless a cheap descriptive forecast is explicitly judged useful for acquisition prioritization.

---

## `LOCAL_TASK_UNRESOLVABLE`

A possible forward task exists but:

```text
p_min >= 0.05
```

or equivalent resolvability failure.

Action:

```text
SKIP Stage 21B
Open Stage 21C
```

---

## `NO_VALID_FORWARD_TASK`

No valid early→later mapping or no usable later target.

Action:

```text
SKIP Stage 21B
Open Stage 21C
```

---

# 12. Expected local geometry is NOT a verdict

Based on the current review, Stage 21A is expected to find something close to:

```text
GSE242423:
    LINK 1
    likely TRUTH 1
    one biological trajectory
    -> TRAJECTORY_PILOT_ONLY / UNRESOLVABLE

GSE165177:
    bulk
    n_donor = 3
    current fate proxy expression-derived
    -> likely UNRESOLVABLE_FOR_PASS

GSE165176:
    experimental CD13/SSEA4 arm metadata exists
    prospective linkage not yet established
    -> audit before deciding
```

These are **predictions recorded before the audit**, not substituted for the audit.

If the files prove them wrong, the audit result wins.

---

# 13. STAGE 21B — conditional cheap forward test

Stage 21B runs only after:

```text
pass_statistically_resolvable = true
```

and a single target is frozen.

No CellFate-Rx retraining occurs.

---

# 14. Freeze one task before fitting

Required fields:

```text
dataset
linkage_level
truth_level

prediction_unit
input_time
outcome_time
prediction_horizon

future_outcome_definition
future_outcome_source
is_outcome_expression_derived

allowed_metadata
gene_panel
split_unit

primary_metric
outer_test
permutation_scheme
```

Once written, none may change because of the observed result.

---

# 15. Target naming rule

This is mandatory.

If:

```text
TRUTH = 1
```

the target must be named in code/results/manuscript notes as something like:

```text
future_expression_fate_surrogate
future_molecular_fate_state
future_signature_state
```

Do **not** name it:

```text
true_fate
future_biological_fate
survival
death
identity_preservation
```

unless an orthogonal outcome directly supports that term.

---

# 16. Preferred outcome hierarchy

## 16.1 First choice

Orthogonally measured later outcome, if available.

Examples:

```text
future experimental reprogramming outcome
later sorted phenotype
later viability/survival outcome
later clone fate
```

---

## 16.2 Second choice

Future expression-derived binary surrogate:

```text
SAFE_SURROGATE
UNSAFE_SURROGATE
```

where the label is explicitly documented as a transformation of `X_(t+k)`.

This is a molecular forecast, not biological fate validation.

---

## 16.3 Three-class secondary analysis

Only if independent-unit class counts make it resolvable.

Never replace a failed binary primary with a three-class result after seeing the outcome.

---

# 17. Models

Stage 21B uses only deliberately simple models.

## M — metadata only

Potential features:

```text
early time
treatment / arm information available at prediction time
dose
donor age
pre-treatment experimental covariates
```

Print the complete feature list.

---

## X — early transcriptome only

Use a frozen gene panel.

Prefer the current CellFate-Rx 2,000-gene panel if coverage permits.

No outcome-aware gene selection.

---

## M+X — metadata + early transcriptome

Primary scientific comparison:

```text
M+X versus M
```

Secondary:

```text
M+X versus X
```

when treatment/metadata variation is sufficient to make that comparison meaningful.

---

# 18. Model family

Binary:

```text
L2-regularized logistic regression
```

Three-class secondary:

```text
multinomial L2 logistic regression
```

No neural network.

Stage 21 asks whether the prospective information exists, not whether the architecture is clever enough to extract it.

---

# 19. Hyperparameter grid

If needed:

```text
C = [0.001, 0.01, 0.1, 1.0, 10.0]
```

Select inside grouped training-only inner CV.

The outer unit must not affect:

- scaling;
- imputation;
- feature selection;
- regularization;
- calibration;
- thresholding.

---

# 20. Splitting

Priority:

```text
held-out donor
held-out independent culture
held-out lineage / clone
held-out independent trajectory
```

Random-cell split is forbidden as primary evidence.

All observations from an independent unit stay entirely on one side.

A single HFF trajectory cannot be converted into donor-level replication by splitting its cells.

---

# 21. Primary metric

Use:

```text
negative log likelihood / log loss
```

Primary quantity:

```text
ΔNLL = NLL(M+X) - NLL(M)
```

Interpretation:

```text
negative -> RNA adds value
zero     -> no incremental value
positive -> RNA hurts
```

Absolute AUROC is not the primary decision statistic.

---

# 22. Secondary metrics

Where defined:

```text
AUROC
PR-AUC
balanced accuracy
Brier score
ECE
```

A secondary metric cannot rescue a failed primary result.

---

# 23. Mandatory molecular-shuffle null

Shuffle `X` while preserving metadata structure as tightly as geometry permits.

Preferred:

```text
shuffle X within treatment × early-time × experiment strata
```

Do not shuffle:

- outcome;
- treatment;
- time;
- donor;
- outer split.

Use at least:

```text
1000 permutations
```

for a final Stage-21B result.

---

# 24. Leakage guards

## 24.1 Same-state leakage

If the input and target are deterministic functions of the same molecular observation:

```text
HARD FAIL
```

---

## 24.2 Temporal leakage

No data measured after the prediction cutoff may enter `X` or metadata.

---

## 24.3 Future-arm leakage

If the later outcome is an experimental arm/category, confirm that this category was **not already known for the early input unit** in a way that makes prediction trivial.

For example, if an early sample is already SSEA4-sorted and the later target is simply "SSEA4 arm", this is not a useful future prediction task.

---

## 24.4 Unit leakage

No donor / culture / clone / lineage defining the outer split may cross train/test.

---

## 24.5 Expression-derived target transparency

If the target is generated by `fate_labels(X_future)`, the result JSON must set:

```text
"outcome_is_expression_derived": true
```

and the console must print:

```text
GROUND TRUTH: EXPRESSION-DERIVED SURROGATE
```

---

# 25. Tests required before execution

## Data-audit tests

Construct synthetic metadata proving:

```text
lineage + orthogonal future assay
    -> LINK3 / TRUTH2

culture follow-up + orthogonal outcome
    -> LINK2 / TRUTH2

future transcriptomic signature only
    -> TRUTH1

one trajectory, many cells
    -> LINK1, independent_n = 1

three donors
    -> p_min = 0.25 under two-sided sign-flip

five donors
    -> p_min = 0.0625

six donors
    -> p_min = 0.03125

cross-sectional data
    -> LINK0

D0 only
    -> LINK0
```

---

## Forward-test tests

Only needed if Stage 21B is licensed.

Must cover:

```text
metadata-only signal
    -> M+X gets no false credit

real early-RNA future signal
    -> M+X improves

random RNA
    -> no systematic gain

restricted shuffled RNA
    -> gain disappears

future feature in X
    -> hard failure

future arm already known at t0
    -> hard failure / task invalid

outer-unit leakage
    -> hard failure

training-fold-only scaling
    -> verified
```

Every verdict branch must be exercised.

---

# 26. Stage 21B verdicts

## `NO_INCREMENTAL_RNA`

Valid resolvable task exists, but:

```text
M+X does not beat M
```

on the primary metric.

Action:

Do not open prospective model development on this dataset.

---

## `SURROGATE_FORWARD_SIGNAL`

Requirements:

```text
TRUTH = 1
resolvable task
M+X beats M
restricted shuffle null rejects
effect not carried by one unit
```

Meaning:

> Early RNA contains information about a later expression-derived molecular fate surrogate beyond metadata.

Action:

Useful as a feasibility lead.

Still open Stage 21C for independent future biological outcome data.

Do **not** call this validated fate prediction.

---

## `ORTHOGONAL_FORWARD_SIGNAL`

Requirements:

```text
TRUTH = 2
LINK >= 2
resolvable task
M+X beats M
restricted shuffle null rejects
effect not carried by one unit
no leakage guard fires
```

Meaning:

> Early RNA contains information about an independently measured later outcome beyond metadata.

Action:

Stage 21 passes the prospective-ground-truth gate.

Stage 22 may open.

---

# 27. No automatic upgrade from surrogate to fate

This sentence should appear in the final Stage-21 result if `TRUTH = 1`:

> **The target is computed from the later transcriptome. This is a genuine forward molecular-state forecast, not an independently verified fate outcome.**

This limitation follows the result into:

```text
CHANGES.md
ARCHITECTURE.md
Stage 25 claim lock
manuscript notes
```

It cannot disappear later because the model score is strong.

---

# 28. STAGE 21C — acquisition handoff

Stage 21C opens when:

```text
NO_VALID_FORWARD_TASK
LOCAL_TASK_UNRESOLVABLE
TRAJECTORY_PILOT_ONLY
LOCAL_SURROGATE_TASK_READY but Paper-1 requires orthogonal fate
```

Stage 21C is **still part of Stage 21**.

This keeps the roadmap numbering unambiguous.

---

# 29. What Stage 21C must ask for

The minimum paper-making dataset is:

```text
early molecular state X
+
proposed / known treatment U
+
later orthogonal outcome Y
+
valid early→later biological linkage
+
enough independent outer units for the frozen test to be resolvable
```

Preferred additions:

```text
multiple donors
multiple treatments or doses
mixed future outcomes within comparable treatment/time strata
lineage / clone information
```

The acquisition requirement is defined by the failed audit field, not by vague "more data."

---

# 30. What counts as orthogonal enough for the paper goal

A later endpoint does not need to be perfect biological truth.

It does need to be **measured independently of the RNA vector used to define the label**.

Examples that can qualify if properly linked:

```text
survival/death from lineage abundance
viability assay
colony formation
imaging
FACS/surface phenotype
experimentally observed reprogramming success/failure
lineage-resolved terminal state
```

The detailed Stage-21C plan must freeze the chosen endpoint before model fitting.

---

# 31. Required Stage-21A result JSON

```json
{
  "stage": 21,
  "phase": "A_data_ground_truth_audit",
  "datasets": [],
  "known_constraints_verified": {},
  "selected_dataset": null,
  "linkage_level": null,
  "truth_level": null,
  "split_unit": null,
  "n_independent_outer_units": null,
  "minimum_attainable_p": null,
  "pass_statistically_resolvable": false,
  "licensed_claim": null,
  "stage21b_licensed": false,
  "stage21c_required": false,
  "verdict": null,
  "src_modified": false
}
```

---

# 32. Conditional Stage-21B result JSON

```json
{
  "stage": 21,
  "phase": "B_forward_test",
  "dataset": "",
  "linkage_level": "",
  "truth_level": "",
  "outcome_is_expression_derived": false,
  "biological_unit": "",
  "split_unit": "",
  "n_independent_outer_units": 0,
  "models": {
    "M": {},
    "X": {},
    "MX": {}
  },
  "primary": {
    "metric": "log_loss",
    "delta_nll_mx_vs_m": null,
    "ci95": [null, null]
  },
  "outer_unit_test": {
    "method": "",
    "p_value": null,
    "minimum_attainable_p": null
  },
  "shuffle_null": {
    "n_permutations": 1000,
    "p_value": null
  },
  "leakage_checks": {},
  "licensed_claim": "",
  "verdict": "",
  "src_modified": false
}
```

---

# 33. Console output — neutral, no pre-selected primary dataset

The first draft accidentally showed GSE165177 as the example primary.

Remove that anchor.

Stage 21A output should look like:

```text
STAGE 21A — PROSPECTIVE GROUND-TRUTH AUDIT
==========================================

dataset       link   truth   independent_n   p_min   resolvable   verdict
-----------   ----   -----   -------------   -----   ----------   -------
<dataset>     ...    ...     ...             ...     ...          ...
<dataset>     ...    ...     ...             ...     ...          ...

ORTHOGONAL OUTCOME AVAILABLE: ...
DECISION-CAPABLE LOCAL TASK:   ...
STAGE 21B LICENSED:            ...
STAGE 21C REQUIRED:            ...

VERDICT: ...
```

Only if Stage 21B is licensed should a model-comparison block print.

---

# 34. Commands

Mandatory audit:

```bash
python experiments/diag_stage21_data_audit.py ^
  "D:\GSE242423" ^
  "D:\Gill" ^
  "D:\GSE165177" ^
  "D:\GSE165178" ^
  "D:\GSE165179" ^
  "D:\GSE113957" ^
  "D:\GSE297234"
```

Then:

```bash
pytest -q
ruff check .
git diff --stat src/
```

Required:

```text
src/ unchanged
```

---

## Stage 21B

Run only if:

```text
stage21b_licensed == true
```

Command:

```bash
python experiments/diag_stage21_forward_fate.py ^
  --audit results/diag_stage21_data_audit_results.json
```

If the audit says Stage 21B is unresolvable, the script must refuse to run unless an explicit descriptive-only override is supplied.

The normal scientific workflow does **not** use that override.

---

# 35. PASS/FAIL rules must not change after audit

The review correctly endorsed the first plan's core evaluation rule.

Keep it.

Example:

```text
time-only AUC = 0.98
M+X AUC       = 0.98
```

is a failure of incremental RNA value.

Likewise, a huge apparent AUC on an expression-derived surrogate does not turn it into independent fate ground truth.

There are now two separate gates:

```text
GROUND-TRUTH GATE
Is Y an independently measured future outcome?

SIGNAL GATE
Does X improve prediction of Y beyond metadata?
```

Both are required for the full prospective-fate claim.

---

# 36. What Stage 21 can prove

Depending on verdict, Stage 21 may establish one of three things.

### A. Data insufficiency

> The existing datasets cannot pose a decision-capable prospective biological-fate problem.

This is a useful acquisition result.

---

### B. Future molecular-state predictability

> Early transcriptomics forecasts a later expression-derived fate surrogate beyond metadata.

This is genuine forecasting.

It is not independently validated biological fate.

---

### C. Prospective orthogonal-outcome predictability

> Early transcriptomics improves prediction of an independently measured later outcome beyond metadata.

This is the result that can open Stage 22.

---

# 37. What Stage 21 cannot prove

Even the strongest local PASS does not automatically establish:

- arbitrary perturbation generalization;
- held-out-treatment generalization;
- clinical use;
- safe rejuvenation;
- accurate future ΔAge;
- working RES;
- single-cell fate if the unit is a culture;
- donor generalization if donors are not held out;
- that deep learning is necessary;
- that CellFate-Rx beats the strongest simple baseline.

Those belong later.

---

# 38. What happens after each outcome

## `LOCAL_ORTHOGONAL_TASK_READY` + `ORTHOGONAL_FORWARD_SIGNAL`

```text
Stage 21 PASS
        ↓
STAGE 22 — Prospective CellFate-Rx
```

Stage 22 may now adapt the fate head to the independently measured future target.

---

## `LOCAL_SURROGATE_TASK_READY` + `SURROGATE_FORWARD_SIGNAL`

```text
useful feasibility lead
        ↓
Stage 21C — acquire/qualify orthogonal prospective outcome data
```

Do not open the paper claim yet.

---

## `TRAJECTORY_PILOT_ONLY`

```text
no decision-capable local replication
        ↓
Stage 21C
```

Do not mistake thousands of cells for thousands of independent trajectories.

---

## `LOCAL_TASK_UNRESOLVABLE`

```text
skip 21B
        ↓
Stage 21C
```

This is specifically designed to avoid spending a week re-deriving a decision already implied by `n`.

---

## `NO_VALID_FORWARD_TASK`

```text
skip 21B
        ↓
Stage 21C
```

---

# 39. Why Stage 21A is still worth running

Even if the expected answer is "acquisition required", Stage 21A remains valuable because it turns a verbal concern into a reproducible artifact.

It will record exactly:

```text
what target exists
how that target was measured
whether it is orthogonal
what early→late linkage exists
what the independent unit is
how many independent units exist
whether a perfect result could meet the statistical bar
what dataset property blocks the claim
```

That information directly specifies what new data must contain.

The expensive modelling phase is no longer automatic.

---

# 40. Why Stage 21B is now conditional

The first plan assumed the cheap logistic-regression test was worth running whenever a forward task could be written.

The review correctly identifies a stronger rule:

> **Do not run a statistical experiment whose geometry cannot change the decision.**

Therefore Stage 21B runs only after the resolvability gate.

If `GSE165177` gives `n=3` independent donors, a conventional donor-level two-sided PASS is impossible under the exact sign-flip logic even if all three donors improve.

If `GSE242423` gives one continuous HFF trajectory, its thousands of cells do not repair the missing independent trajectory count.

Those are data facts, not model failures.

---

# 41. Relationship to the old ΔAge circularity result

The distinction must remain explicit.

## Same-state ΔAge

```text
X_t -> clock(X_t)
```

The target is directly reconstructable from the input.

Circular for prediction claims.

---

## Stage-21 expression-derived forward surrogate

```text
X_t -> fate_labels(X_(t+k))
```

The future target is not present in `X_t`.

This is a valid temporal forecast.

But both input and target live in the transcriptomic modality and the target is still a designed molecular signature.

Therefore it supports **forecasting of a future molecular state**, not independent fate validation.

---

## Desired Paper-1 target

```text
X_t + U -> Y_(t+k)
```

where:

```text
Y_(t+k)
```

comes from a future experimental observation independent of the RNA-derived label.

That is the target Stage 21C must secure if local data cannot provide it.

---

# 42. Claim firewall for later stages

The following exact distinction must survive into Stage 25.

If target is expression-derived:

> **The model forecasts a later transcriptomic fate surrogate.**

If target is orthogonal:

> **The model predicts an independently measured later cellular outcome.**

Only use stronger terms such as:

```text
survival
death
identity preservation
reprogramming success
```

when the target measurement directly warrants them.

---

# 43. Freeze rule

Before execution:

1. commit this Stage-21 plan;
2. commit `test_diag_stage21_data_audit.py`;
3. verify all expected geometry checks are tests, not prose only;
4. verify `src/` clean;
5. run Stage 21A once;
6. record the audit result;
7. compute resolvability before any model fitting;
8. run Stage 21B only if licensed;
9. append the result to `CHANGES.md`;
10. do not rewrite bars or target definitions after seeing results.

---

# 44. One-sentence Stage 21 definition

> **Stage 21 audits whether the data already on disk contain a legitimate, independently grounded and statistically resolvable early-state → later-outcome task; only if that gate passes may a cheap forward molecular-signal test run, otherwise the stage ends by specifying exactly what prospective data must be acquired.**
