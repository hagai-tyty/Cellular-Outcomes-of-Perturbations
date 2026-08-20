# STAGE 21 — Prospective Fate Feasibility on Existing Data

**Status:** PRE-REGISTERED PLAN — do not interpret any result until this document is frozen in git.  
**Primary question:** Can the datasets already held on `D:\` support a legitimate **prospective** fate task, and if so, does an early transcriptomic state contain information about a later fate beyond elapsed time / treatment metadata alone?  
**Primary datasets:** `GSE242423`, then `GSE165177`.  
**Out of scope:** new model architecture, RES, ΔAge optimization, new wet-lab experiments, new public datasets unless this stage reaches the acquisition branch.  
**Current model:** `_s16` remains frozen and untouched.

---

## 0. Why Stage 21 exists

CellFate-Rx was designed to answer a prospective intervention question:

> Given a current molecular state `X` and a proposed perturbation `U`, predict a later cellular outcome.

The current project does **not** yet demonstrate that claim.

What is currently established is narrower:

- fate classification contains real expression-level information **within a timepoint**;
- the clean stratified fate result is AUC `0.917`, permutation `p = 0.0091`;
- that result rests on only **12 safe-vs-unsafe pairs**, from **7 mixed strata across 3 donors**;
- most marginal fate performance is explained by timepoint;
- same-timepoint ΔAge prediction is circular and must not be used as evidence of future prediction;
- the current Gill held-out arm is bulk RNA, so it cannot supply a true per-cell unsafe fraction;
- the current project therefore has a trustworthy instrument but has not yet shown a true `X + U -> future Y` capability.

Before downloading new lineage-tracing datasets or changing CellFate-Rx, Stage 21 asks whether the data already present on disk can support a legitimate forward-looking pilot.

The intended decision is simple:

> **Can we construct a valid future-outcome task from existing data, and does early RNA add signal beyond metadata?**

If no, we stop and acquire the correct data.

If yes, we have earned the right to open the next prospective-model stage.

---

# 1. Scientific question

## 1.1 Primary question

> **Does an early transcriptomic state contain information about a later cellular fate beyond the information already present in elapsed time and other allowed metadata?**

The emphasis is on **later**.

A label derived from the same RNA vector used as input is **not** a prospective target.

A cell measured after the outcome has occurred is **not** a prospective input.

A comparison of unrelated cells at two timepoints is **not** automatically a lineage prediction.

---

## 1.2 Strongest desired form

The strongest valid row would be:

```text
X_early + U
      ↓
independently observed future fate of the same lineage / clone / tracked unit
```

where:

- `X_early` is measured before the outcome;
- `U` is the perturbation, treatment, dose, or exposure state;
- the future fate is measured later;
- the early and later observations are linked by a valid biological unit.

Examples of valid linkage, from strongest to weakest:

1. lineage barcode / clone identity;
2. sister-cell / carbon-copy linkage;
3. matched culture or biological replicate that is explicitly followed forward;
4. population-level trajectory where early state is used to predict a later population outcome.

Stage 21 must report which level is actually supported.

It must **not silently upgrade** a population-level task into a cell-level claim.

---

# 2. Pre-registered hypotheses

## H1 — prospective molecular signal

After accounting for allowed metadata, early transcriptomic state improves prediction of a later fate.

Formally:

```text
Model(metadata + early RNA) > Model(metadata only)
```

on held-out biological units.

---

## H0 — metadata explains the apparent signal

Once time / treatment / donor or other allowed metadata are known, adding early RNA provides no reproducible benefit.

Formally:

```text
Model(metadata + early RNA) ≈ Model(metadata only)
```

or any apparent gain disappears under the molecular-shuffle control.

---

# 3. Stage 21 is a diagnostic stage, not a model-development stage

The following are **forbidden during Stage 21**:

- modifying `src/`;
- retraining `_s16`;
- changing the fate-label definition after seeing results;
- tuning CellFateNet architecture;
- adding new hidden layers / embeddings / losses;
- changing RES;
- changing ΔAge;
- tuning thresholds to improve a verdict;
- using the held-out unit to choose preprocessing or hyperparameters;
- constructing pseudo-lineages by matching cells only because they share a donor and timecourse;
- trying multiple target definitions until one passes.

Stage 21 may add only:

```text
plans/STAGE_21_PROSPECTIVE_FATE_EXISTING_DATA.md
experiments/diag_stage21_data_audit.py
experiments/diag_stage21_forward_fate.py        # only if the audit licenses it
tests/test_diag_stage21_data_audit.py
tests/test_diag_stage21_forward_fate.py         # only if forward test is licensed
results/diag_stage21_data_audit_results.json
results/diag_stage21_forward_fate_results.json  # only if forward test is licensed
```

`src/` must remain byte-unmodified.

---

# 4. Datasets to audit

Stage 21 starts with the datasets already present locally.

## 4.1 `GSE242423` — PRIMARY

Known project role:

- single-cell RNA-seq;
- HFF reprogramming trajectory;
- approximately 42.5k cells;
- 9 timepoints;
- currently the only load-bearing single-cell dataset;
- supplies the large majority of training mass.

Why it is the first candidate:

- it has single-cell resolution;
- it spans a timecourse;
- it already contains the strongest existing fate information;
- no download is needed.

Unknowns Stage 21 must inspect explicitly:

```text
clone_id?
lineage_barcode?
cell lineage across time?
replicate_id?
culture_id?
treatment arms?
dose variation?
explicit future outcome?
sampling structure?
whether cells at different timepoints are biologically linkable?
```

Do not assume any of these exist because the dataset is a timecourse.

---

## 4.2 `GSE165177` — SECONDARY

Known project role:

- bulk RNA;
- donors O1/O2/O3;
- transient reprogramming arm;
- multiple timepoints;
- replicated contemporaneous controls;
- useful for ΔAge analysis;
- not suitable for per-cell `p_unsafe`.

Why it is still worth auditing:

It may support a weaker but valid task such as:

```text
early bulk culture state -> later culture / arm outcome
```

It cannot support:

```text
single-cell RNA -> later fate of that same cell
```

unless the source metadata contain a linkage not previously used.

Stage 21 must therefore label any usable task from `GSE165177` as **population/culture-level**, not cell-level.

---

## 4.3 Datasets excluded from prospective-fate modeling in this stage

### `GSE165178`

Methylation reference / paired anchor.  
Not a prospective single-cell fate dataset.

### `GSE165179`

Methylation twin / reference instrument.  
Not a prospective single-cell fate dataset.

### `GSE113957`

Clock-training cohort.  
Not a reprogramming prospective-fate dataset.

### `GSE297234`

Two donors, D0 only.  
No timecourse, therefore no later fate task.

These may remain available for other project analyses but must not be pulled into Stage 21 to enlarge the sample count artificially.

---

# 5. Stage 21A — Dataset-geometry audit

No model may be fitted until Stage 21A finishes.

## 5.1 Required audit table

For each candidate dataset, output one row containing at least:

```text
dataset
modality
species
cell_line_or_donor_count
biological_replicate_count
n_observations
n_timepoints
treatment_arms
dose_available
time_available
lineage_id_available
clone_id_available
replicate_id_available
culture_id_available
early_rna_available
later_outcome_available
early_to_late_link_type
future_outcome_definition
n_independent_forward_units
n_safe
n_unsafe
n_mixed_strata
strict_cell_level_possible
culture_level_possible
population_level_possible
leakage_risk
verdict
```

---

## 5.2 Valid task levels

The audit must classify each dataset into exactly one strongest supported level.

### LEVEL 3 — `STRICT_LINEAGE`

Valid structure:

```text
early transcriptome
+
known intervention
+
lineage / clone / sister-cell link
->
later independently observed fate
```

This can support a genuine cell/lineage prospective claim.

---

### LEVEL 2 — `CULTURE_FORWARD`

Valid structure:

```text
early culture / replicate transcriptome
+
known intervention
->
later independently observed culture outcome
```

This can support a prospective **population/culture** claim, not a single-cell claim.

---

### LEVEL 1 — `TRAJECTORY_FORWARD`

Valid structure:

```text
early population transcriptome / early cell distribution
->
later population fate distribution
```

No direct lineage mapping.

This is acceptable only as a **pilot feasibility** experiment.

It cannot be the final prospective-paper result by itself.

---

### LEVEL 0 — `INVALID_PROSPECTIVE`

Examples:

- current-state label derived from the same expression vector;
- cells at D4 paired arbitrarily with cells at D10;
- outcome already encoded in the input;
- no valid early-to-late mapping;
- no later outcome;
- no fate variation;
- only D0 data.

No prospective model may be fitted on that dataset.

---

# 6. Stage 21A verdict

The audit returns one of:

```text
STRICT_LINEAGE_AVAILABLE
CULTURE_FORWARD_AVAILABLE
TRAJECTORY_FORWARD_ONLY
NO_VALID_FORWARD_TASK
```

Decision:

```text
STRICT_LINEAGE_AVAILABLE
    -> Stage 21B may run at lineage level.

CULTURE_FORWARD_AVAILABLE
    -> Stage 21B may run at culture level.
       Manuscript language must remain population-level.

TRAJECTORY_FORWARD_ONLY
    -> Stage 21B may run as a pilot only.
       A PASS cannot by itself license the final prospective paper claim.

NO_VALID_FORWARD_TASK
    -> STOP Stage 21.
       Do not modify CellFate-Rx.
       Open a data-acquisition stage for lineage-resolved / second-timecourse data.
```

---

# 7. Stage 21B — Define one forward task

Only one primary task may be selected.

The task must be chosen from the strongest valid structure discovered in Stage 21A.

The exact definition must be written to the result file **before any model is fitted**.

Required frozen fields:

```text
input_time
outcome_time
prediction_horizon
biological_unit
input_features
allowed_metadata
future_outcome_definition
safe_definition
unsafe_definition
split_unit
primary_metric
secondary_metrics
permutation_scheme
```

---

# 8. Preferred target hierarchy

Use the strongest target the data genuinely support.

## First choice — binary future fate

```text
SAFE
UNSAFE
```

where `UNSAFE` may combine:

```text
identity loss OR death
```

only if those labels are independently defined at the later outcome.

Reason:

Stage 21 asks whether any prospective molecular signal exists.

A three-class task wastes power if one or both unsafe subclasses are rare.

---

## Secondary only — three-class fate

```text
identity preserved
identity loss
death
```

Run only if every class has enough independent outcome units.

This cannot replace the primary binary test after seeing the binary result.

---

# 9. GSE242423 fallback pilot if no lineage link exists

If GSE242423 has no lineage/clone mapping but supports a legitimate population trajectory, Stage 21 may run a **trajectory-forward pilot**.

The question becomes:

> Does the molecular state of the earlier HFF population improve prediction of a later fate distribution beyond knowing the earlier timepoint alone?

Example geometry:

```text
D2 molecular state  -> later D4/D7 fate distribution
D4 molecular state  -> later D7/D10 fate distribution
D7 molecular state  -> later D10/D12 fate distribution
...
```

Rules:

1. later cells define the target;
2. no later expression may enter `X_early`;
3. no individual early cell may be described as having a known later fate;
4. the unit of inference is the timepoint / replicate / culture structure supported by metadata;
5. overlapping windows must not be treated as independent if they share the same underlying culture;
6. uncertainty must reflect the number of independent temporal units, not the number of single cells.

This test is explicitly **pilot evidence only**.

---

# 10. GSE165177 fallback pilot

If metadata support a valid early-to-later culture link, the question may be:

> Does an earlier bulk transcriptomic state predict a later culture/arm outcome beyond donor + day + arm metadata?

Do not use `p_unsafe` as the target.

The project already established that a bulk sample cannot express a per-cell unsafe fraction.

Any GSE165177 result must be described as:

```text
culture-level
population-level
sample-level
```

never as single-cell future fate.

---

# 11. Models — intentionally simple

Stage 21 is not testing whether CellFate-Rx beats every baseline.

It is testing whether **prospective molecular information exists at all**.

Use three models.

## M — metadata only

Allowed examples:

```text
early time
treatment identity
dose
donor age
known donor / experiment covariates
```

Do not include anything measured after the prediction time.

---

## X — early transcriptome only

Use the frozen Stage-21 gene panel.

If the existing 2,000-gene CellFate-Rx panel is present in the candidate dataset, use it.

If coverage is incomplete:

- report exact coverage;
- do not pick a new gene set based on outcome association;
- any alternate fixed gene set must be specified before fitting.

---

## M+X — metadata + early transcriptome

This is the primary model.

The main scientific comparison is:

```text
M+X versus M
```

The secondary interaction-usefulness comparison is:

```text
M+X versus X
```

---

# 12. Model family

Primary family:

```text
L2-regularized logistic regression
```

for binary fate.

For three-class secondary analysis:

```text
multinomial L2 logistic regression
```

No neural model in Stage 21.

Reason:

If simple regularized prediction cannot establish incremental prospective signal, architecture development is not justified.

---

# 13. Preprocessing

All preprocessing must be fitted on the training partition only.

Forbidden:

```text
global standardization before split
global feature filtering using all samples
global PCA before split
outcome-aware gene selection
test-fold normalization statistics
```

Allowed:

```text
training-fold standardization
training-fold imputation
fixed pre-existing gene panel
training-only regularization selection
```

If PCA is required because geometry makes raw logistic regression unstable, it must be an explicitly pre-registered secondary sensitivity and fitted inside each training fold.

It cannot replace a failed primary analysis.

---

# 14. Hyperparameter selection

If regularization is tuned, use only this grid:

```text
C = [0.001, 0.01, 0.1, 1.0, 10.0]
```

Select by grouped inner cross-validation on training units only.

The held-out outer unit must never influence:

- `C`;
- scaling;
- imputation;
- feature selection;
- thresholding;
- calibration.

---

# 15. Splitting

Random-cell split is forbidden as the primary evaluation.

Use the strongest independent unit available.

Priority:

```text
1. held-out donor
2. held-out biological replicate / culture
3. held-out clone / lineage
4. held-out temporal block
```

All observations belonging to one biological unit must stay on one side of a split.

If Stage 21A shows only one HFF culture without true independent replicate structure, the result must say so explicitly and cannot be called donor-generalizing.

---

# 16. Primary metric

Use:

```text
negative log likelihood / log loss
```

Primary comparison:

```text
ΔNLL = NLL(M+X) - NLL(M)
```

Interpretation:

```text
ΔNLL < 0  -> early RNA helps
ΔNLL = 0  -> no incremental value
ΔNLL > 0  -> early RNA hurts
```

Why log loss:

- works with probabilistic predictions;
- evaluates probability quality, not only ranking;
- remains defined in cases where one held-out fold has only one class;
- aligns with the eventual CellFate-Rx goal of calibrated outcome probabilities.

---

# 17. Secondary metrics

Where mathematically defined:

```text
AUROC
PR-AUC
balanced accuracy
Brier score
ECE
```

Do not substitute a favorable secondary metric for a failed primary result.

---

# 18. Molecular-shuffle null

This is mandatory.

Purpose:

> Determine whether apparent `M+X` improvement comes from the molecular state or merely from time/treatment structure.

Shuffle transcriptomic profiles while preserving the relevant metadata strata as much as the dataset geometry allows.

Preferred shuffle:

```text
permute X within treatment × early-time × experiment strata
```

If strata are too small, use the nearest valid restricted scheme and record the exact reason.

Do not shuffle:

- future outcome;
- treatment identity;
- time;
- donor label;
- split assignment.

Run at least:

```text
1000 permutations
```

for the final Stage-21 null.

The real `M+X` gain must be stronger than the shuffled-X distribution.

---

# 19. Additional anti-leakage controls

The diagnostic must explicitly test:

## 19.1 Same-state leakage

Confirm that the future label is not computed from the same expression vector used as input.

---

## 19.2 Temporal leakage

Confirm that no gene expression measured after the prediction cutoff enters `X`.

---

## 19.3 Unit leakage

Confirm no clone / donor / replicate identifier appears in both training and test when that identifier defines the split unit.

---

## 19.4 Metadata leakage

Print every feature given to the metadata model.

A field such as:

```text
final_fate
response_status
post_treatment_cluster
survivor_label
```

must trigger a hard error if present in predictors.

---

# 20. Required tests before execution

`tests/test_diag_stage21_data_audit.py` must cover:

```text
dataset with valid lineage links -> STRICT_LINEAGE
dataset with culture links only -> CULTURE_FORWARD
timecourse with no direct link -> TRAJECTORY_FORWARD
cross-sectional dataset -> INVALID_PROSPECTIVE
D0-only dataset -> INVALID_PROSPECTIVE
future-derived input column -> hard failure
```

`tests/test_diag_stage21_forward_fate.py` must cover:

```text
metadata-only synthetic signal
    -> M+X does not receive false credit

RNA-only synthetic future signal
    -> M+X beats M

combined metadata + RNA signal
    -> M+X beats both components when identifiable

random RNA
    -> no systematic incremental benefit

shuffled RNA
    -> molecular gain disappears

test-fold scaling leakage
    -> guard fires

same biological unit in train and test
    -> guard fires

future-outcome feature included in X
    -> guard fires
```

Every verdict branch must be hit by at least one constructed test.

---

# 21. Pre-registered Stage 21 verdicts

## VERDICT A — `NO_VALID_FORWARD_TASK`

Condition:

```text
Stage 21A finds no LEVEL 1/2/3 task.
```

Action:

- stop;
- do not train;
- do not modify `_s16`;
- open a new data-acquisition stage;
- search for lineage-resolved or second dense single-cell timecourse data.

Scientific conclusion:

> Existing local datasets cannot pose the prospective question honestly.

---

## VERDICT B — `FORWARD_TASK_EXISTS_BUT_NO_INCREMENTAL_RNA`

Condition:

A valid task exists, but:

```text
M+X does not improve primary metric over M
```

or the improvement is inconsistent / compatible with zero.

Action:

- do not modify CellFate-Rx;
- do not launch prospective architecture development;
- decide whether to publish the current non-prospective work or acquire stronger prospective data.

Scientific conclusion:

> Existing data do not show that early RNA adds future-fate information beyond metadata.

---

## VERDICT C — `SUGGESTIVE`

Condition:

```text
M+X improves mean ΔNLL
```

but at least one is true:

- uncertainty includes zero;
- improvement is carried by one unit;
- shuffle null is not rejected;
- only LEVEL 1 trajectory-forward geometry exists.

Action:

- no architecture changes yet;
- acquire an independent prospective dataset;
- use result only to justify acquisition.

Scientific conclusion:

> Early molecular state is a plausible prospective signal but not yet established.

---

## VERDICT D — `PASS_PROSPECTIVE_SIGNAL`

Condition:

All must hold:

```text
1. Valid LEVEL 2 or LEVEL 3 forward geometry
2. M+X beats M on the primary metric
3. Improvement is not carried by one pathological unit
4. Real molecular gain exceeds the restricted-shuffle null
5. No leakage guard fires
```

Preferred stronger evidence:

```text
M+X also beats X
```

when treatment / metadata variation is sufficient to test that comparison.

Action:

- freeze Stage 21;
- open Stage 22;
- Stage 22 may adapt CellFate-Rx to true prospective labels;
- `_s16` remains preserved as the old system;
- do not yet claim a final paper result until independent replication.

---

# 22. Important interpretation rule

A large AUC is **not automatically a PASS**.

Examples:

```text
time-only AUC = 0.98
M+X AUC       = 0.98
```

Result:

```text
FAIL prospective molecular contribution
```

because RNA added nothing.

Example:

```text
M NLL   = 0.61
M+X NLL = 0.49
shuffle M+X NLL ≈ 0.61
```

Result:

```text
evidence for genuine molecular contribution
```

The paper-worthy quantity is the **incremental prospective value of X**, not the absolute metric.

---

# 23. Reporting per biological unit

Always print a per-unit table.

Minimum fields:

```text
unit
n_test
safe
unsafe
NLL_M
NLL_X
NLL_MX
delta_NLL_MX_vs_M
AUROC_M
AUROC_X
AUROC_MX
PR_AUC_M
PR_AUC_X
PR_AUC_MX
```

Do not allow an aggregate metric to hide that one donor / clone / replicate carries the entire result.

This is especially important because the current Stage-18 fate signal is heavily concentrated in Y1.

---

# 24. Required result JSON

`results/diag_stage21_data_audit_results.json` should contain:

```json
{
  "stage": 21,
  "phase": "data_audit",
  "datasets": [],
  "selected_dataset": null,
  "selected_task_level": null,
  "selected_forward_task": null,
  "verdict": null,
  "src_modified": false
}
```

If Stage 21B runs, `results/diag_stage21_forward_fate_results.json` should contain at least:

```json
{
  "stage": 21,
  "phase": "forward_fate",
  "dataset": "",
  "task_level": "",
  "biological_unit": "",
  "input_time": "",
  "outcome_time": "",
  "future_outcome_definition": "",
  "split_unit": "",
  "n_independent_units": 0,
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
  "shuffle_null": {
    "n_permutations": 1000,
    "p_value": null
  },
  "per_unit": [],
  "leakage_checks": {},
  "verdict": null,
  "src_modified": false
}
```

---

# 25. Console output

Keep the final console output deliberately plain.

Example:

```text
STAGE 21 — PROSPECTIVE FATE
===========================

DATA AUDIT
GSE242423    TRAJECTORY_FORWARD
GSE165177    CULTURE_FORWARD
GSE165178    INVALID_PROSPECTIVE
GSE165179    INVALID_PROSPECTIVE
GSE113957    INVALID_PROSPECTIVE
GSE297234    INVALID_PROSPECTIVE

PRIMARY TASK
Dataset:             GSE165177
Task level:          CULTURE_FORWARD
Independent units:   12

MODEL COMPARISON
Metadata NLL:        ...
RNA NLL:             ...
Metadata + RNA NLL:  ...
ΔNLL M+X vs M:       ...
95% CI:              ...

SHUFFLE NULL
Permutations:        1000
p:                   ...

VERDICT: ...
```

The numbers above are placeholders only.

---

# 26. Commands

Suggested audit command:

```bash
python experiments/diag_stage21_data_audit.py ^
  "D:\GSE242423" ^
  "D:\GSE165177" ^
  "D:\GSE165178" ^
  "D:\GSE165179" ^
  "D:\GSE113957" ^
  "D:\GSE297234"
```

Only if the audit licenses a forward task:

```bash
python experiments/diag_stage21_forward_fate.py ^
  --audit results/diag_stage21_data_audit_results.json
```

Then:

```bash
pytest -q
ruff check .
git diff --stat src/
```

Required final line:

```text
src/ unchanged
```

If `src/` changed, Stage 21 is invalid until the change is reverted or separately justified and pre-registered.

---

# 27. What Stage 21 does NOT prove

Even a PASS does **not** establish:

- arbitrary perturbation generalization;
- held-out-drug generalization;
- clinical use;
- working RES;
- accurate future ΔAge;
- safe rejuvenation;
- single-cell future fate if only population geometry was available;
- donor generalization if only one HFF line was evaluated;
- that CellFate-Rx itself beats logistic regression.

A PASS proves only:

> **The available prospective geometry contains reproducible future-fate information in early molecular state beyond allowed metadata, sufficient to justify a true prospective CellFate-Rx stage.**

---

# 28. What happens after Stage 21

## If `NO_VALID_FORWARD_TASK`

Open:

```text
STAGE 22 — Prospective Data Acquisition
```

Target:

- lineage-resolved pre-treatment transcriptome + future outcome, or
- second dense single-cell reprogramming timecourse on a different line / donor,
- preferably with multiple perturbations or doses.

---

## If `FORWARD_TASK_EXISTS_BUT_NO_INCREMENTAL_RNA`

Do not build a new prospective network.

Return to publication strategy.

The result means the original prospective hypothesis is not supported by the current data.

---

## If `SUGGESTIVE`

Use Stage 21 only as justification for one independent prospective dataset.

No paper headline yet.

---

## If `PASS_PROSPECTIVE_SIGNAL`

Open:

```text
STAGE 22 — Prospective CellFate-Rx
```

Stage 22 may then compare:

```text
metadata only
X only
X + U logistic regression
simple MLP
prospective CellFate-Rx
```

with true future labels.

Only there does architecture development resume.

---

# 29. Why this is the correct next stage

This stage is designed to prevent the three failure modes that already affected the project:

### 1. Same-state circularity

The ΔAge investigation showed that a model can achieve apparently excellent prediction when the target is mathematically reconstructable from the input.

Stage 21 forbids that geometry.

### 2. Timepoint confounding

The current fate analysis showed that marginal fate metrics can mostly measure elapsed time.

Stage 21 makes metadata-only performance the mandatory null.

### 3. Data geometry mistaken for model quality

The project already showed that bulk held-out data cannot express a per-cell unsafe fraction, and that line/modality shift can dominate apparent failure.

Stage 21 audits geometry before fitting anything.

The stage therefore asks the cheapest meaningful question first:

> **Is there a legitimate prospective signal worth building the next version of CellFate-Rx around?**

If the answer is no, we learn it without spending a retrain.

If the answer is yes, we have a pre-registered foundation for the paper-making research arc.

---

# 30. Freeze rule

Before executing either diagnostic:

1. commit this document;
2. commit all tests for the decision branches;
3. verify `src/` is clean;
4. run the audit once;
5. if the audit licenses Stage 21B, freeze the exact target/split fields in the audit JSON;
6. run the forward test once;
7. record the result additively in `CHANGES.md`;
8. do not alter the pass/fail rule after seeing the result.

---

# 31. One-sentence Stage 21 definition

> **Stage 21 asks whether the data already on disk can support an honest early-state → later-fate prediction task, and whether early RNA adds reproducible information beyond time/treatment metadata before any new CellFate-Rx development is allowed.**
