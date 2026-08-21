# STAGE 23 — Learnability / Interaction Gate

**Status:** DRAFT, PRE-EXECUTION-AUDITED — V1 is historical. Review and commit this V2 plan alone before fitting any model. No Stage-24 model architecture work may begin until Stage 23 is executed, frozen, and reviewed.  
**Plan version:** V2  
**Date opened:** 2026-08-21  
**Depends on:** Stage 22 `STAGE_23_READY`  
**Stage-22 benchmark implementation commit:** `8d6011a`  
**Stage-22 record commit:** `6b6a169`  
**Stage-22 additive record correction:** `9eaddf2`  
**Stage-22 plan commit:** `e97efe1`  
**Repository review base:** `c383bc5`  
**Scope:** establish whether pretreatment transcriptomic state contains prospective signal, whether that signal adds beyond captured clone abundance / treatment-only baselines, and whether a genuine treatment-specific `X × U` interaction is supported. This is a **learnability gate**, not the final production architecture.

---

## 0. Why Stage 23 is split into substages

Stage 23 contains three different scientific questions:

```text
Role A:
  Does pretreatment X predict the later Rewind priming outcome
  beyond trivial prevalence / captured clone-size baselines?

Role B:
  Does pretreatment X predict future treatment response
  beyond treatment U and captured naive clone abundance?

Interaction:
  Does explicit X × U improve held-out prediction
  beyond additive X + U + abundance?
```

Therefore execute:

```text
23A  Protocol / representation freeze
23B  Rewind Role-A learnability gate
23C  WM989 additive state-signal gate
23D  WM989 interaction gate
23E  Negative controls / robustness / leakage audit
23F  Evidence synthesis + Stage-24 gate
```

Every substage is additive. Do not rewrite an earlier substage because a later result is inconvenient.

The committed V1 remains historical. Its independent coding-agent/data-machine audit found one executable blocker (treatment-label casing) plus several pre-execution clarifications/strengthenings. V2 is the additive audited correction. Do not edit V1 in place.

V2 must be reviewed and committed **before** any Stage-23 estimator is fitted. The already-completed independent audit is not a model result and does not replace §3.0: 23A still reruns the executable Stage-22 input/hash audit before pseudobulk construction.

---

# 1. Frozen inputs from Stage 22

## 1.1 Rewind / GSE227151 — Role A

```text
retained cells                  3,905
retained clones                 3,147
primed cells                       42
primed clones                      35
nonprimed cells                  3,863
nonprimed clones                 3,112
expression mapping            3,905 / 3,905
outer split unit                    clone
outer folds                              5
generalization scope          within_R1_clone_heldout
```

Frozen positive clones per outer fold:

```text
7 · 7 · 7 · 7 · 7
```

Outcome:

```text
y_primed = author_top100_iPSC_gDNA_priming
```

It is an operational author-defined future lineage-grounded outcome, not proof that every `nonprimed` clone biologically failed.

Rewind R1 is one true biological replicate. Stage 23 can establish **within-R1 clone-held-out learnability only**.

## 1.2 WM989 / GSE279162 — Role B

```text
post-QC cells                   46,891
assigned cells                  42,771
naive cells                      6,489
naive-observed clones            1,401
clone × treatment rows           8,406
observed-zero rows               6,150   (73.2%)
nonzero rows                     2,256
expression mapping            6,489 / 6,489
outer split unit                    clone
outer folds                              5
```

Treatments:

```text
Dabrafenib
Trametinib
CoCl2
Acid
Cisplatin
Doxorubicin
```

Frozen nuisance warning:

```text
captured naive cells / clone    observed-zero rate
1                               87.1%
2                               82.6%
3-4                             67.9%
5-9                             59.9%
10+                             32.6%
```

Captured clone abundance is therefore a **mandatory competing baseline**. A valid Stage-23 result is that `X` adds little or nothing beyond it.

Frozen expression-feature anchors:

```text
total 10x features              189,656
Gene Expression features        36,601
Custom lineage features        153,055
```

Only the 36,601 Gene Expression rows are eligible for `X`.

### 1.2.1 Audit-derived nuisance / target-dependency facts

The independent pre-execution audit verified two facts that must be carried into Stage 23:

```text
1. The observed-zero confound is keyed strongly to TOTAL captured naive depth.
2. Recomputing author lineage assignment on treated cells only changes
   199 / 39,665 treated assignments = 0.50%.
```

Therefore the WM989 nuisance block `B` must include the explicit total-depth term:

```text
log1p(n_naive_cells)
```

in addition to the three per-naive-sample count terms. This is pre-registered **before fitting** because the frozen Stage-22 diagnostic already showed total captured depth strongly predicts observed-zero. Omitting it would make the nuisance baseline weaker than the measured confound and could bias the incremental `X` comparison in favor of transcriptomics.

The 0.50% assignment difference is an inherited author-pipeline dependency, not pretreatment-Gene-Expression leakage: Stage-22 target counts use post-treatment lineage assignments, and the audit reproduced all targets exactly without using pretreatment Gene Expression values. However, lineage-assignment parameters were estimated jointly and therefore have a small dependence on the presence/count structure of naive cells. Record this as a declared limitation. Do **not** rebuild the frozen Stage-22 target using a treated-only alternative.

## 1.3 Repository-derived Stage-22 preflight facts

The repository review at `c383bc5` verified:

```text
Stage-22 overall                         STAGE_23_READY
all_gates_pass                           true
all ten individual gates                 true
model_fitted                             false
six committed benchmark CSV hashes       match their manifests
two manifest hashes                      match Stage-22 results
WM989 barcode-clustering merge conflicts 0
```

These facts support opening Stage 23. However, Stage 23 must not trust only the serialized `overall` string:

- the current Stage-22 builder derives `overall` from the Role-A verdict alone;
- `all_gates_pass` is computed separately and is not consumed by that derivation;
- the test suite currently preserves that Role-A-only formula;
- Role B is intentionally non-blocking, but Role-A/global blocking gates and Role-B status are not yet represented as separate gate families.

This does **not** invalidate the current Stage-22 verdict because all ten gates are true. Stage-23 preflight must nevertheless fail closed unless the frozen current inputs show:

```text
Role-A verdict is ready
AND all_gates_pass is true
AND every required Role-A/global gate is true
AND model_fitted is false
```

Preserve Role B's verdict/limitations separately; do not let Role B rescue Role A.

Stage-22 G22-2 also has asymmetric evidence strength. Rewind is actually recomputed from future gDNA alone, but the builder currently records the WM989 no-pretreatment-expression conclusion through a literal `True` evidence field. The source path supports that conclusion, but Stage 23 must independently verify that every WM989 target is reconstructed only from post-treatment lineage assignments/counts and that no pretreatment Gene Expression value enters target construction.

## 1.4 Inherited text-hash portability limitation

Stage-22 result tables and manifest hashes are intact. Two provenance hashes are checkout-byte hashes and therefore differ only by text line endings across Windows/Linux:

```text
Stage-22 plan
  recorded Windows/CRLF SHA-256  b26db9fef09359791013d54c7f6e8ea4843d5ccfc3e2823f8794bd67a0bc1ff9
  Linux/LF SHA-256               fd671ad6ea0dff668883b11d2f645436257f0780efe1d9e034413662b34d766a

Stage-22 builder
  recorded Windows/CRLF SHA-256  785b9811cee22fe6b28173dbbe7bd109d519c81cba2592aa4410c6ed7a4cc7b9
  Linux/LF SHA-256               dfd047eff9fb82d0303c24743416d5a3ddc524076976074adceb435cb5d704a7
```

This is an inherited provenance limitation, not a benchmark-content mismatch. Do not rewrite Stage-22 tables to hide it.

For Stage 23, define hashes before implementation:

```text
binary/data artifact hash = SHA-256 of exact bytes
text source/protocol hash = SHA-256 after canonicalizing CRLF and CR to LF
```

Stage-23 determinism must be tested across LF and CRLF representations of the same text, not only by recloning on the development machine.

---

# 2. Global Stage-23 rules

## 2.1 Frozen outer folds

Use the exact Stage-22 `outer_fold` assignments.

Forbidden:

```text
new random seeds
repeated CV until a favorable split appears
rebalancing outer folds after results
dropping difficult folds
creating an alternative main split
```

All hyperparameter choices occur inside the four training folds. The outer test fold is never used for preprocessing, model selection, or hyperparameter choice.

Each observed scientific model produces exactly one frozen OOF prediction per outer-test observation. After those observed OOF predictions are frozen, pre-registered bootstrap/permutation inference in 23B–23E may reuse the held-out outcomes for inference only; it may not feed any result back into representation/model selection.

## 2.2 Independent unit

All inference/resampling uses:

```text
independent unit = clone
```

Never count Rewind cells or WM989 clone×treatment rows as independent biological units.

For WM989, all six treatment rows for one clone stay together.

## 2.3 Feature firewall

### PRIMARY_X

Pretreatment **Gene Expression** features only.

For WM989, the 10x `Custom` lineage/barcode features are explicitly forbidden from `X`; they encode clone identity/lineage assignment and would be direct provenance leakage. Only `feature_type == "Gene Expression"` enters pseudobulk expression.

For every dataset, use the stable 10x feature/gene ID as the canonical feature key rather than gene symbol or row position. If relevant samples do not have identical Gene Expression feature IDs/order, align by the shared feature-ID intersection in deterministic sorted order and report the loss. Never silently align matrices by row number.

### BASELINE / NUISANCE

Rewind:

```text
log1p(n_pretreatment_cells)
n_lanes
```

WM989:

```text
log1p(n_naive_cells)
log1p(n_naive1_cells)
log1p(n_naive2_cells)
log1p(n_naive3_cells)
```

### PROVENANCE-ONLY — never predictive

```text
clone_id
cell_uid
GSM / SampleNum / sample labels
source paths / expression column indices
outcome-rule strings
outer_group
outer_fold
file hashes
```

Treatment identity `U` is allowed only in explicitly treatment-conditioned models.

### TARGET — never in X

```text
Rewind:
  y_primed

WM989:
  n_post_cells
  post_fraction
  post_rank
  post_rank_fraction
  post_tie_size
  detected_post
  outcome_observation_status
  treatment_total_assigned_cells
```

`post_rank*` fields are reporting-only in Stage 23, not training targets.

### Stage-23 firewall override rule

Stage 23 is intentionally allowed to be **stricter** than Stage 22's feature metadata. The Stage-23 classes in this section are authoritative for Stage-23 modeling.

A field that was previously permitted as metadata/nuisance may be reclassified here as target-derived/post-treatment-forbidden without creating a contract failure. Cross-stage firewall tests must enforce:

```text
Stage 23 may tighten eligibility
Stage 23 may not loosen a Stage-22 forbidden/target/provenance field into PRIMARY_X
```

Do not require exact class-label equality across stages. Any tightened reclassification must be enumerated in `stage23_protocol.json` so the difference is explicit rather than accidental.

## 2.4 No representation shopping

23A freezes:

```text
clone expression construction
gene filtering
normalization
PCA candidate dimensions
model families
hyperparameter grids
metrics
bootstrap / permutation rules
gate definitions
```

After the first outer-fold model is fit, those rules are immutable.

If a bug is discovered, correct it additively and rerun the **entire affected comparison set**.

## 2.5 Learned preprocessing is nested and fold-local

There are **two leakage boundaries**:

```text
outer test
inner validation
```

For final outer-fold evaluation, every learned quantity is fitted from the outer-training clones only.

For inner-CV hyperparameter selection, every learned quantity must be **re-fitted inside each inner-training split**. It is not valid to fit a gene filter, scaler, or PCA once on the full outer-training set and then cross-validate only the downstream estimator.

Required inner pipeline:

```text
inner-training clones
-> training-only gene filter
-> training-only gene standardization
-> training-only PCA
-> training-only PC standardization
-> training-only nuisance-feature standardization
-> fit candidate estimator
-> score untouched inner-validation clones
```

After selecting hyperparameters, rebuild the same pipeline from scratch on the **entire outer-training set**, then transform/predict the outer test set once.

Per-clone CP10K/log1p normalization is allowed before these nested splits because it uses only that clone's own pretreatment counts and no across-clone information.

## 2.6 Implementation boundary from the actual repository

Stage 23 is a new experiment-layer learnability gate. Implement it under `experiments/` plus contract tests/results; do not alter `src/` or treat the existing production model as evidence.

The existing `CellFateNet` and `cellfate.evaluation.baselines._LinearBase` are not reused unchanged as Stage-23 estimators: they implement the older cell-level three-class fate/ΔAge interface and fixed baseline settings, not this clone-level prospective nested-CV protocol. Shared low-level utilities may be reused only when their behavior exactly matches this plan and is covered by Stage-23 tests.

---

# 3. Stage 23A — Protocol / representation freeze

**Purpose:** create one fixed clone-level `X_before` representation and freeze evaluation before predictive results.

No outcome model may be fitted until 23A is committed.

## 3.0 Mandatory pre-execution/input audit

Before constructing `X`, execute and record an audit against the committed Stage-22 tables/manifests and the external raw-data roots. This audit is part of 23A but occurs before any outcome is exposed to a fitted estimator.

Required checks:

```text
current Stage-22 counts/schemas/folds match Section 1
all six benchmark CSV hashes match their manifests
both manifest hashes match stage22_prospective_benchmark_results.json
Role-A verdict ready; all_gates_pass true; model_fitted false
the inherited CRLF/LF provenance-hash limitation is recorded, not mistaken for data corruption
every expression_source exists on the data machine
every required raw matrix/features/barcodes file matches the frozen source size/hash
every expression_column_index resolves to the recorded barcode
Gene Expression feature IDs/types are readable for every source sample
WM989 target construction has no dependency on pretreatment Gene Expression values
no Stage-22 artifact, fold, label, or manifest is rewritten
```

Verdict:

```text
STAGE22_INPUTS_AUDITED
STAGE22_INPUTS_BLOCKED
```

If blocked, stop before pseudobulk or model fitting.

## 3.1 Clone pseudobulk

### Rewind

```text
raw pseudobulk count =
    sum raw pretreatment RNA counts across all retained cells of clone
```

The eight Stage-22 ambiguous cells remain excluded.

### WM989

```text
raw pseudobulk count =
    sum raw RNA counts across all assigned Naive1/Naive2/Naive3 cells of clone
```

No treated cell may enter `X_before`.

## 3.2 Frozen expression feature universe

Use only 10x **Gene Expression** features.

WM989:

```text
include: feature_type == "Gene Expression"
exclude: feature_type == "Custom"
exclude: every lineage / LinNNNNN feature
```

Rewind:

```text
use the gene-expression feature table only
```

For each dataset:

```text
canonical feature key = stable feature/gene ID
```

Assert feature-ID uniqueness. If all included samples have identical Gene Expression feature hashes/order, preserve that order. Otherwise, align to the deterministic sorted intersection of feature IDs and report:

```text
features per source sample
shared features retained
features dropped per sample
```

Gene symbols are metadata only and are never used as the join key.

## 3.3 Primary normalization

Per clone:

```text
CP10K_g = 10,000 * count_g / total_clone_counts
X_g     = log1p(CP10K_g)
```

A zero-total clone is a Stage-23A block, not an imputation case.

The repository's `normalize_counts` already applies both CP10K and `log1p` to cell-by-gene inputs. Stage 23 may use a sparse-equivalent implementation, but it must sum **raw** cells into clone pseudobulk first and normalize exactly once. Never sum already-normalized cells and never apply a second `log1p`.

## 3.4 Training-fold gene filter

Within each outer training set:

```text
detected in >= max(5, ceil(0.01 * n_training_clones))
and
non-zero variance
```

Apply that training-derived gene list unchanged to test.

## 3.5 Standardization + PCA

Within the applicable training split only:

```text
standardize retained genes by training mean / SD
fit PCA on training clones
take first K PC scores
standardize those K PC scores again by training mean / SD
```

Freeze PCA implementation as:

```text
sklearn.decomposition.PCA
svd_solver="randomized"
random_state=23023
```

Within a given training split, fit PCA once at the largest mathematically feasible candidate `K` (up to 50) and use prefixes of that same basis for K=10/20/50. Do not refit a different randomized PCA separately for each K.

The second standardization makes L2 regularization comparable across PCs with different eigenvalue scales.

Every continuous nuisance feature (`B`) is likewise standardized using training-only mean / SD before model fitting. Treatment dummy variables are never standardized.

Candidate dimensions:

```text
K = {10, 20, 50}
```

Skip a K only if mathematically impossible for that training matrix. A skipped K must be recorded; no replacement K may be invented.

## 3.6 Inner CV

Rewind:

```text
3-fold StratifiedKFold
shuffle=True
random_state=23023
unit=clone
```

Every inner fold must contain both classes.

Inner selection score:

```text
maximize mean Average Precision (sklearn average_precision_score)
```

WM989:

```text
3-fold GroupKFold
group=clone_id
```

Sort clone IDs deterministically before constructing inner folds.

For C1, the **global benchmark contains 1,401 eligible clones**. Within any outer split, only that split's frozen outer-training clones enter inner CV (approximately 1,120 clones), each represented by all six treatment rows.

For C2, keep the frozen Stage-22 outer-fold assignment but restrict the conditional-abundance dataset **after the outer split is known** to clones/rows with `n_post_cells > 0`. Inner GroupKFold is therefore built from the outer-training clones that have at least one nonzero C2 row, with all of each clone's nonzero treatment rows kept together.

Inner selection scores:

```text
C1 detection:
  minimize mean clone-balanced log loss

C2 conditional abundance:
  minimize mean clone-balanced MAE
```

For C1 every clone has six treatment rows, so ordinary row-average log loss is already clone-balanced.

For C2, clones contribute different numbers of nonzero treatment rows. Define each clone's validation MAE across its own nonzero rows first, then average those clone MAEs. During C2 fitting, each nonzero row receives weight:

```text
1 / number_of_nonzero_treatment_rows_for_that_clone_in_the_training_subset
```

renormalized to mean sample weight 1. This prevents clones observed under many treatments from silently becoming more independent training units.

### Deterministic hyperparameter tie-break

If candidate mean inner scores are equal within `1e-12`, prefer in order:

```text
1. smaller K
2. stronger regularization
   logistic: smaller C
   ridge:    larger alpha
```

For models without `K`, apply rule 2 only.

The inner selection metric is fixed by endpoint; secondary metrics may not choose hyperparameters.

## 3.7 Frozen model grids

Binary logistic:

```text
LogisticRegression
penalty="l2"
solver="liblinear"
C ∈ {0.01, 0.1, 1, 10}
fit_intercept=True
class_weight=None
max_iter=5000
random_state=23023
```

`class_weight` is deliberately frozen to `None` for all logistic comparisons. Rewind imbalance is handled by the **Average Precision** primary metric rather than by changing the effective class prior; this keeps log-loss/Brier probabilities interpretable as unweighted model probabilities.

No SMOTE, synthetic positives, class reweighting, or threshold tuning.

Continuous regression:

```text
Ridge
alpha ∈ {0.1, 1, 10, 100}
fit_intercept=True
```

No boosting, trees, neural nets, wrapper feature selection, or architecture search in Stage 23.

Do not route these comparisons through the existing production `CellFateNet` or the existing fixed cell-level baseline registry. The Stage-23 estimator pipeline must expose its fold-local preprocessing and selected hyperparameters explicitly.

Any convergence warning/error from a candidate estimator is a **protocol failure to investigate**, not permission to silently discard that candidate after seeing its score. If the frozen grid cannot be evaluated as specified, stop the affected substage and record the failure additively.

## 3.8 Frozen treatment coding

Canonical treatment order, matching the frozen Stage-22 table **exactly and case-sensitively**:

```text
Acid
Cisplatin
CoCl2
Dabrafenib
Doxorubicin
Trametinib
```

Use:

```text
reference treatment = Acid
intercept = enabled
five treatment dummy variables for the non-reference treatments
```

The exact same treatment coding is used in W0/W1/W3/W4/W5 for both endpoints.

For W5, interaction terms are constructed only from:

```text
standardized PC score × non-reference treatment dummy
```

The common `X` coefficients therefore represent the reference-treatment state contribution, while interaction terms represent treatment-specific deviations. Predictions, not individual coefficients, are the inferential object.

## 3.9 23A artifacts/tests

Suggested:

```text
results/stage23_protocol.json
results/stage23_rewind_clone_expression_manifest.json
results/stage23_wm989_clone_expression_manifest.json
results/stage23_outer_fold_preprocessing.json
```

Tests:

```text
Stage-22 input audit passes before any estimator fit
Stage-22 `overall` is not trusted without all_gates_pass/required gates
Stage-22 benchmark/manifest artifact hashes verify
WM989 no-pretreatment-expression target dependency is independently checked
canonical text hash is identical for LF and CRLF representations
no treated WM989 cell in X
no excluded Rewind cell in X
one pseudobulk row per clone
WM989 Custom/lineage features absent from X
Gene Expression feature IDs aligned by stable ID, never row number
no target/provenance field in PRIMARY_X
inner-training-only gene filter/scaling/PCA during inner CV
outer-training-only final preprocessing before outer prediction
continuous nuisance features standardized training-only
WM989 nuisance B includes log1p(n_naive_cells) plus the three per-naive-sample terms
Stage-23 stricter firewall reclassification is explicit and never loosens Stage-22 forbidden fields
treatment dummy coding exactly frozen with case-sensitive reference `Acid`
outer folds exactly equal Stage 22
deterministic rerun
```

Verdict:

```text
PROTOCOL_FROZEN
PROTOCOL_BLOCKED
```

If blocked, stop.

---

# 4. Stage 23B — Rewind Role-A learnability

Question:

```text
Does pretreatment X predict future priming beyond prevalence and captured clone size?
```

## 4.1 Models

### R0 — prevalence

Constant outer-training primed prevalence.

### R1 — nuisance only

```text
log1p(n_pretreatment_cells)
n_lanes
```

### R2 — X only

```text
PCA(X)
```

### R3 — X + nuisance

```text
PCA(X)
log1p(n_pretreatment_cells)
n_lanes
```

No other primary Rewind model.

## 4.2 Metrics

Primary:

```text
Average Precision (AP; sklearn `average_precision_score`) at clone grain
```

Secondary:

```text
ROC-AUC
log loss
Brier score
```

Accuracy is reporting-only.

## 4.3 Primary comparisons

Absolute expression signal:

```text
R2 vs R0
```

Primary incremental signal:

```text
R3 vs R1
ΔAP_state = AP(R3) - AP(R1)
```

R2 beating R0 is insufficient if R3 fails to beat R1.

## 4.4 Uncertainty

Stratified clone bootstrap on pooled OOF predictions:

```text
2,000 replicates
seed=23123
resample positive clones with replacement
resample negative clones with replacement
preserve class counts
```

Report:

```text
ΔAP point
95% percentile CI
bootstrap fraction ΔAP <= 0
```

## 4.5 Fold stability

Report `AP(R3)-AP(R1)` for every outer fold.

Each outer fold contains only seven positive clones, so fold-wise AP differences are high-variance diagnostics. They must be reported prominently but **do not create a separate PASS requirement** and must not override the pooled clone-level OOF bootstrap/permutation evidence.

## 4.6 Provisional verdict

`ROLE_A_SIGNAL_PASS` requires:

```text
pooled ΔAP > 0
95% CI lower bound > 0
```

Fold-wise direction remains a robustness diagnostic only.

`ROLE_A_SIGNAL_WEAK`:

```text
pooled ΔAP > 0
but 95% CI includes 0
```

`ROLE_A_SIGNAL_FAIL`:

```text
pooled ΔAP <= 0
```

or invalid leakage/instability.

---

# 5. Stage 23C — WM989 additive state-signal gate

Do not create an arbitrary binary resistance threshold.

Use two descriptive future endpoints.

## 5.1 C1 — future clone detection

```text
detected_post = 1 if n_post_cells > 0
                0 if n_post_cells == 0
```

This means detection in the treated sample, not survival/death/resistance.

Primary:

```text
log loss
```

Secondary:

```text
average precision
ROC-AUC
Brier
```

## 5.2 C2 — conditional future abundance

Only rows with:

```text
n_post_cells > 0
```

Target:

```text
y_abundance = log1p(n_post_cells)
```

Primary:

```text
clone-balanced MAE
```

where MAE is computed within each clone across that clone's nonzero treatment rows, then averaged across clones.

Secondary:

```text
clone-balanced RMSE
per-treatment Spearman
mean of six treatment Spearmans
```

For RMSE, compute each clone's mean squared error over its nonzero rows, average those clone MSEs, then take the square root.

Treatment sample depth is not a predictor; treatment one-hot absorbs treatment-specific mean shifts.

## 5.3 Feature blocks

```text
U = treatment one-hot

B =
  log1p(n_naive_cells)
  log1p(n_naive1_cells)
  log1p(n_naive2_cells)
  log1p(n_naive3_cells)

X = Stage-23A pretreatment pseudobulk PCA
```

## 5.4 Models

For C1 and C2:

```text
W0 = U
W1 = B + U
W2 = X
W3 = X + U
W4 = X + B + U
```

W1 is the load-bearing nuisance baseline.

W4 is the primary additive transcriptomic model.

No `X×U` until 23D.

## 5.5 Primary state comparison

Detection:

```text
ΔLL_state = logloss(W1) - logloss(W4)
```

Conditional abundance:

```text
ΔMAE_state = MAE(W1) - MAE(W4)
```

Positive is better.

Also report W3 vs W0 and W2, but they cannot replace W4 vs W1.

## 5.6 Uncertainty

Clone-cluster bootstrap:

C1 detection:

```text
2,000 replicates
seed=23223
sample 1,401 clones with replacement
carry all six treatment rows for every sampled clone
```

C2 conditional abundance:

```text
2,000 replicates
seed=23224
sample the frozen set of clones with >=1 nonzero treatment row with replacement
carry all nonzero treatment rows for every sampled clone
```

Under the Stage-22 benchmark this C2 clone set is expected to contain **929 clones** (`>=1 treatment` coverage); assert/recompute that count before inference.

Preserve the clone-balanced metric definition above.

Report percentile:

```text
95% CI
97.5% two-sided CI
```

The 97.5% two-sided intervals are the Bonferroni-safe primary intervals for the two endpoints **within this additive-state hypothesis family**.

23C additive-state and 23D interaction are two separately pre-registered scientific hypothesis families. Their verdicts are reported separately; do not combine their p-values or present one as a post-hoc rescue of the other.

## 5.7 Verdict

`ROLE_B_ADDITIVE_PASS` requires:

```text
at least one:
  lower 97.5% CI(ΔLL_state)  > 0
  lower 97.5% CI(ΔMAE_state) > 0

and

the other endpoint is not significantly worse
(its upper 97.5% CI is not < 0)
```

`ROLE_B_ADDITIVE_WEAK`:

```text
one/both point improvements > 0
but neither lower 97.5% bound clears 0
```

`ROLE_B_ADDITIVE_FAIL`:

```text
both point improvements <= 0
```

or significant harm without compensating evidence.

This tests whether an **additive** X contribution adds beyond B+U. A failure here does not preclude a treatment-specific interaction that cancels in the additive average; 23D handles that case explicitly.

---

# 6. Stage 23D — WM989 explicit interaction gate

Question:

```text
Does the contribution of pretreatment state depend on treatment?
```

This is the gate for later state-conditioned treatment ranking.

## 6.1 Interaction features

Let PCA(X)=Z and treatment one-hot=U.

Use only:

```text
Z_j × U_t
```

No full gene-level interaction matrix.

## 6.2 Models

Reference:

```text
W4 = X + B + U
```

Interaction:

```text
W5 = X + B + U + X×U
```

Same:

```text
K={10,20,50}
regularization grids
```

No wider search for W5.

## 6.3 Improvements

Detection:

```text
ΔLL_interaction = logloss(W4) - logloss(W5)
```

Abundance:

```text
ΔMAE_interaction = MAE(W4) - MAE(W5)
```

Also compute the direct **full treatment-conditioned state contribution beyond nuisance**:

```text
Detection:
  ΔLL_full = logloss(W1) - logloss(W5)

Abundance:
  ΔMAE_full = MAE(W1) - MAE(W5)
```

This direct comparison is required because a real interaction may exist even when the additive W4-vs-W1 effect is weak or cancels across treatments.

Use the same 2,000-clone bootstrap and 95% / 97.5% two-sided percentile CIs.

## 6.4 Treatment-level diagnostics

For each treatment report:

```text
detection log-loss improvement
conditional-abundance MAE improvement
```

One favorable treatment is not a broad interaction claim.

## 6.5 Verdicts

`INTERACTION_PASS_MULTI_TREATMENT` requires, on the **same endpoint family**:

```text
1. lower 97.5% CI of W5-vs-W4 interaction improvement > 0
2. lower 97.5% CI of W5-vs-W1 full-state improvement > 0
3. on the other endpoint, neither the W5-vs-W4 nor W5-vs-W1 upper 97.5% CI is < 0
4. W5-vs-W4 improves directionally in >=3/6 treatments
```

This prevents an interaction model from "passing" merely by rearranging error relative to W4 while still failing to beat the load-bearing nuisance baseline W1.

`INTERACTION_LOCAL_ONLY`:

```text
at least one interaction/full-state point improvement is > 0
and there is directional improvement in at least one treatment
but the pre-registered multi-treatment PASS criterion fails
```

This is descriptive/local evidence only and cannot keep Stage-25 broad ranking eligibility.

`INTERACTION_NOT_SUPPORTED`:

```text
W5 fails to improve over W4
```

or reproducibly worsens it.

---

# 7. Stage 23E — Negative controls / robustness / leakage audit

No PASS survives without 23E.

## 7.1 Expression-permutation control

Destroy X↔outcome linkage while preserving captured-abundance structure.

Rewind permutation strata:

```text
n_pretreatment_cells:
  1
  2
  3+
cross with n_lanes where possible
```

If too small, merge deterministically with adjacent size stratum and record it.

WM989 strata preserve both the dominant total-depth effect and naive-library presence pattern:

```text
depth bin:
  1
  2
  3-4
  5-9
  10+

crossed with:
  3-bit naive sample presence pattern
  (Naive1 present?, Naive2 present?, Naive3 present?)
```

If a resulting stratum has fewer than 4 clones, merge it deterministically with the nearest stratum **within the same depth bin** by minimum Hamming distance of the presence pattern; ties break lexicographically. If the depth bin still cannot support permutation, merge with the adjacent depth bin whose midpoint is closest and record the merge.

This negative control preserves the dominant captured-depth/sample-presence structure; it is a secondary robustness test, not a claim of exact conditional randomization on every nuisance count.

The pre-execution audit traced the tiny strata and confirmed the merge rule terminates deterministically. Residual singleton strata can leave at most about 1.1% of outer-test clones fixed under a permutation (and none of the outer-training clones in the audited geometry). Record the actual fixed-clone fraction per fold/permutation family rather than treating a fixed singleton as a failed permutation.

Within each outer fold, permute the **whole clone-level CP10K/log1p Gene Expression profile** as one intact vector before training-derived gene filtering/PCA:

```text
permute training clone profiles among training clones only
permute test clone profiles among test clones only
```

Never permute genes independently within a clone, and never move an expression profile across the outer train/test boundary.

The nested preprocessing/model-selection pipeline is then rerun from the permuted clone-profile mapping.

## 7.2 Permutations

Run permutation testing **only for a claim that is otherwise eligible for PASS** after 23B/23C/23D. Failed/weak claims are recorded as:

```text
PERMUTATION_NOT_REQUIRED_NO_PASS_CANDIDATE
```

For each eligible claim:

```text
200 permutations
base seed=23323
```

Rerun the relevant fitted comparison using the frozen nested-CV algorithm and grids. Do not fix observed-data hyperparameters inside the null runs.

### Runtime rule

The independent audit measured representative PCA operations at roughly sub-second to ~1.5-second scale on the audit machine, so the pre-registered 200-permutation design is computationally substantial but feasible. Runtime is **not** permission to weaken the null procedure after observed results are known.

Forbidden without a new plan version:

```text
reducing 200 permutations because the run is slow
reusing observed-data hyperparameters in null runs
skipping inner preprocessing refits that the permutation changes
stopping permutations early because significance/non-significance looks obvious
```

If the exact frozen protocol is unexpectedly infeasible on the execution machine, stop the affected substage and report the runtime blocker before changing the procedure.

Caching is allowed only when it is mathematically invariant to the permutation.

Because a within-outer-training permutation preserves the **full outer-training profile set**, the final outer-training unsupervised transform may be cached if its inputs are otherwise identical.

Do **not** cache inner-split gene filters/scalers/PCA across permutations: permuting profiles among clone IDs changes which profiles land in each inner-training split, so those transforms must be recomputed.

## 7.3 Permutation statistics

Rewind null:

```text
ΔAP_state
```

Observed Role-A signal must exceed the 95th percentile of permuted ΔAP.

For every permutation-tested claim report the finite-sample empirical tail probability:

```text
p_perm = (1 + number_of_null_statistics >= observed_statistic) / (200 + 1)
```

A permutation gate passes only when:

```text
observed > 95th percentile(null)
and
p_perm <= 0.05
```

WM989 additive nulls:

```text
ΔLL_state
ΔMAE_state
```

Any endpoint used for `ROLE_B_ADDITIVE_PASS` must satisfy the permutation gate above.

WM989 interaction nulls:

```text
ΔLL_interaction
ΔMAE_interaction
ΔLL_full
ΔMAE_full
```

For any endpoint used for `INTERACTION_PASS_MULTI_TREATMENT`, **both** the W5-vs-W4 interaction improvement and the W5-vs-W1 full-state improvement must satisfy the permutation gate above.

## 7.4 Provenance sentinel

Run an explicitly frozen diagnostic sentinel that cannot encode clone identity:

Rewind:

```text
GSM7092515 present?   0/1
GSM7092516 present?   0/1
```

WM989:

```text
Naive1 present?   0/1
Naive2 present?   0/1
Naive3 present?   0/1
```

Use the same endpoint-appropriate logistic/ridge family and nested outer/inner geometry, but **no expression and no clone ID**.

These are provenance-presence flags only; captured counts remain in the scientific nuisance baseline `B`, not in this sentinel.

Never one-hot `clone_id`.

If this sentinel unexpectedly approaches/reproduces the main claimed gain, raise a confounding alert and inspect sample/library structure before 23F.

## 7.5 Optional outcome-shuffle engineering sentinel

A small deterministic training-label shuffle may be run as a **diagnostic only**:

```text
50 shuffles
training side only
evaluate on the unchanged outer test outcomes
```

It does not create or rescue a PASS and has no automatic scientific gate. If it unexpectedly matches the observed claimed gain, raise a leakage alert for manual investigation before 23F.

## 7.6 Outer-test isolation audit

Prove structurally:

```text
outer-test y never enters preprocessing
outer-test X never enters training PCA
outer-test y never enters inner selection
Stage-22 outer_fold unchanged
no treated WM989 expression in X
no target/provenance column in PRIMARY_X
```

Do not hard-code PASS booleans.

## 7.7 Determinism

A fresh-clone rerun must reproduce committed compact Stage-23 artifacts byte-for-byte without rewriting frozen Stage-22 artifacts.

For text protocol/source provenance, use the canonical LF-normalized hash defined in §1.4. Add an explicit contract that the same text encoded with LF versus CRLF produces the same canonical digest. Do not claim cross-platform determinism from a same-machine fresh clone alone.

## 7.8 23E control verdict

Compute:

```text
STRUCTURAL_CONTROLS_PASS =
  outer-test isolation audit passes
  AND feature-firewall audit passes
  AND frozen folds unchanged
  AND fresh-clone determinism passes
```

For each scientific PASS candidate, compute its claim-specific permutation status:

```text
ROLE_A_PERMUTATION_PASS
ROLE_B_ADDITIVE_PERMUTATION_PASS
ROLE_B_INTERACTION_PERMUTATION_PASS
```

A claim may be promoted to final PASS only when:

```text
STRUCTURAL_CONTROLS_PASS
AND its required permutation status == PASS
```

Weak/failed claims do not need permutation testing and remain weak/failed.

The provenance sentinel and optional outcome-shuffle sentinel are diagnostic alerts. If either exposes an actual forbidden-feature/code-path leak, `STRUCTURAL_CONTROLS_PASS = false`; otherwise they do not create a separate arbitrary threshold.

---

# 8. Stage 23F — Evidence synthesis

23F fits no new scientific model. It reads frozen 23B–23E outputs and derives the next gate.

## 8.1 Record structure

Prefer one authoritative:

```text
stage_23_RECORD.md
```

with additive immutable sections:

```text
23A
23B
23C
23D
23E
23F
```

If repository convention strongly favors separate records, use 23A–23F files, but do not maintain conflicting duplicate authorities.

## 8.2 Required final evidence table

Role A:

```text
R0/R1/R2/R3 AP
R3-R1 ΔAP + CI
fold-wise ΔAP
permutation result
verdict
```

Role B additive:

```text
C1 W1 vs W4 log loss
ΔLL + CI
C2 W1 vs W4 MAE
ΔMAE + CI
permutation result
verdict
```

Role B interaction:

```text
C1 W4 vs W5
C2 W4 vs W5
C1/C2 W1 vs W5 full-state comparisons
treatment-wise directions
permutation result
verdict
```

Nuisance-only performance must be prominent, not buried.

---

# 9. Stage-23 gate logic

## 9.1 Role A mandatory

`STAGE_24_READY` requires:

```text
23A = PROTOCOL_FROZEN
23B provisional verdict = ROLE_A_SIGNAL_PASS
23E STRUCTURAL_CONTROLS_PASS = true
23E ROLE_A_PERMUTATION_PASS = true
```

23F then promotes the Role-A verdict to final `ROLE_A_SIGNAL_PASS`.

Role B cannot override Role-A failure.

## 9.2 Role A weak

If:

```text
ROLE_A_SIGNAL_WEAK
```

then:

```text
STAGE_24_HOLD_ROLE_A_WEAK
```

Do not compensate by escalating architecture.

## 9.3 Role A fail

If:

```text
ROLE_A_SIGNAL_FAIL
```

then:

```text
STAGE_24_BLOCKED_ROLE_A
```

Role-B results remain reportable but do not silently replace the mandatory anchor.

## 9.4 Role-B scope if Role A passes

In this section, a Role-B `PASS` means the 23C/23D bootstrap criterion **and** the corresponding 23E permutation gate both passed under `STRUCTURAL_CONTROLS_PASS`.

If:

```text
INTERACTION_PASS_MULTI_TREATMENT
```

then W5 has already been required to beat both W4 and the nuisance baseline W1 on the passing endpoint. Stage 24 may therefore use explicit treatment-conditioned interaction **even if the additive W4-vs-W1 gate was weak/failed because effects canceled across treatments**. Stage 25 state-conditioned ranking remains eligible.

If:

```text
ROLE_B_ADDITIVE_PASS
INTERACTION_NOT_SUPPORTED
```

Stage 24 must keep Role B additive. Do not claim clone-specific treatment ranking from `X×U`. Stage 25 state-conditioned ranking is blocked or redesigned.

If:

```text
INTERACTION_LOCAL_ONLY
```

only a narrow treatment-specific interaction observation is allowed; broad six-treatment ranking remains blocked unless the explicit multi-treatment interaction gate passes.

If:

```text
ROLE_B_ADDITIVE_WEAK or ROLE_B_ADDITIVE_FAIL
AND
INTERACTION_NOT_SUPPORTED
```

Role B is scoped down/removed from the main prospective contribution. This does not block Stage 24 if Role A passed.

---

# 10. What Stage 23 may establish

A PASS may establish:

```text
pretreatment transcriptomic X contains held-out prospective information
beyond explicitly measured captured clone abundance / treatment baselines
within the frozen benchmark geometry
```

An interaction PASS may additionally establish:

```text
the predictive contribution of pretreatment state varies by treatment
in held-out clones across multiple treatments
```

---

# 11. What Stage 23 cannot establish

Even a full PASS does not establish:

```text
external biological-replicate generalization
unseen-treatment generalization
independent-dataset replication
probability calibration
OOD validity
causal treatment effect
lab decision utility
final architecture superiority
```

Rewind remains one biological replicate.

---

# 12. Later-stage boundaries

If `STAGE_24_READY`:

```text
Stage 24 — prospective model
```

Do not move these into Stage 23:

```text
Stage 25 — treatment ranking
Stage 26 — held-out perturbation / unseen-treatment generalization
Stage 27 — independent replication
Stage 28 — calibration / OOD / decision utility
```

---

# 13. Suggested engineering deliverables

```text
experiments/run_stage23_learnability_gate.py
tests/test_stage23_learnability_gate.py

results/stage23_protocol.json
results/stage23_rewind_oof_predictions.csv
results/stage23_rewind_results.json
results/stage23_wm989_detection_oof.csv
results/stage23_wm989_abundance_oof.csv
results/stage23_wm989_results.json
results/stage23_permutation_results.json
results/stage23_gate_results.json

plans/(newer)practical plans/RECORDs/stage_23_RECORD.md
```

Do not commit:

```text
large duplicated expression matrices
raw public data
temporary inner-CV predictions
200 full permutation prediction dumps
```

Commit compact summaries + exact protocol/seeds.

---

# 14. Required automated contracts

## 23A

```text
STAGE22_INPUTS_AUDITED before any estimator fit
Stage-22 overall cross-checked against all_gates_pass and required gates
Stage-22 CSV/manifest hash chain verified
inherited Stage-22 CRLF/LF source-hash mismatch declared, not silently rewritten
WM989 target dependency audit proves no pretreatment Gene Expression enters Y
canonical text hashing is LF/CRLF invariant
Stage-22 fold IDs unchanged
Rewind clone count = 3,147
WM989 naive clone count = 1,401
no treated WM989 cell in X
no ambiguous Rewind cell in X
one pseudobulk vector per clone
WM989 Custom/lineage features excluded from X
feature alignment by stable gene ID
inner-training-only preprocessing during inner CV
outer-training-only refit for outer prediction
PCA solver/seed frozen and prefix basis reused across K
convergence warnings are not silently ignored
```

## 23B

```text
R0/R1/R2/R3 exactly present
Average Precision (`average_precision_score`) is Rewind primary
no accuracy-based gate
one OOF prediction per clone
35 positive clones exactly once in pooled OOF
fold-wise ΔAP is reported but is not a PASS requirement
```

## 23C

```text
C1 detected_post never renamed resistance/survival
C2 scores only n_post_cells > 0
C2 fitting/scoring is clone-balanced
C2 expected eligible clone count = 929 is asserted/recomputed
W0-W4 exactly present
W1 is abundance+treatment baseline
W4-vs-W1 is primary additive comparison
inner metric/tie-break fixed before fitting
```

## 23D

```text
W5 = X+B+U+X×U
`Acid` is the frozen case-sensitive reference treatment with five non-reference dummies
interaction only standardized PCA(X) × frozen treatment dummies
no gene-level full interaction
W4 remains additive reference
W5-vs-W1 full-state comparison present
```

## 23E

```text
permutations preserve outer boundary
permutations preserve clone clustering
Rewind permutations preserve nuisance strata
WM989 permutations preserve depth + naive-sample-presence strata
empirical permutation p computed with +1 correction
forbidden provenance absent
no hard-coded PASS booleans
23E structural/pass logic mechanically derived
```

## 23F

```text
overall gate mechanically derived
Role A mandatory
Role B cannot override Role-A failure
interaction + additive verdicts control Stage-24/25 scope
```

---

# 15. Substage execution sequence

```text
1. preserve/archive committed Stage-23 V1 as historical

2. review and commit this audited Stage-23 V2 plan alone
   do not fit models in the plan commit

3. 23A
   rerun the executable §3.0 input/hash audit
   require STAGE22_INPUTS_AUDITED before pseudobulk construction
   freeze protocol / clone X
   tests + ruff
   append record
   commit

4. 23B
   Rewind learnability
   append record
   commit

5. 23C
   WM989 additive state signal
   append record
   commit

6. 23D
   WM989 interaction
   append record
   commit

7. 23E
   permutation / leakage / determinism
   append record
   commit

8. 23F
   derive final gate
   append final record
   full pytest
   CI-scope ruff
   fresh-clone reproducibility
   commit/push

9. STOP for review
```

Do not begin Stage 24 in the Stage-23 completion change.

---

# 16. Completion report format

```text
commit(s)
tests
ruff
fresh-clone determinism

23A:
  Stage-22 input-audit verdict
  Stage-22 hash-chain verification
  canonical text-hash portability check
  independent WM989 target-dependency audit
  verdict
  pseudobulk dimensions
  Gene Expression feature universe
  retained outer-training genes per fold
  candidate PCA K values
  nested-preprocessing audit

23B Rewind:
  R0/R1/R2/R3 AP
  selected K/C per outer fold for R2/R3
  selected C per outer fold for R1
  R3-R1 ΔAP
  95% CI
  fold-wise sign
  permutation result
  ROLE_A verdict

23C WM989:
  nuisance B definition including total captured depth
  inherited 0.50% joint-assignment dependency declared
  selected hyperparameters per outer fold/model
  detection W1 vs W4
  Δ log loss + CI
  abundance W1 vs W4
  C2 eligible clone count
  Δ MAE + CI
  permutation result
  ROLE_B_ADDITIVE verdict

23D:
  selected hyperparameters per outer fold
  W4 vs W5 detection
  W4 vs W5 abundance
  W1 vs W5 full-state detection/abundance
  treatment-wise directions
  permutation result
  INTERACTION verdict

23E:
  structural-controls verdict
  claim-specific permutation p/percentiles
  provenance sentinel
  optional outcome-shuffle alert
  determinism

23F:
  Stage-24 gate
  Stage-24 Role-B architecture scope
  Stage-25 ranking eligibility
```

Stop for review.
