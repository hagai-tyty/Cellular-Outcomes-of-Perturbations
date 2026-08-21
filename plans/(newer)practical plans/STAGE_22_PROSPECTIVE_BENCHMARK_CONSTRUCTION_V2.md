# STAGE 22 — Prospective Benchmark Construction

**Status:** PLANNED — do not begin Stage 23 until the Stage 22 record is frozen.  
**Plan version:** V2
**Date opened:** 2026-08-21  
**Depends on:** Stage 21D revision 2 (`STAGE_22_READY`)  
**Scope:** data/benchmark construction only — **no model fitting, no architecture changes, no threshold tuning**

---

## 0. Purpose

Stage 22 converts the two Stage-21-qualified prospective datasets into explicit, auditable benchmark tables for the next learnability stage.

The required prospective geometry is:

`X_before + U -> Y_future`

where:

- `X_before` is measured before the perturbation/outcome,
- `U` is the intervention,
- `Y_future` is measured later and is not reconstructed from the same pretreatment RNA vector,
- biological linkage is clone/lineage grounded,
- biological grouping prevents descendants/sister cells from leaking across folds.

Stage 22 does **not** ask whether a model can learn the task. It freezes exactly what the task is.

### 0.1 Pre-start plan freeze

This V2 plan must be committed **before** Stage-22 implementation begins. The plan-freeze commit should contain the plan change only (plus no benchmark outputs or modeling work). Stage-22 manifests then record the frozen plan file hash and reconstruction commit `6c2f2d6`.

Implementation may begin immediately after that dedicated plan-freeze commit, but in a separate change/commit history. Do not edit the frozen plan in response to Stage-22 results; any later correction is additive and explicitly versioned/recorded.

---

## 1. Frozen inputs from Stage 21D

### 1.1 Role A — Rewind / GSE227151

Local root:

```text
D:\GSE227151_Rewind\
```

Primary data:

```text
GSM7092515\
GSM7092516\
filtered10XCells.txt
stepThreeStarcodeShavedReads_BC_10X.txt
stepThreeStarcodeShavedReads_BC_gDNA.txt
GSE227151-GPL18573_series_matrix.txt.gz
GSE227151_family.xml
```

Author code:

```text
D:\GSE227151_Rewind\author_code_zenodo7707418\
```

Load-bearing author scripts:

```text
plotScripts\rewind10X\R1\
  20220921_R1_primedVersusNonPrimedMarkersAndDistribution.R

plotScripts\rewind10X\R1\
  2022.02.14_R1_cellNumberDistributionForPrimedVersusNonPrimed.R
```

Stage-21D reconstruction result:

```text
RECONSTRUCTION_PASS

author-QC cell/sample records      3,921
unique bare cellID strings         3,913
unique clones                      3,149

primed cells                          42
primed clones                         35
nonprimed cells                    3,879
nonprimed clones                   3,114

selected gDNA barcode rows           101
author rule                        slice_max(n=100, with ties)
100th-ranked nUMI cutoff           2,365
```

Important identifier rule:

```text
cell_uid = (SampleNum, cellID)
```

In the current Rewind R1 files there are **0 bare-cellID collisions between the two lanes**. The compound key is still frozen as defensive namespace hygiene because the lanes are separate 10x libraries; it is not evidence that a cross-lane collision was observed.

Important split rule:

```text
outer biological grouping unit = clone
```

because 311 clones span both 10x lanes. A cell-random or lane-random split would leak clone-level outcome information.

Sample mapping remains source-resolved by barcode containment:

```text
SampleNum 1 -> GSM7092516
SampleNum 2 -> GSM7092515
```

Do not infer this mapping from GEO title text.

---

### 1.2 Role B — GSE279162 / WM989

Local root:

```text
D:\GSE279162\
```

Samples:

```text
GSM8562999  Naive1
GSM8563000  Naive2
GSM8563001  Naive3

GSM8563002  Dabrafenib
GSM8563003  Trametinib
GSM8563004  CoCl2
GSM8563005  Acid
GSM8563006  Cisplatin
GSM8563007  Doxorubicin
```

Author code:

```text
D:\GSE279162\author_code_Schaff_manuscript\
```

Primary provenance chain:

```text
How_to_run_code_README.txt
preprocess_GEX.Rmd
preprocess_cDNA_BCs.Rmd
preprocess_gDNA_BCs.Rmd
Find_Markers_Top_Res_lins_in_naive.Rmd
```

Stage-21D reconstruction result:

```text
RECONSTRUCTION_PASS

raw cells                           77,417
post-QC cells                       46,891
assigned lineage cells              42,771
NA lineage cells                     4,120
unique assigned clones               2,215

clones with naive observation        1,401
also seen in >=1 treatment             929
also seen in >=2 treatments             603
seen in all 6 treatments                37
```

Author clone assignment is frozen as:

```text
RNA QC from preprocess_GEX.Rmd
-> Custom features as lineage assay
-> barcode_clustering(cell_lower_limit=100, cor_threshold=0.55)
-> barcode_combine(...)
-> barcoding_posterior(...)
-> barcoding_assignment(difference_val=0.2)
-> assigned_lineage = NA when assigned_posterior < 0.5
```

Do not restore Stage-21D revision-1 dominant-fraction or arbitrary UMI-floor clone calls.

### 1.3 Reuse the Stage-21D reconstruction contract

Stage 22 must consume or directly reuse the **tested Stage-21D revision-2 reconstruction implementation** wherever practical rather than independently re-implementing the author rules a second time.

If code must be refactored for reuse:

```text
author-rule constants
QC thresholds
barcode clustering / posterior assignment
Rewind top-100-with-ties logic
sample mappings
```

must have one source of truth with regression tests against the Stage-21D anchors.

A duplicated implementation that merely happens to reproduce the headline counts is not sufficient if it can silently drift from the source-faithful Stage-21D logic.

### 1.4 Local-path portability

The `D:\\...` paths above describe the current machine only. Benchmark code must accept the dataset roots as explicit CLI arguments (or an equivalent documented environment/config mechanism), e.g.:

```text
--rewind-root <path>
--wm989-root <path>
```

Do not hard-code the current Windows paths into reusable Stage-22 loaders, tests, or benchmark tables. Manifests may record the path used for a run as non-portable provenance, but file identity is established by accession/basename/size/hash.

The Stage-21D diagnostic may keep its existing `D:\\...` module defaults. Its reconstruction functions already accept explicit `base=`/root arguments, so Stage 22 must **always pass the requested roots explicitly** when reusing them. No Stage-21D rewrite is required merely for path portability.

---

## 2. Stage 22 deliverables

Create:

```text
experiments/build_stage22_prospective_benchmarks.py
tests/test_stage22_prospective_benchmarks.py

results/stage22_rewind_benchmark_manifest.json
results/stage22_wm989_benchmark_manifest.json
results/stage22_prospective_benchmark_results.json

results/stage22_rewind_cells.csv
results/stage22_rewind_clones.csv

results/stage22_wm989_cell_assignments.csv
results/stage22_wm989_naive_cells.csv
results/stage22_wm989_clone_treatment.csv
results/stage22_wm989_clones.csv

plans/(newer)practical plans/RECORDs/stage_22_RECORD.md
```

Names may differ slightly if repository conventions require it, but the information content must be equivalent.

Do **not** commit large raw/public matrices. The benchmark outputs should be compact derived tables plus provenance hashes/manifests.

---

# 3. Rewind benchmark

## 3.1 Scientific task

Role A is a fixed-intervention prospective reprogramming task:

```text
pretreatment fibroblast state
+ OSKM/reprogramming intervention
-> future author-defined Rewind priming outcome
```

Operationally, the future outcome is whether a clone barcode is selected by the authors' fixed top-100 iPSC-gDNA abundance rule and therefore labels its linked pretreatment cells as `primed`.

`U` is fixed for the primary benchmark. Rewind therefore tests whether pretreatment molecular state contains prospective information about this later lineage-grounded priming outcome.

Do **not** silently reinterpret `nonprimed` as proven biological death or absolute reprogramming failure. It is the author-defined complement of the top-100 future-gDNA rule and can still be affected by downstream sampling/detection.

Rewind R1 is one true biological replicate. Clone-held-out evaluation within R1 is therefore **within-experiment prospective generalization**, not independent biological-replicate validation. Freeze this limitation in the benchmark manifest and carry it into all Stage-23 claims.

It does **not** test treatment choice or treatment interaction.

---

## 3.2 Rewind cell table

Construct one row per **retained unique pretreatment `cell_uid` after the frozen §3.5 ambiguity exclusion**. The 3,921 Stage-21D rows are source assignment records, not the final cell-table row count.

Required columns:

```text
cell_uid
cellID
SampleNum
gsm
clone_id
nUMI
fracUMI
nLineages

y_primed
outcome_source
outcome_rule
outcome_semantics
biological_replicate
generalization_scope
outer_group
outer_fold
expression_barcode
expression_column_index
expression_source
```

Definitions:

```text
cell_uid = f"{SampleNum}:{cellID}"
clone_id = BC50StarcodeD8
outer_group = clone_id
biological_replicate = "R1"
generalization_scope = "within_R1_clone_heldout"
outcome_semantics = "author_top100_iPSC_gDNA_priming"
y_primed = 1 iff clone_id is selected by the frozen author top-100-with-ties rule
```

`y_primed` must be reproduced only from the author rule:

```text
gDNA arm
-> group by BC50StarcodeD8, SampleNum
-> sum counts
-> slice_max(n=100, with ties)
-> join selected barcodes to filtered10XCells
```

Expected validation:

```text
primed cell/sample records = 42
primed unique clones       = 35
```

No alternative threshold may be substituted.

---

## 3.3 Rewind clone table

Create one row per clone.

Required columns:

```text
clone_id
n_pretreatment_cells
n_lanes
lane_membership
y_primed
outcome_semantics
n_primed_cells
n_nonprimed_cells
outer_group
outer_fold
```

All cells belonging to the same clone must share the same clone-level future outcome.

Assert that no clone receives contradictory `primed` / `nonprimed` labels.

---

## 3.4 Rewind expression linkage

Build a deterministic mapping from `cell_uid` to the corresponding expression column in GSM7092515/GSM7092516.

Do not write a large expression matrix into `results/`.

Instead, freeze the mapping directly in `stage22_rewind_cells.csv` and summarize it in the manifest:

```text
GSM/sample mapping
matrix dimensions
feature hash
barcode hash
cell_uid -> expression barcode
cell_uid -> source column index
cell_uid -> source matrix
```

The benchmark loader used by Stage 23 must be able to reconstruct `X_before` from this mapping without re-deciding any labels.

---

## 3.5 Rewind exclusion / analysis population

The source audit population remains:

```text
3,921 author-QC cell×clone assignment records
3,913 unique cell_uid values
```

There are **8 cell_uid values with two clone assignments each** (`nLineages = 2`), represented by 16 source rows. A single expression cell cannot be assigned to two independently split clone outcomes without creating an ambiguous `outer_fold`.

Therefore the **primary Stage-22 Rewind benchmark pre-registers this exclusion**:

```text
exclude every cell_uid with >1 distinct clone_id
exclusion_reason = "ambiguous_multi_lineage_clone_assignment"
```

All rows for those 8 cell_uids must be removed from the primary benchmark and enumerated in the Rewind manifest/exclusion audit. Do not choose one of the two clones heuristically.

Stage-21D verification established that none of these 8 cells belongs to a primed clone, so the positive anchors remain:

```text
42 primed benchmark cells
35 primed benchmark clones
```

Expected retained unique cell count:

```text
3,913 - 8 = 3,905 cell_uid
```

Recompute and report the post-exclusion negative-cell and unique-clone counts from the data rather than silently carrying the pre-exclusion Stage-21D row counts forward.

Do not apply later UMAP cluster exclusions from figure plotting unless Stage 22 can establish that they were part of the prospective label definition rather than downstream visualization.

If any other exclusion is required, enumerate it separately with its reason and effect on benchmark counts.

---

# 4. GSE279162 / WM989 benchmark

## 4.1 Scientific task

Role B is the treatment-conditioned prospective task:

```text
pretreatment naive state X
+ treatment U
-> future treatment-specific clone response Y
```

This dataset is intended to support later testing of:

```text
U-only
X-only
additive X+U
explicit X×U
```

Stage 22 only constructs the benchmark needed for those comparisons.

---

## 4.2 Freeze the author-assigned cells

Reproduce the author pipeline exactly and create a compact cell-level table.

Required columns:

```text
cell_uid
cell_barcode
sample
gsm
condition
assigned_lineage
assigned_posterior
is_assigned
is_naive
```

Define:

```text
cell_uid = f"{gsm}:{cell_barcode}"
```

Do not assume a 10x barcode string is globally unique across samples. In the current WM989 data, **722 bare 10x barcode strings occur in more than one sample**, so the `gsm:barcode` compound key is required here rather than merely defensive.

Use the post-QC population from the source pipeline. Reproduce lineage clustering/posterior assignment **once on the joint post-QC all-sample object**, as in the author pipeline; do not independently refit barcode clustering or posterior parameters per sample.

For Stage 22 benchmark construction, cells with:

```text
assigned_lineage = NA
```

must remain documented but are excluded from clone-linked prospective rows.

Do not assign them using an alternative heuristic.

---

## 4.3 Naive pretreatment cell table

Create one row per assigned naive cell from Naive1/Naive2/Naive3.

Required columns:

```text
cell_uid
source_naive_sample
gsm
clone_id
assigned_posterior
expression_barcode
expression_column_index
expression_source
outer_group
outer_fold
```

Definitions:

```text
clone_id = assigned_lineage
outer_group = clone_id
```

The three naive libraries constitute pretreatment observations, not three treatment classes.

Pool them only at the scientific `condition = naive` level while preserving source-sample provenance.

### 4.3.1 WM989 expression linkage

For every assigned Naive1/Naive2/Naive3 benchmark cell, freeze a deterministic mapping to its pretreatment expression column:

```text
cell_uid
-> GSM/sample
-> 10x barcode
-> source matrix
-> source column index
```

The WM989 manifest must also freeze the naive expression matrix dimensions plus feature/barcode hashes for all three naive samples.

Stage 23 must be able to reconstruct `X_before` without rerunning lineage assignment or re-deciding which cells belong to a clone.

### 4.3.2 WM989 pretreatment clone table

Create one structural row per clone with at least one assigned naive pretreatment cell.

Required columns:

```text
clone_id
n_naive_cells
n_naive1_cells
n_naive2_cells
n_naive3_cells
naive1_total_assigned_cells
naive2_total_assigned_cells
naive3_total_assigned_cells
naive_pooled_fraction
naive_source_samples
outer_group
outer_fold
```

where:

```text
naive_pooled_fraction =
    n_naive_cells /
    (naive1_total_assigned_cells + naive2_total_assigned_cells + naive3_total_assigned_cells)
```

This table is **pretreatment structure only**. Do not place future treatment outcomes into it as candidate model features. Treatment-specific outcomes remain in `stage22_wm989_clone_treatment.csv`.

Preserve the per-naive-sample counts and denominators rather than only the pooled total so Stage 23 can distinguish captured baseline clone abundance from naive-library depth effects.

---

## 4.4 Clone × treatment outcome table

Create the full prospective matrix for every clone that has a naive pretreatment observation.

Required columns:

```text
clone_id
treatment
n_naive_cells
naive_pooled_fraction
n_post_cells
post_fraction
post_rank
post_rank_fraction
post_tie_size
detected_post
outcome_observation_status
treatment_sample_available
treatment_total_assigned_cells
naive_source_samples
outer_group
outer_fold
```

Treatments:

```text
dabrafenib
trametinib
CoCl2
acid
cisplatin
doxorubicin
```

Definitions:

### Raw outcome

```text
y_count = n_post_cells
```

This is the primary author-supported future response quantity.

### Treatment-normalized abundance

```text
y_fraction = n_post_cells / treatment_total_assigned_cells
```

Freeze it as a derived representation because sequencing/cell recovery depth differs by treatment.

`post_fraction` uses the **full treatment assigned-cell total** as its denominator. Fractions over *all* assigned clones in a treatment sum to one, but the benchmark contains only the 1,401 naive-observed eligible clones, so the benchmark rows generally sum to **less than one**: they sum to the share of treatment-assigned cells attributable to that eligible subset.

The dependence warning still applies because these fractions share a treatment-level denominator and are components of the larger treatment composition. Preserve this limitation explicitly; Stage 23 must not count the 8,406 clone×treatment rows as 8,406 independent biological units.

### Rank

Use one deterministic convention:

```text
post_rank = descending competition rank of n_post_cells
            among the eligible naive-observed benchmark clones only
            (equivalent to pandas rank(method="min", ascending=False))

post_rank_fraction = post_rank / number_of_eligible_naive_clones

post_tie_size = number of clones sharing that n_post_cells value
```

Lower `post_rank` / `post_rank_fraction` means greater observed post-treatment abundance. Zero-count clones therefore tie rather than receiving arbitrary order.

Freeze rank as a representation, not as evidence that ranking is learnable.

---

## 4.5 Zero outcomes must be explicit

For each clone with a valid naive observation, materialize **all six treatment rows**, including treatments in which the clone has zero assigned post-treatment cells.

With the Stage-21D anchor of 1,401 naive-observed clones, the expected structural row count is:

```text
1,401 × 6 = 8,406 clone×treatment rows
```

Treat this as a regression assertion, not as a tuning target.

Every one of the six treatment samples exists, so a zero clone count is an **observed zero count**, not missing outcome data.

Therefore:

```text
n_post_cells = 0
detected_post = false
outcome_observation_status = "observed_zero"
treatment_sample_available = true
```

must be represented explicitly rather than dropping the row.

For `n_post_cells > 0`:

```text
detected_post = true
outcome_observation_status = "observed_nonzero"
treatment_sample_available = true
```

Only a genuinely unavailable treatment sample would use `NA` plus `treatment_sample_available = false`.

This prevents the benchmark from silently becoming “among detected surviving clones, predict abundance.”

Stage 22 must **not relabel `observed_zero` as death, failure, sensitivity, or non-resistance**. A zero can reflect biological non-survival, allocation/sampling, sequencing, or capture limits. The primary table preserves the observed count; biological interpretation of the zero mass remains a declared measurement/censoring limitation for Stage 23.

Report zero rates stratified by `n_naive_cells` as a diagnostic so later work can see whether sparse pretreatment sampling strongly predicts apparent zeros; do not use that diagnostic to alter the benchmark after the fact.

---

## 4.6 Do not freeze a binary resistance threshold in Stage 22

Do not use:

```text
top 5
top 10
top 10%
present / absent at arbitrary count floor
```

as the primary benchmark target.

The source-supported primitive quantity is treatment-specific clone abundance.

Freeze the continuous/count/rank representations first.

If Stage 23 later needs a binary endpoint for a pre-registered diagnostic, it must define that endpoint before fitting and justify it independently of model performance.

## 4.7 Freeze membership, not a post-hoc X representation

Stage 22 freezes which pretreatment cells belong to each clone. It does **not** choose a predictive representation after looking at model results.

For both datasets, the benchmark must permit reconstruction of:

```text
clone_id -> ordered/set membership of pretreatment cell_uids
```

Stage 23 must pre-register **before any model fitting** how `X_before` is represented when a clone has multiple pretreatment cells (for example, clone pseudobulk, an explicitly defined cell-bag encoder, or a cell-level learner with clone-balanced weighting). It may not try several aggregation schemes and keep the one with the best held-out result.

Regardless of representation, clone-level outcomes remain the independent evaluation units.

## 4.8 Feature-eligibility firewall

Stage 22 must classify benchmark columns so Stage 23 cannot accidentally turn provenance or clone-size shortcuts into the claimed transcriptomic predictor.

Freeze three classes:

```text
TARGET:
  y_primed / n_post_cells / post_fraction / post_rank...

PROVENANCE-ONLY:
  cell_uid, clone_id, GSM, SampleNum/sample, source paths,
  expression column indices, outcome_source/outcome_rule,
  outer_group, outer_fold

BASELINE / NUISANCE:
  Rewind captured n_pretreatment_cells and lane counts
  WM989 n_naive_cells and per-naive-sample clone counts
  technical/QC quantities such as nUMI, fracUMI, nLineages,
  assigned_posterior unless explicitly justified otherwise

PRIMARY X:
  pretreatment gene-expression values only
```

The primary transcriptomic claim may not use provenance-only columns.

Clone abundance is especially important because the response is a future clone-abundance/detection outcome. Stage 23 must include abundance-only nuisance baselines and must show that transcriptomic `X` adds value beyond them before claiming state-conditioned predictivity.

Minimum nuisance comparisons to pre-register for Stage 23:

```text
Rewind:
  prevalence-only
  captured pretreatment clone-size baseline
  expression X
  expression X + clone-size nuisance

WM989:
  U-only
  baseline naive clone abundance + U
  expression X
  expression X + U
  expression X + baseline abundance + U
  explicit X×U model with the same baseline abundance control
```

These are comparison requirements, not Stage-22 model fits.

---

# 5. Biological grouping and leakage prevention

## 5.1 Freeze outer folds during Stage 22

Do not leave the primary outer split to Stage 23. Stage 22 must create and commit deterministic clone-level `outer_fold` assignments so model development cannot shop among favorable splits.

Use five outer folds unless an existing repository-wide frozen convention already supersedes this plan.

If no such convention exists:

```text
STAGE22_SPLIT_SEED = 22022
```

### Rewind

Apply the frozen §3.5 ambiguous-cell exclusion **before** constructing the benchmark clone table and outer folds.

Create the split once on the resulting **clone table**, using:

```text
StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=22022
)
```

stratified by clone-level `y_primed`.

Propagate each clone's `outer_fold` to every linked cell.

Assert every fold contains both classes and report positive clones per fold. With only 35 positive clones, this check is mandatory.

Forbidden:

```text
cell-random split
lane-random split
split by GSM alone
repeated reseeding until performance looks good
```

The two lanes are technical/sample provenance, not independent biological replicates.

### WM989

Create the split once on the sorted unique naive-observed clone table, using:

```text
KFold(
    n_splits=5,
    shuffle=True,
    random_state=22022
)
```

Propagate that clone fold to:

```text
all naive cells for clone C
all six treatment rows for clone C
```

A clone appearing under multiple treatments cannot be split so that one treatment arm is in train and another treatment arm for the same clone is in test in the primary state-generalization benchmark.

Held-out-treatment / unseen-perturbation generalization is reserved for **Stage 26** under the roadmap. Stage 23 should use the frozen clone-held-out state-generalization folds and must not introduce a second treatment-holdout geometry unless the roadmap is explicitly revised first.

---

# 6. Benchmark grains

Stage 22 must freeze two grains rather than forcing one table to serve every later experiment.

## 6.1 Cell-grain pretreatment table

Used to obtain `X_before`.

```text
cell -> clone
```

There may be several pretreatment cells per clone.

## 6.2 Clone-grain outcome table

Used for `Y_future`.

```text
clone -> outcome
```

For WM989:

```text
clone × treatment -> outcome
```

The benchmark must make the aggregation relationship explicit.

Do not pretend multiple cells from the same clone are independent outcome observations. If Stage 23 trains at cell grain, it must pre-register how clone multiplicity is weighted so clones with many captured pretreatment cells do not silently contribute more independent outcome weight.

---

# 7. Stage 22 benchmark statistics to report

## 7.1 Rewind

Report:

```text
source author-QC assignment records
source unique cell_uid
source unique clones

ambiguous multi-lineage cell_uid excluded
source rows removed by that exclusion
retained benchmark unique cell_uid
retained benchmark unique clones

bare cellID cross-lane collision count (expected 0)
cells per clone distribution after exclusion
clones spanning both lanes after exclusion
positive cells
negative cells
positive clones
negative clones
class prevalence at cell grain
class prevalence at clone grain
positive/negative clones per outer fold
genes/features in expression matrix
mapped/unmapped benchmark cells
```

Expected source-reproduction anchors include:

```text
42 primed cell/sample records
35 primed clones
```

---

## 7.2 WM989

Report:

```text
post-QC assigned cells by condition
NA lineage cells by condition
unique cell_uid
duplicate bare 10x barcode strings across samples
unique assigned clones
clones with naive observations
naive cells per clone distribution
naive1/2/3 total assigned-cell denominators
naive_pooled_fraction distribution
naive expression mapped/unmapped cells

for each treatment:
  total assigned cells
  clones with nonzero post abundance
  clone abundance distribution
  fraction of naive clones with zero post abundance

clone coverage:
  seen in >=1 treatment
  seen in >=2 treatments
  ...
  seen in all 6 treatments

number of complete clone×6-treatment rows (expected 8,406)
clone counts per outer fold
```

The Stage-21D source-faithful anchors include:

```text
clones with naive observation = 1,401
seen in >=2 treatments        = 603
seen in all 6 treatments      = 37
```

These are regression checks, not tuning targets.

---

# 8. Missingness semantics

Every derived field that can be absent must distinguish:

```text
0
NA / unavailable
not applicable
```

Examples:

- `n_post_cells = 0` means the treatment sample was observed but zero assigned cells for that clone were detected.
- genuinely missing treatment sample data would be `NA`, not zero.
- a Rewind clone has no treatment-varying `U`, so treatment-comparison fields are not applicable.

Do not silently coerce missingness to zero.

---

# 9. Provenance requirements

Each benchmark manifest must include:

```text
source accession
local source path used for this run (non-portable provenance only)
source file names
source file sizes
SHA-256 hashes
author-code file hashes
author rule identifiers
reconstruction commit = 6c2f2d6
plan version / plan file hash
derived table/artifact SHA-256 hashes
```

Do **not** require any file to contain its own hash or the git commit that contains itself; both are circular provenance.

Dataset manifests may hash their derived CSV/table artifacts but must not include their own manifest hash. If a top-level results file hashes the dataset manifests, that top-level file must likewise omit its own hash.

`stage_22_RECORD.md` may record the **implementation commit immediately preceding the record-finalization commit**, or leave final commit identity to git history/the completion report. Do not create an endless commit-self-reference cycle merely to place the final SHA inside the file that is part of that SHA.

Keep committed result/manifests deterministic: do not inject wall-clock timestamps or the mutable repository `HEAD` into content that is expected to reproduce byte-for-byte.

Record the git `HEAD` observed before benchmark generation, the implementation commit, and any run timestamp in `stage_22_RECORD.md` / the completion report instead of the hashed benchmark manifests. If stronger implementation identity is needed inside a manifest, use a deterministic hash of the builder source file rather than mutable repository state.

The source/public data remain outside git.

Derived compact benchmark tables may be committed if repository policy permits.

---

# 10. Frozen benchmark contracts

Add tests that pin at least the following.

## Rewind

```text
source audit reproduces 3,921 assignment rows and 3,913 unique cell_uid
exactly 8 cell_uid have >1 distinct clone_id and are enumerated/excluded
those 8 correspond to 16 source assignment rows
retained benchmark has 3,905 unique cell_uid
retained benchmark cell_uid is unique
cell_uid uses SampleNum + cellID
bare cellID cross-lane collision count is 0 in the current source
no retained cell maps to >1 clone
no retained clone has contradictory outcome
42 primed benchmark cells
35 primed benchmark clones
author top-100-with-ties rule is used
all retained cells from a clone map to one outer_group and one outer_fold
five frozen clone-level outer folds exist
every Rewind fold contains positive and negative clones
Rewind manifest states biological_replicate=R1 and within_R1 generalization scope
nonprimed is stored as an author-defined operational label, not asserted as proven failure
feature-eligibility classes are present
```

## WM989

```text
author QC is applied
Custom features feed the lineage assay
barcode_clustering uses 100 / 0.55
barcoding_assignment uses 0.2
posterior <0.5 becomes NA
no dominant-fraction fallback exists
naive1/2/3 map to pretreatment condition
cell_uid = gsm + 10x barcode and is globally unique
current source has 722 bare 10x barcode strings reused across samples
all naive benchmark cells map to a source expression column or exclusions are explicit
lineage clustering/posterior assignment is run jointly, not refit per sample
every eligible naive clone receives six clone×treatment rows
1,401 eligible clones produce exactly 8,406 structural clone×treatment rows
all six treatment samples are available and have positive assigned-cell denominators
zero post outcomes are retained as observed_zero, not relabeled as failure
rank ties use the frozen competition-rank convention
no Stage-22 binary resistance threshold exists
all rows for a clone share one outer_group and one outer_fold
five frozen clone-level outer folds exist
feature-eligibility classes are present
```

---

# 11. Stage 22 acceptance gates

Stage 22 passes only if all of the following hold.

### G22-1 — Prospective linkage

Every benchmark outcome can be traced from:

```text
pretreatment observation
-> biological clone
-> later independently measured outcome
```

### G22-2 — No label leakage

No target is computed from pretreatment expression values.

### G22-3 — Unique biological grouping

Every benchmark row has an unambiguous `outer_group`.

### G22-4 — Expression resolvability

Every retained pretreatment record used by the benchmark maps to exactly one expression column and exactly one clone/outer fold, or an exclusion is explicitly enumerated and justified. The 8 frozen Rewind multi-lineage cells are the pre-registered ambiguity exclusion, not an expression-mapping failure.

### G22-5 — Outcome completeness

Rewind:

```text
every included clone has a frozen primed/nonprimed outcome
```

WM989:

```text
every eligible naive clone has six explicit treatment rows
```

including zeros.

### G22-6 — Author-rule fidelity

The benchmark reproduces Stage-21D's source-faithful reconstruction rules without heuristic substitution.

### G22-7 — Frozen evaluation geometry

Primary clone-held-out outer folds are materialized during Stage 22 and pass the grouping assertions. Rewind fold construction must preserve positive/negative representation.

### G22-8 — Feature firewall

Target, provenance-only, baseline/nuisance, and primary-expression columns are explicitly classified.

### G22-9 — Claim scope

Rewind's single-biological-replicate limitation and WM989's observed-zero/compositional outcome limitations are present in the manifests/results.

### G22-10 — No modelling

No predictive model is fitted during Stage 22.

---

# 12. Stage 22 verdicts

The diagnostic must compute one of:

```text
BENCHMARK_READY
BENCHMARK_READY_WITH_DECLARED_MISSINGNESS
BENCHMARK_BLOCKED_LINKAGE
BENCHMARK_BLOCKED_EXPRESSION_MAPPING
BENCHMARK_BLOCKED_OUTCOME
BENCHMARK_BLOCKED_LEAKAGE
```

Per-dataset verdicts must be reported separately.

Overall Stage 22:

```text
STAGE_23_READY
```

only if Role A has one of:

```text
BENCHMARK_READY
BENCHMARK_READY_WITH_DECLARED_MISSINGNESS
```

Role B remains high-value but non-blocking according to the Stage-21 V3 rule. If Role B fails benchmark construction, Stage 23 may open on Role A with the treatment-interaction contribution scoped down. The overall result must still preserve Role B's failure/limitation explicitly rather than hiding it behind the non-blocking gate.

Do not hard-code the overall verdict; derive it from per-dataset results.

---

# 13. What Stage 22 does NOT establish

A `BENCHMARK_READY` verdict does **not** mean:

```text
X predicts Y
X adds beyond U
X+U beats X-only or U-only
an X×U interaction exists
treatment ranking is learnable
the model generalizes to unseen treatments
the model is calibrated
the outcome is free from sampling/censoring bias
```

Those are Stage 23+ questions.

Stage 22 establishes only that the prospective benchmark is valid enough to ask them.

In particular, a high WM989 `observed_zero` rate or a strong relationship between captured naive clone abundance and future detection/abundance is **not a Stage-22 failure**. It is exactly why the Stage-23 abundance-only nuisance baselines are mandatory. A scientifically valid Stage-23 result may be that transcriptomic `X` adds little or nothing beyond captured clone size.

---

# 14. Stage 23 handoff contract

If Stage 22 passes, Stage 23 receives frozen tables and must not redefine the benchmark after seeing model performance.

The minimum Stage-23 comparison set remains:

```text
Role A:
  intercept / prevalence
  X-only

Role B:
  U-only
  X-only
  additive X+U
  explicit X×U
```

with clone-held-out evaluation.

Before fitting any Stage-23 model, freeze in the Stage-23 plan:

```text
X representation / clone aggregation rule
expression normalization and feature-selection procedure
clone multiplicity weighting
primary metric(s)
use of baseline clone-abundance nuisance variables
treatment-count normalization choice
handling/interpretation of observed_zero outcomes
```

For Rewind, the metric set must include a class-imbalance-sensitive metric such as **average precision / PR-AUC** at clone grain. Accuracy alone is not an acceptable primary metric for 35 positive versus 3,114 negative clones.

> *Annotation, 2026-08-21 (pre-start, plan text above left as frozen): the illustrative figure
> `3,114` is the PRE-exclusion negative-clone count carried from Stage 21D. After the §3.5
> ambiguity exclusion the benchmark actually holds **35 positive versus 3,112 negative clones**
> (3,147 retained clones), because 2 clones consisted only of the 8 excluded multi-lineage cells.
> The argument for PR-AUC is unchanged. This is an annotation, not a reopening — §3.5's
> instruction to recompute post-exclusion counts from the data remains the binding rule, and the
> executed numbers belong in `stage_22_RECORD.md`.*

These choices may not be selected after comparing held-out performance.

Any learned expression preprocessing (feature selection, scaling parameters, dimensionality reduction, supervised filtering) must be fit on the training portion of each outer fold only. Do not compute a global expression representation using held-out cells and then cross-validate downstream.

Stage 23 is the learnability/interaction gate. It is not part of this plan.

---

# 15. Engineering constraints

```text
src/ unchanged
no training
no architecture changes
no hyperparameter search
no label threshold tuning
no public-data search
no raw public-data commits
no hard-coded D:\\ dataset roots in reusable code
no nondeterministic timestamps in benchmark identity/manifests
no rewriting Stage 21A/21B/21C/21D history
no cell-random outer split
no lane-random Rewind split
no outer-split reseeding after model results are seen
no use of provenance-only columns as predictive features
no divergent duplicate implementation of Stage-21D author rules
```

Any correction discovered during Stage 22 must be recorded additively.

---

# 16. Required completion sequence

```text
1. build Rewind benchmark
2. freeze Rewind clone-level outer folds
3. run Rewind integrity assertions
4. build WM989 benchmark
5. freeze WM989 clone-level outer folds
6. run WM989 integrity assertions
7. write feature-eligibility metadata
8. write manifests/results
9. compute per-dataset verdicts
10. compute overall Stage-23 gate
11. create stage_22_RECORD.md
12. run full pytest suite
13. run CI-scope ruff
14. commit/push
15. STOP for review before Stage 23
```

---

# 17. Completion report format

When Stage 22 is complete, report:

```text
commit
tests
ruff

Rewind:
  verdict
  source rows / source unique cell_uid / source clones
  ambiguous cell_uid excluded / source rows removed
  retained benchmark cells/clones
  positive cells/clones
  positive clones per outer fold
  bare-cellID cross-lane collisions
  expression mapping rate
  outer split unit
  generalization scope

WM989:
  verdict
  naive cells/clones
  unique cell_uid / reused bare barcodes across samples
  clone×treatment rows
  zero-outcome rate by treatment
  zero-outcome rate stratified by n_naive_cells
  clones in >=2 treatments
  clones in all 6
  clones per outer fold
  expression mapping rate
  outer split unit

overall:
  Stage-23 gate
```

Do not start Stage 23 in the same change.
