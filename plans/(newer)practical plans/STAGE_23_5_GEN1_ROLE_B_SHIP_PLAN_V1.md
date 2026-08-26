# STAGE 23.5 — GEN-1 ROLE-B CLAIM REVISION, RANKING PREREGISTRATION, AND STAGE-24 OPENING

**Status:** DRAFT, PRE-EXECUTION. No ranking statistic may be computed and no Stage-24 model/tool work may begin until this plan is independently audited, committed, and recorded by canonical-LF SHA-256.

**Plan version:** V1

**Date opened:** 2026-08-26

**Decision owner:** the user

**Repository:** `hagai-tyty/Cellular-Outcomes-of-Perturbations`

**Repository evidence base:** commit `c52c7ac9a73abcf121624112d52f4f463db06f7a`

**Depends on:** Stage 22 benchmark readiness; final Stage-23 synthesis; Stage-23.2H confirmation record; `STAGE_23_2_HANDOFF_TO_STAGE_24.md`; `STAGE_24_DECISION_BRIEF.md`

**Pre-freeze audit:** performed 2026-08-26 against repository commit
`c52c7ac9a73abcf121624112d52f4f463db06f7a`. Every Stage-22/23/23.2 quantity quoted below was
re-derived from its source artifact rather than accepted from the draft; the verification is
recorded in `RECORDs/stage_23_5_RECORD.md`. The audit found no incorrect number and four
specification defects, all repaired in place before any statistic was computed:

```text
1  the section-8.7 permutation compute budget was unstated      -> 0.2, with the measured precedent
2  the reproduction tolerance was named but never defined       -> 7.1, a number and a stop rule
3  criterion 5 used an unpowered point estimate as a veto while
   simultaneously disclaiming significance testing              -> 8.8 / 8.10, restated as a
                                                                   directional-consistency check
4  the bootstrap CI and the permutation null quantify different
   things and were presented as equivalent                      -> 8.6, asymmetry declared
```

Two clarifications were added at the same time: the fold provenance in 2.2, and the deliberate
direction difference between `delta_RANK` and `delta_TOP1` in 8.3. The six decision criteria and
the forbidden list were hoisted into section 0 so that a reader who stops early still sees them.
No number, threshold, metric, population, model or verdict rule was loosened; item 3 is the only
change to a criterion, and it repairs an internal contradiction rather than lowering a bar.

**Scope:** make the minimum explicit claim-architecture and completion-path revisions needed to ship a scientifically useful Generation-1 paper/tool using only the datasets already in the project; freeze the one new load-bearing ranking analysis before it is inspected; hand a bounded, auditable contract to Stage 24.

---

# 0. Executive decision

Stage 23.5 selects the following Generation-1 path:

```text
ROLE_B_PRIMARY_GEN1
```

Generation 1 will be built around the already-supported WM989 Role-B result:

```text
pretreatment molecular state X
+ candidate known treatment U
        -> independently observed future clonal outcome Y
```

The primary scientific claim is treatment-conditioned prospective prediction within the existing six-condition WM989 system. Rewind Role A remains positive-but-underpowered supporting evidence and is not promoted to a confirmed anchor.

This decision is made because Role B has strong direct evidence for `X × U`, is far separated from its full-refit permutation null, and directly supports the tool being built. It is **not** made because Stage 27 distinguishes Option A from Option B; Stage 27 applies to either route and is revised separately below.

Two roadmap amendments are explicit and must never be collapsed into one:

```text
AMENDMENT 1 — STAGE-24 OPENING RULE
  Stage 24 may open under ROLE_B_PRIMARY_GEN1 after this Stage-23.5 plan,
  audit, digest, and handoff are complete. This is an explicit claim-
  architecture revision. It is not a Stage-23 PASS and it does not change
  any historical Role-A verdict.

AMENDMENT 2 — GEN-1 COMPLETION PATH
  Independent new-system replication is moved from a Generation-1
  publication gate to Generation 2. Stage 26 receives an explicit
  known-treatment-only scoped limit; Stage 27 external replication and
  Stage 28 broad calibration/OOD claims are future work rather than
  prerequisites for the first paper.
```

No additional dataset will be searched, downloaded, qualified, or used for Generation 1.

Generation 1 now has one mandatory completion route:

```text
freeze this plan
  -> build the frozen W5 predictor/tool
  -> execute the preregistered ranking test exactly once
  -> record RANKING_SUPPORTED or RANKING_NOT_SUPPORTED
  -> lock evidence and claims
  -> ship the tool, reproducibility package and manuscript
```

The ranking result changes only whether a validated ranking claim is allowed. It does not decide whether Generation 1 ships. No ranking outcome may reopen dataset search, model selection, metric selection, endpoint selection, an earlier stage, or the roadmap.

## 0.1 The decision criteria, stated up front

The full specification is section 8. It is repeated here because a reader who stops before section 8
must still know what the test is and what cannot be done to it.

```text
STAGE_25_RANKING_SUPPORTED requires ALL of:
  1. observed delta_RANK > 0
  2. lower CI95 bootstrap endpoint > 0
  3. observed delta_RANK > null p95
  4. p_perm <= 0.05 under the 1,000-draw full-refit null
  5. delta_TOP1 >= 0                      directional-consistency check, see 8.8
  6. all integrity, leakage and determinism checks pass

FORBIDDEN after freeze, without exception:
  changing the primary metric, the clone eligibility rule, the weighting,
  the comparator, the null construction, the endpoint, the model,
  the tie handling, the treatment set, or the permutation count;
  early stopping; using C2 or a per-treatment result to rescue C1;
  adding a dataset; any post-result plan revision.
```

Either verdict ships. The full anti-rescue list is section 8.11.

## 0.2 Compute budget — agree this before freezing

The section-8.7 null is a **full-refit** permutation of the WM989 pipeline, and it is expensive.
The cost is stated here rather than discovered during execution:

```text
measured precedent   Stage 23E, WM989 C1, 200 full-refit permutations   335.03 min
                     -> approximately 100 s per draw

this plan            1,000 draws, W4 and W5 refit per draw
                     ~28 h single process
                     ~19-20 h across 3 shards, at the 1.44x effective speedup
                     measured for this pipeline in Stage 23.2H
```

Section 8.7 forbids early stopping, so this cost cannot be trimmed once started. It must be
accepted before the plan is frozen.

**Sharding is required, and each shard MUST write its own cache file.** Stage 23.2H lost a
completed draw to a race when three shards appended to one shared file; the loss was silent and was
caught only by a completeness check. Draw `b` must be seeded from a function of `b` alone so that
shards may be split, interrupted and resumed without changing any value, and the merge step must
assert that every index in `0..n_perm-1` is present before computing a statistic.

---

# 1. Why Stage 23.5 exists

Stage 23 closed with `STAGE_24_BLOCKED_ROLE_A`. Stage 23.2H later produced positive evidence on two untouched Rewind biological replicates but failed its frozen design-power gate. The exact exit remains:

```text
ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
```

Therefore the existing roadmap does not authorize Stage 24. Silently treating the Role-A result as confirmed would be invalid. Silently substituting Role B would also be invalid.

At the same time, Stage 23 established a substantially stronger Role-B result:

```text
ROLE_B_ADDITIVE_PASS
INTERACTION_PASS_MULTI_TREATMENT
C2_INTERACTION_SECONDARY_CONFIRMED
```

The project owner has now imposed two scope constraints:

```text
1. Generation 1 will not use another dataset.
2. Generation 1 must ship a useful positive contribution rather than remain
   indefinitely blocked by attempts to satisfy the original larger claim.
```

Stage 23.5 is the explicit bridge required by those facts. It is a claim-revision and preregistration stage, not a rescue analysis, dataset search, model-development stage, or permission to rewrite prior results.

---

# 2. Immutable historical record

Nothing in Stage 23.5 edits, rescales, re-scores, or reinterprets a historical gate.

## 2.1 Role A — Rewind / GSE227151

Historical Stage-23 result:

```text
verdict                 ROLE_A_SIGNAL_FAIL
delta_AP_state          +0.01050162935116511
p_perm                  0.0845771144278607
null p95                0.014545491038048794
```

Stage-23.2H confirmation result:

```text
cohort                   2,310 clones; 75 positive clones
biological replicates    replicate 2 + replicate 3, ungated samples only
delta_AP                 +0.030504
primary null p95         0.028844
p_perm                   0.03998 at 2,000 permutations
direction                positive in both replicates
transfer                 positive in both directions
frozen power gate        FAIL at recorded 0.64 versus required 0.80
audited power estimate   approximately 0.45; reported, not re-gating
exit                     ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
```

Allowed Role-A wording:

> Pretreatment transcriptional state showed positive, underpowered supporting evidence for the Rewind priming outcome across two independent biological replicates not used to design the correction.

Forbidden Role-A wording:

```text
Role A is confirmed.
The observed delta_AP is an unbiased effect-size estimate.
One replicate carries the signal.
All Rewind replicates share one uniform outcome operationalisation.
```

All usable GSE227151 biological replicates have now influenced diagnosis, confirmation, or the decision record. The Rewind evidence well is closed. No GSE227151 material may be represented as untouched Generation-2 replication.

## 2.2 Role B — WM989 / GSE279162

Frozen cohort:

```text
pretreatment-observed clones       1,401
clone × treatment rows             8,406
treatments                              6
observed-zero rows                 6,150
nonzero rows                       2,256
outer folds                            5   frozen in Stage 22, never re-drawn
outer grouping                     clone
```

The five outer folds, the clone grouping, the outcome construction and every exclusion are the
Stage-22 originals loaded from `results/stage22_wm989_clones.csv` and
`results/stage22_wm989_clone_treatment.csv`. They are not recomputed and not re-seeded.

Primary C1 results:

```text
model                              pooled log loss
W1  B + U                          0.47832325762628236
W4  X + B + U                      0.46895102453285753
W5  X + B + U + X×U                0.45465103771852750

W5 versus W4 interaction gain      0.014299986814330035
bootstrap lower endpoint           0.010103105206688081
permutation null p95               0.0008383045453672399
null draws >= observed             0 / 200
finite-sample p_perm               1 / 201 = 0.004975124378109453

W5 versus W1 full-state gain       0.023672219907754866
bootstrap lower endpoint           0.017499162318074447
permutation null p95               0.0003983520753238711
null draws >= observed             0 / 200
finite-sample p_perm               1 / 201 = 0.004975124378109453
```

The permutation result must be reported as `0/200` exceedances with empirical `p_perm = 1/201`, the minimum attainable value under that run. The observed/null-p95 ratios of approximately 17 and 59 describe separation from the tested null; they are not biological effect-size multipliers and must not be presented as such.

Secondary C2 results:

```text
model                              clone-balanced MAE
W1                                0.5144863827773782
W4                                0.5100798502288132
W5                                0.5007042703677924

W5 versus W4 interaction gain     0.009375579861020777
W5 versus W1 full-state gain      0.013782112409585823
both permutation tests            0/200 exceedances; p_perm = 1/201
```

C2 is endpoint corroboration. It is not independent biological replication.

Treatment-level C1 interaction results:

```text
Acid          +0.0109716227542489     meaningful positive
Cisplatin     +0.0000245928452981     directionally positive, negligible
CoCl2         +0.0189409151730353     meaningful positive
Dabrafenib    +0.0259209299497760     meaningful positive
Doxorubicin   -0.0033224048416363     negative
Trametinib    +0.0332642650052581     meaningful positive
```

Required summary:

> Six conditions were evaluated; four carry meaningful positive interaction, Cisplatin is negligible on C1, and Doxorubicin is negative on both C1 and C2.

Forbidden summary:

> The interaction is supported uniformly across all six treatments.

Captured pretreatment clone abundance remains the dominant predictor:

```text
abundance gain / full state gain = 3.446751843310516
```

This limitation must appear beside, not far away from, the Role-B headline.

## 2.3 Biological-replication boundary

No independent biological replication has been performed for the primary Role-B finding. Clone-held-out folds, two endpoint families, bootstrap resampling, and full-refit permutations do not become biological replication by repetition or wording.

Generation-1 allowed wording:

> The finding is internally supported under clone-held-out evaluation and two endpoint families in one lineage-traced experimental system; independent biological replication remains future work.

---

# 3. Revised Generation-1 claim architecture

## 3.1 Primary claim

The primary Generation-1 claim is:

> Within the existing multi-condition WM989 lineage system, pretreatment Gene Expression contains treatment-specific information about future clonal detection beyond treatment identity and captured pretreatment clone abundance under frozen clone-held-out evaluation.

The claim is deliberately system-bounded and known-treatment-bounded.

## 3.2 Ranking claim and mandatory reporting rule

The preregistered ranking test in §8 is run once and reported regardless of direction.

When every frozen criterion is met, Generation 1 may additionally claim:

> A frozen state×treatment model improves clone-specific ordering of the six observed experimental conditions over a non-interactive additive model.

When any frozen criterion is not met, the manuscript must instead state:

> Clone-specific ordering of the six observed experimental conditions was not supported under the preregistered ranking test.

Both outcomes proceed directly to the same mandatory Generation-1 shipment. No second ranking analysis is authorized.

## 3.3 Supporting claim

Rewind may support only this statement:

> A separately reconstructed reprogramming system showed positive but underpowered evidence that pretreatment transcriptional state carries prospective information about a later lineage outcome.

Rewind cannot validate the WM989 treatment-ranking claim because it does not provide the same multi-treatment task or outcome.

## 3.4 Outcome semantics

WM989 `Y` is an experimentally observed future clonal abundance/detection proxy reconstructed from post-treatment assigned-cell counts. It is not a noise-free measurement of death, sensitivity, resistance, clinical response, or patient benefit.

An observed zero means:

```text
no assigned post-treatment cell was observed for that clone-condition row
```

It must not be silently relabelled as proven death or sensitivity. Consequently the tool uses the neutral labels `future_detection_score`, `future_abundance_score`, and `low-persistence condition` rather than clinical response language.

## 3.5 Claims forbidden in Generation 1

```text
unseen-treatment generalization
cross-cell-line or cross-patient generalization
clinical treatment recommendation
causal treatment-effect estimation
calibrated probability unless calibration is separately frozen and passed
independent biological replication of Role B
uniform benefit across all six conditions
confirmed Role-A prediction
single-cell input equivalence when the model was trained on clone pseudobulk
```

---

# 4. Dataset roles and no-new-data firewall

## 4.1 Primary development/evaluation data

```text
GSE279162 / WM989    Role-B benchmark, model, tool and ranking analysis
```

All Stage-22 clone assignments, outcome construction, five outer folds, feature rules, treatment aliases, and exclusions remain frozen.

## 4.2 Supporting evidence

```text
GSE227151 / Rewind   historical Role-A evidence and limitations only
```

No additional Rewind model fitting is authorized by Stage 23.5.

## 4.3 Prohibited Generation-1 data expansion

Without a separately approved roadmap version, the following are forbidden:

```text
searching for another public dataset
downloading another dataset for performance evaluation
adding a new cell line or biological system
using GSM7092520/21 as if they were untouched confirmation
redefining an inspected dataset as external validation
using source-author marker findings as substitute outcome evidence
```

Literature-only novelty review is permitted because it does not add training or evaluation data. It must be bounded and cannot delay the frozen Generation-1 execution path.

---

# 5. Frozen WM989 model family

## 5.1 Models

Stage 23.5 carries forward the exact Stage-23 models:

```text
W1 = B + U
W4 = X + B + U
W5 = X + B + U + X×U
```

Definitions:

```text
X   pretreatment clone-level Gene Expression representation
B   exact frozen Stage-23 WM989 nuisance block
U   one-hot coding of the six observed conditions
X×U standardized PCA(X) scores crossed with non-reference treatment dummies
```

The exact nuisance-column list, treatment order, reference treatment, PCA grid, regularization grid, inner selection rule, preprocessing, and outer folds must be loaded from and asserted equal to the frozen Stage-23 protocol/results. They must not be retyped from memory.

## 5.2 Generation-1 model decision

W5 is the Generation-1 predictor. Stage 24 is not authorized to conduct an open neural-architecture tournament.

Stage 24 may refactor, serialize, reproduce, test and package W5. It may not replace W5 because a newly tried architecture looks better on the same outer folds. A more complex architecture belongs to Generation 2 unless a new pre-execution protocol creates genuinely protected evidence.

This is an intentional scope decision: Generation 1 values a supported, interpretable model that ships over an unbounded search for a nominally more sophisticated model.

## 5.3 Frozen feature boundary

```text
eligible X features       36,601 Gene Expression features
forbidden X features     153,055 Custom lineage features
representation           clone-level CP10K/log1p pseudobulk followed by
                         training-only filtering/scaling/PCA
```

No lineage feature, outcome-derived feature, post-treatment expression, source-study resistance label, or outer-test-derived transform may enter `X`.

---

# 6. Generation-1 tool contract

## 6.1 Scientific function

The packaged tool implements:

```text
f(X, B, U) -> future-outcome score Y_hat
```

For one starting clone and all six observed conditions:

```text
{f(X, B, U_1), ..., f(X, B, U_6)}
```

The tool always returns the six frozen per-condition outcome scores. A validated within-system ordering is exposed only when the §8 verdict is `STAGE_25_RANKING_SUPPORTED`. Otherwise the tool reports `RANKING_NOT_SUPPORTED` and does not present the score order as validated condition selection. This reporting distinction does not create another development path: the same predictor/tool ships in either case.

## 6.2 Input contract

Supported input forms:

```text
A. preferred raw form
   pretreatment Gene Expression counts for one or more cells
   stable clone identifier for aggregation
   naive-library/sample identity required to reconstruct B

B. pre-aggregated form
   one gene-aligned clone-level expression vector X
   the complete frozen nuisance vector B
```

Expression alone is not equivalent to the evaluated model when `B` is missing. The tool must fail closed or label the result unsupported; it may not impute a convenient default and claim frozen-benchmark performance.

## 6.3 Candidate-condition contract

Generation 1 supports only:

```text
Dabrafenib
Trametinib
CoCl2
Acid
Cisplatin
Doxorubicin
```

Unknown conditions must return `UNSUPPORTED_TREATMENT`; they must not be embedded, nearest-neighbored, or silently mapped to a known condition.

## 6.4 Output contract

Minimum per-condition output:

```text
condition
future_detection_score
model_version
feature_contract_version
support_status
known_limitations
```

Ranking-dependent claim fields:

```text
validated_condition_order  only if STAGE_25_RANKING_SUPPORTED
ranking_status             always; SUPPORTED or NOT_SUPPORTED
future_abundance_score     secondary; must preserve C2 limitations
calibrated_probability     forbidden unless a later frozen calibration gate passes
```

Required support flags include:

```text
SUPPORTED_KNOWN_CONDITION
UNSUPPORTED_TREATMENT
UNSUPPORTED_FEATURE_SCHEMA
MISSING_REQUIRED_NUISANCE
OUT_OF_CONTRACT_INPUT
RANKING_NOT_VALIDATED
```

## 6.5 Minimum software deliverable

Stage 24 must ship, at minimum:

```text
Python prediction API
command-line interface
frozen model artifact
frozen feature/treatment vocabularies
deterministic preprocessing artifact
machine-readable input/output schemas
model card and claim limitations
one example dataset made only from existing permitted benchmark material
unit tests for schema, feature alignment and deterministic scoring
end-to-end reproduction of frozen W5 outer-fold results
```

A web interface is optional and non-gating. Generation 1 does not wait for product polish.

---

# 7. Stage-24 bounded execution contract

Stage 24 becomes a bounded predictor-engineering stage under the Role-B-primary route:

```text
24A  consume Stage-23.5 handoff; freeze engineering plan
24B  reproduce W1/W4/W5 and exact Stage-23 metrics
24C  implement serialization, preprocessing and prediction API
24D  generate one frozen OOF prediction per clone-condition row
24E  verify deterministic scoring and leakage contracts
24F  freeze W5 tool artifacts
24G  hand the frozen predictions/model to Stage 25
```

Stage 24 PASS requires reproduction of the frozen W5 result to the tolerance defined in 7.1, plus
working deterministic tool artifacts. It does not require a new architecture to beat W5.

## 7.1 Reproduction tolerance — frozen

The Stage-23 pipeline is deterministic given its frozen seeds, so the default expectation is exact
reproduction. "Tolerance-declared" is not a licence to accept whatever comes out.

```text
PRIMARY REQUIREMENT
  the regenerated out-of-fold prediction files are BYTE-IDENTICAL to the committed
  Stage-23 artifacts:
      results/stage23_wm989_detection_oof.csv
      results/stage23_wm989_interaction_oof.csv
  -> reproduction PASSES, no tolerance argument needed

FALLBACK, permitted only when byte-identity fails
  every pooled metric (W1, W4, W5 on C1 and C2) agrees with the committed value to
      absolute difference <= 1e-12
  AND the cause of the non-identity is diagnosed and logged as an environment
      difference (library version, BLAS threading, platform), naming the specific cause
  -> reproduction PASSES as TOLERANCE_DECLARED, with the diagnosis recorded

STOP
  any pooled metric differing by more than 1e-12, or non-identity whose cause cannot be
  named, is an INPUT-INTEGRITY STOP. It permits only a correctness repair to the frozen
  W5 implementation. It does not authorize a new model, metric, endpoint, dataset,
  analysis path, or a relaxed tolerance.
```

The 1e-12 bound is chosen to be far tighter than any real modelling difference and loose enough to
absorb last-bit floating-point variation. A discrepancy larger than that is a defect, not noise.

Stage 24 may not inspect the ranking metric defined below. It generates the frozen inputs Stage 25 consumes.

---

# 8. Preregistered Stage-25 ranking test

## 8.1 Scientific question

Among the six observed WM989 conditions, does the frozen interaction model W5 use pretreatment state to improve **clone-specific condition ordering** over the non-interactive additive model W4?

This is the sole load-bearing new capability test for the Generation-1 ranking claim.

## 8.2 Frozen prediction source

Use only the Stage-24-reproduced, frozen outer-fold OOF predictions from the unchanged W5 and W4 algorithms.

Requirements:

```text
one OOF score per clone-condition row
all six rows for a clone assigned to the same outer fold
no score from a model trained on that clone
no ranking-specific retraining
no ranking-specific hyperparameter selection
```

Stage 24 must reproduce the Stage-23 pooled metrics within the frozen tolerance before handing predictions to Stage 25. A reproduction defect permits only a correctness repair to the frozen W5 implementation; it does not authorize a new model, metric, endpoint, dataset, or analysis path.

## 8.3 Primary endpoint

```text
C1 = post-treatment clone detection / observed-zero endpoint
```

Higher model score means greater predicted probability/propensity of future detection. In the
low-persistence selection analysis, lower is therefore treated as preferable. This is
experimental-condition selection, not clinical treatment recommendation.

The two section-8 statistics use this direction differently, and the difference is deliberate:

```text
8.5  delta_RANK   measures DISCRIMINATION. AUC_i asks whether the model scores the
                  conditions where the clone WAS detected above those where it was not.
                  Higher score for a detected condition is correct.

8.8  delta_TOP1   measures SELECTION UTILITY. It picks the condition with the LOWEST
                  predicted detection score -- the one the model expects to leave the
                  clone undetected -- and asks whether that condition was in fact a zero.
```

Both follow from the same scoring direction. Neither reverses the model.

## 8.4 Primary evaluation population

The primary ranking population contains clones with:

```text
at least one C1-positive condition
AND
at least one C1-zero condition
AND
all six frozen OOF scores present for both W4 and W5
```

The Stage-22 manifest implies:

```text
929 clones detected in >=1 condition
37 clones detected in all 6 conditions
expected informative ranking clones = 892
```

Before scoring, mechanically verify this count from the frozen benchmark. A different count is an input-integrity stop, not permission to redefine eligibility.

All-zero clones and all-positive clones are excluded from the primary within-clone AUROC because that metric is undefined for them. Their counts and characteristics remain reported. They may enter the secondary top-choice analysis only if the protocol below defines their outcome without ambiguity; they cannot rescue the primary gate.

## 8.5 Primary metric — equal-clone-weighted within-clone AUROC

For eligible clone `i`, with predicted scores `s_iu` and binary outcomes `y_iu` over six conditions, compute:

```text
AUC_i = mean over every positive/zero condition pair:
          1.0  if score_positive > score_zero
          0.5  if score_positive = score_zero
          0.0  if score_positive < score_zero
```

Then compute:

```text
R(W) = mean_i AUC_i(W)
```

Every eligible clone receives equal weight regardless of how many positive/zero pairs it contributes.

Primary statistic:

```text
delta_RANK = R(W5) - R(W4)
```

W4 is the primary comparator because its additive `X` term cannot create clone-specific `X×U` ordering. W1 is a required secondary comparator:

```text
delta_RANK_FULL = R(W5) - R(W1)
```

`delta_RANK_FULL` is reported but cannot rescue a failed W5-versus-W4 primary gate.

## 8.6 Bootstrap uncertainty

```text
resampling unit       eligible clone
rows kept together   all six condition rows
replicates            2,000
base seed             23501
interval              two-sided 95% percentile interval
```

For each bootstrap draw, resample eligible clones with replacement and recompute `delta_RANK`. No
model is refit in the bootstrap because uncertainty is over the already-frozen OOF comparison.

**Declared asymmetry.** The permutation null of 8.7 refits the entire pipeline; this bootstrap does
not. The interval is therefore *conditional on the fitted models* and quantifies sampling
variability over clones only -- it does not include model-fitting variability, whereas the null
does. This is standard for a frozen-prediction comparison and it is deliberate, but the two
uncertainty statements are not the same kind of quantity and must not be described as if they were.
Report the interval as a clone-resampling interval, not as a full-pipeline interval.

Bootstrap criterion:

```text
lower endpoint of CI95(delta_RANK) > 0
```

## 8.7 Full-refit permutation null

The ranking null must mirror the Stage-23 WM989 expression-permutation geometry rather than inventing a convenient label shuffle.

Frozen WM989 strata:

```text
depth bin:
  1
  2
  3-4
  5-9
  10+

crossed with:
  3-bit naive-sample presence pattern
  (Naive1 present?, Naive2 present?, Naive3 present?)
```

Use the exact deterministic small-stratum merge rule frozen in Stage 23 V2.

For each outer fold and permutation:

```text
permute intact clone-level CP10K/log1p Gene Expression profiles
among outer-training clones within stratum

permute intact clone-level CP10K/log1p Gene Expression profiles
among outer-test clones within stratum

never move a profile across the train/test boundary
never permute genes independently
keep each clone's six-condition outcome vector intact
keep B, U, outcomes and outer-fold identities fixed
```

Rerun the complete relevant pipeline:

```text
training-derived gene filtering
scaling
inner-split PCA
inner hyperparameter selection
W4 fitting
W5 fitting
outer-test prediction
ranking-eligibility assertion
delta_RANK calculation
```

Do not reuse observed-data hyperparameters in null runs. Cache only quantities mathematically invariant to profile reassignment, under the Stage-23 caching rules.

Permutation count and seed:

```text
n_perm       1,000
base seed    23523
no early stopping
```

Draw `b` must be seeded as a function of `b` alone, so that shards are order-independent and
resumable without changing any value. Each shard writes its OWN cache file (0.2). The merge asserts
that all 1,000 indices are present before any statistic is computed; a missing index is an
integrity stop, not a smaller null. Budget approximately 19-20 h across 3 shards (0.2).

Permutation p-value:

```text
p_perm = (1 + number of null delta_RANK >= observed delta_RANK) / 1,001
```

Permutation criterion:

```text
observed delta_RANK > percentile_95(null)
AND
p_perm <= 0.05
```

Report the number of exceedances, empirical p-value, null mean, null SD, null p95, observed-null-p95 margin, actual fixed-clone fraction, and every deterministic stratum merge.

## 8.8 Secondary decision-utility diagnostic

For each primary-eligible clone, select the condition with the **lowest** predicted C1 detection score:

```text
u_star(W) = argmin_u score_iu(W)
```

Tie breaking uses the canonical treatment order in §6.3 and never outcome information.

Define:

```text
LOW_PERSISTENCE_TOP1(W)
  = mean_i 1[y_i,u_star(W) = 0]

delta_TOP1
  = LOW_PERSISTENCE_TOP1(W5) - LOW_PERSISTENCE_TOP1(W4)
```

`delta_TOP1` is reported with a bootstrap interval computed from the **same** clone resampling as
8.6 -- the same 2,000 draws, the same seed, the same eligible clones -- at no additional cost, since
both statistics are functions of the same frozen OOF predictions.

**What kind of criterion this is.** V1 of this plan required `delta_TOP1 > 0` for
`STAGE_25_RANKING_SUPPORTED` while simultaneously stating it was "not independently thresholded for
significance." Those two statements are inconsistent: a bare point estimate with no uncertainty
control was being used to veto a result that had cleared both a bootstrap interval and a 1,000-draw
permutation null. With 892 clones and a coarse `argmin` statistic, `delta_TOP1` can fall below zero
by chance while the ranking improvement is real, so as written the criterion created an
unjustified false-negative path.

It is therefore restated as what it actually is:

```text
CRITERION 5, restated
  delta_TOP1 >= 0        a DIRECTIONAL-CONSISTENCY CHECK, not a significance test

  purpose   a strictly negative top-1 utility alongside a significant positive
            delta_RANK would mean the ranking improvement does not survive contact
            with the decision the ranking is for. That is worth blocking on.
  not       evidence of utility in its own right. It is not powered, carries no
            threshold, and its bootstrap interval is reported but not gated.
```

This remains deliberately conservative and may still produce a false negative. That is accepted
and stated. `delta_TOP1` cannot overturn an unsupported primary `delta_RANK` result in the other
direction: it can only withhold support, never grant it.

This diagnostic does not establish clinical treatment utility because C1 is an observed detection proxy and the six conditions include non-clinical stress contexts.

## 8.9 Secondary/descriptive reporting

Report without allowing rescue:

```text
R(W1), R(W4), R(W5)
delta_RANK_FULL = R(W5)-R(W1)
pairwise condition-ranking accuracy matrix
ranking results by outer fold
ranking results by pretreatment-depth bin
score-tie rate
all-zero and all-positive clone counts
```

C2 remains secondary endpoint corroboration from Stage 23. No C2 ranking analysis is authorized for Generation 1, regardless of the C1 result.

## 8.10 Ranking result and mandatory shipment

`STAGE_25_RANKING_SUPPORTED` requires all of:

```text
1. observed delta_RANK > 0
2. lower CI95 bootstrap endpoint > 0
3. observed delta_RANK > null p95
4. p_perm <= 0.05 under the 1,000-draw full-refit null
5. delta_TOP1 >= 0   (directional-consistency check, 8.8; not a significance test)
6. all integrity, leakage and determinism checks pass
```

If one or more criteria are not met, record:

```text
STAGE_25_RANKING_NOT_SUPPORTED
```

This is a terminal scientific result, not a route to another analysis. It requires all of the following:

```text
the negative/insufficient ranking result is reported in the manuscript
the six per-condition W5 outcome scores remain available in the tool
validated_condition_order is omitted
the tool reports ranking_status = NOT_SUPPORTED
Generation 1 proceeds immediately to evidence lock and shipment
```

There is no third scientific verdict and no Stage-25 result that returns the project to Stage 22, 23, 23.2, 23.5 or 24. Integrity or software defects must be corrected to execute this exact frozen test; they cannot change its scientific contract.

## 8.11 Anti-rescue firewall

After this plan is frozen, the following are forbidden:

```text
changing the primary metric
changing clone eligibility
using micro- instead of equal-clone weighting because it looks better
changing W4 as primary comparator
changing the null construction
reducing permutation count
early stopping permutations
choosing a different endpoint
introducing a new model
changing treatment direction or tie handling
subsetting treatments because Doxorubicin hurts
dropping Cisplatin because its interaction is negligible
using C2 or a per-treatment result to rescue failed C1 ranking
adding another dataset
```

All scientific ambiguities must be resolved before this plan is frozen. After freeze, no post-result plan revision is authorized for Generation 1. A purely mechanical correctness defect may be repaired and logged only when the repair leaves the dataset, population, model, metric, null, threshold and reporting rules unchanged. It may not create an alternative analysis.

---

# 9. Generation-1 paper and tool completion path

The revised finite path is:

```text
STAGE 23.5
  claim revision + ranking preregistration + digest
        -> STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1

STAGE 24
  reproduce, freeze and package W5 predictor/tool
        -> STAGE_24_GEN1_TOOL_READY

STAGE 25
  execute the single frozen ranking analysis
  record exactly one result field:
        STAGE_25_RANKING_SUPPORTED
        or STAGE_25_RANKING_NOT_SUPPORTED
  either value -> GEN1_MANDATORY_SHIP

STAGE 26
  record KNOWN_TREATMENT_ONLY_SCOPED_LIMIT
  no unseen-treatment claim and no rescue experiment

GEN-1 EVIDENCE LOCK
  freeze benchmark, tool, OOF predictions, ranking verdict and limitations

GEN-1 CLAIM LOCK
  freeze abstract-level allowed/forbidden claims

MANUSCRIPT + REPRODUCIBILITY PACKAGE
        -> PREPRINT / SUBMISSION

GENERATION 2, NON-BLOCKING FOR GEN-1
  new-system biological replication
  unseen-treatment transfer
  broad calibration / OOD validation
  more complex architectures
```

Independent biological replication remains scientifically important. Moving it to Generation 2 narrows the Generation-1 claim; it does not manufacture a PASS.

---

# 10. Manuscript contribution contract

Generation 1 must be written around positive, bounded contributions rather than around a universal model that the evidence cannot support.

## Contribution 1 — prospective benchmark

An auditable lineage-grounded benchmark separating:

```text
pretreatment molecular state
candidate condition
independently observed later clonal outcome
```

with clone-held-out splits, outcome-feature firewalls, captured-abundance baselines, and full-refit permutations.

## Contribution 2 — supported treatment-conditioned predictor

In WM989, explicit `X×U` improves known-condition prediction beyond additive `X+U+B` and nuisance+treatment `B+U`, with strong separation from the frozen null and corroboration on C2.

Required limitation in the same presentation:

```text
captured abundance is stronger than state
four conditions carry meaningful interaction
Cisplatin is negligible on C1
Doxorubicin is negative on both endpoints
no independent biological Role-B replication
```

## Contribution 3 — preregistered ranking result

The ranking experiment is always included in the paper. `STAGE_25_RANKING_SUPPORTED` permits a positive ranking claim. `STAGE_25_RANKING_NOT_SUPPORTED` requires a direct statement that clone-specific condition ordering was not supported. The paper and tool ship in either case; the result cannot trigger another experiment.

## Supporting boundary result — Rewind

Rewind demonstrates both prospective promise and the cost of weak outcome geometry: positive direction across two independent replicates, but failed frozen power. It is supporting evidence and a limitation analysis, not the primary success.

## Novelty wording

Until a bounded final prior-art audit is complete, use:

```text
"we introduce"
"we evaluate"
"we show under this benchmark"
```

Do not use:

```text
"the first"
"no previous work"
"unprecedented"
```

without direct, documented support.

---

# 11. Required Stage-23.5 artifacts

Before Stage 24 opens, create and commit:

```text
plans/(newer)practical plans/STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md
results/stage23_5_protocol.json
plans/(newer)practical plans/RECORDs/stage_23_5_RECORD.md
results/stage23_5_handoff_to_stage24.json
```

`stage23_5_protocol.json` must record at minimum:

```text
plan filename
plan version
canonical-LF SHA-256
repository commit containing the frozen plan
source artifact paths and SHA-256 values
two roadmap amendments
dataset firewall
W1/W4/W5 definitions
ranking population rule
primary metric
bootstrap contract
permutation contract
verdict logic
forbidden rescue actions
reproduction tolerance and stop rule
permutation compute budget as accepted
```

The plan cannot contain its own final commit hash without recursion. The companion protocol and Stage-23.5 record are the authoritative digest-bearing freeze artifacts.

The Stage-23.5 handoff must state mechanically:

```text
primary role                       Role B / WM989
Role-A status                      underpowered supporting evidence
Stage-24 model                     frozen W5
Stage-24 opening status            STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1
ranking metric inspected?          false
ranking protocol hash frozen?      true
new datasets authorized?           false
Stage-27 Gen-1 gate?               false
allowed generalization             known conditions, within system
```

---

# 12. Independent audit checklist

An independent audit must verify before commit:

```text
[ ] every quoted Stage-23/23.2 number matches its source artifact
[ ] p_perm is described as 0/200 and 1/201, not as precise tail estimation
[ ] observed/null ratios are not called effect sizes
[ ] four / negligible Cisplatin / negative Doxorubicin wording is preserved
[ ] no independent Role-B biological replication is claimed
[ ] all Rewind evidence is marked consumed for future confirmation purposes
[ ] the two roadmap amendments are separately named
[ ] no new dataset is authorized
[ ] W5/W4/W1 definitions match the frozen Stage-23 implementation
[ ] B and treatment coding are referenced mechanically, not re-invented
[ ] primary ranking population is mechanically reproducible
[ ] whole-clone equal weighting is explicit
[ ] score ties receive 0.5 in AUC_i
[ ] all six treatment rows stay together
[ ] full-refit permutation preserves the outer boundary
[ ] depth × naive-presence strata and merge rules match Stage 23
[ ] 1,000 permutations and no early stopping are enforced
[ ] C2 cannot rescue C1 ranking
[ ] W5 tool and manuscript shipment do not depend on ranking support
[ ] unsupported ranking removes only the validated ordering claim
[ ] both ranking labels lead directly to GEN1_MANDATORY_SHIP
[ ] no result can return the project to an earlier stage
[ ] no calibrated probability or clinical recommendation is promised
[ ] the tool requires B or derives it from supported raw input
[ ] Stage 24 cannot inspect the ranking result
[ ] no Stage-25 statistic exists before plan digest and commit
[ ] the reproduction tolerance is a NUMBER, with a stated stop rule (7.1)
[ ] the permutation compute budget is stated and has been accepted (0.2)
[ ] each permutation shard writes its own cache file and the merge asserts completeness
[ ] delta_TOP1 is labelled a directional-consistency check, not a significance test
[ ] the bootstrap interval is declared conditional on the fitted models (8.6)
[ ] the 892-clone eligibility count is verified mechanically before scoring (8.4)
[ ] section 0 states the six criteria and the forbidden list
```

Every item must be settled during the single pre-freeze audit. Once the plan digest and commit are recorded, no new Generation-1 analysis route or post-result amendment is permitted.

---

# 13. Mechanical Stage-23.5 verdict

`STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1` requires:

```text
1. this plan is independently audited and committed;
2. canonical-LF SHA-256 is recorded in stage23_5_protocol.json;
3. both roadmap amendments are recorded explicitly;
4. the no-new-data firewall is active;
5. the ranking metric, null, threshold and mandatory reporting rule are frozen and uninspected;
6. the Stage-23.5 handoff is complete and internally consistent;
7. no historical Stage-23 or Stage-23.2 verdict is changed.
```

Stage 23.5 does not have a performance PASS. Its successful exit means only that the revised Generation-1 claim architecture and its load-bearing future test are frozen honestly enough for Stage 24 to begin.

---

# 14. One-line Generation-1 contract

> Build and ship an interpretable known-condition `X + U -> Y` predictor from the existing WM989 lineage benchmark; test clone-specific ranking once under a frozen full-refit null; record supported or not supported; then ship without redesign, another dataset, another metric, another model, or a return to an earlier stage. Report Rewind as underpowered support and make no external, clinical, unseen-treatment, calibrated-probability, or independent-replication claim.
