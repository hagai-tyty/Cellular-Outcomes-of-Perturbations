# CellFate-Rx — Current Vision Roadmap v5

**V5 scope:** This revision does not change any executed Stage-21/22/23 result or gate. It only formalizes the post-Stage-23 Role-A resolution path before Stage 24, corrects the unresolved Rewind replicate-status assumption, adds outcome-label reliability as a distinct diagnostic, and requires independent confirmation before a blocked Role-A path can reopen Stage 24.


**Primary scientific target:** Build a genuinely new prospective treatment-conditioned fate model, not merely validate the existing classifier on another dataset.

The intended paper is built around:

```text
pre-intervention molecular state X
+
candidate perturbation U
        ↓
independently measured future outcome Y
```

and, where the experimental design permits:

```text
for the same / comparable starting state:
rank which candidate perturbation is most likely to succeed
```

The old CellFate-Rx audit remains valuable supporting evidence and may become a separate methods paper, but the main research arc now creates a new benchmark, a new prospective model, and a new state-conditioned treatment-ranking capability.

---

## STAGE 21 — Prospective Data Qualification

Audit the local datasets once, then qualify and reconstruct public lineage-resolved prospective data. The stage seeks **two complementary capabilities**: a reprogramming dataset with pre-outcome RNA and an independently measured later outcome (**Role A — required**), plus a multi-perturbation lineage dataset where treatment identity genuinely varies (**Role B — high value, non-blocking**). Stage 21 ends when the Role-A data are reconstructable, statistically usable, leakage-checked, and frozen for benchmark construction; a qualified Role B is folded in on the same terms when the search finds one.

**PASS → STAGE 22**, on either path:

```text
FULL_DATA_PATH   Role A + Role B qualified
CORE_DATA_PATH   Role A qualified, Role B not found inside the frozen search budget
                 -> proceed; the treatment-ranking contribution is SCOPED DOWN
DATA_BLOCKED     Role A not found -> do not build the prospective architecture yet
```

*Corrected 2026-08-21: this paragraph originally read "The stage deliberately seeks two
complementary capabilities … Stage 21 ends only when those data are reconstructable, statistically
usable, leakage-checked, and frozen for benchmark construction", i.e. both roles required. That was
V2's PASS rule. It is superseded by `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §9.8, which
makes Role B non-blocking so a missing multi-perturbation dataset scopes the treatment-ranking
contribution down instead of holding the paper. Executed 21C returned `FULL_DATA_PATH`, so both
roles are in fact qualified — the rule is corrected here regardless, because it governs what a PASS
means, not what happened to be found.*

---

## STAGE 22 — Prospective Benchmark Build

Create a standardized lineage-grounded benchmark from the qualified datasets. Map each experiment into a common `X_before + U -> Y_future` schema, define biologically correct split groups, freeze the outcome ontology, metrics, treatment metadata, leakage rules, and external-test partitions, and reproduce the source-study bookkeeping. The benchmark itself is one contribution: it separates true prospective prediction from post-treatment state classification and prevents random-cell or same-state shortcuts.

**PASS → STAGE 23**

---

## STAGE 23 — Learnability & Interaction Gate

Before building a neural model, prove that the benchmark contains the information needed for the CellFate-Rx thesis. Stage 23 is executed as a frozen sequence of substages: protocol/representation freeze, Role-A Rewind learnability, Role-B additive state signal, explicit `X×U` interaction, permutation/leakage controls, and mechanical evidence synthesis. Compare prevalence/metadata/treatment-only, molecular-state-only, additive `X+U`, and explicit `X×U` baselines under the Stage-22 clone-level holdouts.

Role A remains the mandatory prospective anchor at this stage. Role B is high-value and may establish additive or treatment-specific state signal, but it does not silently rescue a failed Role-A gate. A provisional bootstrap/model-performance PASS is not final until the pre-registered permutation and structural controls survive.

**PASS → STAGE 24**

**`STAGE_24_BLOCKED_ROLE_A` → STAGE 23.2**

---

## STAGE 23.2 — Role-A Resolution / Failure Decomposition

Stage 23.2 is an **interstitial failure-resolution stage between Stage 23 and Stage 24**. It opens only if Stage 23 closes with `STAGE_24_BLOCKED_ROLE_A`. It is not a model-rescue stage, and it cannot rewrite the historical Stage-23 verdict: the failed Stage-23 Role-A permutation result remains failed permanently.

The goal is to determine **why** Role A failed, what can and cannot be learned from the already-inspected Rewind data, and what new evidence would be required before the roadmap may reopen Stage 24.

### Mandatory plan freeze before analysis

Before any decomposition analysis runs, create, independently audit, and commit:

```text
STAGE_23_2_ROLE_A_RESOLUTION_V1.md
```

That plan must freeze:

```text
the diagnostic hypotheses
the exact data each diagnostic may inspect
the experimental unit / grouping rules
null constructions
metrics and decision rules
which analyses are exploratory
what evidence is reserved for confirmation
the final roadmap-resolution logic
```

No 23.2B–23.2G analysis may begin before that plan is frozen.

### Five non-mutually-exclusive explanations

Stage 23.2 must quantify five explanations. They are **not mutually exclusive**; more than one may contribute to the observed failure.

```text
1. NO ROBUST BIOLOGICAL SIGNAL
   Pretreatment molecular state may not contain a reproducible Role-A signal
   for the frozen outcome.

2. MODEL-SELECTION NULL INFLATION
   The larger / different model-selection path for the state model may produce
   positive ΔAP under the full refit permutation null.

3. RESIDUAL DEPTH / SAMPLING STRUCTURE
   The frozen permutation may preserve sampling/depth structure that remains
   partially encoded in permuted expression and is not fully absorbed by the
   nuisance baseline.

4. OUTCOME-LABEL / MEASUREMENT LIMITATION
   The thresholded gDNA-derived Role-A label may be noisy, unstable near the
   cutoff, or an imperfect proxy for the biological event of interest.

5. POWER / EXPERIMENTAL-UNIT LIMITATION
   With 35 positive clones among 3,147 retained clones, the experiment may be
   too weak to distinguish the observed effect from the appropriate null;
   interpretation also depends on the still-unresolved biological-replicate
   structure of the two control GSMs.
```

### Ordered substages

```text
23.2A  resolution protocol + source-design freeze
23.2B  model-selection null decomposition
23.2C  depth / nuisance null decomposition
23.2D  outcome-label reliability
23.2E  power / identifiability analysis
23.2F  diagnostic synthesis + confirmatory protocol freeze
23.2G  independent confirmation / roadmap resolution
```

### 23.2A — Resolution protocol + source-design freeze

In addition to freezing the statistical protocol, 23.2A must resolve the Rewind experimental-unit semantics as far as the source materials permit.

Stages 21–23 established two control GSMs:

```text
GSM7092515   2,030 retained benchmark cells
GSM7092516   1,875 retained benchmark cells
```

They were treated operationally as lanes, and **306 retained clones span both GSMs**. The existing record does not establish whether these GSMs are:

```text
two biological replicates
two technical / 10X lanes from one biological culture
or unresolved from the available source metadata
```

Do not infer biological replication from the number of GSM accessions. Audit the source metadata and author materials. If the design still cannot be established, freeze it explicitly as `REPLICATE_STRUCTURE_UNRESOLVED`; do not silently choose the interpretation that gives the most favorable power calculation.

Until that audit resolves the design, the supported statement is:

```text
35 positive clones among 3,147 retained clones;
biological-replicate structure unresolved at Stage 23.
```

### 23.2B — Model-selection null decomposition

Test whether the positive Role-A permutation-null center can be explained by differences in model-selection flexibility rather than biological signal.

This substage must isolate model-selection effects **without changing the historical Stage-23 test**. Any equalized-search, fixed-complexity, or otherwise simplified comparison is a new diagnostic analysis. It may explain the failure, but it cannot retroactively convert Stage 23 into a PASS.

### 23.2C — Depth / nuisance null decomposition

Test whether the Stage-23 permutation preserves residual sampling/depth information that the nuisance model does not fully absorb.

The diagnostic must distinguish:

```text
expression-derived state information
from
captured-depth / lane / sampling information that survives the null
```

Do not weaken the original permutation simply to make the observed statistic look more exceptional. The Stage-23 null remains authoritative for the historical claim.

### 23.2D — Outcome-label reliability

The Role-A outcome is an **operational gDNA-derived threshold label**, not a noise-free oracle. Stage 23 contains 35 positive clones at 1.11% prevalence, but the negative class is an operational complement rather than independently proven biological failure.

23.2D must audit the measurement rule itself, including:

```text
source-rule reconstruction
rank / count stability around the cutoff
cutoff ties and near-cutoff uncertainty
sensitivity to plausible measurement noise
whether the operational negative class is biologically interpretable
```

Alternative thresholds, soft labels, continuous gDNA quantities, or alternate outcome definitions inspected during Stage 23.2 are **exploratory diagnostics** unless frozen before evaluation on independent evidence. They may show that the original endpoint is measurement-limited; they may not be selected because they make prediction look better.

### 23.2E — Power / identifiability analysis

Power analysis comes **after** 23.2A–23.2D because the appropriate experimental unit, null behavior, and label reliability determine what "power" means.

It must condition on:

```text
35 positive / 3,147 total retained clones
the frozen clone-level outer-fold geometry
the resolved or explicitly unresolved replicate structure
the observed Stage-23 permutation-null distribution
the estimated label reliability / measurement limitation
```

The goal is not to calculate a favorable retrospective power number. The goal is to determine which effect sizes and conclusions are actually identifiable from the existing experiment and what additional independent sample size / biological replication would be needed.

### 23.2F — Diagnostic synthesis + confirmatory protocol freeze

23.2F combines the diagnostics without forcing them into a single-cause story.

Report each mechanism independently as:

```text
SUPPORTED
NOT_SUPPORTED
UNRESOLVED
```

at minimum for:

```text
MODEL_SELECTION_NULL_INFLATION
RESIDUAL_DEPTH_STRUCTURE
OUTCOME_LABEL_LIMITATION
POWER_LIMITATION
ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA
```

Then freeze the next-step requirement **before any new confirmatory evidence is inspected**.

If a corrected model, nuisance representation, outcome definition, or null is proposed after seeing Stage-23/23.2 results, that proposal is exploratory on Rewind. A separate confirmatory protocol must specify what untouched evidence can confirm or reject it.

#### Mandatory Stage-24 handoff contract

23.2F must also create a machine/readable + human-readable Stage-24 readiness dossier before confirmation begins:

```text
STAGE_23_2_HANDOFF_TO_STAGE_24.md
stage23_2_handoff_to_stage24.json
```

The handoff must freeze, at minimum:

```text
ROLE A
  benchmark / outcome version
  exact future-outcome definition
  experimental unit and grouping rule
  replicate-status conclusion or explicit unresolved status
  nuisance variables allowed
  primary metric / decision statistic
  confirmed effect and null-test result
  which evidence was used for diagnosis
  which evidence is reserved for confirmation
  claims allowed / forbidden

ROLE B
  frozen Stage-23 additive + interaction verdicts
  strongest simple baseline to beat (including W5 where applicable)
  exact nuisance block B
  treatment coding
  C1 / C2 endpoint roles
  treatment-level exceptions / limitations that must remain visible

GLOBAL
  feature universe / preprocessing contract
  benchmark version(s)
  split / grouping policy
  datasets already inspected
  datasets reserved for later independent replication
  unresolved limitations
  exact Stage-24 opening condition
```

If Stage 23.2 changes a **material benchmark semantic** — for example the outcome definition, positive/negative ontology, experimental unit, leakage firewall, grouping/split rule, or source reconstruction — Stage 24 may **not** open directly from a successful-looking corrected analysis. The revised benchmark must be versioned and the affected Stage-22/Stage-23 qualification, learnability, structural and permutation gates must be rerun under the revised frozen contract. This prevents Stage 23.2 from bypassing the benchmark and learnability gates by redefining the problem after observing failure.

### 23.2G — Independent confirmation / roadmap resolution

**Diagnosis is not confirmation.**

A blocked Role-A path may reopen Stage 24 only if **both** conditions hold:

```text
1. a pre-registered confirmatory analysis succeeds on evidence that was not
   used to design the Stage-23.2 correction; and
2. STAGE_23_2_HANDOFF_TO_STAGE_24 is complete and benchmark-compatible,
   with no material benchmark change left un-re-gated.
```

Suitable confirmatory evidence may include:

```text
a genuinely independent biological replicate
an independent Rewind-like dataset
a prospectively held-out clone cohort never inspected during Stage 23.2 design
or another independent future-outcome anchor
```

If no suitable untouched evidence is available, Stage 23.2 may finish diagnostically but Stage 24 remains blocked.

Evidence used in Stage 23.2 to reopen Stage 24 is **not** an untouched Stage-27 replication set if it influenced the architecture, protocol, or roadmap decision. Stage 27 must preserve an independent replication test of the eventual frozen model.

### Final Stage-23.2 roadmap statuses

Diagnostic findings may be multi-label, but the roadmap exit must be one of:

```text
ROLE_A_CONFIRMATORY_SUPPORTED
    New, pre-registered, independent confirmatory evidence supports the
    corrected Role-A claim AND the Stage-24 handoff contract is complete
    and benchmark-compatible.
    -> Stage 24 may reopen.

ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
    The failure is explainable or potentially underpowered, but no untouched
    evidence has confirmed a corrected Role-A claim.
    -> Stage 24 remains blocked.

ROLE_A_REDESIGN_REQUIRED
    The current outcome / nuisance / experimental design is not adequate for
    the intended Role-A claim.
    -> freeze a revised benchmark/outcome protocol and obtain new evidence;
       Stage 24 remains blocked.

ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL
    After decomposition, no defensible Role-A signal remains worth preserving
    as the mandatory anchor.
    -> explicitly revise the roadmap and claim architecture before deciding
       whether a different evidence path may open Stage 24.
```

A positive Role-B result from Stage 23 remains valid evidence throughout Stage 23.2, but it cannot silently substitute for Role A. Making Role B the primary anchor would require an explicit roadmap/claim revision.

**`ROLE_A_CONFIRMATORY_SUPPORTED` → STAGE 24**

**`ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE` / `ROLE_A_REDESIGN_REQUIRED` → REMAIN BLOCKED**

**`ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL` → EXPLICIT ROADMAP REVISION BEFORE ANY STAGE-24 DECISION**

---

## STAGE 24 — Build CellFate-Rx Prospective

Stage 24 is the main invention stage and opens only through a valid Stage-23 PASS or a valid Stage-23.2 `ROLE_A_CONFIRMATORY_SUPPORTED` exit with a complete benchmark-compatible handoff contract.

Stage 24 is a **development / architecture-selection stage**, not the independent biological replication stage. The benchmark may be reused under locked grouped/nested evaluation, but because Stage-23 results have already been inspected, Stage-24 performance on those same systems must not be described as untouched external confirmation. Stage 27 preserves that role.

### Mandatory Stage-24 plan freeze

Before fitting any new Stage-24 model, create, independently review, and commit:

```text
STAGE_24_PROSPECTIVE_MODEL_V1.md
```

The Stage-24 plan must consume `STAGE_23_2_HANDOFF_TO_STAGE_24.*` when Stage 23.2 was used, plus the frozen Stage-23 Role-B artifacts. It must freeze:

```text
task / endpoint contract
dataset roles
architecture family and complexity budget
whether tasks are trained separately or with shared representations
how heterogeneous Role-A / Role-B outcomes are handled
input feature universe
perturbation encoding
nuisance handling
interaction module
training loss(es)
outer / inner grouping
hyperparameter search budget
primary complexity gate
secondary metrics
ablation set
seeds / determinism requirements
failure / stop rules
Stage-25 handoff outputs
```

Do not pool biologically different outcome labels into one head merely because both fit the generic `X_before + U -> Y_future` schema. If Role A and Role B use different outcome semantics, Stage 24 must either use explicitly task-specific heads / losses or fit role-specific models under a shared architecture family. That choice must be frozen before model fitting.

### Architecture constraint inherited from Stage 23

Role B already determines whether treatment interaction is justified:

```text
interaction gate supported
    -> Stage 24 may include explicit state × perturbation interaction

interaction gate unsupported
    -> keep the treatment-conditioned component additive or scope it down
```

When the interaction gate is supported, the architecture should expose separate state and perturbation representations before their interaction rather than hiding treatment identity inside an undifferentiated feature vector.

Role A, where treatment does not provide the same multi-treatment variation, must not be used to fabricate evidence for treatment interaction. Its purpose is the prospective pre-state → future-outcome anchor defined by the Stage-23/23.2 contract.

### Ordered substages

```text
24A  task / claim / architecture protocol freeze
24B  benchmark adapters + exact baseline reproduction
24C  architecture implementation + engineering contracts
24D  nested training / model selection on training groups only
24E  locked grouped OOF evaluation on known-treatment benchmark
24F  complexity + ablation gate
24G  model freeze + Stage-25 handoff
```

### 24A — Task / claim / architecture protocol freeze

Freeze the exact Role-A and Role-B tasks that Stage 24 is allowed to model. Explicitly identify:

```text
which task is the mandatory prospective anchor
which Role-B endpoint is primary
which endpoints are secondary
which metrics choose hyperparameters
which metrics are reporting-only
which model outputs will later be used by Stage 25
```

No architecture search may use Stage-25 treatment-ranking results, Stage-26 held-out-treatment results, Stage-27 replication results, or Stage-28 calibration/OOD results.

### 24B — Benchmark adapters + exact baseline reproduction

Before training the new architecture, reproduce the strongest frozen simple baseline(s) from Stage 23 on the Stage-24 evaluation code path. At minimum, where applicable:

```text
Role B:
  W5 = X + B + U + X×U
  C1 detection baseline
  C2 conditional-abundance baseline / secondary endpoint

Role A:
  the exact simple baseline and endpoint frozen by Stage 23.2 confirmation
```

If those baselines cannot be reproduced within the frozen tolerance, stop before fitting the new model.

### 24C — Architecture implementation + engineering contracts

Implement only the architecture family frozen in 24A. Required engineering contracts include:

```text
no target / post-treatment leakage
group identity never enters X
training-only preprocessing
deterministic seed handling
one canonical feature-ID order
treatment coding frozen
nuisance block frozen
explicit task-head semantics
checkpoint / config provenance
no silent fallback to a different architecture
```

### 24D — Nested training / model selection

All preprocessing, early stopping, hyperparameter selection, regularization selection, and architecture choices that remain variable must be made using outer-training data only. Outer-test groups may be scored only after the corresponding outer model is frozen.

The Stage-24 plan must pre-register a bounded search budget. Do not expand the search after seeing outer results.

### 24E — Locked known-treatment evaluation

Generate one frozen OOF prediction per valid biological evaluation unit under the grouped Stage-24 protocol. Compare against the strongest Stage-23 simple baseline using the pre-registered primary metric.

Stage 24 may report known-treatment prospective performance, but it may not claim:

```text
unseen-treatment generalization     -> Stage 26
independent biological replication  -> Stage 27
final calibration / OOD utility     -> Stage 28
```

### 24F — Complexity + ablation gate

The new architecture must **earn its complexity**. Before fitting begins, `STAGE_24_PROSPECTIVE_MODEL_V1.md` must define one mechanical primary complexity rule relative to the strongest simple baseline.

The gate should require a pre-registered improvement on the primary prospective endpoint, with uncertainty / grouped inference appropriate to the benchmark. Secondary improvements may support interpretation but cannot be chosen post hoc to rescue a failed primary gate.

Required ablations should isolate, where applicable:

```text
state encoder contribution
perturbation encoder contribution
state×perturbation interaction contribution
nuisance contribution
task-head / hurdle structure contribution
```

If the neural / complex architecture fails its frozen complexity gate, do not proceed by changing the metric or widening the search. Retain the strongest simple baseline, record `STAGE_24_COMPLEXITY_NOT_EARNED`, and explicitly decide whether the roadmap should stop, redesign, or carry the simple model forward.

### 24G — Model freeze + Stage-25 handoff

On PASS, freeze:

```text
model code / config
weights / seeds
feature and treatment vocabularies
task heads
preprocessing
baseline comparison
OOF predictions
primary / secondary metrics
ablation results
known limitations
exact score used for Stage-25 treatment ranking
```

Stage 25 must consume the frozen Stage-24 model and score without model retraining or result-driven redefinition.

**`STAGE_24_MODEL_PASS` → STAGE 25**

**`STAGE_24_COMPLEXITY_NOT_EARNED` → STOP / EXPLICIT REDESIGN OR SIMPLE-BASELINE ROADMAP DECISION**

---

## STAGE 25 — State-Conditioned Treatment Ranking Challenge

Using the frozen Stage-24 model and score, test the capability that makes the new system useful: for the same or comparable starting molecular state, can it distinguish which candidate treatment is more likely to produce the desired future outcome? Use multi-treatment lineage data to evaluate pairwise treatment ranking, within-state treatment selection, and decision utility against treatment-only and state-only baselines. Do not retrain the model to improve Stage-25 ranking. This should become one of the paper's central figures because it demonstrates that the model has learned `cell state × perturbation`, not merely "good cells" or "good treatments."

**PASS → STAGE 26**

---

## STAGE 26 — Held-Out Perturbation / Generalization Challenge

Where the data make it scientifically resolvable, hold out entire perturbations, treatment classes, experimental batches, or biological contexts and test whether the frozen model transfers. This is a harder claim than prediction under already-seen treatments, so it is separated from Stage 25. A failure here does not erase the known-treatment prospective result, but it determines whether the paper may claim generalization to unseen interventions.

**PASS OR SCOPED LIMIT → STAGE 27**

---

## STAGE 27 — Independent Biological Replication

Freeze the model and core evaluation rules, then test the main prospective result on a separate biological dataset/system not used to choose the architecture, the Stage-23.2 correction, or the decision to reopen Stage 24. Prefer a second lineage-resolved system with an independent future outcome. The goal is to show that the prospective signal and treatment-conditioning result are not peculiar to one experiment. Evidence consumed during Stage 23.2 confirmation is not automatically eligible as this untouched replication set. No model fishing is allowed against the replication set.

**PASS → STAGE 28**

---

## STAGE 28 — Calibration, OOD & Decision Utility

Turn the predictor into a trustworthy decision instrument. Evaluate probability calibration against the actual hard future outcome, uncertainty/risk-coverage, abstention behavior, OOD detection, and whether rejecting uncertain cases improves decision quality. Reuse the lessons from the old Stage-16/17 calibration work, but calibrate only against the prospective ground truth. This stage is where the engineering history becomes a methodological advantage rather than baggage.

**PASS → STAGE 29**

---

## STAGE 29 — Final Evidence Lock

Freeze benchmark versions, model weights, splits, baselines, calibration, OOD rules, metrics, negative controls and manuscript-facing figures. Run the complete final evaluation once under the locked protocol. After this point, result-driven tuning is forbidden; anything new is a separately labelled experiment.

**PASS → STAGE 30**

---

## STAGE 30 — Scientific Claim Lock

Write the authoritative claim document before drafting the manuscript. Separate what is established under known treatments from what generalizes to unseen treatments, distinguish biological outcomes from transcriptomic surrogates, record every limitation, and state exactly which old CellFate-Rx claims were withdrawn. This becomes the claim firewall for every abstract, figure caption and conclusion.

**PASS → STAGE 31**

---

## STAGE 31 — Main Manuscript Build

Build the main paper around the contributions that survive the frozen gates: **(1)** a prospective lineage-grounded benchmark, **(2)** a prospective treatment-conditioned model if Stage 24 is legitimately opened, **(3)** state-conditioned treatment ranking only if the interaction/ranking gates support it, and **(4)** calibrated uncertainty / abstention under biologically correct holdouts where the outcome supports calibration. The old project's circularity, time-confounding, scorecard and calibration failures become supporting evidence explaining why this stricter evaluation framework is necessary.

**PASS → STAGE 32**

---

## STAGE 32 — Internal Review / Submission Lock

Review the manuscript as a hostile external reviewer. Verify every headline against frozen artifacts, confirm that clone/donor/replicate counts are not pseudo-replication, ensure no surrogate is called independent fate, rerun reproducibility checks, and perform the final novelty/prior-art audit. Only corrections are allowed; no performance fishing.

**PASS → SUBMIT**

---

## STAGE 33 — Submission / Preprint

Submit the prospective CellFate-Rx manuscript and release the preprint/reproducibility package if appropriate. Reviewer-requested analyses are recorded as explicit follow-up work rather than silently changing the original evidence base.

**PASS → PUBLICATION / REVISION CYCLE**

---

# PARALLEL TRACK — Methods / Audit Paper

The existing project contains a coherent methodological audit: same-state ΔAge circularity, timepoint confounding, instrument disagreement, RES collapsing because uncertainty exceeds signal, calibration-target mismatch, scorecard bugs, and strong reproducibility discipline. That may support a separate methods/audit manuscript after a focused novelty/venue review. It should run in parallel only if it does not delay Stages 21–27.

---

# AFTER THE PROSPECTIVE PAPER — GENERATION 2 PRODUCT

After the first prospective paper, return to the larger experimental-triage product:

```text
multiple perturbation modalities
richer fate classes
dose optimisation
stronger held-out-treatment transfer
rejuvenation endpoints
future ΔAge with independent anchors
validated RES / integrated utility
```

Paper 1 does not need all of that.

It needs one strong new capability:

> **given a pre-intervention molecular state and candidate treatment, predict and rank independently observed future outcomes under biologically valid prospective evaluation.**

---

# One-line roadmap

```text
STAGE 21
PROSPECTIVE DATA QUALIFICATION
        ↓ PASS

STAGE 22
PROSPECTIVE BENCHMARK BUILD
        ↓ PASS

STAGE 23
LEARNABILITY + X×U INTERACTION GATE
        │
        ├── PASS ───────────────────────────────────────────────┐
        │                                                       │
        └── STAGE_24_BLOCKED_ROLE_A                             │
              ↓                                                 │
            STAGE 23.2                                         │
            ROLE-A RESOLUTION / FAILURE DECOMPOSITION           │
              ↓                                                 │
            FREEZE + AUDIT                                      │
            STAGE_23_2_ROLE_A_RESOLUTION_V1.md                  │
              │                                                 │
              ├── ROLE_A_CONFIRMATORY_SUPPORTED ────────────────┤
              │                                                 │
              ├── ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE         │
              │       ↓                                         │
              │     REMAIN BLOCKED                              │
              │                                                 │
              ├── ROLE_A_REDESIGN_REQUIRED                      │
              │       ↓                                         │
              │     VERSION + RE-GATE BENCHMARK; REMAIN BLOCKED │
              │                                                 │
              └── ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL         │
                      ↓                                         │
                    EXPLICIT ROADMAP / CLAIM REVISION            │
                    BEFORE ANY STAGE-24 DECISION                 │
                                                                │
        ┌───────────────────────────────────────────────────────┘
        ↓
STAGE 24
FREEZE + AUDIT STAGE_24_PROSPECTIVE_MODEL_V1.md
BUILD CELLFATE-RX PROSPECTIVE
        │
        ├── STAGE_24_MODEL_PASS
        │       ↓
        │     STAGE 25
        │
        └── STAGE_24_COMPLEXITY_NOT_EARNED
                ↓
              STOP / EXPLICIT REDESIGN OR
              SIMPLE-BASELINE ROADMAP DECISION

STAGE 25
STATE-CONDITIONED TREATMENT RANKING
        ↓ PASS

STAGE 26
HELD-OUT PERTURBATION / GENERALIZATION
        ↓ PASS OR SCOPED LIMIT

STAGE 27
INDEPENDENT BIOLOGICAL REPLICATION
        ↓ PASS

STAGE 28
CALIBRATION + OOD + DECISION UTILITY
        ↓ PASS

STAGE 29
FINAL EVIDENCE LOCK
        ↓ PASS

STAGE 30
SCIENTIFIC CLAIM LOCK
        ↓ PASS

STAGE 31
MAIN MANUSCRIPT
        ↓ PASS

STAGE 32
INTERNAL REVIEW / SUBMISSION LOCK
        ↓ PASS

STAGE 33
SUBMIT / PREPRINT / PUBLISH
```
