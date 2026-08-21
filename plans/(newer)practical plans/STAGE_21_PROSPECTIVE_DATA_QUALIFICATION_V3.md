# STAGE 21 — Prospective Data Qualification v3

**Status:** PRE-REGISTERED PLAN — 🔒 frozen before the public-dataset search begins.
**Supersedes:** `arcive/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md`, into which the post-21B
execution amendment is now folded. V2 is archived unchanged; nothing in it was edited in place.
**Executed so far:** 21A `f6a0056` · 21B `3b8b644` · 21C `cab27d8` — all frozen. **This file does not reopen them.**
**Current model:** `_s16` remains frozen.
**No model training in this stage. No `src/` changes.**

> Section numbers changed between V2 and V3. The 21A record cites "plan §§2, 5, 6, 11, 19, 20" and
> the 21B record cites the plan as a whole — **those citations resolve against the archived V2**,
> which is why V2 is kept verbatim rather than deleted.

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

## 0.1 What V3 changes relative to V2 — the complete list

V2 was written while the local geometry was still open. 21A and 21B closed it. Exactly three things
change; everything else below is V2's text carried over.

| # | change | why |
|---|---|---|
| 1 | **Phases renumbered** — V2 had 21A local → 21B *public* → 21C download. 21B was actually spent resolving the *local* questions, so public qualification becomes **21C** and download/reconstruct becomes **21D**. | Otherwise two different things are called 21B and the record is ambiguous. |
| 2 | **Role B is no longer blocking.** V2 §15 required `QUALIFIED_ROLE_A` **AND** `QUALIFIED_ROLE_B` to PASS, and §16 sent a single-role result to `PARTIAL_DATA`. V3 replaces that with the `FULL_DATA_PATH` / `CORE_DATA_PATH` / `DATA_BLOCKED` rule in §9.8. | A missing Role B should scope the treatment-ranking contribution down, not hold the paper hostage indefinitely. |
| 3 | **A frozen search budget** (§9.2) — two passes, six named candidate families, then a hard stop. V2 named candidate families but set no stopping rule. | Pre-registering the stopping rule before searching is what stops the search from becoming an open-ended hunt that terminates whenever a result happens to look good. |

Not changed: the linkage classes, the outcome-truth classes, the resolvability arithmetic, the hard
qualification criteria Q1–Q7, the tri-state evidence standard, the leakage guards, the
source-study reproduction gate, and the 21A/21B verdicts.

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

# 2. Stage 21 has four phases

```text
21A — LOCAL DATA AUDIT                          ✅ DONE   (f6a0056)
        ↓
21B — LOCAL SOURCE / DESIGN RESOLUTION          ✅ DONE   (3b8b644)
        ↓
21C — PUBLIC PROSPECTIVE DATASET QUALIFICATION   ✅ DONE   (cab27d8)
        ↓            -> FULL_DATA_PATH
                        Role A = GSE227151 (Rewind)
                        Role B = GSE279162
21D — ACQUIRE + RECONSTRUCT + VERIFY
        ↓
PASS
READY FOR STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

Stage 21 does **not** fit predictive models.

---

# 3. Stage 21A — Local data audit ✅ DONE

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

## 3.1 Result as executed

```text
VERDICT: CULTURE_FORWARD_AVAILABLE
```

Two questions were left explicitly open rather than guessed — `GSE242423`'s lineage status, and
whether `GSE165176`'s SSEA4/CD13 sorting is a *future* outcome — both carried as
`UNKNOWN_REQUIRES_SOURCE_AUDIT`, which is a request for one more file, not a rejection.

Full table, the four parsing bugs it exposed, and the evidence trail: `RECORDs/stage_21A_RECORD.md`.

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

**This section is the pre-registered prediction and is left exactly as it was written before
execution.** What actually happened is recorded separately in §7.1. Do not edit this block.

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

## 7.1 What the audit actually found — added after execution

```text
GSE242423   prediction HELD, and is now proven rather than expected
            0 lineage/clone/CellTag/LARRY/hashtag tokens in 50,648 chars of
            Series Matrix + MINiML; all 24 barcode contexts are 10x vocabulary;
            the extract protocol trypsinises at collection -> destructive sampling
            n_independent_units = 1

GSE165176   prediction UNDERSOLD the metadata and OVERSOLD the geometry
            the sort metadata is better than expected: the SSEA4/CD13 antibody
            sorts are an orthogonal phenotype/state call measured independently
            of RNA -- the only non-RNA readout anywhere in the project
            but the EARLY -> LATER mapping is not merely "not yet established",
            it is absent: 47 of 71 cultures yield BOTH fractions, 118 of 124
            samples are already sorted, no FACS proportions exist, and day 54
            carries only SSEA4 -> no outcome contrast
            effective n = 6 donors, not the 12 that 21A reported
```

The 21A figure of `12` is left frozen in 21A's own result file and record, with the correction
recorded and pinned by test rather than silently overwritten.

---

# 8. Stage 21B — Local source/design resolution ✅ DONE

**Both local prospective routes are closed on evidence.**

```text
GSE242423 -> LINEAGE_ABSENT_PROVEN
GSE165176 -> ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY
```

Neither local dataset can pose a prospective `X_before + U -> Y_future` task. `GSE242423` is one
unlinked, destructively-sampled trajectory. `GSE165176` carries a genuinely orthogonal antibody
readout, but it is a same-timepoint subpopulation split with no proportions, no unsorted early
population beyond day 0, and no terminal outcome variation.

**Nothing here is `UNKNOWN`. No further download changes either answer.** That is what moves the
critical path onto public data. Evidence trail: `RECORDs/stage_21B_RECORD.md`.

---

# 9. Stage 21C — Public prospective dataset qualification

**No modelling. No large reconstruction. No `src/` changes.** The only purpose is to decide which
public datasets genuinely satisfy the required prospective geometry.

A dataset is accepted on experimental geometry, **not** on whether a model later performs well.

## 9.1 The two roles

### ROLE A — core reprogramming prospective dataset (MANDATORY)

Purpose:

> Show that the prospective framework applies to the biological problem CellFate-Rx was originally
> built for.

```text
X_before  +  reprogramming intervention U  ->  independently measured future outcome Y
```

Preferred geometry:

```text
pre-reprogramming fibroblast state
+
OKSM / reprogramming protocol
->
later independently measured reprogramming success / failure
```

Desired: pre-treatment RNA · clone/lineage linkage · orthogonal future success/failure ·
enough independent biological units · reconstructable public files.

### ROLE B — multi-perturbation interaction dataset (HIGH VALUE, NOT BLOCKING)

Purpose:

> Prove that CellFate-Rx can do something more interesting than predict whether a starting state is
> generally "good" or whether a treatment is generally "strong."

```text
X_before  +  multiple treatments U1, U2, …  ->  independent future response Y
```

Ideally the same clone / sister population meets several treatments, so that Stage 25 can later
ask *does treatment preference depend on starting molecular state?*

**If Role B is not found inside the frozen budget, the paper is NOT held.** Scope the
treatment-ranking contribution down and proceed on a strong Role A. This replaces V2 §§15–16.

No candidate is accepted from reputation or abstract alone.

## 9.2 Frozen search budget — exactly two passes

### PASS 1 — six named candidate families

| role | family | accession(s) |
|---|---|---|
| A | Rewind | `GSE227151` / `GSE243933` |
| A | CellTag | `GSE99915` |
| A | CellTag-multi | `GSE216518` / `GSE216521` |
| B | ReSisTrace | `GSE223003` |
| B | multi-treatment melanoma lineage | `GSE279162` |
| B | sequential-treatment / barcode | `GSE253739` |

Accessions and geometry are to be **verified against GEO / the source papers**, not inherited from
any summary in this repository — including summaries written by me.

### PASS 2 — limited expansion

Runs **only for a role still unresolved after Pass 1**. Maximum **3 additional serious candidates
per unresolved role**. A candidate counts toward Pass 2 only if its abstract/methods already make a
plausible case for pre-outcome molecular state **and** real biological linkage **and** a future
outcome.

### After Pass 2

```text
HARD STOP
```

No third pass. No open-ended "one more dataset" loop.

## 9.3 Hard qualification criteria

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

### Additional criterion for the multi-perturbation role

A dataset counts as the main **Role B** dataset only if the perturbation variable contains
scientifically meaningful variation.

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

## 9.4 Cheap rejection cascade — the order those criteria are applied in

Applied in order, so an obvious failure is found before anything large is downloaded.

```text
Q1. Is X measured before Y?                       NO -> reject
Q2. Real biological linkage early X -> later Y?   NO -> reject
Q3. Is Y independently measured, not f(X_future)? NO -> SURROGATE_ONLY
Q4. Can public files reconstruct the X->Y link?   NO -> MISSING_REQUIRED_FILE / reject as primary
Q5. Enough independent biological units?          NO -> PILOT_OR_REPLICATION_ONLY
Q6. Does treatment vary meaningfully?             required for ROLE B
```

**Do not download large raw data to discover a Q1/Q2/Q3 failure.**

## 9.5 Evidence standard — carried forward from 21A/21B

Every important field keeps `value` / `status` / `evidence`, with the status one of:

```text
PRESENT
ABSENT_PROVEN
UNKNOWN_REQUIRES_SOURCE_FILE
```

**"Not found" is never turned into "absent".**

## 9.6 Mandatory qualification table

Per candidate:

```text
dataset · accession · paper · role_A_candidate · role_B_candidate
biological_system · species · cell_type

pre_outcome_rna_available · time_of_X

lineage_or_clone_link · sister_cell_link · linkage_method

future_outcome · outcome_measurement_method · outcome_orthogonal_to_rna

treatments · n_treatments · dose_available · exposure_available
same_clone_or_related_population_across_treatments
multi_treatment_interaction_test_possible

n_independent_experiments · n_independent_clones_or_lineages · n_positive · n_negative
outer_split_unit · minimum_attainable_p · resolvable

processed_files_available · lineage_mapping_available · outcome_mapping_available
join_key_reconstructable · reconstructable

missing_files · evidence · verdict
```

## 9.7 Qualification verdicts

```text
QUALIFIED_ROLE_A · QUALIFIED_ROLE_B · QUALIFIED_BOTH
QUALIFIED_SURROGATE · QUALIFIED_REPLICATION_ONLY
NOT_QUALIFIED · UNKNOWN_REQUIRES_SOURCE_FILE
```

`QUALIFIED_SURROGATE` = excellent forward lineage geometry, but the outcome is still
expression-derived — useful for development, not the independent-fate headline.

**Role A qualifies only with all of:** pre-outcome X · valid biological linkage · orthogonal future
Y · reconstructable public data · resolvable independent units.

**Role B additionally requires** perturbation variation sufficient to later test `X × U`.

## 9.8 Stage 21C decision rule

```text
FULL_DATA_PATH   Role A qualified AND Role B qualified
CORE_DATA_PATH   Role A qualified, Role B not found inside the budget
                 -> proceed; treatment-ranking contribution is SCOPED DOWN
DATA_BLOCKED     Role A not found inside the budget
                 -> do NOT build the prospective architecture yet
```

If Role A qualifies and Role B does not, **do not keep searching**. The main paper continues as a
prospective reprogramming prediction paper; Stage 25 state-conditioned multi-treatment ranking
becomes optional.

## 9.9 Artifacts

```text
experiments/diag_stage21c_public_qualification.py
tests/test_diag_stage21c_public_qualification.py
results/diag_stage21c_public_qualification_results.json
plans/(newer)practical plans/RECORDs/stage_21C_RECORD.md
```

The record stays compact — goal · candidates checked · source files/pages inspected · files
added/modified · bugs/corrections · final qualification table · final verdict · what this proves ·
what it does not prove · next action. **An audit record, not a diary.**

## 9.10 Engineering constraints

```text
src/ byte-unchanged          no predictive modelling      no target tuning
_s16 unchanged               no logistic regression       no outcome fishing
no architecture work         no neural network            no change to 21A/21B frozen verdicts
```

Full relevant tests and CI-scope ruff at the end.

---

# 10. Optional DATA ROLE C — Locked replication candidate

If possible, identify a third dataset during Stage 21 and reserve it.

Do not use it to choose model architecture or hyperparameters.

Its purpose is a later locked external replication.

This is strongly preferred but not required to begin Stage 22.

---

# 11. Stage 21D — Acquire and reconstruct

**Do not automatically download everything when 21C ends.** Report the qualification result first.
Only once the best Role-A (and, if available, Role-B) dataset(s) are frozen does 21D open.

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

# 12. Standard table schema

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

# 13. Source-study reproduction gate

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

# 14. Leakage guards

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

# 15. Stage 21 final PASS

Stage 21 passes when:

```text
1.  Role A is qualified and reconstructed.
2.  Role A uses true pre-outcome molecular state.
3.  Role A has valid biological linkage.
4.  Its primary future outcome is orthogonal to RNA.
5.  It is statistically resolvable for its intended benchmark role.
6.  Source-study reconstruction checks pass.
7.  A standardized frozen table exists.
8.  Leakage guards pass.
9.  src/ is unchanged.
10. The 21C verdict is FULL_DATA_PATH or CORE_DATA_PATH, recorded explicitly.
```

If the verdict is `FULL_DATA_PATH`, criteria 1–8 apply to Role B as well.

Preferred extra:

```text
a third locked replication candidate is identified.
```

**Changed from V2:** V2 required *both* roles reconstructed to PASS. A `CORE_DATA_PATH` PASS is now
legitimate, and the scope-down is declared in the record rather than left implicit.

---

# 16. Stage 21 output

Example:

```text
STAGE 21 — PROSPECTIVE DATA QUALIFICATION
=========================================

LOCAL DATA
final paper geometry available locally: NO   (21A/21B, both routes closed)

ROLE A — REPROGRAMMING
dataset: ...
linkage: ...
outcome: ...
orthogonal: YES
resolvable: YES
reconstruction: PASS

ROLE B — MULTI-PERTURBATION
dataset: ...            (or: NOT FOUND INSIDE BUDGET -> scoped down)
treatments: ...
same-clone multi-U: YES/NO
outcome: ...
orthogonal: YES
resolvable: YES
reconstruction: PASS

ROLE C — LOCKED REPLICATION
candidate: ...

21C VERDICT: FULL_DATA_PATH | CORE_DATA_PATH | DATA_BLOCKED
VERDICT: PASS
NEXT: STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

---

# 17. What Stage 21 proves

A PASS proves only:

> **The project now possesses the datasets required to test a genuinely prospective, treatment-conditioned future-outcome model.**

It does not prove the model works.

---

# 18. What Stage 21 explicitly does NOT do

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

# 19. Next stage

```text
STAGE 22 — PROSPECTIVE BENCHMARK BUILD
```

Stage 22 turns the qualified datasets into a standardized benchmark with frozen tasks, splits, metrics and leakage controls.

Only after that benchmark exists do we ask whether the prospective problem is learnable.

---

# 20. One-sentence definition

> **Stage 21 secures the prospective data a non-trivial CellFate-Rx paper needs — a lineage-grounded
> reprogramming anchor with an independently measured future outcome and, if the frozen search
> budget finds one, a multi-perturbation dataset capable of testing state-dependent treatment
> response — with reconstructable public ground truth in both cases.**
