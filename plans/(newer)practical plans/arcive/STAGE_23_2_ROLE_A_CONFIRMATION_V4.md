# STAGE 23.2 — ROLE-A CONFIRMATION PROTOCOL V4

**Frozen before any untouched confirmation evidence was inspected.**
Stage-23.2 protocol `78edd5d7f9900349925339169a5d5e3e5011fe23e3c7c22608ac98bfe3427bf4`.

**V1, V2 and V3 are preserved unchanged as historical.** None was executed. V4 carries every V3
section forward and adds three final mechanical corrections; the V3 and V2 change logs are retained
below for continuity.

```text
V4 change log

  section 18.2a AMENDED: the 18.1 failure case now says fewer than two qualifying
                     non-R1 BIOLOGICAL REPLICATES, never outcome units.

  sections 16.2, 19.4
                AMENDED: the per-replicate direction gate uses the SAME committed
                         five outer folds restricted to that replicate. Every qualifying
                         replicate must contribute >= 5 positive and >= 5 negative clones,
                         so every replicate-restricted outer test fold has both classes and
                         every outer-training split can run the frozen 3-fold inner CV.

  section 19   AMENDED: separates protocol freeze from label-dependent fold materialisation.
                         The construction rule is frozen before outcome inspection; the fold
                         table and feasibility checks occur only AFTER source-faithful outcome
                         reconstruction and BEFORE any model fitting/performance analysis.

  sections 1-15, 16.1, 16.3-18.1b are otherwise carried over from V3 unchanged.
```

V3's own change log, retained:

```text
V3 change log

  section 11   AMENDED: the positive-count criterion is a property of the COMBINED
               qualifying set, not of each candidate. An individual replicate does
               not need 140 positives.

  section 16.1 AMENDED: outer-fold construction now points at section 19 rather than
               describing the stratification informally.

  section 18.1 AMENDED: the replication gate is stated in BIOLOGICAL REPLICATES, not
               outcome units. Several libraries or selection units belonging to one
               biological replicate count as ONE replicate for this gate.

  section 19   NEW  the exact confirmatory outer-fold construction: fold count,
                    deterministic algorithm and seed, stratification by biological
                    replicate x outcome class, the feasibility conditions checked
                    before any outcome inspection, and the fail-closed rule.

  sections 1-10, 12-18 are otherwise carried over from V2 unchanged.
```

V2's own change log, retained:

V2 added one clarification, in response to a review point that V1 left a real gap: the `>= 140` positive-clone floor was derived
from *within-R1* simulations and V1 did not say what happens when confirmation evidence spans more
than one biological replicate. V1 could therefore have been read as licensing two things it never
intended — pooling separate outcome libraries to reach 140, or treating a single replicate's count
as the whole requirement.

```text
V2 change log

  section 10   amended: the floor is restated as a WITHIN-R1 planning floor that is
               NECESSARY BUT NOT SUFFICIENT, and explicitly does not establish power
               across independent biological replicates

  section 15   NEW  independent biological outcome units; per-unit source-faithful
                    outcome reconstruction; libraries are never pooled before selection

  section 16   NEW  how multiple independent replicates are combined for the
                    confirmatory analysis

  section 17   NEW  how the positive-clone floor applies across replicates

  section 18   NEW  exactly what replicate-level evidence is required for
                    ROLE_A_CONFIRMATORY_SUPPORTED, and therefore for Stage 24 to reopen

  sections 1-9, 11-14 are carried over from V1 unchanged.
```

This protocol may only be executed in 23.2G, and only on evidence that was **not** used to design
the correction. Executing it on the already-inspected Rewind biological replicate R1 is forbidden:
that analysis is exploratory by construction and is already recorded in 23.2C.

---

# 1. Confirmed scientific hypothesis

Under the historical Rewind outcome and the frozen Stage-22/23 evaluation geometry, pretreatment transcriptional state predicts the Role-A outcome beyond a DEPTH-COMPLETE nuisance baseline Bdepth = [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI), log1p(n_detected_GE_features_in_raw_pseudobulk)].

Mechanically derived: RESIDUAL_DEPTH_STRUCTURE is the only SUPPORTED mechanism that implies a correction. MODEL_SELECTION_NULL_INFLATION is UNRESOLVED, so search matching is not indicated; OUTCOME_LABEL_LIMITATION is UNRESOLVED, so V2 §8.7 forbids adopting an alternative outcome. The correction was not chosen by effect size.

# 2. Outcome definition

```text
historical hard label, top-100 gDNA barcodes with ties, unchanged
```

The historical hard label is retained deliberately. `OUTCOME_LABEL_LIMITATION` is UNRESOLVED, not
SUPPORTED, so V2 §8.7 forbids adopting an alternative outcome representation here. If a future
substage establishes the limitation, this protocol must be re-frozen and the material
benchmark-change firewall (V2 §10.7) applies.

# 3. Allowed input X

```text
pretreatment Gene Expression only, 36,601 features
clone pseudobulk: sum RAW counts, then CP10K and log1p exactly once
no lineage-assay feature, no clone identifier, no outcome-derived quantity
```

# 4. Nuisance block

```text
  log1p(n_pretreatment_cells)
  n_lanes
  log1p(total_raw_GE_UMI)
  log1p(n_detected_GE_features_in_raw_pseudobulk)
```

All continuous terms standardised training-only. This is the correction under test: the confirmation
asks whether state predicts outcome **beyond** this depth-complete baseline.

# 5. Model family, grid and preprocessing

```text
family        l2 logistic regression, liblinear, class_weight None
C grid        [0.01, 0.1, 1, 10]
K grid        [10, 20, 50]   (PCA fitted once at max K; K values are prefixes)
preprocessing training-only gene filter, gene scaler, PCA, PC scaler -- refitted inside every
              inner split
inner CV      3-fold StratifiedKFold(shuffle=True, random_state=23023)
```

Unchanged from the historical pipeline. The search path is **not** matched or narrowed:
`MODEL_SELECTION_NULL_INFLATION` is UNRESOLVED, so narrowing it is not an indicated correction.

# 6. Grouping unit

```text
clone. Outer folds grouped at clone level. Cells are never independent outcome replicates.
```

# 7. Primary metric and statistic

```text
metric      average precision at clone grain
statistic   delta_AP = AP(PCA(X) + Bdepth) - AP(Bdepth)
```

# 8. Permutation / null design

```text
whole clone-level expression profiles permuted as intact vectors
within outer-training and within outer-test separately
inside strata: n_pretreatment_cells {1,2,3+} x n_lanes
200 permutations, full nested-CV rerun per draw
p_perm = (1 + #{null >= observed}) / (n_perm + 1)
```

**Known conservatism, declared in advance.** 23.2C established that this null retains technical
similarity between donor and recipient profiles at Spearman ~0.34, which is why its centre is
positive. The design is retained unchanged so the confirmation is comparable to the historical
test; the depth-complete nuisance block is the correction, not a weakened null.

# 9. PASS threshold

```text
observed > null p95   AND   p_perm <= 0.05
```

Both required. No alternative threshold, one-sided relaxation or post-hoc adjustment is permitted.

# 10. Minimum positive-count / design requirement

```text
  scale 1:   35 positive clones -> power 0.290
  scale 2:   70 positive clones -> power 0.520
  scale 4:  140 positive clones -> power 0.940

  requirement: >= 140 positive clones at oracle AUC 0.66
```

the ladder is coarse (35 / 70 / 140 positive clones). The requirement is the smallest TESTED cohort reaching 0.80, not an interpolated value between rungs -- V2 §9.5 forbids extrapolating a precise required N.

A confirmation cohort below this floor may be run, but a null result from it is **not** evidence
against the hypothesis and must be reported as underpowered.

## 10.1 What this floor is, and what it is not — clarified in V2

```text
IS       a planning floor for TOTAL event count, derived from within-R1 simulation
IS NOT   a demonstration of power across independent biological replicates
IS NOT   a per-replicate requirement
IS NOT   a licence to pool outcome libraries, or to alter the frozen
         top-100-with-ties rule, in order to reach 140 positives
```

The 23.2E ladder was produced by resampling biological replicate R1's own covariates with
replacement. Every synthetic cohort on that ladder is still R1. The curve therefore estimates
**within-R1 event-count detectability under the empirical R1 covariate distribution** and nothing
else. Between-replicate variance is not modelled, and at `n_biological_replicates = 1` it is not
estimable at all.

Consequently the true power of a pooled multi-replicate confirmatory design is **unknown**, and may
be lower than the ladder suggests, because a multi-replicate design carries a source of variation
the simulation could not contain. The floor is therefore necessary but **not** sufficient, and
sections 15–18 carry the requirements it does not.

Reaching 140 positives by merging outcome libraries across replicates is specifically forbidden.
Doing so would change the outcome definition — a material benchmark semantic under Stage-23.2 V2
§10.7 — and would additionally destroy the replicate structure that is the entire purpose of
confirmation.

# 11. Source-qualification criteria

Applied to each candidate individually:

```text
pre-intervention molecular measurement present
independently measured later outcome
clone / lineage linkage sufficient for grouped evaluation
raw or processed data sufficient to reconstruct the endpoint
no same-state outcome leakage
outcome-unit structure establishable from source materials      (section 15.2)
a biological replicate INDEPENDENT of R1
```

The last criterion is what distinguishes confirmation from more of the same: R1 is consumed.

## 11.1 The positive-count criterion is a property of the SET — corrected in V3

V2 listed "at least the section-10 positive-count floor" among the per-candidate criteria, which
read as though **each** candidate replicate had to reach 140 positives on its own. That was wrong,
and it contradicted section 17, which already stated the floor applies to the total.

```text
CORRECT
    a candidate qualifies on the seven criteria above, WITHOUT any positive-count test.

    the >= 140 floor is then evaluated ONCE, on the COMBINED set of qualifying
    non-R1 replicates, after source-faithful per-unit reconstruction (sections 15, 17).

WRONG, and removed
    requiring any individual replicate to reach 140 positives
```

An individual replicate contributing, say, 35 positives is a perfectly good qualifying replicate.
What must reach 140 is the sum across the qualifying set. A per-candidate floor would have
disqualified every realistically available replicate and would have imposed a requirement the
23.2E simulations never tested.

The order of operations is fixed and may not be reversed:

```text
1. qualify each candidate on the seven criteria     (no performance quantity read)
2. reconstruct outcomes per unit                    (section 15)
3. sum positives across the qualifying set          (section 17)
4. evaluate the >= 140 floor once, on that sum
```

A candidate may **not** be added to or dropped from the qualifying set on the basis of how it moves
the total, or of any performance quantity. The set is fixed by step 1.

# 12. Search budget for confirmation data

```text
one bounded search, executed only after this protocol is committed
candidates enumerated from the 23.2A reserved ledger first
then a single external search pass
inclusion criteria may NOT be relaxed if results are sparse
```

# 13. Forbidden data — already inspected

```text
GSE227151 biological replicate R1  (GSM7092515, GSM7092516)
the Stage-22 Rewind benchmark and every Stage-23 / 23.2 artifact derived from it
the gDNA table stepThreeStarcodeShavedReads_BC_gDNA.txt
```

Reserved candidates `GSM7092517`-`GSM7092521` have been read at declared-metadata level only and
remain eligible. Whether they carry a reconstructable Role-A outcome is **UNVERIFIED** and must be
established under section 11 before any performance quantity is computed.

# 14. Stage-27 firewall

Any dataset used here becomes development/confirmation evidence and is **not** an untouched
Stage-27 replication set. Stage 27 must preserve an independent biological test of the eventual
frozen Stage-24 model.

---

# 15. Independent biological outcome units and per-unit reconstruction

## 15.1 Definition

```text
independent biological outcome unit
    the outcome assay belonging to ONE biological replicate -- the gDNA library, or the
    set of libraries, that the source study's own selection rule treats as a single
    selection unit for that replicate
```

For R1 this unit was one pooled library: `SampleNum = 3`, 49,554 rows, 1,936 distinct barcodes,
total support `N = 782,826`. That structure is a property of R1's data, **not** an assumption to
carry over. Each candidate replicate's unit structure must be established independently.

## 15.2 Establishing the unit structure — before any outcome value is read

The unit structure for a candidate replicate is determined from declared source metadata, author
code and file organisation only:

```text
which gDNA library or libraries correspond to that biological replicate
whether the author's rule was applied per library or across libraries
the selection-unit key the author's own code groups on
whether any barcode exclusion applies within that unit
```

This determination happens during source qualification (section 11) and completes **before** any
outcome value, barcode count or performance quantity is computed. Reading a library's declared
structure is qualification; reading its counts is inspection.

If the unit structure cannot be established from source materials, the candidate is
**DISQUALIFIED**. It is not adapted, guessed, or assumed to mirror R1.

## 15.3 Reconstruction is per unit, and source-faithful

Within each independent outcome unit, and independently of every other unit, apply the frozen rule
exactly as reconstructed for R1:

```text
group by (barcode, selection-unit key) -> sum support
slice_max(n = 100, with_ties = TRUE)
join to the retained clone set for that replicate
apply the same special-barcode exclusion the source rule applies
```

Each unit yields **its own** positive clone set. The following are forbidden without exception:

```text
pooling libraries from different biological replicates before selection
re-ranking barcodes across units
changing N from 100
changing or relaxing the tie rule
selecting a unit definition because it yields more positives
```

Any of these changes the outcome definition and triggers the Stage-23.2 V2 §10.7 material
benchmark-change firewall: version the benchmark and rerun the affected Stage-22/23 gates first.
None of them may be used to reach the section-10 floor.

## 15.4 Recorded per unit

```text
replicate identity and accession(s)
selection-unit key and its source justification
rows, distinct barcodes, total support N
rank-100 cutoff value and tie size
selected barcode count (>= 100 when ties expand it)
positive clone count
retained clone count
```

---

# 16. Combining independent replicates for the confirmatory analysis

Frozen now, before any reserved matrix or outcome value is inspected.

## 16.1 Primary analysis — pooled, with replicate as a blocking nuisance

```text
1. reconstruct outcomes per unit (section 15). Never pool before this point.
2. pool the resulting clone-level rows across qualifying replicates.
3. replicate identity enters as a NUISANCE covariate, alongside Bdepth.
4. clone-level outer folds, constructed exactly as frozen in section 19.
5. the permutation null permutes whole expression profiles within
   replicate x stratum -- preserving replicate structure under the null exactly as the
   historical null preserved depth structure.
6. everything else -- grid, inner CV, metric, statistic, PASS threshold -- is
   unchanged from sections 5-9.
```

Replicate identity is a nuisance term only. It may **not** be a predictor of interest, and it may
**not** be interacted with `X`: a state × replicate interaction is a different scientific claim and
is out of scope for this confirmation.

## 16.2 Mandatory secondary — per-replicate direction

Run the frozen comparison separately within each qualifying replicate and report, for each:

```text
delta_AP with its bootstrap CI
positive clone count
an explicit underpowered flag when that replicate's positives fall below the
section-10 floor on their own -- which will usually be the case
```

These per-replicate results are **reported and gating on direction only** (section 18.5). They are
not required to reach significance individually; requiring that would impose a per-replicate power
floor the 23.2E simulations never tested.

Each per-replicate comparison uses the **same committed five outer-fold assignments from section 19,
restricted to that replicate**. No replicate-specific folds are redrawn. Therefore every qualifying
replicate must satisfy the per-replicate geometry minimum in section 19.4: at least 5 positive and
5 negative clones. With the frozen round-robin construction, that guarantees every replicate-level
outer test fold contains both classes and every replicate-level outer-training set contains at least
4 positives and 4 negatives, which is sufficient for the unchanged 3-fold inner StratifiedKFold.

## 16.3 Mandatory secondary — cross-replicate transfer

Train on one qualifying replicate and evaluate on a held-out one, both directions where two or more
qualify. Report `ΔAP` per direction.

This is **reported, not gating**. With two replicates it is a single split and far too noisy to
carry a verdict, but a transfer that is strongly negative in both directions while the pooled test
passes is a signal that the pooled result is driven by within-replicate structure, and must be
stated prominently if it occurs.

## 16.4 Forbidden combinations

```text
pooling outcome libraries before per-unit selection            (section 15.3)
treating replicates as exchangeable rows with no blocking term
permuting across replicate boundaries in the null
dropping a qualifying replicate because it weakens the result
adding a replicate after seeing its effect on the pooled statistic
```

The set of qualifying replicates is fixed by section 11 qualification **before** any performance
quantity is computed, and is not revised afterwards.

---

# 17. How the positive-clone floor applies across replicates

```text
the >= 140 floor applies to the TOTAL positive-clone count summed across all
qualifying independent replicates, after per-unit reconstruction.

it is NOT a per-replicate floor.
it is NECESSARY, NOT SUFFICIENT.
```

**Why total and not per-replicate.** The 23.2E ladder varied total event count; it never tested a
design in which each replicate independently carried 140 positives. Imposing a per-replicate floor
would invent a requirement the simulations never examined, and would almost certainly be
unsatisfiable.

**Why necessary but not sufficient.** The ladder is a within-R1 curve (section 10.1). A pooled
multi-replicate cohort carries between-replicate variance that the simulation could not model, so
meeting the total floor does not establish that the confirmatory design has 0.80 power. It
establishes only that the design is not disqualified on event count alone.

**If the total falls below 140.** The analysis may still be run and must still be reported in full.
It **cannot** emit `ROLE_A_CONFIRMATORY_SUPPORTED`. A positive result is recorded as *underpowered
supporting evidence*; a null result is **not** evidence against the hypothesis. The exit is
`ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE` and Stage 24 stays blocked.

**What may not be done to reach 140.**

```text
pooling outcome libraries across replicates
enlarging N beyond 100, or loosening the tie rule
adding a replicate that failed section 11 qualification
counting R1 positives toward the total
```

R1 is consumed by diagnosis and its 35 positives never count toward this floor.

---

# 18. Replicate-level evidence required for `ROLE_A_CONFIRMATORY_SUPPORTED`

All six conditions are required. Each is mechanical.

```text
18.1  >= 2 independent BIOLOGICAL REPLICATES qualify under sections 11 and 15,
      and NONE of them is R1.
      Several libraries or selection units belonging to one biological replicate
      count as ONE replicate for this gate (see 18.1b).

18.2  each qualifying unit's outcome was reconstructed source-faithfully and
      independently, with no pooling before selection and no change to the frozen
      top-100-with-ties rule.

18.3  the TOTAL positive-clone count across qualifying units is >= 140.

18.4  the pooled primary test of section 16.1 passes BOTH frozen gates:
          observed > null p95   AND   p_perm <= 0.05

18.5  delta_AP is POSITIVE in every qualifying replicate analysed separately
      (section 16.2). Individual significance is not required; a negative direction
      in any qualifying replicate blocks the verdict.

18.6  no material benchmark semantic was changed. If reconstruction required one,
      Stage-23.2 V2 §10.7 applies and the benchmark must be versioned and re-gated
      before confirmation, not after.
```

## 18.1b Replicates, not units — corrected in V3

V2 phrased this gate in *outcome units*. That was too weak: a single biological replicate whose
gDNA happened to be split across two libraries would have presented as two units and satisfied a
gate whose entire purpose is to establish biological replication.

```text
counting rule for gate 18.1

    the unit of counting is the BIOLOGICAL REPLICATE
    several libraries / selection units within one replicate  ->  ONE replicate
    a replicate contributes at most 1 toward the ">= 2" requirement
```

Sections 15 and 16 continue to operate on **outcome units**, because reconstruction and pooling are
per-library operations. Only the replication gate counts replicates. The two are deliberately
different, and the qualification record must state, for every qualifying replicate, which units it
comprises.

If the biological-replicate identity of a candidate's libraries cannot be established from source
materials, those libraries cannot be counted toward gate 18.1 at all — the same disqualification
rule as section 15.2, applied to replicate identity rather than unit structure.

## 18.1a Why two replicates, and why neither may be R1

Stage-23.2 V2 §9.5.2 clears `BIOLOGICAL_REPLICATION_LIMITATION` only when the claim is supported by
two or more independent biological replicates — which is exactly why gate 18.1 counts replicates
rather than units. R1 cannot be one of them, for two independent reasons:

```text
1. R1 is consumed. It designed the correction, so V2 §11.1 and the Stage-27 firewall
   exclude it from confirmation evidence.

2. R1 does not support the claim. The corrected same-data diagnostic on R1 came out
   NEGATIVE (O11 +0.00872 vs q95 +0.00835, p_diag 0.0547). A replicate whose corrected
   analysis failed cannot be counted as one of the replicates supporting the claim.
```

Reason 2 is the stronger one and is easy to overlook: it is not merely that R1 is procedurally
disqualified, but that its own corrected result did not support the hypothesis.

## 18.2a What each outcome means

```text
all six satisfied
    ROLE_A_CONFIRMATORY_SUPPORTED
    BIOLOGICAL_REPLICATION_LIMITATION -> NOT_SUPPORTED
    Stage 24 may reopen, subject to a complete benchmark-compatible handoff

18.3 fails (total positives < 140), everything else satisfied
    ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
    the result is recorded as underpowered supporting evidence
    Stage 24 stays blocked

18.1 fails (fewer than two qualifying non-R1 biological replicates)
    ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
    BIOLOGICAL_REPLICATION_LIMITATION remains SUPPORTED
    Stage 24 stays blocked -- a single replicate cannot clear it

18.4 or 18.5 fails
    ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
    Stage 24 stays blocked. A failed confirmation at an adequate event count is
    informative and must be reported as such, but per section 10 it is not
    evidence of absence unless the design met the floor.

18.2 or 18.6 fails
    ROLE_A_REDESIGN_REQUIRED
    the benchmark or the reconstruction is inadequate for the intended claim;
    version, re-gate, and return
```

## 18.3a Declared feasibility risk

This is recorded now, before any reserved matrix is inspected, so that it cannot later look like a
post-hoc excuse.

R1 yielded **35 positive clones from 3,147 retained clones**. If comparable replicates yield a
comparable rate, reaching a **total of 140 non-R1 positives** would require roughly four R1-sized
biological replicates. The 23.2A reserved ledger declares three non-R1 `hiFT` samples, of which two
(`GSM7092520`, `GSM7092521`) are sorted for cycling status and therefore may fail the
"same scientific claim" requirement of section 11 outright.

It is therefore a realistic and legitimate outcome that **the reserved evidence cannot satisfy this
protocol**, and that 23.2G exits at `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE` with Stage 24 still
blocked. That is a finding about what the available data can support, not a failure of the stage,
and it must not be worked around by relaxing sections 15, 17 or 18.

Whether the reserved replicates carry a reconstructable Role-A outcome at all remains
**UNVERIFIED**. Nothing in this section is based on inspecting them.

---

# 19. Confirmatory outer-fold construction — rule frozen before outcome inspection; folds materialized after label reconstruction

The historical Role-A folds were frozen in Stage 22 and simply reused by every later substage.
Confirmation evidence has no such pre-existing folds. The **construction rule** is therefore frozen
here before any untouched outcome value is inspected, but the fold table itself necessarily depends
on reconstructed outcome class `y`. It is materialized only after source-faithful outcome
reconstruction and before any model fitting or performance analysis.

## 19.0 Frozen execution order

```text
1. qualify each candidate on section-11 source criteria using declared metadata,
   author code and file organisation; fix the qualifying replicate set
2. reconstruct each qualifying outcome unit source-faithfully under section 15
3. evaluate the section-17 total-positive floor and the section-19.4 geometry conditions
4. construct and commit the deterministic outer-fold table + SHA-256
5. only then begin expression-model fitting, null/permutation runs or any performance analysis
```

Outcome labels are therefore permitted at steps 2-4 because they are required to reconstruct the
frozen endpoint and stratified folds. They may not be used to alter the qualifying set, model/grid,
threshold, nuisance block, outcome rule or fold algorithm. No expression value or model performance
quantity is consulted before step 5.

## 19.1 Geometry — unchanged from the historical pipeline

```text
outer folds        5
grain              clone
inner CV           3-fold StratifiedKFold(shuffle=True, random_state=23023)
```

Identical to sections 5 and 6. Nothing about the nested-CV geometry is re-tuned for confirmation.

## 19.2 Deterministic construction

```text
seed                       23451     (SEED_CONFIRMATION_FOLDS, frozen here)
stratification cells       (biological replicate) x (outcome class y)

algorithm, per stratification cell:
    1. take the cell's clone ids, sorted lexicographically      (order-independent input)
    2. permute them with numpy default_rng(23451)               (one generator, cells
                                                                 visited in sorted cell order)
    3. deal the permuted ids round-robin into outer folds 0..4
```

Dealing round-robin rather than slicing guarantees each fold receives `floor(n/5)` or `ceil(n/5)`
members of every cell, so both biological replicate and outcome class are balanced across folds by
construction rather than by chance.

The construction consumes only clone ids, biological-replicate identity and the reconstructed
outcome class. It reads **no** expression value and **no** performance quantity, and it is computed
once, before any model is fitted.

## 19.3 Frozen before use

```text
results/stage23_2/stage23_2g_confirmation_folds.csv     clone_id, replicate, y, outer_fold
```

committed together with its SHA-256 **before** any model is fitted. Every confirmatory statistic
pins that digest. A fold table that does not match the committed digest invalidates the run.

## 19.4 Feasibility conditions — checked after source-faithful outcome reconstruction and before model fitting/performance analysis

F1-F3 are the mathematical geometry conditions. F4 is the separate confirmation event-count floor
checked at the same pre-performance checkpoint.

```text
F1  every qualifying biological replicate contributes >= 5 positive AND >= 5 negative clones
    (the minimum needed for the same five committed outer folds, restricted to that replicate,
     to place both classes in every replicate-level outer test fold)

F2  every outer fold receives >= 1 positive and >= 1 negative clone BOTH:
      - in the pooled confirmation cohort, and
      - when restricted to each qualifying biological replicate
    (otherwise the corresponding AP / delta_AP contrast is not defined as frozen)

F3  every outer-TRAINING set contains >= 3 positive AND >= 3 negative clones BOTH:
      - in the pooled confirmation cohort, and
      - when restricted to each qualifying biological replicate
    (otherwise StratifiedKFold(n_splits=3) is not constructible for the frozen inner CV)

F4  total positive clones across the qualifying set >= 140          (section 17)
```

With the section-19.2 round-robin assignment, F1 implies at least one positive and one negative from
each qualifying replicate in every outer test fold, leaving at least four of each class in every
replicate-restricted outer-training set. F2 and F3 are still asserted explicitly as contracts rather
than inferred silently.

F1-F3 are checked on reconstructed labels after section-15 outcome reconstruction and before any
model fitting or performance quantity. If F1-F3 fail, the geometry is infeasible and section 19.5
applies. F4 does **not** determine mathematical constructibility: if F4 alone fails, section 17's
underpowered-confirmation rule applies and the analysis may still run, but it cannot emit
`ROLE_A_CONFIRMATORY_SUPPORTED`.

## 19.5 Fail-closed rule

If F1, F2 or F3 fails after source-faithful outcome reconstruction, the frozen geometry is not constructible for the available evidence:

```text
emit    CONFIRMATION_GEOMETRY_INFEASIBLE
exit    ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
Stage 24 stays blocked
```

**Forbidden responses to an infeasible geometry:**

```text
reducing the outer-fold count below 5
changing the inner CV scheme, its fold count or its seed
merging biological replicates to manufacture a workable stratification
dropping a qualifying replicate to make the geometry fit
pooling outcome libraries                                  (section 15.3)
re-drawing folds with a different seed
constructing separate replicate-specific folds instead of restricting the committed five-fold table
```

The geometry is part of what is being confirmed. A confirmation run under a geometry chosen to
accommodate the data would not be a confirmation of the historical result.

## 19.6 No post-hoc revision

```text
the fold count, the CV scheme, the seed, the stratification cells and the rule of
restricting the committed fold table for per-replicate analyses are fixed by this section
and may NOT be changed after any performance quantity has been seen.
```

If a fault in the construction is discovered after performance is observed, the correct response is
to record it, stop, and issue a new protocol version — not to re-draw folds inside the current run.
This is the same discipline sections 15 and 17 apply to the outcome rule and the floor.
