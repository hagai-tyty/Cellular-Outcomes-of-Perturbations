# stage_23_5_RECORD — Gen-1 Role-B ship plan: pre-freeze audit

## Goal
Audit `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md` before it is frozen, verify every quantity it quotes
against its source artifact, and repair any specification defect while no statistic has yet been
computed.

## Status
**FROZEN.** The plan remains **V1** — it is the first version the project believes will actually
ship, so it is frozen rather than superseded.

```text
  plan canonical-LF sha256   8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48
  section-12 checklist       38 / 38 PASS, walked once in full
  compute budget             ACCEPTED by the decision owner before freeze
  Stage 24                   OPEN under STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1
```

No Stage-25 ranking statistic exists. The 1,000-draw permutation run has **not** been started.

## Provenance
The plan was drafted externally and imported. It was **not** accepted as given: every number was
re-derived from the repository at commit `c52c7ac9a73abcf121624112d52f4f463db06f7a`.

## Verification — 28 of 28 quantities match exactly

```text
  COUNTS
    clones                     plan     1401   repo     1401   MATCH
    clone x treatment rows     plan     8406   repo     8406   MATCH
    nonzero rows               plan     2256   repo     2256   MATCH
    observed-zero rows         plan     6150   repo     6150   MATCH
    treatments                 plan        6   repo        6   MATCH
    outer folds                plan        5   repo        5   MATCH
    detected in >=1 condition  plan      929   repo      929   MATCH
    detected in all 6          plan       37   repo       37   MATCH
    RANKING-ELIGIBLE           plan      892   repo      892   MATCH
  
  METRICS
    W1 C1 pooled log loss      MATCH   0.47832325762628236
    W4 C1 pooled log loss      MATCH   0.46895102453285753
    W5 C1 pooled log loss      MATCH   0.4546510377185275
    W1 C2 clone-balanced MAE   MATCH   0.5144863827773782
    W4 C2 clone-balanced MAE   MATCH   0.5100798502288132
    W5 C2 clone-balanced MAE   MATCH   0.5007042703677924
    C2 W5 vs W4                MATCH   0.009375579861020777
    C2 W5 vs W1                MATCH   0.013782112409585823
    C1 W5 vs W4                MATCH   0.014299986814330035
    C1 W5 vs W1                MATCH   0.023672219907754866
    C1 W5vW4 null p95          MATCH   0.0008383045453672399
    C1 W5vW1 null p95          MATCH   0.0003983520753238711
    abundance / state ratio    MATCH   3.446751843310516
  
  PER-TREATMENT C1 INTERACTION
    Acid                       plan +0.0109716227542   repo +0.0109716227542   MATCH
    Cisplatin                  plan +0.0000245928453   repo +0.0000245928453   MATCH
    CoCl2                      plan +0.0189409151730   repo +0.0189409151730   MATCH
    Dabrafenib                 plan +0.0259209299498   repo +0.0259209299498   MATCH
    Doxorubicin                plan -0.0033224048416   repo -0.0033224048416   MATCH
    Trametinib                 plan +0.0332642650053   repo +0.0332642650053   MATCH
```

Sources: `results/stage22_wm989_clones.csv`, `results/stage22_wm989_clone_treatment.csv`,
`results/stage23_wm989_interaction_results.json`, `results/stage23_final_synthesis.json`.

The 892-clone ranking population was recomputed from the benchmark rather than taken from the
plan's arithmetic. `929 - 37 = 892` holds because eligibility is "at least one detected AND at
least one zero", and the 37 all-six clones are exactly those with no zero. Confirmed directly:
`((positives >= 1) & (positives < 6)).sum() == 892`.

Distribution of detected conditions per clone, for the record:

```text
  0 positives   472 clones      excluded, within-clone AUROC undefined
  1 positives   326
  2 positives   215
  3 positives   185
  4 positives   107
  5 positives    59
  6 positives    37 clones      excluded, within-clone AUROC undefined
                ---
  eligible      892
```

## Stage numbering — checked, no collision

The plan places the ranking test at Stage 25 and the scoped generalization limit at Stage 26. The
live roadmap already defines Stage 25 as *"State-Conditioned Treatment Ranking Challenge"* and
Stage 26 as *"Held-Out Perturbation / Generalization Challenge"*. The plan slots into the existing
numbering rather than inventing new stages.

## Four specification defects found and repaired

No incorrect number was found. Four defects were, all in the specification rather than the
evidence, and all repaired in place before any statistic was computed.

### 1. The permutation compute budget was unstated

Section 8.7 calls for 1,000 **full-refit** permutations and forbids early stopping, with no
indication of cost. The measured precedent makes it substantial:

```text
  Stage 23E, WM989 C1, 200 full-refit permutations    335.03 min   -> ~100 s/draw
  this plan, 1,000 draws with W4 and W5 refit each
      ~28 h single process
      ~19-20 h across 3 shards at the 1.44x effective speedup measured in Stage 23.2H
```

Added as §0.2, on the face of the document, to be accepted before freeze rather than discovered at
hour six. §0.2 and §8.7 additionally require that **each shard write its own cache file** and that
the merge assert all 1,000 indices are present before computing a statistic — Stage 23.2H lost a
completed draw to a race between three shards appending to one shared file, silently, caught only
by a completeness check.

### 2. The reproduction tolerance was named but never defined

§7 said "exact or tolerance-declared reproduction" and §8.2 said "within the frozen tolerance". No
tolerance appeared anywhere in the document. That is the same class of gap that let a broken power
instrument survive undetected through Stage 23.2H.

Added as §7.1: byte-identical OOF files as the primary requirement, a 1e-12 absolute bound on
pooled metrics as a fallback that also demands the cause of non-identity be named, and an explicit
input-integrity STOP otherwise.

### 3. Criterion 5 used an unpowered point estimate as a veto

§8.10 required `delta_TOP1 > 0` for `STAGE_25_RANKING_SUPPORTED`, while §8.8 stated the same
quantity was "not independently thresholded for significance". Those cannot both hold: a bare point
estimate with no uncertainty control was being allowed to veto a result that had cleared both a
bootstrap interval and a 1,000-draw permutation null. With 892 clones and a coarse `argmin`
statistic, `delta_TOP1` can fall below zero by chance while the ranking improvement is real.

Restated in §8.8 and §8.10 as an explicit **directional-consistency check** (`delta_TOP1 >= 0`),
with its bootstrap interval reported from the same clone resampling at no additional cost. It
remains deliberately conservative and may still produce a false negative; that is now stated rather
than concealed. It can withhold support, never grant it.

**This is the only change to a decision criterion.** It repairs an internal contradiction; it does
not lower a bar.

### 4. The bootstrap interval and the permutation null were presented as equivalent

The null of §8.7 refits the entire pipeline; the bootstrap of §8.6 does not. The interval is
therefore conditional on the fitted models and covers clone-sampling variability only. Standard
practice for a frozen-prediction comparison, but the two are different kinds of quantity. Declared
in §8.6.

## Two clarifications added

- **§2.2** now states that the five outer folds are the Stage-22 originals, loaded and never
  re-drawn, at the point where the cohort numbers appear.
- **§8.3** now spells out why `delta_RANK` and `delta_TOP1` use the scoring direction differently —
  discrimination versus selection utility — since the two read as contradictory otherwise.

The six decision criteria and the forbidden-actions list were hoisted into §0.1 so a reader who
stops before §8 still sees them. §11 and §12 were extended to cover the new frozen quantities.

## What was NOT changed

```text
  no number, threshold, metric, population, model, comparator, endpoint or verdict rule loosened
  no historical Stage-23 or Stage-23.2 verdict touched
  no new dataset authorized
  the two roadmap amendments remain separately named, not collapsed
  the plan remains V1 and remains DRAFT, PRE-EXECUTION
```

## Open items before freeze

```text
  the §0.2 compute budget must be explicitly accepted by the decision owner
  results/stage23_5_protocol.json must be created with the canonical-LF SHA-256
  results/stage23_5_handoff_to_stage24.json must be created
  the §12 checklist must be walked once, in full, and its result recorded
```

Until those are done the plan is audited but not frozen, and Stage 24 may not open.


---

# FREEZE — two corrections applied after the audit, then frozen

The decision owner accepted the §0.2 compute budget (1,000 full-refit permutations, three
independent shards, ~19–20 h, no early stopping) and required two corrections before freeze. Both
were applied; nothing else was reopened.

## Correction 1 — §7.1 is now a row-level reproduction gate

The pooled-metric fallback is **removed entirely**. Pooled metrics can agree to any precision while
individual rows differ in compensating directions, and the §8 ranking test never reads a pooled
metric — it is a function of within-clone orderings and nothing else. The fallback is replaced by
four requirements, all of which must hold:

```text
  R1  8,406 rows in each file; identical ordered (clone_id, treatment, outer_fold, y) tuples
  R2  |regenerated - frozen| <= 1e-12 for EVERY prediction cell of EVERY model column,
      checked cell by cell; a maximum, mean or aggregate is not a substitute
  R3  for all 1,401 clones, and for W4 and W5 independently, sign(s_iu - s_iv) unchanged for
      every within-clone condition pair, INCLUDING tie structure
  R4  the cause of non-identity named as a specific environment difference;
      "floating point" alone is not a cause
```

R3 is the load-bearing addition. A reproduction that preserved every pooled metric while flipping a
single near-tie would not have reproduced the input the ranking test consumes.

## Correction 2 — the `delta_TOP1` change is logged as a relaxation

The earlier record described the `delta_TOP1` change as repairing an internal contradiction "rather
than lowering a bar." That was not the whole truth, and §8.8 now separates the two components:

```text
  RECLASSIFICATION, not a relaxation
    calling it a directional-consistency check rather than a significance criterion.
    Removes a contradiction — an unpowered point estimate was vetoing a result that had
    cleared a bootstrap interval and a 1,000-draw null — and moves no boundary.

  RELAXATION, limited but real
    delta_TOP1 > 0  ->  delta_TOP1 >= 0
    Admits exactly one additional outcome: delta_TOP1 exactly zero, meaning W5 and W4
    select the same-quality top-1 condition on average. The statistic is a difference of
    two means of binary indicators over 892 clones, so it lives on a grid of 1/892 and
    exact zero is attainable, not a measure-zero curiosity.
```

Untouched: `delta_RANK > 0`, the bootstrap lower endpoint, the null-p95 comparison, `p_perm <= 0.05`
at 1,000 full-refit draws, and the primary metric, comparator, population, weighting, endpoint and
null construction. **Every primary threshold stands exactly as drafted.** The relaxation is confined
to the secondary consistency boundary and was made before any statistic existed.

## Section-12 checklist — 38 / 38

Walked once, in full, mechanically where mechanical was possible. Every quoted number was
re-derived from source; every textual commitment was asserted present in the frozen plan; the
absence of any Stage-25 artifact was checked on disk. Two extra items were added beyond the plan's
own list: that the stage numbering does not collide with the live roadmap, and that no historical
verdict was altered. The full item-by-item result is embedded in
`results/stage23_5_protocol.json` under `audit.items`.

One check failed on first run and was traced to the **check**, not the plan: an assertion had
assumed a line wrap that does not exist inside a fenced block. The assertion was corrected to be
whitespace-insensitive and the checklist re-run clean. The plan was not edited to satisfy it.

## Freeze artifacts

```text
  results/stage23_5_protocol.json              digest-bearing freeze artifact
  results/stage23_5_handoff_to_stage24.json    mechanical handoff
```

The plan cannot contain its own digest without recursion, so the protocol is authoritative for it.

## Stage 24 is open

```text
  STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1
```

Stage 24 is a bounded predictor-engineering stage: reproduce W1/W4/W5 under the §7.1 gate, build
and freeze the W5 tool, emit one frozen OOF prediction per clone-condition row, hand them to Stage
25. It may not inspect the ranking metric, replace W5 on the same folds, or add data.

## Not started

The §8.7 permutation run has **not** been started. It is ~19–20 h across three shards and is Stage
25's, not Stage 24's.
