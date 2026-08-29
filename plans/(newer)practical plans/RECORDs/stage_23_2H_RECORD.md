# stage_23_2H_RECORD — Role-A independent confirmation on biological replicates 2 and 3

## Goal
Execute `STAGE_23_2_ROLE_A_CONFIRMATION_V5.md` end to end and resolve the Stage-24 opening question.

## Why V5 existed at all

Stage 24 opens on one condition, from the roadmap:

> a pre-registered confirmatory analysis succeeds on evidence that was not used to design the
> Stage-23.2 correction, **and** the handoff contract is complete and benchmark-compatible.

Two things stood between us and attempting it, and V4 could not resolve either.

**The blocker of record was wrong.** 23.2G concluded `QUALIFYING_SET_EMPTY_FROM_FROZEN_SEARCH_SPACE`
— replicates 2 and 3 had a pretreatment transcriptome and no later outcome at all. That was a true
answer to the wrong question: the qualification enumerated GEO supplementary files, and those
replicates' outcome materials were never deposited in GEO. `stage_23_2G_RECORD.md` stands unedited;
`stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md` records the correction.

**The real blocker was V4 §18.3** — `>= 140` total positive clones, against 76 available. Plus a
second, structural problem: V4 §15.3 demanded a *source-faithful* reconstruction and in the same
breath mandated `slice_max(n = 100, with_ties = TRUE)`. That is coherent only while R1 is the only
replicate, because top-100-with-ties **is** R1's source rule. It is not replicate 2's and not
replicate 3's.

## What V5 changed, and why each change is not a relaxation

```text
§6   R1's top-100-with-ties  ->  each replicate's OWN source rule
§7   the author spike-in indexing bug: corrected primary, author-faithful sensitivity
§9   the imported >= 140 floor  ->  a MEASURED power gate on the realized geometry
§4   eligibility frozen to ungated samples only
§8   a new, separately versioned cohort with its own folds and strata
```

**§6.** V4 froze R1's *instance* of a principle in the place where the *principle* belonged, at a
time when the two were indistinguishable. Every anti-gaming prohibition in §15.3 survives verbatim:
no pooling across replicates, no re-ranking across units, no selecting a unit definition because it
yields more positives. And the rules were fixed by reproducing the authors' serialized objects
exactly, before any statistic existed — they cannot have been tuned to a result.

**§9.** V4 §10.1 states in its own words that the 23.2E ladder measures *"within-R1 event-count
detectability under the empirical R1 covariate distribution and nothing else."* The confirmation
cohort is not R1:

```text
                       clones   positives   prevalence
  R1  (Stage 22)         3147          35        1.11%
  confirmation cohort    2310          75        3.25%
```

Power at fixed oracle AUC depends on prevalence and covariate structure, not on positive count
alone. V5 keeps the instrument and changes only the cohort it is pointed at. The gate can fail, and
nothing downstream is conditioned on it passing.

## Inputs
- `STAGE_23_2_ROLE_A_CONFIRMATION_V5.md` (live); V1–V4 archived unedited in `arcive/`
- `D:\GSE227151_Rewind\NEW DATA\` — 15 GEX files, 2 Figshare metadata workbooks, 3 author outcome
  objects (see that folder's `MANIFEST.md` for sha256 of every file)
- `R2/` and `R3/` `filtered10XCells.txt`, the shared `stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt`
- author code `plotScripts/rewind10X/R2`, `R3`

## Files added
- `experiments/run_stage23_2h_confirmation.py`
- `tests/test_stage23_2h_confirmation.py`
- `plans/(newer)practical plans/STAGE_23_2_ROLE_A_CONFIRMATION_V5.md`
- `results/stage23_2h/` — benchmark, representation, bdepth, power, confirmation, verdict

## Files modified
- `tests/test_stage23_2_role_a_resolution.py` — two path resolvers so V4's archival does not break
  contracts that pin it
- `experiments/run_stage23_2_role_a_resolution.py` — same, for 23.2G's V4 digest lookup

## What did NOT change
- No Stage-22 file. R1's outcome rule, clone set and folds are untouched.
- No Stage-23 gate re-run or re-interpreted. The historical Role-A FAIL remains permanent.
- `stage_23_2G_RECORD.md` untouched, including its now-superseded conclusion.

---

## 23.2H-A — benchmark construction

### The reconstructions reproduce the authors' own objects exactly

This is the load-bearing check. A rule we invented would be worthless; a rule that regenerates the
authors' serialized output bit-for-bit is theirs.

```text
  R2  primedCellsInd.rds
        lineages          26 author  ==  26 reconstructed    sets identical
        cells             79 author  ==  79 reconstructed    sets identical
        SampleNum         {S4, S5} only   -> replicate 2, confirming the map a third time
        spike-in coef     LSD1_4A rel.err 3.5e-16 vs the author's stored nUMINorm/nUMI
                          LSD1_4B rel.err 1.5e-16

  R3  primedCellIDList.rds elements 6, 15, 24  (top-200 per unit)
        FS_1   27 author  ==  27 reconstructed    identical
        FS_2   75 author  ==  75 reconstructed    identical
        FS_3   83 author  ==  83 reconstructed    identical
```

R2's rule was recovered, not found: V4 recorded it as unavailable from shipped code, and no producer
script exists. It is `min(nUMINorm_A, nUMINorm_B)` over the inner join, restricted to lineages with
a linked 10X cell, top 26 — ranks 1–26 all in the author set, 27+ none, with a clean 60.22 → 33.34
boundary gap.

### The author spike-in bug is real, material, and confined to units 2 and 3

```text
  corrected vs author-coefficient selection, per unit, of 200 lineages each
    FS_1   200 / 200 identical      the one unit the author code scales correctly
    FS_2   186 / 200
    FS_3   170 / 200
```

That is the signature the bug predicts exactly: `lmr[[i]]`/`lmr[[i+1]]` with the 1..3 loop counter
is correct for the first pair and wrong for the other two.

### The frozen cohort

```text
  ELIGIBLE (ungated, linked, post-ambiguity-exclusion)      cells   clones   prevalence
    replicate 2   S4+S5   GSM7092517/18                      3480     1827      1.42%
    replicate 3   S1      GSM7092519                           598      483     10.14%
                                                    total    4078     2310      3.25%

  POSITIVES, PRIMARY arm (corrected coefficients)          cells   clones
    replicate 2   R2_MIN_PAIRED_TOP26_V1                       79       26
    replicate 3   R3_MAX_PAIRED_TOP200_UNION_V1                60       49
                                                    total     139       75

  POSITIVES, SENSITIVITY arm (author coefficients)         cells   clones
    replicate 3   R3_MAX_PAIRED_TOP200_UNION_AUTHORBUG         61       50
                                                    total     140       76
```

Folds: 5, clone level, seed 23511, stratified on `(biological_replicate, y_primed)`. Every fold
carries both replicates. Sorted samples S2/S3 never enter — asserted by contract on the cells table.

### A guard fired, and it caught a real error

V5 §8's frozen expectation table originally carried `rep 3 = 61 cells / 50 clones` for both arms,
because the numbers had been transcribed from `primedCellIDList.rds` — which is the **sensitivity**
arm's output, not the primary's. 23.2H-A halted on the first run rather than proceeding with a
mismatch. The primary arm's true count is 60/49. V5 §8 now states both arms separately and records
the correction inline. No statistic had been computed at that point.

---

## 23.2H-B — representation

```text
  clone pseudobulk      2310 clones x 36601 genes
  normalisation         sum raw counts per clone -> CP10K -> log1p, applied exactly once
  Bdepth                4 columns, outcome-free, recomputed from the raw matrices
  cross-check           detected features == the normalised matrix's non-zero pattern   PASS
  all clones            total raw UMI > 0                                               PASS
  runtime               1.3 min
```

---

## 23.2H-C — the V5 §9 power gate

Recorded before any confirmatory statistic existed, as §10 requires.

```text
  cohort                 2310 clones, 75 positives, prevalence 3.25%
  positives per fold     16 / 14 / 15 / 15 / 15
  beta calibration       0.62066  ->  median oracle AUC 0.6596   (target 0.66)
  null      n=200        p95 = 0.007374   mean = 0.000183   median oracle AUC 0.5006
  alternative n=100      mean delta_AP = 0.012811          median oracle AUC 0.6585

  MEASURED POWER         0.64
  threshold              0.80
  GATE 18.3              FAILS          verdict DESIGN_UNDERPOWERED
```

### What the corrected measurement actually bought

The V4 rule and the V5 rule reach the same verdict here, and it is worth being exact about why that
does **not** make V5 pointless.

```text
  V4 §18.3   >= 140 positive clones            75 available    FAIL, on an imported number
  V5 §9      measured power >= 0.80            0.64 measured   FAIL, on this cohort's own number
```

V4's own ladder would have implied roughly 0.52 at this event count, because it was built on R1's
1.11% prevalence. The realized figure is **0.64** — materially better, and in the direction V5
predicted, because this cohort carries ~3× R1's prevalence. So the category error V5 identified was
real and the correction moved the number. It simply did not move it far enough.

That is the outcome a pre-registered gate is supposed to be able to produce. The threshold was
frozen in the protocol and in code before the measurement, the failure branch was written out in
advance (§9.3), and nothing downstream was conditioned on the gate passing.

### Consequence, per V5 §9.3 and V4 §17

The confirmatory analysis is still run and still reported in full. A positive result is recorded as
**underpowered supporting evidence**; a null result is **not** evidence against the hypothesis.
`ROLE_A_CONFIRMATORY_SUPPORTED` cannot be emitted.

Runtime 2.9 h across two parallel shards (200 null + 100 alternative), plus 0.6 min to merge.

## 23.2H-D — the confirmatory analysis

Run in full despite the failed power gate, as V5 §9.3 requires. Both V5 §7 arms, 200 permutations
each, 115.3 min per arm.

### Primary arm — corrected spike-in coefficients

```text
POOLED
  AP  R1  Bdepth only          0.069013
  AP  R3  state + Bdepth       0.099517
  delta_AP                    +0.030504

  null p95   0.028101      null mean 0.007744      5 of 200 null draws >= observed
  p_perm     0.0299        exceeds_null_p95 TRUE

  GATE 18.4  pooled dual gate            PASSES

PER REPLICATE                                    (V4 §16.2, direction-gating only)
  rep 2   1827 clones   26 pos   delta_AP +0.058251   CI95 [+0.003038, +0.163776]
  rep 3    483 clones   49 pos   delta_AP +0.023291   CI95 [-0.035922, +0.088588]

  GATE 18.5  positive in every replicate  PASSES

CROSS-REPLICATE TRANSFER                         (V4 §16.3, reported not gating)
  2 -> 3   delta_AP +0.019420
  3 -> 2   delta_AP +0.028963
```

Both transfer directions are positive, which is the reassuring case: V4 §16.3 flags a *negative*
transfer in both directions alongside a passing pooled test as the signature of a within-replicate
structure artifact. That signature is absent.

### Sensitivity arm — author spike-in coefficients

```text
  delta_AP   +0.021686      null p95 0.009472      0 of 200 null draws >= observed
  p_perm      0.0050        exceeds_null_p95 TRUE

  rep 2   26 pos   delta_AP +0.004289   CI95 [-0.069588, +0.062345]
  rep 3   50 pos   delta_AP +0.056267   CI95 [+0.004417, +0.127957]

  GATE 18.4  PASSES        GATE 18.5  PASSES
```

The V5 §7 bug decision does not drive the verdict. Both arms pass both statistical gates, agree on
direction in every replicate, and agree that the pooled effect is positive.

### An instability the two arms expose

```text
                     primary          sensitivity
  rep 2 delta_AP     +0.058251        +0.004289
  rep 3 delta_AP     +0.023291        +0.056267
  dominant replicate     2                 3
```

The arms differ by **one positive clone** — yet they disagree about which replicate carries the
larger effect, and about which replicate's bootstrap CI excludes zero. This is not a finding about
either replicate; it is what a 0.64-power design does to a per-replicate decomposition. **No claim
may be made about where the signal lives.** The pooled direction and the pooled gate are the only
stable quantities here, and even those must carry the underpowered caveat.

---

## 23.2H-E — verdict

```text
  18.1  >= 2 independent non-R1 biological replicates      PASS
  18.2  source-faithful reconstruction                     PASS
  18.3  measured design power >= 0.80                      FAIL   (0.64)
  18.4  pooled dual gate                                   PASS
  18.5  positive in every replicate                        PASS
  18.6  no un-re-gated material benchmark change           PASS

  EXIT      ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
  STAGE 24  BLOCKED
```

Five of six gates pass. The one that fails is the one that cannot be argued around: the design did
not have the power to make a confirmatory claim, and V5 §9.3 fixed that consequence in advance.

---

## Bugs found

- **The published Rewind R3 spike-in normalisation is mis-indexed.** Detailed above and in V5 §7.
  Verified against the authors' own serialized intermediate, not inferred from the code alone.
- **The frozen R1 loader path is broken.** `D:\GSE227151_Rewind\` no longer holds
  `filtered10XCells.txt`, `stepThreeStarcodeShavedReads_BC_10X.txt` or `..._BC_gDNA.txt` at its
  root; they are in `r1\`. `diag_stage21d_public_reconstruction.py::rewind_required` reads them from
  the root. Not repaired here — it touches frozen benchmark inputs and 23.2H does not use them.
- **Eight V4-pinned contracts would have broken on archival.** `test_..._pins_the_live_confirmation_protocol`
  and the module's 23.2G digest lookup both hard-coded V4's live path. Both now resolve through the
  same live-or-archived helper the rest of the module already used.

## Tests
- 32 Stage-23.2H contracts, 0 skipped · 110 Stage-23.2 contracts still pass after V4's archival
- ruff clean on both modules and both test files

## Scientific interpretation

**Proves:** the corrected Role-A hypothesis — that pretreatment transcriptional state predicts the
Rewind priming outcome beyond a depth-complete nuisance baseline — is **supported on two
independent biological replicates that were never used to design the correction**. The pooled
`ΔAP` of `+0.0305` exceeds its permutation null's 95th percentile with `p_perm = 0.0299`, the
direction is positive in both replicates, both cross-replicate transfer directions are positive,
and the result survives the author-coefficient sensitivity arm at `p_perm = 0.0050`.

That is a real change in the evidential position. Stage 23's Role-A failure (`ΔAP +0.0105`,
`p_perm 0.0846`) has now been followed by a positive, pre-registered result on untouched evidence
under the mechanically indicated correction.

**Does NOT prove, and must not be reported as:**

- **A confirmed Role-A claim.** Measured power is 0.64 against a pre-registered 0.80 threshold. V5
  §9.3 and V4 §17 both fixed this consequence before the number existed: a positive result at this
  power is *underpowered supporting evidence*, not confirmation. `ROLE_A_CONFIRMATORY_SUPPORTED`
  cannot be emitted and Stage 24 does not open.
- **An effect size.** An underpowered design that clears significance is selected on significance,
  so `+0.0305` is more likely to overstate the truth than to understate it. It is a direction, not
  a magnitude.
- **Anything about which replicate carries the signal.** A one-clone difference between arms flips
  the dominant replicate and flips which bootstrap CI excludes zero.
- **A uniform endpoint across replicates.** Three replicates, three source-defined
  operationalisations (V5 §6.2). This is a property of the source study, not a choice we made, and
  it stands wherever this result is reported.

**What would resolve it.** The gate failed on event count under this geometry, not on mechanism.
Reaching 0.80 needs more positive clones from replicates that qualify under §11 — which the
existing GSE227151 material cannot supply, since replicates 2 and 3 are now consumed and the sorted
samples are a different population. That is a new-evidence requirement, which is exactly what the
`ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE` exit means.

**Stage-27 cost.** Replicates 2 and 3 are consumed by this confirmation and are no longer available
as an untouched replication set (V5 §12).

## Next action
Stage 24 remains **BLOCKED**. The roadmap exit is `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE`. The
decision now facing the project is a scope decision, not a mechanical one, and it is the user's:
either acquire a qualifying Role-A replicate large enough to clear 0.80, or revise the roadmap's
claim architecture — the vision's `ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL` branch explicitly
contemplates deciding whether a different evidence path may open Stage 24, and Role B remains a
valid but non-substitutable positive.

---

# ADDENDUM — extending the permutation null to n = 2000

## Why

Gate 18.4's first criterion is `observed > null p95`. At n = 200 that p95 rests on roughly ten
order statistics, and the two arms — which differ by **one clone** — returned p95 values of
`0.028101` and `0.009472`, a 3× spread. That is far more likely to be Monte Carlo error in the
*estimate* than instability in the *design*, and it is the single largest unforced weakness in the
result as it stands. It is fixable with compute alone: no new data, no threshold moved, no
parameter changed.

More permutations do not change the estimand. They estimate the same p-value and the same p95 more
precisely, and the outcome can just as easily weaken the result as strengthen it.

## Pre-commitment, made BEFORE any draw beyond 200 was inspected

```text
n_perm                2000, both arms, fixed now
supersession          the n=2000 result REPLACES the n=200 result WHICHEVER SIDE IT LANDS ON.
                      If p_perm rises above 0.05, that is the reported result and gate 18.4 FAILS.
no optional stopping  no intermediate value is inspected. n is not revised after looking.
                      The cache makes resuming cheap, which is exactly why this clause exists.
preservation          the n=200 results are archived at
                        stage23_2h_confirmation_n200.json
                        stage23_2h_confirmation_authorbug_n200.json
                      and pinned by contract, so the superseded numbers cannot quietly vanish.
unchanged             every other parameter: cohort, folds, strata, seeds, statistic, thresholds.
```

The numbers on the table at the moment this commitment was made:

```text
  PRIMARY      delta_AP +0.030504   null p95 0.028101   5/200 exceed   p_perm 0.0299
  SENSITIVITY  delta_AP +0.021686   null p95 0.009472   0/200 exceed   p_perm 0.0050
```

Note that the power gate (18.3) is **not** affected by this and remains FAILED at 0.64. Extending
the null cannot open Stage 24; it can only establish whether the marginal pooled result is real.

## Engine

Draw `b` is seeded from `default_rng(SEED_PERMUTATION + b)` and nothing else, so its value is
independent of how many draws preceded it. Shards may be split arbitrarily, run in any order,
interrupted and resumed. The cache stores the raw null AP rather than the delta, so a shard needs
no observed statistic and can run before one exists.

Verified by `--stage 23.2h-smoke` on 18 real draws against a real sequential ground truth:

```text
  sequential run produced every draw                          PASS  6/6
  3 shards run OUT OF ORDER are bit-identical to sequential   PASS  max abs diff 0.000e+00
  interrupt-and-resume reproduces sequential exactly          PASS  6/6
  mixed-protocol cache is refused                             PASS  guard fired
  the n=200 result is preserved on disk                       PASS
```

Measured 26.2 s/draw single-process unloaded, 33.8 s/draw observed under load.

## Amendment to the pre-commitment — the sensitivity arm stays at n = 200

Decided by the user **before the primary arm's n=2000 result existed**, so the choice cannot be
outcome-driven. Recorded here rather than applied silently, because it departs from the
"both arms" clause written above.

```text
  PRIMARY      n = 2000    proceeds as committed
  SENSITIVITY  n =  200    unchanged

  reason       cost. The sensitivity arm is not marginal -- 0 of 200 null draws reach the
               observed statistic, p_perm 0.0050, with a margin of +0.012214 against a null p95
               of 0.009472. Extra precision there cannot change any gate or any decision, whereas
               the primary arm's 5-of-200 margin genuinely can.
  cost avoided ~9.7 h of compute
  status when  the primary arm stood at 1546 / 2000 draws with no intermediate value inspected
  decided
```

The asymmetry is recorded as a limitation: the two arms are reported at different permutation
counts, and the sensitivity arm's `p_perm = 0.0050` carries the coarser 1/201 resolution.

The primary arm was allowed to finish because it is the marginal one and was 2.2 h from completion
with 1546 draws already computed; stopping it would have discarded seven hours of finished work and
left the fragile number in place.

---

# AUDIT FINDING — the power study rejects against a different null than the test

Found during a post-hoc audit of the whole stage, **after** gate 18.3 had already been recorded as
failed. It makes that failure more severe, not less, and it weakens the case for opening Stage 24.

## The two nulls are not the same distribution

```text
  POWER STUDY null   delta = AP(R3 | y_synth) - AP(R1 | y_synth)
                     both models refit on a synthetic y drawn at beta = 0,
                     positives placed uniformly at random within each outer fold
                     mean +0.000183   sd 0.005522   p95 0.007374

  TEST null          delta = AP(R3 | PERMUTED profiles, real y) - AP(R1 | real y, OBSERVED)
                     R1 pinned at its observed value -- it is expression-free, so the permutation
                     provably cannot move it. Profiles permuted within
                     (size|n_lanes|replicate) x train/test side.
                     mean +0.007744   sd 0.010241   p95 0.028101

  rejection thresholds differ by 3.81x
```

The test null is centred positive and wider because R1 is pinned while R3 is re-selected over
`(K, C)` on every permuted draw, so the null absorbs the model-selection advantage. That makes the
**test** conservative, which is correct and desirable. It also means the power study's threshold is
not the threshold the test actually applies.

## What this does to the recorded power

```text
  power against the POWER-STUDY null    0.640   <- the number recorded for gate 18.3
  power against the TEST's own threshold 0.080  <- audit recomputation
```

## What can and cannot be concluded

**Can:** `0.64` overstates the power of the test as actually run. The direction is unambiguous —
the real rejection threshold is 3.81x higher than the one the power study used.

**Cannot:** that `0.080` is the correct figure. That recomputation mixes constructions too: the
alternative draws refit R1 on the synthetic y, whereas the test null pins R1 at its observed value.
A correct estimate requires nesting a full permutation null inside every alternative draw —
100 x 200 = 20,000 nested fits, roughly four days of compute. Not feasible, and not attempted.

## The estimator is NOT the problem

Checked directly, because if the two stages had used different pipelines the power number would
have described a different estimator entirely:

```text
  _observed()             delta_AP +0.0305040385     used by 23.2H-D, the test
  S232._delta_ap_once()   delta_AP +0.0305040385     used by 23.2H-C, the power study
  identical: True   |diff| 0.000e+00
```

The discrepancy is purely in the null construction.

## Provenance and ownership

The mismatch is inherited from 23.2E / V2 §9, which built its ladder against the same
label-randomisation null. It is not something 23.2H introduced. But V5 §9.2 asserts that *"the
measurement instrument is unchanged. Only the cohort it is applied to changes"* — and that assertion
was made without auditing the instrument. The cohort was re-derived carefully; the null the
instrument uses was taken on trust. That is the error, and it is this stage's.

## Consequence for the decision

Gate 18.3 fails either way, so no verdict changes. What changes is the strength of the argument for
amending the opening rule: if the design's real power against its own test is nearer 0.08 than 0.64,
then the positive predictive value of the observed significant result falls below 0.5 at any
plausible prior, and "underpowered supporting evidence" becomes a thin basis for opening even a
development stage.

**This finding must accompany the Role-A result wherever it is reported, including in the paper's
limitations.**

---

# THE n = 2000 PRIMARY RESULT

```text
  delta_AP            +0.030504      (unchanged -- this is the observed statistic)
  null p95             0.028844      mean +0.007952   sd 0.011318
  margin              +0.001660      5.4% of the statistic
  null >= observed     79 / 2000
  p_perm               0.03998

  GATE 18.4  observed > null p95 AND p_perm <= 0.05     PASSES
  GATE 18.5  positive in every replicate                PASSES
    rep 2   delta_AP +0.058251   CI95 [+0.003038, +0.163776]
    rep 3   delta_AP +0.023291   CI95 [-0.035922, +0.088588]
```

The extension was worth doing: `p_perm` moved from `0.0299` to `0.0400`. Still under 0.05, so the
gate holds, but appreciably closer to the boundary than the 200-draw estimate suggested. The
pre-commitment is honoured — this replaces the n=200 result, and it happens to have moved in the
less favourable direction.

## A correction to the audit finding above

The audit hypothesised that the 3.81x gap between the two arms' null p95 values (`0.028101` primary
vs `0.009472` sensitivity) was *"very likely Monte Carlo noise in the p95 estimate, since 200
permutations means the p95 rests on about ten order statistics."*

**That hypothesis is now disproven by the primary arm's own data.** Going from 200 to 2000 draws
moved the primary's null p95 by less than 3%:

```text
  primary   null p95 at n=200    0.028101
  primary   null p95 at n=2000   0.028844      +2.6%
```

The estimate was already stable at 200. So the between-arm gap is **not** estimation noise. The two
arms differ by one positive clone, and because folds are stratified on `(replicate, y)`, that one
clone reshuffles the fold assignment — and the null distribution turns out to be genuinely sensitive
to which clones land in which fold.

That is a **worse** finding than the one it replaces. Estimation noise would have been fixable with
compute. Fold sensitivity of the null at this event count is a property of the design, and it sits
alongside the flipping dominant replicate as evidence that the confirmation cohort is too small to
support fine-grained conclusions. It does not invalidate the pooled test — that test is valid
whatever the null's shape, because the null is constructed under the exchangeability the test
assumes — but it does mean the *margin* of the pooled result should not be read as precise.

The claim that extra permutations would settle the between-arm discrepancy was wrong, and it was
mine. The permutations settled a different and still useful question: whether `p_perm = 0.0299` was
a 200-draw artifact. It was, partly — the better estimate is `0.0400`.

---

# THE POWER QUESTION, ANSWERED PROPERLY

The audit above established that the recorded power measured against the wrong null but could not
say by how much, offering only a range of `0.08` to `0.64`. That range was too wide to decide on, so
the measurement was redone with the null built the way the **test** builds it.

```text
  method   2 synthetic y's at oracle AUC 0.66, real per-fold event counts
           150 profile-permutation draws each, R1 pinned at that y's own observed value
           resulting p95 applied to the 100 alternative draws already cached
           300 draws, ~70 min across 3 shards

                                 null p95      corrected power
    synthetic y0                 0.010834             0.55
    synthetic y1                 0.014090             0.40
    mean threshold               0.012462             0.45

    threshold spread                            1.30x
    recorded instrument (label-randomisation)   p95 0.007374  ->  0.640
```

**Corrected power at the pre-registered alternative: 0.45.** Still fails the 0.80 threshold. Gate
18.3's FAILED verdict, recorded at 0.64 under the frozen instrument, stands unchanged; 0.45 is
recorded beside it and nothing is re-scored.

## Two of my own earlier claims are corrected by this

**The `0.08` figure was too pessimistic.** It borrowed the *real* data's null threshold
(`0.028844`), but under a synthetic alternative the correctly-constructed threshold is only
`0.011`-`0.014`. The truth is much nearer the optimistic end of the range I gave.

**The null is not badly fold-sensitive.** Two synthetic cohorts with matched per-fold event counts
give thresholds `1.30x` apart, not the `3.9x` seen between the primary and sensitivity arms. So the
earlier conclusion — that the between-arm gap proved a design-level fold sensitivity — overreached.
Something the synthetic draws do not reproduce is responsible, most plausibly the item below.

## A further flaw, found in the same pass and NOT corrected

`_assign()` places synthetic positives uniformly at random within each outer fold, **ignoring
replicate**. The real cohort is severely imbalanced -- 10.14% prevalence in replicate 3 against
1.42% in replicate 2 -- so the synthetic cohorts do not reproduce the real one's structure. This is
inherited from 23.2E, where a single replicate made it moot, and it is very likely why the real
data's null (`p95 0.028844`) is so much wider than the synthetic alternatives' (`0.011`-`0.014`).

`0.45` is therefore a better-constructed approximation, not an exact figure. It is the number to
report. It is not a number to defend to three decimal places, and correcting `_assign()` to
stratify on replicate is the obvious next improvement if the power question is ever reopened.

## Decision impact

```text
   prior     power 0.08     power 0.45     power 0.64
    0.10           0.15           0.50           0.59
    0.25           0.35           0.75           0.81
    0.50           0.62           0.90           0.93
```

Positive predictive value of the observed result at alpha = 0.05. The decision-relevant band narrows
from `[0.35, 0.81]` to `[0.75, 0.81]` at a sceptical 0.25 prior.

## An operational defect worth recording

One of the 300 audit draws (index 10) was lost to a race: three shards appending to a single cache
file. No malformed lines -- the write simply vanished. It was detected by the completeness check,
recomputed, and verified. The 2000-draw permutation cache was independently verified complete.

The design should give each shard its own file rather than relying on append atomicity across
processes. Not changed here, because doing so would invalidate the existing caches; recorded as a
known defect for any future sharded run.

---

## Addendum 2026-08-28 — the run-order check could never survive a clone

`test_the_power_gate_was_recorded_before_the_confirmatory_statistics` asserted V5 §10's ordering
rule from **file modification times**:

```python
assert POWER.stat().st_mtime <= CONFIRM.stat().st_mtime
```

**Git does not preserve mtime.** On any fresh checkout both files carry the checkout instant in
arbitrary order, so the assertion was a coin flip on the filesystem it happened to run on. It
passed on the long-lived working tree that wrote the files and failed on every CI run — it was the
reason CI had been red continuously, across commits that had nothing to do with Stage 23.2.

**Run order is not recoverable from the committed artifacts.** Neither JSON carries a timestamp,
and only the confirmation records a git commit. Rather than assert something the repository cannot
support, the contract now checks what the ordering rule was there to protect:

```text
  the power gate is recorded, verdict DESIGN_UNDERPOWERED
  gate_18_3_measured_power is False, power 0.64 < threshold 0.80
      -> a gate that FAILED cannot have been tuned to pass
  the confirmatory run is pinned to the frozen V5 protocol, so it did not run
      against a revision
```

The run order itself remains documented in this record — 23.2H-C above 23.2H-D — and that is a
**record of the claim, not proof of it.** Nothing in the Stage-23.2 result changes; only an
untestable assertion was replaced with testable ones, and the limitation is now stated rather than
papered over.
