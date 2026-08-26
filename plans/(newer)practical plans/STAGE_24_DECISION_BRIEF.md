# STAGE 24 — DECISION BRIEF

**Status:** FINAL. All results in.
**Decision owner:** the user. This brief does not decide anything; it lays out what is known, what
each option costs, and what I would do.

---

# 1. The question in one paragraph

Stage 23 ended with Role A failing its permutation gate, and the roadmap blocked Stage 24 until a
pre-registered confirmatory analysis succeeded on evidence not used to design the correction.
Stage 23.2H ran that confirmation on biological replicates 2 and 3 under protocol V5. **It found a
positive, replicated signal, and it failed one gate: measured design power.** The decision is
whether a development stage may open on positive-but-underpowered evidence, or whether Role A is
set aside and the project's claim architecture is revised.

---

# 2. The variables

Every quantity the decision turns on, defined once.

## 2.1 The statistic

```text
AP                  average precision of out-of-fold predictions, pooled over the 5 outer folds
R1  (baseline)      Bdepth-only model:  [log1p(n_pretreatment_cells), n_lanes,
                    log1p(total_raw_GE_UMI), log1p(n_detected_GE_features), is_replicate_3]
R3  (state model)   the same nuisance block PLUS K expression PCs
delta_AP            AP(R3) - AP(R1).  The quantity of interest: what transcriptional STATE adds
                    beyond technical depth and replicate identity.
```

`delta_AP` is the same statistic Stage 23 used. Nothing about it changed.

## 2.2 The null

```text
permutation         whole clone expression profiles are permuted WITHIN
                    (size{1,2,3+} x n_lanes x biological_replicate) x (train|test side)
what it preserves   depth structure, replicate structure, and the train/test boundary
p_perm              (1 + #{null >= observed}) / (n_perm + 1)
n_perm              PRIMARY 2000   SENSITIVITY 200   (see §3.5 for why they differ)
```

Permuting inside replicate is deliberate: it puts any replicate-identifying batch signal into the
null as well as the observed statistic, so the test cannot be passed by that signal.

## 2.3 The gates

```text
18.1   >= 2 independent biological replicates, none of them R1
18.2   each replicate's outcome reconstructed from ITS OWN source rule, verified by
       reproducing the authors' serialized objects exactly
18.3   measured design power >= 0.80 at the pre-registered alternative (oracle AUC 0.66)
18.4   observed delta_AP > null p95   AND   p_perm <= 0.05
18.5   delta_AP positive in EVERY qualifying replicate analysed separately
18.6   no material benchmark semantic changed without re-gating
```

## 2.4 Power, and what it is not

```text
measured power      P(the pipeline detects a signal of oracle AUC 0.66 | that signal is real)
                    estimated by simulation on THIS cohort's real X, Bdepth, folds and event counts
the alternative     0.66, chosen in advance by 23.2E to match R1's observed R3 ROC-AUC of 0.6628
threshold           0.80, frozen before the measurement
```

Power is a property of the **design**, not of the result. It does not tell you whether the observed
effect is real; it tells you how often this design would have found one of that size.

---

# 3. The results

## 3.1 The cohort

```text
                                              cells   clones   positives   prevalence
  replicate 2   S4+S5   GSM7092517/18          3480     1827          26        1.42%
  replicate 3   S1      GSM7092519               598      483          49       10.14%
                                      total    4078     2310          75        3.25%
```

Ungated samples only. The sorted samples `GSM7092520/21` are excluded on population grounds — they
select on proliferation speed and are therefore a different pre-state population from the Stage-22
R1 benchmark.

## 3.2 The reconstructions are the authors', not ours

This is what makes 18.2 meaningful.

```text
  R2  primedCellsInd.rds        26/26 lineages, 79/79 cells, sets IDENTICAL
                                spike-in coefficients agree to 1e-16 with the author's stored values
  R3  primedCellIDList.rds      27/75/83 cells per unit, all three sets IDENTICAL
```

R2's rule had been recorded as unrecoverable (no producer script ships). It was recovered:
`min(normA, normB)` over the inner join, restricted to lineages with a linked 10X cell, top 26.

We also found a real bug in the authors' published R3 normalisation — a spike-in coefficient
mis-index that leaves unit 1 correct and corrupts units 2 and 3. Both arms are run: corrected
coefficients primary, author coefficients as a declared sensitivity.

## 3.3 Power — gate 18.3

```text
  cohort                  2310 clones, 75 positives, prevalence 3.25%
  beta calibration        0.6207  ->  median oracle AUC 0.6596  (target 0.66)
  null       n=200        p95 0.007374
  alternative n=100       mean delta_AP 0.012811

  MEASURED POWER at AUC 0.66      0.64        threshold 0.80        GATE 18.3 FAILS
```

### The power curve — reported, never gating

```text
  oracle AUC 0.66   power 0.640     <- the pre-registered GATE point
  oracle AUC 0.70   power 0.900
  oracle AUC 0.75   power 0.970
```

The curve is steep exactly where it matters: 0.66 -> 0.70 takes power from 0.64 to 0.90. The design
sits just below the threshold-crossing region, adequately powered for effects of roughly AUC >= 0.68
and upward.

**But every number in that curve is measured against the wrong null — see §3.8, and §3.10 for how far off.**

The curve exists to convert *"underpowered"* into *"adequately powered for effects of AUC >= X."*
It may not be used to evaluate the gate at any point other than 0.66 — Addendum 1 forbids that
explicitly, along with reading power off the curve at the observed effect size (the observed-power
fallacy).

### Why V4's rule and V5's rule both fail, and why V5 was still worth writing

```text
  V4 §18.3   >= 140 positive clones     75 available    FAIL, on a number imported from R1
  V5 §9      measured power >= 0.80     0.64 measured   FAIL, on this cohort's own number
```

V4's ladder would have implied roughly 0.52 at this event count, because it was built on R1's 1.11%
prevalence. The realized figure is 0.64 — better, in the direction V5 predicted, because this cohort
carries ~3x that prevalence. The category error was real. It did not move the number far enough.

## 3.4 The confirmatory result — gates 18.4 and 18.5

### Primary arm, corrected coefficients, n_perm = 2000  -- THE REPORTED RESULT

```text
  AP  R1 baseline        0.069013
  AP  R3 state+baseline  0.099517
  delta_AP              +0.030504
  null p95               0.028844      margin +0.001660   (5.4% of the statistic)
  null mean             +0.007952      sd 0.011318
  79 of 2000 null draws >= observed    p_perm 0.03998

  GATE 18.4   observed > null p95 AND p_perm <= 0.05      PASSES
```

Extending from 200 to 2000 draws moved `p_perm` from `0.0299` to `0.0400`. The gate holds, but the
result is closer to the boundary than the 200-draw estimate implied. This honours the pre-commitment
recorded before the extension: the n=2000 figure replaces the n=200 figure, and it happened to move
against the result.

### Primary arm at n_perm = 200 — the superseded predecessor, preserved

```text
  AP  R1 baseline        0.069013
  AP  R3 state+baseline  0.099517
  delta_AP              +0.030504
  null p95               0.028101      margin +0.002403   (8% of the statistic)
  5 of 200 null draws >= observed      p_perm 0.0299
```

### Per replicate — direction-gating only (V4 §16.2)

```text
  replicate 2   1827 clones   26 pos   delta_AP +0.058251   CI95 [+0.003038, +0.163776]
  replicate 3    483 clones   49 pos   delta_AP +0.023291   CI95 [-0.035922, +0.088588]
```

### Cross-replicate transfer — reported, not gating (V4 §16.3)

```text
  2 -> 3   delta_AP +0.019420
  3 -> 2   delta_AP +0.028963
```

Both positive. The failure signature V4 §16.3 warns about — negative transfer in both directions
alongside a passing pooled test, indicating a within-replicate artifact — is **absent**.

### Sensitivity arm, author coefficients, n_perm = 200

```text
  delta_AP              +0.021686
  null p95               0.009472      margin +0.012214
  0 of 200 null draws >= observed      p_perm 0.0050

  replicate 2   26 pos   delta_AP +0.004289   CI95 [-0.069588, +0.062345]
  replicate 3   50 pos   delta_AP +0.056267   CI95 [+0.004417, +0.127957]
```

The bug decision does not drive the verdict: both arms pass both statistical gates and agree on
direction in every replicate.

## 3.5 Why the two arms have different permutation counts

The primary arm was extended to 2000 because its 5-of-200 margin genuinely could flip. The
sensitivity arm stayed at 200 because it is not marginal — 0 of 200 draws reach the observed value.
The decision was made **before** the primary arm's 2000-draw result existed, so it cannot be
outcome-driven, and it is recorded in `stage_23_2H_RECORD.md` with that timestamp.

## 3.6 An instability worth knowing about

```text
                     primary          sensitivity
  rep 2 delta_AP     +0.058251        +0.004289
  rep 3 delta_AP     +0.023291        +0.056267
  dominant replicate     2                 3
```

The arms differ by **one positive clone**, and they disagree about which replicate carries the
larger effect and which bootstrap CI excludes zero. This is what 0.64 power does to a per-replicate
decomposition. **No claim may be made about where the signal lives.**

## 3.7 Gate summary

```text
  18.1  two independent non-R1 replicates          PASS
  18.2  source-faithful reconstruction             PASS
  18.3  measured design power >= 0.80              FAIL   0.64
  18.4  pooled dual gate                           PASS   p_perm 0.03998 at n=2000
  18.5  positive in every replicate                PASS
  18.6  benchmark compatibility                    PASS
```

## 3.8 AUDIT FINDING — the power study rejects against a different null than the test

Found after gate 18.3 had already been recorded as failed. It makes the failure **more** severe and
weakens the case for opening Stage 24.

```text
  POWER STUDY null   delta = AP(R3 | y_synth) - AP(R1 | y_synth)     both models refit on
                     a synthetic y drawn at beta = 0
                     mean +0.000183   sd 0.005522   p95 0.007374

  TEST null          delta = AP(R3 | PERMUTED profiles, real y) - AP(R1 observed)
                     R1 PINNED -- it is expression-free, so the permutation cannot move it
                     mean +0.007952   sd 0.011318   p95 0.028844

                                             rejection thresholds differ by 3.91x
```

The test null is centred positive and wider because R1 is pinned while R3 is re-selected over
`(K, C)` on every permuted draw, so the null absorbs the selection advantage. That makes the **test**
conservative, which is correct. It also means the power study's threshold is not the test's.

### The whole curve, recomputed against the test's own threshold

```text
   oracle AUC     power vs power-null     power vs TEST null     status
        0.66                    0.640                  0.080     GATE POINT
        0.70                    0.900                  0.210     reported only
        0.75                    0.970                  0.560     reported only
```

**Neither column was correct**, and the audit was therefore carried out properly rather than left
as a range. See §3.10.

**Not the estimator's fault.** Checked directly: `_observed()` and the power study's
`_delta_ap_once()` return `+0.0305040385` on identical inputs, `|diff| 0.000e+00`. The discrepancy
is purely the null construction.

**Ownership.** The mismatch is inherited from 23.2E / V2 §9. But V5 §9.2 asserts *"the measurement
instrument is unchanged; only the cohort changes"* — an assertion made without auditing the
instrument. The cohort was re-derived carefully; the null was taken on trust. That error is this
stage's.

## 3.9 The null p95 was NOT estimation noise -- but see 3.10

An earlier hypothesis — that the 3.9x gap between the two arms' null p95 values was Monte Carlo
error at 200 draws — is disproven by the primary arm's own extension:

```text
  primary null p95 at n=200     0.028101
  primary null p95 at n=2000    0.028844      +2.6%
```

The estimate was already stable at 200. The between-arm gap is therefore genuine sensitivity of the
null to fold assignment: the arms differ by one positive clone, folds are stratified on
`(replicate, y)`, and that one clone reshuffles them.

This looked like a property of the design rather than of the estimate. **§3.10 partly walks that
back too**: under matched conditions — two synthetic cohorts with the same per-fold event counts —
the null thresholds differ by only 1.30x, not 3.9x. So the null is *not* wildly fold-sensitive in
general. The large primary-vs-sensitivity gap comes from something the two synthetic draws do not
reproduce, most plausibly the real cohort's severe replicate-prevalence imbalance (§3.10).

Either way it does not invalidate the pooled test — that test is valid whatever the null's shape —
but the *margin* should still not be read as precise.

## 3.10 THE POWER QUESTION, ANSWERED PROPERLY

§3.8 established that the recorded power measured against the wrong null but could not say by how
much. That gap was closed by building the null the way the **test** builds it — profile permutation
with R1 pinned — on synthetic cohorts drawn under the pre-registered alternative.

```text
  method   2 synthetic y's at oracle AUC 0.66, real per-fold event counts
           150 profile-permutation draws each, R1 pinned at that y's own observed value
           the resulting p95 applied to the 100 alternative draws already cached
           300 draws, ~70 min

                                 null p95      corrected power
    synthetic y0                 0.010834             0.55
    synthetic y1                 0.014090             0.40
    mean threshold               0.012462             0.45

    threshold spread between the two y's        1.30x
    recorded instrument (label-randomisation)   p95 0.007374  ->  power 0.640
```

**Corrected power at the pre-registered alternative: 0.45.**

The truth sits between the two flawed estimates of §3.8 and much nearer the optimistic one. The
`0.08` recomputation was too pessimistic: it borrowed the *real* data's null threshold (`0.028844`),
but under a synthetic alternative the correctly-constructed threshold is only `0.011`–`0.014`.

### What this does to the decision

```text
   prior     power 0.08     power 0.45     power 0.64
    0.10           0.15           0.50           0.59
    0.25           0.35           0.75           0.81
    0.50           0.62           0.90           0.93
```

Positive predictive value of the observed result, at alpha = 0.05. The decision-relevant band
narrows from `[0.35, 0.81]` to roughly **`[0.75, 0.81]`** at a sceptical 0.25 prior. That is a range
one can act on.

### A further flaw, found in the same pass and NOT corrected

`_assign()` places synthetic positives uniformly at random within each outer fold, **ignoring
replicate**. The real cohort is severely imbalanced — 10.14% prevalence in replicate 3 against 1.42%
in replicate 2 — so the synthetic cohorts do not reproduce the real one's structure. This is
inherited from 23.2E, where a single replicate made it moot.

Consequently `0.45` is a **better-constructed approximation, not an exact figure**. It is the number
to report. It is not a number to defend to three decimal places.

### Status

Reported, not gating. Gate 18.3 was evaluated and recorded under the frozen V5 instrument and its
FAILED verdict stands unchanged at 0.64. The corrected figure is recorded beside it. Nothing is
re-scored, and 0.45 still fails the 0.80 threshold.

---

# 4. What is and is not established

**Established.** Pretreatment transcriptional state predicts the Rewind priming outcome beyond a
depth-complete technical baseline, on two biological replicates that were never used to design the
correction, under outcome rules that reproduce the authors' own objects exactly, with the direction
positive in both replicates and in both transfer directions.

**Not established, and not reportable:**

```text
  a confirmed Role-A claim              power 0.64 against a pre-registered 0.80
  an effect size                        a significant result from an underpowered design is
                                        selected on significance and overstates
  which replicate carries the signal    flips on a one-clone difference
  a uniform endpoint across replicates  three replicates, three source-defined operationalisations
```

---

# 5. The decision

## Option A — amend the opening rule, open Stage 24 as a development stage

Add a roadmap exit `ROLE_A_SUPPORTED_UNDERPOWERED`: gates 18.1/18.2/18.4/18.5/18.6 pass, 18.3 fails,
and Stage 24 may open **as a development stage only**, with the limitations propagating.

**What this means concretely:**

```text
  Stage 24        opens. Architecture selection and model development proceed.
  gate 18.3       STAYS FAILED, on the record, permanently. Nothing is lowered or re-scored.
  the anchor      Role A is cited as UNDERPOWERED SUPPORTING EVIDENCE, never as "confirmed"
  the paper       the Role-A claim is reportable with its limitation; no effect size is quoted
  Stage 27        remains a HARD confirmation gate, unchanged, on an UNUSED system
  cost            replicates 2 and 3 are consumed and INELIGIBLE for Stage 27 -- that stage now
                  needs a genuinely different lineage-resolved system
  risk            if Role A is noise (~19% at a 0.25 prior), Stage 24 develops against a weak
                  anchor. Stage 27 is what catches that.
```

**Why it is defensible:** the roadmap already defines Stage 24 as *"a development / architecture-
selection stage, not the independent biological replication stage… Stage 27 preserves that role."*
Stage 24 never made a confirmatory claim. Opening it on underpowered evidence therefore does not
launder that evidence into a confirmation — the confirmation burden stays where it always was.

**Why it is not the same as lowering the threshold:** lowering 0.80 would assert *"this design was
adequately powered."* That is false. The amendment asserts *"we choose to develop on underpowered
supporting evidence, and every downstream claim carries that caveat."* That is true and auditable.

## Option B — keep Stage 24 blocked, revise the claim architecture

Accept `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE`, and make Role B the primary anchor through an
explicit roadmap revision.

**What this means concretely:**

```text
  Stage 24        stays blocked on the Role-A path
  Role A          reported as a diagnosed failure plus an underpowered positive; not an anchor
  Role B          becomes primary. ROLE_B_ADDITIVE_PASS and INTERACTION_PASS_MULTI_TREATMENT
                  already hold, so this is not a hypothetical
  the paper       narrower. The prospective benchmark (contribution 1) survives intact; the
                  treatment-conditioned model becomes Role-B-only
  cost            a replicated positive result is set aside despite passing 5 of 6 gates
  Stage 27        unchanged
```

## Option C — acquire more Role-A evidence

**Rejected, on your reasoning and I agree.** A new dataset would not carry the thing that makes this
one trustworthy: reconstructions verified bit-for-bit against the authors' own objects, a
mechanically proven replicate map, and a documented bug in their normalisation. A better number we
understand less is a worse position than a marginal number we understand completely.

---

# 6. Recommendation

## 6.1 Option A — open Stage 24 as a development stage

I withdrew my first recommendation when the audit showed the power figure could not be trusted, and
said the band had to be narrowed before anyone decided. §3.10 narrowed it. The recommendation now
stands on the measured number rather than on the broken one.

**Corrected power 0.45 at the pre-registered alternative**, which puts the positive predictive value
of the observed result at about **0.75 at a sceptical 0.25 prior**, and 0.90 at an even one.

```text
  five of six gates pass, and the failing one is a self-imposed planning standard
  p_perm 0.040 from a valid permutation test at n = 2000
  direction positive in BOTH replicates and BOTH cross-replicate transfer directions
  outcome rules reproduce the authors' own serialized objects exactly
  Stage 24 is DEVELOPMENT by the roadmap's own definition; Stage 27 is the real gate
  Stage 27 remains fully intact on a system used nowhere upstream
```

The asymmetry that decides it: a false positive here costs wasted development effort, and Stage 27
catches it. Declining a result that passed 5 of 6 gates and replicates in direction across two
independent biological replicates costs a genuine finding.

## 6.2 What must ride along, without exception

```text
  gate 18.3 STAYS FAILED on the record, at 0.64 under the frozen instrument, with the
    corrected 0.45 recorded beside it. Neither is re-scored.
  the Role-A anchor is UNDERPOWERED SUPPORTING EVIDENCE, never "confirmed"
  no effect size is quoted -- +0.0305 is a direction, not a magnitude
  no claim about which replicate carries the signal -- it flips on one clone
  the endpoint is NOT uniform across replicates: three replicates, three source-defined rules
  the power instrument's two known flaws (wrong null; synthetic positives ignore replicate)
    are reported wherever the power figure is
  replicates 2 and 3 are CONSUMED and ineligible for Stage 27
```

## 6.3 If you would rather not

**Option B remains defensible and I would not argue against it.** A 0.75 posterior at a sceptical
prior is a reasonable thing to build on and also a reasonable thing to decline. The difference
between the two options is appetite for risk on a development stage, not a disagreement about what
the evidence says.

**What I will not do under either option:** move the 0.80 threshold, evaluate gate 18.3 at any
alternative other than 0.66, quote `+0.0305` as an effect estimate, claim which replicate carries
the signal, or describe the Role-A result as confirmed.
