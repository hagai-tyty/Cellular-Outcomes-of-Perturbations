# STAGE 23.2 — ROLE-A INDEPENDENT CONFIRMATION PROTOCOL, V5

**Supersedes V4.** V1–V4 remain historical in `arcive/` and are not edited. V5 exists because V4
was written under a factual premise that has since been disproved, and because two of its clauses
are mutually unsatisfiable against the evidence that actually exists.

**Frozen before execution.** Nothing in this document may be changed after any confirmatory
statistic is computed. The three decisions V5 makes — the outcome-rule resolution (§6), the
spike-in-bug resolution (§7) and the power gate (§9) — are frozen here, in advance, with their
justifications, precisely so that they cannot be chosen later to suit a result.

---

# 1. Why V5 exists

## 1.1 V4's factual premise was wrong

V4 was executed at step 1 (23.2G source qualification) and returned
`QUALIFYING_SET_EMPTY_FROM_FROZEN_SEARCH_SPACE`, on the finding that biological replicates 2 and 3
carry a pretreatment transcriptome and **no later outcome of any kind**. `stage_23_2G_RECORD.md`
records that verdict and stands unedited.

The finding was wrong, and it was wrong for a locatable reason: the qualification enumerated GEO
supplementary files, and the outcome materials for replicates 2 and 3 were never deposited in GEO.
They are public, but they live in the authors' Dropbox and on Figshare. `GSE227151_family.xml`
cannot see them, so a GEO-complete search returned a true answer to the wrong question.

The correction is recorded in `RECORDs/stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md`, which re-derives
every load-bearing fact from the primary files rather than from the accompanying third-party manual.

## 1.2 V4 §15.3 and §18.2 are unsatisfiable as written

§15.3 requires the reconstruction to be **source-faithful** and, in the same breath, mandates
`slice_max(n = 100, with_ties = TRUE)` while forbidding any change to `N` or the tie rule. That was
coherent only while R1 was the sole replicate, because top-100-with-ties **is** R1's source rule.
It is not replicate 2's and it is not replicate 3's. Once a second replicate enters, "source-faithful"
and "R1's parameters" name different things and cannot both be obeyed. §6 resolves this.

## 1.3 V4 §18.3 imports a floor derived under a different geometry

§18.3 requires ≥140 total positive clones. The available same-claim total is 76. But §18.3's own
authority, V4 §10.1, states that the 23.2E ladder estimates *"within-R1 event-count detectability
under the empirical R1 covariate distribution and nothing else."* The confirmation cohort is not
R1 and does not share its geometry:

```text
                       clones   positives   prevalence
  R1  (Stage 22)         3147          35        1.11%
  confirmation cohort    2310          76        3.29%
```

Power at fixed oracle AUC is a function of prevalence, cohort size and covariate structure — not of
positive count alone. Transplanting R1's rung onto a cohort with ~3× the prevalence is a category
error in both directions: it could as easily be too lenient as too strict. §9 replaces the imported
constant with the same measurement, made on the cohort that will actually be analysed.

**This is not the removal of a gate.** It is the same gate, evaluated against the right
distribution, pre-registered before any confirmatory statistic exists, and fully capable of failing.

---

# 2. Carried forward from V4, unchanged

The following are inherited verbatim and are **not** reopened by V5:

```text
§9    the corrected confirmatory hypothesis (depth_complete_nuisance_control)
§10   the 23.2E power ladder as a HISTORICAL within-R1 result
§11   the seven source-qualification criteria
§12   the search budget
§13   the forbidden / already-inspected data list
§14   the Stage-27 firewall
§15.1 the definition of an independent biological outcome unit
§15.2 unit structure must be established from source materials before outcome values are read
§15.4 what is recorded per unit
§16   the pooled primary test and per-replicate secondary test
§18.1 / 18.1b  >= 2 independent non-R1 BIOLOGICAL REPLICATES, counted as replicates not units
§18.4 the frozen dual gate: observed > null p95 AND p_perm <= 0.05
§18.5 delta_AP positive in every qualifying replicate analysed separately
§18.6 no un-re-gated material benchmark change
```

The corrected hypothesis is restated here for completeness because everything else refers to it:

> **depth_complete_nuisance_control** — Under the historical Rewind outcome and the frozen
> Stage-22/23 evaluation geometry, pretreatment transcriptional state predicts the Role-A outcome
> beyond a depth-complete nuisance baseline
> `Bdepth = [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI), log1p(n_detected_GE_features_in_raw_pseudobulk)]`.

---

# 3. The confirmation evidence, as verified

Every quantity below was re-derived from primary files. Nothing is taken from the third-party
`Rewind_File_Connection_and_Usage_Manual.docx`, which is treated as an unverified secondary source.

## 3.1 The replicate map is mechanically proven

Each `cellID` in each author lineage table was looked up in each GEX sample's barcode vector:

```text
             GEX S1     S2     S3     S4     S5
  table S1   100.0%   0.1%   0.0%   0.1%   0.1%
  table S2     0.2% 100.0%   0.0%   0.0%   0.0%
  table S3     0.0%   0.3% 100.0%   0.0%   0.0%
  table S4     0.1%   0.1%   0.0% 100.0%   0.1%
  table S5     0.1%   0.0%   0.0%   0.1% 100.0%
```

100% on the diagonal, ~0.1% off it, for both the R2- and R3-folder tables independently.

```text
  S1 -> GSM7092519   biol rep 3   ungated control
  S2 -> GSM7092520   biol rep 3   sorted, fast-proliferating
  S3 -> GSM7092521   biol rep 3   sorted, slow-proliferating
  S4 -> GSM7092517   biol rep 2   ungated control
  S5 -> GSM7092518   biol rep 2   ungated control
```

Independently corroborated: `primedCellsInd.rds` carries `SampleNum ∈ {S4, S5}` and nothing else.

37 raw 16-mers occur in more than one GEX sample. **Every join must carry the sample key.**

## 3.2 The outcome-unit structure (V4 §15.2, now satisfied)

From `metadata_fileShare_gDNABarcodes.xlsx`, and confirmed against the eight gDNA library keys
present in the shared master:

```text
  biological replicate 3    FS_1A/FS_1B, FS_2A/FS_2B, FS_3A/FS_3B      3 units, each an A/B split
  biological replicate 2    LSD1_4A / LSD1_4B                          1 unit, two treatment arms
```

Replicate 2's two libraries are **two treatment arms, not a technical split**: one control, one
LSD1-inhibited. The workbook contradicts itself about which is which (row name and treatment column
disagree in one row). That ambiguity is **inert** here, because the selection statistic actually used
is symmetric in the two libraries (§5.1). It is recorded, not resolved, and it is not guessed.

## 3.3 The folders are not the replicates

The R2- and R3-named `stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt` files are byte-identical
(MD5 `1c39b05a21da8219301ab228a81e72c1`), and both `filtered10XCells.txt` tables span S1–S5. They
are two analysis filters over one shared master, not two replicate extracts. Replicate identity
comes from `SampleNum` (§3.1), never from folder location. The shared master counts **once**.

---

# 4. Frozen eligibility

```text
  ELIGIBLE (primary confirmation)
    biological replicate 2    S4 (GSM7092517) + S5 (GSM7092518)    ungated control
    biological replicate 3    S1 (GSM7092519)                      ungated control

  EXCLUDED
    biological replicate 3    S2 (GSM7092520)  sorted fast
    biological replicate 3    S3 (GSM7092521)  sorted slow
```

S2 and S3 are excluded because sorting on proliferation speed changes the **pre-state population**,
and the Role-A claim is about an unselected pretreatment state. The Stage-22 R1 benchmark is ungated;
admitting sorted cells would compare a different population against it.

They are excluded **before** any statistic is computed, on population grounds alone. They may not be
admitted later. If a future protocol wishes to use them it must declare a separate claim, and V4
§10.7 / §18.6 apply.

## 4.1 Linkage table assignment — frozen

```text
  biological replicate 2   R2/filtered10XCells.txt   restricted to SampleNum in {S4, S5}
  biological replicate 3   R3/filtered10XCells.txt   restricted to SampleNum == S1
```

Each replicate uses the author list **against which its own outcome object was constructed**. This
is a provenance rule, not a yield rule; it is fixed by which list the author's own code consumed,
and it may not be swapped for the other list.

## 4.2 Ambiguity exclusion — inherited from Stage 22 §3.5

Any `(SampleNum, cellID)` mapping to more than one clone is excluded, with the same audit record
Stage 22 emits. Clone tables and folds are built on the post-exclusion population.

---

# 5. Frozen outcome rules, per replicate

Both rules below reproduce the corresponding author object **exactly**. That exact reproduction is
the evidence that the rule is the author's and not ours.

## 5.1 Biological replicate 2 — `R2_MIN_PAIRED_TOP26_V1`

V4 recorded R2's rule as unavailable from shipped code. It was recovered by reconstruction and
verified against `primedCellsInd.rds`:

```text
  1. spike-in normalise LSD1_4A and LSD1_4B
       coef = lm(c(20000, 5000) ~ 0 + nUMI_spike)   per library, author form
  2. keep lineages detected in BOTH libraries                inner join   5035 -> 636
  3. keep lineages with at least one linked 10X cell                       636 -> 92
  4. rank by  min(nUMINorm_A, nUMINorm_B)  descending
  5. take the top 26
```

Verification:

```text
  re-derived coef LSD1_4A   0.007452741453   author stored nUMINorm.x / nUMI.x   0.007452741453
  re-derived coef LSD1_4B   0.04555769882    author stored nUMINorm.y / nUMI.y   0.04555769882
  ranks 1-26   all in the author set
  ranks 27+    none in the author set
  boundary     60.22 -> 33.34
  result       26 lineages, 79 cells, SampleNum in {S4, S5}   == primedCellsInd.rds exactly
```

`foldchange` in the object is `log2((nUMINorm_B + 1) / (nUMINorm_A + 1))`, ranges −3.62 to +4.84
with no threshold, and is therefore **not** the selection criterion. It is a recorded column.

Step 5 is a fixed count of 26, equivalently any threshold in `(33.34, 60.22]`. The two cannot be
distinguished from the shipped object and produce identical labels, so the ambiguity does not
affect the reconstruction. Recorded, not resolved.

`min(A, B)` is symmetric, which is why §3.2's control/LSD1i label ambiguity cannot affect the labels.

## 5.2 Biological replicate 3 — `R3_MAX_PAIRED_TOP200_UNION_V1`

Reconstructed from `20221021_R3_identifyingPrimedCellsByCutoff.R`:

```text
  1. spike-in normalise each of FS_1A/1B, FS_2A/2B, FS_3A/3B
  2. within each of the three units: full join A/B, nUMIMax = max(nUMINorm_A, nUMINorm_B)
  3. slice_max(nUMIMax, n = 200, with_ties = FALSE)          per unit, independently
  4. map selected lineages back to that replicate's linked cells
  5. union the three per-unit positive sets
```

`N = 200` is fixed by the authors' own downstream scripts, which set `i = 6` against
`cutoffList <- c(10,25,50,100,150,200,250,500,1000)`. The single later section using `i = 7`
(top-250) is not the primary definition.

Per-unit reconstruction and union — **not** re-ranking across units — is what the author code does
(`c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+9]]), unlist(primedCellIDList[[i+18]]))`)
and is what V4 §15.3 requires.

---

# 6. Resolution of the V4 §15.3 rule conflict

```text
  V5 RULE

  Each qualifying biological replicate is labelled by ITS OWN source-defined outcome rule,
  frozen from author code and author objects before execution, and declared in §5.

  R1's top-100-with-ties is R1's source rule. It is not a universal constant and it is not
  imposed on other replicates.
```

**Why this is the correct reading and not a relaxation.** §15.3's own heading is
*"Reconstruction is per unit, and source-faithful"*, and its forbidden list is aimed squarely at
outcome-shopping: pooling across replicates, re-ranking across units, and *"selecting a unit
definition because it yields more positives"*. Every one of those prohibitions survives V5 intact.
What V4 did was freeze R1's **instance** of the principle in the place where the **principle**
belonged, at a time when R1 was the only replicate and the two were indistinguishable.

**Anti-gaming constraints, retained in full.** The rules in §5 were fixed by reproducing the authors'
serialized objects exactly, before any confirmatory statistic was computed. They cannot have been
tuned to a result because no result existed. Additionally:

```text
  FORBIDDEN, without exception
    changing N, the tie rule, the ranking statistic or the join for any replicate after this
      document is committed
    substituting one replicate's rule into another replicate
    pooling outcome libraries across biological replicates before selection
    re-ranking barcodes across units
    selecting any rule, unit definition or cutoff because of how it moves positives or performance
    admitting the sorted samples S2/S3
```

## 6.1 Why this does not trigger V4 §10.7 / §18.6 re-gating

§10.7 protects the **existing frozen benchmark** from being redefined after failure. Nothing in the
Stage-22 Rewind benchmark changes:

```text
  R1's outcome rule            unchanged (top-100 with ties)
  R1's clone set and folds     unchanged
  Stage 22 / Stage 23 gates    unchanged, not re-run, not re-interpreted
  the historical Stage-23 Role-A FAIL   permanent, as always
```

V5 defines how **new, previously unlabelled replicates** are labelled from their own sources. It
adds a cohort; it does not edit one. The confirmation cohort is written as a **new, separately
versioned benchmark artifact** (`stage23_2h_confirmation_*`), never by mutating a Stage-22 file.

## 6.2 The declared limitation this creates

Replicates 1, 2 and 3 measure the same biological concept — a lineage's future reprogramming
priming — through three different source-defined operationalisations:

```text
  R1   top 100 on summed raw gDNA UMI, ties INCLUDED, one pooled library
  R2   top 26  on min(normA, normB), inner join, spike-in normalised, one two-arm unit
  R3   top 200 on max(normA, normB), ties EXCLUDED, spike-in normalised, three A/B units
```

This is a property of the source study, not a choice available to us. It must be reported as a
standing limitation wherever the confirmation result is reported, and it may not be described as
a uniform endpoint. Prevalence heterogeneity between the confirmation replicates (1.4% vs 10.4%)
is part of the same limitation and must be reported alongside any pooled statistic.

---

# 7. Resolution of the author spike-in indexing bug

`20221021_R3_identifyingPrimedCellsByCutoff.R` loops `i in 1:3` over `sampleList <- c(1,3,5)` but
scales with `lmr[[i]]` / `lmr[[i+1]]` instead of `lmr[[sampleList[i]]]` / `lmr[[sampleList[i]+1]]`.
Coefficients 5 and 6 are never used. Verified by reading the shipped `overalapTableList.rds` and
matching stored `nUMINorm / nUMI` ratios to 10 significant figures:

```text
             observed ratio        correct coef      actually used
  FS_1A      0.006731403886          0.0067314      lmr[[1]]  CORRECT
  FS_1B      0.01470233947           0.0147023      lmr[[2]]  CORRECT
  FS_2A      0.01470233947           0.0092321      lmr[[2]]  WRONG, is FS_1B's
  FS_2B      0.009232116866          0.0081997      lmr[[3]]  WRONG, is FS_2A's
  FS_3A      0.009232116866          0.0134957      lmr[[3]]  WRONG, is FS_2A's
  FS_3B      0.008199715901          0.0166534      lmr[[4]]  WRONG, is FS_2B's
```

The bug is not neutral: selection ranks on `max(nUMINorm_A, nUMINorm_B)`, the two sides receive
different scale factors, and mis-assigning them changes top-200 membership.

```text
  V5 DECISION, frozen before execution

  PRIMARY      corrected coefficients   R3_MAX_PAIRED_TOP200_UNION_V1
  SENSITIVITY  author coefficients      R3_MAX_PAIRED_TOP200_UNION_AUTHORBUG
```

**Why corrected is primary.** The confirmation tests a biological hypothesis, not the reproducibility
of a figure. `sampleList` exists in the author's own code for the sole purpose of indexing the right
library, which makes the intended rule unambiguous; `lmr[[i]]` is an indexing slip, not a modelling
choice. Carrying a demonstrably wrong per-library scale factor would inject a known technical error
into the outcome definition.

**Both are run and both are reported.** The sensitivity arm reproduces `primedCellIDList.rds`
exactly, which is also how the reconstruction is validated. The primary arm may **not** be swapped
for the sensitivity arm after results are seen, in either direction, and disagreement between them
is reported as a limitation rather than resolved by preference.

---

# 8. The confirmation cohort

```text
  ELIGIBLE (ungated, barcode-linked, post-ambiguity-exclusion)     cells   clones
    biological replicate 2   S4+S5   GSM7092517/18                  3480     1827
    biological replicate 3   S1      GSM7092519                       598      483
                                                            total   4078     2310
```

The eligible cohort is identical in both §7 arms: eligibility and linkage are outcome-free, so the
spike-in decision cannot move them. Positives differ by arm, and are frozen separately:

```text
  POSITIVES                                                       cells   clones
  PRIMARY      R2_MIN_PAIRED_TOP26_V1                    rep 2       79       26
               R3_MAX_PAIRED_TOP200_UNION_V1             rep 3       60       49
                                                         total      139       75

  SENSITIVITY  R2_MIN_PAIRED_TOP26_V1                    rep 2       79       26
               R3_MAX_PAIRED_TOP200_UNION_AUTHORBUG      rep 3       61       50
                                                         total      140       76
```

Counts are recorded here as the frozen expectation. The executable stage recomputes them from source
and **must** reproduce the numbers for the arm it is running, or halt.

> **Pre-execution correction, recorded rather than silently applied.** This table first carried
> `rep 3 = 61 cells / 50 clones` for *both* arms, because the figures had been transcribed from the
> authors' `primedCellIDList.rds` — which is the **sensitivity** arm's output, not the primary's.
> The executable §8 guard caught it on the first run of 23.2H-A and halted, before any statistic
> existed. The primary arm's true count is 60 cells / 49 clones: correcting FS_2's and FS_3's
> spike-in coefficients moves 14 and 30 of their 200 selected lineages respectively (FS_1 is
> unaffected, as expected — it is the one unit the author code scales correctly), and the net effect
> on rep 3's linked, ungated, S1-restricted clone set is one clone. Nothing else changed, and no
> confirmatory statistic had been computed at the time of the correction.

## 8.1 Fold construction

Outer folds are drawn fresh for this new cohort — R1's Stage-22 folds are not reused and not
touched — at **clone level**, `N_OUTER = 5`, stratified on `(biological_replicate, y_primed)` so
that every fold carries both replicates and a comparable event count.

```text
  seed   STAGE_23_2H_SPLIT_SEED = 23511
  unit   clone
  rule   StratifiedKFold(n_splits=5, shuffle=True, random_state=23511) over
         sorted(clone_id), stratified on f"{replicate}|{y_primed}"
```

Stratifying on `y` is legitimate here and does not repeat V3's error: this is a **cohort-construction**
step performed once, exactly as Stage 22 built R1's folds with `StratifiedKFold` on `y_primed`. What
V3 got wrong was requiring a fold table to be frozen *before outcome values were inspected* while
simultaneously stratifying on them. V5 does not claim that: outcomes are reconstructed first (§5),
folds are drawn from them second, and no model has been fitted at either point.

## 8.2 Replicate identity as a blocking nuisance covariate

V4 §16.1 point 3 is carried forward and made explicit here, because this cohort makes it
load-bearing in a way R1 alone never did. The nuisance design is Bdepth **plus a replicate
indicator**:

```text
  [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI),
   log1p(n_detected_GE_features), is_biological_replicate_3]
```

Replicate identity is a nuisance term only. It may not be a predictor of interest and it may not be
interacted with `X` — a state × replicate interaction is a different scientific claim (V4 §16.1).

**Why it matters here.** The two replicates have very different prevalence — 1.4% in replicate 2
against 10.1% in replicate 3. Expression carries replicate identity through ordinary batch
structure. Without the blocking term, a model that merely recognised which replicate a clone came
from would predict its outcome better than chance, and the pooled `ΔAP` would credit that to
transcriptional state. The blocking term puts that information in the **baseline**, where it
belongs, so `ΔAP` measures only what state adds beyond it.

## 8.3 Permutation strata

Inherited from Stage 23 §7.1, extended by replicate:

```text
  stratum = f"{size}|{n_lanes}|{replicate}"     size in {1, 2, 3+} by n_pretreatment_cells
```

Adding the replicate term is required, not optional, and it does two jobs. It prevents a profile
from crossing a replicate boundary, which would destroy the very independence the confirmation
exists to establish. And because a permuted profile stays inside its own replicate, whatever
replicate-identifying signal expression carries is present in the **null** as well as in the
observed statistic — so the test cannot be passed by that signal even if the §8.2 blocking term
failed to absorb all of it. The two protections are deliberately redundant.

---

# 9. The power gate — V4 §18.3 replaced

## 9.1 What is replaced, and by what

```text
  V4 §18.3    total positive clones >= 140
  V5 §9       the confirmation design must achieve power >= 0.80 at oracle AUC 0.66
              UNDER ITS OWN REALIZED GEOMETRY, measured by the 23.2E machinery
              applied to the actual confirmation cohort
```

The measurement instrument is unchanged. Only the cohort it is applied to changes — from R1
resampled with replacement, to the cohort that will actually be analysed.

## 9.2 Method — frozen

Identical in form to 23.2E, re-pointed at the confirmation cohort:

```text
  z            label-free synthetic direction: PC1 of the confirmation X residualised on Bdepth,
               deterministically oriented (largest |loading|, ties -> smallest feature index)
  positives    assigned by weighted sampling WITHOUT replacement inside each outer fold,
               at the cohort's REAL per-fold positive counts
  beta         bisection-calibrated so the oracle z has median AUC 0.66
  statistic    delta_AP = AP(R3) - AP(R1) from the frozen nested pipeline
  null         positives assigned at beta = 0
  power        fraction of alternative draws whose delta_AP exceeds the null p95
  n_sims       200 null, 100 alternative        (23.2E's allocation, unchanged)
  seeds        SEED_COVARIATE 23540, SEED_NULL 23541, SEED_ALT 23542, SEED_BETA 23543
```

`z` is a simulation generator only. It never scores the real outcome and never enters the fitted
evaluation pipeline. This is 23.2E's constraint, retained verbatim.

## 9.3 The decision rule — frozen before the measurement

```text
  power >= 0.80    §18.3 SATISFIED. The design is not disqualified on detectability.
  power <  0.80    §18.3 FAILS.
```

On failure, V4 §17's consequences apply unchanged: the confirmatory analysis is still run and still
reported in full; a positive result is recorded as **underpowered supporting evidence** and a null
result is **not** evidence against the hypothesis; the exit is `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE`
and Stage 24 stays blocked.

**The gate may fail, and no step below is conditioned on its passing.** The power study is executed
and its value recorded before the confirmatory permutation null is computed, so that the threshold
cannot be revisited afterwards.

## 9.4 What may still not be done to satisfy §9

```text
  enlarging the cohort by admitting the sorted samples S2/S3
  switching a replicate to the other author linkage list
  changing N, the tie rule or the ranking statistic in any §5 rule
  adding a replicate that failed §11 qualification
  counting R1 positives toward anything
  re-running the power study with different seeds and keeping the better number
  reporting SUPPORTED on a design whose measured power is below 0.80
```

---

# 10. Execution order — frozen

```text
  23.2H-A   build the confirmation benchmark
            reconstruct §5 outcomes, apply §4 eligibility, build cells/clones/folds/strata,
            reproduce the §8 counts or HALT

  23.2H-B   representation
            clone pseudobulk (raw sum -> CP10K -> log1p, exactly once) and Bdepth,
            using the frozen Stage-23 constructors

  23.2H-C   power  (§9)
            record power BEFORE any confirmatory statistic exists

  23.2H-D   confirmatory analysis
            observed delta_AP pooled and per replicate, then the permutation null
            (200 permutations, strata per §8.2), then the §18.4 dual gate

  23.2H-E   verdict and handoff
            evaluate §11, write the exit, update STAGE_23_2_HANDOFF_TO_STAGE_24
```

A stage may not be reordered. 23.2H-C may not be run after 23.2H-D.

---

# 11. Gates for `ROLE_A_CONFIRMATORY_SUPPORTED`

All six required. 18.1, 18.4, 18.5, 18.6 are V4's, unchanged.

```text
  18.1   >= 2 independent non-R1 BIOLOGICAL REPLICATES qualify under §11 and §15.2
  18.2   each qualifying unit reconstructed source-faithfully and independently under §5/§6,
         with no pooling before selection and no post-hoc rule change
  18.3   measured design power >= 0.80 at oracle AUC 0.66 under the realized geometry   (§9)
  18.4   the pooled primary test passes BOTH frozen gates:
             observed delta_AP > null p95   AND   p_perm <= 0.05
  18.5   delta_AP POSITIVE in every qualifying replicate analysed separately
  18.6   no material benchmark semantic changed without re-gating                        (§6.1)
```

## 11.1 Exits

The roadmap's four exits are unchanged; V5 adds none.

```text
  all six gates pass                    ROLE_A_CONFIRMATORY_SUPPORTED   -> Stage 24 may open
  18.1/18.2/18.6 pass, 18.3 fails,
    18.4/18.5 pass                      ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
                                        reported as underpowered supporting evidence
  18.4 or 18.5 fails                    ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
  the design is shown inadequate for
    the intended claim                  ROLE_A_REDESIGN_REQUIRED
  no defensible signal remains          ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL
```

A `ROLE_A_CONFIRMATORY_SUPPORTED` exit additionally requires the Stage-24 handoff contract to be
complete and benchmark-compatible, per the roadmap.

---

# 12. Stage-27 firewall — restated

Replicates 2 and 3 are consumed by this confirmation. They are **not** available as the untouched
Stage-27 replication set. Stage 27 must preserve an independent biological test of the eventual
frozen Stage-24 model, on a system not used here.

---

# 13. Record of what V5 changes

```text
  §6   V4 §15.3's "top-100 with ties" replaced by each replicate's own source rule
  §7   the author spike-in indexing bug: corrected primary, author-faithful sensitivity
  §9   V4 §18.3's imported >= 140 floor replaced by a measured design-specific power gate
  §4   eligibility frozen to ungated samples only
  §8   a new, separately versioned confirmation cohort with its own folds and strata
```

Everything else is V4. V4 remains readable, unedited, in `arcive/`.
