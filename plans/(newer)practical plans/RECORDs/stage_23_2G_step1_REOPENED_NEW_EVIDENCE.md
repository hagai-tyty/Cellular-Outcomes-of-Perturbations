# stage_23_2G step 1 — REOPENED on new evidence (independent verification)

**This file does not modify `stage_23_2G_RECORD.md`.** That record stands as written. It reached
`QUALIFYING_SET_EMPTY_FROM_FROZEN_SEARCH_SPACE` on the premise that biological replicates 2 and 3
have **no later outcome measurement of any kind**. That premise is now known to be **false**. It was
false because the qualification searched GEO only, and the outcome materials for replicates 2 and 3
were never deposited in GEO — they live in the authors' public Dropbox and on Figshare. The original
record is preserved unchanged; this file records what is now believed instead.

## Provenance of the new evidence

User-supplied drop at `D:\GSE227151_Rewind\NEW DATA\` (reorganised, see `MANIFEST.md` there),
accompanied by `Rewind_File_Connection_and_Usage_Manual.docx`, an AI-written guide. **Nothing below
is taken from that manual.** Every claim was re-derived from the primary files. Where this
verification and the manual disagree, both are stated.

## What was verified, and how

### 1. The five reserved GEX samples are present and structurally sound

```text
GSM7092517_2_4_control    barcodes 4468   features 36601   mtx 36601 x 4468   24,373,399 nnz
GSM7092518_2_5_control    barcodes 5157   features 36601   mtx 36601 x 5157   25,457,555 nnz
GSM7092519_3_1_control    barcodes 2722   features 36601   mtx 36601 x 2722   14,094,638 nnz
GSM7092520_3_2_fast       barcodes 3188   features 36601   mtx 36601 x 3188   17,003,651 nnz
GSM7092521_3_3_slow       barcodes 2422   features 36601   mtx 36601 x 2422   12,246,116 nnz
```

All 15 files pass `gzip -t`. Every matrix header agrees with its own barcodes and features files.
Common 36,601-feature space.

### 2. The SampleNum → GSM → biological replicate map is now MECHANICALLY PROVEN

Not inferred from folder names or titles. Every `cellID` in each author lineage table was looked up
in each GEX sample's barcode vector:

```text
             GEX S1     S2     S3     S4     S5
  table S1   100.0%   0.1%   0.0%   0.1%   0.1%
  table S2     0.2% 100.0%   0.0%   0.0%   0.0%
  table S3     0.0%   0.3% 100.0%   0.0%   0.0%
  table S4     0.1%   0.1%   0.0% 100.0%   0.1%
  table S5     0.1%   0.0%   0.0%   0.1% 100.0%
```

100% on the diagonal, ~0.1% off it, for **both** the R2- and R3-folder tables. Combined with the
`metadata_fileShare_cDNABarcodes.xlsx` library descriptions, the frozen map is:

```text
S1 -> GSM7092519  biol rep 3  ungated control
S2 -> GSM7092520  biol rep 3  sorted, fast-proliferating
S3 -> GSM7092521  biol rep 3  sorted, slow-proliferating
S4 -> GSM7092517  biol rep 2  ungated control
S5 -> GSM7092518  biol rep 2  ungated control
```

37 raw 16-mers occur in more than one GEX sample, so every join must carry the sample key.

### 3. The gDNA library → biological replicate map, from `metadata_fileShare_gDNABarcodes.xlsx`

```text
R3_1A / R3_1B   "biol rep 3, split 1A / 1B"   FS-1A / FS-1B
R3_2A / R3_2B   "biol rep 3, split 2A / 2B"   FS-2A / FS-2B
R3_3A / R3_3B   "biol rep 3, split 3A / 3B"   FS-3A / FS-3B
R2_control      "biol rep 2, control"          LSD1-4A     treatment column: "LSD1 inhibition"
R2_LSD1i        "biol rep 2, LSD1 inhibition"  LSD1-4B     treatment column: "none"
```

The shared master carries exactly these eight gDNA libraries and no others. This closes the
section-15.2 outcome-unit-structure criterion for both replicates.

**Replicate 2's two libraries are two treatment arms, not a technical A/B split.** The workbook
contradicts itself about which is which (row name vs treatment column are swapped in one of them).
That ambiguity is **inert** for the outcome definition actually used, because the selection
statistic is symmetric in the two libraries — see item 4.

### 4. Replicate 2's outcome rule was RECOVERED (the manual states it is unrecoverable)

`primedCellsInd.rds`: 79 rows, 79 distinct cellIDs, **26 distinct lineage barcodes**, `SampleNum`
values `{S4, S5}` only — i.e. exactly biological replicate 2, confirming the map a third time.

Spike-in calibration was re-derived from the master with the authors' own two spike barcodes and
their `lm(c(20000, 5000) ~ 0 + nUMI)` form:

```text
LSD1_4A coef 0.007452741453     author stored nUMINorm.x / nUMI.x = 0.007452741453
LSD1_4B coef 0.04555769882      author stored nUMINorm.y / nUMI.y = 0.04555769882
```

Agreement to 10 significant figures fixes `.x = LSD1_4A`, `.y = LSD1_4B` and proves the calibration
is reproduced exactly. `foldchange` is `log2((nUMINorm.y + 1) / (nUMINorm.x + 1))`; it ranges
-3.62 to +4.84 with no threshold, so it is a recorded column, **not** the selection criterion.

The selection rule reproduces the author's 26-lineage set **exactly**:

```text
1. spike-in normalise LSD1_4A and LSD1_4B
2. keep lineages detected in BOTH libraries          (inner join: 5035 -> 636)
3. keep lineages having at least one linked 10X cell (636 -> 92)
4. rank by  min(nUMINorm_A, nUMINorm_B)   descending
5. take the top 26
```

Ranks 1–26 are all in the author set and ranks 27+ are all outside it, with a clean 60.22 → 33.34
gap at the boundary. Steps 1–4 are certain; step 5 is a fixed count of 26 or equivalently any
threshold in `(33.34, 60.22]` — the two cannot be distinguished from the shipped object alone, and
they produce identical labels, so the ambiguity does not affect the reconstruction.

`min(A, B)` is symmetric, which is why the control/LSD1i label swap does not matter.

### 5. The author's published R3 object carries a SPIKE-IN INDEXING BUG

`20221021_R3_identifyingPrimedCellsByCutoff.R` iterates `i in 1:3` over `sampleList <- c(1,3,5)` but
scales with `lmr[[i]]` and `lmr[[i+1]]` instead of `lmr[[sampleList[i]]]` and
`lmr[[sampleList[i]+1]]`. Coefficients 5 and 6 are never used.

Read directly out of the shipped `overalapTableList.rds` and matched to 10 significant figures:

```text
             observed nUMINorm/nUMI     correct coef      coef actually used
  FS_1A          0.006731403886          0.0067314      lmr[[1]]  CORRECT
  FS_1B          0.01470233947           0.0147023      lmr[[2]]  CORRECT
  FS_2A          0.01470233947           0.0092321      lmr[[2]]  WRONG (is FS_1B's)
  FS_2B          0.009232116866          0.0081997      lmr[[3]]  WRONG (is FS_2A's)
  FS_3A          0.009232116866          0.0134957      lmr[[3]]  WRONG (is FS_2A's)
  FS_3B          0.008199715901          0.0166534      lmr[[4]]  WRONG (is FS_2B's)
```

The bug is **not neutral**: selection ranks on `max(nUMINorm.x, nUMINorm.y)`, and the two libraries
in a pair receive different scale factors, so mis-assigning them changes which lineages enter the
top-N. Reproducing the authors' published labels requires reproducing the bug; correcting it
produces a different label set. This is a V5 decision and it is not mentioned in the manual.

### 6. Cohorts and available positives

R3's top-200 sets are `primedCellIDList` elements 6, 15, 24 (ladder length 9, pair-major order),
carrying 27 / 75 / 83 cells; downstream author scripts union the three.

```text
ELIGIBLE (ungated, barcode-linked)        cells   clones
  rep 2   S4+S5   GSM7092517/18            3480     1827
  rep 3   S1      GSM7092519                598      483

POSITIVES under each replicate's own source rule
  rep 2   author primedCellsInd, S4+S5        79       26
  rep 3   author top-200 union, S1 only       61       50
```

## Status against the V4 gates

```text
11    source qualification          BOTH rep 2 and rep 3 now PASS all seven criteria
15.2  outcome-unit structure        ESTABLISHED for both, from source materials
18.1  >= 2 independent non-R1
      biological replicates         PASSES -- 2 (rep 2 = 1 unit; rep 3 = 3 units)
18.3  >= 140 total positive clones  FAILS
18.2  source-faithful under the
      frozen top-100-with-ties rule CONFLICTED -- see below
```

### 18.3 arithmetic

```text
  V4 §15.3 literal top-100    ungated only   26 + 32 =  58   FAIL
  R3 source-faithful top-200  ungated only   26 + 50 =  76   FAIL
  R3 source-faithful top-200  incl. sorted   26 + 94 = 120   FAIL, and breaks same-claim
```

No same-claim configuration reaches 140. Per V4 §17 the analysis may still be run and must be
reported in full, but it cannot emit `ROLE_A_CONFIRMATORY_SUPPORTED`; the exit stays
`ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE` and Stage 24 stays blocked on that gate alone.

### 18.2 / 15.3 conflict

V4 §15.3 mandates `slice_max(n = 100, with_ties = TRUE)` per unit and forbids changing N or the tie
rule. Neither replicate's source rule is that rule:

```text
  R1  top 100, ties INCLUDED, one pooled library, raw summed UMI
  R2  top 26 on min(normA, normB), inner join, spike-in normalised
  R3  top 200 on max(normA, normB), ties EXCLUDED, spike-in normalised, 3 A/B pairs
```

"Source-faithful" and "the frozen top-100-with-ties rule" cannot both be satisfied. This is a
protocol contradiction, not a data problem, and it requires a V5 decision.

## Bugs / breakage found outside the new data

`D:\GSE227151_Rewind\` no longer holds `filtered10XCells.txt`,
`stepThreeStarcodeShavedReads_BC_10X.txt` or `stepThreeStarcodeShavedReads_BC_gDNA.txt` at its root;
they were moved into `r1\`. `experiments/diag_stage21d_public_reconstruction.py::rewind_required`
reads them from the root, so the frozen R1 loader path is currently broken. Not repaired here, since
it touches frozen benchmark inputs.

## What did NOT change

- No Stage-23 or Stage-23.2 artifact was modified. `stage_23_2G_RECORD.md` is untouched.
- No model was fitted. No performance quantity was computed for any reserved candidate.
- No `X_before` matrix was loaded beyond reading barcode/feature/header dimensions.
- 23.2G remains STOPPED for review. Steps 2–5 of the V4 §19.0 order were not attempted.
