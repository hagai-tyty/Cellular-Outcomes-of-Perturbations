# stage_21D_ RECORD

## Goal
Acquire and reconstruct the two datasets Stage 21C qualified, without modelling: build
`pretreatment expression -> cell -> clone -> future outcome` for `GSE227151` (Rewind, Role A) and
`GSE279162` (WM989 six-treatment, Role B), and decide whether Stage 22 can open.

## Inputs
- `D:\GSE227151_Rewind\` — `filtered10XCells.txt`, `stepThreeStarcodeShavedReads_BC_10X.txt`,
  `stepThreeStarcodeShavedReads_BC_gDNA.txt`, `GSE227151_family.xml`, series matrix,
  `GSM7092515` + `GSM7092516` 10x triplets — **11 required paths, all present**
- `D:\GSE279162\` — 9 samples × (barcodes / features / matrix) + family XML + series matrix —
  **29 required paths, all present**
- frozen commits: 21A `f6a0056`, 21B `3b8b644`, 21C `cab27d8`, roadmap correction `6e8a362`
- plan: `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §11 (21D), §13 (reproduction gate), §14
  (leakage guards)
- model: `_s16` frozen, never loaded
- no download performed — every byte was already on `D:\`; nothing raw entered git

## Files added
- `experiments/diag_stage21d_public_reconstruction.py`
- `tests/test_diag_stage21d_public_reconstruction.py` (27 tests)
- `results/diag_stage21d_public_reconstruction_results.json`
- `results/stage21d_rewind_clone_table.tsv` (3,149 clones, 194 KB)
- `results/stage21d_gse279162_clone_table.tsv` (4,018 naive clones × 9 samples, 126 KB)

## Files modified
- none — **additive only**. 21A / 21B / 21C results and records untouched (pinned by test)

## What changed
- Both prospective linkages are now reconstructed from the actual files, not from any summary
- Provenance for all 40 raw inputs is committed as path + size + SHA-256; the ~1.1 GB of matrices
  stays outside git
- Tri-state rule carried over unchanged: `PRESENT` / `ABSENT_PROVEN` /
  `UNKNOWN_REQUIRES_SOURCE_FILE`
- Verdict is derived from the findings, not asserted: an `UNKNOWN` on the outcome rule *forces*
  `RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE`, and a test pins that coupling

## What did NOT change
- `src/` unchanged (verified: `git diff --name-only src/` empty)
- model unchanged — nothing fitted; no sklearn / torch / statsmodels imported (asserted by AST walk,
  not by string search, so the docstring may name them)
- labels unchanged
- no threshold was moved to reproduce a published number

## Tests
- 1670 passed
- ruff clean (CI scope: `src/ tests/ scripts/ plan_tests/`)
- every headline number below was re-derived independently in the shell (`sort`/`comm`/`awk`,
  no pandas) and matched exactly

## Result

**OVERALL: `STAGE_22_PENDING_OUTCOME_RULE`** — Stage 22 does **not** open yet.

### A. `GSE227151` (Rewind, Role A) → `RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE`

The linkage is complete and the label is threshold-free; only the author's exact *definition* is
missing.

```text
expression cells                    13,665   (GSM7092515 6,569 + GSM7092516 7,096)
cells with author-QC clone           3,913   (28.6%)  from filtered10XCells.txt, 3,921 rows
unique clones                        3,149
positive (primed) clones                82
negative clones                      3,067   positive rate 2.6%
positive cells                         102
negative cells                       3,819
ambiguous rows (nLineages > 1)          16
sample/replicate structure           2 10x lanes of ONE biological replicate ("biol rep 1")
```

- `primed` := the clone's barcode is recovered from the post-reprogramming gDNA arm (`SampleNum 3`)
- the gDNA table is **pre-collapse**: 49,554 rows over 1,936 distinct barcodes, so counts are summed
  per clone
- summed gDNA reads are strongly **bimodal** — 1,309 clones at 1 read, 153 at 2, then a colony mode
  with 357 clones ≥ 500
- **every one of the 82 primed clones carries ≥ 562 summed reads.** None sits in the noise mode, so
  the primed set is *identical* for any read floor in `[1, 562]` — **no threshold had to be chosen**
- 311 clones appear in **both** 10x lanes → a lane-wise split would leak clones; the outer split
  unit must be the clone
- per lane: `SampleNum 1` 50 cells / 42 clones · `SampleNum 2` 52 cells / 46 clones

**Source-study reproduction gate: FAILED, and recorded rather than chased.** The published figure of
42 primed cells is not reproduced by *any* cell-level reading (pooled 102, lane 50, lane 52). It
coincides with the *clone* count of one lane (42). A coincidence at this resolution does not
establish the rule, so it was recorded and **not adopted**.

Unresolved (`UNKNOWN_REQUIRES_SOURCE_FILE`): is `primed` per clone or per cell · are the two lanes
pooled or is one the anchor · was any gDNA floor applied (immaterial over `[1, 562]`, but unstated).
The blocking input is named: `D:\GSE227151_Rewind\author_code_zenodo7707418\` **exists but is
empty**, and GEO's Data-Processing block documents only the cellranger scRNA path — it says nothing
about barcode calling. This is pending a code drop, **not** pending more data.

### B. `GSE279162` (WM989, Role B) → `RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE`

The design is exactly what Role B needs, but every count is floor-conditional.

```text
feature structure          189,656 rows = 36,601 genes + 153,055 'Custom' LinNNNN features
                           the clone id is a FEATURE ROW inside the same matrix, not a side file
treatments                 Dabrafenib Trametinib CoCl2 Acid Cisplatin Doxorubicin
pretreatment samples       Naive1 Naive2 Naive3
naive clones (floor 1)     4,018      (floor 2: 1,969 · floor 5: 1,825 · floor 50: 781)
clones in >= 2 treatments  1,500
clones in all 6            250
naive cells per clone      median 1, max 262; 1,645 clones have >= 2
```

- design confirmed from GEO Overall-Design: barcoded WM989 cells doubled → a subset taken as the
  untreated naive arm → remainder split six ways (4 wk dabrafenib / trametinib / CoCl2 / acid, or
  2 wk cisplatin / doxorubicin + 2 wk holiday) → re-sequenced. **X precedes U precedes Y by
  construction.**
- positives per treatment at floor 1: Acid 1,357 · Doxorubicin 1,349 · Cisplatin 1,126 ·
  CoCl2 834 · Dabrafenib 779 · Trametinib 736
- **per-cell clone assignment is not clean**: the dominant barcode carries a median of only ~50% of
  a cell's lineage UMIs, and the median cell shows 6 lineage features. A per-cell call requires an
  explicit dominance + UMI rule
- unlike the Rewind gDNA arm there is **no empty region** separating signal from noise, so the naive
  pool falls 4,018 → 781 across the floor sweep. No floor was selected

Unresolved (`UNKNOWN_REQUIRES_SOURCE_FILE`): per-cell minimum UMI and minimum dominant fraction ·
the clone-presence floor for "surviving" · whether Naive1/2/3 pool or act as replicates.

## Bugs found
1. **`SampleNum` is transposed relative to the GEO titles.** GEO calls `GSM7092515` "biol rep 1,
   sample 1", but the barcode tables' `SampleNum 1` is `GSM7092516`. Resolved by cellID containment
   against each GSM's own `barcodes.tsv` — 1878/1878 (100%) vs 8/1878 (0.4%), and 2035/2035 vs
   3/2035 — which is unambiguous. Trusting either label would have attached every cell's expression
   to the wrong lane and silently corrupted the whole Role-A table. Pinned by test
2. **The two intersection figures carried into this stage do not reproduce.** The mapping direction
   is confirmed exactly, but not the counts. Measured on disk: `GSM7092516` has **7,096** 10x
   barcodes and `GSM7092515` has **6,569**; the clone-assigned intersections are **1,878/1,878** and
   **2,035/2,035**. The figures `6479/6479` and `6189/6189` match nothing in these files. They are
   left as stated in the incoming brief; this record carries the measured values
3. **The Rewind author-code directory is empty.** `author_code_zenodo7707418\` is present with zero
   files, which is precisely why the outcome rule is `UNKNOWN` rather than resolved. Named as a
   specific missing input rather than folded into a vague "cannot determine"
4. A house convention was missed on the first pass: every results-writer must define
   `_RESULTS = Path(__file__).resolve().parents[N] / "results"`. `tests/test_results_paths.py`
   caught it before commit

## Scientific interpretation
**Proves:** both qualified datasets carry a real, reconstructable `X_before + U -> Y_future`
geometry, verified from the files rather than from any summary. Rewind's primed/nonprimed label is
**threshold-free** over `[1, 562]` — an unusually clean property for a presence-based outcome, and
the strongest single result of this stage. `GSE279162` genuinely supports the interaction test it
was qualified for: 1,500 clones are re-observed under two or more treatments and 250 under all six,
so `does treatment preference depend on starting state?` is answerable in principle.

**Does NOT prove:** that either outcome column is the authors'. Both stages of the label — which
cells belong to a clone, and which clones count as surviving — rest on rules that the public files
do not state. For Rewind that gap is narrow and enumerable (three readings, all sharing the same
82-clone set); for `GSE279162` it is wide, because the counts move by 5× across the floor sweep.
Nor does it prove resolvability at the replicate level: only **one** biological replicate of Rewind
is locally reconstructable, so clone-level inference is available but replicate-level generalisation
is not. And it proves nothing about learnability — no model was fitted.

## Next action
Resolve the two outcome rules before Stage 22 opens. For Rewind that means obtaining the Zenodo
`7707418` code drop (the local directory is empty) or the paper's barcode-calling methods; for
`GSE279162` it means the authors' clone-calling and survival thresholds. Both are *rule* gaps, not
data gaps — no further dataset search is warranted, and 21C's `FULL_DATA_PATH` verdict stands
unchanged.

---

# ADDENDUM — 2026-08-21 — source-code resolution (Stage 21D revision 2)

**Everything above this line is the original Stage-21D record, frozen at `30ca7f0`, and is left
exactly as written.** It was correct on the evidence available to that run. This addendum records
what changed when the author code arrived.

## Original verdict, and why it was reasonable then

```text
GSE227151  RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE
GSE279162  RECONSTRUCTABLE_PENDING_EXACT_OUTCOME_RULE
OVERALL    STAGE_22_PENDING_OUTCOME_RULE
```

The required author-code directories were not available to that run.
`D:\GSE227151_Rewind\author_code_zenodo7707418\` existed but was **empty**, and GEO ships
`GSE279162` as a raw feature-barcode matrix whose Data-Processing block states only that barcode
reads were linked to 10x cell barcodes — not how a clone is called. Both outcome rules were
genuinely unknown, and the original run recorded them as `UNKNOWN` instead of inventing them.

## New evidence

The complete author code for both datasets was subsequently placed on `D:\`. It stays there; the
repository carries paths and SHA-256 digests only, for all seven scripts.

```text
D:\GSE227151_Rewind\author_code_zenodo7707418\plotScripts\rewind10X\R1\
    20220921_R1_primedVersusNonPrimedMarkersAndDistribution.R
    2022.02.14_R1_cellNumberDistributionForPrimedVersusNonPrimed.R
D:\GSE279162\author_code_Schaff_manuscript\
    How_to_run_code_README.txt · preprocess_GEX.Rmd · preprocess_cDNA_BCs.Rmd
    preprocess_gDNA_BCs.Rmd · Find_Markers_Top_Res_lins_in_naive.Rmd
```

## Rewind resolution → `RECONSTRUCTION_PASS`

The exact rule, transcribed from both R1 scripts:

```text
filter(cellID == "dummy")
group_by(BC50StarcodeD8, SampleNum) -> summarise(nUMI = sum(UMI))
slice_max(order_by = nUMI, n = 100)                  # dplyr default with_ties = TRUE
inner_join(filtered10XCells, by = "BC50StarcodeD8")  -> primed
anti_join (filtered10XCells, by = "BC50StarcodeD8")  -> nonprimed
```

**Equivalence check first, as required.** Our separate `stepThreeStarcodeShavedReads_BC_gDNA.txt`
*is* the `cellID == "dummy"` arm: all 49,554 rows carry `cellID = "dummy"` (so the author's filter
is inert) and `SampleNum` is uniformly 3. **One real schema difference: the per-row count column is
named `counts` here and `UMI` in the author script.** The semantics were not assumed equal — they
are confirmed downstream by the 42-cell reproduction.

```text
grouped barcodes                     1,936
slice_max cutoff (100th largest)     nUMI = 2,365
barcodes tied at the cutoff          2      -> 101 selected, not 100 (with_ties = TRUE)
selected nUMI range                  2,365 .. 18,735

primed cells                         42     <- reproduces the published 42
primed unique clones                 35
nonprimed cells                      3,879
nonprimed unique clones              3,114
by lane   SampleNum 1 (GSM7092516)   24 cells / 21 clones
          SampleNum 2 (GSM7092515)   18 cells / 17 clones
```

Nothing was tuned. `n = 100` is fixed by the author script; **42 was asserted after the rule was
implemented, never targeted.** The tie was checked both ways and the cell count is 42 either way,
which is reported rather than relied on.

**The excluded barcode `ATTCTAGTTGTAGTACGAGTAGCACATGTTCTACGTGGAGGACGAGAACG` is ABSENT from
`filtered10XCells.txt` and from the gDNA arm, so the second script's exclusion is INERT** — because
the upstream QC that produced `filtered10XCells.txt` had already removed exactly it. Both author
scripts give identical primed sets.

## GSE279162 resolution → `RECONSTRUCTION_PASS`

Executed in the README's order. Condition-specific RNA QC from `preprocess_GEX.Rmd`; the `Custom`
features become the `lineage` assay; `naive1/2/3` collapse to one `naive` condition — which answers
one of the original record's three open questions outright. Then `preprocess_cDNA_BCs.Rmd`:
drop zero-count lineages → `barcode_clustering(cell_lower_limit = 100, cor_threshold = 0.55)` →
`barcode_combine` → `barcoding_posterior` → `barcoding_assignment(difference_val = 0.2)` →
`assigned_lineage = NA` where `assigned_posterior < 0.5`.

```text
raw cells                    77,417
post-QC cells                46,891
lineages after rowSums > 0   14,883
clustering                   604 candidates -> 153 correlated pairs -> 90 clusters (204 lineages)
                             author stop("Merging happening") fired 0 times
lineages after combine       14,769
assigned cells               42,771      NA 4,120  (4,007 of them by the 0.5 posterior floor)
unique assigned clones       2,215

per condition (post-QC / assigned / NA / clones)
  naive          7,226 / 6,489 /   737 / 1,401
  Acid           8,004 / 7,246 /   758 /   950
  Cisplatin      7,390 / 6,815 /   575 /   346
  CoCl2          5,173 / 5,134 /    39 /   410
  Dabrafenib     6,556 / 6,164 /   392 /   509
  Trametinib     8,711 / 8,091 /   620 /   425
  Doxorubicin    3,831 / 2,832 /   999 /   662

clones with a naive observation      1,401
  re-observed in >= 1 treatment        929
  re-observed in >= 2 treatments       603
  re-observed in all 6                  37
naive cells per clone                median 2, max 88; 816 clones have >= 2
```

Future outcome, exactly as `Find_Markers_Top_Res_lins_in_naive.Rmd` builds it —
`table(assigned_lineage[OG_condition == condition])`:

```text
Y(clone, treatment) = post-treatment assigned-cell abundance, rankable within each treatment

              clones present / absent   assigned cells   median / max abundance
Dabrafenib          321 / 1,080             4,411            2 /   476
Trametinib          274 / 1,127             4,471            2 /   776
CoCl2               319 / 1,082             4,326            2 / 3,290
Acid                665 /   736             6,823            3 /   464
Cisplatin           272 / 1,129             6,717            1 / 2,791
Doxorubicin         405 /   996             2,282            2 /   577
```

The script's `num_lin = 5` is figure-specific and was **not** adopted as a Stage-22 target; a test
asserts no such binding exists in our code.

## What this corrects in the original record

| quantity | original (revision 1) | author rule (revision 2) |
|---|---|---|
| Rewind primed clones / cells | 82 / 102 | **35 / 42** |
| Rewind rule | present in the gDNA arm at any read count | top-100 barcodes by summed count |
| GSE279162 naive clones | 4,018 (UMI floor 1) | **1,401** |
| GSE279162 clones in ≥2 treatments | 1,500 | **603** |
| GSE279162 clones in all 6 | 250 | **37** |
| overall gate | `STAGE_22_PENDING_OUTCOME_RULE` | **`STAGE_22_READY`** |

The original record noted that Rewind's presence-based label was invariant for any read floor in
`[1, 562]`. That remains true — and remains beside the point, because the authors do not use a
floor. **Being invariant to the wrong knob is not the same as being right**, which is exactly why
the original run refused to call it a PASS. Revision 1's exploratory dominant-fraction / UMI-floor
sweep for `GSE279162` is **superseded and must not be used as the production clone call**; a test
asserts its numbers do not come back.

One correction to an over-general statement above: the original record said the author QC step
"drops 1,054 of 4,975 rows (21.2%)" without explaining it. It is precisely one hyper-abundant
barcode — `filtered10XCells.txt` **is** `stepThreeStarcodeShavedReads_BC_10X.txt` minus all 1,054
rows of `ATTCTAGT…GAGAACG`, verified as an exact set identity.

## Findings from the original record that still stand, unchanged

- Rewind `SampleNum` must be resolved by **cellID containment**, not by the GEO title text:
  `SampleNum 1 -> GSM7092516` (1878/1878), `SampleNum 2 -> GSM7092515` (2035/2035). Still
  transposed relative to GEO's "sample N" wording.
- 311 clones span both Rewind lanes, so a lane-wise outer split would leak. **The clone remains the
  biological grouping unit.**
- The `6479/6189` figures are still unsupported by the files and were **not** restored. Measured:
  `GSM7092516` 7,096 barcodes, `GSM7092515` 6,569; clone-assigned 1,878 and 2,035.
- Only **one** biological replicate of Rewind is locally reconstructable. Resolving the outcome rule
  does not change that, and it still limits replicate-level generalisation.

## What is still NOT proven

`GSE279162` has **no published count in the five provided scripts that can serve as an independent
reproduction check** — `preprocess_GEX.Rmd` writes `cellsPreandPostFilt.xlsx`, and that output is
not in the archive. The pipeline executes cleanly end to end, so this is a **missing validation
anchor, not a reproduction failure**; it is carried as `UNKNOWN_REQUIRES_SOURCE_FILE` in the results
file so nobody later mistakes "ran without error" for "checked against the authors' numbers". Rewind
has such an anchor and passes it exactly.

Three deviations from the authors' object are declared rather than implied: `min.cells = 1` is
subsumed by the merged `rowSums > 0` filter the author script applies anyway; normalisation, PCA,
UMAP and clustering are not run because they feed visualisation and never `assigned_lineage`; and
the sample merge order is fixed alphabetically, which can only matter for `which.max` ties, and a
tie yields a zero top-two difference that the 0.2 rule sends to `NA` regardless.

Learnability is still untested. No model was fitted.

## Verification

- **1684 passed**, ruff clean (CI scope) — the 21D test file grew from 27 to 41 tests
- `src/` unchanged; 21A / 21B / 21C results and records untouched; Stage 21C stays frozen
- the Rewind headline (42 cells / 35 clones, cutoff 2,365, 101 selected barcodes) was re-derived
  independently in the shell with `sort`/`comm`/`awk`, and the `Naive1` post-QC count (1,595 of
  2,356) was re-derived independently with `awk`; both matched
- `results/diag_stage21d_public_reconstruction_results.json` carries the revision-1 state under
  `supersedes` and a `previous_verdict` per dataset — the old state is preserved, not deleted

## Next action

Stage 22 — Prospective Benchmark Build may open. **Not started in this change.**
