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
