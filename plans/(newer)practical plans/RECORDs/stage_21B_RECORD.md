# stage_21B_ RECORD

## Goal
Resolve the two questions Stage 21A left open, from source metadata only: does `GSE242423` contain
any real lineage/clone linkage, and can `GSE165176`'s SSEA4/CD13 sorting define a legitimate
**future** culture-level outcome rather than a contemporaneous phenotype. Fit nothing.

## Inputs
- `D:\GSE242423\GSE242423_series_matrix.txt.gz` and `GSE242423_family.xml.tgz` (both newly
  downloaded for this stage) plus the existing barcode/matrix files
- `D:\Gill\GSE165176_series_matrix.txt.gz` — 124 samples
- frozen commits: 21A at `f6a0056`, CI fix `392b69d`; executed at `3b8b644`
- plan: `(newer)practical plans/arcive/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md`
  — path updated 2026-08-21 when V2 was archived and superseded by
  `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md`; the archived file is byte-unchanged
- model: `_s16` frozen, not loaded

## Files added
- `experiments/diag_stage21b_source_design.py`
- `tests/test_diag_stage21b_source_design.py` (25 tests)
- `results/diag_stage21b_source_design_results.json`

## Files modified
- none — **additive only**. `results/diag_stage21_data_audit_results.json` was not touched
  (verified: `git diff` empty on that file; its verdict is still `CULTURE_FORWARD_AVAILABLE`)

## What changed
- Both 21A pending questions are closed on evidence
- 21A's `n_independent_units = 12` for `GSE165176` is **superseded by 6** — verified, not inherited
- Same tri-state rule carried over: `PRESENT` / `ABSENT_PROVEN` / `UNKNOWN_REQUIRES_SOURCE_FILE`

## What did NOT change
- `src/` unchanged (verified)
- model unchanged — nothing fitted, no sklearn/torch imported (asserted by test)
- labels unchanged
- Stage 21A's result file and verdict unchanged

## Tests
- 1639 passed
- ruff clean (CI scope)

## Result

**A. `GSE242423` → `LINEAGE_ABSENT_PROVEN`**
- **0** occurrences of clone / lineage / CellTag / LARRY / hashtag / HTO / multiplex / CMO /
  sister / barcoded across **50,648 chars** of Series Matrix + MINiML
- all **24** distinct `barcode` contexts are 10x matrix vocabulary or `barcodes.tsv.gz` filenames
- per-sample characteristic tags are exactly `["cell type", "genotype"]`
- 20 supplementary files, none a lineage map; no replicate token in any of 9 titles
- extract protocol **trypsinises at collection** → each day is a destructive sample
- `n_independent_units = 1` continuously-cultured trajectory. The ~42,500 cells are **not**
  independent trajectories

**B. `GSE165176` → `ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY`**
- the sort is an **orthogonal phenotype/state call measured independently of RNA**: CD13 → *"Failing to reprogram
  fibroblast"*, SSEA4 → *"Reprogramming fibroblast"* / *"iPSC"*, unsorted → *"Dermal fibroblast"*
- **47 of 71** donor × day × experiment cultures yield **both** fractions → marker identity is a
  within-culture subpopulation label, not a culture outcome
- **118 of 124** samples are already sorted; only **6** unsorted samples exist, all day 0 → a
  sorted early input predicting a later marker is leakage
- **no** FACS proportions/percentages/frequencies → `early RNA -> later %SSEA4+` cannot be built
- **no** colony count, efficiency, survival or terminal assay
- at day 54 all 6 donors appear with **only SSEA4** → no outcome contrast
- Exp1 = days [7,9,11,13,15]; Exp2 = days [0,11,21,29,34,40,47,54]; overlap **day 11 only** →
  time blocks, not replicates → **effective n = 6 donors**

## Bugs found
1. `parse_gill_titles` dropped unparseable titles, so row order no longer aligned with the
   per-sample characteristics arrays — the marker→cell-type mapping could have gone silently
   off-by-N. Fixed by carrying the original index
2. Stage 21A's `n_independent_units = 12` (6 donors × 2 experiments) overstates the unit count.
   21A's figure is left frozen in its own file; a test pins both numbers so the correction is
   auditable rather than a silent overwrite

## Scientific interpretation
**Proves:** neither local dataset can pose a prospective `X_before + U -> Y_future` task.
`GSE242423` is one unlinked, destructively-sampled trajectory. `GSE165176` carries a genuinely
orthogonal antibody **phenotype/state call measured independently of RNA**, but it is a
same-timepoint subpopulation split with no proportions, no unsorted early population beyond day 0,
and no terminal outcome variation.

**Does NOT prove:** that the sorting data are worthless — the SSEA4/CD13 call remains the only
non-RNA phenotype readout in the project and stays available for *contemporaneous* analyses. It also
says nothing about public datasets; only that manufacturing a task from local data would require
inventing geometry the experiments never had.

Nothing here is `UNKNOWN`. No further download changes either answer.

## Next action
Both local routes are closed, so per the plan the next action is **Stage 21C — Public
Prospective Dataset Qualification** — not another attempt on local data. Before searching,
pre-register a fixed search budget and stopping rule (finite candidate set or defined number of
passes), with Role-B scoped down rather than allowed to hold the paper hostage.

*Corrected 2026-08-21 (terminology): every place that called the SSEA4/CD13 sort a "fate call" or
"fate readout" — the Result section and both halves of the Scientific interpretation — now reads "orthogonal
phenotype/state call measured independently of RNA". The sort is orthogonal to RNA, which is
the real finding; calling it a *fate* call overstates it, because 21B's own evidence shows it
is contemporaneous, not a future outcome. The evidence strings in
`results/diag_stage21b_source_design_results.json` are left frozen as executed.*

*Corrected 2026-08-21 (next action): this originally read "Stage 21B public Role-A / Role-B dataset
qualification", written before the 21C/21D split existed. Keeping it would have left two different
things called 21B. Sequence now frozen in
`STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §2 — briefly held in a standalone
`STAGE_21_EXECUTION_AMENDMENT_AFTER_21B.md`, since folded into V3 so there is one plan file, not
a plan plus a patch.*
