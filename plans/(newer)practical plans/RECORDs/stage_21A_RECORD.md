# stage_21A_ RECORD

## Goal
Audit the geometry of every dataset already on disk and decide whether any of them can pose a
prospective `X_before + U -> Y_future` task. Classify each into the strongest defensible task
level. Fit nothing.

## Inputs
- datasets: `GSE242423`, `GSE165176` (mirrored at `D:\Gill`), `GSE165177`, `GSE165178`,
  `GSE165179`, `GSE113957`, `GSE297234` — local metadata and filenames only, no expression values
- frozen commits: plan `8e7f8ff` (pre-registration), archive `8ed6288`; executed at `f6a0056`,
  CI follow-up `392b69d`
- plan: `(newer)practical plans/arcive/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md`
  §§2, 5, 6, 11, 19, 20 — path updated 2026-08-21 when V2 was archived and superseded by
  `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md`. Section numbers here are V2's; the
  archived file is byte-unchanged, so they still resolve
- model: `_s16` frozen, not loaded

## Files added
- `experiments/diag_stage21_data_audit.py`
- `tests/test_diag_stage21_data_audit.py` (23 tests)
- `results/diag_stage21_data_audit_results.json`

## Files modified
- `tests/test_diag_stage12_rebuild_verdict.py` — read a plan by fixed path; broke when `plans/`
  was reorganised. Now locates by filename anywhere under `plans/`
- `tests/test_diag_stage16_safety_floor.py` — same fix
- `tests/test_stage16_recalibration_verified.py` (`392b69d`) — see **Bugs found** #5

## What changed
- Tri-state audit: `PRESENT` / `ABSENT_PROVEN` / `UNKNOWN_REQUIRES_SOURCE_AUDIT`
- Every field is a `Finding(value, status, evidence)` naming the file and what was read
- A level is ruled out only when evidence proves it out; an `UNKNOWN` suffixes the verdict
  `_PENDING_SOURCE_AUDIT` — a request for one more download, not a rejection
- Dataset location accepts aliases (`GSE165176` lives at `D:\Gill`)

## What did NOT change
- `src/` unchanged (verified: `git diff --name-only src/` empty)
- model unchanged — `_s16` never loaded, nothing fitted
- labels unchanged
- no expression values read

## Tests
- 1610 passed
- ruff clean (CI scope: `src/ tests/ scripts/ plan_tests/`)

## Result
**VERDICT: `CULTURE_FORWARD_AVAILABLE`**

| dataset | level | independent units |
|---|---|---|
| GSE242423 | `TRAJECTORY_FORWARD_PENDING_SOURCE_AUDIT` | 1 |
| **GSE165176** | **`CULTURE_FORWARD`** | **12** ⚠ *revised to 6 by 21B* |
| GSE165177 | `TRAJECTORY_FORWARD` | 6 |
| GSE165178 | `CULTURE_FORWARD` | 4 |
| GSE165179 | `TRAJECTORY_FORWARD` | 6 |
| GSE113957 | `INVALID_PROSPECTIVE_PENDING_SOURCE_AUDIT` | — |
| GSE297234 | `INVALID_PROSPECTIVE_PENDING_SOURCE_AUDIT` | 1 timepoint |

Key finding: `GSE165176` titles carry `SSEA4` / `CD13` — **antibody surface-marker sorts**, an
identity readout measured independently of the RNA vector. The project had never used it; `y_cls`
has always come from `fate_labels()` scoring marker programs on the cell's own expression.

`GSE242423` keeps LEVEL 3 open: its barcodes are provably plain 10x (`ABSENT_PROVEN` for a tag in
the barcode *file*), but no series matrix was on disk, so whether the *study* used lineage tracing
was `UNKNOWN`.

## Bugs found
1. Day tokens appear as `d11` **and** `13days`; missing the second made `GSE165177` look like it
   had no timecourse, driving a false `INVALID`
2. `GM00731_D0` + `GM23815_D0` is two donors at **one** timepoint, not two timepoints — counting
   raw labels turned a D0-only corpus into a fake timecourse
3. An `UNKNOWN n_timepoints` silently drove `INVALID` instead of `PENDING` — a parsing failure was
   rejecting a dataset
4. Blanket-asserting "no independent outcome" hid the SSEA4/CD13 sorts entirely
5. (`392b69d`) Two tests in `test_stage16_recalibration_verified.py` read gitignored
   `cellfate_loocv_*/bundle/` while their guard only checked committed JSON — **this was the red X
   on CI**. Reproduced in a clean clone, fixed with a `needs_folds` marker on just those two

## Scientific interpretation
**Proves:** local corpora were classified from evidence, with every verdict traceable to a named
file. One dataset (`GSE165176`) carries a non-RNA identity readout.

**Does NOT prove:** that any prospective task exists. `CULTURE_FORWARD` here means the audit's
criteria were met on file-level evidence — the audit explicitly recorded
`outcome_is_contemporaneous = True`, i.e. the sort happens at collection and is not yet a *future*
outcome. Nothing about `GSE242423`'s study design was settled.

## Next action
Stage 21B — resolve the two open questions: `GSE242423`'s lineage status from the newly downloaded
source metadata, and whether `GSE165176`'s sorting can define a genuine future culture-level
outcome.
