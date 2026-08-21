# stage_22_ RECORD

## Goal
Convert the two Stage-21D-reconstructed datasets into frozen, auditable prospective benchmark
tables — `X_before + U -> Y_future` with clone-grounded linkage, deterministic clone-level outer
folds, an expression-column mapping, and a feature-eligibility firewall — and decide whether Stage
23 may open. No model is fitted.

## Inputs
- `D:\GSE227151_Rewind\` (11 required paths) and `D:\GSE279162\` (29 required paths), both already
  local; no download performed
- author code: `author_code_zenodo7707418\` (2 R1 scripts) and `author_code_Schaff_manuscript\`
  (5 primary files) — referenced by path + SHA-256, never committed
- frozen commits: 21D revision 2 `6c2f2d6`; plan `e97efe1`, annotated `f2d6c86`
- plan: `STAGE_22_PROSPECTIVE_BENCHMARK_CONSTRUCTION_V2.md`
  (`b26db9fe…` at build time — recorded in both manifests)
- model: `_s16` frozen, never loaded

## Files added
- `experiments/build_stage22_prospective_benchmarks.py`
- `tests/test_stage22_prospective_benchmarks.py` (44 tests)
- `results/stage22_prospective_benchmark_results.json`
- `results/stage22_rewind_benchmark_manifest.json`, `results/stage22_wm989_benchmark_manifest.json`
- `results/stage22_rewind_cells.csv` (3,905 rows) · `results/stage22_rewind_clones.csv` (3,147)
- `results/stage22_wm989_cell_assignments.csv` (46,891) ·
  `results/stage22_wm989_naive_cells.csv` (6,489) ·
  `results/stage22_wm989_clones.csv` (1,401) ·
  `results/stage22_wm989_clone_treatment.csv` (8,406)

## Files modified
- none — **additive only**. 21A / 21B / 21C / 21D results and records untouched (pinned by test)

## What changed
- Both benchmarks are built by **importing** the Stage-21D author-rule implementation, not by
  re-implementing it. A test asserts the builder defines none of `slice_max_with_ties`,
  `qc_and_lineage`, `barcode_clustering`, `barcode_combine`,
  `barcoding_posterior_and_assignment` itself (plan §1.3)
- Clone-level outer folds are **frozen here**, once, over sorted clone keys — not left to Stage 23
- Dataset roots are CLI arguments; no `D:\` literal exists in the builder (plan §1.4)
- Manifests carry no timestamp and no repository `HEAD`; the whole artifact set was verified
  **byte-identical across two consecutive runs**

## What did NOT change
- `src/` unchanged (verified: `git diff --name-only src/` empty)
- model unchanged — nothing fitted; `sklearn` is imported for `KFold`/`StratifiedKFold` only, and
  the no-modelling gate proves it from the syntax tree
- labels unchanged — the author rules are byte-for-byte the Stage-21D ones
- no threshold tuned, no binary resistance endpoint frozen

## Tests
- 1732 passed
- ruff clean (CI scope: `src/ tests/ scripts/ plan_tests/`)
- expression mapping (3,905 Rewind + 6,489 WM989 cells), fold consistency across all five tables,
  and the Rewind label were each re-derived independently outside the builder and matched

## Result

**OVERALL: `STAGE_23_READY`** — all ten gates pass, derived from the two per-dataset verdicts.

### A. `GSE227151` (Rewind, Role A) → `BENCHMARK_READY_WITH_DECLARED_MISSINGNESS`

```text
source assignment records            3,921
source unique cell_uid               3,913
source unique clones                 3,149
bare cellID cross-lane collisions        0

ambiguous cell_uid excluded              8   (16 source rows, reason
                                              ambiguous_multi_lineage_clone_assignment)
retained benchmark cell_uid          3,905
retained benchmark clones            3,147
clones lost entirely to exclusion        2

primed cells / clones               42 / 35
nonprimed cells / clones         3,863 / 3,112
prevalence            cell 1.08%   clone 1.11%
clones spanning both lanes             306
expression mapped                3,905 / 3,905
```

Folds (`StratifiedKFold`, 5 splits, seed 22022, on the post-exclusion clone table):

```text
fold   clones   positive clones   positive cells
  0      630           7                9
  1      630           7                8
  2      629           7                8
  3      629           7                8
  4      629           7                9
```

- `y_primed` is the author top-100-with-ties gDNA cut: cutoff `nUMI = 2,365`, 2 barcodes tied at
  it, **101** barcodes selected
- `cell_uid = SampleNum:cellID`; `SampleNum 1 -> GSM7092516`, `2 -> GSM7092515`, re-derived by
  containment, never from the GEO title text
- cells per clone: 2,584 clones have 1 cell, 416 have 2, then a tail to 11

### B. `GSE279162` (WM989, Role B) → `BENCHMARK_READY_WITH_DECLARED_MISSINGNESS`

```text
raw cells                           77,417
post-QC cells                       46,891      (author per-condition RNA QC)
assigned cells                      42,771
NA-lineage cells                     4,120      documented, never reassigned
unique cell_uid                     46,891
reused bare barcodes  raw 722  ·  post-QC 248

unique assigned clones               2,215
clones with a naive observation      1,401
clone x treatment rows               8,406      = 1,401 x 6
observed_zero rows                   6,150      73.2%
naive expression mapped          6,489 / 6,489

coverage   >=1 treatment 929  ·  >=2 treatments 603  ·  all 6 treatments 37
naive cells per clone   median 2, max 88; 816 clones have >=2
folds (KFold, seed 22022)   281 / 280 / 280 / 280 / 280 clones
```

Per treatment, over the 1,401 eligible clones:

```text
             assigned cells   clones nonzero   zero rate   benchmark fraction sum   max abundance
Dabrafenib        6,164            321          77.1%             0.716                 476
Trametinib        8,091            274          80.4%             0.553                 776
CoCl2             5,134            319          77.2%             0.843               3,290
Acid              7,246            665          52.5%             0.942                 464
Cisplatin         6,815            272          80.6%             0.986               2,791
Doxorubicin       2,832            405          71.1%             0.806                 577
```

`Y(clone, treatment) = post-treatment assigned-cell abundance`, plus `post_fraction`,
`post_rank` (min competition rank over the eligible clones), `post_rank_fraction` and
`post_tie_size`. **No binary resistance threshold was frozen**, and the script's figure-specific
`num_lin = 5` was not turned into a target — a test asserts no such binding exists in our code.

### The abundance confound, measured

```text
naive cells per clone      clones     observed-zero rate
        1                    585            87.1%
        2                    242            82.6%
      3-4                    218            67.9%
      5-9                    189            59.9%
      10+                    167            32.6%
```

Captured clone size predicts the outcome hard — a 2.7× swing. This is recorded as a declared
limitation, not corrected, and it is why the Stage-23 abundance-only nuisance baselines are
mandatory rather than a formality.

## Bugs found
1. **Two gates were hard-coded `True` in the first draft** (`G22-2` no-label-leakage, `G22-10`
   no-modelling). Both are now *computed*: the Rewind target is recomputed from the gDNA arm alone
   without opening an expression file and compared with what was written, and the modelling check
   walks the builder's own syntax tree. Asserting a gate is exactly the failure mode Stage 13
   caught in the scorecard
2. **The no-modelling gate failed on itself.** Its first version searched the source text for the
   literal `.fit(`, which matched its own docstring. Rewritten to inspect `ast.Call` nodes instead.
   It only surfaced *because* the gate was computed rather than asserted
3. **`722` reused barcodes is the RAW figure, not the benchmark's.** Over the post-QC population it
   is **248**. Both are now recorded with a note on which population each measures, rather than
   picking whichever matched the plan's prose
4. Missed the repo convention that every results-writer defines `_RESULTS` — caught by
   `tests/test_results_paths.py` before commit, as in Stage 21D

## Scientific interpretation
**Proves:** both prospective tasks are now frozen artifacts rather than descriptions. Every
benchmark row traces pretreatment cell → clone → independently measured later outcome; every clone
carries exactly one `outer_group` and one `outer_fold` across every table it appears in; every
eligible pretreatment cell resolves to a specific expression column; and the whole artifact set
reproduces byte-for-byte. Stage 23 can no longer redefine the task after seeing performance.

The pre-registered ambiguity exclusion was the right call and the plan's insistence on
*recomputing* post-exclusion counts was vindicated: 2 clones consisted only of excluded cells, so
the retained clone count is **3,147, not 3,149**. Had the plan frozen the Stage-21D figure as a
retained-count assertion it would have failed on the first run.

**Does NOT prove:** anything about learnability. Not that `X` predicts `Y`, not that `X` adds
beyond `U`, not that an `X × U` interaction exists, not that the outcome is free of sampling or
censoring bias. Three limitations are frozen into the manifests rather than left to be discovered:
Rewind is **one** biological replicate, so clone-held-out evaluation is within-experiment
generalisation and `nonprimed` is an operational complement rather than proven failure; WM989's
`post_fraction` shares a treatment-level denominator, so the 8,406 rows are **not** 8,406
independent units; and the observed-zero mass is a measurement/censoring question, not a biological
label.

A legitimate Stage-23 outcome is that transcriptomic `X` adds little or nothing beyond captured
clone size. The benchmark is built so that result would be visible rather than hidden.

## Next action
Stage 23 — learnability / interaction gate. **Not started.** Before any model is fitted, Stage 23
must freeze in its own plan: the `X` representation and clone-aggregation rule, expression
normalisation and feature selection (fit on training folds only), clone-multiplicity weighting,
primary metrics (PR-AUC / average precision at clone grain for Rewind, given 35 positive vs
**3,112** negative clones), the baseline clone-abundance nuisance variables, treatment-count
normalisation, and the interpretation of `observed_zero`.
