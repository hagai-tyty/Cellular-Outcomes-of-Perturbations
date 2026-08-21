# stage_23_ RECORD

Authoritative record for Stage 23. Sections are **additive and immutable**: 23A is written now,
23B–23F are appended as they execute. An earlier section is never rewritten because a later result
is inconvenient; corrections are appended with a dated note.

---

# 23A — Protocol / representation freeze

## Goal
Run the pre-registered §3.0 input audit against the frozen Stage-22 benchmark and the external raw
matrices, build one fixed clone-level `X_before`, and freeze every evaluation choice before any
outcome model exists. **No model is fitted in 23A.**

## Inputs
- plan: `STAGE_23_LEARNABILITY_INTERACTION_GATE_V2.md` (V1 archived under `arcive/` after the
  independent pre-execution audit)
- frozen Stage-22 benchmark commit `8d6011a`; record `6b6a169`; correction `9eaddf2`
- `D:\GSE227151_Rewind\` (11 source + 2 author-code files) · `D:\GSE279162\` (29 + 5)
- model: `_s16` frozen, never loaded

## Files added
- `experiments/run_stage23_learnability_gate.py`
- `tests/test_stage23_learnability_gate.py` (23 tests)
- `results/stage23_protocol.json`
- `results/stage23_rewind_clone_expression_manifest.json`
- `results/stage23_wm989_clone_expression_manifest.json`
- `results/stage23_outer_fold_preprocessing.json`

## Files modified
- `tests/test_ci_portability.py` — the new Stage-23 test module is registered with the
  CI-condition guard so it is re-run with every local dataset root absent

## What changed
- The §3.0 audit is **executable code, not prose**: 38 checks, and pseudobulk construction is
  gated behind `STAGE22_INPUTS_AUDITED`
- `X_before` exists as a frozen artefact with a content hash; the matrix itself stays in the
  gitignored `_cc_cache/` and never enters git
- Every later-substage choice — K, grids, seeds, treatment coding, tie-break, metrics, bootstrap
  and permutation rules — is written into `stage23_protocol.json` now

## What did NOT change
- `src/` unchanged · no model fitted · no Stage-22 artifact, fold, label or manifest rewritten
- the frozen Stage-22 CRLF/LF provenance mismatch is **declared, not repaired** (V2 §1.4)
- the frozen WM989 target is **not** rebuilt from a treated-only alternative (V2 §1.2.1)

## Tests
- 1857 passed · ruff clean (CI scope)
- 121 passed under `CELLFATE_NO_LOCAL_DATA=1` (the CI condition) for the Stage-23 + portability
  modules; the pseudobulk cache is absent there and those tests skip rather than fail

## Result

**VERDICT: `PROTOCOL_FROZEN`** — input audit `STAGE22_INPUTS_AUDITED` (38/38 checks).

### Input audit

```text
six benchmark CSV hashes match their manifests          OK
two manifest hashes match stage22 results               OK
results file omits its own hash                         OK
Role-A verdict ready                                    BENCHMARK_READY_WITH_DECLARED_MISSINGNESS
all_gates_pass / every individual gate / model_fitted    true / 10-of-10 true / false
preflight derived independently of the `overall` string  OK
40 external source + 7 author-code files                present, byte-identical
every expression_column_index resolves                  0 mismatches over 10,394 cells
36,601 Gene Expression features per sample              OK   (Custom block 153,055 where present)
all 11 samples share ONE gene feature-ID list           OK
WM989 targets rebuilt from treated cells only           8,406 rows, max abs diff 0
canonical text hash is LF/CRLF invariant                OK
```

Two inherited limitations are recorded rather than silently fixed:

- **Stage-22 gate derivation.** Stage 22 derives `overall` from the Role-A verdict without
  consuming `all_gates_pass`. Confirmed harmless here because all ten gates are true — and Stage 23
  fails closed on its own derived condition regardless of the serialized string.
- **CRLF/LF provenance.** The Stage-22 plan and builder digests are checkout-byte hashes and differ
  between Windows and Linux for identical content. Stage 23 uses a canonical LF-normalised digest
  for text; a test asserts LF, CRLF and CR inputs all produce one hash.

### Clone `X_before`

```text
                 clones   genes    nnz          density   cells/clone (min/med/max)
GSE227151         3,147   36,601   15,880,583   0.1379    1 / 1 / 11
GSE279162         1,401   36,601   10,697,195   0.2086    1 / 2 / 88
```

Built by summing **raw** pretreatment counts per clone, then CP10K and `log1p` **exactly once**.
Verified by inverting the transform: `expm1` row sums are 10,000.000 for every clone in both
matrices, and one 11-cell Rewind clone was recomputed from the raw `.mtx` independently, agreeing
to `4.4e-16`.

Only pretreatment lanes were read — the two Rewind `*_control_*` matrices and the three WM989
`*_Naive*` matrices. No treated column was opened, and `X` is exactly 36,601 columns wide, so none
of WM989's 153,055 `Custom` lineage features can be present.

### Outer-training gene filter (descriptive)

```text
              fold 0   fold 1   fold 2   fold 3   fold 4    detection floor
GSE227151     13,610   13,589   13,627   13,606   13,607    26
GSE279162     18,385   18,411   18,270   18,227   18,361    12
```

`max_feasible_K >= 50` in every fold, so **no candidate K may be skipped**. This table is
descriptive only: per V2 §2.5 the filter is refitted inside every inner-training split.

## Bugs found
1. **The audit surfaced nothing new** — every V2 §1 anchor reproduced. That is the expected
   outcome given the pre-execution audit, and is recorded as a pass rather than as an absence
2. A first draft of the module wrote the pseudobulk into `results/`. Corrected before commit: a
   3,147 × 36,601 matrix is not a compact benchmark artefact, and V2 §13 forbids it. It now lives
   in the gitignored `_cc_cache/`, with its content hash committed so a rebuild stays verifiable

## Scientific interpretation
**Proves:** the Stage-22 benchmark is intact and externally consistent with the raw data at the
moment Stage 23 opened, and one fixed clone-level representation now exists whose construction is
independently checkable from the committed hashes. Every evaluation choice that could otherwise be
made after seeing a result is frozen.

**Does NOT prove:** anything about learnability. No model has been fitted, no outcome has been
exposed to an estimator, and no comparison has been run. The observed-zero/abundance confound
recorded in Stage 22 is untouched by this substage and remains the load-bearing risk for Role B —
which is why V2 §1.2.1 added `log1p(n_naive_cells)` to the nuisance block before any fitting.

## Next action
23B — Rewind Role-A learnability. Not started.
