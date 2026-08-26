# stage_24_RECORD — Gen-1 Role-B predictor engineering

## Goal
Execute Stage 24 under the frozen Stage-23.5 contract: reproduce the frozen W5 result, package it,
and hand Stage 25 a set of out-of-fold predictions it can rank. This is a bounded engineering
stage, not an architecture search.

## Authority
`STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, recorded in
`results/stage23_5_protocol.json`. Opening status `STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1`.

## Inputs
- `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md` @ `8da16fca...bced48`, plus
  `results/stage23_5_protocol.json` and `results/stage23_5_handoff_to_stage24.json`
- the frozen Stage-23 WM989 artifacts, read-only: `stage23_wm989_detection_oof.csv`,
  `stage23_wm989_abundance_oof.csv`, `stage23_wm989_interaction_oof.csv`,
  `stage23_wm989_interaction_abundance_oof.csv`, `stage23_wm989_results.json`,
  `stage23_wm989_interaction_results.json`
- `results/stage22_wm989_clones.csv` for the nuisance block
- `_cc_cache/stage23/GSE279162_pseudobulk.npz`, the frozen clone-level representation
- `D:/GSE279162` for the feature-schema contract

## Files added
- `experiments/run_stage24_gen1_tool.py`
- `src/cellfate/gen1_predictor.py` — the shipped API
- `tests/test_stage24_gen1_tool.py`, `tests/test_gen1_predictor.py`
- `results/stage24/` — engineering plan, reproduction report, serialization report, artifact
  metadata, and `repro/` holding the regenerated files the gate compared against

## Files modified
- `tests/test_stage23_learnability_gate.py` — generalised the interstitial exclusion
- `.gitignore` — the 44 MB tool artifact
- no Stage-22, Stage-23, Stage-23.2 or Stage-23.5 artifact was modified

## What changed
- Stage 24 opened and 24A-24C executed
- W5 serialized into a loadable artifact with a public prediction API

## What did NOT change
- `src/cellfate/` gained a new module but no existing module was altered
- the W5 model, its features, its grid and its selection rule are the frozen Stage-23 ones
- every frozen Stage-23 artifact still hash-matches `git HEAD`, asserted by contract
- no Stage-25 statistic exists; the ranking metric has not been inspected

## Progress

```text
  24A  consume handoff, freeze engineering plan     DONE
  24B  reproduce W1/W4/W5 under the §7.1 gate       DONE   BYTE_IDENTICAL
  24C  serialization, preprocessing, prediction API  DONE   SERIALIZED_AND_EQUIVALENT
  24D  frozen OOF per clone-condition row           DONE   8,406 rows, 892 eligible
  24E  deterministic scoring + leakage contracts     DONE   7/7 checks
  24F  freeze W5 tool artifacts                      DONE   all §6.5 deliverables
  24G  hand frozen predictions/model to Stage 25     DONE   STAGE_24_GEN1_TOOL_READY
```

---

## 24A — the freeze it consumed

A handoff that does not match the plan it claims to come from is not a handoff. All 12 checks pass:

```text
  plan digest matches protocol            plan digest matches handoff
  plan status FROZEN                      Stage 24 open
  ranking metric NOT inspected            ranking protocol hash frozen
  no new datasets authorized              Stage 27 not a Gen-1 gate
  audit fully passed (38/38)              compute budget accepted
  no source-artifact drift                no ranking artifact exists on disk
```

The drift check re-hashes every source artifact the protocol pinned. The ranking check globs
`results/` for any file whose name contains "rank" — Stage 24 generates Stage 25's inputs and must
not be able to see the answer.

## 24B — reproduction

```text
  C1_W0toW4   BYTE_IDENTICAL   R1 pass  R2 pass  R3 pass
  C2_W0toW4   BYTE_IDENTICAL   R1 pass  R2 pass  R3 pass
  C1_W5       BYTE_IDENTICAL   R1 pass  R2 pass  R3 pass

  23C reproduced in 1.29 min   23D in 1.00 min
  R2 worst absolute difference across every prediction cell   0.0
  R3 within-clone orderings verified                          1,401 clones, W4 and W5
```

The §7.1 fallback never engaged. Byte-identity satisfied the primary requirement outright, so no
tolerance argument was needed and R4 has nothing to name.

### How the frozen artifacts stayed frozen

`run_23c` and `run_23d` write to module-level path constants — including
`results/stage23_wm989_interaction_oof.csv`, the very file 24B must compare against. Running them
unmodified would have overwritten the comparison target with the comparison.

24B rebinds those five constants to `results/stage24/repro/` for the duration of the call. The
frozen files become the target and are never written to; the models are the frozen implementation
**called**, not re-typed (plan §5.1). All six frozen artifacts still hash-match `git HEAD`, and that
is asserted by contract rather than claimed.

## Bug found — in my own gate, caught by the artifact contradicting itself

The first 24B run reported `C2_W0toW4` as `BYTE_IDENTICAL` **and** `R1=False`. Those cannot both be
true, and the contradiction was the tell.

```text
  cause    EXPECTED_ROWS was a single constant, 8406.
           C1 (detection) scores every clone x condition row      8,406
           C2 (abundance) is defined only where a clone was
              detected, so it carries the nonzero rows            2,256
           R1 also short-circuited the key check behind the row-count check,
           so a count mismatch masked the key result instead of reporting it.

  impact   none here -- byte-identity carried the verdict. But a reproduction that
           was merely tolerance-clean rather than identical would have hit a
           SPURIOUS INPUT_INTEGRITY_STOP on C2 and blocked Stage 24 for no reason.

  fixed    per-endpoint row counts; the three R1 sub-checks computed independently;
           and a gate_self_consistent assertion that RAISES if a byte-identical file
           ever fails a sub-gate again.
```

The last fix is the one that matters: it converts this class of defect from something spotted by eye
into something the code refuses to proceed past.

## 24C — serialization, preprocessing and the prediction API

The frozen builder computes W5 but never returns a fitted object: `expression_block` returns
transformed arrays and `_fit_logistic` returns predictions. A tool needs the learned state itself,
so 24C rebuilds it from the builder's own helpers -- `_frozen_pipeline_cache`,
`standardize_train_only`, `treatment_dummies`, `interaction_block` -- and then proves the rebuild
is faithful the only way that means anything.

```text
  ARTIFACT vs the frozen pred_W5 column
    rows compared                8,406
    tolerance                    1e-12
    max absolute difference      4.996e-16
    rows over tolerance          0
    verdict                      SERIALIZED_AND_EQUIVALENT

  SHIPPED API vs the frozen pred_W5 column   (40 clones, via src/cellfate/gen1_predictor.py)
    max absolute difference      2.220e-16
```

Nothing was taken on the resemblance of the code. If the artifact could not regenerate the frozen
column the module raises and 24C fails.

### Two components are shipped, and they answer different questions

```text
  fold0..fold4   one per outer fold, each with its own gene filter, PCA, scalers, (K, C) and
                 coefficients. These reproduce the frozen OOF and are what Stage 25 consumes.
                 A benchmark clone is scored by the single fold model that did not train on it.
  deployment     the same W5 specification and the same inner-CV selection rule applied once to
                 ALL 1,401 clones, for scoring a NEW clone that belongs to no fold.
                 K=50, C=0.1, 18,290 genes retained, 309 design columns.
```

The deployment component is **packaging, not a new model** -- same features, same grid, same
selection rule. It is recorded as NOT validated on held-out data: its performance is *estimated* by
the frozen out-of-fold result, which came from the fold components. The plan does not specify which
model a tool should ship, so this is an engineering decision taken inside Stage 24's remit and
recorded rather than slid past. The §8 ranking contract is untouched.

### The API refuses three things, and the refusals are the point

`src/cellfate/gen1_predictor.py` implements the §6.2-§6.4 contract:

```text
  MISSING_REQUIRED_NUISANCE     B is part of the evaluated model. Expression alone is not
                                equivalent to it, so a missing, short or non-finite nuisance
                                vector returns no score. It is never imputed to a default.
  UNSUPPORTED_TREATMENT         unknown conditions are never embedded, nearest-neighboured or
                                mapped onto a known condition.
  validated_condition_order     withheld unless a Stage-25 verdict file records
                                RANKING_SUPPORTED. The six scores are always returned; their
                                ORDER carries no validated meaning until then.
```

A tool that quietly imputed `B`, or mapped an unknown drug onto a known one, would report
frozen-benchmark numbers for something the benchmark never evaluated.

The 44 MB artifact is gitignored -- five fold components plus deployment, each carrying a PCA basis
over ~18k retained genes. It rebuilds in ~0.5 min from committed inputs, its sha256 is recorded in
the committed `stage24_w5_artifact.json`, and for the reproducibility package it ships as a release
asset.

## 24D — the table Stage 25 consumes

```text
  rows                    8,406      one per clone x condition
  clones                  1,401      six rows each, one outer fold each
  columns                 clone_id, treatment, outer_fold, y, pred_W1, pred_W4, pred_W5,
                          detected_post, ranking_eligible
  ranking-eligible        892        >=1 C1-positive AND >=1 C1-zero  (plan §8.4)
  excluded                472 all-zero, 37 all-positive -- within-clone AUROC undefined
  integrity checks        8 / 8
```

`y` was asserted equal to `detected_post` rather than assumed. The 892 count is verified here
because §8.4 requires it verified **before** scoring, and 24D is the last point at which Stage 24
touches the table.

**No ranking statistic was computed.** Stage 24 is forbidden from inspecting the ranking metric, so
24D emits inputs and asserts their integrity — no AUC, no `delta_RANK`, no top-1 quantity. A
contract scans the artifact for those names to keep it that way.

## 24E — determinism and leakage

```text
  deterministic within a session                 True    60 clones
  deterministic across independent loads         True
  reproduces the frozen out-of-fold column       max |diff| 2.498e-16

  fold isolation      every fold component verified DISJOINT from the clones it scores,
                      read from the artifact's own recorded training set rather than
                      inferred from the fact that the OOF reproduces
  deployment          trained on all 1,401 clones, as declared
  no outcome array    the artifact carries no array of the outcome's length
  feature space       36,601 Gene Expression features; the 153,055 WM989 Custom lineage
                      features cannot be addressed by a filter indexing into a GE space
```

That fold-isolation check is why 24C started recording each component's training clone set: proving
a model never saw what it scores is worth more than inferring it.

## 24F — the §6.5 deliverables, frozen and hashed

```text
  Python prediction API        src/cellfate/gen1_predictor.py
  command-line interface       src/cellfate/gen1_cli.py
  frozen model artifact        stage24_w5_artifact.npz  (gitignored, 44 MB, rebuilds in ~0.5 min)
  frozen vocabularies          treatments, nuisance order, feature contract
  preprocessing artifact       gene filter + PCA basis + scalers, per component
  machine-readable schemas     tool/io_schema.json
  model card                   tool/MODEL_CARD.md
  example dataset              tool/example_clones.csv, from permitted benchmark material only
  unit tests                   36 contracts across two modules
  end-to-end reproduction      24B BYTE_IDENTICAL, 24C equivalence 5e-16
```

The CLI distinguishes outcomes in its **exit code**: `0` every condition scored, `2` at least one
refused, `3` input unreadable. A refusal that exits `0` would let a caller treat a missing score as
a real one.

## 24G — handoff

```text
  STAGE_24_GEN1_TOOL_READY

  table            results/stage24/stage24_oof_for_stage25.csv, hashed in the handoff
  model            stage24_w5_artifact.npz + its metadata
  population       892 eligible clones, verified mechanically
  ranking metric   NOT inspected by Stage 24
  ranking stat     NOT computed by Stage 24
```

The handoff also names what Stage 25 may not do: change the metric, population, weighting,
comparator, endpoint or null; reduce the permutation count or stop early; use C2 or a per-treatment
result to rescue a failed C1 ranking; add a dataset; or revise the plan after seeing a result.

## Tests
- 12 Stage-24 contracts + 13 predictor contracts, 0 skipped
- **Mutation-tested**, all five caught: flipping `gate_self_consistent`, changing C2's row
  expectation to 8406, cutting R3's clone count to 900, injecting 3 over-tolerance cells, and
  marking one within-clone ordering as changed
- `pytest` 2108 passed, 1 skipped · ruff clean (CI scope)

## Bugs found outside this stage
Two existing CI invariants caught real omissions of mine, both fixed:

1. `test_the_determinism_set_covers_every_committed_stage23_artifact` excluded interstitials with a
   hardcoded `startswith("stage23_2")`, which stopped covering the moment Stage 23.5 wrote its first
   artifact. Generalised to `stage23_<digit>`, plus an assertion that the exclusion has not
   swallowed a real Stage-23 artifact.
2. `test_results_paths` requires every writer to define `_RESULTS` in a literal `__file__`-relative
   form. The new module did not. Conformed.

Neither was a CI defect. Both were the checks doing exactly what they exist for.

## Result

```text
  24A   handoff integrity                 12 / 12 checks
  24B   reproduction                      BYTE_IDENTICAL on all three files
          R2 worst absolute difference     0.0 over every prediction cell
          R3 within-clone orderings        1,401 clones, W4 and W5, unchanged
  24C   serialization                     SERIALIZED_AND_EQUIVALENT
          artifact vs frozen pred_W5       max |diff| 4.996e-16 over 8,406 rows
          shipped API vs frozen            max |diff| 2.220e-16 over 40 clones
          deployment component             K=50, C=0.1, 18,290 genes, 309 columns
  24D   handoff table                     8,406 rows, 892 eligible, 8/8 checks
  24E   determinism + leakage             7/7 checks, every fold component isolated
  24F   §6.5 deliverables                 all present and hashed
  24G   handoff                           STAGE_24_GEN1_TOOL_READY
```

## Scientific interpretation

**Proves:** the frozen Stage-23 W5 result is exactly reproducible from the committed code and data —
byte-for-byte, including every out-of-fold prediction and every within-clone condition ordering the
Stage-25 ranking test will consume. Stage 25's inputs are therefore the frozen inputs, not an
approximation of them.

**Does NOT prove:** anything about the ranking claim. No Stage-25 statistic exists, and Stage 24 is
forbidden from inspecting the ranking metric. Reproducibility is a property of the pipeline, not
evidence for the hypothesis.

## Next action
**Stage 25.** Compute `delta_RANK` from the frozen table only, verify the 892-clone population
before scoring, and run the 1,000-draw full-refit null (~19-20 h across three shards, per-shard
cache files, completeness assertion, no early stopping). Record `STAGE_25_RANKING_SUPPORTED` or
`STAGE_25_RANKING_NOT_SUPPORTED` once; either verdict proceeds to `GEN1_MANDATORY_SHIP`.

That permutation run has **not** been started.
