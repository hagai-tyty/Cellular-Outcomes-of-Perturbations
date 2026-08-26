# stage_24_RECORD — Gen-1 Role-B predictor engineering

## Goal
Execute Stage 24 under the frozen Stage-23.5 contract: reproduce the frozen W5 result, package it,
and hand Stage 25 a set of out-of-fold predictions it can rank. This is a bounded engineering
stage, not an architecture search.

## Authority
`STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, recorded in
`results/stage23_5_protocol.json`. Opening status `STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1`.

## Progress

```text
  24A  consume handoff, freeze engineering plan     DONE
  24B  reproduce W1/W4/W5 under the §7.1 gate       DONE   BYTE_IDENTICAL
  24C  serialization, preprocessing, prediction API  DONE   SERIALIZED_AND_EQUIVALENT
  24D  frozen OOF per clone-condition row            --
  24E  deterministic scoring + leakage contracts     --
  24F  freeze W5 tool artifacts                      --
  24G  hand frozen predictions/model to Stage 25     --
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

## Scientific interpretation

**Proves:** the frozen Stage-23 W5 result is exactly reproducible from the committed code and data —
byte-for-byte, including every out-of-fold prediction and every within-clone condition ordering the
Stage-25 ranking test will consume. Stage 25's inputs are therefore the frozen inputs, not an
approximation of them.

**Does NOT prove:** anything about the ranking claim. No Stage-25 statistic exists, and Stage 24 is
forbidden from inspecting the ranking metric. Reproducibility is a property of the pipeline, not
evidence for the hypothesis.

## Next action
24C — serialization, preprocessing and the prediction API. The §8.7 permutation run (~19-20 h across
three shards) belongs to Stage 25 and has **not** been started.
