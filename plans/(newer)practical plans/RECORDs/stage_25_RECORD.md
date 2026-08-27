# stage_25_RECORD — the preregistered clone-specific ranking test

## Goal
Execute §8 of the frozen Stage-23.5 plan exactly once: does the frozen interaction model W5 use
pretreatment state to improve clone-specific ordering of the six observed WM989 conditions over the
non-interactive additive model W4?

## Inputs
- `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md` @ canonical-LF `8da16fca...bced48`, FROZEN
- `results/stage24_handoff_to_stage25.json` — `STAGE_24_GEN1_TOOL_READY`
- `results/stage24/stage24_oof_for_stage25.csv` — the frozen out-of-fold table, hash-checked
- `results/stage22_wm989_clones.csv` — nuisance block and depth bins
- `_cc_cache/stage23/GSE279162_pseudobulk.npz` — the frozen clone representation, for the null

## Files added
- `experiments/run_stage25_ranking.py`
- `tests/test_stage25_ranking.py`
- `results/stage25/stage25a_observed.json`, `stage25_smoke.json`, `stage25_verdict.json`

## Files modified
- `.gitignore` — the null shard cache

## What changed
- The one load-bearing new capability test in Generation 1 was run, once, and recorded.

## What did NOT change
- No Stage-22, Stage-23, Stage-23.2, Stage-23.5 or Stage-24 artifact was touched
- No model was refitted on observed data, no hyperparameter re-selected, no threshold moved
- The plan digest and the out-of-fold table hash were asserted before a number was read

---

## Result

```text
  STAGE_25_RANKING_SUPPORTED          all six §8.10 criteria pass

  eligible clones            892      verified mechanically before scoring
  R(W1)  0.692654
  R(W4)  0.692176
  R(W5)  0.743781

  delta_RANK                +0.051605
  bootstrap CI95            [+0.037197, +0.065571]        lower endpoint > 0
  null p95                   0.008672
  margin over null p95      +0.042933
  null draws >= observed     0 / 1000
  p_perm                     1/1001 = 0.000999
  delta_TOP1                +0.115471   CI95 [+0.082960, +0.145740]
```

### The margin is the story, not the p-value

```text
  observed is 6.0x the null p95
  observed is 11.8 null SDs above the null mean
  the largest of 1,000 full-refit null draws is 0.013722 -- the observed value
    exceeds EVERY null draw
```

`p_perm = 0.000999` is the **floor** of a 1,000-permutation test and must be reported as
`p < 0.001 (0 of 1,000)`, not as a precise tail estimate — the same discipline applied to Role B's
`0/200`. The number that carries the weight is the separation, not the p-value: nothing the null
produced came close.

### It holds everywhere it was broken down

```text
  BY OUTER FOLD                          BY PRETREATMENT DEPTH
    fold 0   183 clones   +0.0435            1     277 clones   +0.0535
    fold 1   184 clones   +0.0548            2     143 clones   +0.0528
    fold 2   176 clones   +0.0658          3-4     169 clones   +0.0314
    fold 3   173 clones   +0.0506          5-9     165 clones   +0.0462
    fold 4   176 clones   +0.0435          10+     138 clones   +0.0779

  score-tie rate (W5)   0.0000    -- the 0.5 tie rule never had to fire
```

Positive in all five folds and all five depth bins. These are reported, not gating; §8.9 forbids
them rescuing anything, and here they are not asked to.

### The comparator choice was vindicated by the data

```text
  R(W1)  0.692654        nuisance + treatment
  R(W4)  0.692176        + additive X
  R(W5)  0.743781        + explicit X x U
```

`R(W4)` sits **below** `R(W1)` by 0.0005. The additive `X` term contributes nothing to *ordering*,
which is exactly why §8.5 named W4 the primary comparator: its additive term cannot create
clone-specific `X×U` ordering. The entire ranking gain is the interaction. `delta_RANK_FULL`
(+0.051127) is essentially identical to `delta_RANK` (+0.051605) for the same reason.

### delta_TOP1

```text
  LOW_PERSISTENCE_TOP1(W5)   0.828475
  LOW_PERSISTENCE_TOP1(W4)   0.713004
  delta_TOP1                +0.115471   CI95 [+0.082960, +0.145740]
```

Choosing the condition with the lowest predicted detection score finds a genuine zero for 82.8% of
eligible clones under W5 against 71.3% under W4. This was a **directional-consistency check**, not a
significance test (§8.8) — it could only withhold support, never grant it. It did not withhold.

## Tests
- 15 Stage-25 contracts, 0 skipped
- smoke test 7/7 before any long compute was committed
- full suite green

## Bugs found
- none in Stage 25. The 25C incomplete-null guard was verified firing at 855/1000 draws before the
  run finished, so its refusal path is exercised rather than assumed.

## Runtime

```text
  smoke                 66 s/draw single-process, unloaded
  real run             115 s/draw per shard across three shards
  total                ~10.7 h wall, against the 19-20 h accepted in §0.2
```

The ~1.7× per-shard penalty is the memory-bandwidth contention §0.2 predicted; three shards buy
about 1.7×, not 3×. The measurement came in under the accepted budget and nothing in the contract
changed because of it.

## Scientific interpretation

**Proves:** in the WM989 lineage system, under clone-held-out evaluation, an explicit
state×treatment interaction improves clone-specific ordering of six observed experimental
conditions over a non-interactive additive model. The improvement is `+0.0516` in equal-clone-
weighted within-clone AUROC, it exceeds all 1,000 full-refit permutation draws, and it is positive
in every outer fold and every pretreatment-depth stratum. The test was preregistered in full — the
metric, population, weighting, comparator, endpoint, null construction, permutation count and
verdict logic were all fixed before any of these numbers existed.

**Does NOT prove:**
- **that this generalises beyond WM989.** One lineage-traced system, six known conditions. No
  independent biological replication has been performed; clone-held-out folds and two endpoint
  families are not replication.
- **anything about unseen treatments.** The six conditions are the vocabulary; the tool returns
  `UNSUPPORTED_TREATMENT` outside it.
- **clinical utility.** C1 is an observed detection proxy, not death, sensitivity, resistance or
  patient response, and the six conditions include non-clinical stress contexts.
- **that state dominates.** Captured pretreatment clone abundance remains ~3.45× the whole state
  contribution. The ordering is abundance first, then treatment-specific state.
- **uniform benefit across conditions.** Four of six carry meaningful interaction; Cisplatin is
  negligible on C1 and Doxorubicin is negative on both endpoints.
- **anything about Role A.** It remains positive-but-underpowered supporting evidence with gate
  18.3 FAILED at 0.64 (audited ~0.45).

## Consequence for the tool

`validated_condition_order` is now exposed — and only because a real verdict file exists. Verified
end to end: loading the artifact without the verdict returns `ranking_status NOT_SUPPORTED` and no
order; loading it with the verdict returns `SUPPORTED` and an order. **The six scores are identical
either way.** The verdict unlocks a claim, not a computation.

## Next action
`GEN1_MANDATORY_SHIP`. Stage 26 records `KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`, then evidence lock,
claim lock, manuscript and reproducibility package. Independent new-system replication remains
Generation 2 and is not a Generation-1 gate.

No second ranking analysis is authorized. This is a terminal scientific result.
