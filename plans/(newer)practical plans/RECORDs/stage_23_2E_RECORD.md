# stage_23_2E_ RECORD

## Goal
Determine which effect sizes the historical Role-A pipeline can actually detect at different
rare-positive counts, and separate that question from the one no simulation can answer — whether
the evidence rests on a single biological replicate.

## Inputs
- 23.2A frozen protocol `78edd5d7f9900349…`; `Bdepth` frozen
- label-free synthetic direction `z`: clone-level `X` residualised on `Bdepth`, one PCA, PC1
  deterministically oriented and standardised. `y_primed` is never touched
- frozen seeds: covariate resample 23440, null 23441, alternative 23442, beta calibration 23443

## Files added
- `results/stage23_2/stage23_2_power_identifiability.json`
- `_cc_cache/stage23_2/power_shards/*.json` (9 shards, gitignored)

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2e`, `--stage 23.2e-merge`
- `tests/test_stage23_2_role_a_resolution.py` — 11 further contracts (66 total)

## What changed
- Nothing outside the additive 23.2E artifacts. No Stage-23 or earlier-substage file touched

## What did NOT change
- `src/` unchanged · no simulation count reduced, no scale dropped, no early stopping
- the historical Role-A verdict is untouched

## Tests
- 66 passed · ruff clean

## Result

**`WITHIN_R1_EVENT_COUNT_LIMITATION = SUPPORTED`**
**`BIOLOGICAL_REPLICATION_LIMITATION = SUPPORTED`**

1,200 full nested-pipeline runs, all at the frozen counts.

```text
scale   N        pos/fold   q95_null     power @ AUC 0.66   power @ AUC 0.70
  1     3,147        7      +0.00644          0.290              0.310
  2     6,294       14      +0.00466          0.520              0.740
  4    12,588       28      +0.00227          0.940              1.000
```

Beta calibration landed on target at every scale (achieved median oracle AUC 0.658–0.674 against
0.66; 0.694–0.702 against 0.70).

**At the historical cohort size, a genuine state signal of oracle AUC 0.66 would be detected 29% of
the time.** Doubling the events to 70 raises that to 52%; quadrupling to 140 raises it to 94%. The
frozen rule — power below 0.50 at scale 1, reaching 0.80 or better at scale 2 or 4 — returns
`SUPPORTED`.

Two features of the ladder are worth recording beyond the headline:

- **The null tightens faster than the signal grows.** `q95_null` falls from `+0.00644` to `+0.00227`
  across the ladder while the alternative mean `ΔAP` stays roughly flat (`+0.0056 → +0.0051` at
  0.66). Power rises mainly because the null narrows, not because the estimator gets better at
  recovering the signal.
- **The effect-size gradient is weak at scale 1 and steep later.** Going 0.66 → 0.70 buys only
  0.290 → 0.310 at 35 positives, but 0.520 → 0.740 at 70. At the historical event count the
  procedure is close to insensitive to how strong the signal actually is.

### The second status is a design fact, not an estimate

`BIOLOGICAL_REPLICATION_LIMITATION = SUPPORTED` because the claim rests on
`n_biological_replicates = 1`. Per V2 §9.5.2 it admits only `SUPPORTED` or `NOT_SUPPORTED` —
`UNRESOLVED` is not an available value, because the replicate count is known rather than inferred —
and it can become `NOT_SUPPORTED` only when a Role-A claim is supported by two or more independent
biological replicates. A contract asserts that a power curve reaching 1.000 does **not** clear it.

### Bias direction is not claimed

Scales 2 and 4 are empirical resamples of biological replicate R1. Repeated covariate profiles
reduce effective covariate diversity and therefore make the projected curve an approximation whose
bias direction is **not guaranteed**. The study estimates within-R1 event-count detectability under
the empirical R1 covariate distribution; it does not estimate power gained from additional
independent biological replicates. A contract fails if the words "optimistic" or "conservative"
appear anywhere in the artifact.

### Runtime — and a large execution failure worth recording

```text
shard                 runtime      note
scale1_null            482.2 min   3-way contended
scale1_alt_auc66       227.9 min   3-way contended
scale1_alt_auc70       313.5 min   3-way contended
scale2_null            791.8 min   3-way contended
scale2_alt_auc66       111.9 min   2 processes, --threads 4
scale2_alt_auc70       115.0 min   2 processes, --threads 4
scale4_null            354.4 min   2 processes, --threads 4, resumed from 5
scale4_alt_auc66       229.1 min   2 processes, --threads 4
scale4_alt_auc70       229.1 min   2 processes, --threads 4
```

The first four shards ran roughly **four times slower per simulation** than the last five for the
same work. `OMP_NUM_THREADS` was unset, so three concurrent processes each spawned eight OpenBLAS
threads on an eight-core machine — twenty-four threads contending for eight cores. Aggregate
throughput fell *below* running the shards one at a time. Two processes with `--threads 4` fixed
it: `scale2_alt_auc66` ran at 1.11 min/sim where the contended `scale2_null` had managed 3.96.

## Bugs found

1. **My runtime estimates were wrong three times, by up to 5×** — 13 h, then 17 h, then 25 h,
   against an actual ~5 h once the threading was fixed. Every estimate came from solo probes that
   never reproduced the contention of the real run. The lesson recorded here: measure the
   configuration you will actually run, not a convenient one.
2. **Thread oversubscription** as above — a genuine engineering fault, not merely a bad estimate.
   Parallelism made the job slower than sequential.
3. **V2 §18 resume was not implemented** until it was needed. When the contended run had to be
   stopped, three partial shards were unrecoverable. Resume was then added (execution-neutral,
   verified bit-identical on five already-completed simulations) and immediately proved itself:
   `scale4_null` resumed from 5 banked simulations.
4. **A cosmetic key-format inconsistency** in the merge: `str(0.70)` yields `"0.7"` while
   `str(0.66)` yields `"0.66"`. Normalised to two decimal places and re-merged. No value changed

## Scientific interpretation

**Proves:** the historical Role-A experiment was underpowered for the effect sizes it was asked to
detect. With 35 positive clones the pipeline recovers a true AUC-0.66 signal 29% of the time, so a
null result at that event count carries little evidence against a real effect of that size. It also
establishes, separately, that the entire body of Role-A evidence rests on one biological replicate.

**Does NOT prove:**
- **That a Role-A signal exists.** The study injects a synthetic signal and asks whether the
  procedure would find it. It says nothing about whether Rewind contains one.
- **That more clones would settle the question.** Scales 2 and 4 are resamples of R1, so the curve
  answers "more events from this same biological context", not "more biology". The replication
  limitation is untouched by every point on the ladder.
- **That AUC 0.66 is the true effect size.** It was chosen in advance to approximate the historical
  R3 discrimination, and the 0.70 curve is planning sensitivity only.
- **Anything about the corrected same-data analysis.** 23.2C's `p_diag = 0.0547` is a separate,
  exploratory result and is not evidence that the corrected effect is real.

## Next action
23.2F — diagnostic synthesis and confirmatory protocol freeze.
