# stage_23_2B_ RECORD

## Goal
Measure how much of the positive Stage-23 permutation-null centre is attributable to `R3`'s broader
`K × C` model-selection path, using the historical permutation mappings so the comparison is paired
draw-for-draw.

## Inputs
- 23.2A frozen protocol `78edd5d7f9900349…`, `STAGE_23_2_PROTOCOL_FROZEN`
- `results/stage23_2/stage23_2_historical_null_d00.json` — the committed 200 historical `D00` values
- mapping set `6cff09bd5423b7f3…` (cache-only), regenerated from the frozen seed
- frozen Stage-22 outer folds, `y_primed`, and the historical `B0` nuisance block

## Files added
- `results/stage23_2/stage23_2_model_selection_decomposition.json`
- `_cc_cache/stage23_2/decomposition_B0.npz` (per-draw arrays, gitignored)

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2b`; 23.2A untouched
- `tests/test_stage23_2_role_a_resolution.py` — 10 further contracts (36 total)

## What changed
- Nothing in Stage 23, Stage 22, or 23.2A. Only additive 23.2B artifacts were written

## What did NOT change
- `src/` unchanged · the historical `D00` array is **read**, never recomputed for a statistic
- the historical observed `ΔAP = +0.01050` is unchanged and was not reinterpreted

## Tests
- 36 passed · ruff clean

## Result

**STATUS: `MODEL_SELECTION_NULL_INFLATION = UNRESOLVED`** — runtime 89.6 min.

### Engine fidelity, checked before any statistic was read

```text
cell 00 recomputed from the shared-scores engine
    reproduces the committed D00 array           EXACTLY (max abs diff 0.0)
expression-free reference R1
    reproduces the historical AP                 EXACTLY (0.010346763671910921)
```

The engine scores all twelve `(K, C)` candidates once per draw and then applies each selection rule
as a different `argmax` over the same scores. That is an identity, not an approximation — and the
`cell 00` check above is what proves it, since `hist12` must land on the historical value bitwise.

### Primary paired comparison

```text
μ00   historical full K x C search        +0.00350
μ10   no-K-selection (3-arm equal mean)   +0.00299

selection_shift = mean(D00_j - D10_j)     +0.00051
                  95% CI                  [-0.00011, +0.00122]
                  fraction of draws > 0    0.470

fraction of the null mean explained by search    14.5%
```

The CI includes zero, so the status is `UNRESOLVED` by the frozen rule.

**The substantive number is the one the status does not capture: 85% of the positive null centre
survives with no K selection at all.** Removing the entire `K` selection path moves the null mean
from `+0.00350` to `+0.00299`.

The `fraction of draws > 0 = 0.470` is worth recording alongside the positive mean: per draw, the
sign of `S_j` is close to a coin flip, so the positive mean is carried by a skewed minority of
draws rather than by a consistent per-draw advantage.

### Fixed-K arms (design change 1) — reported as diagnostics

```text
arm       paired mean      95% CI
K = 10      +0.00301   [+0.00233, +0.00379]
K = 20      +0.00264   [+0.00215, +0.00320]
K = 50      +0.00331   [+0.00266, +0.00404]

arm dispersion  0.00067
```

Dispersion is small relative to the arms themselves, so the equal-weight reference is not hiding a
wide spread. This is the check V1's single `K = 20` arm could not have provided — and note that
`K = 20`, the arm V1 would have used alone, is the *lowest* of the three, which would have inflated
`selection_shift` by roughly 0.0004 relative to the balanced reference.

### Search-width ladder — unconditional in V2

```text
 4-candidate  (fixed K, mean of arms)   +0.002988
 8-candidate  (K ∈ {10,20})             +0.003069
12-candidate  (K ∈ {10,20,50})          +0.003495     = historical

monotone increase: TRUE
```

Monotone, which is descriptive support that wider search does raise the null centre — while the
primary CI says the total effect is not distinguishable from zero at 200 draws.

### Observed-data sensitivity (V2 §6.5, effect attribution only)

```text
rule              observed ΔAP
hist12 (frozen)      +0.01050
K = 10               +0.01013
K = 20               +0.00915
K = 50               +0.00957
no-K-selection       +0.00962
ladder8              +0.00755
```

The null-centred separation is almost unchanged by removing K selection:

```text
historical      0.01050 - 0.00350 = 0.00700
no-K-selection  0.00962 - 0.00299 = 0.00663
```

**Removing model selection does not rescue Role A.** Both the observed statistic and the null centre
fall by similar small amounts, so the gap between them barely moves. No p-value was recomputed and
the historical failure is untouched.

## Bugs found
- None in this substage. The two fidelity checks (`cell 00` bitwise, `R1` bitwise) were built in
  precisely so that a silent engine divergence would surface as a failure rather than as a plausible
  number, and both passed on the first run

## Scientific interpretation

**Proves:** the historical comparison can be re-run draw-for-draw against a no-K-selection reference
that spans the same grid, and doing so leaves most of the positive null centre in place. Model
selection contributes something — the ladder is monotone and the mean shift is positive — but at 200
draws its total contribution is not distinguishable from zero, and it is at most about a seventh of
the null mean.

**Does NOT prove:**
- **That model selection contributes nothing.** `UNRESOLVED` is not `NOT_SUPPORTED`. The monotone
  ladder and the positive mean both point the same way; the sample size does not resolve it.
- **That the remaining null centre is biological.** 23.2B removes one candidate mechanism from
  contention as the *dominant* explanation. What produces the other 85% is 23.2C's question, and the
  realized strata recorded in 23.2A — 82.1% of clones in a single cell — are the reason residual
  depth structure is the live hypothesis.
- **Anything about Role A's status.** The historical verdict is permanent. The observed sensitivity
  above is attribution, not a corrected test.

## Next action
23.2C — residual depth / nuisance decomposition. Started.
