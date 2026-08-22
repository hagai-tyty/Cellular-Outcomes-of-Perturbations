# stage_23_2C_ RECORD

## Goal
Ask how much of the positive Stage-23 permutation-null centre remains because the permuted
expression profile still carries technical depth / sparsity structure that the historical nuisance
block `B0` does not represent.

## Inputs
- 23.2A frozen protocol `78edd5d7f9900349…`; `Bdepth` frozen for all 3,147 clones
- the same 200 historical permutation mappings used by 23.2B (paired draw-for-draw)
- 23.2B's `D10` array, read from cache, so the 2×2 shares one construction

## Files added
- `results/stage23_2/stage23_2_depth_decomposition.json`
- `_cc_cache/stage23_2/decomposition_Bdepth.npz`

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2c`
- `tests/test_stage23_2_role_a_resolution.py` — 9 further contracts (55 total)

## What changed
- Nothing outside the additive 23.2C artifact

## What did NOT change
- `src/` unchanged · Stage 23 untouched · the historical observed `ΔAP` is preserved as cell `O00`
- **`Bdepth` was not promoted to anything.** It is a diagnostic block, per V2 §7.7

## Tests
- 55 passed · ruff clean

## Result

**STATUS: `RESIDUAL_DEPTH_STRUCTURE = SUPPORTED`** — runtime 98.9 min.

### The 2×2 of permutation-null centres

```text
                          B0            Bdepth
full K x C search     μ00 +0.00350    μ01 +0.00233
no-K-selection        μ10 +0.00299    μ11 +0.00161

both corrections together remove +0.00188  =  53.9% of the historical null centre
```

```text
contrast                      mean       95% CI                 frac draws > 0
depth   (μ00 - μ01)        +0.00117   [+0.00057, +0.00180]           0.860
selection (μ00 - μ10)      +0.00051   [-0.00010, +0.00124]           0.470
factor interaction         -0.00021   [-0.00075, +0.00036]           0.510
```

**Residual depth is the dominant mechanism — more than twice the selection contribution — and it is
the only one of the two whose CI excludes zero.** The factor interaction is indistinguishable from
zero, so the two mechanisms are approximately additive; neither masks the other.

### The mechanism is confirmed without touching an outcome

Spearman correlation between the donor profile's technical properties and the recipient clone's,
across all 200 historical mappings:

```text
donor nonzero-gene count  vs recipient log1p(total UMI)     median +0.3372  95% [+0.3357, +0.3391]
donor nonzero-gene count  vs recipient log1p(detected)      median +0.3440  95% [+0.3421, +0.3455]
donor total UMI           vs recipient total UMI            median +0.3306  95% [+0.3296, +0.3331]
```

All three exclude zero in the positive direction. **A permuted profile still resembles its recipient
technically at ρ ≈ 0.34.** That is the concrete reason the null centre is positive: the frozen strata
(`n_pretreatment_cells {1,2,3+} × n_lanes`) are too coarse to break the association — 82.1% of
clones sit in the single `1|1` cell, inside which depth varies freely and is preserved by the
shuffle. This is exactly what the Stage-23 design *intended* — the null deliberately preserves
abundance structure — so the finding is not that Stage 23E was wrong, but that its conservatism has
a measurable size.

### Observed 2×2

```text
O00  historical full search + B0        +0.01050
O10  no-K-selection        + B0         +0.00962
O01  full search           + Bdepth     +0.00694
O11  no-K-selection        + Bdepth     +0.00872
```

### The corrected same-data diagnostic — `NEGATIVE`, and only just

```text
O11        +0.00872
q95_11     +0.00835        O11 exceeds it
p_diag_11   0.0547         but this is > 0.05

CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC = NEGATIVE   (both gates are required)
```

It is worth being precise about *why* the corrected p-value improves at all, because the obvious
reading is wrong. The correction does **not** enlarge the effect:

```text
                       observed    null centre    separation    null p95
historical             +0.01050      +0.00350       +0.00700     +0.01455
corrected (cell 11)    +0.00872      +0.00161       +0.00711     +0.00835
```

The null-centred separation is essentially unchanged (`+0.00700 → +0.00711`). What changes is the
**spread** of the null: its 95th percentile falls from `+0.01455` to `+0.00835`. Removing the two
mechanisms makes the null tighter, not the signal bigger — and that alone moves `p` from `0.0846`
to `0.0547`, which still does not clear the frozen `0.05`.

### Lane-composition sensitivity — not executed

23.2A returned `WITHIN_R1_STRUCTURE_UNRESOLVED`, and V2 §7.5 permits the per-sample cell-count
block only under `WITHIN_R1_TECHNICAL_LANES`. It was therefore skipped, and a contract asserts no
`n_cells_GSM…` term appears anywhere in the artifact.

### Fixed-K arms under `Bdepth`

```text
K = 10  +0.00138      K = 20  +0.00151      K = 50  +0.00195      dispersion 0.00056
```

Same construction as 23.2B, which is what makes the interaction term interpretable.

## Bugs found
- None in this substage. Two of my own **test predicates** needed correcting during review: one
  computed a JSON blob it never used, and one iterated a name it never referenced — both dead code
  in checks that were otherwise sound

## Scientific interpretation

**Proves:** the positive Stage-23 permutation-null centre is substantially a *methodological*
artefact, and the larger share of it is residual technical structure rather than model-selection
breadth. Together the two mechanisms account for 53.9% of it. The depth mechanism is established at
the outcome level (CI excludes zero) *and* mechanistically, by an outcome-free measurement showing
permuted profiles retain ρ ≈ 0.34 technical similarity to the clones that received them.

**Does NOT prove:**
- **That Role A has a signal.** Correcting both mechanisms leaves `p_diag = 0.0547`, which fails
  its own frozen gate — and even a `POSITIVE` here would have been exploratory on already-inspected
  Rewind data, unable to emit `ROLE_A_CONFIRMATORY_SUPPORTED`. The historical `ROLE_A_SIGNAL_FAIL`
  is untouched and permanent.
- **That the effect got bigger.** It did not. The corrected p-value improves because the null
  narrows, not because the separation grows — `+0.00700` versus `+0.00711`.
- **That `Bdepth` is the right production nuisance block.** It was chosen to diagnose this null.
  V2 §7.7 forbids promoting it into Stage 24 without a frozen Stage-24 plan.
- **That the remaining 46% of the null centre is understood.** Two mechanisms were tested; the
  residue is unexplained, and label limitation (23.2D, `UNRESOLVED`) and power (23.2E) remain open.

## Next action
23.2E — power / identifiability. A runtime feasibility measurement is required first, per V2 §9.4.
