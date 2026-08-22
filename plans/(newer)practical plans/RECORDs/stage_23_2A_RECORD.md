# stage_23_2A_ RECORD

## Goal
Prove the Role-A failure decomposition is executable before any diagnostic model is fitted: verify
the closed Stage-23 state and every frozen anchor, settle the within-R1 sample structure as far as
source evidence permits, recover the exact historical permutation basis, and freeze an immutable
protocol identity. **23.2A fits no diagnostic predictor.**

## Inputs
- plan: `STAGE_23_2_ROLE_A_RESOLUTION_V2.md` (V1 audited, preserved at `c06fc98`, archived)
- Stage-23 closure anchored at `2e04ccf`; artifacts verified by hash, not by trust
- `D:\GSE227151_Rewind\` — 2 control GSMs, 3 barcode tables, GEO series matrix + family XML,
  author code
- frozen Stage-22 benchmark: 3,147 clones, 35 positives, 5 outer folds

## Files added
- `experiments/run_stage23_2_role_a_resolution.py`
- `tests/test_stage23_2_role_a_resolution.py` (26 contracts)
- `results/stage23_2/stage23_2_protocol.json`
- `results/stage23_2/stage23_2_source_design.json`
- `results/stage23_2/stage23_2_reserved_confirmation_candidates.json`
- `results/stage23_2/stage23_2_historical_null_d00.json`
- `results/stage23_2/stage23_2_bdepth.csv`
- `results/stage23_2/stage23_2a_results.json`

## Files modified
- `tests/test_ci_portability.py` — the new module joins the CI-condition guard

## What changed
- The protocol identity is a **canonical-JSON digest of the scientific surface only** (V2 §4.1),
  with git commit, source hashes, dependency versions and platform recorded *beside* it (§4.2)
- The 200 historical `D00` values are now a committed artifact; the mapping table stays cache-only
  with a committed digest
- `Bdepth` exists as a frozen, outcome-free table for all 3,147 clones

## What did NOT change
- `src/` unchanged · no Stage-23 artifact rewritten · no Stage-22 label, fold or mapping re-derived
- no diagnostic predictor fitted · no reserved-candidate matrix downloaded or opened

## Tests
- 26 passed · ruff clean
- Mutation-tested: flipping the within-R1 status without its consequence flag, flipping the replay
  reproduction flag, and corrupting one replayed `D00` value broke four separate contracts

## Result

**VERDICT: `STAGE_23_2_PROTOCOL_FROZEN`** — all ten V2 §5.7 gates pass. Runtime 97.0 min.

```text
stage23_2_protocol_sha256   78edd5d7f9900349925339169a5d5e3e5011fe23e3c7c22608ac98bfe3427bf4
```

### Preflight — 36 checks, 0 failed

Every §2 anchor reproduced: 3,147 / 35 / 3,112 clones, 7 positives per fold, R0–R3 pooled AP
(0.01112 / 0.01035 / 0.01923 / 0.02085), R3 ROC-AUC 0.6628, ΔAP `+0.01050` with 95% CI
`[+0.00397, +0.02258]`, and the full historical null. Gate confirmed `STAGE_24_BLOCKED_ROLE_A`,
Role-A `ROLE_A_SIGNAL_FAIL`, both Role-B verdicts unchanged, `STRUCTURAL_CONTROLS_PASS` true, and
the closure record still declaring Stage 23 formally closed with the permanence rule intact.

Two V2 corrections confirmed against the repository rather than assumed:

```text
F3  results/stage23_permutation_results.json holds SUMMARY statistics only;
    the per-permutation array was never committed
F4  there is no separate historical bootstrap artifact -- the block lives inside
    stage23_rewind_results.json (replicates 2000, seed 23123)
```

### The historical permutation basis is exactly recoverable

All 200 draws regenerated from `default_rng(23323 + b)` through the frozen `permute_within`, once
per outer fold in fold order. Every mapping verified train→train, test→test, within-stratum and
bijective.

```text
replayed vs committed        replayed          committed         match
mean                   +0.003495157024   +0.003495157024        yes
sd                     +0.006305127914   +0.006305127914        yes
p95                    +0.014545491038   +0.014545491038        yes
max                    +0.051437072593   +0.051437072593        yes
min                    -0.001046290766   -0.001046290766        yes
null >= observed                    16                 16        yes
p_perm                  0.084577114428     0.084577114428        yes

bitwise identical to the Stage-23 cache: TRUE (all 200 values)
```

```text
mapping-set sha256   6cff09bd5423b7f38b9cba12b88050e411fcc6b66c1456de7f58e1f5f243e74f
mapping rows         3,147,000   (200 draws x 5 folds x 3,147 clones), cache-only
D00 values sha256    fcba622b9d4c483e...   committed
mean fixed clones per fold   9.96 / 9.70 / 9.79 / 10.50 / 9.88
```

### Realized permutation strata — five non-empty cells

```text
1|1   2,584      2|1   220      2|2   196      3+|2   110      3+|1   37
```

`1|2` cannot exist: one pretreatment cell cannot span two lanes. **82.1% of clones sit in a single
stratum**, which is the structural reason 23.2C's residual-depth question is live.

### Within-R1 source design — `WITHIN_R1_STRUCTURE_UNRESOLVED`

The biological-replicate question is settled and was not reopened: the benchmark records
`biological_replicate = R1` for all 3,905 cells with `generalization_scope =
within_R1_clone_heldout`, and GEO titles read *"biol rep 1, sample 1"* and *"biol rep 1, sample 2"*.
**One biological replicate, two samples.**

The finer question is not resolvable from the local materials:

```text
declared characteristics       IDENTICAL across both GSMs
                               (tissue, cell line, cell type, genotype, treatment)
metadata declares lane split   NO
metadata declares separate
  culture / harvest            NO
distinct BioSample accessions  YES   <- non-discriminating
distinct SRA experiments       YES   <- non-discriminating
per-sample GEM loading         YES   <- non-discriminating
per-sample 10x indexing        YES   <- non-discriminating
```

The last four are true of a single suspension split across two 10X channels *and* of two
separately handled libraries, so they were deliberately excluded from the decision rule.

**Consequence: V2 §7.5 lane-composition sensitivity is FORBIDDEN in 23.2C.** Under an unresolved
status, per-sample cell counts could absorb biological-unit structure.

### Sample-numbering conflict (F2), recorded not repaired

```text
GEO title      GSM7092515 = "sample 1"     GSM7092516 = "sample 2"
file naming    GSM7092515_1_2_control_*    GSM7092516_1_1_control_*
benchmark      SampleNum 2 -> GSM7092515   SampleNum 1 -> GSM7092516
```

The V2 tie-break applies: author `SampleNum` / file naming wins, because the benchmark and the
author's own barcode tables are keyed on it. The frozen Stage-22 mapping was **not** re-derived; a
contract asserts it still matches the live benchmark.

### gDNA source rule reproduced exactly

```text
grouping           (BC50StarcodeD8, SampleNum) -> sum -> slice_max(n=100, with_ties=TRUE)
SampleNum in gDNA  {3} only -- grouping by it is a NO-OP (M4)
support column     `counts`, not `nUMI` (M4)
rows / barcodes    49,554 rows, 1,936 distinct barcodes, total N = 782,826
rank-100 cutoff    2365, tie size 2  ->  101 barcodes selected
positive clones    35   EXACT match to the frozen label
```

### `Bdepth` — frozen, outcome-free, complete

```text
total_raw_GE_UMI                       min 504    median 24,558    max 283,308    sum 90,850,570
n_detected_GE_features_in_pseudobulk   min 380    median  4,969    max  10,506
```

The detected-feature count matches the frozen 23A matrix's non-zero pattern **exactly for all
3,147 clones**, which independently validates the raw re-accumulation against the normalised cache.

### Reserved confirmation ledger — metadata only

13 declared samples: 2 used by Stage 23, **3 reserved candidates** (`GSM7092517`/`GSM7092518` =
biol rep 2, `GSM7092519` = biol rep 3), 2 reserved but sorted for cycling status (different
design), 6 outcome-side iPS samples. No matrix was downloaded or opened; every
`matching_future_outcome_declared` is `null` and stays `UNVERIFIED` until 23.2F freezes the
confirmation protocol. A contract asserts the builder never names a reserved accession.

## Bugs found
Three, all in my own 23.2A code, all caught before any result was trusted:

1. **GEO parser read one line of four.** `_series_field` returned only the first matching line, but
   `!Sample_extract_protocol_ch1` spans several — half the protocol was invisible. The status came
   out `UNRESOLVED` either way, but on incomplete evidence.
2. **The derivation rule over-claimed.** It would have emitted `WITHIN_R1_SEPARATE_LIBRARIES` from
   distinct BioSample/SRA accessions plus per-sample GEM loading and indexing — none of which
   distinguish the two candidate designs. Rule tightened to require genuinely discriminating
   metadata; the non-discriminating evidence is still recorded, just not used to decide.
3. **The preflight made 23.2A non-re-runnable.** It forbade *any* pre-existing Stage-23.2 artifact,
   including 23.2A's own, so the substage could never be re-executed — which would have made a
   determinism check impossible. Scoped to later-substage artifacts only, per V2 §5.1's intent.

## Scientific interpretation

**Proves:** the decomposition experiment is executable exactly as frozen. The historical null is
not merely re-samplable but *exactly* recoverable — all 200 values reproduce bitwise — so 23.2B and
23.2C can be genuinely paired against the historical comparison rather than against a fresh sample
of it, and the Monte-Carlo noise a paired design removes really is removed. The frozen label
reconstructs exactly from source, and `Bdepth` is computable for every clone without touching an
outcome.

**Does NOT prove:** anything about why Role A failed. No diagnostic model has been fitted and no
mechanism has been tested. In particular the `WITHIN_R1_STRUCTURE_UNRESOLVED` finding is a
statement about the available metadata, not evidence that the two samples are equivalent — and it
*removes* a diagnostic (§7.5) rather than adding one. The biological-replicate count remains 1,
which no substage of 23.2 can change.

## Next action
23.2B — model-selection null decomposition. Started.
