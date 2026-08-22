# stage_23_2D_ RECORD

## Goal
Ask whether the frozen Role-A outcome — an operational top-gDNA-barcode threshold label — is stable
enough to carry a hard 35-positive / 3,112-negative classification claim. This substage studies the
**measurement**, not predictive performance, and fits nothing.

## Inputs
- 23.2A frozen protocol `78edd5d7f9900349…`
- `stepThreeStarcodeShavedReads_BC_gDNA.txt` — 49,554 rows, 1,936 barcodes, `SampleNum = 3` only
- frozen Stage-22 clone table and `y_primed`

## Files added
- `results/stage23_2/stage23_2_label_reliability.json`

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2d`; earlier substages untouched
- `tests/test_stage23_2_role_a_resolution.py` — 10 further contracts (46 total)

## What changed
- Nothing outside the additive 23.2D artifact

## What did NOT change
- `src/` unchanged · no predictor fitted · no alternate label became a prediction target
- the frozen N=100 outcome remains the historical target

## Tests
- 46 passed · ruff clean

## Result

**STATUS: `OUTCOME_LABEL_LIMITATION = UNRESOLVED`** — runtime 0.02 min.

### This is the case V2 design change 3 was written for

All four of V1's `NOT_SUPPORTED` criteria are met:

```text
mean frozen-positive retention    0.9751   >= 0.90   PASS
positives with P(selected) < 0.80      2   <= 3      PASS
Jaccard(top90,  top100)           0.9429   >= 0.90   PASS
Jaccard(top110, top100)           1.0000   >= 0.90   PASS
```

**Under V1's rule this substage would have returned `NOT_SUPPORTED`** — the ledger would have
recorded "the outcome label is sound" on the strength of two diagnostics that observe only
sequencing-count noise and cutoff position. V2 requires, in addition, independent outcome-assay
replication of the same clones. The Rewind materials contain none, so the status is `UNRESOLVED`
and `NOT_SUPPORTED` is formally unreachable here. A contract asserts exactly this, and would fail
if the gate were ever relaxed.

### Multinomial sampling stability (V2 §8.3)

5,000 conditional resamples, seed 23431, over the single pooled selection unit
(`N = 782,826` counts across 1,936 barcodes).

```text
mean frozen-positive retention        0.9751
median                                (see artifact)
minimum                               0.4322
positives with P(selected) < 0.50          1 / 35
positives with P(selected) < 0.80          2 / 35

positive clones per draw     mean 34.49   sd 0.85   range 32-38
Jaccard vs the frozen set    mean 0.9653   p05 0.9167
expected intruding retained clones per draw   0.364
```

### The set is stable; its boundary is not

```text
rank    counts   tie   gap to prev   gap to next
  98     2380     1        14            12
  99     2368     1        12             3
 100     2365     2         3             0
 101     2365     2         0             1
 102     2364     1         1            14
```

**A single gDNA count separates a positive from a negative.** Rank 102 sits at 2,364 against the
2,365 cutoff. The rank-100 clone is retained in only 43.2% of resamples — a coin flip — and it is
the same tie that turns "top 100" into 101 selected barcodes.

The honest summary is that these two facts are not in tension: the label is stable *as a set*
(96.5% Jaccard) while being arbitrary *at its margin*. Which of those matters depends on whether a
downstream claim leans on the marginal clones, and with only 35 positives a single boundary clone
is 2.9% of the positive class.

### Cutoff sensitivity (V2 §8.4)

```text
        barcodes  positives  Jaccard   lost  gained
top80         80         32   0.9143      3       0
top90         90         33   0.9429      2       0
top100       101         35   1.0000      0       0
top110       110         35   1.0000      0       0
top120       120         37   0.9459      0       2
```

`top110` selects nine more barcodes than `top100` yet yields the identical clone set — those
barcodes are not among the 3,147 retained clones. The frozen N=100 outcome is unchanged and no
alternate N was selected for any reason.

### Removed by V2

`cross_gsm_gdna_concordance` is recorded as `REMOVED IN V2`. The gDNA table carries `SampleNum = 3`
on every row, so per-GSM outcome support is not identifiable; a contract asserts neither GSM
accession appears anywhere in this artifact.

## Bugs found
- None in the analysis. One of my own **test predicates** was wrong: a check for predictive
  quantities scanned for the substring `predict`, which matched the artifact's own
  `no_predictive_model_fitted` declaration — the check was flagging the very statement it existed
  to verify. Replaced with a scan for actual metric names over the non-declarative fields

## Scientific interpretation

**Proves:** under the two noise sources that can be modelled from the available materials — resampling
the pooled gDNA library, and moving the cutoff rank — the frozen 35-positive set is reproduced
closely, and the source rule reconstructs it exactly. It also proves the label's boundary is decided
by differences of one count, and that one frozen positive is a coin flip under resampling.

**Does NOT prove:**
- **That the label is a reliable measurement of the biological event.** That is precisely what these
  diagnostics cannot show, which is why `NOT_SUPPORTED` is unreachable. Colony-level biological
  sampling, PCR duplication and assay-to-assay variation are unobserved, and the multinomial model
  is a *lower bound* on instability.
- **That label noise is not contributing to the Role-A failure.** `UNRESOLVED` means undetermined,
  not absent.
- **That an alternate threshold would be better.** No alternate N was evaluated for performance, and
  none may be, per V2 §3.7.

## Next action
23.2C is running; 23.2E — power / identifiability — follows.
