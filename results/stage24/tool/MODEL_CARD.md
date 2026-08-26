# CellFate-Rx Gen-1 — Model Card

**Model** `gen1-w5-c1-v1` — W5 = `X + B + U + X*U`
**Feature contract** `wm989-ge-36601-v1`
**Endpoint** C1, post-treatment clone detection
**System** WM989 (GSE279162), 1401 lineage-traced clones, six conditions

## What it does
For one starting clone it returns a `future_detection_score` for each of the six observed
conditions: the model's propensity that the clone is still detected after that condition.

## What it does NOT do
```text
  known conditions only; UNSUPPORTED_TREATMENT for anything outside the six
  requires the complete frozen nuisance block B; it may not be imputed
  trained on clone-level pseudobulk, so a single cell is not an equivalent input
  captured pretreatment clone abundance remains ~3.45x the state contribution
  no calibrated probability; the score is not a calibrated risk
  no independent biological replication of the Role-B finding
  ranking is NOT validated until Stage 25 records RANKING_SUPPORTED
```

## Validation
```text
  frozen out-of-fold reproduction   max |diff| vs the Stage-23 frozen column   2.498e-16
  determinism, same session         True
  determinism, across loads         True
  fold isolation                    every fold component verified disjoint from its test clones
  eligible ranking clones           892 of 1401
```

The **deployment** component is packaging, not validation: it is the same specification and the same
selection rule fitted once on all clones, and it is **not** validated on held-out data. Its
performance is estimated by the frozen out-of-fold result, which came from the fold components.

## Ranking
`ranking_status` is `NOT_SUPPORTED` until Stage 25 records `STAGE_25_RANKING_SUPPORTED` under its
pre-registered test. Until then the six scores are returned but their **order is not a validated
condition ranking** and `validated_condition_order` is withheld.

## Intended use
Research use on WM989-like lineage-traced pretreatment data, for the six observed experimental
conditions. **Not** a clinical tool, **not** a treatment recommendation, **not** a calibrated
probability, and **not** applicable to unseen treatments, other cell lines, or patients.
