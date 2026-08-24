# STAGE 23.2 — ROLE-A CONFIRMATION PROTOCOL V1

**Frozen at 23.2F, before any untouched confirmation evidence was inspected.**
Stage-23.2 protocol `78edd5d7f9900349925339169a5d5e3e5011fe23e3c7c22608ac98bfe3427bf4`.

This protocol may only be executed in 23.2G, and only on evidence that was **not** used to design
the correction. Executing it on the already-inspected Rewind biological replicate R1 is forbidden:
that analysis is exploratory by construction and is already recorded in 23.2C.

---

# 1. Confirmed scientific hypothesis

Under the historical Rewind outcome and the frozen Stage-22/23 evaluation geometry, pretreatment transcriptional state predicts the Role-A outcome beyond a DEPTH-COMPLETE nuisance baseline Bdepth = [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI), log1p(n_detected_GE_features_in_raw_pseudobulk)].

Mechanically derived: RESIDUAL_DEPTH_STRUCTURE is the only SUPPORTED mechanism that implies a correction. MODEL_SELECTION_NULL_INFLATION is UNRESOLVED, so search matching is not indicated; OUTCOME_LABEL_LIMITATION is UNRESOLVED, so V2 §8.7 forbids adopting an alternative outcome. The correction was not chosen by effect size.

# 2. Outcome definition

```text
historical hard label, top-100 gDNA barcodes with ties, unchanged
```

The historical hard label is retained deliberately. `OUTCOME_LABEL_LIMITATION` is UNRESOLVED, not
SUPPORTED, so V2 §8.7 forbids adopting an alternative outcome representation here. If a future
substage establishes the limitation, this protocol must be re-frozen and the material
benchmark-change firewall (V2 §10.7) applies.

# 3. Allowed input X

```text
pretreatment Gene Expression only, 36,601 features
clone pseudobulk: sum RAW counts, then CP10K and log1p exactly once
no lineage-assay feature, no clone identifier, no outcome-derived quantity
```

# 4. Nuisance block

```text
  log1p(n_pretreatment_cells)
  n_lanes
  log1p(total_raw_GE_UMI)
  log1p(n_detected_GE_features_in_raw_pseudobulk)
```

All continuous terms standardised training-only. This is the correction under test: the confirmation
asks whether state predicts outcome **beyond** this depth-complete baseline.

# 5. Model family, grid and preprocessing

```text
family        l2 logistic regression, liblinear, class_weight None
C grid        [0.01, 0.1, 1, 10]
K grid        [10, 20, 50]   (PCA fitted once at max K; K values are prefixes)
preprocessing training-only gene filter, gene scaler, PCA, PC scaler -- refitted inside every
              inner split
inner CV      3-fold StratifiedKFold(shuffle=True, random_state=23023)
```

Unchanged from the historical pipeline. The search path is **not** matched or narrowed:
`MODEL_SELECTION_NULL_INFLATION` is UNRESOLVED, so narrowing it is not an indicated correction.

# 6. Grouping unit

```text
clone. Outer folds grouped at clone level. Cells are never independent outcome replicates.
```

# 7. Primary metric and statistic

```text
metric      average precision at clone grain
statistic   delta_AP = AP(PCA(X) + Bdepth) - AP(Bdepth)
```

# 8. Permutation / null design

```text
whole clone-level expression profiles permuted as intact vectors
within outer-training and within outer-test separately
inside strata: n_pretreatment_cells {1,2,3+} x n_lanes
200 permutations, full nested-CV rerun per draw
p_perm = (1 + #{null >= observed}) / (n_perm + 1)
```

**Known conservatism, declared in advance.** 23.2C established that this null retains technical
similarity between donor and recipient profiles at Spearman ~0.34, which is why its centre is
positive. The design is retained unchanged so the confirmation is comparable to the historical
test; the depth-complete nuisance block is the correction, not a weakened null.

# 9. PASS threshold

```text
observed > null p95   AND   p_perm <= 0.05
```

Both required. No alternative threshold, one-sided relaxation or post-hoc adjustment is permitted.

# 10. Minimum positive-count / design requirement

```text
  scale 1:   35 positive clones -> power 0.290
  scale 2:   70 positive clones -> power 0.520
  scale 4:  140 positive clones -> power 0.940

  requirement: >= 140 positive clones at oracle AUC 0.66
```

the ladder is coarse (35 / 70 / 140 positive clones). The requirement is the smallest TESTED cohort reaching 0.80, not an interpolated value between rungs -- V2 §9.5 forbids extrapolating a precise required N.

A confirmation cohort below this floor may be run, but a null result from it is **not** evidence
against the hypothesis and must be reported as underpowered.

# 11. Source-qualification criteria

```text
pre-intervention molecular measurement present
independently measured later outcome
clone / lineage linkage sufficient for grouped evaluation
raw or processed data sufficient to reconstruct the endpoint
no same-state outcome leakage
at least the section-10 positive-count floor
a biological replicate INDEPENDENT of R1
```

The last criterion is what distinguishes confirmation from more of the same: R1 is consumed.

# 12. Search budget for confirmation data

```text
one bounded search, executed only after this protocol is committed
candidates enumerated from the 23.2A reserved ledger first
then a single external search pass
inclusion criteria may NOT be relaxed if results are sparse
```

# 13. Forbidden data — already inspected

```text
GSE227151 biological replicate R1  (GSM7092515, GSM7092516)
the Stage-22 Rewind benchmark and every Stage-23 / 23.2 artifact derived from it
the gDNA table stepThreeStarcodeShavedReads_BC_gDNA.txt
```

Reserved candidates `GSM7092517`-`GSM7092521` have been read at declared-metadata level only and
remain eligible. Whether they carry a reconstructable Role-A outcome is **UNVERIFIED** and must be
established under section 11 before any performance quantity is computed.

# 14. Stage-27 firewall

Any dataset used here becomes development/confirmation evidence and is **not** an untouched
Stage-27 replication set. Stage 27 must preserve an independent biological test of the eventual
frozen Stage-24 model.
