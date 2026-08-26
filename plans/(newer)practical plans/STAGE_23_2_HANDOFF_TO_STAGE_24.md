# STAGE 23.2 → STAGE 24 HANDOFF

**Status: NOT Stage-24-ready.** 23.2H has now been executed under confirmation protocol V5 and
returned `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE`. Five of the six §18 gates pass; the measured
design-power gate does not. Stage 24 remains BLOCKED.

> **23.2G's `QUALIFYING_SET_EMPTY` verdict is superseded.** It rested on the premise that
> biological replicates 2 and 3 carry no later outcome. That premise was false: their outcome
> materials were never deposited in GEO, and a GEO-complete search therefore returned a true answer
> to the wrong question. `stage_23_2G_RECORD.md` stands unedited; see
> `RECORDs/stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md` and `RECORDs/stage_23_2H_RECORD.md`.

## 23.2H result — independent confirmation on replicates 2 and 3

```text
  cohort            2310 clones, 4078 cells, 75 positive clones, prevalence 3.25%
                    rep 2 = GSM7092517/18 (S4+S5)    rep 3 = GSM7092519 (S1)    ungated only

  PRIMARY ARM       delta_AP  +0.030504     null p95 0.028101     p_perm 0.0299
    rep 2           delta_AP  +0.058251     CI95 [+0.003038, +0.163776]
    rep 3           delta_AP  +0.023291     CI95 [-0.035922, +0.088588]
    transfer        2->3 +0.019420          3->2 +0.028963

  SENSITIVITY ARM   delta_AP  +0.021686     null p95 0.009472     p_perm 0.0050
    (author spike-in coefficients; agrees on the dual gate and on direction)

  18.1 replicates      PASS        18.4 pooled dual gate     PASS
  18.2 source-faithful PASS        18.5 direction per rep    PASS
  18.3 power >= 0.80   FAIL 0.64   18.6 benchmark compat     PASS

  EXIT  ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE       STAGE 24  BLOCKED
```

### Claims this adds

**Allowed:**
- the corrected Role-A hypothesis is **supported as underpowered evidence** on two independent
  biological replicates not used to design the correction
- the direction is positive in both replicates and in both cross-replicate transfer directions
- the result is not an artifact of the author spike-in indexing bug

**Forbidden:**
- Role A is confirmed, or Stage 24 may open
- `+0.0305` is an effect estimate — an underpowered design selected on significance overstates
- any statement about which replicate carries the signal (the two arms disagree, on a one-clone
  difference)
- describing replicates 1/2/3 as sharing a uniform endpoint

### Evidence firewall — updated

```text
  consumed by diagnosis      GSM7092515, GSM7092516, ..._BC_gDNA.txt, the Stage-22 Rewind benchmark
  consumed by CONFIRMATION   GSM7092517, GSM7092518, GSM7092519  (+ the R2/R3 author outcome objects)
  still reserved             GSM7092520, GSM7092521  -- sorted, a DIFFERENT pre-state population,
                             disqualified from this claim rather than saved for it
```

Replicates 2 and 3 are no longer available to Stage 27. Stage 27 must preserve an independent
biological test on a system not used here.

---

## Historical record below — the pre-23.2H state, retained unedited

Stage-23.2 protocol `78edd5d7f9900349925339169a5d5e3e5011fe23e3c7c22608ac98bfe3427bf4`

## Role A

```text
historical Stage-23 verdict   ROLE_A_SIGNAL_FAIL   (permanent)
  observed delta_AP           +0.01050
  p_perm                      0.0846

diagnostic ledger
  MODEL_SELECTION_NULL_INFLATION             UNRESOLVED
  RESIDUAL_DEPTH_STRUCTURE                   SUPPORTED
  OUTCOME_LABEL_LIMITATION                   UNRESOLVED
  WITHIN_R1_EVENT_COUNT_LIMITATION           SUPPORTED
  BIOLOGICAL_REPLICATION_LIMITATION          SUPPORTED
  ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA   UNRESOLVED
```

### Corrected confirmatory hypothesis

**depth_complete_nuisance_control** — Under the historical Rewind outcome and the frozen Stage-22/23 evaluation geometry, pretreatment transcriptional state predicts the Role-A outcome beyond a DEPTH-COMPLETE nuisance baseline Bdepth = [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI), log1p(n_detected_GE_features_in_raw_pseudobulk)].

*Why this and nothing else:* RESIDUAL_DEPTH_STRUCTURE is the only SUPPORTED mechanism that implies a correction. MODEL_SELECTION_NULL_INFLATION is UNRESOLVED, so search matching is not indicated; OUTCOME_LABEL_LIMITATION is UNRESOLVED, so V2 §8.7 forbids adopting an alternative outcome. The correction was not chosen by effect size.

### Minimum design requirement for confirmation

```text
  scale 1:   35 positive clones -> power 0.290 at oracle AUC 0.66
  scale 2:   70 positive clones -> power 0.520 at oracle AUC 0.66
  scale 4:  140 positive clones -> power 0.940 at oracle AUC 0.66
  requirement: >= 140 positive clones (smallest TESTED cohort reaching 0.80)
```

the ladder is coarse (35 / 70 / 140 positive clones). The requirement is the smallest TESTED cohort reaching 0.80, not an interpolated value between rungs -- V2 §9.5 forbids extrapolating a precise required N.

### Claims

**Allowed:**
- Stage 23's Role-A permutation failure is partly explained by residual technical depth structure surviving the abundance-preserving permutation
- the historical experiment was underpowered for an AUC-0.66 signal at 35 positives
- a depth-complete nuisance control is the single mechanically indicated correction

**Forbidden:**
- Role A has a demonstrated signal
- the corrected same-data analysis confirms anything
- Rewind has no biological signal
- more clones from R1 would resolve the biological-replication limitation

### Evidence firewall

```text
consumed by diagnosis   GSM7092515, GSM7092516, stepThreeStarcodeShavedReads_BC_gDNA.txt, the full Stage-22 Rewind benchmark
reserved for confirmation   GSM7092517, GSM7092518, GSM7092519, GSM7092520, GSM7092521
```

Reserved accessions have been read at declared-metadata level only. No matrix has been downloaded, and no outcome or performance quantity computed for any of them.

## Role B — frozen, carried forward unchanged

```text
additive          ROLE_B_ADDITIVE_PASS
interaction       INTERACTION_PASS_MULTI_TREATMENT
C2 secondary      C2_INTERACTION_SECONDARY_CONFIRMED
baseline to beat  W5 on C1, pooled log loss 0.45465
```

Treatment-level limitations that must remain visible:
- Doxorubicin: W5 is WORSE than W4 on BOTH endpoints (-0.00332 C1, -0.00346 C2)
- Cisplatin C1 improvement is +0.00002, numerically negligible
- 'multi-treatment' means four treatments carry the interaction, not six
- captured abundance remains the dominant predictor (~3.45x the state contribution)

## Global

```text
feature universe   pretreatment Gene Expression only, 36,601 features; WM989's 153,055 Custom lineage features permanently excluded
split policy       clone-level outer folds frozen in Stage 22; never re-drawn
unresolved         MODEL_SELECTION_NULL_INFLATION, OUTCOME_LABEL_LIMITATION, ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA
```

### Exact Stage-24 opening rule

Stage 24 opens only on ROLE_A_CONFIRMATORY_SUPPORTED from 23.2G/23.2H -- a pre-registered confirmatory analysis succeeding on evidence not used to design the correction -- together with a complete, benchmark-compatible handoff. Role B positives cannot substitute for it without an explicit roadmap revision.
