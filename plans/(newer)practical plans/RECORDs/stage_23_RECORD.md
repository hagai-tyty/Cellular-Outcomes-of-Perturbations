# stage_23_ RECORD

Authoritative record for Stage 23. Sections are **additive and immutable**: 23A is written now,
23B–23F are appended as they execute. An earlier section is never rewritten because a later result
is inconvenient; corrections are appended with a dated note.

---

# 23A — Protocol / representation freeze

## Goal
Run the pre-registered §3.0 input audit against the frozen Stage-22 benchmark and the external raw
matrices, build one fixed clone-level `X_before`, and freeze every evaluation choice before any
outcome model exists. **No model is fitted in 23A.**

## Inputs
- plan: `STAGE_23_LEARNABILITY_INTERACTION_GATE_V2.md` (V1 archived under `arcive/` after the
  independent pre-execution audit)
- frozen Stage-22 benchmark commit `8d6011a`; record `6b6a169`; correction `9eaddf2`
- `D:\GSE227151_Rewind\` (11 source + 2 author-code files) · `D:\GSE279162\` (29 + 5)
- model: `_s16` frozen, never loaded

## Files added
- `experiments/run_stage23_learnability_gate.py`
- `tests/test_stage23_learnability_gate.py` (23 tests)
- `results/stage23_protocol.json`
- `results/stage23_rewind_clone_expression_manifest.json`
- `results/stage23_wm989_clone_expression_manifest.json`
- `results/stage23_outer_fold_preprocessing.json`

## Files modified
- `tests/test_ci_portability.py` — the new Stage-23 test module is registered with the
  CI-condition guard so it is re-run with every local dataset root absent

## What changed
- The §3.0 audit is **executable code, not prose**: 38 checks, and pseudobulk construction is
  gated behind `STAGE22_INPUTS_AUDITED`
- `X_before` exists as a frozen artefact with a content hash; the matrix itself stays in the
  gitignored `_cc_cache/` and never enters git
- Every later-substage choice — K, grids, seeds, treatment coding, tie-break, metrics, bootstrap
  and permutation rules — is written into `stage23_protocol.json` now

## What did NOT change
- `src/` unchanged · no model fitted · no Stage-22 artifact, fold, label or manifest rewritten
- the frozen Stage-22 CRLF/LF provenance mismatch is **declared, not repaired** (V2 §1.4)
- the frozen WM989 target is **not** rebuilt from a treated-only alternative (V2 §1.2.1)

## Tests
- 1857 passed · ruff clean (CI scope)
- 121 passed under `CELLFATE_NO_LOCAL_DATA=1` (the CI condition) for the Stage-23 + portability
  modules; the pseudobulk cache is absent there and those tests skip rather than fail

## Result

**VERDICT: `PROTOCOL_FROZEN`** — input audit `STAGE22_INPUTS_AUDITED` (38/38 checks).

### Input audit

```text
six benchmark CSV hashes match their manifests          OK
two manifest hashes match stage22 results               OK
results file omits its own hash                         OK
Role-A verdict ready                                    BENCHMARK_READY_WITH_DECLARED_MISSINGNESS
all_gates_pass / every individual gate / model_fitted    true / 10-of-10 true / false
preflight derived independently of the `overall` string  OK
40 external source + 7 author-code files                present, byte-identical
every expression_column_index resolves                  0 mismatches over 10,394 cells
36,601 Gene Expression features per sample              OK   (Custom block 153,055 where present)
all 11 samples share ONE gene feature-ID list           OK
WM989 targets rebuilt from treated cells only           8,406 rows, max abs diff 0
canonical text hash is LF/CRLF invariant                OK
```

Two inherited limitations are recorded rather than silently fixed:

- **Stage-22 gate derivation.** Stage 22 derives `overall` from the Role-A verdict without
  consuming `all_gates_pass`. Confirmed harmless here because all ten gates are true — and Stage 23
  fails closed on its own derived condition regardless of the serialized string.
- **CRLF/LF provenance.** The Stage-22 plan and builder digests are checkout-byte hashes and differ
  between Windows and Linux for identical content. Stage 23 uses a canonical LF-normalised digest
  for text; a test asserts LF, CRLF and CR inputs all produce one hash.

### Clone `X_before`

```text
                 clones   genes    nnz          density   cells/clone (min/med/max)
GSE227151         3,147   36,601   15,880,583   0.1379    1 / 1 / 11
GSE279162         1,401   36,601   10,697,195   0.2086    1 / 2 / 88
```

Built by summing **raw** pretreatment counts per clone, then CP10K and `log1p` **exactly once**.
Verified by inverting the transform: `expm1` row sums are 10,000.000 for every clone in both
matrices, and one 11-cell Rewind clone was recomputed from the raw `.mtx` independently, agreeing
to `4.4e-16`.

Only pretreatment lanes were read — the two Rewind `*_control_*` matrices and the three WM989
`*_Naive*` matrices. No treated column was opened, and `X` is exactly 36,601 columns wide, so none
of WM989's 153,055 `Custom` lineage features can be present.

### Outer-training gene filter (descriptive)

```text
              fold 0   fold 1   fold 2   fold 3   fold 4    detection floor
GSE227151     13,610   13,589   13,627   13,606   13,607    26
GSE279162     18,385   18,411   18,270   18,227   18,361    12
```

`max_feasible_K >= 50` in every fold, so **no candidate K may be skipped**. This table is
descriptive only: per V2 §2.5 the filter is refitted inside every inner-training split.

## Bugs found
1. **The audit surfaced nothing new** — every V2 §1 anchor reproduced. That is the expected
   outcome given the pre-execution audit, and is recorded as a pass rather than as an absence
2. A first draft of the module wrote the pseudobulk into `results/`. Corrected before commit: a
   3,147 × 36,601 matrix is not a compact benchmark artefact, and V2 §13 forbids it. It now lives
   in the gitignored `_cc_cache/`, with its content hash committed so a rebuild stays verifiable

## Scientific interpretation
**Proves:** the Stage-22 benchmark is intact and externally consistent with the raw data at the
moment Stage 23 opened, and one fixed clone-level representation now exists whose construction is
independently checkable from the committed hashes. Every evaluation choice that could otherwise be
made after seeing a result is frozen.

**Does NOT prove:** anything about learnability. No model has been fitted, no outcome has been
exposed to an estimator, and no comparison has been run. The observed-zero/abundance confound
recorded in Stage 22 is untouched by this substage and remains the load-bearing risk for Role B —
which is why V2 §1.2.1 added `log1p(n_naive_cells)` to the nuisance block before any fitting.

## Next action
23B — Rewind Role-A learnability. Not started.

---

# 23B — Rewind Role-A learnability

## Goal
Ask whether pretreatment transcriptomic state predicts the later author-defined Rewind priming
outcome **beyond prevalence and captured clone size**, under the protocol frozen in 23A. This is
the first substage that fits an estimator.

## Inputs
- `results/stage23_protocol.json` (23A, `PROTOCOL_FROZEN`) — referenced by SHA-256 in the results
- clone `X_before` from the 23A cache, re-verified against its committed content hash before use
- frozen Stage-22 outer folds, `y_primed`, and the two nuisance columns
- no new seed, no new split, no new grid

## Files added
- `results/stage23_rewind_oof_predictions.csv` (3,147 rows — one frozen OOF row per clone)
- `results/stage23_rewind_results.json`

## Files modified
- `experiments/run_stage23_learnability_gate.py` — `--stage 23b` added; 23A untouched
- `tests/test_stage23_learnability_gate.py` — 13 further contracts (36 total)

## What changed
- Four models fitted exactly as pre-registered: `R0` prevalence, `R1` nuisance-only, `R2` X-only,
  `R3` X + nuisance
- Every learned quantity is refitted inside each inner-training split; the outer test fold is
  transformed and predicted exactly once
- Hyperparameters selected by mean inner Average Precision, with the frozen tie-break

## What did NOT change
- `src/` unchanged · outer folds identical to Stage 22 · no Stage-21/22 artifact rewritten
- 23A's protocol, manifests and cached `X` untouched
- no threshold tuning, no class reweighting, no SMOTE, no model outside the frozen grids

## Tests
- 1873 passed · ruff clean (CI scope)
- **0 convergence warnings** across every candidate fit — V2 §3.7 treats one as a protocol failure
  to investigate, not as permission to drop a candidate

## Result

**PROVISIONAL VERDICT: `ROLE_A_SIGNAL_PASS`** — provisional until 23E structural controls and
`ROLE_A_PERMUTATION_PASS`.

### FROZEN 23B HEADLINE — accepted 2026-08-22, not to be re-run or re-tuned

The 23B protocol, models, folds, preprocessing, hyperparameters, metrics and thresholds are
**frozen exactly as reported**. Six facts carry forward together; quoting the verdict without them
would misrepresent the result.

```text
1. dAP_state = AP(R3) - AP(R1) = +0.01050   95% CI [+0.00397, +0.02258]
2. fold-wise R3-R1 positive in 5/5 folds    DIAGNOSTIC ONLY, not a gate
3. absolute signal is WEAK                  AP(R3) = 0.02085 , ROC-AUC = 0.6628
4. R3 does NOT improve log loss or Brier over the prevalence model R0
       log loss  R3 0.06209  vs  R0 0.06109      (R3 worse)
       Brier     R3 0.011078 vs  R0 0.010998     (R3 worse)
5. R1 is WORSE than R0 (AP 0.01035 vs 0.01112, ROC-AUC 0.4747, below chance),
   so the Rewind nuisance baseline is poor and dAP_state is partly "R1 is weak"
6. FINAL Role-A PASS still requires the 23E permutation gate
       ROLE_A_PERMUTATION_PASS + STRUCTURAL_CONTROLS_PASS
```

Point 4 deserves emphasis: the state signal shows up as **ranking** improvement (AP, ROC-AUC), not
as better-calibrated probabilities. A model that ranks slightly better while scoring marginally
worse on log loss and Brier than a constant is a weak-signal signature, not a usable classifier.

### Pooled out-of-fold metrics (3,147 clones, 35 positive)

```text
model                      AP     ROC-AUC    log loss     Brier
R0  prevalence        0.01112      0.4998     0.06109  0.010998
R1  nuisance only     0.01035      0.4747     0.10181  0.015630
R2  X only            0.01923      0.6043     0.07493  0.012592
R3  X + nuisance      0.02085      0.6628     0.06209  0.011078
```

### Selected hyperparameters per outer fold

```text
fold    R1 C     R2 K/C     R3 K/C    inner AP  R1       R2       R3
0       0.01     10/0.01    50/10              0.01167  0.03094  0.03015
1       0.01     10/0.1     10/0.1             0.01065  0.03270  0.03216
2       0.01     10/0.1     10/0.1             0.01067  0.03729  0.03024
3       1        20/10      10/0.1             0.01093  0.04895  0.04863
4       1        10/1       10/1               0.01080  0.05441  0.05358
```

`K = 10` is selected in 8 of 10 X-model fits — the tie-break prefers the smaller K, and the larger
bases did not earn their place.

### Primary inference

```text
ΔAP_state  = AP(R3) - AP(R1) = +0.01050
             95% CI [+0.00397, +0.02258]     P(Δ <= 0) = 0.0015
ΔAP_abs    = AP(R2) - AP(R0) = +0.00812
             95% CI [+0.00191, +0.04393]     P(Δ <= 0) = 0.0010
```

Stratified clone bootstrap, 2,000 replicates, seed 23123, positive and negative clones resampled
separately so class counts are preserved.

### Fold diagnostics (reported, not a gate)

```text
fold      R0       R1       R2       R3    ΔAP R3-R1   ΔAP R2-R0   pos   genes
0    0.01111  0.00998  0.03890  0.02487     +0.01489    +0.02779     7  13,610
1    0.01111  0.01061  0.03100  0.02820     +0.01759    +0.01989     7  13,589
2    0.01113  0.01050  0.02523  0.02810     +0.01760    +0.01410     7  13,627
3    0.01113  0.01083  0.02388  0.05479     +0.04396    +0.01275     7  13,606
4    0.01113  0.01166  0.02443  0.02471     +0.01304    +0.01330     7  13,607
```

Both deltas are positive in 5/5 folds. Per V2 §4.5 this is a high-variance diagnostic on seven
positive clones per fold and creates no separate PASS requirement — it is consistent with the
pooled result rather than independent evidence for it.

## Bugs found
1. Three defects in my own new tests, all caught before commit: the "23A fits nothing" contract was
   written against the whole module and broke once 23B legitimately fitted estimators (now scoped
   to the 23A functions by AST); a tolerance compared a 6-dp rounded field at `1e-9`; and two
   assertions ended in `or True`, which made them unconditionally true. The `or True` pair is the
   worst of the three — a test that cannot fail is worse than no test, because it reads as coverage

## Scientific interpretation
**Proves:** under the frozen geometry, pretreatment expression carries held-out prospective
information about the later priming outcome that the captured-clone-size baseline does not. The
lower bound of the pre-registered interval clears zero on both the incremental comparison
(`R3` vs `R1`) and the absolute one (`R2` vs `R0`), and the direction is consistent across all five
folds.

**Three things this does not say, which matter more than the verdict:**

- **The absolute signal is very weak.** `AP = 0.021` against a prevalence of `0.011` — a model
  roughly twice as good as guessing, on 35 positives. ROC-AUC `0.663`. This clears a learnability
  gate; it is nowhere near a usable predictor, and the Stage-24 scope should be set by the effect
  size, not by the word PASS.
- **`R1` is worse than `R0`.** Nuisance-only AP `0.01035` is *below* prevalence and its ROC-AUC is
  `0.4747` — below chance. Captured clone size is not merely uninformative for Rewind priming, it
  is faintly anti-predictive out-of-fold. So `ΔAP_state` is partly "R1 is a poor baseline" rather
  than purely "X is good". The descriptive arithmetic `AP(R3) - AP(R0) = +0.0097` is of similar
  magnitude, so the conclusion does not appear to be an artifact of a weak `R1` — but that
  comparison was **not pre-registered**, carries no interval, and must not be promoted into
  evidence. It is recorded so the asymmetry is visible.
  This is the opposite of the WM989 picture, where captured depth is a strong competitor.
- **One biological replicate.** Rewind R1 remains a single experiment, so this is within-R1
  clone-held-out generalisation and nothing wider.

The verdict is provisional. It is not final until 23E shows the structural controls pass and the
observed `ΔAP_state` beats its permutation null.

## Next action
23C — WM989 additive state-signal gate. Not started.

---

# 23C — WM989 additive state-signal gate

## Goal
Ask whether pretreatment `X` adds beyond treatment `U` and captured naive clone abundance `B`, on
two descriptive future endpoints. **No interaction terms** — those are 23D.

## Inputs
- 23A protocol (`PROTOCOL_FROZEN`), referenced by SHA-256; clone `X` re-verified against its
  committed content hash before use
- frozen Stage-22 WM989 outer folds, `n_post_cells`, and the four nuisance columns
- `B = log1p(n_naive_cells) + log1p(n_naive1_cells) + log1p(n_naive2_cells) + log1p(n_naive3_cells)`
- treatment coding frozen: canonical order `Acid, Cisplatin, CoCl2, Dabrafenib, Doxorubicin,
  Trametinib`, reference `Acid`, five non-reference dummies, never standardized

## Files added
- `results/stage23_wm989_detection_oof.csv` (8,406 rows) ·
  `results/stage23_wm989_abundance_oof.csv` (2,256 rows) ·
  `results/stage23_wm989_results.json`

## Files modified
- `experiments/run_stage23_learnability_gate.py` — `--stage 23c`; 23A/23B untouched
- `tests/test_stage23_learnability_gate.py` — 15 further contracts (51 total)

## Tests
- 1887 passed · ruff clean (CI scope) · **0 convergence warnings**

## Result

**VERDICT: `ROLE_B_ADDITIVE_PASS`** — provisional until 23E structural controls and
`ROLE_B_ADDITIVE_PERMUTATION_PASS`.

The PASS rests entirely on **C1 detection**. C2 conditional abundance is not significant, and is
not significantly worse — which is exactly the pre-registered condition for the other endpoint.

### C1 — future clone detection (8,406 rows, 1,401 clones)

```text
model                  log loss        AP     ROC-AUC      Brier
W0  U                   0.55992   0.36895     0.60973    0.18716
W1  B + U               0.47832   0.61365     0.76410    0.15356
W2  X                   0.54772   0.40072     0.68358    0.18290
W3  X + U               0.52408   0.49944     0.71531    0.17226
W4  X + B + U           0.46895   0.63190     0.77983    0.14998

ΔLL_state = logloss(W1) - logloss(W4) = +0.00937
    95%   CI [+0.00583, +0.01259]
    97.5% CI [+0.00516, +0.01315]      P(Δ ≤ 0) = 0.0000
```

### C2 — conditional future abundance (2,256 nonzero rows, 929 clones)

```text
model      clone-balanced MAE   clone-balanced RMSE   mean treatment Spearman
W0  U                 0.58423               0.82872                  -0.10189
W1  B + U             0.51449               0.76305                  +0.36283
W2  X                 0.59963               0.83582                  +0.12475
W3  X + U             0.58048               0.82046                  +0.10850
W4  X + B + U         0.51008               0.75131                  +0.39987

ΔMAE_state = MAE(W1) - MAE(W4) = +0.00441
    95%   CI [-0.00184, +0.01048]
    97.5% CI [-0.00247, +0.01111]      P(Δ ≤ 0) = 0.0830
```

Per-treatment Spearman on pooled OOF (secondary):

```text
model      Acid  Cisplatin   CoCl2  Dabrafenib  Doxorubicin  Trametinib
W0      -0.0465    -0.1172 -0.0829     -0.1013      -0.0917     -0.1717
W1      +0.6162    +0.4658 +0.3491     +0.1902      +0.4150     +0.1407
W4      +0.5892    +0.4761 +0.3469     +0.2901      +0.4235     +0.2734
```

### Selected hyperparameters per outer fold (`K / C` for C1, `K / alpha` for C2)

```text
C1   fold      W0        W1        W2        W3        W4
     0      None/10   None/10     50/1      50/1      50/1
     1      None/10    None/1     50/1      50/1      50/1
     2      None/10   None/10    50/10     50/10      50/1
     3      None/10    None/1    50/10      50/1      50/1
     4      None/10   None/10     50/1      50/1      50/1

C2   fold      W0        W1        W2        W3        W4
     0     None/0.1    None/1   10/0.1    10/0.1     20/10
     1       None/1   None/10   10/100    10/0.1     50/10
     2     None/0.1    None/1   10/0.1    10/0.1      20/1
     3     None/0.1   None/10   10/0.1    10/0.1     20/10
     4     None/0.1   None/10   10/0.1    10/0.1     50/10
```

Every C1 X-model selected `K = 50`; C2 selected `K = 10` for the X-only models and 20–50 once `B`
was present. Retained outer-training genes: 18,227–18,411 (C1), 19,026–19,256 (C2).

## Bugs found
1. **A duplicate test name.** The 23C hyperparameter contract was written with the same function
   name as the 23B one, so Python would have silently replaced the 23B test and it would never
   have run again. Ruff's `F811` caught it; renamed, and a check confirms all 51 test names are
   now unique. A shadowed test is worse than a missing one — it still shows up as coverage
2. No deviations from the frozen protocol, no failed assumptions, and no convergence warnings

## Deviations / clarifications
- V2 §3.6 restricts the C2 dataset to nonzero rows **after** the outer split is known. It does not
  spell out which clone set the C2 gene filter / scaler / PCA is fitted on. Implemented as: fitted
  on the **C2-eligible outer-training clones** — the actual training population — never on clones
  outside the analysis set. Recorded here because it is a reading of the plan rather than a
  quotation of it
- Expression transforms are fitted on **unique clones**, not on repeated clone×treatment rows, so a
  clone observed under six treatments does not get six times the weight in the PCA basis

## Scientific interpretation
**Proves:** on the detection endpoint, pretreatment state adds a small but reliably positive amount
beyond treatment identity and captured naive abundance. The lower 97.5% bound clears zero
(`+0.00516`) and no bootstrap replicate out of 2,000 put Δ at or below zero.

**The result that matters more than the verdict: `B` is doing most of the work.**

```text
C1 log-loss gain from adding B to U      0.55992 -> 0.47832   = 0.0816
C1 log-loss gain from adding X to B+U    0.47832 -> 0.46895   = 0.0094
```

Captured clone abundance buys roughly **nine times** what pretreatment expression adds on top of
it. This is precisely the confound Stage 22 measured (observed-zero rate 87.1% → 32.6% across
capture depth), and it is why V2 §1.2.1 made the total-depth term mandatory before any fitting. A
Role-B PASS here means "X adds a little beyond a strong abundance baseline", not "X predicts
treatment response".

**Does NOT prove:**

- **Anything on conditional abundance.** `ΔMAE_state` is `+0.0044` with a 97.5% interval spanning
  zero and `P(Δ ≤ 0) = 0.083`. Under the pre-registered rule this endpoint contributes only the
  "not significantly worse" half of the PASS. It must not be reported as a positive finding.
- **That `X` is useful alone.** `W2` (X only) is *worse* than `W0` (treatment only) on C2 MAE
  (`0.59963` vs `0.58423`), and its mean treatment Spearman is `+0.125` against `W1`'s `+0.363`.
  Pretreatment expression without the abundance baseline is a poor predictor here.
- **That the interaction question is answered.** No `X × U` term was fitted. The C2 per-treatment
  Spearman table shows `W4` improving over `W1` mainly on Dabrafenib (`+0.190 → +0.290`) and
  Trametinib (`+0.141 → +0.273`) while being flat or slightly lower elsewhere — a *hint* that the
  state contribution may be treatment-dependent, which is exactly what 23D exists to test properly.
  It is recorded as a secondary observation and carries no inferential weight here.

The verdict is provisional until 23E.

## Next action
23D — WM989 explicit interaction gate. Not started.
