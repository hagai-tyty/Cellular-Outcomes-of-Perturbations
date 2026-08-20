# Stage 16 — the safety floor rejects cells that are demonstrably safe

**Status:** PLAN, then EXECUTE. Read-only diagnostic; `src/` is **not** touched by this stage.
A positive result buys the right to propose a Change, separately pre-registered.

> ### 🟩 STATUS AT CLOSE — EXECUTED, DIAGNOSIS CORRECTED, THEN VERIFIED 2026-08-18
>
> ⚠️ **§16.8's proposed fix was NOT the one shipped.** Pre-build diagnostics found that Platt was ALREADY applied at inference and that the real defect was the TARGET it was fitted against (soft `y_cls` vs the hard class every consumer reads). See `CHANGES.md` 2026-08-18 'STAGE 16 DIAGNOSIS CORRECTED'. The fix landed in `training/train_model.py`, and was **verified on recalibrated `_s16` folds**: sensitivity 0.275 -> 0.670 with specificity UNCHANGED at 0.929 -- the feared trade did not occur.
>
> *Commit: `e85111e / 9ff570c`. The status line above is the pre-execution claim and is left unedited, per the append-never-rewrite rule.*

---

## 16.0 The observation, and what it can and cannot buy

Stage 15 found that in five of six folds **more than half the held-out cells are
`REJECTED_UNSAFE`** (`S < τ_safe − 3w = 0.76`), while `fate_prauc` is 0.965–0.992 and
`fate_roc` 0.968–0.983. A head that ranks that well while its probabilities sit that low is the
signature of miscalibration, not of danger.

**What this CANNOT buy, stated first so nothing is oversold:** it cannot make RES non-zero. Stage
15 established `R_eff = 0` for **119 of 119** cells, and that gate is independent of the safety
gate. Even if every unsafe rejection were repaired, `g(R_eff) = 0` still zeroes the product.

**What it CAN buy:** whether the safety head is usable as a decision instrument at all. Fate
classification is the half of the project that is still alive, and `REJECTED_UNSAFE` is its
user-facing verdict. If that verdict is wrong on the majority of cells, the fate head's excellent
PR-AUC is not currently reaching any decision.

## 16.1 THE GATE — measured BEFORE this plan was written

A high rejection rate is not evidence of anything if the cells are genuinely unsafe. That had to
be checked first, with the power to kill the hypothesis outright. It did not.

| fold | n | TRUE safe | true safe % | rejected @0.76 | rej % | median S (true safe) | median S (true unsafe) |
|---|---|---|---|---|---|---|---|
| N2 | 19 | 19 | **100.0 %** | **15** | 78.9 % | 0.604 | — |
| N3 | 20 | 15 | 75.0 % | 12 | 60.0 % | 0.832 | 0.178 |
| O1 | 21 | 17 | 81.0 % | 19 | 90.5 % | 0.704 | 0.170 |
| O2 | 20 | 16 | 80.0 % | 19 | 95.0 % | 0.648 | 0.064 |
| Y1 | 18 | 8 | 44.4 % | 10 | 55.6 % | 0.785 | 0.621 |
| Y2 | 21 | 16 | 76.2 % | 15 | 71.4 % | 0.737 | 0.241 |

**N2 is decisive: 19 of 19 cells are truly safe and 15 are rejected as unsafe.** There is no
possible reading of that fold in which the rejections are correct.

The separation is excellent — median S is ~0.17 for truly-unsafe cells against ~0.70 for
truly-safe. The head *knows*. But the safe class sits near **0.70 instead of near 1.0**, and the
0.76 threshold cuts straight through the middle of it.

*Honesty note:* this gate was measured before the plan existed, so it is **not** pre-registered.
Everything in §16.3 was written before being run. (Same convention as Stage 14 §14.3.)

## 16.2 The candidate mechanisms

| | mechanism | status |
|---|---|---|
| **H1** | **Miscalibration** — probabilities systematically compressed toward the middle | open |
| **H2** | **The cells really are unsafe** | **KILLED by §16.1.** N2 is 100 % safe and 79 % rejected |
| **H3** | **Threshold placement** — 0.76 is simply the wrong bar for a well-ranked head | open |
| **H4** | **Class-prior shift** — the model trains on HFF (D0→iPSC, mixed) and is scored on a Gill donor that is 75–100 % safe. A posterior learned under one prior is biased under another | open |

H1, H3 and H4 are **different diseases with different cures**, and they are separable:

- If **H3** alone: the head is calibrated, and simply lowering the bar fixes it. ECE would already
  be low. It is not (0.249–0.266), so H3 alone is unlikely — but the *size* of the threshold error
  still needs measuring, because it bounds what any fix can achieve.
- If **H1**: a calibrator fitted on held-out data repairs it.
- If **H4**: the calibrator must be fitted on a cohort with the *target* prior, and one fitted on
  `calib` (which is HFF, not the held-out donor) will **not** transfer.

**H4 is the one that would be easy to miss**, and it is the reason the test below fits the
calibrator two ways rather than one.

## 16.3 The tests — PRE-REGISTERED

`experiments/diag_stage16_safety_floor.py`, read-only.

**Pooled, not per-fold.** n ≈ 20 cells per fold is too few for a per-fold rate, but every cell is
predicted by a model that never saw it, so pooling across folds is the more correct LOOCV estimate
— the identical argument `scorecard.pooled_fate_ece` already relies on. Pooled n = 119. Per-fold
numbers are reported alongside so a single-fold artefact cannot hide inside the pool.

**T1 — The cost, stated as a confusion at the shipped threshold.** Count FALSE REJECTIONS (truly
safe, `S < 0.76`) and false approvals (truly unsafe, `S ≥ 0.76`). This is the quantity the whole
stage is about; everything else explains it.

**T2 — Is S under-confident on the safe class?** Report the distribution of S among truly-safe
cells. A calibrated head should put them near 1.0. Report the fraction of truly-safe cells whose S
falls below the threshold.

**T3 — Does calibration repair it? Fitted TWO ways.**
  - **(a) deployable:** Platt fitted on `calib`, applied to `test`. This is what
    `scorecard._platt` already does and the only version that could ship.
  - **(b) oracle:** Platt fitted on `test` itself. **Not deployable, and labelled as such
    everywhere.** It is an upper bound: it says what calibration could achieve if the transfer
    problem did not exist.
  - The **gap between (a) and (b) is the measurement of H4.**

**T4 — Where is the optimal threshold?** Sweep the threshold and find the one maximising balanced
accuracy on the pooled set (oracle). Report the distance from 0.76. This bounds H3.

**T5 — Prior shift, measured directly.** Safe fraction in `train`, in `calib`, and in `test`. If
the train/calib prior differs materially from the test prior, H4 has a mechanism and not just a
symptom.

## 16.4 PRE-REGISTERED READING

| T3(a) deployable Platt | T3(b) oracle Platt | conclusion |
|---|---|---|
| false rejections **fall materially** | falls | **H1 — plain miscalibration.** A deployable calibrator fixes it; propose it as a Change |
| **does not fall** | **falls** | **H4 — prior/cohort shift.** Calibration is possible but does not transfer from `calib` to a held-out donor. The fix is a prior correction or a donor-matched calibrator, NOT Platt-on-calib |
| does not fall | does not fall | **H3 or something structural.** Calibration is not the lever; report the optimal threshold from T4 and treat the bar as a policy question |

**Bar for "materially":** false rejections must fall by **more than half** on the pooled set
(≥ 50 % reduction), because the current rate is so high that a marginal improvement changes no
decision. Stated before the run.

**A guard against the obvious trap:** a calibrator that simply pushes every probability up would
"fix" false rejections by destroying the false-approval rate. **Both** directions are reported, and
any drop in false rejections must not be bought with a rise in false approvals beyond the number of
truly-unsafe cells available. Balanced accuracy is reported so the trade is visible in one number.

## 16.5 What this stage will NOT do

- **No `src/` change.** Diagnosis only.
- **No claim that RES becomes non-zero.** §16.0. The `R_eff` gate is independent and still closed.
- **No re-tuning of `τ_safe` to make the numbers look better.** If T4 says the threshold is
  misplaced, that is *reported*, and any change to a shipped safety threshold is a separate,
  pre-registered decision with its own bar — a safety floor is a policy, not a hyperparameter to
  fit.
- **No use of the oracle calibrator as a result.** It is an upper bound and is labelled as one
  everywhere it appears.

## 16.6 Verification

| item | how |
|---|---|
| the gate can kill | the true-safe fraction is computed and reported per fold, not assumed |
| no leak in (a) | the deployable Platt sees only `calib`; a test asserts the fit set and the score set are disjoint |
| the oracle is marked | a test asserts every oracle field carries `oracle` in its name |
| the trap is closed | false approvals reported alongside false rejections; balanced accuracy reported |
| pooling is honest | per-fold numbers reported beside the pooled ones |
| record | `CHANGES.md` with the verdict, including if it is negative |
| suite | full `pytest` green + `ruff` clean |

---

# ANNOTATION — added 2026-08-17, AFTER the run

*Everything above is the plan as written before execution and is left unedited.*

## 16.7 RESULT — verdict **H1**, and the deployable calibrator clears the bar

Pooled n = 119 (91 truly safe, 28 truly unsafe). **70.3 % of truly-safe cells fall below the
0.76 bar.** Median S: true-safe **0.704**, true-unsafe **0.217** — the head separates cleanly; its
safe class simply sits under the bar.

| arm | false rej | drop | false appr | sens | spec | bal acc |
|---|---|---|---|---|---|---|
| raw (as shipped) | **64** | — | 2 | **0.297** | 0.929 | 0.613 |
| Platt on `calib` **[DEPLOYABLE]** | **26** | **−59.4 %** | 5 | **0.714** | 0.821 | 0.768 |
| Platt on `test` *[ORACLE]* | 27 | −57.8 % | 0 | 0.703 | 1.000 | 0.852 |

**The headline: shipped sensitivity is 0.297.** The safety gate approves under a third of
genuinely safe cells, while the head behind it scores PR-AUC 0.965–0.992. The head's quality is
not reaching the decision.

**H4 is refuted.** The deployable arm matches the oracle on false rejections (26 vs 27), so
calibration *does* transfer from `calib` to a held-out donor. Better: the deployable calibrator is
fitted on all six folds, whereas the oracle is **undefined on N2** (19/19 safe → unidentifiable
boundary → passes through). The deployable arm is the stronger one in practice.

**The trap is closed.** False approvals rise only 2 → 5 while false rejections fall 64 → 26, and
balanced accuracy rises 0.613 → 0.768. This is a favourable trade, not a shifted operating point.

**H3 also contributes.** The oracle-best threshold on *raw* scores is **0.495** against the shipped
**0.76** (balanced accuracy 0.832). So H1 and H3 are two views of one fact: the probabilities are
compressed low. Recalibrating them is the principled repair; lowering the safety bar would be
fitting a safety policy to data, and §16.5 rules it out.

**Post-hoc, not pre-registered:** the false-rejection rate tracks how far a fold's test prior sits
from the training prior (train 53.2 % safe; test 44.4–100 %). Spearman **+0.771**, n = 6, against a
two-sided critical ρ of 0.886 — **suggestive, not significant.** Y1, the fold whose prior is
closest to training (44.4 %), has the lowest false-rejection rate (2 of 8); N2, the furthest
(100 %), has among the highest (15 of 19). Recorded as a mechanism candidate, not a finding.

## 16.8 PRE-REGISTRATION of the Change this licenses

**Proposed:** apply the Platt calibrator — fitted on `calib`, as `scorecard._platt` already does —
to `S` **before** the safety gate in `compute_res_batch` / `compute_res`, behind a flag that is
**off by default** until the Change is accepted.

**Target metric:** pooled safety **sensitivity**, bar ≥ 0.60 (measured 0.714; bar set below the
observed value deliberately, so the Change is judged on clearing a threshold rather than on
reproducing one run).

**Guards, any one of which blocks acceptance:**

| guard | bar |
|---|---|
| pooled false approvals | must not exceed **7** (measured 5; raw is 2) |
| pooled balanced accuracy | must **rise** vs raw |
| `fate_prauc`, `fate_roc` | **unchanged** — Platt is monotone, so ranking is invariant by construction. Any movement means the implementation is wrong, not that it helped |
| RES | must stay **0** — the `R_eff` gate is independent and still closed. A non-zero RES here would mean something else moved |

**Explicitly NOT part of this Change:** `τ_safe` and `w` are untouched. The safety floor is a
policy; this Change fixes the probabilities fed to it, not the floor itself.

**Not to be done:** accepting on the oracle arm, or re-fitting the calibrator on test. The oracle
is an upper bound and is undefined on N2.
