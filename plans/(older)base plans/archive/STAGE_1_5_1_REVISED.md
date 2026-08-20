# STAGE 1.5.1 — **REVISED**: fix the ΔAge control, not the clock

**Supersedes** `STAGE_1_5_1_CLOCK_PRECISION.md` (V1), `STAGE_1_5_1_NEW_CHANGES.md` (the review) and
`STAGE_1_5_1_NEW_V2.md` (V2) as the execution plan. All three are left **byte-unmodified**.

**Status:** the diagnosis below is **EXECUTED** (2026-07-26) on real Gill data — every number is
re-runnable via `experiments/diag_gill_replication.py`. The steps in §4 are PLANNED.
`git diff --stat src/` is empty.

> **Why this document exists.** V1, the review and V2 all shared one unexamined premise: *the clock
> is too imprecise*. That premise is **false**. Tested directly, our current clock reproduces Gill
> et al.'s published ~30 yr rejuvenation to within 2 years. The failure was in **our analysis and our
> ΔAge label definition**, not in the instrument. V1/V2 would have spent ~4 h of GPU and a clock
> refit fixing something that was never broken — and, per V2 §0, could have shipped a *worse* clock
> while reporting success.

---

## 1. The finding

Gill's dataset labels every sample by outcome: **`Reprogramming fibroblast` (65)** — cells that are
successfully reprogramming — versus **`Failing to reprogram fibroblast` (47)**, plus 6 iPSC and 6
day-0 dermal fibroblasts. Applying our **current, unmodified clock**, contrasting each donor's day-0
baseline against Gill's reported optimum (days 10–13):

| arm | day 0 → peak |
|---|---|
| **Responders** (reprogramming) | **+8.2 yr** |
| **Non-responders** (failing to reprogram) | **+36.5 yr** |
| **Paired difference (responder − non-responder)** | **−28.3 yr**, 95% CI **[−51.6, −5.0]**, **5/6 donors negative** |

**Gill reported ~−30 yr. We measure −28.3 yr. That is a replication.**

Two things follow, and both are measurements rather than interpretations:

**(a) The identity artefact is real, large, and now quantified.** Cells that were exposed to OSKM and
**failed to reprogram** read **+36.5 years older after 11 days**. That is biologically impossible; it
is the clock reading loss of fibroblast identity as ageing. It is exactly the mechanism a reviewer
described (MET / somatic-identity collapse) and exactly what our own §9 H3 attribution had already
hinted at — its top drivers were **KRT7** (epithelial marker), **SULF1/TGFBI** (ECM), **COTL1**
(cytoskeleton), **GREM1/DKK1** (BMP/Wnt).

**(b) A contemporaneous non-responder control cancels it.** Both arms undergo the same OSKM exposure
and the same identity change; only one rejuvenates. Differencing them removes the artefact and leaves
the biology.

**Time course — it peaks where Gill says it does:**

| day | responder − non-responder gap | 95% CI |
|---|---|---|
| 7 | −9.7 | [−19.2, −0.2] |
| 9 | +30.0 | [−10.5, +70.5] ← noisy outlier |
| 11 | −23.8 | [−54.2, +6.5] |
| **13** | **−37.3** | **[−51.4, −23.1]** ← tightest, Gill's optimum |
| 15 | −40.6 | [−68.2, −12.9] |
| 21 | −23.1 | [−43.5, −2.6] |
| 29 | −22.3 | [−46.7, +2.1] |

---

## 2. What went wrong in our analysis — two errors, both ours

| # | Error | Consequence |
|---|---|---|
| **A1** | **Cell-type pooling.** `diag_e1_trajectory.py:179` excludes only `iPSC`. The **47 non-responders were averaged in with the 65 responders** — 42% of samples that by definition cannot rejuvenate | Diluted the effect toward zero. Pooled, the same contrast reads **+22.4 yr** instead of −28.3 |
| **A2** | **Wrong control in the label itself.** `sources.py:471` (Gill) and `:607` (HFF) define `is_control` as **day 0**, so `delta_age` measures against a baseline that has *not* undergone the identity change. The +36.5 yr artefact stays in every label | ΔAge is contaminated at source — this is upstream of every model, metric and stage |
| **A3** | **Wrong statistic for the shape.** E1/E1b used a **monotonic Spearman**; the effect is a **dip** that recovers after OSKM withdrawal (~day 13). A fall-then-rise rank-correlates to ≈0 | Produced "NO_TREND" and "WRONG_DIRECTION" from data containing a significant effect |

**Every escalation in the Stage 1.5 arc traces to these three.** M1's failure is separately real
(neonatal out-of-range, §6), but E1, E1b, D2 and the entire "ΔAge target is unvalidated" conclusion
were artefacts of our own method.

---

## 3. Why this fix helps — quantified, not asserted

The bars in V1/V2 were built on an effect size of **~11 yr** (`MASTER_PLAN` §5b-ter). That number was
measured **with the contaminated day-0 control** — the artefact was cancelling most of the signal.

| | effect | label noise (measured, T7) | SNR |
|---|---|---|---|
| current labels (day-0 control) | **11.4 yr** | 17.9 yr | **0.64** |
| corrected labels (non-responder control) | **28–37 yr** | 17.9 yr | **1.6–2.1** |
| corrected + condition-level averaging (§5b-ter, k≈10) | 28–37 yr | ~5.7 yr | **≫2** |

**The fix works by restoring the effect size, not by reducing noise** — which is why it succeeds
where the clock route could not. V2 §5 measured that even a *perfect* clock (mean MAE 4.0) reaches
only SNR 1.88 at k=1; this reaches 1.6–2.1 **with the clock we already have, at zero refit cost**,
and clears 2 comfortably once combined with the aggregation `MASTER_PLAN` §5b-ter already prescribes.

---

## 4. The plan, with the measurement for every step

**Discipline:** one change at a time; predictions pre-registered in the lab notebook before each
step; `FRAGILE` reported within 0.5 yr of any boundary; annotate-never-rewrite.

### Step 1 — Harden the finding *(read-only, hours)*

Already measured; three robustness checks remain, all cheap:

| Check | Bar |
|---|---|
| leave-one-donor-out on the paired statistic | separation stays negative for **all 6** leave-outs (guards against the N3 outlier carrying it) |
| alternative peak windows (11–13, 13–15, 10–15) | separation negative and within **±10 yr** of −28.3 across all windows |
| non-responders vs their *own* earlier timepoints | shows **no** rejuvenation — confirms they are a valid negative control, not merely a different population |

**How we measure the help:** this step adds no capability; it establishes whether the rest is worth
doing. **If the LOO check flips sign on any donor, STOP** — the effect is not robust at n=6.

### Step 2 — Define the control for **both** sources *(the real design work)*

⚠️ **This is the hard part, and it is where this plan can fail.** Gill is trivial; HFF is not.

| source | cells | responder/non-responder available? |
|---|---|---|
| **Gill** (bulk, 6 donors) | ~75 | ✅ **yes** — CD13/SSEA4 sorting is already parsed into `ctype` |
| **HFF / GSE242423** (single-cell) | **33,613 (79% of all age labels)** | ❌ **no** — unsorted |

Three options for HFF, to be decided by the pre-registered test below, not by preference:

- **H-a — marker-based split.** Score each cell for pluripotency/MET markers at its timepoint and
  split responders from non-responders computationally. Most faithful; needs its own validation.
- **H-b — mask HFF's ΔAge** (`age_mask=False`), keeping HFF for the fate head only and letting Gill
  carry the age signal. Simple and safe, but discards 79% of age labels.
- **H-c — leave HFF unchanged**, correct only Gill. Cheapest, but leaves 79% of labels contaminated
  and mixes two label definitions in one target — **likely the worst of the three**.

**Pre-registered decision test:** apply H-a's marker split to HFF and measure the same
responder-minus-non-responder separation. **If it reproduces a clearly negative separation, adopt
H-a. If it does not, adopt H-b** (mask), because a contaminated label is worse than a missing one.

**How we measure the help:** the separation statistic on HFF, with the same CI treatment as Gill.

### Step 3 — Implement *(small, surgical)*

Touches `sources.py` (control definition per source) and possibly `aging.py`. **No model, training,
calibration or inference code changes.**

**How we measure the help — a pipeline invariant that is unit-testable:** rebuild a single fold and
assert that the **built `y_age` reproduces the analysis**: Gill responders at day 13 must sit
**≈ −37 yr** relative to their matched controls, within the measured CI. If the artefact analysis and
the pipeline disagree, the implementation is wrong — this test catches it *before* any retrain.

### Step 4 — Rebuild, retrain, rescore *(~4 h GPU)*

Pre-registered predictions, so the result is interpretable either way:

| metric | prediction | why |
|---|---|---|
| `dage_mae_model` | **improves** | the labels carry a real signal instead of an artefact |
| `conformal_width` | **narrows** | label noise falls relative to effect |
| `res_approvals` | may become **non-zero** | R_eff = −(μ + z·σ) can finally be positive |
| four guards | **will move** | `y_age` changed — **by construction, not by defect** |

⚠️ **The `+0.000` guard streak ends here, and Stage 1's PARTIAL verdict does not carry over.** Stated
in advance so it is not misread as a regression.

### Step 5 — Revalidate with the corrected labels

Re-run the existing suite. **Specific pre-registered expectations:**

| test | expectation |
|---|---|
| E1b (reprogramming-phase trend) | the WRONG_DIRECTION verdict should resolve — it was measuring the artefact |
| E1 / D2 | re-run **split by responder status**, and with a contrast rather than a monotonic trend |
| **M1** | **should still fail.** It is out-of-range extrapolation (§6), a separate and unfixed issue |
| §9 H2 in-range tracking | should still pass — nothing here touches the clock |

---

## 5. What this plan does NOT fix — stated so it is not oversold

1. **The neonatal failure (R5) stands.** GSE113957 has **0 samples below age 1** (measured, T6), so
   N2, N3 and HFF remain out of the clock's fitted range. Their **absolute** ages are unusable. This
   is a *data* limit no control change can address. It does not block ΔAge, which is a difference.
2. **The compression stands.** Predicted-vs-true slope is **0.717**, so ΔAge magnitudes are ~28%
   understated. Rankings are unaffected; absolute claims should carry the caveat.
3. **The CI is wide** — [−51.6, −5.0] at n=6 donors. The effect is significant, not precise.
4. **The non-responder arm may itself partially rejuvenate**, which would make −28.3 a *conservative*
   underestimate. Acceptable direction, but it means the true effect could be larger.
5. **Deployment needs an internal control.** Scoring a new culture requires identifying
   non-responders within it — natural for single-cell, unsolved for bulk.

---

## 6. Disposition of the earlier plans

| document | status |
|---|---|
| `STAGE_1_5_1_CLOCK_PRECISION.md` (V1) | **PARKED.** Its premise (precision is the bottleneck) is refuted. Its bars were unreachable (V2 §4) and its lead candidate falsified (V2 §6) |
| `STAGE_1_5_1_NEW_CHANGES.md` (review) | its R4 refutation and C3 elimination stand and are correct; its P3 was wrong (V2 §1.3) |
| `STAGE_1_5_1_NEW_V2.md` (V2) | **PARKED**, but its measurements stand — especially §5 (a perfect clock still reaches only SNR 1.88), which is *why* the clock route was never going to be sufficient |

**Revisit the clock only if Step 4 shows the corrected labels are still too noisy** — and then for the
*compression*, which averaging cannot fix, not for precision.

---

## 7. Reproduction

```bash
python experiments/diag_gill_replication.py "D:\Gill"   # -> diag_gill_replication_results.json
```

| quantity | value |
|---|---|
| responders, day 0 → peak | +8.2 yr |
| non-responders, day 0 → peak | **+36.5 yr** (biologically impossible ⇒ artefact) |
| **paired separation** | **−28.3 yr, 95% CI [−51.6, −5.0], 5/6 negative** |
| pooled (as E1/E1b/D2 did) | +22.4 yr — the dilution that hid the effect |
| peak of the time course | **day 13–15** (−37.3 / −40.6) |
| Gill et al. published | ~−30 yr at days 10–13 |

**Source:** [Gill et al. 2022, *eLife* 71624](https://elifesciences.org/articles/71624) — including
their own report that *"existing transcription clocks failed to accurately predict the age of our
negative control samples"*, the same failure we hit.
