# THE PATH — what is settled, what is open, and the one decision left

**Written 2026-08-16 from measured results only.** Every number here is recorded in `CHANGES.md`
with the script that produced it. Nothing below is a proposal to re-measure.

---

## 1. Settled. Do not spend anything more on these.

| question | answer | evidence |
|---|---|---|
| Can we predict same-timepoint ΔAge from expression? | **The question is void — it is circular** | 1 956 of the 2 000 panel genes carry clock weights; the clock's own weights reconstruct the label at **ρ 0.96–0.97** |
| Can early ΔAge predict late ΔAge? | **No** | partial **−0.064** after donor age; donor age alone gives **0.971** and is free |
| Can early *expression* predict the late residual? | **No** — fragile | **3 of 9** configurations, effectively 1 of 9 |
| Is there a forward Δt signal (Stage 3a)? | **Unanswerable at that geometry** | regime A **UNRESOLVABLE at every alpha**, raw and logit |
| Would more per-timepoint replication fix it? | **No** | regime E **P0-void** — `p_unsafe` is a fraction of *cells*; a bulk sample is already an average |
| Did C-7 fix the labels? | **Yes, largely** | fold spread **16.67 → 3.69 yr**; magnitudes onto §5.14's predicted −8.196 |
| Can ΔAge be verified to 10 %? | **No — unverifiable, not merely hard** | 10 % of a 12.66 yr SD is **1.36 yr**; the two references disagree by **7.30 yr MAE** |

---

## 2. 🔑 The decision nobody has framed correctly

§5.13 posed it as **"adopt `top100`, or don't."** That is not the choice. **Something has to ship,
and the incumbent is `raw`.** The real choice is **which of two instruments ships**, and on the
like-for-like table `top100` **dominates `raw` on every axis measured**:

| | vs multi-tissue | vs skin & blood | beats predicting **nothing**? |
|---|---:|---:|---|
| **`raw` — what ships today** | **22.69** | **25.49** | ❌ **NO, on both** — the floor is 11.71 / 9.89 |
| **`top100`** | **7.15** ✅ at the instrument floor | **11.27** | ✅ on multi-tissue |

> **The clock currently in production is roughly twice as bad as a predictor that outputs zero, on
> both reference clocks.** `top100` is better than it by **15.5 yr** and **14.2 yr** respectively.

**Both fail the SPLIT rule, so SPLIT does not separate them.** It rejects both equally. The rule was
written to stop a favourable reference being cherry-picked — it was never meant to decide *which
failing instrument to keep*, and applied that way it silently keeps the worse one.

**This is not narrowing the estimand.** 1.5.2 refuted that and this analysis independently agrees
with 1.5.2. The claim stays "agreement with both references". `top100` still does not meet it. The
question is only whether the project ships the instrument that misses by 7.15/11.27 or the one that
misses by 22.69/25.49.

**That decision is the user's.** It is the only genuine judgement call left in the ΔAge arc.

---

## 3. The work, in order

### Now — free, no retrain

Nothing. **ΔAge accuracy is fully measured.** Any further same-timepoint modelling buys panel
fidelity (§1), and every variant in the ledger family has been scored against the floor.

### Next — ONE retrain, and it answers four open questions at once

C-7 §5 already requires a full guard re-report on adoption, and §5.13 freed that budget by
cancelling step 4. **One run settles all of:**

1. **RES** — currently **unknown**, not "unchanged". The newest scorecard predates C-7
2. **Every Stage 1 guard under clean labels**, over 6 folds
3. **Whether the clock choice moves anything downstream** — run it with whichever clock §2 decides
4. **Whether C-7's label fix changes the fate head** — expected: no, it consumes no ΔAge

**Pre-register before it runs:** C-7 flag ON, clock frozen by §2's decision, all guards reported
before/after, snapshot and rollback. **Do not bundle any other change into it.**

### After — and only after

| | |
|---|---|
| **Integrate GSE165177** | 95 samples, on disk, in no training config, **replicated contemporaneous controls** — the one dataset where `n = 1` zero-points cannot occur. Its own retrain, own snapshot |
| **C-2** | **after Stage 6 only.** Enabling it masks ~99.8 % of the corpus |
| **Stage 6 acquisition** | the ask is now **specific**: a single-cell time course with several donors. Bulk cannot express `p_unsafe` — regime E proved that, so "more bulk" is not the requirement |

---

## 4. What the project can claim today, without qualification

| | status |
|---|---|
| **Fate classification** | ✅ `fate_roc` **0.983**, `fate_prauc` 0.992 — untouched by every ΔAge problem |
| **Within-donor ranking** | ✅ Spearman **0.925–0.983**, every fold |
| **ΔAge as a measurement** | ✅ with `top100`: **at the instrument floor on multi-tissue** (7.15 vs 7.30, CI spans 0), spread preserved (SD ratio 0.98), ordering **ρ 0.810** — *better than the two gold standards order each other* (0.613) |
| **ΔAge as a prediction target** | ❌ circular at ρ 0.96 |
| **Forward / stopping-time tool** | ❌ not supported by this data |

**The methodological findings stand on their own and are publishable independently:** a production
transcriptomic clock that **loses to predicting nothing on two independent references**; a **−14.10
yr density bias removable by sparsification**; a **degenerate GEO column** that inflated 99.7 % of a
corpus threefold through a variance-floor path; and an instrument floor that **bounds what any RNA
readout can be validated to.**

---

## 5. The one thing that would change the picture

**A third independent reference.** Every ceiling here comes from having exactly two methylation
clocks that agree with each other at **ρ 0.613 / MAE 7.30 yr**. With a third, the shared component
becomes identifiable and "agreement with both" stops being a bar whose difficulty is set by the
references' mutual disagreement rather than by the RNA readout.

**That is a data question, not a modelling one**, and it is the only acquisition on this page that
would move the ΔAge claim itself rather than the tool.
