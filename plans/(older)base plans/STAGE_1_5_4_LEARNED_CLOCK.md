# STAGE 1.5.4 — Can a model *learn* age from RNA? (the question M-2a never asked)

**Status:** ✅ **EXECUTED 2026-08-03. Verdict: NOT LEARNABLE (SPLIT on all three families).**
See §9 for results and §10 for three method corrections made *while blind to the real result*.

**Implements:** the last free route to repairing ΔAge — see §1.
**Depends on:** nothing unrun. Both datasets are already on disk.
**Blocking for:** nothing. It runs beside Stage 1.5.3 step 6 and does not touch it.

**Scope:** 1 new script, 1 test file, 0 lines changed in `src/`. **No label moves in this stage under
any outcome** — it measures whether better labels are *obtainable*, and producing them would be a
separate, separately-gated change.

---

## 1. Why this stage exists — the distinction M-2a did not draw

Stage 1.5.2's M-2a asked:

> *Does the **existing Fleischer clock's scalar output** track methylation age on reprogramming cells?*

Answer: **no.** ρ_partial **0.267** (skin & blood) / **0.516** (multi-tissue) → SPLIT ⇒ NOT CALIBRATABLE.

**That is a fact about one clock, not about the transcriptome.** The Fleischer clock was:

* fitted on **quiescent dermal fibroblasts**, not reprogramming cells;
* fitted against **chronological** age, never against **methylation** age;
* a dense RidgeCV over **33,155 genes from 133 samples**, with no feature selection.

So its failure says that *this* mapping does not transfer. It does **not** establish that the
information is absent from the transcriptome. **Nobody has asked the second question**, and it is
the last route to repairing ΔAge that costs nothing:

> **Can a model *trained* on the transcriptome predict methylation age on these cells?**

If yes, ΔAge labels become obtainable at scale for every cell we hold. If no, the RNA route is
closed on evidence rather than by inference from one clock, and Stage 6's spend is justified rather
than assumed.

**Honest prior: this probably fails.** M-2a's asymmetry finding (RNA↔multi-tissue 0.516 vs
RNA↔skin&blood 0.267, while the two methylation clocks agree with each other at **0.568**) is
evidence that RNA in these cells tracks something clock-specific — most likely identity, given
`corr(age_rna, pluripotency) = −0.62`. And ~45 training samples against ~20,000 genes is the same
p ≫ n regime that produced the Fleischer clock's problems. It is worth running because it is free
and decisive, not because it is likely.

---

## 2. Data — already on disk, nothing to acquire

| | |
|---|---|
| RNA | **GSE165177** — 95 samples, transient reprogramming |
| methylation | **GSE165179** — 96 samples, the same experiment |
| joined | **68 conditions**, donors **O1 (22), O2 (23), O3 (23)**, days 10/13/15/17, 6 arms |
| donor ages | O1 = 53, O2 = 53, O3 = 38 |

This is the same join M-2a used, reproduced on a second machine with every ρ bit-identical.

---

## 3. Design

**Target.** `age_meth` — the clock's **linear predictor**, not years. This is the natural regression
target and it is what M-2a's rows already carry. *(Note for readers: the field is named `age_meth`
but holds the lp; `anti_trafo` maps lp → years, and is linear above lp = 0.)* Both Horvath clocks
are targets, and **both are always reported**, whichever way they fall.

**Split — leave-one-donor-out across O1 / O2 / O3.** Everything is fitted inside the fold:
standardisation, any feature selection, and the regularisation strength. **The held-out donor
influences nothing.** With 3 donors this means ~45 training samples per fold.

**Estimand — pooled ρ_partial at n = 68.** Out-of-fold predictions are collected for all 68
conditions and scored *once*, pooled.

> **Why pooled and not per-fold:** the training split and the estimand are different choices.
> 3-donor LODO is what prevents leakage; it is **not** what the statistic is computed on. §5b's
> record already settles which geometry works here — pooled ρ_partial at n = 68 is **RESOLVABLE at
> 0.9940**, while every smaller geometry tried (n = 22, n = 11, n = 9) came back UNRESOLVABLE. Using
> the pooled estimand also makes the result **directly comparable to the Fleischer clock's 0.267 /
> 0.516**, because it is the identical metric on the identical conditions.

**Model families — all pre-registered, all reported, none selected on the answer.**

| | family | what a pass would mean |
|---|---|---|
| **F1** | ridge on the full transcriptome | the information is in the transcriptome somewhere |
| **F2** | ridge restricted to the Fleischer clock's gene set | the **genes** were right and only the **weights** were wrong |
| **F3** | PCA(k) → ridge | the signal is low-dimensional and dense |

Regularisation strength is chosen by **inner** cross-validation on the training donors only.
Reporting all three prevents the failure where one family is tried, then another, until something
passes.

---

## 4. The confound, and the two guards that make the result mean anything

**The confound is the same one M-2a faced, and it is fatal if unhandled.** Both RNA and methylation
move with reprogramming progress: methylation age falls 24–28 yr along that axis, and
`corr(age_rna, pluripotency) = −0.62`. **A model that learns only "how far along is this cell" would
score high on a raw correlation while carrying no age information at all.**

**Guard 1 — ρ_partial is the only pass criterion.** Pluripotency is partialled out using the
existing `OSKM_PLURIPOTENCY` signature, reused verbatim so it cannot be tuned for this stage.
ρ_all and ρ_within are computed and reported but are **never** a pass criterion — the same rule
M-2a operated under.

**Guard 2 — a label-shuffle null, which M-2a did not need and this stage does.** With ~20,000
features and ~45 training samples, a pipeline can manufacture correlation from nothing. So the whole
LODO procedure is re-run with the **training** targets shuffled inside each fold.

> **If the shuffled control does not collapse to ρ_partial ≈ 0, the pipeline is broken and no
> positive result from it may be believed.** This is a gate on the *method*, checked before the real
> result is read, not a sensitivity analysis reported beside it.

---

## 5. Bars — frozen before any model is fitted

| criterion | bar | why this value |
|---|---|---|
| **Primary** — pooled ρ_partial, out-of-fold | **≥ 0.50 on BOTH clocks** | identical to M-2a's registered bar, so the two stages are directly comparable. RESOLVABLE at n = 68 (**0.9940**) |
| **SPLIT** | one clock passes, one fails | counts as a **failure**, per M-2a §7 |
| **Guard 2** — shuffled-label null | pooled ρ_partial **≤ 0.20** on both clocks | a correct pipeline returns ≈ 0; 0.20 is slack for n = 68 sampling noise |
| **Secondary** *(descriptive, never a pass)* | beats Fleischer's 0.267 / 0.516 | reported for every family |

**No new resolvability simulation is required for the primary**, because the bar, the metric and the
geometry are *the same ones already registered and simulated* for M-2a — `stage_1_5_2_resolvability_results.json`,
row "M-2a rho_partial (ACTUAL)", n = 68, pass rate 0.9940. Re-using a registered bar is deliberate:
it removes any freedom to pick a friendlier one for a stage whose result I have already predicted
will be negative.

---

## 6. What each outcome licenses

| result | verdict | licenses |
|---|---|---|
| ρ_partial ≥ 0.50 on both clocks, shuffle collapses | **LEARNABLE** | a separate, separately-gated stage to *produce* labels. **Not** a label change here |
| SPLIT, or one clock only | **NOT LEARNABLE** | closing the RNA route on evidence. Stage 6's spend is then justified, not assumed |
| both < 0.50 | **NOT LEARNABLE** | as above, more strongly |
| shuffle does **not** collapse | **VOID** | the pipeline is broken. No result from this run may be reported at all |

**In three of four outcomes this stage produces no label change**, and that is the expected case. Its
value is that it converts *"the RNA route is probably closed"* into *"the RNA route is closed, and
here is the measurement"* — which is what Stage 6 needs in order to be a decision rather than a hope.

---

## 7. What this stage cannot do, stated up front

* **It cannot fix HFF.** Even a perfect RNA→methylation-age model gives HFF an accurate ΔAge of
  ≈ 0 — a neonatal line has almost no chronological age to remove. HFF's problem is biological, not
  instrumental.
* **It cannot escape 3 donors.** The training split has three, of which two share an age (53, 53, 38).
  A pass would need confirming on more donors before any label is produced from it.
* **It changes no label under any outcome.** Producing labels would be a separate stage with its own
  snapshot, guards and rollback.

---

## 8. Artefacts and verification

| file | role |
|---|---|
| `experiments/diag_learned_clock.py` | the measurement; read-only; pure logic separated from I/O |
| `tests/test_diag_learned_clock.py` | unit tests with no repo data present, per the pattern of the five existing `diag_*` scripts |
| `results/diag_learned_clock_results.json` | full output including per-family, per-clock, per-fold values |

**Verification of the stage itself:** `git diff --stat src/` **empty**; full suite green; the join
re-verified at 68 conditions before any model is fitted (§10 step 1 of 1.5.2's pattern — state the
shape before the statistic); and the shuffled-label gate passed before the real result is read.

```bash
python experiments/diag_learned_clock.py "<GSE165177 dir>" "<GSE165179 dir>"
```


---

## 9. ✅ RESULTS — executed 2026-08-03

**Verdict: NOT LEARNABLE.** All three model families return **SPLIT**, which §6 counts as a failure.

`n = 68` conditions × 35,720 genes, donors O1 (22) / O2 (23) / O3 (23), LODO, everything fitted
in-fold.

| ρ_partial (plu **and** donor removed) | skin & blood | multi-tissue | verdict |
|---|---:|---:|---|
| **Fleischer clock** *(baseline, same estimand)* | **0.309** | **0.517** | — |
| learned — full transcriptome | 0.277 | **0.627** | SPLIT |
| learned — clock genes only | 0.247 | **0.604** | SPLIT |
| learned — PCA(10) → ridge | 0.386 | **0.579** | SPLIT |

**Every one of the six beat its own permutation null** (null medians −0.06 to −0.11, q95 +0.09 to
+0.28). So the models *do* learn something real from the transcriptome — it simply is not age.

### 9.1 The finding that matters: the asymmetry survives retraining

M-2a's asymmetry was the strongest evidence that the RNA clock tracks something clock-specific
rather than age. **It reproduces here with models trained from scratch:**

| | sb | mt | ratio |
|---|---:|---:|---:|
| Fleischer | 0.309 | 0.517 | 1.67× |
| learned, full | 0.277 | 0.627 | **2.26×** |
| learned, PCA | 0.386 | 0.579 | 1.50× |

Training on the transcriptome **improves agreement with multi-tissue** (+0.06 to +0.11 over
Fleischer) and **does not improve agreement with skin & blood** (0.247–0.386 against 0.309, two of
three *worse*). A measurement of the shared age signal cannot behave that way: the two methylation
clocks agree with each other, so anything tracking what they share must track **both** at similar
strength.

> **The asymmetry is therefore a property of the DATA, not of the Fleischer clock.** Retraining
> does not remove it. That closes the last cheap route: it is not that we had the wrong clock — RNA
> in these cells does not carry the shared age signal.

### 9.2 What this licenses

**Does license:** closing the RNA route on measurement rather than inference. Stage 6's spend is now
justified rather than assumed — the cheap alternative was tried and it failed.

**Does not license:** any label change (none was possible under any outcome, §7); any claim about
HFF, which is neonatal and unfixable by better instrumentation; any claim beyond 3 donors.

---

## 10. ⚠️ Three method corrections, all made BEFORE any real ρ was read

G2 exists precisely so the method can be repaired without touching the result. It fired, and it was
right to. Recorded rather than quietly fixed.

**1. The LODO mean-reversion artefact — the first run was VOID.** With training labels shuffled,
ρ_partial came back at **−0.45 / −0.36**, not ≈ 0. Cause: a fold trained without donor *d* predicts
*d* with roughly the mean of the *other* donors, so the prediction is **anti-correlated with the
donor mean by construction**. A model that had learned nothing produced a large correlation.
**Fix:** partial out donor as well as pluripotency. That also sharpens the question to the one this
data can answer — within-donor tracking — since between-donor variance is 3 points at ages 53/53/38.
Pinned by `test_lodo_mean_reversion_is_removed_by_partialling_donor`, which reproduces it
synthetically.

**2. My own G2 bar was UNRESOLVABLE — the §5b failure, in this stage.** §5 registered
`|ρ| ≤ 0.20` on the worst of 6 comparisons, with "0.20 is slack for n = 68 sampling noise" asserted
and **never simulated**. At n = 68 a correct pipeline has SD(ρ) ≈ 0.122, so P(|ρ| > 0.20) ≈ 0.102
per comparison and **≈ 0.474 across six** — it would fire on a sound pipeline almost half the time.
**Fix:** replaced with a **permutation null** (20 draws per family × clock; the real value must beat
its own 95th percentile). This self-calibrates and needs no chosen constant.

**3. A PCA implementation bug.** Components were fitted on centred training data and then applied
to *uncentred* matrices, mixing the training mean into every score. Visible in G2 as a residual
artefact in the `pca` family only (−0.171 / −0.261 while the others sat near zero).

**All three were found and fixed with the real ρ values unread.** The bar that decides the verdict
— ρ_partial ≥ 0.50 on both clocks — was **not** touched.
