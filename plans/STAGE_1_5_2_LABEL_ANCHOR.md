# STAGE 1.5.2 — Is the transcriptomic clock CALIBRATABLE? (paired RNA ↔ methylation)

> ## ✅ **STATUS 2026-07-31: EXECUTED AND CLOSED. Answer — NO.**
>
> *§10 says "the status line at the top changes only when it has [run]". It has. **Everything below
> the original status line is the pre-registration exactly as written before the runs**; the results
> are §11–§16, appended, with nothing rewritten.*
>
> | | |
> |---|---|
> | **Verdict** | **NOT CALIBRATABLE** (§7 row 4). **Phase 2 does not run.** |
> | **M-2a** | SPLIT — ρ_partial **+0.267** / **+0.516** against a bar of 0.50 (§11) |
> | **§11's falsification check** | run, and **passed** — the anchor is sound, so the verdict is **accepted** (§12-R) |
> | **M-2b** | AGREE_FRAGILE at exactly 7/11 — but **0/3 at the one discriminating timepoint** (§14) |
> | **M-2c** | **not run**, gated on M-2a — and §12-R found a second, independent reason it would have been meaningless |
> | **G-a / G-b** | **closed** — the only `src/` changes in this stage, both record-only, ΔAge bit-identical (§13) |
> | **G-c step 1** | run — **RUN_STEP_2**, and it **refuted §0's own evidence** (§15) |
> | **What is still open** | exactly one thing: **G-c step 2**, which needs a retrain and belongs to whatever stage next rebuilds (§16) |
> | **Re-audit** | **§17** — every number re-verified against its JSON; **§11's per-arm *reading* is corrected** (the verdict is not) |
>
> **Read §16 first** — it is the one-page answer to every question this document asked.
> **Then §17**, which corrects how §11's per-arm table may be read and reports the strongest single
> result in the stage: in the arm where the two methylation clocks agree at **+0.936 — the sharpest
> reference of any arm — the RNA clock is negatively correlated.**

**Original status line, left as written:**
🔵 **PRE-REGISTERED — NOT YET EXECUTED, AND GATED.** Bars are fixed in §6 *before*
GSE165178 is opened. Nothing in this document reports a result.
**🔒 Does not start until the two prerequisites in §0 are met** — baseline-count visibility and
donor age wired into `src/`. Both are small; neither needs new data.

**Implements:** the step that no existing stage owns — see §1.
**Depends on:** **§0's G-a and G-b**. And `STAGE_1_5_1_REV_FINAL.md` (executed) for the clocks, the
transform, the pairing code and the negative-control result it rests on.
**Blocking for:** Stage 2's *premise* (§1), and any absolute ΔAge claim in Stage 5.
**Not blocking for:** Stage 3. The fate head does not consume ΔAge.

**Scope, Phase 1:** 1 new script, 1 new test file, 1 public download. **`src/` untouched.**
**Scope, Phase 2:** gated on Phase 1 (§8); touches labels; full snapshot and rollback required.

> **This document is additive.** No existing plan file is edited by it. Per the standing rule, only
> worked/done items are ever appended to a plan — and this file is marked NOT EXECUTED throughout
> until results exist.

---

## 0. 🔒 GATE — two small prerequisites, and one thing that is deliberately NOT a prerequisite

> ### ⚠️ CORRECTED 2026-07-30 (same day) — the first version of this gate was built on a false premise
>
> It read *"this stage does not start until D2 and D3 are closed"* and described both as unmeasured.
> **They were measured on 2026-07-24** by `experiments/diag_zero_point.py`; I asserted they were open
> without checking whether the diagnostic had already answered them. Worse, **D2's scientific half
> cannot be closed by analysis at all**, so the gate as written was *unsatisfiable* and would have
> blocked this stage indefinitely. Corrected below. The original wording is in git history.

**What `diag_zero_point.py` actually returned** (2026-07-24, `diag_zero_point_results.json`):

| | the question it asked | result |
|---|---|---|
| **M1** (D3's question) | does the clock read age on this data? | 🔴 **FAIL** — extreme contrast **11.8 yr** against a bar of **20.2 yr** on a true 53-yr gap, at power **0.996** |
| **M2** (D1's question) | is the zero-point's cross-batch structure driving the offset? | ✅ **NO_BATCH_EFFECT** — **−2.99 yr, 95% CI [−13.12, +7.14]**, n=12 |
| **M3** (D2's question) | is the per-donor offset real biology, or `n=1` baseline noise? | ⚠️ **INDETERMINATE** |

**M3 is the one that matters here, and it is decisive about what is *not* knowable.** Observed
offset SD **16.4 yr** against **12.3 yr** expected from a single unreplicated baseline ⇒ the baseline
explains **56% of the variance, 95% CI [9%, 100%]**, leaving 10.9 yr SD for biology + batch + model.

> **That CI spans almost the whole range. D2 is not unmeasured — it is measured and unresolvable at
> n = 6.** No further analysis closes it; it needs more donors. **So D2's scientific half is NOT a
> gate on this stage, and must never be written as one.**

Its recorded decision was **ESCALATE** — *"the clock does not separate the age extremes on this
data, so ΔAge's target is unvalidated… Stage 2's premise is void as stated."* **That is the strongest
existing statement of why this stage exists**, and M1's failure is precisely what the methylation
anchor is meant to resolve.

### The gate, restated to only what is actually closeable

| # | prerequisite | why it gates | cost |
|---|---|---|---|
| **G-a** | **`_control_baseline` records baseline count and composition.** `aging.py:81-90` averages whatever controls exist and reports neither; Stage 1.5 made `n=0` visible, **`n=1` is still silent** | M-2b's RNA-side contrast inherits that baseline. It must be *visible in the output* which donors rest on `n=1`, or the result cannot be interpreted | small code change + test |
| **G-b** | **Donor chronological age parsed in `src/`** (GEO declares it: N2/N3 = 0, Y1 = 29, Y2 = 35, O1/O2 = 53) | independent ground truth this stage would otherwise download methylation to approximate. *Unwired*, not unknown — `REV FINAL` §2 **used** these values as its guard and got **MAE 4.0 / 4.4 yr** | small, near-free |

**Explicitly NOT gates:**

* **D2's scientific question** — unresolvable at n=6 (above). Carried as a stated limitation, not a
  blocker.
* **D1** — answered, `NO_BATCH_EFFECT`. See §9-R3.
* **D3's scientific question** — answered: **M1 FAIL**. Only the *wiring* (G-b) remains.

**On G-b, do not oversell it:** donor age is a per-donor **constant**. It cannot measure
rejuvenation *within* a donor, so it does **not** replace methylation or make this stage
unnecessary — it anchors the **absolute** calibration question only.

**Gate condition:** G-a and G-b are implemented, tested, and recorded in `CHANGES.md` and the lab
notebook. Both are small and neither depends on new data. Until then this document stays
🔵 NOT EXECUTED.

> ## 🆕 G-c ADDED 2026-07-31 — decide HFF's `age_mask`. **Gates PHASE 2 ONLY, not Phase 1.**
>
> *Additive; no existing gate, table or bar is modified. Placement is deliberate: Phase 1 is
> measurement with `src/` untouched and must **not** be blocked by a decision that needs a retrain —
> that would repeat the over-broad gate this section already had to correct once.*
>
> ### The gap this closes
>
> **`age_mask` appears nowhere in this document or in `STAGE_1_5_1_REV_FINAL.md`** (verified by
> grep, 2026-07-31). So under **every** branch of §7 — including `CALIBRATABLE` — the model keeps
> training on **33,613 HFF ΔAge labels (99.78% of the training split)** that are produced by a clock
> now known to be **doubly out of scope for HFF**:
>
> | | |
> |---|---|
> | **out of age range** | HFF is a **neonatal** line; GSE113957 contains **0 samples below age 1**, and HFF's D0 baseline reads **84.5 yr** |
> | **out of domain** | reprogramming cells — established by `REV FINAL` §1, the finding this whole stage rests on |
>
> §2's ledger states the volume honestly ("≈99.8% stays unanchored either way") but **draws no
> consequence from it.** This gate draws it.
>
> ### The evidence already in hand
>
> `REV FINAL` established, on methylation, what real rejuvenation looks like in this system. The HFF
> RNA labels do not look like it:
>
> | | dose-response | monotone? | source |
> |---|---|---|---|
> | **methylation ground truth** (Gill intermediates) | **−3.30 / −3.15 yr/day**, ρ −0.885 / −0.842, p ≤ 0.0006 | **yes** | `REV FINAL` §4.4(b) |
> | **HFF RNA labels** (D0–D14 pseudobulk) | **−0.36 yr/day**, ρ **−0.214** | **no** | `diag_d2_replication_results.json` |
>
> ≈**9× weaker and non-monotone.** ⚠️ **Not decisive on its own** — different protocol (standard vs
> MPTR), different line, different modality — which is exactly why this is a gate with a test, not a
> conclusion.
>
> ### G-c — the decision, pre-registered
>
> **Question:** are HFF's ΔAge labels informative, or is the ΔAge head learning artefact from
> 33,613 contaminated labels while 75 usable ones are drowned out?
>
> **Step 1 (cheap, no retrain).** Test whether HFF's ΔAge labels carry the rejuvenation signature
> methylation established: Spearman(ΔAge, reprogramming day) over the HFF trajectory, with the
> iPSC endpoint excluded (a cell-type change, per the standing rule), and the dose-response slope
> reported beside the methylation figures above.
>
> **Step 2 (decisive, and it must ride in Phase 2's rebuild).** If Step 1 is not clearly positive,
> compare `age_mask=True` vs `age_mask=False` for HFF in **one** retrain, on the existing scorecard.
>
> | Step 1 result | action |
> |---|---|
> | HFF shows the signature (monotone, slope within ~2× of methylation's) | **keep** HFF labels; record the check |
> | ambiguous | run Step 2; decide on the scorecard, pre-registering the metric before the run |
> | no signature (current evidence points here) | **mask HFF's ΔAge in Phase 2's rebuild** and state the consequence plainly: the age head trains on ~75 labels, which may be too few — that is a finding, not a failure |
>
> ### Why this gates Phase 2 specifically
>
> Phase 2 already changes labels and pays for a rebuild + retrain. **Deciding HFF separately would
> cost a second rebuild** and would mean shipping a run that repairs 0.2% of labels while knowingly
> leaving 99.8% contaminated. **Both label decisions must ride in the same rebuild.**
>
> **What G-c does NOT do:** it does not block Phase 1, does not touch `src/` by itself, does not
> presuppose the answer, and does not claim the HFF labels are worthless — only that the question is
> unowned, cheap to ask, and larger in scope than the one Phase 2 answers.

---

## 1. Why this stage exists

Stage 1.5.1 proved the transcriptomic clock is **out of domain** on reprogramming cells, and proved
it with an identity-matched methylation design. It changed **nothing**: `src/` untouched, zero
labels moved. That was correct for a diagnosis, but it leaves a hole in the stage graph:

```
1.5.1   "the ΔAge labels are produced by an instrument that fails here"   ✅ done
  ???   "here is whether that instrument can be repaired"                 ← NO STAGE OWNS THIS
Stage 2 "correct per-donor offsets ON those labels"                       ⏸ blocked on exactly that
```

Stage 2's own annotation states the blocker verbatim: *"re-measure the per-donor level shift on
corrected labels before spending — **that was never done, because the labels were never
corrected**."* Nothing in Stages 2, 4 or 6 corrects them:

| stage | why it is not the owner |
|---|---|
| **1.5.1** | Closed and validated; scope is measurement-only; its §6.2 explicitly **withdrew** this test as ill-defined *on GSE165179*. |
| **Stage 2** | **Consumes** ΔAge labels. Downstream of the repair, not the repair. |
| **Stage 4** | Depends on Stage 3. Validates tool recommendations. |
| **Stage 6** | Owns *acquisition*. GSE165178 passes its hard gate cleanly (sorted arms with a per-sample outcome label = its ✅ ideal row), but Stage 6 is about **adding training data**, not anchoring existing labels. |

**This stage is that missing step, and only that step.**

---

## 2. The honest ledger — read this before §3

Stated up front so the stage cannot be oversold later, in the form 1.5.1 §8.2 should have used from
the start:

| | |
|---|---|
| does it settle whether the RNA clock tracks methylation age? | **yes — that is the whole point** |
| does it repair Gill's ΔAge labels? | **only if Phase 1 passes**, and only Gill's |
| how many training labels does that touch? | **≈75 of 33,688 — about 0.2%** |
| does it repair HFF? | **no.** No methylation exists for HFF in any public series |
| what fraction stays unanchored either way? | **≈99.8%** |

**So the value here is not label volume — it is the instrument verdict.** Whether the RNA clock is
calibratable at all is what Stage 2's premise and Stage 5's claims actually rest on, and it is
currently unmeasured. A negative result is as useful as a positive one and costs the same.

**What must NOT be claimed on the back of this stage:** that ΔAge is fixed; that HFF is anchored;
that rejuvenation is demonstrated in our labels. §7 pins each outcome to exactly what it licenses.

---

## 3. The data and the join

**GSE165178** `[Sendai array]` — Illumina methylation on the **same samples** as GSE165176, the RNA
series the model trains on. Both are SubSeries of SuperSeries **GSE165180**.

Verified in `STAGE_1_5_1_REV_FINAL.md` §8.3 against real sample titles:

* 22 methylation samples, titled `{donor}_{day}_{marker}` (e.g. `Y2_d11_SSEA4`);
* our RNA titles are that same key plus a batch suffix (`Y2_d11_SSEA4_Sendai_Exp1`);
* **22/22 join on `donor_day_marker`, zero unmatched**;
* donors **O1, O2, Y1, Y2** (4 of 6 — the missing N2/N3 are neonatal, outside the clock's fitted
  range anyway); days **9, 11, 15**;
* the sort marker **is** the arm label: `CD13` → *Failing to reprogram fibroblast*,
  `SSEA4` → *Reprogramming fibroblast*.

**Known geometry constraint, stated before the run.** 4 donors × 3 days × 2 markers = 24 cells in
the grid; 22 exist, so the grid is nearly full but **there is no day-0 arm and no untreated arm in
GSE165178.** Consequences, pre-committed:

1. The only internal contrast available is **SSEA4 vs CD13**.
2. CD13 is a *treated non-responder*, not an untreated control. 1.5.1 showed non-responders are
   **inert on methylation** (+0.5 / −2.4 yr) — but that was measured in the **transient** arm
   (GSE165179), and GSE165178 is the **Sendai** (continuous OSKM) arm. **The inertness result does
   not automatically transfer.** It is an assumption here, and §9-R1 says what happens if it fails.
3. **The exact day grid and arm counts are re-verified on download** (§10 step 1) before any
   statistic is computed — ground rule §11, state the shape before the statistic.

> ## 🆕 ADDED 2026-07-31 — **GSE165177 is a second, larger paired set, and it fixes §3's geometry constraint**
>
> *Additive; §3 above is unmodified. Found by checking the SuperSeries rather than assuming
> GSE165178 was the only pairing.*
>
> **GSE165180 has four SubSeries, and we hold two.** The plan names GSE165178. But **GSE165177
> `[Transient_RNAseq]` is the RNA half of the experiment whose methylation we already hold and have
> already validated** (GSE165179 — `REV FINAL` §3's contrasts A/B/C, which reproduce exactly on
> re-run).
>
> | | samples | pairs with | untreated arm? | status |
> |---|---|---|---|---|
> | **GSE165177** `[Transient_RNAseq]` | **95** | **GSE165179 — already held & validated** | ✅ **yes** (`negative_control`) | ⬇ **download** |
> | GSE165178 `[Sendai array]` | 22 | GSE165176 (our training RNA) | ❌ **none** (§3 above) | ⬇ download |
>
> **Title-format match, verified:** GEO lists GSE165177 titles as
> `O1_transiently_reprogrammed_17days_exp1`, `O3_failed_to_transiently_reprogram_15days_exp2`. Our
> local GSE165179 titles read `O1_failed_to_transiently_reprogram_15days_exp1` — **same
> `{donor}_{arm}_{N}days_{exp}` key, same arm vocabulary, same donors (O1/O2/O3), same day grid
> (10/13/15/17).** *(From GEO metadata; the exact join is re-verified on download per §10 step 1.)*
>
> ### Why this materially improves the stage
>
> 1. **It removes §3's stated geometry constraint.** GSE165178 has *"no day-0 arm and no untreated
>    arm"*, forcing the only contrast to be SSEA4-vs-CD13 and requiring the **assumption** that
>    non-responder inertness transfers from the transient to the Sendai arm (§3 point 2, §9-R1).
>    GSE165177/165179 carries a **real `negative_control` arm**, so that assumption is not needed —
>    the very contrast REV FINAL validated.
> 2. **≈95 pairs instead of 22.** §6 anticipates that **ρ_within may return UNRESOLVABLE at n=11**
>    and pre-commits a fallback. At this sample size that risk largely disappears, and the fallback
>    may not be needed.
> 3. **The methylation half is already measured, validated, and re-run.** No new methylation
>    analysis is required — only the RNA side is missing.
>
> ### The two are complementary, not alternatives
>
> | series | answers |
> |---|---|
> | **GSE165177** | **the verdict** — *is the RNA clock calibratable against methylation?* (§2 calls this "the value here") |
> | **GSE165178** | **the label attachment** — the only series that joins the samples the model actually trains on |
>
> If only one is downloaded, **take GSE165177**: the stage's primary product is the verdict, and this
> is the better instrument for it.
>
> ### Files to fetch — minimum set only
>
> The measurement script reads **exactly two files** per series
> (`diag_methylation_anchor.py:383`); the `signal_intensities` matrix is **never opened** and the
> IDAT `RAW.tar` is not needed.
>
> | series | file | size |
> |---|---|---|
> | GSE165177 | `GSE165177_series_matrix.txt.gz` | small |
> | GSE165177 | `GSE165177_Log2_RPM_Transient_reprogramming.txt.gz` | 3.2 MB |
> | GSE165177 | `GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz` | 12.1 MB |
> | GSE165178 | `GSE165178_series_matrix.txt.gz` | small |
> | GSE165178 | `GSE165178_Matrix_processed_sendai.txt.gz` | 142.5 MB |
>
> **Do NOT fetch** `*_signal_intensities_*.txt.gz` (188.5 MB) or `GSE165178_RAW.tar` (466.7 MB) —
> unused. *(Note: the held `GSE165179_Matrix_signal_intensities_transient.txt.gz`, 807 MB, is also
> unused and can be deleted.)*
>
> ### One implementation caveat, flagged now
>
> GSE165177's titles use the **transient** vocabulary (`O1_transiently_reprogrammed_17days_exp1`),
> **not** the Sendai format `GillReprogrammingSource` parses (`N2_d11_CD13_Sendai_Exp1`), and its
> `Log2 RPM` matrix is **split across two files**. Expect a small dedicated loader rather than reuse
> of the production source. The clock application itself is unchanged.

---

## 4. The confound that decides the whole design

**This is the hole that would sink a naive version of this stage, so it is closed here rather than
flagged.**

The obvious test — "does `age_rna` correlate with `age_meth` across the 22 samples?" — is
**not sufficient**, and a high correlation would prove nothing on its own.

Reason: 1.5.1 measured `corr(age_rna, pluripotency) = −0.62`. Methylation age **also** falls sharply
during reprogramming (−24 to −27 yr, 1.5.1 §3). So **both modalities move with reprogramming
progress.** A clock that detects only "is this cell reprogramming?" — carrying no age information
whatever — would still produce a strong `age_rna` ↔ `age_meth` correlation across a sample set whose
dominant axis is exactly that. This is the same identity artefact that produced the +36.5 yr
reading, re-entering through the back door.

**Therefore the headline correlation is reported as descriptive only and is NOT a bar.** Three
readings are computed, and the decisive ones are the confound-free ones:

| reading | n | role |
|---|---|---|
| **ρ_all** — across all samples | 22 | descriptive only. **Never a pass criterion.** |
| **ρ_within** — within each arm separately | 11 + 11 | **decisive.** Within CD13 the cells are not reprogramming, so a surviving correlation cannot be identity |
| **ρ_partial** — partialling out the pluripotency score | 22 | **decisive.** Removes the shared axis directly |

The pluripotency score is the existing `OSKM_PLURIPOTENCY` signature share
([diag_clock_validity.py:69](experiments/diag_clock_validity.py:69)) — reused, not reinvented, so
it cannot be tuned for this stage.

**Stage 4 G1 compliance (negative control).** CD13 *is* the negative-control arm: cells that got
OSKM and did not reprogram. Every statistic in §5 is reported on it. Per G1, **a claim whose
negative control moves with it is not validated** — and here that has teeth: if ρ_within(CD13) is
the *only* thing that holds, the agreement is real age signal; if ρ_all is high while both ρ_within
and ρ_partial collapse, the agreement is identity and the stage returns a **negative** verdict.

---

## 5. The three measurements

Run in order. **Each gates the next** — M-2c is meaningless if M-2a fails, because no monotone
calibration can repair a clock that is not tracking the target at all.

### M-2a — Does the RNA clock track methylation age? *(the instrument question)*

Spearman ρ between `age_rna` and `age_meth` per sample, in the three readings of §4, on **both**
Horvath clocks (skin & blood 2018, multi-tissue 2013), both always reported whichever way they fall.

Absolute ages, no control needed. This is the control-free core of the stage.

### M-2b — Do the two modalities agree on the CONTRAST? *(the ΔAge-shaped question)*

Per (donor, day): `Δ = mean(SSEA4) − mean(CD13)` computed independently on each modality; replicates
**averaged, not treated as independent** (1.5.1's unit-of-analysis rule, and the fix for the pairing
bug that once silently dropped 6 of 9 pairs — `pair_by_donor_day` is reused verbatim, with its
regression test).

Reported: sign agreement across pairs, and ρ(Δ_rna, Δ_meth).

**Expected direction, pre-committed:** methylation says reprogramming cells are **younger**
(1.5.1: −24 to −27 yr). RNA has been reading treated cells **older** (+36.5 yr on non-responders).
**Disagreement is the live hypothesis, and it is a decisive result, not a failure of the stage.**

### M-2c — Can a correction be learned? *(gated on M-2a)*

Fit a **monotone** map `age_meth ≈ f(age_rna)` — affine first (2 parameters, honest at n=22),
isotonic reported as sensitivity only. Validate **leave-one-donor-out across the 4 donors** —
never leave-one-sample-out, which would leak donor identity.

Reported: LODO MAE, against the RNA clock's uncorrected `cv_mae = 12.27 yr` as the reference point.

---

## 6. Pre-registered bars, and their resolvability

Ground rule §5b: **a bar that a correct system fails is not a bar.** It is checked forward, before
the run, via `audit_metrics.bar_verdict(sim_null, bar, lower_is_better)` against
`MIN_PASS_RATE = 0.95`.

**Procedure, fixed now and executed before GSE165178 is opened:**

1. Simulate each metric at the geometry it will be graded on, under a system that **meets the intent
   exactly** (recipes below).
2. Read `verdict`. If **RESOLVABLE** at the proposed bar, that bar is frozen. If **UNRESOLVABLE**,
   take `usable_bar`, or change the geometry, or drop the criterion — per §5b, **now, not after**.
3. Write the resulting numbers into this file, add each to `tests/test_bars_resolvable.py`
   (a bar with no resolvability test is **not considered pre-registered**), and commit **before**
   the download is touched.

| criterion | geometry | intent simulated | proposed bar |
|---|---|---|---|
| **M-2a** ρ_within | n=11 per arm | bivariate normal, ρ_true = **0.70** — the level at which monotone calibration is worth attempting | ρ ≥ 0.50 |
| **M-2a** ρ_partial | n=22, 1 covariate | same, ρ_true = 0.70 | ρ ≥ 0.50 |
| **M-2b** sign agreement | up to 11 pairs | Bernoulli(0.85) per pair | ≥ 8/11 |
| **M-2c** LODO MAE | 4 donor folds | residual SD = 5 yr (≈ the 4.0/4.4 yr methylation guard spread) | ≤ 8.0 yr |

**Two things are anticipated now rather than discovered later:**

* **n=11 is small, and ρ_within may well come back UNRESOLVABLE.** Pre-committed fallback, so it
  cannot be chosen after seeing data: **if ρ_within is UNRESOLVABLE at n=11, the decisive criterion
  becomes ρ_partial at n=22**, and ρ_within is demoted to descriptive. The fallback is fixed here,
  in advance, and does not change what a pass means (§7).
* **M-2c is gated on M-2a.** If M-2a fails, M-2c is **not run and not reported** — fitting a
  calibration to a clock that is not tracking the target would manufacture a number with no meaning.

**Both clocks must agree for a pass.** A criterion met on one clock and not the other is recorded as
**SPLIT**, which is a failure for the purpose of §7, not a pass.

> ## 🔒 EXECUTED 2026-07-31 — bars frozen at the ACTUAL geometry, before any statistic
>
> *§6's procedure, run by `experiments/stage_1_5_2_resolvability.py` →
> `stage_1_5_2_resolvability_results.json`. **Titles only** — no expression or beta value was read,
> so no bar could be tuned to the answer. `src/` untouched.*
>
> **Actual geometry** (GSE165177 × GSE165179, verified on download): 90 joined pairs → **68 unique
> (donor, arm, day) conditions**; donors **O1, O2, O3** (**3**, not the 4 the LODO bar assumed);
> days 10/13/15/17; smallest arm **9** conditions (`transiently_reprogrammed`).
>
> | bar | geometry | pass rate | verdict |
> |---|---|---|---|
> | M-2a ρ_within | registered n=11 | 78.1% | ❌ UNRESOLVABLE (usable 0.236) |
> | M-2a ρ_within | actual n=9 | 75.8% | ❌ UNRESOLVABLE (usable 0.167) |
> | **M-2a ρ_partial** | **registered n=22** | **92.3%** | ❌ **UNRESOLVABLE (usable 0.457)** |
> | **M-2a ρ_partial** | **actual n=68** | **99.4%** | ✅ **RESOLVABLE** |
> | M-2b sign agreement ≥8/11 | 11 pairs | 93.0% | ❌ UNRESOLVABLE (usable **7**/11) |
> | M-2c LODO MAE ≤8.0 | 4 folds (registered) | 100% | ✅ RESOLVABLE |
> | M-2c LODO MAE ≤8.0 | **3 folds (actual)** | 100% | ✅ RESOLVABLE |
>
> ### 🔴 The finding that matters: the pre-committed fallback would also have failed
>
> §6 anticipated ρ_within might be UNRESOLVABLE and fixed a fallback — *"the decisive criterion
> becomes ρ_partial at n=22."* **Measured, that fallback is itself UNRESOLVABLE at 92.3%.** On the
> GSE165178-only geometry the stage would have had **no valid decisive criterion at all**, and would
> have produced a verdict that ground rule §5b forbids relying on.
>
> **GSE165177 is what rescues it:** at the actual n=68, ρ_partial passes at **99.4%**. This is a
> measured vindication of adding that series, not an argument for it.
>
> ### Bars as frozen (these govern §7 from here)
>
> | criterion | status |
> |---|---|
> | **M-2a ρ_partial ≥ 0.50 at n=68** | ✅ **DECISIVE** — the criterion §7 is graded on |
> | M-2a ρ_within | **demoted to descriptive** — the §6 fallback fires as written. Reported, never gated on |
> | **M-2b sign agreement** | bar **moves ≥8/11 → ≥7/11** (its `usable_bar`), per §5b: move it now, not after |
> | **M-2c LODO MAE ≤ 8.0 yr** | ✅ frozen, **3 folds** — the reduced donor count does not weaken it |
>
> **Unchanged:** both clocks must agree or it is SPLIT; M-2c stays gated on M-2a.
>
> ⚠️ **Honest note on M-2b:** moving a bar to its `usable_bar` makes it easier to pass. It is done
> here because §5b requires an unresolvable bar to be moved *before* the run — but it is a **weaker**
> test than registered, and any M-2b pass must be reported with that caveat attached.

---

## 7. What each outcome licenses — decided before the run

| M-2a | M-2b | verdict | what it licenses | what it does NOT license |
|---|---|---|---|---|
| pass | agree | **CALIBRATABLE** | Phase 2 (§8) for Gill's ≈75 cells | anything about HFF |
| pass | disagree | **TRACKS-BUT-BIASED** | Phase 2 with the sign/scale correction from M-2c, *and* a mandatory re-derivation of every 1.5.1 contrast under it | treating the old labels as merely noisy |
| ρ_all high, ρ_within **and** ρ_partial collapse | either | **IDENTITY, NOT AGE** | closing the RNA-clock route permanently; ΔAge comes from methylation where it exists, Gill-only | any RNA-derived ΔAge claim, including existing ones |
| fail | either | **NOT CALIBRATABLE** | same as above | Phase 2. It does not run |

**In three of four outcomes Phase 2 does not happen, and that is a real result.** The two negative
verdicts are the strongest thing this stage can produce: they retire a route the project has already
spent four failed attempts on, and they do it on paired ground truth rather than argument.

---

## 8. Phase 2 — the label change, gated

**Runs only on CALIBRATABLE or TRACKS-BUT-BIASED.** This is deliberately a separate phase with its
own risk class, because it is the one place in this project where the **training target itself
moves**: every guard, every scorecard number and every downstream stage moves with it.

Requirements, all mandatory:

1. **One change only** (ground rule §2). Applying the calibration is the whole change; nothing else
   in the same commit.
2. **Full snapshot before, rollback path verified** — not assumed, exercised.
3. **Every Stage 1 guard re-run** and reported before/after. A guard that moves is a finding to
   pre-register, not something to absorb.
4. **The HFF question is answered explicitly, not by silence.** Applying a Gill-learned correction
   to HFF is a *separate, additionally justified* decision. Default is **do not apply it**: HFF is a
   different cell system and there is no ground truth to validate transfer against. If it is applied,
   the justification and the plausibility evidence go in this file first.
5. **1.5.1's contrasts are re-derived** under the corrected labels and any change reported. 1.5.1's
   methylation results are independent of this and cannot move; its *RNA-side* statements can.

---

## 9. Risks and failure modes — with the response fixed in advance

**R1 — CD13 is not inert in the Sendai arm.** The inertness result (+0.5/−2.4) is from the transient
arm. If Sendai non-responders have genuinely drifted, `Δ` in M-2b is measured against a moving
baseline. **Response:** M-2a is *unaffected* (it uses absolute ages, no control), so the instrument
verdict survives regardless. M-2b is reported with the caveat and is **not** used alone for any
verdict. Additionally, CD13's absolute methylation ages are compared to donor chronological age —
a direct check of drift that costs nothing.

**R2 — 4 donors is few.** LODO across 4 folds is weak. **Response:** the §6 resolvability check is
what decides whether M-2c's bar is even meaningful at that geometry; if UNRESOLVABLE, M-2c is
reported descriptively and cannot trigger Phase 2 on its own.

**R3 — the batch suffix hides a real batch effect.** RNA titles carry `Exp1`/`Exp2`; methylation may
not.

> **⚠️ Corrected 2026-07-30 — this was already measured, and an earlier draft of this file re-raised
> it as if it were novel.** It is **finding D1**, and Stage 1.5's follow-up settled it: the paired
> Exp1−Exp2 offset is **−2.99 yr, 95% CI [−13.12, +7.14], n = 12 → `NO_BATCH_EFFECT`**. Re-deriving
> a closed finding as a fresh "risk" is exactly the drift §0 exists to stop, so the measured number
> is cited here instead.

**Response:** RNA replicates are averaged per `donor_day_marker` before joining — the same rule
1.5.1 uses — and the exp1/exp2 spread is still reported, because **the D1 null must not be
over-read**: its CI half-width (~10 yr) excludes a *large* batch effect, not a meaningful one.
D1 also remains **structurally** true (all baselines `Exp2`, ~50% of samples `Exp1`), so the spread
is reported as a visible quantity rather than assumed negligible.

**R4 — probe coverage differs from GSE165179.** 1.5.1 got 100% / 94.6% CpG coverage there.
**Response:** coverage is recomputed and reported per clock; below 90% the clock is reported as
degraded, matching 1.5.1's handling.

**R5 — I am wrong about the join.** It was verified on titles, and titles can mislead. **Response:**
§10 step 1 re-verifies 22/22 on the downloaded matrices before anything else runs, and the script
**aborts** rather than proceeding on a partial join.

---

## 10. Execution, verification, artefacts

**Step 1 — shape before statistic.** Load GSE165178. Print and record: sample count, donor roster,
day grid, arm counts, probe count, CpG coverage per clock, and the join result against GSE165176.
**Abort if the join is not 22/22.** No statistic is computed in this step.

**Step 2 — freeze the bars.** Run the §6 resolvability simulations, write the frozen numbers into
§6, add tests to `tests/test_bars_resolvable.py`, commit.

**Step 3 — run M-2a → M-2b → M-2c**, respecting the gates.

**Step 4 — record.** `CHANGES.md`, and pre-registration + result in
`experiments/DELTAAGE_LAB_NOTEBOOK.md` with predictions written **before** step 3.

```bash
python experiments/diag_label_anchor.py "<dir containing GSE165178>"
```

**Artefacts:** `experiments/diag_label_anchor.py`, `tests/test_diag_label_anchor.py`,
`diag_label_anchor_results.json`, new rows in `tests/test_bars_resolvable.py`.

**Verification of the stage itself:**

* `git diff --stat src/` is **empty** for Phase 1 — the guarantee that no guard can move
* full suite green (currently 455 passing)
* every pure function unit-tested with no repo data present, per the pattern of the four existing
  `diag_*` scripts
* `pair_by_donor_day` reused **verbatim** from `diag_methylation_anchor.py`, keeping its regression
  test — the replicate-dropping bug is not to be reintroduced by a copy

---

## 11. What would falsify this stage's own conclusion

Stated so a reviewer does not have to construct it:

* **A positive verdict is falsified** if the correlation is shown to be carried by the arm axis after
  all — so ρ_within and ρ_partial are the criteria, and ρ_all is explicitly barred from being one.
* **A negative verdict is falsified** if the methylation ages themselves are unreliable here — so the
  clocks are checked against donor chronological age on the CD13 arm (R1) before any negative
  verdict is accepted.
* **Both are falsified** if the join is wrong — so it is re-verified on the actual matrices, with an
  abort, before anything else runs (§10 step 1).

**Nothing in this file may be reported as a result until it has run.** The status line at the top
changes only when it has.

---

# 11. ✅ M-2a EXECUTED 2026-07-31 — **SPLIT ⇒ NOT CALIBRATABLE. Phase 2 does not run.**

`python experiments/diag_m2a_calibratability.py "D:\GSE165177" "D:\GSE165179"` →
`diag_m2a_calibratability_results.json`. Bars were frozen **before** this ran (§6 block above),
from titles only. **`src/` untouched.**

| | ρ_all *(descriptive)* | **ρ_partial** *(DECISIVE, bar ≥ 0.50)* | verdict |
|---|---|---|---|
| Horvath **skin & blood** 2018 | +0.444 | **+0.267** | ❌ **FAIL** |
| Horvath **multi-tissue** 2013 | +0.690 | **+0.516** | ⚠️ pass by **0.016** |

**§6: "A criterion met on one clock and not the other is recorded as SPLIT, which is a failure for
the purpose of §7."** → **M-2a FAILS** → §7 row 4: **NOT CALIBRATABLE. Phase 2 does not run.**

⚠️ **The one "pass" is FRAGILE and must not be quoted as a pass.** +0.516 clears 0.50 by **0.016**.
This project has been burned twice by exactly this (E1b cleared its bound by 0.009, D2 by 0.014, in
opposite directions — `REV FINAL` records the lesson). At n=68 a margin of 0.016 is noise.

### The confound control did real work — §4 was right to insist on it

| clock | ρ_all → ρ_partial | drop |
|---|---|---|
| skin & blood | 0.444 → 0.267 | **−40%** |
| multi-tissue | 0.690 → 0.516 | **−25%** |

A quarter to a *half* of the raw agreement is the shared reprogramming axis, not age. Had §4 not
pre-committed to partialling out pluripotency, ρ_all = +0.69 would have been reported as a
comfortable pass — and it would have been the identity artefact re-entering "through the back door",
exactly as §4 predicted.

### The within-arm pattern — descriptive, but consistent across 6 arms × 2 clocks

*(ρ_within was demoted to descriptive by the §6 freeze — UNRESOLVABLE at n=9–12 — so no single
number here is reliable. The pattern across all twelve is worth recording anyway.)*

| arm | skin & blood | multi-tissue | reprogramming? |
|---|---|---|---|
| `failed_to_transiently_reprogram` | +0.322 | **+0.538** | no |
| `negative_control` | +0.413 | +0.385 | no |
| `negative_control_intermediate` | +0.091 | +0.371 | no |
| `failing_..._intermediate` | +0.140 | +0.084 | no (failed) |
| **`transiently_reprogrammed`** | +0.267 | **+0.033** | **yes** |
| **`transient_reprogramming_intermediate`** | **−0.136** | **−0.191** | **yes** |

**The RNA clock tracks methylation age in cells that are NOT reprogramming, and stops tracking — or
inverts — in exactly the cells that are.** That is the same in-domain/out-of-domain boundary
`REV FINAL` §1 established, now visible in the calibration itself.

### Pipeline validated before the negative result was accepted

The methylation side reproduces `REV FINAL` contrast A independently:
`negative_control_intermediate` +1.487 vs `transient_reprogramming_intermediate` +0.348 →
**1.139 × 21 = 23.9 yr**, against REV FINAL's **−24.1 yr**. The failure is real, not a wiring bug.

### What this licenses (§7, decided before the run)

| ✅ licensed | ❌ not licensed |
|---|---|
| **Closing the RNA-clock route.** Five attempts have now failed: refit, precision, control-swap, statistical fixes, and calibration | Phase 2. **It does not run** |
| ΔAge from **methylation where it exists** (Gill's 3 donors, `REV FINAL` §3) | Any RNA-derived ΔAge claim — **including the existing labels** |
| **G-c becomes the live question** (§0): if the labels cannot be repaired, the decision is whether to keep training on them | Treating the current labels as "merely noisy" |

**Not run: M-2b and M-2c.** M-2b is the ΔAge-shaped question and is gated on **G-a**; M-2c is gated
on M-2a, which failed — §6: *"fitting a calibration to a clock that is not tracking the target would
manufacture a number with no meaning."*

---

# 12. 🔵 PRE-REGISTERED 2026-07-31 — §11's falsification check on M-2a. **NOT YET RUN.**

> **§11 sets a precondition that §11's own result did not meet.** It says: *"A negative verdict is
> falsified if the methylation ages themselves are unreliable here — so the clocks are checked
> against donor chronological age on the CD13 arm (R1) before any negative verdict is accepted."*
>
> **That check had not been run when M-2a was recorded.** So the SPLIT ⇒ NOT CALIBRATABLE verdict
> above is *recorded but not yet accepted*, and this section is the missing precondition. Written
> and committed **before** the measurement, so no bar and no decision rule can be chosen after
> seeing the answer.

Script: `experiments/diag_r1_anchor_reliability.py`. Read-only, `src/` untouched. Phase 1
(pre-registration) has run and reads no beta value; phase 2 is the measurement.

### The four checks

| | question | geometry | intercept-dependent? |
|---|---|---|---|
| **R1a** | LODO chronological-age recovery — intercept from the **other** donors, predict the held-out one | 3 donors, 1 day-0 sample each | no (nothing grades its own intercept) |
| **R1b** | the real 15-yr gap (O1/O2 = 53, O3 = 38) on the **untreated** `Negative control` arm | 21 samples, 3 donors | **no** — `age_i − age_j = 21·(lp_i − lp_j)`, REV FINAL §4.3 |
| **R1c** | drift of the treated **non-responder** arms vs their own day-0 (§9-R1's own words) | 33 samples | no (within-donor difference) |
| **R1d** | do the **two Horvath clocks agree with each other**, under M-2a's own criterion? | **68 conditions** | no |

### Bars, frozen before the data was opened

| bar | pass rate for a correct clock | verdict | action taken |
|---|---|---|---|
| R1a LODO MAE ≤ **5.0** yr | 76.7% | ❌ UNRESOLVABLE | bar **moved to 7.17 yr** (§5b: now, not after) |
| R1b \|gap − 15\| ≤ **7.07** yr | 87.5% | ❌ UNRESOLVABLE | bar **moved to 9.08 yr** |
| **R1d ρ_partial ≥ 0.50 at n=68** | **99.4%** | ✅ **RESOLVABLE** | **reused verbatim from §6 — not re-derived** |

Null: a clock with published Horvath accuracy (MAE 3.0 yr), error modelled as a **donor-level
offset with no averaging benefit** — deliberately pessimistic, since an independent-error model
would shrink the SE by √7 and make every bar look easily resolvable.

### 🔴 The problem with R1a and R1b, stated before they are read

This series has **3 donors carrying only 2 distinct ages.** Both chronological-age checks came back
UNRESOLVABLE and had to be loosened. **Loosening a bar that gates a negative verdict makes that
verdict harder to falsify** — a bias in favour of the result already recorded. So neither R1a nor
R1b may carry this decision alone, and that is fixed here rather than discovered later.

### Why R1d is the decisive one

R1d asks the same question — *are the methylation readings reliable on these samples?* — as: **do
the two independent Horvath clocks agree with each other**, over the same 68 conditions, under the
same pluripotency partialling, against **the identical ρ_partial ≥ 0.50 bar that M-2a was graded
on**. No new bar is introduced.

That symmetry is the point. RNA↔methylation and methylation↔methylation are scored by the same
criterion on the same samples, so **if one passes and the other fails, the failure is localised to
an instrument rather than to the data.** The CpG overlap between the two clocks is reported
alongside, because agreement driven by shared probes would be partly trivial.

### Decision rule, pre-committed

| R1d | R1a / R1b | outcome |
|---|---|---|
| **FAIL** | either | **M-2a's negative verdict is WITHDRAWN.** The anchor cannot arbitrate the RNA clock if it cannot arbitrate itself |
| PASS | both pass | **M-2a's negative verdict is ACCEPTED.** §11's falsification condition is not met |
| PASS | either fails | **ACCEPTED WITH CAVEAT** — the well-powered check supports the anchor, the under-powered ones do not, and the tension is recorded rather than resolved by picking the convenient one |

---

## 12-R. ✅ EXECUTED 2026-07-31 — **the anchor holds. M-2a's verdict is ACCEPTED.**

`python experiments/diag_r1_anchor_reliability.py --run "D:\GSE165179" "D:\GSE165177"` →
`diag_r1_anchor_reliability_results.json`. **`src/` untouched.**

| check | Horvath skin & blood 2018 | Horvath multi-tissue 2013 | bar | |
|---|---|---|---|---|
| **R4** CpG coverage | **100.0%** (391) | **94.6%** (353) | ≥ 90% | ✅ both OK |
| **R1a** LODO age MAE | **6.03 yr** | **6.63 yr** | ≤ 7.17 | ✅ pass |
| **R1b** intercept-free 15-yr gap | **+10.39** (\|err\| 4.61) | **+6.48** (\|err\| 8.52) | \|err\| ≤ 9.08 | ✅ pass |
| **R1d** inter-clock ρ_partial | *(one number for the pair)* | **+0.568** | ≥ 0.50 | ✅ pass |

**⇒ §11's falsification condition is NOT met. M-2a's SPLIT ⇒ NOT CALIBRATABLE is ACCEPTED.**

Note both R1a and R1b clear only the **loosened** bars — at the originally proposed 5.0 yr, R1a
would have failed on both clocks. That loosening was frozen and committed *before* the run (§12
above), but it is a weaker test than proposed and is reported as such.

### R1c — the negative control is inert in absolute terms too

Drift of the treated non-responders against **their own day-0 fibroblast** (within-donor, so
intercept-free): **−0.76 / −0.24 yr** (skin & blood), **+2.96 / −2.56 yr** (multi-tissue). All four
are inside the clock's own error. §9-R1's worry — that non-responders had silently drifted — is
measured and does not hold here. This corroborates `REV FINAL`'s +0.5 / −2.4 yr from a different
direction: that was a *contrast against the negative-control arm*, this is an *absolute* comparison
to the donor's own starting point.

### 🔴 The finding that was not pre-registered, and that qualifies what the verdict means

R1d was included to prove the anchor could arbitrate. It does — but its **value** is the more
interesting number:

| pairing | ρ_partial | as a fraction of the meth↔meth ceiling |
|---|---|---|
| **methylation ↔ methylation** (the ceiling) | **+0.568** | — |
| RNA ↔ Horvath multi-tissue | +0.516 | **91.0%** |
| RNA ↔ Horvath skin & blood | +0.267 | 47.1% |

**Two independent methylation clocks, reading the same DNA from the same samples and sharing only
60 CpGs (17% of the smaller panel), agree with each other at ρ_partial = 0.568.** The bar was 0.50.
The ceiling clears it by 0.068.

**M-2a's bars were simulated against a null with ρ_true = 0.70. Nothing on this data reaches 0.70 —
including methylation against methylation.**

Two readings, and both are defensible:

* **Supporting the recorded verdict:** RNA ↔ skin & blood reaches 47% of the ceiling. Against that
  reference the RNA clock genuinely fails.
* **Qualifying it:** RNA ↔ multi-tissue reaches **91%** of the ceiling. Against *that* reference the
  RNA clock does nearly as well as a second methylation clock does. So the SPLIT is partly a
  property of **disagreement between the two methylation references**, not of the RNA clock alone.

> ### 🆕 ADDED on review — the ASYMMETRY is itself the evidence, and it settles the two readings
>
> *Added 2026-08-01. The table above was recorded honestly, including the reading that qualifies
> the verdict — but the inference it licenses was never drawn, and it runs in the verdict's
> favour. A reviewer who reaches the "91% of ceiling" line and stops there will reach the wrong
> conclusion.*
>
> **These are not two RNA clocks. They are ONE RNA clock (Fleischer) measured against two
> methylation references** — and those two references agree with each other at **+0.568**.
>
> | | ρ_partial |
> |---|---:|
> | Horvath-mt ↔ Horvath-sb — *the two references, to each other* | **+0.568** |
> | Fleischer RNA ↔ Horvath-mt | +0.516 |
> | Fleischer RNA ↔ Horvath-sb | **+0.267** |
>
> **If the RNA clock were tracking true biological age, it could not do this.** The two
> methylation clocks are largely measuring one shared quantity — that is what ρ = 0.568 between
> them means, on only 60 shared CpGs. Anything that tracks that shared quantity must correlate
> with **both** of them at broadly similar strength, because each reference's own reliability
> bounds how well *any* third measurement can agree with it.
>
> Instead the RNA clock agrees with one at **0.52** and the other at **0.27** — a **2× asymmetry
> against references that agree with each other at 0.57.** A measurement of the common signal
> cannot be twice as correlated with one noisy proxy of that signal as with an equally noisy
> second proxy. **The asymmetry is the signature of tracking something clock-specific rather than
> age**, which is exactly §1's diagnosis arrived at from an independent direction.
>
> **So the "qualifying" reading does not survive.** RNA↔multi-tissue reaching 91% of the ceiling
> is not evidence the RNA clock nearly works; **paired with RNA↔skin-&-blood at 47%, it is
> evidence it does not.** A verdict of SPLIT understates this — the correct reading is that the
> RNA clock fails a test the two methylation clocks pass against each other.
>
> **Recorded as an inference from already-published numbers, not a new measurement.** It changes
> no verdict — M-2a's SPLIT ⇒ NOT CALIBRATABLE stands exactly as recorded — and nothing here was
> computed after the fact to support it. What changes is only that the strongest objection to the
> verdict now has its answer written down beside it.
>
> ### 🆕 2026-08-01 — the argument above, made arithmetic
>
> *Added on cross-review. The "2× asymmetry" reasoning is right and can be stated exactly rather
> than qualitatively, which makes it checkable.*
>
> If a **single shared age factor** were all three instruments were measuring, each correlation
> would be the product of two loadings: `ρ(i,j) = λᵢ·λⱼ`. Solving the three measured values:
>
> ```
> λ_RNA  = sqrt(0.267 × 0.516 / 0.568) = 0.493
> λ_sb   = 0.267 / 0.493               = 0.542
> λ_mt   = 0.516 / 0.493               = 1.048   ← a correlation loading CANNOT exceed 1
> ```
>
> **The three numbers are not jointly consistent with one common factor.** Fitting them requires
> the multi-tissue clock to correlate with the shared signal at **1.048**, which is not a possible
> value. There is no assignment of loadings that reproduces this pattern.
>
> ⚠️ **Do not over-read it.** These are *Spearman partial* correlations at n = 68, and 1.048 exceeds
> the boundary by only 5% — comfortably inside sampling noise, so this is **at** the edge of what a
> single common factor permits, not provably past it. The honest statement is that the data sit on
> the boundary of the simplest model that could have rescued the RNA clock, and give it no room.
> That is corroboration of the paragraph above, not a proof on its own.

**What both readings agree on, and what actually matters:** a monotone calibration fitted here would
be fitted against a reference that is only ~0.57 self-consistent after the reprogramming axis is
removed. **M-2c would have been meaningless even if M-2a had passed** — which independently
vindicates §6's decision to gate it, on a ground §6 never anticipated.

### Why the ceiling is low — visible in R1a's folds

| held out | true age | predicted |
|---|---|---|
| O1 | 53 | **44.0** |
| O2 | 53 | **58.5** |
| O3 | 38 | 41.6 |

**Two donors of identical chronological age read 44.0 and 58.5** — a 14.5 yr spread on skin & blood
(43.7 vs 52.3 on multi-tissue). Donor-level methylation-clock error on this data is roughly ±7 yr,
which is the same order as the 15-yr age contrast the series contains. That is the direct cause of
the low ceiling, and it is a property of **3 donors**, not of either instrument.

### What this changes, and what it does not

| ✅ | ❌ |
|---|---|
| M-2a's verdict is now **accepted**, not merely recorded — §11's precondition is satisfied | It does **not** upgrade the verdict: NOT CALIBRATABLE still means Phase 2 does not run |
| The methylation anchor is corroborated for **contrasts** (R1c; `REV FINAL` §3 stands unchanged) | It does **not** license absolute methylation ages at n=3 donors — ±7 yr donor error |
| §6's gate on M-2c is vindicated for a second, independent reason | It does **not** rescue any RNA-derived ΔAge claim |

---

# 13. ✅ G-a and G-b CLOSED 2026-07-31 — §0's gate condition is met

Both were "small code changes" in §0 and are now implemented, tested and recorded. **These are the
only `src/` changes in this stage**; every measurement above ran with `src/` untouched.

| gate | what shipped | where |
|---|---|---|
| **G-a** | `_control_baseline` records per line: `n_control`, `n_cells`, `source` (`controls` / `self_fallback`), `unreplicated`, and the composition of the baseline **vs the whole line**. Persisted to `dataset_summary.json`; `verify_stage1_5.py` gains an unreplicated / cross-batch column | `data/aging.py`, `data/build_dataset.py`, `verify_stage1_5.py` |
| **G-b** | donor chronological age parsed (**both** GEO spellings: `donor age` and `donor age (years)`), plus `batch` derived from the title suffix — the thing D1 says nothing recorded. Both ride in `obs` as metadata | `data/sources.py` |

### The hard guard held

*"ΔAge values must come out **bit-identical** before/after. It records, it does not compute."*
Asserted in a unit test **and** re-checked on all six real Gill donors: `np.array_equal`, not
`allclose`.

### What G-a printed the first time it ran — D1 and D2, visible at last

```
N2..Y2:  n_control=1 / 19–21 cells   baseline batch=['Exp2'], line spans ['Exp1','Exp2']
```

**All six donors rest on an unreplicated baseline, and all six baselines are `Exp2` while every
donor spans both batches.** That is exactly D1 and D2 — previously reconstructible only by hand,
now emitted by the pipeline itself. G-a's purpose was never to fix them; it was to stop them being
silent, and it does.

### Two decisions worth stating

* **The new flags are reported BESIDE the Stage 1.5 verdict, never folded into it.** That PASS
  means one specific thing — the no-control fallback did not fire — and four runs are recorded
  against it. Redefining it retroactively would invalidate that record.
* **`donor_age` is metadata, never a model input.** It is a per-donor constant; it anchors
  *absolute* calibration only and cannot measure rejuvenation within a donor. The deployed request
  schema forbids extra fields, and a test pins that it stays that way.

### ⚠️ One correction to §0's own reasoning

§0 justifies gating M-2b on G-a with *"M-2b's RNA-side contrast inherits that baseline."*
**It does not.** M-2b as specified in §5 is `Δ = mean(SSEA4) − mean(CD13)` — a direct difference
between two treated arms, with no vehicle-control baseline anywhere in it. The gate's stated
rationale is wrong. G-a was implemented anyway (it is a §0 gate condition for the stage, and it is
worth having on its own merits), and M-2b ran after it, so nothing was skipped — but the reason
recorded in §0 should not be relied on by a future reader. *Per the standing rule, §0 is left as
written; this is the correction beside it.*

---

# 14. ✅ M-2b EXECUTED 2026-07-31 — **AGREE_FRAGILE, and the agreement is an artefact of the day axis**

`python experiments/diag_m2b_contrast_agreement.py --run "D:\GSE165178" "D:\Gill"` →
`diag_m2b_contrast_agreement_results.json`. `src/` untouched by this measurement.

**§10 step 1 passed first, on the actual matrices rather than titles:** join **22/22**, donors
O1/O2/Y1/Y2, days 9/11/15, 11 CD13 + 11 SSEA4, coverage **100.0% / 94.6%**. §9-R5's abort was armed
and did not fire.

| clock | sign agreement (bar ≥ 7/11) | ρ(Δ_rna, Δ_meth) | mean Δ RNA | mean Δ METH |
|---|---|---|---|---|
| Horvath skin & blood | **7/11** | +0.645 | −22.19 yr | −26.85 yr |
| Horvath multi-tissue | **7/11** | +0.491 | −22.19 yr | −26.59 yr |

⚠️ **Both land EXACTLY on the bar**, and the bar itself was already loosened from 8/11 to 7/11 by
the §6 freeze. One pair flipping fails it. Reported as `AGREE_FRAGILE`, not `AGREE`.

> ### 🆕 ADDED on review — the audit trail proving the bar was DERIVED, not chosen
>
> *Added 2026-08-01. "You loosened the bar and then landed exactly on it" is the obvious challenge
> to the paragraph above, and it is a fair one to raise. The answer is fully in the artefacts, but
> it was split across two files and never stated, so it could not be checked without reconstructing
> it by hand. Stated here instead.*
>
> | | | source |
> |---|---|---|
> | resolvability simulation ran | **13:11:39** | `stage_1_5_2_resolvability_results.json` |
> | M-2b ran | **13:53:13** — **42 minutes later** | `diag_m2b_contrast_agreement_results.json` |
> | registered bar 8/11 | **UNRESOLVABLE**, pass rate **0.9297** vs the 0.95 floor | resolvability |
> | `usable_bar` **computed** by `audit_metrics` | **7.0** | resolvability |
> | bar actually used by M-2b | **7** — *identical to the computed value* | M-2b |
>
> **The 7/11 bar is not a number someone picked. It is the output of §5b's `usable_bar`, computed
> from a simulated null before the data was touched, and frozen 42 minutes before the run.** §5b's
> instruction on an unresolvable bar is exactly *"move the threshold to `usable_bar`, or change the
> geometry, or drop the criterion — but do it now, not after a run wears the failure."* That is what
> happened, in that order, and the timestamps prove it.
>
> **What remains true and is not defended away:** the result landed *exactly* on the bar, so one pair
> flipping changes the label. That is why it is recorded as `AGREE_FRAGILE` and why §14's conclusion
> rests on the **0/3 at the discriminating timepoint**, not on the 7/11. **The fragility is real; the
> goalpost-moving is not.**

### 🔴 The pooled number hides the whole story. Split by day:

| day | agreement | mean Δ RNA | mean Δ METH (s&b) |
|---|---|---|---|
| **9** | **0 / 3** | **+38.82 yr** | **−2.97 yr** |
| 11 | 3 / 4 | −38.37 yr | −3.55 yr |
| **15** | **4 / 4** | −51.76 yr | **−68.06 yr** |

**Perfectly graded, and it inverts the headline.** At **day 15**, where methylation reports a huge
real effect (−68 yr), the two agree 4/4 — but at that magnitude *any* instrument that responds to
reprogramming at all gets the sign right. At **day 9**, the one timepoint that discriminates —
methylation says **nothing has happened yet** (−2.97 yr) — the RNA clock reports **+38.82 years of
ageing**, and agrees **0 of 3**.

**+38.82 yr against `REV FINAL` §1's +36.5 yr.** That is the identity artefact, reproduced to
within 2.3 years, on **the very samples the model trains on**, against paired ground truth from the
same cells. §4's confound warning applies to M-2b exactly as it does to M-2a: agreement
concentrated where both effects are large is agreement about the **day axis**, not about age.

### What M-2b changes: nothing about the verdict, something about the evidence

§7 row 4 was pre-committed: **M-2a fail ⇒ NOT CALIBRATABLE, whatever M-2b says.** That holds.

But M-2b was expected to *disagree* — §5: *"Disagreement is the live hypothesis."* It technically
agreed, at the bar, and **the pre-registered expectation was wrong as stated.** The finer reading is
that the pre-registered expectation was right about the *mechanism* and wrong about the *statistic*:
the modalities disagree precisely where it matters and agree where nothing could distinguish them,
so a pooled sign test was the wrong instrument for the question. That is recorded as a miss, not
reframed as a hit.

**M-2c remains NOT RUN** — gated on M-2a, which failed. §12-R added a second, independent reason:
the methylation reference is only ~0.57 self-consistent here, so a calibration fitted to it would
have no meaning either.

---

# 15. ✅ G-c STEP 1 EXECUTED 2026-07-31 — **RUN_STEP_2, and §0's evidence for G-c was wrong**

`python experiments/diag_gc_hff_signature.py --run runs/cellfate_loocv_O1` →
`diag_gc_hff_signature_results.json`. Bars frozen and committed **before** any label was read
(`d57fcd7`); resolvability **98.9% RESOLVABLE**, so no §5b move was needed. `src/` untouched.

### 🔴 §0's G-c evidence table does not survive measurement

§0 justified G-c with this row, and predicted the answer would be "no signature":

| | dose-response | monotone? | source |
|---|---|---|---|
| **§0's claim: HFF RNA labels** | **−0.36 yr/day**, ρ **−0.214** | **no** | `diag_d2_replication_results.json` |
| **measured on the actual labels** | **−1.526 yr/day**, ρ **−0.905** | **ρ says yes** | this run |

**ρ = −0.905 is stronger than methylation's own −0.885 / −0.842.** §0's ρ = −0.214 is off by a
factor of four, and its "≈9× weaker and non-monotone" verdict does not hold.

**Why the two disagree, and why this run is the right one for G-c's question.**
`diag_d2_replication` measured a **pseudobulk of 2000 sampled cells per timepoint** on **absolute
predicted age**. G-c asks about the ΔAge *labels the model trains on*: per-cell `y_age`,
**control-relative**, **after cell-cycle deconfounding**. Those are different quantities, and the
difference is not small. §0 cited the wrong one.

To §0's credit, it flagged its own evidence as *"⚠️ Not decisive on its own … which is exactly why
this is a gate with a test, not a conclusion."* **The test was the right call and it changed the
answer.** *Per the standing rule §0 is left as written; this is the correction beside it.*

### The measurement

| day | n cells | mean ΔAge | SEM |
|---:|---:|---:|---:|
| 0 | 4 983 | 0.000 | 0.183 |
| 2 | 3 947 | **+3.854** | 0.224 |
| 4 | 4 760 | **+3.538** | 0.214 |
| 6 | 4 681 | −5.798 | 0.244 |
| 8 | 4 732 | −3.978 | 0.226 |
| 10 | 4 919 | −6.311 | 0.193 |
| 12 | 4 872 | −8.227 | 0.285 |
| 14 | 4 799 | **−24.023** | 0.268 |

| criterion | measured | bar | |
|---|---|---|---|
| ρ_timepoint (monotonicity) | **−0.905** | ≤ −0.50 | ✅ **PASS** |
| slope (within ~2× of methylation) | **−1.526** yr/day | [−6.45, −1.61] | ❌ **FAIL by 0.084** |

**⇒ G-c step 1: `RUN_STEP_2`** — one criterion holds, one does not. That is the pre-registered
"ambiguous" row, and it routes to the `age_mask=True` vs `age_mask=False` retrain comparison.

⚠️ **A fourth hairline margin.** The slope misses the band edge by **0.084 yr/day**. Had the band
been "within 2×" of the *shallower* methylation figure (−3.15/2 = −1.575) instead of the mean, it
would have missed by 0.049. The verdict is right on the boundary and is reported as such.

### Robustness — descriptive, added after the first run

Leave-one-timepoint-out:

* **ρ is robust**: spans **[−0.964, −0.857]** across all eight folds. The monotone trend is not
  carried by any single point.
* **the slope is not**: spans **[−1.923, −0.938]**, and dropping **day 14 alone halves it** to
  −0.938 — comfortably outside the band in the *weaker* direction.

Day 14 is the last timepoint before the iPSC endpoint the standing rule already excludes as a
cell-type change. **So the magnitude of HFF's apparent rejuvenation depends chiefly on the point
closest to that identity change**, while its monotonicity does not. That is precisely the ambiguity
step 2 exists to resolve, and it is now a concrete reason rather than a shrug.

### The early rise, and its agreement with M-2b

Days 2–4 read **+3.9 / +3.5 yr** — the labels say the cells get *older* first. §14 found the same
sign at the discriminating early timepoint in the Gill data (day 9: RNA **+38.8 yr** while
methylation said −3.0). Two independent datasets, two different cell systems, same direction:
**the RNA clock reports early reprogramming as ageing.** The magnitudes differ by an order of
magnitude, so this is a shared direction, not a shared effect size — but it is the same failure
mode appearing twice.

### Label volume, measured

| | cells | share |
|---|---:|---|
| HFF | **42 481** | **99.71%** |
| non-HFF (Gill, all 6 donors) | 124 | 0.29% |

Consistent with §2's "≈75 of 33,688" once the train-split restriction is applied — §2 counts the
**training split**, this counts **all** cells in the built run. The ratio is the same either way,
and it is the ratio that matters: **masking HFF leaves the age head with of order 10² labels.**

### What this changes

| ✅ | ❌ |
|---|---|
| G-c step 1 is **executed**, and its answer is **not** the one §0 predicted | It does **not** show the HFF labels are good — one of two criteria failed |
| G-c step 2 (the retrain comparison) is **licensed and now necessary**, with a specific hypothesis: does the trend survive without day 14's magnitude? | It does **not** license keeping the labels by default, nor masking them by default |
| §0's cited evidence is **corrected by measurement** | It does **not** touch §7 or the M-2a verdict — G-c gates Phase 2 only, and Phase 2 does not run |

---

# 16. ✅ STAGE CLOSED 2026-07-31 — every question this document asked, answered

## 16.1 The questions, and where each is answered

| # | the question | answer | where |
|---|---|---|---|
| **M-2a** | does the RNA clock track methylation age? | ❌ **SPLIT** — ρ_partial +0.267 / +0.516 vs bar 0.50 | §11 |
| **M-2b** | do the modalities agree on the ΔAge-shaped contrast? | ⚠️ **AGREE_FRAGILE** at exactly 7/11 — and **0/3** where it counts | §14 |
| **M-2c** | can a correction be learned? | 🚫 **NOT RUN** — gated on M-2a. Two independent reasons | §11, §12-R |
| **§11 falsification** | is a negative verdict safe to accept? | ✅ **yes** — the anchor passes all four checks | §12-R |
| **§9-R1** | are the non-responders inert, or drifting? | ✅ **inert**: −0.76 / −0.24 / +2.96 / −2.56 yr vs their own day-0 | §12-R |
| **§9-R4** | is CpG coverage adequate? | ✅ **100.0% / 94.6%** on both series | §12-R, §14 |
| **§9-R5** | is the join real, or a title artefact? | ✅ **22/22** on the actual matrices; abort armed, did not fire | §14 |
| **§10 step 1** | shape before statistic | ✅ done on both pairings | §6, §14 |
| **G-a** | is the baseline visible? | ✅ **closed** — and it printed D1 + D2 on first run | §13 |
| **G-b** | is donor age wired? | ✅ **closed** — both GEO spellings, plus `batch` | §13 |
| **G-c step 1** | do HFF's labels carry the signature? | ⚠️ **RUN_STEP_2** — ρ passes, slope fails by 0.084 | §15 |
| **G-c step 2** | `age_mask` True vs False on the scorecard | ⏳ **THE ONE OPEN ITEM** | §16.4 |

## 16.2 What this stage licenses, and what it does not

| ✅ licensed | ❌ not licensed |
|---|---|
| **Closing the RNA-clock route.** Five attempts have now failed: refit, precision, control-swap, statistical fixes, calibration | **Phase 2.** It does not run |
| ΔAge from **methylation contrasts** where they exist (`REV FINAL` §3, corroborated absolutely by §12-R's R1c) | Any **RNA-derived** ΔAge claim, including the labels currently in the training set |
| Stating that Stage 2's premise **remains void** — the labels were never corrected, because they could not be | Treating the current labels as "merely noisy" |
| G-c step 2 as a **necessary** experiment with a specific hypothesis | Deciding `age_mask` without it, in either direction |

## 16.3 What went wrong in this document, recorded rather than quietly fixed

Four of this stage's own claims did not survive its own tests. Each is annotated in place; none is
rewritten.

| where | the claim | what measurement said |
|---|---|---|
| **§6** | the pre-committed fallback (ρ_partial at n=22) would rescue an unresolvable ρ_within | **it was itself UNRESOLVABLE at 92.3%.** Only GSE165177's n=68 made a valid verdict possible |
| **§5** | "disagreement is the live hypothesis" for M-2b | it **agreed**, at the bar. The mechanism was right, the pooled sign statistic was the wrong instrument |
| **§0** | G-c's evidence: HFF ρ **−0.214**, non-monotone | measured **ρ −0.905** — §0 cited a pseudobulk of absolute ages, not the labels |
| **§0** | M-2b is gated on G-a because it "inherits that baseline" | **it does not** — M-2b is a difference of two treated arms, with no vehicle baseline in it |

**And one in the execution, not the plan:** the first pass of §12's script graded R1a against the
*proposed* 5.0 yr rather than the *committed* 7.17. Both numbers are reported in §12-R.

**Four hairline margins** now sit in this project's record: E1b 0.009, D2 0.014, M-2a **0.016**,
G-c **0.084**. That is not bad luck — it is what happens when bars are set near the resolution of
the instrument. §12-R's ceiling finding says why: **nothing on this data reaches ρ = 0.70, including
methylation against methylation.**

## 16.4 The one open item — G-c step 2

> ### 🆕 2026-07-31 — **`STAGE_1_5_3_EXECUTE.md` now owns the code side of this.**
>
> *Additive; nothing in §16.4 is modified.* Writing that stage surfaced something §16.4 did not
> know: **G-c step 2 cannot be run today.** `age_mask` is a function of `source` alone
> (`aging.py:219`), and Gill and HFF both report `source = "reprogramming"` — so there is no
> expressible policy that masks HFF and keeps Gill. 1.5.3 C-1 is that blocking change, and it also
> itemises five others 1.5.2 implies, including the finding that **RES is multiplied by a
> ΔAge-derived term with no way to report that the term is unvalidated** (`res.py:38-41`).

**Not run here, and deliberately.** It needs a rebuild + retrain, and this stage's entire Phase 1
guarantee is `git diff --stat src/` staying empty for the measurements. Handing a retrain to a
measurement-only phase is how scope creep starts.

What the next stage that rebuilds must do, pre-specified so it cannot be chosen after the fact:

1. Compare `age_mask=True` vs `age_mask=False` for HFF in **one** retrain, on the existing scorecard.
2. **Pre-register the metric before the run** (§0's G-c wording), through `audit_metrics.bar_verdict`.
3. Carry §15's specific hypothesis: HFF's ρ is robust (spans [−0.964, −0.857] under
   leave-one-timepoint-out) but its **slope halves without day 14**, the point nearest the excluded
   iPSC identity change. If masking helps, that is why.
4. State the consequence either way: **masking leaves the age head of order 10² labels.** Too few is
   a finding, not a failure.

## 16.5 Verification of the stage itself (§10's own checklist)

| requirement | status |
|---|---|
| `git diff --stat src/` empty for Phase 1 | ✅ every measurement ran with `src/` untouched; the only `src/` changes are G-a and G-b (§13), which are the §0 **gate**, not Phase 1 |
| full suite green | ✅ **537 passing**, 0 failing (was 455 at the start of this stage — **+82**) |
| every pure function unit-tested with no repo data present | ✅ 4 new test files, plus 4 new bar rows in `tests/test_bars_resolvable.py` |
| `pair_by_donor_day` reused verbatim, keeping its regression test | ✅ §12-R reuses `diag_methylation_anchor.py` unmodified |
| every bar carries an entry in `tests/test_bars_resolvable.py` | ✅ added for §12 and G-c; a bar without one is not pre-registered |
| ruff clean | ✅ on `src/`, `scripts/` and every file this stage touched. *(12 pre-existing errors remain in three older `experiments/` files and one older test — untouched, not introduced here.)* |
| recorded in `CHANGES.md` + the lab notebook | ✅ |

**Deviation from §10's named artefacts, recorded.** §10 anticipated one script,
`experiments/diag_label_anchor.py`. The stage shipped **five**, because the measurements turned out
to be separable and separately gated:

| planned | actual | why |
|---|---|---|
| `diag_label_anchor.py` | `stage_1_5_2_resolvability.py` | §6's bar freeze, run first and alone so no measurement could leak into it |
| | `diag_m2a_calibratability.py` | M-2a, on GSE165177 × GSE165179 |
| | `diag_r1_anchor_reliability.py` | §11/§9-R1/§9-R4, which §10 never allocated a script to |
| | `diag_m2b_contrast_agreement.py` | §10 step 1 + M-2b, on GSE165178 × GSE165176 |
| | `diag_gc_hff_signature.py` | G-c step 1, added to this plan after §10 was written |

Each has a matching `tests/test_*.py` and a `*_results.json`.

---

# 17. 🔍 RE-AUDIT 2026-07-31 — **§11's per-arm reading was wrong, and the real finding is stronger**

> *Additive. Nothing above is modified. Every load-bearing number in §11–§16 was re-checked against
> its JSON artefact and all reproduce exactly; this section reports the one thing that did not
> survive re-reading. The verdict does not move — §7 was decided on ρ_partial, not on this table.*

## 17.1 The defect

§11 reported RNA↔methylation **per arm** and read it as:

> *"The RNA clock tracks methylation age in cells that are NOT reprogramming, and stops tracking —
> or inverts — in exactly the cells that are. That is the same in-domain/out-of-domain boundary
> `REV FINAL` §1 established, now visible in the calibration itself."*

**That table has a numerator and no denominator.** A low RNA↔methylation correlation inside an arm
means one of two completely different things, and §11 did not distinguish them:

* the methylation reference is **sharp** there and the RNA clock disagrees with it → the RNA clock
  is failing, and that is the strongest evidence this stage can produce;
* the methylation reference is **blunt** there — the two Horvath clocks do not even agree with each
  other → **nothing** can be concluded about the RNA clock in that arm, in either direction.

§12-R measured the ceiling **pooled** (ρ_partial +0.568 over all 68 conditions) and §11 measured the
numerator **per arm**. Neither was wrong; putting them side by side was never done.

`python experiments/diag_m2a_per_arm_ceiling.py` → `diag_m2a_per_arm_ceiling_results.json`.
Re-cut of the numbers M-2a already wrote — nothing re-measured, no raw data reopened.

## 17.2 The denominator, per arm

| arm | n | **meth↔meth** *(the ceiling)* | RNA (mean of both clocks) | % of ceiling | reprogramming? |
|---|---:|---:|---:|---:|:--|
| **`transient_reprogramming_intermediate`** | 11 | **+0.936** | **−0.164** | −17% | **YES** |
| `negative_control` | 12 | +0.860 | +0.399 | 46% | no |
| `failing_to_transiently_reprogram_intermediate` | 12 | +0.762 | +0.112 | 15% | no |
| `negative_control_intermediate` | 12 | +0.671 | +0.231 | 34% | no ⚠️ *too blunt* |
| `failed_to_transiently_reprogram` | 12 | +0.566 | +0.430 | 76% | no ⚠️ *too blunt* |
| **`transiently_reprogrammed`** | 9 | **+0.233** | +0.150 | 64% | **YES** ⚠️ *too blunt* |

*"Too blunt" = meth↔meth below 0.70, the level M-2a's own null assumed a real agreement looks like.*
*Unpartialled within-arm Spearman; not numerically comparable to §12-R's pooled partialled +0.568.
Conditioning on arm already removes most of the between-arm axis the partialling existed to handle.*

**Only 3 of 6 arms have a reference sharp enough to arbitrate anything.**

## 17.3 What actually follows — three corrections to §11

**① §11's headline sentence is too clean, and is withdrawn as stated.**
`failing_to_transiently_reprogram_intermediate` is a **non-reprogramming** arm with a **sharp**
reference (+0.762), and the RNA clock reads **+0.112** there — **15% of the ceiling.** It does not
"track" in that arm. So the failure is *not* confined to reprogramming cells, and the clean
in-domain/out-of-domain boundary §11 claimed to see is not in this table.

**② §11 counted an uninterpretable arm as evidence.** `transiently_reprogrammed` was listed with
+0.267 / +0.033 and read as the clock "stopping". Its meth↔methylation ceiling is **+0.233** — the
*lowest of all six arms*. The two references barely agree with each other there, so **that row
supports nothing** and should never have been read.

**③ The one row that does hold is far stronger than §11 made it look — and §11 buried it.**
In `transient_reprogramming_intermediate` the two methylation clocks agree at **+0.936, the sharpest
reference of any arm in the study**, and the RNA clock is **negatively correlated (−0.164)**.

> **Where the ground truth is at its most reliable, the transcriptomic clock runs backwards.**

That is not a power problem, not a blunt-reference artefact, and not the pooled ρ_partial that
§12-R's ceiling finding qualified. It is the single most damning result in the stage, and §11
presented it as one cell in a six-row table of equals.

## 17.4 Why this changes no verdict

| | |
|---|---|
| **§7 was decided on ρ_partial at n=68**, which is unaffected — this is a different cut of the same rows | ✅ **NOT CALIBRATABLE stands** |
| every arm's n is **9–12**, which the §6 freeze established is **UNRESOLVABLE** for a ρ bar | ✅ nothing here is a criterion, by construction |
| §11 already labelled the per-arm table "descriptive" and "no single number here is reliable" | ✅ the correction is to the **reading**, not to a result |

**But "descriptive" was doing a lot of work.** §11 labelled the table descriptive and then drew a
structural conclusion from it in the next sentence. That is the defect worth naming: *a caveat does
not license a claim.*

## 17.5 It also sharpens §12-R's ceiling finding

§12-R reported one pooled ceiling, ρ_partial **+0.568**, and concluded that no instrument on this
data reaches the ρ_true = 0.70 the bars assumed. **The per-arm cut shows that pooled number is an
average over a 4× range — +0.233 to +0.936.**

So the correct statement is not *"methylation is uniformly mediocre here"* but:

> **The methylation reference is excellent in some cell states and near-useless in others, and the
> pooled figure hides which.** A calibration fitted across all conditions — which is exactly what
> M-2c would have been — would have been fitted against a reference whose reliability varies 4× with
> the very variable being calibrated.

That is a **third** independent reason M-2c would have been meaningless, and a stronger one than
§12-R's: not merely that the reference is imprecise, but that its precision is **confounded with the
axis under study.**

## 17.6 Everything else in §11–§16 re-verified

Checked against the JSON artefacts, not against the prose:

| | |
|---|---|
| §6 geometry and all 7 bar verdicts | ✅ exact |
| §11 M-2a: ρ_all +0.444/+0.690, ρ_partial +0.267/+0.516, SPLIT, 18 pluripotency genes | ✅ exact |
| §12-R: coverage 1.000/0.946, LODO 6.03/6.63, gap +10.39/+6.48, R1d +0.568 on 60 shared CpGs (17.0%), ceiling 47.1%/91.0% | ✅ exact |
| §14 M-2b: 22/22 join, 7/11 both clocks, per-day 0/3, 3/4, 4/4 | ✅ exact |
| §15 G-c: ρ −0.905, slope −1.526, 8 timepoints, 37 693 cells, LOO ranges | ✅ exact |
| §16.5's test count | ⚠️ **537 at close; 555 now** — REV FINAL §11's donor-identity work added 18 after this stage closed. Not a discrepancy, but the figure is a snapshot |
| §16.5's "12 pre-existing lint errors" | ✅ still 12, in the same four files, none introduced here |
