# STAGE 1.5.2 — Is the transcriptomic clock CALIBRATABLE? (paired RNA ↔ methylation)

**Status:** 🔵 **PRE-REGISTERED — NOT YET EXECUTED, AND GATED.** Bars are fixed in §6 *before*
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
