# STAGE 1.5.1 — REV FINAL: anchor ΔAge to an identity-independent measurement

**Supersedes as the execution plan:** `STAGE_1_5_1_CLOCK_PRECISION.md` (V1),
`STAGE_1_5_1_NEW_CHANGES.md` (review of V1), `STAGE_1_5_1_NEW_V2.md` (V2),
`STAGE_1_5_1_REVISED.md` (V3) and `STAGE_1_5_1_REVISED_REVIEW.md` (review of V3).
**All five are left byte-unmodified** and remain the audit trail.

**Depends on:** Stage 1 closed (PARTIAL); Stage 1.5 §9/§10 and the V3 review, all executed.
**Blocking for:** Stage 2, and every quantitative rejuvenation claim in Stage 5.
**Scope:** acquire one public dataset; compute methylation age; one decisive comparison. Steps 1–2
are **read-only** — `git diff --stat src/` must be empty until Step 4.

**Status:** ✅ **EXECUTED** (2026-07-26). §0–§9 are the pre-registration, left as written. §10 amends them against the real data; **§11 is the close-out: results, plus EIGHT logical holes this plan contained, each named and closed.** Read §11 first if you are reviewing.

---

## 0. Why this document exists

Four prior plans proposed four different fixes — refit the clock (V1), fix its precision (V2),
redefine the control to the non-responder arm (V3), fix the pooling/statistic errors (the review).
**All four were tested. All four failed**, and the reason is now established rather than suspected:

> **The transcriptomic clock is correctly built and correctly applied — but it is out of domain on
> reprogramming cells, and no RNA-only analysis can fix that, because every RNA route to "age" runs
> through the same clock.**

**Scoping the claim is not an option** (user decision, 2026-07-26). So the instrument problem must
be *solved*, and that requires a measurement of age that does not come from this clock.

---

## 1. What is settled — every line measured, not argued

### 1.1 The clock is sound *in its domain*

| test | result |
|---|---|
| GSE113957, its own domain, through our production path | **MAE 0.77 yr, ρ 0.99**, 100% weighted gene coverage |
| Gill day-0 adults, in the clock's fitted range [1, 96] | **+18.0 yr contrast for a 21 yr true gap, ρ +0.60** |

**It is not broken, not mis-applied, and not imprecise in the way V1/V2 assumed.**

### 1.2 It is out of domain on reprogramming cells

| evidence | value |
|---|---|
| non-responders after 11 days of OSKM | **+36.5 yr** — biologically impossible |
| mean responder age, days 7 → 9 → 11 | **78.7 → 101.3 → 77.2** — a +23 yr swing in 2 days, reversed in 2 more |
| what the clock is actually reading | corr(age, **pluripotency**) = **−0.62**; corr(age, fibroblast identity) = **+0.62** |
| the two reprogramming datasets | **opposite signs** (Gill +0.205, HFF −0.214) |
| Gill et al. themselves | *"existing transcription clocks failed to accurately predict the age of our negative control samples"* — which is why they built their own |

### 1.3 All three candidate label definitions fail

| definition | result |
|---|---|
| day-0 control, as built | responders **+8.2 yr**, CI [−20.1, +36.5] — no effect at any window |
| day-0 control, A1+A3 fixed | **same** — the two errors were real but are not the cause |
| non-responder control (V3) | the effect **is** the control arm: baseline-free, responders **−1.1 yr** (ns) while non-responders **+17.5 yr** (CI [+6.0, +28.9], 0/6 dissenting); and a **−9.7 yr gap already exists at day 7** |

The baseline-free responder CI **[−14.7, +12.5] excludes a Gill-scale −30 yr effect.** An effect
≲15 yr is not excluded — but nothing in RNA can resolve it, because the instrument has ~20 yr
sample-level noise on these cells.

### 1.4 The day-0 baseline is itself noise-dominated

`corr(day-0 baseline age, responder ΔAge) = −0.986`. Responder ages cluster (66–86, sd ≈ 7) while
the **n=1** baselines scatter (36–99, sd ≈ 21). ΔAge-vs-day-0 is largely the baseline error with its
sign flipped. **This is a consequence of clock imprecision (±12.27 yr on a single sample), so a more
precise instrument shrinks it directly.**

---

## 2. What this stage does

Acquire **GSE165179** — Gill's own multi-omic companion series: **96 Illumina MethylationEPIC
samples, same experiment and same donors** as the RNA data we already hold (GSE165176) — and use it
as the identity-independent anchor.

**Why methylation and not a fourth RNA dataset.** More RNA inherits the same out-of-domain
instrument. Methylation is a different molecular layer whose clocks are independently validated,
**it is what Gill used to establish the ~30 yr claim**, and — decisively — it is ~4× more precise:

| instrument | error | ΔAge noise (√2·err) | SNR vs a ~30 yr effect |
|---|---|---|---|
| Fleischer transcriptomic (current) | 12.27 yr | 17.4 yr | **1.7** |
| Horvath skin & blood (methylation) | **≈ 2.5–3.5 yr** | ≈ 3.5–5.0 yr | **≈ 6–9** |

That precision gain also removes the §1.4 problem: with a ~3 yr clock the n=1 day-0 baseline stops
dominating.

**Cost: a public download.** No wet lab, no new samples, no GPU.

---

## 3. THE BARS — set now, before any data is touched

### 3.0 Bar-resolvability (ground rule §5b) — checked *first*

With 6 donors and a methylation-clock error of ≈3.5 yr, the SE of a paired per-donor ΔAge is
`3.5·√2/√6 ≈ 2.0 yr`, so `t = 2.571` puts the minimum detectable effect at **≈ 5.2 yr**.

| effect | detected? |
|---|---|
| −30 yr (Gill's claim) | ✅ overwhelmingly |
| −15 yr | ✅ |
| **−5 yr** | ⚠ borderline — **report FRAGILE** |

**M-1 is RESOLVABLE.** This must be re-confirmed with the *actual* measured spread before the
verdict is read, and registered in `tests/test_bars_resolvable.py`.

### 3.1 M-1 (PRIMARY) — does methylation show rejuvenation in responders?

Per donor: methylation age of responder samples in the pre-registered window (**days 10–13**,
Gill's optimum, fixed before looking) minus that donor's own day-0 methylation age. Paired across
6 donors.

| verdict | condition | meaning |
|---|---|---|
| **PASS** | 95% CI excludes 0 **and** is negative | rejuvenation is real and measurable — **ΔAge has a valid anchor** |
| **NULL** | CI includes 0 | no rejuvenation detectable even with the good instrument |
| **CONTRADICTS** | CI excludes 0 and is positive | cells read older on methylation too |

### 3.2 M-2 — does the transcriptomic ΔAge agree with the methylation ΔAge?

Paired per-sample (and per-donor) comparison on samples present in both series.

| verdict | condition | consequence |
|---|---|---|
| **AGREES** | ρ ≥ 0.7 **and** the RNA/meth slope ∈ [0.5, 2.0] | the RNA clock is usable on these cells after calibration → **Step 3a** |
| **CALIBRATABLE** | ρ ≥ 0.7 but slope outside that range | usable *after* rescaling → **Step 3a** with an explicit correction |
| **DISAGREES** | ρ < 0.7 | RNA cannot measure this quantity here → **Step 3b** |

### 3.3 M-3 — the negative control (ground rule §10), and the direct adjudication

Non-responders, same window, under methylation.

This settles a live contradiction: **we measure +36.5 yr; Gill reports non-responders showed *"a
moderate reduction in transcription age."*** Methylation decides which is right.

| result | reading |
|---|---|
| non-responders ≈ 0 or mildly negative | **confirms the +36.5 yr was a transcriptomic artefact** — §1.2 is vindicated |
| non-responders strongly positive | they genuinely deteriorate; the V3 control was measuring something real after all |

**M-3 must be reported with M-1 in the same breath**, per §10. A claim whose negative control moves
with it is not validated.

### 3.4 Guards — all must hold

| # | Guard |
|---|---|
| G1 | sample pairing between GSE165179 and GSE165176 is **explicit and verified** (donor + day + arm), not assumed |
| G2 | the methylation clock reproduces **known chronological age on the day-0 samples** (the six donors' true ages are 0, 0, 29, 35, 53, 53) — if it cannot, it is not an anchor either |
| G3 | any verdict within **0.5 yr** of a boundary is reported **FRAGILE** |
| G4 | `git diff --stat src/` **empty** through Steps 1–3 |

**G2 is the load-bearing guard.** It is the same in-domain check that vindicated the transcriptomic
clock in §1.1, applied to the new instrument before we trust it. Note the neonatal caveat carries
over: donors at age 0 sit below most clocks' fitted range, so G2 is judged on the **adult** donors.

---

## 4. Steps

### Step 0 — Verify the data pairs *(minutes, read-only)*
Confirm GSE165179's 96 samples map to GSE165176's 124 by **donor + day + arm**. Record how many
pair. **If fewer than 4 donors have paired day-0 + peak-window responder samples, STOP** — M-1 is
not estimable and the plan must change rather than proceed on a fraction.

### Step 1 — Compute methylation ages *(hours, read-only)*
Implement **Horvath skin & blood (2018)** — fitted on fibroblasts, the right clock for this tissue —
from its published coefficients, plus **Horvath 2013 pan-tissue** as a cross-check. Apply to the EPIC
beta matrix through a documented, unit-tested path (`experiments/diag_methylation_anchor.py`).
**Run G2 before anything else.**

### Step 2 — The three measurements *(read-only)*
M-1, M-2, M-3 against §3's bars. **Predictions pre-registered in the lab notebook before this runs.**

### Step 3 — The decision fork *(no code until this resolves)*

| M-1 | M-2 | Action |
|---|---|---|
| PASS | AGREES / CALIBRATABLE | **3a — calibrate.** Fit the RNA→true-age correction on paired samples, validate leave-one-donor-out, apply to all data including HFF. **This keeps the quantitative claim and covers the 79% of labels methylation cannot reach.** |
| PASS | DISAGREES | **3b — re-target.** Methylation becomes the ΔAge source for Gill; HFF's age labels are unvalidated and its role must be decided explicitly (it keeps the fate head regardless) |
| NULL | any | **3c — the effect is not there at ~5 yr resolution.** A major, publishable finding: it would mean the rejuvenation signal Gill reported is not reproducible on these samples with an independent instrument. Escalates to Stage 4/5 |
| CONTRADICTS | any | **3d — stop and re-examine.** Two independent instruments agreeing that cells get older would overturn the project premise; treat as a bug hunt first (ground rule §6) |

### Step 4 — Implement, rebuild, rescore *(~4 h GPU)* — only on 3a or 3b
Pre-registered predictions, and the standing warning: **`y_age` changes, so the four guards move by
construction, and Stage 1's PARTIAL verdict does not carry over.**

### Step 5 — Revalidate
Re-run the existing suite (§9 clock validity, E1/E1b **split by responder status with a window
contrast**, D2). M1 should still fail — it is out-of-range extrapolation (§5.1), untouched by this.

---

## 5. What this does NOT fix

1. **The neonatal out-of-range limit stands.** GSE113957 has no samples below age 1, so N2/N3/HFF
   absolute ages remain unusable. A *data* limit, not a method one.
2. **HFF has no methylation.** 79% of age labels cannot be directly anchored. Step 3a addresses this
   by calibration; if we land on 3b, HFF's age labels are a live problem, not a solved one.
3. **n = 6 donors.** Better instrument, same donor count. Precision improves; generalisation does not.
4. **Methylation clocks are not oracles.** They have their own tissue and batch sensitivities. This
   buys *independent* evidence, which we currently have none of — not certainty.
5. **It does not make reprogramming intermediates in-domain for the RNA clock.** It measures how
   wrong that clock is there, which is what enables a correction.

---

## 6. Disposition of the earlier plans

| document | status |
|---|---|
| V1 `CLOCK_PRECISION` | **PARKED.** Premise (imprecision is the bottleneck) refuted — the clock reproduces at 0.77 yr in-domain. Its R4 was factually wrong (no cross-sample scaler exists) |
| review of V1 (`NEW_CHANGES`) | R4 refutation and C3 elimination **stand** |
| V2 (`NEW_V2`) | **PARKED**; its §5 finding — a *perfect* clock still reaches only SNR 1.88 — stands and is *why* the clock route was never sufficient |
| V3 (`REVISED`) | **PARKED as an execution plan.** Its A1/A3 discoveries **stand and are adopted** (§5 below); its non-responder control is refuted (baseline-free, the effect is the control arm) |
| review of V3 | its §5–§6 measurements are the evidence base for this document |

**Adopted from V3 unconditionally, independent of this plan:** stop pooling non-responders into the
treatment arm (A1); never test a dip with a monotonic statistic (A3). Both are already ground rules
§10/§11.

---

## 7. Discipline

- **Pre-register before running.** Predictions for Steps 1 and 2 in the notebook first.
- **One change at a time.** Only the age *anchor* moves. No model, calibration or RES changes ride along.
- **Report the negative control (§10)** with every effect claim — M-3 beside M-1, always.
- **State the shape before the statistic (§11)** — window contrast, not a monotonic trend.
- **Report FRAGILE** within 0.5 yr of any boundary.
- **Record failures** with the same prominence as successes. A NULL at M-1 is a real result.
- **Annotate, never rewrite.** All five prior plans stay byte-unmodified.

---

## 8. Deliverables

| Artefact | Purpose |
|---|---|
| `experiments/diag_methylation_anchor.py` (new, read-only) | pairing check (G2/Step 0), Horvath clock implementation, M-1/M-2/M-3. Writes `diag_methylation_anchor_results.json` |
| `tests/test_diag_methylation_anchor.py` (new) | every verdict branch, per the `verify_1a` lesson; plus the Horvath transform against published test vectors |
| `tests/test_bars_resolvable.py` (extended) | M-1's bar registered, per ground rule §5b |
| notebook + `CHANGES.md` | pre-registration, then results, then the §3 decision |

---

## 9. What to download

**GEO accession `GSE165179`** — *"Multi-omic rejuvenation of human cells by maturation phase
transient reprogramming [Transient array]"*, Illumina MethylationEPIC (GPL21145), 96 samples.

| file | why |
|---|---|
| the **series matrix** (`GSE165179_series_matrix.txt.gz`) | sample titles → donor / day / arm, for pairing (Step 0) |
| the **processed beta-value matrix** (supplementary) | the methylation levels the clock consumes |

⚠ **Check the supplementary file's format before running** — as with GSE113957, the first file GEO
offers may be the wrong representation. Beta values (0–1 per CpG, rows = `cg…` probe IDs) are what
is needed; if only IDATs are provided, that requires a preprocessing step and Step 1 must budget for it.

---

# 10. AMENDMENT (2026-07-26) — the data is in, and it is better than §2–§4 assumed

§0–§9 are left as written (pre-registration integrity). This section corrects them against the
**actual** contents of GSE165179, now downloaded and inspected. Nothing has been run beyond parsing
the metadata and one row of the beta matrix.

## 10.1 Format — a non-issue

| worry | reality |
|---|---|
| "structured very differently — a huge problem?" | **No.** The processed matrix is **comma-separated** (not tab) with a `Detection Pval` column interleaved after every sample. Both are two lines of parsing. |
| beta values or IDATs? | ✅ **Beta values.** Probe IDs `cg…`, values in [0, 1] (row 1 spans 0.0118–0.0264). Exactly what a Horvath clock consumes — **no preprocessing needed**, so §9's IDAT contingency does not apply. |
| do the columns match the series matrix? | ✅ **96 sample columns, 96 samples.** |

## 10.2 🔴 What §2–§4 got wrong

| § | assumed | actually |
|---|---|---|
| §2 | "same donors" as the RNA data | **3 donors: O1 (53), O2 (53), O3 (38).** The RNA set has 6. Overlap is **O1 and O2 only** |
| §3.1 | ΔAge = responder vs **own day-0** | the dataset carries a **proper untreated negative control at every timepoint** — a far better comparator (10.3) |
| §3.0 | resolvability for **6 donors** | recomputed below — the unit is the **(donor, day) pair**, not the donor |
| **§4 Step 0** | *"if fewer than 4 donors … STOP"* | **this gate is WRONG and would have wrongly halted the stage.** 3 donors, but **9 identity-matched pairs** — ample. The gate counted the wrong unit. |

**Step 0's gate is replaced by:** proceed if **≥ 6 identity-matched (donor, day) pairs** exist for
M-1. Measured: **9**.

## 10.3 What the dataset contains — and why it solves the core problem

96 samples, days 10/13/15/17, donors O1/O2/O3:

| type | n | identity |
|---|---|---|
| **Transiently reprogrammed fibroblast** | 13 | **fibroblast** — went through MPTR and *returned* |
| **Negative control fibroblast** | 21 | **fibroblast** — untreated, cultured in parallel |
| **Failed to transiently reprogram fibroblast** | 21 | **fibroblast** |
| Transient reprogramming *intermediate* | 12 | mid-reprogramming (identity changed) |
| Negative control *intermediate* | 12 | mid-reprogramming control |
| Failing to reprogram *intermediate* | 12 | mid-reprogramming |
| Fibroblast (day 0) | 3 | baseline |
| iPSC | 2 | pluripotent |

> **This is the identity-matched design that was impossible in the RNA data.**
> `Transiently reprogrammed fibroblast` vs `Negative control fibroblast` compares **two fibroblast
> populations** — same identity, same culture time, differing only in treatment. It removes the
> confound that broke every RNA analysis (`corr(age, pluripotency) = −0.62`) **by design rather than
> by adjustment**.

**And it supplies the negative control the RNA data does not have.** The review §2.4 verified that
GSE165176 contains **no untreated sample at any day > 0**. GSE165179 has **21**. This is Gill's
actual comparison, now available to us.

**Bonus:** the *intermediate* samples measure the identity artefact directly on methylation — the
same cells mid-reprogramming — a clean read on whether the +36.5 yr RNA reading is
instrument-specific.

## 10.4 Corrected design and resolvability

Unit of analysis: the **(donor, day) pair**, matched between arms.

| measurement | comparison | pairs | SE | **min detectable effect** |
|---|---|---|---|---|
| **M-1** | transiently reprogrammed *vs* negative control | **9** | 1.41 yr | **3.3 yr** |
| **M-3** | failed to reprogram *vs* negative control | **12** | 1.22 yr | **2.7 yr** |

(assuming a 3.0 yr methylation-clock error; to be re-confirmed against the measured spread.)

**Both RESOLVABLE with room to spare** — a Gill-scale −30 yr effect is detectable ~10× over, and
even −5 yr clears the bar. **This is the first time in the entire Stage 1.5 arc that the instrument
is comfortably sharper than the effect.**

## 10.5 What is weaker than §2 hoped — stated plainly

- **M-2 (RNA ↔ methylation agreement) is thin.** Only **O1 and O2** appear in both datasets, and the
  day grids differ (RNA 7/9/11/13/15/21/29 vs methylation 10/13/15/17), overlapping at **13 and 15**.
  M-2 therefore rests on ~2 donors × 2 days. **Report as indicative; never use as a gate.** §3.2's
  ρ ≥ 0.7 threshold must be read with that n in mind.
- **Consequence for HFF.** §4 Step 3a assumed M-2 could license calibrating the RNA clock for the
  79% of labels methylation cannot reach. With ~4 pairs that licence is weak. **Step 3a is downgraded
  from "calibrate and apply broadly" to "estimate a correction and report its uncertainty
  honestly"**; if the correction is not well determined, **Step 3b** (methylation as the target for
  Gill, HFF's age labels flagged) is the honest landing.
- **The two datasets are companion experiments, not the same samples** — different donors, day grids
  and sorting. M-1/M-3 stand alone; only M-2 needs the pairing.

## 10.6 Net

**M-1 and M-3 — the questions that decide whether ΔAge has a valid anchor — are fully answerable, on
an identity-matched design, with an instrument sharper than the effect.** That is a better position
than any point in Stage 1.5 so far.

**M-2 — whether the RNA clock can be rescued for HFF — is under-powered and must not be oversold.**

Next action unchanged in spirit, corrected in detail: implement the Horvath clock (**guard G2**
first — now judged on donors aged 53/53/38, all comfortably in range, so G2 is a fair test), then
run M-1 and M-3.

---

# 11. EXECUTED — results, and **eight holes in §0–§10 of this plan**, each closed

§0–§10 are left as written. This section is the honest close-out: what was run, what it found, and
**every logical hole the execution exposed in my own pre-registration**. Written so a reviewer can
audit it without me.

## 11.1 Status of every registered measurement

| ID | registered as | status |
|---|---|---|
| **G2** | clock reproduces known chronological age | ✅ run — REPRODUCES on both clocks, **but see Hole 4** |
| **M-1** | "responder vs control" rejuvenation | ✅ run — **but the plan's definition was ambiguous; see Holes 1–3** |
| **M-3** | negative control | ✅ run — **inert on both clocks**, the cleanest result here |
| **M-2** | RNA ↔ methylation agreement | ❌ **NOT RUN — and it is ill-defined; see Hole 8** |

## 11.2 The results

Two clocks (Gill used several and reported both of these rejuvenated), three contrasts, all
identity-matched within a (donor, reprogramming-length) row:

| contrast | Horvath skin & blood | Horvath multi-tissue |
|---|---|---|
| **intermediates** (cells *still* reprogramming) | **−24.1** [−31.1, −17.0] REJUVENATION | **−27.5** [−33.7, −21.4] REJUVENATION |
| **returned fibroblasts** (MPTR — Gill's claim) | −5.8 [−19.5, +7.9] NO_EFFECT | −9.4 [−18.3, −0.5] REJUVENATION_**FRAGILE** |
| **failed to reprogram** (NEGATIVE CONTROL) | **+0.5** [−2.3, +3.2] | **−2.4** [−5.7, +0.8] |

By reprogramming length, the returned fibroblasts peak at **13 days** on both clocks
(−14.1 / −18.4) and diminish at 15–17 — **exactly the shape Gill describe**. The intermediates run
the opposite way (−14 → −36 monotonically), which is simply "closer to iPSC".

**Gill's ~30 yr is reproduced on the intermediates. The negative control is inert on both clocks, so
the transcriptomic +36.5 yr artefact is dead twice over.**

## 11.3 The eight holes — each stated, then closed

### 🔴 Hole 1 — the decision fork (§4 Step 3) **cannot be evaluated**

It is a table of `M-1 × M-2`. But M-2 was never run (Hole 8), and "M-1" turns out to name **two
different contrasts with opposite verdicts** (intermediates REJUVENATION, returned fibroblasts
NO_EFFECT/FRAGILE). **As written, the plan cannot tell anyone what to do next.**

**Closed by 11.4**, which restates the fork on quantities that actually exist.

### 🔴 Hole 2 — the intermediate arm was never a registered measurement

§10.3 calls it a *"bonus"*. It became the headline. **That transition must be declared, not glossed:**
the intermediates were pre-registered only as a read on the identity artefact, and they are where the
effect was found. **Anyone reviewing this should treat the intermediate result as
hypothesis-generating-then-confirmed-on-a-second-clock, not as a pre-registered primary endpoint.**
It is honest because (a) it was named in the plan before running, (b) both clocks agree, (c) the
negative control is inert, and (d) the day-profile distinguishes it from the returned-fibroblast
contrast on an *a priori* basis. It would be dishonest to present it as the registered primary.

### 🔴 Hole 3 — §3.1's definition of M-1 is stale and would mislead

§3.1 says *"minus that donor's own **day-0** methylation age"* over *"**6** donors"*. Neither is what
was done or what should be done: the dataset carries **matched untreated negative controls**, which
are strictly better, and there are **3 donors / 9–12 (donor, day) pairs**. §10.2 flagged this but did
not rewrite §3.1's verdict table, so a reader following §3.1 literally would run the wrong analysis.
**Superseded: the comparator is the matched negative control, and the unit is the (donor, day) pair.**

### 🔴 Hole 4 — G2 is partly self-fulfilling, which §3.4 did not anticipate

The coefficient tables ship **no intercept row**, so an intercept was derived from the three
known-age day-0 samples. G2's **MAE is therefore partly circular** and "REPRODUCES at MAE 4.0/4.4 vs
a 5.0 tolerance" is **not** strong evidence on its own.

**Closed:** G2 must be judged on (i) the **spread** of the implied intercept across donors (5.3 yr
skin & blood, 6.4 yr multi-tissue — mediocre), and (ii) the **intercept sweep**, which is what
actually carries the result: over −0.60 to +0.70 the intermediates stay REJUVENATION and the negative
control stays inert, while **the returned-fibroblast verdict flips** between NO_EFFECT and
REJUVENATION_FRAGILE. Any future use of this clock must report the sweep.

### 🔴 Hole 5 — Step 3's every "PASS" row depends on M-2, which does not exist

So Steps 3–4 were unreachable as written. Closed by 11.4.

### 🟠 Hole 6 — "same donors" (§2) is false

3 donors (O1=53, O2=53, O3=38) vs the RNA set's 6; overlap **{O1, O2}** *by label*, and **it is not
verified that these are the same physical donors** — only that both are labelled O1/O2 and aged 53.
§10.2 corrected the count but not this identity assumption. **Treat donor identity across the two
series as unverified.**

### 🟠 Hole 7 — no bar was ever set for the actual scientific question

The live question is now *"how much rejuvenation survives the return to fibroblast identity?"* — the
returned-fibroblast contrast. The plan has **no pre-registered bar for it**, so the −5.8 / −9.4
result cannot be graded. **It must be pre-registered before any further analysis of that contrast**,
per ground rule §5b, and its resolvability checked at n=9.

### 🔴 Hole 8 — M-2 is not merely under-powered, it is **ill-defined**

§3.2 assumed the two datasets could be paired sample-to-sample. Measured, they cannot:

| | RNA (GSE165176) | methylation (GSE165179) |
|---|---|---|
| arms | `Reprogramming fibroblast`, `Failing to reprogram`, `iPSC`, `Dermal fibroblast` | `Transiently reprogrammed fibroblast`, **`Transient reprogramming intermediate`**, `Negative control` (×2), `Failed`(×2), `Fibroblast`, `iPSC` |
| untreated control at day > 0 | ❌ **none** | ✅ 21 + 12 |
| days | 7, 9, 11, 13, 15, 21, 29, 34, 40, 47 | 10, 13, 15, 17 |

Three independent blockers:
1. **The RNA labels do not distinguish "still reprogramming" from "returned to fibroblast."** The
   methylation data does, and the two give **−24 vs −5.8** — so the RNA arm cannot be mapped to a
   methylation arm without choosing the answer.
2. **The contrasts differ.** RNA can only reference day 0 (no untreated control exists); methylation
   references matched untreated controls. Correlating them compares two different quantities.
3. Overlap is **2 donors × 2 days**.

**Closed: M-2 as specified is withdrawn.** It cannot be executed honestly on these two datasets.
Consequence: **§4 Step 3a's promise to "calibrate the RNA clock and apply it to HFF" has no
evidential basis and must not be attempted on this data.**

## 11.4 The decision fork, restated so it can actually be evaluated

Replaces §4 Step 3. Keyed on quantities that exist.

| condition | met? | action |
|---|---|---|
| negative control inert on both clocks | ✅ **yes** (+0.5, −2.4) | the design is valid; the RNA +36.5 artefact is **closed** |
| rejuvenation detectable on an identity-matched methylation contrast | ✅ **yes**, intermediates −24.1 / −27.5, both clocks, intercept-robust | **ΔAge has a valid anchor. The target is not the problem — the RNA instrument is.** |
| rejuvenation *retained* after return to fibroblast identity | ⚠️ **unresolved** — right sign, Gill's exact day-shape, but ns on one clock, FRAGILE on the other, and intercept-sensitive | **needs a pre-registered bar + more donors. Do not claim it.** |
| RNA ↔ methylation calibration feasible | ❌ **no** — ill-defined (Hole 8) | **Step 3a is off the table.** HFF's age labels cannot be anchored from this data |

**Therefore the licensed next step is neither 3a nor 3b as written.** It is:

> **Use methylation as the ΔAge source where methylation exists (Gill's 3 donors), and treat HFF's
> age labels as unanchored** — with the fate head unaffected. Any attempt to extend the age target to
> HFF requires methylation data for HFF, which is a Stage 6 acquisition, not an analysis choice.

## 11.5 Discipline audit — pre-registered vs discovered

Stated plainly so a reviewer can judge the epistemics rather than reconstruct them:

| item | status |
|---|---|
| G2, M-1, M-3, their verdict tables, the FRAGILE rule | **pre-registered** (§3, before any data) |
| the (donor, day) pair as the unit; ≥6 pairs gate | **amended §10.4 after seeing metadata only**, before any age was computed |
| running **both** clocks | **not** pre-registered — adopted from Gill's published methods after reading them. Both are always reported |
| the **intermediate** contrast | named as a "bonus" in §10.3 pre-run; **promoted to headline after the fact** (Hole 2) |
| the intercept sweep | **not** pre-registered — added because the missing intercept was discovered mid-run |
| replicate averaging | **bug fix** — the first run silently dropped 6 of 9 pairs |

**Nothing was selected on its p-value**: all three contrasts and both clocks are reported in every
table, including the one that undercuts the retention claim.

## 11.6 What a reviewer should check hardest

1. **Is promoting the intermediate contrast legitimate?** (Hole 2.) My argument: named pre-run, two
   independent clocks agree, negative control inert, and the day-profile separates it from the
   returned-fibroblast contrast on a priori grounds. A reviewer may reasonably want this
   re-registered and re-run on new donors before it enters a manuscript.
2. **Is the implied intercept acceptable?** (Hole 4.) The sweep says the two robust conclusions do
   not depend on it. Someone should nevertheless obtain the **published** Horvath intercepts and
   confirm.
3. **Is donor identity across GSE165176/165179 real?** (Hole 6.) Assumed from labels, never verified.
4. **Is the retention question worth n=9?** (Hole 7.) I think it needs more donors; the plan should
   say so before anyone spends GPU on Step 4.
