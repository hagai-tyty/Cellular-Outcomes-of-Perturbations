# STAGE 1.5.1 — REV FINAL: anchor ΔAge to an identity-independent measurement

**Supersedes as the execution plan:** `STAGE_1_5_1_CLOCK_PRECISION.md` (V1),
`STAGE_1_5_1_NEW_CHANGES.md` (review of V1), `STAGE_1_5_1_NEW_V2.md` (V2),
`STAGE_1_5_1_REVISED.md` (V3) and `STAGE_1_5_1_REVISED_REVIEW.md` (review of V3).
**All five are left byte-unmodified** and remain the audit trail.

**Depends on:** Stage 1 closed (PARTIAL); Stage 1.5 §9/§10 and the V3 review, all executed.
**Blocking for:** Stage 2, and every quantitative rejuvenation claim in Stage 5.
**Scope:** acquire one public dataset; compute methylation age; one decisive comparison. Steps 1–2
are **read-only** — `git diff --stat src/` must be empty until Step 4.

**Status:** PLANNED. This document is the pre-registration. Nothing in §4 has been executed.

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
