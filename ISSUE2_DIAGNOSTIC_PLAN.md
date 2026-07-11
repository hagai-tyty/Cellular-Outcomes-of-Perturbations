# Issue 2 — Diagnostic Plan: Why can't the model predict ΔAge better than a simple linear baseline?

**Status:** planning (no tests built yet).
**Goal of this document:** decide *what* to test and *how*, so the experiments pinpoint the
**specific cause** of the ΔAge problem — not just produce more numbers.

---

## 1. The precise problem (state it exactly)

From the leave-one-donor-out cross-validation:

| metric | model | ridge (simple linear) |
|---|---|---|
| ΔAge MAE (years) | 14.3 ± 9.7 | 14.1 ± 8.0 |
| Ranking (Spearman) | 0.69 ± 0.10 | — |

**The problem is NOT "the model is broken."** The problem is:

> On ΔAge *magnitude*, our 256-dim neural-net ensemble ties a trivial linear formula.
> A big model tying a tiny model means we are **not extracting more signal than a linear
> fit does** — either because the signal isn't there (data-limited) or because the model
> can't capture it (model-limited). **We do not yet know which.**

The purpose of this plan is to determine **which**, with evidence, and only then decide
whether there is anything to "fix."

**Definition of "fixed":** either (a) we prove the tie is *expected* given the data
(nothing to fix — report honestly), or (b) we find a *specific, demonstrated* cause with a
*specific fix* that measurably improves ΔAge MAE without overfitting the test set.

---

## 2. Guardrails (the rules we will NOT break)

These matter more than the tests themselves. Violating them produces fake results.

1. **Run each experiment ONCE, read it, do not tune against it.** The moment we tweak a
   knob, re-run, and keep the better number, we are training on the test set. Every rung is
   a *diagnosis*, not an *optimization*.
2. **Every rung has a positive control** — a version the model MUST pass. If it fails the
   thing it should trivially pass, the bug is in our code/architecture, and we stop and fix
   that before going further.
3. **Change ONE variable per rung.** If two things change and performance drops, we've
   learned nothing about which one caused it.
4. **The key metric is the model-minus-ridge gap**, not absolute MAE. We are diagnosing
   *why the model doesn't beat the linear baseline*, so the comparison is the measurement.
5. **No claims from a single fold.** Where donor structure matters, use the leave-one-out
   machinery, report mean ± std.

---

## 3. Candidate causes (the hypotheses we're distinguishing between)

Each is a concrete, testable explanation for the tie. The ladder is designed so that each
rung rules one or more of these in or out.

| # | Hypothesis | If true, the fix is… |
|---|---|---|
| H1 | **Architecture/code bug** — the model can't learn a clean ΔAge mapping even on easy data | fix the model/code (real bug) |
| H2 | **Sample-limited** — ~118 bulk points is too few for a neural net to beat linear | more data (v2); not a model fix |
| H3 | **Label/clock noise** — the ΔAge target itself is noisy (clock is imperfect) | denoise, or a better/different clock |
| H4 | **Signal is linear** — the true ΔAge↔expression map is ~linear, so ridge is already optimal | nothing to fix; report ridge is near-optimal |
| H5 | **Modality gap** — bulk↔single-cell difference specifically degrades ΔAge | improve harmonization |
| H6 | **Donor heterogeneity** — donors differ so much that leave-one-out can't transfer | more donors / donor-adaptation |
| H7 | **Feature/target mismatch** — the 2000-gene model input doesn't carry the clock's genes, so the model literally can't see what it must predict | put clock genes in the panel (**easy fix**) |

**Note on H7:** this is the one hypothesis that is both *plausible* and *cheaply fixable*.
The clock (Fleischer) reads a specific gene set; the model input is the 2000-HVG panel. If
those don't overlap well, the model is being asked to predict something from features that
don't contain it. **We should check this first — it's a 10-minute check, not a full rung.**

---

## 4. The method: a synthetic "ground-truth ΔAge" difficulty ladder

**Core idea:** build fake datasets where we *know the true ΔAge analytically* (because we
plant it), start from trivially easy, and add ONE difficulty per rung. Measure where a
known-good model first fails. **The rung where it breaks is the cause.**

Why synthetic: on real data we can never separate "model is bad" from "signal isn't there,"
because we don't know the true ΔAge — the clock is our only (imperfect) estimate. On
synthetic data we *define* the true ΔAge, so any failure to recover it is unambiguous.

**What we control per rung:** sample count, label noise, signal linearity/complexity,
modality gap, donor heterogeneity, feature/target overlap.

**What we measure each rung:**
- ΔAge MAE (absolute error, years)
- R² / correlation (fraction of ΔAge variance the model explains)
- **model-minus-ridge gap** (the diagnostic metric)
- (where relevant) leave-one-donor-out Spearman

---

## 5. The ladder (rung by rung)

Each rung lists: **what changes**, **what we measure**, and **what each outcome tells us**.

### Rung 0 — Positive control (must pass)
- **Setup:** synthetic data, true ΔAge = a clean **nonlinear** function of a few input
  genes, huge sample size (~40k), no noise, single modality, clock genes fully in the input.
- **Measure:** model MAE ≈ 0, R² ≈ 1, and model **beats** ridge (because the signal is
  nonlinear and the model has capacity + data).
- **Outcome logic:**
  - PASS → architecture is sound; H1 ruled out. Continue.
  - **FAIL → H1 confirmed (real code/architecture bug). STOP and fix this first.**

### Rung 1 — Cut the sample size (tests H2)
- **What changes vs Rung 0:** drop samples from ~40k to ~118 (matching Gill). Nothing else.
- **Measure:** does the model-minus-ridge gap collapse toward 0?
- **Outcome logic:**
  - Gap collapses (model now ties ridge) → **H2 confirmed: sample-limited.** This is the
    predicted result. Fix = more data, not a better model.
  - Gap stays large (model still beats ridge at n=118) → H2 ruled out; sample size is not
    the bottleneck. Continue.

### Rung 2 — Add label/clock noise (tests H3)
- **What changes vs Rung 0:** keep large n and the clean signal, but add Gaussian noise to
  the true ΔAge target at increasing levels (e.g. σ = 0.5, 1, 2, 4 "years").
- **Measure:** at what noise level does the model-minus-ridge gap vanish?
- **Outcome logic:**
  - Gap vanishes only at high noise → clock noise *alone* isn't the explanation (real clock
    noise would have to be shown to be that high).
  - Gap vanishes at low noise → **H3 plausible: ΔAge is noise-limited.** Fix = better clock.

### Rung 3 — Linear vs nonlinear signal (tests H4 — the most important rung)
- **What changes:** two variants at large n, no noise:
  - **(3a) LINEAR** true ΔAge = a linear combination of genes.
  - **(3b) NONLINEAR** true ΔAge = interactions/thresholds between genes.
- **Measure:** model-minus-ridge gap in each.
- **Outcome logic:**
  - 3a: model ties ridge (**by design** — ridge is optimal for linear signal).
  - 3b: model **beats** ridge (it can capture nonlinearity ridge can't).
  - **This tells us what a tie *means* on real data:** if real ΔAge behaves like 3a, the
    tie is *expected and correct* (H4 confirmed — nothing to fix). If real ΔAge is like 3b
    but the model still ties → the model is failing to capture nonlinearity (points back to
    H1/H2). *This is the rung that decides whether "fixing" is even the right frame.*

### Rung 4 — Modality gap (tests H5)
- **What changes:** two synthetic "datasets" — one single-cell-like (many noisy cells), one
  bulk-like (few averaged samples) — sharing the *same* underlying ΔAge signal but on
  different scales. Run with harmonization ON and OFF.
- **Measure:** ΔAge MAE with/without harmonization; does the cross-modality case degrade
  vs a single-modality case with the same total signal?
- **Outcome logic:**
  - Cross-modality degrades even with harmonization ON → **H5: the modality gap is a real
    ΔAge bottleneck** the current harmonization doesn't fully solve. Fix = better harmonization.
  - No degradation → H5 ruled out; harmonization is doing its job for ΔAge.

### Rung 5 — Donor heterogeneity (tests H6)
- **What changes:** create K synthetic "donors" where the ΔAge *response* varies per donor
  (some rejuvenate, some age — mimicking N2 vs N3). Then leave-one-donor-out.
- **Measure:** leave-one-out MAE/Spearman vs a homogeneous-donor control.
- **Outcome logic:**
  - Heterogeneous case collapses under leave-one-out → **H6: donor heterogeneity limits
    generalization** (matches the real N2/N3 opposite-sign observation). Fix = more donors.
  - Holds up → H6 ruled out.

### Rung 6 — Feature/target gene overlap (tests H7 — check this FIRST as a cheap pre-step)
- **What changes:** make the clock genes progressively **absent** from the model's input
  panel (100% overlap → 50% → 0%).
- **Measure:** ΔAge MAE as overlap drops.
- **Outcome logic:**
  - MAE degrades as overlap drops → **H7: the panel must include clock genes.** Fix = add
    them to the panel (**easy, real, likely worth doing regardless**).
  - No effect → H7 ruled out.
- **Pre-step (do before building the full ladder):** just *check the actual overlap* between
  the current 2000-HVG panel and the Fleischer clock genes on the real data. If it's low,
  we may have found a real, cheap fix before running a single synthetic experiment.

---

## 6. The decision tree (diagnosis → fix)

Run in this order; the **first** rung that breaks a known-good model is the dominant cause.

```
Pre-step: clock-gene ∩ panel overlap on real data
   └─ low overlap?  → likely H7. Try the easy fix (add clock genes), re-measure ONCE.

Rung 0 (positive control)
   ├─ FAIL → H1 (code/architecture bug) → fix the model. [STOP HERE]
   └─ PASS → architecture sound, continue

Rung 1 (cut n to 118)
   └─ gap collapses → H2 (sample-limited) → fix = more data (v2), OR accept + report honestly

Rung 2 (add noise)
   └─ gap vanishes at realistic noise → H3 (clock noise) → fix = better/denoised clock

Rung 3 (linear vs nonlinear)  ← decides if "fixing" is even the right frame
   ├─ real ΔAge ≈ linear (3a-like) → H4 → NOTHING TO FIX; ridge is near-optimal, report it
   └─ real ΔAge nonlinear but model still ties → back to H1/H2

Rung 4 (modality gap)
   └─ degrades with harmonization ON → H5 → fix = better harmonization

Rung 5 (donor heterogeneity)
   └─ leave-one-out collapses → H6 → fix = more donors / donor-adaptation

Rung 6 (feature overlap)
   └─ degrades as overlap drops → H7 → fix = clock genes in the panel
```

Multiple causes can contribute. The **order of collapse** tells us the *dominant* one — the
one worth fixing first.

---

## 7. What we need to build

1. **A controllable synthetic ΔAge generator.** The existing `SyntheticSource` plants a
   *fate* signal but not a controllable *ΔAge*. We need a generator where the true ΔAge is a
   **known function** of input genes, with knobs for:
   - sample count
   - label noise level
   - signal linearity (linear vs nonlinear)
   - number of "datasets" (modalities) and their scale gap
   - number of "donors" and per-donor response variation
   - overlap between clock genes and the model-input panel
2. **A thin experiment runner** that, for each rung, builds the synthetic data, trains the
   model + ridge, and prints the one diagnostic line (MAE, R², model-minus-ridge gap). Reuse
   the existing training/eval machinery; only the data generation is new.
3. **A results table** (one row per rung) so the decision tree can be read at a glance.

**Reuse, don't rebuild:** the model, training loop, ridge baseline, and eval metrics already
exist and are tested. This is a *data-generation + orchestration* task, not new modeling.

---

## 8. Honest prediction (so we can be held to it)

**My bet: Rung 1 breaks it — the problem is sample size (H2), possibly compounded by H6
(donor heterogeneity).** I expect Rung 0 to pass (architecture is fine), and the model to
start tying ridge the moment we cut to ~118 samples. If so, the honest conclusion is:
*"ΔAge magnitude is data-limited, not model-limited; the fix is more donors, and for the
preprint we report the tie honestly."*

**But I could be wrong, and that's the point.** Three ways I could be wrong, each of which
changes the plan:
- If **Rung 0 fails**, there's a real bug and I'm wrong that the architecture is sound.
- If the **pre-step** shows low clock/panel overlap, there's a cheap real fix (H7) and the
  tie was partly an own-goal.
- If **Rung 3** shows real ΔAge is linear (H4), then the tie is *correct* and there was
  never anything to fix — which would save us from chasing a non-problem.

The ladder exists to rule my prediction in or out with evidence, not to confirm it.

---

## 9. Open questions to resolve before building

1. Do we run the cheap **pre-step (clock ∩ panel overlap)** first? (Recommended: yes — it
   could find a real fix or rule out H7 in minutes.)
2. What sample sizes define "large" (Rung 0) vs "small" (Rung 1)? (Proposed: 40k vs 118 to
   match real HFF vs Gill.)
3. Which nonlinearity do we use for Rung 3b — interactions, thresholds, or both? (Affects how
   strong the "model should beat ridge" signal is.)
4. Do we need Rungs 4–6 if Rung 1 already explains the tie? (Possibly not — if H2 is
   confirmed and sufficient, later rungs become confirmation rather than necessity.)
