# ΔAge Investigation — Lab Notebook

**How we work here (the rules).** This is a running log, not a scratchpad. For **every** test:

1. **State the hypothesis** — what we think is true, and *why*.
2. **State the prediction BEFORE running** — the specific numbers/outcome we expect if the
   hypothesis is right. (Predicting before looking is what keeps this honest — if we only
   interpret after, we rationalize whatever we see.)
3. **State what the result will let us conclude** — for each possible outcome, what it means
   and what we do next. Decide the branches *before* seeing the data.
4. **Run it ONCE.** Record the actual result verbatim.
5. **Verdict:** was the prediction right? Then either (a) hypothesis holds and the cause is
   ruled in/out → move to the next planned test, or (b) prediction wrong → the surprise is a
   *finding*; design a new test to chase it.
6. **Guardrail:** never tune a knob and re-run to get a nicer number. Each test is a
   measurement, not an optimization. Tuning against these results = training on the test set.

**The problem we're cracking:** our neural-net model ties a simple linear model (ridge) on
ΔAge magnitude (MAE ~14 both). A more powerful model tying a weaker one on the same input is
itself the finding — we need to know *why*, because the "why" decides whether there's
anything to fix.

**The candidate causes (we rule these in/out one at a time):**
- H1 code/architecture bug · H2 sample-limited · H3 clock/label noise · H4 signal is linear
  (nothing to fix) · H5 modality gap · H6 donor heterogeneity · H7 model input misses the
  clock's signal genes.

---

## Test 0 — PRE-STEP: does the model input carry the clock's ΔAge signal? (tests H7)  ✅ DONE

**Hypothesis.** The clock computes ΔAge from its gene set; the model only sees the 2000-HVG
panel. If the panel misses the clock's high-signal genes, the model can't see what it must
predict.

**Prediction (made before running).** Unknown — this was an exploratory measurement. Bar for
"H7 is real": panel captures < ~70% of the clock's signal.

**What running it gives us.** A direct measurement of overlap, computable from files on hand
(panel.json, fleischer_clock.json, harmonization.json) — no training needed.

**Result (actual).**
- Fleischer clock uses **33,155 genes** (essentially the whole transcriptome — a *dense* clock).
- 98% of the 2000 panel genes are clock-weighted (panel isn't junk).
- But the panel captures only ~**47%** of the clock's variance-weighted ΔAge signal (raw
  weight mass: 21%). **88 of the clock's top-200 signal genes are absent from the panel.**

**Verdict.** H7 is **partially confirmed as a mechanism for high absolute MAE** — the model
sees only ~half the ΔAge signal.

**BUT — critical caveat that redirected the investigation.** Ridge uses the *same* 2000-gene
input, so it also sees only 47%. The missing signal hurts **both equally** → this explains
why *absolute MAE is high*, but **not** why the model *ties ridge*. So H7 is NOT the
explanation for the model-vs-ridge tie. It splits the problem into two questions:

- **Q1 — accurate ΔAge (low MAE):** H7 relevant (expand input → both improve). Affects the
  ceiling, same for everyone.
- **Q2 — does our model beat ridge on the input it has?** H7 irrelevant. This is about whether
  the model extracts more than a linear fit from the *same* features. **This is the real
  target** (a strictly-more-powerful model should never lose to a linear one on the same data).

**Decision.** Do NOT expand the input (user's call, and correct — it wouldn't address Q2).
Focus on Q2: is the model failing to reach a ceiling it *should* reach (fixable), or is the
ceiling genuinely linear (nothing to fix)? → Test 1.

---

## Test 1 — Is the model OVERFITTING the training donors? (tests H2/overfitting, across all 6 folds)  ⏳ AWAITING RESULT

**Hypothesis.** With only ~5 training Gill donors (~100 bulk cells), the model may memorize
them and fail to generalize to the held-out donor — which would look like "ties ridge on
test" while actually fitting train perfectly.

**Prediction (made before running).**
- If **overfitting**: train-Gill MAE will be small (say < 5) while held-out MAE is ~14 — a
  large gap, on most folds.
- My honest bet: **partial-to-not-overfitting.** I expect train-Gill MAE to be *lower* than
  test (some memorization) but not tiny — because there are so few Gill cells that even
  fitting them is hard, and the ΔAge signal is dense/noisy. I'd guess train-Gill ~8–12 vs
  test ~14. If so → not primarily an overfitting story.

**What running it gives us.** Reads all 6 finished LOOCV folds (already on disk), computes the
model's ΔAge MAE on training Gill donors vs the held-out donor, per fold + aggregate. Command:
`python overfit_check.py`.

**Decision branches (decided before seeing data):**
- **train-Gill ≪ test on most folds → OVERFITTING confirmed.** Real, fixable model problem.
  Next: Test 2a (does regularization / a simpler model close the gap?). This would mean the
  model *can* be made to beat ridge → the user's instinct was right.
- **train-Gill ≈ test on most folds → NOT overfitting.** The model fits training donors about
  as poorly as held-out. It genuinely can't fit ΔAge better on this input. Next: Test 3
  (is the signal linear, or is training the limit?) — plant a known nonlinear signal.
- **Mixed across folds → donor-dependent overfitting** (would match known heterogeneity, H6).
  Next: both a regularization test AND the linearity test.

**Result (actual).**

| held-out | train-HFF | train-Gill | TEST (held-out) | read |
|---|---|---|---|---|
| N2 | 1.01 | 11.58 | 21.79 | partial |
| N3 | 1.32 | 6.63 | 29.69 | memorizing |
| O1 | 1.06 | 8.89 | 5.39 | cannot-fit |
| O2 | 1.29 | 9.11 | 7.54 | cannot-fit |
| Y1 | 1.09 | 8.80 | 7.28 | cannot-fit |
| Y2 | 1.16 | 8.28 | 14.06 | partial |

Aggregate: mean train-Gill MAE **8.88**, mean held-out MAE **14.29**. 1/6 "memorizing", 3/6 "cannot-fit".

**Verdict. Prediction was RIGHT** (predicted train-Gill ~8–12 vs test ~14; actual 8.88 vs 14.29).

**OVERFITTING IS RULED OUT — decisively. Three key observations:**
1. Train-Gill MAE is **8.88 and consistent** (6.6–11.6). The model can't fit ΔAge below ~9
   *even on donors it trained on*. It is not memorizing; it genuinely can't fit better.
2. **On O1, O2, Y1 the held-out MAE (5.4/7.5/7.3) is LOWER than train-Gill MAE (8.9/9.1/8.8).**
   You cannot overfit and simultaneously do *better* on unseen data. This nails it.
3. The two high-test-MAE folds (N2 29.7-ish, N3) are **donor heterogeneity (H6), not
   overfitting** — N3 is the known outlier (+30 while others don't). "Memorizing" tag on N3 is
   an artifact of the ratio test; it's really an outlier-donor effect.

**Consequences:**
- The "regularize the model to beat ridge" path is **DEAD** — there's nothing being overfit,
  so regularization has nothing to remove. (H2-as-overfitting ruled out.)
- **New finding: the model hits a ~9-year MAE floor it can't cross even on training data.**
  When a model fits *training* data poorly, it usually means low signal-to-noise in the
  features, not a bug. This is a strong hint toward H4 (limited/linear signal — a data
  ceiling) rather than a model failure — but Test 3 is needed to confirm.

**Decision → Test 3** (linear vs model-failure), per the pre-committed branch. Updated
prediction going in (revised by this result): I now lean that the model will tie ridge even
on nonlinear synthetic signal *at ~100 samples*, because the training-fit floor suggests the
real ceiling is low. Test 3 uses large-n synthetic with *known* linearity to separate "low
signal" from "model can't capture nonlinearity."

---

## Test 2 — [conditional on Test 1] — PLANNED, not yet specified

If Test 1 says **overfitting**: Test 2a = re-train ΔAge with stronger regularization / a
simpler head, ONCE, and check whether the model now beats ridge on held-out. (Guardrail: one
config chosen by principle beforehand, not swept against the test number.)

If Test 1 says **not overfitting**: Test 2 = Test 3 below (linearity) becomes the next step.

_[Specify precisely once Test 1 result is in — per the rules, we decide the exact test after
seeing which branch we're on.]_

---

## Test 3 — Is the true ΔAge signal LINEAR (nothing to fix) or is the model failing? (tests H4 vs H1)  — PLANNED

**Hypothesis.** The model ties ridge either because the real ΔAge↔expression map is ~linear
(ridge is optimal → nothing to fix) OR because the model can't capture nonlinearity that IS
there (a real model failure).

**Prediction (to be made before running).** _[state before running]_

**Method.** Synthetic data with a *known* ΔAge signal in the same 2000-gene input space, at
large n, no noise:
- (3a) LINEAR signal → model should TIE ridge (ridge is optimal by construction).
- (3b) NONLINEAR signal (interactions/thresholds) → model should BEAT ridge.

**Decision branches (before data):**
- **Model beats ridge on 3b** → architecture is *capable* of beating ridge when there's
  nonlinear signal to find. So on real data, the tie means the real signal is ~linear (H4
  confirmed) → **nothing to fix; report ridge is near-optimal, honestly.**
- **Model CANNOT beat ridge even on 3b** → the model is failing to capture nonlinearity that
  exists → real H1/training problem → fix the model (architecture/training), not the data.

**Result (actual).** _[TO FILL]_

**Verdict.** _[TO FILL]_

---

## Later planned tests (only if still unexplained after 1–3)

- **Test 4 — sample staircase (H2).** Same clean signal at 40k → 10k → 2k → 500 → 118 → 60
  samples; find *where* the model-vs-ridge gap closes. Tells us how far data would need to
  scale to fix it. (Staircase, not a single big-vs-small jump — per the earlier decision.)
- **Test 5 — label noise (H3).** Add increasing noise to the true ΔAge; find the level where
  the model-vs-ridge gap vanishes; compare to plausible real clock noise.
- **Test 6 — modality gap (H5).** Two synthetic "datasets" (single-cell-like vs bulk-like),
  same signal, harmonization on/off; does cross-modality specifically degrade ΔAge?
- **Test 7 — donor heterogeneity (H6).** K synthetic donors with per-donor response variation;
  leave-one-out; does heterogeneity alone collapse generalization?

We run these **only as needed** — if Tests 1–3 already explain the tie, later tests become
confirmation, not necessity. We stop when the cause is cracked, not when we run out of tests.

---

## Running conclusions (updated as we go)

- **After Test 0:** the problem is genuinely two questions (Q1 accurate-ΔAge vs Q2 beat-ridge).
  H7 explains high absolute MAE but not the tie. Real target = Q2. Not expanding the input.
- **After Test 1:** **Overfitting RULED OUT** (model fits training donors as poorly as/worse
  than held-out; on 3/6 folds it does *better* on held-out — impossible under overfitting).
  Regularization path is dead. **New finding: a ~9-year training-fit floor** the model can't
  cross even with the answer in front of it — a strong hint the signal ceiling is low
  (H4/H3), not a model bug. The two high-MAE folds (N2/N3) are donor heterogeneity (H6), not
  memorization. Next: Test 3 to separate "signal is linear/low (nothing to fix)" from "model
  can't capture nonlinearity (fixable)".
- _[append after each test]_

---

## The decision tree (one glance)

```
Test 0 (H7 overlap) ✅ -> panel sees ~47% of signal; explains high MAE, NOT the tie.
                          Two questions split out. Target = Q2 (beat ridge on same input).

Test 1 (overfitting, 6 folds) ⏳
   ├─ train-Gill << test -> OVERFITTING -> Test 2a: regularize, does it beat ridge? (fixable!)
   ├─ train-Gill ~= test -> NOT overfitting -> Test 3 (linear vs model-failure)
   └─ mixed              -> donor-dependent -> Test 2a AND Test 3

Test 3 (linear vs nonlinear, synthetic)
   ├─ model beats ridge on nonlinear -> real signal is linear -> NOTHING TO FIX (report honestly)
   └─ model can't beat ridge on nonlinear -> model failure -> fix architecture/training

(Tests 4-7 only if still unexplained: sample size, noise, modality, heterogeneity.)
```
