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

### Test 3 RESULT (run once, refined with smooth-nonlinearity + train-MAE)

**Prediction going in (revised after Test 1):** "model will tie ridge even on nonlinear at
this sample size." **This prediction was WRONG** — the model beat ridge decisively on smooth
nonlinearity. Being wrong here is the finding: the architecture is more capable than I bet.

| signal | predict-mean | ridge MAE | model MAE | model-train | read |
|---|---|---|---|---|---|
| linear | 11.84 | **0.00** | 1.33 | 0.73 | ridge optimal, model close |
| smooth | 12.62 | 9.09 | **2.31** | 0.88 | **model WINS (+6.78)** |
| products | 11.77 | 11.25 | 11.33 | 1.31 | tie (memorizes train, can't generalize) |
| noise | 12.04 | 12.38 | 13.20 | 1.39 | both fail (control passes) |

**Reproduced independently on the user's GPU hardware** (numbers match the sandbox run within
training noise: smooth gap +6.78 vs +6.64; model-smooth MAE 2.31 vs 2.45). The result is
real and reproducible, not a one-machine fluke.

**Positive controls all pass:** linear→ridge optimal; noise→both at floor (harness invents no
signal). **Key result: model WINS on smooth nonlinearity (2.31 vs 9.09).**

**Verdict — the architecture is NOT broken.** It captures smooth nonlinearity and beats ridge
by a wide margin when such signal exists. The products-tie is the well-known universal MLP
limitation (memorizes train 1.37, fails test 11.28), not specific to this model.

### THE STRUCTURAL INSIGHT THAT CLOSES THE CASE

The Fleischer clock is a **LinearClock: age = w·x + b**. Therefore:

    ΔAge = clock(x_pert) − clock(x_ctrl) = w·(x_pert − x_ctrl)   -- LINEAR in expression.

Every transform in between (harmonization Z-score, Gill projection) is **affine**. So the
entire ΔAge target is **linear in the model's input, by construction.** There is no nonlinear
signal to capture → ridge is optimal → the model can only tie it. Test 3 proves the model
*would* beat ridge if nonlinearity existed (smooth +6.64); it doesn't on real data because
**the linear clock guarantees a linear target.** **H4 CONFIRMED, two independent ways**
(architecture-is-capable proof + structural linearity of the clock).

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
- **After Test 3: CRACKED.** Model beats ridge on smooth nonlinearity (+6.64) → architecture
  is fine. And the clock is **linear**, so ΔAge is **linear in expression by construction** →
  ridge is the optimal tool → the model tying ridge is the **CORRECT, EXPECTED** result, NOT
  a failure. **H4 confirmed. There is nothing to "fix" about the model's ΔAge prediction.**
- **After Test 5.0 (reopened by user):** the apparent "worse than ridge" (14.29 vs 14.05) is
  **statistically NOT distinguishable from zero** (paired 95% CI [−2.16, +2.64], model wins
  2/6 folds). The model is **statistically TIED** with ridge — now demonstrated with a paired
  test, not assumed. This closes the "why worse, not tied" gap: there is no real deficit.
- **After Test 6 (empirical close, user ran on real data):** NO nonlinear model beats ridge on
  real ΔAge (ridge 14.05 is best; trees 15.63, forest 17.81, kernel 24.68 — all worse, and
  worsening with flexibility). This is the signature of a genuinely linear signal and rules
  out the last open worry (nonlinear per-donor-offset inference). **Real ΔAge IS linearly
  predictable from this input — confirmed empirically, not by structural assumption.** The
  neural-net tie with ridge is CORRECT.

## FINAL DIAGNOSIS (the whole picture, consistent)

The "ΔAge problem" was **two real effects that are both understood, neither of which is a
model bug:**

1. **The model ties ridge on ΔAge magnitude** — because a **linear clock** makes ΔAge a
   **linear** function of expression. A linear model is therefore optimal; a neural net can at
   best match it. **This is correct behavior, not a failure. Do not try to "beat ridge" here.**
2. **Absolute ΔAge MAE is ~9–14 years** — because the model input (2000 HVG) captures only
   ~47% of the clock's genes (Test 0), and ΔAge = w·(x_pert − x_control) has a per-donor
   control offset the model must infer (donor heterogeneity, Test 1). Both limit accuracy for
   *everyone* (model and ridge alike).

**Implications:**
- For the preprint: report honestly that ΔAge magnitude is well-served by a linear predictor
  (the clock is linear), and the model's contribution is the **joint calibrated fate +
  ranking + uncertainty**, not beating ridge on ΔAge. The ranking result (Spearman 0.69±0.10)
  stands as the real ΔAge-related contribution.
- If lower ΔAge MAE is ever wanted (Q1): the levers are **more clock genes in the input**
  (H7) and **more donors** (H6) — NOT a fancier model. And a more powerful ΔAge model is
  pointless *unless the clock itself becomes nonlinear* (e.g. a nonlinear single-cell clock),
  in which case Test 3 shows the architecture is ready to exploit it.
- **The investigation is complete.** Tests 4–7 (sample staircase, noise, modality,
  heterogeneity) are no longer needed to explain the tie — it's explained. They would only
  quantify the *absolute-MAE* ceiling, which is a separate, lower-priority question.

---

## Test 4 — "feed training X back, does it reproduce Y?" as a coverage diagnostic (user's idea)  ✅ DONE

**Hypothesis.** Train the real model on (X→Y), feed the same X back, check it reproduces Y —
done as a diagnostic by hiding signal genes from the input. If reproduction degrades as we
hide signal genes, imperfect reproduction on real data is the missing-genes effect, not a bug.

**Prediction (before running).** 100% coverage → near-perfect reproduction; reproduction
degrades as coverage drops toward the predict-mean floor.

**Result (actual).** Reproduction stayed near-perfect at ALL coverages: 100%→0.93, 75%→1.11,
47%→1.50, 25%→1.28 (floor 12.82). **Prediction WRONG** — no degradation.

**Verdict — the surprise IS the finding.** On *training* data the model memorizes Y through
the **760 noise genes as per-row fingerprints**, so it reproduces training Y even when the
signal genes are hidden. Consequences:
- **The model reproduces its training data → the training path works, NOT broken.** (The
  user's "good sign" — confirmed.)
- **Reproduction-on-training is a WEAK test:** it passes via memorization even without seeing
  the signal, so it can't show the missing-genes effect. This is *why* we judge on held-out
  data, not training reproduction.
- **Explains a real number:** synthetic train-repro ~1 vs real train-Gill ~9 (from Test 1) is
  the **HFF-domination** effect — in the real pipeline the ~100 Gill cells are 0.24% of
  training, drowned by 42k HFF, so they're NOT memorized; here all 8k rows are one
  distribution and all get memorized. Consistent, not contradictory.
- To SEE the missing-genes effect cleanly, measure **held-out error vs coverage**
  (generalization), not training reproduction. (Optional Test 4b; the diagnosis already holds
  from Tests 0/1/3.)

Diagnosis unchanged and reinforced: nothing broken; the ΔAge tie is the linear-clock result.

## Test 4.1 — HELD-OUT error vs gene coverage (the clean version of Test 4)  ✅ DONE

**Hypothesis.** Measuring *generalization* (not training reproduction, which memorization
masks), held-out MAE should degrade as we hide signal genes — the direct signature of the
missing-genes effect.

**Prediction (before running).** 100% → low test MAE; degrades smoothly toward the
predict-mean floor as coverage drops; ~near floor by 25%.

**Result (actual).** **Prediction RIGHT.** model TEST MAE: 100%→2.99, 75%→9.16, 47%→12.06,
25%→12.95 (floor 12.80). Ridge tracks alongside (9.31/10.53/11.82/12.51). At 100% the model
crushes ridge (2.99 vs 9.31) — architecture confirmed again.

**Verdict — missing-genes effect DIRECTLY confirmed.** Held-out error rises toward the floor
as signal genes are hidden: you cannot predict unseen rows from genes you can't see. Both
model and ridge climb together — the ceiling is the *visible signal*, not the estimator.

**Honest nuance that completes the picture.** At 47% synthetic coverage the model is ~at the
floor (12.06), but the *real* model at ~47% gets ~9–14 — above ideal, below floor. Reason:
**synthetic genes are independent** (hiding one loses it entirely), whereas **real genes are
correlated**, so the visible 47% *partially proxies* the hidden 53%. Synthetic therefore
*overestimates* the degradation; real is milder. This is exactly why real ΔAge MAE (~9–14)
sits **between** the full-coverage ideal (~3) and the predict-mean floor (~13). All consistent.

---

## Test 5 — REOPENED: why is the model *worse* than ridge, not tied? (new question)  ⏳

**Why reopened (user caught this).** The "closed" conclusion said the model *ties* ridge
because ΔAge is linear. But the real LOOCV numbers are **model 14.29 ± 9.67 vs ridge
14.05 ± 8.03 — the model is WORSE by 0.24**, not equal. The linearity story explains why the
model doesn't *beat* ridge; it does NOT explain why the model is *worse*. A model that can
represent any linear function ridge can should never *lose* to it on a linear target. So the
0.24 deficit is a genuinely unexplained observation. Not closed.

### Test 5.0 — is the 0.24 deficit REAL, or noise? (do this first)

**Hypothesis.** With std ±8–9 across only 6 folds, a 0.24 mean difference may not be
statistically distinguishable from zero. Per-fold the model won O1 (5.39 vs 8.25) and the
aggregate is nearly equal — smells like noise + donor effects, not a systematic deficit.

**Prediction (before running).** The paired difference (model − ridge) across the 6 folds is
NOT significantly different from 0 (CI includes 0; wins/losses mixed). → honest statement
becomes "statistically tied," and there is no deficit to explain.

**Method.** Read per-fold model & ridge reg_MAE from all 6 `cellfate_loocv_*/reports/
holdout.json`, compute the paired per-fold difference, mean, a paired test, and win/loss
count. Uses data already on disk. Run once.

**Decision branches (before data):**
- **Difference indistinguishable from 0 (mixed wins)** → 0.24 is noise → "statistically tied"
  → the linearity story stands, investigation closes honestly. No Test 5.1 needed.
- **Difference is consistent & real (model loses on most folds systematically)** → a genuine
  deficit → proceed to Test 5.1 to find the cause.

**Result (actual).** _[TO FILL — user runs test5_ridge_gap.py]_

**Verdict.** _[TO FILL]_

### Test 5.0 RESULT (user ran it — real per-fold data)

| fold | model MAE | ridge MAE | model−ridge | who wins |
|---|---|---|---|---|
| N2 | 21.79 | 21.59 | +0.21 | ridge |
| N3 | 29.69 | 26.27 | +3.42 | ridge |
| O1 | 5.39 | 8.25 | −2.87 | model |
| O2 | 7.54 | 9.31 | −1.77 | model |
| Y1 | 7.28 | 6.44 | +0.83 | ridge |
| Y2 | 14.06 | 12.45 | +1.61 | ridge |

mean(model − ridge) = **+0.24**, std of diffs = 2.28, **95% CI = [−2.16, +2.64]** (includes 0),
paired t = +0.26. Model wins 2/6, ridge wins 4/6.

**Verdict — prediction RIGHT: the 0.24 is NOISE.** The 95% CI comfortably includes 0 and the
per-fold diffs swing widely (model wins O1 by 2.87, loses N3 by 3.42). The difference is not
statistically distinguishable from zero. **Honest statement: the model is STATISTICALLY TIED
with ridge on ΔAge MAE** — exactly what a linear target predicts. No systematic deficit
exists, so **Test 5.1 is not needed** (per the pre-committed branch). This turns the earlier
hand-wave ("ties ridge") into a demonstrated, reviewer-proof claim.

### Test 5.1 — [only if 5.0 says the deficit is REAL] is it the multi-task tradeoff?

**Hypothesis.** The model optimizes fate + ΔAge *jointly*; ridge optimizes ΔAge *alone*. The
shared representation may sacrifice a little ΔAge accuracy to serve fate classification — so a
ΔAge-only model should match ridge while the multi-task model is slightly worse.

**Prediction (before running).** _[state before running]_

**Method.** Train the model on ΔAge only (drop/zero the fate loss), compare its held-out
ΔAge MAE to ridge and to the multi-task model. Gap closes → it's the fate/age tradeoff (a
real, explainable *cost*, not a bug). Gap persists → optimization gap (→ Test 5.1.1: train
longer / lower LR / less dropout on the age head, chosen a priori, once).

**Result (actual).** _[TO FILL]_

---

## Test 6 — can ANY model beat ridge on the REAL ΔAge? (the empirical test, not assumed)  ⏳

**Why this test (user caught the real hole).** Tests 3/4.1 used *synthetic* data and a
*structural* argument (clock is linear → ΔAge linear → ridge optimal). But that argument has
a gap: the model predicts ΔAge from x_pert *without* being given x_ctrl, so it must **infer
the per-donor control offset** from x_pert — and inferring which donor / what baseline need
NOT be linear. So it is *possible* a more powerful model beats ridge on real ΔAge by better
offset inference. We never measured this on real data. The stakes:
- **If nothing beats ridge on real ΔAge** → ridge is genuinely at the ceiling; the neural net
  tying it is correct; its ΔAge contribution is honestly "matches the linear optimum." Closed
  empirically.
- **If something beats ridge** → real exploitable structure exists → our net only tying ridge
  is a REAL underperformance to diagnose (optimization / multi-task / architecture).

**Method.** On the REAL Gill x_pert → ΔAge pairs, under the SAME leave-one-donor-out protocol,
compare off-the-shelf, deterministic, no-GPU models:
1. **ridge** (linear incumbent)
2. **gradient-boosted trees** + **random forest** (nonlinearity + interactions)
3. **kernel ridge (RBF)** (smooth nonlinearity)
Report per-fold + aggregate MAE, and whether any nonlinear model beats ridge by more than
noise (reuse the paired logic from Test 5.0).

**Prediction (before running, on record).** I lean that trees/kernels will **NOT**
meaningfully beat ridge on real ΔAge — because the per-donor offset, though nonlinear in
principle, is hard to infer from only ~5 training donors (too few to learn donor structure).
But I genuinely don't know; that's why we run it. If a tree/kernel DOES beat ridge, I'm wrong
and there's a real lead.

**Decision branches (before data):**
- **No nonlinear model beats ridge (within noise)** → ridge is at the ceiling; ΔAge is
  linearly predictable from this input; neural-net tie is correct → **closed empirically.**
- **A nonlinear model beats ridge (beyond noise)** → real structure exists → the neural net
  underperforms → REOPEN with a concrete target (Test 6.1: why doesn't the net capture what a
  tree can? optimization vs multi-task vs architecture).

**Result (actual).** _[TO FILL — user runs test6_beat_ridge.py]_

**Verdict.** _[TO FILL]_

### Test 6 RESULT (user ran it on REAL data)

| fold | ridge | boosted_trees | random_forest | kernel_rbf |
|---|---|---|---|---|
| N2 | 21.59 | 24.83 | 24.91 | 20.06 |
| N3 | 26.27 | 22.21 | 24.30 | 33.12 |
| O1 | 8.25 | 10.58 | 13.71 | 21.24 |
| O2 | 9.31 | 7.64 | 9.25 | 19.80 |
| Y1 | 6.44 | 8.94 | 12.04 | 20.30 |
| Y2 | 12.45 | 19.60 | 22.66 | 33.53 |

Aggregate MAE: **ridge 14.05** (best), boosted_trees 15.63, random_forest 17.81,
kernel_rbf 24.68. Paired vs ridge: boosted trees +1.58 (tied, noise), random forest +3.76
(tied, noise), kernel +10.63 (significantly WORSE).

**Verdict — prediction RIGHT: NO nonlinear model beats ridge on real ΔAge.** Ridge is the
best of all four. Crucially, **performance DEGRADES with model flexibility** (linear best →
trees worse → forest worse → kernel worst). This is the textbook signature of a **genuinely
linear signal**: extra flexibility can't find nonlinear structure that isn't there, so it just
overfits noise and gets worse. If hidden nonlinear structure existed (e.g. in the per-donor
offset inference), at least one flexible model would have beaten ridge — the OPPOSITE happened.

**This EMPIRICALLY closes the gap the structural argument left open.** The worry that the
per-donor control offset (ΔAge = w·(x_pert − x_ctrl), x_ctrl unknown to the model) might be
nonlinearly exploitable is now ruled out: trees and kernels had every chance to exploit it and
found nothing. Real ΔAge IS linearly predictable from this input. The neural net tying ridge is
CORRECT — not an underperformance — now demonstrated on real data, not assumed.

---

## Test 7 — does the model's RANKING beat ranking-by-ridge-ΔAge? (does the model earn its keep on ranking?)  ⏳

**Why (user's question).** We proved ridge matches the model on ΔAge *magnitude* (Test 6). But
the model's headline is RANKING (Spearman 0.69) via its RES score (combined fate + ΔAge). The
honest question: does that RES ranking actually beat simply **sorting perturbations by ridge's
predicted ΔAge**? If not, the RES machinery isn't earning its keep for ranking either.

**The subtlety (state it so the test is fair).** "Ranking quality" = how well a predicted
ordering matches the TRUE ordering by actual ΔAge (which perturbation rejuvenates most). Ridge
predicts ΔAge directly, so "sort by ridge ΔAge" is a STRONG baseline — it's ordering by an
estimate of the exact quantity. That's the right bar: if the model's RES can't beat sorting by
a linear ΔAge prediction, RES adds nothing for ranking.

**Method.** Under the SAME leave-one-donor-out protocol, on the held-out donor's cells:
- **Model ranking:** rank by the model's RES score → Spearman vs true ΔAge.
- **Ridge ranking:** rank by ridge's predicted ΔAge → Spearman vs true ΔAge.
- (Also try ridge-ΔAge directly as the score, and the model's own ΔAge as a score, to
  separate "RES vs ΔAge-sort" from "model vs ridge".)
Report per-fold Spearman for each, aggregate mean ± std, and the paired difference
(model − ridge) across the 6 folds with a 95% CI (same paired logic as Test 5/6).

**Prediction (before running, on record).** My honest lean: **model_RES will be roughly TIED
with ranking-by-ridge-ΔAge**, because ΔAge is the dominant ranking signal and ridge predicts
it as well as the model (Test 6). I expect model_RES ≈ model_dAge ≈ ridge_dAge, all near the
~0.69 aggregate. If model_RES clearly beats ridge_dAge, the fate/uncertainty in RES adds
ranking value; if it clearly loses, RES is over-engineered for ranking. Prediction: tied (CI
includes 0).

**Decision branches (before data):**
- **Model RES ranking beats ridge-ΔAge ranking (beyond noise)** → the model earns its keep on
  ranking; RES's fate+ΔAge combination genuinely helps ordering. Real contribution confirmed.
- **Tied** → ranking is driven by ΔAge, which ridge predicts equally well → the model's ranking
  is NOT a unique contribution over a linear baseline (honest, important finding).
- **Model loses** → RES actively hurts ranking vs a simple ΔAge sort → reconsider RES for ranking.

**Result (actual).** _[TO FILL — user runs test7_ranking.py]_

**Verdict.** _[TO FILL]_

---

**Result (actual).** _[TO FILL — user runs test7_ranking.py]_

**Verdict.** _[TO FILL]_

### Test 7 RESULT (user ran it)

| fold | model_RES | model_dAge | ridge_dAge |
|---|---|---|---|
| N2 | +0.742 | +0.910 | +0.957 |
| N3 | +0.804 | +0.909 | +0.925 |
| O1 | +0.684 | +0.990 | +0.960 |
| O2 | +0.507 | +0.970 | +0.952 |
| Y1 | +0.706 | +0.960 | +0.951 |
| Y2 | +0.674 | +0.947 | +0.983 |

Aggregate Spearman: **model_RES 0.686 ± 0.091**, **model_dAge 0.948 ± 0.030**,
**ridge_dAge 0.955 ± 0.017**. Paired (model_RES − ridge_dAge) = −0.268, 95% CI
[−0.381, −0.155]. model_RES loses on 6/6 folds.

**Verdict — prediction WRONG (I said tied; it's clearly worse). BIG finding.** Two facts:
1. **ridge_dAge (0.955) beats model_RES (0.686) on ALL 6 folds** — a simple linear ΔAge sort
   ranks far better than the RES score.
2. **model_dAge (0.948) ≈ ridge_dAge** — ranking by the model's OWN ΔAge is also ~0.95. So the
   model's ΔAge predictions rank beautifully; **it is RES that destroys the ranking**, dragging
   0.95 down to 0.69.

So the "0.69 headline" is the RES ranking — the WORST of the three. The model's predictions are
fine; the RES *score* is the problem. RES = φ·S^k·g·exp(−λ·P_loss) multiplies ΔAge by fate
terms, and that scrambles the clean ΔAge ordering.

**CRUCIAL caveat before concluding RES is "bad".** Test 7 scores every ranking against **true
ΔAge only** (quality = −ΔAge). But RES is deliberately NOT trying to rank by ΔAge alone — it
ranks by **SAFE rejuvenation** (rejuvenating AND safe), down-weighting cells that rejuvenate
but are unsafe. So Test 7 may be penalizing RES for doing its job — measuring it against a
target it was explicitly designed NOT to optimize. Two competing explanations, must
distinguish: (a) RES genuinely adds noise and is bad for ranking; (b) RES optimizes a
different, better objective (safe rejuvenation) and Test 7 used the wrong ground truth. → Test 7.1.

### Test 7.1 — score the rankings against SAFE REJUVENATION (RES's actual objective)  ⏳

**Hypothesis.** RES ranks by *safe rejuvenation*, not raw ΔAge. If we score the three rankings
against a ground truth that rewards cells that are BOTH rejuvenating AND safe (not just
most-negative ΔAge), RES may win — because it deliberately penalizes unsafe-but-rejuvenating
cells that a pure ΔAge sort ranks too high. If RES still loses even against safe-rejuvenation,
RES genuinely doesn't work and we should rank by ΔAge.

**Prediction (before running).** Genuinely uncertain, but my honest lean: **RES will STILL lose (or at best tie) even against safe-rejuvenation.** Two reasons: (1) the held-out Gill cells may be mostly one fate class (little "unsafe" to exploit → degenerates to Test 7), and (2) the model's OUT-OF-DONOR fate calibration was poor in LOOCV (ECE 0.26), so RES's fate weighting likely injects noise rather than aligning with true safety. If RES loses here too → "RES doesn't earn its keep for ranking; rank by ΔAge." If RES wins → it does its job and Test 7 used the wrong target. (I've been wrong repeatedly here, so low confidence.)

**Method.** Under leave-one-donor-out, define a "safe-rejuvenation quality" ground truth on the
held-out cells — e.g. reward = rejuvenation magnitude but ZEROED/penalized for cells whose true
fate is loss/death (unsafe). Rank the same three scores (model_RES, model_dAge, ridge_dAge)
against this target; report per-fold + aggregate Spearman + paired (model_RES − best_ΔAge_sort).
Try a couple of reasonable safe-rejuvenation definitions to avoid gaming one.

**Decision branches (before data):**
- **RES beats the ΔAge sorts against safe-rejuvenation** → RES does its job; Test 7 was unfair
  (wrong target); report the *right* metric. Model earns its keep on the objective it targets.
- **RES still loses even against safe-rejuvenation** → RES genuinely doesn't help ranking →
  honest finding: rank by ΔAge, reconsider/redesign or drop RES.

### Test 7.1 RESULT (user ran it)

Fate composition of held-out age-valid cells (unsafe = loss/death): N2 0/21, N3 3/21,
O1 4/21, O2 5/21, Y1 8/19, Y2 5/21 — so most folds DO have unsafe cells (safety can matter).

| target | model_RES | model_dAge | ridge_dAge | paired (RES−ridge) 95% CI |
|---|---|---|---|---|
| GATED | −0.005 | 0.292 | 0.295 | −0.300 [−0.567,−0.034] RES WORSE |
| PENALIZED | 0.137 | 0.414 | 0.414 | −0.277 [−0.525,−0.030] RES WORSE |

precision@5 (top-5 truly safe & rejuvenating): RES 0.20, model_dAge 0.27, ridge_dAge 0.30.

**Verdict — prediction RIGHT: RES loses even against safe-rejuvenation, its OWN objective.**
The Test 7 caveat is resolved: RES isn't losing due to a wrong target; it genuinely doesn't
help ranking. **But the per-fold pattern PINPOINTS the cause:** RES quality tracks INVERSELY
with the number of unsafe cells —
  N2(0 unsafe): +0.742 · N3(3): +0.676 · O1(4): −0.257 · O2(5): −0.431 · Y1(8): −0.534.
The MORE unsafe cells a fold has, the WORSE RES does. If the model's fate predictions were
right, RES should get BETTER with more unsafe cells to correctly penalize. The opposite means
**the model's fate predictions are WRONG on held-out donors — RES down-weights the wrong
cells**, scrambling the ranking exactly when safety matters. Consistent with the LOOCV
out-of-donor calibration failure (ECE 0.26).

**Conclusions (two, both evidenced):**
1. **RES does NOT earn its keep for ranking** — worse than a plain ΔAge sort against BOTH raw
   ΔAge (Test 7) and safe-rejuvenation (Test 7.1), on ~every fold, by Spearman AND precision@5.
   The RES machinery SUBTRACTS value. Honest action: for ranking, rank by ΔAge (model or ridge,
   ~equal); reconsider/redesign RES.
2. **Root cause = the fate predictions, not the RES formula per se.** RES is only as good as the
   fate probabilities feeding it, and those don't generalize across donors. This PREVIEWS the
   next question (does fate classification work?) — early signal: not well enough out-of-donor
   to be useful in RES. → run the fate-classification test next.

**Reframes the headline honestly:** the defensible ranking number is ranking-by-ΔAge (~0.95 vs
true ΔAge, ~0.30–0.41 vs safe-rejuvenation), NOT the RES 0.69. The model's *value* narrows to:
calibrated-in-distribution fate + ΔAge that matches linear + honest uncertainty — with RES and
out-of-donor fate flagged as not-yet-working.

---

## Test 7.2 — RES formula vs plain ΔAge, SAME ΔAge input (clean isolation, user's design)  ⏳

**Why (user).** Test 7 had a confound: model_RES used the MODEL's ΔAge while ridge_dAge used
RIDGE's ΔAge, so two things varied (the RES formula AND ΔAge quality). Fix: hold ΔAge constant.
Feed the SAME ΔAge (ridge's — the best predictor, Test 6, ranks ~0.95) into both:
  A = rank by that ΔAge directly.
  B = rank by RES(that same ΔAge, model fate S/P_loss, uncertainty, OOD).
Now the ONLY difference is the RES formula. If B loses to A, it is unambiguously the RES
formula hurting — not a ΔAge-quality difference.

**Prediction (before running, on record).** B (RES) loses to A — the RES formula, fed even a
GOOD ΔAge, still hurts ranking, because the model's out-of-donor fate predictions inject noise
(matches Test 7.1's pattern: RES worse where more unsafe cells). If B ties/wins, that FLIPS the
earlier verdict: RES is fine and Test 7's loss was just the model's worse ΔAge feeding it.

**Method.** Per held-out donor, build A = ridge ΔAge, B = RES(ridge ΔAge + model fate +
uncertainty + OOD). Score both against (i) true ΔAge and (ii) safe-rejuvenation (gated). Report
per-fold + aggregate Spearman + paired (B − A) 95% CI. Same cells, same ΔAge, one variable.

**Decision branches (before data):**
- **B < A (RES hurts even with good ΔAge)** → the RES FORMULA degrades ranking → for ranking,
  use ΔAge directly; RES needs redesign (or is for a different purpose than ranking).
- **B ≈ A** → RES neither helps nor hurts ranking → it is redundant for ranking; ΔAge suffices.
- **B > A (RES helps when fed good ΔAge)** → RES formula is sound; Test 7's loss was the model's
  worse ΔAge, not RES → pair RES with the better ΔAge and re-evaluate. (Would flip the verdict.)

### Test 7.2 RESULT (user ran it)

Same ΔAge (ridge) fed to both; only the RES transform differs.

| target | A = ridge ΔAge | B = RES(ridge ΔAge) | paired (B−A) 95% CI | verdict |
|---|---|---|---|---|
| vs TRUE ΔAge | 0.955 | 0.654 | −0.300 [−0.473,−0.128] | **RES WORSE (6/6 folds)** |
| vs SAFE-REJUV | 0.295 | 0.116 | −0.179 [−0.420,+0.061] | tied (noise); RES wins N3/Y1 only |

**Verdict — prediction RIGHT: the RES FORMULA ITSELF degrades ranking.** With ΔAge held
identical, wrapping it in RES drops Spearman-vs-true-ΔAge from 0.955 to 0.654 on ALL 6 folds
(CI excludes 0). Clean isolation Test 7 could not do: UNAMBIGUOUSLY the RES transform, not the
ΔAge feeding it. Against safe-rejuvenation RES is a wash (CI includes 0, wins N3/Y1), so RES
adds no value on either target and actively hurts pure-ΔAge ranking.

**RES ranking question CLOSED (Tests 7 + 7.1 + 7.2 converge):** RES multiplies ΔAge by fate
terms (S, P_loss); those fate predictions are unreliable out-of-donor (ECE 0.26), so the
multiplication SCRAMBLES a clean ΔAge ordering. RES takes good information and makes it worse.
**Honest action: for ranking, use ΔAge directly (Spearman ~0.95 vs true ΔAge); do NOT use RES
for ranking.** Salvage depends on whether fate is recalibratable (Test 8.2).

---

## TESTS 8-12 — does the model earn its keep, do embeddings help, and where do failures come from?

**Numbering rule:** a genuinely different question gets the next whole integer; a zoom-in on an existing question gets a decimal under it. Map to subcommands: Test 8 = `fate_baseline` (8.1 = `indist_vs_donor`, 8.2 = `fate_cal_disc`); Test 9 = `string_dage` (9.1 = `string_fate`); Test 10 = full-transcriptome [pending]; Test 11 = `input_ablation`; Test 12 = `per_donor`.

  ⏳

**Context / honest guardrails baked in:**
- **Data-scale law (critical):** more-flexible models need MORE data to beat linear. Beating
  ridge on huge data (GenBio 2026, 4 lines × 2000 perts) does NOT imply beating it on 6 donors
  — the opposite is typical (flexibility overfits noise at small n; Test 6 showed exactly this).
  So NONE of these tests are assumed wins; each is a real question.
- **STRING unreachable from build env** → using **gene2vec** (200-dim, co-expression-based,
  external/frozen, covers 94% of panel) as a proxy for the prior-knowledge-embedding principle.
  Not identical to STRING (protein-interaction); if a signal appears, confirm with STRING later.
- **Epistemic-extrapolation warning:** the clock was trained on NATURAL aging; reprogramming is
  EXTREME rejuvenation outside its training domain — read all clock-based results with that caveat.

### Test 8 — does the model's FATE classification beat a simple baseline (logistic regression)?
**Hypothesis.** Fate (safe/loss/death) may be where the NN earns its keep (classification,
possibly pathway/interaction-dependent) — or it may only tie logistic regression, esp.
out-of-donor (LOOCV ECE was 0.26). **Prediction:** model ~ties logistic regression on held-out
donors (maybe wins in-distribution); no decisive out-of-donor win. **Branches:** model beats
LogReg beyond noise → fate is a real contribution; tie → fate head not earning its keep vs linear.

### Test 10 — can RIDGE predict the FULL transcriptome response of reprogramming? (feasibility)  [PENDING target def]
**Hypothesis.** The GenBio-style task (predict all genes' change) is where nonlinearity CAN help
— but first, can linear ridge even do it, as the baseline? **Prediction:** ridge captures the
bulk linear response (decent global R²) but misses the sparse/nonlinear DEGs. **Branches:** sets
the baseline that 9.3 (embeddings) must beat. NOTE: target definition needs care (predict
perturbed-vs-control expression change); flag if data structure doesn't support clean pairing.

### Test 9 — do gene embeddings (gene2vec) help on ΔAge? (ridge-raw vs ridge-emb vs MLP-emb)
**Hypothesis.** ΔAge is LINEAR (proven). Embeddings encode gene INTERACTIONS, which a linear
target ignores — so embeddings likely add nothing; a linear projection X@E is strictly lossy vs
raw X, and MLP-on-emb risks overfitting at n≈100. **Prediction (strong):** NO improvement on
ΔAge; ridge-raw ≥ both embedding variants. **Branches:** if embeddings DON'T help → confirms
ΔAge signal is purely linear, interaction structure irrelevant (clean finding). If they DO →
surprising, investigate.

### Test 9.1 — do gene embeddings help FATE? (logreg-raw vs logreg-emb vs MLP-emb)  [zoom-in on Test 9]
**Hypothesis.** Fate is classification, plausibly interaction-dependent — the place embeddings
COULD help if anywhere. **Prediction (uncertain, low confidence):** slight help at best;
out-of-donor generalization likely still limited by 6 donors. **Branches:** embeddings improve
held-out fate beyond noise → real, field-consistent lever → pursue STRING specifically; no help
→ fate limited by data/donors, not representation.

### Test 10.1 — do gene embeddings help FULL transcriptome response? (zoom-in on Test 10)  [PENDING]
**Hypothesis.** This is the task where the 2026 literature shows embeddings win — AT SCALE.
**Prediction:** on 6 donors, embeddings give little/no gain over ridge (data-scale law); a gain
here would be the strongest pro-embedding signal. **Branches:** as 9.1/9.2. Depends on 8.2 target.

**Results (actual).** _[TO FILL — user runs test_suite.py <subcommand>]_

---

## FAILURE-ISOLATION TESTS — isolate ONE factor at a time, across all 6 folds  ⏳

**Principle:** the earlier tests MEASURE gaps; these PINPOINT causes. Each holds everything
fixed and varies exactly ONE thing, so a moving number names its own cause. Every test reports
per-fold + aggregate across all 6 LODO folds and states explicitly WHAT HELPED / WHAT DIDN'T.

### Test 8.1 — IN-DISTRIBUTION vs OUT-OF-DONOR (zoom-in: fate/ΔAge failure = fitting or generalization?)
**Vary only: what is held out.** Compare model & baselines on (a) held-out CELLS from training
donors (in-distribution) vs (b) held-out DONOR (LODO), for BOTH fate and ΔAge.
**Question:** is the failure in *fitting* (bad even in-distribution) or *generalization* (fine
in-dist, breaks across donors)? **Prediction:** ΔAge ~same both ways (linear, no gen. gap);
fate GOOD in-distribution but COLLAPSES out-of-donor (the ECE 0.26 story). **What it pinpoints:**
if fate is in-dist-good/out-of-donor-bad → the problem is DONOR SHIFT, not the fate head design
→ fix = more donors / donor adaptation, not architecture.

### Test 11 — INPUT ABLATION (does the model use cell×perturbation together?)  [distinct question]
**Vary only: what the model is given.** For fate & ΔAge, compare a baseline trained on
{cell-state only} vs {perturbation only} vs {both}. **Question:** does using BOTH beat using
either alone? (If not, the "joint" modeling adds nothing.) **Prediction:** perturbation is ~fixed
(OSKM) so u_only is weak; x_only ~ x+u → the model rides cell state, perturbation adds little
here. **What it pinpoints:** whether the state×perturbation interaction — the model's core thesis
— is actually load-bearing on THIS data, per fold.

### Test 8.2 — CALIBRATION vs DISCRIMINATION for fate (zoom-in: which kind of fate failure?)
**Vary only: the metric axis.** For fate, separately measure DISCRIMINATION (PR-AUC / ROC — does
it rank safe vs unsafe correctly?) and CALIBRATION (ECE / reliability — are the probabilities
right?), in-dist and out-of-donor. **Question:** does the fate head rank correctly but report
miscalibrated probabilities (fixable by recalibration), or is the ranking itself bad (needs a
better classifier)? **Prediction:** decent discrimination, poor calibration out-of-donor →
recalibration (temperature/Platt per donor) is the cheap fix; RES could then use recalibrated
fate. **What it pinpoints:** whether RES's failure (Test 7.1) is fixable by fate recalibration
without touching the model.

### Test 12 — PER-DONOR JACKKNIFE (are aggregates robust or hostage to specific donors?)  [distinct question]
**Vary only: which donor is excluded from the aggregate.** Recompute the headline metrics
dropping one donor at a time (jackknife over the 6). **Question:** are the poor aggregates driven
by 1-2 outlier donors (e.g. N3, +30 ΔAge, 0 unsafe cells) or broadly true? **Prediction:** N3 (and
maybe N2) dominate the bad numbers; dropping them changes the story. **What it pinpoints:** whether
conclusions are robust or hostage to specific donors — essential before any claim in the writeup.

**Cross-cutting reporting rule:** every test prints a "WHAT HELPED / WHAT DIDN'T" line per factor,
and a per-fold table so no aggregate hides a fold-level effect. Guardrail unchanged: each run once,
nothing tuned against results.

**Results (actual).** _[TO FILL — user runs test_suite.py <subcommand>]_

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
