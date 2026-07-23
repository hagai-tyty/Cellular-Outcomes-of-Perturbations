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

## Test 1 — Is the model OVERFITTING the training donors? (tests H2/overfitting, across all 6 folds)  ✅ DONE

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

## Test 2 — [conditional on Test 1] — N/A, NOT NEEDED (Test 1 ruled out overfitting → regularization path dead)

If Test 1 says **overfitting**: Test 2a = re-train ΔAge with stronger regularization / a
simpler head, ONCE, and check whether the model now beats ridge on held-out. (Guardrail: one
config chosen by principle beforehand, not swept against the test number.)

If Test 1 says **not overfitting**: Test 2 = Test 3 below (linearity) becomes the next step.

_[Specify precisely once Test 1 result is in — per the rules, we decide the exact test after
seeing which branch we're on.]_

---

## Test 3 — Is the true ΔAge signal LINEAR (nothing to fix) or is the model failing? (tests H4 vs H1)  ✅ DONE — see RESULT below

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

**Result (actual).** → recorded below in "Test 3 RESULT".

**Verdict.** → see "Test 3 RESULT" below.

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

> ⚠️ **Scope note (2026-07-13).** The substance below is correct and has held up — ΔAge is linear
> (reconfirmed by Tests 5.0, 6, and 9). But it "closes the case" only for ΔAge **magnitude**; the
> investigation went on to ranking/RES, fate, and input. For current status read **"Where the
> investigation stands"** + the decision tree.

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

> ⚠️ **SUPERSEDED (pre-Test-5 plan) — history only.** The actual Tests 4–12 pursued different
> questions (coverage diagnostic → ridge-deficit → beat-ridge → ranking/RES → fate discrimination &
> calibration → input ablation → per-donor jackknife → embeddings). The "Tests 4–7" sketched below
> (staircase / noise / modality / heterogeneity) were **not** the path taken. See the decision tree.

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

> ⚠️ **Partial log — current only through Test 6.** Does not include Tests 7–12. For the complete
> current picture see **"Where the investigation stands"** + the decision tree.

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

> ⚠️ **SUPERSEDED — do not read as current.** Written after Test 3, before the ranking/fate/input
> tests. It calls the investigation "complete" (it wasn't) and cites the RES ranking (Spearman 0.69)
> as the model's contribution — both later **overturned**: Tests 7/7.1/7.2 showed RES *hurts* ranking
> (rank by ΔAge, ~0.95), and Test 8.2 showed fate is *mis*calibrated. Current status: **"Where the
> investigation stands"** + the decision tree.

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

## Test 5 — REOPENED: why is the model *worse* than ridge, not tied? (new question)  ✅ RESOLVED (5.0 done; 5.1 not needed)

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

**Result (actual).** → recorded below in "Test 5.0 RESULT".

**Verdict.** → see "Test 5.0 RESULT" below.

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

**Result (actual).** → NOT RUN. Test 5.0 showed the deficit is noise, so 5.1 was unnecessary (pre-committed branch).

---

## Test 6 — can ANY model beat ridge on the REAL ΔAge? (the empirical test, not assumed)  ✅ DONE — see RESULT below

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

**Result (actual).** → recorded below in "Test 6 RESULT".

**Verdict.** → see "Test 6 RESULT" below.

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

## Test 7 — does the model's RANKING beat ranking-by-ridge-ΔAge? (does the model earn its keep on ranking?)  ✅ DONE — see RESULT below

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

**Result (actual).** → recorded below in "Test 7 RESULT".

**Verdict.** → see "Test 7 RESULT" below.

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

### Test 7.1 — score the rankings against SAFE REJUVENATION (RES's actual objective)  ✅ DONE — see RESULT below

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

## Test 7.2 — RES formula vs plain ΔAge, SAME ΔAge input (clean isolation, user's design)  ✅ DONE — see RESULT below

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

  **Status: MIXED** — Tests 8, 8.1, 8.2, 11, 12 are DONE (see their RESULT sections below); Tests 9, 9.1, 10, 10.1 remain PENDING.

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

**Failure-isolation principle (applies to Tests 8.1, 8.2, 11, 12).** The earlier tests MEASURE
gaps; these PINPOINT causes. Each holds everything fixed and varies exactly ONE thing, so a moving
number names its own cause. Every test reports per-fold + aggregate across all 6 LODO folds, prints
a "WHAT HELPED / WHAT DIDN'T" line per factor and a per-fold table so no aggregate hides a
fold-level effect. Guardrail: each run once, nothing tuned against results.

### Test 8 — does the model's FATE classification beat a simple baseline (logistic regression)?  ✅ DONE
**Hypothesis.** Fate (safe/loss/death) may be where the NN earns its keep (classification,
possibly pathway/interaction-dependent) — or it may only tie logistic regression, esp.
out-of-donor (LOOCV ECE was 0.26). **Prediction:** model ~ties logistic regression on held-out
donors (maybe wins in-distribution); no decisive out-of-donor win. **Branches:** model beats
LogReg beyond noise → fate is a real contribution; tie → fate head not earning its keep vs linear.

#### Test 8 — RESULT (fate: model vs logistic regression) — read past the "tied" verdict.

| fold | model PR-AUC | logreg PR-AUC | margin |
|---|---|---|---|
| N3 | 0.997 | 0.997 | tie |
| O1 | 1.000 | 0.968 | +0.032 |
| O2 | 1.000 | 0.984 | +0.016 |
| Y1 | 0.961 | 0.636 | **+0.325 (model wins big)** |
| Y2 | 1.000 | 0.964 | +0.036 |

Aggregate model 0.992 vs logreg 0.910 (+0.082). Paired CI [−0.087,+0.251] includes 0 → script
says "tied". **But reading the folds: the model WINS or TIES on all 5, never loses.** The
non-significance is the WIDE CI on n=5 with 4 folds SATURATED (~1.0, task easy there), not a
genuine tie.

**The Y1 fold is the tell:** on the HARD fold, logreg COLLAPSES to 0.636 while the neural net
holds at 0.961. When fate is easy both saturate; when it gets hard, the LINEAR model breaks and
the net does NOT — the signature of the net capturing nonlinear structure logreg can't. This is
the OPPOSITE of the ΔAge/ranking story (where simple baselines matched/beat the model everywhere).

**Honest call:** NOT a statistically significant win (5 folds, saturation, ~21 pts/fold), so cannot
claim "model significantly beats logreg on fate." BUT the direction is genuinely favorable — wins/
ties every fold, never loses, clearly beats logreg on the one hard fold. **Fate discrimination is
the model's STRONGEST claim to earning its keep, and a linear baseline does NOT cleanly match it**
(unlike ΔAge). Caveat: rides heavily on Y1; more donors would confirm or wash it out.

### Test 8.1 — IN-DISTRIBUTION vs OUT-OF-DONOR (zoom-in: fate/ΔAge failure = fitting or generalization?)  ✅ DONE
**Vary only: what is held out.** Compare model & baselines on (a) held-out CELLS from training
donors (in-distribution) vs (b) held-out DONOR (LODO), for BOTH fate and ΔAge.
**Question:** is the failure in *fitting* (bad even in-distribution) or *generalization* (fine
in-dist, breaks across donors)? **Prediction:** ΔAge ~same both ways (linear, no gen. gap);
fate GOOD in-distribution but COLLAPSES out-of-donor (the ECE 0.26 story). **What it pinpoints:**
if fate is in-dist-good/out-of-donor-bad → the problem is DONOR SHIFT, not the fate head design
→ fix = more donors / donor adaptation, not architecture.

#### Test 8.1 — RESULT (in-dist vs out-of-donor) — read from raw numbers. CLARIFYING.

| fold | fate PR in-dist | ΔAge MAE in-dist | fate PR out-donor | ΔAge MAE out-donor |
|---|---|---|---|---|
| N2 | 0.940 | 3.21 | n/a | 21.79 |
| N3 | 0.931 | 4.36 | 0.997 | 29.69 |
| O1 | 0.929 | 4.53 | 1.000 | 5.39 |
| O2 | 0.933 | 4.44 | 1.000 | 7.54 |
| Y1 | 0.934 | 4.00 | 0.961 | 7.28 |
| Y2 | 0.935 | 4.46 | 1.000 | 14.06 |

**CORRECTION to earlier claims (user forced this by checking data scale):**
1. **Fate DISCRIMINATION is genuinely GOOD** — in-dist PR-AUC 0.929–0.940 on THOUSANDS of held-out
   HFF cells (not 21). Steady ~0.93. And it does NOT collapse out-of-donor (0.96–1.00). So the
   earlier "fate fails out-of-donor" narrative was WRONG. Fate discriminates safe-vs-unsafe well,
   confirmed at real sample volume. (The out-of-donor 1.0 is mildly inflated by ~21 pts, but the
   in-dist 0.93 is trustworthy.)
2. **ΔAge FITTING is genuinely GOOD** — in-dist MAE ~4 years, tight (3.2–4.5) across all folds. The
   model predicts ΔAge well when test cells resemble training. The "~14 MAE" we kept quoting is an
   ARTIFACT of averaging in donor-generalization failure on 2 atypical donors.

**The clean fitting-vs-generalization split (what this test was for):**
- Fitting: FINE (ΔAge ~4 in-dist; fate ~0.93 in-dist).
- Generalization: the real limit — out-of-donor ΔAge is ~5–7 for TYPICAL donors (O1/O2/Y1) but
  22–30 for ATYPICAL ones (N2 zero-unsafe, N3 +30 outlier). So the model is GOOD on typical
  held-out donors and BREAKS on the 2 weird ones. Donor-heterogeneity (H6), quantified.

**Reframe of the whole picture (more accurate + more positive):** the model WORKS in-distribution
(fate discrimination ~0.93, ΔAge ~4); its weakness is generalizing to ATYPICAL donors from a
6-donor cohort. That is a data-scale/heterogeneity story, not "the model can't predict." Unchanged:
ridge matches model on ΔAge (linear target); RES hurts ranking (but now known NOT to be a fate-
discrimination problem → must be fate CALIBRATION → Test 8.2).

### Test 8.2 — CALIBRATION vs DISCRIMINATION for fate (zoom-in: which kind of fate failure?)  ✅ DONE
**Vary only: the metric axis.** For fate, separately measure DISCRIMINATION (PR-AUC / ROC — does
it rank safe vs unsafe correctly?) and CALIBRATION (ECE / reliability — are the probabilities
right?), in-dist and out-of-donor. **Question:** does the fate head rank correctly but report
miscalibrated probabilities (fixable by recalibration), or is the ranking itself bad (needs a
better classifier)? **Prediction:** decent discrimination, poor calibration out-of-donor →
recalibration (temperature/Platt per donor) is the cheap fix; RES could then use recalibrated
fate. **What it pinpoints:** whether RES's failure (Test 7.1) is fixable by fate recalibration
without touching the model.

#### Test 8.2 — RESULT (fate discrimination vs calibration) — diagnosis CONFIRMED, fixable kind.

| fold | ROC-AUC (rank) | ECE raw (calib) | ECE recal | reduction |
|---|---|---|---|---|
| N3 | 0.981 | 0.275 | 0.145 | 47% |
| O1 | 1.000 | 0.316 | 0.147 | 53% |
| O2 | 1.000 | 0.271 | 0.099 | 63% |
| Y1 | 0.932 | 0.271 | 0.243 | 10% (stubborn) |
| Y2 | 1.000 | 0.270 | 0.132 | 51% |

**Confirmed (as Test 8.1 predicted): fate RANKS well but is MISCALIBRATED.** ROC-AUC 0.93–1.00
(near-perfect discrimination) but ECE 0.27–0.32 (probabilities ~30% off). The two are decoupled →
the classic "good discrimination, bad calibration" pattern = the CHEAP-to-fix kind (recalibration),
NOT the "need a better classifier" kind.

**Platt recalibration roughly HALVES ECE on 4/5 folds** (0.28→0.13 typical; Y1 barely moves).
Doesn't reach the 0.05 bar, but a trivial post-hoc fix that touches nothing in the model cuts
calibration error ~half. Caveat: ~21 out-of-donor pts/fold → noisy; Y1 stubborn = not magic.

**UNLOCKS a real path to salvage RES (first in a while):** RES hurt ranking (7.2) by multiplying
ΔAge by fate probabilities; those probs are miscalibrated (8.2); recalibration halves the error.
So **recalibrate fate BEFORE feeding RES** might stop RES scrambling the ranking. PLAUSIBLE, not
proven (recal only halved error; must re-test 7.2 with recalibrated fate). Concrete next step:
re-run the RES-isolation with Platt-recalibrated fate probabilities.

### Test 9 — do gene embeddings (gene2vec) help on ΔAge? (ridge-raw vs ridge-emb vs MLP-emb)  ✅ DONE
**Hypothesis.** ΔAge is LINEAR (proven). Embeddings encode gene INTERACTIONS, which a linear
target ignores — so embeddings likely add nothing; a linear projection X@E is strictly lossy vs
raw X, and MLP-on-emb risks overfitting at n≈100. **Prediction (strong):** NO improvement on
ΔAge; ridge-raw ≥ both embedding variants. **Branches:** if embeddings DON'T help → confirms
ΔAge signal is purely linear, interaction structure irrelevant (clean finding). If they DO →
surprising, investigate.

#### Test 9 — RESULT (gene2vec embeddings on ΔAge) — prediction RIGHT: embeddings don't help.

ΔAge MAE (lower better), per fold:

| fold | ridge_raw | ridge_emb | mlp_emb |
|---|---|---|---|
| N2 | 21.80 | 21.48 | 23.90 |
| N3 | 26.59 | 25.20 | 23.81 |
| O1 | 8.03 | 12.25 | 7.85 |
| O2 | 8.90 | 8.67 | 9.49 |
| Y1 | 6.54 | 12.42 | 5.97 |
| Y2 | 12.85 | 17.91 | 14.24 |

Aggregate MAE: ridge_raw 14.12 · ridge_emb 16.32 · mlp_emb 14.21. Paired vs ridge_raw:
ridge_emb +2.20 [−1.15, +5.55] (tied/worse); mlp_emb +0.09 [−1.71, +1.89] (tied).

**Verdict — prediction RIGHT: embeddings add NOTHING to ΔAge.** The linear projection X@E is
strictly lossy (ridge_emb worse on O1/Y1/Y2); MLP-on-emb only matches raw ridge. This CONFIRMS —
a third independent way after Tests 3 and 6 — that the ΔAge signal is **purely linear**; the
interaction structure embeddings encode is irrelevant to a linear target. No lever here.

### Test 9.1 — do gene embeddings help FATE? (logreg-raw vs logreg-emb vs MLP-emb)  [zoom-in on Test 9]  ✅ DONE
**Hypothesis.** Fate is classification, plausibly interaction-dependent — the place embeddings
COULD help if anywhere. **Prediction (uncertain, low confidence):** slight help at best;
out-of-donor generalization likely still limited by 6 donors. **Branches:** embeddings improve
held-out fate beyond noise → real, field-consistent lever → pursue STRING specifically; no help
→ fate limited by data/donors, not representation.

#### Test 9.1 — RESULT (gene2vec embeddings on FATE) — at most a marginal, non-significant bump.

safe-class PR-AUC (higher better), per fold:

| fold | logreg_raw | logreg_emb | mlp_emb |
|---|---|---|---|
| N3 | 0.997 | 0.943 | 1.000 |
| O1 | 0.971 | 0.880 | 0.997 |
| O2 | 0.984 | 0.903 | 1.000 |
| Y1 | 0.639 | 0.703 | 0.707 |
| Y2 | 0.964 | 0.920 | 1.000 |

Aggregate PR-AUC: logreg_raw 0.911 · logreg_emb 0.870 · mlp_emb 0.941. Paired vs logreg_raw:
logreg_emb −0.041 [−0.118, +0.035] (tied); mlp_emb +0.030 [−0.001, +0.060] (tied, at the edge).

**Verdict — embeddings do NOT clearly help fate.** The linear embedding (logreg_emb) hurts; the
nonlinear MLP-on-emb is nominally best (0.941) and nudges the hard Y1 fold up (0.707 vs 0.639),
but its CI lower bound is −0.001 — a whisker short of significance on n=5. Honest read: **at most
a marginal, unconfirmed bump; representation is NOT the fate lever — data/donors are** (per the
pre-committed branch). Not worth chasing STRING on this data; revisit only with more donors.

### Test 10 — can RIDGE predict the FULL transcriptome response of reprogramming? (feasibility)  [PENDING target def]
**Hypothesis.** The GenBio-style task (predict all genes' change) is where nonlinearity CAN help
— but first, can linear ridge even do it, as the baseline? **Prediction:** ridge captures the
bulk linear response (decent global R²) but misses the sparse/nonlinear DEGs. **Branches:** sets
the baseline that 9.3 (embeddings) must beat. NOTE: target definition needs care (predict
perturbed-vs-control expression change); flag if data structure doesn't support clean pairing.

### Test 10.1 — do gene embeddings help FULL transcriptome response? (zoom-in on Test 10)  [PENDING]
**Hypothesis.** This is the task where the 2026 literature shows embeddings win — AT SCALE.
**Prediction:** on 6 donors, embeddings give little/no gain over ridge (data-scale law); a gain
here would be the strongest pro-embedding signal. **Branches:** as 9.1/9.2. Depends on 8.2 target.

**Results for Tests 10 / 10.1:** _[TO FILL — pending]_

### Test 11 — INPUT ABLATION (does the model use cell×perturbation together?)  [distinct question]  ✅ DONE
**Vary only: what the model is given.** For fate & ΔAge, compare a baseline trained on
{cell-state only} vs {perturbation only} vs {both}. **Question:** does using BOTH beat using
either alone? (If not, the "joint" modeling adds nothing.) **Prediction:** perturbation is ~fixed
(OSKM) so u_only is weak; x_only ~ x+u → the model rides cell state, perturbation adds little
here. **What it pinpoints:** whether the state×perturbation interaction — the model's core thesis
— is actually load-bearing on THIS data, per fold.
> ⚠️ **Premise later FALSIFIED (2026-07-13).** The prediction above assumed "perturbation is ~fixed."
> It is not: `dose_time = [log10(dose), log(time_h)]` varies richly (≈12 timepoints d0/Fib–d54;
> fitted scaler std ≈ [1.29, 3.26], both nonzero). The pre-registered guess is kept for the record;
> see the **CORRECTION** in the RESULT below and **Test 11.1**.

#### Test 11 — RESULT (input ablation) — x+u ≈ x_only, u_only fails for ΔAge; interaction not load-bearing *as measured*.  [⚠ INTERPRETATION CORRECTED 2026-07-13 — see CORRECTION below]

ΔAge MAE — x_only / u_only / x+u (per fold):
N2 21.80/16.97/21.59 · N3 26.59/46.76/26.27 · O1 8.03/36.12/8.25 · O2 8.90/32.50/9.31 · Y1 6.54/29.88/6.44 · Y2 12.85/51.13/12.45.
→ x+u ≈ x_only (diffs 0.1–0.4 = noise); u_only bad for ΔAge (17–51).

Fate PR-AUC — x_only / u_only / x+u (per fold):
N3 0.997/0.984/0.997 · O1 0.971/1.000/0.968 · O2 0.984/0.923/0.984 · Y1 0.639/0.831/0.639 · Y2 0.964/0.988/0.964.
→ x+u = x_only every fold — BUT u_only is GOOD for fate (0.83–1.00), and on Y1 u_only 0.831 BEATS
x_only 0.639. So the perturbation carries real FATE signal that x+u fails to pick up.
(6× LinAlgWarning, rcond ~9e-8 — one per fold; cause explained in CORRECTION.)

**Reading (as measured): adding u to x does not beat x alone (x+u ≈ x_only), for BOTH targets.**
On its face this undercuts the "jointly models cell×perturbation" framing — but see the CORRECTION:
for fate, u *alone* is strong, so x+u ≈ x_only is partly a linear-concat dilution artifact, not proof.

> ### ⚠️ CORRECTION (2026-07-13) — the "perturbation is constant" explanation was WRONG.
>
> **RETRACTED claim:** that u_only failed / x+u ≈ x_only *because* "the perturbation is ~CONSTANT"
> and "the perturbation doesn't VARY / is UNTESTED here." **This is false**, verified against the
> actual bundles in this repo:
> - `dose_time = [log10(dose), log(time_h)]` (see `local_runners/diag_harmonize.py`).
> - The raw Gill data spans **~12 timepoints** (Fib/d0, d7, d9, d11, d13, d15, d21, d29, d34, d40,
>   d47, d54) × 2 markers (CD13/SSEA4) × 2 experiments — a dense time course, not a constant.
> - The fitted `dose_time` scaler std is **nonzero on BOTH dims across every fold**: dose-dim
>   std ≈ 1.29, **time-dim std ≈ 3.26** (`runs/*/scalers.json`). Time is captured, not collapsed.
> - Test 11's `u_only` = `[fingerprint, dose_time]` (`test_suite.py`) **already includes** this
>   varying time, yet still failed on ΔAge.
>
> **What the LinAlgWarning (rcond ~9e-8, one per fold) actually is — and isn't.** It is NOT collapsed
> time. It is the **constant fingerprint**: every reprogrammed cell got the same cocktail (OSKM), so
> the `fp` columns of `u` are identical cell-to-cell (and `x_std` shows some panel genes are constant
> too — min std 0). That collinearity is what ill-conditions the ridge matrix. Under StandardScaler
> the constant `fp` columns zero out, so u_only is effectively ridge on the **2** varying `dose_time`
> columns → predicting ΔAge from [dose, time] alone, which is genuinely weak because the same day
> gives different ΔAge across the 6 donors (time needs cell-state context). Decisive-check outcome:
> `dose_time` is **NOT** collapsed (std [1.29, 3.26]) → the "pipeline threw the timing away"
> hypothesis is **ruled out**.
>
> **What still stands (the measurements):** x+u ≈ x_only, and u_only is terrible for ΔAge — those
> numbers above are real. The interaction adds nothing *as measured*.
>
> **The corrected reading:** since the perturbation (time) DOES vary and IS available to the model,
> "constant perturbation" cannot be the cause. The real, still-open question is WHY the varying
> time doesn't help ΔAge in this ablation. Two live hypotheses, not yet decided:
> - **(A) Redundant** — cell state already encodes time (a day-20 cell *looks* day-20), so ΔAge is
>   predictable from state and time-as-a-separate-input is redundant. Benign; no modeling bug.
> - **(B) Under-used** — the varying time carries ΔAge signal that state does NOT fully contain,
>   but the model (and/or the linear u_only ablation) underweights it. A real, fixable modeling issue.
>
> **New (2026-07-13, from the full console output): fate breaks the pattern — and it implicates the
> LINEAR TEST, not necessarily the net.** For fate, u_only is NOT weak: O1 1.000, Y2 0.988, and Y1
> **0.831 vs x_only's 0.639** — so the perturbation carries real signal. Yet x+u = x_only every fold
> (Y1 stays 0.639). Reason: concatenating **2** `dose_time` columns with **2000** genes lets a linear
> model drown the low-dimensional perturbation signal. So "interaction not load-bearing" is, at least
> in part, a **feature-dilution artifact of the linear concat** — NOT evidence the signal is unusable,
> and NOT a verdict on the trained net (which has its own perturbation encoder and is never called by
> Test 11). Whether the **net** uses the signal is decided by **Test 11.1 part 3** (hold the cell
> fixed, sweep the time input into the net).
>
> **What is genuinely unchanged:** the model's core thesis — predicting outcomes ACROSS *different*
> perturbations (factor combinations / doses) — remains untested here, because this dataset varies
> **time under one cocktail (OSKM)**, not the cocktail itself. Only the false "nothing varies"
> mechanism is retracted; the data–capability scoping stands. **Test 11.1 (below) decides A vs B.**

### Test 11.1 — ISOLATE THE TIME DIMENSION (does time carry ΔAge signal beyond cell state, and does the model use it?)  [DESIGNED — pending data to run]

**Context.** Follows the Test 11 correction (2026-07-13): time varies richly and is encoded in
`dose_time`, yet `u_only` failed and `x+u ≈ x_only` for ΔAge. This test decides between
**(A) time is redundant with cell state** (benign) and **(B) time carries ΔAge signal the model
under-uses** (a real, fixable modeling issue).

**Hypothesis (and why).** Leaning slightly toward (A): in a reprogramming time-course a cell's
expression state is largely a readout of how far along it is, so "time" and "state" are heavily
correlated and state should absorb most of time's ΔAge signal. This is NOT confident — (B) is fully
live, because the model could be riding state out of convenience while ignoring a real,
partially-independent time signal.

**Prediction BEFORE running (honest: no confident direction).** I do not have a confident numeric
prediction — that is the point of running it. The pre-registered *branches*:
- If time explains ≈0 of the ΔAge residual-after-state AND ΔAge is ~flat vs time within fixed state
  → **(A) redundant.** Expect partial-R²(resid ~ time) < ~0.03.
- If time explains meaningful residual ΔAge variance (partial-R² ≳ 0.1) → state does NOT fully
  contain time → **(B)**, pending the model-sensitivity check.

**Method (three parts; per-fold, train-fit / test-measured, nothing tuned on results).**
1. **Variance decomposition (model-free).** R²(ΔAge ~ time-only), R²(ΔAge ~ state X), and the
   **partial** R²(residual-of-ΔAge-after-X ~ time). The partial term is decisive: it is the ΔAge
   signal in time that state does NOT already carry.
2. **State-controlled time effect.** Within narrow cell-state neighborhoods (kNN in expression
   space), test whether ΔAge still varies with time. Nonzero within-neighborhood slope ⇒ time
   carries state-independent signal.
3. **Model sensitivity (does the TRAINED model use time).** On held-out cells, hold X fixed and
   sweep the `dose_time` time-input across its observed range; measure how much the model's ΔAge
   output moves. Flat ⇒ model ignores time (supports B); responsive ⇒ model uses it, and then the
   `u_only` failure is a property of the *linear* ablation, not the model.

**What each outcome concludes.**
- **(A) redundant** (parts 1–2 ≈ 0): Test 11's "interaction not load-bearing" is benign and
  substantively correct — the model rightly rides the richer signal (state), which already encodes
  time. No modeling bug; nothing to fix beyond honest scoping.
- **(B) under-used** (parts 1–2 show residual time signal, part 3 flat): real issue — the model
  wastes an available, informative signal. Action: strengthen the perturbation/time encoder.
- **Ablation artifact** (part 3 responsive despite u_only failing): the u_only=terrible result
  reflects the *linear* ablation, not the model ignoring time. Re-state Test 11 accordingly.

**Status: DESIGNED, NOT YET RUN.** Blocked on data — the repo ships trained bundles (`runs/…`) but
**not** the per-cell dataset (X, `dose_time`, ΔAge labels). To run this I need either the built
dataset shards / `manifest.parquet` for the Gill/HFF folds, or the raw Gill counts + config to
rebuild via the pipeline. Bundles for part 3 are already present.

**Results (actual).** _[TO FILL — once data is provided]_

### Test 12 — PER-DONOR JACKKNIFE (are aggregates robust or hostage to specific donors?)  [distinct question]  ✅ DONE
**Vary only: which donor is excluded from the aggregate.** Recompute the headline metrics
dropping one donor at a time (jackknife over the 6). **Question:** are the poor aggregates driven
by 1-2 outlier donors (e.g. N3, +30 ΔAge, 0 unsafe cells) or broadly true? **Prediction:** N3 (and
maybe N2) dominate the bad numbers; dropping them changes the story. **What it pinpoints:** whether
conclusions are robust or hostage to specific donors — essential before any claim in the writeup.

#### Test 12 — RESULT (per-donor jackknife) — read from RAW per-fold numbers, not the script verdict

Per-fold: ΔAge MAE N2=21.79 N3=29.69 O1=5.39 O2=7.54 Y1=7.28 Y2=14.06 (ALL=14.29).
Fate PR-AUC N2=n/a(0 unsafe) N3=0.997 O1=1.000 O2=1.000 Y1=0.961 Y2=1.000.

**ΔAge MAE is bimodal, inflated by 2 atypical donors.** Good folds cluster 5.4–7.5; bad folds
N2=21.8, N3=29.7 are 3–4× worse. Dropping N3 → 11.21 (−22% off the 14.29 headline); dropping N2
→ 12.79; dropping a GOOD fold makes it worse (drop O1 → 16.07). So the "MAE~14" story is really
"~6 for typical donors, ~25 for the 2 atypical ones (N2 zero-unsafe, N3 +30 outlier)." Quantified
donor-heterogeneity (H6): remove the 2 weird donors and MAE ≈ 10, not 14.

**SURPRISE — fate DISCRIMINATION is near-perfect out-of-donor (PR-AUC ~0.96–1.00 every fold).**
This PARTIALLY CONTRADICTS the earlier "fate fails out-of-donor" claim. High PR-AUC means the
classifier RANKS safe-vs-unsafe almost perfectly on held-out donors — discrimination is NOT the
problem. (Caveat: ~21 cells/fold with few unsafe → PR-AUC=1.0 is easy to hit and may overstate;
the PATTERN of strong discrimination is still real.) This means the fate issue (ECE 0.26) is most
likely CALIBRATION, not discrimination — fate ranks right but reports miscalibrated probabilities.
Strong preview of Test 8.2; and it means RES's failure (multiplying by fate) is plausibly fixable
by RECALIBRATING fate, since the fate ranking it needs is actually good.

**Process note (user correction):** read effect sizes from raw per-fold numbers; the script's
one-line verdict uses a hardcoded threshold and can mislead. Applies to all tests.

---


### Test 7.4 — does recalibration help RES on a NON-rank metric? (7.3 used rank-invariant Spearman)  ⏳

**Why (user caught a real flaw in 7.3).** 7.3 measured ranking via **Spearman**, which is
**rank-invariant**. Platt recalibration is a **monotonic** transform — it rescales probabilities
but preserves their order — so RES's *rank* correlation is nearly unchanged **by construction**
(C−B was literally +0.000 on 5/6 folds, +0.014 on N3). That is NOT evidence recalibration is
useless; it means 7.3 tested the one axis recalibration cannot move. Recalibration's real benefit
would appear in **threshold / absolute-score** quality, which Spearman ignores.

**Hypothesis.** Recalibrating fate makes RES a better *actionable score* (better precision at a
decision threshold, and better-calibrated as a score) even though it can't change rank order.

**Prediction (before running, on record).** C > B on the threshold/calibration metrics (recal
helps where it structurally can) — vindicating the concern. Whether C beats A (plain ΔAge sort)
on these metrics is genuinely open; low confidence.

**Method.** Same A / B / C (ridge ΔAge / RES-raw / RES-recal), same held-out cells, but score with
metrics that are NOT rank-invariant:
- **precision@k and recall@k** for retrieving truly safe-AND-rejuvenating cells (top-k by score);
- **score calibration** of RES itself where meaningful;
- report per-fold + aggregate + paired (C−B) and (C−A) on each.

**Decision branches (before data):**
- **C > B on threshold metrics** → recalibration DOES help RES as an actionable score; 7.3's null
  was a metric artifact → report RES quality on calibrated-score terms, not rank.
- **C ≈ B even here** → recalibration genuinely doesn't help RES on any axis → 7.3's conclusion
  stands and is now robust.
- **C ≥ A** → recalibrated RES matches/beats a plain ΔAge sort on the actionable metric → RES
  salvageable for thresholded triage (even if not for ranking).

**Result (actual).** Ran on real data. Only **4/6 folds usable** (N3, Y2 degenerate — no
variation in "truly safe AND rejuvenating" among their test cells).

SCORE CALIBRATION error (lower better):

| fold | A ridge | B RES raw | C RES recal | C−B |
|---|---|---|---|---|
| N2 | 0.319 | 0.762 | 0.744 | −0.018 (helped) |
| O1 | 0.208 | 0.095 | 0.088 | −0.007 (helped) |
| O2 | 0.105 | 0.143 | 0.143 | +0.000 (nothing) |
| Y1 | 0.249 | 0.104 | 0.132 | **+0.028 (HURT)** |

paired C−B = +0.001, 95% CI [−0.031, +0.032] (n=4) → **tied**.

PRECISION/RECALL at the calib cutoff: **identical B vs C on every fold**, C−B = +0.000,
CI [+0.000, +0.000] (n=3).

**DESIGN ERROR IN THIS TEST (own it).** The precision/recall half is **structurally void** — the
same flaw as 7.3 in a new costume. I defined the cutoff as a **calib-set QUANTILE**, but a
quantile is **monotone-equivariant**: Platt shifts the scores and the 70th percentile shifts
identically, so **exactly the same cells are flagged** by construction. Verified explicitly. The
zero-width CI [+0.000,+0.000] is the fingerprint of a forced result, not a measurement. A genuine
test needs a **fixed absolute value** cutoff on the native RES scale (→ Test 7.4.1).

**Verdict — prediction WRONG again (I predicted C > B on value metrics; it is tied).** On the ONE
valid half (score calibration, verified value-sensitive), recalibration moves RES in **both
directions**: helps on N2/O1, nothing on O2, **hurts on Y1 by more than it helped anywhere**. That
is noise, not a small consistent benefit.

**This directly tests the user's hypothesis** ("recal always changed the result for good, just at
small scale"). On a metric that *can* see monotonic changes, **it does not always help**. The
always-≥0 pattern in 7.3 was itself an artifact of the rank-blind metric.

**Status of the RES question:** two different metric families (rank in 7.3, value-calibration in
7.4) both show no systematic recalibration benefit → the "drop RES for ranking" conclusion is
broader-based than before. **Caveats:** n=4 only, and half of 7.4 was void — not bulletproof.
**Side finding:** RES calibration is wildly fold-dependent (0.095 on O1 vs **0.762** on N2) — RES
is not uniformly miscalibrated, it is catastrophically miscalibrated on specific donors.

### Test 7.4.1 — PATCH: fixed ABSOLUTE cutoff + the APPROVED gate  ⏳

**Why.** 7.4's decision cutoff was a calib QUANTILE — monotone-equivariant, so Platt could not
change which cells were flagged (verified: identical sets). RES is documented as living in
**[0, 1)**, so a FIXED numeric cutoff is meaningful and IS sensitive to recalibration (verified
synthetically: 5–9 cell decision flips at cutoffs 0.05–0.50, vs 0 flips for the quantile).

**Extra lever found in the source:** `compute_res_batch` also returns a **status**, and
`REJECTED_UNSAFE` fires on the ABSOLUTE test `S < tau_safe - 3*w`. Platt changes S, so it moves
this gate — this is the product's real accept/reject decision, not a proxy.

**Method.** Same isolation (ΔAge held = ridge; only fate calibration differs).
B = RES(raw fate), C = RES(recal fate). Report at fixed cutoffs {0.05, 0.10, 0.20, 0.30, 0.50}:
cells flagged, precision/recall for truly safe-AND-rejuvenating cells; plus the same for the
APPROVED gate. **Built-in sensitivity check** (7.4's lesson): the script reports how many
cell-level decisions actually FLIP between B and C, and refuses to report "tied" if the answer
is zero — it declares itself blind instead.

**Prediction (before running, on record).** The sets WILL differ (recal shifts S, so the fixed
gate moves), but **precision will be tied** — recalibration moves the decision boundary without
improving the decision, because the damage is the multiplicative formula, not its inputs (the
7.3/7.4 pattern). Low-moderate confidence; I have been wrong on 7.3 and 7.4 predictions.

**Scope limit (state up front):** 7.3 already showed recalibrated RES ranks far BELOW a plain
ΔAge sort (C−A = −0.298, CI excludes 0). This test can only change *how much* RES
underperforms, not whether it does.

**Result (actual).** Ran on real data. 4/6 folds usable (N3, Y2 skipped — no safe-rejuv
variation). Sensitivity check PASSED: 18 cell-level decision flips → the metric can see recal.

**THE HEADLINE IS THE SCORE-RANGE TABLE:**

| fold | B RES range (raw) | C RES range (recal) | good cells |
|---|---|---|---|
| N2 | [0.000, 0.000] | [0.000, 0.232] | 16/21 |
| O1 | [0.000, 0.000] | [0.000, 0.121] | 2/21 |
| O2 | [0.000, 0.000] | [0.000, 0.000] | 3/21 |
| Y1 | [0.000, 0.007] | [0.000, 0.535] | 2/19 |

**RAW RES IS DEGENERATE — collapsed to ~0 for every cell on every fold.**

APPROVED gate (the product's real accept/reject decision):

| fold | n approved B | n approved C | prec C | rec C |
|---|---|---|---|---|
| N2 | **0/21** | 3/21 | 1.00 | 0.19 |
| O1 | **0/21** | 1/21 | 1.00 | 0.50 |
| O2 | **0/21** | 0/21 | — | 0.00 |
| Y1 | **0/19** | 2/19 | 0.50 | 0.50 |

**Raw RES approves ZERO cells on every fold** — as shipped it rejects 100% of candidates
out-of-donor. Recalibrated RES approves 6 cells total, 5 of them truly safe-and-rejuvenating.

**SCRIPT BUG (own it):** every `paired C-B precision ... (n=0) -> tied` line is an ARTIFACT.
Precision for B is n/a everywhere because B flagged nothing, so there were zero pairs; with n=0
the code falls through to "tied". It should read "cannot compare — B flagged nothing." Ignore
those verdict lines; the score-range and APPROVED tables carry the information.

**Verdict — prediction WRONG (3rd time); the user was RIGHT, and more specifically than claimed.**
User's words: "it always changed the result for good — just the scale isn't good." That is exactly
the finding: **the SCALE was the whole problem.** Raw RES sits at ~1e-4 or below; recalibration
lifts it into a usable [0, 0.5] range.

**MECHANISM (deduced, not guessed).** RES = Φ(S)·S^k·g(R_eff)·exp(−λ·P_loss), with
Φ(S)=sigmoid((S−τ_safe)/w). Recalibration changes ONLY S and P_loss — it never touches µ_age or
σ_age, so g(R_eff) is IDENTICAL in both arms. Since C reaches 0.535, g must be > 0. Therefore the
collapse in B comes entirely from the SAFETY terms: **raw out-of-donor S sits far below τ_safe,
and the double penalty (sigmoid floor × S^k) annihilates the score.** This ties directly to Test
8.2 (fate systematically under-confident out-of-donor, ECE 0.28): **RES's safety gate is
mis-tuned relative to the model's out-of-donor calibration.**

**REVISED RES CONCLUSION (supersedes the 7.3/7.4 reading):**
- **As a GATE/decision score: raw RES is broken (approves nothing); recalibration is not optional,
  it is REQUIRED for RES to function at all.** This is a real, actionable fix.
- **As a RANKING: unchanged — RES still loses to a plain ΔAge sort** (7.3: C−A = −0.298, CI
  excludes 0). Ranking and gating are different questions; recalibration fixes the second, not the
  first.

**CAVEATS (do not overclaim):** n=4 folds; approvals are 1–3 cells per fold (precision 1.00 on 3
cells is weak evidence); recall is LOW (C finds only 3 of N2's 16 good cells); O2 stays collapsed
even after recal. This shows RES becomes *functional*, not that it becomes *good*.

### Test 7.4.2 — MEASURE the mechanism: decompose RES into its four factors  ⏳

**Why.** 7.4.1's mechanism (raw S below tau_safe → safety floor annihilates RES) was **deduced**,
not measured. This measures it by decomposing
`RES = Phi(S) * S^k * g(R_eff) * exp(-lam*P_loss)` and reporting each factor's magnitude, raw vs
recalibrated. Phi and S^k change under recalibration; **g is IDENTICAL in both arms** (it depends
only on mu_age/sigma_age), so whichever factor is ~0 in the raw arm names the defect.

**Four candidate diagnoses (pre-registered):**
1. **Phi(S) ~ 0** → out-of-donor S below tau_safe; narrow safety floor kills RES. CALIBRATION defect.
2. **g(R_eff) ~ 0** → sigma_age so large the upper age bound is never negative; never *confidently*
   rejuvenating. UNCERTAINTY defect — recalibration could not have helped (but it did, so unlikely).
3. **in_dist mostly False** → the OOD gate rejects everything. OOD defect.
4. **exp(-lam*P_loss) ~ 0** → risk penalty dominates (only possible if lam > 0).

Also reports the **STATUS breakdown** (APPROVED / REJECTED_OOD / REJECTED_UNSAFE /
REJECTED_NO_REJUVENATION) — the rejection reason named directly.

**Pre-flight checks done (sandbox):**
- The decomposition reproduces `compute_res_batch` **exactly** (verified numerically).
- With defaults `tau_safe=0.85, w=0.03, k=2.0`, the floor is razor-sharp:
  S=0.60 → Phi·S^k = **8.65e-05** (displays as "RES = 0.000"); S=0.85 → 0.361.
  **A 0.25 shift in S changes RES by ~4000x.** So mild under-confidence is sufficient to
  annihilate RES — the 7.4.1 mechanism is numerically plausible before we even run it.
- **`lam` defaults to 0.0** → `exp(-lam*P_loss) = 1` always → the P_loss term is **INERT**, so
  recalibrating P_loss in 7.3/7.4/7.4.1 was a **no-op by construction**. The test reports whether
  this holds in the actual bundles. (If so: only S ever mattered.)

**Prediction (on record).** Diagnosis 1: Phi(S) median ~1e-3 or smaller in the raw arm while
g stays healthy (>0.1), most raw cells REJECTED_UNSAFE, and lam == 0 confirming the P_loss no-op.
Moderate confidence — but note this prediction is *entailed* by 7.4.1's result (recal helped, and
recal cannot touch g), so it is closer to a consistency check than an open question.

**Result (actual).** Ran on real data, all 6 folds.
`RES PARAMS: tau_safe=0.85 w=0.03 k=2.0 kappa=5.0 z_conf=1.0 lam=0.0` → **lam=0 CONFIRMED**:
the P_loss term is inert, so recalibrating P_loss in 7.3/7.4/7.4.1 was a **no-op by
construction**. Only S ever mattered.

| fold | Phi B | Phi C | S^k B | S^k C | **g (both)** | med mu_age | med sigma | med R_eff |
|---|---|---|---|---|---|---|---|---|
| N2 | 5.00e-04 | 1.97e-01 | 0.387 | 0.653 | **0.000** | **+8.76** | 2.23 | 0.00 |
| N3 | 9.72e-03 | 8.95e-01 | 0.506 | 0.836 | **0.000** | **+4.60** | 2.54 | 0.00 |
| O1 | 4.60e-03 | 7.78e-01 | 0.474 | 0.788 | **0.000** | **+18.14** | 2.69 | 0.00 |
| O2 | 1.20e-02 | 8.95e-01 | 0.515 | 0.836 | **0.000** | **+20.44** | 2.39 | 0.00 |
| Y1 | 2.09e-03 | 6.45e-01 | 0.442 | 0.753 | **0.000** | **+11.69** | 2.62 | 0.00 |
| Y2 | 6.26e-02 | 9.45e-01 | 0.591 | 0.875 | **0.000** | **+21.12** | 2.14 | 0.00 |

**Verdict — prediction WRONG (4th time). My 7.4.1 mechanism ("safety floor is the killer") was
INCOMPLETE.** The real structure is TWO STACKED BLOCKERS:

1. **PRIMARY — no predicted rejuvenation.** `med mu_age is POSITIVE on every fold` (+4.6 to
   +21 yr). `R_eff = max(0, -(mu + z*sigma))` needs mu < -2.4 for ANY credit; we are 7–23 years
   away. So **g = 0 and RES is zero BEFORE the safety floor gets a vote.** Recalibration
   **cannot** touch this (g depends only on age terms).
2. **SECONDARY — the safety floor.** Phi is 5e-4–6e-2 raw, 0.20–0.95 recalibrated. Recalibration
   DOES fix this — which is why 7.4.1 saw a *few* approvals (the minority of cells that do have
   R_eff > 0; medians hide that tail).

**The STATUS table proves the handoff.** Status order is OOD → UNSAFE → NO_REJUVENATION →
APPROVED. N3 raw: 15 UNSAFE, 0 no-rejuv. N3 recal: 4 UNSAFE, **6 no-rejuv**. Cells did not become
approved — they **moved from "rejected as unsafe" to "rejected as not rejuvenating."**
Recalibration lifted them past gate 1, exposing gate 2 behind it.

**THE UNIFYING HYPOTHESIS (new, important).** `mu_age` positive while ranking is excellent
(Spearman 0.95) is the signature of a **systematic OFFSET**, not noise. One offset would explain
three separate results simultaneously:
- ΔAge MAE ~14 yr — an offset directly inflates MAE;
- ranking works (0.95) — rank correlation is offset-INVARIANT;
- RES is dead — RES uses the ABSOLUTE value, so an offset zeroes the rejuvenation credit.

**COMPETING EXPLANATION — must distinguish before acting.** Is positive mu_age model BIAS or
biological TRUTH?
- **N2 says BIAS:** 16/21 cells are truly safe-and-rejuvenating (true ΔAge < 0) yet ridge
  predicts +8.76.
- **O1 may say TRUTH:** only 2/21 are truly safe-and-rejuvenating, so a positive prediction may
  be CORRECT there.
If the cells genuinely are not rejuvenating by the Fleischer clock, **RES is correctly approving
nothing** and the defect is in the DATA/clock (cf. the epistemic-extrapolation warning: clocks
trained on natural aging, applied to reprogramming), not in RES. → **Test 7.4.3** must compare the
TRUE ΔAge distribution against the predicted one, per fold.

**SIDE FINDING (OOD preview).** The OOD gate flags 34/124 held-out cells (~27%) — 4–10 per fold.
Unexamined until now; relevant to the pending OOD-detector test.

### Test 7.4.3 — BIAS or TRUTH? (the decisive discriminator)  ⏳

**Why.** 7.4.2 proved RES dies because predicted ΔAge is positive (g=0), but could not say WHY.
The two causes need OPPOSITE fixes: model **BIAS** (fix the offset → RES revives) vs biological
**TRUTH** (cells really do not rejuvenate by this clock → RES is correctly rejecting; the defect
is the DATA/CLOCK and fixing RES is the wrong target).

**Method (3 parts).**
1. TRUE vs PREDICTED ΔAge distributions per fold + the offset. Large positive offset WITH good
   rank correlation = bias signature.
2. **ORACLE RES** — feed RES the **TRUE** ΔAge (model's sigma, recalibrated fate, real OOD flags).
   **This is decisive.**
3. **OFFSET-CORRECTED RES** — subtract each fold's offset, estimated on **CALIB only** (never on
   the held-out donor → leak-free; will UNDER-correct if the offset is donor-specific).

**Discriminator verified in sandbox** (both scenarios give PRED=0, matching the real data, but
the oracle column separates them):
- BIAS scenario (truth rejuvenates, prediction offset +12): PRED=0, CORRECTED=19, **ORACLE=19**.
- TRUTH scenario (cells genuinely age): PRED=0, CORRECTED=0, **ORACLE=0**.

**Prediction (on record).** **MIXED, leaning BIAS on N2 and TRUTH on O1/O2/Y2.** Rationale: N2 has
16/21 truly safe-and-rejuvenating cells yet a +8.76 prediction (bias); O1 has only 2/21 truly
safe-and-rejuvenating (truth). Expect oracle approvals clearly > predicted on N2, and ≈ 0 on the
high-mu folds. Moderate confidence — my last four predictions in this thread were wrong.

**Decision branches.**
- **ORACLE >> PRED** → BIAS-dominated → RES logic sound, ΔAge prediction is the bottleneck →
  fixing the offset revives the product AND explains MAE~14 with ranking~0.95.
- **ORACLE ≈ PRED and TRUE-eligible fraction low** → TRUTH-dominated → RES correctly approves
  nothing; the honest finding is that **this dataset does not exhibit clock-measurable
  rejuvenation**, so no scoring fix can help. This would be the single most consequential result
  in the whole investigation — it would mean the product claim needs different DATA, not a better
  model.
- **MIXED** → read per-fold: high `frac true<0` + low `frac pred<0` = bias; both low = truth.

**Result (actual).** Ran on real data, all 6 folds.

| fold | med TRUE | med PRED | **error** | offset(calib) | frac true<0 | spearman | appr PRED | appr ORACLE |
|---|---|---|---|---|---|---|---|---|
| N2 | −11.35 | +8.76 | **+20.11** | −0.02 | 0.76 | +0.957 | 3/21 | 4/21 |
| N3 | +29.00 | +4.60 | **−24.40** | −0.01 | 0.00 | +0.925 | **7/21** | **0/21** |
| O1 | +12.43 | +18.14 | **+5.71** | −0.06 | 0.24 | +0.960 | 1/21 | 2/21 |
| O2 | +7.40 | +20.44 | **+13.04** | +0.00 | 0.24 | +0.952 | 0/21 | 3/21 |
| Y1 | +15.96 | +11.69 | **−4.27** | −0.06 | 0.32 | +0.951 | 2/19 | 2/19 |
| Y2 | +29.92 | +21.12 | **−8.80** | −0.05 | 0.14 | +0.983 | 1/21 | 0/21 |

**totals: PRED=14, CORRECTED=14, ORACLE=11. Median calib offset = −0.03 yr.
Mean TRUE-eligible fraction = 0.27.**

**Verdict — my "systematic offset" hypothesis is WRONG as stated, and the truth is worse.**

**(1) There is NO global bias.** The calib-estimated offset is −0.03 yr — essentially zero. The
predictor is **unbiased in-distribution**. Hence CORRECTED = PRED exactly (14=14) and MAE is
unchanged (21.59→21.61). Part 3 was a null **by construction**: there was nothing to correct.

**(2) The error is DONOR-SPECIFIC, huge, and SIGNED BOTH WAYS — so it cancels.**
Per-donor error (med PRED − med TRUE): N2 **+20.11**, N3 **−24.40**, O1 +5.71, O2 +13.04,
Y1 −4.27, Y2 −8.80.
**Mean = +0.23 yr (cancels). Mean ABSOLUTE = 12.72 yr. Std = 14.71 yr.**
This is hypothesis **H6 (donor heterogeneity) quantified precisely** — not a bias that can be
subtracted, but a per-donor level shift that is **unknowable from calib** (calib is
in-distribution, where the model is fine). No global correction can ever fix this.

**(3) ORACLE (11) is LOWER than PRED (14) — the model OVER-approves.** Feeding RES the TRUE ΔAge
yields FEWER approvals than feeding it predictions. So RES's approvals are not merely sparse,
they are **wrong in composition**.

**(4) N3 is a FALSE-POSITIVE FACTORY — the safety-relevant failure.** True median ΔAge **+29.00**,
**frac true<0 = 0.00, TRUE-eligible = 0.00** — not a single cell genuinely rejuvenates. Yet the
model predicts +4.60 and **approves 7 cells**. The system would recommend 7 conditions that
actually AGE the cells. (And N3's Spearman is still 0.925 — it ranks correctly while being
catastrophically wrong about the level.)

**(5) NOT truth-dominated: real rejuvenation IS present.** Mean TRUE-eligible = 0.27; N2 is 0.76.
So the clock does register genuine rejuvenation in this data — the earlier "maybe nothing
rejuvenates" worry is ruled out.

**THE UNIFYING PICTURE (replaces the offset hypothesis).** The model has **excellent
within-donor ordering (Spearman 0.925–0.983 on every fold) but badly shifted per-donor absolute
levels (±12.7 yr typical).** That single fact explains everything observed:
- ranking ~0.95 ✓ (rank is invariant to a per-donor level shift)
- ΔAge MAE ~14 ✓ (the level shift IS the error)
- RES dead / wrong ✓ (RES thresholds are ABSOLUTE, so a shifted level makes them meaningless)

**ACTIONABLE CONSEQUENCE.** The honest product statement: the model is usable for
**within-donor RANKING** of conditions (which is the realistic use case — rank conditions for one
patient), and is **NOT reliable for absolute-threshold decisions out-of-donor** (RES approval).
The fix direction is **per-donor calibration**: with a few labelled reference samples from a new
donor, the level shift becomes estimable and absolute decisions become possible. Without them,
report ranks, not approvals.

**Prediction scorecard:** I predicted "MIXED, bias on N2, truth on O1/O2/Y2" — roughly right on
the per-fold pattern (first decent prediction in a while), but I **missed the two biggest facts**:
that the errors CANCEL (no global offset) and that the model **over-approves** (oracle < pred).

## WHERE THE INVESTIGATION STANDS — canonical current summary (updated 2026-07-13)

*Single source of truth for the project's state. The earlier blocks ("FINAL DIAGNOSIS", "Later
planned tests", "THE STRUCTURAL INSIGHT THAT CLOSES THE CASE") predate the Test 5 reopening and are
**superseded** by this section + the decision tree — read those as history, not current status.*

### What the model is, and what it must prove
CellFate-Rx maps (cell state X, perturbation u = fingerprint + dose_time) → 3-way fate (safe / loss /
death) + ΔBiological-age + uncertainty, and folds them into a Rejuvenation Efficacy Score (RES) for
safe-reprogramming triage. Its reason to exist is **perturbation response** — predict how the outcome
*changes* as you vary the perturbation. Validation data = HFF + Gill: an **OSKM reprogramming
time-course over 6 donors** (N2, N3, O1, O2, Y1, Y2).

### Established (with confidence)

**1. ΔAge magnitude — SETTLED; the NN tying ridge is correct, not a failure.** The Fleischer clock is
linear (age = w·x + b) and every transform is affine, so ΔAge = w·(x_pert − x_ctrl) is linear in
expression *by construction*. Ridge is optimal; the NN can only match it — and Test 3 proves the
architecture *would* beat ridge if nonlinearity existed (+6.64 on synthetic). Confirmed three ways:
statistically tied not worse (Test 5.0, paired CI [−2.16, +2.64]); no nonlinear model beats ridge on
real ΔAge (Test 6 — trees 15.63 / forest 17.81 / kernel 24.68 all worse than ridge 14.05, flexibility
*hurts*); and the clock's structural linearity. In-dist MAE ~4 yr (Test 8.1); the ~14 aggregate is
inflated by 2 atypical donors — drop N2/N3 → ~10, typical-donor MAE ≈ 6 (Test 12). High *absolute*
MAE is explained by input coverage (panel sees only ~47% of the clock's signal; 88 of its top-200
genes absent — Test 0) and the per-donor control offset (Test 1) — both limit ridge equally. Lower MAE
would come from **more clock genes + more donors**, not a fancier model. Gene embeddings don't help
either (Test 9: ridge_emb worse, mlp_emb ties) — a third confirmation, after Tests 3 and 6, that ΔAge
is purely linear.

**2. Fate — the model's ONE clean edge over linear; discrimination real, calibration the fixable gap.**
Fate discrimination is GOOD and trustworthy: PR-AUC ~0.93 in-dist on *thousands* of held-out HFF cells,
holding out-of-donor 0.96–1.00 (Tests 8.1, 12). It is the single place the NN beats a linear baseline —
decisively on the hard Y1 fold (0.961 vs logreg 0.636; Test 8). But the probabilities are
**miscalibrated** (ECE ~0.28), decoupled from the good ranking = the cheap-to-fix kind; Platt recal
~halves ECE on 4/5 folds (Test 8.2). Calibration, not discrimination, is the fate problem. A richer
representation doesn't rescue it either: gene embeddings give at most a marginal, non-significant bump
(Test 9.1 — mlp_emb +0.030 with CI touching 0) — the fate limit is data/donors, not representation.

**3. Ranking / RES — RES does NOT earn its keep (one salvage path left).** The model's ΔAge predictions
rank rejuvenation beautifully (Spearman ~0.95 vs true ΔAge, ≈ ridge). The RES *score* destroys that,
dragging 0.95 → 0.69 (Test 7 — loses 6/6 vs a plain ΔAge sort). It loses even against
**safe-rejuvenation, its own objective** (Test 7.1), and with ΔAge held identical the RES *formula
itself* drops Spearman 0.955 → 0.654 on all 6 folds (Test 7.2, clean isolation). Cause is pinned:
RES multiplies ΔAge by fate probabilities that are miscalibrated out-of-donor, so it down-weights the
wrong cells (RES quality tracks *inversely* with the number of unsafe cells). → For ranking today,
**use ΔAge directly**; RES as formulated subtracts value. Only untested lever: recalibrate fate first
(Test 7.3).

**4. Generalization — the real practical limit is donor heterogeneity, not capability.** Fitting is
fine (in-dist ΔAge ~4, fate ~0.93); the failure is generalizing to *atypical* donors from a 6-donor
cohort — N2 (zero unsafe cells) and N3 (+30 ΔAge outlier) dominate every bad aggregate; typical donors
generalize well (Tests 8.1, 12). Fix = more donors / donor adaptation, not architecture.

**5. Input / perturbation — the CORE THESIS is UNTESTED here (the pivotal finding).** Test 11 (a
*linear* ablation that never touches the net) found x+u ≈ x_only. Corrected reading (2026-07-13): NOT
because the perturbation is constant — dose_time varies richly (12 timepoints; scaler std [1.29, 3.26]).
x+u ≈ x_only is partly a linear-concat dilution artifact (2 dose_time columns drowned among 2000 genes),
and u_only carries *real fate signal*. So Test 11 cannot indict the trained net. More fundamentally,
**this dataset varies only TIME under ONE cocktail (OSKM)** — never the cocktail itself. The model's
entire reason to exist (outcomes across *different* perturbations) is therefore **untestable on this
data**. The validation set was the wrong data to prove the core claim — the single most important
strategic result of this investigation.

### The bottom line
On this data the model *has* delivered: a ΔAge predictor that matches the linear optimum, a fate
classifier with a genuine (if not yet statistically significant) discrimination edge over linear, and
honest uncertainty/OOD machinery — fate discrimination being its one clean win. It has *not* delivered
its headline: RES actively hurts ranking, and the cross-perturbation thesis is untested because the
data can't pose it. Net: a sound, honestly-self-reporting system whose **core claim remains unproven
for want of the right data.**

### What needs to be done (prioritized)

**Tier 1 — close the open questions on THIS data (cheap; needs the per-cell data re-provided).**
- **Test 11.1** — does the trained net actually *use* the real, varying time signal, or waste it?
  Decides a real modeling bug (B) vs benign redundancy (A). Needs shards / `manifest.parquet`, or raw
  Gill counts + config to rebuild. (Bundles for the model-sensitivity part are already in the repo.)
- **Test 7.3** — RES with Platt-recalibrated fate: does fixing calibration (8.2) stop RES scrambling
  the ranking? The one path that could make RES earn its keep. (Uses existing bundles.)
- **Test 10 / 10.1** — full-transcriptome response; first needs a clean perturbed-vs-control target
  definition — flag if the data can't support clean pairing.

**Tier 2 — implement the fixes the diagnosis already points to.**
- Bake **per-donor fate recalibration** (Platt/temperature) into the pipeline, then feed RES from it.
- If lower ΔAge MAE is a goal: expand the panel toward the clock's high-signal genes + add donors.

**Tier 3 — the real validation (to put the actual thesis on trial).**
- Get data with **varied perturbations** (multiple factor cocktails / doses / drugs), not OSKM-only —
  the field's perturbation datasets (GenBio 2026 / Ahlmann-Eltze) have this — and re-run the acceptance
  gates (beats_all_baselines, ece_ok, coverage_ok, ranking_ok) on it. Only then are the cell×perturbation
  claim and RES actually on trial.
- Consider a **nonlinear single-cell aging clock**: it would make ΔAge nonlinear, giving the
  architecture (proven capable in Test 3) something a linear model cannot match.

---

## Test 13 — ΔAge TRAJECTORY SHAPE, and is the per-donor error time-varying?  ⏳ PLANNED

**Provenance (important).** An external review of the master plan asserted that the MPTR protocol
produces a **biphasic** trajectory (early transcriptomic stress reading as **+30 to +50 years**,
then late-phase rejuvenation), and that "fixing" the per-donor level shift would mathematically
erase that biology.

**That claim is NOT established anywhere in this notebook.** We have never analysed ΔAge as a
function of time. Partial consistency exists — Test 7.4.3 measured true median ΔAge of **+29.00
(N3)** and **+29.92 (Y2)**, and MPTR is the correct protocol for GSE165176 — but "biphasic"
requires showing the curve **comes back down**, which we have not measured, and "+50" appears
nowhere in our data. So it enters as a **hypothesis with a test attached**, not as a premise.

**Two questions, one test.**

**Q1 (biology).** Is the true ΔAge trajectory **biphasic** (rise then fall) or monotonic?

**Q2 (methodology — the one that actually gates the plan).** Test 16 (per-donor calibration)
assumes the per-donor error is a **constant scalar offset**. Is it? If the error instead **varies
with time**, a scalar correction is **misspecified** and would distort trajectories — subtracting
real biological signal rather than model error. **This matters regardless of whether Q1 is true.**

**Method (per donor, on its held-out time-course; time = `dose_time[:,1]` = log time_h).**
- **Part 1** — TRUE ΔAge vs time: fit linear + quadratic; a significant **negative** quadratic
  coefficient with an **interior vertex** = hump = biphasic.
- **Part 2** — PREDICTED ΔAge vs time, same fits: does the model reproduce the shape or flatten it?
- **Part 3** — **RESIDUAL (pred − true) vs time. THE TEST-16 GATE.**

**Prediction (on record, before running).**
- **Q1:** **MIXED, not a clean universal hump.** Donor heterogeneity is extreme (N2 is 76%
  truly-rejuvenating; N3 is 0%), and ~21 cells over ~12 timepoints is low-powered. I expect a hump
  on some donors at most.
- **Q2:** **the residual IS time-dependent on ≥2 folds — and CURVED rather than sloped.**
  Reasoning: within-donor Spearman is 0.925–0.983, high but **not 1.0**. A perfectly constant
  offset preserves rank *exactly* (Spearman = 1.0), so the shortfall implies genuine shape
  distortion. If the model flattens a humped trajectory, the residual humps.
  → **I expect Test 16 to need redesign.** Moderate confidence.

**Decision branches (set before data).**

| Part 3 result | Consequence for Test 16 |
|---|---|
| no significant slope **or** curve on ≥5/6 folds | **Test 16 as designed is VALID** — scalar per-donor correction is well-specified |
| significant slope **or** curve on ≥2 folds | **REDESIGN Test 16** — time-aware correction, or calibrate only on *matched timepoints*; never a single global scalar |

| Part 1 result | Consequence |
|---|---|
| hump on ≥half the folds | biphasic supported **on our data** → RES's absolute per-cell gating is *biologically* misspecified (it scores trajectory points in isolation) → state this in the writeup |
| hump on 0 folds | **do not build the plan on the biphasic claim** (low power: "not demonstrated" ≠ "disproven") |
| some folds only | donor-specific, not a universal protocol effect |

**Pre-flight verification (sandbox, before running on real data).**
- Shape detector correctly identified all four synthetic cases: hump (`quad_t = −9.96`),
  monotonic-down, monotonic-up, and no-trend.
- **Gate design flaw found and fixed:** the first version tested only the residual's **linear
  slope**. A **curved** residual (`slope_t = −0.05`, `curve_t = −7.65`) was wrongly reported as
  "~constant" — and that is *exactly* the shape produced if the model flattens a biphasic
  trajectory. The gate now fires on **either** a significant slope **or** a significant curve.

**Caveat carried in the output:** ~21 cells per donor. These fits are low-powered; a
non-significant quadratic means **"not demonstrated"**, not "monotonic proven".

**Result (actual).** _[TO FILL — user runs `python tests_13_16.py`]_

**Verdict.** _[TO FILL]_

---

## Test 14 — CONFORMAL INTERVAL VALIDATION  ⏳ PLANNED (never measured before)

**Why.** RES consumes `sigma_age` inside `R_eff = max(0, -(mu + z*sigma))`, and Test 7.4.2 showed
`R_eff = 0` on every fold. If the uncertainty is miscalibrated, that is an **independent upstream
cause** of RES's collapse, separate from the level shift. We have never measured it.

**Method.** Per fold: coverage = fraction of TRUE ΔAge inside `[mu - q, mu + q]`, compared with
the nominal conformal level. Also interval width, and width relative to the data's own spread
(a ratio >> 1 means the interval is wider than the variation it is meant to describe).

**Prediction (on record).** **Intervals will UNDER-cover out-of-donor**, and substantially.
Mechanism: `q` is fitted on the **calib** split, which is *in-distribution*, where ΔAge MAE is
~4 yr (Test 8.1). It is then applied to held-out donors where MAE is ~14 yr (Test 5/6). An
interval sized for a 4-year error cannot cover a 14-year error. **Moderate-high confidence** —
this one follows from numbers already measured.

**Decision branches.** UNDER-covers → uncertainty is a *second* independent cause of RES failure;
fix before judging RES further. CALIBRATED → uncertainty ruled out; the level shift stands alone.
OVER-covers → honest but uninformative; check the width/2sd ratio.

---

## Test 15 — OOD DETECTOR VALIDATION  ⏳ PLANNED (never measured before)

**Why.** `res = where(in_dist, res, 0)` — the OOD gate **zeroes RES outright**, and Test 7.4.2
found it flags **34/124 (~27%)** of held-out cells. Never validated. If it fires on cells that are
*not* actually erroneous, it is a third independent cause of RES's collapse.

**Method.** Per fold: OOD rate; mean \|error\| for flagged vs kept cells; and
**AUC(error → flagged)** — does the prediction error actually predict which cells get flagged?
AUC 0.5 = flags are random with respect to error.

**Prediction (on record).** **AUC ≈ 0.5–0.6 — uninformative to mildly informative.** Reasoning:
every held-out donor is out-of-distribution in some sense, so flagging exactly 27% looks like a
threshold artifact rather than a discovered property. **Low confidence** — genuinely unmeasured.

**Decision branches.** AUC > 0.6 → detector works, zeroing RES on those cells is defensible.
AUC < 0.45 → **misleading**: it discards cells that are not the erroneous ones → fix or disable
the OOD gate in the RES path. ≈ 0.5 → discards ~27% of cells for no measurable benefit.

---

## Test 16 — PER-DONOR CALIBRATION FEASIBILITY  ⏳ PLANNED  ***the gate for the main fix***

**Why.** Test 7.4.3 established a per-donor level shift (±12.7 yr) that cancels on average and is
invisible from calib. The only way to estimate it for a new donor is from **labelled reference
cells of that donor**. This test asks how many are needed — which decides whether the main fix is
a practical protocol change or an unusable one.

**Method.** Per fold, sweep k = 1, 3, 5, 10 with **40 random draws** each. Estimate the offset from
k reference cells, correct the **remaining** cells, and evaluate only on those — **reference cells
are always excluded from evaluation, so there is no leakage.** Two variants:
- **SCALAR** — one global offset from the k cells.
- **MATCHED** — per-cell offset from the reference cell nearest in **time** (handles a
  time-varying error, i.e. the case Test 13 Part 3 is testing for).

**PRE-REGISTERED CRITERION (from MASTER_PLAN §7b, fixed before running):**
\|shift\| reduced **≥50% on ≥4/6 folds at k ≤ 5** → PASS → implement.
Only at k ≥ 10 → BORDERLINE → document the cost, decide on practicality.
Not even at k = 10 → **FAIL → STOP**; report as a within-donor ranker.

**Prediction (on record).** **PASS at k = 3–5 on the scalar criterion, BUT the MATCHED variant
will be notably better** — indicating the error is time-varying and the implemented fix should be
time-aware rather than a single scalar. This is consistent with my Test 13 prediction (curved
residual). **Moderate confidence.**

**Pre-flight verification (sandbox).** The k-sweep logic was checked against two synthetic
regimes: with a genuinely CONSTANT offset, scalar reduced \|shift\| by **94–97%** and matched added
nothing; with a TIME-VARYING error, scalar plateaued at **60–74%** (4.6–7.2 yr residual) while
matched reached **1.3–2.4 yr**. So the scalar-vs-matched gap is a working detector of
misspecification, not decoration.

---

## Tests 13–16 — ONE RUN

All four are in `tests_13_16.py`; `python tests_13_16.py` runs the whole pre-change battery and
prints a FINAL SUMMARY mapping the results onto the master-plan decisions. **This is the last
testing step before code changes begin.**

---

## The decision tree (one glance) — updated 2026-07-13 (actual path + outcomes)

```
ΔAge MAGNITUDE — is the tie-with-ridge a failure or the ceiling?
Test 0   (H7 overlap)                        ✅ panel sees ~47% of signal -> explains high MAE, NOT the tie
Test 1   (overfitting? 6 folds)              ✅ NOT overfitting (train-Gill 8.88 ≈ test 14.29; O1/O2/Y1 test < train)
   └─ "regularize to beat ridge" path DEAD -> Test 2 = N/A
Test 3   (linear vs nonlinear, synthetic)    ✅ real ΔAge signal is ~LINEAR
Test 5.0 (is the 0.24 deficit real?)         ✅ NOISE -> model STATISTICALLY TIED with ridge (5.1 not needed)
Test 6   (can ANY model beat ridge, real?)   ✅ NO — trees/kernels worse; flexibility hurts
Test 9   (gene2vec embeddings on ΔAge)       ✅ no help (emb worse / mlp ties) — 3rd confirmation of linearity
   => ΔAge is linearly predictable; the NN tying ridge is CORRECT, not underperformance.   ► CLOSED

RANKING / RES — does RES earn its keep?
Test 7   (RES vs sort-by-ridge-ΔAge)         ✅ RES LOSES
Test 7.1 (RES vs safe-rejuvenation)          ✅ RES still LOSES (its own objective); cause = wrong out-of-donor fate
Test 7.2 (RES formula, ΔAge held constant)   ✅ the RES FORMULA itself degrades ranking
   => rank by ΔAge directly (~0.95); do NOT use RES for ranking.   ► CLOSED (salvage hinges on fate recal -> 7.3)

FATE — is fate where the model earns its keep?
Test 8   (model vs logreg)                   ✅ wins/ties every fold, wins big on hard Y1 -> the model's best claim
Test 8.1 (in-dist vs out-of-donor)           ✅ fitting FINE (~0.93 / ΔAge ~4); limit = generalization to atypical donors
Test 8.2 (discrimination vs calibration)     ✅ ranks well but MISCALIBRATED (ECE ~0.28); Platt recal ~halves it
Test 12  (per-donor jackknife)               ✅ aggregates hostage to 2 atypical donors (N2/N3); typical-donor MAE ≈ 6
Test 9.1 (gene2vec embeddings on fate)       ✅ marginal & non-sig -> representation is NOT the lever; data/donors are
   => fate DISCRIMINATION is real; CALIBRATION is the fixable gap -> unlocks a possible RES salvage

INPUT — does the model use cell × perturbation together?
Test 11  (input ablation)                    ✅ x+u ≈ x_only; u_only bad for ΔAge (but STRONG for fate)
   NOT a data problem (dose_time varies); "interaction not load-bearing" is partly a linear-concat dilution artifact
   └─ does the trained NET use the (real) time signal? -> Test 11.1

OPEN: Test 11.1 (time isolation, needs data) · Test 10 / 10.1 (full transcriptome) ·
      Test 7.3 (recalibrated-RES retest)
```

---
---

# PART II — CODE CHANGES BEGIN HERE

**Everything above this line is a MEASUREMENT.** Tests 0–18 changed no behaviour; they recorded
what the system already did. Everything below **modifies the system**, and is therefore held to
the change protocol in `plans/REF_GROUND_RULES.md` §2: snapshot, one change, snapshot, compare,
accept only on a paired CI that excludes zero.

**Baseline frozen 2026-07-20** (`scorecard.py snapshot --tag baseline`). Every number below is
measured against it.

---

## STAGE 1 — Cross-donor calibration (Change A)  ⏳ IMPLEMENTED, NOT YET RUN

### Hypothesis

**One architectural mistake, four manifestations.** Every calibration parameter in the bundle is
fitted on data from donors the model trained alongside, then applied to a held-out donor whose
error regime is completely different:

| Parameter | Fitted on | Symptom out-of-donor | Test |
|---|---|---|---|
| `temperature` | `val` split | fate ECE **0.281** | T8.2 |
| `conformal.q` | `calib` residuals | coverage **0.401 vs 0.90**; **0.000** on N2/N3 | T14 |
| `sigma_age` | ensemble spread, never calibrated | **~2.4 yr** vs true error **~14 yr** | T7.4.2 |
| `ood` reference | `train` trunk features | AUC **0.47** ≈ chance | T15 |

If the hypothesis is right, refitting the first three on **cross-donor** statistics (inner
leave-one-donor-out within the training donors) should move all three toward honesty **without
touching discrimination or point predictions** — because nothing about the model changes, only
what the calibrators were shown.

### Prediction, made BEFORE running

| Metric | Baseline | Predicted after | Basis |
|---|---|---|---|
| `conformal_coverage` | 0.401 | **0.85–0.95** | arithmetic on T14 |
| `conformal_q` | 8.86 | **~30–40 yr** | P90 ≈ 2.5–3.0× mean error (14.29) |
| `conformal_width` (= 2q) | 17.72 | **~70–86 yr** | 2× the above |
| `fate_ece` | 0.281 | **~0.13, bar ≲0.17** | measured, T8.2 (Platt) |
| `sigma_scale` | — (1.0) | **~5–6** | 14 yr error / 2.4 yr spread |
| `ood_rate` | 0.273 | **unchanged** | OOD refit deliberately not implemented |
| `res_approvals` | 3 (oracle 0) | **0** | honest σ ⇒ `R_eff` = 0, MASTER_PLAN §5b-ter |
| `dage_mae_model` | 14.291 | **unchanged** | same weights |
| `rank_model_dage` | 0.948 | **unchanged** | same weights |
| `fate_prauc` | 0.992 | **unchanged** | temperature preserves argmax |

**Sub-stage 1a is predicted to be BIT-IDENTICAL, not merely "noise."** It appends a column that
the network never sees (`forward` takes x, u, dose_time only); indices 0–5 are unchanged, and
adding a tensor consumes no RNG, so shuffle order and weights are identical.

### What each outcome lets us conclude — branches fixed before seeing data

| Outcome | Conclusion | Next |
|---|---|---|
| Coverage 0.85–0.95 **and** ECE drops ≥40%, guards all `noise` | **Hypothesis CONFIRMED.** In-distribution calibration was the root cause. The generalizable methodological finding stands | Stage 2 decision (reference cells?) |
| Coverage overshoots >0.95 | `q` inflated by the N2/N3 outliers — **expected, record it.** Do NOT tune `q` down; that is fitting the test | **FAIL** — a failed target, not a qualified success; then a new pre-registered bar correcting the CV bias (see ruling §1 below) |
| Coverage still ≈0.40 | cross-donor residuals are not representative — most likely too few inner donors | check `xdonor_n_donors == 5` before anything else |
| ECE unchanged | logits pooled from too few donors, or the defect is not calibration | diagnose separately; the three refits are independent |
| **Any of `fate_prauc`, `fate_roc`, `rank_model_dage`, `dage_mae_model` REGRESSES** | the change reached something it must not | **REVERT.** A bug, not a trade-off |
| 1a moves any metric at all | the column edit did more than add a column | **REVERT and re-audit** |

**`rank_res` and `res_approvals` are NOT guards.** RES is expected to move — see below.

### Two ambiguities in the plan, resolved BEFORE the run (user ruling, 2026-07-20)

The plan contradicts itself on one bar and is silent on a near-miss. Both were decided in advance,
because deciding after seeing the numbers is how every change comes to look like an improvement
(ground rules §5).

**1. Coverage above 0.95 → FAIL.** §3 sets the bar at 0.85–0.95; §1b.4 calls overshoot "expected"
and forbids tuning `q` down. Those cannot both govern. **Ruling: the bar is the bar.** Coverage of
0.97 is recorded as a failed target, not a qualified success.

> **This is likely to bite, and the reason is structural.** `q` is fitted on residuals from inner
> models trained on **4** donors, then applied to a deployed model trained on **5**. Fewer training
> donors ⇒ weaker generalization ⇒ larger residuals ⇒ `q` biased high. This is the standard
> pessimistic bias of cross-validation (CV estimates the error of a model trained on *n(k−1)/k*
> samples, not *n*), compounded by N2/N3 inflating the P90. **Predicted before running: overshoot
> is more likely than landing inside the window.**
>
> If it fails this way, the honest response is a **new test with a new pre-registered bar** —
> correcting the donor-count bias explicitly, which is a principled fix. Shrinking `q` until
> coverage lands in the window is **fitting the test** and is forbidden by §1b.4.

**2. `fate_ece` between 0.17 and 0.22 → FAIL, then fix.** A real improvement that misses the ≥40%
bar is recorded as a **failed target**, not accepted retroactively. The fix is then a *separate
change with its own snapshot and its own bar* — most likely a Platt calibrator, which already
demonstrated **0.153** on this data (baseline `fate_ece_platt`), versus a single temperature
scalar which is strictly less flexible.

**§3 states the three refits are independent**, so a `fate_ece` failure does not invalidate the
coverage and `sigma_scale` results; they are adopted or rejected on their own evidence.

### Sharper prediction on the guards than the plan requires

The plan asks for `noise`. The deployed ensemble trains **before** `crossdonor_stats`, with the
same seeds and data, and `set_global_seed` enables `cudnn.deterministic` and
`use_deterministic_algorithms`. So the deployed weights should reproduce **bit-for-bit**:

| Guard | Predicted |
|---|---|
| `dage_mae_model` | **exactly 14.291** |
| `rank_model_dage` | **exactly 0.948** |
| `level_shift_model` | **exactly** the baseline per-donor values |
| `ood_rate` | **exactly 0.273** |

**Any movement at all in these — even in the third decimal — means the change reached something it
must not.** That is a strictly harder test than "the CI includes zero", and it should be applied.

### The RES prediction, stated in advance so it cannot be read as a regression

`sigma_scale` widens `sigma_age` ~5×, and `R_eff = max(0, −(mu + z·σ))` consumes σ. Per-cell RES
should therefore approve **nothing**, and that is the **correct** result: honest per-cell
uncertainty (~19 yr) exceeds the real effect (~11 yr). This is `MASTER_PLAN` §5b-ter's central
arithmetic playing out, and it is a finding — per-cell confident rejuvenation is unreachable at
this data scale. The RES verdict itself stays deferred to Change C (Stage 4).

### What was built

- **1a** — donor column (7th tensor) sourced from the shard's `cell_line`; two positional unpacks
  in `training/train.py` converted to indexed access.
- **1b** — `training/xdonor_calib.py`: inner-LODO over training donors, pooling residuals, logits
  and ensemble spread; `temperature`, `q` and `sigma_scale` fitted on those. `sigma_scale`
  persisted in `ConformalParams` (defaulted, so pre-existing bundles still load) and applied in
  `Predictor`.
- **Rollback** — `TrainConfig.xdonor_calibration = False` restores every in-distribution path.

**Four defects were found in the plan and fixed; one plan step was found to be unimplementable and
deliberately skipped.** All five are recorded in `plans/STAGE_1_DEVIATIONS.md`. The two that
affect interpretation:

1. **The inner-LODO as specified leaked** — it passed the held-out donor as the early-stopping
   monitor, which would have made the residuals best-case and understated `q`. Fixed.
2. **The OOD refit is not implementable** — the plan pools trunk features across independently
   seeded inner models, whose latent bases differ by arbitrary rotation, while `OODDetector`
   compares the *deployed* model's features. Pooling them makes the Mahalanobis distance
   meaningless. Left unchanged; `ood_rate` is predicted not to move, and the gate's fate is
   deferred to the Stage 3d disable decision (`STAGE_1` §1b.4 anticipates exactly this).

### Load-bearing precondition, unverified at time of writing

**Does `cell_line` carry donor identity, at donor granularity?** `python verify_1a.py` answers it.
Expect **exactly 5** distinct donors in a LOOCV training split.

- fewer → smaller inner-LODO pool, noisier fit
- **more → `cell_line` is finer than donor** (donor × timepoint, say). Holding out such a group is
  *not* holding out a donor; the same donor's cells stay in training, residuals understate true
  cross-donor error, and `q` comes out too small. **This failure looks like success.** Stop and
  inspect the values before proceeding.

### RESULT (actual) — PENDING

> Nothing has been executed. The implementation was written on a machine with no Python, no
> dataset shards and no venv; not even an import check has run. Fill this in from the first real
> run, verbatim, whatever it says.

```
python verify_1a.py                    ->  [paste]
python -m pytest tests/ -q             ->  [paste]
python scorecard.py snapshot --tag A_xdonor
python scorecard.py compare baseline A_xdonor   ->  [paste]
```

| Metric | Baseline | Predicted | Actual | Verdict |
|---|---|---|---|---|
| `conformal_coverage` | 0.401 | 0.85–0.95 | | |
| `fate_ece` | 0.281 | ≲0.17 | | |
| `conformal_width` | 17.72 | ~70–86 | | |
| `sigma_scale` | 1.0 | ~5–6 | | |
| `ood_rate` | 0.273 | unchanged | | |
| `dage_mae_model` | 14.291 | noise | | |
| `rank_model_dage` | 0.948 | noise | | |
| `fate_prauc` | 0.992 | noise | | |

### RESULT — RUN 1 (2026-07-21) — ⛔ **INVALID. The experiment did not test the hypothesis.**

Ran on D: (RTX 2050, CUDA), 6 folds, 212 min. **The run is void for a reason found in its own
diagnostics: `cell_line` is not donor.**

```
DONOR VOCAB: {'HFF': 0, 'N3': 1, 'O1': 2, 'O2': 3, 'Y1': 4, 'Y2': 5, 'N2': 6}
cells per donor (N2 fold train): HFF=33613, N3=14, O1=16, O2=18, Y1=13, Y2=14
```

The training split is the **GSE242423 HFF reprogramming corpus (33,613 cells) merged with the six
Gill donors (~14 cells each)**. `cell_line` labels both, so the inner-LODO rotated over HFF as if
it were a seventh donor:

```
- xdonor.fold  donor=0 (HFF)  n_train=75      n_held=33613
- xdonor.fold  donor=1 (N3)   n_train=33674   n_held=14
  ... (four more, ~14 cells each)
- xdonor.done  n_donors=6  n_residuals=33688
```

**The HFF fold trained on 75 cells** — val_loss **33.0** against the deployed model's **5.3** — and
because it is also the largest fold it contributed **33,613 of 33,688 pooled residuals (99.8%)**.

> **So `q` and `sigma_scale` were calibrated against data starvation, not donor shift.** The
> quantity Stage 1 set out to measure was never measured.

**The smoking gun is `sigma_scale`, which should be similar across folds:**

| fold | N2 | N3 | O1 | O2 | Y1 | **Y2** |
|---|---|---|---|---|---|---|
| `sigma_scale` | 6.28 | 22.56 | 11.40 | 16.40 | 11.79 | **74.45** |
| `q` (yr) | 19.73 | 40.23 | 37.53 | 41.33 | 34.65 | 36.87 |

A **12× spread**. Y2's 74.45 implies a median ensemble spread of **0.50 yr** against a P90 residual
of 36.9 — five members trained on 75 cells agreeing closely with each other while collectively
catastrophically wrong. That is the ensemble-under-shift failure in its purest form, but it is
*our own construction*, not a property of the data.

### Scorecard, run 1 — recorded for completeness, interprets nothing

| Role | Metric | Baseline | Run 1 | Verdict vs bar |
|---|---|---|---|---|
| TARGET | `conformal_coverage` | 0.401 | **0.873** | ACCEPT, CI [+0.171,+0.774] — **but see below** |
| TARGET | `fate_ece` | 0.281 | **0.227** | **FAIL** — `noise`, 19% drop vs the ≥40% bar |
| GUARD | `dage_mae_model` | 14.291 | 14.291 | ✅ **bit-identical**, +0.000 every fold |
| GUARD | `rank_model_dage` | 0.948 | 0.948 | ✅ **bit-identical** |
| GUARD | `fate_prauc` | 0.992 | 0.988 | ✅ noise |
| GUARD | `fate_roc` | 0.983 | 0.978 | ✅ noise |
| watch | `conformal_width` | 17.72 | **70.12** | as predicted (70–86) |
| watch | `ood_rate` | 0.273 | 0.273 | ✅ unchanged, as predicted |
| watch | `res_approvals` | 3 | **0** | over-approval gap → 0, as predicted |

**The coverage "ACCEPT" is an averaging artifact and must not be reported as a pass.** Per fold:

```
N2 0.381 | N3 0.857 | O1 1.000 | O2 1.000 | Y1 1.000 | Y2 1.000
```

Four folds cover **everything** (q far exceeds their error) and N2 covers **38%** (q=19.7 below its
21.8 MAE). The 0.873 mean is two opposite failures cancelling. Under the pre-registered ruling
(>0.95 → FAIL), **four of six folds fail**.

### Two things the run DID establish

1. **The guards behaved exactly as predicted**, including the sharper bit-identical prediction:
   `dage_mae_model` and `rank_model_dage` moved by **+0.000 on every fold**. Stage 1 provably does
   not touch the model — only the calibration layer.
2. **`fate_prauc`/`fate_roc` moved slightly (0.992→0.988) and this is CORRECT, not a leak.** `S` is
   `softmax(logits/T)[:,0]`, and for 3-class softmax the ordering of one class's probability across
   cells is **not** temperature-invariant (the normaliser depends on all three logits). Temperature
   legitimately changed, so a small rank change follows. Verified analytically with a counterexample.

### Root cause of the invalidity — a defect in the verification, not just the calibration

`verify_1a.py` **detected this and printed the warning verbatim** —

> `!! MORE than the expected 5; saw 6. THIS IS THE DANGEROUS DIRECTION.`

— and then **graded the run `PASS`**, because the verdict logic only escalated to `STOP` on
*too few* donors. The operator followed a PASS. **This cost 3.5 hours of GPU time and a void
experiment.** The check existed, fired, and was ignored by its own scoring rule.

### Fixes applied (not yet run)

| Where | Fix |
|---|---|
| `xdonor_calib.py` | `MIN_INNER_TRAIN_FRAC = 0.5` — skip any inner fold whose held-out donor leaves <50% of the training split. Such a fold measures data starvation, not donor shift. Raises if fewer than 2 usable folds survive |
| `verify_1a.py` | `STOP` (not PASS-with-warning) when a donor holds >50% of a training split, **or** when the donor count is anything other than the expected 5 |
| `tests/test_training.py` | two tests: a 90%-dominant donor must be skipped and must not reach the residual pool; a 95/5 split must raise rather than calibrate off one donor |

**Bars are unchanged.** This is ground rule §6 — *"when a result surprises you, the default
assumption is a bug in the test"* — not a retroactive threshold move. Run 2 will pool ~75 honest
Gill-donor residuals (5 donors × ~15 cells) instead of 33,613 HFF ones.

### Predictions for RUN 2, recorded before it runs

| Metric | Run 1 | Run 2 predicted | Why |
|---|---|---|---|
| `xdonor_n_donors` | 6 | **5** | HFF skipped |
| `xdonor_n_residuals` | 33,688 | **~75** | Gill donors only |
| `sigma_scale` spread | 6.3–74.5 | **narrower, ~3–8** | no starved fold |
| `conformal_coverage` | 0.873 (0.38–1.00) | **less saturated** | q from honest residuals |
| guards | bit-identical | **bit-identical** | unchanged |

> **Honest caveat, stated now.** With only ~15 training cells per Gill donor, the inner-LODO pool
> is ~75 residuals from 5 donors. That is enough for a 90% split-conformal quantile
> (`ceil(76×0.9)`=69th of 75), but the per-donor heterogeneity already measured (MAE 5.4 on O1 vs
> 29.7 on N3) means **a single global `q` may be unable to hit 0.85–0.95 on every fold regardless
> of how correctly it is fitted.** If run 2 still shows saturated folds, that is a finding about
> donor heterogeneity — not a further bug — and the response is a per-donor or conditional
> interval, pre-registered separately.

### That caveat is now QUANTIFIED, before run 2 — `experiments/q_power_analysis.py`

The caveat above was a hunch. It has been simulated from the baseline's own per-fold MAEs, on
the pool geometry run 2 will actually have (5 donors × ~14 cells ≈ 70 residuals). No GPU, no
bundles, reproducible.

**Result 1 — `q` is very noisy at this pool size.** Over 4,000 resamples: median **36.2 yr**,
90% range **[23.7, 48.2]** — a **68% spread from sampling alone**, before any modelling error.

**Result 2 — a single global `q` cannot serve these donors.** Coverage that the median `q`
delivers to each held-out donor:

| fold | MAE | coverage | |
|---|---|---|---|
| O1 | 5.39 | **1.000** | saturated |
| Y1 | 7.28 | **1.000** | saturated |
| O2 | 7.54 | **1.000** | saturated |
| Y2 | 14.06 | **0.960** | saturated |
| N2 | 21.79 | **0.808** | under-covers |
| N3 | 29.69 | **0.666** | under-covers |

**Aggregate mean 0.906 — inside the 0.85–0.95 window — while `0/6` individual folds are.**

> **This is the single most important thing to know before reading run 2.** Donor error scales
> differ **5.5×**, so one scalar `q` mathematically cannot hold every fold at 90%. The aggregate
> can land in the window purely by averaging saturation against under-coverage — which is
> **exactly what run 1 showed** (mean 0.873 from 0.381 / 0.857 / 1.000 / 1.000 / 1.000 / 1.000).
>
> The model is conservative: half-normal residuals give `P90/mean = 2.07`, while `MASTER_PLAN`
> §5d puts our heavy-tailed mixture at **2.67**. The real `q` is therefore *larger*, and
> saturation *more* likely, than simulated.

**Consequence for the pre-registered bar.** A run-2 miss on `conformal_coverage` is now expected
on structural grounds, and would be **a finding about donor heterogeneity, not evidence that the
calibration code is wrong**. The two must not be confused. The legitimate responses are:

1. accept the failure and report it — a single global interval is not appropriate for donors
   whose errors differ 5.5×; or
2. run a **new** test with a **new** pre-registered bar, for a per-donor or conditional interval.

**Never** shrink `q` until coverage lands in the window. `STAGE_1` §1b.4 already forbids it, and
the reason is now numerical rather than rhetorical: with a 68% sampling spread, *any* target
coverage can be hit by choosing a `q` after seeing the result.

**Also worth watching:** the aggregate `conformal_coverage` is a mean over folds, so it can read
"in range" while nothing is. **Read the per-fold row, not the mean.**

### VERDICT — RUN 1: **INVALID, RE-RUN REQUIRED**

Not a failure of the hypothesis. The hypothesis was never tested: the calibration set was 99.8%
data-starvation residuals. Both targets are void; both guards passed and are informative.

---

### RESULT — RUN 2 (2026-07-22) — ✅ **VALID.** Scored: 1 target passed, 1 REGRESSED.

Ran on D: (RTX 2050), 6 folds, 210 min. The run-1 defect is gone, and it is visible in the logs
rather than asserted:

```
inner-LODO: SKIPPING donor 0 -- holding it out leaves 75 of 33688 cells (0.2%, below the 50% floor)
xdonor.done  n_donors=5  n_skipped=1  n_residuals=103  residuals_per_donor={21,21,21,19,21}
```

| | run 1 | run 2 |
|---|---|---|
| HFF | rotated as a donor, supplied **99.8%** of residuals | **skipped** |
| inner-model val_loss | HFF fold **33.0** vs deployed 5.2 | every fold **5.17–5.27**, matching deployed |
| cells per donor | 13–18 (train only) | **21** (train+val+calib, all unseen by the inner model) |

Every inner model is now a genuine proxy for the one that ships. That is what makes the run valid.

#### Scored against the pre-registered bars (§3)

| role | metric | bar | result | |
|---|---|---|---|---|
| GUARD | `dage_mae_model` | noise | 14.291 → 14.291 | ✅ |
| GUARD | `rank_model_dage` | noise | 0.948 → 0.948 | ✅ |
| GUARD | `fate_prauc` | noise | 0.992 → 0.992 | ✅ |
| GUARD | `fate_roc` | noise | 0.983 → 0.983 | ✅ |
| GUARD | `ood_rate` | unchanged (deviation A2) | 0.273 → 0.273 | ✅ |
| GUARD | `level_shift_model` | noise | −5.713 → −5.713 | ✅ |
| **TARGET** | `conformal_coverage` | 0.85–0.95 | 0.401 → **0.889**, ACCEPT | ✅ |
| **TARGET** | `fate_ece` | ACCEPT + ≥40% drop | 0.281 → **0.364**, CI [+0.058,+0.109] | ❌ **REGRESSION** |

**All six guards came back with `max |diff| = 0.00e+00` on every fold** — not merely "noise" but
bit-identical, which was the sharper prediction recorded above. Stage 1 provably touches
calibration and nothing else.

#### Prediction scorecard (all recorded before the run)

| predicted | actual | |
|---|---|---|
| guards bit-identical | 0.00e+00 on all six | ✅ |
| coverage aggregate lands in range | 0.889 | ✅ |
| coverage per-fold bimodal | 1.000 ×5, N3 0.333 | ✅ |
| `fate_ece` misses the 40% bar | it **regressed** — worse than predicted | ⚠️ |
| `sigma_scale` ≈ 5–6 | **9.9–18.6** | ❌ my arithmetic: I divided by MEAN error; the formula uses P90 |
| temperature > 1 (softening) | **0.755–0.849** (sharpening) | ❌ predicted from the synthetic rehearsal, not from data |

#### Finding 1 — ONE donor sets `q` for the entire study

| fold | own MAE | q | q/MAE | coverage | N3 in pool? |
|---|---|---|---|---|---|
| N2 | 21.79 | 33.76 | 1.55 | 1.000 | yes |
| **N3** | **29.69** | **24.39** | **0.82** | **0.333** | **no** |
| O1 | 5.39 | 34.64 | 6.43 | 1.000 | yes |
| O2 | 7.54 | 36.27 | 4.81 | 1.000 | yes |
| Y1 | 7.28 | 34.41 | 4.73 | 1.000 | yes |
| Y2 | 14.06 | 34.19 | 2.43 | 1.000 | yes |

`q` is 33.8–36.3 on every fold where N3 is in the calibration pool and **24.4** on the one fold
where it is not. The pooled P90 lands inside N3's residuals, so **N3's offset alone determines the
interval everyone else receives** — and LOOCV then removes it from its own pool, dropping `q`
below N3's own error. `q/MAE` spans **0.82 → 6.43**.

#### Finding 2 — the residuals are SHIFTED, not scattered

N2's MAE is 21.79, yet **all 21** of its cells fall inside q = 33.76. Under the half-normal my
power analysis assumed, roughly a third should have exceeded it. So residuals cluster around a
per-donor **offset** rather than spreading from zero — T7.4.3's level shift appearing directly in
the coverage numbers, and making coverage nearly a **step function**: `q` either clears a donor's
offset (→1.000) or it does not (→0.333).

> **The "heavy-tailed mixture" framing in the power analysis above is RETRACTED.** The error shape
> is offset-dominated, not heavy-tailed. That points at **Stage 2**, which compresses per-donor
> MAE from [5.4 … 29.7] to [4.3 … 10.0] (T16) — a 5.5× spread down to 2.3×.

#### On coverage passing

0.889 against a nominal 0.90 is split conformal's **marginal** guarantee working. The per-fold
spread is **conditional** coverage, provably unachievable distribution-free (Barber, Candès,
Ramdas & Tibshirani 2021). The aggregate is a mean over folds and can read "in range" while no
fold is — so **read the per-fold row**. Recorded as a property of the method, not a bug.

#### Finding 3 — why `fate_ece` regressed: four quantities, no two the same

| stage | quantity |
|---|---|
| `calibrate.py:_nll` optimised | **multi-class** NLL |
| `metrics.py:ece` reported | **top-1** confidence ECE |
| `scorecard.py:_ece` grades | **binary** ECE on `P(safe)` |
| `res.py` + `STAGE_3` §0.1 consume | **`S` = `P(safe)`, `P_loss`** |

Plus a fit/apply mismatch: temperature is fitted on `ensemble_logits` (mean of member logits) but
applied per-member and then averaged — `softmax(mean(lg)/T) ≠ mean(softmax(lg/T))` by Jensen.

Cross-donor multi-class NLL chose T ≈ 0.755–0.849 while the held-out folds want ≤0.54; every fold
worsened by +0.053…+0.108. The cross-donor top-1 ECE *improved* at the same time (0.269→0.217) —
not a contradiction, a different metric.

### VERDICT — RUN 2: **PARTIAL PASS.** Stage 1 does not pass as a whole.

§3: *"Accept only if the TARGET metric says ACCEPT and no GUARD says REGRESSION."* A target
regressed. §3 also makes the refits independent, so:

| refit | verdict |
|---|---|
| conformal `q` | **ADOPTED** — 0.401 → 0.889, with the per-fold caveat on record |
| `sigma_scale` | **ADOPTED** — drove RES approvals 3 → 0 and the over-approval gap to zero, as predicted |
| temperature | **REJECTED** — regression; the in-distribution fit was better on the graded metric |
| OOD | not attempted (deviation A2 — not implementable as specified) |

---

## Change A″ — calibrate `P(safe)`, the quantity the product ships  ⏳ IMPLEMENTED, NOT YET RUN

### Hypothesis

The fate calibrator optimises a quantity nothing uses. `res.py` consumes `S` and `P_loss`,
`STAGE_3` §0.1 needs a risk threshold on `P(unsafe)`, and the scorecard grades binary ECE on
`P(safe)` — while `fit_temperature` minimises multi-class NLL. Calibrating what is actually
shipped and graded should recover the bar.

**This is what the plan's own bar already assumes.** T8.2's table is, cell for cell, the
scorecard's `fate_ece` and `fate_ece_platt` columns — and its "ECE recal" is **Platt fitted on
the calib split**. `MASTER_PLAN` §5a names the defective quantity as "`S`, `P_loss`" with
"**YES — Platt halves it**". So `STAGE_1`'s ≲0.17 bar was derived from an in-distribution-fitted
Platt; §1b.2's `fit_temperature(xstats…)` is the line that never matched §2's own expected effect.

### Is this a departure from the cross-donor principle? **No.**

The principle says *calibrate on data whose error regime matches deployment*. Measured:

| quantity | in-distribution | out-of-donor | premise met? |
|---|---|---|---|
| ΔAge error | ~4 yr | ~14 yr | **yes** → `q`/`sigma_scale` use the pool alone |
| fate discrimination | 0.929–0.940 | **0.96–1.00** (T8.1) | **no** — no degradation |
| fate calibration | — | calib-fitted Platt halves out-of-donor ECE on 4/5 folds (T8.2) | **no** — it transfers |

So the in-distribution split *qualifies* for fate, and there is 43× more of it. The calibrator is
fitted on the **union** (~4,593 cells). Fitting on the 103-cell pool alone would discard 97.8% of
the available data for a 2-parameter fit — and because cells within a donor share that donor's
offset, the pool's **effective** n is nearer 5 than 103.

**The principle is tested, not assumed:** the strict pool-only Platt is fitted on every run and
reported (`xdonor_only_platt_a/b`, `xdonor_only_safe_ece_insample`, `shipped_safe_ece_on_pool`),
never shipped.

### Prediction for RUN 3, recorded before it runs

| metric | run 2 | predicted | why |
|---|---|---|---|
| `fate_ece` | 0.364 | **≈0.15–0.17** | in-dist Platt measures 0.153 (T8.2); the union adds the pool |
| `conformal_coverage` | 0.889 | **≈0.889, unchanged** | the calibrator does not enter `q` |
| `fate_prauc` / `fate_roc` | 0.992 / 0.983 | **unchanged** | positive slope ⇒ monotone, so no reordering. Audit found two ways it could still MERGE cells (EPS clamp, float32 output cast) and both are fixed |
| `dage_mae_model`, `rank_model_dage` | unchanged | **bit-identical** | same weights |

**Bar unchanged: ACCEPT + ≥40% drop (≤0.169).** Not weakened because the specification was wrong.
A corrected calibrator that still misses is a real result about calibrator capacity at n≈103.

> **If `fate_ece` lands just above 0.169** — accept the failure, or pre-register a new bar for a
> different calibrator. Do **not** re-read the number and rationalise. That rule is what made
> run 2's regression informative instead of embarrassing.

### RESULT — RUN 3 (2026-07-23): **PARTIAL** — coverage passes, `fate_ece` accepts but misses the drop

Executed on the data machine. 222 tests pass; 6/6 folds retrained in **229.0 min**; snapshot
`B_fatecal`; compared against both `baseline` and `A_xdonor`.

**Scored against the bars in `STAGE_1_CALIBRATION.md` §3 (lines 402–407):**

| Role | Metric | Bar | baseline → B_fatecal | Verdict |
|---|---|---|---|---|
| TARGET | `conformal_coverage` | reach 0.85–0.95 | 0.401 → **0.889** ACCEPT | ✅ **PASS** |
| TARGET | `fate_ece` | ACCEPT + ≥40% drop (≤0.169) | 0.281 → **0.249** ACCEPT, **−11.0%** | ❌ **MISS** |
| GUARD | `fate_prauc` | noise | 0.992 → 0.992 (+0.000) | ✅ |
| GUARD | `fate_roc` | noise | 0.983 → 0.983 (+0.000) | ✅ |
| GUARD | `rank_model_dage` | noise | 0.948 → 0.948 (+0.000) | ✅ |
| GUARD | `dage_mae_model` | noise | 14.291 → 14.291 (+0.000) | ✅ |

**All four guards bit-identical (+0.000, CI [+0.000,+0.000]) for the third consecutive run.**
Stage 1 provably does not touch the model. `interval_width` 17.717 → 65.888 reads REGRESSION but
is **not a guard** and widening was the pre-registered consequence of an honest `q`.

`A_xdonor → B_fatecal`: `fate_ece` 0.364 → 0.249 ACCEPT, `fate_ece_platt` 0.161 → 0.140 ACCEPT,
coverage unchanged at 0.889. The union fit **did** repair run 2's regression — it just did not
repair enough.

#### The prediction was falsified, and the reason is precise

Predicted `fate_ece` **≈0.15–0.17**; measured **0.249**. The stated reasoning was *"in-dist Platt
measures 0.153 (T8.2)"*. **That equated two different quantities.**

`scorecard.py:189` computes `fate_ece_platt = _ece(_platt(S_cal, …, S), st)`, and `S` at
`scorecard.py:157` comes from `est.rows(...)` → `predictor._raw_batch`, which has **already
applied the bundle's calibration** (`predictor.py:170`). So `fate_ece_platt` is not a standalone
in-distribution Platt — it is a **SECOND calibration layer stacked on top of whatever the bundle
ships**, fitted on the calib split.

That single fact explains the whole table:

| snapshot | bundle ships | `fate_ece` (bundle alone) | `fate_ece_platt` (+ stacked layer) |
|---|---|---|---|
| `baseline` | temperature | 0.281 | 0.153 |
| `A_xdonor` | cross-donor temperature | 0.364 | 0.161 |
| `B_fatecal` | union Platt | **0.249** | **0.140** |

**The stacked layer lands at 0.140–0.161 regardless of what the bundle does.** It, not the
bundle's calibrator, was doing the work in every T8.2 number. The 0.153 was never available to a
single-layer bundle calibrator, so a bar derived from it was never within reach of the thing being
built. This is the same class of error as the earlier T8.2 retraction, one level deeper.

#### The bar itself is fair — checked before blaming it

`fate_ece` is measured on 19–21 held-out cells with 10 bins (~2 cells/bin), so the estimator is
biased upward and the bar could have been below its resolution. Simulated a **perfectly
calibrated** model (`y ~ Bernoulli(p)`, so any ECE is pure estimator bias) at run-3's geometry:

| regime | n | median ECE | P(≥0.17) |
|---|---|---|---|
| confident (matches PR-AUC 0.992) | 21 | 0.075 | 1.7% |
| confident, **mean of 5 folds as scored** | 5×21 | **0.078**, 90% range [0.057, 0.105] | **0.0%** |

The floor is **0.078**; the bar 0.169 sits at ~2× it. **The bar is attainable and stands.**
0.249 is a real miss, not a measurement artefact.

#### Why the union fit under-delivered

`fate.calibrated` reports `total=4509 in_dist=4406 xdonor=103` — the cross-donor pool is
**2.28%** of the fitting data. The fitted slopes show it was drowned:

| | slope `a` across the 6 folds | mean |
|---|---|---|
| shipped (union) | 2.574, 2.599, 2.584, 2.621, 2.595, 2.622 | **2.599** |
| pool-only (diagnostic, never shipped) | 1.533, 1.238, 1.404, 1.417, 1.542, 1.144 | **1.380** |

The union slope is ~1.9× the pool-only slope and is **tight to ±0.024 across folds** — the
signature of a fit determined by the 4406 shared in-distribution rows, not by the 103 rows that
differ per fold. The union is, to three digits, the in-distribution fit. This is the deviation
from *"calibrate on data whose error regime matches deployment"* that was flagged when it was
made, and it cost the target.

#### Recorded, not yet explained

A synthetic probe of the two calibrator families (`sigmoid(a·logit p + b)` as shipped, vs
`LogisticRegression` on raw `p` as `scorecard._platt` uses) **failed to reproduce** the observed
gap — it made logistic-on-`p` *worse*, not better. The leading hypothesis is that logistic-on-`p`
is **bounded** (it cannot exceed `sigmoid(w+c)`, the calib empirical rate at high `p`) whereas
logit-Platt with `a>1` drives saturated inputs to exactly 1.0, so only the former can pull the top
ECE bin down. **The mechanism is not confirmed** and is not relied on below.

#### Next step is a measurement, not a change

`train_model.py:240-243` already computed `xdonor_only_safe_ece_insample` and
`shipped_safe_ece_on_pool` per fold and wrote them to `bundle/metrics.json`; the console printed
only the slopes. Those numbers score the pool-only calibrator against the shipped one **on the
same pool, at zero compute cost**, and no further change is pre-registered until they are read.

---

## RUN 3 POST-MORTEM (2026-07-23) — the diagnostics came back. Three of my claims are wrong.

Offline analysis of `diag_dump/` (6 folds; pool, calib and test arrays, raw and calibrated).
The pipeline reproduces the graded `fate_ece` from raw probabilities to **0.00e+00**, so every
number below is the same quantity `scorecard.py` grades.

### RETRACTION 1 — "the bar is fair and attainable". It is not: it sits BELOW the estimator floor.

I computed the floor from `beta(12,1)` (most p > 0.9), inferring saturation from fate PR-AUC
0.992, and reported a floor of 0.078 with the bar at ~2× it. **The dump shows no saturation
whatever**: on test, P(safe) spans 0.09–0.88 with **0.0% above 0.99**; on the pool, max is 0.65.
Near-perfect *ranking* does not imply saturated *probabilities* — that inference was simply wrong.

Recomputing the floor from the **actual** test P(safe) vectors (`y ~ Bernoulli(p)`, so any ECE is
pure estimator bias):

| | floor median | floor p90 | observed | percentile |
|---|---|---|---|---|
| mean over the 5 scored folds | **0.183** | 0.225 | 0.249 | 99.5% |

**A perfectly calibrated model clears the 0.169 bar only 26.9% of the time at this geometry.**
The bar is below what n≈21 with 10 bins can resolve. The miss stands as a pre-registered fact,
but it carries almost no information about calibration quality.

The observed 0.249 still sits at the **99.5th percentile** of the perfect-calibration
distribution, so there *is* real miscalibration — just not 0.249 worth of it.

### RETRACTION 2 — "the union fit cost the target". It did not: it is the best candidate tried.

Scoring all three pre-registered candidates on the graded folds, each against the floor **of its
own output** (raw ECE is not comparable across candidates — see Retraction 3):

| candidate | ECE | own floor | **excess** | vs bar |
|---|---|---|---|---|
| identity (no calibration) | 0.364 | 0.172 | **+0.192** | miss |
| **union Platt [SHIPPED]** | 0.249 | 0.178 | **+0.071** | miss |
| pool-only Platt [the PRINCIPLE] | 0.308 | 0.164 | **+0.144** | miss |

**Reverting to the cross-donor principle would have been twice as bad.** The 2.28%-drowning
argument correctly described the *fit*, but the conclusion drawn from it — that the deviation
cost the target — is refuted. Recorded rather than deleted.

On excess, the shipped calibrator removes **63%** of the miscalibration present with no
calibration at all (+0.192 → +0.071). The effect Stage 1 was built to produce is there; the
metric it was graded on cannot show it.

### RETRACTION 3 — raw ECE rewards SHARPENING, so it cannot rank calibrators

Refitting Platt on four donors' held-out cells and applying to the fifth appeared to take ECE
0.249 → 0.103, apparently *beating* the 0.179 floor — impossible for an honest estimate. The
refit sharpens (a = 3.4–5.7), and sharper probabilities occupy more extreme bins where Bernoulli
variance is smaller, which **lowers the floor**. Against its own floor the refit scores +0.035,
not −0.076.

**Sharpening alone accounts for 0.110 of the 0.146 apparent gain — 75%.** A calibrator can move
`fate_ece` toward the bar by sharpening while getting no better calibrated. This is a defect in
the target metric, not a strategy: it is the one way this bar could be "landed" dishonestly, and
it is recorded here so that route is closed.

### What the pool said, and why it was misleading

| family (LODO **within** the pool) | ECE |
|---|---|
| identity — no calibration at all | **0.130 — already "passes" 0.169** |
| logit-Platt [shipped] | 0.067 |
| logistic-on-p [T8.2 family] | 0.101 |
| isotonic | 0.059 |

On the pool everything passes, including doing nothing. **The pool is not a proxy for the graded
donor.** Base rates: calib **0.514**, pool **0.64**, test **0.754** (N2 1.000, Y1 0.579). The
deployed calibrator is fitted for a 0.51-safe world and graded on a 0.75-safe one — the residual
is *label shift*, which no calibrator fitted on source data can correct.

Answers to the four questions: **Q1** family is fine, logit-Platt beats logistic-on-p 0.067 vs
0.101; **Q2** union beats pool-only, +0.071 vs +0.144; **Q3** capacity is ample on the pool
(isotonic 0.059) but no pre-registered candidate clears the bar on test; **Q4** ICC **0.342**,
n_eff **13.7 of 103** — donor offsets dominate, so the pool is ~14 independent points.

### Where the residual actually lives

Per-fold excess for the other-donor refit: N3 +0.044, O1 −0.023, O2 −0.014, Y2 0.000, and
**Y1 +0.165**. Y1's base rate is 0.579 against 0.76–0.86 for the rest. Donor calibration
transfers between donors that resemble each other and fails on the one that does not — the same
donor heterogeneity that leaves N3's conformal coverage at 0.333.

**That is Stage 2's subject, not Stage 1's.** No further calibrator change is pre-registered:
the family is right, the fitting set is right, and the remaining error is not a calibration
problem.

### INSTRUMENT AUDIT (2026-07-23) — "probably can't resolve it" was checked, not assumed

Challenged on asserting the bar was unresolvable rather than demonstrating it. Measured
(`audit_metrics.py`): simulate a system that satisfies the intent EXACTLY and ask how often the
criterion reports it as passing.

| criterion | null median | **pass rate for a system that IS correct** |
|---|---|---|
| `fate_ece ≤ 0.169`, **as graded** (mean of per-fold, n≈21 × 5) | 0.183 | **26.9%** |
| `fate_ece ≤ 0.169`, **pooled** (n = 103) | 0.090 | **99.6%** |
| `conformal_coverage ∈ [0.85, 0.95]`, pooled marginal (n = 124) | 0.903 | **93.0%** |

**Confirmed: as graded, a perfectly calibrated model fails the bar 73% of the time.** The bar
would have to be ≥ **0.225** for a perfect model to pass 95% of the time. The criterion is
measuring the sample size, not the model.

**And the fix works.** Pooling held-out cells across folds instead of averaging per-fold ECEs
takes the pass rate for a correct system from 26.9% to **99.6%**. Pooling is the more correct
LOOCV estimate — every cell is still predicted by a model that never saw it. Under it,
run 3 scores **0.211 against the 0.169 bar: still a MISS**, now at the 100th percentile of the
null, i.e. unambiguously real. Repairing the instrument does not hand Stage 1 a pass; it converts
"cannot tell" into "genuinely short, by a measurable amount".

**`conformal_coverage` is sound as written** — 93.0% pass rate for a correctly-90% system, and
the pooled marginal rate is exactly what a conformal guarantee promises. Stage 1's coverage PASS
survives audit.

#### Correction made during the audit

The first version of the guard analysis reported a minimum detectable effect from each metric's
**own fold-to-fold SD** — e.g. `dage_mae_model` SD 9.67, "MDE 10.15 years", which would have
meant the guard was nearly blind. **That was wrong.** The paired test is built on the
DIFFERENCES, where the metric's own spread cancels. The sensitivity is

> minimum detectable **mean** effect = k × SD(**effect** across folds),  k = 1.05 (6 folds), 1.24 (5)

so a **uniform** change is caught at **any** magnitude — which is precisely why Stage 1's guards
reading +0.000 with CI [+0.000, +0.000] is strong evidence rather than luck. Calibrated against
the one real non-zero change measured so far (A_xdonor → B_fatecal on `fate_ece`: mean 0.115,
CI half-width 0.0275 ⇒ SD(effect) 0.0221), that change sat at **4.2×** the detection threshold.

**The real blind spot** is a change that helps some folds and hurts others: it can be large in
the mean and still read as noise. Stage 1 could not hit it — the model was untouched — but every
change from Stage 2 on will be heterogeneous across donors by construction. Guard verdicts must
be read with the per-fold column, not the aggregate alone.

#### Not reachable from this dump

`dage_mae_ridge`, `level_shift_ridge`, `rank_ridge_dage`, `rank_res`, `res_*` need ridge
predictions and RES recomputed. **None is a TARGET or a GUARD**, so no acceptance criterion
depends on them; auditing them needs a dump extension, not a rerun.

### STAGE 1 RE-SCORED (2026-07-23) — repaired estimator, unchanged verdict

`scorecard.py` now reports the calibration target in its resolvable form: pooled over all
held-out cells, with the estimator floor and the excess above it. Verified against `diag_dump/`
to reproduce the graded metric to **0.00e+00** before being trusted.

| | per-fold **[as graded]** | pooled **[repaired]** |
|---|---|---|
| `fate_ece` | 0.249 | **0.211** |
| floor (perfect model) | 0.179 | **0.091** |
| excess | +0.071 | **+0.121** |
| a *correct* system passes 0.169 | **26.9%** | **99.6%** |
| verdict | MISS — uninterpretable | **MISS — real, 100th pctile of the null** |

**The verdict is unchanged, and that is what makes the repair legitimate.** Had it converted a
miss into a pass it would have been goalpost-moving; instead it converts "cannot tell" into
"genuinely short, by a measurable amount".

Note the direction: pooled ECE is *lower* (0.211 vs 0.249) but pooled **excess is higher**
(+0.121 vs +0.071). Per-fold averaging lets each donor's offset be absorbed into its own fold's
bins; pooling forces donors with different base rates into shared bins, where the label shift
shows up as what it is. The pooled number is the harsher and more honest one.

### FINAL VERDICT — STAGE 1: **PARTIAL**

| Role | Metric | Bar | Result | |
|---|---|---|---|---|
| TARGET | `conformal_coverage` | 0.85–0.95 | 0.401 → **0.889** pooled marginal | ✅ **PASS**, audited (93.0% pass rate for a correctly-90% system) |
| TARGET | `fate_ece` | ≤ 0.169 | **0.211** pooled | ❌ **MISS**, real and measurable |
| GUARD ×4 | PR-AUC, ROC, rank, MAE | noise | **+0.000** bit-identical ×3 runs | ✅ |

**What Stage 1 delivered:** conformal coverage went from 0.401 to 0.889 — the headline defect
the stage existed to fix, and it is solved. `q` and both `sigma_scale` factors are calibrated on
cross-donor statistics. The model is provably untouched.

**What it did not:** `fate_ece` does not reach its bar, and the diagnostics say why — the
residual is **donor-level label shift** (base rates: calib 0.514, pool 0.64, test 0.754),
concentrated in **Y1** (0.579 against 0.76–0.86 elsewhere). No calibrator fitted on source data
can correct an unknown target prior. Of the pre-registered candidates the shipped union Platt is
the best (excess +0.071 per-fold vs +0.144 pool-only vs +0.192 uncalibrated).

**Carried into Stage 2, not papered over:** more donors is the only thing that addresses a prior
that varies donor to donor. Stage 2's k≈3 reference cells per donor is exactly that intervention,
so `fate_ece` becomes a Stage 2 acceptance metric rather than an unresolved Stage 1 failure.

**Carried as a measurement warning:** every Stage 2 change touches the model and will be
heterogeneous across donors. A guard's paired CI detects a mean effect only above ≈1.05 × the
effect's own fold-to-fold SD, so a change that helps some donors and hurts others can be large
and still read "noise". **Guard verdicts must be read with the per-fold column from here on.**
