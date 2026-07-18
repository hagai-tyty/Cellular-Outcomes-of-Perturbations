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

**Result (actual).** _[TO FILL — user runs test7_4_res_threshold.py]_

**Verdict.** _[TO FILL]_

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
