# CellFate-Rx — Master Plan (v4)

**What changed from v2, and why.** v2 was reviewed and found to contain real errors, not just
pessimistic tone:

1. It said **"ranking is broken."** That is **false.** Ranking scores **0.948–0.955** (Test 7) —
   one of the strongest results in the project. Only the *RES score used as a ranker* scored
   0.686. v2 collapsed "RES underperforms as a ranker" into "ranking is broken".
2. It **demoted RES prematurely.** Tests 14/15/16 subsequently showed RES fails because of **three
   upstream defects**, all fixable, two with confirmed fixes. **RES has never been evaluated with
   correct inputs.** No verdict is possible until they are fixed.
3. It **buried a positive finding from Test 11**: `u_only` predicts fate at **0.83–1.00**, and on
   Y1 **beats** `x_only` (0.831 vs 0.639). The perturbation channel carries real signal.
4. It repeatedly framed **calibration/engineering defects as fundamental limits.**

**Governing rule for this document:** never write "broken" without naming *which component, on
which metric, and whether the cause is fixable*. Distinguish **capability** (ranking) from
**implementation** (the RES formula).

---

# PART I — WHERE WE ACTUALLY STAND

## 1. Executive summary

1. **Ranking works, and works well.** Within-donor Spearman **0.925–0.983** across all six folds.
   This is a real, defensible capability.
2. **Fate discrimination works out-of-donor** (PR-AUC 0.96–1.00), and the perturbation channel
   alone carries fate signal (Test 11).
3. **ΔAge is at the linear optimum.** It ties ridge because the clock is linear; nothing beats
   ridge here. That is correct behaviour, not failure.
4. **Four components have identified, fixable calibration defects** — per-donor level shift,
   conformal intervals, fate probabilities, OOD gate. Fixes are known for all four, demonstrated
   for two.
5. **RES's verdict is DEFERRED.** It sits downstream of all four defects. Its current failure is
   fully explained by them, and it has never been run with corrected inputs.
6. **The one genuine, code-immune limit** is the data: a single OSKM cocktail means the
   cross-*perturbation* thesis cannot be posed here.

## 2. Component status — four categories, not "broken"

| Component | Status | Evidence |
|---|---|---|
| **Within-donor ranking** | ✅ **WORKS** | Spearman **0.925–0.983**, every fold (T7, T7.4.3) |
| **ΔAge vs baselines** | ✅ **WORKS (at optimum)** | ties ridge (T5); nothing beats ridge (T6); architecture proven capable (T3) |
| **Fate discrimination** | ✅ **WORKS** | PR-AUC 0.929–0.940 in-dist, **0.96–1.00 out-of-donor** (T8.1); holds on hard fold Y1 (0.961 vs logreg 0.636) |
| **Perturbation channel** | ✅ **CARRIES SIGNAL** | `u_only` fate PR-AUC **0.83–1.00**; on Y1 **beats** `x_only` (T11) |
| **Harmonization** | ✅ **WORKS** | parameter-free, leak-safe, unit-tested |
| **ΔAge absolute level** | 🔧 **FIXABLE — fix demonstrated** | ±12.7 yr per-donor shift (T7.4.3); T16: **k=3 reference cells → MAE 14.3→7.1 (−50%)** |
| **Conformal intervals** | 🔧 **FIXABLE — cause identified** | coverage **0.40 vs 0.90** (T14); `q` fitted in-distribution (MAE≈4) applied out-of-donor (MAE≈14) |
| **Fate calibration** | 🔧 **FIXABLE — fix demonstrated** | ECE 0.28; Platt halves it to ~0.13 (T8.2) |
| **OOD gate** | 🔧 **FIXABLE or DISABLE** | AUC **0.47** ≈ chance (T15); zeroes RES for no measurable benefit |
| **RES score** | ⏸️ **VERDICT DEFERRED** | fails on all measured axes — **but every input it consumes is defective.** Never run with corrected inputs. See §4 |
| **Cross-perturbation thesis** | 🚫 **DATA-LIMITED** | one cocktail; time-only variation (T11, T11.1). No code change helps |

## 3. What the ranking result actually says

Because v2 got this wrong, stating it precisely:

| What was ranked | Spearman vs true ΔAge | Reading |
|---|---|---|
| model ΔAge | **0.948** | ranking capability: excellent |
| ridge ΔAge | **0.955** | ranking capability: excellent |
| RES score | 0.686 | **the RES formula degrades a good ranking** |

**The model ranks correctly. The RES transform was throwing that away** — and we now know why
(it multiplies a good ΔAge by three defective signals). **"Ranking is broken" was never true and
must not appear in any writeup.**

---

# PART II — DECISIONS

## 4. RES: verdict deferred, with a defined re-test

**Why v2's "demote RES" was wrong.** RES = `Φ(S)·S^k·g(R_eff)·exp(−λ·P_loss)`. Every input is
currently defective:

| RES input | Defect | Fixable? |
|---|---|---|
| `mu_age` (via `R_eff`) | per-donor level shift ±12.7 yr → `R_eff = 0` everywhere (T7.4.2) | **YES — demonstrated**, k=3 (T16) |
| `sigma_age` (via `R_eff`) | intervals cover 0.40 vs 0.90 → sigma understated (T14) | **YES — cause identified** |
| `S`, `P_loss` | ECE 0.28 out-of-donor (T8.2) | **YES — Platt halves it** |
| `in_dist` | AUC 0.47 ≈ chance, zeroes RES arbitrarily (T15) | **YES — fix or disable** |
| `lam` | `lam = 0` makes the `P_loss` term **inert** (T7.4.2) | config, trivially |

**A score built on four defective inputs cannot be judged on its formula.** The honest position:
**RES is untested, not refuted.**

**Defined re-test — Test 7.6 (after fixes):** rerun Tests 7, 7.1 and 7.4.1 with corrected
`mu`, `sigma`, fate probabilities, and OOD disabled.

| Outcome | Decision |
|---|---|
| `rank_res` within **0.03** of the ΔAge sort **and** over-approval gap ≤ 0.5 | **RES is vindicated** — keep as the headline score |
| within **0.05**, gap ≤ 1.0 | keep as a **secondary safety filter**, rank by ΔAge |
| still ≥0.10 below the ΔAge sort | **then** demote — and only then is the formula itself implicated |

**Independent of the re-test, one structural finding stands (T13):** the trajectory is **biphasic
on 4/6 donors**, so scoring each cell *in isolation against absolute thresholds* presumes
monotonic rejuvenation. Even a vindicated RES should be reconsidered as a **trajectory-aware**
score. That is a design improvement, not a condemnation.

## 5. What to fix — ONE root cause, then one protocol change

### 5a. The root cause (found by auditing `train_model.py`)

**Every calibration parameter in the bundle is fitted on data from the TRAINING donors, then
applied to a held-out donor whose error regime is completely different.** Three lines:

| Parameter | Fitted on (source line) | Out-of-donor symptom | Test |
|---|---|---|---|
| `temperature` | `fit_temperature(cal_logits, cal_target)` — **val** split | fate ECE **0.28** | T8.2 |
| `conformal.q` | `fit_conformal(abs_res, …)` — **calib** residuals (MAE≈4 yr) | coverage **0.40 vs 0.90** | T14 |
| `ood` reference | `fit_ood(train_feats)` — **train** trunk features | AUC **0.47** ≈ chance | T15 |

The per-donor level shift is the same family: the model is **unbiased in-distribution**
(calib offset −0.03) and shifted only out-of-donor (T7.4.3).

> **This is not four defects. It is one architectural mistake with four manifestations:
> the system calibrates itself against donors it has already seen.**

That reframing matters for the writeup: "we calibrated in-distribution and deployed out-of-donor"
is a *generalizable methodological finding*, not a list of bugs.

### 5b. FIX A — cross-donor calibration (pure code, no new data)

> Refit **all three** calibration parameters on **cross-donor** statistics: run an inner
> leave-one-donor-out *within the training donors*, pool the resulting out-of-donor logits,
> residuals and trunk features, and fit `temperature`, `q` and the OOD reference on **those**.

- Touches three lines in `train_model.py` (116, 130, 134) plus a shared helper producing the
  inner-LODO statistics.
- Requires **no new experimental data** — it only changes which split the calibrators see.
- `fit_conformal` already applies the finite-sample `n/(n+1)` correction, so no extra work there.

**Expected effect (falsifiable):** `q` rises from **8.9** toward **~30–40 yr** on the *uncorrected*
model; coverage 0.40 → ~0.90; fate ECE 0.28 → ~0.13; OOD AUC moves off 0.47 or the gate is dropped.
**Wider intervals are correct, not a regression** — out-of-donor ΔAge genuinely carries that
uncertainty, and saying so is a finding.

### 5b-bis. CRITICAL CORRECTION — `sigma_age` is NOT the conformal `q`

An earlier draft of this plan assumed refitting the conformal quantile would fix `R_eff`.
**It would not.** They are different parameters:

| quantity | source | value | used by |
|---|---|---|---|
| `sigma_age` | `age.std(0)` — **ensemble / MC-dropout spread** | ~**2.4 yr** | **`R_eff`**, i.e. RES |
| `q` | `fit_conformal(...)` — conformal quantile | ~**8.9 yr** | the reported interval, and Test 14 |

`R_eff = max(0, −(mu + z_conf·sigma_age))` uses the **ensemble spread**, not `q`. So refitting `q`
alone fixes the reported intervals and leaves RES exactly as broken.

**And the ensemble spread is the more overconfident of the two:** members disagree by ~2.4 yr while
the true out-of-donor error is ~14 yr. The ensemble agrees with itself while being collectively
wrong — the classic failure of ensemble uncertainty under distribution shift.

**Therefore Change A must also address `sigma_age`**, either by rescaling it with a cross-donor
factor, or — cleaner — by making `R_eff` use the calibrated bound directly:
`R_eff = max(0, −(mu + q))`. That matches the stated intent ("the upper age bound must be
negative"); the ensemble spread was only ever a proxy for it.

### 5b-ter. The consequence — and it is the most important finding in this review

Once `R_eff` uses an **honest** uncertainty, the arithmetic is unforgiving:

| uncertainty used | mu | u | R_eff | g |
|---|---|---|---|---|
| ensemble spread (miscalibrated) | −11.0 | 2.4 | 8.6 | 0.63 |
| honest, uncorrected model | −11.0 | ~39 | **0.0** | 0.00 |
| honest, after level correction | −11.0 | ~19 | **0.0** | 0.00 |

**With honest per-cell uncertainty (~19 yr) and a real effect (~11 yr), R_eff = 0.**
**Per-cell "confident rejuvenation" is unachievable at this data scale** — not because of a bug,
but because the uncertainty genuinely exceeds the effect size. **RES's current design asks a
question the data cannot answer.**

**But the same arithmetic points at the fix.** Uncertainty on a *mean* shrinks by √n:

| n cells | SE of mean (q=17–21) | effect −11.35 |
|---|---|---|
| 5 | 7.6–9.4 yr | detectable |
| 10 | 5.4–6.6 yr | detectable |
| 21 | 3.7–4.6 yr | **comfortably detectable** |

> **RES should score CONDITIONS (populations of cells), not individual cells.** That also matches
> the real use case — you rank reprogramming *conditions* for a patient, not single cells. This is
> a **design change with a quantitative justification**, and it converts an unanswerable question
> into an answerable one.

**This revises the Test 7.6 expectation (§4):** per-cell RES will very likely *still* approve
nothing after changes A and B, and that will be the **correct** result. The real test is whether a
**condition-level** RES beats a condition-level ΔAge sort.

### 5c. FIX B — per-donor level correction (protocol change, needs reference cells)

> From **k ≈ 3 cells of the new donor with known true ΔAge**, estimate
> `d = median(pred − true)` and subtract it.

**Conditional rule (concrete).** Compute the estimate's standard error
`SE ≈ 1.253·sd/√k` and **apply the correction only if `|d| > 2·SE`.**
*(Numerically verified: at k=3 this asymptotic formula **overstates** the true SE of a median by
~8% — i.e. it errs toward *not* correcting, which is the safe direction. SE of the **mean** is
exact (`sd/√k`) if you prefer the simpler estimator; the median matches what T16 measured.)*
At k=3 this declines O1's 0.6 yr shift — the donor the correction damaged — while accepting N2's
15.0 and N3's 28.3.

**FIX B IS NOT A PURE CODE FIX.** It requires clock readings on control *and* perturbed samples
for k≈3 cells of every new donor. Small and concrete, but it *is* an experimental ask; without it
the model stays a within-donor ranker.

### 5d. Ordering — and why it is A → B → A′

Fix B changes the error scale that Fix A's `q` depends on, so calibration must be refitted after
it:

| step | state | mean \|error\| | required q (P90 ≈ 2.5–3.0× mean for our heavy-tailed mixture) |
|---|---|---|---|
| now | in-distribution calibration | 14.3 | **8.9** (far too small) |
| after **A** | cross-donor calibration | 14.3 | ~35–43 |
| after **B** | level corrected | 6.9 | ~17–21 |
| after **A′** | recalibrated on corrected residuals | 6.9 | **~17–21** ✓ |

*(P90/mean verified numerically: 2.07 for normal, 2.31 for exponential, **2.67** for a heavy-tailed
mixture like ours — O1 MAE 5.4 vs N3 29.7. Hence a range, not a point estimate.)*

**Why not B first?** B alone leaves `sigma` understated while `mu` is corrected — the T14 coupling
that makes RES ~2× too permissive. **A first is the safe direction** (over-wide intervals →
RES too *conservative*, never too permissive).

| # | Change | Kind | Status |
|---|---|---|---|
| **A** | cross-donor calibration (temperature + q + OOD) | code only | root cause identified; 3 lines |
| **B + A′** | per-donor level correction **bundled with** a calibration refit | protocol + code | B demonstrated (T16); bundled because B changes the error scale A′ depends on |
| **C** | Test 7.6 — re-evaluate RES with A and B applied | evaluation | the deferred verdict (§4) |

## 6. What is genuinely NOT fixable by code

Short list, and it is short on purpose:

- **Cross-perturbation generalization.** One cocktail. Needs new data (§11).
- **Proving time-redundancy is benign** (T11.1) — state and time are collinear here.
- **Statistical power on the fate edge** — 4/6 folds saturate near 1.0; needs more donors.

Everything else on the defect list has an identified fix.

---

# PART III — WORK REMAINING

## 7. Testing is complete

Tests 13–16 closed the last blocking questions. **No further diagnostics are required before
changing code.**

**But two tests gate PUBLICATION, not code.** Both can change what may be claimed:

| Test | Question | Why it gates a headline claim |
|---|---|---|
| **Test 19 — the second clock** | `clock_fit.py` fits a **linear** clock, but the published Fleischer clock is a **nonlinear LDA ensemble**. Test 6's "nothing beats ridge" is therefore **conditional on a linear target we generated ourselves** | A reviewer will call this close to circular. Implement a nonlinear `AgingClock` subclass, regenerate ΔAge, rerun Test 6. **Ridge still wins → the claim strengthens. A flexible model wins → the finding becomes "deep ties linear *when the clock is linear*", which is more novel. Either outcome improves the paper** |
| **Test 8.3 — the Y1 probe** | The fate-beats-linear claim rests on **one fold**: four saturate near PR-AUC 1.0, the whole signal is Y1 (0.961 vs 0.636) | If Y1 is an artifact the claim must be withdrawn. **Better we find it than a reviewer** |

Cross-state generalization and full-transcriptome remain genuinely optional and gate nothing.

## 8. The change protocol — `scorecard.py`

```bash
python scorecard.py snapshot --tag baseline    # BEFORE any change
# ... exactly ONE change ...
python scorecard.py snapshot --tag f1_conformal
python scorecard.py compare baseline f1_conformal
```

`compare` reports, per metric, the **paired per-fold difference and 95% CI**, and returns
**ACCEPT (better)** / **noise (CI incl. 0)** / **REGRESSION** / **(context)**.

> A change is accepted only if the **TARGET** metric says ACCEPT **and** no **GUARD** metric says
> REGRESSION.

## 9. Pre-registered criteria per change

| Change | Run | TARGET | Threshold | GUARDS (must not REGRESS) | If it fails |
|---|---|---|---|---|---|
| **A** cross-donor calibration | `compare baseline A_xdonor` | `conformal_coverage` **and** `fate_ece` | coverage reach **0.85–0.95** (from 0.40); ECE **ACCEPT** + ≥40% drop (0.28 → ≲0.17) | `fate_prauc`, `fate_roc`, `rank_model_dage`, `dage_mae_model` must all be **noise** — A changes only calibration, so if discrimination or ΔAge moves, the implementation touched something it should not | inspect which of the three refits failed; they are independent and can be adopted separately |
| **B + A′** level correction + refit | `compare A_xdonor B_percalib` | `dage_mae_model` | **ACCEPT** + ≥25% drop (T16 predicts ~50%) | `rank_model_dage` **noise or ACCEPT, never REGRESSION** (a level shift is rank-invariant — if ranking moves, the implementation is wrong); no fold may worsen >20%; `conformal_coverage` stays 0.85–0.95; `conformal_width` should **fall** (~35–43 → ~17–21 yr) | revert the bundle; keep A; report as a within-donor ranker with honest wide intervals |
| **C** Test 7.6 (RES verdict) | rerun T7 / 7.1 / 7.4.1 | `rank_res` | see §4 table | over-approval gap must not exceed baseline | only *then* is demoting RES justified |

**A note on `conformal_width`.** It is expected to **rise** under A (8.9 → ~35–43 yr) and **fall**
under B+A′ (→ ~17–21 yr). A rise under A is **correct behaviour**, not a regression — the current
narrow interval is the defect. The scorecard marks `conformal_width` "lower is better", so read
this row against the plan, not the arrow.

## 10. Roadmap

**Phase 0 — freeze.** `scorecard.py snapshot --tag baseline`.

**Phase 1 — Change A: cross-donor calibration.** Pure code, no new data, fixes three symptoms at
once (coverage, fate ECE, OOD). Safe direction: it can only make RES more conservative. This is
also the natural first change because it *validates the scorecard pipeline* on something whose
expected effect is known.

**Phase 2 — Change B + A′: per-donor level correction, bundled with a calibration refit.**
Largest single expected gain (MAE −50%). Requires k≈3 reference cells per donor — confirm that
experimental ask is acceptable **before** implementing. Verify per-fold that no donor regressed
>20% and that `rank_model_dage` did not move.

**Phase 3 — Change C: Test 7.6, the deferred RES verdict.** Only here is RES judged, on corrected
inputs, against the §4 thresholds.

**Phase 4 — REQUIRED before publication: the Y1 investigation.** *(v3 listed this as "optional"
while the risk register rated it High impact — a contradiction, now resolved in favour of the risk
register.)* The fate-beats-linear claim is a headline result and **rests almost entirely on one
fold**: four folds saturate near PR-AUC 1.0, and the whole signal is Y1 (0.961 vs logreg 0.636).
If Y1 is a quirk, the claim evaporates. **Do not publish the fate edge without understanding why
Y1 is hard.** Working hypothesis (from T7.1): Y1 carries the most transcriptomic stress — 8 unsafe
cells of 19, the most of any fold — and the linear model cannot see through that noise.

**Phase 5 — genuinely optional.** Cross-state generalization; trajectory-aware RES redesign if
T13's biphasic finding warrants it.

**Phase 6 — consolidate and write up** (§12).

---

# PART IV — BEYOND THE CODE

## 11. Data strategy

| Limit | Consequence | Fix |
|---|---|---|
| **One cocktail (OSKM), time-varying only** | cross-perturbation thesis **unposed** | a dataset varying the **cocktail** (≥5 regimes) |
| **6 donors** | per-donor shift unlearnable without reference cells; fate edge rests on one fold | **≥15 donors** |

Priority asks for the next dataset: (1) multiple distinct perturbations; (2) more donors;
(3) reference/control samples per donor — which would make the level correction automatic; (4) varied time at fixed
state.

**Note:** T16 already shows **k≈3 labelled reference cells per donor** suffices for the level fix.
That is a small, concrete experimental ask, not a redesign.

## 12. Publication strategy

**Lead with what works — it is more than v2 implied:**
- **Within-donor ranking, Spearman 0.925–0.983** across six held-out donors.
- **Fate discrimination out-of-donor, PR-AUC 0.96–1.00**, holding where a linear baseline
  collapses (Y1: 0.961 vs 0.636) — reported as suggestive (n=5, 4 folds saturated).
- **A parameter-free, leak-safe cross-modality harmonization method.**
- **A rigorous null result with a structural explanation:** deep models match but do not beat
  linear baselines for clock-ΔAge on small multi-donor data — because the clock is linear.
  Consistent with *Nature Methods* 2025.
- **A quantified account of per-donor calibration**: ±12.7 yr level shift, and **k≈3 reference
  cells recover ~50% of the error** — a practical protocol recommendation.
- **A biphasic trajectory finding (4/6 donors)** and its consequence: absolute per-cell scoring is
  mis-specified for a non-monotonic process.

**State plainly, without overclaiming:**
- Cross-perturbation generalization is **untested** — the data cannot pose it.
- Absolute cross-donor ΔAge requires reference cells.
- Uncertainty and OOD components required recalibration; report before/after.

**Do not write:** "ranking is broken" (false), "RES doesn't work" (untested with correct inputs),
or any claim of cross-perturbation generalization.

## 13. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Fate edge is a Y1 artifact | Medium | High | investigate Y1 before publishing |
| the level correction hurts well-calibrated donors (O1, Y1) | **Confirmed** | Medium | conditional application; report per-fold |
| the conformal refit over-widens intervals into uselessness | Medium | Medium | width/2sd guard in §9 |
| RES stays poor even after changes A and B | Medium | Low | that is what change C is for — a real verdict, not a guess |
| Reviewer runs RES vs ΔAge-sort | High | High if unreported | report it, with the input-defect explanation |
| **Over-correcting into "the model is worthless"** | **Realized in v2** | High | this document is the correction |

---

## 14. What exactly each change buys

**Read this before committing effort to any of them.** Each entry states the concrete gain, the
evidence class behind the number, what it does *not* buy, and what it costs.

**Evidence classes:** **[M]** measured in a test · **[A]** derived from arithmetic on measured
quantities · **[U]** genuinely unknown until run.

---

### Change A — cross-donor calibration *(temperature + `q` + `sigma_age` + OOD)*

| Gain | From | To | Class |
|---|---|---|---|
| Conformal coverage | **0.40** (0.00 on N2/N3) | **~0.90** | [A] |
| Fate ECE | **0.28** | **~0.13** | **[M]** T8.2 |
| Cells wrongly discarded by OOD | **27%** (48% on N2) | 0% if the gate is dropped, or a justified rate if it becomes informative | **[M]** T15 |
| `sigma_age` honesty | ensemble spread 2.4 yr vs true error 14 yr | calibrated | [A] |

**What it buys: TRUSTWORTHINESS, not accuracy.** After A the model is no better at predicting —
it becomes honest about what it doesn't know. For a safety-critical decision (which reprogramming
conditions are safe), **a miscalibrated confidence is worse than no confidence**: right now the
intervals cover **zero** cells on two donors while claiming 90%.

**What it does NOT buy:** ΔAge accuracy (MAE stays 14.3), ranking (0.948 unchanged), or RES
approvals — it will likely make RES *more* conservative, which is correct.

**Cost:** pure code. Three lines in `train_model.py` plus a shared inner-LODO helper. **No new
experimental data.**

---

### Change B + A′ — per-donor level correction *(bundled with a calibration refit)*

| Gain | From | To | Class |
|---|---|---|---|
| ΔAge MAE (aggregate) | **14.3** | **~7.0 (−50%)** | **[M]** T16 |
| N2 MAE | 21.8 | **7.1** | **[M]** |
| N3 MAE | 29.7 | **10.0** | **[M]** |
| Per-donor level shift | ±12.7 yr | within the estimator's noise | **[M]** |
| Interval width (after A′) | ~35–43 yr | **~17–21 yr** | [A] |

**What it buys: absolute ΔAge becomes usable out-of-donor for the first time.** Today a prediction
can be wrong by ±12.7 yr in a donor-specific direction, so *no* absolute-threshold decision is
defensible. This is the single largest measured gain available.

**What it does NOT buy:** ranking (a level shift is rank-invariant — ranking *must* stay at 0.948,
and if it moves the implementation is wrong), and it does **not** by itself make per-cell RES work
(§5b-ter).

**Cost — this is the one with a real price:**
- **k ≈ 3 cells per new donor with known true ΔAge** — clock readings on control *and* perturbed
  samples. A protocol change, not a code change.
- **Hurts 2 of 6 donors** unless applied conditionally (O1: 5.4 → 6.1; Y1: 7.3 → 8.5), because
  correcting a donor that needs no correction only injects noise.

---

### Condition-level RES *(the redesign §5b-ter argues for)*

| Gain | From | To | Class |
|---|---|---|---|
| Uncertainty on the scored quantity | ±17–21 yr (per cell) | **±3.7–4.6 yr** (mean of 21 cells) | [A] |
| Detectability of an ~11 yr effect | **impossible** (u > effect) | **comfortable** (u ≈ ⅓ of effect) | [A] |

**What it buys: it makes the question answerable at all.** Per-cell confident rejuvenation is
arithmetically out of reach here — honest uncertainty (~19 yr) exceeds the real effect (~11 yr).
Averaging over a condition's cells shrinks uncertainty by √n and puts the effect comfortably
inside detection range. It also matches the actual use case: you rank **conditions** for a
patient, not individual cells.

**What it does NOT buy:** any improvement in the underlying predictions — it changes *what is
scored*, not how well it is predicted.

**Cost:** a genuine redesign of the RES aggregation layer. Should be treated as a **new hypothesis
with its own pre-registered test**, not a tweak.

---

### Change C — Test 7.6 *(re-evaluate RES on corrected inputs)*

**Buys a defensible verdict instead of a guess.** RES currently fails, but every input it consumes
is defective, so its formula has never actually been tested. All three outcomes are worth having:

| Outcome | What you gain |
|---|---|
| RES matches the ΔAge sort | the headline score is recovered and justified |
| RES needs condition-level aggregation | a specific, quantitatively motivated redesign |
| RES still fails on corrected inputs | demotion **with evidence** — a claim that survives review |

**Cost:** evaluation only, no code change. **[U]** — the outcome is genuinely unknown, which is
precisely why it must be run rather than assumed.

---

### Y1 investigation *(required before publishing the fate claim)*

**Buys protection for the single strongest positive result.** Fate discrimination beating a linear
baseline is the model's clearest win — and four of six folds saturate near PR-AUC 1.0, so the whole
signal is Y1 (0.961 vs logreg 0.636). Either the claim survives scrutiny and becomes publishable,
or it is an artifact and **we find out before a reviewer does**. **[U]**

**Cost:** one analysis, no code change.

---

### The honest counterfactual — what doing NOTHING already gets you

Stated so the fixes are judged against the real alternative, not against zero:

- **Within-donor ranking, Spearman 0.925–0.983** across six held-out donors — already publishable.
- **Fate discrimination out-of-donor, PR-AUC 0.96–1.00** — already publishable (with the Y1 caveat).
- **A parameter-free, leak-safe harmonization method** — already publishable.
- **A rigorous null result with a structural explanation** (linear clock ⇒ linear target ⇒ ridge
  optimal) — already publishable, and consistent with *Nature Methods* 2025.

> **The fixes make the system more useful and more honest. They are not what makes it valid.**
> Change A is the one closest to mandatory, because shipping intervals that cover 0% while
> claiming 90% is a correctness problem, not a polish problem.

---

### What none of this fixes

| Limitation | Why no change helps |
|---|---|
| **Cross-perturbation generalization** | one OSKM cocktail — the question cannot be *posed*, let alone answered |
| **Statistical power on the fate edge** | 4/6 folds saturate; needs more donors, not more code |
| **Proving time-redundancy is benign** | state and time are collinear in this dataset |

**All three need data, not engineering.** That remains the highest-leverage investment in the
project (§11).

---

## 15. The forward model — scope, and the identifiability wall

### 15a. It is an addition, not a rewrite

`CellFateNet.forward` already returns the trunk latent `z`:
`return self.cls_head(z), self.age_head(z), z`. A forward model is a module **on** `z`:

```
now:   (x, u) → z → {fate, ΔAge}
add:   (z, Δt) → z'          ← new dynamics module
       z' → {fate, ΔAge}     ← reuses the existing heads, unchanged
```

One new module, trained on the same bundles. The original design already contemplated a
WorldModel. **Architecturally this is modest.**

### 15b. But architecture was never the blocker — identifiability is

**Test 11.1's finding changes sign depending on what is being built:**

| Building a… | State/time collinearity is… |
|---|---|
| **query** model ("what if I ran 5 more days?") | **fatal** — the model reads time off state, so the time knob does nothing (0.035 yr) |
| **forward** model ("given this state, what comes next?") | **the mechanism** — state encodes trajectory position, which is exactly what you propagate |

The property that broke the query model is what makes the forward model well-posed.

### 15c. What is learnable from the existing data, and what is not

| Question | Learnable? | Why |
|---|---|---|
| How does the trajectory evolve in time? | ✅ | 12 timepoints × 6 donors; state encodes position |
| What happens if I stop at day 13 vs day 20? | ✅ | a time question along the measured trajectory |
| What is the optimal **dose**? | ❌ | **one dose** — not identifiable by any model |
| Which **cocktail** is best? | ❌ | one cocktail |

> **"Best dose and time" splits in half. Time is answerable. Dose is not** — that is an
> identifiability wall, not a modelling gap, and it needs a **dose × time grid** to cross.

### 15d. The time half is the dataset's own headline question

Gill's data is **MPTR — Maturation Phase Transient Reprogramming**. The premise is *transient*:
reprogram, then withdraw. **The central question is when to stop** — long enough to rejuvenate,
short enough to avoid identity loss. That is precisely the trade-off between the two heads the
model already has.

### 15e. Hard constraints on any implementation

- **~126 age-valid samples.** A generative model over 2000 genes is hopeless. Work in **latent
  space** (`z` is 256-d; reduce to ~10–20 components first).
- **Population-level, not cell-level.** Sampling is destructive — the same cell is never seen
  twice. This is distribution-at-t → distribution-at-t+Δ, the PRESCIENT / Waddington-OT problem
  class. Established, but genuinely hard at n=126.
- **Scope and name it honestly:** a *transient-reprogramming stopping-time model*, **not** a dose
  optimizer.

### 15f. What it would deliver

1. **Optimal stopping time** — sweep Δt forward; find where P(safe) is still high and ΔAge most negative.
2. **Early abort** — predict the day-20 outcome from the day-5 state; kill bad runs early.
3. **Trajectory-level scoring** — pairs with the condition-level RES insight (§5b-ter): uncertainty
   ±4 yr on a trajectory vs ±19 yr on a cell.
4. **Counterfactual safety** — "if I withdrew now, what is P(loss of identity)?"

**Do not promise dose–response.** Without a dose × time grid it is unidentifiable, and claiming it
would be the one genuinely indefensible thing in the whole project.
