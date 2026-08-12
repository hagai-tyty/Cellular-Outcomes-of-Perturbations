# STAGE 3 — The stopping-time tool

**Implements:** `MASTER_PLAN.md` §15 (forward model) **and** §5b-ter (condition-level RES) — the
same object, see §0.3.
**Depends on:** Stage 1 **required**; Stage 2 optional.
**Scope:** 3 new modules (~600 lines), 4 internal sub-stages.

> ## 🆕 ADDED 2026-07-26 — the tool needs an INTERNAL CONTROL at inference time
>
> *Additive note; nothing below is modified. Added after `STAGE_1_5_1_REVISED.md`.*
>
> **A ΔAge computed against a day-0 baseline is not trustworthy on reprogramming intermediates** — it
> carries an identity-change artefact measured at **+36.5 yr over 11 days**, enough to invert the
> sign of a real ~30 yr effect. The fix (a contemporaneous non-responder control) is straightforward
> for *training* data, but the tool must also do it **at inference**, and that is an open design
> question this stage owns:
>
> - **Single-cell input (expected case):** score the culture's non-reprogramming cells as the internal
>   control and report ΔAge as responders **minus** that control. Natural, and it also gives the
>   condition-level aggregation §5b-ter already requires.
> - **Bulk input:** **no internal control exists.** Either require a paired non-responder sample, or
>   refuse to report absolute ΔAge and fall back to ranking only. **Decide this explicitly** — a bulk
>   ΔAge with no control is exactly the number that reads +36.5 yr on cells that did nothing.
>
> **Add to the sub-stage that defines the inference contract:** the input schema must carry (or let
> the tool derive) a responder/non-responder split, and the tool should **decline to report ΔAge**
> when it cannot — consistent with this project's "can say no" design.
>
> **Unaffected:** the fate/safety head needs no control and is not touched by this. Test 18's STOP on
> the forward Δt signal is a separate, still-open question.

| Sub-stage | Produces | Blocking for |
|---|---|---|
| **3a Gate** | GO / WEAK GO / STOP | everything below |
| **3b Data** | `forward_pairs.py` | 3c |
| **3c Training** | `train_forward.py` | 3d |
| **3d Decision** | `stopping.py` | Stage 4 |

---

## 0. Framing

### 0.1 What is being built

A researcher submits a transcriptome of their culture mid-protocol and receives, for each
candidate withdrawal time: predicted ΔAge with honest error bars, a calibrated P(identity loss),
and a recommended withdrawal window.

### 0.2 Why this product survives what the data cannot support

It asks a **relative** question along the **one axis this data varies**. Comparing "day 15 vs day
21 **for the same donor**" cancels the ±12.7 yr per-donor level shift. **No dose axis, no
cross-perturbation generalization** — the two things the data cannot pose.

### 0.3 The convergence

`MASTER_PLAN.md` §5b-ter argues RES must move from per-cell to **condition-level** scoring. The
arithmetic that forces it — `R_eff = max(0, −(mu + u))`, credit requires the *upper* bound negative:

| uncertainty used | mu | u | R_eff | g |
|---|---|---|---|---|
| ensemble spread (miscalibrated) | −11.0 | 2.4 | **8.6** | 0.63 |
| honest, uncorrected model | −11.0 | ~39 | **0.0** | 0.00 |
| honest, after level correction | −11.0 | ~19 | **0.0** | 0.00 |

**With honest per-cell uncertainty, R_eff = 0** — uncertainty (~19 yr) exceeds the real effect
(N2's true median ΔAge is **−11.35 yr**). But uncertainty on a *mean* shrinks by √n:

| n cells | SE of mean (q = 17–21) | vs effect −11.35 |
|---|---|---|
| 5 | **7.6–9.4 yr** | marginal |
| 10 | 5.4–6.6 yr | detectable |
| 21 | **3.7–4.6 yr** | **comfortably detectable** |

**That redesign and this tool are the same object.** Built once, in 3d.

### 0.4 Hard constraint discovered in the code

```python
# src/cellfate/common/constants.py:66
N_DOSE_TIME: int = 2               # [log10(dose_uM), log(time_h)]
# src/cellfate/common/schemas.py:76  -- a validator RAISES on any other length
```

So the forward model gets its **own bundle**, reusing the width-2 tensor but reinterpreting
column 1:

| | column 0 | column 1 |
|---|---|---|
| scoring bundle (existing, untouched) | log10 dose | **log absolute time** |
| forward bundle (new) | log10 dose | **log Δt** |

A `mode` field records which, **asserted at load** (3d.6).

> ⚠️ **BLOCKER FOUND IN AUDIT.** `BundleMeta` is declared `model_config = ConfigDict(extra="forbid")`
> (`schemas.py:269`). **You cannot simply set `bundle_meta.mode`** — pydantic will reject the extra
> key and the bundle will fail to write. The field must be **declared** first:
>
> ```python
> # src/cellfate/common/schemas.py, inside class BundleMeta
> mode: str = "scoring"          # "scoring" | "forward"; defaulted so OLD BUNDLES STILL LOAD
> ```
>
> A declared field with a default is compatible with `extra="forbid"`, and every existing bundle
> keeps validating. **Do this before writing any forward bundle.**

---

# SUB-STAGE 3a — The gate

**Run before writing any tool code:**

```powershell
python experiments/test18_forward_gate.py
```

**What it tests.** Test 11.1 showed the current model ignores its time input entirely (ΔAge shifts
**0.035 yr** across a full sweep) because time is redundant with state along one trajectory. The
hope is that training on `(t_i → t_j)` **pairs** fixes this: with pairs the same starting state
appears with *different* Δt and *different* targets, so Δt is no longer redundant.

**That is a hypothesis.** The gate checks it with **ridge** — because this project has established
repeatedly (T3, T5, T6, T9) that at this scale ridge matches or beats flexible models. **If ridge
cannot find a Δt signal, a neural net will not either.**

| Verdict | Condition | Action |
|---|---|---|
| **GO** | Δt improves prediction beyond noise **and** the sweep moves >2 yr | build 3b–3d |
| **WEAK GO** | sweep moves but Δt does not beat state-only | build, tempered |
| **STOP** | neither | **do not write tool code.** Ship the scoring model; go to Stage 5 |

> ## 🔬 **3a DIAGNOSIS — 2026-08-12. A forward time signal IS present; 3a's estimator could not have found it, and its bar could not have registered it.**
>
> *Additive. The withdrawal banner and the result box below are unchanged apart from two dated
> annotations inside them.* Artefacts: `experiments/stage3a_diagnose.py`,
> `results/stage3a_diagnose_results.json`, `tests/test_stage3a_diagnose.py` (20 tests).
> **READ-ONLY** — no retrain, no rebuild, `src/` untouched. `test18_forward_gate.py` is imported
> **unmodified** and its own `timepoint_table` / `build_pairs` / `feats` / `paired_ci` are the
> primitives, so what is diagnosed is the estimator 3a actually ran.
>
> **Anchored first.** All ten of 3a's Part A/C numbers reproduce to the digit (Y1's 311.47 and
> 7.589 included). The fast ridge path used by §D4 — features frozen so only the target varies —
> is checked against `sklearn.linear_model.Ridge` on the real targets before use: **max|Δ| =
> 3.9e−10**, and pinned as a property in the test file.
>
> ### D0 — the finding. A model-free predictor using only t_j beats every arm 3a ran
>
> *Not pre-registered; added once the target's shape became visible.* Part C's target is the mean
> of **1.6–1.8 binary cells per timepoint**, and **63/70 = 90 % of its values sit exactly at 0 or
> 1**. Per donor it is close to a monotone **step in time**:
>
> | donor | timepoints | cells/tp | unsafe fraction by timepoint | at 0 or 1 | std |
> |---|---|---|---|---|---|
> | **N2** | 11 | 1.7 | `0 0 0 0 0 0 0 0 0 0 0` | 11/11 | **0.000** |
> | N3 | 12 | 1.7 | `0 0 0 .5 0 0 0 0 0 1 1 1` | 11/12 | 0.431 |
> | O1 | 12 | 1.8 | `0 0 0 0 0 0 0 0 1 1 1 1` | 12/12 | 0.471 |
> | O2 | 12 | 1.7 | `0 0 0 0 0 0 0 0 1 1 1 1` | 12/12 | 0.471 |
> | Y1 | 11 | 1.6 | `0 0 .5 .75 .5 .5 .5 0 1 1 1` | 6/11 | 0.376 |
> | Y2 | 12 | 1.8 | `0 0 0 0 0 0 .5 0 1 1 1 1` | 11/12 | 0.462 |
>
> If the target is essentially "when does this culture turn", the forward question has a
> **model-free ceiling**: predict the held-out donor's value at `t_j` from the **other** donors'
> value at that same `t_j`. That uses only `t_j = t_i + Δt` — no genes, no fitting, and nothing
> from the held-out donor (pinned by a leakage test). It is the most forward-in-time information
> Δt can carry at this geometry.
>
> | Part C, the 5 folds 3a graded | MAE |
> |---|---|
> | pooled-mean baseline | 0.493 |
> | 3a's **state + Δt, raw** (the arm the STOP was read from) | **2.010** |
> | 3a's state only | 0.409 |
> | state + Δt, bounded (logit) | 0.352 |
> | **oracle on `t_j` alone** | **0.157** |
>
> **paired oracle − pooled = −0.336, 95 % CI [−0.450, −0.221], n = 5 → `t_j` HELPS**, decisively;
> the CI is nowhere near 0. The oracle beats 3a's best arm by **2.24×** and its state-only by
> **2.6×**. Per fold the gain is +0.329, +0.400, +0.400, +0.180, +0.369 — **every fold, same
> direction.**
>
> The same holds for Part A's target: ΔAge oracle **14.208** against a pooled-mean 23.615 and 3a's
> state-only 20.67, paired **−9.407, CI [−16.332, −2.482] → `t_j` HELPS**.
>
> **So a forward time signal exists in this corpus, on the exact geometry 3a graded**, and 3a
> reported that it does not.
>
> #### The limit on that claim, stated as plainly as the claim
>
> The five graded folds exclude **N2**, and N2 is the one donor that **never becomes unsafe** —
> flat 0 across all 11 timepoints. Include it and the oracle's advantage falls to **−0.183, CI
> [−0.411, +0.046] → tied**, with N2 the only negative fold (**−0.240**): the shared time course
> predicts a transition that N2 never has. N2 is absent from 3a's grading for an unrelated reason
> (C-7 rule 4 masks its **ΔAge**, and `timepoint_table` filters on the age mask before computing a
> target that does not depend on ΔAge at all), and it is separately ungradable because a
> zero-variance target is `test18`'s own documented skip. **So: established on the graded
> geometry; not established across all six donors.** The point estimate favours `t_j` either way.
> Resolving N2's status is a precondition of any re-run.
>
> ### D1 — the divergence is real, and BOTH proposed mechanisms were wrong
>
> The audit proposed that Y1's Δt sits outside the training folds' support. **It does not.** All
> five folds share the *identical* Δt range **[0.14, 11.77]** with **0 held-out pairs outside
> training support** — Y1's missing timepoint is **interior**, so only its pair count is smaller
> (55 vs 66). My own pre-registered guess (**P2**: the Δt block carries it) is **also refuted** —
> Y1's Δt coefficient (2.089) and mean |Δt contribution| (1.634) sit inside the other folds' range
> (1.745–1.832, 1.356–1.443).
>
> **It is the gene block.**
>
> | held-out | \|gene\| state-only | \|gene\| state+Δt | ratio | cos(w_gene) | mean\|z_gene\| train | mean\|z_gene\| held-out | pred outside [0,1] |
> |---|---|---|---|---|---|---|---|
> | N3 | 0.134 | 1.563 | 11.7× | +0.475 | 0.604 | 0.503 | **77 %** |
> | O1 | 0.147 | 1.440 | 9.8× | +0.485 | 0.604 | 0.455 | 33 % |
> | O2 | 0.162 | 1.513 | 9.4× | +0.473 | 0.604 | 0.488 | 62 % |
> | **Y1** | 0.373 | **7.614** | **20.4×** | +0.517 | 0.671 | **2.546** | **100 %** |
> | Y2 | 0.139 | 1.074 | 7.7× | +0.477 | 0.605 | 0.448 | 33 % |
>
> Adding the two Δt columns re-solves the ridge and **rotates the gene-weight vector by ~60°**
> (cos ≈ 0.47–0.52), amplifying the gene block's contribution **7.7×–20.4× on every fold**. Y1's
> held-out state is the one that sits far outside the training scaler's support — mean |z| **2.546
> against 0.671** in training, while every other fold's held-out mean |z| is *below* its own
> training value. Rotated weights × out-of-support state = predictions of **5.44 to 10.92** on a
> target bounded by 1. **And it is not only Y1** — the raw estimator emits out-of-range predictions
> on **every fold**: 33 %, 33 %, 62 %, 77 %, 100 %.
>
> ### D2 — bounding the predictor repairs the instrument; it does not, by itself, find the signal
>
> | Part C, 5 folds | mean-only | state | state+Δt | paired mean | 95 % CI | width |
> |---|---|---|---|---|---|---|
> | **raw** (as 3a ran it) | 0.446 | 0.409 | 2.010 | +1.601 | [−2.270, +5.472] | **7.742** |
> | **clip to [0,1]** | 0.446 | 0.386 | 0.358 | −0.028 | [−0.105, +0.050] | **0.155** |
> | **logit link** | 0.446 | 0.380 | 0.352 | −0.027 | [−0.109, +0.054] | **0.163** |
>
> Y1's state+Δt MAE goes **7.589 → 0.314**; the CI **narrows 50×** and the point estimate flips
> sign. All three still read *tied* — so the bound fixes the instrument and does **not** change
> this run's reading. Read against D0 that is the informative part: the signal is there, and a
> bounded *ridge on 2000 genes plus two Δt columns* still does not reach it. The failure is the
> **model class and the geometry**, not the corpus.
>
> **P4 failed and is recorded as failed:** the state-only arm is *not* broken on the graded
> geometry — it beats mean-only 0.409 vs 0.446, and 0.380 vs 0.446 once bounded. There was a
> working comparator. (Without the age mask the raw state arm *is* worse than mean-only, 0.564 vs
> 0.461, and only the bounded versions beat it.)
>
> ### D3 — Part B's identical swing, closed arithmetically
>
> SD of the swing across folds = **2.5e−14**; the analytic value `w_Δt·Δz(Δt) + w_Δt²·Δz(Δt²)` =
> **−269.12592**; **max |observed − analytic| = 5.7e−14**. The expression contains **no x₀ term**,
> so the five rows cannot disagree — one fit printed five times. **P3 confirmed**, and pinned as a
> theorem in `tests/test_stage3a_diagnose.py` rather than as an observation.
>
> Re-run correctly (one LODO fit per fold, each swept over its own Δt range) the swings finally
> differ — **−241.91, −247.74, −264.16, −267.33, −274.71**, SD 13.80 yr — and remain
> **nonphysical**. Part B never "passed" the >2 yr clause; it was out of range, exactly as Part C's
> 7.589 is.
>
> ### D4 — `bar_verdict` at the real geometry: the bar could not have registered the signal
>
> §5b's check, run at last, on the real features, the real Δt and the real fold/pair structure —
> only the **target** is simulated, so a Δt effect is present **by construction**. ρ is the share
> of the simulated signal carried by Δt (ρ = 1 is a *pure function of Δt*). 3a's rule is graded
> verbatim: PASS iff the paired 95 % CI upper end < 0. 2000 trials/cell, `MIN_PASS_RATE = 0.95`.
>
> | ρ | σ | pass rate, **raw** (what 3a used) | pass rate, **logit** |
> |---|---|---|---|
> | 0.00 | 0.05 / 0.001 | 0.000 / 0.000 | 0.000 / 0.000 |
> | 0.25 | 0.05 / 0.001 | 0.000 / 0.000 | 0.000 / 0.000 |
> | 0.50 | 0.05 / 0.001 | 0.000 / 0.000 | 0.000 / 0.000 |
> | 0.75 | 0.05 / 0.001 | 0.273 / **1.000** | 0.055 / 0.000 |
> | 1.00 | 0.05 / 0.001 | **0.000** / **0.000** | 0.447 / **1.000** |
>
> 1. **The raw estimator is not a detector.** Its pass rate is **non-monotone in the true effect**
>    — 0.000 at ρ = 0.5, 1.000 at ρ = 0.75, back to **0.000 at ρ = 1.0**, where the target *is* Δt
>    with near-zero noise. A statistic whose verdict is not a function of the effect it is meant to
>    measure cannot support a verdict in either direction.
> 2. **Even repaired, the bar needs ρ ≈ 1.** The logit version clears 0.95 only at ρ = 1.0 with
>    σ = 0.001; at ρ = 0.75 it reaches 0.055, at ρ = 1.0 / σ = 0.05 it reaches 0.447 (0.823 at six
>    folds). At ρ = 0 both are 0.000, so there is no false-positive problem: this is **specificity
>    without sensitivity**.
>
> At 5 folds × 55–66 pairs × 11–12 timepoints, this bar can only register a Δt effect that explains
> **essentially all** of the forward signal. Per the accepted item verbatim: *"If a correct system
> cannot clear the bar at that scale, then the dataset cannot answer the question, and THAT is the
> finding, not STOP."* D0 sharpens it — the dataset **does** carry the answer; **this bar** could
> not have registered it.
>
> ### Pre-registered expectations, graded (written before the run, in the script's docstring)
>
> | | expectation | held? |
> |---|---|---|
> | **P1** | Y1's Δt is nested inside the training range, not outside | ✅ **YES** — 0 pairs outside |
> | **P2** | the Δt block carries the divergence | ❌ **NO** — ratio 1.13; the **gene** block does |
> | **P3** | Part B's swing is identical and equals the analytic value | ✅ **YES** — SD 2.5e−14 |
> | **P4** | the state-only arm is also broken (worse than mean-only) | ❌ **NO** — 0.409 vs 0.446 |
> | **P5** | the bar is UNRESOLVABLE even at pure Δt | ✅ **YES** for raw (0.000) |
>
> Two of five failed. Both are recorded as failures, and every section above is written from what
> was measured rather than from what was expected.
>
> ### What this changes
>
> | | |
> |---|---|
> | 3a's STOP | **stays withdrawn**, and the ground has shifted: it is not only unsupported, it points the wrong way on the graded geometry |
> | *"the forward signal may well be absent"* (withdrawal banner) | **superseded.** D0 measures it present and strong on those five folds; annotated in place, not deleted |
> | the audit's Y1-extrapolation mechanism | ❌ **refuted** — identical Δt support, 0 pairs outside |
> | my own P2 (the Δt block) | ❌ **refuted** — the gene block, rotated ~60° by adding the Δt columns |
> | Part B | **must not be read at all** — constant by arithmetic before the fix, nonphysical after |
> | Part C's prediction | **must be bounded.** Raw least squares on a fraction is misspecified twice: unbounded output, and a target that is 90 % at 0 or 1 |
> | 3a's bar | **UNRESOLVABLE as written.** No 3a verdict may be taken again until a resolvable bar is registered under §5b |
> | 3b / 3c / 3d | **still unwritten.** D0 is an existence proof, not a model — nothing here licenses building them |
> | Stage 5 | **still not entered on this basis** |
>
> ### What a valid 3a re-run needs, before it is run
>
> 1. A **bounded** predictor (logit link; `clip` is the floor) — D2.
> 2. A **resolvable bar** registered under §5b at the real geometry, with `tests/test_bars_resolvable.py`
>    updated. D4 says the current rule is not one. Pooling timepoints is the mitigation test 18's
>    own cells-per-timepoint warning already named.
> 3. **N2 resolved** — decoupled from the age mask, and a stated rule for a donor whose target has
>    zero variance.
> 4. A **model class that can reach the ceiling D0 measures.** A ridge over 2000 genes plus two Δt
>    columns cannot; the oracle needs only `t_j`.
>
> **What is NOT claimed.** That the tool is buildable, that 3b–3d should be written, or that the
> oracle is a model — it needs `t_j` and other donors' measurements at `t_j`, so it is a **ceiling
> and an existence proof**, not a predictor for an unseen culture. And it is established on five
> donors, not six.

> ## ⛔ **3a's STOP IS WITHDRAWN — 2026-08-08, later the same day.** The run below is NOT VALID and must not be acted on.
>
> *The result box is left exactly as it was written. Nothing in it is deleted — it is the record of a verdict that should not have been taken.*
>
> **Three defects, each verified here before accepting the challenge that raised them:**
>
> 1. **Part C's decisive statistic is produced by a DIVERGING FIT, not by absent signal.** The target is `unsafe.mean()` — a FRACTION, bounded `[0,1]` by construction (`test18_forward_gate.py:85`). Y1 reports **MAE 7.589**, so the model is emitting values around **8** on a quantity that cannot exceed 1. Out-of-range output is not evidence about forward safety.
> 2. **Y1 drives the entire verdict.** With Y1: mean +1.601, 95 % CI **[−2.836, +6.038]**, width **7.742**. Without it: mean **+0.211**, CI **[−0.337, +0.758]**, width **1.095**. *(**CORRECTED 2026-08-12** — the two intervals as written are wrong. Both were recomputed by hand with the `T_CRIT` entry for one fold too few; the means and the with-Y1 width are right, the endpoints are not. Measured by `stage3a_diagnose.py`: with Y1 **[−2.270, +5.472]**, width 7.742 — matching the recorded run exactly; without Y1 **[−0.194, +0.615]**, width **0.810**. No conclusion drawn from them changes. Left standing above rather than edited.)* **Y1 contributes 89.5 % of the mean.** And Y1 is structurally different — **11 timepoints and 55 pairs** against 12 and 66 for every other fold.
> 3. **Part B's per-fold column is CONSTANT BY CONSTRUCTION, and my own explanation of it was wrong.** I flagged the identical −269.13 as suspicious but attributed it to *"the ridge extrapolating on Δt²"*. The real cause is at `test18_forward_gate.py:224-228`: Part B fits **ONE global model on all donors' pairs** and varies only `x0` per row, so the swing is `β·Δ(features)` with identical `β` and an identical Δt range — **identical for every donor by arithmetic.** It is one fit printed five times, not five fits agreeing.
>
> **And no `bar_verdict` was ever run for 3a** — `REF_GROUND_RULES.md` §5b requires showing a system that meets intent clears the bar ≥ 95 % of the time **before** the run. With n = 5 and a CI spanning 7.7 units on a `[0,1]` target, the honest verdict is **UNRESOLVABLE**, not *"tied"*.
>
> **Why this matters more than an ordinary correction:** STOP is the one **terminal** verdict in this project — *"do not write tool code; ship the scoring model; go to Stage 5."* My script's own text says *"a NEGATIVE result is decisive"*, and that rule holds only when the negative comes from the **data**. A negative produced by a diverging fit on one structurally-odd fold is decisive about nothing.
>
> **What is NOT withdrawn:** the qualitative direction. Excluding Y1, Δt still fails to help (mean +0.211, CI includes 0). The forward signal may well be absent — **this run did not establish it.** *(**SUPERSEDED 2026-08-12 — and in the direction that favours the tool, not the STOP.** §D0 of the diagnosis box above measures a model-free ceiling: on the very five folds 3a graded, a predictor using **only `t_j = t_i + Δt`** reaches unsafe-fraction MAE **0.157** against 3a's best arm at 0.352 and its state-only at 0.409, paired **−0.336, CI [−0.450, −0.221]**. A forward time signal is **present and strong** on that geometry. The sentence above is left standing as what was believed on 2026-08-08; it should no longer be relied on. The limit is in §D0 — add donor N2, which never becomes unsafe, and the six-fold CI includes 0.)* 3b/3c/3d stay unwritten pending a valid re-run, and Stage 5 is NOT entered on this basis.
>
> Diagnosis and re-run: `experiments/stage3a_diagnose.py`.
>
> ---
>
> ## 🚨 **3a RESULT — 2026-08-08: STOP.** Run on C-7-clean labels AND on the contaminated ones; they agree.
>
> *Additive. Nothing above is changed.* `experiments/stage3a_forward_gate.py` imported
> `test18_forward_gate.py` **unmodified** and called its own `main()` with `resolve_root`
> redirected. Full record: `experiments/DELTAAGE_LAB_NOTEBOOK.md`, "STAGE 3a RUN".
>
> **Run on two arms, because this branch is terminal.** §5.14 of `STAGE_1_5_6_SPARSE_CLOCK.md`
> proved arm A's `y_age` is inflated ~3× in five of six folds by one degenerate GEO column, and a
> terminal decision must not rest on that. `_c7` (C-7 gate on, degenerate control rejected, donor
> N2's ΔAge masked) is **operative**; `_armA` is context.
>
> | | `_c7` (operative) | `_armA` |
> |---|---|---|
> | Part A ΔAge MAE, state only | **20.67** | 25.53 |
> | Part A state + Δt | 86.77 | 90.54 |
> | Part C unsafe-frac MAE, state only | **0.409** | 0.386 |
> | Part C state + Δt | 2.010 | 0.983 |
> | Part C paired | +1.601, CI [−2.270, +5.472] | +0.598, CI [−0.821, +2.016] |
> | **verdict** | **STOP** | **STOP** |
>
> Cleaning the labels improved the state-only fit by **19 %** (25.53 → 20.67) and **did not change
> the verdict.**
>
> **Δt does not merely fail to help — it makes the fit worse**, by ~65 MAE in Part A and 0.6–1.6
> in Part C. Both CIs include zero, so both are formally *tied*, but the point estimate is
> consistently the wrong way. That is §0.4's own diagnosis: along one trajectory, time is
> **redundant with state**.
>
> ⚠️ **Part B must NOT be read as passing the ">2 yr sweep" clause.** Its swing is 255–269 yr,
> and in `_armA` the per-fold swing is **identical to the digit (−255.76) in all six folds** — a
> real per-fold response cannot be identical across six different held-out donors. That is the
> ridge extrapolating on Δt², not signal. **Part C is decisive and test 18 says so itself.**
>
> **Consequence, per the table above: 3b, 3c and 3d are NOT written.** This routes to *"ship the
> scoring model; go to Stage 5."* The **fate head is untouched** and works (ROC 0.983, PR-AUC
> 0.992) — that is what shipping the scoring model means.
>
> **What it does not decide:** that forward prediction is impossible in principle. Test 18's own
> caveat governs — *"a screen, not a proof … a NEGATIVE result is decisive"* — decisive **for this
> corpus**, at ~12 timepoints × 6 donors. The two options §3a names, **(a)** a ΔAge-trajectory
> readout with no safety recommendation and **(b)** data with more unsafe-cell variation (Stage 6),
> are scope decisions and are **recorded, not taken.**

**A STOP is a real result** — this dataset cannot support forward prediction, which is worth
knowing and worth reporting.

---

# SUB-STAGE 3b — The data layer

## 3b.1 New file: `src/cellfate/data/forward_pairs.py`

```python
"""Forward (t_i -> t_j) pairs for the stopping-time model.

Sampling is destructive - the same cell is never observed twice - so a cell-to-cell pairing does
not exist. The unit is a POPULATION SNAPSHOT: mean expression at t_i paired with the mean outcome
at t_j.

From 12 timepoints this yields up to 66 ordered pairs per donor (~396 across six donors), versus
the 6 examples a day-0 -> endpoint pairing would give.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.data import gather_split
from cellfate.common.constants import LOSS_IDX, DEATH_IDX

MIN_CELLS_PER_TP = 1          # raise to 3 if the data allows
REGIME = "holdout"


def _resolve_root(name: str) -> str:
    """Fold folders may sit at the repo root or under runs/ - check both."""
    from pathlib import Path
    for base in (".", "runs", ".."):
        if (Path(base) / name).exists():
            return str(Path(base) / name)
    return name


@dataclass(frozen=True)
class ForwardPair:
    donor: str
    x_i: np.ndarray          # (G,) population-mean expression at t_i
    t_i: float               # log time_h of the source timepoint
    dt: float                # log-gap to the target timepoint  <-- THE NEW SIGNAL
    y_age_j: float           # mean TRUE ΔAge at t_j        (age-valid cells only)
    y_frac_unsafe_j: float   # fraction loss/death at t_j   (all cells)
    n_i: int
    n_j: int


def _timepoint_rows(donor: str, regime: str = REGIME):
    """Per-(donor, timepoint) population summaries."""
    te = gather_split(ArtifactPaths.of(_resolve_root(f"cellfate_loocv_{donor}")), regime, "test")
    t_all = np.asarray(te.dose_time[:, 1], float)
    X_all = np.asarray(te.X, float)
    y_all = np.asarray(te.y_age, float)
    cls   = te.y_cls.astype(int)
    am    = te.mask.astype(bool)

    rows = []
    for tp in np.unique(np.round(t_all, 6)):
        sel = np.isclose(t_all, tp)
        if sel.sum() < MIN_CELLS_PER_TP:
            continue
        sel_age = sel & am                       # NOTE: two different masks - do not conflate
        unsafe = (cls[sel] == LOSS_IDX) | (cls[sel] == DEATH_IDX)
        rows.append({
            "t": float(tp),
            "x": X_all[sel].mean(0),             # ALL cells -> the state
            "y_age": float(y_all[sel_age].mean()) if sel_age.any() else np.nan,
            "frac_unsafe": float(unsafe.mean()),
            "n": int(sel.sum()),
            "n_age": int(sel_age.sum()),
        })
    return sorted(rows, key=lambda r: r["t"])


def build_forward_pairs(donors, regime: str = REGIME) -> list[ForwardPair]:
    out = []
    for d in donors:
        rows = _timepoint_rows(d, regime)
        if len(rows) < 3:
            continue
        for i, ri in enumerate(rows):
            for rj in rows[i + 1:]:              # ordered: t_j > t_i, guaranteed by the sort
                if not np.isfinite(rj["y_age"]):
                    continue                     # target has no age-valid cells
                out.append(ForwardPair(
                    donor=d, x_i=ri["x"], t_i=ri["t"], dt=rj["t"] - ri["t"],
                    y_age_j=rj["y_age"], y_frac_unsafe_j=rj["frac_unsafe"],
                    n_i=ri["n"], n_j=rj["n"]))
    return out


def pairs_to_arrays(pairs):
    """-> (X, dose_time, y_age, y_unsafe, donor).  dose_time[:,1] = log Δt, NOT absolute time."""
    X  = np.vstack([p.x_i for p in pairs]).astype(np.float32)
    dt = np.zeros((len(pairs), 2), np.float32)
    dt[:, 1] = [p.dt for p in pairs]             # column 1 REINTERPRETED - see §0.4
    return (X, dt,
            np.array([p.y_age_j for p in pairs], np.float32),
            np.array([p.y_frac_unsafe_j for p in pairs], np.float32),
            np.array([p.donor for p in pairs]))
```

## 3b.2 The leakage rule — the single most important constraint in this stage

> **Split by DONOR. Never split by pair.**

Pairs within a donor share timepoints: `(t1→t5)` and `(t2→t5)` share the **target**; `(t1→t5)` and
`(t1→t9)` share the **source**. A random pair split therefore puts the *same measurement* on both
sides and will produce a beautiful, meaningless result.

```python
# CORRECT
train = [p for p in pairs if p.donor != held_out]
test  = [p for p in pairs if p.donor == held_out]

# CATASTROPHIC - never do this
train, test = train_test_split(pairs, test_size=0.2)     # shares timepoints across the split
```

**Stage 4 re-verifies this independently.** It is the most likely explanation for a suspiciously
good result.

## 3b.3 Sanity checks the builder must print

```python
def audit(pairs):
    import collections
    per = collections.Counter(p.donor for p in pairs)
    dts = np.array([p.dt for p in pairs])
    assert (dts > 0).all(),            "ordering bug: non-positive Δt"
    assert dts.std() > 0,              "all Δt identical - Δt cannot be learned"
    assert all(v >= 3 for v in per.values()), f"donor with <3 pairs: {per}"
    print(f"pairs={len(pairs)} donors={len(per)} "
          f"Δt min/med/max = {dts.min():.2f}/{np.median(dts):.2f}/{dts.max():.2f}")
    print("per donor:", dict(per))
```

| Check | Failure means |
|---|---|
| `dt > 0` everywhere | ordering bug |
| `dt` spans a usable range | if all gaps are similar, Δt cannot be learned — **stop** |
| ≥3 pairs per donor | that donor cannot be held out |
| no pair crosses donors | leakage |

## 3b.4 The weighting decision — make it before coding

A pair built from 1 cell at `t_i` is far noisier than one from 20. Options: (a) unweighted,
(b) weight by `min(n_i, n_j)`, (c) drop pairs below a threshold.

**Recommendation: (c) with a low threshold, then unweighted.** Weighting adds a hyperparameter this
dataset is too small to tune honestly — and this project has repeatedly found added flexibility
hurts at n≈100.

## 3b.5 Acceptance for 3b

- [ ] `build_forward_pairs(DONORS)` returns a non-empty list
- [ ] `audit()` passes every assertion
- [ ] Δt range is reported and non-degenerate
- [ ] **no training code is written until this holds**

---

# SUB-STAGE 3c — Training

## 3c.1 New file: `src/cellfate/training/train_forward.py`

Mirrors `train_model.py` with **three** differences. Everything else — `CellFateNet`, encoders,
heads, scalers, the bundle writer — is **reused unchanged**.

| # | Difference | Why |
|---|---|---|
| 1 | input is 3b pairs; `dose_time[:,1] = dt` | the forward signal |
| 2 | targets are `y_age_j` and `y_frac_unsafe_j` | population outcomes |
| 3 | calibration via Stage 1's `crossdonor_stats` | never in-distribution |

## 3c.2 The ordering rule

```
   train model on pairs
        ↓
   inner-LODO over training donors  ->  pooled cross-donor residuals / logits / features
        ↓
   fit temperature, q, sigma_scale, OOD on THOSE
        ↓
   write bundle   (bundle_meta.mode = "forward")
```

> **Calibration is fitted last, to the model that will ship.** Reusing the existing bundle's
> calibration, or fitting before the final training run, reproduces exactly the defect Test 14
> measured (coverage 0.40 vs 0.90).

## 3c.3 The head decision

The existing fate head is 3-class trained on per-cell labels. 3b targets are **population
fractions**.

| Option | Implication |
|---|---|
| (a) regress the unsafe fraction — single sigmoid | simplest; **loses the loss-vs-death distinction** |
| **(b) keep 3-class, train on soft labels** (class proportions) | preserves it; needs soft-target cross-entropy |

**Recommendation: (b).** Losing identity and dying are different failures with different
consequences, and soft-target cross-entropy is a two-line change:

```python
# hard labels:  F.cross_entropy(logits, targets.argmax(1))
# soft labels:  -(targets * F.log_softmax(logits, dim=1)).sum(1).mean()
```

## 3c.4 What must NOT change

- `N_DOSE_TIME = 2` — hard contract with a validator that raises
- the existing scoring bundle — untouched, still measured by `scorecard.py`
- `CellFateNet` — same architecture. This project's own evidence (T3, T6, T9) is that added
  capacity **hurts** at this scale; what changes is what it trains on, not how big it is

## 3c.5 Expected effects — record before running

| Metric | Now | After 3c | Basis |
|---|---|---|---|
| **`dt_response`** | **0.035 yr** | **> 2 yr** | the 3a gate threshold |
| forward coverage | — | 0.85–0.95 | Stage 1's method |
| forward ΔAge MAE | — | comparable to per-cell MAE | — |

> **If `dt_response` stays near zero after training on pairs, 3c has FAILED** even if the losses
> look fine — the model has again learned to read time off state. **Check this before 3d.**

---

# SUB-STAGE 3d — The decision layer

## 3d.1 New file: `src/cellfate/inference/stopping.py`

```python
"""Stopping-time recommendations from the forward bundle.

Answers: "given my culture as it is now, when should I withdraw?" -- a RELATIVE question along
the one axis this data varies, which is why it survives the ±12.7 yr per-donor level shift.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Sequence
import numpy as np

from cellfate.inference import Predictor


@dataclass
class StoppingOption:
    dt_hours: float
    withdraw_day: float
    delta_age: float
    delta_age_lo: float      # calibrated, POPULATION-level
    delta_age_hi: float
    p_unsafe: float          # calibrated
    in_range: bool           # False => Δt outside the trained range


@dataclass
class StoppingReport:
    options: list[StoppingOption]
    recommendation: StoppingOption | None
    reason: str
    warnings: list[str] = field(default_factory=list)


def recommend_stopping(X_now, bundle_root, risk_threshold=0.20,
                       dt_grid_hours: Sequence[float] | None = None,
                       dt_range=None, donor_offset=None) -> StoppingReport:
    pred = Predictor(bundle_root)

    # --- 1. ASSERT the right bundle. Loading the scoring bundle here silently produces
    #        nonsense, because column 1 means absolute time there and Δt here.
    #  AUDIT FIX: the attribute is `pred.meta` (predictor.py:67), NOT `pred.bundle_meta`.
    #  Using getattr(pred, "bundle_meta", {}) would silently yield {} -> mode=None -> the
    #  assert would fire on EVERY bundle, blocking the tool entirely. ---
    mode = getattr(pred.meta, "mode", "scoring")
    if mode != "forward":
        raise ValueError(f"stopping requires a FORWARD bundle, got mode={mode!r}")

    X_now = np.atleast_2d(np.asarray(X_now, np.float32))
    n_cells = X_now.shape[0]
    warnings: list[str] = []

    # --- 2. population mean: the tool scores CONDITIONS, not cells (§0.3) ---
    x_pop = X_now.mean(0, keepdims=True)

    grid = list(dt_grid_hours) if dt_grid_hours is not None else _default_grid(dt_range)
    options = []
    for dt in grid:
        dose_time = np.array([[0.0, np.log(dt)]], np.float32)   # col 1 = log Δt
        # AUDIT FIX: Predictor has no `.arch`; the architecture lives on the members.
        # NOTE: arch holds BOTH "n_fp" (fingerprint width) and "n_pert" (actual input width;
        # network.py:43 sets n_pert = len(TF_VOCAB) for tf-kind, else n_fp). Use n_pert -
        # "n_fp" is wrong for TF cocktails, which is what this dataset uses.
        fp = np.zeros((1, pred.members[0].arch["n_pert"]), np.float32)
        row = pred.predict_encoded(x_pop, fp, dose_time)[0]

        mu = float(row["mu_age"])
        if donor_offset is not None and donor_offset.applied:
            mu -= donor_offset.d                                # Stage 2, optional

        # --- 3. population interval: sigma/sqrt(n) is a LOWER BOUND (cells are correlated) ---
        half = pred.q / np.sqrt(max(n_cells, 1))
        options.append(StoppingOption(
            dt_hours=dt, withdraw_day=_to_day(dt),
            delta_age=mu, delta_age_lo=mu - half, delta_age_hi=mu + half,
            p_unsafe=float(1.0 - row["S"]),
            in_range=(dt_range is None or dt_range[0] <= np.log(dt) <= dt_range[1]),
        ))

    # --- 4. recommendation rule ---
    ok = [o for o in options if o.p_unsafe <= risk_threshold and o.in_range]
    rec = min(ok, key=lambda o: o.delta_age) if ok else None
    reason = ("most rejuvenation with P(unsafe) below threshold" if rec
              else f"no withdrawal time meets P(unsafe) <= {risk_threshold:.0%}")

    # --- 5. mandatory warnings ---
    warnings.append("intervals are population-level and assume independent cells "
                    "(LOWER BOUND on true uncertainty)")
    if donor_offset is None or not donor_offset.applied:
        warnings.append("absolute ΔAge carries a ±12.7 yr per-donor offset; comparisons WITHIN "
                        "this donor are reliable, absolute values are not")
    if n_cells < 10:
        warnings.append(f"only {n_cells} cells - uncertainty estimate unreliable below ~10")
    if any(not o.in_range for o in options):
        warnings.append("some options are outside the trained Δt range and were not recommended")
    if rec is None:
        warnings.append("NO SAFE WITHDRAWAL WINDOW FOUND")

> ## 🆕 ADDED 2026-07-31 — the recommendation rule sorts on a label Stage 1.5.2 could not validate
>
> *Additive; §4's rule and the `StoppingOption` dataclass above are unmodified.*
>
> `rec = min(ok, key=lambda o: o.delta_age)` makes **ΔAge the sole ordering key** of this tool's
> output. `STAGE_1_5_2_LABEL_ANCHOR.md` returned **NOT CALIBRATABLE**: the transcriptomic clock does
> not track methylation age (ρ_partial +0.267 / +0.516 against a pre-frozen 0.50 bar), and §12-R
> confirmed the anchor itself is sound, so the failure is the RNA clock's.
>
> **What this does and does not invalidate:**
>
> | | |
> |---|---|
> | ❌ **absolute** ΔAge as a quantity to compare across donors | the number is not validated in years |
> | ✅ **within-donor ordering** of withdrawal times | `rank_model_dage` is **0.91–0.99** across all six folds (`scorecard/baseline.json`), and 1.5.2 §17 found the RNA clock reaches 91% of the meth↔meth ceiling against one of the two references |
>
> **So the rule survives as a *ranking* rule and fails as a *reporting* rule** — which is precisely
> the distinction §0.3's internal-control note above already asks this stage to make. The two should
> be resolved together, not separately.
>
> **Concretely, for 3d:** `StoppingOption` (line 385) needs a validity field beside `delta_age`, and
> `_to_day`'s output should be presented as an ordered shortlist rather than a ΔAge in years, unless
> and until the label basis changes. `STAGE_1_5_3_EXECUTE.md` C-4 costs the equivalent change at the
> inference boundary; do not solve it twice with two conventions.

    return StoppingReport(options, rec, reason, warnings)
```

## 3d.2 Refuse-to-extrapolate

`in_range = dt_min_observed <= log(dt) <= dt_max_observed`, taken from 3b.

Out-of-range options are **computed and shown but never recommended**, and flagged. The model has
no information beyond the observed gaps, and **silently extrapolating a safety prediction is the
single most dangerous thing this tool could do.**

## 3d.3 The √n caveat is not decoration

`sigma/√n` assumes independent cells. **Cells from one culture are correlated**, so this is a
**lower bound** on true uncertainty. It is stated in the report, printed to the user, and repeated
in the manuscript limitations.

## 3d.4 Mandatory warnings

| Condition | Warning |
|---|---|
| always | population intervals are a lower bound |
| no donor offset applied | absolute ΔAge carries ±12.7 yr; **within-donor comparisons are the reliable output** |
| `n_cells < 10` | uncertainty unreliable |
| any option out of range | not recommended, shown only |
| OOD flag fires | "this culture is unlike the training trajectories" |
| no option meets threshold | **no safe window found** |

**`recommendation = None` is a correct, useful answer** — it means this culture has no safe
window, which is exactly what a researcher needs to hear when true.

## 3d.5 Output the researcher sees

```
Culture: donor D, sampled day 10, 847 cells

  withdraw day 12   ΔAge  −6 [−10, −2]   P(unsafe)  8%
  withdraw day 15   ΔAge  −9 [−13, −5]   P(unsafe) 14%
  withdraw day 18   ΔAge −11 [−15, −7]   P(unsafe) 29%   ← above threshold
  withdraw day 21   ΔAge −12 [−16, −8]   P(unsafe) 47%   ← above threshold

  RECOMMENDATION: withdraw day 15
    most rejuvenation with P(unsafe) below your 20% threshold

  ⚠ intervals are population-level, assume independent cells (lower bound)
  ⚠ ΔAge comparisons within this donor are reliable; absolute values carry ±12.7 yr
```

## 3d.6 Training/serving skew — the silent failure

The live sample **must** be preprocessed exactly as 3b did: same gene panel, same scaler, same
population-mean step. **A mismatch here will not announce itself** — it produces confident wrong
answers.

**AUDIT FIX:** `Predictor` has no `.panel`. The panel identity lives in the metadata, and the
Predictor already cross-checks it against the scalers at construction (`predictor.py:74`). Use:

```python
# gene count: the scalers carry the panel width
assert X_now.shape[1] == len(pred.scalers.params.x_mean), \
    f"expected {len(pred.scalers.params.x_mean)} genes, got {X_now.shape[1]}"
# panel identity: meta vs scalers is already asserted inside Predictor.__init__,
# so a successful load already guarantees they agree. Re-assert only against YOUR data:
assert forward_pairs_panel_hash == pred.meta.gene_panel_hash, \
    "the pairs were built on a different gene panel than the bundle was trained on"
```

## 3d.7 Explicit non-goals

- **No per-cell recommendations** — per-cell uncertainty exceeds the effect size
- **No dose recommendations** — identifiability wall, one dose in the data
- **No cross-donor absolute claims** without reference cells

## 3d.8 Failure modes

| Symptom | Cause | Action |
|---|---|---|
| `ValueError: requires a FORWARD bundle` | wrong bundle loaded | correct — this assert exists to prevent silent nonsense |
| every option flagged out-of-range | `dt_range` not passed from 3b | pass it; do not disable the guard |
| `recommendation` always `None` | risk threshold too strict, or fate uncalibrated | check Stage 1 acceptance passed |
| ΔAge identical across all Δt | `dt_response` ≈ 0 → **3c failed** | return to 3c; do not ship |
| intervals absurdly narrow | `q` from the scoring bundle, not the forward one | check `bundle_meta.mode` |

## 3d.9 Acceptance for Stage 3

```powershell
python scorecard.py snapshot --tag C_forward
python scorecard.py compare B_percalib C_forward
```

| Role | Metric | Bar |
|---|---|---|
| **TARGET** | `dt_response` | **> 2 yr** (from 0.035) |
| **TARGET** | `forward_coverage` | 0.85–0.95 |
| **GUARD** | the existing scoring metrics | must be **noise** — the forward bundle is a sibling, not a replacement |

## 3d.10 Done when

`recommend_stopping` runs on a held-out donor's day-N sample, returns a report whose intervals are
calibrated, and every §3d.4 warning fires when it should.
