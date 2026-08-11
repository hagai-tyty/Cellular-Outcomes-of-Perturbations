# STAGE 3a — the STOP verdict is NOT SAFE TO ACT ON

**Status:** 🔴 **AUDIT 2026-08-11. The gate RAN and returned STOP. The STOP is produced by a
diverging fit, not by an absence of signal. `src/` untouched by this audit; nothing withdrawn.**

**Why this matters more than a normal disagreement:** STOP is **terminal** — *"do not write tool
code. Ship the scoring model; go to Stage 5."* It ends the forward-prediction arm of the project.

---

## 1. First, what DID get fixed — this is real and it is large

C-7 was implemented and adopted on the dataset. The labels measurably improved:

| fold | day-14 ΔAge **before** (arm A) | **after C-7** |
|---|---:|---:|
| N2 | −7.352 | −7.337 |
| N3 | −22.121 | **−5.628** |
| O1 | −24.023 | **−6.514** |
| O2 | −22.891 | **−6.674** |
| Y1 | −22.049 | **−4.606** |
| Y2 | −23.869 | **−8.292** |
| **fold spread** | **16.671 yr** | **3.686 yr** |

**The fold instability fell 78 %**, and the magnitudes collapsed from ~−22 to ~−6.5 — landing on
§5.14's pre-registered prediction of **−8.196**. `n_cells` 42 605 → **42 600** (exactly the five
rejected columns) and N2's **19** cells are masked in **every** fold (rule 4 firing as specified).

**This is a genuine fix, not a reclassification.** It also found a real bug on the way — the gate
flag never reached injected sources, so the gate was **inert** when first wired (`e6fc183`).

---

## 2. 🔴 But Stage 3a's STOP rests on a model that is diverging

### The decisive fact: an impossible number

Part C's target is built at `test18_forward_gate.py:85`:

```python
unsafe = ((cls == LOSS_IDX) | (cls == DEATH_IDX)).astype(float)   # then .mean() per timepoint
```

**`p_unsafe` is a FRACTION. It cannot leave [0, 1].** Yet:

| held-out | pairs | state only | state + Δt |
|---|---:|---:|---:|
| N3 | 66 | 0.437 | 0.997 |
| O1 | 66 | 0.411 | 0.362 |
| O2 | 66 | 0.408 | 0.592 |
| **Y1** | **55** | **0.424** | **7.589** |
| Y2 | 66 | 0.362 | 0.509 |

> **An MAE of 7.589 on a quantity bounded by 1 means the model is predicting values around 8.**
> That is not "no forward safety signal". That is a regression producing out-of-range output.

### It is not isolated

| symptom | reading |
|---|---|
| Part A, Y1: ΔAge MAE **22.84 → 311.47** by adding one feature | a standardized ridge at `alpha=1.0` cannot degrade 14× on signal grounds |
| Part B: swing **−269.13 yr for ALL FIVE folds, identical to 2 dp** | five independent leave-one-donor-out fits cannot agree to 2 dp by chance. The Δt coefficient is not being learned from each fold's data |
| Part B magnitude: **269 yr** of aging swing | against a `> 2 yr` threshold. It does not "pass" — it is nonphysical |
| the paired CI: mean **+1.601**, 95 % **[−2.270, +5.472]**, n = 5 | **driven entirely by Y1.** Without it: +0.049, −0.184, −0.560, −0.147 — a completely different picture |

**Y1 is the fold with 11 timepoints, not 12** (no d29) and 55 pairs, not 66. Held out, its Δt
distribution sits outside the training folds', the scaler maps it to extreme z-scores, and the
linear model extrapolates. **That is a data-shape artefact, not a property of forward prediction.**

### And the bar was never shown resolvable

`REF_GROUND_RULES.md` §5b requires simulating a system that meets intent **before** the run and
showing it clears the bar at ≥ 95 %. **There is no `bar_verdict` for Stage 3a** in
`results/stage3a_forward_gate_results.json`. With **n = 5** and a CI spanning **7.7 units** on a
[0,1] target, the honest verdict is **UNRESOLVABLE**, not "tied".

> The script's own caveat says *"A NEGATIVE result is decisive."* **That holds only when the
> negative comes from the data.** A negative produced by a diverging fit is not decisive about
> anything.

---

## 3. What is actually established, and what is not

| | |
|---|---|
| ΔAge labels are materially better | ✅ **established** — 78 % less fold spread, contamination removed |
| the C-7 gate works and is adopted on the dataset | ✅ established — A1–A4 pass |
| **RES improved** | ❌ **UNKNOWN — not measured.** RES is a *model* output and **no retrain has happened.** The newest scorecard is `gc2_D_stratshuffle_hff_s0.json` (2026-08-07), which **predates C-7** |
| every Stage 1 guard under the new labels | ❌ **not re-reported.** C-7 §5 requires it on adoption |
| **there is no forward Δt signal** | ❌ **NOT established.** The fit diverges on the fold that dominates the CI |

**"Stage 1.5.6 closed" is fair. "Ready for Stage 5" is not.**

---

## 4. What to do instead — cheap, and it does not re-open anything

1. **Diagnose Y1.** Print the held-out Δt range against the training range per fold, and the fitted
   Δt coefficient. If Y1's Δt lies outside the training support, the fold is being extrapolated and
   must be reported as such rather than averaged in.
2. **Bound the Part C prediction to [0, 1].** The target is a fraction by construction; a predictor
   that leaves the range is misspecified. `clip` is the minimum; a logit link is the honest fix.
3. **Explain the identical −269.13 swing** before Part B is read at all.
4. **Run `bar_verdict` for 3a** — simulate a system with a real Δt effect at this geometry (5 folds,
   55–66 pairs, 11–12 timepoints) and confirm it clears the bar at ≥ 95 %. **If it cannot, the
   dataset cannot answer the question and that — not STOP — is the finding.**

**None of this requires a retrain and none of it re-opens 1.5.6.** It is the same discipline applied
to 1.5.6's own bars, applied to the one verdict in this project that is irreversible.
