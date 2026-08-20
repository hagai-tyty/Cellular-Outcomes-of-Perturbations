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

---

# STAGE 3a — the diagnosis is right; "the blocker is more data" is NOT

**Status:** 🟡 **AUDIT-2, 2026-08-12.** The STOP withdrawal and the resolvability study are **correct
and I agree with them.** The *conclusion drawn from them* is wrong: **the dataset the gate needs is
already on disk.** `src/` untouched, nothing withdrawn.

---

## 1. What is right, including one place I was wrong

| claim | verdict |
|---|---|
| Stage 3a's STOP came from a diverging fit, not the data | ✅ **confirmed** — `partC_frac_pred_outside_unit = 0.773`. **77 % of Part C's predictions leave [0, 1]** on a target that is a fraction by construction |
| the gate was graded on ~0.3 % of the cells | ✅ Gill's ~100 bulk columns out of 42 600 |
| **regime A — the geometry 3a actually ran on — is UNRESOLVABLE** | ✅ **at every alpha, raw and logit; pass rates 0.000–0.017 against `min_pass` 0.95.** A correct system could never have returned GO. This is the finding, and it is a strong one |
| the estimator, not the data, produced the STOP | ✅ 2000-trial null, 4 regimes, 4 alphas, both links — proper §5b work |

### ❌ And one thing I got wrong

My audit blamed **Δt extrapolation** on Y1 (11 timepoints, not 12). **That is refuted by measurement:**

```
heldout_pairs_outside_train_dt_support = 0
z_dt_heldout_absmax = z_dt_train_absmax = 2.4427      # identical
```

The extrapolation is in the **gene block**, not Δt: `z_gene_heldout_absmax = 39.45` against
`z_gene_train_absmax = 7.24`. **My mechanism was wrong; the conclusion it supported — that the fit
diverges — was right for a different reason.** Recorded rather than quietly dropped.

---

## 2. 🔴 Where I disagree: the blocking dataset is on disk and has never been tried

The diagnosis says the binding constraint is per-timepoint replication and a line/modality shift.
**Both are true of `gill_bulk`. Neither is true of `GSE165177`**, which sits in
`C:/Users/hagay/Desktop/GSE165177`, is **in no training config**, and has been item #2 on the
standing shortlist — *"95 adult in-range methylation-paired samples, in no training config, free"* —
since before this arc began.

### Measured from the raw matrices, just now

| | `gill_bulk` (Sendai) — what 3a ran on | **`GSE165177` (transient)** |
|---|---|---|
| donors | 6 | **3** (O1, O2, O3) |
| timepoints | 12 | **4** (10, 13, 15, 17 d) |
| **samples per (donor, timepoint)** | **≈ 1.7** | **7.0 – 8.0** |
| **controls** | **1 per donor, day 0 only — 6 total** | **2–3 per donor PER TIMEPOINT — 33 total** |
| condition arms | reprogramming time course | **6**, incl. an explicit `negative_control` arm |

### Why the control column is the one that matters

**Every failure this arc has fought traces to `gill_bulk` having one unreplicated day-0 control per
donor:** G-a's "`n = 1` visible" gate; C-7 and `N2_Fib` (one bad column corrupting five of six
folds); Group D's silent self-centring; Group E's per-chunk control absence; §4.7's 16.67 yr fold
instability (`σ_ref` from five single samples); and rule 4 / B2′ existing at all.

> **`GSE165177` has replicated, contemporaneous controls at every timepoint. Not one of those
> failure modes can arise in it.**

### And it directly fixes the mechanism the diagnosis names

`p_unsafe` is a fraction estimated per (donor, timepoint). At **1.7** samples it can only take
**{0, 0.5, 1}** — which is why `saturated_at_bounds = 63` of 70 profile values. **At 7–8 samples it
moves in eighths.** That is the difference between a regressable target and a saturated one.

---

## 3. The honest cost — this is a trade, not a free upgrade

**It is not strictly better, and the write-up must not say it is:**

| | `gill_bulk` | `GSE165177` |
|---|---:|---:|
| LOO folds | **6** | **3** |
| ordered forward pairs per donor | **66** (12 tp) | **6** (4 tp) |
| total pairs | **319** | **≈ 18** |

**Six pairs per donor against sixty-six is a large loss**, and 3 folds widen every CI. The gain is
4–5× the per-timepoint replication, real contemporaneous controls, and no modality shift — 319
*noisy, highly-correlated* pairs against ~18 *well-estimated* ones.

**Nobody has evaluated that trade, and the diagnosis does not mention it.** "More data is the
blocker" skips the question rather than answering it.

---

## 4. What to do — free, and the machinery already exists

**Add `GSE165177` as regime E to `experiments/stage3a_bis_resolvability.py`** — the same 2000-trial
null, the same alphas, raw and logit — at its real geometry: 3 folds, ~6 pairs/donor, 7–8 samples
per timepoint, controls at every timepoint.

| outcome | reading |
|---|---|
| **regime E RESOLVABLE** | the gate can be run for real on data we already hold. **No acquisition, no retrain** |
| **regime E UNRESOLVABLE** | *then* the acquisition claim is established — and it will be established with a number instead of asserted |

Either way the answer costs one script run and settles it. **A second option worth costing in the
same pass:** with 6 condition arms, forward pairs can also be built **across arms at matched
timepoints** (`negative_control` → `failed` → `transiently_reprogrammed`) rather than only
`t_i → t_j`. That is a different and possibly richer design than the one 3a assumed, and it does not
depend on the timepoint count at all.

**Do not go to Stage 5, and do not open Stage 6 acquisition, until regime E has been run.**
