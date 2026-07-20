# REFERENCE — Program architecture

**Read after `00_START_HERE.md`.** That file tells you what to type; this one tells you why the
stages are ordered as they are, and what each one buys.

Covers `MASTER_PLAN.md` §1, §2, §3, §5d, §14, §15.

---

## 1. Component status — where the project stands

**Four categories, not a broken/working binary.**

| Component | Status | Evidence |
|---|---|---|
| **Within-donor ranking** | ✅ **WORKS** | Spearman **0.925–0.983**, every fold (T7, T7.4.3) |
| **ΔAge vs baselines** | ✅ **WORKS (at optimum)** | ties ridge (T5); nothing beats ridge (T6); architecture proven capable (T3) |
| **Fate discrimination** | ✅ **WORKS** | PR-AUC 0.929–0.940 in-dist, **0.96–1.00 out-of-donor** (T8.1); holds on Y1 (0.961 vs logreg 0.636) |
| **Perturbation channel** | ✅ **CARRIES SIGNAL** | `u_only` fate PR-AUC **0.83–1.00**; on Y1 `u_only` **0.831 beats** `x_only` **0.639** (T11) |
| **Harmonization** | ✅ **WORKS** | parameter-free, leak-safe, unit-tested |
| **ΔAge absolute level** | 🔧 **FIXABLE — demonstrated** | ±12.7 yr per-donor shift (T7.4.3); k=3 cells → MAE 14.3→7.1 (T16) |
| **Conformal intervals** | 🔧 **FIXABLE — cause known** | coverage **0.40 vs 0.90** (T14) |
| **Fate calibration** | 🔧 **FIXABLE — demonstrated** | ECE 0.28 → ~0.13 with Platt (T8.2) |
| **OOD gate** | 🔧 **FIXABLE or DISABLE** | AUC **0.47** ≈ chance (T15) |
| **RES score** | ⏸️ **VERDICT DEFERRED** | every input it consumes is defective; never tested on corrected inputs |
| **Cross-perturbation thesis** | 🚫 **DATA-LIMITED** | one cocktail (T11, T11.1). No code change helps |

## 2. What the ranking result actually says

**This must not be misstated, and an earlier draft of the plan did misstate it.**

| What was ranked | Spearman vs true ΔAge | Reading |
|---|---|---|
| model ΔAge | **0.948** | ranking capability: excellent |
| ridge ΔAge | **0.955** | ranking capability: excellent |
| RES score | 0.686 | **the RES formula degrades a good ranking** |

**The model ranks correctly.** The RES transform was throwing that away — because it multiplies a
good ΔAge by three defective signals. **"Ranking is broken" was never true and must not appear in
any writeup.**

## 3. Dependency graph

```
                 STAGE 1 — CALIBRATION  (Change A)
                 fixes temperature, q, sigma_age, OOD
                 ┌──────────────┴──────────────┐
                 │                             │
                 │            STAGE 3a — GATE (Test 18)
                 │                 │  STOP if no forward Δt signal
                 │                 ▼
     STAGE 2 — LEVEL CORRECTION   STAGE 3b–3d — THE TOOL
     (Change B + A′)               │  data → training → decision layer
     absolute ΔAge only            │  = condition-level RES
                 │                 │
                 └──────┬──────────┘
                        ▼
              STAGE 4 — VALIDATION  (incl. Change C, the RES verdict)
                        ▼
              STAGE 5 — PUBLICATION  (second clock, Y1, writeup)
```

**Stage 1 is the foundation for everything, including the tool.** The tool needs a calibrated
`P(unsafe)` for its risk threshold, honest `q` for its error bars, corrected `sigma_age` for
`R_eff`, and a working OOD flag for its warning. **All four are Change A.** Building the tool first
would ship a recommender whose stated 90% intervals cover 40%.

## 4. Ordering — and why it is A → B → A′

Stage 2 changes the error scale that Stage 1's `q` depends on:

| step | mean \|error\| | required q (P90 ≈ 2.5–3.0× mean) |
|---|---|---|
| now (in-distribution calibration) | 14.3 | **8.9** — far too small |
| after **A** | 14.3 | ~35–43 |
| after **B** | 6.9 | ~17–21 |
| after **A′** (refit) | 6.9 | **~17–21** ✓ |

**Why not B first?** B alone leaves `sigma` understated while `mu` is corrected — the coupling that
makes the safety score ~2× more permissive. **A first is the safe direction**: over-wide intervals
make RES *conservative*, never permissive.

## 5. What each stage buys

**Evidence classes:** **[M]** measured · **[A]** arithmetic on measured quantities · **[U]** unknown.

### Stage 1 — calibration

| Gain | From | To | Class |
|---|---|---|---|
| conformal coverage | 0.40 (0.00 on N2/N3) | ~0.90 | [A] |
| fate ECE | 0.28 | ~0.13 | **[M]** T8.2 |
| cells wrongly discarded by OOD | 27% (48% on N2) | 0% or justified | **[M]** T15 |

**Buys TRUSTWORTHINESS, not accuracy.** MAE stays 14.3, ranking stays 0.948. But for a safety
decision, **a miscalibrated confidence is worse than no confidence** — the intervals currently
cover *zero* cells on two donors while claiming 90%. **Cost: pure code, no new data.**

### Stage 2 — level correction

| Gain | From | To | Class |
|---|---|---|---|
| ΔAge MAE | 14.3 | **~7.0 (−50%)** | **[M]** T16 |
| N3 MAE | 29.7 | 10.0 | **[M]** |
| interval width (after A′) | ~35–43 yr | ~17–21 yr | [A] |

**Buys usable absolute ΔAge out-of-donor** — the largest measured gain. Does **not** change ranking
(rank-invariant by construction). **Cost: k≈3 reference cells per donor — a protocol change** —
and it *hurts* O1 and Y1 unless applied conditionally.

### Stage 3 — the tool

| Gain | From | To | Class |
|---|---|---|---|
| uncertainty on the scored quantity | ±17–21 yr (per cell) | **±3.7–4.6 yr** (mean of 21) | [A] |
| detecting an ~11 yr effect | **impossible** | **comfortable** | [A] |

**Buys the ability to ask the question at all.** Per-cell confident rejuvenation is arithmetically
out of reach; averaging over a condition puts the effect in range. **Cost: a genuine build.**

### Stages 4 and 5

Buy **verdicts instead of guesses** — the deferred RES decision, and whether the linearity and fate
claims survive. **[U]** by nature; that is why they must be run rather than assumed.

## 6. The honest counterfactual — what doing nothing already gets you

- **Within-donor ranking, Spearman 0.925–0.983** — publishable today
- **Fate discrimination out-of-donor, 0.96–1.00** — publishable (with the Y1 caveat)
- **A parameter-free, leak-safe harmonization method** — publishable
- **A rigorous null result with a structural explanation** — publishable

> **The stages make the system more useful and more honest. They are not what makes it valid.**
> Stage 1 is closest to mandatory: shipping intervals that cover 0% while claiming 90% is a
> correctness problem, not a polish problem.

## 7. Program-level risks

| Risk | Stage | Mitigation |
|---|---|---|
| Building the tool on broken calibration | 3 | Stage 1 is a hard prerequisite, not a parallel track |
| Leakage via pair splitting | 3b | split by donor; Stage 4 re-verifies independently |
| Calibration fitted to a superseded model | 3c | fit last, on the deployed model |
| The two bundles confused at inference | 3d | `bundle_meta.mode` assert |
| Suspiciously good validation | 4 | 6/6 triggers a leakage audit, not celebration |
| Linearity claim rests on a clock we fitted ourselves | 5 | second-clock gate before publishing |
| Fate claim rests on one fold (Y1) | 5 | Y1 gate before publishing |
| F2 hurts already-calibrated donors | 2 | conditional rule; report per-fold |
| **Over-correcting into "the model is worthless"** | all | §1 and §2 of this document exist to prevent it |

## 8. Three honest exits

- **Stage 3a STOP** → no forward signal. Ship the scoring model; skip to Stage 5.
- **Stage 4 fails** → honest but not useful. Ship a calibrated readout, not a recommender.
- **Stage 5 Gate A overturns linearity** → the null result changes, and the paper gets *more*
  interesting.

**Each exit still leaves a publishable result. That is by design.**
