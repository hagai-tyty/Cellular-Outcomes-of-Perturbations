# P3 — PRE-REGISTRATION. Does molecular progress beat calendar day? Committed BEFORE the run.

**Scope: HFF only** (`GSE242423`, 9 timepoints, ~42k cells). One line. No transfer test here —
that is P5, and it is gated on this.

**The idea under test.** Two lines can run the same course at different speeds, so *calendar day*
cannot transfer while *trajectory position* might. This also reinterprets an existing result: Test
11.1 and Test 18 found time redundant with state along one trajectory and filed it as failure;
under this model that is the expected observation.

---

## 1. Why this is not a trivial comparison in either direction

Within one line, day and progress look redundant — but they differ in one way that matters:

> **All cells at day 8 share a day. They do not share a molecular progress.** Reprogramming is
> asynchronous, so at any timepoint some cells are far along and others have barely started.
> Progress can explain **within-timepoint** variation in risk that day, being constant across the
> timepoint, cannot.

That also fixes the metric. **Evaluated per cell inside a held-out timepoint, `day` predicts a
constant** — so AUC would hand progress an automatic win, and a per-timepoint metric would hand
day one. Neither is informative.

**Primary metric: LOG-LOSS on held-out cells** — a proper scoring rule. A constant predictor that
names the held-out timepoint's true risk exactly scores *well*. Progress only wins if per-cell
variation is real **and** predictable. Brier score reported alongside.

---

## 2. Design — leakage-proofed

**Leave-one-timepoint-out**, 9 folds. For each held-out timepoint, **everything** is fit on the
training timepoints only:

1. **Scaling and dimensionality reduction** — fit on training cells.
2. **The progress coordinate** — a ridge from training-cell expression onto training day.
   `progress(cell) = predicted day`. Continuous and per-cell even though trained on discrete days,
   because every cell's expression differs.
3. **The risk model** — fit on training cells, from whichever predictor the arm uses.

Held-out cells are then **projected into that frozen coordinate** and scored. Nothing derived from
the held-out timepoint — not the scaler, not the reduction, not the progress fit, not the risk
model — touches any of it. *If any of them saw the held-out day, this is not out-of-time
validation and the run is void.*

### Arms

| arm | predictor | note |
|---|---|---|
| **day only** | calendar day | constant within a held-out timepoint; interpolates the training day→risk curve |
| **progress only** | per-cell progress | the arm under test |
| **day + progress** | both | |
| *(context)* full expression → risk | all 2000 genes | **not part of the gate.** Reported as a ceiling so a weak result can be read against what the data supports at all |

### Heads — run separately, never combined

Per P4, pre-specified: **no composite endpoint.**

- **identity-loss head:** `y = 1` if the cell's hard label is LOSS
- **apoptosis head:** `y = 1` if the cell's hard label is DEATH

---

## 3. 🔒 PRE-REGISTERED OUTCOMES

Unit is the **fold** (9 held-out timepoints), paired across arms, `t(0.975, df=8) = 2.306`.

| # | result | verdict |
|---|---|---|
| **G1** | paired (progress − day) log-loss CI **excludes 0, negative** | **PROGRESS BEATS DAY** |
| **G2** | CI **includes 0** | **TIES** |
| **G3** | CI **excludes 0, positive** | **PROGRESS LOSES** |

### The gate — fixed now, not adjustable afterwards

- **P5 becomes eligible** if **G1** fires on **at least one head**. The two heads are distinct
  biology and the tool needs whichever is predictable.
- **STOP** if both heads return G2 or G3. **P5 does not happen.**

**This gate is not modified based on what the result looks like.** If progress fails within the
single line it was designed on, it does not get a second chance on another line.

### Reported but NOT part of the gate

- The context arm. A large gap between it and every gated arm means the risk is predictable but
  neither coordinate captures it — informative, and not a pass.
- Per-fold results for the two extreme timepoints (earliest, latest), where the day arm must
  **extrapolate** rather than interpolate. Flagged so a win driven only by those is visible.
- The iPSC timepoint, a terminal state rather than a point on the reprogramming course.

---

## 4. Declared limits

1. **One line.** Nothing here speaks to transfer; that is exactly what P5 would test.
2. **Progress is defined as predicted-day-from-expression.** That is one construction of
   "trajectory position", not the only one. A negative result licenses "this construction fails",
   **not** "no progress coordinate can work" — and G3 explicitly routes to reconsidering the
   construction before concluding anything.
3. **Cells are not independent** within a timepoint. The fold, not the cell, is the unit for every
   interval.
4. **Labels are the pipeline's own marker-based calls**, inheriting whatever they inherit.
5. No retrain, no harmonizer change, `src/` untouched.

## 5. Recording

`results/p3_progress_results.json`; write-up to the work order, `CHANGES.md` and the notebook;
unit tests in `tests/test_p3_progress.py`. Every outcome graded as written, including failures.

---

## 6. RESULT — 2026-08-14. **G2 TIES on both heads. STOP. P5 does not happen.**

*Graded against §3 exactly as written. Artefacts: `experiments/p3_progress.py`,
`results/p3_progress_results.json`, `tests/test_p3_progress.py` (13 tests).*
42,481 HFF cells, 9 timepoints, leave-one-timepoint-out.

### Mean held-out log-loss (lower is better)

| head | day | **progress** | day+prog | day_interp | *context* |
|---|---|---|---|---|---|
| identity loss | 0.6573 | **0.6569** | 0.6553 | **0.6096** | *1.2419* |
| apoptosis | 0.1399 | **0.1191** | 0.1268 | **0.1170** | *0.3769* |

| head | paired (progress − day) | 95 % CI | verdict |
|---|---|---|---|
| identity loss | −0.0004 | [−0.0937, +0.0930] | **G2 TIES** |
| apoptosis | −0.0208 | [−0.0690, +0.0274] | **G2 TIES** |

### 🔒 The gate, applied as written

**Both heads G2 → STOP. P5 does not happen.** The gate was fixed before the run and is not
adjusted to the result.

### The pre-flagged checks, and neither rescues it

- **The terminal iPSC state was flagged in advance** as a fold that could distort the verdict.
  Dropping it: identity loss −0.0339 [−0.0948, +0.0271], apoptosis −0.0236 [−0.0791, +0.0320].
  **Both still include zero.** The tie is not an artefact of that timepoint.
- **`day_interp` — the strongest fair day arm — BEATS progress on both heads** (0.6096 vs 0.6569;
  0.1170 vs 0.1191). Given its best form, day is not merely tied with progress but slightly ahead.

### What the context arm revealed, which was not the plan

It was included as a **ceiling**. It is a **floor**: full 2000-gene expression scores **1.2419**
and **0.3769**, roughly **twice as bad** as a single scalar coordinate. Out of time, a naive
high-dimensional model is badly calibrated, and **none of these predictors is strong in absolute
terms** — identity-loss log-loss ~0.66 is close to what predicting a constant base rate achieves.

**The honest summary is not "day wins". It is that at this held-out geometry, none of these
coordinates predicts per-cell risk well, and progress adds nothing over day.**

### What this does and does not license

| | |
|---|---|
| **P5** (`GSE221739` cross-line test) | ❌ **does not happen.** The gate is not relaxed |
| the reframing *as an idea* | **not refuted in general** — see the limit below |
| **this construction** of progress | ❌ **fails.** Predicted-day-from-expression carries no risk information beyond day itself |
| Test 11.1 / Test 18's "time is redundant with state" | **corroborated**, now per-cell, out-of-time, and on a proper scoring rule |

**The limit declared in §4.2 governs the reading.** Progress here is
*predicted-day-from-expression* — one construction of trajectory position, and one that is
regularised toward day by design. An unsupervised construction (diffusion pseudotime, an
independent trajectory fit) is **untested** and this result does not speak to it.

**But the gate does not care.** §3 routes G2 on both heads to STOP, and a different construction
would need its own pre-registration and its own justification for why it should be tried after
this one failed — not an extension of this one.

### What is NOT claimed

That per-cell risk is unpredictable in principle — the context arm's poor showing is a statement
about naive high-dimensional models out of time, not about the biology. That day is a *good*
predictor; it is merely not worse. Or that any of this transfers to another line, which is exactly
the question P5 would have asked and now will not.
