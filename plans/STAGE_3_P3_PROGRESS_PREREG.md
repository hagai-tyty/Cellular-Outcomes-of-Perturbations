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
