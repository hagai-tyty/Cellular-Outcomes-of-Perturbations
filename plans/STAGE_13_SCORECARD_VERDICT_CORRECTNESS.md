# Stage 13 — the scorecard's verdicts are computed on the wrong quantity

**Status:** PLAN, then EXECUTE. Changes `scorecard.py` only (aggregation + verdict + printing).
**Does NOT change `measure_fold`, any stored snapshot, or any model.** See §13.7.

---

## 13.1 Why this is Stage 13 and not something to defer

`scorecard.py` is the instrument that decides whether a change to CellFate-Rx is **accepted or
reverted**. Every Change in `CHANGES.md` since 2026-07-19 was judged by it. A defect here does not
degrade a result — it **inverts a decision**.

Two defects were found while reading the C-7 comparison and recorded in `CHANGES.md` as *"Neither
is fixed yet; both need their own change."* This is that change. Below, both are re-measured from
the snapshots on disk rather than taken from that note.

**The central discovery of this plan** (not in the original note): defect 1 has **three distinct
faces**, and the one that matters most was not the one recorded.

---

## 13.2 Defect A — signed metrics judged as if they were magnitudes

`METRICS` marks `level_shift_model` and `level_shift_ridge` as direction `"abs"`, and
`scorecard.py:378` states the intent outright: *"judge |level shift|, not signed"*. The code does
not do that. `abs()` is applied to **the aggregate**, once, at display time — and to nothing else.

### A1 — the aggregate cancels instead of accumulating *(the worst face; previously unrecorded)*

`_agg` averages the **signed** per-fold shifts, then `abs()` is taken of that mean. Level shift is a
**per-donor** bias whose sign varies by donor, so the average measures **cancellation across the
donor panel**, not the size of the error.

Measured on the snapshots (`scorecard/gc2_A_keep_hff.json`, `scorecard/c7_A_keep_hff.json`):

| metric | snapshot | printed `|mean(signed)|` | true `mean(|shift|)` | understated by |
|---|---|---|---|---|
| `level_shift_model` | gc2_A | 5.713 | **13.120** | 2.3× |
| `level_shift_model` | c7_A | 4.830 | **9.621** | 2.0× |
| `level_shift_ridge` | gc2_A | **0.230** | **12.723** | **55.2×** |
| `level_shift_ridge` | c7_A | 0.644 | 11.162 | 17.3× |

The ridge row prints **0.230** for a quantity whose true mean magnitude is **12.72 yr**.

This is not an obscure row. **±12.7 yr per-donor level shift is the founding measurement of
Stage 2** (`MASTER_PLAN.md:81, 387, 452`) and the entire justification for buying k≈3 reference
cells. The scorecard has been printing that exact quantity as ≈ 0 — i.e. as *"there is no level
shift"* — in the one table the project uses to decide things.

`MASTER_PLAN.md:452` records the shift as `±12.7 yr`. `mean(|shift|)` = **12.72**. The correct
statistic reproduces the project's own headline number; the printed one destroys it.

### A2 — the paired difference is computed on signed values

`_paired` (line 377) is called **before** the `abs()` on line 378-379 and reads raw fold values, so
the CI that drives the verdict is built on **signed** differences for a metric declared `"abs"`.

### A3 — `_verdict` then applies `better_is_down` to that signed difference

`_verdict` sets `better_is_down = direction in ("lower", "abs")` and tests the sign of the CI. On a
signed quantity this is meaningless: moving −28 → −22 is an **improvement** in magnitude and reads
as an **increase**.

**The consequence, measured:**

| metric | current output | correct output |
|---|---|---|
| `level_shift_model` | `+5.030` CI `[+1.218,+8.842]` → **REGRESSION** | `−3.118` CI `[−9.100,+2.865]` → **noise** |
| `level_shift_ridge` | `+4.389` CI `[+2.805,+5.972]` → **REGRESSION** | `−0.084` CI `[−5.399,+5.231]` → **noise** |

Both rows reported a **decisive regression**. The correct verdict for both is **no detectable
change, with the point estimate in the improving direction.** The sign of the effect is flipped.

This has misled at least twice in the recorded history — `CHANGES.md:962` (C-7) and
`CHANGES.md:7778` (2026-08-03, *"A reader trusting that row would draw the opposite conclusion"*).

---

## 13.3 Defect B — the two columns can be means over different folds

`_agg(A)` and `_agg(B)` each average over whatever folds are non-null **in their own snapshot**,
while `_paired` correctly restricts to folds present and valid in **both**. The three numbers on one
row can therefore describe three different populations.

In the C-7 comparison N2 errors out in `c7_A` (*"too few age-valid cells"*), so column A is 6 folds
and column B is 5. **13 of 18 metrics are affected.** The clearest case:

| metric | col A (n=6) | col B (n=5) | implied diff | reported mean diff (n=5) |
|---|---|---|---|---|
| `dage_mae_model` | 14.291 | 15.713 | +1.422 | **+2.922** |
| `dage_mae_ridge` | 14.052 | 11.300 | −2.751 | **−1.245** |

The reader sees `+1.422` and a verdict computed from `+2.922` — the verdict-driving number is more
than **double** the visible one, and nothing on the row says why.

The same defect sits in the RES over-approval block (lines 401-407), which calls `_agg` on each
snapshot independently.

---

## 13.4 The fix

### F1 — aggregate `"abs"` metrics as magnitudes, and keep the signed mean as its own row

Per-fold cells stay **signed** — the direction of a donor's shift is real information and must not
be hidden. What changes is the aggregate and the paired statistic:

    aggregate     mean(|shift_d|)                 <- the criterion (matches the "|.| lower better" docstring)
    paired diff   mean(|shift_d^B| - |shift_d^A|)  <- the accept/reject statistic
    verdict       lower-is-better on the above     <- now applied to a genuine magnitude

**And print both.** `mean(signed)` is not garbage — it answers a *different, also useful* question
("is there a **global** offset, or do donors cancel?"). Replacing one number with another loses
that. So an `"abs"` metric gets its magnitude row (judged) **plus** a `signed mean` context row
(never judged). This adds information instead of trading one blind spot for another.

*Rejected alternative:* changing `measure_fold` to store `abs(...)`. That destroys the per-donor
sign permanently, breaks `experiments/diag_zero_point.py:447` (which reads the signed per-donor
shifts and is correct to), and makes every snapshot on disk unreadable in its own terms.

### F2 — comparison columns must use the paired fold set

Compute both aggregates over exactly the folds `_paired` uses, and state which folds those are,
naming any donor dropped and why.

This buys a **checkable invariant**, which is what makes the fix testable rather than merely
plausible:

> **`col_B − col_A == mean diff`, to floating-point tolerance, for every metric and every
> direction — including `"abs"`.**

That identity is false in the current code for 13 of 18 rows and must hold for all of them after.
It is the single strongest test in this stage: it cannot pass unless F1 and F2 are *both* right,
because for an `"abs"` metric it only holds if the columns and the paired statistic are taking the
magnitude at the same point in the computation.

### F3 — the snapshot table's `mean` column

`_print_snapshot` has the same A1 flaw. Its `mean` column for an `"abs"` metric becomes
`mean(|.|)`, explicitly labelled so no reader has to guess which one they are looking at.

---

## 13.5 Retro-verdicts — what this stage can do that Stage 12 could not

Stage 12's defect was in **written artefacts**, so fixing it could not repair the past without a
rebuild. This defect is in **aggregation and verdict only**; every per-fold number on disk is
correct and signed.

**Therefore every past comparison can be re-judged correctly, right now, with no rebuild and no
retrain.** This stage produces a retro-verdict table over all snapshot pairs in `scorecard/`,
recording where the corrected rule changes a verdict.

Per the standing rule, **no past record is rewritten.** The original entries stay exactly as
written; the retro-table is recorded beside them as *what we now believe*.

---

## 13.6 Why the correction is not itself a judgement call

A fix to a decision rule can smuggle in a new bias. Three reasons this one does not:

1. **It implements the stated intent, it does not change it.** `"judge |level shift|, not signed"`
   was written in the source before this stage; the code simply never did it.
2. **It reproduces the project's own independently-derived number.** `mean(|shift|)` = 12.72 vs the
   `±12.7 yr` recorded in `MASTER_PLAN.md` from Test 7.4.3. The printed 0.230 does not.
3. **It is not directionally convenient.** It converts two `REGRESSION`s into `noise` — but it also
   raises the reported level shift from 0.230 to 12.72, which makes the model look **worse**, and
   restores `dage_mae_model`'s true `+2.922` degradation in place of a visible `+1.422`. It is not
   a fix that flatters the project.

---

## 13.7 What this stage explicitly does NOT do

- **No model change, no rebuild, no retrain, no re-snapshot.** Aggregation and printing only.
- **`measure_fold` is untouched** — stored per-fold values keep their sign.
- **No stored snapshot is modified or rewritten.** All nine remain readable and are re-judged in
  place.
- **No past record is edited.** Corrections are appended beside the originals (standing rule).
- **The `"neutral"` and `"lower"`/`"higher"` metrics keep their current semantics.** Only the
  fold-set alignment (F2) touches them, and that is a pure correctness fix.

---

## 13.8 Verification

| item | how |
|---|---|
| A1 — magnitude aggregation | test: signed shifts that cancel (`+10, −10`) must aggregate to 10, not 0 |
| A1 — real data | test pinning `mean(\|level_shift_ridge\|)` = 12.72 on `gc2_A`, vs the 0.230 the old path gave |
| A2/A3 — verdict on magnitudes | test: `−28 → −22` across folds reads as improvement, not increase |
| A3 — the exact false REGRESSION | regression test on the two real snapshots: both rows must read `noise`, not `REGRESSION` |
| B — fold alignment | test: a snapshot pair where B is missing a fold must aggregate both columns over the common set |
| **the invariant** | test asserting `col_B − col_A == mean diff` for **every** metric, on the real snapshot pair and on synthetic pairs including `"abs"` |
| signed row retained | test that the signed mean is still reported for `"abs"` metrics |
| no measurement change | test that `measure_fold`'s stored value is still signed |
| nothing else moved | full `pytest` green, `ruff check` clean |
| record | `CHANGES.md` entry + retro-verdict table; originals annotated, not rewritten |
