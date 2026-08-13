# REGIME E — PRE-REGISTRATION. Written and committed BEFORE the run.

**Status:** 🔒 pre-registration. Committed before `experiments/stage3a_regime_e.py` was executed.
**Decides:** whether Stage 3a's forward gate is gradeable on `GSE165177`, which is already on disk
and in no training config — i.e. whether the data-acquisition ask in
`DATA_REQUIREMENT_SECOND_TIMECOURSE.md` is needed at all.
**Prompted by:** the other machine's `AUDIT-2` (`22de08a`), whose central claim I verified against
the raw GEO matrix and accepted.

---

## 1. Why this run exists

3a-bis measured that the held-out `gill_bulk` geometry is UNRESOLVABLE at every effect size — a
correct system could never have returned GO. I concluded "acquire a second dense time course."
**That conclusion skipped a dataset already in hand.**

| | `gill_bulk` — what 3a ran on | **`GSE165177`** |
|---|---|---|
| donors | 6 | 3 (O1, O2, O3) |
| donor ages | 0–53, several outside the clock's `[1, 96]` | **53, 53, 38 — all adult, all in range** |
| samples per (donor, timepoint) | ≈ **1.7** | **6–9** |
| contemporaneous negative controls | **6 total** (1/donor, day 0 only) | **33** (2–3 per donor *per timepoint*) |
| usable time positions | 12 | 5 — day 0, 10, 13, 15, 17 |
| ordered forward pairs per donor | 66 | **10** |

The mechanism 3a-bis blamed is `gill_bulk`-specific: at 1.7 samples a fraction can only be
{0, 0.5, 1}, hence 63 of 70 profile values pinned at the bounds. At 6–9 samples it moves in sixths
to ninths. **Different target. Never measured.**

---

## 2. What will be run

`experiments/stage3a_regime_e.py`, **READ-ONLY** — no build, no retrain, `src/` untouched.

- **Source:** `D:\GSE165177` raw `Log2_RPM` matrices (24 samples `exp1`, 71 samples `exp2`, 95
  total) + `GSE165177_series_matrix.txt.gz` for donor / day / arm / age.
- **Fate labels:** the pipeline's own `cellfate.data.labels.fate_labels` — pluripotency,
  somatic-identity and apoptosis marker programs, z-scored **against each donor's own
  `negative_control` samples**. Hard label = argmax; `unsafe = LOSS or DEATH`, exactly as
  `test18_forward_gate.py:85`.
- **Target:** unsafe fraction over the **non-control** samples at each (donor, day).
- **Trajectory unit:** **donor** (3). The only defensible independence unit here — arms within a
  donor share the culture.
- **Null:** identical to 3a-bis — HFF's measured unsafe curve scaled by `alpha`, observed at the
  **real per-(donor, day) sample counts** via `Binomial(n, p)/n`.
- **Bar:** 3a's own rule, graded verbatim — **PASS iff the paired 95 % CI upper end < 0**.
  `MIN_PASS_RATE = 0.95`, 2000 trials per cell, `alpha ∈ {0, 0.25, 0.5, 1.0}`, raw and logit.

### The 2×2 that attributes any failure

3 folds is a severe power cost — the paired CI uses `t(0.975, df=2) = 4.303` against 2.776 at 6.
So fold count and replication are **again** confounded unless separated deliberately:

| | at `gill_bulk`'s ≈2 samples/tp | at `GSE165177`'s real 4–9 samples/tp |
|---|---|---|
| **3 folds** (real) | E-a | **E (the operative cell)** |
| **6 folds** (counterfactual) | E-c | E-d |

---

## 3. 🔒 PRE-REGISTERED OUTCOMES — what each result MEANS and what it TRIGGERS

Graded on the **logit** estimator at the **operative cell** (3 folds, real sample counts), because
3a-bis established the raw estimator is not a detector.

| # | Result | Reading | Action it triggers |
|---|---|---|---|
| **E1** | pass ≥ 0.95 at **α = 1.0 AND α = 0.5** | **RESOLVABLE.** The gate can be graded on data already held, for effects down to half HFF's amplitude. | **The acquisition ask is WITHDRAWN.** Proceed to a graded 3a run on `GSE165177`, pre-registered separately. No expert conversation needed. |
| **E2** | pass ≥ 0.95 at **α = 1.0 only** | **MARGINAL.** Gradeable only for a full-amplitude effect. | Run the graded 3a, but pre-register in advance that **a null result is uninformative about anything smaller than α = 1** and may not be read as absence of signal. Acquisition ask stays ON HOLD, not withdrawn. |
| **E3** | pass < 0.95 at **both α = 1.0 and α = 0.5** | **UNRESOLVABLE.** Even replicated contemporaneous controls at 6–9 samples/timepoint cannot carry the bar. | **The acquisition ask is ESTABLISHED — with a number, not an assertion.** The expert question sharpens to: *"replicated controls at 6–9 samples/timepoint across 3 adult donors still are not enough; what would be?"* |
| **E4** | pass > 0.05 at **α = 0** | **THE NULL IS BROKEN.** A no-effect world must not pass. | **Discard the entire run**, fix the simulation, re-register. No conclusion may be drawn from E1–E3 in this case. |

### Attribution, read only if E3 fires

| # | Result | Reading | Action |
|---|---|---|---|
| **A1** | E-d (6 folds, real n) ≥ 0.95 while E (3 folds, real n) < 0.95 | **FOLD COUNT binds.** Replication is already sufficient. | The ask becomes **more donors**, not more cells per timepoint — a materially cheaper and more findable request. Rewrite the acquisition spec accordingly. |
| **A2** | E (3 folds, real n) ≥ 0.95 while E-a (3 folds, ≈2/tp) < 0.95 | **REPLICATION binds**, and `GSE165177` already clears it. | Confirms the 3a-bis mechanism and E1/E2 governs. |
| **A3** | both E-a and E-d < 0.95 | **BOTH bind.** Neither fix alone suffices. | The ask needs **both** more donors and more replication — state both numbers to the expert. |

### Precondition — checked and reported BEFORE the null is read

| # | Condition | Consequence |
|---|---|---|
| **P0** | the observed unsafe fraction has **zero variance across timepoints in ≥ 2 of 3 donors** | **REGIME E IS VOID.** A target with no time variation cannot be forward-predicted regardless of what the null says, and `test18`'s own rule already skips such a fold. Report and stop; do not read E1–E4. |

---

## 4. Declared limits of this run — stated before seeing any number

1. **Safety target only.** No ΔAge, no clock, no harmonizer. This measures the geometry of the
   `p_unsafe` question, which 3a itself calls *"THE DECISIVE ONE"*. It says nothing about ΔAge.
2. **No cross-dataset training.** `GSE165177` is full-transcriptome in its own gene space; the
   `gill`/HFF fold bundles are in the 2000-gene panel space. Joining them needs a gene-space
   mapping that is its own piece of work. Regime E therefore trains **within `GSE165177`**
   (2 donors) and holds out the third. This is the honest floor, and it **understates** what a
   joined training set could do.
3. **Therefore E3 does not by itself prove acquisition is required** — it proves it for the
   within-`GSE165177` training geometry. If E3 fires, the gene-space join is the next thing to
   cost **before** any acquisition, and that is recorded here in advance.
4. **`exp1` / `exp2` batch structure** (24 vs 71 samples) is present and is **not** corrected for.
   It will be reported per fold so it cannot be mistaken for signal.
5. **The audit's "no modality shift" claim is not tested here.** It is only true of a bulk-trained
   variant; testing it needs the gene-space join in limit 2.

---

## 5. Recording

Results to `results/stage3a_regime_e_results.json`; write-up appended to `plans/STAGE_3_TOOL.md`,
`CHANGES.md` and `experiments/DELTAAGE_LAB_NOTEBOOK.md`; unit tests in
`tests/test_stage3a_regime_e.py`. Outcomes graded against §3 **as written above**, including any
that fail.

---

## 6. RESULT — 2026-08-12. **P0 FIRED. Regime E is VOID**, and the reason generalises.

*Graded against §3 exactly as written above. Artefacts: `experiments/stage3a_regime_e.py`,
`results/stage3a_regime_e_results.json`, `tests/test_stage3a_regime_e.py` (26 tests).*
**READ-ONLY** — no build, no retrain, `src/` untouched.

93 of 95 samples loaded (the 2 excluded are the `iPSC 13` / `iPSC 21` lines, which the series
files under donor `iPSC`, day 51 — not a point on any donor's trajectory), 35,720 genes, donors
O1/O2/O3 aged 53/53/38, 33 contemporaneous negative controls, batches `exp1` 24 / `exp2` 69.

### The pre-registered precondition fired

| donor | age | day 0 | day 10 | day 13 | day 15 | day 17 | SD across timepoints |
|---|---|---|---|---|---|---|---|
| O1 | 53 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 | 0.100 |
| O2 | 53 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |
| O3 | 38 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |

**2 of 3 donors have zero variation → P0 fires → E1–E4 must not be read**, per §3. They have not
been read, and the null was not run.

### 🔬 The reason, which matters far more than the verdict

**The untreated day-0 fibroblasts label `loss` — P(loss) = 0.966, 0.876, 0.730.** They are the
starting material; somatic identity is what they *are*. That cannot be biology, so it was chased
down before P0 was reported:

| arm | n | P(safe) | P(loss) | P(death) | hard labels |
|---|---|---|---|---|---|
| **day0_fibroblast** | 3 | 0.071 | **0.857** | 0.071 | **loss = 3** |
| failed_to_transiently_reprogram | 21 | 0.019 | 0.962 | 0.019 | loss = 21 |
| failing_to_transiently_reprogram | 12 | 0.000 | 0.998 | 0.002 | loss = 12 |
| **negative_control** | 33 | **0.699** | 0.163 | 0.138 | safe = 30, loss = 3 |
| transient_reprogramming | 11 | 0.000 | 1.000 | 0.000 | loss = 11 |
| transiently_reprogrammed | 13 | 0.059 | 0.613 | 0.328 | loss = 9, death = 3, safe = 1 |

`fate_labels` z-scores each program against the `is_control` samples. Here those are fibroblasts
cultured 10–17 days, so **anything that differs from them lands on the unsafe side** — including
the day-0 starting material. The split produced is **control vs non-control, not a time course.**

### The structural finding — it is not about this dataset

> **`p_unsafe` is a fraction OF CELLS. A bulk RNA-seq sample is ALREADY a population average, so a
> hard label per sample collapses that fraction to 0 or 1 before it can be counted, and the
> "fraction" becomes a fraction of SAMPLES.**

This explains, with one mechanism, two things previously attributed to sample size:

- **`gill_bulk`'s 63 of 70 profile values pinned at the bounds.** Blamed on 1.7 samples per
  timepoint. Wrong: at 4–6 samples per timepoint with real contemporaneous controls, `GSE165177`
  is *more* saturated, not less — 11 of 12 cells at exactly 1.000.
- **Why more bulk replication cannot help.** The quantity is not expressible in bulk at any
  replication. This is a property of the measurement, not of the sample count.

**It also means regime E could never have succeeded**, and the 2×2 attribution (A1/A2/A3) is moot:
fold count and replication were never the binding constraint for the *safety* target.

### What this does and does not settle

| | |
|---|---|
| the acquisition ask in `DATA_REQUIREMENT_SECOND_TIMECOURSE.md` | **RE-ESTABLISHED, and on much better grounds than I first gave.** Not "we need more replication" but "`p_unsafe` requires single-cell resolution, and GSE242423 is the only single-cell dataset we have" |
| requirement **H2** of that spec (single-cell, not bulk) | **promoted from a requirement to the whole point.** It was listed alongside seven others; it is the one that binds |
| AUDIT-2's central challenge — *"you skipped a dataset on disk"* | **was right to make, and testing it was the correct call.** It cost one script run and replaced an assertion with a mechanism |
| AUDIT-2's expectation that `GSE165177` would carry the gate | **refuted**, for a reason neither of us had identified |
| `GSE165177`'s value to the project | **NOT dismissed.** Its 4–6× replication, 33 contemporaneous controls and three adult in-range donors (53/53/38) are real advantages — for **ΔAge**, which is continuous per sample and does not have this problem at all. That is a separate, open, and promising question |
| my own §2 pair count (5 positions, 10 pairs/donor) | **held** — but only after fixing a bug that silently dropped it to 4 and 6; see below |

### ⚠️ A bug this run caught in itself, recorded because nothing failed

The first implementation matched only `{donor}_{arm}_{N}days_{exp}`. `GSE165177` names its day-0
fibroblasts **`O1 Fib`** — space-separated, no `days` token — so **all three were dropped
silently**, cutting every trajectory's first timepoint and the pair count from 10 per donor to 6.
No error, no warning; the run simply answered a smaller question and would have agreed with
AUDIT-2's "4 timepoints, 6 pairs" for the wrong reason. Fixed, and pinned by
`tests/test_stage3a_regime_e.py::test_parses_the_day_zero_fibroblast_that_the_first_implementation_dropped`.

**Including day 0 did not rescue P0** — the day-0 samples label `loss` too, which is what exposed
the mechanism above. The bug mattered for the pair count and for the diagnosis, not for the verdict.

### Declared limits — unchanged from §4, all still apply

Safety target only; trains within `GSE165177`; `exp1`/`exp2` batch not corrected; the "no modality
shift" claim untested. **Limit 3 is now moot for a better reason:** the gene-space join was listed
as the thing to cost before acquisition, but joining gene spaces cannot make a bulk sample express
a per-cell fraction. **Do not spend on the join for the safety target.** It remains open and
worthwhile for ΔAge.

---

## 6. RESULT — 2026-08-12. **P0 FIRED. Regime E is VOID**, and the reason generalises.

*Graded against §3 exactly as written above. Artefacts: `experiments/stage3a_regime_e.py`,
`results/stage3a_regime_e_results.json`, `tests/test_stage3a_regime_e.py` (26 tests).*
**READ-ONLY** — no build, no retrain, `src/` untouched.

93 of 95 samples loaded (the 2 excluded are the `iPSC 13` / `iPSC 21` lines, which the series
files under donor `iPSC`, day 51 — not a point on any donor's trajectory), 35,720 genes, donors
O1/O2/O3 aged 53/53/38, 33 contemporaneous negative controls, batches `exp1` 24 / `exp2` 69.

### The pre-registered precondition fired

| donor | age | day 0 | day 10 | day 13 | day 15 | day 17 | SD across timepoints |
|---|---|---|---|---|---|---|---|
| O1 | 53 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 | 0.100 |
| O2 | 53 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |
| O3 | 38 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |

**2 of 3 donors have zero variation → P0 fires → E1–E4 must not be read**, per §3. They have not
been read, and the null was not run.

### 🔬 The reason, which matters far more than the verdict

**The untreated day-0 fibroblasts label `loss` — P(loss) = 0.966, 0.876, 0.730.** They are the
starting material; somatic identity is what they *are*. That cannot be biology, so it was chased
down before P0 was reported:

| arm | n | P(safe) | P(loss) | P(death) | hard labels |
|---|---|---|---|---|---|
| **day0_fibroblast** | 3 | 0.071 | **0.857** | 0.071 | **loss = 3** |
| failed_to_transiently_reprogram | 21 | 0.019 | 0.962 | 0.019 | loss = 21 |
| failing_to_transiently_reprogram | 12 | 0.000 | 0.998 | 0.002 | loss = 12 |
| **negative_control** | 33 | **0.699** | 0.163 | 0.138 | safe = 30, loss = 3 |
| transient_reprogramming | 11 | 0.000 | 1.000 | 0.000 | loss = 11 |
| transiently_reprogrammed | 13 | 0.059 | 0.613 | 0.328 | loss = 9, death = 3, safe = 1 |

`fate_labels` z-scores each program against the `is_control` samples. Here those are fibroblasts
cultured 10–17 days, so **anything that differs from them lands on the unsafe side** — including
the day-0 starting material. The split produced is **control vs non-control, not a time course.**

### The structural finding — it is not about this dataset

> **`p_unsafe` is a fraction OF CELLS. A bulk RNA-seq sample is ALREADY a population average, so a
> hard label per sample collapses that fraction to 0 or 1 before it can be counted, and the
> "fraction" becomes a fraction of SAMPLES.**

This explains, with one mechanism, two things previously attributed to sample size:

- **`gill_bulk`'s 63 of 70 profile values pinned at the bounds.** Blamed on 1.7 samples per
  timepoint. Wrong: at 4–6 samples per timepoint with real contemporaneous controls, `GSE165177`
  is *more* saturated, not less — 11 of 12 cells at exactly 1.000.
- **Why more bulk replication cannot help.** The quantity is not expressible in bulk at any
  replication. This is a property of the measurement, not of the sample count.

**It also means regime E could never have succeeded**, and the 2×2 attribution (A1/A2/A3) is moot:
fold count and replication were never the binding constraint for the *safety* target.

### What this does and does not settle

| | |
|---|---|
| the acquisition ask in `DATA_REQUIREMENT_SECOND_TIMECOURSE.md` | **RE-ESTABLISHED, and on much better grounds than I first gave.** Not "we need more replication" but "`p_unsafe` requires single-cell resolution, and GSE242423 is the only single-cell dataset we have" |
| requirement **H2** of that spec (single-cell, not bulk) | **promoted from a requirement to the whole point.** It was listed alongside seven others; it is the one that binds |
| AUDIT-2's central challenge — *"you skipped a dataset on disk"* | **was right to make, and testing it was the correct call.** It cost one script run and replaced an assertion with a mechanism |
| AUDIT-2's expectation that `GSE165177` would carry the gate | **refuted**, for a reason neither of us had identified |
| `GSE165177`'s value to the project | **NOT dismissed.** Its 4–6× replication, 33 contemporaneous controls and three adult in-range donors (53/53/38) are real advantages — for **ΔAge**, which is continuous per sample and does not have this problem at all. That is a separate, open, and promising question |
| my own §2 pair count (5 positions, 10 pairs/donor) | **held** — but only after fixing a bug that silently dropped it to 4 and 6; see below |

### ⚠️ A bug this run caught in itself, recorded because nothing failed

The first implementation matched only `{donor}_{arm}_{N}days_{exp}`. `GSE165177` names its day-0
fibroblasts **`O1 Fib`** — space-separated, no `days` token — so **all three were dropped
silently**, cutting every trajectory's first timepoint and the pair count from 10 per donor to 6.
No error, no warning; the run simply answered a smaller question and would have agreed with
AUDIT-2's "4 timepoints, 6 pairs" for the wrong reason. Fixed, and pinned by
`tests/test_stage3a_regime_e.py::test_parses_the_day_zero_fibroblast_that_the_first_implementation_dropped`.

**Including day 0 did not rescue P0** — the day-0 samples label `loss` too, which is what exposed
the mechanism above. The bug mattered for the pair count and for the diagnosis, not for the verdict.

### Declared limits — unchanged from §4, all still apply

Safety target only; trains within `GSE165177`; `exp1`/`exp2` batch not corrected; the "no modality
shift" claim untested. **Limit 3 is now moot for a better reason:** the gene-space join was listed
as the thing to cost before acquisition, but joining gene spaces cannot make a bulk sample express
a per-cell fraction. **Do not spend on the join for the safety target.** It remains open and
worthwhile for ΔAge.

---

## 6. RESULT — 2026-08-12. **P0 FIRED. Regime E is VOID**, and the reason generalises.

*Graded against §3 exactly as written above. Artefacts: `experiments/stage3a_regime_e.py`,
`results/stage3a_regime_e_results.json`, `tests/test_stage3a_regime_e.py` (26 tests).*
**READ-ONLY** — no build, no retrain, `src/` untouched.

93 of 95 samples loaded (the 2 excluded are the `iPSC 13` / `iPSC 21` lines, which the series
files under donor `iPSC`, day 51 — not a point on any donor's trajectory), 35,720 genes, donors
O1/O2/O3 aged 53/53/38, 33 contemporaneous negative controls, batches `exp1` 24 / `exp2` 69.

### The pre-registered precondition fired

| donor | age | day 0 | day 10 | day 13 | day 15 | day 17 | SD across timepoints |
|---|---|---|---|---|---|---|---|
| O1 | 53 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 | 0.100 |
| O2 | 53 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |
| O3 | 38 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **0.000** |

**2 of 3 donors have zero variation → P0 fires → E1–E4 must not be read**, per §3. They have not
been read, and the null was not run.

### 🔬 The reason, which matters far more than the verdict

**The untreated day-0 fibroblasts label `loss` — P(loss) = 0.966, 0.876, 0.730.** They are the
starting material; somatic identity is what they *are*. That cannot be biology, so it was chased
down before P0 was reported:

| arm | n | P(safe) | P(loss) | P(death) | hard labels |
|---|---|---|---|---|---|
| **day0_fibroblast** | 3 | 0.071 | **0.857** | 0.071 | **loss = 3** |
| failed_to_transiently_reprogram | 21 | 0.019 | 0.962 | 0.019 | loss = 21 |
| failing_to_transiently_reprogram | 12 | 0.000 | 0.998 | 0.002 | loss = 12 |
| **negative_control** | 33 | **0.699** | 0.163 | 0.138 | safe = 30, loss = 3 |
| transient_reprogramming | 11 | 0.000 | 1.000 | 0.000 | loss = 11 |
| transiently_reprogrammed | 13 | 0.059 | 0.613 | 0.328 | loss = 9, death = 3, safe = 1 |

`fate_labels` z-scores each program against the `is_control` samples. Here those are fibroblasts
cultured 10–17 days, so **anything that differs from them lands on the unsafe side** — including
the day-0 starting material. The split produced is **control vs non-control, not a time course.**

### The structural finding — it is not about this dataset

> **`p_unsafe` is a fraction OF CELLS. A bulk RNA-seq sample is ALREADY a population average, so a
> hard label per sample collapses that fraction to 0 or 1 before it can be counted, and the
> "fraction" becomes a fraction of SAMPLES.**

This explains, with one mechanism, two things previously attributed to sample size:

- **`gill_bulk`'s 63 of 70 profile values pinned at the bounds.** Blamed on 1.7 samples per
  timepoint. Wrong: at 4–6 samples per timepoint with real contemporaneous controls, `GSE165177`
  is *more* saturated, not less — 11 of 12 cells at exactly 1.000.
- **Why more bulk replication cannot help.** The quantity is not expressible in bulk at any
  replication. This is a property of the measurement, not of the sample count.

**It also means regime E could never have succeeded**, and the 2×2 attribution (A1/A2/A3) is moot:
fold count and replication were never the binding constraint for the *safety* target.

### What this does and does not settle

| | |
|---|---|
| the acquisition ask in `DATA_REQUIREMENT_SECOND_TIMECOURSE.md` | **RE-ESTABLISHED, and on much better grounds than I first gave.** Not "we need more replication" but "`p_unsafe` requires single-cell resolution, and GSE242423 is the only single-cell dataset we have" |
| requirement **H2** of that spec (single-cell, not bulk) | **promoted from a requirement to the whole point.** It was listed alongside seven others; it is the one that binds |
| AUDIT-2's central challenge — *"you skipped a dataset on disk"* | **was right to make, and testing it was the correct call.** It cost one script run and replaced an assertion with a mechanism |
| AUDIT-2's expectation that `GSE165177` would carry the gate | **refuted**, for a reason neither of us had identified |
| `GSE165177`'s value to the project | **NOT dismissed.** Its 4–6× replication, 33 contemporaneous controls and three adult in-range donors (53/53/38) are real advantages — for **ΔAge**, which is continuous per sample and does not have this problem at all. That is a separate, open, and promising question |
| my own §2 pair count (5 positions, 10 pairs/donor) | **held** — but only after fixing a bug that silently dropped it to 4 and 6; see below |

### ⚠️ A bug this run caught in itself, recorded because nothing failed

The first implementation matched only `{donor}_{arm}_{N}days_{exp}`. `GSE165177` names its day-0
fibroblasts **`O1 Fib`** — space-separated, no `days` token — so **all three were dropped
silently**, cutting every trajectory's first timepoint and the pair count from 10 per donor to 6.
No error, no warning; the run simply answered a smaller question and would have agreed with
AUDIT-2's "4 timepoints, 6 pairs" for the wrong reason. Fixed, and pinned by
`tests/test_stage3a_regime_e.py::test_parses_the_day_zero_fibroblast_that_the_first_implementation_dropped`.

**Including day 0 did not rescue P0** — the day-0 samples label `loss` too, which is what exposed
the mechanism above. The bug mattered for the pair count and for the diagnosis, not for the verdict.

### Declared limits — unchanged from §4, all still apply

Safety target only; trains within `GSE165177`; `exp1`/`exp2` batch not corrected; the "no modality
shift" claim untested. **Limit 3 is now moot for a better reason:** the gene-space join was listed
as the thing to cost before acquisition, but joining gene spaces cannot make a bulk sample express
a per-cell fraction. **Do not spend on the join for the safety target.** It remains open and
worthwhile for ΔAge.
