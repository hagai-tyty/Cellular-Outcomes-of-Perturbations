# CHANGES

Running log of every modification to this repository, newest first.

**Convention.** One entry per stage or task. Every entry states **what** changed, **why**, and
**whether it has been executed**. Nothing is marked verified until it has actually run on the data
machine — "the code looks right" is not verification and is recorded as such.

Files added by the user (`scorecard.py`, `test18_forward_gate.py`, `plans/*` except the deviation
log, `experiments/score + test 18.docx`) are noted where relevant but are not entries here.

---

## 2026-08-18 - STAGE 18: the fate head IS partly reading a clock -- and there is real signal underneath, on 12 pairs

**Status:** MEASURED, read-only. `experiments/diag_stage18_fate_beyond_day.py`,
`results/diag_stage18_fate_beyond_day_results_s16.json`, 17 tests. Run on `_s16` (the recalibrated
folds). **`src/` untouched.**

`fate_prauc` 0.958 is the only strong number the project has. `dose_time` is a MODEL INPUT and it
encodes the timepoint, so this asks how much of that survives once the clock is taken away.

### The structural fact that makes every marginal metric suspect

On the held-out donors, **fate is very nearly a function of timepoint**:

| fold | timepoints | carrying >1 class | time-only PR-AUC | model PR-AUC |
|---|---|---|---|---|
| N2 | 11 | **0** | n/a (19/19 safe) | n/a |
| N3 | 12 | 1 | 0.892 | 1.000 |
| **O1** | 12 | **0** | **1.000** | 0.994 |
| **O2** | 12 | **0** | **1.000** | 1.000 |
| Y1 | 11 | **5** | 0.660 | 0.796 |
| Y2 | 12 | 1 | 0.988 | 1.000 |

**Only 7 of 70 timepoints carry more than one class.** On O1 and O2 the timepoint ALONE reaches
PR-AUC **1.000** -- a lookup table on the hour, using no genes at all, is unbeatable there. Any
model metric on those folds is measuring a calendar.

### The decisive test: within-timepoint concordance

Inside a timepoint `dose_time` is constant and cannot help, so anything the model gets right there
came from expression. Over every (safe, unsafe) pair drawn from the SAME timepoint of the SAME
donor:

> **stratified AUC 0.917 over 12 pairs, permutation p = 0.0091**
> (scores shuffled WITHIN strata, 20,000 draws -- a global shuffle would also destroy the
> between-timepoint structure the model is not being credited for, and would be trivially beaten)

Per fold: N3 1.000 (4 pairs), Y1 0.857 (7 pairs), Y2 1.000 (1 pair). N2/O1/O2 contribute **zero**
pairs -- they have no within-timepoint contrast to offer.

### The verdict, both halves stated

**The marginal 0.93-0.96 is very largely the clock.** That is now measured, not suspected.

**But there IS signal underneath it.** 11 of 12 same-timepoint pairs ranked correctly is
significant against a proper within-stratum null. The fate head is not only reading a calendar.

**And the entire evidence base for that is 12 pairs, from 7 strata, across 3 donors.** It is a
real result on a very thin base, and the honest phrasing is exactly that. This is not a number to
build a headline on until more mixed-timepoint data exists.

### The fold that carries it is the one that looks worst

Y1 supplies **7 of the 12 pairs** and has **5 of the 7 mixed timepoints** -- and it is also the
fold with the WORST marginal PR-AUC (0.796 against 1.000 elsewhere). Those are the same fact seen
twice: Y1 is the only donor whose fate does not track its clock, so it is simultaneously the
hardest fold for a clock-reader and the only fold that can test one.

**Consequence for the LOOCV design:** the aggregate `fate_prauc` is dominated by folds that cannot
discriminate the hypothesis. A future evaluation should weight, or at least report, folds by how
much within-timepoint contrast they actually contain.

### What this does NOT say

It does not say the fate head is worthless -- it says its headline metric is inflated by an input
it was handed. It does not settle whether the within-timepoint signal generalises: 12 pairs cannot.
And it does not touch ΔAge, which was already established as circular for a different reason
(`diag_clock_circularity`, rho 0.96-0.99).

---

## 2026-08-18 - STAGE 16 VERIFIED ON REAL ARTEFACTS: recalibrated, evaluated, and the feared safety loss did NOT occur

**Status:** the hard-label calibrator is now **implemented, tested, AND empirically validated on
recalibrated artefacts**. New: `local_runners/recalibrate_folds.py`, folds
`cellfate_loocv_*_s16`, snapshot `scorecard/c7t_stage16.json`,
`results/stage16_recalibration_results.json`,
`results/diag_stage16_safety_floor_results_{s12,s16}.json`, 12 tests.

The previous entry ended with the fix in `src/` and the honest caveat that no artefact had been
recalibrated. **A code fix plus unit tests is not evidence that shipped behaviour changed.** This
closes that gap.

### How it was recalibrated without retraining

Platt is fitted during training, so a rerun of training would have cost ~5 h. Only that one step
needed redoing. `recalibrate_folds.py` reproduces the shipped calibration block exactly --
`fit_platt_binary(concat(calib P(safe), xdonor P(safe)), concat(hard, hard))` -- taking the
cross-donor pool from the bundle's own `xdonor_stats.npz` (the same 100 cells training used) and
recovering the raw calib probabilities by **inverting** the shipped Platt, which is exact and
asserted to 1e-9 before anything is fitted.

`_s16` is produced by HARDLINKING `_s12` (so 242 MB of shards per fold cost nothing) and replacing
only `bundle/temperature.json`, unlinked first so the original is untouched. **`_s12` and the
`c7t_stage12` snapshot taken from it remain valid**, and a test asserts both fold sets still carry
their own distinct coefficients.

Slopes roughly doubled, as the composed-coefficient prediction said they would: **2.54-2.69 ->
4.82-5.09**.

### THE VERIFICATION -- held-out safety evaluation, same instrument, both fold sets

Pooled over 119 held-out cells (91 truly safe / 28 truly unsafe):

| | `_s12` soft-fit | `_s16` hard-fit | **I predicted** |
|---|---|---|---|
| safety sensitivity | 0.275 | **0.670** | 0.714 |
| **specificity** | 0.929 | **0.929** | 0.821 |
| balanced accuracy | 0.602 | **0.799** | 0.768 |
| false rejections | 66 | **30** | 26 |
| **false approvals** | **2** | **2** | 5 |
| median S, truly safe | 0.703 | **0.874** | -- |
| % truly-safe below the 0.76 bar | 72.5 % | **33.0 %** | -- |

**THE FEARED TRADE DID NOT HAPPEN.** The decision to ship was taken on an expected specificity
drop of **0.929 -> 0.821** and false approvals **2 -> 5**. On the real artefacts specificity is
**unchanged** and false approvals are **unchanged at 2**. The safety posture did not loosen at
all; sensitivity simply more than doubled. **Every fold improved or held; none got worse.**

The earlier prediction came from a Platt fitted on calib alone and STACKED on the shipped one. The
actual recalibration fits on calib + the cross-donor pool from the raw probabilities -- a different
and, as it turns out, strictly better estimator. Predicting the wrong trade-off was the cost of
not having done this run first.

*(`_s12`'s baseline reads 0.275 / 66 where the original `_c7t` Stage 16 run read 0.297 / 64. Those
are different fold builds -- `_s12` is post-Stage-12 -- so the small difference is expected.)*

### Every pre-registered guard, checked on the scorecard

`compare c7t_stage12 c7t_stage16`:

| guard | result |
|---|---|
| `fate_prauc` | 0.958 -> 0.958, **diff exactly 0.000** |
| `fate_roc` | 0.954 -> 0.954, **diff exactly 0.000** |
| every ΔAge metric | **bit-identical** |
| `conformal_coverage`, `conformal_width`, `ood_rate` | **bit-identical** |
| RES | **still zero** (approvals 0 on every fold) |
| **`fate_ece`** | **0.276 -> 0.182, CI [-0.156, -0.033], 5/0\* unanimous -> ACCEPT (better)** |

**`fate_ece` is the only metric whose aggregate moved, and it moved the right way, decisively.**
Rank preservation is confirmed to 1e-12 on real data rather than argued from monotonicity.

Two details worth recording rather than rounding away:

* **Y1's `res_max` is floating-point DUST, not a score** -- 8.1e-11 before, 1.2e-12 after. That is
  the same residue behind the retracted "Spearman 0.40 over RES" headline, so the test asserts it
  as dust (`< 1e-9`, approvals 0) instead of asserting an exact zero that was never true.
* **`fate_ece_platt` moves per fold but its aggregate does not** (0.180 -> 0.180, CI includes 0).
  That is the expected signature of a second calibrator that had been partly compensating for the
  defect this fix removes.

### Status, stated precisely

**Stage 16 is deployed and empirically validated on the `_s16` folds.** The `train_model.py` fix
makes every FUTURE training run correct by construction, and `recalibrate_folds.py` brings
existing bundles up to date without retraining.

**It is not retroactively true of every artefact on disk.** `_c7t`, `_s12` and any other bundle
not passed through recalibration still carry soft-fitted coefficients and still behave as before.

---

## 2026-08-18 - STAGE 14 + STAGE 16 BUILT: calibrated dAge at the reporting boundary, and a fate calibrator fitted on the right target

**Status:** EXECUTED, both. Decisions taken by the user after the trade-offs were laid out.
`src/` changes: `inference/dage_calibration.py` (new), `inference/schema.py`,
`inference/service.py`, `training/train_model.py`. 30 tests. Full suite 1558 green, ruff clean.

The user asked for the plans to be improved and for diagnostics to run **before** building. That
sequencing changed Stage 16's fix entirely -- see the diagnosis-correction entry above.

### STAGE 14 -- k_var = 0.5991, at the reporting boundary, alongside raw

**Where:** `build_response` in `service.py`, and nowhere else. A test enumerates every `src/` file
importing `dage_calibration` and requires the list to be exactly `["service.py"]`.

**Why not upstream.** `res.py`'s `kappa = 5.0` is a rejuvenation half-saturation **in years**, so
calibrating before RES would silently reinterpret it -- the same class of defect as the
`huber_delta` knee that ruled out rescaling the training target. RES still sees raw
`mu_age`/`sigma_age`, and a test slices the `compute_res` call to prove it. The scorecard and
every evaluation path stay raw too, so **a 1/k drop in dAge MAE cannot appear anywhere** and be
mistaken for an improvement.

**Which factor.** `k_var = 0.5991` over `k_LS = 0.3637`. `k_LS` reaches MAE 6.78 against the
7.30 yr instrument floor -- the quotable headline -- but its SD ratio is 0.597: it wins partly by
under-reporting magnitude 40%. For a reporting transform the objective is unbiased magnitude, so
10.68 ships as the honest number. Both constants carry provenance and are pinned against the
recorded Stage 11 LODO fits so neither can drift.

**Alongside, never instead.** `delta_age_mean` / `_interval` / `epistemic_std` keep their raw
values. Four new fields (`delta_age_calibrated`, `_interval_calibrated`,
`epistemic_std_calibrated`, `delta_age_calibration_k`) default to `None`, so every existing
Response stays valid. `k` was fitted on O1/O2/O3 of the transient arm and that cohort is
**disjoint** from the training one, so a caller must be able to see both numbers. The caveat says
UNTESTED and DISJOINT in those words, and the pre-existing "do not read the number as years"
warning is **not** weakened.

**On the interval:** scaling both endpoints preserves the nominal coverage exactly, because the
RNA-clock label scales with them. It does **not** establish coverage against methylation truth,
which has never been measured. Stated in the caveat rather than implied by the units.

Reported dAge magnitudes now fall ~1.7x and interval width goes 63.4 -> 38.0 yr. **That is a
change of units, not an improvement**, and nothing in the scorecard will suggest otherwise.

### STAGE 16 -- fit the fate calibrator on the HARD class

`train_model.py` fitted Platt on `cal_target[:, SAFE_IDX]`, the **soft** probability in `y_cls`.
The calibrator was excellent at that (ECE 0.009-0.013 on calib) and wrong for every consumer,
all of which read `S` as P(hard class = safe).

The fix is one expression -- `_hard_safe(t) = (argmax(t) == SAFE_IDX)` -- applied to the shipped
fit **and** to the cross-donor diagnostic beside it, since two calibrators fitted against
different targets cannot be compared and that block exists to compare them.

**What was deliberately NOT done:**

* **No second calibrator stacked at inference.** That was the original 16.8 proposal. Because
  `platt_safe(p,a,b) = sigmoid(a*logit(p)+b)`, two Platts compose exactly into one
  (`a1*a2`, `a2*b1+b2`) -- so stacking would have worked numerically while leaving the root cause
  in place. A test pins the composition identity, and another pins that inference performs
  exactly one calibration step.
* **`tau_safe` was not moved.** The empirically optimal threshold on raw scores is 0.495 against
  a shipped 0.85, and lowering the bar to suit a soft-scale `S` would be fitting a safety policy
  to data. A test pins 0.85.

### FORWARD-ONLY, and this matters

Platt is fitted during training. **Bundles already on disk -- including the `_s12` folds built
today -- still carry soft-fitted coefficients.** The fix takes effect on the next training run;
it does not retro-correct anything, and a green suite must not be read as "the shipped folds are
now calibrated". Realising the measured gain (sensitivity 0.297 -> ~0.714) requires re-running
training, or a recalibration pass that re-fits the coefficients from the existing members.

### The safety posture, stated plainly

The user chose to make this the **default** rather than ship it behind a flag, with the trade in
front of them: the gate today catches **92.9%** of unsafe cells and approves **29.7%** of safe
ones; hard-calibrated it catches **82.1%** and approves **71.4%**. In counts, false rejections
64 -> 26 against false approvals 2 -> 5. That is a deliberate loosening of a safety gate, taken
knowingly, and it is recorded here so it is never mistaken for an incidental side effect.

Nothing about RES changes: it stays exactly 0, for the independent `R_eff` reason Stage 15
established.

---

## 2026-08-18 - STAGE 16 DIAGNOSIS CORRECTED: S is calibrated to the SOFT label; the gate wants a HARD one

**Status:** MEASURED, read-only. Found during the pre-build diagnostics the user asked for before
implementing 16.8. **`src/` untouched.** This **corrects the mechanism recorded in the Stage 16
entry of 2026-08-17**, which stands as written. The empirical findings there are unchanged; the
explanation is now wrong and is superseded here.

### What Stage 16 got right, and what it got wrong

**Right, and unchanged:** 91 of 119 held-out cells are truly safe; 70.3% of them fall below the
0.76 bar; shipped sensitivity is 0.297; a hard-label-fitted correction lifts it to 0.714.

**Wrong:** the diagnosis "H1 - plain miscalibration". Two things I did not check before concluding:

1. **Platt is ALREADY applied at inference.** `predictor.py:176` calls `apply_platt` whenever the
   bundle carries coefficients, and every fold ships them (`platt_a` 2.54-2.69, `platt_b`
   0.38-0.52). So the `S` Stage 16 measured was **already calibrated**, and its "deployable
   Platt" arm was a SECOND Platt stacked on the first. Because
   `platt_safe(p,a,b) = sigmoid(a*logit(p)+b)`, two Platts compose exactly into one
   (`a = a1*a2`, `b = a2*b1 + b2`) -- so that arm was really measuring *a different single
   calibrator*, not "calibration vs none".

2. **The shipped calibrator is fitted against the SOFT class target, not the hard label.**
   `train_model.py:207-215` fits on `cal_target[:, SAFE_IDX]` -- the soft distribution stored in
   `y_cls` -- while the safety gate, `scorecard.fate_ece`, and Stage 16's own analysis all score
   it against `argmax(y_cls)`.

### The measurement that settles it

On `calib`, per fold, shipped `S` against the two targets:

| fold | mean S | soft mean | hard mean | **ECE vs SOFT** | **ECE vs HARD** |
|---|---|---|---|---|---|
| N3 | 0.502 | 0.500 | 0.540 | **0.009** | 0.110 |
| O1 | 0.500 | 0.500 | 0.540 | **0.012** | 0.106 |
| O2 | 0.500 | 0.500 | 0.540 | **0.011** | 0.110 |
| Y1 | 0.501 | 0.500 | 0.540 | **0.009** | 0.109 |
| Y2 | 0.502 | 0.500 | 0.540 | **0.013** | 0.113 |

**The calibrator is doing its job almost perfectly -- against the soft target.** Mean `S` tracks
the soft mean to three decimals and misses the hard base rate by exactly the soft-vs-hard gap.

And the causal link to the rejections, on the held-out cells themselves:

> **46 of 91 hard-SAFE cells carry a SOFT safe target below the 0.76 bar.**

A perfectly soft-calibrated `S` is therefore *expected* to sit under the gate for **half** of the
genuinely safe cells. That is not a calibration failure; it is the gate and the estimate answering
two different questions.

### The corrected diagnosis

Not "the head is miscalibrated" (it is well calibrated), not prior/cohort shift (already refuted
by the oracle/deployable agreement). It is a **TARGET MISMATCH**:

* `S` estimates **P(soft label = safe)**.
* `tau_safe = 0.85` is a statement about **P(identity is actually preserved)** -- a hard event.
* The two differ by enough that half the truly-safe cells legitimately fall below the bar.

### What this changes about the fix

The 16.8 proposal -- "apply a Platt calibrator to `S` before the safety gate" -- would have shipped
a **second calibrator to compensate for the first being fitted against the wrong target**. It would
have worked numerically and hidden the root cause.

The correct fix is to fit the shipped calibrator on the **hard** label, so `S` means what every
consumer already assumes. Because Platts compose, this is expressible as **corrected coefficients
on a single calibrator**, not a stack.

Deliberately NOT chosen: moving `tau_safe` to suit a soft-scale `S`. That is fitting a safety
policy to data, and the Stage 16 plan ruled it out in advance (16.5).

### Open question this raises, recorded not answered

Which target is *right* for a safety gate is a genuine question, not a formality. The soft label is
the more honest description of an ambiguous cell; the hard label is what "identity preserved" means
operationally. This work adopts **hard**, because `tau_safe` is written as a statement about the
event, and because every existing consumer already reads `S` that way. The alternative -- keep `S`
soft and define the gate on the soft scale -- is coherent but would require re-deriving `tau_safe`
from evidence rather than inheriting it.

---

## 2026-08-18 - STAGE 17: the scorecard rewarded over-coverage, and hid whether the folds agreed

**Status:** EXECUTED. `plans/STAGE_17_COVERAGE_DIRECTION_AND_FOLD_CONSISTENCY.md`. Changes
`scorecard.py` (metric direction, aggregation, printing) only. **No model change, no rebuild, no
snapshot rewritten, no verdict rule changed.** 22 tests
(`tests/test_stage17_coverage_direction.py`).

Both defects were found while executing Stage 12 12.9 and deliberately left unpatched until that
measurement finished -- changing the instrument mid-measurement invalidates the measurement.

### D1 - `conformal_coverage` was judged "higher is better"

Coverage is **target-seeking**: it should approach `conformal_level` (0.90, stored per fold), not
climb. Coverage 1.000 is not a triumph, it means the intervals are too wide -- which is precisely
why `conformal_width` sits at **63-81 years**.

**Under the old direction, a change that simply widened every interval until nothing escaped would
have scored ACCEPT.** A test now constructs exactly that and requires it to read REGRESSION.

Not hypothetical: across the committed snapshots **4-5 of 6 folds are OVER-covering**. Mean
coverage looks respectable (0.85-0.92) while the mean *distance from nominal* is 0.10-0.18.

`conformal_coverage` is now direction `"target"`, judged on `|coverage - conformal_level|` per
fold. The target is read **per fold**, never hard-coded -- 0.90 is data. The raw coverage is kept
as a `(context, never judged)` row, because the distance alone cannot say **which side** of nominal
a fold sits on, and that is the actionable half.

**Retro pass: 4 of 21 distinct coverage verdicts change.**

### D2 - the table never showed whether the folds AGREED

`scorecard.py`'s own header says: *"check the per-fold column before trusting an aggregate
verdict."* **There was no per-fold column.** Only `dage_mae_model` ever printed one; the other 17
metrics printed a mean and a CI and nothing else. The instrument asked the reader to perform a
check it did not make available.

Every row now carries a `b/w` tally -- folds better / worse -- with `*` marking a unanimous
direction across >= 4 folds.

### The two fixes are load-bearing TOGETHER, which the retro pass shows

`baseline -> c7_A_keep_hff` flips `ACCEPT (better)` -> `noise` under the corrected statistic. Taken
alone that looks like a loss of information: coverage really did improve, 0.401 -> 0.923. But the
tally on the same row reads **5/0\*** -- unanimous. The old rule got that pair right *by accident*
(it was under-covering, so "higher" happened to align); the new statistic is correct but
underpowered at n=5; **the tally carries the signal neither one preserves alone.**

### CORRECTION owed to the Stage 12 entry (2026-08-18, which stands as written)

That entry recorded `conformal_width` **-6.97 yr [-19.93, +5.98]** as *"directional, NOT
significant, not claimed"*. **That remains true and the verdict is unchanged.** What the table
could not show at the time, and now can:

> the change was **unanimous across all five folds (5/0)**.

A 5/5 unanimous direction is a near miss, not a shrug. It is still **not significant** by the
pre-registered rule, Stage 12's verdict -- the pre-registered null -- is **unchanged**, and nothing
is claimed from it. What changes is only how strong the unclaimed directional hint was. Across all
21 distinct comparisons this is the **only** unanimous-but-noise width row.

### The guard that keeps the tally honest

A direction tally is a hair's breadth from a second significance test, and a second test on the
same data is how a project talks itself into findings the first test rejected.

- **It never produces a verdict.** `_verdict` still takes only `(direction, md, lo, hi)`; a test
  asserts that signature and that the tally cannot reach it. The accept/reject rule is unchanged.
- **It prints counts, never a p-value.** With 5 paired folds the smallest achievable two-sided
  sign-test p is **2/2^5 = 0.0625** -- a unanimous result can *never* clear 0.05 at this n.
  Printing it would invite the exact misreading the tally exists to prevent. Tests assert no
  p-value machinery is present.
- **It is symmetric**: unanimous regressions are flagged as loudly as unanimous improvements.
- **Ties are not agreement.** A row of ties is not unanimous, and unanimity needs >= 4 folds.

### Deliberately NOT changed

`conformal_width` stays `("lower", ...)`. Narrower **is** genuinely better at equal coverage, and
coverage is judged separately. Recorded as a decision, with a test pinning it, rather than left as
an oversight.

### A tripwire that fired

Stage 12's helper carried a test asserting `METRICS["conformal_coverage"][0] == "higher"` with the
note *"if the scorecard is ever fixed, this should be REVISITED rather than silently left
passing."* It fired on this Change. Revisited into a stronger assertion: the helper and the
scorecard are now two independent implementations of the same statistic and must **agree per fold
on real data**.

1528 tests pass, ruff clean.

---

## 2026-08-18 - STAGE 12 CLOSED: the rebuild ran, the split changed exactly as predicted, the metrics did not

**Status:** EXECUTED, §12.9 discharged, **Stage 12 CLOSED**. Six folds rebuilt + retrained (~5 h),
`scorecard/c7t_stage12.json`, `experiments/diag_stage12_rebuild_verdict.py`,
`results/diag_stage12_rebuild_verdict_results.json`. **No `src/` change in this entry.**

### The comparison is genuinely one change — verified before launch

`_c7t` was built 2026-08-15 02:01-05:59, **after** the gate fix `18b0c49` (01:30). `git log -- src/`
shows exactly two commits since 2026-08-14: that gate fix and Stage 12 (`72c0981`). So `src/`
differs between baseline and rebuild **by Stage 12 alone**.

Confirmed again by the run itself: the C-7 exact-match guard fired and passed — same 5 rejected Gill
samples, same 19 masked ΔAge labels, 42,600 cells, 42,581 labelled. **The labels are identical; only
the key differs.**

### Pre-flight before spending the compute (the arm-D lesson)

One fold, dataset-only, `MAX_CELLS=200`: cell_ids 1868/1868 unique, split map 1868/1868 entries,
guard silent, key format `reprogramming:HFF:b0:0`. Clear to launch.

### The split changed exactly as predicted — every number

`diag_stage12_split_effect` predicted these from artefacts alone, with no rebuild. Against the real
rebuild:

| quantity | predicted | actual |
|---|---|---|
| train / val / calib | 34,079 / 4,325 / 4,177 | **34,079 / 4,325 / 4,177** |
| D0 share train / val / calib | 11.8 / 11.7 / 11.4 % | **11.755 / 11.699 / 11.396 %** |
| distinct cell_ids | 42,600 | **42,600** (was 1,100) |
| split-map entries | 42,600 | **42,600** (was 1,100) |

Fixed as a side effect: `dataset_summary.json`'s `split_sizes` used to read
`{train: 850, calib: 115, val: 116, test: 19}` — **map entries, not cells**. It now reports cells.

### THE RESULT — the pre-registered null, on every metric

Target `|conformal_coverage − 0.90|` paired over 5 folds: **+0.0095, 95 % CI [−0.0169, +0.0360]**,
includes 0 → **NO DETECTABLE MOVE**, which is §12.9's pre-registered outcome #2.

**All 18 rows of `compare` read `noise (CI incl. 0)`.** Guards clean (`fate_prauc` −0.007,
`fate_roc` −0.013, `rank_model_dage` −0.001). RES still exactly 0, as Stage 15 predicted and
Stage 16 said in advance it would be.

§12.9 said what to do with this outcome: *"record it, do not re-run looking for a win."* Recorded.
**Stage 12 is CORRECT-BUT-INERT on model metrics.**

### The honest limit on that null

`conformal_coverage` is a fraction of ~20 held-out cells, so it is **quantised in steps of ~0.05
per fold**. Per fold: N3 19/20 → 19/20, O1 20/21 → 20/21, O2 20/20 → 20/20, Y1 18/18 → 18/18,
Y2 **15/21 → 14/21**. **The entire observed change is one cell in one fold**, and two folds sit
saturated at 1.000 where movement toward nominal is impossible without narrowing the intervals.

The null is real, but it is a null **at this resolution**. It must not be quoted as "Stage 12
provably changes nothing".

*Directional, NOT significant, not claimed:* interval width −6.97 yr [−19.93, +5.98]; OOD flag rate
0.451 → 0.343. Both point the way a better-composed calib set would; neither clears the bar.

### A metric the scorecard still judges the wrong way

§12.9's target needed a statistic `scorecard.py` does not compute. `conformal_coverage` is
registered `("higher", ...)`, but coverage is **target-seeking**: 1.000 is not better than 0.900, it
means the intervals are too wide — which is exactly why `conformal_width` sits at 63-81 years and
why O2 and Y1 both read 1.000. So `diag_stage12_rebuild_verdict.py` judges `|coverage − nominal|`,
the same shape Stage 13 established for `level_shift`. **The scorecard row itself is left alone and
this is logged as a known defect**, not silently patched mid-Change.

### Two defects found while running this

1. **`run_loocv.py` advertises a tag that would destroy the baseline.** It derives the snapshot tag
   from *arm + gate only*, so re-running the same arm prints the same tag as the first run — here,
   `c7_A_keep_hff`, the very comparator. Same class as the defect its own comment at lines 84-86
   records and fixes for the `gc2_` → `c7_` case, which does not cover re-running an arm. Recorded,
   not patched (the runner is not part of this Change); a background task is filed.
2. **The Stage 13 retro-audit's scope was a directory glob.** Writing `c7t_stage12.json` silently
   enlarged it and broke three pinned counts hours after Stage 13 shipped. Fixed properly rather
   than by updating the numbers: the retro pass audits comparisons **judged by the broken rule**,
   which is a closed historical set of nine snapshots. A post-Stage-13 snapshot was never judged by
   the old rule, so counting it would invent comparisons that never happened. The set is now frozen
   and two tests pin the scope.

### Status

**Stage 12 is CLOSED.** The defect was real and is fixed; the split composition changed exactly as
measured in advance; the model-metric effect is a null at the available resolution. The fix stands
on **correctness** — a unique key is unambiguously right — which is what §12.6 said in advance it
would have to stand on.

1506 tests pass, ruff clean.

---

## 2026-08-17 - STAGE 16: the safety gate rejects 70% of demonstrably safe cells, and calibration fixes it

**Status:** MEASURED, hypothesis CONFIRMED, Change pre-registered but **NOT executed**.
`plans/STAGE_16_SAFETY_FLOOR_MISCALIBRATION.md`, `experiments/diag_stage16_safety_floor.py`,
read-only, `results/diag_stage16_safety_floor_results.json`, 35 tests. **`src/` untouched.**

Stage 15 logged as an untested hypothesis that the `REJECTED_UNSAFE` rate looked like
miscalibration rather than danger. Tested now.

### The gate that could have killed it

A high rejection rate proves nothing if the cells are genuinely unsafe, so that was measured
first, with the power to end the stage. **It did not.** Pooled over six folds: **91 of 119
held-out cells are TRUE SAFE.** N2 is decisive — **19 of 19 truly safe, 15 rejected as unsafe.**
There is no reading of that fold in which the rejections are correct.

*(This gate was measured before the plan was written and is therefore not pre-registered; §16.3's
tests all were. Same convention as Stage 14 §14.3.)*

### The cost, and the repair

Median S is **0.704** for truly-safe cells and **0.217** for truly-unsafe — the head separates
cleanly, and its safe class simply sits below the 0.76 bar. **70.3 % of truly-safe cells fall
below it.**

| arm | false rej | drop | false appr | sens | spec | bal acc |
|---|---|---|---|---|---|---|
| raw (as shipped) | **64** | — | 2 | **0.297** | 0.929 | 0.613 |
| Platt on `calib` **[DEPLOYABLE]** | **26** | **−59.4 %** | 5 | **0.714** | 0.821 | 0.768 |
| Platt on `test` *[ORACLE, not deployable]* | 27 | −57.8 % | 0 | 0.703 | 1.000 | 0.852 |

Pre-registered bar was a ≥50 % drop in false rejections. **The deployable arm clears it at 59.4 %.**

**THE HEADLINE: shipped safety sensitivity is 0.297.** The gate approves under a third of
genuinely safe cells while the head behind it scores PR-AUC 0.965–0.992. **The fate head's quality
is not reaching the decision it exists to make.**

### H4 (prior/cohort shift) is REFUTED — and the deployable arm is the stronger one

The two-way fit was the point of the design: if calibration worked but did not transfer, the
oracle would repair it and the deployable arm would not. **The deployable arm matches the oracle
(26 vs 27 false rejections)** — calibration transfers from `calib` to a held-out donor.

Better than that: the oracle is **undefined on N2** (19/19 safe → unidentifiable Platt boundary →
passes through unchanged, which is why N2's oracle column equals its raw column). The deployable
calibrator is fitted on all six folds, because `calib` always carries both classes. **The
"deployable" arm is not a compromise; it is the better estimator here.**

### The trap was checked, not assumed

A calibrator could "fix" false rejections by shoving every probability upward. False approvals rise
only **2 → 5** while false rejections fall **64 → 26**, and balanced accuracy rises **0.613 →
0.768**. A favourable trade, not a shifted operating point.

### H3 also contributes, and is the same fact seen sideways

The oracle-best threshold on **raw** scores is **0.495** against the shipped **0.76** (balanced
accuracy 0.832). So "the probabilities are too low" and "the bar is too high" are one phenomenon.
Recalibrating the probabilities is the principled repair; **lowering the safety bar would be
fitting a safety policy to data**, and the plan rules it out in advance (§16.5).

### Post-hoc, explicitly not pre-registered

Training prior is **53.2 % safe**; test priors run 44.4–100 %. The false-rejection rate tracks the
gap between them: Spearman **+0.771**, n = 6, against a two-sided critical ρ of **0.886** —
**suggestive, NOT significant.** Y1 (prior closest to training, 44.4 %) has the lowest
false-rejection rate, 2 of 8; N2 (furthest, 100 %) among the highest, 15 of 19. Recorded as a
mechanism candidate, not a finding.

### What this does NOT buy

**RES stays 0.** Stage 15 established `R_eff = 0` for 119 of 119 cells, and that gate is
independent of the safety gate. This was stated in the plan before the run (§16.0) so no one reads
a 59 % repair as movement on the headline metric. What it buys is a safety verdict that can be
trusted — the fate head is the half of the project that is still alive.

### The Change this licenses (pre-registered, not executed)

Apply the `calib`-fitted Platt calibrator to `S` before the safety gate, behind a default-off flag.
Target: pooled sensitivity ≥ 0.60 (measured 0.714). Guards: false approvals ≤ 7; balanced accuracy
must rise; `fate_prauc`/`fate_roc` **unchanged** (Platt is monotone, so ranking is invariant by
construction — any movement means the implementation is wrong); RES must stay 0. `τ_safe` and `w`
are untouched — the floor is a policy, and this fixes the probabilities fed to it, not the floor.
Full text: `plans/STAGE_16_SAFETY_FLOOR_MISCALIBRATION.md` §16.8.

---

## 2026-08-17 - STAGE 15: RES is zero because the model has no confidence, and it is over-determined

**Status:** MEASURED, question CLOSED. `experiments/diag_stage15_res_zero.py`, read-only (inference
on ~20 held-out cells x 6 folds), `results/diag_stage15_res_zero_results.json`, 19 tests. **`src/`
untouched — this is a diagnosis, not a fix.**

`res_median`, `res_max` and `res_approvals` have read 0.000 in every snapshot ever taken and the
cause was never established. It is now, by attribution rather than inference.

### Which factor is zero — all four checked, on all 119 cells

`RES = φ(S) · S^k · g(R_eff) · exp(−λ·P_loss)`

| factor | result |
|---|---|
| `φ(S)` = sigmoid | **> 0 for all 119 cells** — a sigmoid has no zero |
| `S^k` | **> 0 for all 119** |
| `exp(−λ·P_loss)` | **exactly 1.000** — `lam: 0.0` ships, so this factor is inert |
| **`g(R_eff)`** | **zero for 119 of 119** |

So the zero is **`g(R_eff)` alone**, exhaustively, not "probably".

### Why R_eff = 0: σ dwarfs μ

`R_eff = max(0, −(μ_age + z_conf·σ_age))` — credit only for *confident* rejuvenation.

| fold | μ min | \|μ\| med | σ med | σ/\|μ\| | min(μ+zσ) | z needed |
|---|---|---|---|---|---|---|
| N2 | −19.85 | 14.75 | 44.46 | 3.01 | +10.72 | 0.416 |
| N3 | −19.62 | 15.10 | 30.04 | 1.99 | **+2.00** | **0.898** |
| O1 | −12.48 | 11.56 | 52.41 | 4.53 | +27.32 | 0.235 |
| O2 | −12.57 | 11.31 | 26.08 | 2.31 | +7.64 | 0.455 |
| Y1 | −19.17 | 19.04 | 82.57 | 4.34 | +5.40 | 0.350 |
| Y2 | −16.13 | 14.75 | 31.23 | 2.12 | +13.96 | 0.381 |

**The model's uncertainty is 2.0–4.5× its signal.** No cell in any fold comes within 2 years of
the credit threshold. **RES is not broken — it is working exactly as designed, and correctly
reporting that there is no confident rejuvenation to credit.**

The σ values are not spuriously inflated: conformal coverage is 0.714–1.000 against a nominal 0.90,
so the intervals are honest-to-conservative. The model genuinely does not know.

**Re-tuning the gate would not rescue this.** `z needed` — the largest `z_conf` at which *any*
cell would qualify — is 0.235–0.898 against a shipped 1.0. Even `z_conf = 0.9` would light up
**one cell in one fold**. This is not a near miss for the system.

### Second finding: RES ≡ 0 is OVER-DETERMINED — three gates, all closed

| gate | cells hit |
|---|---|
| `R_eff = 0` (no confident rejuvenation) | **119 of 119** |
| `REJECTED_UNSAFE` (S < τ−3w = 0.76) | 11–12 of ~20 in five of six folds |
| `REJECTED_OOD` | 6–7 per fold; **16 of 18 for Y1** |

**Fixing any one would still leave RES at zero.** This is structural, not a bug.

### Third finding: the status field understates the failure it reports

`compute_res_batch` applies precedence OOD → UNSAFE → NO_REJUVENATION, so only cells passing the
first two gates get labelled with the third. The status counts show **1–3 cells per fold** as
`REJECTED_NO_REJUVENATION` when in fact **100 % of cells** have `R_eff = 0`. A reader trusting
those counts would conclude rejuvenation was a minor issue. Pinned by a test so it cannot mislead
later.

### Observation, with its hypothesis marked as untested

Over half of held-out cells are `REJECTED_UNSAFE` in five of six folds, yet `fate_prauc` is
0.965–0.992 — the head **ranks** well while its probabilities sit low (median S 0.60–0.72 against
a 0.76 floor). That is consistent with **miscalibration** rather than genuinely unsafe cells, and
it connects to the pooled-ECE work. **Not tested here; recorded as a hypothesis.**

### What this closes and what it opens

**Closed:** why RES is zero. Fully attributed.
**Opened:** RES cannot become non-zero until σ_age falls below |μ_age| for at least some cells.
That is a statement about the model's uncertainty, not about the RES formula — and it is the same
wall as the ΔAge work: the signal does not exceed the noise. No change to `src/` is proposed,
because there is no defect here to fix.

---

## 2026-08-17 - STAGE 14 PRE-FLIGHT: rescaling the ΔAge target is NOT a units change, and I had that backwards

**Status:** MEASURED + PRE-REGISTERED. `experiments/diag_stage14_calibration_equivariance.py`,
read-only, trains nothing, `results/diag_stage14_calibration_equivariance_results.json`, 14 tests.
Plan: `plans/STAGE_14_ADOPT_CALIBRATED_DAGE.md`. **`src/` untouched. Nothing adopted yet.**

Stage 11 found the dense clock was never broken, only mis-scaled. The obvious next move is to
adopt `y_age -> k * y_age`. This ran the checks first, and one of them came back the opposite way
from what I predicted in writing.

### E1 — the linear path is exactly equivariant (as expected)

Ridge refitted on `k*y` gives predictions equal to `k` times the original to **1e-12** on all five
available folds: `MAE(k*y) = k*MAE(y)` exactly, `Δρ = 0.00e+00` exactly. **For the linear path,
calibration buys units and nothing else** — measured on real folds, not argued from algebra.

### E2 — the neural path is NOT, and my prediction was wrong

**I expected** `huber_delta = 2.0` yr to sit far below a ΔAge SD of 13-23, making essentially every
residual outside the knee, the loss effectively L1, and a rescale a near-pure units change.

**Measured** (ridge in-sample training residuals, 5 folds): median |residual| is **1.36-2.40 yr** —
*comparable to* the knee, not far beyond it. Fraction of residuals inside the quadratic region:

| fold | now | after x k |
|---|---|---|
| N3 | 44.9% | **87.2%** |
| O1 | 42.9% | **85.0%** |
| O2 | 43.2% | **85.4%** |
| Y1 | 66.6% | **97.3%** |
| Y2 | 42.9% | **85.1%** |

The loss today is a genuine Huber mix; after rescaling it becomes **85-97% quadratic — effectively
plain MSE.** So rescaling the training target silently converts a robust loss into an
outlier-sensitive one *at the same time* as changing units, and the two would be inseparable in the
result. **A two-change experiment wearing the costume of a one-change experiment.**

*Limitation stated in the plan:* the percentages come from ridge in-sample residuals standing in
for the network's; the exact values depend on the network's fit. The direction is monotone and a
test pins that, so only the magnitude is indicative.

### What this changes about the recommendation

Rescaling the **training target** is now **NOT RECOMMENDED**. It needs a rebuild + retrain + full
re-score, restarts the guard record, changes the loss regime as a side effect (and would require
rescaling `huber_delta` by the same `k` to avoid confounding), and buys nothing the cheap option
does not — because the model learns the same function up to scale.

**Recommended instead: calibrate at the REPORTING boundary.** Apply `k` to `mu_age`, the conformal
half-width, and any ΔAge quoted in years. No rebuild, no retrain, no guard-record restart, trivially
reversible, and it achieves the entire actual goal — absolute ΔAge claims and interval widths become
honest. Since ρ is rank-invariant, every ranking result stands unchanged and needs no re-audit.

### Two calls recorded in advance

**Which `k`: the variance-matched 0.599, not the least-squares 0.364.** `k_LS = ρ·SD(truth)/SD(pred)`
is strictly smaller and wins on MAE by *under-reporting magnitude* (SD ratio 0.597 — it shrinks the
spread 40%). That is the shrinkage trap this project already hit. For a reporting transform the
objective is unbiased magnitude, not minimal MAE, so 10.68 is the honest number and 6.78 is the
flattering one.

**A 63% drop in ΔAge MAE must never be reported as an improvement.** It is a change of units. The
`compare` table will print a large ACCEPT and that verdict is meaningless here. A guard table of the
exact per-fold values a pure rescale must produce is recorded in the plan (§14.5), with
`rank_model_dage` required to be bit-identical — any deviation means the change did more than
convert units.

### The transfer problem, stated rather than buried

`k` was fitted on donors **O1/O2/O3 of the transient arm only** — the only rows in
`results/dage_ledger.csv` carrying `TRUTH_meth_dage_mt` (68 of 90; the Sendai cohort carries none).
The cohort with methylation truth and the cohort the model trains on are **disjoint**. Stage 11
§11.4 forbade claiming transfer and that stands: calibrated ΔAge is to be reported **alongside**
raw, never silently in place of it, carrying the caveat that transfer to Sendai/HFF is untested.
§14.7 records what would license it.

---

## 2026-08-17 - STAGE 12 EFFECT: the split-composition half, measured exactly with no rebuild

**Status:** MEASURED. `experiments/diag_stage12_split_effect.py`, read-only,
`results/diag_stage12_split_effect_results.json`, 17 tests
(`tests/test_diag_stage12_split_effect.py`). **No rebuild, no retrain, no `src/` change.**

Stage 12 fixed the colliding `cell_id` and deliberately claimed **nothing** about the size of the
effect. That claim splits in two, and one half needs no compute at all:

- **the split map itself** — measurable now, exactly, from artefacts already on disk;
- **the effect on model metrics** — still needs a rebuild + retrain + re-score. **Not claimed
  here.**

A built fold stores everything the split depends on: `manifest.parquet` carries `cell_id` (the old
colliding key) alongside `shard_id` and `row_idx` (which together *are* the fixed key), in build
order. So both split maps can be re-derived with the real `holdout_split`, the real seed, and the
real rows.

### The canary, because a re-derivation is worth nothing unless it is faithful

The script reproduces the **old** map first and requires it to equal `splits/holdout.json`
**exactly**, aborting if not. It does: **1100 stored entries, 1100 rebuilt, identical.** So the
numbers below are measured against the build that actually ran, not simulated. (A test corrupts
the stored map and confirms the run aborts — a check that cannot fail is not a check.)

### The measurement

**42,600 cells carried 1,100 distinct ids.** D0 — the control anchor — is **4,988 cells sharing
117 ids**, so the control timepoint's split was decided **117 times, not 4,988.**

| split | n (old) | D0% (old) | n (new) | D0% (new) | shift |
|---|---|---|---|---|---|
| train | 33,686 | 11.9% | 34,079 | 11.8% | −0.1% |
| val | 4,405 | **13.3%** | 4,325 | 11.7% | −1.6% |
| calib | 4,490 | **9.0%** | 4,177 | **11.4%** | **+2.4%** |

**The finding:** under the colliding key the three splits scatter across **4.3 points**; under the
fixed key they converge to within **0.4 points** of each other and of the population rate (11.7%) —
a **>10× reduction in spread**. That is the signature of the sample size changing from ~117
decisions to 42,600, and it identifies the mechanism rather than merely correlating with it.

`calib` — the split conformal intervals are computed on, and the one most depleted of the control
timepoint — is the split that gains most from the fix (+2.4 points).

**33.6% of all cells (14,334) land in a different split** under the fixed key. The split map is
materially different, which is why the model-metric consequence is worth a rebuild rather than
being waved off as negligible.

### Correction to the Stage 12 record (the original entry stands as written)

Stage 12 recorded *"D0 occupies indices 0–111 exclusively"* and put the D0 decision count at
**112**. Measured across the union of all chunks it is **117** — 112 was the count in shard `b0`
alone, and chunks differ slightly in size, so D0 spans indices 0–116 over the union. The
load-bearing numbers are unchanged and reproduce exactly: **calib 9.0% vs val 13.3%**, and the
effective split size being ~10² rather than ~10⁴.

### Not claimed

The effect on any model metric. That requires the rebuild + retrain + re-score, pre-registered as
its own Change (see the annotation appended to `plans/STAGE_12_CELL_ID_UNIQUENESS.md`).

---

## 2026-08-17 - STAGE 13: the scorecard judged the wrong quantity, and it favoured the shuffle controls

**Status:** EXECUTED. `plans/STAGE_13_SCORECARD_VERDICT_CORRECTNESS.md`. Changes `scorecard.py`
**aggregation, verdict and printing only**. `measure_fold` untouched, **no snapshot rewritten, no
rebuild, no retrain, no model touched**. New: `experiments/diag_stage13_retro_verdicts.py`,
`results/diag_stage13_retro_verdicts_results.json`, 42 tests
(`tests/test_stage13_scorecard_verdicts.py`, `tests/test_diag_stage13_retro_verdicts.py`).

This discharges the two defects recorded on 2026-08-15 as *"Neither is fixed yet; both need their
own change"* (this file, §"Two scorecard defects found while reading the C-7 comparison"). **Both
of those notes stand as written and are not edited.** Re-measuring them from the snapshots found
that defect 1 has **three** faces, and that the note named the least important one.

### A1 — the aggregate cancelled instead of accumulating (NOT previously recorded)

`abs()` was applied to the **mean**, so for a per-donor bias whose sign varies by donor the column
measured how far the **panel cancels**, not how large the error is.

| snapshot | metric | printed | true `mean(\|.\|)` | understated |
|---|---|---|---|---|
| `gc2_A_keep_hff` | `level_shift_ridge` | **0.230** | **12.723** | **55.2x** |
| `gc2_B_mask_hff` | `level_shift_ridge` | 0.873 | 16.299 | 18.7x |
| `gc2_A_keep_hff` | `level_shift_model` | 5.713 | 13.120 | 2.3x |

**±12.7 yr per-donor level shift is the founding measurement of Stage 2** (`MASTER_PLAN.md:81,
387, 452`) and the whole justification for k≈3 reference cells. The scorecard printed that exact
quantity as **0.230** — i.e. as *"there is no level shift"*. The corrected statistic, 12.72,
reproduces the project's own independently-derived number from Test 7.4.3; the printed one erased
it.

### A2/A3 — the verdict was computed on signed differences

`_paired` ran before the display-time `abs()`, so the CI was built on signed values, and
`_verdict`'s `better_is_down` then read `-28 -> -22` (a 6 yr **improvement** in magnitude) as an
increase. On the C-7 comparison:

| metric | as printed | correct |
|---|---|---|
| `level_shift_model` | `+5.030` [+1.218,+8.842] **REGRESSION** | `-3.118` [-9.100,+2.865] **noise** |
| `level_shift_ridge` | `+4.389` [+2.805,+5.972] **REGRESSION** | `-0.084` [-5.399,+5.231] **noise** |

### B — the two columns could average different fold sets

`_agg` aggregated each snapshot over whatever folds were valid **in itself**, while `_paired` used
the intersection. N2 errors out in `c7_A` only, so **13 of 18 rows printed a 6-fold mean beside a
5-fold one**. `dage_mae_model` showed `14.291 -> 15.713` (a visible +1.42) with a verdict driven by
**+2.922** — more than double the visible number, with nothing on the row saying why.

### THE FINDING: the defect was biased, not merely noisy

Because the per-fold data on disk is intact, **every past comparison could be re-judged with no
rebuild** — the thing Stage 12 could not do. Over the 9 committed snapshots (5 distinct on these
metrics; `baseline == A_xdonor == B_fatecal == B_fatecal_pooled == gc2_A_keep_hff`):

**12 of 20 distinct verdicts change — and 8 of the 12 are comparisons in which a SHUFFLE CONTROL
scored `ACCEPT (better)`.**

    arm A -> arm C (full label shuffle)   level_shift_model   ACCEPT (better)  ->  REGRESSION
    arm A -> arm D (stratified shuffle)   level_shift_model   ACCEPT (better)  ->  noise

A shuffle destroys the donor structure the level shift measures; it cannot improve it. Mean
|level shift| in the raw data: **A 13.12 -> D 16.25 -> C 23.24**, monotone worsening exactly as a
negative control should. The old rule called that an improvement because the signed mean fell
(-5.71 -> -20.28) and `better_is_down` rewards falling. **The error had a systematic direction: it
flattered destroyed-label controls.**

A1 and A3 are **different mechanisms** and neither explains the other — A1 is about the printed
column (cancellation), A3 about the verdict (`better_is_down` on a signed quantity). A test
constructs a panel where A1 is inert (ratio 1.0) and A3 still inverts the verdict in both
directions.

Also recorded: **8 unchanged verdicts does not mean 8 sound rows.** 5 of the 8 had the *sign of
their point estimate flip* (e.g. `-5.388 -> +1.610`) and still read `noise` — the same verdict for
the opposite reason.

### The fix

- `"abs"` metrics aggregate as `mean(|per-fold|)` and pair on `|B_d| - |A_d|`. **Per-fold cells
  stay signed** — the direction of a donor's shift is real information, is read by
  `experiments/diag_zero_point.py:447`, and taking `abs()` at measurement time would destroy it
  permanently and make every snapshot unreadable in its own terms.
- The signed mean is **kept as its own row**, marked `(context, never judged)`. It answers a
  different, still-useful question (is there a *global* offset, or do donors cancel?), and
  replacing one number with another would have lost it.
- Comparison columns are averaged over the **paired fold set**, and dropped donors are named.
- The RES over-approval block uses one fold set for all four terms.
- `_print_snapshot`'s `mean` column for an `"abs"` metric is now `mean(|.|)` and is labelled `|.|`.

**The invariant that makes this testable:** `col_B - col_A == mean diff`, to floating point, for
**every** metric and direction. False for 13 of 18 rows before; holds for 18 of 18 after. For an
`"abs"` metric it can only hold if the columns and the paired statistic take the magnitude at the
same point in the computation — so one assertion pins both fixes at once.

### Why this is not a judgement call dressed as a fix

It implements intent already written in the source (`"judge |level shift|, not signed"`,
`scorecard.py:378`); it reproduces the project's own 12.7 yr number; and it is **not
directionally convenient** — it converts two `REGRESSION`s into `noise`, but it also raises the
reported level shift from 0.230 to 12.72 (worse), restores `dage_mae_model`'s true `+2.922`
degradation in place of a visible `+1.422` (worse), and turns the shuffle controls from
improvements into regressions.

### Not claimed

No model result changes — this stage re-judges, it does not re-measure. Which *decisions* in the
record were actually taken on a flipped verdict has not been audited; the retro-table is recorded
so that audit is now possible. The 2026-08-03 note (this file, *"A reader trusting that row would
draw the opposite conclusion"*) shows the row misled at least twice before it was fixed.

---

## 2026-08-16 - ALL NINE VARIANTS vs the instrument floor, with a shrinkage control

**Status:** Measured. `experiments/diag_instrument_floor.py`, read-only, ledger only. **`src/`
untouched, no label moved.** Bar amended **before** the run to require beating a constant-zero
predictor.

Floor MAE 7.30 (mt-sb). Truth SD mt 12.66 / sb 13.55. **Constant-zero control: MAE 11.71 (mt),
9.89 (sb).**

### Four findings

**1. The shipped dense clock LOSES to predicting nothing, on BOTH references** - 22.69 against
11.71 (mt) and 25.49 against 9.89 (sb), roughly 2x worse than outputting 0. §0 said this for one
clock and one statistic; it now holds on both, like-for-like.

**2. top100 is the only pass, and it is not an artefact** - MAE 7.15 vs floor 7.30 (CI spans 0),
**SD ratio 0.98**, **rho +0.81**.

**3. The shrinkage control earned its place on the first run.** `ranknorm` beats constant-zero on
**both** clocks and its delta-floor CI **spans zero on both** - on headline numbers, the best
all-rounder in the family. **SD ratio 0.30 / 0.28, rho +0.14.** A collapsed predictor. Without the
SD and rho diagnostics it would have been reported as a second candidate.

**4. THE SKIN-AND-BLOOD FAILURE IS FULLY EXPLAINED BY THAT CLOCK'S DISAGREEMENT WITH MULTI-TISSUE.**

| | Spearman |
|---|---|
| **mt vs sb (reference vs reference)** | **+0.613** |
| top100 vs mt | **+0.810** |
| top100 vs sb | +0.450 |
| predicted if sb reaches the RNA readout only through mt (0.810 x 0.613) | **+0.497** |
| observed | **+0.450**, gap **-0.046** |

**top100 orders multi-tissue (0.810) better than skin and blood orders multi-tissue (0.613)** - the
RNA readout agrees with one gold standard more closely than the two gold standards agree with each
other. **And its sb agreement is within 0.05 of pure mediation through mt, so no residual failure
remains to attribute to the RNA side.**

Family-wise: rho against mt spans **+0.14 to +0.81** (0.67 wide); rho against sb spans **+0.29 to
+0.48** (0.19 wide). **Nine variants ranging from excellent to useless on mt move sb's ordering by
0.19.** Whatever sb measures beyond mt is invariant to everything tried on the RNA side.

### What it licenses

**Not** narrowing the estimand - 1.5.2 refuted that on factor loadings (lambda_mt = 1.048 > 1), and
this **independently agrees with 1.5.2**: rho(mt,sb) = 0.613 beside rho(rna,mt) = 0.810 is exactly
what "three instruments not consistent with one common factor" looks like.

**Does** establish a bounded negative: **the sb-specific component is not accessible from RNA in
this dataset across the entire nine-variant family.** A measured limit, not a preference.

**§5.13's rejection stands.** What changes is the reason on record: not "the RNA readout is too
weak" but **"the two references agree with each other at rho 0.613, and the RNA tracks one at
0.810."**

---


## 2026-08-16 - THE LIKE-FOR-LIKE COMPARISON: top100 sits ON the instrument floor for multi-tissue

**Status:** Measured. `experiments/diag_instrument_floor.py`, read-only, reads only the ledger.
**`src/` untouched, no label moved.** Bar pre-registered in the docstring before the numbers.

44 non-control conditions carrying **both** methylation truths - every pairing on the **same rows,
same statistic, paired per condition**, which is exactly what §2 refused to assert without.

### The floor

**methylation vs methylation (mt - sb): MAE 7.30 yr, RMS 9.45**, mean signed -2.79. Methylation dAge
SD: mt 12.66, sb 13.55 - so **the floor is 54% of the truth's own spread**.

Corroboration: `dage_meth_concordance` got RMS 9.07 on 9 control groups; this gets 9.45 on 44
conditions. Different row sets, same answer.

### The result

| variant | vs | MAE | dMAE vs floor | 95% CI | |
|---|---|---:|---:|---|---|
| **raw - the SHIPPED dense clock** | mt | **22.69** | **+15.39** | [+12.00,+18.93] | FAIL |
| raw | sb | 25.49 | +18.19 | [+13.48,+22.83] | FAIL |
| **top100** | **mt** | **7.15** | **-0.16** | **[-2.81,+2.38]** | **PASS** |
| top100 | sb | 11.27 | +3.97 | [+1.78,+6.20] | FAIL |
| top500 / top2000 | both | 16.3-24.9 | +9 to +18 | | FAIL |

### Three things now established that were not

1. **The shipped dense clock is 3x outside the instrument floor** - 22.69 against 7.30. Measured on
   a principled scale rather than against zero.
2. **top100 sits ON the floor for multi-tissue** - 7.15 vs 7.30, **CI spanning zero**, statistically
   indistinguishable from the two gold standards' own disagreement. The strongest positive claim
   available about dAge in this project.
3. **The sb/mt SPLIT is quantified for the first time**: top100 misses skin & blood by **+3.97 yr
   beyond the floor, CI excluding zero**. §5.13's rejection of the sparse clock is vindicated with a
   number rather than a rule.

### What it settles about the 10% goal

10% of a 12.66 yr SD is **MAE <= 1.36 yr**. The floor is **7.30 yr - 5.4x larger**. **You cannot
verify agreement to 1.36 yr with two rulers that disagree by 7.30 yr.** Not hard; unverifiable with
these instruments, now measured rather than argued.

**Reachable restatement, half of it already met:** *dAge agreement at or inside methylation's own
self-disagreement.* Met on multi-tissue; not met on skin & blood.

### Not established

Two instruments agreeing does not make either correct - 1.5.2's factor-loading arithmetic
(lambda_mt = 1.048 > 1) showed the three are not jointly consistent with one age factor. It does not
rescue same-timepoint dAge PREDICTION, circular at rho 0.96. And it is **transient arm only, n = 44**
- no Sendai condition carries both truths, so `__POOLED__` is the same rows and is not an
independent check.

---


## 2026-08-16 - AUDIT-3: the work is sound, the model did not get worse, and "10% error" is mis-specified

**Status:** Read-only audit. **`src/` untouched, nothing withdrawn.**

### Three findings decide everything, and all three are his

**A. CIRCULARITY - the most valuable result of the arc.** Of the model's 2,000-gene panel, **1,956
carry Fleischer clock weights**, retaining **21.4%** of the clock's absolute weight mass. The
clock's own weights on that panel reconstruct the label at **rho 0.96-0.97**; ridge reproduces it at
**0.96-0.99** on all five evaluable folds. **The label is a linear functional of the input.** So the
5.84 yr MAE measures panel fidelity, not prediction. **A 10% target here is achievable and
meaningless** - add clock genes to the panel and it falls toward zero having learned nothing.

**B. "early dAge -> late dAge" is already dead.** Early->late Pearson 0.741, but **donor age ->
late is 0.931 (Spearman 0.971)**, and early->late **partialling out donor age is -0.064**. Donor age
does all the work and is known at t=0. The follow-up (early EXPRESSION -> late residual after age)
is **FRAGILE: 3 of 9, effectively 1 of 9**.

**C. Regime E fired P0, and my GSE165177 recommendation was WRONG.** I argued 7-8 samples per
timepoint would fix `p_unsafe`'s saturation. It cannot, for a structural reason: *"p_unsafe is a
fraction of CELLS; a bulk sample is already a population average, so a per-sample hard label
collapses the fraction to 0/1."* Measured `unsafe_sd_by_donor`: O1 0.10, O2 **0.00**, O3 **0.00**.
I confused replication at the SAMPLE level with resolution at the CELL level. **He is right.** What
GSE165177 does still fix - replicated contemporaneous controls - stands, and is why
`dage_gse165177` reproduced Gill 2022.

### The ceiling, and it is not the model

`dage_meth_concordance`: **inter-clock RMS between Horvath skin&blood and Horvath multi-tissue on
the same samples = 9.07 yr.** Typical HFF day-14 dAge after C-7 is ~-6.5 yr. **The two reference
instruments disagree with each other by more than the effect being measured.** A 10% target on
-6.5 yr is +/-0.65 yr. **No architecture, dataset or loss moves that.**

§1's RNA-vs-methylation MAE of **5.36 yr** should be read against 9.07, not against zero - on its
face the RNA readout agrees with one methylation clock better than the two clocks agree with each
other. **NOT yet established like-for-like** (MAE on 68 conditions vs RMS on 9 control groups -
§0's ERROR 2 was exactly this mistake), **but one cheap analysis from being established, and it is
the highest-value number left.**

### Did the model get worse? No.

`diag_target_path` shows median compression 0.826 -> 0.534 after C-7. That reads as worse and is
not: **C-7 removed a contaminant that inflated HFF's labels ~3x** (fold spread 16.67 -> 3.69 yr), so
a model partly fitting an artefact now has less artefact to fit. His own audit says it: *"N3's
'improvement' is not the model improving."* **Degrading against a contaminated baseline is the
expected sign of a real correction.**

### What to do - dAge STAYS; what changes is which claim is defensible

1. **Stop optimising same-timepoint dAge accuracy.** Circular at rho 0.96. Biggest change of
   direction available.
2. **Establish the instrument-floor comparison like-for-like** - RNA vs each methylation clock and
   clock vs clock, same conditions, same statistic, same pairing. If RNA sits inside the 9.07 yr
   envelope the honest headline is *"agreement at the limit of methylation's own reproducibility"* -
   a strong claim, not a scoped-down one, and possibly already true.
3. **Drop early-dAge -> late-dAge.** Partial -0.064.
4. **Keep what works:** `fate_roc` 0.983, within-donor ranking 0.925-0.983.

**On the 10% goal:** against the RNA clock's own output it is trivial (circular); against methylation
it is below the references' mutual disagreement. Neither is clearable by a correct system - which is
§5b applied, this once, to a project goal rather than a test. **The reachable version of the same
ambition: dAge agreement with methylation at or inside 9.07 yr RMS, demonstrated like-for-like on
both clocks.** §1's 5.36 suggests it may already be met on one. Item 2 settles it, cheaply.

---


## 2026-08-17 - STAGE 12: `cell_id` was not unique, and the split map was keyed on it. FIXED in src/

**Status:** Executed. Pre-registered in `plans/STAGE_12_CELL_ID_UNIQUENESS.md`.
**This stage CHANGES `src/`** -- the first in this arc that does. New:
`tests/test_stage12_cell_id_uniqueness.py` (11 tests). Suite green (1336), ruff clean.

### The defect

`cell_id` was `f"{source}:{cell_line}:{index_within_chunk}"` -- **the chunk was not in the key**.
HFF is planned as 45 chunks, so `reprogramming:HFF:0` existed 45 times: **42,481 cells carried 981
distinct ids**, and `splits/holdout.json` held **1,100 entries for 42,600 cells**.

`make_splits` keys the split on cell_id. So ONE index-slot decision was applied to all 45 shards.

### The harm -- measured, not assumed

Each shard holds all 9 timepoints, and within a shard **D0 occupies indices 0-111 exclusively**
(the other 8 days interleave from index 112). D0's split assignment was therefore decided by
**~112 index-slots**, not ~4,700 cells:

| split | share of D0 |
|---|---|
| calib | **9.0 %** |
| train | 11.9 % |
| val | **13.3 %** |

`calib` is depleted of the control timepoint by ~23 % relative to `val` -- and **calib is what the
conformal intervals are computed on**. The ±4 % spread is exactly the sampling noise of n=112,
which is what identifies the mechanism rather than merely correlating with it.

Effective split sample size: **981, not 42,481**. For D0: **112**.

### The fix

`CellChunk.id` is already guaranteed globally unique (`chunking.py:33` raises on a collision) and
was simply not used. Both construction sites now key on it:

    before  f"{source}:{cell_line}:{i}"   ->  reprogramming:HFF:0      (x45)
    after   f"{chunk_id}:{i}"             ->  reprogramming:HFF:b0:0   (unique)

**BOTH sites were wrong** (`sources.py:208` synthetic, `sources.py:363` reprogramming); fixing one
would have left the other colliding, and a test pins both.

### The guard

A build-time assertion before `make_splits` -- the consumer of the bad key -- raising with the
collision count. This defect survived because **nothing ever checked**.

**Build-time only, deliberately.** Folds already on disk carry colliding ids; a read-time assertion
would make every recorded artefact unloadable. Tests pin that the read paths carry no such check.

### What this does NOT do

**No rebuild, no re-score, no retrain. No recorded result is revised.** Existing folds remain
readable -- their shards and splits were written together and stay mutually consistent. The fix is
**forward-only**: it changes what the NEXT build writes.

**The size of the effect on model metrics is UNKNOWN and is not claimed.** Quantifying it needs a
rebuild under the fix and a paired comparison -- its own Change, separately pre-registered.

Validated end to end: four existing suites (`test_correctness`, `test_evaluation`, `test_inference`,
`test_tf_encoder`) call `build_dataset.run`, so the new guard executes on real builds. It does not
fire -- the ids are now unique.

---

## 2026-08-17 - STAGE 11: the dense clock was never broken -- it was MIS-SCALED. One parameter closes the gap

**Status:** Executed, READ-ONLY. Pre-registered in `plans/STAGE_11_DAGE_SCALE_CALIBRATION.md`
BEFORE the run. New: `experiments/diag_stage11_scale.py`,
`tests/test_diag_stage11_scale.py` (17 tests), `results/diag_stage11_scale_results.json`.
Suite green (1325). **`src/` untouched.**

### The result

`k` fitted LEAVE-ONE-DONOR-OUT against methylation ΔAge, 44 conditions, 3 donors.

| | uncalibrated | | least squares | | **variance-matched** | | |
|---|---|---|---|---|---|---|---|
| variant | MAE | SD | MAE | SD | **MAE** | **SD** | rho |
| raw | 22.69 | 1.66 | 6.78 | 0.60 | **10.68** | **0.99** | 0.770 |
| top100 | 7.15 | 0.98 | 6.17 | 0.77 | **7.50** | **0.99** | 0.810 |
| **top500** | 16.27 | 1.74 | 5.62 | 0.75 | **7.28** | **0.98** | 0.757 |
| top2000 | 21.62 | 1.78 | 6.45 | 0.65 | 9.41 | 0.99 | 0.771 |
| resid_pluri | 13.00 | 1.14 | 10.25 | 0.63 | 12.39 | 1.04 | 0.354 |

**Verdict: SCALE IS THE PROBLEM.** `raw` goes 22.69 -> 6.78 on a ONE-PARAMETER rescale, closing
**102 %** of its 15.5 yr gap to top100. `k` = 0.37 / 0.33 / 0.39 across donors, spread **1.19x**,
well inside the 2x transferability bar.

### What this reframes

**top100 was never a better instrument -- it was an accidentally well-scaled one.** Its
uncalibrated SD ratio is 0.98 while raw's is 1.66, which is the entire source of its apparent
15.5 yr advantage. Once both are calibrated, the sparse variants land 7.3-7.5 and raw lands 10.7:
sparsity still buys ~3 yr, not 15.5.

The dense clock's ORDERING was correct all along (rho 0.770). It reported the right story in the
wrong units, roughly 2.7x too large.

### The trade, stated rather than buried

Least squares reaches MAE 6.78 -- **below the 7.30 methylation floor** -- but at SD ratio **0.60**.
It gets there by SHRINKING: `k_LS = rho * SD(truth)/SD(pred)`, so with imperfect correlation it
deliberately under-shoots. That is the same trade that flattered `ranknorm` and `resid_pluri`.

A number reported as "these cells got N years younger" must therefore use the **variance-matched**
calibration, where raw reaches **10.68 at SD 0.99** and **top500 reaches 7.28 at SD 0.98 -- at the
methylation floor with magnitude preserved.** Any MAE below the floor should be treated as
shrinkage, not accuracy.

### Not settled

n = 3 donors for the LODO calibration; transient arm only (no Sendai condition carries both
methylation truths, so there is no independent cohort to confirm `k` on); whether `k` transfers to
a new dataset is untestable here and is not claimed. And a correct scale does NOT rescue
same-timepoint ΔAge PREDICTION, which is circular (rho 0.96-0.99) whatever the target's units.

`src/` is not changed by this stage. The most this licenses is PROPOSING a calibrated-ΔAge Change,
separately pre-registered (plan 11.0, same rule as Stage 10).

### ⚠️ CORRECTION APPENDED 2026-08-20 — the headline above overstates the result

*The entry above stands as written; this qualification is appended, not substituted.*

**The title claims "one parameter closes the gap". It does not close it.** The load-bearing
numbers, restated:

* Variance-preserving calibration takes raw MAE **22.69 → 10.68**. That is a large reduction and it
  is real.
* **It does not reach the 7.30 yr methylation-vs-methylation floor.** A gap of ~3.4 yr remains.
* The **6.78** figure that does beat the floor is the *least-squares* `k`, which wins partly by
  under-reporting magnitude (SD ratio 0.597). It is **not** the factor that shipped — Stage 14
  chose the variance-matched `k_var = 0.5991` precisely to avoid that shrinkage.

**The accurate summary is therefore:** *scale is a **major component** of dense-clock error, not the
whole of it.* Variance-preserving calibration reduces raw MAE 22.69 → 10.68 but does not reach the
7.30 methylation floor, and **transfer of the reporting calibration outside its methylation cohort
remains untested**.

The "Not settled" paragraph above already recorded the n = 3 / transient-arm-only / untested-transfer
caveats. What is corrected here is the **headline**, which read as if the gap had been closed.
`ARCHITECTURE.md` §13 has been updated to the wording above.


---

## 2026-08-17 - STAGE 10: pluripotency is MEDIATION, not contamination. The removal recommendation is WITHDRAWN

**Status:** Executed, READ-ONLY. Pre-registered in
`plans/STAGE_10_PLURIPOTENCY_CONTAMINATION_OR_MEDIATION.md` BEFORE the run. New:
`experiments/diag_stage10_pluri.py`, `tests/test_diag_stage10_pluri.py` (14 tests),
`results/diag_stage10_pluri_results.json`. Suite green (1304), ruff clean. **`src/` untouched.**

### The recommendation being withdrawn

The previous session concluded "pluripotency has to come out of the ΔAge readout" from a single
number (`resid_pluri` 22.69 -> 13.00 MAE vs methylation). **That was wrong**, and the objection
raised against it was right: OSKM INDUCES pluripotency, so pluripotency induction may BE the
mechanism of rejuvenation, in which case regressing it out deletes the signal.

### What is actually being removed (10.1)

The signature is 5 genes carrying **0.0005% of the clock's |w| mass**, ranked 17,829-26,826 of
33,155. **The clock barely reads them.** So `resid_pluri` does not remove a direct reading -- it
removes the component of ΔAge that CO-VARIES with a pluripotency score, a far more invasive
operation. Two of the five (`POU5F1`=OCT4, `SOX2`) are **OSKM transgenes**, so the score partly
measures vector dose. `top100` excludes all five by construction (all rank > 17,000), so truncation
and residualisation are NOT the same operation.

### Three tests, all MEDIATION

| test | result | reading |
|---|---|---|
| **A** transient-vs-control gap after residualising | shrinks **58%** (O1 0.67, O2 0.58, O3 0.46) | bar >50% -> **MEDIATION**: removing it destroys the outcome signal |
| **B** pluripotency among CONTROLS | signature is **exactly constant, sd = 0** | the genes are OFF without OSKM, so there is NO baseline covariation for contamination to act through -> **MEDIATION** |
| **C** agreement with methylation | rho **0.770 -> 0.354** | methylation cannot see RNA pluripotency, yet removing it HALVES agreement -> the component carried real signal -> **MEDIATION** |

**Verdict: MEDIATION, 3 of 3.**

Test B is worth stating plainly: the five genes are not expressed in untreated fibroblasts at all,
so the score has zero variance until OSKM is delivered. A covariate that does not exist without the
treatment cannot be a baseline confound. (The first version of Test B used ΔAge and returned NaN --
ΔAge is DEFINED relative to the controls, so controls carry ~zero ΔAge by construction. Fixed to
raw clock age; the NaN then turned out to be the answer rather than a failure.)

### Why the original number was misleading -- and it is the shrinkage trap again

| variant | MAE | SD ratio | rho vs methylation |
|---|---|---|---|
| raw | 22.69 | 1.66 | **0.770** |
| resid_pluri | 13.00 | 1.14 | **0.354** |
| ranknorm | 10.15 | 0.30 | 0.138 |
| **top100** | **7.15** | **0.98** | **0.810** |

`resid_pluri` improved MAE by **halving the ordering**. That is the same failure mode the
shrinkage control caught for `ranknorm` -- but it slipped through, because an SD ratio of 1.14 does
not LOOK collapsed. The tell is rho, not SD.

**So `raw`'s problem was never pluripotency. It was SCALE** -- right ordering (0.770), 66%
over-magnitude. `top100` fixes the scale (0.98) while IMPROVING the ordering (0.810), and it does so
without deleting the pluripotency-mediated signal. That is why top100 is the right variant and
`resid_pluri` is not.

### Standing

The earlier recommendation is withdrawn, not amended. No `src/` change is licensed by this stage
under any outcome (plan 10.5), and none was made.

**Not settled:** partial mediation cannot be excluded -- if ΔAge = a(true rejuvenation) +
b(spurious pluripotency) and the true effect is itself partly pluripotency-driven, all three tests
can read MEDIATION while some spurious component remains. And with OCT4/SOX2 in the signature,
"pluripotency" and "vector dose" are not separable in this data. n = 3 donors.

---

## 2026-08-16 - THREE TESTS (+ a Phase 0 nobody asked for): every one lands NEGATIVE, and two of my own hypotheses died

**Status:** Executed, READ-ONLY, pre-registered in `plans/THREE_TESTS_PREREG.md` BEFORE any run.
New: `experiments/diag_clock_difference_capacity.py`, `diag_phase1_top100.py`,
`diag_phase2_harmonized_transfer.py`, `diag_phase3_within_donor_forward.py`,
`tests/test_three_phases.py` (12 tests), four results JSONs. Suite green (1286), ruff clean.

### PHASE 0 -- which clock, and the two hypotheses it killed (both mine)

The brief was "redo everything on top100". That hides an assumption, and testing it changed what
every later phase used.

**Part A was INVALID and is recorded as such.** Scoring the shipped clock on GSE113957 gave
**MAE 0.13 yr, r = 1.000** against a published `cv_mae` of **12.27** -- a 94x gap. The clock was
FITTED on those 133 donors; that is memorisation, not capacity. The pre-registration flagged this
cohort as "optimistic for raw by construction" -- right in direction, wrong by two orders of
magnitude. Kept, not deleted: it quantifies how overfit the shipped dense clock is.

**Part B (refit inside CV folds, truncate the refit, score out of fold) -- MECHANISM NOT SUPPORTED.**
My hypothesis was that dense noise-weights cancel in absolute age but compound in a DIFFERENCE, so
sparsity should help ΔAge. It does not:

| level | MAE_abs | MAE_diff | (alpha=10) |
|---|---|---|---|
| full | **12.36** | **17.80** | reproduces the published cv_mae 12.27 -- the refit is sound |
| top2000 | 12.67 | 18.19 | |
| top100 | 20.24 | 27.73 | clearly worse on BOTH |

Absolute and difference winners are identical in 4/4 alphas. Two by-products worth keeping:
`MAE_diff ~ sqrt(2) x MAE_abs` throughout (17.80 vs 17.48 predicted), so **the earlier ΔAge noise
estimate was right** -- errors compound as independent; and the refit reproducing 12.27 validates
the procedure.

**Part C -- my replacement hypothesis died too.** "Truncation removes cohort-specific overfitting,
so it helps OUT of cohort." On Gill's day-0 fibroblasts (n=5, ages 0-53): raw rho **0.872**,
top2000 0.872, top500 0.616, **top100 0.616**. Truncation does not help out of cohort either.
(GSE165177 rho 0.000 at every level -- 38/53/53, no dynamic range, exactly as pre-registered.)

**What actually explains top100's instrument-floor win: the PERTURBATION.** From the ledger, 44
reprogramming conditions vs multi-tissue methylation:

| variant | MAE | SD ratio | removes |
|---|---|---|---|
| raw | 22.69 | 1.66 | nothing |
| resid_cc | 23.65 | 1.67 | cell cycle -- **no help** |
| **resid_pluri** | **13.00** | **1.14** | **pluripotency** |
| top2000 -> top500 -> top100 | 21.62 -> 16.27 -> **7.15** | 1.78 -> 1.74 -> **0.98** | progressively more |

Regressing out pluripotency alone nearly halves the error and removes most of the 66% magnitude
inflation; cell cycle does nothing. Combined with Parts B and C -- where truncation HURTS on resting
samples -- the mechanism is that **the dense clock reads the reprogramming programme on top of age**.

**Licensed rule, replacing "use top100 everywhere": dense for RESTING samples, top100 for PERTURBED
ones.**

### PHASE 1 -- both ΔAge forward tests on top100. UNCHANGED, and one gets WORSE

- **early -> late partial given donor age: -0.064 -> -0.443.** A large move, still short of the
  df=3 bar (|r| > 0.878). Verdict **UNCHANGED**, magnitude reported as the plan required.
  Notably the dense clock's headline correlations largely evaporate on the clean instrument
  (`early_cd13 ~ late`: +0.830 -> **-0.063**), consistent with them having been contamination.
- **late residual from early expression: SIGNAL -> NULL, 0/5 alphas** (LOO rho 0.657 vs null p95
  0.771). The 2026-08-15 "first positive result" was an artifact of the contaminated clock. The
  robustness sweep already said FRAGILE; a better instrument now kills the baseline outright. Two
  independent lines, same answer.

### PHASE 2 -- transfer through the project's own Harmonizer. DOES NOT WORK

Gill, n=5: best MAE **14.84**, Spearman **0.667**, **2/5 alphas** -- short of a majority.
Baseline (predict the training median) is 17.60.

**The Harmonizer is no better than the crude z-score** it was meant to replace (z-score got 3/5 at
the same 14.84 / 0.667). So the earlier "transductive correction" caveat was not the issue: the
project's own cross-cohort machinery does not rescue transfer either. GSE165177 excluded from the
verdict as pre-registered (no dynamic range).

### PHASE 3 -- the within-donor forward test. RAW SIGNAL, then NULL on control

GSE165177's design holds donor age constant by construction: 17 (donor, arm) trajectories over
d10->d17, 3 donor clusters. Precondition met (median |move| 3.60). Raw result: **LODO Spearman
0.684, beats the permutation null 5/5 alphas** -- and it also beats a **structure-preserving**
(within-donor) null 5/5, so it is not merely reading donor offsets.

**Then the persistence control killed it:**

```
spearman(early, late) directly     : +0.971
early SCORE alone, 1 feature, LODO : +0.968
full 1903-gene expression, LODO    : +0.684
```

The trajectory barely moves across the window, so "predict late from early" is predicting a nearly
unchanged number -- and **1,903 genes do WORSE than one number**. FINAL: **NULL**. Without that
control this would have been reported as a forward result with the donor-age confound removed.

### The honest summary of the four

Every phase landed negative, two of my own mechanistic hypotheses were refuted by my own tests, and
the one thing that survived is a **mechanism**: ΔAge under reprogramming is contaminated by
pluripotency, and both truncation and explicit residualisation reduce it. That is the finding
worth keeping.

---

## 2026-08-16 - FATE AUDIT: the head is PARTLY riding day, but not only day -- and my earlier claim was wrong

**Status:** Executed, READ-ONLY, computed from `scorecard/c7_A_keep_hff.json` (`_fate_S`/`_fate_y`)
and the built shards. No new build.

**Correcting my own claim first.** The 2026-08-15 entry stating "the fate label is a DAY THRESHOLD"
was generalised from ONE donor's shard (O1), which happens to have a clean split. It does not hold
across the cohort and should not have been stated as a general finding.

Tested properly -- if the label were a day threshold and `dose_time` is a model input, then DAY
ALONE should score as well as the model:

| fold | n | model AUC | day-only AUC | label separation |
|---|---|---|---|---|
| N3 | 20 | 1.000 | **0.760** | unsafe from d11 -- overlaps |
| O1 | 21 | 1.000 | **1.000** | clean split d29/d34 |
| O2 | 20 | 1.000 | **1.000** | clean split d29/d34 |
| Y1 | 18 | 0.838 | **0.669** | unsafe from d9 -- overlaps |
| Y2 | 21 | 1.000 | **0.969** | unsafe from d21 |

Mean model 0.968, mean day-only 0.880.

**Verdict: mixed, and materially better than the ΔAge head.** In O1 and O2 the label IS a clean day
threshold and the model adds nothing over a calendar. In N3, Y1 and Y2 the label genuinely overlaps
days and the model beats day-only by 0.24, 0.17 and 0.03. So the headline number is INFLATED by the
day correlation, but there is real separation underneath that a day threshold cannot produce.

**This is NOT the ΔAge circularity repeating.** There, the target is a linear readout of the input
at rho 0.96-0.99. Here the label is not recoverable from the input by construction, and the margin
over day-only is genuine signal.

Caveats: n=18-21 per fold, all Gill bulk; AUC 1.000 over ~20 points with ~4 positives is easy to
reach; and this pits a 33k-cell-trained model against a single feature.

---

## 2026-08-16 - OUT-OF-COHORT TRANSFER: NOT established. Raw transfer fails catastrophically

**Status:** Executed, READ-ONLY, PRE-REGISTERED. New: `experiments/diag_age_transfer.py`,
`tests/test_diag_age_transfer.py` (8 tests), `results/diag_age_transfer_results.json`.
Full suite green (1258), ruff clean.

**This qualifies the entry below.** The capacity result (MAE 11.95, n=133) was measured on the
cohort the clock was fitted on, and the open question was whether it TRANSFERS. Train on
GSE113957's 133 donors; predict the chronological age of day-0 fibroblasts in cohorts the model
never saw. Day-0 only -- a reprogramming sample's miss would be uninterpretable.

### Result

| cohort | n | ages | treatment | best MAE | Spearman | verdict |
|---|---|---|---|---|---|---|
| GSE165176 Gill | 5 | 0,29,35,53,53 | raw | 20.97 (worst 118.48) | -0.36 | **does not** |
| GSE165176 Gill | 5 | | zscore | **14.84** | +0.667 | TRANSFERS (3/5) |
| GSE165177 | 3 | 38,53,53 | raw | 12.86 | **-0.866** | does not |
| GSE165177 | 3 | | zscore | 10.48 | **-0.866** | does not |

**The solid finding is the failure.** RAW cross-cohort transfer is catastrophic -- MAE 70-118 yr
on Gill at low alpha, with NEGATIVE correlation. A batch effect between GSE113957 (raw counts) and
the Log2-RPM cohorts dominates the age signal completely.

### Why the one PASS should not be leaned on

The pre-registered rule fires for Gill/zscore, but three things weaken it:

1. **n=5**, and Spearman 0.667 needs 1.0 to reach p<0.05 at that n.
2. **MAE 14.84 against a baseline of 17.60** -- predicting the TRAINING median. A 16% improvement.
3. **The z-score is TRANSDUCTIVE.** It standardises features using all five test samples, so it
   needs the whole test cohort in hand. For a single new sample it is not available. Not a label
   leak, but not deployable either.

### GSE165177 is uninformative, not contrary evidence

Its donors span 38-53 yr -- a **15 yr range against a ~12 yr instrument error**. There is no
dynamic range to rank, and with n=3 and two tied ages Spearman is close to meaningless. The -0.866
should NOT be read as evidence against transfer; the cohort cannot test it either way.

### Standing answer to "we just need more donors"

**Not established.** Within a cohort the representation reads age well (n=133, MAE 11.95, r 0.846).
Across cohorts it does not, without a transductive correction whose supporting evidence is n=5 and
marginal. The requirement is therefore not simply "more donors" -- it is **more donors PLUS a
cross-cohort normalisation that works on a single sample**, and neither is demonstrated here.

Two silent-failure bugs found and pinned: GSE165177 ships its samples across two matrices and the
day-0 fibroblasts are only in `part2`; and its series matrix titles them `O1_Fib` while the
expression header says `O1 Fib`, so the age join missed and the whole cohort vanished to an empty
frame with no error.

---

## 2026-08-15 - AGE CAPACITY at n=133: the representation CAN learn age. n was the problem, not the method

> **[ANNOTATED 2026-08-16 -- original left standing.]** **QUALIFIED.** The transfer test above shows
> this result does NOT cross cohorts without a transductive correction: raw transfer to Gill's day-0
> fibroblasts gives MAE 70-118 yr with negative correlation. The capacity claim stands exactly as
> written -- the representation CAN carry age, measured within cohort with donors held out -- but the
> sentence "n was the problem, not the method" is too strong. Within-cohort n was *a* problem;
> cross-cohort batch effect is a second one, and it is not fixed by more donors.

**Status:** Executed, READ-ONLY, PRE-REGISTERED. New: `experiments/diag_age_capacity.py`,
`tests/test_diag_age_capacity.py` (9 tests), `results/diag_age_capacity_results.json`.
Full suite green (1242), ruff clean.

**This overturns a conclusion stated repeatedly earlier today.** Every negative result was reported
against 6 reprogramming donors, and the standing framing was "the constraint is n". That framing
was never tested, and it conflated two different questions. GSE113957 -- 143 dermal fibroblasts,
declared ages 1-96 -- was on disk the whole time and answers the capacity question at 24x the n.

### Result: CARRIES AGE, 5 of 5 alphas, both feature sets, both cohorts

| cohort | features | best MAE | ratio to baseline | Pearson |
|---|---|---|---|---|
| Normal, n=133 | pipeline panel | **11.95 yr** | 0.47 | **0.846** |
| Normal, n=133 | top-2000 HVG | 12.29 yr | 0.48 | 0.842 |
| GPL18573 only, n=120 | pipeline panel | 13.29 yr | 0.48 | 0.824 |
| GPL18573 only, n=120 | top-2000 HVG | 13.09 yr | 0.48 | 0.837 |

Mean-baseline MAE (predict the median age) is 25.37 yr. **The pipeline's own gene panel, with
plain ridge and donors held out, reaches MAE 11.95 yr against the published `cv_mae` of 12.27** --
matching a purpose-built clock. It holds on a single platform, so it is not a batch artifact, and
the pipeline's panel does as well as a panel fit locally.

### NOT circular, and that distinction is load-bearing

The target is the **GEO-declared chronological age of the donor** -- metadata that no transform of
the expression produced. `diag_clock_circularity` found the ΔAge regression predicts a linear
readout of its own input; that failure mode cannot arise here. Pinned by a test asserting the
script never imports a clock.

### The honest caveat

GSE113957 is the cohort the Fleischer clock was FITTED on, so age signal is guaranteed present.
This is therefore a CAPACITY result -- "the representation can carry age when age signal is there"
-- and NOT an out-of-cohort generalisation result. It also measures CHRONOLOGICAL age in resting
fibroblasts, which is a much larger signal than reprogramming-induced ΔAge.

HGPS (progeria, n=10) excluded from the primary: it ages abnormally fast and would inflate the
result. Its samples carry no parseable age in GEO, so they drop out regardless.

### What this changes

The day's failures were about (a) circularity in how ΔAge was POSED, and (b) n=6 for the
reprogramming question. **They were not about the representation being incapable.** With adequate
n and a non-circular target, this pipeline reads age at published-clock accuracy. The reprogramming
question remains limited by donor count -- but that is now demonstrably a data-collection problem
rather than a method problem, which is the distinction that was missing.

---

## 2026-08-15 - RNA/methylation concordance: direction agrees, magnitude does not (r ~ 0.87, n=6 cells)

**Status:** Computed from the recorded P1 join (`results/dage_meth_concordance_results.json`), no
new run. GSE165177 (RNA) x GSE165179 (methylation), the same donors, arms and days.

| clock | Pearson | Spearman | slope RNA~meth |
|---|---|---|---|
| Horvath skin & blood 2018 | +0.864 | +0.771 | 1.74 |
| Horvath multi-tissue 2013 | +0.869 | +0.657 | 1.98 |

An INDEPENDENT molecular assay agrees with the RNA ΔAge on direction and ordering, which is real
evidence the target is not pure noise.

**But the evidence is thinner than the correlation suggests.** The joined cells are 3 donors x 2
arms = 6 points, and the correlation is largely carried by one large between-arm separation
(transient intermediates ~ -65 yr vs transient fibroblasts ~ -25 yr in RNA). Within-donor Spearman
across arms is undefined -- only 2 arms per donor, and rank correlation needs 3. So this
establishes "both assays see intermediates rejuvenating more than fibroblasts", not quantitative
agreement at donor resolution.

The slope of ~1.7-2.0 means **RNA reports roughly twice the magnitude methylation does** -- a
calibration discrepancy consistent with the previously recorded "RNA 2.5-2.9x larger", and
unresolved.

---

## 2026-08-15 - ROBUSTNESS SWEEP: the positive result is FRAGILE (3 of 9, and really 1 of 9)

**Status:** Executed, READ-ONLY, post-hoc, PRE-REGISTERED bar. New:
`experiments/diag_residual_robustness.py`, `tests/test_diag_residual_robustness.py` (8 tests),
`results/diag_residual_robustness_results.json`. Full suite green (1229).

**This demotes the entry below.** The robustness checks named when that result was recorded were
run, and the result did not survive them.

### Result: FRAGILE

| run | features | median rho | alphas | verdict |
|---|---|---|---|---|
| BASELINE 7-29d, all markers, panel | 1903 | 0.943 | 5/5 | **SIGNAL** |
| window 7-15d | 1903 | 0.943 | 5/5 | SIGNAL |
| window 7-21d | 1903 | 0.943 | 5/5 | SIGNAL |
| window 11-29d | 1903 | 0.429 | 0/5 | null |
| marker CD13 only | 1903 | 0.257 | 0/5 | null |
| marker SSEA4 only | 1903 | 0.771 | 0/5 | null |
| features clock genes | 18928 | 0.429 | 0/5 | null |
| features top-500 variable | 500 | 0.771 | 1/5 | null |
| features random-500 | 500 | 0.029 | 0/5 | null |

3 of 9 against a pre-registered bar of 6. **And it is weaker than 3 of 9 looks:** the three passing
runs all return EXACTLY 0.943 -- the identical ranking of the same six donors. They are not three
independent confirmations, they are one configuration appearing three times. Effectively **one of
nine** configurations passes.

### What the failures say

- **Dropping days 7-9 kills it** (11-29d: 0.429). The whole effect lives in the two earliest
  timepoints. Meanwhile dropping days 21-29 changes nothing (7-15d identical). A signal that
  depends on two samples per donor is not a signal at n=6.
- **CD13-only collapses to 0.257** -- and CD13 was the marker the earlier clock-age correlation
  identified as most informative (+0.83). No contradiction, because that analysis predicted RAW
  late age (donor-age mediated) while this predicts the RESIDUAL, but it does mean the narrative
  "CD13 carries the signal" does NOT transfer to the residual task.
- **The clock's own genes do not carry it** (18,928 features, 0.429). Only the 2,000-HVG panel does.

### The one informative positive inside a negative verdict

**random-500 gives 0.029** -- essentially nothing. So IF there is any signal, it is specific to a
gene set and is NOT merely "donors whose early expression resembles each other end up alike." That
conditional is worth keeping; it just has nothing to condition on yet.

### Standing verdict

The headline result should NOT be relied on. It holds for one combination -- all markers, the HVG
panel, and days 7-9 included -- and vanishes under every single-axis change tested. At n=6, with
Spearman granularity so coarse that 0.943 is one adjacent transposition from perfect, that pattern
is what a lucky configuration looks like.

It is not affirmatively DISPROVEN either: a real effect could be carried by early timepoints and by
HVGs specifically, and n=6 has no power to tell those apart. The way to find out is unchanged and
is now the only route -- **more donors, ideally age-matched pairs**, which is where the residual
lives.

---

## 2026-08-15 - FIRST POSITIVE RESULT: early EXPRESSION predicts the late residual after donor age

**Status:** Executed, READ-ONLY, PRE-REGISTERED. New: `experiments/diag_residual_expression.py`,
`tests/test_diag_residual_expression.py` (11 tests), `results/diag_residual_expression_results.json`.
Full suite green (1217), ruff clean.

> **[ANNOTATED 2026-08-15, same day -- original left standing.]** **DEMOTED.** The three robustness
> checks this entry names as untested were run immediately afterwards (entry above) and the result
> **did not survive them**: 3 of 9 configurations against a pre-registered bar of 6, and the three
> passers return the identical ranking, so effectively 1 of 9. The effect requires all markers, the
> HVG panel, AND days 7-9; every single-axis change removes it. The measurement below is correct as
> executed and the procedure still passes its own controls -- what was wrong was calling it a
> "strong lead". It is one configuration out of nine at n=6. Not disproven, but not to be relied on.

### Result

Leave-one-donor-out, 6 donors, 1,903 panel genes, donor-age fit REFIT INSIDE EACH FOLD:

| run | LOO Spearman | null p95 | percentile | alphas passing |
|---|---|---|---|---|
| POSITIVE CONTROL: raw late age | 0.771 | 0.600 | 96.9 | 4 of 5 |
| **TEST: late residual \| donor age** | **0.943** | 0.829 | **98.6-98.9** | **5 of 5** |

**SIGNAL** on the pre-registered rule. Early expression predicts the late outcome after donor
chronological age is removed -- which the early clock-age SUMMARY could not do (partial -0.064,
entry below). The 2,000 dimensions the summary discards were carrying it.

The residuals being predicted, and why they are not donor age:

| donor | donor age | residual |
|---|---|---|
| N2 / **N3** | 0 / 0 | -0.90 / **+4.26** |
| O1 / **O2** | 53 / 53 | -0.47 / **+5.06** |
| Y1 / Y2 | 29 / 35 | -5.87 / -2.09 |

In BOTH age-matched pairs the same donor is higher. That is the earlier "2 of 2, sign test p=0.5,
meaningless but unconfounded" observation -- now recovered by a model that never saw the held-out
donor, in the one place donor age cannot be the explanation.

### Why the procedure is trusted

- **The positive control passes.** Without it a positive test would be uninterpretable.
- **The donor-age fit is refit per fold.** Residualising globally first would let a donor help
  define its own residual. Pinned by a test that constructs an exactly-age-determined target and
  requires NaN.
- **No analytic p-value.** At n=6 with 1,903 features none is credible; significance comes from a
  2,000-draw permutation null running the identical procedure.
- **The null captures the LOO artifacts.** At alpha=1e4 the control's null sits NEGATIVE (p95
  -0.657): with a mean-predicting model, leave-one-out makes the train mean anti-correlated with
  the held-out value by construction. The null reproduces that, so the comparison stays honest --
  and that is why the control "fails" at that one alpha.
- **All five alphas reported, none selected.** The effect holds across four orders of magnitude.

### A bug this found in its own instrument, and it is the RES bug again

`test_the_donor_age_fit_is_refit_per_fold_not_once_globally` failed on first run. With a target
exactly determined by age the residuals are ~1e-14 floating-point dust, and the `std == 0` guard
tests EXACT zero -- so it computed a Spearman over numerical residue. **Identically the defect that
turned RES values maxing at 1.6e-4 into "Ranking generalizes: Spearman 0.40."** Fixed with a
guard relative to the scale of the quantity the residual came from.

Also corrected before the run: the primal ridge solve (2000x2000 per fold, ~120,000 solves, does
not finish) was replaced by the equivalent dual form (5x5), pinned against an explicit primal
computation because it is a change to the MATH.

### What is NOT claimed

**n = 6.** Spearman 0.943 is one adjacent transposition from perfect; one donor moving one rank
changes it materially. Nominal p ~ 0.014, single dataset, after many tests today. This is a
STRONG LEAD, not an established result, and it needs donors that were not used to find it.

Untested robustness that should come next: a different early window, CD13-only vs all-marker
input, and a different feature set. If the effect is fragile to those, it is noise.

---

## 2026-08-15 - EARLY->LATE FORWARD TEST: the forward signal is DONOR CHRONOLOGICAL AGE (partial r = -0.064)

**Status:** Executed, READ-ONLY. New: `experiments/diag_early_late_forward.py`,
`tests/test_diag_early_late_forward.py` (14 tests), `results/diag_early_late_forward_results.json`.
Suite green, ruff clean.

### Why this was tried

The same-timepoint ΔAge regression is circular (entry below). An early->late formulation escapes
that by construction: the target is measured from a DIFFERENT, later sample, which the input
cannot contain. This was the user's proposal and it is the right shape for the project's goal --
read the early trajectory, call the late outcome.

### It escaped circularity, and then died on a second confound

On RAW clock ages (no control subtraction, so the shared-zero-point artifact cannot apply), with
N2 restored because a within-donor comparison needs no control, **n = 6**:

| relation | Pearson |
|---|---|
| early CD13 -> late | +0.830 |
| **donor age -> late** | **+0.931** |
| donor age -> early CD13 | +0.902 |
| **early CD13 -> late, donor age removed** | **-0.064** |

The partial is **-0.064**. Once donor chronological age is known, the early assay adds NOTHING.
The textbook mediation signature -- and `test_a_fully_mediating_covariate_drives_the_partial_to_zero`
simulates exactly this case and requires the partial below 0.2, so the instrument was validated
against the answer it returned.

The reverse partial, `donor_age ~ late | early_cd13` = +0.759, does not clear the df=3 threshold
(0.878) either; at n=6 nothing here reaches significance on its own.

### The timing hypothesis PASSED its pre-registered bar -- and is still not a result

Rule, written into the module BEFORE it was computed: `T_day` = the first measured day below the
donor's own early/late midpoint from which every later day stays below. "Sustained" is
load-bearing -- the early single samples swing 40-50 yr between adjacent days.

  N2 34, N3 34, Y1 34, Y2 40, O1 40, O2 47  |  **spearman(donor_age, T_day) = +0.906 > 0.886**

Three reasons it is a LEAD, all recorded before the run: the hypothesis was formed BY EYE after
seeing the trajectories (post-hoc); n=6; and `T_day` takes only 3 distinct values here, so one
donor moving one timepoint would move it materially. And it is still donor age doing the
predicting.

### The one donor-age-free scrap

Gill has two same-age pairs. In BOTH, the donor with the higher early CD13 age had the higher late
age: N2 76.1->42.2 vs N3 81.2->47.4 (both age 0); O1 109.0->65.8 vs O2 120.8->71.4 (both age 53).
2 of 2, right direction. Sign test p = 0.5 -- **statistically meaningless**, recorded only because
it is the sole unconfounded evidence available and it does not point the wrong way.

### The constraint this identifies

**Not a modelling problem.** With 6 donors spanning ages 0-53, donor age dominates the variance
and cannot be separated from reprogramming signal at df=3. No architecture changes that. The
design that would break it is **more donors at MATCHED chronological ages** -- hold age constant,
let outcome vary, which is what the two same-age pairs already gesture at. That is a sharper form
of `plans/DATA_REQUIREMENT_SECOND_TIMECOURSE.md`, currently ON HOLD.

Still untested and cheap: the late residual (after donor age) against early **expression** rather
than the early clock-age summary, which discards 2,000 dimensions. Power at n=6 is nearly nil.

---

## 2026-08-15 - CIRCULARITY: the ΔAge regression predicts a linear readout of its own input

**Status:** Executed, READ-ONLY. New: `experiments/diag_clock_circularity.py`,
`tests/test_diag_clock_circularity.py` (16 tests), `results/diag_clock_circularity_results.json`.

### The result: CIRCULAR on all five evaluable C-7 folds

`clock_panel(X) = X @ w` -- the Fleischer clock's OWN weights applied to the model's OWN 2,000-gene
input.

| fold | T1 clock->label | T2 ridge->clock | T3 ridge->label |
|---|---|---|---|
| N3 | 0.9643 | 0.9875 | 0.9653 |
| O1 | 0.9743 | 0.9788 | 0.9898 |
| O2 | 0.9704 | 0.9825 | 0.9945 |
| Y1 | 0.9723 | 0.9603 | 0.9930 |
| Y2 | 0.9596 | 0.9817 | 0.9938 |

The label is reconstructable from the model's input by the clock at rho ~ 0.96-0.97, and ridge
reproduces that readout at rho ~ 0.96-0.99. `xonly~clk` matches `T2` to four decimals, consistent
with the separate finding that the perturbation contributes nothing.

**Ridge is not predicting age. It is re-deriving a linear functional of the vector it was handed.**
Its good MAE (5.84 yr on Y1 against a target SD of 27.5) measures how well the panel preserves the
clock's direction -- a compression diagnostic, not a prediction. There is no forward step in it.

### Why the verdict is trustworthy

- **The test can say no.** Pre-C-7 it returned NOT CIRCULAR (N3) and LABEL-RECOVERABLE (O1, O2, Y1).
- **Not circular by construction**, checked before running: the label comes from the FULL ~33k-gene
  profile (`build_dataset.py:342-346`), the panel retains only **21.35%** of the clock's |w| mass,
  and `X` is harmonized while ΔAge is deconfounded and re-anchored. Any of those could have broken
  the correspondence; none did.
- **Consistency carries it, not one fold.** At n~20 a single rho of 0.97 has a 95% CI reaching near
  0.92; five of five in 0.96-0.97 is the evidence.

### Uncomfortable: C-7 made the labels cleaner AND more circular

Pre-C-7 was mixed (1 not circular, 3 label-recoverable, 2 circular); after removing the degenerate
control all five are circular. Part of what made the task look non-trivial was the corruption.

### Alternative reading NOT excluded

The label may track a dominant expression axis rather than the clock specifically. It does not
change the verdict -- either way the target is determined by the input -- but it would change the
mechanism. Separating them needs a PC decomposition.

### Scope

Does NOT invalidate the clock, and does NOT touch C-7, which rests on the degenerate column.
It DOES invalidate every ΔAge MAE in this project as a measure of predictive skill, including the
C-7 vs pre-C-7 comparison recorded below.

---

## 2026-08-15 - RES IS DEGENERATE, the fate label is a DAY THRESHOLD, and two scorecard defects

**Status:** Measured from the frozen snapshots and built folds. READ-ONLY. Recorded here because
each retracts or reframes a previously reported number.

### RES is identically zero, and the "Spearman 0.40" headline was numerical residue

Per-fold `res_max` from the snapshots:

- **`c7_A_keep_hff`: exactly 0.0 on all five folds.**
- **`gc2_A_keep_hff`: 1.6e-4 (N2), 7.6e-5 (O1), 1.2e-6 (Y2), exactly 0.0 on N3/O2/Y1.**

The three folds with a non-zero `res_max` are EXACTLY the three that produced a Spearman. So the
pre-C-7 headline -- *"Ranking generalizes across held-out donors: Spearman 0.40 +/- 0.20 (n=3)"* --
was a rank correlation over RES values whose maximum was 1.6e-4, 7.6e-5 and 1.2e-6. That is
floating-point residue, not signal.

**Now believed:** RES was effectively constant-zero all along; C-7 removed the last crumbs that let
the metric compute at all. `res_approvals` and `res_approvals_oracle` are 0 in BOTH snapshots --
the gate approves nothing, ever. Independently corroborates Stage 3a's "the estimator is not a
detector". The earlier claim is left standing in its own entry, annotated.

### The fate head's PR-AUC 1.0 is a day threshold

Per-donor fate labels against reprogramming day (Gill bulk, O1 shown; the identical zero-positions
12-15 across O1/O2/Y2 show the same structure in every donor):

  days 0, 7, 9, 11, 13, 15, 21, 29 -> **safe** (all 17 samples)
  days 34, 40, 47, 54              -> **loss** (all 4 samples)

A clean step function with no overlap -- and `dose_time`, which encodes day, is a direct model
INPUT. PR-AUC 1.0 is obtainable by learning "day > 31 -> unsafe". The label is biologically
sensible (later reprogramming = more identity loss) but statistically redundant with a feature the
model is handed. **Not a prediction**: knowing it is day 54 of OSKM requires no model.

### Two scorecard defects found while reading the C-7 comparison

1. **Signed metrics are judged lower-is-better.** `level_shift` ranges -28.3 to +15.0 and the
   scorecard averages the SIGNED values, so C-7 moving all five folds by a consistent +5.03 --
   TOWARD zero, since most started negative -- was flagged `REGRESSION`. Mean |level shift| actually
   IMPROVED 12.74 -> 9.62 (ridge: 11.25 -> 11.16, flat). **Neither was a regression.**
2. **Aggregate columns can span different fold counts.** `gc2` ΔAge MAE averages 6 folds, `c7`
   averages 5 (N2 has no ΔAge), so "14.291 -> 15.713" is not a before/after pair. Verified by
   arithmetic; the paired diff (+2.922 over the 5 common folds) is the correct number and the
   table reports it correctly.

Neither is fixed yet; both need their own change.

---

## 2026-08-15 - TARGET-PATH AUDIT: nothing normalises ΔAge, and N3's "improvement" is not the model improving

**Status:** Executed, READ-ONLY (inference on ~20 cells/fold). New:
`experiments/diag_target_path.py`, `tests/test_diag_target_path.py` (14 tests),
`results/diag_target_path_results.json`. Full suite green (1163).

> **[ANNOTATED 2026-08-15, later the same day -- original left standing.]** The scale-mismatch
> result below is not retracted, but its IMPORTANCE is now much smaller than stated. The
> circularity entry above shows the ΔAge regression predicts a linear readout of its own input
> (rho 0.96-0.99), so "compression of the prediction against the truth" is a statement about how
> faithfully a linear model reproduces a clock readout across a scale change -- not about
> predictive skill. The H-SUPPORTED verdict stands on its own numbers; what it MEANS is narrower
> than this entry implies.

### The target path, read from source

- `scalers.json` holds `x_mean`, `x_std`, `dt_mean`, `dt_std`, `proliferation_coef` -- and **no
  age scaler**.
- `training/dataset.py:58`: `ya = np.where(am, arr["y_age"], 0.0)` -- **raw years**, while `X` and
  `dose_time` both go through `scalers.transform_*`.
- `models/losses.py:58`: `F.huber_loss(age_pred[m], age_true[m], delta=2.0)` -- an **unweighted**
  mean over masked cells, `delta` fixed in **years**.

So there is no target normalisation and no per-source weighting. **Whatever scale HFF's labels
have is the scale the age head learns**, and HFF is 99.8% of the age-valid training cells
(33,613 vs 60-77 Gill). A second consequence: `delta=2.0` is absolute, so halving HFF's spread
moved delta from 0.105 to 0.228 of the target SD -- the loss's shape relative to the data changed
without anyone touching the loss.

### Scale-mismatch test: H-SUPPORTED by the pre-registered rule, but read the caveat

compression = SD(pred)/SD(true) on the held-out donor. Median **0.826 -> 0.534**; more compressed
under C-7 in **3 of 5** folds (O1 0.826->0.437, O2 0.880->0.428, Y2 0.905->0.534); N3 (0.609->0.687)
and Y1 (0.801->0.825) went the other way.

The rule fired, but its majority clause is nearly powerless: with 5 folds, P(>=3 in one direction |
no effect) = 0.5, a coin flip. **The magnitude carries this, not the count**, and the effect is
heterogeneous. That caveat is stated in the code beside the rule, not appended after the result.

### The finding that changes an earlier reading

Decomposing MAE into bias (|mean_pred - mean_true|) and spread:

| fold | arm | bias | MAE | bias/MAE |
|---|---|---|---|---|
| N3 | pre-C-7 | -29.69 | 29.69 | **1.00** |
| N3 | C-7 | -21.29 | 21.29 | **1.00** |
| O1 | C-7 | -0.39 | 12.73 | 0.03 |
| O2 | C-7 | +1.37 | 13.78 | 0.10 |

**N3's MAE is 100% bias in BOTH arms.** The model predicts ~0 (0.372, then 0.665) against a truth
of +30.07, then +21.96. So N3's headline "improvement" (MAE 29.70 -> 21.29, the one fold the
scorecard marked `+ better`) is **not the model improving** -- it is the target moving 8.1 yr
closer to a flat near-zero prediction the model makes regardless. The model carries no ΔAge signal
on N3 in either arm.

Conversely O1 and O2 under C-7 have almost no bias, so their doubled MAE is **spread** -- which is
where the compression result actually bites.

### Also found, NOT caused by C-7, not yet assessed

`splits/holdout.json` has **1,100 entries for 42,600 cells**. Splits are keyed on `cell_id`, which
indexes WITHIN a chunk, so all ~45 shards share one index->split map and train/val/calib is
assigned over ~981 HFF index-slots rather than 42,481 cells. Proportions come out right (79.1%
train), which is why it is invisible. Whether it biases val/calib depends on whether position
within a chunk correlates with timepoint -- **untested, and it needs its own owner.**

### Not claimed

That normalising the target would make the model correct, or that the pre-C-7 target was the right
one. C-7's justification is the degenerate control and none of this bears on it.

---

## 2026-08-15 - PAIRED TARGET AUDIT: C-7 changed HFF's ΔAge target NONLINEARLY, and the N2 fold proves the cause

**Status:** Executed, READ-ONLY, no rebuild. New: `experiments/diag_target_shift.py`,
`tests/test_diag_target_shift.py` (17 tests), `results/diag_target_shift_results.json`.
Full suite green (1145), ruff clean.

> **[ANNOTATED 2026-08-15, later the same day -- original left standing.]** Everything measured
> here about the LABEL is unaffected: the target genuinely moved nonlinearly, HFF's spread genuinely
> halved, and the N2-fold internal control genuinely pins the cause to one column. What has changed
> is the framing of the MODEL comparison this was serving. The circularity entry above shows the
> ΔAge regression is a linear readout of its own input, so the "the model got worse / ridge got
> better" question that motivated this audit was never measuring predictive skill on either side.
> The label findings stand; the motivation was mis-posed.

### The question

The C-7 retrain's ΔAge MAE worsened in 4 of 5 folds while ridge's improved. But C-7 changes the
LABELS, so pre/post MAE are errors against two different targets. This pairs the label sets
directly, per cell, within fold. Readings A/OFFSET, B/SCALE, C/NONLINEAR were **pre-registered**
as module constants before running (`SLOPE_TOL=0.05`, `R2_FLOOR=0.90`, `TIME_TOL=2.0`).

### Result: C, and not marginally

HFF (42,481 cells/fold), pre-C-7 -> C-7: **slope 0.26-0.44** (1.0 would be a pure offset),
**r2 0.50-0.80** (below the 0.90 floor), and the per-timepoint mean shift **spreads 22.8-28.0 yr**.
The shift ramps monotonically along the trajectory: ~0 at D0, ~-2 at D2/D4, +2 to +6 by D10/D12,
**+16 at D14, +20 to +25 at D21**. Independent corroboration: step 3c recorded HFF day-14 moving
+18.558 yr when the degenerate control leaves; this audit measures +15.6 to +17.5 at D14 by a
different route.

**The target's spread halved.** HFF SD 19.13 -> 8.77 (excluding the N2 fold). Gill's did NOT:
25.78 -> 27.65, slightly UP, and Gill's own change is nearly a pure offset (slope 1.03-1.09,
r2 0.98, shift -6 to -7 yr).

### The N2 fold is an internal control, and it attributes the cause exactly

Fold N2 shows **slope 1.000, r2 1.000, mean shift 0.012 +/- 0.022 yr, timepoint spread 0.02 yr** --
no target change at all. In that fold N2 is the HELD-OUT donor, so its degenerate control was
already excluded from the harmonizer pre-C-7 (run log: "1 control(s) excluded" before, "0" after,
the gate having already removed it). The other four rejected columns are treatment samples, and
the harmonizer fits on CONTROLS only.

So the entire target change is attributable to **one column** -- `N2_Fib_Sendai_Exp2` leaving the
harmonizer's control set. Not to the other four rejections, not to anything else C-7 does.

### Consequence for the model comparison

Train/test scale ratio (Gill SD / HFF SD) went **1.35 -> 3.15**. The network now trains on a
target with half the dynamic range and is evaluated on one with the same range as before, so
compressed predictions carry a guaranteed MAE penalty. That is a stated mechanism consistent with
the observed direction, **not** a demonstrated cause -- separating it needs its own test.

### A diagnostic bug, caught by its own canary

The first version joined on `cell_id`, which is NOT a key: it indexes WITHIN a chunk, so
`reprogramming:HFF:0` occurs 45 times (1,100 unique ids over 42,600 rows). The join exploded
42,600 -> 1,843,299 rows and produced a plausible-looking table of pure noise (slope 0.011-0.034,
r2 0.001). It was caught by a check built in for free -- `_c7` vs `_c7t` labels, which are the same
config on the same data and must differ by exactly 0, read **73.77**. Those numbers were discarded
and never reported as a finding.

Corrected pairing: positional within each (shard, timepoint) group, and ONLY where the group has
equal size in both builds -- a Gill donor that lost a column has every later index shifted, so
those 9 cells are dropped rather than guessed at. Every pair is then verified on `cell_line` and
`time`, which the harmonizer never touches, and a mismatch raises. `_c7` vs `_c7t` now reads
**0.000e+00** in all six folds, which is also runbook §3 check 3 discharged.

### Not claimed

Which target is scientifically correct. C-7's justification is the degenerate control itself
(`integrity.py`: library 1.03e8 vs the 1e6 that RPM means, dynamic range 1.74 log2 vs 9.00-15.26
for every admitted column), and nothing measured here bears on it. A worse fit to cleaner labels
is not an argument for dirtier labels.

---

## 2026-08-15 - C-7 RETRAIN RUN 1 INVALID: the gate reported ON and did nothing (third time)

**Status:** Run 1 executed, VOIDED, root cause found and fixed, full suite green (1041 tests).
Re-run pending. `src/cellfate/data/sources.py` and `local_runners/run_multi_local.py` changed;
`tests/test_c7_reaches_the_retrain.py` extended to 17 tests.

### What happened

The six-fold arm-A retrain ran to completion with `CELLFATE_BULK_GATE=1`. The header printed
`[C-7] bulk_integrity_gate = ON`. **The labels it produced were pre-C-7.**

| | required | run 1 produced |
|---|---|---|
| cells | 42,600 | **42,605** |
| ΔAge labels masked | 19, all N2, `no_control_baseline` | **0** |
| Gill columns kept | 119 | **124** |

`n_age_labeled` came out **equal to** `n_samples` in all six folds. That is the unambiguous
signature of "nothing masked", independent of how cells are counted — the specific alternative
reading that was raised and checked before the run was voided.

### Root cause - the same defect class, a third distinct incarnation

1. (`e6fc183`) the flag was wired into `build_sources`, which `run` skips when sources are injected.
2. (`ebb2a95`) `run_multi_local.py` built its `DataConfig` without the field at all.
3. **This one:** the flag arrived correctly and was applied to the right sources — **too late.**

`run_multi_local.py` called `gill.plan()` to list donors ~30 lines before the `DataConfig` existed.
`plan()` reads and **caches** the matrix, and the gate's screen lives inside that read (`_load`,
which early-returns on its cache). `apply_source_flags` then set a flag that read `True` over an
already-cached, unscreened 124-column matrix.

Confirmed directly rather than inferred: constructing the source and setting the flag **before**
`plan()` rejects 5 columns; setting it **after** rejects 0. Same source, same file, same flag value.

The five columns are `N2_Fib_Sendai_Exp2`, `N2_d21_CD13_Sendai_Exp2`, `N3_d21_SSEA4_Sendai_Exp2`,
`O2_d9_SSEA4_Sendai_Exp1`, `Y1_d7_CD13_Sendai_Exp1`. The first is **N2's only control** — losing it
is what leaves N2 without a baseline and masks its remaining 19 ΔAge labels under rule 4.

### Fixes

- `bulk_integrity_gate` is now a **property** on `GillReprogrammingSource`. Its setter drops the
  cached read when the value **changes**, and is idempotent when it does not (`apply_source_flags`
  runs twice by design). The flag now means the same thing regardless of when it is set.
- `run_multi_local.py` sets the gate on the source **before** `plan()`, so the donor list is also
  derived from the gated corpus.
- **The check moved into the run, as an EQUALITY.** All three failures were "flag on, nothing
  happened", each caught only by inspecting artefacts after the compute was spent. The runner now
  **aborts** unless the build matches the frozen C-7 artefact exactly: the 5 named rejected columns,
  `n_samples` 42,600, `n_age_labeled` 42,581, masked `{("N2","no_control_baseline"): 19}`.

  The first version of this guard was directional (`n_age_labeled < n_samples`) and was **too
  weak** — externally caught. It would have accepted 42,605 cells with 1 masked label, and 42,605
  with 19 masked; both are wrong, and this defect's whole signature is a plausible-looking silent
  discrepancy. The comparison is now a pure function (`c7_mismatches`), unit-tested with no build,
  and validated end-to-end: it accepts all six frozen `_c7` folds and rejects all six invalid `_c7t`
  folds with the precise per-invariant diagnosis. The invariant was verified **identical in all six**
  frozen folds before being frozen as constants.

  Recorded for the next person: if a deliberate QC/`MAX_CELLS`/corpus change moves these numbers,
  re-freeze the constants from a verified build. Widening the check to make it pass reintroduces
  the defect.

### Not claimed

Nothing about C-7's effect on model performance. Run 1's metrics are pre-C-7 arm A and are **not** a
C-7 result; they were not snapshotted and no comparison was drawn from them. The frozen `_c7`
dataset-only folds are unaffected and remain correct (42,600 / 19 masked, re-verified here), so
every recorded analysis built on them stands.

---

## 2026-08-14 - DOUBLE-log1p BUG: every age in the dAge run was wrong; plus the clock compresses age ~3.4x

**Status:** Bug found, fixed, both runs redone, records corrected. **READ-ONLY**, `src/` untouched.
New: `experiments/clock_gse297234.py`, `tests/test_clock_gse297234.py` (13 tests),
`results/clock_gse297234_results.json`, `plans/STAGE_1_5_8_CLOCK_ON_GSE297234_PREREG.md`.

### The bug - mine, caught by my own unit test

`normalize_counts` applies **CP10k AND log1p** itself (`normalize.py:29`). Both scripts written in
this arc wrapped it in a **second** `np.log1p`. Every other script in the project calls it bare -
`build_dataset.py:184`, and `clock_fit.py:61` says so outright. Found while hand-computing an
expected pseudobulk value for a unit test; the run itself could never have caught it, because
`log1p` is **monotone** - it compressed every magnitude while preserving every ordering, so the
result still reproduced a published direction AND a published timing optimum on corrupted data.

### What changed in the dAge run (magnitudes ~2.4x larger)

| quantity | as published | corrected |
|---|---|---|
| pooled untreated age | 94.1 (bias +46.1) | **106.8 (bias +58.8)** |
| culture drift (controls - day 0) | +17.5 | **+32.6** |
| **dAge(transient), donor-clustered** | -17.89 [-26.52, -9.25] | **-42.45 [-67.39, -17.51]** |
| transient - failed, donor-clustered | -9.58 [-19.29, +0.12] | **-24.19 [-49.53, +1.14]** |
| M-E4 control SD | 5.04 | **8.43** |
| M-E5 control-arm batch shift | -8.52 | **-14.33 [-17.43, -11.23]** |

**Survives:** M-E0 CLEAN (operates pre-normalisation); M-E1 FAIL-CALIBRATION (worse); M-E3
REPRODUCED (donor CI still excludes 0); the D10/D13 optimum (D10 -46.19, D13 -43.02 - robust
because log1p is monotone); M-E5's arm-stratified control shift.

**Changed:** M-E4's branch - 5.04 was the "< half cv_mae" branch, 8.43 is the INCONCLUSIVE band, so
the bar whose inference I withdrew no longer fires at all.

**REVERSES:** the 2026-08-13 claim that "our -17.9 sits BELOW Gill's ~30 headline". At **-42.45** we
are **ABOVE** it, ~1.4x. The "suspiciously large" failure mode M-E3 watched for is now LIVE and must
be asked about, not dismissed.

### GSE297234: two donors 74 years apart (22 and 96)

Pre-registered before the run. GM23815 (22) -> **84.1** (bias +62.1); GM00731 (96) -> **106.1**
(bias +10.1). N1 ordering **CORRECT**. N2 **FAIL-CALIBRATION**. N3 **SOURCE IS NOT THE DRIVER** -
mean bias +36.1 vs GSE165177's +30, inside the +-10 band. N5 coverage 62.8%/92.7% vs 57.1%/89.2%,
so coverage does not confound N3.

**Two findings.** (1) **The Coriell hypothesis is dead, and it was mine.** P-N3 predicted a smaller
bias because these are Coriell lines like Fleischer's own training stock while GSE165177 used Lonza.
The bias is not smaller. Supplier drops out of the acquisition spec and ComBat-style harmonisation -
what Gill actually did - is the only route left. (2) **The clock COMPRESSES age rather than merely
shifting it: slope 0.297**, so a real 74-year gap renders as 22 years. An additive intercept
correction cannot fix that. Independent support: pseudobulk 0.297 and per-cell 0.307 now agree
closely, where before the fix they disagreed badly (0.160 vs 0.307).

**Open question, recorded not answered:** if sensitivity to CHRONOLOGICAL age is ~0.3 out of domain,
is sensitivity to reprogramming-induced change attenuated too? Between-donor discrimination and a
within-dataset same-day contrast are different quantities, so this run cannot say - but it is now a
live question about the magnitude of every dAge this project reports.

---

## 2026-08-13 - GILL 2022 VERIFIED: our +30 yr bias has a documented precedent AND a documented fix

**Status:** External claims checked against the paper (eLife 2022;11:e71624) before use, per the
standing rule. **READ-ONLY.** Recorded as section 8 of
`plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md`.

All five asserted claims verify. The two that matter:

**1. Gill ComBat-corrected the Fleischer reference BEFORE training.** Verbatim: *"trained a
transcription age-predictor using random forest regression on published fibroblast RNA-seq data
from donors aged 1-94 years old that was batch corrected to our transient reprogramming data set"*,
using `combat` from `sva`. Their clock's median absolute error is **12.57 yr against our cv_mae
12.27** - so the two clocks have essentially the same PRECISION. The difference between us is
**calibration, not resolution**: they harmonised the reference to their target before fitting; we
applied Fleischer's weights unchanged across a study boundary. **Our +30 yr floor is the documented
consequence of the one step they took and we did not.**

**2. Our data independently reproduces their TIMING result, which we had never checked.** Gill:
*"10 or 13 days may be the optimum for transcriptional rejuvenation."* Our mean dAge(transient) by
day: **D10 -22.65 (strongest), D13 -17.90, D17 -16.59, D15 -14.40**. Day 10 then day 13, recovered
with an UNCORRECTED clock. Per donor O1 and O3 peak at D10, O2 at D15. A second independent axis of
agreement, outside the pre-registration.

Magnitudes: Gill report ~30 yr (custom clock), ~20 yr (BiT), ~10 yr (Sarkar's data) vs negative
controls. Ours is -17.9 - inside that spread and BELOW their headline, so the "suspiciously large"
failure mode M-E3 was watching for does not fire.

**New candidate work item, not yet costed or pre-registered:** ComBat-harmonise GSE113957 to the
target dataset and refit, as Gill did. NOT claimed to fix our absolute ages - it changes the clock
from a frozen external artefact into a per-target fit, which every recorded dAge would have to be
re-derived against. Recorded as an option, not a decision.

---

## 2026-08-13 - CORRECTION: M-E3 was pseudoreplicated; the paired contrast does not survive

**Status:** External critique, verified before acceptance, then fixed. **READ-ONLY.** The §6 result
text is left standing; the correction is §7 of
`plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md`.

M-E3's CI was computed over 12 (donor, day) cells. Those come from **3 donors**, 4 timepoints each,
so they are not independent and the unit a generalisation claim rests on is the DONOR.

| quantity | 12 cells | 3 DONORS | excludes 0 at donor level? |
|---|---|---|---|
| dAge(transient) | -17.88 [-21.53, -14.24] | **-17.89 [-26.52, -9.25]** | **YES** |
| transient - failed | -9.58 [-13.17, -6.00] | **-9.58 [-19.29, +0.12]** | **NO** |

Per-donor: O1 -18.32/-10.62, O2 -21.13/-12.86, O3 -14.21/-5.26.

**M-E3's headline verdict STANDS** - dAge(transient) still excludes zero at n=3 with t(df=2)=4.303,
so Gill 2022's direction is still reproduced. **The paired transient-failed contrast is
DOWNGRADED** to "direction consistent in 3/3 donors, interval includes 0 at n=3" and may not be
quoted as an established effect.

**A second bug found while checking the first, and it is mine:** `ci()` fell back to the NORMAL
quantile 1.96 for any df > 10. M-E3 runs at n=12 (df=11, t=2.201), so the recorded 12-cell CI was
~12% too narrow. Fixed to fall through to `scipy.stats.t`. No verdict changes.

The general lesson, twice now in this project: build the interval on the unit the claim generalises
over. 3a-bis pooled two causes; M-E3 pooled within-donor cells. Both were caught externally, not by
the bar.

---

## 2026-08-12 - dAge ON GSE165177: Gill 2022's rejuvenation REPRODUCED; absolute age FAILS; two of my own bars were badly posed

**Status:** Run and recorded, graded against `plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md`
(committed BEFORE the script existed). **READ-ONLY - no build, no retrain, `src/` untouched.** New:
`experiments/dage_gse165177.py`, `tests/test_dage_gse165177.py` (18 tests),
`results/dage_gse165177_results.json`. 93 samples x 35,720 genes.

Normalisation was the main hazard and is handled by reusing the pipeline's own path for the
identically-formatted gill_bulk (`sources.py:506`): `2**log2 - 1` -> `normalize_counts(1e4)` ->
`log1p`. The clock declares log1p_cp10k and the file ships Log2 RPM; feeding it straight in would
have inflated every age by ~1.44x.

### M-E0 CLEAN. 0 of 93 columns rejected by the C-7 integrity gate

The same gate rejected five gill_bulk columns, one of which (N2_Fib_Sendai_Exp2) corrupted five of
six folds. **By the project's own standard GSE165177 is cleaner than the dataset every recorded
dAge was computed on.**

### M-E1 FAIL-CALIBRATION. Absolute age is unusable here

O1 (53) -> 89.2, O2 (53) -> 97.3, O3 (38) -> 95.9; pooled **94.1 vs a true mean of 48.0, |d| = 46.1
yr, nearly 4x one cv_mae**. Gene coverage does NOT excuse it: 57.1% of clock genes are present but
that is **89.2% of the clock's total |weight|**. The 53-vs-38 contrast was pre-registered as not
gated and is reported as indicative only - and it points the wrong way (O3, the youngest, reads
second-highest).

### M-E2 The contemporaneous-control dAge, computed for the first time in this project

`clock(sample) - mean(clock(controls of the SAME donor at the SAME day))`. failed -6.2 to -11.0 yr;
transient -3.1 to -24.3 yr. gill_bulk is structurally incapable of this - its only baseline is one
day-0 sample per donor, cross-batch for ~half the data.

### M-E3 REPRODUCED. Transient reprogramming rejuvenates, by our own clock

dAge(transient) vs its own contemporaneous control: **-17.88 yr, 95% CI [-21.13, -14.64]**, n=12
cells. Paired transient - failed: **-9.58 yr, CI [-12.77, -6.39]**. 11 of 12 (donor, day) cells
negative. Gill 2022's central claim recovered on data that was in no training config.

**Why M-E3 succeeds while M-E1 fails, and it is the load-bearing point:** dAge is a DIFFERENCE, so
the clock's +46 yr bias and every missing-gene term appear in the treated sample and its control
alike and cancel exactly. Absolute age needs the clock to be accurate; dAge needs it only to be
consistent. First direct evidence that the control-relative design does the job it was chosen for.

### M-E4 The bar fired, but the inference I attached to it does not follow - WITHDRAWN

Pooled within-(donor, day) control SD = **5.04 yr** against cv_mae 12.27, i.e. the "< half"
branch, whose pre-registered reading was "the per-donor offset is MORE LIKELY REAL BIOLOGY; first
evidence for Stage 2's premise". **That reading is wrong and this same run disproves it.** cv_mae
is cross-validated error against TRUE age across 133 donors; 5.04 is replicate scatter within one
condition. A clock can be reproducible and badly biased - M-E1 shows exactly that, 5 yr scatter
alongside a +46 yr bias. **Reproducibility is not accuracy.** The bar fired as written; the
conclusion is withdrawn. Stage 2's premise is NOT strengthened. What stands is narrower and still
useful: the zero-point is reproducible to ~5 yr once controls are contemporaneous and replicated.

### M-E5 The pooled bar said "sub-error"; the pooled statistic was the wrong one

Pooled exp1-exp2 = -2.91 yr, CI [-6.35, +0.54]. Stratified by arm (not pre-registered): **control
-8.52 yr, CI [-10.13, -6.92]**; failed -4.29 [-7.52, -1.06]; transient +7.59 [-3.60, +18.77].
Control and transient move in OPPOSITE directions, so the pooled mean cancels an effect real in
both. The control CI excludes zero: **the zero-point itself shifts -8.52 yr between batches, 0.69
cv_mae**, and because dAge is measured against that control an arm-dependent offset does NOT
cancel. **D1 is confirmed material for any cross-batch comparison** - my pooled bar could not have
seen it.

### Net

GSE165177 is **valuable for dAge and now demonstrated to be** - exactly the opposite of its verdict
for p_unsafe, and precisely the split regime E predicted. The control-relative design is
**vindicated**; the clock's absolute age on bulk fibroblasts is **unusable**; D1 is **confirmed**;
D2 is **quantified** at 5.04 yr. **Not claimed:** that the clock is calibrated here, that 53-vs-38
resolves, or that these values compare to the project's harmonised figures - no harmonizer was
applied, by design.

---

## 2026-08-12 - REGIME E: P0 fired, and p_unsafe turns out not to be expressible in bulk at all

**Status:** Run and recorded, graded against `plans/STAGE_3A_REGIME_E_PREREG.md` (committed
81602fe BEFORE the script existed). **READ-ONLY - no build, no retrain, `src/` untouched, no 3a
verdict taken.** New: `experiments/stage3a_regime_e.py`, `tests/test_stage3a_regime_e.py` (26
tests), `results/stage3a_regime_e_results.json`.

### The pre-registered precondition fired

93 of 95 GSE165177 samples loaded (2 excluded iPSC lines), 35,720 genes, donors O1/O2/O3 aged
53/53/38, 33 contemporaneous controls. Unsafe fraction by (donor, day): **1.000 everywhere except
O1 at day 10 (0.750)**. SD across timepoints: O1 0.100, O2 **0.000**, O3 **0.000**. Two of three
donors flat -> **P0 fires -> E1-E4 must not be read**, and they were not; the null was not run.

### Why - the part that matters more than the verdict

The **untreated day-0 fibroblasts label `loss`**, P(loss) = 0.966 / 0.876 / 0.730. They are the
starting material, so that cannot be biology. `fate_labels` z-scores each program against the
`is_control` samples, which here are fibroblasts cultured 10-17 days; anything differing from that
reference lands on the unsafe side. The split produced is **control vs non-control, not a time
course** - negative controls P(safe) 0.699, every other arm P(loss) 0.61-1.00.

### The structural finding

**`p_unsafe` is a fraction OF CELLS. A bulk sample is ALREADY a population average, so a hard
label per sample collapses the fraction to 0/1 before it can be counted, and the "fraction"
becomes a fraction of SAMPLES.**

One mechanism explains two things previously blamed on sample size: gill_bulk's 63/70 values
pinned at the bounds, and why more bulk replication cannot help. GSE165177 has 4-6x the
replication and real contemporaneous controls and is **more** saturated, not less - 11 of 12 cells
at exactly 1.000. Regime E could never have succeeded, and the A1/A2/A3 attribution 2x2 is moot:
fold count and replication were never the binding constraint for the safety target.

### What it settles

The acquisition ask is **RE-ESTABLISHED on much better grounds**: not "more replication" but
"`p_unsafe` requires single-cell resolution, and GSE242423 is the only single-cell dataset we
have". Requirement H2 of the spec (single-cell, not bulk) is **promoted from one of eight
requirements to the whole point**.

AUDIT-2's challenge was **right to make** - testing it cost one script run and replaced an
assertion with a mechanism. Its expectation that GSE165177 would carry the gate is **refuted**, for
a reason neither of us had identified. **GSE165177 is NOT dismissed**: its replication, 33
contemporaneous controls and three adult in-range donors are real advantages for **dAge**, which is
continuous per sample and does not have this problem. Separate and still open.

### A bug this run caught in itself

The first implementation matched only `{donor}_{arm}_{N}days_{exp}`. GSE165177 names its day-0
fibroblasts **`O1 Fib`** - space-separated, no `days` token - so all three were **dropped
silently**, cutting every trajectory's first timepoint and the pair count from 10 per donor to 6.
Nothing raised; the run just answered a smaller question, and would have agreed with AUDIT-2's
"4 timepoints, 6 pairs" for the wrong reason. Fixed and pinned by a named regression test.
Including day 0 did **not** rescue P0 - those samples label `loss` too, which is what exposed the
mechanism.

Also recorded: the gene-space join listed as the next thing to cost before acquisition is **moot
for the safety target** - joining gene spaces cannot make a bulk sample express a per-cell
fraction. It stays open and worthwhile for dAge.

---

## 2026-08-12 - AUDIT-2: the 3a diagnosis is right; "the blocker is more data" is not

**Status:** Audit. **`src/` untouched, nothing withdrawn.**

### What is right, including one place I was wrong

The STOP withdrawal and the resolvability study are correct and I agree with both.
`partC_frac_pred_outside_unit = 0.773` - **77% of Part C's predictions leave [0,1]** on a target that
is a fraction by construction. And **regime A, the geometry 3a actually ran on, is UNRESOLVABLE at
every alpha, raw and logit** (pass rates 0.000-0.017 against min_pass 0.95). A correct system could
never have returned GO. That is a strong finding, properly done - 2000 trials, 4 regimes, both links.

**I was wrong about the mechanism.** My audit blamed dt extrapolation on Y1. Measured:
`heldout_pairs_outside_train_dt_support = 0` and `z_dt_heldout_absmax = z_dt_train_absmax = 2.4427`,
identical. The extrapolation is in the **gene block** - `z_gene_heldout_absmax = 39.45` against
`z_gene_train_absmax = 7.24`. Right conclusion, wrong reason. Recorded, not dropped.

### Where I disagree - the blocking dataset is on disk

The diagnosis names per-timepoint replication and a line/modality shift as binding. Both are true of
`gill_bulk`; **neither is true of `GSE165177`**, which is on disk, in **no training config**, and has
been shortlist item #2 - *"95 adult in-range methylation-paired samples, in no training config,
free"* - since before this arc began. Measured from the raw matrices:

| | gill_bulk (what 3a ran on) | **GSE165177** |
|---|---|---|
| donors | 6 | **3** |
| timepoints | 12 | **4** |
| **samples per (donor, timepoint)** | **~1.7** | **7.0-8.0** |
| **controls** | **1 per donor, day 0 only (6 total)** | **2-3 per donor PER TIMEPOINT (33 total)** |
| condition arms | time course | **6**, incl. explicit `negative_control` |

**Every failure this arc has fought traces to gill_bulk's one unreplicated day-0 control per
donor** - G-a, C-7/N2_Fib, Group D's self-centring, Group E, §4.7's 16.67 yr instability, rule 4 and
B2'. **GSE165177 has replicated contemporaneous controls at every timepoint; none of those failure
modes can arise in it.**

And it fixes the named mechanism directly: `p_unsafe` from 1.7 samples can only be {0, 0.5, 1},
which is why `saturated_at_bounds = 63` of 70. From 7-8 samples it moves in eighths.

### The honest cost - a trade, not a free upgrade

3 LOO folds not 6; **~6 ordered forward pairs per donor not 66** (4 timepoints vs 12); ~18 total
pairs against 319. That is a large loss in pairs against 4-5x the replication and real controls -
319 noisy correlated pairs versus ~18 well-estimated ones. **Nobody has evaluated that trade**, and
"more data is the blocker" skips it rather than answering it.

### What to do - free, machinery already exists

Add **GSE165177 as regime E** to `experiments/stage3a_bis_resolvability.py` at its real geometry.
RESOLVABLE means the gate can run for real on data we already hold, with no acquisition and no
retrain. UNRESOLVABLE means the acquisition claim is **established with a number instead of
asserted**. Worth costing in the same pass: with 6 condition arms, forward pairs can be built
**across arms at matched timepoints**, a design that does not depend on the timepoint count at all.

**Do not go to Stage 5, and do not open Stage 6 acquisition, until regime E has been run.**

---


## 2026-08-12 - CORRECTION to 3a-bis: the binding constraint is the LINE/MODALITY shift, not cells per timepoint

**Status:** Run and recorded. READ-ONLY. Found by checking my own claim before registering a bar on
it. The 3a-bis entry below is **left standing**; the correction is recorded here and as a box inside
the "3a-bis" section of `plans/STAGE_3_TOOL.md`. New: two registered bars in
`tests/test_bars_resolvable.py`. 976 tests pass.

### What was wrong

The 3a-bis write-up attributed regime A's failure to the **held-out cells per timepoint**. Regimes A
and B differ in **two** ways at once, not one: held-out precision AND held-out distribution - an HFF
pseudo-replicate is the same line and the same single-cell modality as the training data, while a
Gill donor is a **different line measured in bulk**. A-vs-B cannot separate them, so the attribution
was a confound.

It surfaced from an idealised check while drafting the registered bar: give both arms their optimal
predictors and binomial noise at 2 cells/timepoint does NOT stop the comparison resolving. Target
imprecision alone therefore could not be the whole story.

### The 2x2 that settles it

C = the within-HFF trajectory scored at 2 cells/tp; D = the Gill donor scored at 472, a
counterfactual Gill can never have, included only to separate the causes. Pass rate at alpha=1,
logit:

| held-out trajectory | at 2 cells/tp | at 472 cells/tp |
|---|---|---|
| **HFF pseudo-replicate** | **C: 0.965** RESOLVABLE | **B: 1.000** RESOLVABLE |
| **Gill donor** | **A: 0.000** | **D: 0.000** |

**The binding constraint is the LINE/MODALITY SHIFT.** Handing the Gill geometry 472 cells per
timepoint changes nothing - 0.000 at every alpha. Handing the HFF geometry only 2 cells still
resolves at full amplitude.

**What precision does control is the SENSITIVITY FLOOR**, and that is real: at alpha = 0.5 the
within-HFF geometry drops from **0.990 to 0.429** when scored on 2 cells instead of 470.

**Limit on the correction:** every Gill donor is bulk and HFF is single-cell, so line and modality
are **perfectly confounded** in this corpus. The 2x2 separates precision from the shift, and no
further.

### Registered bars (REF_GROUND_RULES 5b)

Two entries added to `REGISTERED_BARS`, both on the PRECISION axis, which is the only one that is
synthesisable self-contained: held-out at ~470 cells/tp on HFF's full measured curve -> RESOLVABLE;
held-out at 2 cells/tp at half amplitude -> UNRESOLVABLE. The second carries an explicit note that
it does NOT explain the Gill geometry's failure, so the correction cannot be silently re-lost.

### Downstream

Stage 6's requirement is unchanged in substance - a second dense single-cell time course on a
different line - but for a different reason: a held-out line **of the same modality** is the only
way to test transfer at all, not more cells per timepoint. Still no 3a verdict taken; 3b/3c/3d stay
unwritten.

---

## 2026-08-12 - STAGE 3a-bis: 3a was graded on 0.3% of the cells, and as designed could never have returned GO

**Status:** Run and recorded. **READ-ONLY - no retrain, no rebuild, `src/` untouched. NO 3a verdict
is taken.** New files `experiments/stage3a_bis_resolvability.py` and
`tests/test_stage3a_bis_resolvability.py` (18 tests, pass); new artefact
`results/stage3a_bis_resolvability_results.json`. Recorded as the "3a-bis" box in
`plans/STAGE_3_TOOL.md`. This is the REF_GROUND_RULES 5b precondition that must precede any bar.

### What 3a was reading

`test18_forward_gate.py:74` builds every row from `gather_split(..., REGIME, "test")`, and under
the `holdout` regime the test split is **the held-out Gill donor and nothing else**. So 3a ran on
18-21 bulk samples per fold at **1.7 cells per timepoint** - while the SAME bundle's train split
holds **33,613 HFF single cells over 9 timepoints at ~3,735 cells per timepoint**, whose unsafe
fraction runs **0.0835 -> 0.9996** with a per-timepoint SE of **0.006** (SD across timepoints over
typical SE = **42:1**). The gate was decided on **~0.3%** of the available cells, and the 99.7% it
ignored are the precise ones.

### The check

The simulated truth is **HFF's own measured curve** rather than an invented effect size:
`p(t_j) = g_bar + alpha*(g(t_j) - g_bar)`, observed at the **real cell counts**
`u_j ~ Binomial(n_j, p)/n_j`, so 1-2 cells and ~470 cells enter as the noise they actually are.
Under this truth `state+dt` can reach t_j and `state` alone cannot, so a working test MUST detect
it. 3a's rule graded verbatim, 2000 trials/cell, MIN_PASS_RATE 0.95.

**Both regimes share almost the same training set** - 643 vs 679 pairs, both containing every HFF
pseudo-replicate - so the only thing that differs is how precisely the HELD-OUT trajectory is
measured. A controlled comparison, not two experiments.

| alpha | B: held out on ~470 cells/tp (raw/logit) | A: held out on 1-2 cells/tp (raw/logit) |
|---|---|---|
| 0.00 | 0.038 / 0.017 | 0.006 / 0.011 |
| 0.25 | 0.468 / 0.337 | 0.005 / 0.017 |
| 0.50 | **0.976 / 0.990** | 0.001 / 0.005 |
| 1.00 | **1.000 / 1.000** | **0.000 / 0.000** |

### The finding

**Regime A - holding out a Gill donor - is UNRESOLVABLE at every amplitude**, including alpha = 1
where the truth is HFF's full measured curve (pass rate **0.000**), and it FALLS as the effect grows
(0.011 -> 0.017 -> 0.005 -> 0.000). **So Stage 3a as designed could not have returned GO for any
signal whatever** - its STOP was not merely unsupported, the test was structurally incapable of any
other answer.

**Regime B is RESOLVABLE from alpha = 0.5 upward** (0.990, then 1.000): the corpus can recover its
own forward curve when the held-out target is measured on ~470 cells rather than 1.7. **The binding
constraint is the HELD-OUT cells per timepoint, not the training-set size** - regime A has MORE
training pairs than B and still reads 0.000.

At alpha = 0 every cell is <= 0.038, so there is no false-positive problem: the bar is specific and
blind. Even in the good regime **alpha = 0.25 fails (0.337)** - the detectable effect must be at
least about half HFF's measured amplitude.

### Pre-registered: 4 of 4 held

Q1 regime B resolvable at alpha 1 **YES** (1.000). Q2 regime A raw not resolvable at alpha 1
**YES** (0.000). Q3 held-out cells bind rather than training size **YES**. Q4 no false positive at
alpha 0 **YES** (max 0.038).

### What this changes

3a's STOP **stays withdrawn** and is now known to have been **unfalsifiable**. Section 3a's closing
line - *"this dataset cannot support forward prediction"* - is **wrong as stated**: the dataset
supports it where it is measured densely, and it is the Gill BULK arm that cannot serve as a
held-out target at any effect size. Test 18's cells-per-timepoint warning is **promoted from a
warning to the binding constraint, quantified**. Stage 6 is likewise quantified: what is needed is
a **second dense single-cell time course on a DIFFERENT line**, not more donors measured at 1-2
cells.

**Not claimed:** that the tool is buildable. Regime B holds out a **pseudo-replicate of the same
culture**, so it shows the forward curve is recoverable, never that it transfers to a new line. The
product question - cross-line transfer - **remains unanswerable with this corpus**, now as a
measured statement rather than an inference from a broken run. 3b/3c/3d stay unwritten.

---

## 2026-08-12 - STAGE 3a DIAGNOSED: a forward time signal IS present; the estimator could not find it and the bar could not register it

**Status:** Run and recorded. **READ-ONLY - no retrain, no rebuild, `src/` untouched.** New files
`experiments/stage3a_diagnose.py` and `tests/test_stage3a_diagnose.py` (20 tests, pass); new
artefact `results/stage3a_diagnose_results.json`. Recorded as the "3a DIAGNOSIS" box in
`plans/STAGE_3_TOOL.md` and as "STAGE 3a WITHDRAWN AND DIAGNOSED" in the lab notebook. All four
items accepted when 3a's STOP was withdrawn (81b70ad) are run here, plus one section (D0) added
after the shape of the target became visible and labelled as not pre-registered.

`test18_forward_gate.py` is imported UNMODIFIED and its own primitives are used, so what is
diagnosed is the estimator 3a actually ran. All ten of 3a's Part A/C numbers reproduce to the digit.
The fast ridge path used by D4 is validated against `sklearn.linear_model.Ridge` first
(max|d| = 3.9e-10) and pinned as a property in the test file.

### D0 - the finding: a model-free predictor using only t_j beats every arm 3a ran

Part C's target is the mean of **1.6-1.8 binary cells per timepoint**, and **63/70 = 90% of its
values are exactly 0 or 1** - per donor, close to a monotone step in time. So the forward question
has a model-free ceiling: predict the held-out donor's value at t_j from the OTHER donors at that
same t_j. No genes, no fitting, nothing from the held-out donor (pinned by a leakage test).

On the five folds 3a graded: pooled-mean baseline 0.493, **3a's state+dt raw 2.010** (the arm the
STOP was read from), 3a's state-only 0.409, bounded logit 0.352, **oracle on t_j alone 0.157**.
**Paired -0.336, 95% CI [-0.450, -0.221], n=5 -> t_j HELPS**, every fold the same direction; the
oracle beats 3a's best arm by 2.24x. Same for dAge: oracle **14.208** vs pooled 23.615, paired
**-9.407, CI [-16.332, -2.482]**.

**The limit, stated as plainly as the claim.** Those five folds exclude **N2**, the one donor that
never becomes unsafe (flat 0 across 11 timepoints). Include it and the advantage falls to **-0.183,
CI [-0.411, +0.046] -> tied**, N2 the only negative fold. N2 is absent from 3a's grading for an
unrelated reason: `timepoint_table` filters on the AGE mask before computing a target that does not
depend on dAge at all, so C-7's rule 4 removed it from the SAFETY analysis too. **Established on the
graded geometry; not established across all six donors.**

### D1 - the divergence is real; BOTH proposed mechanisms are wrong

The audit's proposal - Y1's dt lies outside training support - **does not hold**: all five folds
share the identical dt range [0.14, 11.77] with **0 held-out pairs outside**, Y1's missing timepoint
being interior. My own pre-registered guess (P2, the dt block) is **also refuted**: Y1's dt
coefficient (2.089) and mean |dt contribution| (1.634) sit inside the other folds' range.

It is the **gene** block. Adding the two dt columns re-solves the ridge and rotates the gene-weight
vector ~60 degrees (cos 0.47-0.52), amplifying the gene block **7.7x-20.4x on every fold**. Y1's
held-out state is the one far outside the training scaler's support - mean |z| **2.546 vs 0.671**.
Rotated weights times out-of-support state gives predictions of **5.44 to 10.92** on a target
bounded by 1. Not only Y1: raw Part C is out of range on **every fold** (33%, 33%, 62%, 77%, 100%).

### D2 - bounding repairs the instrument, not the answer

Part C, 5 folds: raw +1.601 CI [-2.270, +5.472] width **7.742**; clip -0.028 CI [-0.105, +0.050]
width **0.155**; logit -0.027 CI [-0.109, +0.054] width **0.163**. Y1 goes **7.589 -> 0.314**, the
CI narrows **50x**, the point estimate flips sign - and all three still read "tied". Read against
D0 that is the informative part: the signal is there and a bounded ridge over 2000 genes plus two
dt columns still does not reach it. The failure is the **model class and the geometry**, not the
corpus. Pre-registered P4 **failed**: state-only was not broken (0.409 vs a 0.446 mean-only
baseline; 0.380 once bounded).

### D3 - Part B closed arithmetically

SD of the swing across folds **2.5e-14**; analytic value **-269.12592**; max |observed - analytic|
**5.7e-14**. The expression has **no x0 term**, so the five rows cannot disagree. Re-run correctly
(one LODO fit per fold over its own dt range) the swings finally differ - -241.91, -247.74, -264.16,
-267.33, -274.71, SD 13.80 yr - and remain **nonphysical**. Part B never "passed" the >2 yr clause.

### D4 - `bar_verdict` at the real geometry (REF_GROUND_RULES 5b), never run for 3a

Real features, real dt, real fold/pair structure; only the TARGET is simulated so a dt effect is
present by construction. 2000 trials/cell, MIN_PASS_RATE 0.95, 3a's rule graded verbatim.

**The raw estimator is not a detector**: its pass rate is non-monotone in the true effect - 0.000 at
rho 0.5, **1.000** at rho 0.75, back to **0.000 at rho 1.0** where the target IS dt with sigma
0.001. Even repaired, the logit version clears 0.95 only at rho = 1.0 with sigma = 0.001 (0.055 at
rho 0.75; 0.447 at rho 1.0/sigma 0.05, 0.823 at six folds). At rho = 0 both are 0.000 -
**specificity without sensitivity**.

### Pre-registered, graded: 3 of 5 held

P1 dt nested inside training **YES** (0 outside). P2 dt block carries the divergence **NO** (the
gene block does). P3 Part B identical and analytic **YES**. P4 state-only also broken **NO**. P5 bar
unresolvable even at pure dt **YES** for raw. Both failures recorded as failures.

### What this changes

3a's STOP **stays withdrawn**, and the ground has shifted - on the graded geometry it points the
wrong way. The withdrawal banner's *"the forward signal may well be absent"* is **superseded** and
annotated in place, not deleted. Part B **must not be read at all**. Part C's prediction **must be
bounded**. **No 3a verdict may be taken again until a resolvable bar is registered under 5b** with
`tests/test_bars_resolvable.py` updated. 3b/3c/3d **stay unwritten**; Stage 5 is **still not entered
on this basis**.

**Not claimed:** that the tool is buildable, or that the oracle is a model - it needs t_j and other
donors' measurements at t_j, so it is a **ceiling and an existence proof**. A valid re-run needs
four things, listed in the plan: a bounded predictor, a resolvable bar, N2's status resolved, and a
model class that can reach the ceiling D0 measures.

### Correction to the withdrawal banner (81b70ad)

The two CIs in item 2 of that banner were recomputed by hand with the `T_CRIT` entry for one fold
too few. Means and the with-Y1 width are right; the endpoints are not. Measured: with Y1
**[-2.270, +5.472]** (width 7.742, matching the recorded run); without Y1 **[-0.194, +0.615]**
(width **0.810**, not 1.095). Corrected as a dated note; the original text is left standing. No
conclusion drawn from them changes.

---

## 2026-08-11 - AUDIT: C-7's fix is REAL and large; Stage 3a's STOP is NOT safe to act on

**Status:** Audit. **`src/` untouched, nothing withdrawn.** New file
`plans/STAGE_3A_VERDICT_AUDIT.md`.

### What genuinely got fixed - this is not a reclassification

C-7 implemented and adopted on the dataset. HFF day-14 dAge per fold, before -> after:
N2 -7.352 -> -7.337, N3 -22.121 -> **-5.628**, O1 -24.023 -> **-6.514**, O2 -22.891 -> **-6.674**,
Y1 -22.049 -> **-4.606**, Y2 -23.869 -> **-8.292**.

**Fold spread 16.671 -> 3.686 yr, a 78% reduction**, magnitudes landing on §5.14's pre-registered
-8.196. n_cells 42605 -> 42600 (exactly the five rejected columns); N2's 19 cells masked in every
fold. A real bug was found on the way - the gate flag never reached injected sources, so the gate
was inert when first wired (e6fc183).

### But Stage 3a's STOP rests on a diverging fit

`p_unsafe` is built at `test18_forward_gate.py:85` as a **fraction** - `((cls==LOSS)|(cls==DEATH))
.astype(float)` then `.mean()`. **It cannot leave [0,1].** Yet Part C's Y1 fold reports state+dt
**MAE 7.589**. A model predicting ~8 for a quantity bounded by 1 is producing out-of-range output.

Corroborating: Part A's Y1 goes **22.84 -> 311.47** by adding one feature under a standardized ridge
at alpha=1.0; Part B's swing is **-269.13 for all five folds, identical to 2 dp**, which five
independent LOO fits cannot do by chance; Part B's magnitude is **269 yr**, nonphysical against a
>2 yr threshold; and the paired CI (mean +1.601, 95% [-2.270, +5.472], n=5) is **driven entirely by
Y1** - without it the gains are +0.049, -0.184, -0.560, -0.147.

**Y1 has 11 timepoints not 12** and 55 pairs not 66, so held out its dt distribution sits outside
the training support and the linear model extrapolates. **A data-shape artefact, not a property of
forward prediction.**

**And no resolvability check.** §5b requires simulating a correct system before the run; there is no
`bar_verdict` for 3a. With n=5 and a CI spanning 7.7 units on a [0,1] target the honest verdict is
**UNRESOLVABLE**, not "tied". The script's own caveat - *"a NEGATIVE result is decisive"* - holds
only when the negative comes from the data.

### What is NOT established

**RES improved: UNKNOWN, not measured.** RES is a model output and **no retrain has happened** - the
newest scorecard, `gc2_D_stratshuffle_hff_s0.json`, predates C-7. Every Stage 1 guard under the new
labels is likewise not re-reported, though C-7 §5 requires it on adoption.

**"Stage 1.5.6 closed" is fair. "Ready for Stage 5" is not.**

### Recommended, none of it needing a retrain

Diagnose Y1's dt support; bound Part C's prediction to [0,1] (clip at minimum, logit link honestly);
explain the identical -269.13 swing before Part B is read; and run `bar_verdict` for 3a at this
geometry. **If a correct system cannot clear the bar with 5 folds and 55-66 pairs, then the dataset
cannot answer the question - and that, not STOP, is the finding.**

---


## 2026-08-08 - C-7's five "open decisions" are not decisions: four dissolve, one has a third answer

**Status:** Recorded. **Still NOT implemented. `src/` untouched, no label moved, §11 unmodified.**

The spec closed on five items marked "yours to decide", with recommendations but no data. A project
whose discipline is *measure, don't decide* should not close a change that way - and it does not
have to. Every answer below is read off the code with the line cited.

| | asked | answer | basis |
|---|---|---|---|
| **A1** | gate in `fetch` or the build loop | **fetch** | `src.fetch` is called at **three** sites - `build_dataset.py:170, 289, 345`. Build-loop gating means three edit sites and **missing :345 leaves the degenerate column in `sigma_ref`**, which is the whole defect. Source gating is one site, covers three |
| **A2** | assert against single-cell sources | **question DISSOLVES** | under A1 a single-cell source never calls the gate. And G1 is defined on **RPM** while `GSE242423SingleCellSource` yields **raw UMI counts** (~1e3-1e4/cell; `normalize_counts` runs later), so G1 would reject every cell by construction of the units |
| **B1** | which census pass | **dedicated - and the flagged cost is not real** | rule 4 fires only where a line can lose a control; only the gate removes controls; the gate is bulk-only. **So the census is bulk-only** - `gill_bulk`, one 8 MB file, 124 columns. There is no "second full corpus read". (`fit_harmonizer` is also conditional, `:383` - reusing it would leave non-harmonized builds ungated) |
| **C1** | rule 4 first or second | **before `donor_out_of_clock_range`** | `age_mask_reason` is a **persisted parquet column** (`io.py:139, 265`) and `schemas.py:57-59` says the order is meaningful. Free today because C-2 is off, but **N2 is donor age 0, so once C-2 activates it matches BOTH rules** and the order decides the recorded reason. "No zero-point exists" is *undefined*; "outside the clock's range" is *out-of-validity* - undefined is stronger. Decide it now while it costs nothing |
| **D1** | thread the mask or assert at the call site | **NEITHER - the census already records it** | `aging.py:121` already writes `"source": "controls" if ctrl.any() else "self_fallback"` per line. **That is Stage 1.5.2's G-a gate; B2' is the same census one field further** - one assertion, at the single place the fallback occurs, not two places it is consumed |

**D1's real work item**, smaller than either option offered: of the two `_control_baseline` call
sites only `delta_age` passes a census (`aging.py:301`); `aging.py:245` does not. That site must
pass one or its fallback stays invisible.

Recorded because the pattern matters more than the items: **"I left these to you" on a change this
mechanical is a signal to go and look, not a signal to choose.**

---


## 2026-08-08 - C-7 decides: option (c), and B2 does not have to block it

**Status:** Decision recorded. **Still NOT implemented. `src/` untouched, no label moved.**

### Option (c) adopted; (a) superseded

**Reject the degenerate control, MASK N2's dAge, keep the donor and the fold.** §3's (a) conflated
three separable decisions and answered all three with the harshest available answer:

* should the degenerate control enter the harmonizer? **no**, definitively (§5.14)
* should N2's 21 dAge labels survive? **no** - zero-point **98.65 yr** for a donor of true age **0**
* should N2's **cells** survive? **yes** - the fate head consumes no dAge and runs at `fate_roc`
  **0.983**, *"untouched by every dAge problem"*. Dropping the donor destroys working fate data to
  fix a broken age label

Verified independently - the clock on each day-0 fibroblast: N2 **98.65** (true 0), O2 79.50, O1
79.12, Y1 64.92, Y2 57.66, N3 36.44. N2 sits **+35.12** above the other five's mean of 63.53.

*Context, not a new defect:* the clock over-predicts **every** donor by +22 to +36 yr, and **N3 is
also age 0 yet reads 36.44** - it cannot separate a neonate from a 35-year-old. Absolute age, where
the intercept does not cancel (§0 ERROR 1), so known behaviour. Recorded because it bears on C-2.

**(c) keeps LOOCV at SIX folds** - §5's "re-report every guard over 5 folds" does not apply, §4.7's
record stays comparable, and donors are not spent.

### The B2 collision dissolves

The objection is real: rejecting N2's control leaves line N2 with zero controls and
`_control_baseline` self-centres. **But B2's purpose was never to forbid reaching the fallback - it
was to forbid the fallback producing a label that is KEPT**, the `age_label_policy` fail-open. So:

> **B2' - no line may reach the fallback AND retain its dAge label.**
> `assert not (fell_back and not masked)`.

### Rule 4 should be GENERAL, not a donor special case

> **Rule 4 - a `cell_line` with zero admissible controls has no zero-point, so its dAge is
> undefined and is masked.**

Keyed on **data integrity**, not identity. No donor name. And it fires on **exactly** the condition
that triggers the fallback, which makes B2' automatic rather than colliding. It also closes Stage
1.5's **Group D** defect properly for every future dataset, not just N2.

### The subtlety that decides implementability

`_control_baseline` falls back when a line has no controls **in this chunk**, and dAge is computed
per chunk. Two distinguishable cases: **no controls at all** (rule 4 / B2' - mask) versus **controls
exist but none in this chunk** (Stage 1.5 **Group E**, separately owned). **B2' must test the first
and not the second**, or it fires on Group E's case and blocks C-7 for the wrong reason. The
predicate is global per `cell_line`, evaluated in the same pre-pass `fit_harmonizer` already runs.

### Sequencing

**Rule 4 ships WITH C-7, not after.** The gate alone creates the orphaned line - reject the control
without rule 4 and line N2 has no zero-point and no mask, the exact window B2 forbids. C-7 is
already a `src/` change, so adding rule 4 costs nothing extra and removes a state the system must
never be in.

**Still open:** drop vs **re-quantify** from SRA SRP302546 - now *optional* rather than urgent, since
under (c) the donor and fold survive either way and re-quantification becomes an upgrade that would
restore N2's dAge labels. And the ten degenerate non-control columns remain in scope, unresolved.

---


## 2026-08-08 - GEO checked: the defect IS in the deposit. C-7 is needed.

**Status:** C-7's one open question, resolved. **`src/` untouched, no label moved, no gate
implemented.**

C-7 §6 flagged one thing worth checking before implementation - whether the degeneracy is GEO's
deposit or our read. Checked.

### The file we hold is the file GEO serves

GEO lists exactly one supplementary file for GSE165176: `GSE165176_Log2_RPM_Sendai_reprogramming
.txt.gz`, **8.0 Mb**, and **124 samples** (GSM5027507-GSM5027630). Ours is 8 337 920 B = 7.95 MiB
with 124 columns. The `(1)` in our filename is the browser's, not a different file.

### Not damaged, and not our parse

`gzip -t` intact; 35 806 lines; **136 fields on every single line** (12 annotation + 124 samples),
so no column can shift. Values below were pulled with **awk**, not pandas.

### What the raw text contains

`N2_Fib_Sendai_Exp2` (column 33) takes **four distinct strings** across 35 805 genes: `11.489547`
x **35 690 (99.7%)**, `11.64155` x 106, `12.64155` x 8, `13.226513` x 1. The first eight rows read
11.489547 every time while sound `O1_Fib` varies normally. **`Y1_d7_CD13_Sendai_Exp1` has TWO
distinct values, the modal one covering 100.0%.**

**Not a transcriptome, not a parsing artefact. It is what was deposited. C-7 is needed, not
redundant.**

### A route that does not cost a donor

The GEO record points at **SRA SRP302546** for raw reads, so **`N2_Fib` could be re-quantified
rather than dropped.** That matters: C-7 §3's recommendation costs a whole LOOCV fold, reaches C-2
(N2 is donor age 0) and reaches §4.7 (whose 16.67 yr spread is defined over six folds including
N2's). One FASTQ restores N2's zero-point at no cost in donors, folds or guard re-reports - and it
is far cheaper than Stage 6 acquisition. **Recorded as an option, not adopted:** re-quantifying one
sample through a different pipeline than the other 123 introduces a batch term to be checked, not
assumed. The gate is unaffected either way.

### One gate candidate tried and REJECTED

**Distinct-value count does not reproduce the five.** The five flags carry 2, 4, 5, 5 and **27**
distinct values; the other 119 carry **22**-693. `N2_d21_CD13` (flagged) has 27 while sound
`O2_d40` has 22 and `O2_d34` has 26 - **the populations overlap**, so distinct-count would flag two
sound columns before reaching the fifth degenerate one. G1 and G2 separate cleanly; this does not.
Recorded so it is not retried.

---


## 2026-08-08 - C-7 pre-registered: bulk sample integrity gate (NOT implemented)

**Status:** Pre-registered only. **`src/` untouched, no label moved, no existing plan edited.**
New file `plans/CHANGE_C7_BULK_SAMPLE_INTEGRITY.md`.

Takes the next free **change ID** rather than a section number, because two machines pushed a
`## 5.8` concurrently and neither saw the other's. C-1..C-6 are Stage 1.5.3's; this is C-7.

### The gate - two conditions, both justified by UNITS

An earlier proposal thresholded `mean - min` at 1/5 of the cohort median. **Rejected** (§5.10): it
cuts a continuous distribution 8% from its neighbour, so on a new cohort it flags or misses
arbitrarily. Replaced by:

* **G1 library** - the matrix is Reads Per Million, so a sound column's linear values must sum to
  ~1e6 **by definition**. Accept `[1e5, 1e7]`, a decade either side. The 5 degenerate columns sit at
  1.694e+07 - 2.148e+09; the other 119 at 2.859e+05 - 3.880e+06. Margins 2.58x below, 1.69x above.
* **G2 dynamic range** - any real transcriptome spans orders of magnitude. Require log2 range >= 8
  (256-fold). The 5 span 0.15-7.26; the other 119 span 9.00-15.26. **No overlap.**

Each alone flags exactly the same 5 with 0 false positives. Kept as two because they fail
differently - G1 catches a mis-scaled library, G2 a collapsed distribution.

**Does NOT catch** seven further columns that look poor on `mean-min` but pass both, including
`Y1_Fib` - whose library and range are normal and whose downgrade to "not established" stands.
Recorded as open.

### The consequence that makes this more than one assertion

**`N2_Fib` is N2's ONLY control.** Rejecting it leaves N2 with none, and `aging.py:88` then
self-centres - forcing N2's mean dAge toward 0 silently. **That is the exact behaviour Stage 1.5's
Group D pinned as a defect.** Rejecting the sample without deciding the donor trades a known-bad
control for a silent fallback, which is worse because the first is visible.

Recommended: **reject the sample AND the donor** - a donor with no sound control cannot carry
control-relative dAge. Corpus goes 124 -> 100 Gill columns, **6 donors -> 5**.

Which reaches further: **LOOCV goes 6 folds to 5**, so every Stage 1 guard and step 6's MDE must be
re-reported; it reaches **C-2** (N2 is donor age 0, so C-2's "two neonatal donors" becomes one plus
HFF); and it reaches **§4.7**, whose 16.67 yr spread is *defined over the six folds including N2's*.

> **Sequencing, and it is not negotiable: C-7 is WRITTEN now and ADOPTED after 3b and 3c report.**
> Writing costs nothing and cannot be undone by their results; adopting first would delete their
> evidence rather than answer their question.

### Bars

B1 separation (exactly 5 flagged, 0 of 119); B2 **no silent fallback** - a donor losing its last
control must **raise**, per the `age_label_policy` fail-open precedent; B3 the gate can fail, both
branches executing in tests; B4 **bit-identical when off** - ships off, enabled by its own run
exactly as C-2 did. All four are **deterministic** classifications on a fixed matrix, so
`bar_verdict` records *resolvability N/A* rather than simulating a null - claiming a power
calculation here would be theatre.

### Not licensed

C-7 does **not** re-measure §§4.5/4.6/1c/1d, which were computed on a contaminated `sigma_gill`
(§5.11) - that costs one HFF stream and belongs to 3b/3c. Nothing is withdrawn. And it does not
decide whether the defect is GEO's deposit or our read of it - **worth ten minutes against the GEO
supplementary file before the gate is implemented**, since if it is our read C-7 becomes redundant
rather than wrong.

---


## 2026-08-08 - 5.10 accepted; and N2_Fib is inside every number 1.5.6 has produced

**Status:** Recorded. **`src/` untouched, no label moved, §5.9 and §5.10 unmodified.**

### 5.10's corrections verified, and both errors are mine

Next SOUND column is `N3_d11_SSEA4` at **0.2967** (rank 13); rank 12 is `Y1_Fib` at 0.2745, flagged.
**Margin 1.08x, no gap** - §5.9 claimed 3.4x by comparing `Y1_Fib` against the six CONTROLS while
the screen ran on 124 COLUMNS. Two populations, one number, the same error species as A5. `Y1_Fib`
**downgraded**: library 1.51e+06 and range 14.43 both sit inside the sound population. The
five-column / one-control / 4.37x library separation reproduces exactly. **Gate moves to the library
tell.**

### THE EXPOSURE nobody had stated

`diag_harmonization_gain.py:112` selects **all six** `_Fib_` controls, `N2_Fib` included. So
`sigma_gill` was inflated by the degenerate column in **every** measurement built on it: the gain
**2.152** (1c), the top-100 gain **2.769** and the whole k-sweep (1d), and HFF day-14 **-21.43 /
-29.70**. `diag_harmonizer_refit_sparse.py` inherits it (written, never run).

Same contamination in the recorded builds: `N2_Fib` is `is_control` for `gill_bulk`, so it enters
`sigma_ref` in **five of six folds** - including **O1**, July's -24.02 reference and the anchor
`d_f/d_O1` is normalised against throughout 5.3-5.10.

**So 4.6's headline - "sparsification makes the gain worse, 2.769 vs 2.152" - was computed on a
contaminated sigma_gill.** Both share the contamination so the comparison may survive; the absolute
values do not, and whether the ordering survives is a measurement nobody has made. **Nothing
withdrawn** - what is recorded is that the inputs are known-contaminated and the re-measurement
costs one HFF stream, the same one step 3c needs.

Not touched: §1's MAE 16.61 -> 5.36 (transient arm, different matrix), and §5.6's G0, which is a
bit-exact reproduction of the SHIPPED harmonizer - a fidelity gate, not a correctness one.

---


## 2026-08-08 - 5.7 confirmed, and the defect is larger: 12 columns and TWO controls

**Status:** Verified independently. **`src/` untouched, no label moved, §5.7 unmodified.** Census
over the raw GSE165176 matrix before any transform - no run, no HFF stream.

### A correction I owe on 5.6

§5.7 is right that **the containment interval cannot falsify.** `d_hat_f/d_hat_O1` is a
`delta*r*w`-weighted average of rho, the clock's weights are near-balanced in sign (2648 +,
2589 -), so the average is **not** bounded by [min rho, max rho]. T2's bound was valid because d is
affine in one scalar; T3's is not. `diag_t3_sigma_gill_leverage.py`'s docstring says exactly this -
and its output then prints an "inside?" column reading YES and a verdict reading T3_STILL_LIVE. **A
stated limit that the output contradicts is a defect, not a caveat.** The Spearman +1.000 is the
real evidence in §5.6; the interval is not evidence at all.

### 5.7's root cause reproduces exactly

N2_Fib_Sendai_Exp2: mean-min **0.0008**, range **1.74**, **99.7%** of genes at the column min,
linear RPM sum **1.03e+08**. Sound controls: mean-min 0.96-1.56, range 13.1-14.5, RPM sum
1.48-1.66e+06. These are reads per million, so a sound column must sum to ~1e6 - the 68x library is
the mechanical tell.

**A screen that does NOT work, recorded so nobody retries it:** "fraction of genes at the column
min" flags **all 124 columns** at >50%, because ~60% zero-inflation is normal here and every
zero-count gene lands on the log2 floor. `mean - min` discriminates; %at-min alone does not.

### The defect is larger than 5.7 states

Screening all 124 at `mean-min < 1/5 x cohort median` (median 1.4196) flags **12 columns**, and
**two are controls**:

| sample | mean-min | range | role |
|---|---:|---:|---|
| **Y1_d7_CD13_Sendai_Exp1** | **0.0000** | **0.15** | entirely constant, library 2.15e+09 |
| N3_d21_SSEA4_Sendai_Exp2 | 0.0004 | 2.47 | |
| O2_d9_SSEA4_Sendai_Exp1 | 0.0005 | 2.15 | |
| **N2_Fib_Sendai_Exp2** | **0.0008** | 1.74 | **CONTROL** - 5.7's finding |
| N2_d21_CD13_Sendai_Exp2 | 0.0142 | 7.26 | |
| O2_d40 / O2_d34 / O1_d34 | 0.025-0.090 | 9.1-9.8 | |
| N2_d11_CD13_Exp2 / Y2_d34 / O1_d11_CD13_Exp2 | 0.097-0.103 | 9.0-10.1 | |
| **Y1_Fib_Sendai_Exp2** | **0.2745** | 14.43 | **CONTROL - a second one** |

Clean separation: every flagged column is below 0.284, the next sound one is 0.964 - a 3.4x gap.

### The second control explains what 5.5 left dangling

§5.5 flagged Y1's floor ratio as "genuinely anomalous, 34% low" with normal labels and left it
unexplained. **Y1_Fib is the second defective control.** And §5.6's own table carried the signature:
the two folds with the lowest |w|-weighted rho are **N2 (0.870)** and **Y1 (0.915)** - exactly the
two folds that REMOVE a defective control, ranked in the same order as the severity of the defect.

**So §5.6's "N2 is the atypical donor" reading is superseded. It was never the donor, and never
donor age 0 - it is the sample.** Spearman +1.000 was reading contamination, not biology.

### Ten degenerate NON-control columns are in the build as Gill training labels

`gill_bulk` is a training source; its 124 columns carry age labels into the corpus.
`Y1_d7_CD13_Sendai_Exp1` is entirely constant across ~20k genes. **`apply_qc` passes all of them** -
its gates are `min_genes` and `max_mito_frac`, designed for single cells, and a constant bulk column
clears both trivially. **This does not reach §1's headline** (MAE 16.61 -> 5.36 is the *transient*
arm, a different matrix); it reaches anything scored on the Sendai arm.

### The gate this argues for

**`mean - min` in log2 space per bulk sample, floored at 1/5 of the cohort median.** Flags 12/124
with a 3.4x margin and no false positives. Equivalently: linear RPM must sum to ~1e6. Stage 1.5.2's
G-a made `n = 1` **visible**; nothing checks whether that `n = 1` is **sound**. One assertion wide.

**Not established:** that removing these columns reproduces d/d_O1 = 0.306. That is step 3c.

---


## 2026-08-08 - T3 survives the test that killed T2; my averaging argument is dead

**Status:** Measured. **`src/` untouched, no label moved, §5.3/§5.5 unmodified.** New read-only
script `experiments/diag_t3_sigma_gill_leverage.py` - **no HFF stream needed**, Gill's six control
samples and the shipped O1 harmonizer are enough.

### My §5.4 prediction was wrong and the precheck killed it

§5.4 predicted T2 > T3. §5.5 measured the floor ratios: **anti-aligned**. N2, the fold that
collapses, is 1.4% off O1; Y1's ratio is 34% low with entirely normal labels. Max-leverage then
eliminated T2 outright - **-24.35 predicted vs -7.35 actual, a 17.00 yr miss on a 16.671 yr
spread.** Withdrawn. It cost two numbers per fold and died on the first fold it touched.

### G0 passes BIT-EXACTLY

Recomputing `sigma_gill` from the raw Gill matrix with **O1 held out**, floored, against
`runs/cellfate_multi/harmonization.json`: **2664 unclamped genes, median relative error 0.0000,
p90 0.0000, 5328/5328 genes aligned.**

Three things follow, none assumed: the shipped artifact **is** the O1 fold (proven by
reconstruction, not inferred from gene count); `sources.py:417`'s control definition is the
pipeline's; and the Gill side of the reconstruction is validated end to end.

### T3-only counterfactual - T3 is NOT eliminated

| fold | observed | d_f/d_O1 | rho range on clock genes | \|w\|-mean | inside? |
|---|---:|---:|---|---:|:---:|
| **N2** | **-7.35** | **0.306** | [0.108, 1.766] | **0.870** | yes |
| Y1 | -22.05 | 0.918 | [0.542, 1.926] | 0.915 | yes |
| N3 | -22.12 | 0.921 | [0.144, 2.021] | 0.987 | yes |
| O2 | -22.89 | 0.953 | [0.426, 1.999] | 0.995 | yes |
| Y2 | -23.87 | 0.994 | [0.427, 2.007] | 0.999 | yes |
| O1 | -24.02 | 1.000 | [1.000, 1.000] | 1.000 | yes |

Every fold's observation lies inside T3's containment interval; N2's reaches to **0.108** against
the 0.306 it needs.

**Spearman(|w|-weighted rho, observed d_f/d_O1) = +1.000 over 6 folds.** T3's leverage rank-orders
all six folds exactly as the labels are. Magnitude is under-predicted (0.870 vs 0.306), which is
exactly what §5.2 A1 says a |w|-weighted mean must do since it drops delta. The ordering statistic
is scale-free, so A1 does not reach it. **Right channel, wrong estimator.**

### Why the averaging argument was wrong

Dropping one of five controls removes **one entire donor profile**, so every gene's sigma responds
to a **single shared latent** - which donor was dropped. The perturbation is coherent across genes
in exactly the way I attributed to the floor. It also explains the **N2/N3 asymmetry** that donor
age alone could not (both are age 0): in the N3 fold, N2 is still in the control set and still
carries the spread.

### Consequence for §5.5's inversion

§5.5 argued the residue points outside harmonization because T3 should average out. **That premise
is now measured false.** The live hypothesis is T3 - the five-sample sigma_ref estimate - whose
leakage-free remedies are shrinkage or a reference change, and whose replication remedy is
donor-blocked. **NARROWS the field and shows the ordering is right; does NOT establish magnitude**,
which needs delta and therefore the HFF stream.

---


## 2026-08-08 - Audit of 5.3, and the variance floor turns out to BE the transform's centre

**Status:** Recorded. **`src/` untouched, no label moved, §5.3 unmodified.** Nothing here required a
run - every number comes from the code or from `runs/cellfate_multi/harmonization.json`, already on
disk.

### 5.3's structural fact is right, and stronger than its own citation

§5.3 credits `Harmonizer.fit`. The real guarantee is `fit_harmonizer` (`build_dataset.py:352-354`):
controls are selected by `is_ctrl & ~cell_line.isin(heldout)` - *"decidable from cell_line alone,
**before the full split**"*. The train/val/calib split never reaches the harmonizer, so HFF's
control set is **bit-identical** across folds, not merely "the same donors".

**Sharpening:** HFF's *admissible set* is fold-invariant too, and `genes_G` is the intersection. So
`G^(f)` moves through Gill's side alone - **T1, T2 and T3 are three channels out of one five-sample
estimate**, not three independent sources.

`F` is well-conditioned - §5.3 redefines `spread` as max-min in years, so §5.1's divide-by-median
degeneracy does not carry over. Checked, no issue.

### G0 compares two different quantities

`diag_harmonizer_refit_sparse.py` selects **all six** `_Fib_` controls with no fold exclusion, so its
regime A is an all-six-donor fit, **not fold O1** (which uses five). A 0.5 yr gate between them can
fail for reasons unrelated to fidelity. Extending it per fold fixes the quantity but forfeits
independence - both sides then run one code path.

**The independent reference already exists.** `runs/cellfate_multi/harmonization.json` **is the O1
fold's shipped harmonizer** - confirmed twice: 5328 genes matches §4.7's O1 column, and the held-out
donor has 21 split entries where Y1 has 19. It carries the pipeline's own per-gene mu and post-floor
sigma, so reconstructing against it validates T1/T2/T3 **gene by gene against ground truth**.

### THE MEASUREMENT - the floor is not a second-order term

| | |
|---|---|
| floor, gill_bulk / hff_sc | **0.15821** / **0.42388** |
| **floor ratio** | **0.3732** |
| clamped at the floor | **2664/5328 = 50.0%** in each dataset (mechanical - floor = median sigma) |
| clamped in **BOTH** | **1848 = 34.7%**, ratio **exactly 0.3732**, min = max |
| **median ratio, all 5328 genes** | **0.3732 - the median gene's ratio IS the floor constant** |
| median ratio, clamped in neither | 0.5335 |

**More than a third of the genes carry one identical ratio, and it is the median of the whole
distribution.** The floor sets the transform's central tendency rather than trimming its tail.

### Why this reorders 5.3's ladder before it runs

A shift in `floor_gill/floor_hff` between folds moves **1848 genes in lockstep** - coherent and
non-averaging. T3's per-gene sigma noise is large but *independent across genes*, so it largely
averages out in a weighted sum over thousands. **Predicts T2 > T3, and it is checkable from each
fold's two floor scalars alone, before any reconstruction is written.**

That raises the stakes on §5.3's asymmetry point: T1/T2 have leakage-free fixes, T3 largely does
not. If T2 also dominates, the instability is **the cheapest to fix and not donor-blocked** - moving
it off Stage 6's critical path. **Stated as a prediction with a mechanism, not a result.**

### Checked and clean

`split_sizes = {852, 115, 117, 21}` against `n_samples = 42605` reads like 41,500 unassigned cells.
It is not - those are **split-map entries** (1105: HFF 981, ~21 per Gill donor, Y1 19). No defect.
Recorded so nobody re-checks it.

---


## 2026-08-08 - CORRECTION: A5 withdrawn, A2 downgraded

**Status:** Corrected after review by the second machine. **`src/` untouched, no label moved. A1-A5
are left in place and marked, not deleted** - the wrong version stays visible as the record, per the
same convention being applied to §5.1.

### A5 - WITHDRAWN. My misreading.

§4.7 says *"the 99.8% quoted elsewhere in this document is the earlier **33,688-cell corpus**"*.
33,688 is the **corpus** size; §4's 33,613 is **HFF's** count. Checked: 33,613 / 33,688 =
**99.777% = 99.8%**. Consistent. I read a sentence about a corpus as a second claim about HFF's
count. There was no contradiction.

### A2 - DOWNGRADED. The elimination claim was an extrapolation.

A2 argued the deconfounder cannot carry the fold spread because step 1b bounded S3+S4 at +2.11 yr.
**Step 1b ran on the unharmonized path** - its own docstring opens *"Harmonization is OFF (no config
sets `harmonize: true`)"*, written before §4.3 established `HARMONIZE = True` in the real build. So
+2.11 yr is the deconfounder's contribution to a ΔAge ~2.15x smaller than the one it actually sees,
and S3's real inputs are harmonized values. Using it to bound **fold-to-fold variance** in a
harmonized build is an extrapolation, not a proof.

**It is the same species of error as ruling out harmonization from the absence of `harmonize: true`
in the YAML - committed against the very measurement that corrected it.**

### What survives, and why the correction strengthens it

The elimination framing falls; the **instrument** does not, and it is now better motivated. If the
deconfounder's fold behaviour is **unknown** rather than bounded, a statistic that presumes
harmonization is the mechanism is the wrong tool, while a **reconstruction** separates them by
construction: where harmonization-only reconstruction tracks `d_f`, harmonization is attributed;
**the residue where it fails to track is where the deconfounder or anything else lives.** A
measurement, not an assumption that the answer is known.

A1 (the statistic), A3 (the leaky remedy), A4 (the variance floor / gene set) stand unchanged.

### Repo state

`HEAD == origin/main == acf6ef0`, 0 ahead / 0 behind - **nothing to pull.** The second machine has
not pulled `acf6ef0` yet; that is his side, not a divergence.

---


## 2026-08-07 - Audit of the incoming step 3b, and C-2 is not free

**Status:** Recorded. **`src/` untouched, no label moved, no existing plan section edited.**

Pulled 9 commits from the other machine (arm D, the cross-machine reproduction, and 1.5.6 step 3b).
Audited step 3b and re-checked the standing shortlist against the code and the raw data.

### Step 3b - 4 defects and 1 inconsistency (`STAGE_1_5_6_SPARSE_CLOCK.md` §5.2)

| | |
|---|---|
| **A1** | `G_f` is a |w|-weighted mean of ratios - the family §4.5 disproved, because it drops `δ_g`. And the exact gain makes `d_f/G_f` fold-invariant **by algebra**, so `R = 0` always. Proxy measures the proxy; exact measures nothing |
| **A2** | The question is settled by **elimination** - only `σ_gill`, the gene set, and the deconfounder can move across folds, and step 1b bounded the deconfounder at +2.11 yr. The right instrument is a per-fold **reconstruction**, which cannot return "six folds cannot decide it" |
| **A3** | The ATTRIBUTED branch's remedy - *"fit-once-and-freeze"* - **reintroduces donor leakage**, the exact thing `harmonize.py:59-60` and Group C exist to prevent. Four non-leaking candidates recorded instead |
| **A4** | §4.7's own table lists 5026-5402 admissible genes per fold, then attributes everything to `σ_ref`. But `harmonize.py:112` floors σ at `median(sigma)` - a **set-level** statistic - so the gene set moves the floor for **every gene at once**. N2 has the fewest genes and is the fold that collapses. Untested |
| **A5** | §4 says 33,613 HFF labels; §4.7 says 33,688. One is wrong |

Verified and **correct**: five-single-controls (confirmed from the raw matrix - exactly 6 `_Fib_`
day-0 baselines among 124 columns), the 16.67 yr arithmetic, and `_clock_range`, which looks like a
fail-open but is not - `LinearClock.from_json` lifts `meta.age_range` and returns (1.0, 96.0).

### Two items our own 1b-1d work left unrecorded

* **U1** - `diag_harmonization_gain.py` picks Gill's controls with a regex whose failure mode is
  *unparseable -> control*. 118 of 124 names carry `_dNN_`, **none** carry `_d0_`, and the 6 that do
  not parse are the fibroblast baselines. **1c/1d used the right rows by accident.** The script is
  left **unmodified** so those results stay reproducible; the defect is recorded in the plan.
* **U2** - §4.6's option 2 can only move the **variance floor** and the **admissible mask**; `mu_g`
  and `sigma_g` are per-gene and do not move. **So option 2 and step 3b are the same lever**, and
  running them separately spends two efforts on one question.

### C-2 is not free - it masks 99.8% of the age labels (`00_START_HERE.md`)

The shortlist has carried *"Turn on C-2 - free, worth more than 13 donors"* for weeks. Checked:

* clock `meta.age_range` = **[1.0, 96.0]**
* Gill donors **N2 = 0, N3 = 0**, O1 53, O2 53, Y1 29, Y2 35
* **HFF `DONOR_AGE_YEARS = 0.0`** (`sources.py:557`), and **C-3 shipped**, so it is stamped
  (`sources.py:731`)
* rule 3 excludes on `donor_age < lo` (`aging.py:200`)

**Enabling it masks all 42,481 HFF cells plus N2 and N3.** The justifying comment at
`build_dataset.py:123` counts *"30 of the 75 non-HFF labels"* - correct when written, before C-3
reached HFF. **No bug: the flag is off and the guard held. The shortlist's cost estimate is what is
wrong.**

C-2 firing on HFF may well be **correct** - every HFF label is the clock extrapolated below its own
fitted range, which is a candidate explanation for why HFF has resisted 1.5.2 through 1.5.6. But
that makes C-2 **downstream of acquisition**, inverting its position in the shortlist. Recorded with
a proposed amendment to the order; **the order itself is left unedited pending that decision.**

Also noted: the fold that collapses (N2, -7.35) is a donor-age-0 donor - but N3 is age 0 and does
not collapse, so this is a **lead, not an explanation**.

### Shortlist status

| # | item | status |
|---|---|---|
| 1 | C-2 | 🔴 not free - after acquisition, not before |
| 2 | GSE165177 | ✅ free, already ordered |
| 3 | Stratified shuffle | ✅ **DONE** - arm D landed |
| 4 | Acquire | ✅ Stage 6 - and now C-2's prerequisite |

Not audited line by line: arm D's result and the reproduction claims. Plan text and the C-2 /
harmonizer facts only.

---


## 2026-08-08 (later) - STAGE 3a RUN: STOP, on clean labels and contaminated ones alike

**Status:** RUN, read-only. Goal reached: C-7 implemented, adopted on the dataset, and Stage 3a
executed on the resulting clean labels.

Run on TWO arms because 3a's STOP branch is TERMINAL - "do not write tool code, ship the scoring
model, go to Stage 5" must not rest on labels 5.14 proved were inflated ~3x in five of six folds.
`_c7` operative, `_armA` context. If they had disagreed, the disagreement would itself have been
the finding.

    Part A dAge MAE state-only     _c7 20.67   _armA 25.53
    Part A state + Dt              _c7 86.77   _armA 90.54
    Part C unsafe-frac state-only  _c7 0.409   _armA 0.386
    Part C state + Dt              _c7 2.010   _armA 0.983
    Part C paired   _c7 +1.601 CI [-2.270,+5.472]   _armA +0.598 CI [-0.821,+2.016]
    VERDICT                        _c7 STOP    _armA STOP

Cleaning the labels improved the state-only fit by 19 percent and did NOT change the verdict.

Dt does not merely fail to help, it makes the fit WORSE - about 65 MAE in Part A and 0.6 to 1.6 in
Part C. Both CIs include zero so both are formally tied, but the point estimate is consistently the
wrong way. That is section 0.4's own diagnosis: along one trajectory, time is redundant with state.

**PART B IS NOT EVIDENCE and must not be read as passing the 2-yr sweep clause.** Its swing is 255
to 269 yr, and in `_armA` the per-fold swing is identical to the digit, -255.76, in all six folds. A
real per-fold forward response cannot be identical across six different held-out donors. That is the
ridge extrapolating on the Dt-squared term, not signal. Part C is decisive and test 18 says so in
its own output.

`_c7` runs n=5 in Part A because N2's dAge is masked by C-7 rule 4, exactly as intended.

CONSEQUENCE: 3b, 3c and 3d are NOT written. This routes to shipping the scoring model and going to
Stage 5. The fate head is untouched and works - ROC 0.983, PR-AUC 0.992 - and that is what shipping
the scoring model means.

NOT DECIDED: that forward prediction is impossible in principle. Test 18's own caveat governs - a
screen, not a proof - and it is decisive for THIS corpus, at roughly 12 timepoints across 6 donors.
The two options section 3a names, a dAge-trajectory readout with no safety recommendation, or data
with more unsafe-cell variation via Stage 6, are scope decisions. Recorded, not taken.

---

## 2026-08-08 (later) - C-7 adopted on the DATASET: A1-A4 all pass, 5.14's prediction confirmed

**Status:** Six `_c7` folds built (6380 s), dataset only - no training, no bundle. `src/`
unchanged since e6fc183. Nothing adopted for trained-model consumers.

### A bug caught first, by checking the artefact rather than the wiring

The first build produced `cellfate_loocv_N2_c7` with `n_samples=42605` - byte-identical to arm A
- and N2's baseline census still reading `n_control=1, source=controls`. **The gate had done
nothing.** `run()` does `sources = sources if sources is not None else build_sources(cfg)`, and
the flag was set only inside `build_sources`; every real caller INJECTS sources
(`run_multi_local.py`, the C-7 driver, every test with a synthetic source), so `build_sources`
never ran and the flag never reached the source.

It would not have failed loudly. It would have produced six "C-7" folds identical to arm A and a
Stage 3a verdict reported as taken on clean labels. Fixed with `apply_source_flags(cfg, sources)`
called from BOTH `build_sources` and `run`, plus a test pinning exactly this. (e6fc183)

### The gate then fired identically in all six folds

    rejected: N2_Fib_Sendai_Exp2, N2_d21_CD13_Sendai_Exp2, N3_d21_SSEA4_Sendai_Exp2,
              O2_d9_SSEA4_Sendai_Exp1, Y1_d7_CD13_Sendai_Exp1
    lines_without_controls: ['N2']
    n_samples: 42605 -> 42600

  fold  HFF day-14 _c7   HFF day-14 _armA
  N2         -7.337            -7.352
  N3         -5.628           -22.121
  O1         -6.514           -24.023
  O2         -6.674           -22.891
  Y1         -4.606           -22.049
  Y2         -8.292           -23.869

  A1  N2's dAge fully masked, every fold, reason `no_control_baseline`   PASS
  A2  donor and fold survive -- 6 folds, N2's cells present              PASS
  A3  mean |day14 - (-8.196)| = 1.719 yr against 5.14's reconstruction   PASS
  A4  spread 16.671 -> 3.686 yr, a 78% reduction                         PASS

### The strongest evidence is one nobody designed

**The N2 fold barely moved: -7.352 -> -7.337.** It is the one fold whose harmonizer ALREADY
excluded the degenerate control, so C-7 had almost nothing to remove there. The other five, which
all contained it, moved from -22..-24 down to -4.6..-8.3 - **they moved to join N2**. The fold
that looked anomalous was the clean one all along, exactly as 5.7 and 5.14 predicted, now
measured on built labels rather than reconstructed.

### What it does NOT establish

That ~-6.5 is the CORRECT HFF dAge. It is the UNCONTAMINATED one; whether it is right is still
what arms C and D narrowed and did not close. It does sit close to the raw single-dataset value
(-8.5 to -10.6), consistent with 5.14's implication that the "harmonization gain" was
substantially this defect.

**And the residual is not zero:** 3.686 yr of spread remains. C-7 removed 78%; something still
moves HFF's labels a few years across folds, and that is not investigated here.

N2 now contributes 19 samples not 21 (two of its columns rejected), dAge masked, so it feeds the
FATE head only - precisely what option (c) was chosen to preserve.

---

## 2026-08-08 (later) - C-7 IMPLEMENTED (components A-D), ships OFF

**Status:** IMPLEMENTED, flag OFF. No label moved, no build touched, nothing adopted.
25 new tests, full suite green, ruff clean.

### The four components, landed together

Any three of them leave the system in a state C-7's own bars forbid, so they ship as one change.

**A. `src/cellfate/data/integrity.py`** - pure, no I/O, no config. G1 (linear RPM sum in
[1e5,1e7], because the matrix IS RPM and must sum to ~1e6 by definition of the units) AND G2
(log2 dynamic range >= 8, i.e. >= 256-fold). Both unit-justified, not fitted to this cohort's
quantiles - the mean-min-at-1/5-median proposal was rejected in 5.10 for cutting a continuous
distribution 8% from its neighbour. Verified against the real matrix AND the committed
124-column census: exactly 5 rejected, 0 false positives, each condition independently
rejecting all five.

**B. `GillReprogrammingSource`** gates in `_load()`, the single place the matrix is read, so
`plan()` and all THREE `src.fetch` call sites (build_dataset.py:170, :289, :345) are covered by
one edit. Missing :345 would leave the degenerate column inside `sigma_ref`, which is the entire
defect being fixed. Plus `lines_without_controls()`, GLOBAL by construction because `_load`
already holds every column - so there is no second corpus pass, and the chunk-local case cannot
be confused for the global one.

**C. `age_label_policy` rule 4 `no_control_baseline`** - a `cell_line` with zero admissible
controls has no zero-point, so its dAge is UNDEFINED and is masked. Keyed on the condition, not
on identity, so it works for every future dataset and closes Stage 1.5 Group D generally. Placed
AFTER `cancer_source` and BEFORE `donor_out_of_clock_range`: N2 is donor_age 0 so it matches both
once C-2 activates, and `age_mask_reason` is a persisted parquet column (io.py:139, :265), so the
order decides what is written to the shard. "No zero-point exists" is undefined; "outside the
fitted range" is out-of-validity. Undefined is stronger and is what gets recorded.

**D. `assert_no_unmasked_fallback`** - B2', reading the census `_control_baseline` already writes
(`"source": "self_fallback"`, gate G-a). And `recenter_on_control_arrays` now ACCEPTS a census:
it is `_control_baseline`'s SECOND call site (the S4 re-centring) and it passed none, so its
fallback was invisible. A B2' guarding only `delta_age` would have passed while S4 silently
self-centred the same orphaned line.

**ONE FLAG** (`DataConfig.bulk_integrity_gate`), set centrally in `build_sources`, so the gate
and rule 4 cannot be switched independently - the gate alone would strip the control and leave
the dAge unmasked, which is the state B2' exists to forbid.

### A design error the test suite caught, before any rebuild

I first made B2' UNCONDITIONAL. That broke
`test_the_silent_no_control_fallback_self_centres_a_line_to_zero`, which pins today's Group E
behaviour and says in terms that changing it must be "a deliberate, reviewed act". Making the
assertion unconditional WAS such a change, made without review. B2' is now gated on the flag, so
the flag-off path is untouched. **That is B4 doing its job at unit-test speed rather than after a
5-hour rebuild.**

### Verified end to end on the real Gill matrix

  gate OFF -> 124 samples, 6 donors, nothing rejected, empty census
  gate ON  -> 119 samples, 6 donors STILL, rejects exactly the five, and
              lines_without_controls == {'N2'}

Option (c) working as designed: **the donor and the fold survive**; only N2's zero-point goes,
and rule 4 masks its 21 labels. C-7 section 5's "re-report over 5 folds" stays corrected to 6.

### Tests

12 gate tests (B1 separation on the recorded 124-column cohort, margins, each condition
independently rejecting, B3a/B3b both branches) + 13 rule-4/B2' tests (B3c, **B3d**, C1 ordering,
B2'a/b/c including the S4 site, B4 defaults). **B3d is the load-bearing one**: a chunk with no
controls whose line HAS them globally must NOT trip rule 4, or C-7 blocks on Stage 1.5 Group E
for the wrong reason.

One test expectation was wrong and the code right: a constant column at a HIGH value fails G1
first, because 36k genes at 11.489547 gives a library of 1.03e+08 - N2_Fib's actual signature.
It now has its own test, and the G2 test pins the library in-band so it grades G2 specifically.

### NOT done

Adoption. The flag is off; enabling it is a separate pre-registered run with its own snapshot and
a full Stage 1 guard re-report over 6 folds. That is the retrain 5.13 freed by cancelling 1.5.6
step 4.

---

## 2026-08-08 (later) - C-7 section 9 verified: (c) agreed, B2' is an AMENDMENT, and the predicate's host is wrong

**Status:** VERIFIED and RECORDED. No run. `src/` untouched, no build touched, no label moved.
C-7 section 10 appended (83 insertions, 0 deletions; section 9 unmodified).

### Verified and agreed

The donor error table is exact - recomputed here from the clock on each day-0 control: N2 +98.65,
N3 +36.44, Y1 +35.92, O2 +26.50, O1 +26.12, Y2 +22.66. And the framing is stronger than my
"+35.12 above the mean": the clock is biased HIGH on every donor (+22.66 to +36.44), and **N2 is
2.71x the next worst**. fate_roc 0.983 "untouched by every dAge problem" confirmed.

**Option (c) agreed.** **Rule 4 general rather than donor-named agreed - and their version is better
than mine.** I keyed on IDENTITY (masked_cell_lines); they key on the CONDITION (zero admissible
controls => no zero-point => dAge undefined), which fires exactly where the fallback fires, works
for every future dataset, and closes Stage 1.5 Group D generally.

**The two-fallback distinction is correct and sharp:** _control_baseline falls back "when a line has
no controls IN THIS CHUNK", and Stage 1.5 Group E is the chunk-local case, so B2' must test the
GLOBAL predicate or it fires on Group E and blocks C-7 for the wrong reason. **Rule 4 ships WITH
C-7** - the gate alone creates the orphaned line. **SRA becomes an optional upgrade** under (c).

### Bar discipline: B2' is an AMENDMENT and should be labelled one

Section 9 says the B2 collision "dissolves rather than blocks". The substance is right; the framing
understates it. B2 as pre-registered says a donor losing its last control "must RAISE". B2' replaces
that with "may fall back if masked" - a change to the bar's TEST, not a reading of it, though
faithful to the bar's own TITLE ("no silent fallback"). Under 5b a bar may be amended before the run
with the reason recorded, which is exactly this situation, so the amendment is legitimate and
agreed. It should be recorded as "B2 -> B2', amended 2026-08-08, reason: the original conflated the
mechanism (raise) with the invariant (no unmasked fallback label)". This project has been bitten
four times by bars that moved without the move being labelled.

### BLOCKING CORRECTION: the predicate cannot live in fit_harmonizer's pre-pass

Section 9 says the global predicate is "decidable in the pre-pass fit_harmonizer already runs". It
is not, for two independent reasons, both in the code:

1. `harmonizer = fit_harmonizer(cfg, work) if cfg.harmonize else None` (build_dataset.py:383). The
   pre-pass runs ONLY when harmonization is on. A rule-4 predicate hosted there would silently not
   exist in any harmonize=False build - and those are real (arm B/C/D probes, any single-dataset
   build). A data-integrity invariant that evaporates when a flag is off is a guard that cannot
   fire, which is the exact defect class this project keeps catching.
2. `controls.setdefault(str(ds), ...)` pools per DATASET_ID, not per cell_line. The loop has
   obs["cell_line"] in hand so the tally is easy to add, but "already runs" overstates it: the loop
   runs, the tally does not.

Where it belongs: `work = plan_all(sources)` and `load_or_fit_panel(cfg, work)` both run
UNCONDITIONALLY (build_dataset.py:379-382). The global per-cell_line control census belongs in an
unconditional pass over `work`, independent of fit_harmonizer - which also makes it available to
harmonize=False builds and to G-a's existing baseline_census, whose job is already "what each dAge
zero-point actually rests on".

**Blocking for the implementation, not for the decision.** Option (c), rule 4's general form, the
two-fallback distinction and the ship-together sequencing all stand. Only the predicate's HOST
changes.

---

## 2026-08-08 (later) - POST-CLOSURE: order changes to C-7 BEFORE Stage 3a; one overreach corrected

**Status:** RECORDED. No run. `src/` untouched, no build touched, no label moved. Section 5.16
added; section 5.14's one overreach flagged in place (3 lines changed, not deleted).

### 1. My overreach, accepted

5.14 attributed the +2.73 yr baseline gap to S3+S4 "near the magnitude step 1b measured (+2.11 yr)".
**Step 1b ran on the UNHARMONIZED path** - invoking it is the same extrapolation A2 was downgraded
for, repeated two sections after recording why it was wrong. Corrected in place. Nothing downstream
moves: the +2.73 is a scale sanity-check and the verdict rests on the Delta column.

### 2. Scalar-vs-reconstruction discrepancy RESOLVED, and it confirms A1

The second machine ran the O1-fold leave-one-control-out on the Gill side alone: |w|-weighted ratio
N2 0.904, Y1 0.795, N3/O2/Y2 1.006/1.009/1.030. That reproduces 5.7's -11.8%/-20.1% and puts Y1
AHEAD of N2 - the opposite order to 3c. **Not tension: A1 confirming itself.** The scalar drops
delta, and dAge is a near-cancelling signed sum (2648 positive weights, 2589 negative), so a ~10%
per-gene ratio change yields a 69% swing in the sum. 5.10's mixed-sign argument and this are the
same fact from opposite directions.

**And the sensitivity is not generic, which is why 3c survives it:** if it were mere numerical
fragility every drop would swing. Four of five move 0.32-2.98 yr; only N2's moves 18.558. The
near-cancellation AMPLIFIES a real specific perturbation; it does not manufacture one.

### 3. THE ORDER CHANGES: C-7 before Stage 3a. My proposal is withdrawn.

I proposed Stage 3a next, run twice (six donors, and again with N2 excluded) so the defect could not
decide it. **Withdrawn - it was wrong, and worse than I understood.**

test18_forward_gate.py targets "mean TRUE dAge at t_j", read from fold folders and POOLED ACROSS
FOLDS - the quantity 5.14 proved is inflated ~3x in five of six folds. And 3a's STOP branch is
TERMINAL: "do not write tool code. Ship the scoring model; go to Stage 5." A genuine Dt signal could
sit under a 16.67 yr between-fold label artefact and read as STOP.

**And the two-arm design makes it worse, not better.** Excluding donor N2 removes N2's ROWS, but
N2's CONTROL still sits in the harmonizer of every fold that does not hold N2 out. So it drops the
N2 fold - THE ONLY FOLD WHOSE HARMONIZER IS CLEAN - and keeps the five contaminated ones. Backwards.

C-7's blocker has cleared: its section 3 reads "C-7 is WRITTEN now and ADOPTED after 3b and 3c
report"; 3c reported ATTRIBUTED and 3b is unnecessary (5.14). And 5.13 cancelled step 4, freeing the
retrain budget adoption needs. **Order: C-7 (adopt) -> then Stage 3a on clean labels.**

### 4. A THIRD donor option, neither side had costed

C-7 section 3 poses two: (a) drop donor N2, (b) re-quantify N2_Fib from SRA SRP302546. Option (a)
conflates three separable decisions - should the degenerate CONTROL enter the harmonizer (no,
definitively), should N2's own 21 dAge LABELS survive (no - its zero-point reads 98.65 yr from a
constant vector, highest of six for a donor of age 0, an unrecoverable +35.12 yr offset), and should
N2's CELLS survive at all (not obviously - the fate head consumes no dAge, works at ROC 0.983, and
donors are THE binding constraint; dropping N2 also removes one of only two age-0 donors).

**Option (c): reject the control and mask N2's dAge labels, keep the donor and the fold.** Removes
the contaminant everywhere, discards only labels already known garbage, keeps SIX folds - so C-7's
guard re-report stays at 6 rather than the 5 its section 5 currently assumes.

Cost stated honestly: age_label_policy keys on source, masked_datasets and donor_age
(aging.py:135-176) - none can address a single donor by cell_line. Option (c) needs a fourth rule
(masked_cell_lines), small and following C-1's own precedent, but it IS a src/ change and therefore
its own Change.

**And it collides with B2 - which is the design question.** Rejecting N2's control leaves cell-line
N2 with ZERO controls, and _control_baseline "falls back to the line's own mean when a line has no
controls in this chunk" - the silent self-centring B2 exists to forbid. So under (c) B2 must fire
unless masking also prevents the baseline being requested for that line. Whether it does is an
implementation question to settle BEFORE (c) is chosen; C-7 section 5 already routes a B2 failure to
"blocking. The donor-level decision is not optional; fix that first."

Recorded as an option to cost, not adopted. The choice among (a), (b), (c) belongs to C-7.

---

## 2026-08-08 (later) - STAGE 1.5.6 CLOSED: step 3c ATTRIBUTED, step 4 cancelled, no retrain spent

**Status:** CLOSED. **No retrain spent.** `src/` untouched, no build touched, no label moved,
`configs/clocks/` untouched. Sections 5.14 (3c result) and 5.15 (closure) added; the status line
appended to, not replaced. 126 insertions / 0 deletions in the plan.

### Step 3c RUN -> ATTRIBUTED

G0 first: the five-control fit reproduces the shipped sigma_gill on all 2664 unclamped genes at
median relative error 6.14e-09. Scale confirmed - baseline d_hat -26.755 against the recorded
d_O1 -24.023, a +2.73 yr difference which is the S3+S4 contribution (step 1b measured +2.11). HFF
stream 5782 day-14 cells, 5981 day-0.

  dropped   d_hat      Delta
  N2        -8.196   +18.558
  N3       -27.174    -0.419
  O2       -28.321    -1.566
  Y1       -27.078    -0.323
  Y2       -29.735    -2.980

The four healthy drops all move it the SAME way, more negative, by 0.32-2.98 yr. N2's drop moves it
the OPPOSITE way, +18.558 yr, at 6.23x the largest healthy drop. Not the biggest of a family - not a
member of it.

B1 outlier 6.23x (bar >=2x) PASS. B2 gap closed 1.113 (bar >=0.70) PASS. B3 direction +18.558 PASS.
Section 5.8 stated in advance that exact reproduction of the N2 fold's -7.352 was NOT predicted; it
nearly happens anyway (-8.196).

### The consequence that reaches furthest

Removing one column drops the harmonized magnitude 69.4%. Beside step 1d's own figures: unharmonized
direct dense -9.96; harmonized dense -21.43 ("gain 2.152"); harmonized WITHOUT the degenerate
control -8.196. **With the contaminant gone the harmonized value lands on the unharmonized one** -
residual gain ~0.82, essentially none. So sections 4.3-4.6's "harmonization gain x2.152" is
substantially NOT a property of the transform but one degenerate control inflating sigma_gill.
Strong implication, not a closed result: 1c/1d floored over the clock genes while this floors over
the pipeline's genes_G (defect U1), so the figures are not exactly comparable.

### Step 2 DECIDED, and it cancelled step 4

Estimand stays BOTH clocks -> the sparse clock fails (MAE 8.79 vs <=8, sign 0.41-0.68 vs >=0.80, at
every k) and is NOT adopted. Narrowing to multi-tissue is unavailable: Stage 1.5.2 refuted it with
arithmetic (lambda_mt = 1.048, and a correlation loading cannot exceed 1). Steps 3, 4, 5 do not run.
**Step 4's retrain CANCELLED, not deferred.** Step 3b is UNNECESSARY - the mechanism its ladder was
designed to find is now named and measured.

### What leaves the stage

1. The clock's -14.10 yr density bias: real, removable at k~100 (MAE 16.61 -> 5.36), CONFINED to
   multi-tissue dAge on the transient arm. Not adopted.
2. N2_Fib_Sendai_Exp2 is a degenerate GEO column and the carrier of HFF's fold instability,
   inflating the harmonized magnitude ~3x in five of six folds. **The N2 fold is the clean one; the
   other five, including O1 = July's -24.02 reference, carry inflated labels.**

### Not established

That -8.2 is the CORRECT HFF dAge. This establishes that -24 is CONTAMINATED. Whether the
uncontaminated value is right is what arms C and D narrowed and did not close. **No recorded result
is withdrawn** - arms A-D, step 6 and July's reference stand, now with a named defect inside them.

### Superseded leads, annotated in 00_START_HERE rather than deleted

Both candidate explanations for N2 are answered and neither was the mechanism: donor age 0 (N3 is
also 0 and does not collapse) and fewest admissible genes / the variance-floor lever (section 5.5
eliminated the floor at maximum leverage, F = 1.398; and Y1, whose floor ratio IS anomalous at 34%
low, has entirely normal labels). It is the SAMPLE.

### Handed forward with owners

C-7 owns the fix (exclude or re-quantify from SRA SRP302546) and the bulk-integrity gate, plus the
ten degenerate non-control columns. C-2's activation moves to AFTER Stage 6 (masks ~99.8% of the
corpus). "Is HFF's dAge CORRECT age?" remains open.

ruff clean; tests/test_results_paths.py 178 pass - it caught the new script's _RESULTS spelling
before commit.

---

## 2026-08-08 (later) - STEP 2 DECIDED: estimand stays BOTH clocks, sparse clock FAILS, step 4 CANCELLED

**Status:** DECIDED on recorded evidence. No run, no retrain. **`src/` untouched, no build touched,
no label moved, `configs/clocks/` untouched.** Section 5.13 added; the step 2/3/4/5 rows marked.

### The decision

**The estimand is dAge agreement with BOTH methylation references. The sparse clock does not meet it
and is NOT adopted. Steps 3, 4 and 5 do not run.**

### Why it cannot narrow to multi-tissue

Section 5.12 showed step 2's bar fails on skin & blood by both clauses at every k, leaving only one
escape: narrow the estimand. **Stage 1.5.2 already examined that move and refuted it with
arithmetic** (added 2026-08-01). If one shared age factor were all three instruments measured, each
correlation is a product of loadings:

  lambda_RNA = sqrt(0.267 x 0.516 / 0.568) = 0.493
  lambda_sb  = 0.267 / 0.493               = 0.542
  lambda_mt  = 0.516 / 0.493               = 1.048   <- a correlation loading CANNOT exceed 1

The three numbers are not jointly consistent with one common factor. 1.5.2's own words: "RNA<->
multi-tissue reaching 91% of the ceiling is not evidence the RNA clock nearly works; paired with
RNA<->skin-and-blood at 47%, it is evidence it does not." And its standing rule: "both clocks must
agree or it is SPLIT." Adopting the narrow reading inside 1.5.6 would overturn a closed verdict by
convenience rather than evidence.

### Why no further data can change it

The inventory (STAGE_6_NEW_DATA_REV section 2) holds TWO methylation instruments and no third -
GSE165178 (Sendai, 22/22) and GSE165179 (transient, 68 conditions). Both are already used: section
1's headline is multi-tissue on the transient arm, section 3's failure is skin & blood on the same.
There is no unused reference clock and no in-range donor set that re-opens the question.

### What it does NOT do

Does not withdraw section 1's finding. The -14.10 yr density bias is real; removing it at k~100
takes MAE 16.61 -> 5.36 with spread preserved (ratio 1.04) and LODO generalisation intact, and
section 0's "the dense clock is worse than predicting a constant" (16.61 vs an 8.45 zero-floor)
stands. It CONFINES that finding to multi-tissue dAge on the transient arm - one clock, one arm, one
estimand - which is what section 0's "what survives, precisely" table already said.

### Step dispositions

- 2 DECIDED. Bar not registered: a bar no correct system can clear is UNRESOLVABLE under 5b.
- 3 NOT DONE deliberately. Shipping fleischer_clock_top100.json would build an artefact for a
  candidate rejected at step 2. "Free" is not "warranted".
- 4 **CANCELLED, not deferred.** The rebuild existed to measure the sparse clock end to end; with
  the candidate rejected there is nothing for it to decide. **One retrain not spent.**
- 5 MOOT - no label change to decide.
- 3b / 3c UNAFFECTED, and reframed: they were never about the clock. They measure why HFF's labels
  move 16.67 yr across folds, which reaches arms A-D, step 6 and July's reference regardless.

Section 4.6's four-way option table resolves permanently to "do neither" - its own listed honest
position - because the option it was choosing among is not adopted.

### Re-homed on closure

- C-2's activation (STAGE_1_5_3_EXECUTE assigns it to 1.5.6) -> AFTER Stage 6. 00_START_HERE records
  that enabling it masks ~99.8% of the corpus; it is downstream of acquisition.
- The clock's -14.10 yr density bias -> recorded and closed as a measured property of the Fleischer
  clock. Not fixed here, and not fixable without changing the estimand.

---

## 2026-08-08 (later) - Step 2's bar is already measured to FAIL; the "both arms" ambiguity resolves itself

**Status:** RECORDED, no run. **`src/` untouched, no build touched, no label moved.** Section 5.12
added; step 2's table row flagged. 73 insertions / 1 deletion - the deletion is that row.

Raised by the second machine, verified here against sections 0 and 3 before recording.

### Both clauses of step 2's bar are already measured on skin & blood

Step 2: "MAE <= 8 yr AND sign agreement >= 0.80 vs methylation, on BOTH arms, k fixed at 100."

  MAE            bar <= 8      measured 8.79 at k=100   (section 0, ERROR 3)
  sign agreement bar >= 0.80   measured 0.41-0.68       (section 3, "sometimes below chance")

Section 3 forecloses the escape: "Horvath skin & blood does not come right at any k." At k=50 the
MAE clause passes (6.69) while sign agreement still fails. Under REF_GROUND_RULES 5b a bar no
correct system can clear is UNRESOLVABLE and must move BEFORE the run. Registering it as written
would be the 5b violation this project has now caught four times - knowingly, against numbers
already on the page.

### The ambiguity resolves itself - one reading is impossible

"Both arms" could mean the protocol arms (transient/Sendai) or the reference clocks (skin&blood/
multi-tissue). **It cannot mean the protocol arms:** section 1.3 records that GSE165178 has no
untreated control, so dAge CANNOT BE FORMED on the Sendai arm - it is scored on absolute age. A dAge
bar "on both protocol arms" is impossible by construction, not merely unmet. Therefore "both arms"
means both methylation clocks, and under that reading the bar fails on skin & blood by both clauses
at every k.

### Narrowing the estimand is not a free escape

The obvious repair - restate the bar as multi-tissue only - is bigger than a bar edit, because
section 3 already ruled: "This does not clear M-2a's SPLIT rule. A variant passing on one clock and
not the other is still a SPLIT." Narrowing re-opens Stage 1.5.2's M-2a verdict, which is closed. So
step 2's decision reaches OUT of this stage and cannot be granted inside it.

### Step 2's real job, and three honest outcomes

Not "write down a threshold" - a SCOPE decision. (1) Estimand stays BOTH clocks: the bar is right
and the sparse clock already fails it, so steps 3/4/5 do not happen and section 1's finding is
recorded as real but confined to multi-tissue dAge on the transient arm - which is what section 0's
"what survives, precisely" table already says. (2) Estimand narrows to MULTI-TISSUE: defensible only
if M-2a's SPLIT rule is explicitly amended in Stage 1.5.2 with its own justification; step 2 blocks
on that. (3) Something else entirely (e.g. bias removal rather than agreement): must be stated and
pre-registered before step 3. Recorded, not decided.

### Consequence for the working order

Steps 2 and 3 were listed as free and independent of 3b/3c. **Step 3 still is** - writing
fleischer_clock_top100.json ships a new file and switches nothing. **Step 2 is not**: it carries a
scope decision reaching into Stage 1.5.2, and its bar must not be registered in its current form.
The order 3c -> 3b is unaffected; step 4 stays blocked on both.

### Agreed from the same exchange, verified

- C-7 (plans/CHANGE_C7_BULK_SAMPLE_INTEGRITY.md) written and NOT adopted until 3b/3c report. The
  reasoning is right and worth restating: adopting it first drops donor N2, which deletes the fold
  section 4.7's 16.67 yr spread is DEFINED OVER - it would remove the question rather than answer it.
  Taking a separate change ID rather than a section number is also the right call after the 5.8
  collision.
- The GEO deposit is the source, not our read: same file GEO serves, gzip intact, 136 fields on all
  35,806 lines so no column can shift, awk agrees with pandas. N2_Fib takes FOUR distinct strings
  across 35,805 genes with 11.489547 covering 99.7%.
- The SRA route (SRP302546) is a genuinely better option than dropping the donor - it restores N2's
  zero-point at no cost in donors, folds, or guard re-reports. Their caveat is right and should stay
  attached: one sample re-quantified through a different pipeline than the other 123 introduces a
  batch term to CHECK, not assume.
- Rejected gate candidate recorded so neither machine retries it: DISTINCT-VALUE COUNT does not
  separate. The five flagged carry 2, 4, 5, 5 and 27 distinct values; the other 119 carry 22-693 -
  N2_d21_CD13 (flagged, 27) overlaps sound O2_d40 (22). Library and dynamic range separate cleanly;
  this does not.

---

## 2026-08-08 (later) - Section 5.9's separation claim does not reproduce; Y1_Fib downgraded; numbering collision fixed

**Status:** VERIFIED against the raw matrix. **`src/` untouched, no build touched, no label moved.**
140 insertions / 1 deletion in the plan - the single deletion is a duplicate section header.

### Numbering collision fixed first

Two machines pushed a `## 5.8` concurrently (78fd8a9 = step 3c's pre-registration, 09bc61f = the
audit of 5.7) and neither saw the other's. The audit is renumbered **5.8 -> 5.9**. **Only the header
number changed; every word of its body is as written.** The step-3b and step-3c table rows point at
5.8, which now resolves unambiguously.

### What reproduces from 5.9

12 of 124 columns flagged at mean-min < 1/5 x cohort median (median 1.4196, threshold 0.2839) - same
twelve. The %at-min screen is useless (124/124 exceed 50%), a genuinely useful negative result.
N2_Fib to the digit. The ten degenerate NON-control columns are live gill_bulk training labels;
Y1_d7_CD13_Sendai_Exp1 is entirely constant (range 0.15, library 2.15e+09) and apply_qc passes it
because min_genes and max_mito_frac are single-cell gates a constant bulk column clears trivially.
And 5.9's framing of the containment interval as a DEFECT rather than a caveat is the right call.

### What does NOT reproduce

5.9: "Every flagged column is below 0.284; the next sound one is 0.964. A 3.4x gap, no overlap."

**The next sound column is N3_d11_SSEA4 at 0.2967, not 0.964. The margin is 1.08x and there is no
gap** - the threshold falls between two adjacent values 8% apart. 0.964 is O1_Fib, the next sound
CONTROL. The screen ran on 124 columns; the separation was computed against 6. Different populations.

### Y1_Fib as a second defective control - NOT ESTABLISHED

Library 1.51e+06 (sound population 8.47e+05 - 2.29e+06) and log2 range 14.43 (sound 13.1-14.5) are
both entirely normal. Only mean-min is low, and it is lowest of six controls - but n=6, and against
124 columns it is unremarkable. Its max is the HIGHEST of any control (16.600), which reads more
like a library dominated by a few transcripts than a constant column. Possibly a quality concern;
not the defect N2_Fib has.

### What survives: FIVE columns, ONE control

Flagged by BOTH mechanical tells (mean-min < 0.015 AND library >= 1.69e+07): Y1_d7_CD13 2.15e+09,
N3_d21_SSEA4 2.38e+08, O2_d9_SSEA4 1.63e+08, **N2_Fib 1.03e+08 (CONTROL)**, N2_d21_CD13 1.69e+07.
Next column by library is 3.88e+06 - a **4.4x gap**. That is the real separation; it is just not
where the mean-min threshold put it. Section 5.7 unchanged.

### 5.9's supersession of 5.6 held at PARTIAL

It requires Y1_Fib to be defective, and rests on Y1 (0.918) vs N3 (0.921) in observed d/d_O1 - a
0.003 margin, the same fragile pair flagged when 5.6's Spearman was checked. WHAT STANDS: 5.6's
"N2 is the atypical donor / donor age 0" reading IS superseded; N2_Fib is a sample defect. WHAT DOES
NOT: that Y1 is the second instance. Recorded as "N2_Fib explains N2; Y1 remains unexplained".

### The gate moved to the library tell

5.9's proposed mean-min gate is "validated on this cohort - 12/124, 3.4x margin, no false
positives". That validation does not hold, and a threshold cutting a continuous distribution 8% from
its neighbour will flag or miss arbitrarily on a new cohort. Use 5.9's own stronger half instead:
**assert each bulk sample's linear RPM sum lies within a stated band of 1e6** - justified by the
matrix's own units rather than a quantile of this cohort, and separating the five by 4.4x. One
supporting figure corrected: "sound columns land at 1.48-1.66e6" is the six CONTROLS' range; over
all 124 the sound range is 8.47e+05 - 2.29e+06, so the band must be set from the full population.

### Step 3c UNCHANGED, deliberately

5.9's covering note asks 3c to test removing the TWO defective controls jointly. Declined: (1) 3c
already answers it - its design is leave-one-CONTROL-out over ALL FIVE O1-fold controls
individually, and Y1_Fib is one of them, so if it matters 3c returns it as a second outlier, B1's
">= 2x the second largest" fails, and the run routes to PARTIAL by measurement; (2) joint removal
would destroy the separability that makes the question answerable and would bake in an unestablished
premise. No change to 5.8 and none needed.

Carried: 5.9's reconciliation request stands and 3c settles it - 5.7's leave-one-out figures
(-11.8% N2, -20.1% Y1) are over all genes unweighted while 5.6's rho is |w|-weighted on clock genes
and orders them the other way. Different quantities; 3c computes actual dAge and adjudicates.

---

## 2026-08-08 (later) - Step 3c added: does removing the degenerate control reproduce the spread?

**Status:** PLAN ONLY. Nothing executed. **`src/` untouched, no build touched, no label moved.**
107 insertions / 1 deletion in the plan - the single deletion is step 3b's row, repointed to note
it runs only if 3c does not settle it.

Turns section 5.7's "what is NOT established" into an actual pre-registered step rather than a
paragraph of open questions.

### Why a new step and not a rung of 3b's ladder

5.7 established by direct measurement that N2_Fib_Sendai_Exp2 is nearly a constant vector and that
it enters sigma_gill in every fold that does not hold N2 out. That is a NAMED candidate. 3b's ladder
was designed when the candidate was a diffuse property of a five-sample estimate. Testing a named
contaminant is cheaper and strictly more decisive, so 3c runs first and may make 3b unnecessary.

What 5.7 did not settle: whether removing that column reproduces the magnitude. A scalar sigma
argument does not deliver it - leave-one-donor-out moves the |w|-weighted sigma_gill -11.8% dropping
N2 but -20.1% dropping Y1, and Y1's labels are normal.

### The test

On the O1 fold (July's reference, d_O1 = -24.023, and one of the five folds that INCLUDE the
contaminant), refit the harmonizer five times, each dropping one of its five controls, everything
else held exactly, and recompute HFF's day-14 dAge each time. MIN_REPLICATES = 3, so four controls
remain legal in every arm. Needs the HFF stream, because delta is required for a magnitude and 5.6
established the mixed-sign weights (2648 +, 2589 -) make every delta-free shortcut unbounded - the
same machinery 3b needs, used once instead of eighteen times.

**The built-in negative control is the point.** Dropping ANY control changes sigma_gill. The claim
is not that N2's removal moves the number but that it is an OUTLIER among the five, and the four
healthy drops are measured in the same run by the same code on the same fold.

### Bar

- B1 outlier (primary): |Delta_N2| is the largest of the five AND >= 2x the second largest.
- B2 magnitude (primary): gap closed A = Delta_N2 / (d_N2 - d_O1) >= 0.70, denominator +16.671 yr.
- B3 direction (gate): Delta_N2 > 0. Removing the contaminant must move O1's dAge TOWARD zero;
  the wrong sign falsifies the mechanism outright.
- bar_verdict resolvability run BEFORE the measurement, bar moves to usable_bar first if
  UNRESOLVABLE (REF_GROUND_RULES 5b).

Stated in advance so it cannot be read as a shortfall later: "O1 minus N2's control" is NOT the N2
fold - the N2 fold holds out N2 and therefore includes O1's control, on a different admissible set.
Exact reproduction of -7.352 is not predicted and is not the bar.

Branches: ATTRIBUTED (3b becomes unnecessary) / PARTIAL (3b runs on the residue) / GENERIC (the
contaminant is not the carrier, 3b runs as written) / FALSIFIED (wrong sign - record it, do not
rescue it).

### Why this matters more than the ladder

5.3's remedy table ranked fixes by leakage-safety, with replication donor-blocked behind
STAGE_6_NEW_DATA_REV section 3 / D2. A degenerate INPUT SAMPLE sits outside that table: excluding or
repairing one bad column never touches the held-out donor (leakage-free), needs no new donors (not
donor-blocked), and is a data-ingest fix rather than a change to the estimator. **If 3c attributes,
the instability is the cheapest thing on the page to fix and comes off Stage 6's critical path.**

### Carried with 3c but not part of its measurement

- 3c.2 - is the defect GEO's deposited matrix or our read of it? Decides whether the fix is exclude
  or re-read. This stage, reported beside 3c.
- 3c.3 - apply_qc passed a control constant to 1.74 log2 with a library 68x the cohort. A guard is a
  src/ change and therefore its own Change with its own bar. NOT this stage; named so it is not lost.
- 3c.4 - five further degenerate Gill TREATMENT samples. NOT this stage; they do not enter
  sigma_gill but they do enter their own dAge and any prior Gill analysis.

Includes a falsifiability self-test: synthetic controls with no degenerate sample must return
GENERIC; one with a planted constant column must return ATTRIBUTED. Both branches execute in tests.

---

## 2026-08-08 (later) - ROOT CAUSE: a Gill CONTROL sample is degenerate in the raw GEO matrix

**Status:** RUN, read-only, raw GEO only. **`src/` untouched, no build touched, no label moved.**
Added `experiments/diag_gill_control_integrity.py` + results JSON; section 5.7 added to
STAGE_1_5_6_SPARSE_CLOCK.md (5.3/5.5/5.6 all untouched).

Found while verifying the second machine's 5.6. Sections 5.5 and 5.6 both computed statistics OVER
sigma_gill; neither asked whether the six control samples it is fitted on are sound.

### N2_Fib_Sendai_Exp2 is nearly a constant vector

Raw Log2 RPM before any transform: min = median = mean = 11.490, max 13.227. Dynamic range **1.74
log2 units** where every other control spans 13-15, and a mean sitting 0.0008 above its own floor -
which real RNA-seq cannot do. Implied library after 2**x - 1: **1.03e+08** vs ~1.5e+06 for all five
others (68x). log1p-CP10k profile SD 0.011 vs ~0.58. Rank agreement with the other five controls
**0.096** vs 0.50-0.69.

Six of 124 Gill columns are degenerate; exactly one is a control (N2's). All 20 of N2's other
samples are normal - one bad sample, not a bad donor. The other five are treatment samples
(N2_d21_CD13, N3_d21_SSEA4, O2_d40_SSEA4, O2_d9_SSEA4, Y1_d7_CD13; Y1_d7_CD13's range is 0.152).

### Why it does not stay in its own donor

The day-0 _Fib_ sample is is_control (sources.py:417), so it is BOTH N2's entire dAge zero-point in
every fold AND one of the controls sigma_gill is fitted on - and sigma_gill/sigma_hff is the gain
applied to HFF's labels, 99.7% of the corpus, in every fold that does not hold N2 out. A
near-constant column inflates sigma_gill; removing it deflates it. Exactly one fold removes it -
N2's - and that is the fold reading -7.35 against the others' -22 to -24.

**The fold that looked anomalous is the one whose harmonizer is clean.** The five agreeing folds
agree because they share a contaminant, including O1 - July's -24.02 reference. Agreement across
folds is not corroboration when the folds share the defect.

### Established vs not

Established by direct measurement: the six degenerate columns, N2's control among them, its 0.096
rank agreement, N2's other 20 samples normal. NOT established: that it accounts for the whole
16.67 yr spread - leave-one-out on sigma_gill moves the |w|-weighted mean -11.8% dropping N2 but
-20.1% dropping Y1, and Y1 is normal. The reconstruction is still required. Also open: whether the
defect is GEO's deposit or our read of it.

### Gate gap

apply_qc runs on every fetched chunk and this column passed it. Stage 1.5's G-a made n=1 controls
VISIBLE; nothing checks whether that n=1 is SOUND.

### On 5.6, verified independently

G0 passes BIT-EXACTLY (median rel err 6.1e-09, p90 5.6e-08, 5328/5328 aligned) - settles the
artifact identity by reconstruction and validates the Gill side end to end. The ordering statistic
survives stripping the O1 anchor and N2 (n=4, rho +1.000, exact p 0.042), though it hinges on the
Y1/N3 pair separated by 0.003 in observation. ONE CORRECTION: the containment interval cannot
falsify - d_hat_f/d_hat_O1 is a c_g-weighted average with c_g = delta*r*w, and the clock's weights
are near-balanced in sign (2648 +, 2589 -), so the average is NOT bounded by [min rho, max rho].
Their script says so; the summary reads it as a passed test. T2's bound WAS valid (d affine in one
scalar), which is why T2 could be and was eliminated. "T3 not eliminated" comes from a test that
could not eliminate it. The ordering is the real evidence.

Nothing is withdrawn on the strength of this entry. It names a defect and its reach.

---

## 2026-08-08 - Precheck RUN: the variance floor is NOT the carrier (section 5.5)

**Status:** RUN. **`src/` untouched, no build touched, read-only.** Added
`experiments/diag_fold_floor_precheck.py` + its results JSON; section 5.5 added to
STAGE_1_5_6_SPARSE_CLOCK.md (108 insertions, 0 deletions - 5.3 and 5.4 both untouched).

The second machine's section 5.4 measured that the variance floor is the transform's CENTRE, not a
tail-trim: floor_gill 0.15821, floor_hff 0.42388, ratio 0.3732, and 1848 of 5328 genes (34.7%)
clamped in BOTH datasets carry that one ratio to 12 dp - which is also the median ratio of the whole
distribution. From that it predicted **T2 > T3**, and noted the prediction is checkable from two
scalars per fold. Reproduced their measurement exactly, then ran the check.

### The prediction is NOT supported

Floor ratio by fold: N2 0.3784, N3 0.3674, O1 0.3732, O2 0.3676, Y1 0.2447, Y2 0.3725.
**Anti-aligned with the anomaly.** N2 - the fold that collapses to -7.35 - has a floor ratio 1.4%
off O1's. Y1 - whose floor ratio IS anomalous, 34% low - has entirely normal labels (-22.05).
Spearman(R, day-14) = -0.14.

Maximum-leverage bound: `d = sum(delta*r*w)` is affine in R_f with slope sum_B(delta*w) over the
clamped-both set, so `d_f = d_O1 * R_f/R_O1` is T2's ceiling. T2: F = 1.398, explains -39.8%, worst
miss 17.00 yr on a 16.671 yr spread -> **ELIMINATED at any leverage fraction**. T1 (clock-weight
coverage): F = 0.957, explains +4.3% -> **NOT A CARRIER**; N2 in fact carries more top-100 clock
genes than O1 (49 vs 48).

**Scope, stated in the script and the section:** R_f and C_f are SCALARS. This eliminates the
lockstep-constant and total-coverage channels. Floor effects acting through WHICH genes clamp, and
mask effects acting through gene IDENTITY, are not scalar and are NOT tested. The ladder is
narrowed, not closed.

### The inversion

Within harmonization only T3 (per-gene sigma_gill) survives - and T3 is the term section 5.4's own
averaging argument says should largely cancel over ~5000 weighted genes. If T1 and T2 are out by
measurement and T3 is expected to average out, the residue points OUTSIDE harmonization, at the
deconfounder, which A2's downgrade left unbounded. Argued, not shown: dropping 1 of 5 samples can
move a per-gene sigma a long way and those moves are not independent of gene identity. T3 must be
measured, not reasoned away. The reconstruction is now more necessary, not less.

### Corrections to 5.3, carried in 5.5 rather than by editing it

- **G0 replaced.** It gated against `diag_harmonizer_refit_sparse.py` regime A, which selects every
  `_Fib_` sample across ALL SIX donors while fold O1 fits on five - two different quantities. Now
  gates against `runs/cellfate_multi/harmonization.json`, the O1 fold's SHIPPED harmonizer
  (confirmed independently: its 5328 genes are unique to O1 among the six folds). Validates inputs
  gene by gene against ground truth instead of one scalar to 0.5 yr.
- **Citation corrected** to `fit_harmonizer` in build_dataset.py - `keep = is_ctrl & not_test`,
  decidable from cell_line alone before the split, so HFF's control set is BIT-IDENTICAL across
  folds, not merely the same donors.
- **Sharpening added:** `admissible[ds]` is per dataset from that dataset's own pooled controls
  (harmonize.py:87-88), so HFF's admissible set is fold-invariant too and G^(f) moves through Gill
  alone. T1/T2/T3 are three channels out of ONE five-sample estimate, not three sources.
- **Ladder reordered:** the reconstruction should spend its effort on T3 and the residue.
  NOT ATTRIBUTED is now the expected branch on current evidence.

### Housekeeping verified, recorded so nobody re-checks it

`split_sizes {852,115,117,21}` against `n_samples 42605` is not 41,500 unassigned cells - those are
split-map ENTRIES (1105 total: HFF 981, ~21 per Gill donor, Y1 19); every fold's holdout.json has a
1105-entry map. No defect. And `loocv_results/{folds,summary}.json` are tracked and committed on
main; the locally modified copies are run output, not a missing commit.

Ran ruff (clean) and tests/test_results_paths.py (171 pass) on the new script.

---

## 2026-08-08 - Step 3b rewritten as a RECONSTRUCTION (section 5.3 supersedes 5.1)

**Status:** PLAN ONLY. Nothing executed. **`src/` untouched, no build touched, no label moved.**
186 insertions / 1 deletion in the plan doc - the single deletion is the step-3b table row being
repointed and re-described. **Section 5.1 is left fully visible** under a supersession banner, the
same convention section 5.2 applied to its own withdrawn A5.

Written after the second machine's audit (section 5.2) and the correction round on it. Three of its
five findings stand and are acted on here; A5 was withdrawn by its author as a misreading and A2 was
downgraded from "settled by elimination" to "recommends the instrument".

### What was wrong with 5.1, verified independently before acting

- **A1 (fatal).** `G_f` was a |w|-weighted MEAN OF RATIOS, dropping the per-gene deltas - exactly
  the statistic section 4.5 disproved (median ratio 0.608 against net gain x2.152). And the exact
  gain does not rescue it: with `G_f = sum(d*r*w)/sum(d*w)`, `d_f/G_f == sum(d*w)`, fold-invariant
  BY ALGEBRA, so R = 0 for any data. Proxy measures the proxy; exact form measures nothing.
- **A3.** 5.1's ATTRIBUTED remedy was "fit-once-and-freeze", which REINTRODUCES LEAKAGE - every Gill
  donor has exactly one control, so no subset is in-train for all six folds. Struck.
- **A4.** 4.7 attributed the spread to sigma_ref alone and never mentioned the VARIANCE FLOOR.
  Verified at harmonize.py:112-113: `floor = median(sigma)` over the ADMISSIBLE genes - a set-level
  statistic - then `sigma = max(sigma, floor)`, which clamps roughly HALF the genes. Changing the
  gene set moves sigma for genes whose own sigma did not move. N2 has the fewest genes (5026) and is
  the fold that collapses.

### What section 5.3 replaces it with

A **reconstruction**: recompute `d_f` from each fold's own harmonizer inputs and check it reproduces
the recorded `d_f`. Where it tracks, harmonization is attributed; the residue where it fails to
track is where the deconfounder, or anything else, lives. Not algebraically degenerate - `d_hat` is
computed from raw inputs, never from `d`.

**The structural fact that sharpens it:** HFF's control cells are FOLD-INVARIANT (every fold holds
out a Gill donor). So `sigma_raw^hff` is the same number per gene in all six folds, and HFF's side
of the ratio can move through exactly two channels - the admissible mask and the floor. Any
HFF-side variation IS the mask or the floor.

**Three named terms, carried explicitly:** T1 mask (which genes enter), T2 variance floor (A4's
set-level median), T3 sigma_gill (the five-single-controls mechanism, 4.7's original claim).
Reported as an ablation ladder, with T1 and T2 stated as NON-ORTHOGONAL so the write-up cannot
present a clean variance decomposition it does not have.

**Gate G0 (fidelity):** the reconstruction at fold O1 must agree with
`diag_harmonizer_refit_sparse.py` regime A to <= 0.5 yr - two independent implementations of the
same quantity. Deliberately NOT a check against `d_O1`, which would bake in an assumed S3+S4
magnitude that A2's downgrade says is unknown.

**Primary metric:** `F = spread(d_f - d_hat_f) / spread(d_f)`, spread(d_f) = 16.671 yr.
Bar **F <= 0.25**, lower-is-better, `bar_verdict` resolvability run BEFORE the measurement.
Branches: ATTRIBUTED / PARTIAL / NOT ATTRIBUTED, with step 4 blocked on the first two and
proceeding with a stated nuisance term on the third.

### Two things folded in, as requested

1. **Section 4.6 option 2 is the same lever as A4.** `mu_g` (:110) and `sigma_g` (:111) are
   per-gene, so restricting the gene set moves neither; only the floor (:112) and the admissible
   mask (:88, :91) move. So "re-fit on the sparse gene space" is mechanically a change to the
   variance floor and nothing else - the lever the folds vary incidentally and option 2 varies
   deliberately. Carried as term T2. One measurement, not two.
2. **`experiments/diag_harmonizer_refit_sparse.py` is REUSED, not discarded** (304 lines, NOT YET
   RUN). It already selects Gill's controls by sources.py:417's definition with an asserted count
   and computes the three floor regimes with reconciliation to 1c/1d. 3b extends it per fold.

### The remedy table is now leakage-aware, and the asymmetry matters

T1 and T2 have **leakage-free** fixes (a floor or mask that is not a function of the current fold's
gene set never touches the held-out donor). T3 largely does not - shrinkage or a ref_dataset change
are the only non-leaking options, and more replicates is donor-blocked (STAGE_6_NEW_DATA_REV
section 3 / D2). **Which term carries the variance therefore decides whether this is cheap to fix or
data-blocked.** That is itself a reason to run the ladder.

### U1 carried forward, and a precision limit it implies

`diag_harmonization_gain.py` selected Gill's controls with a regex whose `else 0.0` default makes any
UNPARSEABLE name a control. It happened to be right - 118 of 124 names carry `_dNN_`, none carries
`_d0_`, and the 6 that fail to parse are exactly the day-0 dermal fibroblasts. **1c/1d got the right
6 rows by accident.** The script is left unmodified so 1c/1d stay reproducible.

Consequence recorded, nothing withdrawn: 1c/1d floored at the median over the CLOCK genes while the
pipeline floors over the full admissible space and applies an expression floor the clock-gene set
never had. **So 2.152 and 2.769 are near-pipeline gains, not the pipeline's own.** Regime A
reconciles them; 3b reports it.

---

## 2026-08-07 (later) - Stage 1.5.6: step 3b added, and it GATES step 4

**Status:** PLAN ONLY. Nothing executed. **`src/` untouched, no build touched, no label moved.**

Written after the matched HFF reproduction (entry above) found that HFF's day-14 dAge swings
**16.67 yr across the six LOOCV folds** - N2 reads -7.35 against a median of -22.51, a 3.1x
compression - even though HFF is never the held-out line and supplies 99.7% of the age-labelled
corpus.

### Why it belongs in 1.5.6 and NOT in a new stage

An earlier draft of this proposed a new Stage 1.5.7. That was wrong and is withdrawn. The check
belongs inside 1.5.6 for a concrete reason: **step 4 is a rebuild + LOOCV comparing sparse vs dense
clock across those same six folds.** Run as written, its paired CI carries a 16.67 yr nuisance term
and the retrain would be spent on a confounded comparison. So 3b is not a new question - it is a
gate step 4 already needed and nobody knew it. It is also the same subject: steps 1c and 1d already
measure the harmonization gain (2.152 dense, 2.769 sparse).

### The mechanism is predicted by this project's own closed form

Stage 1.5's audit replaced the false "batch-immune by construction" claim with
`dAge = sum_g delta_g * sigma_ref,g / (sigma_d,g + EPS) * w_g`. `Harmonizer.fit` takes training
control cells only -- "already excludes the held-out donor" (harmonize.py:59-60) -- with
`ref_dataset = "gill_bulk"`, and every Gill donor's zero-point is ONE unreplicated control (audit
section 5.2). So `sigma_ref` is estimated from **five single control samples**, and which five
changes every fold. HFF's entire label scale is multiplied by a factor estimated from five cells.

A mean-of-ratios does NOT explain it (N2 0.425 and Y1 0.410 are similar; only N2 collapses), which
is expected from section 4.5 - the gain is per-GENE. The right statistic is the **clock-weighted**
gain, and nobody has computed it per fold. That is step 3b.

### Added to plans/STAGE_1_5_6_SPARSE_CLOCK.md (140 insertions, 0 deletions)

- a dated pointer box above section 0
- **section 4.7** - the measured fold instability and its consequence for step 4
- a **3b** row in the section 5 plan table, between steps 3 and 4
- **section 5.1** - the full pre-registration: statistic `G_f` (clock-weighted gain), primary metric
  `R` (residual spread ratio), bar **R <= 0.35** graded through `bar_verdict` with a resolvability
  simulation required BEFORE the run, three fixed decision branches, an explicit
  what-this-does-not-license list, and a falsifiability self-test (the verify_1a lesson)
- two artefact rows in section 7

### The decision branches, fixed in advance

ATTRIBUTED -> step 4 is BLOCKED; fixing the fit protocol becomes its own Change with its own bar.
NOT ATTRIBUTED -> step 4 proceeds, the spread carried as a stated nuisance term, new owner needed.
UNRESOLVABLE -> step 4 proceeds, instability handed to STAGE_6_NEW_DATA_REV section 3 beside D2 -
the same n=1 control problem, already donor-blocked.

### Correction carried in the same breath

G-c step 2 is **not** an open retrain. It ran as step 6 - the snapshots are literally
`gc2_A_keep_hff.json` / `gc2_B_mask_hff.json` - and returned INCONCLUSIVE (MDE 5.045 vs delta*
3.572). It does not need a retrain; it needs lower variance, which is the same question 3b asks.
Step 3b explicitly does NOT re-open step 6; that would need its own pre-registration.

---

## 2026-08-07 (later) - Matched HFF reproduction: EXACT bit-for-bit; and HFF's labels are NOT fold-stable

**Status:** RUN and VERIFIED. **`src/` untouched, no build touched, read-only.** The July reference
`results/diag_gc_hff_signature_results.json` was preserved byte-exact (the July script writes that
fixed path, so each run's output was moved aside and the file restored from git; verified clean).

Follow-up to the entry above, which left the HFF half of the user's ask open. User: "now run the HFF
one on matched data."

### Correction to the previous entry's withdrawal reasoning - recorded, not rewritten

The entry above withdrew the HFF comparison as "1-dataset reference vs 2-dataset build". The
withdrawal was right; **the reason was not the operative cause.** There are TWO HFF references:

- `diag_pipeline_decompose_results.json` - raw GSE242423, ONE dataset, harmonization off, day-14
  -8.5 to -10.6. The user's point applies here and this comparison stays withdrawn.
- `diag_gc_hff_signature_results.json` - BUILT SHARDS of a TWO-dataset harmonized LOOCV build,
  day-14 **-24.02**. This is the number actually being compared against, and it was never
  single-dataset.

**The real cause was a FOLD mismatch:** arm A's `N2` fold read against a reference produced from the
`O1` fold (`diag_gc_hff_signature.py` defaults to `runs/cellfate_loocv_O1`). Which fold is read
turns out to matter enormously - see R2.

### Added

- `experiments/repro_hff_signature_armA.py` - imports `load_hff` and `trajectory_stats` from the
  July script UNMODIFIED so the arithmetic is identical on both sides; only the run dir changes.
  Pre-registered R1 (matched fold exact) and R2 (six folds agree within 2.0 yr / 0.30 yr per day).
  Deliberately never writes the July reference path.
- `results/repro_hff_signature_armA_results.json`,
  `results/diag_gc_hff_signature_armA_O1_results.json`.

### R1 - PASS. The matched fold reproduces July BIT-FOR-BIT

n_cells 37693, rho_timepoint -0.9047619047619048, slope -1.5255573306808494, rho_percell
-0.4160726187605165, slope_percell -1.5055164919911033, 5/7 descending, and `days` / `mean_dage` /
`sem_dage` / `n_per_day` elementwise identical. Also exact: `label_volume`, the `verdict` block, and
**all eight leave-one-timepoint-out folds**. Full float64 equality, not rounding agreement.

The G-c step-1 verdict is unchanged and re-confirmed on the current build: rho -0.905 PASSES the
<= -0.50 bar, slope -1.526 yr/day FAILS the [-6.45, -1.61] band, ambiguous => RUN_STEP_2.

### R2 - FAIL. A new finding: HFF's dAge labels depend on which donor is held out

| fold | day-14 dAge | rho | slope | descending | harmonized genes |
|---|---|---|---|---|---|
| **N2** | **-7.352** | -0.8095 | **-0.4885** | 3/7 | 5026 |
| N3 | -22.121 | -0.7857 | -1.2649 | 4/7 | 5258 |
| O1 | -24.023 | -0.9048 | -1.5256 | 5/7 | 5328 |
| O2 | -22.891 | -0.9048 | -1.4351 | 5/7 | 5304 |
| Y1 | -22.049 | -0.9048 | -1.3995 | 5/7 | 5402 |
| Y2 | -23.869 | -0.9048 | -1.4737 | 5/7 | 5305 |

day-14 spread **16.671 yr** against a 2.0 yr tolerance; slope spread 1.037 against 0.30. N2 is a
**3.1x compression** off a median of -22.506.

HFF supplies **42481 of 42605** age-labelled cells (99.7%) and is never the held-out line - only a
Gill donor is, worth ~21 cells. The training target for 99.7% of the corpus should not move when
0.05% of cells are withheld. It moves by 3x.

**Candidate mechanism, NOT established:** harmonization is refit per fold and the fitted gene set
varies (5026-5402; N2 fewest). The `gill_bulk` reference side is small enough that dropping a donor
may perturb its per-gene mu/sigma, and step 1c/1d already showed the gain is per-gene and largest
where the clock's weights sit. But N2 and Y1 have similar `gill_bulk` profiles and only N2's labels
collapse, so a single-variable story does not fit.

### What this does and does not change

Arms A/B/C/D all aggregate across these six folds, so the instability sits inside every step-6
number. It is a **candidate** source of the between-fold variance that made step 6 inconclusive
(SD 4.808, MDE 5.045 vs delta* 3.572) - a hypothesis, not a claim; nothing here measures its
contribution. **No previously recorded result is withdrawn on the strength of this.**

It does not touch the RES/ranking reproduction above: per-fold `model_dAge` and `ridge_dAge` were
identical to July on all six folds, so the *ranking* is robust to this instability even though the
label *magnitudes* are not.

### New open items

- Why does the N2 fold compress HFF's dAge 3x? Per-fold harmonization refit is the suspect.
- Does the fold instability explain step 6's SD 4.808?

Notebook: 118 insertions, **0 deletions** - the earlier entry was annotated, never rewritten.

---

## 2026-08-07 - Cross-machine reproduction: arm A reproduces July's ranking EXACTLY; RES differs exactly as pre-registered

**Status:** RUN and VERIFIED on the data machine. **`src/` untouched, no build touched, read-only.**

The user asked for a direct check: the other machine ran ΔAge/RES experiments in July, before Stage
1.5.6; re-run them on arm A and see whether they match. "A simple test that will tell us if it's
actually working or not."

### Half the ask was withdrawn before it ran - the user caught the flaw

An earlier attempt compared arm A's HFF ΔAge trajectory against `diag_gc_hff_signature`. The user
stopped it: *"arm A is from 2 data sets while the tests run on 1 data set - of course it would be
different results."* Correct. Those diagnostics run the clock on **raw single-dataset** input; arm A
is a **two-dataset harmonized** build. **The comparison is withdrawn as invalid - it measured
nothing**, and is recorded that way in the notebook rather than quietly dropped. The HFF half now
needs a dataset-matched re-run and is logged as open.

The RES/ranking half has no such confound: it runs on the **built** LOOCV artefacts, and both sides
are the same build. Verified from metadata first, reading the deleted July builds out of git at
`f353526^`: n_samples 42605, n_shards 51, splits 852/115/117/21, `gene_panel_hash`
`783f269a214aa972`, label_distribution 22635/1095/18875 - **identical on both sides**. Only
`baseline_census` differs, because it did not exist in July (added by G-a).

### Added

- `experiments/repro_test7_res_armA.py` - imports the three July scripts **unmodified** and only
  redirects `resolve_root` at the `_armA` roots. Pre-registered bar (P1-a/b/c primary, P2-a/b/c
  secondary, outcomes O1-O4) written into the header before any arm-A number was read.
- `experiments/diag_res_degenerate_armA.py` - separates the two mechanisms that can make RES a
  constant: the OOD gate zeroing every cell, versus `R_eff = max(0, -(mu + z*sigma))` flooring to 0.
- `results/repro_test7_res_armA_results.json`, `results/diag_res_degenerate_armA_results.json`.

### Result 1 - the ΔAge ranking reproduces EXACTLY, 6/6 folds, three decimals

`model_dAge` +0.910 / +0.909 / +0.990 / +0.970 / +0.960 / +0.947 and `ridge_dAge` +0.957 / +0.925 /
+0.960 / +0.952 / +0.951 / +0.983 - **every value identical to July**. Aggregates 0.948 / 0.955,
delta -0.000 on both. Test 7.1 exact too (gated 0.292/0.295, penalized 0.414/0.414, precision@5
0.27/0.30, same unsafe counts). Test 7.2's A arm exact at 0.955.

**This survived a different machine, a full rebuild, and Stages 1.5.3-1.5.6.**

### Result 2 - RES went degenerate, and Stage 1 Change A predicted it in writing

`model_RES` is a **constant** (zero variance, Spearman undefined) on N3/O2/Y1, and lower on the
other three (N2 0.742->0.222, O1 0.684->0.369, Y2 0.674->0.609). The diagnostic pins the mechanism:
on the constant folds `R_eff = 0` for **100%** of cells because `mu + z*sigma >= 0` throughout -
**not** the OOD gate. `sigma_age` went ~2.4 yr (July) -> **37.35 yr** (arm A), which is `sigma_scale`
from Stage 1 Change A doing exactly what it was built to do. Change A's own pre-registration says
per-cell RES *"should therefore approve nothing, and that is the correct result."*

Both Change-A invariance guards held exactly: `rank_model_dage` **0.948** (predicted "exactly
0.948"), `ood_rate` **0.2732** (predicted "exactly 0.273") - independent confirmation from a test
written for another purpose.

Direction of the July finding is intact and stronger. On the three folds where RES is defined,
paired (RES - ridge) = **-0.567 [-1.019, -0.115]**, excludes 0. Test 7.2, which isolates the RES
*formula* with ΔAge held constant, gives **-0.529 [-0.812, -0.246]** against July's -0.300 [-0.473,
-0.128].

### Two things recorded against myself

1. **The bar fired O3 "DOES NOT REPRODUCE" and that verdict is left standing as it fired.** I wrote
   it without first checking Stage 1 Change A's own RES prediction, which lives in the same notebook
   and says RES should approve nothing. Had I read it, P1-b/P2-b would have been registered as
   *expected to fail* and the outcome would have been O1. The record is annotated, not rewritten.
2. **`sigma_scale` came in off-prediction:** predicted ~5-6, measured 9.89-18.64 (mean 13.33). The
   prediction used mean-error / mean-spread (14/2.4); `sigma_scale_factor`
   (`xdonor_calib.py:399-403`) uses **P90(|residual|) / median(spread)**. The implementation matches
   its own documented arithmetic. Corollary worth carrying forward: with sigma scaled to a P90
   half-width and `z_conf = 1.0`, `mu +- 1*sigma` is a **90%** interval - `sigma_age` is not a 1-sigma
   and must not be read as one.

Also noted, deliberately not fixed: `test7_ranking.paired_ci` does not NaN-filter while 7.1 and 7.2
do, which is why P1-b's CI came back undefined rather than signed. The July scripts are the record
being reproduced and were not edited.

### Scope - what this does NOT establish

REPRODUCTION, not validation. It shows arm A ranks by ΔAge as well as July did and that RES still
degrades that ranking. It says **nothing** about whether the ΔAge labels are *correct age* - the
question arms C and D narrowed and did not close.

The July entries in `DELTAAGE_LAB_NOTEBOOK.md` were annotated with forward pointers and otherwise
left byte-for-byte intact (188 insertions, **0 deletions**).

---

## 2026-08-06 - Step 1d: the sparse clock and harmonization INTERACT ADVERSARIALLY on HFF

**Status:** MEASURED, and it runs AGAINST the sparse clock. **`src/` untouched, no label moved.**

Step 1c closed with a warning rather than a result: the gain is applied **per gene**, so a clock that
changes *which* genes carry ΔAge has its own gain, and neither §1's Gill-side number nor 1c's 2.152
transfers to the combination. Step 1d measures it instead of assuming it, sweeping k over the same
cells with the pipeline's own arithmetic and the same median variance floor.

| k | direct | harmonized | **gain** | median sigma-ratio on kept genes |
|---:|---:|---:|---:|---:|
| 50 | -8.22 | -28.17 | **3.429** | 0.838 |
| **100** | -10.73 | **-29.70** | **2.769** | 0.836 |
| 150 | -11.68 | -30.97 | 2.651 | 0.780 |
| 300 | -15.45 | -32.99 | 2.135 | 0.690 |
| 1000 | -12.78 | -24.83 | 1.944 | 0.608 |
| **all 33,155** | -9.96 | **-21.43** | **2.152** | 0.608 |

### The finding

**The sparse clock's gain is HIGHER than the dense clock's - 2.769 against 2.152 - and it rises
monotonically as k falls.** The mechanism is in the last column: the clock's largest-|weight| genes
sit precisely where `sigma_gill / sigma_hff` is largest (**0.836** among the top 100 against **0.608**
over all genes). **Sparsification concentrates the clock onto exactly the genes harmonization
amplifies most.**

Under the pipeline as it actually runs, `top100` gives HFF day-14 **-29.70 yr** - **further from
plausible than the dense clock's -21.43.**

### What this does to §1

> **§1's finding does not transfer to HFF. It inverts.** Sparsifying removes a -14 yr bias on Gill
> and *increases* the magnitude on HFF. Gill is the harmonization reference, so its gain is ~1 by
> construction and this effect is invisible there. The two changes do not compose - they compound in
> opposite directions on the two datasets.

§1's status is downgraded from "MEASURED, validated leave-one-donor-out" to **"MEASURED and
CONFINED"**. The leave-one-donor-out validation stands exactly as measured; what does not stand is
the extrapolation from it to HFF.

### The plan consequence

**A single clock cannot be adopted globally on this evidence.** The same change improves Gill and
degrades HFF, and **HFF is 99.8 % of the age labels.** Recorded in §4.6 with the four options and
why the honest position is to adopt nothing until step 4 measures the combination end-to-end.
**Step 4's rebuild was never optional; it is now the gate rather than a formality.**

This is why 1d was registered before it was run - it could only ever have confirmed or overturned
§1, and it overturned it for the dataset that carries the labels.

---


## 2026-08-04 — Step 1c: the harmonization gain is MEASURED at 2.152, and it REWEIGHTS rather than rescales

**Status:** ✅ Confirmed. **843 tests pass, ruff clean, `src/` untouched, no label moved.**

Step 1b eliminated the deconfounder (+2.11 yr, wrong direction). Step 1c measures the remaining
candidate using the pipeline's own arithmetic rather than an approximation: `transform` is
`(x − mu_d)/(sigma_d + EPS)` and `project_to_clock` is `x_scaled·sigma_ref + mu_ref`, and because
ΔAge is a **difference** both `mu` terms cancel exactly, leaving
`Σ_g w_g·(x_pert − x_ctrl)·sigma_gill/(sigma_hff + EPS)`. The variance floor is applied exactly as
`harmonize.py:112` does it (σ floored at the **median**) — omitting it would manufacture the very
gain being measured.

| | |
|---|---|
| HFF day-14 ΔAge, direct | **−9.96 yr** |
| HFF day-14 ΔAge, harmonized | **−21.43 yr** |
| **measured gain** | **2.152** — predicted **2.26** ✅ |
| recorded shard `y_age` | −24.02 |

Harmonization closes **~11.5 of the ~14 yr**. The residual ≈2.6 yr is within what differing cell
subsets (50,241 vs 37,693) and a deconfounder refit on harmonized data would move.

### 🔴 The finding that is worse than "a gain"

**Median σ ratio 0.608, mean 0.560** — for most genes `sigma_gill < sigma_hff`, which would *shrink*
ΔAge. The net effect is nonetheless **×2.152**.

> **Harmonization is not rescaling ΔAge — it is REWEIGHTING it.** The ratio is per gene, and the
> clock's heavy-weight genes sit where `sigma_gill/sigma_hff` is large. **A majority of genes are
> damped while an amplified minority carries the clock.**

Stronger than Stage 1.5 Group B's closed form implies on its face: Group B said `sigma_d` survives as
a per-dataset **gain**; measured, it survives as a per-**gene** reweighting whose net effect on this
clock is ×2.15, with the median gene pulling the opposite way.

### Consequence — 1.5.6's two effects do not compose additively

Sparsifying to the top-100 weights **changes which genes carry ΔAge**, and the gain is per gene. So
**a sparse clock has a different harmonization gain from the dense one**, and neither §1's number
(measured with Gill as reference, where the gain is ≈1 by construction) nor the 2.152 transfers to
the combination. Step 4's rebuild must measure the sparse clock's gain on HFF directly — added as
step **1d**.

**§1's Gill-side result cannot be extrapolated to HFF at all** — not because it is wrong, but because
Gill is the reference dataset and is the one place this effect is invisible.

---


## 2026-08-04 — Step 1b: the 13.4 yr is NOT the deconfounder. It is the harmonization gain

**Status:** ✅ Executed. **`src/` untouched, no label moved.** It refuted my own attribution.

### What I claimed, and what the measurement says

1.5.6 §4.2 said *"about half of HFF's ΔAge magnitude is contributed by pipeline processing"*.
**False.** Decomposing the chain step by step:

| step | day-0 | day-6 | day-14 | contributes at day 14 |
|---|---|---|---|---|
| S1 clock, absolute age | 78.65 | 78.30 | 68.04 | — |
| S2 control-relative (ΔAge) | −0.00 | −0.35 | **−10.62** | — |
| S3 cell-cycle deconfounded | −0.32 | +0.04 | −8.83 | **+1.79** |
| S4 re-centred = `y_age` | 0.00 | +0.36 | **−8.51** | +0.32 |

**Deconfounding and re-centring contribute +2.11 yr, and they move ΔAge TOWARD zero.** They cannot
explain a −24.02 shard value; the gap is still **15.51 yr**.

### The real source — and I had ruled it out on bad evidence

**Harmonization was ON in the actual build.** I checked `configs/data/*.yaml`, found no
`harmonize: true`, and concluded it was off. The build is driven by the runner:
`local_runners/run_multi_local.py:161` — `harmonize=HARMONIZE, harmonize_ref_dataset="gill_bulk"`.

**And Stage 1.5 §2 Group B derived the mechanism nine days ago:**

```
ΔAge = Σ_g (x_pert,g − x_ctrl,g) · sigma_ref,g / (sigma_d,g + EPS) · w_g
```

`sigma_d` does **not** cancel — it survives as a **per-dataset multiplicative gain**, and HFF carries
`sigma_gill / sigma_hff`. That is exactly why Group B recorded "batch-immune by construction" as an
overstatement: ΔAge is immune to *additive* batch effects, **not to scale ones**.

**Implied gain: 24.02 / 10.62 ≈ 2.26** — the shape Group B predicts.

**Attribution by elimination plus a matching closed form, NOT a measurement.** The confirming test is
one number and is now step **1c**: compute `Σ|w_g| · sigma_gill,g / sigma_hff,g` over the clock's
genes and check it lands near 2.26.

### Why this changes the sparse-clock plan

| | size | status |
|---|---|---|
| clock density bias (Gill, ΔAge) | **−14.10 yr** | ✅ measured, LODO-validated |
| deconfound + re-centre | **+2.11 yr** | ✅ measured — small, wrong direction |
| **harmonization gain** | **≈ ×2.26** | ⚠️ attributed, not measured |

**A gain and a bias compose differently.** Sparsifying changes the weighted sum; harmonization then
*multiplies* it. So **the sparse clock must be evaluated with harmonization ON** — which §1's
Gill-side work did not do, because **Gill is the reference dataset and its own gain is ≈ 1, making
the effect invisible there.** That is a genuine limitation of §1's result as applied to HFF, and it
is why step 4's rebuild cannot be skipped.

---


## 2026-08-04 — Plans 1 → 1.5.5 CLOSED, and a self-audit of 1.5.6 found three errors in it

**Status:** ✅ Closures applied; 1.5.6 corrected. **835 tests pass, ruff clean, `src/` untouched.**

### Closures — two were stale in ways that mattered

| plan | was | now |
|---|---|---|
| **Stage 1** | carried as **PARTIAL / open**, "Change A″ awaiting a third LOOCV run" | ✅ **CLOSED.** A″ **is** in `src/` (`platt_safe`, `fit_platt_binary`) and **was** scored — `scorecard/B_fatecal.json`, 2026-07-22 22:08: `conformal_coverage` 0.4006 → **0.8889**, `fate_ece` 0.2806 → **0.2491**, `fate_ece_platt` 0.1535 → **0.1399**. The regression that made it PARTIAL came from `A_xdonor` = **run 1, recorded INVALID** (HFF rotated as a donor) |
| **Stage 1.5** | EXECUTED, three findings | ✅ **CLOSED.** D1 measured and downgraded; D2 measured **INDETERMINATE** (56 %, CI [9 %, 100 %]) and data-blocked with a sized owner, its **code half shipped as G-a**; D3 answered (M1 FAIL) with the **wiring shipped as G-b** |
| **Stage 1.5.1** | §11 gave three items owners | ✅ **CLOSED.** Two are data-blocked and sized; the third (**HFF's `age_mask`**) was **delivered** by 1.5.3 step 6 |
| **Stage 1.5.2** | already closed | ✅ unchanged |
| **Stage 1.5.3** | **"steps 5–7 open"** and step 6 marked **"🛑 BLOCKED — STOPPED before running"** | ✅ **CLOSED.** Step 6 **ran three times** (confounded, clean rerun, arm C) — `results/STEP6_REPORT.md` and three `scorecard/gc2_*.json` snapshots prove it. The table had been wrong for two days |
| **Stage 1.5.4 / 1.5.5** | executed | ✅ unchanged |

### The 1.5.6 self-audit — I checked my own headline before letting it stand

**The core result survives the hardest test available.** A MAE gain is fake if the predictor merely
shrinks toward zero. It does not:

* `top100` SD **11.82** vs truth's **11.41** — ratio **1.04**, spread *matched*, not collapsed
* raw SD **22.99** — the dense clock **over-disperses 2.02×** on top of its −14 yr bias
* **a constant-zero predictor scores MAE 8.45. `top100` scores 5.36. The full clock scores 16.61 —
  worse than predicting nothing at all.**

That last line is *stronger* than what §1 originally claimed.

**Three errors found and corrected:**

1. **§1.3 "the Sendai arm agrees independently" — WITHDRAWN.** Sendai is scored on *absolute* age
   where the intercept does not cancel. Sparsifying drags every prediction toward **b0 = 72.43**:
   raw mean 98.86 → top100 **67.71** (4.73 from the intercept), while truth sits at **28.29**, and
   top100's SD collapses to **6.50** against truth's 14.75. It moved toward the intercept, not toward
   truth. **Not corroboration.**
2. **"below methylation's own ±7 yr error" — WITHDRAWN.** The ±7 is a **donor-level absolute-age**
   error; MAE 5.36 is a **condition-level ΔAge** error. Not like-for-like.
3. **Skin & blood was understated.** §3 said ordering "stays poor". Measured: `top100` scores MAE
   **8.79** against a constant-zero floor of **6.83** — **it loses to predicting nothing.** The spread
   is right and the ordering is wrong, which is worse than abstaining.

**Net: the finding is real but confined** — multi-tissue ΔAge, transient arm, one estimand. Not the
two-arm agreement the document originally claimed.

---


## 2026-08-04 — STAGE 1.5.6: the clock's DENSITY is the defect. Sparsifying it cuts ΔAge error 3x

**Status:** ✅ Measured and validated leave-one-donor-out. **`src/` untouched, no label moved.**
Falsifier (HFF) running; §5 step 1 of the plan.

### The finding

The Fleischer clock is a dense RidgeCV over **33,155 genes fitted from 133 samples**. Restricting it
to its **~100 largest-|weight| genes**, changing nothing else:

| | MAE vs methylation | bias | ρ | sign agreement |
|---|---|---|---|---|
| **full clock (33,155 genes)** | **16.61 yr** | **−14.10** | +0.703 | 0.62 |
| **top-100 genes** | **5.36 yr** | **−1.61** | **+0.835** | **0.94** |

*(68 conditions, Gill transient arm, Horvath multi-tissue as truth, ΔAge vs ΔAge, replicates averaged)*

**MAE 5.36 yr is below the reference instrument's own donor-level error of ±7 yr.** The RNA clock now
agrees with methylation about as closely as methylation agrees with itself.

### The mechanism, not a tuned number

| k | 20 | 50 | **100** | 150 | 300 | 1000 | all |
|---|---|---|---|---|---|---|---|
| MAE | 6.19 | 5.99 | **5.36** | 5.42 | 9.99 | 12.33 | **16.61** |
| bias | +5.06 | +3.22 | **−1.61** | −1.99 | −7.69 | −10.15 | **−14.10** |

**The bias crosses zero exactly where MAE bottoms out.** Thousands of near-zero weights each
contribute a little drift; summed over 33,155 genes they become a **−14 yr systematic offset**.
Dropping them removes the offset rather than shrinking noise.

**Leave-one-donor-out:** k chosen on two donors, scored on the third — held-out MAE **6.70 / 6.84 /
5.45** against the full clock's **16.59 / 16.48 / 16.75**. The full clock is worse than every sparse
variant for **every donor individually**, so this is not selection. The Sendai arm independently puts
its MAE minimum at k=50 (28.52 vs 65.63) — different protocol, partly different donors, same region.

### Why nine months of correlation tests never saw it

**Spearman is shift- and scale-free, so a uniform −14 yr offset is invisible to it.** M-2a, 1.5.4 and
this stage's own first sweep all scored ρ_partial and all reported "weak but present". MAE and bias
against a real instrument are what exposed it. The lesson generalises: **a correlation is the wrong
summary for a quantity whose units are the claim.**

### What is NOT fixed

**Horvath skin & blood does not come right at any k.** MAE improves (17.84 → 6.69) but ordering stays
poor throughout (ρ ≤ 0.43, sign 0.41–0.68, sometimes below chance). The sb/mt asymmetry from 1.5.2
and 1.5.4 **survives**, so this does not clear M-2a's SPLIT rule. What changed is knowing *what kind*
of failure it is: biased on one axis, mis-ordered on the other.

### Two errors of my own, both found and corrected here

1. **Pseudo-replication.** The first two-arm run scored 30 Sendai rows from 22 independent
   methylation samples and 90 transient rows from 68 conditions — exp1/exp2 left unaveraged. Fixed:
   every modality is now collapsed to one row per condition **before** scoring.
2. **The wrong summary statistic.** I ran an entire 9-variant sweep on ρ_partial and reported "no
   variant helps". That sweep was blind to the thing that mattered by construction.

### Deliverables

`results/DAGE_LEDGER.md` — 68 per-condition rows, each with **TRUTH** (methylation ΔAge),
**EXPECTED** (what the pipeline should produce), **ACTUAL** (what it produced) and **ERROR** in
years. Plus `results/dage_ledger.csv` (90 × 60) and the full k-sweep JSON. Plan in
`plans/STAGE_1_5_6_SPARSE_CLOCK.md`.

### Why this plausibly reaches HFF — and the falsifier now running

**The −14 yr bias is a property of the clock's weights, not of any dataset**, so it is applied to
HFF's 33,613 labels too. Pre-registered prediction: HFF's day-14 ΔAge should move from the recorded
**−24.0 yr** to roughly **−10**, with the trajectory shape surviving (ρ ≈ −0.9). **If the shape
collapses instead, the hypothesis is wrong and this stage stops.**

### ⚠️ The falsifier ran, and it VOIDED its own prediction — plus found something larger

Predicted: HFF's day-14 ΔAge moves −24.0 → ≈ −10 under the sparse clock, shape preserved.

| k | ρ(day, ΔAge) | day-14 |
|---|---|---|
| top50 | −0.905 | −8.36 |
| **top100** | −0.881 | **−10.72** |
| **all 33,155** | −0.857 | **−10.62** |

**The shape survived** (the actual falsification condition). **But the predicted shift did not
happen, because the full clock in this run already reads −10.62, not −24.0.**

**My script reported CONFIRMED. That was wrong**, and it is the same error class I flagged in others
this week: it checked the result against the prediction **without checking that the baseline the
prediction rested on still held**. Corrected to `BASELINE_NOT_REPRODUCED`.

### 🔑 The 13.4-year gap — the real finding

| day | pipeline `y_age` (built shards) | clock applied directly |
|---|---|---|
| 2 | **+3.85** | −0.52 |
| 6 | −5.80 | −0.35 |
| 12 | −8.23 | −1.26 |
| **14** | **−24.02** | **−10.62** |

G-c step 1 read **built shards** — the pipeline's `y_age`, which adds harmonization (Gill
Projection), cell-cycle deconfounding and control re-centring. This run applies the clock directly.

> **About half of HFF's apparent ΔAge magnitude comes from pipeline processing, not from the clock.**
> Nobody has audited that 13.4 yr.

A sparse clock addresses a −14 yr bias **in the clock**; it cannot address a +13.4 yr contribution
from **downstream processing**. **Both are the same size, and only one has been measured** — so
auditing the pipeline's contribution is now step 1b and outranks writing the sparse-clock config.

Also unexplained: the direct route reads **+3.58 yr at day 4** — cells four days into reprogramming
reading three and a half years OLDER, in both routes. Same class of impossibility as the +36.5 yr
non-responder reading, surviving in HFF.

---


## 2026-08-03 — Stage 6 rewritten against the data actually held, and the ask is sized: 2 donors, not 15

**Status:** ✅ Written. `plans/STAGE_6_NEW_DATA_REV.md`. The original is **byte-unmodified**; its
species blocker and `input_ablation` gate are carried forward unchanged because both still hold.

The original names **2 datasets**. Six are held, and four of the questions it was written to answer
have since been closed by measurement — executed as written it would buy data to answer questions
the data already answers.

### The sizing nobody had done

Every unresolved statistic here is **donor**-limited. Derived from measured instrument errors:

| open question | measured spread | donors needed |
|---|---|---|
| step 6's arm comparison, σ known | SD 4.808 | 10 |
| step 6's arm comparison, **σ-robust** | SD 4.808 | **19** |
| contrast B retention, multi-tissue | SD 11.6 | 17 |
| contrast B retention, skin & blood | SD 18.2 | **38** |
| M3 — offset real or `n=1` noise? | σ CI factor 3.93 at n=6 | ~20 |

σ-robust is the bar the step-6 rerun correctly **refused** to claim at n=4, where the χ² interval put
the MDE at 6.704 against Δ\* 3.572.

### 🔑 The finding that changes what to buy

Dropping the two donors outside the clock's fitted range collapses the per-fold difference SD
**4.808 → 1.130** — a factor of 4.3, on 2 of 6 folds. Re-running the sizing on that:

| | σ-robust donors needed |
|---|---|
| all donors | **19** |
| **in-range only** | **6** |

**Four in-range donors are already held. The gap is 2, not 15.**

**Fixing the instrument's domain is worth more than tripling the donor count.** C-2
(`enforce_clock_age_range`) is already built and shipped inert; this is the arithmetic that justifies
turning it on. That reframes the acquisition from "a cohort" to "two donors with paired methylation".

### One section deliberately left open

**§6 — what to acquire — is GATED on the stratified shuffle.** If HFF's residual structure is
systematic artefact, HFF must be *replaced* and the target is volume. If it is real signal, HFF is an
asset and the target is the 2 donors §4 sizes. Writing a target now would be guessing at the one
number this stage exists to get right.

### Ordering, and why it is not negotiable

1. turn on **C-2** — free, code already written, worth more than 13 donors
2. integrate **GSE165177** — 95 adult in-range methylation-paired samples, in no training config, free
3. the **stratified shuffle** — settles §6's gate, free
4. **then** acquire

Steps 1–3 are free and must precede 4, because each changes what 4 should buy.

Also recorded, because a plan that only lists what is buyable is dishonest: **no acquisition fixes
HFF** (one cell line — settled by methylation *on HFF*, not cohort size), **the RNA route stays
closed**, and **contrast B on skin & blood needs 38 donors**, which is likely out of reach and is
stated as a limit rather than carried as a plan.

---


## 2026-08-03 — STAGE 1.5.5: HFF's ΔAge is NOT an identity readout, and not a depth readout either

**Status:** ✅ Executed and recorded. `src/` untouched, **no label moved**, CI lint clean.

Arm C established that HFF's labels carry *"consistent exploitable structure of unknown provenance"*
and could not say what it was, because shuffling destroys real signal and systematic artefact alike.
This measures the provenance directly.

**The question:** within a single timepoint, how much of a cell's ΔAge is explained by identity?
Between timepoints ΔAge and pluripotency move together trivially, so a pooled correlation is
uninformative **by construction**. Holding day constant is the whole design.

**Why it could resolve when nothing else has:** every dead end in this project is *donor*-limited —
3 donors, 6 folds, MDE 1.049 × SD, M3's CI of [9 %, 100 %]. HFF is the one place that is not:
**44,473 cells, ~5,800 per timepoint, SE(ρ) ≈ 0.013.** Two orders of magnitude more power than
anything else in the data.

### Verdict: NOT_DOMINATED — 0 of 8 timepoints reach the pre-registered 0.50 bar

| day | ρ(age, pluri) | ρ(age, somatic) | R² identity | R² technical | **R² both** |
|---|---:|---:|---:|---:|---:|
| 0 | +0.079 | +0.171 | 0.033 | 0.086 | 0.100 |
| 2 | −0.079 | +0.265 | 0.082 | 0.075 | 0.123 |
| 4 | +0.036 | +0.111 | 0.034 | 0.030 | 0.054 |
| 6 | −0.158 | +0.281 | 0.081 | 0.008 | 0.083 |
| 8 | −0.130 | +0.133 | 0.020 | 0.007 | 0.026 |
| 10 | −0.128 | +0.251 | 0.076 | 0.030 | 0.114 |
| 12 | −0.270 | +0.397 | 0.162 | 0.013 | 0.167 |
| 14 | −0.194 | +0.256 | 0.068 | 0.010 | 0.112 |

**Identity and technical covariates together explain 2.6–16.7 % — so 83–97 % is neither.**

**The design's own confound check, in one number:** pooled ρ = **−0.216** against within-timepoint
values near zero. Pooling would have manufactured exactly the artefact the test exists to detect,
which is why it is reported as descriptive and never graded.

### What is now measured rather than argued

| candidate explanation for HFF's ΔAge | status |
|---|---|
| identity artefact (+36.5 yr mechanism, `corr(age,pluri) = −0.62`) | ❌ **rejected**, R² ≤ 0.16 |
| technical / sequencing depth | ❌ **rejected**, R² ≤ 0.09 |
| clock noise | ⚠️ still live |
| real biological signal | ⚠️ still live |

**This does not establish the labels are age.** It is the first positive characterisation they have
ever had, and it removes the two cheapest ways to dismiss them.

### Limitations that bound the claim, recorded because they cut against it

1. **Identity is proxied by 22 genes** (18 pluripotency + 4 somatic). The strong form of the identity
   hypothesis is "position on the reprogramming manifold", which a small marker set under-measures.
   A richer axis could explain materially more.
2. **R² is linear on ranks** — monotone-robust, but blind to a non-monotone relationship.

Both push the same direction: **the measured identity share is a floor, not a ceiling.** The honest
reading is "identity does not dominate", not "identity is absent".

### Two implementation notes

* **The first run did not fail slowly, it thrashed.** `cells_per_run=None` densifies every cell into
  one chunk — ~48,000 × 36,601 × 4 bytes ≈ **7 GB**. The source's own docstring says that parameter
  "bounds peak RAM"; set to 4,000. Batching cannot change a result here because every statistic is
  computed per timepoint *after* the chunks are concatenated.
* **I reintroduced `N802`** — capitals in a test function name — the exact defect that held CI red
  from 2026-07-26 and that I diagnosed and fixed two days ago. Caught by ruff before commit this
  time, which is the argument for the lint gate being where it is.

---


## 2026-08-03 — STAGE 1.5.4 EXECUTED: NOT LEARNABLE. The last free route to repairing ΔAge is closed

**Status:** ✅ Executed and recorded. **787 tests pass**, CI lint clean, **`src/` untouched**, no
label moved (none could, under any outcome).

M-2a asked whether the **existing Fleischer clock's output** tracks methylation age → no. That is a
fact about one clock, not about the transcriptome. 1.5.4 asked the question nobody had: **can a model
TRAINED on the transcriptome predict methylation age on these cells?**

**Verdict: NOT LEARNABLE — SPLIT on all three model families.**

| ρ_partial (pluripotency **and** donor removed) | skin & blood | multi-tissue | |
|---|---:|---:|---|
| **Fleischer clock** *(baseline, identical estimand)* | **0.309** | **0.517** | — |
| learned — full transcriptome | 0.277 | **0.627** | SPLIT |
| learned — clock genes only | 0.247 | **0.604** | SPLIT |
| learned — PCA(10) → ridge | 0.386 | **0.579** | SPLIT |

All six beat their own permutation null, so the models **do** learn something real from RNA. It is
not age.

### The finding: the asymmetry survives retraining

Training improves agreement with multi-tissue (+0.06 to +0.11 over Fleischer) and **does not**
improve agreement with skin & blood (two of three families are *worse* than the untrained clock).
The ratio widens rather than closing — 1.67× for Fleischer, **2.26×** for the full-transcriptome
model.

**A measurement of the shared age signal cannot behave that way.** The two methylation clocks agree
with each other, so anything tracking what they share must track both at similar strength.
**The asymmetry is a property of the DATA, not of the Fleischer clock** — it was never "we had the
wrong clock." RNA in these cells does not carry the shared age signal, and retraining cannot
manufacture it.

### Three method corrections, ALL made before any real ρ was read

G2 exists so the method can be repaired without touching the result. It fired three times and was
right each time.

1. **The first run was VOID.** Shuffled training labels gave ρ_partial **−0.45 / −0.36**, not ≈ 0.
   Cause: a LODO fold predicts the held-out donor with roughly the mean of the *others*, so the
   prediction is **anti-correlated with the donor mean by construction** — a large correlation from
   a model that learned nothing. **Fix:** partial out donor as well as pluripotency. Reproduced
   synthetically in `test_lodo_mean_reversion_is_removed_by_partialling_donor`.
2. **My own G2 bar was UNRESOLVABLE — the §5b failure inside my own stage.** §5 registered
   `|ρ| ≤ 0.20` on the worst of 6 comparisons with "slack for sampling noise" **asserted and never
   simulated**. At n = 68, SD(ρ) ≈ 0.122 ⇒ P(|ρ| > 0.20) ≈ 0.102 per comparison and **≈ 0.474 across
   six**: it fires on a sound pipeline almost half the time. **Fix:** replaced with a self-calibrating
   **permutation null** (20 draws per family × clock; the real value must beat its own q95).
3. **A PCA bug** — components fitted on centred data, applied to uncentred matrices. Visible in G2
   as a `pca`-only residual artefact.

**The bar that decides the verdict — ρ_partial ≥ 0.50 on both clocks — was not touched.** It is
M-2a's registered bar, re-used deliberately so no friendlier one could be chosen for a stage whose
result was predicted negative in advance.

### One correction to my own reporting

The script initially compared against M-2a's published 0.267 / 0.516, which was **not**
apples-to-apples: 1.5.4's estimand partials out donor and M-2a's does not. Recomputed the Fleischer
baseline under the identical estimand (**0.309 / 0.517**) before reading anything into the gap.

### What it licenses

**Does:** closing the RNA route on **measurement** rather than inference. Stage 6's spend is now
justified — the cheap alternative was tried and it failed.

**Does not:** any label change; any claim about HFF (neonatal, unfixable by better instrumentation);
anything beyond 3 donors.

---


## 2026-08-02 — The order of work from here, written into `00_START_HERE.md`

**Status:** ✅ Documentation only. Additive section; the 2026-07-22 status table is left byte-intact.

Asked to write the plan and running order into a file. **Put it in `00_START_HERE.md` rather than a
new document**, because that file's own opening line is *"The running order for the whole project"* —
a second document claiming to be the order is exactly the drift this project keeps having to
correct. Its "Where you are right now" table was last touched **2026-07-22**, so it predates the
entire Stage 1.5 → 1.5.3 arc; that table is preserved as the record of its moment and the new
section sits beside it.

### What the section records

**Where the project actually stands**, with every figure re-verified from `scorecard/baseline.json`
before it was written down: fate head **works** (`fate_roc` 0.983, `fate_prauc` 0.992); ΔAge
*relative* consistent but **not externally validated** (`rank_model_dage` 0.948 is the model
reproducing its own labels); ΔAge *absolute* **broken** (14.29 yr against a 12.27 yr instrument,
**SNR ≈ 0.90**; coverage 0.401 vs nominal 0.9); and the product **dead** — `res_median` **0.000 on
every fold**, because `R_eff` needs `mu_age < −30` against a measured ~−11.

**The order:** (1) finish 1.5.3 by running step 6; (2) **Stage 1.5.4** — can a model *learn* age
from RNA, the question M-2a never asked; (3) integrate **GSE165177**; (4) rewrite Stage 6; (5)
acquire. **Steps 1–3 are free and none of them fixes ΔAge at scale — only #5 does.** 1–4 exist so
that #5 buys the right thing.

**The number nobody has computed:** how many **donors** the age arm needs. Every unresolved statistic
here is donor-limited — M3's *[9 %, 100 %]* CI, contrast B's ≈16-vs-9 pairs, step 6's MDE of
`1.049 × SD` at 6 folds. **GSE165177 triples the labels and adds one donor**, which is precisely the
trap to avoid when sizing Stage 6.

**Why HFF may not be fixable at all:** a neonatal line has almost no chronological age to remove, so
even a perfect instrument reads its true ΔAge as ≈ 0. Better measurement on HFF buys an accurate
zero. That is why the age arm needs **adult-donor** data — and why "B better" is the *expected*
outcome of step 6 rather than a surprise.

### One deliberate omission

**Stage 1.5.4 and the Stage 6 rewrite are named as work, not linked as documents**, and no row was
added for them to the file table. `STAGE_6_NEW_DATA.md` once carried an acceptance gate naming a test
that had never been written — a gate that could therefore never fail. Listing a plan file before it
exists is the same mistake, and the section says so explicitly.

---


## 2026-08-02 — Fix the late-crash encoding bug: scripts finished their work, then died before writing it

**Status:** ✅ Fixed in 4 forward-path scripts and verified by re-running the one that failed.
**756 passed**, CI lint clean.

Found while checking whether step 6 could run locally. `experiments/diag_m2a_calibratability.py`
computed everything, printed its SPLIT verdict, and then **raised `UnicodeEncodeError` on a `Δ`
character under codepage `cp1255`** — *before* writing its results JSON. The recorded result
survived only because the crash beat the write.

**The failure mode is the dangerous kind: late.** All the compute succeeds, the verdict reaches the
screen, and the artefact is silently never persisted. On a long diagnostic that is an entire run
thrown away, and nothing in the output says so.

**22 scripts print characters this codepage cannot encode.** Fixed the four on the forward path:

| script | why it matters |
|---|---|
| `plan_tests/verify_age_mask_identical.py` | the step-1 **bit-identity gate** |
| `scorecard.py` | what **step 6** compares its arms with |
| `experiments/diag_m2a_calibratability.py` | the one that actually lost its output |
| `experiments/diag_gc_hff_signature.py` | current G-c step 1 diagnostic |

`plan_tests/verify_stage1_5.py` **already carried this guard** (line 270) with the same cp1255
rationale, so the convention was copied rather than invented. Extended to `stderr` as well, because
a traceback whose source line contains one of those characters fails the same way — which is exactly
what happened to the detector script written to find this.

The remaining 18 are historical one-off experiment scripts (`test3_linearity.py`, `test7_*`,
`test18_forward_gate.py`, …) that are already run and recorded and sit on no forward path. Left alone
deliberately rather than swept up.

### Verified by re-running the failure

`diag_m2a_calibratability.py` now exits **0**, renders `Δ` and `§` correctly, and **writes its
JSON** (timestamp moved from 2026-07-31T13:18:41 to 2026-08-02T06:36:54).

**And it reproduced bit-for-bit on a different machine:** every `rho_all` / `rho_partial` identical
to the data-machine run (0.4445 / 0.2671 skin & blood, 0.6902 / 0.5163 multi-tissue), verdict SPLIT.
So `results/diag_m2a_calibratability_results.json` changes only in its timestamp in this commit —
the numbers are unchanged, and the re-run is now an independent cross-machine reproduction of M-2a
rather than a single recorded result.

*(`scorecard.py` carries 2 pre-existing ruff errors at lines 232 and 266 — `UP017`, `E702`. Not
introduced here, not on CI's lint path, left alone.)*

---


## 2026-08-02 — Step 6's bar REGISTERED, and `k = 4` pinned. Suite verified locally at 756 passing

**Status:** ✅ Both gaps from the deep review are closed. **Full suite run locally for the first
time this session — 756 passed**, CI lint clean.

### GAP 1 closed — the arm comparison now has a bar

Step 6 decides whether **99.7 % of the project's age labels** are kept or discarded, and its
criterion had no registered bar. Every bar for this stage graded a *mechanism* (B1/B2, A1/A2/A3).
Ground rule §5b: *"a bar with no such test is not considered pre-registered."*

**`plan_tests/register_gc_step2_bar.py`** → `results/register_gc_step2_bar_results.json`, plus **3
rows** in `tests/test_bars_resolvable.py`.

**Δ\* = 3.57 yr, derived not chosen:** Stage 2 §12 already registers *"≥ 25 % drop in
`dage_mae_model`"* as its TARGET; applied to the 14.29 yr recorded baseline. Using the project's own
existing threshold for the same metric avoids inventing one.

| SD(per-fold difference) | MDE | P(detect Δ\*) | verdict |
|---|---|---|---|
| 0.5 | 0.52 | 1.0000 | RESOLVABLE |
| **1.0** | **1.05** | **1.0000** | **RESOLVABLE** |
| 2.0 | 2.10 | 0.9338 | UNRESOLVABLE |
| 3.0 | 3.15 | 0.6476 | UNRESOLVABLE |
| 5.0 | 5.25 | 0.2955 | UNRESOLVABLE |
| 13.7 *(arms independent)* | 14.38 | 0.0752 | almost pure noise |

False-positive rate at a true effect of 0: **0.0508** — the CI is honest, it is only weak.

**Δ\* is detectable at ≥ 95 % only if the arms track each other to within ~1 yr per fold**, on a
metric whose baseline already ranges 5.39 → 29.69 across folds. That is demanding and **is not known
to hold.**

**The reading rule is pre-registered because the NULL is the dangerous outcome.** "B better" is
self-limiting. *"CI includes 0"* read as "HFF's labels contribute nothing, discard them" would throw
away 99.7 % of the labels **on a null that may simply be underpowered.** So: CI includes 0 **with
MDE > Δ\*** is **INCONCLUSIVE and licenses nothing** — explicitly not a licence to discard. The run
must report its observed SD and MDE beside the effect, because which row of the outcome table
applies depends on the MDE.

### GAP 2 closed — `age_window_k = 4` pinned into step 6

5c ships inert at `k = 1`, and 1 means OFF. Step 6's gate said only *"5c must have shipped"*, and no
command set `k`. **Run as written, both arms would have used k = 1, arm B would be starved, and
problem #1 from the readiness audit would return silently** — the confound 5c exists to remove,
reintroduced by a default. `k = 4` is B2's registered value, not a new choice.

### Found by running the tests rather than reading them

The new script defined `_RESULTS` as `ROOT / "results"`. `test_results_paths.py` checks that form
**by regex** and cannot follow the indirection, so it failed. Spelled out to
`Path(__file__).resolve().parents[1] / "results"`. Worth recording: the convention test earned its
keep on the first new file written after it.

### Verification capability improved

`numpy`, `scipy`, `pyarrow`, `pandas`, `torch` (CPU) and the `[dev]` extra are now installed here,
so the suite runs locally. Previous entries in this file carried an explicit
*"suite not re-run, asserted not verified"* caveat; that no longer applies from this entry onward.

*(2 ruff errors remain in `experiments/` — `F541`, `F401`. Not in CI's lint path and out of scope.)*

---


## 2026-08-02 — Deep review of steps 1-5: work is sound, but step 6 is NOT ready (two gaps)

**Status:** ✅ 3 fixes applied and pushed. **2 gaps found in step 6's readiness — not yet fixed.**

Reviewed the other machine's steps 1-5 across logic, bars, coding and flow, recomputing the
load-bearing numbers rather than trusting them.

### Verified correct (checked, not assumed)

* **The core guarantee holds.** `verify_age_mask_identical` -> **IDENTICAL, max_abs_delta 0.0** over
  7 chunks -- and it carries a **self-test** proving it can detect a 1-ULP change, a mask flip and a
  reason appearing. A gate whose only exercised path says PASS is not a gate; this one is not that.
* **Every switch ships inert:** `AGE_MASKED_DATASETS` empty, `enforce_clock_age_range=False`,
  `age_window_k=1`.
* **C-5's bar arithmetic reproduces exactly.** p = 75/33688 = 0.002226 -> 1.140 cells/batch (claimed
  1.14); P(empty) = 0.3195 exact vs 0.3199 Poisson (claimed ~32%); B1 = 0.6805 (claimed 68.9%);
  B2 status quo = 0.0286 (claimed 2.9%). `k = 4` halves the per-update SE **exactly** (1/sqrt 4).
* **The bar DISCRIMINATES and the script exits non-zero if it stops doing so** -- stronger than §5b
  requires.
* **The plan's own recommendation was overturned by measurement** (Option 1 -> Option 2) and
  withdrawn openly rather than quietly replaced.
* **Mutation-tested guard:** they re-injected the exact fixed-W bug the readiness audit found,
  confirmed the guard fails, then restored.
* `huber_age_window` reduces once over concatenated cells, **not** a mean of batch means.
* `_AgeWindow` **re-forwards** buffered cells rather than storing stale activations.
* C-4's `age_ok` defaults to **False** when provenance is absent -- the conservative direction.
* Cleanup deleted **no `.py` or `.md`** -- only zips, a cache and a notebook.
* All three of my earlier fixes survived, including the fail-open `raise`.

### Fixed in this commit

1. **`zip(..., strict=True)`** in `experiments/verify_rev_final_4_4.py` (2 sites, mine). An
   unchecked `zip` truncates silently -- the same shape as the census collision.
2. **That script wrote to the repo root**, which the tidy-up had moved to `results/`. Repointed both
   its read and its write to `_RESULTS`, and cleared its 7 `E701`s with a lookup table. **All four
   checks still reproduce byte-for-byte** (V1 -24.05/-27.55, V2 -1.13/-3.62, V3 rho -0.885/-0.842).
3. **A hole in `tests/test_results_paths.py`.** Its `_RESULTS` check began
   `if "_RESULTS" not in t: pytest.skip("reads results but does not write any")` -- an **assumption,
   not a check**. The script in (2) mentioned a `*_results.json`, defined no `_RESULTS`, and wrote to
   root; it was skipped under a message asserting it did not write. The next run would have dropped a
   stray JSON into the root and turned `test_no_results_json_is_left_in_the_repo_root` red -- a
   latent CI failure the file existed to prevent. The skip is now conditional on the script
   containing no write call. Simulated over all 23 writers: none trips it.

### 🔴 GAP 1 — step 6's decision has NO registered bar

Registered bars: B1/B2 (C-5), A1/A2/A3 (C-5c), and Stage 1.5.2's. **All grade mechanisms.** The
comparison that step 6 actually decides on -- arm A vs arm B on `dage_mae_model`, paired across 6
donor folds -- has **no `bar_verdict` row and no resolvability check**.

Ground rule §5b: *"a bar with no such test is not considered pre-registered."* By the project's own
standard, **step 6's criterion is not pre-registered.**

It matters here more than usual. `sensitivity_multiplier(6)` gives **MDE = 1.050 x SD(per-fold
difference)**, and that SD has never been measured. Baseline `dage_mae_model` already ranges
5.39 -> 29.69 across folds (SD 9.67). If the paired difference is anywhere near as heterogeneous, a
real effect would read as noise -- **exactly the §5b failure that bit Stage 1 twice on `fate_ece`
and that Stage 1.5.2 caught for M-2a.** This is the step that decides whether 99.7% of the age
labels are discarded; it should not be the one bar nobody checked.

### 🔴 GAP 2 — step 6 never pins `age_window_k = 4`

5c ships inert at `age_window_k = 1`, and 1 means OFF. Step 6's gate row says only *"5c must have
shipped"*. **Nothing instructs the operator to set `k = 4` in both arms**, and no command in PART E
sets it. Run as written, both arms use k = 1, arm B is starved again, and **problem #1 from the
readiness audit returns silently** -- the confound 5c was created to remove. The value is derived
(B2's registered `k`), so this is a one-line pin, not a decision.

---


## 2026-08-01 — CI GREEN on 3.11 and 3.12. This closes the "not verified" caveat on three entries

**Status:** ✅ Verified by execution, not by inspection.

The first CI run after the lint fix (`84800fc`) is **green on both Python 3.11 and 3.12**. Because
`ruff` had been aborting the job before `pytest` since **2026-07-26**, this is the first time the
suite has actually executed in CI in roughly a week -- so green establishes considerably more than a
routine pass:

| what had never been CI-verified until now | status |
|---|---|
| the whole Stage 1.5.2 / 1.5.3 arc from the other machine (~19 commits) | ✅ passes |
| the `src/` changes for gates **G-a** and **G-b** (`aging.py`, `sources.py`, `build_dataset.py`, `data/__init__.py`) | ✅ passes |
| **`test_delta_age_is_bit_identical_with_and_without_the_census`** -- the hard guard that G-a records and does not compute | ✅ passes |
| the census key-collision fix (`chunk_id::line`) | ✅ passes |
| the `verify_stage1_5.py` ragged-row fix | ✅ passes |
| the **three regression tests added in `1380cb2`**, which had never run anywhere | ✅ pass |
| the suite on **two interpreters**, not just the data machine's | ✅ passes |

**This retroactively closes an explicit caveat carried by the three preceding entries in this file**,
each of which recorded that the suite could not be run locally (no `numpy`/`torch`/`pytest` on this
machine) and was therefore *asserted-not-verified*. It is now verified, by execution, on a clean
machine, twice.

**What it still does not establish:** that any *label* is correct. CI proves the code does what its
tests say, including that ΔAge is bit-identical across the G-a change. It says nothing about whether
the ΔAge values themselves are right -- that is the question Stage 1.5.2 answered in the negative and
Stage 1.5.3 exists to act on.

---


## 2026-08-01 — CI red X diagnosed: it was the LINT step, failing since 2026-07-26

**Status:** ✅ Fixed and verified locally. Two renames, no logic touched.

"Tests 3.11 and 3.12" are the **Python version matrix** in `.github/workflows/ci.yml`, not test
names -- so the red X was the whole CI job failing on both interpreters, not two specific tests.

**The failure was `ruff`, not `pytest`.** CI runs `ruff check src/ tests/ scripts/` **before**
`pytest -q`, so the lint step was aborting the job and **pytest has not executed in CI since the
lint break was introduced.** Reproduced locally with the exact CI command:

```
N802 Function name `test_census_keys_must_survive_one_cell_line_spanning_MANY_chunks` should be lowercase
N802 Function name `test_pairing_KEEPS_exp1_exp2_replicates_for_averaging` should be lowercase
Found 2 errors.
```

| offending name | introduced in | when |
|---|---|---|
| `test_pairing_KEEPS_...` | **`c7199d6`** | 2026-07-26 -- the original break |
| `test_census_keys_..._MANY_chunks` | **`1380cb2`** | 2026-08-01 -- **mine**, added to an already-red build |

Both used capitals inside a function name for emphasis. `pyproject.toml` deliberately selects the
`N` (pep8-naming) rules and ignores five specific codes, **each with a written justification**
(`N812`, `N818`, `N803`, `N806`, `N815`). **N802 is not among them**, so lowercase function names
are the project's intended standard.

**Fixed by renaming both to lowercase**, not by adding `N802` to the ignore list. Adding it would
have widened a deliberately narrow, individually justified exception list in order to keep a
stylistic flourish -- and "I wanted to shout in a function name" does not belong beside the reasons
already there. The emphasis was already carried by both docstrings, so nothing was lost.

Neither name is referenced in any `.md`, so the rename creates no cross-document drift.

**Verified:** `python -m ruff check src/ tests/ scripts/` -> **All checks passed!**

### ⚠️ What this does NOT establish

**Whether `pytest` passes.** It has been unreachable in CI behind the lint failure since
2026-07-26, and this machine has no `numpy`/`torch`/`pytest` to run it. Clearing the lint gate means
the test step will now execute for the first time in weeks, and **it may surface failures that were
simply never reached** -- including the three regression tests added in `1380cb2`, which have still
never run anywhere. If the X persists after this, the cause is a genuine test failure and the CI log
will finally name it.

*(`experiments/` carries 11 ruff errors, but CI does not lint that path, so they are not the cause
and are out of scope here.)*

---


## 2026-08-01 (addendum) — the 6th item: M-2b's bar was DERIVED, and now the proof is written down

**Status:** ✅ Applied. Documentation only; no code, no labels, no verdicts.

The review produced **six** items — 3 bugs and 3 attack points — and the previous commit actioned
**five**. The sixth (attack C, "M-2b passed by exactly zero margin") was judged to need no fix
because §14 already discloses it as `AGREE_FRAGILE`. **That judgement was half right.** The
*fragility* is disclosed. What was NOT written down is the answer to the sharper form of the
challenge: *"you loosened the bar from 8/11 to 7/11 and then landed exactly on it."*

Checked, and the answer is clean:

| | | source |
|---|---|---|
| resolvability simulation ran | **13:11:39** | `stage_1_5_2_resolvability_results.json` |
| M-2b ran | **13:53:13** — 42 minutes later | `diag_m2b_contrast_agreement_results.json` |
| registered bar 8/11 | **UNRESOLVABLE**, pass rate 0.9297 vs the 0.95 floor | resolvability |
| `usable_bar` **computed** by `audit_metrics` | **7.0** | resolvability |
| bar actually used | **7** — identical to the computed value | M-2b |

**The 7/11 bar is the output of §5b's `usable_bar`, computed from a simulated null before the data
was touched and frozen 42 minutes before the run** — not a number chosen to fit a result. §5b's
instruction on an unresolvable bar is *"move the threshold to `usable_bar` ... but do it now, not
after a run wears the failure"*, and the timestamps show that is the order it happened in.

This was fully present in the artefacts but split across two JSON files and never stated, so it
could not be checked without reconstructing it by hand. Recorded in §14 — the same class of gap as
the ceiling asymmetry: **the defence existed in the data and nobody had written it down.**

**Not defended away:** the result landed exactly on the bar, so one pair flipping changes the label.
That is why §14's conclusion rests on the **0/3 at the discriminating timepoint**, not on the 7/11.
The fragility is real; the goalpost-moving is not.

---


## 2026-08-01 — Line-by-line review of the Stage 1.5.2 / 1.5.3 work: 3 bugs fixed, 2 documentation gaps closed

**Status:** ✅ Applied. **No label moves, no verdict changes.** Two shipped-code bugs fixed with
regression tests; one bug fixed in code that has not been written yet; two documentation gaps closed.

Reviewed every change pulled from the other machine — `src/` diff line by line, and the load-bearing
numbers recomputed from the raw artefacts rather than taken on trust.

### What was verified as correct

* **`aging.py`'s census is genuinely additive.** The baseline arithmetic is textually unchanged, and
  `test_delta_age_is_bit_identical_with_and_without_the_census` uses `np.array_equal` (not
  `allclose`) *and* asserts the census was non-empty, so it cannot pass vacuously.
* **`sources.py`'s `_maybe_float` returns `None`, never `0.0`** — correct, because 0 is a real age
  here (N2/N3 are neonatal) and a silent default would be indistinguishable from it.
* **`verify_stage1_5.py` keeps the new G-a warnings OUT of `status`** — right call: four runs are
  recorded against the Stage 1.5 PASS, and folding a new condition into it would retroactively
  redefine what those PASSes meant.
* **M-2a's verdict logic** correctly demoted ρ_within and used ρ_partial — the fallback that was
  pre-registered in §6.
* **The resolvability work is the strongest part, and it corrected me.** My §6 fallback (ρ_partial at
  n=22) was **itself UNRESOLVABLE at 0.9233**. It was caught, and fixed by changing the *geometry*
  (n=68 via GSE165177) rather than the bar — §5b applied exactly as written.
* **M-2b** — registered bar 8/11 was unresolvable, moved to the §5b `usable_bar` of 7/11, landed at
  exactly 7/11, correctly marked FRAGILE, and the pass flagged as a day-axis artefact.
* **C-2's evidence recomputed from `scorecard/baseline.json` and matches exactly:** `dage_mae_model`
  3.01x, `dage_mae_ridge` 2.63x, `conformal_coverage` **0.000** on both out-of-range donors,
  `rank_model_dage` 0.910 vs 0.967.
* **GSE165178's join re-confirmed** at 22/22, 0 unmatched.
* **The changed `diag_methylation_anchor_results.json`** differs only in `utc` and `data_dir` — the
  re-run reproduced **identically**.

### Bug 1 (shipped) — the baseline census silently discarded 44 of HFF's 45 chunks

`build_dataset.py` merged each chunk's census with `baseline_census.update(chunk_census)`, keyed on
`cell_line`. **`cell_line` is not unique across chunks:** `verify_stage1_5_results.json` records
**HFF in 45 of them**. Every chunk overwrote the previous one, so the manifest kept **one** record —
for the dataset carrying ~99.8% of the age labels. A baseline problem in any chunk but the last was
invisible, which is precisely what gate G-a exists to prevent. It also undercut C-3, whose HFF
metadata fix is verified *through* this census.

**Fixed** by keying on `f"{chunk_id}::{line}"` and stamping `chunk_id` / `cell_line` into the record.
Demonstrated: three chunks of one line collapse to 1 record under the old key and keep all 3 (and all
4 warnings) under the new one. All 15 existing tests were single-chunk, so nothing caught it.

### Bug 2 (shipped) — `verify_stage1_5.py` crashed on any errored chunk

G-a widened the census table to six columns; the error branch still appended **five**. `render_table`
indexes `row[i] for i in range(len(headers))`, so a short row raises `IndexError` — crashing the
renderer on exactly the path `scan_build` goes out of its way to survive (*"recorded per chunk, never
aborts the scan"*). **Fixed**; header and both append sites now verified at 6 cells.

### Bug 3 (not yet written) — `age_label_policy` failed open

C-1's planned helper read `if masked_datasets and "dataset_id" in obs.columns:`. If a withholding
policy is switched on and the column it needs is absent, the silent outcome is to **keep labels that
were meant to be withheld** — the unsafe direction, and invisible. Not hypothetical: C-3 records that
G-b reached Gill and never reached HFF, so `donor_age` is missing on HFF today. **The spec now
raises**, with two tests added, and distinguishes a missing *column* (policy inapplicable — error)
from a missing *value* (recorded absence — never acted on).

### Gap 1 — C-2's range evidence does not explain Y2

The in-range mean coverage of 0.601 hides a split: O1 0.810, O2 0.667, Y1 0.737, **Y2 0.190**. Three
folds have broken coverage and the range criterion identifies **two**. Y2 is comfortably inside
`[1, 96]` and still fails. The claim as worded ("the two worst folds", "the range field is
informative") is correct, but it is **not a complete account of the calibration failure** and must
not be presented as one. Recorded in C-2.

### Gap 2 — the ceiling asymmetry argues FOR the verdict, and nobody had written it down

§12-R honestly recorded both readings of the ceiling, including the one that *qualifies* the verdict
(RNA↔multi-tissue reaches 91% of the meth↔meth ceiling). The inference it licenses was never drawn,
and it runs the other way.

**These are not two RNA clocks — they are ONE RNA clock against two methylation references**, and
those references agree with each other at **+0.568**:

| | ρ_partial |
|---|---|
| Horvath-mt ↔ Horvath-sb (the two references, to each other) | **+0.568** |
| Fleischer RNA ↔ Horvath-mt | +0.516 |
| Fleischer RNA ↔ Horvath-sb | **+0.267** |

Anything genuinely tracking the shared signal must correlate with **both** references at broadly
similar strength, since each reference's own reliability bounds how well any third measurement can
agree with it. A **2x asymmetry against references that agree with each other at 0.57** is the
signature of tracking something **clock-specific rather than age** — §1's diagnosis reached from an
independent direction. **So the qualifying reading does not survive**, and SPLIT understates the
result. Recorded as an inference from already-published numbers; nothing was computed after the fact.

### Not verified

The full suite was **not re-run** — this machine has no `numpy`/`pytest`. The two shipped fixes were
syntax-checked, the table widths verified programmatically, and `census_warnings` exercised directly
on the new record shape. **The three new regression tests in `tests/test_baseline_census.py` have not
been executed** and must be run on the data machine.

---


## 2026-07-30 (correction, same day) — Stage 1.5.2's gate was unsatisfiable; D2/D3 had already been measured

**Status:** ✅ Corrected. Markdown only; `src/` untouched. Supersedes the gate wording committed in
`75c331e` a few hours earlier.

Asked to confirm that 1.5.1 REV FINAL and 1.5.2 were both "ready and pushed". Checking rather than
confirming exposed a defect in the gate **I wrote today**.

### The defect

`STAGE_1_5_2_LABEL_ANCHOR.md` §0 read *"this stage does not start until D2 and D3 are closed"* and
described both findings as unmeasured. **Both were measured on 2026-07-24** by
`experiments/diag_zero_point.py`. I asserted they were open without checking whether the diagnostic
had already answered them — the same failure as the §9-R3 re-derivation of D1, one day later.

**And the gate was worse than merely redundant: it was unsatisfiable.** D2's scientific half cannot
be closed by analysis, so as written it would have blocked 1.5.2 permanently.

### What `diag_zero_point.py` actually returned

| | question | result |
|---|---|---|
| **M1** (D3's question) | does the clock read age on this data? | 🔴 **FAIL** — extreme contrast **11.8 yr** vs a **20.2 yr** bar on a true 53-yr gap, power **0.996** |
| **M2** (D1's question) | is the cross-batch zero-point driving the offset? | ✅ **NO_BATCH_EFFECT** — −2.99 yr, 95% CI [−13.12, +7.14], n=12 |
| **M3** (D2's question) | per-donor offset: real biology or `n=1` baseline noise? | ⚠️ **INDETERMINATE** |

**M3 in detail:** observed offset SD **16.4 yr** vs **12.3 yr** expected from a single unreplicated
baseline ⇒ the baseline explains **56% of the variance, 95% CI [9%, 100%]**, leaving 10.9 yr SD for
biology + batch + model. **That CI spans nearly the entire range — D2 is measured and unresolvable
at n = 6.** More donors would close it; more analysis will not.

Its recorded decision was **ESCALATE**: *"the clock does not separate the age extremes on this data,
so ΔAge's target is unvalidated… Stage 2's premise is void as stated. Do not proceed to Phase 2/3."*

### The corrected gate

Narrowed to the two halves that are actually closeable, neither needing new data:

* **G-a** — `_control_baseline` (`aging.py:81-90`) must record baseline **count and composition**.
  Stage 1.5 made `n=0` visible; **`n=1` is still silent**, and M-2b's RNA-side contrast inherits
  whatever that baseline does, so which donors rest on `n=1` must be visible in the output.
* **G-b** — donor chronological age parsed in `src/` (GEO declares N2/N3=0, Y1=29, Y2=35, O1/O2=53).
  *Unwired*, not unknown: `REV FINAL` §2 already **used** these values as its guard at MAE 4.0/4.4 yr.

**Explicitly recorded as NOT gates:** D2's scientific question (unresolvable at n=6 — carried as a
stated limitation), D1 (answered), and D3's scientific question (answered — M1 FAIL; only the wiring
remains).

**The relationship also runs the other way, and the earlier framing had it backwards.** M1's failure
is precisely what a methylation anchor resolves, so **1.5.2 is the response to that ESCALATE, not
something queued behind it.**

### Also corrected

`STAGE_1_5_1_REV_FINAL.md` §6's box (added earlier today) listed D2 and D3 as flatly "🔴 OPEN" and
repeated the "gated behind D2 and D3" framing. Both rows now carry the measured results and the
open/closed split, and the sequencing note states why D2's scientific half must never be written as
a gate.

**Header and `Depends on` lines in 1.5.2 updated to match.** Grep confirms no stale "gated behind
D2 and D3" wording survives outside the two self-corrections that quote it.

---


## 2026-07-30 — Adversarial audit of REV FINAL and STAGE 2 before external review; 5 defects fixed

**Status:** ✅ Executed and verified. **1 new artefact script + 2 annotated plan files.**
`src/` untouched (`git diff --stat src/` empty). Suite not re-run — no existing `.py` modified; the
new script is standalone pure-stdlib and was run directly, output below.

Asked to make both documents bulletproof before external critique. Audited for a **different**
failure mode than last time: not "is each sentence true" but **is the claim reproducible, is the
ledger complete, and does it survive being quoted out of context.**

### REV FINAL — three defects

**F1 — §4.3's "the intercept cancels EXACTLY" never admitted its own violations.** `anti_trafo` is
linear only at age >= 20; below that it is exponential and the intercept does **not** cancel. This is
why §3 and §4.3 disagree (−24.1 vs −24.5, −27.5 vs −28.3) — and the document never said so. A
reviewer comparing the two tables would have hit it immediately. **Measured:** 4 of 66 predicted ages
per clock fall below 20, all deeply-rejuvenated day-15/17 intermediates.

**The fix makes it a strength.** Where **zero** pairs violate the condition, the two forms agree to
**exactly 0.00** — contrast C on both clocks, contrast B on multi-tissue. Zero violations => zero
difference, every time, which *demonstrates* the algebra rather than asserting it. Max deviation is
**0.79 yr against an effect of −24 to −28**, and every deviation is **negative**, so the §3 headline
values are the **conservative** ones. Claim restated precisely instead of unconditionally.

**F2 — §4.4 had NO artefact, and it is the load-bearing defensive section.** §4.4 answers the one
challenge a reviewer is most likely to press (contrast A's post-hoc promotion, §7). Neither of its
two checks is produced by anything: `diag_methylation_anchor.py`'s `CONTRASTS` has **no
failing-intermediate arm and no dose-response**, and §9 nonetheless described the results JSON as
"full output". So the numbers were right but **unreproducible**.

Closed with **`experiments/verify_rev_final_4_4.py`** — an independent re-derivation from the raw
555 MB beta matrix, **pure stdlib, no numpy, no shared code** with the measurement script, so
agreement is corroboration rather than the same code answering twice. It reproduces a **known**
value first (V1) before its new numbers are trusted. **All four checks reproduce:**

| check | recomputed | document |
|---|---|---|
| V1 contrast A *(pipeline validation)* | −24.05 [−31.12, −16.98] / −27.55 [−33.69, −21.40] | −24.1 / −27.5 |
| V2 §4.4(a) failing intermediates | **−1.13 [−2.75, +0.49] / −3.62 [−5.07, −2.16]** | −1.1 / −3.6 |
| V3 §4.4(b) dose-response | **rho −0.885 p 0.0001 / −0.842 p 0.0006** | −0.885 / −0.842 |
| V3 slope | **−3.30 / −3.15 yr/day** | −3.30 / −3.15 |

Also found: **§4.4(b)'s slopes are the intercept-free form** (derived-intercept gives −3.10 / −2.77)
and the table never said which convention it used. Now stated. rho and p are identical either way.

**F3 — two sentences false or overclaimed when quoted out of context.** §8.1's *"There is no join
key, so a methylation age cannot be attached to any cell the model trains on"* is true of GSE165179
only — **GSE165178 joins 22/22** — and would have been quoted against §8.3. §4.5's *"ΔAge has a valid
anchor"* reads as "the labels are anchored", which §8.2 explicitly denies. Both scoped.

### STAGE 2 — two defects, both bearing on the wet-lab spend

**S1 — §1 and §4 quote different models and the file never says so.** §1's table is the **ridge**
baseline's shifts; §4's is the **model's**. Every donor disagrees (O1: +5.72 vs +0.64). Recorded in
`STAGE_1_DEVIATIONS.md` §C1, but a reader of Stage 2 sees only two contradictory tables. Mapping
table added; §4's are the ones to use.

**S2 — §2's headline benefit predates §4's rule, and nobody re-measured.** §2 reports
**14.3 → 6.9 (−52%)**; §4 then says the same T16 run *"helps 4 donors and hurts 2"*, so **§2's figure
is the UNCONDITIONAL correction applied to everyone.** §4's `|d| > 2·SE` rule exists to suppress some
of those corrections. The benefit of the stage **as specified** is therefore between −50% and zero
and is **currently unknown**.

**S3 — whether the rule fires at k = 3 was never computed, and it decides the spend.** Substituting
k = 3 into `SE ≈ 1.253·s/√k` gives `fires <=> |d| > 1.447·s`, where `s` is the within-donor sd of
`pred − true`. With `s ≈ 1.253 × 6.9 ≈ 8.65` (inferred from §2's corrected MAE, **not measured**),
the threshold is **12.5 yr** and the rule **fires for only 3 of 6 donors — N3, Y2, N2**. It correctly
declines both donors T16 damaged (O1, Y1), but **also declines O2, which T16 helped**; capturing O2
needs `k > 6.28·s²/d²` ≈ **11 cells**, not 3. So **"k = 3 minimum" is the number at which the
*unconditional* correction passed — it is not established for the *conditional* one this document
specifies.** All arithmetic re-verified numerically.

**Pre-registered before any spend:** measure `s` directly (needs no new cells); re-measure the
benefit with the rule active at k = 3 and k = 5; **grade §12's ≥25% TARGET bar on the conditional
number**, since that is what would ship. If it misses: raise k, or take §0's stated fallback of a
within-donor ranker — **not** relax the bar (ground rule §5), and **not** silently revert to the
unconditional correction that damaged O1 and Y1.

**Nothing was rewritten in either document.** All changes are additive annotation boxes plus two
scoping corrections to sentences that were false as written.

---


## 2026-07-30 — Two open findings were being carried silently; 1.5.2 gated behind them

**Status:** ✅ Annotations applied. **Markdown only — no `.py` file touched, `src/` untouched**
(`git diff --stat src/` empty). Suite not re-run: this machine has no `numpy`/`pytest`, and nothing
that could move it was modified.

Challenged on whether the original Stage 1.5 problem was ever actually solved, or whether stages
were being stacked instead of closed. Checked the record rather than answering from memory. **The
challenge was half right, and the half that was right matters more.**

### What WAS solved — the premise that this was all drift is wrong

Stage 1.5 asked one question: **is the ±12.7 yr per-donor offset an artefact of the silent
zero-point fallback at `aging.py:88`?** It was executed, pre-registered, and answered:

> **51 of 51 chunks carry >=1 vehicle control. The fallback never fired.**

It also delivered the 21 tests that four plan documents already *claimed* existed while **no test
imported `harmonize.py`** — `STAGE_6`'s acceptance gate had named a test that could never fail, and
`STAGE_5` had promised a reviewer a proof nobody wrote. And it corrected two overstatements by
measurement: "batch-immune by construction" is false (a per-dataset multiplicative gain survives),
and intercept cancellation is numerical, not bit-identical.

Nor did `REV FINAL` ever conclude "no fix needed" — its own §8.2 ledger reads *"are the ΔAge labels
now fixed? **no** — and this stage cannot fix them."*

### What was NOT solved, and was being carried silently

Stage 1.5 surfaced three findings. **D1 was measured and downgraded** (paired Exp1−Exp2 offset
**−2.99 yr, 95% CI [−13.12, +7.14], n=12, `NO_BATCH_EFFECT`** — structurally true but not
demonstrated to drive the offset; the ~10 yr CI half-width excludes a *large* effect, not a
meaningful one). **D2 and D3 are still OPEN**, have a fix plan recorded in
`STAGE_1_5_HARMONIZATION_AUDIT.md` §5, and **nothing was executed**:

* **D2** — every Gill donor's zero-point rests on **one unreplicated control sample**.
* **D3** — donor chronological age is **parsed nowhere in `src/`** (re-grepped: still zero hits)
  though GEO declares it (N2/N3=0, Y1=29, Y2=35, O1/O2=53).

### Two failures of my own, recorded rather than quietly patched

1. **`REV FINAL` §6 ("what this stage does not establish") omitted D2 and D3.** Every statement in
   that section is true; the **ledger of open work was incomplete**. A reviewer reading it cold
   would close it believing the harmonization arc was finished. I had validated whether each
   *sentence* was true, not whether the *list* was complete — and had called the document hole-free
   on that basis.
2. **`STAGE_1_5_2` §9-R3 re-raised D1 as a novel risk** when a measured answer already existed.
   Re-deriving a closed finding as a fresh risk is the same drift in miniature.

### Applied

| file | change |
|---|---|
| `plans/STAGE_1_5_1_REV_FINAL.md` §6 | **Additive box only; the body list is byte-unmodified.** Records D1/D2/D3 with status, states that the omission was the defect, and explains **why this stage's conclusions survive D2 by design** — every §3 contrast is a *paired arm comparison* (same donor, same day, methylation) that **never touches the RNA day-0 baseline**. Notes that immunity is not a fix, and that with D1 measured small and the clock convicted in §1, **`n=1` is now one of only two live explanations for the ±12.7 yr offset Stage 2 is premised on**. Also records that **D3 is *unwired*, not *unknown*** — this stage's own guard **used** those donor ages and returned **MAE 4.0 / 4.4 yr**, so they parse and are accurate; only `src/` ignores them |
| `plans/STAGE_1_5_2_LABEL_ANCHOR.md` | New **§0 GATE**: does not start until D2 and D3 are closed. Header status and `Depends on` updated. §9-R3 corrected to cite D1's measured −2.99 yr, with an explicit "do not over-read the null" |

**Sequencing rationale, recorded so it can be challenged:** D2 is cheaper than 1.5.2, needs no
download, and bears directly on **Stage 2's premise** — a bigger question than whether the RNA clock
is calibratable. D3 is nearly free and supplies ground truth that 1.5.2 would otherwise download
methylation to approximate. **D3's limit is stated in §0 so it is not oversold:** donor age is a
per-donor *constant*, so it cannot measure rejuvenation *within* a donor and does **not** make
1.5.2 unnecessary — it anchors the absolute-calibration question only.

**Next:** execute D2 and D3 (read-only metadata work — no downloads, no label changes).

---


## 2026-07-30 — Stage 1.5.2 written: the missing stage between "labels are wrong" and "correct them"

**Status:** 🔵 **PRE-REGISTERED, NOT EXECUTED.** New plan file only. No code, no data, no labels
touched. `src/` untouched — verified with `git diff --stat src/` (empty); the only changes are this
entry and one new `plans/` file. **The suite was NOT re-run:** this machine has no `numpy`/`pytest`
(they live on the data machine), and since **no `.py` file was modified** there is nothing here that
could move it. Per the convention at the top of this file, that is recorded as unverified rather
than asserted — an earlier draft of this entry claimed "455 tests still pass", which had not been
checked.

Asked which stage the newly-unblocked RNA↔methylation agreement test belongs to. Checking the stage
graph showed the honest answer is **none of them** — there is a real gap:

```
1.5.1   "the labels come from an instrument that fails here"   done
  ???   "here is whether that instrument can be repaired"      NO OWNER
Stage 2 "correct per-donor offsets ON those labels"            blocked on exactly that
```

`STAGE_2_LEVEL_CORRECTION.md`'s own annotation names the blocker: *"re-measure the per-donor level
shift on corrected labels before spending — that was never done, because the labels were never
corrected."* Verified against each candidate owner: 1.5.1 is closed and measurement-only (its §6.2
withdrew this test as ill-defined **on GSE165179**); Stage 2 consumes the labels; Stage 4 depends on
Stage 3; Stage 6 owns acquisition, not label anchoring. So a new file was added rather than
reopening a validated one — reopening a doc stamped EXECUTED to slip in a change that moves the
training target is precisely the silent-target-change the ground rules exist to prevent.

**Added:** `plans/STAGE_1_5_2_LABEL_ANCHOR.md`. Additive only; no existing plan file edited.

**The design point that decides the stage, closed in §4 rather than flagged.** The obvious test —
"does age_rna correlate with age_meth?" — is **not sufficient and would prove nothing alone**.
1.5.1 measured `corr(age_rna, pluripotency) = -0.62`, and methylation age also falls sharply during
reprogramming (-24 to -27 yr). Both modalities move with reprogramming progress, so a clock carrying
**zero** age information would still correlate strongly across a sample set whose dominant axis is
exactly that — the +36.5 yr identity artefact re-entering through the back door. The headline
correlation is therefore **barred from being a pass criterion**; the decisive readings are
within-arm (CD13 only, where cells are not reprogramming) and partialled on the existing
`OSKM_PLURIPOTENCY` signature, reused so it cannot be tuned for this stage.

**Geometry constraints recorded before the run, not discovered after:** GSE165178 has 4 donors x 3
days x 2 markers = 24 grid cells, 22 exist, and **no day-0 or untreated arm**. So the only internal
contrast is SSEA4 vs CD13, and CD13 is a *treated non-responder*, not an untreated control. 1.5.1's
inertness result (+0.5/-2.4) came from the **transient** arm; GSE165178 is the **Sendai** arm, so
that result **does not automatically transfer** and is carried as an explicit assumption with a
fixed response (§9-R1).

**Two errors in my own recollection, caught while writing and corrected in the file:**

1. I was about to cite ground rules "§10 negative controls / §11 shape-before-statistic". Those
   sections **do not exist** — `REF_GROUND_RULES.md` ends at §6. The negative-control and shape
   gates are **G1/G2 in `STAGE_4_VALIDATION.md`**. Cited correctly.
2. The value of this stage is **not** label volume, and the file says so up front (§2): Gill is
   **~75 of 33,688 age labels, about 0.2%**. HFF's ~99.8% stays unanchored either way because no
   public methylation exists for it. The deliverable is the **instrument verdict**, which is what
   Stage 2's premise and Stage 5's claims actually rest on.

**Pre-registration discipline (ground rule §5b).** Bars are stated with their intent, geometry and
resolvability recipe, to be frozen by `audit_metrics.bar_verdict` **before GSE165178 is opened**,
with each added to `tests/test_bars_resolvable.py` — a bar with no resolvability test is not
considered pre-registered. Anticipated in advance rather than discovered later: **n=11 per arm may
come back UNRESOLVABLE**, so the fallback (ρ_partial at n=22 becomes decisive) is fixed now so it
cannot be chosen after seeing data. M-2c is **gated** on M-2a — fitting a calibration to a clock
that is not tracking the target would manufacture a meaningless number.

**Three of the four pre-registered outcomes do not produce a label change**, and the file says that
is a real result: the two negative verdicts retire a route the project has already spent four failed
attempts on, on paired ground truth rather than argument.

**Phase 2 (the actual label change) is gated and separated** because it is the one place the
training target moves: one change only, snapshot with an *exercised* rollback, every Stage 1 guard
re-run and reported, and applying a Gill-learned correction to HFF **defaults to NO** — different
cell system, no ground truth to validate transfer against.

**Nothing here is executed.** The status line changes only when it has run.

---


## 2026-07-25 — **Stage 1.5.1 planned: clock precision (option B).** PLAN ONLY, nothing executed

**Status:** 📋 **PLAN ONLY** — `plans/STAGE_1_5_1_CLOCK_PRECISION.md`. No code written, no fit run,
`git diff --stat src/` empty.

**Why (the one number).** Stage 1.5 ran five measurements whose results only cohere one way: the
clock's own cross-validated error is **±12.27 yr**, and *every* effect the project measures is the
same size — per-donor offset **±12.7 yr**, D2's entire D0→D14 spread **13.1 yr**, the rejuvenation
effect to grade **~11 yr**. **SNR ≈ 1.** That is why nothing replicated, why E1b (+0.205) and D2
(−0.214) flipped sign on independent data, and why both verdicts were decided by hundredths.

Not a modelling, calibration or target-definition problem — a **precision problem in the
instrument**, upstream of all of them. ΔAge is the training label, so its noise is a hard floor on
`sigma_age`, conformal width, RES, and every quantitative claim.

**Root-cause candidates**, read from the artefact and `clock_fit.py`: dense `RidgeCV` over **33,155
genes from 133 samples** with no feature selection (R1); the compression signature is visible —
`cv_pearson 0.837` against `cv_mae 12.27`, and every donor reads high near the 72.4 intercept (R2); a
dense clock is fragile to the 57% gene-coverage gap (R3); possible CV optimism from scaling before
the split (R4); and out-of-range at age 0 is a **data** limit no refit can fix (R5).

**Bars set now, derived from the science.** ΔAge is a difference of two clock readings (noise
≈ √2·cv_mae), so an 11 yr effect at ≥2σ requires **`cv_mae ≤ 4.0 yr` = PASS**; 4–6 MARGINAL
(ranking only); **>6 FAIL** → pre-registered fallback (more clock training data / restrict claims to
large effects / ranking-only framing). Bars fixed *before* any fit so they cannot be tuned to the
result. Per §10's lesson, any verdict within **0.5 yr** of a boundary is reported **FRAGILE**.

**Structure:** cheap measurement gates expensive computation — audit the current fit under leak-free
nested CV, evaluate four pre-specified candidates (ElasticNet / fold-internal gene filtering / slope
recalibration / dense-ridge control) on one identical harness, and only then rebuild + retrain.

**Stated in advance:** changing the clock changes `y_age`, so the four-run `+0.000` guard streak ends
**by construction, not by defect**, and Stage 1's PARTIAL verdict does not automatically carry over.

---

## 2026-07-25 (§10 D2 EXECUTED) — **E1b does not replicate.** No reliable age trend either way

**Status:** ✅ Run on the data machine. Pre-registration committed *before* the run (`5360c24`).
Predicted *replicates* (~65%); **falsified**. `git diff --stat src/` empty.

GSE242423 (Kundaje lab, **single-cell**, different donor, different protocol — shares nothing with
Gill but the clock), 8 timepoints D0–D14: **rho −0.214** vs Gill's E1b **+0.205** → `CONTRADICTS`.

**⚠️ Fragile by 0.014.** Pre-committed boundary −0.20, cleared by 0.014, trajectory non-monotonic.
Verdict honoured as pre-registered, but it means **"failed to replicate"**, not "opposite effect
demonstrated".

**Establishes:** E1b's age-rise does not replicate → the **E1b escalation largely dissolves**, and
per the pre-registered branch **D1 does not run** (its spec stays locked, unused).

**Does NOT establish that ΔAge is valid.** Gill +0.205 and GSE242423 −0.214 are two *weak* effects of
*opposite sign*; with E1's NO_TREND (−0.064) the defensible claim is that **the clock reads no
reliable age trend during reprogramming, in either direction** — noise-dominated, not backwards.

**Two incidental corroborations:** HFF is a **neonatal** line yet reads **84.5 yr** at D0,
independently reproducing M1's ~+80 yr over-prediction on age-0 donors (different lab, different
modality); and iPSC reads **63.4**, ~20 yr below every fibroblast timepoint — the identity axis again.

**Methodological lesson:** E1b cleared its bound by 0.009 and D2 by 0.014, in opposite directions.
**Boundary-crossing verdicts at n≈6–8 are unstable.** Future bars should require a margin, or report
`FRAGILE` within a stated distance of the boundary. Three of my predictions in this arc are now
falsified; the consistent error is expecting cleaner signals than n≈6–8 can deliver.

**Also recorded (prior commit `5360c24`):** my earlier ruling-out of the OOD-gating idea via T15's
`AUC 0.47` is **withdrawn** — T15 measured `AUC(error → flagged)`, a different and finer task. The
real reason the shipped detector cannot implement it is structural: `train_model.py:291` fits the
reference on `train_ds`, which already contains the whole D0→iPSC trajectory, so it can never flag
it. It measures distance from the *model's* training data, not the *clock's* fitted domain.

---

## 2026-07-24 (Phase 1 EXECUTED) — **M1 FAILED. The clock does not read age on this data. ESCALATE.**

**Status:** ✅ **RUN on the data machine.** `python experiments/diag_zero_point.py "D:\Gill"` →
`diag_zero_point_results.json`. **The pre-registered prediction was FALSIFIED** — it predicted
`PHASE_2_AND_3` with M1 clearing. `git diff --stat src/` empty. Full record in the lab notebook
under *RESULT — PHASE 1*; summary in `plans/STAGE_1_5_HARMONIZATION_AUDIT.md` §7.

| Measurement | Verdict |
|---|---|
| **M1** clock vs chronological age | ❌ **FAIL** — extreme contrast **11.8 yr** vs bar **20.2** (true gap 53 yr) |
| **M2** Exp1/Exp2 batch offset | ⚠️ `NOT_ESTIMABLE` — but that verdict is a **stub**, see below |
| **M3** share of offset variance from one baseline | ⏳ **INDETERMINATE** as predicted — 56%, CI [9%, 100%] |

Per-donor: N2 (age **0**) → predicted **98.7**, i.e. older than both 53-year-olds; N3 (age 0) →
36.4. **Two donors of identical age read 62 yr apart.**

**Consequence (pre-registered branch):** ΔAge's target is **unvalidated**. This reaches past Stage
1.5 into **Stage 4**, and **Stage 2's premise is void as stated**. Phases 2–4 blocked.

**The failure is structured, which is the lead.** O1/O2 (both 53) agree to **0.4 yr**, and across the
four *adult* donors the old-vs-young separation is ≈18 yr against a true 21 yr gap. The catastrophe
is confined to the **neonatal** donors — and `fleischer_clock.json` was fit on **adult** dermal
fibroblasts (GSE113957), so age 0 is extrapolation outside its fitted domain. All six donors are
also over-predicted (+22.7 to +98.7). Hypothesis for the escalation: usable on adults, invalid on
neonates — which would leave two of six LOOCV folds with an unvalidated target.

### Two defects found while reviewing the run

1. **M2's verdict is a stub, and its claim is false.** `diag_zero_point.py:326` calls
   `m2_verdict([])`, so it always emits *"no matched (donor, day, marker) pairs … fix option (a) is
   impossible."* **Matched pairs demonstrably exist** — `N2_d11_CD13_Sendai_Exp1` and `…_Exp2` are
   both in the series matrix. The true statement is narrower: the *pipeline's `obs`* discards batch
   identity (D1), but the diagnostic can parse the titles and does not. This is the pre-registered
   *"M2 estimable after all"* branch. It does not change the ACTION (M1 short-circuits `decide()`).
2. **The rebuild is currently broken.** `run_multi_local.py:53` points `CLOCK` at
   `local_runners/configs/clocks/fleischer_clock.json`, **which does not exist**. `build_clock` fails
   loud (correct design), so a rebuild aborts at the clock step — the "we can always harmonize
   again" fallback does not currently work. Only `configs/clocks/fleischer_clock.json` is tracked.

Both are being fixed next, then Phase 1 is re-run and re-recorded.

### Both fixed, Phase 1 re-run — M1 unchanged, **D1 downgraded**

- **M2 now measures.** `parse_title()` / `group_matched_pairs()` read `(donor, day, marker, Exp)`
  from the series-matrix titles — the only place batch identity survives — plus 8 branch tests.
  Result: **12 matched pairs**, offset **−2.99 yr, 95% CI [−13.12, +7.14] → `NO_BATCH_EFFECT`.**
  The stub's claim is disproven by measurement.
- **`run_multi_local.py:53`** `CLOCK` now resolves to the tracked
  `configs/clocks/fleischer_clock.json`; **a rebuild is possible again.**

**This corrects finding D1, against the earlier claim.** The cross-batch zero-point was recorded as
"a real defect"; measured, the Exp1↔Exp2 term is **not distinguishable from zero**. D1 stays
structurally true (all baselines Exp2, ~50% of samples Exp1) but is **not demonstrated to drive** the
±12.7 yr offset — Phase 3 option (a) would have little to remove. **Not over-read:** the CI
half-width (~10 yr) is the same order as the offset, so this excludes a *large* batch effect, not a
meaningful one.

**Verdict unchanged: M1 still FAILS, ACTION remains ESCALATE**, Phases 2–4 blocked. The live
explanations are now the clock's validity and the `n=1` baseline.

### Escalation scoped (plan doc §8) — the severity is provisional, not settled

**M1 tested ABSOLUTE age; the model trains on ΔAge, which is control-relative.** For `age = w·x+b`,
ΔAge = `w·(x_pert − x_base)`, so the intercept, any additive per-donor baseline offset, and **every
gene Gill is missing cancel**. Measured this session: Gill covers 57% of the clock's genes / **89% of
its weight mass**, so **10.8% reads as zero** — an absolute-age error that vanishes in ΔAge. So M1's
failure proves the clock's *absolute* readings invalid on this data; it does **not** prove ΔAge's
target invalid — a separate, unmeasured question. §8.3 pre-registers the tests that settle it, first
of which (E1) is a within-donor age-*trajectory* check, the only one that bears on ΔAge. The failure
is also structured (O1/O2 both 53 agree to 0.4 yr; only the age-0 neonates N2/N3 blow up, below the
clock's ~1–94 yr fitted range), so it may localise to 2 of 6 folds. This is a self-correction of my
own framing, not of the M1 result — same shape as the D1 downgrade. Handoff: Stage 4 validation.

### E1 EXECUTED — NO_TREND. Prediction falsified; escalation now supported on both axes

`experiments/diag_e1_trajectory.py` (+13 branch tests, pre-registered and committed BEFORE the run
in `3a81cb6`). Predicted PASS (moderate); result **`NO_TREND`**.

- Primary (iPSC excluded): mean per-donor Spearman(age, day) **−0.064**, 95% CI **[−0.232, +0.104]**,
  4/6 negative, every |rho| ≤ 0.28. Adults-only also NO_TREND (−0.055).
- With-iPSC PASS (−0.179) **does not count** — carried by pluripotent endpoints (cell-type change,
  not aging), which E1 excluded on purpose.

M1 (absolute age) and E1 (within-donor change) now **agree**: on this data the frozen clock does not
demonstrably read the aging axis, absolute or relative. The §8.4 NO_TREND branch fires — **the deep
escalation stands**: ΔAge's target is unvalidated (Stage 4 / Stage 5); Stage 2's premise remains
void as stated.

Two caveats, so the null is neither over-read nor explained away: (1) null at n=6, but the per-donor
rhos are weak and sign-inconsistent, so low power is not the story; (2) the monotonic metric may be
mis-specified for Gill's **transient** (MPTR) protocol (OSKM withdrawn ~day 13 → non-monotonic
trajectory) — a limitation of my pre-registration, not grounds to dismiss the result. Next step, to
pre-register before running: **E1b** over the reprogramming phase only (days 0→~15). Stated guard:
not a retry until something passes — a null E1b plus E1 is strong evidence against ΔAge validity.
Until then ΔAge's rejuvenation signal is **NOT validated**. `src/` untouched.

### E1b EXECUTED — WRONG_DIRECTION. Escalation hardens; diagnostics stop

Pre-registered and committed before the run (`15ad575`, cutoff `REPROG_PHASE_DAY_MAX = 15.0` chosen
from the protocol, not the ages). Predicted ~45% PASS; result **`WRONG_DIRECTION`**.

- E1b (reprogramming phase, day ≤ 15): mean per-donor Spearman(age, day) **+0.205**, 95% CI
  **[+0.009, +0.401]**, 5/6 donors positive. In the OSKM window where cells should rejuvenate, the
  clock reads them getting **older**. Weak (CI lower bound at +0.009) but the wrong sign, robustly.

All four tests now agree the clock does not read the aging axis on this data: **M1 FAIL** (absolute),
**E1 NO_TREND** (full trajectory), **E1b WRONG_DIRECTION** (reprogramming phase); the only PASS is
with-iPSC, which is the fibroblast→iPSC *identity* axis, not aging. Coherent read: the clock tracks
identity (iPSC = young) but not rejuvenation during reprogramming, where it runs backwards. ΔAge —
computed mostly on non-iPSC reprogramming cells — is **not a validated rejuvenation target here.**

This is upstream of the whole model (ΔAge is its target), so it reaches into **Stage 4 / Stage 5**,
not just Stage 2 — the most consequential finding of the Stage 1.5 arc.

**Diagnostics stop.** Two trajectory tests were pre-registered and both failed; a third metric tweak
would be fishing and is **not** proposed. Next is a Stage 4 decision — is the frozen Fleischer clock
a valid ΔAge source for OSKM-reprogramming cells at all, and if not, what is a valid rejuvenation
target on this data? Caveats kept in view: n=6, bulk, fibroblast clock out of domain on reprogramming
cells — the finding is about this clock on this data, not about reprogramming biology. `src/`
untouched.

Also fixed while here: a pre-existing `N802` in `tests/test_diag_zero_point.py` that would have
failed CI's `ruff check src/ tests/ scripts/`, and two dead imports in the diagnostic. **CI lint is
still red from 11 other pre-existing errors elsewhere in `tests/` (e.g. `test_verify_1a.py:108`) —
not touched here.**

---

## 2026-07-26 (answered) — Why we did not match Gill: we DO, on the intermediates. ΔAge has a valid anchor.

**Status:** ✅ Executed. `experiments/diag_methylation_anchor.py` extended to both Horvath clocks and
all three contrasts; 28 tests; 455 pass; `src/` untouched.

Read Gill's methylation methods instead of speculating. Three facts: their comparison was
transiently reprogrammed **vs negative control fibroblasts** (what we did ✓); their optimum was
**13 days** (where our largest effect sits ✓); and they used **several clocks**, reporting
**multi-tissue *and* skin & blood both rejuvenated** — the ~30 yr being the **median across clocks**,
not one clock (✗ we had run one). So both clocks were run and the **intermediate** arm added.

| contrast | skin & blood | multi-tissue |
|---|---|---|
| **INTERMEDIATES** (still reprogramming) | **−24.1** [−31.1, −17.0] **REJUVENATION** | **−27.5** [−33.7, −21.4] **REJUVENATION** |
| transiently reprogrammed fibroblasts (MPTR) | −5.8 [−19.5, +7.9] NO_EFFECT | −9.4 [−18.3, −0.5] REJUVENATION_**FRAGILE** |
| **failed to reprogram (NEGATIVE CONTROL)** | **+0.5** [−2.3, +3.2] | **−2.4** [−5.7, +0.8] |

**The answer: Gill's ~30 yr is reproduced on the cells still in the reprogramming phase** — −24.1 and
−27.5 yr, both clocks, both highly significant, 12 identity-matched pairs. **We had been measuring
the wrong arm.**

**The negative control is inert on both clocks** (+0.5, −2.4; both CIs contain zero), so the design
is valid and the transcriptomic **+36.5 yr artefact is dead twice over**.

**The returned fibroblasts match Gill's SHAPE exactly** — by reprogramming length, skin & blood
−2.7 / **−14.1** / −0.6 / −5.1 and multi-tissue −13.4 / **−18.4** / −5.6 / −0.6 for 10/13/15/17 d:
**maximal at 13 days, diminished at 15–17**, precisely as Gill describe. The intermediates show the
*opposite* trend (−14 → −36 monotonically), which is simply "closer to iPSC" — the day-profile
distinguishes two different quantities.

**Coherent reading:** rejuvenation during reprogramming is large (−24 to −27.5 yr); retention after
return to fibroblast identity is partial and peaks at 13 days. First result in the arc where sign,
magnitude, shape *and* the negative control all agree.

**What must NOT be claimed:** the MPTR-fibroblast retention is **not solid** — skin & blood says
NO_EFFECT, multi-tissue says FRAGILE (CI bound −0.5), and the **intercept sweep flips it**. The
intermediates and the negative control are intercept-robust; this is not. Also: G2 is weak (MAE
4.0/4.4 vs a 5.0 tolerance on 3 known-age samples, implied intercept, 5.3/6.4 yr spread);
multi-tissue coverage is 94.6%; n = 9–12 pairs from 3 donors.

**Settles for the project:** ΔAge **has a valid anchor**; the transcriptomic clock's failure is fully
localised to the instrument (the biology is real, the RNA clock cannot see it); the +36.5 yr artefact
is closed. **Open:** how much rejuvenation survives the return to fibroblast identity.

---

## 2026-07-26 (unblocked) — The paired dataset exists: GSE165178 anchors our RNA labels

**Status:** ✅ Verified against real sample titles. No code changed; 455 tests pass. Plan §8 corrected.

Asked what the next step is if 1.5.1 passes review. Checking that found **my own previous answer was
wrong twice over**, and the correction unblocks the route I had declared closed.

Both of our series are SubSeries of SuperSeries **GSE165180**, which has **four** parts:

| accession | contents | have it? |
|---|---|---|
| GSE165176 | `[Sendai_RNAseq]` — the RNA we train on | yes |
| GSE165177 | `[Transient_RNAseq]` | no |
| **GSE165178** | **`[Sendai array]` — methylation on the SAME Sendai samples** | **no — get this** |
| GSE165179 | `[Transient array]` — §3's results | yes |

**GSE165178 pairs to our training data sample-for-sample.** Verified on the real titles, not assumed:

- 22 methylation samples titled `{donor}_{day}_{marker}` (e.g. `Y2_d11_SSEA4`);
- our RNA titles are the same key plus a batch suffix (`Y2_d11_SSEA4_Sendai_Exp1`);
- **22/22 join on `donor_day_marker`, zero unmatched**;
- donors **O1, O2, Y1, Y2** — 4 of our 6, the two missing being the neonatal N2/N3 that sit outside
  the clock's fitted range anyway; days 9/11/15;
- and the **sort marker IS the arm label** in our data: `CD13` → *Failing to reprogram fibroblast*
  (47), `SSEA4` → *Reprogramming fibroblast* (65). The arm assignment transfers unambiguously.

**What I had said, and why it was wrong.** §6.2 withdrew the RNA↔methylation agreement test (M-2) as
*"ill-defined"* — arms unmappable, overlap 2 donors × 2 days. §8.3 then said the remaining work needed
*"new profiling, since no public series pairs methylation to GSE165176."* **Both statements are true
of GSE165179 and false in general.** I checked only the series I had rather than the SuperSeries, and
generalised from it. The zero-overlap finding stands for GSE165179; it is not a property of the study.

**What this unlocks:** M-2 becomes well-defined and adequately powered (22 paired samples, 4 donors,
arms mapped); Gill's RNA labels gain a direct methylation anchor; and the calibration route — old
Step 3a, previously written off — is testable again. Whether a Gill-trained correction generalises to
HFF remains open, since HFF is a different cell system with no methylation anywhere.

**Next action:** download **GSE165178** (series matrix + processed beta matrix; check the format
first, as with the others), pre-register bars including resolvability at n=22/4 donors, then run M-2:
*does the transcriptomic ΔAge agree with the methylation ΔAge on the same samples?* Agreement ⇒ the
RNA clock is calibratable and ΔAge is recoverable for Gill. Disagreement ⇒ localises exactly where it
fails, against paired ground truth. Either answer is decisive and needs no new experiments.
**GSE165177** is worth taking at the same time — it pairs with GSE165179 and extends the comparison
to the transient arm.

---

## 2026-07-26 (correction) — REV FINAL §8 was not executable; the two series share zero samples

**Status:** ✅ Corrected. No code changed; 455 tests pass.

Asked whether anything in `STAGE_1_5_1_REV_FINAL.md` remained to execute. Checking that exposed an
error in its own §8, which instructed *"use methylation as the ΔAge source where methylation
exists."* **It is not executable.** Measured directly:

```
GSE165176 (RNA)   124 samples, e.g. N2_d11_CD13_Sendai_Exp1
GSE165179 (meth)   96 samples, e.g. O1_negative_control_15days_exp1
sample-title overlap: 0
```

The two series are **separate experiments** — no shared samples, different donor rosters
(N2/N3/Y1/Y2/O1/O2 vs O1/O2/O3), different day grids (7–47 vs 10–17), different arm vocabularies.
**There is no join key, so a methylation age cannot be attached to any cell the model trains on.**

**Consequence:** the project's ΔAge labels are **unchanged** by that stage and remain RNA-derived
from a clock the same stage proved is out of domain on reprogramming cells. The stage delivered
*knowledge* — the +36.5 yr artefact is closed, rejuvenation is real (−24 to −28 yr), the failure is
localised to the instrument — but **not labels**, and it cannot produce them.

§6 item 2 carried the same gap: it implied Gill's RNA samples could be anchored and only HFF could
not. With zero overlap, **neither** can. Both sections corrected, and §8 now carries an explicit
ledger of what the stage did and did not deliver.

**Executable state:** nothing remains in 1.5.1. **Stage 2** may proceed on the both-hypotheses
justification in its annotation; **Stage 3** depends on Stage 1 (required) and Stage 2 (optional) and
is not gated by any of this. Both open questions need **data**, not code: more donors with paired
methylation for the retention question (≈16 pairs), and methylation on the samples we actually train
on to anchor the labels — the latter meaning **new profiling**, since no public series pairs
methylation to GSE165176 or to HFF.

### ⚠️ Known intermittent test failure — recorded, not dismissed

A single test has now failed on **2 of ~15** full-suite runs, passing on every other run including
**5 consecutive** clean runs immediately after. The failing test's name was **not captured** either
time, so it is unidentified. Most likely a Windows temp-file lock (this repository has hit that
before), but that is **unverified**.

**Why this is recorded rather than waved off:** this project has already had one "flake" that turned
out to be a real batch-size-dependent defect. Anyone seeing a red suite should capture the test name
(`pytest -p no:warnings --tb=short` and keep the output) rather than immediately re-running, since
re-running is what has destroyed the evidence twice.

---

## 2026-07-26 (closed) — The three weak points in REV FINAL are now RESOLVED, not flagged

**Status:** ✅ All three challengeable points closed by measurement. 455 tests pass; `src/` untouched.
Plan is `plans/STAGE_1_5_1_REV_FINAL.md` (438 lines).

### 1. The derived intercept — **resolved algebraically** (§4.3)

Horvath's transform is linear above age 20 (`anti_trafo(x) = 21x + 20`), so for any pair of samples
both predicting >20 yr:

```
age_t − age_c = [21(lp_t+k)+20] − [21(lp_c+k)+20] = 21·(lp_t − lp_c)
```

**The intercept cancels exactly.** Every contrast here is a difference, so none depends on it.
Recomputed intercept-free: A **−24.5 / −28.3**, B **−6.0 / −9.4**, C **+0.5 / −2.4** — matching the
derived-intercept values to under a year. **This is no longer a robustness argument; it is algebra.**
The missing intercept row now matters only for *absolute* ages, which the document does not use.

### 2. Contrast A's post-hoc promotion — **corroborated by two things that could not have selected it** (§4.4)

**(a) An internal negative control specific to A.** Ran the previously unused
`Failing to transiently reprogram intermediate` arm against the same comparator:

| | A: transient-reprog intermediates | A-control: **failing** intermediates | paired A − A-control |
|---|---|---|---|
| skin & blood | −24.5 [−32.2, −16.7] | **−1.1** [−2.7, +0.5] | −23.3 [−31.1, −15.5] |
| multi-tissue | −28.3 [−35.4, −21.2] | **−3.6** [−5.1, −2.2] | −24.7 [−31.3, −18.1] |

Same OSKM exposure, same culture, same batch, same timepoints — but **failed**: ≈0 to −3.6 yr.
Succeeded: −24 to −28. **Rules out OSKM exposure per se, batch and culture duration.** (The small
real effect in failing cells matches Gill's own note that reprogramming-factor expression alone
rejuvenates some aspects.)

**(b) A dose-response.** Spearman(reprogramming length, effect) = **−0.885 (p=0.0001)** and
**−0.842 (p=0.0006)**, slope ≈ **−3.2 yr/day**, both clocks. **A contrast selected post hoc from
noise does not produce a monotonic dose-response at p<0.001.**

A now rests on four independent legs: two clocks, an internal negative control, a dose-response, and
the intercept-free formulation. The post-hoc promotion stays disclosed in §7 — it is simply no longer
the only support.

### 3. Contrast B's power — **corrected, and it was wrong in our favour's opposite direction** (§5)

An earlier statement in this record said *"MDE ≈13.7 yr, ~17 pairs needed, so n=9 is hopeless."*
**That used the skin & blood spread alone and applied it to both clocks.** Computed per clock:

| | sd | MDE at n=9 | observed | detectable? |
|---|---|---|---|---|
| skin & blood | 18.2 | **14.0** | −6.0 | ❌ |
| multi-tissue | 11.6 | **8.9** | −9.4 | ⚠️ **just barely** |

**On multi-tissue n=9 is already adequate** — which is the real reason one clock reaches significance
and the other does not. The two clocks do not disagree about the effect; they differ in precision.
Verified that this is a genuine instrument property, not luck on one contrast: multi-tissue is tighter
on **4 of 5** quantities including the untreated-control ages where no effect exists (**5.2 vs 6.6**),
though not uniformly (negative-control pairs favour skin & blood, 4.3 vs 5.1). Neither clock is
declared the winner; both are reported throughout.

Also checked: pooling across reprogramming lengths is **well specified** — day-heterogeneity p =
**0.852 / 0.255**, so averaging is not averaging over a varying effect.

**Honest state of B:** both clocks are consistent with a **real but small retention effect of ≈−6 to
−9 yr**. Neither excludes it; one detects it marginally. **Not established, and not dismissed** — the
question sits at its resolution boundary, and ≈16 pairs would settle it on both clocks.

### Also fixed

§10.4 and §10.6 contradicted the body after these changes (they still carried the superseded power
claim and the "mitigated by a sweep" framing). Both rewritten, with the supersession recorded in
place rather than silently edited.

---

## 2026-07-26 (executed) — Methylation anchor RUN: +36.5 yr artefact confirmed dead; no rejuvenation detected

**Status:** ✅ Executed on GSE165179. `experiments/diag_methylation_anchor.py` + 28 tests; `src/`
untouched. Horvath clocks exported to `configs/clocks/` from `biolearn` (Biomarkers of Aging
Consortium) and verified against the publications — 2013 multi-tissue **353** CpGs, 2018
**skin & blood 391** CpGs. **100% probe coverage** on all 96 samples.

**M-3 (the negative control) is the headline.** Failed-to-reprogram fibroblasts vs their matched
untreated control: **+0.5 yr, CI [−2.3, +3.2]** over 12 identity-matched pairs — *indistinguishable
from untreated cells*. The transcriptomic clock read the same cells at **+36.5 yr**. That artefact
is now **definitively dead**, measured against a real control with a sharp instrument rather than
argued about.

**M-1: no significant rejuvenation.** Transiently reprogrammed vs untreated control: **−5.8 yr,
CI [−19.5, +7.9]**, 5/9 pairs negative. By reprogramming-phase length: 10 d −2.7, **13 d −14.1
(Gill's optimum, n=2)**, 15 d −0.6, 17 d −5.1. The sign is right and largest at the pre-registered
optimum, but underpowered. **Neither a demonstration nor a clean refutation.**

**Bug found and fixed mid-run:** the first pairing required a unique sample per (donor, day, arm),
which silently discarded **6 of 9** M-1 pairs — GSE165179 runs every condition as `exp1` **and**
`exp2`. The first run therefore reported only 3 day-10 pairs. Replicates are now averaged, with a
regression test naming the defect.

**Robustness:** the coefficient tables carry no intercept row, so one was implied from the three
known-age day-0 samples — making G2 partly self-fulfilling. Swept the intercept from −0.60 to +0.70:
**M-1 stays −5.2 to −6.0 and M-3 stays +0.5 throughout.** The conclusions do not depend on it,
because most predictions sit above age 20 where Horvath's transform is linear and a constant cancels
in a difference.

**Honest limits:** G2 passes but barely (3 known-age samples, MAE 4.0 against a 5.0 tolerance, 5.3 yr
implied-intercept spread); we do **not** reproduce Gill's ~30 yr and the reason is not established;
n=9 pairs from 3 donors; a ~15 yr effect is not excluded. The treated arm is far more heterogeneous
(age sd 14.9) than the control arm (7.1) — and since M-3's CI is only ±2.8 yr on the same
instrument, that spread is most likely real biological variability, not measurement noise.

---

## 2026-07-26 — STAGE 1.5.1 REV FINAL written: anchor ΔAge to methylation (GSE165179)

**Status:** PLAN ONLY, nothing executed. `plans/STAGE_1_5_1_REV_FINAL.md`. All five prior 1.5.1
documents left byte-unmodified (verified).

Four fixes were proposed across V1/V2/V3/review; all four were tested and all four failed. The
reason is now established rather than suspected: **the transcriptomic clock is correctly built and
correctly applied, but out of domain on reprogramming cells — and no RNA-only analysis can fix that,
because every RNA route to "age" runs through that same clock.** Scoping the claim is not an option
(user decision), so the instrument problem must be solved.

**The plan:** acquire **GSE165179** — Gill's own multi-omic companion, 96 Illumina MethylationEPIC
samples, *same experiment and donors* as our GSE165176 RNA data — and use it as the
identity-independent anchor. It is a public download: no wet lab, no new samples, no GPU.

**Why methylation rather than a fourth RNA dataset:** different molecular layer, independently
validated clocks, it is what Gill used for the ~30 yr claim, and it is ~4× more precise —
Horvath skin & blood ≈3 yr vs Fleischer 12.27 yr, taking ΔAge SNR from **1.7 to ≈7**. That precision
also dissolves the n=1 day-0 baseline problem (`corr(baseline, ΔAge) = −0.986`), which is a
consequence of clock imprecision.

**Three pre-registered measurements:** M-1 does methylation show rejuvenation in responders
(resolvability checked: MDE ≈5.2 yr at n=6, so a −30 yr effect is overwhelmingly detectable);
M-2 does transcriptomic ΔAge agree with it (decides whether the RNA clock can be *calibrated* rather
than replaced — this is what covers the 79% of labels HFF holds, which methylation cannot reach);
M-3 the negative control (§10), which directly adjudicates our +36.5 yr against Gill's reported
*"moderate reduction"* in the same cells.

**Load-bearing guard G2:** the methylation clock must first reproduce known chronological age on the
six day-0 samples — the same in-domain check that vindicated the transcriptomic clock. If it cannot,
it is not an anchor either.

**Four-way decision fork** pre-registered, including the outcome that would be bad news
(M-1 NULL = the effect is not there at ~5 yr resolution, a publishable finding escalating to
Stage 4/5) and the one that would mean a bug hunt (CONTRADICTS).

**Adopted from V3 unconditionally:** A1 (stop pooling non-responders) and A3 (never test a dip with
a monotonic statistic) — already ground rules §10/§11.

---

## 2026-07-26 (stress test) — self-correction; conclusion strengthened; the anchor exists (GSE165179)

**Status:** ✅ Stress-tested the previous entry after being asked "are we completely sure?". The
conclusion survives; **its justification did not and is corrected.** Detail in
`plans/STAGE_1_5_1_REVISED_REVIEW.md` §6–§7.

**🔴 Self-correction.** The decomposition `−28.3 = +8.2 − 36.5` used two terms measured against the
**same n=1 day-0 baseline** — the noisiest quantity in the dataset. Measured:
`corr(day-0 baseline, responder ΔAge) = −0.986`, because responder ages cluster tightly (66–86,
sd≈7) while baselines scatter (36–99, sd≈21). That decomposition was baseline-contaminated and did
not support the claim I drew from it. *(It also cuts in the plan's favour: `R − F` is baseline-free,
a real advantage the review had not credited.)*

**Baseline-free re-test — same conclusion, now properly supported.** Within-arm, day 7 → peak
window, using no day-0 sample:

| arm | change | 95% CI | donors negative |
|---|---|---|---|
| responders | **−1.1 yr** | [−14.7, +12.5] | 4/6 |
| non-responders | **+17.5 yr** | **[+6.0, +28.9]** | **0/6** |

Only the control arm moves. And the gap is **already −9.7 yr at day 7** (CI [−18.5, −0.9]), widening
to −28.3 — the widening (≈18.6) is accounted for by the non-responder rise (17.5). So the contrast
measures a **standing population difference plus the control arm deteriorating.**

**Stronger than before:** the baseline-free CI [−14.7, +12.5] **excludes a Gill-scale −30 yr effect**
in responders. The earlier "56% power" caveat applied to the day-0-referenced test only. An effect
≲15 yr is still not excluded.

**Instrument noise, for the record:** mean responder age by day — 7→78.7, **9→101.3**, 11→77.2,
13→78.3. A +23 yr swing in two days, reversed in two more. The only clean monotone signal in the
series is the late approach to iPSC (34→47 d): the identity axis again.

**Scoping the claim is not an option (user decision), so the instrument must be fixed — and the
anchor exists.** **GSE165179**: Gill's own multi-omic companion, **96 Illumina MethylationEPIC
samples**, same experiment and donors as our GSE165176 RNA data. Methylation clocks (Horvath
skin & blood 2018, fitted on fibroblasts) do not read pluripotency the way a transcriptomic ridge
does, and methylation is how Gill established the ~30 yr claim.

**Next step is a data acquisition, not a code change:** get GSE165179 and ask first whether
methylation and transcriptomic ΔAge agree on the six donors' responder arms. Agreement vindicates
ΔAge and resolves the escalation with the target intact; disagreement localises exactly where the
transcriptomic clock fails. Links to Stage 6.

---

## 2026-07-26 (A1/A3 re-run) — labels not rescued; the "-28.3 rejuvenation" is entirely the control arm

**Status:** ✅ Executed on real Gill data. `experiments/diag_e1_corrected.py` + 21 tests; 405 tests
pass; `src/` untouched. Detail in `plans/STAGE_1_5_1_REVISED_REVIEW.md` §5.

Ran the one test both plans agree on: fix only the **undisputed** errors from
`STAGE_1_5_1_REVISED.md` — A1 (stop pooling the 47 non-responders into the treatment arm) and A3
(window contrast instead of a monotonic Spearman on a dip) — while **keeping the current day-0
control**, so the result is independent of the disputed control swap.

**Responders vs their own day-0: NO_EFFECT at every window** (10–13 d: **+8.2 yr**, CI [−20.1,
+36.5]; also 7–9, 13–15, 15–21, 21–29). Leave-one-donor-out STABLE; every point estimate positive.
Non-responders are significantly **AGEING** from day 10 on (+36.5 … +44.6).

**The decomposition is the finding:** `−28.3 = +8.2 − 36.5`. **100% of the "rejuvenation" comes from
the CONTROL arm rising, 0% from the treatment arm falling.** This is arithmetic — immune to sample
size and to the identity-adjustment argument. Redefining `is_control` to the non-responder arm would
define ΔAge as "how much less the reference inflates."

**Power, stated honestly:** at n=6 (sd 27.0) power for a Gill-scale −30 yr effect is **56%**, so the
null is underpowered, *not* proof no rejuvenation exists. What is informative is the direction — the
estimate is +8.2, not negative-but-short. The decomposition carries no such caveat.

**All three candidate label definitions now fail** (day-0 as built; day-0 with A1/A3 fixed;
non-responder control). This is the honest end of the transcriptomic-only route on this dataset —
Gill used **methylation**, and their paper states existing transcription clocks "failed to accurately
predict the age of our negative control samples."

**Next step is not a label change:** either an identity-independent anchor (second modality → Stage 6)
or scoping the quantitative rejuvenation claim. The fate/safety head (PR-AUC 0.99) is untouched.

---

## 2026-07-26 (review, tested) — R4 refuted BY RUNNING; Step 1 effectively complete; C3 eliminated

**Status:** ✅ Tests run on the real GSE113957 (~3 s compute). Recorded in **`plans/STAGE_1_5_1_NEW_CHANGES.md`** (a companion file — the original
`STAGE_1_5_1_CLOCK_PRECISION.md` is left **byte-identical**). No code changed.

R4 was challenged rather than accepted, so it was tested three ways instead of argued:

| Test | Result |
|---|---|
| does any cross-sample statistic exist? | `normalize_counts(X[:5])` vs `normalize_counts(X)[:5]` → **max diff 0.0**; no `StandardScaler` anywhere. **The scaler R4 describes does not exist** |
| is there a *different* leak (group leakage from repeated donors)? | **143 samples, 143 unique cell ids, 0 repeats** — none |
| reproduce the CV directly | **12.67 yr / ρ 0.841** vs the artefact's 12.27 / 0.837 — **it reproduces** |

**`cv_mae = 12.27` is honest.** R4 is wrong in mechanism *and* conclusion. The error is
understandable — scaler-before-split is *the* classic leak, and the code reads
`Xn = normalize_counts(X)` right above `cross_val_predict` — but **normalisation ≠
standardisation**: per-sample library size uses one row's own total; per-gene standardisation uses
a column statistic across samples. Right instinct, wrong identification. Removing R4 removes the
one route by which the SNR problem could have been *overstated*.

**The same run answered two more Step 1 items and killed a candidate:**

- **R2 confirmed:** predicted-vs-true slope on held-out folds = **0.717** (bar S1 wants 0.85–1.15).
- **R1 confirmed:** `alpha` = 0.272, near the *bottom* of `logspace(-1,4)` → penalty barely binding;
  with in-sample 0.77 vs CV 12.67 (**16×**), this is memorisation, not mild over-regularisation.
- **C3 eliminated:** slope recalibration (fitted out-of-fold) gives **12.78 yr — 1% *worse***.
  Rescaling a memorising model's slope does not remove its error. Step 2 should run C1/C2 vs C4 only.

Also noted: the artefact says `n_samples = 133` but GSE113957 has **143** — 10 samples were excluded
when the clock was fit, and which ten is not recorded. Minor, but it means my earlier "in-sample"
reproduction covered 133 of 143.

---

## 2026-07-26 (review) — Reviewed Stage 1.5.1: diagnosis endorsed, R4 refuted, OOD route withdrawn

**Status:** ✅ Review only, no code changed. Recorded in **`plans/STAGE_1_5_1_NEW_CHANGES.md`**; the original plan doc is left byte-identical.

Pulled 3 commits (D2 replication + the 1.5.1 plan) and checked the plan against the tree.

**Endorsed:** the SNR≈1 diagnosis is correct and the pre-registered bars are sound (`√2·cv_mae ≤ 5.5
⇒ cv_mae ≤ 3.9` recomputed and correct; artefact confirmed dense ridge, 33,155 non-zero weights).

**❌ R4 is factually wrong.** It claims `cross_val_predict` runs on data "standardised before the
split". **`clock_fit.py` has no cross-sample scaler at all** — the only transform is
`normalize_counts`, which is per-row (library size), and `RidgeCV` is refit inside each fold. So
`cv_mae = 12.27` is already leak-free and Step 1's "re-measure leak-free" will find nothing. This
*removes the hope* that the SNR problem was overstated. Step 1 keeps its other two items
(error-by-decile, predicted-vs-true slope); the in-fold guard stays correct for the new C1/C2
candidates, where feature selection genuinely must be inside the fold.

**New evidence strengthening R1:** the §9 reproduction gave **0.77 yr in-sample** vs **12.27 yr CV**
— a **16× gap**. That is memorisation, not mild over-regularisation, and it is the most direct
evidence in the record that the ridge penalty is not binding at 33k features / 133 samples. It also
predicts C3 (slope recalibration alone) will underperform — rescaling a memorising model's slope
does not remove its error.

**OOD-detector route withdrawn** (I proposed it before §10; recorded rather than dropped). Three
independent failures: (1) the detector is a Gaussian over the *model's* latent fitted on `train_ds`,
which **contains** the reprogramming intermediates — they are in-distribution by construction and
would never be flagged; (2) its measured AUC is **0.47** (chance), already documented at
`train_model.py:288-290`; (3) §10's D2 showed the trajectory sign **flips** between datasets
(+0.205 vs −0.214), so "reprogramming is out-of-domain" is not a stable property — it is noise at
SNR≈1. Gating would also flag the entire use case, making it option C in disguise.

---

## 2026-07-25 (§9 reproduction) — Gold check RUN on GSE113957: clock reproduces age (0.77 yr, ρ0.99). H1 refuted.

**Status:** ✅ Ran locally against the real NCBI GSE113957 files. `_load_known_age_fibroblasts`
rewritten for the NCBI layout + unit-tested; 391 tests pass; `src/` untouched.

The `run_reproduction` gold check now has its number. On 143 dermal fibroblasts (ages 1–96), the
frozen clock through the production path (`normalize_counts` → `LinearClock`):

| | |
|---|---|
| MAE | **0.77 yr** |
| Spearman / Pearson (pred vs age) | **+0.99 / +0.99** |
| weighted gene coverage | **100%** (33,155/33,155) |
| verdict | **REPRODUCES** |

**H1 (mis-applied) is definitively refuted** — the pipeline applies the clock correctly. *Honest
scope:* 0.77 yr is in-sample (clock fit on GSE113957), so it confirms application correctness, not
generalization; generalization is carried by H2's **out-of-sample** Gill result (+18/21 yr, ρ+0.60).
Application-correct **and** generalizes to held-out fibroblasts. The M1/E1/E1b escalation is fully
explained by out-of-range (neonatal, age 0) + out-of-domain (reprogramming) inputs. **ΔAge stays.**

Loader work: reads NCBI raw counts (GeneID × GSM), joins GeneID→Symbol via the annotation table,
dedups duplicate symbols by highest total (matching the clock's `dedup: highest_expressed` and
guaranteeing unique symbols so `predict_age` can't double-count), and merges GSM→age across both
platform series matrices. New pure helpers (`parse_age_value`, `series_gsm_to_age`,
`dedup_symbols_highest_total`) unit-tested, plus an end-to-end loader test on synthetic NCBI files.

*(Run against the four files in Downloads; move them to `D:\GSE113957\` and re-run the full
diagnostic on the data machine to regenerate the complete results JSON with the reproduction block.)*

---

## 2026-07-25 (§9 EXECUTED) — Clock validity scored: the escalation was OVER-READ. Clock is in-domain OK; ΔAge stays.

**Status:** ✅ Run on `D:\Gill` (reproduction check skipped — GSE113957 absent). Scored against the
pre-registration; full record in the notebook (*RESULT — §9*) and plan §9.5. `src/` untouched.

**The §7–§8 escalation ("ΔAge target unvalidated; clock can't read age") is refuted.**

| Check | Result |
|---|---|
| **H2 in-range tracking** | ✅ **TRACKS** — +18.0 yr for a 21 yr gap, Spearman +0.60. M1 failed by anchoring on the two age-0 neonatal donors, below the clock's [1,96] range (N2 read 98.7) |
| H1 coverage | DEGRADED — 57% of genes but **89.2% of weight** (the fibroblast genes that matter are all present) |
| H1 intercept dominance | MOVES — predictions SD 21.4 yr, not collapsed onto the intercept |
| H1 CP10k denominator | STABLE — 4.3 yr |
| H3 attribution | DIFFUSE — OSKM 0.007%, cell-cycle 0.65%, senescence 2.7%; **none of my categories explain it** |
| ACTION | **IN_DOMAIN_OK_INVESTIGATE_REPROGRAMMING** |

**H3 reframed by the exact shares:** the +20.1 yr reprogramming drift is the residual of ≈150 yr of
positive vs ≈130 yr of negative gene contributions — the clock summing a transcriptome in upheaval,
i.e. **extrapolating outside its training distribution (out-of-domain by nature)**, not a
marker-gene confound. My H3 prediction (OSKM/cell-cycle) was wrong; recorded. This is exactly what
the model's existing OOD detector is for.

**Standing conclusion:** the clock reads fibroblast age; ΔAge's instrument is valid in-domain. On
reprogramming intermediates it is out-of-domain — an interpretation limit, not a broken target. Fix
options A/B/C/D **not triggered**; C (retreat to fate) off the table. **ΔAge stays.**

**Next (priority):** (1) GSE113957 reproduction — turns the n=4 H2 into an n=133 powered
validation; (2) domain-aware ΔAge via the existing OOD detector.

---

## 2026-07-25 (§9) — Clock-validity diagnostic: is the clock broken, mis-applied, or out-of-domain?

**Status:** ✅ Written and tested (**384 tests**, 28 new). ⏳ **Not yet run on the data machine** —
needs `D:\Gill` (and optionally `D:\GSE113957` for the gold reproduction check). `src/` untouched.

### Why

The §7–§8 escalation ("ΔAge target unvalidated; diagnostics stop") **reaches past its evidence**.
An independent re-check found three confounds that each produce M1/E1/E1b's failures *without* the
clock being wrong about aging, and none were ruled out:

- **H1 (mis-applied):** `predict_age` sums only over genes present (`weights.get(g, 0)`). The clock
  has **33,155 genes**; if the Gill matrix misses many, most of the model is silently dropped and
  predictions collapse toward the 72.4 intercept — which is exactly where the M1 ages (36–99)
  cluster.
- **H2 (out-of-range):** M1's "young" anchor was the two **neonatal donors (age 0)**, below the
  clock's fitted range [1, 96]; N2 read 98.7. Among **in-range adults** the day-0 contrast is
  **~18 yr for a ~21 yr true gap** — the clock tracking in-domain fibroblast age.
- **H3 (out-of-domain):** days 0–15 of OSKM are cells leaving fibroblast identity; a "reads older"
  signal driven by pluripotency/cell-cycle genes is cell-STATE, not aging.

Also corrected: the `with-iPSC` config **PASSES** (6/6 donors, p=0.0295) — the clock has real
signal; E1b is marginal (**p=0.0445**); and E1 is underpowered (needs ρ≈−0.4, the transient effect
gives ρ≈−0.1) so its null is uninformative. All four fix options assume the instrument is broken —
this settles that first.

### `experiments/diag_clock_validity.py` (new) — read-only, four independent axes

| Check | Verdicts | Decides |
|---|---|---|
| H1 gene coverage | OK / DEGRADED / CRIPPLED (by fraction of \|weight\|, not gene count) | is the clock fully applied |
| H1 own-domain reproduction | REPRODUCES / DEGRADED / BROKEN / SKIPPED | is it applied correctly on known-age fibroblasts |
| H2 in-range tracking | TRACKS_IN_RANGE / NO_IN_RANGE_TRACKING (neonatal excluded, median split) | does it track age it was fit to read |
| H3 directional attribution | OUT_OF_DOMAIN_CONFOUND / AGING_GENES_DRIVE_IT / DIFFUSE | is the reprogramming reversal cell-state |

`decide()` folds them: CRIPPLED/BROKEN → **FIX_APPLICATION** (recoverable, ΔAge stays as-is);
in-range TRACKS + reprogramming CONFOUNDED → **TARGET_RECOVERABLE_DOMAIN_FIX**; clean application +
no in-range tracking → **GENUINE_CLOCK_LIMITATION** (the first point at which A/B/D are earned, not
assumed). Bars pre-registered (§5b); predictions recorded in the notebook before running.

Supporting checks: intercept dominance, and CP10k-denominator sensitivity (predictions normalised
over the full data gene set vs the clock-overlap set only).

### Tests — `tests/test_diag_clock_validity.py` (28)

Every verdict branch, plus hand-worked coverage math (weight ≠ gene count), the in-range median
split that excludes the out-of-range neonatal donors (the M1 error), attribution counting only the
positive "age rises" contribution, and the full `decide()` table including application-fix
priority. Found one bug while writing them: the in-range contrast used single min/max donors
(dropping Y2); switched to a median group split — more robust at n=4.

---

## 2026-07-24 (Phase 1) — Stage 1.5 Phase 1 written: zero-point diagnostics, bars pre-registered

**Status:** ✅ Written and tested (**332 tests**, 28 new + 1 new registry bar). ⏳ **Not yet run on
the data machine** — it needs `D:\Gill`. `src/` untouched. Status ledger with ✅ markers added to
the head of `STAGE_1_5_HARMONIZATION_AUDIT.md`.

### `experiments/diag_zero_point.py` (new) — read-only, decides whether Phase 3 is needed at all

| Measurement | Question |
|---|---|
| **M1** | does the frozen clock read chronological age on *this* data? (GEO donor ages 0,0,29,35,53,53) |
| **M2** | is there an Exp1/Exp2 batch effect? (all six baselines are Exp2 — finding D1) |
| **M3** | how much of the ±12.7 yr offset *variance* could one unreplicated baseline explain? |

All repo-data imports are confined to `baseline_ages()`, so the decision logic is data-free and
unit-tested. Uses the production normalisation (`normalize_counts` → log1p-CP10k), the space the
frozen clock was fitted in, so predicted ages are comparable to its own CV MAE.

### Bars pre-registered and resolvability-checked BEFORE running (§5b, tightening T1)

| | Bar | Null | Correct system passes | |
|---|---|---|---|---|
| M1 | contrast ≥ **20.2 yr** | a clock reading **nothing** | **99.6%** | ✅ RESOLVABLE |
| M2 | paired 95% CI excludes 0 | no batch effect | depends on pair count | ⚠ CONDITIONAL |
| M3 | share ≥ 50% or < 25% | — | — | ❌ **expected UNRESOLVABLE at n=6** |

**M1's bar is deliberately not "contrast > 0"** — a clock reading pure noise clears that half the
time. It is set at `z₀.₉₅ × SE` under the null (20.2 yr), which a correct clock still clears 99.6%
of the time. The 29-vs-35 middle contrast is **not gated**: at 6 yr it is half the clock's error.
M1 is now an entry in `tests/test_bars_resolvable.py` (registry extended with a `higher` kind).

**M3 is pre-registered as expected-unresolvable**, which is the point of checking bars forward: the
point estimate is **56%** (observed offset SD 16.4 yr vs 12.27 yr from a single baseline, residual
10.9 yr), but the χ² CI at 6 donors spans ~**[9%, 100%]**. Saying so *now* stops a 56% point
estimate being mistaken for a finding later.

### Tests — `tests/test_diag_zero_point.py` (new, 28)

Every branch every function can emit, including the ones we hope not to see (**M1 FAIL → ESCALATE**,
which reaches past Stage 1.5 into Stage 4) and the one we expect (**M3 INDETERMINATE at n=6**,
asserted on the real level-shift values). Also pins that M1 does **not** gate on the underpowered
middle contrast, and that a Phase 3 decision states it reopens **both** Stage 1 targets (T4).

**Predicted outcome, recorded in the lab notebook:** `PHASE_2_AND_3` with **no** Phase 3 lead —
M1 clears, M2 is NOT_ESTIMABLE (D1 discards batch identity, so option (a) is unavailable until
Phase 2), M3 indeterminate. Useful precisely because it says *instrument first*: the data needed to
choose a Phase 3 option does not exist yet.

---

## 2026-07-24 (review) — Independently reviewed the Stage 1.5 work; verified it, then tightened its plan

**Status:** ✅ Review complete, 303 tests pass, `src/` untouched. Recorded as
`STAGE_1_5_HARMONIZATION_AUDIT.md` **§6**. No code changed by this entry.

The Stage 1.5 tests, verifier, Group E run and fix plan (5 commits) were re-checked **against the
tree rather than taken on trust**.

**Every checkable claim held:** clock `cv_mae_years` = **12.2688** (claimed 12.27); `donor age`
genuinely unused (0 grep hits; `_parse_series` reads only `days of reprogramming` + `cell type`);
Exp1/Exp2 identity genuinely discarded (docstring only); `git diff --stat src/` **empty** across all
5 commits; **21/21** new tests and **303** total pass; Group E **51/51** chunks with the fallback
never fired, **covering all six LOOCV donors** so the PASS is not vacuous; every Gill donor carries
**exactly 1 control**.

**The tests are not decorative — they were mutation-tested.** Four defects were injected into `src/`
(variance floor removed; control branch killed; `sigma_ref` dropped from the Gill Projection;
`_align` made positional) and each was caught by the right test. `src/` restored after each.

**§5 corrected the Stage 1.5 doc, and was right to.** §2 Group A had specified intercept
cancellation as *bit-identical*; it is not — `(age+b) − mean(age_ctrl+b)` re-rounds, so the
cancellation is numerical (~1e-14), not symbolic. Reproduced independently.

**One concern raised, then dismissed by checking:** the verifier counts controls per *chunk* while
production groups per *cell_line within* a chunk. Every source emits one chunk per cell line by
construction (`sources.py:364/459/507`), so the check is exactly equivalent — not a defect, though
the invariant is unasserted (tightening T5).

**Five gaps closed in the fix plan (§6.2):**

| | Gap | Tightening |
|---|---|---|
| T1 | Phase 1's measurements had implicit bars and skipped `bar_verdict` — violating the §5b ground rule adopted the day before | each of M1–M3 gets a pre-registered bar, a resolvability check and an entry in `tests/test_bars_resolvable.py` **before** running. M1's power stated: SE(diff) 12.27 yr vs a 53 yr contrast ≈ 4.3σ |
| T2 | M3 was measured but appeared in no decision | given its own decision table; it is the quantity that should *size* Phase 3 |
| T3 | option (c) is largely **redundant** — `sigma_scale_factor` already fits `sigma_age` to residuals that contain the baseline error | keep only if made **per-donor**; as written it should be struck |
| T4 | option (a) needs matched Exp1/Exp2 pairs that may not exist; and Phase 3 changes `y_age`, so it reopens **both Stage 1 targets**, not just the guards | M2 must report pair counts first, "(a) impossible" is a permitted outcome, and Stage 1's PARTIAL verdict must be declared re-openable **before** Phase 3 |
| T5 | the gate's chunk↔line assumption is unasserted | group by `raw.obs["cell_line"]` (one line, with Phase 2) |

**Verdict:** the work is sound and the plan is now concrete. **Phase 1 is the correct next action** —
read-only, cheap, and able to escalate past this whole stage if M1 shows the clock does not read age
on this data.

---

## 2026-07-24 — Stage 1.5: the harmonization claim made true, and the ΔAge zero-point gate built

**Status:** ✅ **FULLY EXECUTED on the data machine.** Groups A–D: 21 new tests pass, full suite
**303**, ruff clean. Group E: **PASS — 51/51 chunks carry ≥1 control; the `aging.py:88` fallback
never fired.** `src/` **untouched** (`git diff --stat src/` empty), so no guard can have moved.
Predictions were pre-registered in the lab notebook *before* the run and were confirmed.

### Group E result, and the finding it surfaced

| Source | Chunks | Controls per chunk |
|---|---|---|
| GSE242423 HFF | 45 stratified batches | **111–112** of ~980 cells |
| Gill | 6 donor chunks | **exactly 1** of 19–21 cells |

**Ruled out:** the ±12.7 yr per-donor offset is not an artefact of the self-centring fallback.
**Surfaced, and still open:** every Gill donor's zero-point rests on **one unreplicated control
sample**, so any error in that single day-0 measurement propagates 1:1 as a per-donor additive
offset — the same shape as the effect Stage 2 is premised on, and not distinguishable from it by
anything measured so far. Read with deviation **C1** (the ±12.7 is the *ridge* shift; the model's
mean shift is −5.71, 95% CI [−22.9, +11.5], including zero), the Stage 2 premise is weaker than
"established biology". A finding, not a defect to patch here — and exactly what Stage 2's k≈3
reference cells per donor would address.

**Why.** Four plan documents assert cross-modality harmonization is "unit-tested" with "intercept
cancellation **proven**" (`MASTER_PLAN.md:48`, `STAGE_5_PUBLICATION.md:127`,
`STAGE_6_NEW_DATA.md:143`) — and **no test imported `harmonize.py`**. `STAGE_6`'s acceptance gate
therefore named a test that could never fail, and `STAGE_5` promised a reviewer a proof that was
never written. This stage makes the existing claim *true*, not weaker.

| File | Change |
|---|---|
| `tests/test_harmonize.py` (new, +21) | **A** intercept / `mu_d` / `mu_ref` cancel; additive batch offset immune. **B** the exact closed form. **C** fit leak-safety, variance floor, sorted-intersection gene space, `MIN_REPLICATES` / unknown-dataset / missing-reference raise, `_align`, JSON round-trip. **D** per-line zero-point **and the silent fallback pinned**. **E** every branch of `decide_verdict` |
| `verify_stage1_5.py` (new) | the runnable gate. A **pure** `decide_verdict()` separated from all I/O (the `verify_1a` lesson — a decision function whose only exercised path says PASS is not a gate), plus a read-only replay that censuses vehicle controls per chunk and writes `verify_stage1_5_results.json` |

### Two overstatements the tests corrected

1. **"batch-immune by construction"** (`harmonize.py:9`) is false as written. ΔAge is immune to
   *additive* batch effects but carries a per-dataset multiplicative **gain**:
   `ΔAge = Σ_g δ_g · sigma_ref,g / (sigma_d,g + EPS) · w_g`, now pinned as a closed form. The same
   raw δ gives a *different* ΔAge in a dataset with different spread — measured, not argued.
2. **"intercept cancellation is bit-identical"** (plan, Group A) is not exact. The cancellation is
   *numerical*, not symbolic — `age + b` then subtracting a control mean re-rounds. Immune to
   ~1e-12; `np.array_equal` fails. Found by writing the test the plan asked for.

### Defects found in my own draft before commit, not after

`sys.modules` dataclass-load crash under `importlib` (collection error); an inverted
gene-intersection fixture that asserted the wrong answer; the over-strict `array_equal` above; one
`UP017`. Recorded because "the tests passed" is only meaningful if the first draft did not.

**Deliberately NOT done:** the wording fixes to `harmonize.py`'s docstring and the two
reviewer-facing rows (`STAGE_5:127`, `STAGE_6:143`) are **proposed, not applied** — plan §4 makes
that the user's call, not this stage's.

### Fix plan recorded (PLAN ONLY — nothing executed)

Following the Group E census into the Gill metadata produced three findings, and a fix plan is
recorded in **`plans/STAGE_1_5_HARMONIZATION_AUDIT.md` §5** — appended after the original
pre-registration (§0–§4 left exactly as written, never substituted), so the plan as written and
what actually happened stay auditable side by side in one file. **Not** a new document. The lab
notebook carries only a pointer to it, so the two cannot drift.

- **D1 — the zero-point is cross-batch.** All six baselines are `*_Fib_Sendai_`**`Exp2`**, while
  ~**50%** of every donor's treatment samples are **Exp1** (10 per donor). Half of `y_age` is
  therefore `age(Exp1) − age(Exp2 baseline)` — a batch term inside the target's *definition*.
- **D2 — baseline replication is invisible.** `_control_baseline` records neither count nor
  composition; Stage 1.5 made `n=0` visible, `n=1` is still silent.
- **D3 — `donor age` is parsed nowhere** (grep: zero hits) though GEO declares it
  (N2/N3=0, Y1=29, Y2=35, O1/O2=53) — the only ground truth able to test whether the clock reads
  age on this data.

**The number that makes it urgent:** the clock's own metadata carries `cv_mae_years = 12.27`, and
the per-donor offset Stage 2 exists to correct is ±12.7 (ridge) / 13.12 (model). The offset is the
size of **one** clock measurement's error — and each donor's zero-point **is** one clock
measurement. Not proof it is noise; proof the two are currently indistinguishable.

Plan is sequenced measurement-first (M1 clock-vs-chronological-age, M2 Exp1/Exp2 batch effect, M3
bound the noise share) with pre-registered branches, so the cheap measurements decide whether the
rebuild-and-re-score change is needed at all. Explicitly left alone: the ΔAge definition, the
clock's weights, Stage 1's calibration and its four-run `+0.000` guard record, the Exp1 samples,
and every prior record.

---

## 2026-07-23 — Made "audit the bar before the run" a ground rule, not a lesson learned twice

**Status:** ✅ Written and tested (282 tests). The transferable win from the Stage 1 scoring saga:
coverage and `fate_ece` were both audited *after* they misfired. This turns that into a forward
habit — every new acceptance bar is checked for **resolvability before it is pre-registered**.

| File | Change |
|---|---|
| `plans/REF_GROUND_RULES.md` | new **§5b** — a pre-set bar (§5) must also be RESOLVABLE: simulate a system that meets the intent EXACTLY at the grading geometry and confirm it passes ≥ 95% *before* registering the bar. Cites both Stage 1 cases (fate_ece 26.9% → pool → 99.6%; coverage 93% confirmed). No existing rule renumbered. |
| `audit_metrics.py` | new `bar_verdict(null, bar, …)` → **RESOLVABLE / UNRESOLVABLE** against `MIN_PASS_RATE = 0.95`; docstring section on forward use. `resolvability()` was already the reusable core. |
| `tests/test_bars_resolvable.py` (new, +10) | one entry per registered TARGET bar, asserting a correct system's pass rate matches its required verdict. Includes the **retired** per-fold `fate_ece` bar asserting it stays UNRESOLVABLE — the lesson made executable — and one assertion that pooling flips the same bar's verdict. Adding a bar means adding an entry here; a bar with no entry is, by rule, not pre-registered. |

Bug caught while writing the tests: my first `higher_is_better` case expected RESOLVABLE at a
90% pass rate — but 90% < `MIN_PASS_RATE`, so UNRESOLVABLE was correct. The code was right; the
test expectation was wrong. Fixed the test.

**This does not touch any run, bundle, or scorecard column** — it is process + one helper + tests.

---

## 2026-07-23 (latest) — Wrote the Stage 1.5 plan doc (harmonization & ΔAge zero-point audit)

**Status:** ✅ Plan document committed. ⏳ The stage itself is **not run** — this is its
pre-registration, nothing under `src/` or `tests/` is touched yet.

Stage 1.5 existed only as an out-of-repo plan file; now `plans/STAGE_1_5_HARMONIZATION_AUDIT.md`
records it in the repo. It is a **measurement-only** stage (0 lines change in `src/`) that sits
between Stage 1 (closed) and Stage 2. It exists because four plan docs assert harmonization is
"unit-tested" / "intercept cancellation proven" (`MASTER_PLAN:48`, `REF_ARCHITECTURE:20`,
`STAGE_5:127`, `STAGE_6:143`) while **no test exercises `harmonize.py`** — the Stage 6 gate names
a test that does not exist. Reading the module surfaced two concrete facts the audit pins:

1. ΔAge cancels additive batch effects (`mu_d`, `mu_ref`, clock intercept) but carries a
   per-dataset **scale gain** `sigma_ref/(sigma_d+EPS)` — so "batch-immune by construction" is an
   overstatement, and Group B asserts the exact invariant instead.
2. `_control_baseline` has a **silent fallback** ([aging.py:88](src/cellfate/data/aging.py)): a
   donor in a chunk with no vehicle controls is self-centred, forcing its mean ΔAge toward 0.
   Whether this fired on the real build is what distinguishes the ±12.7 yr per-donor offset being
   real biology (Stage 2's premise) from an artefact — Group E checks it directly.

The doc specifies the test groups (A: the promised intercept proof; B: the true scale invariant;
C: fit/leak-safety; D: the ΔAge zero-point incl. the fallback; E: real-data replay of `plan_all`)
and a `verify_stage1_5.py` gate mirroring `verify_1a.py`. No existing plan doc was modified —
additive, in the style of `STAGE_1_DEVIATIONS.md`.

---

## 2026-07-23 (latest) — Repaired the calibration target and re-scored Stage 1 against it

**Status:** ✅ **Re-run on the data machine** (`rescore_results.zip`, commit `0003ff8`). 273
tests pass there. The live `scorecard.py snapshot --tag B_fatecal_pooled` printed the pooled
block **ECE 0.211 / floor 0.091 / excess +0.121 / 100th pctile** — identical to the offline
prediction from `diag_dump/` to **0.00e+00**. Guards vs the pre-repair `B_fatecal` snapshot:
`max|Δ| = 0.00e+00` on all four, so the additive scorecard change did not perturb any measured
value. `baseline` (pre-repair snapshot) correctly reports pooled ECE `n/a`.

### Why

`fate_ece` is graded as the mean of per-fold ECEs over ~21 held-out cells in 10 bins. Measured
(`audit_metrics.py`): a **perfectly calibrated** model scores 0.183 and clears the 0.169 bar only
**26.9%** of the time. The criterion was measuring the sample size, not the model. Pooling the
held-out cells across folds — the more correct LOOCV estimate, since every cell is still
predicted by a model that never saw it — raises that to **99.6%**.

### `scorecard.py` (the user's file; additive only, no existing metric changed)

| Change | Detail |
|---|---|
| `measure_fold` stores `_fate_S` / `_fate_y` | raw per-fold safe probabilities and labels. Underscored, so `METRICS`-driven tables ignore them |
| new `pooled_fate_ece(folds)` | pooled ECE + **floor** + **excess** + null percentile. Returns `None` for snapshots predating it, so `baseline.json` still loads |
| `_print_snapshot` / `cmd_compare` | print the pooled block; compare shows both snapshots' raw ECE **and** excess |
| `cmd_compare` header | states that the paired CI's sensitivity comes from the **consistency** of a change across folds, not the metric's own spread, and that a heterogeneous change can be large in the mean and still read as noise |

**`floor`** is the median ECE a perfectly calibrated model with that exact probability vector
would score (`y ~ Bernoulli(p)`, so all of it is estimator bias). **`excess = ece − floor`** is
the only quantity comparable across calibrators: raw ECE also moves when a calibrator merely
*sharpens*, because sharper probabilities sit in extreme bins where the floor is lower. On run 3,
**75%** of one apparent improvement was exactly that.

### Stage 1 re-scored

| | per-fold **[as graded]** | pooled **[repaired]** |
|---|---|---|
| `fate_ece` | 0.249 | **0.211** |
| floor | 0.179 | **0.091** |
| excess | +0.071 | **+0.121** |
| pass rate for a *correct* system | **26.9%** | **99.6%** |
| vs bar 0.169 | MISS (uninterpretable) | **MISS (real, 100th pctile of null)** |

**The verdict does not change, which is the point** — repairing the instrument could not have
been goalpost-moving, because Stage 1 fails either way. What changes is that the failure is now
*interpretable*: at 100% of the null it is unambiguously real, not an artefact of n≈21.

**Stage 1 final: PARTIAL.** `conformal_coverage` PASS (0.889 pooled marginal; audited at 93.0%
pass rate for a correctly-90% system). `fate_ece` MISS. Four guards +0.000, bit-identical, three
runs running.

### Tooling added this session (all read-only w.r.t. runs; logged here late — the changelog rule was missed on the first three)

| File | Purpose |
|---|---|
| `dump_pool_diag.py` (+9 tests) | reads back `xdonor_only_safe_ece_insample` / `shipped_safe_ece_on_pool`, computed by run 3 and printed nowhere |
| `dump_diag_bundle.py` (+8 tests) | packages pool + calib + test arrays, raw **and** calibrated, into a ~2 MB sendable dump so calibrators can be refitted offline instead of by retraining |
| `diag_calibrators.py` (+11 tests) | compares calibrator families by leave-one-donor-out **within** the pool; reports ICC / effective n |
| `audit_metrics.py` (+12 tests) | asks of every criterion: how often does a system that satisfies the intent EXACTLY get reported as passing |
| `tests/test_scorecard_pooled.py` (+9) | pins the repair, above all that `excess` calls a purely sharpened model **worse** |

Two defects found by writing those tests: `donor_ids_from_counts` must refuse to reconstruct pool
donor labels when residual and fate row counts disagree (it returns `None` rather than guessing);
and a boundary bug where `0.250 - 0.230 = 0.019999999999999990` reported a gain of exactly the
threshold as below it.

---

## 2026-07-23 (later) — Diagnostics read. Three of yesterday's claims retracted; the bar is below the estimator floor

**Status:** ✅ Analysed `diag_dump/` from the data machine. Pipeline reproduces the graded
`fate_ece` from raw probabilities to **0.00e+00**. Full detail in the lab notebook under
*RUN 3 POST-MORTEM*. **No source changed.**

| Retracted | Replaced by |
|---|---|
| "the bar is fair and attainable, ~2× the 0.078 floor" | Floor recomputed on the **actual** P(safe) vectors is **0.183**. A perfectly calibrated model clears 0.169 only **26.9%** of the time. The bar is below what n≈21 × 10 bins can resolve. |
| "the union fit cost the target; revert to the pool-only principle" | Union **excess +0.071** vs pool-only **+0.144** vs identity **+0.192**. The principle would have been twice as bad. The shipped calibrator is the best candidate tried. |
| "P(safe) saturates, so the top ECE bin cannot move" | **0.0%** of test rows exceed 0.99; P(safe) spans 0.09–0.88. Near-perfect *ranking* (PR-AUC 0.992) does not imply saturated *probabilities* — that inference was wrong, and the family hypothesis built on it is dead. |

**The metric rewards sharpening.** An other-donor refit appeared to take ECE 0.249 → 0.103,
seemingly beating its own 0.179 floor — impossible. Sharpening (a = 3.4–5.7) moves probabilities
into extreme bins where Bernoulli variance is smaller, **lowering the floor**; 0.110 of the 0.146
apparent gain (75%) is that artefact. Recorded so the one dishonest route to "landing" the bar is
closed explicitly.

**Excess over own floor is the comparable quantity.** By it, Stage 1 removes **63%** of the
miscalibration present with no calibration at all (+0.192 → +0.071) — the effect the stage was
built to produce, on a metric that can show it.

**Where the residual lives:** base rates are calib 0.514, pool 0.64, test 0.754. The calibrator is
fitted for a 0.51-safe world and graded on a 0.75-safe one; that is *label shift*, uncorrectable
from source data. Per-fold, the failure concentrates on **Y1** (base rate 0.579 vs 0.76–0.86
elsewhere) — the same donor heterogeneity behind N3's 0.333 coverage. **Stage 2's subject.**

**No further calibrator change is pre-registered.** Family right, fitting set right, residual not
a calibration problem.

---

## 2026-07-23 — RUN 3 executed and scored: PARTIAL. `fate_ece` misses; the bar it was set from was measuring a stacked calibrator

**Status:** ✅ Run on the data machine (229.0 min, 6/6 folds, 222 tests pass). Scored, logged in
`experiments/DELTAAGE_LAB_NOTEBOOK.md` under *RESULT — RUN 3*. **No code changed by this entry.**

**Verdict against `STAGE_1_CALIBRATION.md` §3:** 5 of 6 criteria met.

| Role | Metric | Bar | Result | |
|---|---|---|---|---|
| TARGET | `conformal_coverage` | 0.85–0.95 | 0.401 → **0.889** | ✅ |
| TARGET | `fate_ece` | ACCEPT + ≥40% drop (≤0.169) | 0.281 → **0.249** (−11.0%) | ❌ |
| GUARD ×4 | `fate_prauc`, `fate_roc`, `rank_model_dage`, `dage_mae_model` | noise | all **+0.000** | ✅ |

Guards bit-identical for the third consecutive run — Stage 1 provably does not touch the model.
`interval_width` 17.7 → 65.9 reads REGRESSION but is not a guard; widening is the pre-registered
consequence of an honest `q`.

### The finding: `fate_ece_platt` is a stacked layer, not an alternative calibrator

`scorecard.py:189` fits its Platt on `S_cal` and applies it to `S` — and `S`
(`scorecard.py:157`) is the **predictor's output**, which already has the bundle's calibration
applied (`predictor.py:170`). So `fate_ece_platt` measures **bundle calibration + a second
calib-fitted layer**, not a standalone in-distribution Platt.

It lands at 0.140–0.161 in all three snapshots regardless of what the bundle ships (baseline
0.153, A_xdonor 0.161, B_fatecal 0.140). **The second layer was doing the work in every T8.2
number.** The run-3 prediction of ≈0.15–0.17 was derived from 0.153 as though a single-layer
bundle calibrator could reach it. It could not. Prediction falsified; the reason is a
specification error on my side, recorded rather than re-rationalised.

### The bar was checked before being blamed, and it holds

`fate_ece` is estimated on 19–21 cells over 10 bins, so estimator bias could in principle have put
0.169 below its resolution. Simulating a perfectly calibrated model (`y ~ Bernoulli(p)`) at run-3's
geometry gives a floor of **0.078** (90% range [0.057, 0.105] for the 5-fold mean, `P(≥0.17)=0.0%`).
The bar sits at ~2× the floor. **It is attainable; 0.249 is a real miss.** The bar is not moved.

### Why the union fit under-delivered

`total=4509 in_dist=4406 xdonor=103` → the cross-donor pool is **2.28%** of the fit. Shipped slope
`a` = 2.599 ± 0.024 across folds; the pool-only diagnostic slope = 1.380, ranging 1.144–1.542. The
shipped slope being ~1.9× larger *and* far tighter across folds is the signature of a fit
determined by the 4406 rows the folds share, not the 103 that differ. The union is the
in-distribution fit to three digits — the deviation from *"calibrate on the deployment regime"*
that was flagged when it was made, and it cost the target.

### Not explained

A synthetic probe of the two calibrator families failed to reproduce the observed gap (it made
`LogisticRegression`-on-raw-`p` *worse*). The boundedness hypothesis — logistic-on-`p` cannot
exceed `sigmoid(w+c)` while logit-Platt drives saturated inputs to exactly 1.0 — is unconfirmed
and nothing below depends on it.

### Reporting gap found (cosmetic, not fixed yet)

`retrain_stage1.py:249` prints `ECE pre`/`ECE post` from `xdonor_ece_before_temp`/`_after_temp`,
which apply `softmax(logits / temperature)`. Stage 1 sets `temperature = 1.0` whenever Platt is
fitted, so those two columns are now **identical by construction** — which is exactly what run 3
printed (0.269/0.269, 0.294/0.294, …). Not a calibration bug; the summary table is showing a
guaranteed no-op and hiding `xdonor_safe_ece_before`/`_after`, the binary figure that matches what
`scorecard.py` grades. Fix belongs with the next change, not on its own.

---

## 2026-07-22 — Full audit of the session's code; one real guard bug found and fixed

**Status:** ✅ 221 tests, smoke 34/34. Everything committed and pushed.

A line-by-line audit of everything changed this session, run against live code rather than by
re-reading it. Most of it confirmed what was claimed; one thing did not.

### 🐛 The bug: calibration could move the rank GUARDS

Platt is monotone, so it can never *reorder* cells — but it can **merge** them, and a merged pair
changes a rank metric. Two mechanisms, both measured:

| mechanism | effect |
|---|---|
| `EPS` clamp at `1e-6` | collapsed **4 of 8** float32-representable values near 1.0 onto one number (float32's ulp there is ~6e-08, so the clamp was coarser than the input's own resolution) |
| casting calibrated probs back to **float32** in `_summaries` | merged values the map left distinct. At slope 20: **PR-AUC 1.000 → 0.941, ROC-AUC → 0.966** |

`_PLATT_BOUNDS` permits a slope up to `1e2`, so a steep fit is reachable on real data. Had one
occurred, `fate_prauc` would have shown a **REGRESSION** — a Stage 1 guard — and the correct
response under §3 is to *revert*. We would have reverted a working change because of a rounding
artefact.

**Fixed:** `EPS` → `1e-9` (two orders below the float32 ulp, so every representable input except
exact 0/1 survives distinct), a numerically stable sigmoid for the wider logit range this admits,
and `_summaries` no longer downcasts — `_rows` converts to Python floats and `res.py` upcasts to
float64 anyway, so nothing downstream wanted the narrower type. Guards now hold at slopes 2, 8,
20 and 100; `test_calibration_does_not_move_the_rank_guards_even_at_a_steep_slope` pins it.

**Claims corrected.** Four places said Platt makes the rank guards "mathematically invariant" or
"bit-identical". That was too strong — monotone means *no reordering*, not *no merging*. All four
now say what is true, in `CHANGES.md`, the lab notebook, `smoke_stage1.py` and
`common/calibration.py`.

### Verified, not assumed

| check | method | result |
|---|---|---|
| biology untouched | `git diff 18d7e69..HEAD -- src/cellfate/data/ models/ evaluation/` | **empty** — clock, harmonization, fate labels, ΔAge targets, network all unchanged |
| column binding | indices 0–5 vs pre-session | `X_I…AM_I` still 0–5, `DONOR_I` appended |
| donor never a feature | `forward(x, u, dose_time)`; grep for `DONOR_I` | only in grouping logic |
| **row alignment** | rebuilt a dataset, compared every donor code against the shard's `cell_line` | **144/144 rows match** |
| Platt recovers miscalibration | 3× sharpen, +1.8 bias, and both | recovered `a`,`b` within 0.02 of the true inverse; mean\|p−p_true\| ≈ 0.002 |
| simplex invariants | saturated / zero / uniform input | finite, rows sum to 1, in [0,1], loss:death ratio preserved |
| schema guards | negative slope, half-specified pair | both rejected |
| back-compat | legacy `TemperatureParams` / `ConformalParams` | load unchanged, `sigma_scale` 1.0, both modes allowed |
| xstats round-trip | save → load | all seven arrays plus both dicts |

### Scope check on real bundles

Retrained the six rehearsal folds with the current code and compared against the same folds
trained *before* Change A″:

```
conformal_q  (N2)  0.47744181752204895  ->  0.47744181752204895
sigma_scale  (N2)  7.795770789209797    ->  7.795770789209797
temperature        1.498                ->  1.0   (Platt replaces it)
```

**Bit-identical** — the calibrator change provably does not reach `q` or `sigma_scale`. This is
the same check to run on the real data when run 3 lands.

### Held-out comparison (synthetic, 3 folds × 10 cells — weak, directional only)

| | mean ECE on a truly held-out donor |
|---|---|
| no calibration | 0.161 |
| **cross-donor temperature** (what run 2 shipped) | **0.190 ← worst** |
| in-distribution temperature | 0.160 |
| pool-only Platt | 0.172 |
| **union Platt (shipped)** | **0.153 ← best** |

Cross-donor temperature being worst independently reproduces run 2's regression on data it was
never fitted to. The synthetic setup does not reproduce the real miscalibration magnitude
(baseline 0.281 there vs ~0.16 here), so this is **directional support, not a prediction** that
run 3 clears the bar.

---

## 2026-07-21 — Stage 1 run 1 was INVALID; bulk-corpus guard added

**Status:** ✅ **Fixes written, NOT yet run.** Run 1 executed fully (6 folds, 212 min) and is void.

**What happened.** `cell_line` is not donor. The training split merges the **GSE242423 HFF corpus
(33,613 cells)** with the **six Gill donors (~14 cells each)**, and both are labelled by
`cell_line` — so the inner-LODO rotated over HFF as a seventh donor. Holding HFF out left a model
trained on **75 cells** (val_loss 33.0 vs the deployed 5.3), and because that fold is also the
largest it contributed **33,613 of 33,688 pooled residuals (99.8%)**. `q` and `sigma_scale` were
therefore calibrated against data starvation, not donor shift.

The tell: `sigma_scale` ranged **6.28 to 74.45** across folds for a quantity that should be
similar. Y2's 74.45 implies a median ensemble spread of 0.50 yr against a P90 residual of 36.9.

**My defect, not just the plan's.** `verify_1a.py` *detected this and printed the warning* — "MORE
than the expected 5; saw 6. THIS IS THE DANGEROUS DIRECTION" — and then **graded the run `PASS`**,
because the verdict logic only escalated to STOP on *too few* donors. The operator followed a PASS.
Cost: 3.5 h of GPU time and a void experiment. A check that fires and is then overruled by its own
scoring rule is worse than no check.

| File | Fix |
|---|---|
| `src/cellfate/training/xdonor_calib.py` | `MIN_INNER_TRAIN_FRAC = 0.5` — skip any inner fold whose held-out donor leaves <50% of the training split; raise if <2 usable folds survive |
| `verify_1a.py` | `STOP` when any donor holds >50% of a training split, **or** when the donor count differs from the expected 5. Both were previously PASS-with-warning |
| `tests/test_training.py` | two regression tests: a 90%-dominant donor must be skipped and must not reach the residual pool; a 95/5 split must raise |

**Bars unchanged** — this is ground rule §6 ("the default assumption is a bug in the test"), not a
retroactive threshold move. Run 1 numbers, per-fold coverage, and run-2 predictions are recorded in
the lab notebook.

**What run 1 did establish:** the guards behaved exactly as predicted, including the sharper
bit-identical prediction — `dage_mae_model` and `rank_model_dage` moved **+0.000 on every fold**.
Stage 1 provably does not touch the model. `fate_prauc` moved 0.992→0.988, which is *correct*: `S`
is `softmax(logits/T)[:,0]` and 3-class softmax is not rank-preserving in one class under a
temperature change.

---

## 2026-07-22 — The "flaky" test was real: batch-size float sensitivity, now pinned

**Status:** ✅ Fixed and verified — **7 consecutive clean full-suite runs (220 tests)** against a
check that failed 2-of-3 before the fix.

The previous entry logged a transient two-test failure and attributed it to a Windows file lock.
**That was wrong.** Chasing it properly found a real numerical property.

### Finding it

Rather than hope it recurred, I replaced the guess with a stronger check —
`test_batch_size_does_not_change_any_row`, which sweeps several batch sizes instead of comparing
only batch-of-5 against singletons. It failed **immediately and repeatedly**, converting a
1-in-N flake into a deterministic signal.

### The cause — upstream of this change, and not a defect

Measured on a trained bundle:

```
RAW ensemble probability (no calibration)   max |batch24 - single| = 8.9e-08
after Platt (slope a ~ 8)                                          = 5.0e-07   (5.5x)
sigma_age (multiplied by sigma_scale ~12)                          = 1.2e-06
```

torch selects different CPU kernels for different batch sizes, so identical rows differ in the
last float32 ulp **before any of this code runs**. Two shipped factors then amplify it: Platt
works in logit space so it multiplies by roughly its slope, and `sigma_age` is scaled by
`sigma_scale`. Both magnitudes are numerically irrelevant.

**The defect was the assertion, not the arithmetic.** `test_batch_and_single_agree` asserted
`model_dump() == model_dump()` — bit-exact float equality, a guarantee torch never made. It
passed by luck; the amplification exhausted the luck.

### The fix

Agreement is now asserted to a **relative** tolerance (`rel_tol=1e-4`, `abs_tol=1e-7`), not an
absolute one. Absolute was tried first at `1e-6` and **still failed** — on `sigma_age`, whose
scale and amplification differ from a probability's. An absolute bound would need re-tuning
whenever a fitted parameter moves, which is how tests rot. Relative does not: float32 carries
~1.2e-07 relative precision, amplification is capped by the Platt slope bound (1e2) and
`sigma_scale`, so ~1e-5 is the ceiling and 1e-4 leaves an order of magnitude.

This keeps every defect the test exists for — misaligned rows, leaked state, bad indexing all
move values by O(0.1–1) **relative**, four orders above the bound.

Also added `test_platt_clip_bounds_the_logit_blowup`: `P(safe)` values that round to exactly 1.0
in float32 would give an infinite logit and a NaN probability, and this model saturates there
routinely. The `EPS` clamp is load-bearing, and now documented as such in
`common/calibration.py` along with the amplification-scales-with-slope property.

---

## 2026-07-22 — Stage 1 run 2 scored; Change A″ calibrates `P(safe)`

**Status:** ✅ Run 2 **executed and scored**. Change A″ written and tested locally (218 tests,
smoke 32/32); the real-data run is pending.

### Run 2 result

| role | metric | bar | result |
|---|---|---|---|
| GUARD ×6 | `dage_mae_model`, `rank_model_dage`, `fate_prauc`, `fate_roc`, `ood_rate`, `level_shift_model` | noise | **max abs diff 0.00e+00 on every fold** ✅ |
| TARGET | `conformal_coverage` | 0.85–0.95 | 0.401 → **0.889**, ACCEPT ✅ |
| TARGET | `fate_ece` | ACCEPT + ≥40% drop | 0.281 → **0.364** ❌ **REGRESSION** |

Per §3's independence clause `q` and `sigma_scale` are adopted; only the fate calibrator changes.

**What run 2 established about coverage** (recorded, not "fixed" — it is a property, not a bug):
`q` = 33.8/34.6/36.3/34.4/34.2 on every fold where N3 sits in the pool, and **24.4** on the one
fold where it does not. **N3's error offset alone sets the interval for the whole study**, and
LOOCV removes it from its own pool — hence 0.333 there. `q/MAE` spans 0.82 → 6.43. N2's MAE is
21.79 yet all 21 of its cells fall inside q=33.76, so residuals cluster around a per-donor
**offset** rather than scattering — T7.4.3's level shift, which is Stage 2's target. The 0.889
aggregate is split conformal's **marginal** guarantee; per-fold is **conditional** coverage,
provably unachievable distribution-free (Barber, Candès, Ramdas & Tibshirani 2021).

### Why `fate_ece` regressed — four quantities, no two the same

| stage | quantity |
|---|---|
| `calibrate.py:_nll` optimised | multi-class NLL |
| `metrics.py:ece` reported | top-1 confidence ECE |
| `scorecard.py:_ece` grades | **binary ECE on `P(safe)`** |
| `res.py` + `STAGE_3` §0.1 consume | **`S` = `P(safe)`, `P_loss`** |

Plus a fit/apply mismatch: temperature is fitted on `ensemble_logits` (mean of member logits) but
applied per-member then averaged — `softmax(mean(lg)/T)` ≠ `mean(softmax(lg/T))` by Jensen.

**The plan already pointed here.** `MASTER_PLAN` §5a names the defective quantity as
"`S`, `P_loss`" and records "**YES — Platt halves it**" (T8.2); `REF_ARCHITECTURE`:23 reads
"ECE 0.28 → ~0.13 **with Platt**". `STAGE_1`'s ≲0.17 bar is derived from that Platt measurement —
while §1b.2 specified `fit_temperature`. Change A″ resolves that inconsistency in favour of the
plan's own evidence.

### The change

**Fitted on ALL held-out cells, not just the cross-donor pool.** My first cut fitted Platt on the
cross-donor pool alone (~103 cells) and would have missed the bar:

| | mean `fate_ece` | drop | |
|---|---|---|---|
| in-dist temperature (baseline) | 0.281 | — | |
| cross-donor temperature (run 2) | 0.364 | −30% | REGRESSION |
| **cross-donor Platt** (first cut) | **≈0.199** | ~29% | **misses** the 0.169 bar |
| in-dist Platt (`fate_ece_platt`) | 0.153 | +45.3% | ACCEPT |

Decomposed: the **family** change (temperature → Platt) is worth **−45%**; the **fitting-data**
change (in-distribution → cross-donor) costs **+30%**. The first cut fixed the family and kept
the data restriction that run 2 had already measured as harmful.

So the calibrator is fitted on the **union** — calib/val split **∪** cross-donor pool (~4,593
cells). Restricting to the pool means fitting 2 parameters on 103 cells while discarding 4,490.

**RETRACTED: this is NOT a departure from the cross-donor principle.** An earlier version of this
entry called it one. Checking `T8.2` in the lab notebook shows otherwise — its table is, cell for
cell, the scorecard's own columns:

| fold | T8.2 "ECE raw" | `fate_ece` | T8.2 "ECE recal" | `fate_ece_platt` |
|---|---|---|---|---|
| N3 | 0.275 | 0.275 | 0.145 | 0.145 |
| O1 | 0.316 | 0.316 | 0.147 | 0.147 |
| O2 | 0.271 | 0.271 | 0.099 | 0.099 |
| Y1 | 0.271 | 0.271 | 0.243 | 0.243 |
| Y2 | 0.270 | 0.270 | 0.132 | 0.132 |

T8.2's "recal" is **Platt fitted on the calib split**. So `STAGE_1`'s ≲0.17 bar was itself derived
from an in-distribution-fitted Platt. Holding the calibrator to a bar measured with a method we
refused to use would be incoherent; §1b.2's `fit_temperature(xstats...)` is the line that never
matched §2's own expected effect.

The principle says *calibrate on data whose error regime matches deployment*. Its premise is
measured and decisive for ΔAge (~4 yr in-distribution vs ~14 yr out-of-donor) and **not met for
fate**: discrimination is 0.929–0.940 in-distribution against **0.96–1.00 out-of-donor** (T8.1,
no degradation), and a calib-fitted Platt **halves out-of-donor ECE on 4 of 5 folds** (T8.2 — it
transfers). So the in-distribution split *qualifies* for fate, and there is 43× more of it.

**And the principle is now tested rather than assumed.** The strict pool-only Platt is fitted on
every run and reported as a diagnostic — never shipped — via `xdonor_only_platt_a/b`,
`xdonor_only_n`, `xdonor_only_safe_ece_insample` and `shipped_safe_ece_on_pool`. On the synthetic
geometry the shipped (all-data) fit scores **0.103** on the cross-donor pool against the pool-only
fit's **0.109 in-sample** — the union wins on the pool's own data even though the pool-only fit is
being graded on exactly what it was fitted to.

`fate_calib_n` in `metrics.json` records the split (`total` / `in_dist` / `xdonor`) so the
composition of the fit is auditable rather than implied.

| file | change |
|---|---|
| `src/cellfate/common/calibration.py` **(new)** | `platt_safe` / `apply_platt`. In `common` because both layers need it and **`inference` must not import `training`** — an invariant my first draft broke |
| `training/train.py` | `ensemble_probs` — the shared helper, so the calib split and the cross-donor pool cannot be computed two different ways |
| `training/calibrate.py` | `fit_platt_binary(p_safe, y_safe)` — 2-param Platt on safe-vs-rest log-loss, slope constrained **positive** so the map is rank-preserving. Same guards as `fit_temperature` (identity fallback, never-worse-than-identity). `fit_temperature` kept as fallback |
| `training/xdonor_calib.py` | `probs_mean` — the ensemble-averaged probability, byte-for-byte `Predictor`'s `pbar`, so fit and application see the same quantity. `save_xstats`/`load_xstats` persist the pool |
| `common/schemas.py` | `TemperatureParams` gains `platt_a`/`platt_b` (defaulted `None`), validated as a pair with a positive slope. **`SCHEMA_VERSION` again not bumped** |
| `training/train_model.py` | fits Platt, leaves `temperature = 1.0` (one calibrator, not two stacked), persists xstats, reports `xdonor_safe_ece_before/after` — the metric the scorecard grades |
| `inference/predictor.py` | applies Platt to `pbar`; loss/death ratio preserved so `P_loss` stays meaningful to RES |

**Persisting the pool is the enabler:** `crossdonor_stats` costs ~35 min/fold and its output was
discarded, so every calibration experiment cost another 3.5 h. Future calibrators are now a
seconds-long offline refit — with the standing rule that selection uses **that pool only**, never
the held-out folds.

### Bar unchanged

`fate_ece` must still say ACCEPT with a **≥40% drop** (0.281 → ≤0.169). Not weakened because the
specification was wrong. Guards must stay bit-identical; Platt's positive slope makes
`fate_prauc`/`fate_roc` stable -- monotone, so it never REORDERS cells. It can still MERGE
them, which a rank metric would feel; both merge paths (the EPS clamp and a float32 output cast)
were found in audit and fixed, and a test now pins the guards at slopes up to the 1e2 bound.

On synthetic data the graded metric moves the right way — binary `P(safe)` ECE **0.176 → 0.080**
on the cross-donor pool — but that is indicative only, not evidence about the real folds.

### One test I had to fix

`test_platt_recovers_a_miscaled_and_a_BIASED_p_safe` initially "sharpened" a score that was never
calibrated, so there was no correct slope to recover and it failed for the wrong reason. Rebuilt
from `y ~ Bernoulli(sigmoid(z))`, so the true inverse is known: it now asserts a ≈ 1/3 for a 3×
over-sharpening and b ≈ −1.8 for a +1.8 bias, **and** that no pure slope can fix the biased case —
which is precisely the failure a temperature cannot address.

---

## 2026-07-21 — Dress rehearsal on the real layout; two more defects found

**Status:** ✅ **RUN.** 211 tests pass. The three Stage 1 scripts were executed end-to-end
against synthetic `cellfate_loocv_*` folds built to mirror the production layout: a bulk corpus
at **94.4%** of the training split (real HFF: 99.8%) plus six donors held out one at a time.

### `verify_1a.py` — correct on the real geometry

```
6 labels -> 5 usable ;  BULK_L0=840(SKIP), DONOR_L1..L5=10 each
VERDICT: PASS -- exactly 5 usable training donors per fold (['BULK_L0'] skipped as bulk corpora)
```

### `retrain_stage1.py` — the skip fires where it matters

```
SKIPPING donor 0 -- holding it out leaves 50 of 890 training cells (5.6%, below the 50% floor)
xdonor.done  n_donors=5  n_residuals=50  residuals_per_donor={1:10, 2:10, 3:10, 4:10, 5:10}
temperature 1.498 | q 0.477 | sigma_scale 7.796
```

**Temperature came out 1.498 — above 1, i.e. SOFTENING.** Run 1 produced 0.28–0.60 (sharpening),
because the pool was 99.8% HFF. Softening is the direction theory predicts for a model that is
over-confident out-of-donor, so the fix moves this quantity the way it should.

### Defect 1 — one missing bundle destroyed the whole snapshot (`scorecard.py`)

`measure_fold` wraps the split loading in `try/except` and returns `{"_error": ...}` per fold —
but `Predictor(root)` sat **outside** that block. A single fold with a missing, incomplete or
schema-mismatched bundle raised out of `cmd_snapshot` and **discarded every fold already
measured**. A 6-fold retrain that dies at hour 3, or a deliberate partial retrain, would cost all
the surviving results. Bundle loading is now inside the same error contract.

*(This is in the user's file, changed because the fold-level `_error` contract already existed —
the call had simply landed on the wrong side of it.)*

### Defect 2 — the gate's decision table had only ever run its PASS branch

Every STOP/FAIL path in `verify_1a.py` lived inside `main()`, reachable only by constructing a
whole dataset. That is precisely how run 1 proceeded: the one branch that ever executed was the
one that said PASS. Extracted `bulk_and_usable()` and `decide_verdict()` as pure functions and
added `tests/test_verify_1a.py` — 12 tests driving **every** branch, including:

- the run-1 geometry (corpus present → PASS, and the corpus is **named**)
- `cell_line` finer-grained than donor → STOP
- too few donors surviving the skip → STOP
- folds disagreeing on donor count → STOP
- a corpus is skipped across 51%–99% dominance, not just the extreme

The last test pins the **known gap**: a donor at 49% is kept (holding it out leaves 51%, above
the floor) yet supplies ~49% of the pooled residuals, tripping neither the skip nor the >50% pool
warning. Whether 50% is the right floor is a threshold decision — the test exists so changing it
is deliberate rather than accidental.

*Writing that test also caught an error in the test itself: I first asserted 51% was not skipped,
when it is. The boundary is now asserted in both directions.*

---

## 2026-07-21 — **EXECUTED.** Python installed locally; 199 tests + 26 smoke checks pass

**Status:** ✅ **RUN, not just written.** This supersedes every "IMPLEMENTED, NEVER EXECUTED"
caveat below for the unit tests and the smoke test. The *real-data* Stage 1 run is still pending.

Installed Python 3.11.9 (winget) and a venv at `C:\cfv` — short path deliberately: torch's nested
license directories exceed Windows `MAX_PATH` from this repo's depth, and the install fails with
`WinError 206`. torch is the CPU wheel from the PyTorch index.

### What running it immediately caught — a total blocker

```
TypeError: non-default argument 'feats' follows default argument
```

`XDonorStats.residuals_per_donor` was added *before* `feats`, and a defaulted dataclass field
cannot precede a non-defaulted one. **The package did not import at all.** Every claim in the
preceding entries — reviewed three times, "lint clean", "syntax verified" — was made against code
that could not be loaded.

Fixed by moving the field last, with a comment naming the constraint.

### Then one stale test

`test_predictor_refuses_a_mode_the_bundle_was_never_calibrated_for` set `sigma_scale_mc = 1.0`
and expected a raise — the *old* value-inference contract, written before status moved to
`sigma_calibrated_modes`. Updated to the new contract, and extended with the two cases the old
form could not express: (b) a calibrated mode whose factor clamped to 1.0 must **still load**,
and (c) a legacy bundle must behave exactly as before.

### Results

```
tests/          199 passed
smoke_stage1.py  26/26 checks, 10s
```

Selected smoke output, on the run-1 geometry:

| | |
|---|---|
| bulk corpus skipped | `SKIPPING donor 0 -- leaves 96 of 216 cells (44.4%, below the 50% floor)` |
| donors rotated | 6, corpus excluded |
| residual pool | `{1:16, 2:16, 3:16, 4:16, 5:16, 6:16}` — balanced, corpus contributes **nothing** |
| per-mode factors | ensemble **4.22**, mc_dropout **2.62** — distinct, each from its own spread |
| degenerate temperature | correctly refused (T=1.0 instead of a collapse to the 0.01 bound) |
| **reproducibility** | sigma_scale, q and temperature **identical** across two runs |

That last row **measures** the claim that mc_dropout's dropout passes don't disturb training
reproducibility — previously argued from "train_member re-seeds", never tested.

> **Honest limit:** the synthetic corpus is 55.6% of the training split; the real one (HFF) is
> 99.8%. The mechanism is exercised, but at a milder ratio than production. A donor sitting just
> under the 50% floor would be neither skipped nor flagged by the >50% pool warning — a real gap
> in the threshold design, not covered by this test.

---

## 2026-07-21 — End-to-end smoke test, and the bug writing it exposed

**Status:** ✅ Written, not run. `smoke_stage1.py` at repo root, CPU, ~2 min.

**Why the existing tests could never have caught run 1's failure.** Every test fixture uses
**balanced** synthetic sources — two cell lines, equal cells each. The real dataset is one bulk
corpus (HFF, 33,613 cells) plus six tiny donors (~14 each). The bug lived entirely in that
*geometry*, so it was invisible to the suite by construction.

`smoke_stage1.py` builds a dataset with the same shape — `BULK_L0` ~300 cells, `DONOR_L0..5` ~20
each — and runs build → train → calibrate → bundle → predict, asserting every Stage 1 invariant.
It would have caught the bulk-corpus rotation, a silent fallback to in-distribution calibration,
an uncalibrated inference mode, and a lopsided residual pool. It also **measures** the claim that
the mc_dropout passes don't disturb training reproducibility — previously argued, never tested.

### The bug it exposed before it even ran — **a factor of 1.0 is ambiguous**

Tracing the script by hand, `Predictor(mode="mc_dropout")` would have **raised on a correctly
calibrated bundle**. My guard inferred calibration status from the factor's *value*:

```python
if self.sigma_scale == 1.0 and max(ens, mc) != 1.0:   # WRONG
```

But `sigma_scale_factor` is **clamped at 1.0**, so 1.0 means *either* "measured, and the spread
was already adequate" *or* "never measured". Conflating them refuses to serve a bundle whose
spread simply needed no widening — and on well-fit data that is the normal case, not an edge one.

Fixed by recording status explicitly: `ConformalParams.sigma_calibrated_modes` (defaulted to
`[]`, so legacy bundles keep their old behaviour) plus `is_calibrated_for(mode)`. The guard now
reads the list instead of guessing from a number.

This is the third bug in a row found by *constructing the adversarial case* rather than
re-reading the code — worth weighting when judging how much confidence a review pass deserves.

---

## 2026-07-21 — Code audit: three defects Stage 1b newly exposes

**Status:** ✅ Written, not run. Code only — no test was altered to accommodate any of these.

Stage 1b shrinks the calibration pool from ~4,400 in-distribution cells to **~75 cross-donor
cells** (5 Gill donors × ~15) once HFF is skipped. Several things that were safe at the old scale
are not at the new one.

### 1. `fit_temperature` could ship a maximally overconfident T — **real bug, fixed**

Temperature is **unidentifiable on single-class data**: NLL falls monotonically as T → 0, because
"always this class, with certainty" is optimal. The optimiser runs to the lower bound (`1e-2`),
and the existing *"never worse than T=1"* guard **passes** — the fit genuinely is better on that
data — so `T = 0.01` ships and every fate probability saturates.

Unreachable before: the old pool was ~4,400 HFF cells with ample class variation. Reachable now:
~75 Gill cells whose unsafe fraction ranges 0/21 to 8/19 per donor, so a pool that is nearly all
one class is a real possibility.

Fixed in `calibrate.py` (the method's own property, so it protects every caller):
`has_class_variation()` requires ≥2 classes carrying ≥1% of the mass, else return T=1.0 with a
warning. Uncalibrated beats confidently wrong.

### 2. A lopsided residual pool is invisible — **fixed (diagnostic)**

`q` is a *quantile of the pooled residuals*, so a donor owning most of the pool sets it almost
alone. That is exactly how run 1 failed (HFF: 99.8%), and the >50% bulk-corpus skip only catches
the extreme. `XDonorStats.residuals_per_donor` now records the composition, it reaches
`metrics.json`, and `crossdonor_stats` warns when any donor exceeds 50% of the pool.

### 3. `sigma_scale` is multiplicative, so it fixes magnitude but not SHAPE — **measured, not
silently fixed**

A cell the ensemble happens to agree on keeps a near-zero sigma even after a 6× scaling. RES
consumes sigma via `R_eff = max(0, −(mu + z·σ))`, so that cell is scored as if its ΔAge were
near-certain and can be **APPROVED** on that basis — while its true out-of-donor error is ~`q`.
**That is the permissive direction, the dangerous one.**

`MASTER_PLAN` §5b-bis anticipated this and offered `R_eff = max(0, −(mu + q))` as the *"cleaner"*
alternative; `STAGE_1` specified the rescaling instead. Changing RES is a scored behaviour with a
deferred verdict (Change C, Stage 4), so this is **deliberately not fixed here**. Instead
`metrics.json` now reports `xdonor_sigma_over_q_p10/p50/p90` and
`xdonor_sigma_under_half_q_frac`, so the size of the gap is measured and the choice can be made
on evidence rather than argument.

### Also

`mc_dropout_spread`'s `DataLoader` is now explicitly `shuffle=False` — the caller indexes the
result with the age mask, so a future edit flipping that default would misalign spreads with
residuals **silently**.

---

## 2026-07-21 — mc_dropout is now actually calibrated (the guard was right)

**Status:** ✅ Written, not run.

**Two wrong answers before the right one.** The `ConfigError` on `Predictor(mode="mc_dropout")`
was not a bug in the guard — it was the guard correctly reporting that **the code had never
calibrated that mode**. My first two responses both dodged that:

1. an `xfail(strict=True)` on the failing test — silencing the alarm;
2. downgrading the raise to *drop the factor and warn* — making the alarm quieter, and rewriting
   the test to assert the quieter behaviour. That is fitting the test. The justification offered
   ("mc_dropout was uncalibrated before Stage 1 too") defends a new bug with an older one, and
   contradicts `REF_ARCHITECTURE` §5: *a miscalibrated confidence is worse than no confidence.*

**The actual job the code wasn't doing:** produce a `sigma_scale` for mc_dropout. It is cheap —
the inner-LODO has already trained the members, so it is T extra forward passes on ~15 held-out
cells per fold.

| File | Change |
|---|---|
| `xdonor_calib.py` | new `mc_dropout_spread()` mirrors `Predictor._raw_batch`'s mc branch exactly (dropout-only train mode, ONE tiled forward, `std(0, unbiased=False)`); `XDonorStats` gains `sigma_pred_mc`; `sigma_scale_factor(..., mode=)` selects the matching spread |
| `schemas.py` | `ConformalParams` gains `sigma_scale_mc` (defaulted, so old bundles still load) plus `scale_for(mode)` |
| `train_model.py` | fits **both** factors from the same held-out rows; `TrainConfig.mc_dropout_T = 50` matches `Predictor`'s default; `assert_mode_matches` deleted — obsolete once every mode has its own factor |
| `predictor.py` | selects the factor for its mode; **raises** if the bundle was calibrated but not for that mode |

**The guard survives, narrowed:** it now fires only when a bundle genuinely lacks the requested
mode's factor (e.g. a run-1 bundle). It no longer fires on every Stage-1 bundle, because every
Stage-1 bundle now has both. The `xfail` is gone and
`test_mc_dropout_is_single_batched_call` is back to its original form — passing because the
underlying defect is fixed, not because the test was loosened.

New tests: both modes carry distinct, >1.0 factors end-to-end; each factor scales *its own*
spread to the same honest width; a bundle missing one mode's factor still raises.

**Also:** `retrain_stage1.py` now sets `CUBLAS_WORKSPACE_CONFIG=:4096:8` before importing torch.
Run 1 printed torch's warning that cuBLAS GEMMs are nondeterministic on CUDA ≥ 10.2 without it.
The guards came back bit-identical anyway, but that was luck — and "bit-identical" is the sharpest
evidence we have that Stage 1 leaves the model untouched.

---

## 2026-07-20 — Follow-up task: per-mode sigma_scale for mc_dropout

**Status:** ⏳ **blocked on Stage 1 score** — the xfail marker is in place (`tests/test_inference.py`).

**What:** `test_mc_dropout_is_single_batched_call` is marked `xfail` (strict) with a placeholder reason,
because mc_dropout mode now requires its own `sigma_scale` calibration. Currently only the ensemble
spread is calibrated (xdonor produces a factor ~5–6× for ensemble). The raw mc_dropout spread is a
different magnitude (T-pass jitter vs 5-member disagreement), so it needs its own inner-LODO pass to
measure and scale.

**Why now is blocked:** Implementing this edits `xdonor_calib.py` / `train_model.py` / `predictor.py`
— the exact code being measured in Stage 1. Adding the calibration mid-experiment would contaminate
the result (one change → measure, vs. two changes → whose fault?). So it's blocked until after
`scorecard.py compare baseline A_xdonor` returns a clean result, and then it becomes the next task.

**Implementation sketch:** In `train_model.py`, after the ensemble `sigma_scale` calibration, run a
*parallel* inner-LODO measuring mc_dropout spread instead, fit a separate factor, store both
`sigma_scale` and `sigma_scale_alt` in the bundle with their modes, and have `Predictor` pick the
right one. The schema change is additive (defaults to 1.0) so all existing bundles keep loading.

**Tracking:** The strict `xfail` will force removal of the marker the moment this lands and tests
start passing — it cannot be forgotten.

---

## 2026-07-20 — Tooling: JSON output + UTF-8 console fix for the Stage 1 scripts

**Status:** ⏳ **Patched; execution in progress.** The UTF-8 fix is **confirmed working** — the first
live run of `verify_1a.py` on the data machine printed the `—` in its header instead of crashing,
which is the exact code path that failed before. The `verify_1a_results.json` write has not yet been
confirmed (the run was still in its load phase when this was recorded).

**Why.** The first real execution of the Stage 1 CLIs surfaced a blocker the "never executed"
implementation could not have caught: this machine's console codepage is **cp1255 (Hebrew)**, which
cannot encode the box-drawing characters in `render_table` (or `Δ`). Every script that prints one of
those tables raised `UnicodeEncodeError` at the first table and aborted mid-run. (Found when two
copies of `verify_1a.py` ran at once; the captured crash pointed at `cp1255.py`, "position 0–63" —
the table's top border, which is entirely box-drawing.) The user also asked for `verify_1a`'s result
to be saved to a file, as JSON, rather than only printed.

| File | Change | Why |
|---|---|---|
| `verify_1a.py` | writes **`verify_1a_results.json`** — per-fold checks plus a machine-readable `verdict.status` (`PASS` / `STOP` / `FAIL` / `CANNOT_VERIFY`), assembled and saved **before** any console table; plus `sys.stdout.reconfigure(encoding="utf-8", errors="replace")` | the verdict must survive a console that cannot print it, and the run must not die on a `print` |
| `retrain_stage1.py` | same UTF-8 `reconfigure` guard (it already wrote `retrain_stage1_results.json` per fold) | a stray non-ASCII print must not kill a multi-hour training run |
| `scorecard.py` | **deliberately untouched** | it works and owns `scorecard/baseline.json`; it already writes its snapshot JSON before printing, so its data survives a console crash. It still needs `$env:PYTHONUTF8 = "1"` for console output, since its `compare` subcommand only prints |

**Operational note.** Set `$env:PYTHONUTF8 = "1"` once per PowerShell session: the two patched
scripts no longer require it, but the untouched `scorecard.py` still does for its console tables.
Result files stay JSON (not Markdown) per the user's request — `compare` reads them as JSON.

---

## 2026-07-20 — Stage 1: cross-donor calibration (Change A)

**Plan:** `plans/STAGE_1_CALIBRATION.md` · **Deviations:** `plans/STAGE_1_DEVIATIONS.md`
**Status:** ⚠️ **IMPLEMENTED, NEVER EXECUTED.** Written on a machine with no Python, no `D:`
drive and no dataset shards. Not even an import check has run. Every claim below is from reading
code, not from running it.

**Goal:** every calibration parameter was fitted on donors the model trained alongside, then
applied to a held-out donor with a completely different error regime — one architectural mistake
with four manifestations (fate ECE 0.281, conformal coverage 0.401, `sigma_age` 2.4 yr vs ~14 yr
true error, OOD AUC 0.47). Refit them on inner leave-one-donor-out statistics instead.

### Source — sub-stage 1a (donor labels)

| File | Change | Why |
|---|---|---|
| `src/cellfate/training/dataset.py` | `DONOR_I` added as a 7th tensor column, sourced from the shard's `cell_line`; `DONOR_VOCAB` + `_donor_code()` give stable integer codes; column appended in **both** return paths, including the empty-split branch | inner-LODO is impossible without donor identity in the training tensors |
| `src/cellfate/training/train.py` | two positional unpacks (`_eval_loss`, `train_member`) converted to indexed access via `X_I…AM_I` | `for x, fp, dt, yc, ya, am in dl` breaks the moment a 7th column exists |

**Indices 0–5 are unchanged**, so the first six columns bind exactly as the old positional unpack
did. The donor column is never fed to the network (`forward` takes `x, u, dose_time` only), and
adding a tensor consumes no RNG — so **1a is expected to be bit-identical**, not merely "noise".

### Source — sub-stage 1b (cross-donor calibration)

| File | Change | Why |
|---|---|---|
| `src/cellfate/training/xdonor_calib.py` **(new)** | `crossdonor_stats()` runs inner-LODO over training donors, pooling out-of-donor residuals, logits and ensemble spread; `sigma_scale_factor()` derives the multiplier that makes the spread match reality; `n_train_donors()` exposes the precondition | produces statistics from the regime deployment actually faces |
| `src/cellfate/training/train_model.py` | `temperature`, `q` and `sigma_scale` now fit on those statistics, each with a logged in-distribution fallback; `ResParams` construction moved above the calibration block; `TrainConfig` gains `xdonor_calibration` and `inference_mode`; `_report` records `xdonor_*` diagnostics and cross-donor ECE | the actual fix |
| `src/cellfate/training/conformal.py` | `fit_conformal()` accepts `sigma_scale` / `sigma_scale_mode` and passes them into `ConformalParams` | keeps construction in one place rather than mutating the object afterwards |
| `src/cellfate/common/schemas.py` | `ConformalParams` gains `sigma_scale: float = 1.0` and `sigma_scale_mode: str = "ensemble"`, both validated | `sigma_age` needs its own rescaling; `q` alone does not reach RES |
| `src/cellfate/inference/predictor.py` | reads `sigma_scale`, applies it to `sigma_age`, and raises `ConfigError` if a non-unit factor meets a mode it was not calibrated for | applying an ensemble-calibrated factor to MC-dropout spread calibrates the wrong quantity, silently |
| `src/cellfate/training/__init__.py` | exports `crossdonor_stats`, `sigma_scale_factor`, `n_train_donors`, `XDonorStats` | — |

**`SCHEMA_VERSION` deliberately NOT bumped.** Both new fields are additive and defaulted, and
`Predictor` raises on a version mismatch — a bump would make every bundle in `runs/` fail to load.

### Tests

| File | Change |
|---|---|
| `tests/test_correctness.py` | the existing loader test now also asserts 7 columns, integer dtype, length match, 2 donor codes, and that the empty-split branch returns 7 |
| `tests/test_training.py` | `_toy_dataset` grew a donor column (`n_donors` arg); **8 new tests** — sigma factor widens / never shrinks / handles empty and `z_conf=0`; `crossdonor_stats` refuses a single donor; pre-Stage-1b `conformal.json` still loads; bundle records cross-donor provenance; `sigma_scale` reaches the Predictor; the mode guard fires and does not false-positive on a unit factor |

### Supporting files

| File | Purpose |
|---|---|
| `retrain_stage1.py` **(new)** | **Required before any Stage 1 snapshot.** `scorecard.py` does not train — it loads each fold's existing `bundle/` via `Predictor(root)`, and every Stage 1 change is in the training path. Snapshotting without retraining measures the OLD bundles and shows no change, which would read as "Stage 1 did nothing" when Stage 1 never ran. This retrains the six LOOCV folds **in place**, reusing shards/scalers/splits and redoing only train → calibrate → bundle. Uses `run_multi_local.py`'s exact hyperparameters so the comparison stays one-change. Backs up each `bundle/` to `bundle_pre_stage1/` first; `--donors N2` smoke-tests one fold; `--no-xdonor` produces the 1a-only snapshot |
| `verify_1a.py` **(new)** | Answers the precondition that gates 1b: does `cell_line` distinguish donors, at **donor granularity**? Prints raw `cell_line` values, per-fold donor counts and column counts. **Expect exactly 5** in a LOOCV training split — flags both too few and too many |
| `plans/STAGE_1_DEVIATIONS.md` **(new)** | Every departure from the plan, with reasoning |
| `experiments/DELTAAGE_LAB_NOTEBOOK.md` | Appended the Stage 1 entry, pre-registered: hypothesis, predictions and decision branches written **before** any numbers exist. Also marks the boundary where the project moves from measurement to modification |

### Post-implementation audit (same day)

A full review pass over every changed file, since none of it can be executed here. It found one
real bug and three doc/consistency defects:

| Finding | Fix |
|---|---|
| **The mode guard was defeated by its own input.** `sigma_scale_mode` was written from `cfg.inference_mode` — the label the *caller declared* — while `sigma_pred` is always the spread across ensemble members. Setting `inference_mode="mc_dropout"` would have stamped an ensemble-derived factor with an mc_dropout label, and the load-time guard, finding a matching label, would have waved it through | `SIGMA_SCALE_MODE` is now a module constant, so the label always describes what was **computed**; `assert_mode_matches()` implements the plan's §3 write-time check as a real `ConfigError` and is unit-tested. Found by checking implementation against the plan line by line — the plan asked for this assert and substituting the runtime guard alone opened the hole |
| **`tests/test_training.py` asserted an invariant Stage 1b breaks.** `fit_temperature` only promises "never worse than T=1" **on its fitting split** — which is now the cross-donor pool, not calib. The in-distribution NLL is now free to rise, and *should*: baseline T=0.542 (sharpening, because the model is under-confident in-distribution) while out-of-donor needs T>1 (softening). One scalar cannot serve both | `_report` gains `xdonor_nll_before/after_temp`; the test now asserts on whichever split the temperature was actually fitted on |
| `xdonor_calib.py`'s module docstring still claimed the OOD reference is fitted on these statistics — contradicting the implementation | corrected to "three of the four, not four", with the reason |
| `DONOR_VOCAB` code values depend on first-seen order, which is safe **only** because every pooled statistic is order-invariant — an undocumented constraint a future change could break | documented as a requirement at the definition, naming the test that would catch it |
| `schemas.py` had a 101-char line | wrapped (cosmetic; `E501` is in the repo's ruff ignore list, so it was never a CI failure) |

A second pass, line by line against `STAGE_1_CALIBRATION.md`, closed three more gaps:

| §  | Gap | Fix |
|---|---|---|
| 1a.2 | the plan's snippet prints `sorted(arr.keys())`; `verify_1a.py` did not | added |
| 1a.5 | the plan prints per-donor **counts** (`torch.bincount`); `verify_1a.py` reported only the donor set. A donor with a handful of cells makes its inner-LODO fold nearly useless and the pooled calibration quietly inherits that | added, with a "thin donor" flag below 20 cells |
| 1b.2 | Edit 2's `if/elif` has **no `else`** — as written, `temperature` is unassigned when xdonor logits *and* both in-distribution splits are empty (`NameError`). The pre-Stage-1 code had that branch | kept the original `TemperatureParams(temperature=1.0)` fallback |

Also verified: every call site of the changed signatures (`fit_conformal`, `load_split_tensors`,
`_report`, `ConformalParams`) is backward-compatible; every `Predictor()` construction in the repo
uses the default `mode="ensemble"`, so the new mode guard cannot fire spuriously; `run()` stays
reproducible; and the four ruff rules that are active (`F`, `I`, `B`, `N`) are satisfied.

### Plan defects found and fixed

1. **The inner-LODO leaked.** §1b.1 passes the held-out donor as `train_ensemble`'s monitoring
   split, so each inner model would early-stop on the very donor whose residuals then fit `q`.
   Residuals would be best-case, understating exactly what Stage 1 exists to widen. **Fixed:**
   pass the outer val split with that donor removed.
2. **The OOD refit is not implementable.** §1b.2 Edit 4 pools trunk features across
   independently-seeded inner models, whose latent bases differ by arbitrary rotation, while
   `OODDetector` compares the *deployed* model's features. **Not done** — sanctioned by §3 (the
   refits are independent) and §1b.4 (disable the gate rather than chase it). **`ood_rate` should
   not move.**
3. **`cfg.inference_mode` did not exist** — the §3 assert would have raised `AttributeError`.
   Added, plus a load-time guard in `Predictor`.
4. **`ResParams` was constructed twice.** Fixed.

### Expected effects — read against these, not the scorecard arrows

| Metric | Baseline | Expected |
|---|---|---|
| `conformal_coverage` | 0.401 | **0.85–0.95** (target) |
| `fate_ece` | 0.281 | **≲0.17** (target) |
| `conformal_width` (= **2q**) | 17.72 | **~70–86** — rising is correct |
| `sigma_scale` | 1.0 | **~5–6** |
| `ood_rate` | 0.273 | **unchanged** (see defect 2) |
| `res_approvals` | 3 (oracle 0) | **0** — the predicted correct result, not a regression |
| `dage_mae_model`, `rank_model_dage`, `fate_prauc`, `fate_roc` | — | **noise** (the four guards) |

### Pre-registered rulings (2026-07-20, before the run)

The plan contradicts itself on one bar and is silent on a near-miss. Both decided in advance:

- **Coverage > 0.95 → FAIL.** §3's bar wins over §1b.4's "overshoot is expected". Overshoot is
  *predicted* to be likely: `q` is fitted on inner models trained on 4 donors and applied to a
  deployed model trained on 5, the standard pessimistic bias of cross-validation, compounded by
  N2/N3 inflating the P90. If it fails this way the response is a **new test with a new bar**
  correcting that bias — never shrinking `q` until coverage fits, which is fitting the test.
- **`fate_ece` in 0.17–0.22 → FAIL, then fix separately.** Likely fix is a Platt calibrator
  (already 0.153 on this data) rather than a single temperature scalar. §3 makes the three refits
  independent, so this does not invalidate the coverage or `sigma_scale` results.
- **Guards must be *identical*, not "noise".** The deployed ensemble trains before
  `crossdonor_stats` with the same seeds, and `set_global_seed` enables deterministic cuDNN — so
  `dage_mae_model` should read **exactly 14.291**, `rank_model_dage` **exactly 0.948**, `ood_rate`
  **exactly 0.273**. Any movement means the change reached something it must not.

### To verify

```powershell
python verify_1a.py                        # 1. gates everything: exactly 5 donors per fold?
python -m pytest tests/ -q                 # 2. 198 + 9 new
python retrain_stage1.py --donors N2       # 3. ONE fold first — confirm it runs, check the cost
python retrain_stage1.py                   # 4. all six  (~6x the usual training time)
python scorecard.py snapshot --tag A_xdonor
python scorecard.py compare baseline A_xdonor
```

**Step 3/4 are not optional.** `scorecard.py` reads each fold's existing `bundle/`; without a
retrain it measures the pre-Stage-1 model and reports no change.

---

## 2026-07-20 — Baseline analysis (no code changed)

Read `scorecard.py` and `test18_forward_gate.py` output (user-supplied, `experiments/score + test
18.docx`). Findings recorded for the project record:

- **Baseline confirms every number the plans predicted** — MAE 14.291, rank 0.948/0.955/0.686,
  ECE 0.281, coverage 0.401 (0.000 on N2/N3), OOD 0.273.
- **Test 18 returns STOP.** Part C (forward unsafe-fraction, the decisive one) is tied. Two
  supporting observations: Part B is structurally void — its swing is identical on all six folds,
  which a linear model in `[x, dt, dt²]` guarantees by construction regardless of signal — and
  Parts A and C are numerically blown up (Y1's unsafe-fraction MAE of 2.928 on a target bounded in
  [0,1]). The STOP is probably right but the null is not clean.
- **1.8 cells per timepoint.** Per-timepoint SE 12.9–15.9 yr against an 11.35 yr effect — exceeds
  it on every donor. This breaks the ±3.7–4.6 yr arithmetic in `MASTER_PLAN` §5b-ter, which
  assumed 21 cells *at one timepoint*.
- **Three bookkeeping errors in the plan docs** (details in `STAGE_1_DEVIATIONS.md` §C): the
  "±12.7 yr" figure quoted throughout is the **ridge** baseline's shift, not the model's (the
  model's is 13.12, mean −5.71); `conformal_width` is 2q, not q; the RES over-approval figure is
  3 vs 0 here, not 14 vs 11.

  > **Retracted the same day:** I initially inferred from the −5.71 mean that "part of the model's
  > shift is global, so a free global correction is available." **Wrong.** With n=6 and sd 16.39,
  > SE is 6.69 and the 95% CI is [−22.9, +11.5] — it includes zero. The mislabelling in the plan
  > is real; the inference I drew from it was not, and it is withdrawn. Recorded rather than
  > quietly deleted, because a retracted claim is part of the record.

**Not fixed in the source plan documents** — flagged for a decision, since the first would
otherwise reach the manuscript.

---

## 2026-07-31 — Stage 1.5.2 executed and closed. Answer: **the clock is NOT calibratable.**

Plan: `plans/STAGE_1_5_2_LABEL_ANCHOR.md` (§11–§16 are the results; §0–§10 are the
pre-registration, unchanged). Downloads used: GSE165177, GSE165178 (+ the held GSE165179,
GSE165176).

### The verdict

**M-2a: SPLIT ⇒ NOT CALIBRATABLE ⇒ Phase 2 does not run.** ρ_partial **+0.267** (skin & blood)
and **+0.516** (multi-tissue) against a pre-frozen bar of 0.50. §6: a criterion met on one clock
and not the other is a failure, not a pass.

**This closes the RNA-clock route.** Five repair attempts have now failed: refit, precision,
control-swap, statistical fixes, and calibration.

### §11's falsification check was missing, and is now run

§11 requires "the clocks are checked against donor chronological age before any negative verdict
is accepted." That had not been done when M-2a was recorded — so the verdict was recorded but not
accepted. Four checks, bars frozen and committed first (`5e61147`):

| | | bar | |
|---|---|---|---|
| R4 CpG coverage | 100.0% / 94.6% | ≥ 90% | ✅ |
| R1a LODO age MAE | 6.03 / 6.63 yr | ≤ 7.17 | ✅ |
| R1b intercept-free 15-yr gap | +10.39 / +6.48 | \|err\| ≤ 9.08 | ✅ |
| R1d meth↔meth ρ_partial | **+0.568** | ≥ 0.50 | ✅ |

**The verdict is accepted.** R1c also settles §9-R1 directly: non-responders drift −0.76 / −0.24 /
+2.96 / −2.56 yr against **their own day-0**, inside clock error.

### The finding that qualifies all of it

**Two methylation clocks sharing only 60 CpGs (17%) agree with each other at ρ_partial +0.568** —
clearing the same 0.50 bar by 0.068. RNA vs multi-tissue reaches **91%** of that ceiling; RNA vs
skin & blood reaches 47%. **M-2a's bars assumed ρ_true = 0.70 and nothing here reaches it,
methylation included.** So M-2c would have been meaningless even had M-2a passed — §6's gate on it
is vindicated on a ground §6 never anticipated. R1a's folds show why: two donors of identical age
53 read **44.0** and **58.5**.

### M-2b: the pooled number inverts once split by day

7/11 on both clocks — **exactly** on a bar already loosened from 8/11. Recorded as AGREE_FRAGILE.

```
day  9:  0/3 agree   RNA +38.82   METH  -2.97
day 11:  3/4 agree   RNA -38.37   METH  -3.55
day 15:  4/4 agree   RNA -51.76   METH -68.06
```

At day 15 methylation reports −68 yr and any instrument responding to reprogramming gets the sign
right. At **day 9**, the one timepoint that discriminates — methylation says nothing has happened —
the RNA clock reports **+38.82 years of ageing** and agrees **0 of 3**. Against `REV FINAL` §1's
**+36.5 yr**: the identity artefact reproduced to within 2.3 years, on the very samples the model
trains on, against paired ground truth from the same cells.

§5 pre-committed that disagreement was the live hypothesis. It technically agreed, so **the
pre-registered expectation was wrong as stated** — recorded as a miss.

### `src/` changes — gates G-a and G-b, the only ones in this stage

* **G-a** — `_control_baseline` now records per line: `n_control`, `n_cells`, `source`
  (`controls` / `self_fallback`), `unreplicated`, and the composition of the baseline **vs the
  whole line**. Persisted to `dataset_summary.json`; `verify_stage1_5.py` gains an
  unreplicated / cross-batch column.
* **G-b** — donor chronological age parsed (**both** GEO spellings), plus `batch` from the title
  suffix — the thing D1 says nothing recorded.

**On the real data these immediately print what was previously invisible:** all six Gill donors
rest on **n=1** baselines, and all six baselines are **Exp2** while every donor spans **both**
batches. That is D1 and D2, emitted by the pipeline instead of reconstructed by hand.

**Hard guard held:** ΔAge is **bit-identical** with and without the census — `np.array_equal`, in a
unit test *and* re-checked on all six real donors. The new flags are reported **beside** the Stage
1.5 verdict, never folded into it: that PASS means one specific thing and four runs are recorded
against it.

### G-c step 1 — and it refutes §0's own evidence for G-c

§0 predicted "no signature" from `diag_d2_replication`'s −0.36 yr/day, ρ −0.214. Measured on the
**actual per-cell `y_age` labels**: slope **−1.526** yr/day, ρ_timepoint **−0.905** — *stronger*
than methylation's own −0.885 / −0.842, and off from §0's figure by 4×. The two disagree because
§0 cited a **pseudobulk of absolute predicted age**, not the control-relative post-deconfounding
labels the model trains on.

**Verdict: RUN_STEP_2** (the pre-registered "ambiguous" row) — ρ passes, slope misses the band edge
by **0.084**. Leave-one-timepoint-out gives step 2 a concrete hypothesis: ρ is robust
([−0.964, −0.857]) but **dropping day 14 alone halves the slope** to −0.938, and day 14 is the last
point before the iPSC endpoint already excluded as a cell-type change.

### Four hairline margins now on record

E1b **0.009**, D2 **0.014**, M-2a **0.016**, G-c **0.084**. Not bad luck — what happens when bars
sit near the resolution of the instrument, which the ρ = 0.568 ceiling explains.

### Still open — exactly one thing

**G-c step 2**: `age_mask=True` vs `False` for HFF in one retrain, on the existing scorecard, metric
pre-registered through `audit_metrics.bar_verdict` before the run. Not done here because it needs a
rebuild and this stage's Phase 1 guarantee is `src/` untouched for every measurement. Masking leaves
the age head of order 10² labels — too few is a finding, not a failure.

### To verify

```powershell
python -m pytest tests/ -q                                  # 537 passing (was 455)
python experiments/diag_r1_anchor_reliability.py --run "D:\GSE165179" "D:\GSE165177"
python experiments/diag_m2b_contrast_agreement.py --run "D:\GSE165178" "D:\Gill"
python experiments/diag_gc_hff_signature.py --run runs/cellfate_loocv_O1
python verify_stage1_5.py "D:\GSE242423" "D:\Gill"          # now shows the G-a baseline column
```

---

## 2026-07-31 (later) — Stage 1.5.1 REV FINAL closed out, and Stage 1.5.2 re-audited

Two pieces of work: give every open question in `STAGE_1_5_1_REV_FINAL.md` an owner, and re-check
Stage 1.5.2 now that it is closed. Both produced findings. `src/` untouched by either.

### REV FINAL §6.3 is ANSWERED — the donors *are* the same people

§10.6 listed *"O1/O2 are the same physical donors"* as **❌ not verifiable**. That is true of the
**metadata** and false of the **data**: methylation carries a genotype fingerprint, and both
GSE165178 and GSE165179 are arrays.

The roster asymmetry is the control — GSE165178 has O1/O2/**Y1/Y2**, GSE165179 has O1/O2/**O3**, so
Y1 and Y2 *cannot* match and measure what a spurious match looks like:

| query (Sendai) | O1 | O2 | O3 | best | margin |
|---|---|---|---|---|---|
| **O1** | **0.9619** | 0.8416 | 0.4272 | **O1** ✅ | **0.1203** |
| **O2** | 0.7719 | **0.9755** | 0.3925 | **O2** ✅ | **0.2036** |
| Y1 *(none)* | 0.7382 | 0.6754 | 0.5897 | — | 0.0628 |
| Y2 *(none)* | 0.7033 | 0.6529 | 0.5817 | — | 0.0504 |

Both pre-registered conditions met — correct assignment **and** margin separation. The second one
matters: a panel with no identity signal gets both right **10.9%** of the time. **⇒ `SAME_DONORS`.**

**The route there is the more interesting half.** The run **aborted twice** before the assignment was
ever computed:

| | panel | cross-arm stability (bar ≥ 0.95) | |
|---|---|---|---|
| attempt 1 | top 5000 by between-donor F | 0.821 / 0.942 / 0.966 | ❌ aborted |
| attempt 2 | 419 **trimodal** genotype-shaped probes | 0.938 / 0.985 / 0.990 | ❌ aborted |

Then I audited the bar itself — **I had set 0.95 by assertion, the same §5b violation this project
has caught four times.** Simulated with array noise from GSE165179's own exp1/exp2 replicates: a
**perfect** panel scores 0.9681 and clears 0.95 **100%** of the time. **The bar was fair, so it was
not moved.**

**The diagnosis is a finding.** The panel is stable against OSKM exposure (failed arm **0.990 /
0.994 / 0.995**) and moves *only* in cells that **succeeded** — most in O1, whose two reprogrammed
samples are day 10 and day **17**, the deepest. Global demethylation during successful reprogramming
reaches even genotype-shaped CpGs, in proportion to depth. That corroborates REV FINAL §4.2 from a
direction nothing was looking in.

**Nothing depended on this, and that is now checked**: every contrast in 1.5.1 and 1.5.2 is within a
single experiment, so no result crosses the Sendai/transient boundary on a donor label.

### Every other REV FINAL question now has an owner (§11)

| answered | |
|---|---|
| §6.2 **Gill's** labels | ❌ **NO** — Stage 1.5.2 M-2a. §8.3's "the calibration route is back on the table" is superseded |
| §6.5 absolute methylation ages | **quantified**: ±7 yr donor-level (§12-R) |
| §10.7's uncaptured flake | **closed** — not reproduced in ~15 suite runs since |
| §10.6 "Gill's ~30 yr" | **won't fix** — a claim about their paper; nothing depends on it |

**Genuinely open, three items, all owned:** §5's retention and HFF's labels need **Stage 6**
(more donors / new data); HFF's `age_mask` needs **1.5.2 G-c step 2** (one retrain, no new data).

⚠️ One thing 1.5.2 changed about §5: the −6 to −9 yr retention effect is **the same size as the
±7 yr between-donor error of the instrument measuring it.** More donors help the pairing; they will
not make the instrument sharper. Stage 6 should size for that, not just for n.

### Stage 1.5.2 §17 — the re-audit found §11's per-arm *reading* was wrong

Every load-bearing number in §11–§16 re-verified against its JSON. All exact. One thing did not
survive re-reading.

§11 reported RNA↔methylation per arm and concluded the clock *"tracks in cells that are NOT
reprogramming and stops — or inverts — in exactly the cells that are."* **That table has a numerator
and no denominator.** Adding it:

| arm | n | **meth↔meth** | RNA | |
|---|---:|---:|---:|---|
| **`transient_reprogramming_intermediate`** | 11 | **+0.936** | **−0.164** | REPROG |
| `negative_control` | 12 | +0.860 | +0.399 | |
| `failing_..._intermediate` | 12 | +0.762 | +0.112 | |
| `negative_control_intermediate` | 12 | +0.671 | +0.231 | ⚠️ too blunt |
| `failed_to_transiently_reprogram` | 12 | +0.566 | +0.430 | ⚠️ too blunt |
| **`transiently_reprogrammed`** | 9 | **+0.233** | +0.150 | ⚠️ too blunt |

**Only 3 of 6 arms have a reference sharp enough to arbitrate anything.** Three corrections:

1. **§11's headline is withdrawn as stated.** `failing_..._intermediate` is a **non-reprogramming**
   arm with a **sharp** reference (+0.762) where the RNA clock reads **+0.112** — 15% of ceiling. The
   failure is not confined to reprogramming cells.
2. **§11 counted an uninterpretable arm as evidence** — `transiently_reprogrammed` has the *lowest*
   ceiling of all six.
3. **The row that does hold is far stronger than §11 made it look, and §11 buried it:** where the two
   methylation clocks agree at **+0.936 — the sharpest reference in the study — the RNA clock is
   negatively correlated.** Where the ground truth is most reliable, the transcriptomic clock runs
   backwards.

**No verdict moves** — §7 was decided on ρ_partial at n=68, and every arm's n=9–12 was frozen as
UNRESOLVABLE by §6. The defect is that §11 labelled the table "descriptive" and then drew a
structural conclusion from it in the next sentence. **A caveat does not license a claim.**

It also sharpens §12-R: the pooled ceiling +0.568 is an average over a **4× range** (+0.233 to
+0.936), so **the reference's precision is confounded with the axis under study** — a third and
stronger reason M-2c would have been meaningless.

### To verify

```powershell
python -m pytest tests/ -q                                   # 564 passing (was 537)
python experiments/diag_donor_identity.py --run "D:\GSE165178" "D:\GSE165179"
python experiments/diag_m2a_per_arm_ceiling.py
```

---

## 2026-08-01 — Working tree tidied. No result, label or verdict altered.

Housekeeping only, recorded because the standing rule is that everything is recorded. **567 tests
pass and the CI lint command passes after every step below.** `src/` behaviour is untouched.

### Root: 40 files → 12

| what | why |
|---|---|
| **untracked `gene2vec_cache.txt` (55 MB)** | it is a **download cache**, not an artefact — `experiments/test_suite.py:64` re-fetches it on demand. Gitignored; the file stays on disk |
| **deleted 5 `.zip` files (49 MB)** | verified byte-for-byte that each is an **exact duplicate** of a directory that is *also* tracked (10/10, 7/7, and 13/13 files present unpacked). `*.zip` added to `.gitignore` so it cannot recur |
| **moved 19 `*_results.json` → `results/`** | every one was referenced by code, so this was a repoint, not a move: 18 writers now resolve a `_RESULTS` constant, 6 test files and 2 cross-reading scripts follow |
| **deleted `demo.ipynb`** | superseded |

**The check that matters for the results move:** `pytest -rs` reported **no skips**. Those tests
`pytest.skip` when their results file is missing, so "no skips" is what proves they found the new
location rather than passing vacuously.

### Scripts sorted by what they are

| moved to | files |
|---|---|
| **`experiments/`** | `test18_forward_gate.py` — an exploratory numbered test, joining `test5_ridge_gap.py` and the rest |
| **`plan_tests/` (new)** | `verify_1a.py`, `verify_stage1_5.py`, `smoke_stage1.py` — the **per-stage verification gates**, i.e. the scripts a plan names as the thing that decides whether that stage passed. With `HOW_TO_RUN.md` |
| **stayed at root** | `scorecard.py`, `retrain_stage1.py`, `audit_metrics.py` (imported by 6 files), and the three `diag_*`/`dump_*` diagnostics `tests/` imports |

**Four breakages the move caused, found before shipping rather than after:**

1. `verify_1a.py` and `verify_stage1_5.py` both resolved `results/` as `__file__.parent`, which after
   the move pointed at `plan_tests/results`. → `parents[1]`.
2. `verify_stage1_5.py` resolved `local_runners/` the same way. → `parents[1]`.
3. **`tests/test_harmonize.py` and `tests/test_verify_1a.py` load these scripts by PATH** via
   `spec_from_file_location`, so a grep for `import X` missed them entirely. This turned the suite
   red mid-way and is the reason the tests were run *before* committing.
4. `tests/test_baseline_census.py` imports `verify_stage1_5` by name → `plan_tests/` added to its
   `sys.path`.

### Stage 1.5.1 drafts ARCHIVED, not deleted

The five superseded drafts moved to `plans/archive/`. **They were not deleted, and that is
deliberate:** `STAGE_1_5_1_REV_FINAL.md` §10.7 records as a *verified check* that they are
"byte-unmodified", which only means something if they are readable next to the final document — and
they are cited by **nine** other files, `STAGE_1_5_1_REVISED.md` twelve times.

All five SHA-256 hashes verified identical across the move; `git mv` used so history follows.
`plans/archive/README.md` explains what each draft was and what superseded it.

### CI lint scope widened

`plan_tests/` added to `ruff check src/ tests/ scripts/`. Moving `verify_1a.py` into a linted
directory surfaced **two pre-existing errors** (`F841` dead local, `UP017` `timezone.utc`) — it had
never been linted, because root was never in scope. Both fixed.

**`experiments/` deliberately left OUT of scope:** it carries 12 pre-existing errors in older
scripts, and cleaning those is its own change with its own diff, not something to smuggle into a
tidy-up.

### Documentation follow-through

* `plans/00_START_HERE.md` gains a **"where things live"** map, and its two runnable `test18`
  references now point at `experiments/`.
* `plans/STAGE_1_5_3_EXECUTE.md`: the lint command in PART E and §6 widened to include
  `plan_tests/`, and the step-1 guard script `verify_age_mask_identical.py` reassigned from
  `experiments/` to `plan_tests/` — it is a per-stage gate, which is what that folder is for.
* **Historical command lines in `CHANGES.md` and the lab notebook are left exactly as written.**
  They record what was actually run at the time, which is the point of them.

### One correction this surfaced, unrelated to the tidy-up

The review commits earlier the same day added 8 lines to `build_dataset.py` around line 313. That
shifted every citation below it by **+6**, and `STAGE_1_5_3_EXECUTE.md` cited the cell-cycle
deconfounder block four times. Corrected: `445-451` → **`449-457`**, `456` → **`462`**,
`457-460` → **`463-466`**. All 38 of that document's `src/` citations were then re-verified against
the files by content, not just by range.

### 🔴 A bug the tidy-up itself introduced, found and fixed the same day

The results-file move repointed 18 writers with a regex, `Path("x.json")` -> `_RESULTS / "x.json"`.
**That was wrong in 20 places across 16 files**, because `.` binds tighter than `/` in Python:

```python
_RESULTS / "x.json".write_text(...)      # calls .write_text on the STRING -> AttributeError
(_RESULTS / "x.json").write_text(...)    # correct
```

**No existing test could catch it.** The unit tests exercise the pure functions and read the
recorded JSON; none of them calls `main()`, so all 567 passed against code that could not write its
own output. It surfaced only when a writer was actually executed as part of the pre-flight check
for Stage 1.5.3.

All 20 fixed. **`tests/test_results_paths.py` added** so the class of bug cannot pass again: it
statically checks every writer for the missing parentheses, for bare CWD-relative
`Path("x_results.json")`, for a `_RESULTS` constant that is `__file__`-relative at the right depth,
and that no `*_results.json` is left in the repo root. Verified the guard works by reintroducing the
bug and watching it fail.

**Verified afterwards by running writers end to end**, including from a different working
directory, to confirm the paths are `__file__`-relative in fact and not just in intent. The two
regenerated artefacts were then **restored to their committed versions**: they differed only in
`utc` and in `set`-iteration order, and `STAGE_1_5_2_LABEL_ANCHOR.md` §14 cites
`13:11:39` as evidence that the bar was frozen 42 minutes *before* M-2b ran. Overwriting that
timestamp would have destroyed the provenance it proves.

### Also corrected in the same pass

* Stale usage strings inside the moved scripts — `python verify_stage1_5.py` etc. still printed the
  old path in their `--help` text and in the "source data not found" message a user actually sees.
* `STAGE_1_5_3_EXECUTE.md`: the lint command widened to include `plan_tests/`, and the step-1 guard
  script `verify_age_mask_identical.py` reassigned from `experiments/` to `plan_tests/` — it is a
  per-stage gate, which is what that folder is for.

### Pre-flight sweep for Stage 1.5.3 — dangling references

Swept every `python <path>.py` command in every markdown file against the filesystem. Seven were
stale after the reorganisation, all in operator-facing DO plans, all repointed:
`STAGE_3_TOOL.md` ×1, `STAGE_6_NEW_DATA.md` ×3, `00_START_HERE.md` ×2, `REF_DATA_STRATEGY.md` ×1.

Four references remain to files that do not exist, and **all four are correct as written**:

| reference | why it is fine |
|---|---|
| `validate_stopping.py`, `test19_second_clock.py` | Stage 4/5 scripts those stages specify but nobody has written. Now marked ⚠️ in `00_START_HERE.md`'s command table so an operator is not surprised |
| `experiments/diag_label_anchor.py` | the name §10 planned; §16.5 already records that the stage shipped five differently-named scripts instead |
| `plan_tests/verify_age_mask_identical.py` | Stage 1.5.3's step-1 guard. The plan says explicitly that writing it *is* step 1 |

Historical command lines in `CHANGES.md`, `DELTAAGE_LAB_NOTEBOOK.md` and `STAGE_1_DEVIATIONS.md`
were deliberately **not** touched — they record what was actually run at the time.

### ✅ The `test_evaluation` order dependence — FIXED, not just disclosed

Earlier the same day this was reopened, characterised, given an owner, and judged not to block
Stage 1.5.3. All of that was true, and **"does not block" is not "no issue"** — the fix was a few
lines, so it was done.

**Root cause:** `evaluate()` ran inside `test_evaluate_writes_reports_and_wellformed_gates`, and
**three** tests read the `reports/cell_line.json` it produced. Two of them only worked if pytest
happened to run the writer first.

**Fix:** report generation extracted into a module-scoped fixture `eval_reports` returning
`(reports_dir, gates)`. Tests now depend on the fixture rather than on each other. **No assertion
changed**; the reports are still built exactly once per module.

| check | before | after |
|---|---|---|
| the 3 tests run **individually** | ❌ 2 of 3 failed, deterministically | ✅ all 3 pass |
| full suite, 4 consecutive runs | 1 failure in ~5 | ✅ 645 passed, 1 skipped, ×4 |

**Not overclaimed:** the *intermittent* half has not recurred in four clean runs, which is
precisely the evidence that proved too weak when this was first closed on "~15 runs, no failures".
What is established is that its most likely amplifier is gone, and that a future recurrence would
be a real fixture/tmpdir question rather than an artefact of test ordering.

---

## 2026-08-01 — **Stage 1.5.3 steps 1–4 EXECUTED.** No label moved.

`plans/STAGE_1_5_3_EXECUTE.md` steps 1–4. **676 tests pass** (was 645), ruff clean, and the
bit-identity gate reads **max|Δ| = 0.00e+00** after every step.

### The gate came first, and it self-tests that it can fail

`plan_tests/verify_age_mask_identical.py` was written and its baseline captured **before any
`src/` edit** — the only moment that can be done honestly. It compares ΔAge and `age_mask` by
`np.array_equal` and SHA-256 of the raw float64 bytes, never a tolerance.

**A gate that cannot fail is not a gate** (the `verify_1a` lesson). So every run first injects
three faults into a copy of its own baseline and aborts unless all three are caught:

| injected fault | caught |
|---|---|
| one ULP on a single ΔAge value | ✅ |
| one flipped `age_mask` bit, ΔAge untouched | ✅ |
| a reason string appearing while the policies are off | ✅ |
| *(control)* an unchanged copy must PASS | ✅ |

**Geometry:** all six Gill donors + one 1800-cell HFF chunk = **7 chunks, 1944 cells**.

### What shipped

| step | change | gate |
|---|---|---|
| **1** | **C-6** `age_mask_reason` through `Sample`, `ManifestRow`, both parquet schemas, `assemble_samples`. **C-3** HFF stamps `DONOR_AGE_YEARS = 0.0` + empty `batch` | IDENTICAL, 0.0 |
| **2** | **C-1** `AGE_MASKED_DATASETS` + the pure `age_label_policy()`; `delta_age` returns a 3-tuple | IDENTICAL, 0.0 |
| **3** | **C-2** `LinearClock.age_range` carried from `meta`; `DataConfig.enforce_clock_age_range = False` | IDENTICAL, 0.0 |
| **4** | **C-4 option (a)** `AgeProvenance` + two defaulted `Response` fields + a warning list; **PART B.2's 7 annotations** to 6 plans | `res.py` untouched; **zero deletions** in `plans/` |

### The blocking capability, demonstrated on real data

```
THE BLOCKING CAPABILITY -- one chunk, both datasets, same `source`:
   hff_sc     age_mask=False reason=dataset_policy
   hff_sc     age_mask=False reason=dataset_policy
   gill_bulk  age_mask=True  reason=None
   gill_bulk  age_mask=True  reason=None
```

**G-c step 2 is now runnable. It was not, before this stage** — `age_mask` keyed on `source`
alone and both reprogramming sources report `"reprogramming"`.

Also verified live: the clock now reports `age_range = (1.0, 96.0)`, and switching C-2 on masks
all 21 cells of the neonatal donor N2 with reason `donor_out_of_clock_range` — while the default
path leaves every one of them untouched.

### 🔴 A deviation from the plan, and why

**C-6 in the plan chose the STRICT migration** ("require the column, and rebuild"), reasoning that
C-1/C-2 move labels and force a rebuild anyway.

**That reasoning does not hold for steps 1–4.** Both policies ship with their flags **off**, so no
label moves and no rebuild happens. Requiring the column would break every committed shard in
`runs/` — read by `training/dataset.py`, `evaluation/data.py`, `inference/service.py` and three
runners — **for zero benefit**. The plan's own caveat says exactly this: it *"must not ship in a
release that does not already rebuild."*

So `shard_to_numpy` reads the column **tolerantly**, with the reasoning at the call site. Step 6's
rebuild is where it may be tightened.

### One assertion changed in four steps, and it is called out rather than buried

`tests/test_inference.py` asserted `(warning is not None) == (status == REJECTED_OOD)` — that
`warning` existed for exactly one reason. **C-4 deliberately adds a second**: the ΔAge label class
can be unvalidated on a perfectly in-distribution query, and `OODDetector` (a latent Mahalanobis
test) cannot express that. The biconditional became the implications that are actually true, which
is **strictly stronger** in the direction that matters — OOD must still always warn — plus the new
one. **Every other assertion in the stage is untouched**, including
`test_delta_age_masks_cancer_sources`, where only the tuple unpacking widened.

### Defects of my own, caught before commit

| | |
|---|---|
| `io.py` uses a **relative** import, so my absolute-form edit silently no-oped and `load_age_provenance` raised `NameError` in nine tests | |
| my `predictor.py` import edit broke a parenthesised multi-line import | |
| a test asserted `"not calibratable"` against a note reading `"NOT calibratable"` | |
| an empty-table fixture exercised a numpy reshape edge case instead of the tolerance it was meant to test | |
| a missing `ValidationError` import; two unsorted import blocks; one `N814` | |

### What is NOT done, and is not supposed to be

**Steps 5, 6, 7 remain open by design.** Step 5 is C-5's design plus its bar; **step 6 is G-c
step 2**, the retrain, which is the first thing in this whole stage that moves a label; step 7 is
whatever C-4 option (c) becomes at Stage 3. **`AGE_MASKED_DATASETS` is still empty and
`enforce_clock_age_range` is still `False`** — the capability exists, and using it is a separate,
pre-registered decision.

---

## 2026-08-02 — Stage 1.5.3 **step 5**: C-5's bar registered, and it overturned the recommendation

`python plan_tests/register_c5_bar.py` -> `results/register_c5_bar_results.json`. No `src/` file
touched, no label moved, no retrain. **699 tests pass**, ruff clean.

### The bar had to grade the mechanism, not the outcome

Step 5's gate is *"bar RESOLVABLE **before any retrain**"*, which rules out `dage_mae_model` -- that
needs step 6's run. So the bar measures what the mechanism delivers per optimiser update:

| | | |
|---|---|---|
| **B1** | P(update contributes **any** age gradient) | ≥ 0.95 |
| **B2** | P(that gradient uses **≥ 4 cells**) | ≥ 0.95 |

**B1 alone would have been too easy.** C-5's diagnosis is not only the 32 % empty batches, it is
also that the survivors carry *"a Huber loss over one or two cells"*. `k = 4` is the smallest value
that halves the per-update standard error against a single cell (SE ∝ 1/√m).

### The result

| candidate | mean cells/update | B1 | B2 | |
|---|---:|---:|---:|---|
| status quo (uniform shuffling) | 1.15 | 68.9 % | 2.9 % | ❌ FAIL |
| Option 3 — pin `s_age` only | 1.14 | 68.4 % | 2.8 % | ❌ FAIL |
| **Option 2 — accumulate, W = 8** | **9.13** | **100 %** | **98.2 %** | ✅ **PASS** |
| Option 1 — sampler, w = 7.1 | 7.97 | 100 % | 96.2 % | ✅ PASS |

**Resolvable:** the dense regime (today, before masking) clears both at 100 %.
**Discriminating:** the bar separates the candidates, and the script **exits non-zero** if it ever
stops doing so — a bar everything passes decides nothing.

### 🔴 The plan recommended Option 1. The measurement says Option 2.

1. **Option 2 scores higher on the harder bar** — 98.2 % vs 96.2 % on B2.
2. **Option 2 costs the fate task nothing.** Option 1 needs `w = 7.1`, oversampling the 75 age cells
   **7.0×** (0.223 % → 1.563 % of every batch, a **1.34 %** shift in the fate head's training mix).
   C-5 called that *"not free"*; this is the number, and Option 2's is zero because it changes no
   sampling at all.

Option 1's only advantage was simplicity, and it buys that by putting Stage 1's `+0.000`
bit-identical guard record at risk for no measured gain.

**Option 3 is dead, and now provably so:** it is `weight=1, accumulate=1` — *identical to the status
quo by construction*. Pinning `s_age` does nothing about occupancy.

### What the bar cannot settle

Whether the age head actually **learns** from 75 labels is `dage_mae_model` at step 6, and no
simulation answers it. The fate guards must still read "noise" there — with Option 2 there is no
resampling to disturb them, which is precisely why it is the safer choice.

Registered as 6 rows in `tests/test_bars_resolvable.py`, with 12 unit tests on the pure functions
including closed-form checks: the uniform mean reproduces C-5's 1.14, and the empty-batch rate
matches both the exact binomial `(1−p)^512` and the plan's `e^−1.14` estimate.

---

## 2026-08-02 — Stage 1.5.3 **step 5b**: deeper tests before committing to C-5's option

`python plan_tests/c5_deeper_tests.py` -> `results/c5_deeper_tests_results.json`. READ-ONLY: no
`src/` file touched, no label moved, no training. **721 tests pass** (18 new), ruff clean.

*(Correction to the step-5 entry above: its committed state is **703** passing, not 699 — the figure
was written mid-step and four more tests landed in the same commit. Left as written per the
annotate-never-rewrite rule.)*

### Why, when step 5 had already chosen

B1/B2 grade **occupancy** — does an update get an age gradient, over how many cells. Choosing a
design on that alone is choosing on the one axis that happened to get measured. Seven axes it cannot
see were tested (D1–D7), plus a fourth **hybrid** candidate so the comparison was not forced between
two extremes. **Two of the seven changed the reading, and one of my own step-5 claims was weaker than
I had stated it.**

| candidate | eff cells | dup | **grad upd** | cover | donor | **reps/ep** | fate churn |
|---|---:|---:|---:|---:|---:|---:|---:|
| status quo (shuffle) | 1.14 | 1.00 | 2 660 | 98.8 % | 1.01 | 0.99 | 0.0 % |
| Option 1 — sampler w = 7.1 | 7.59 | 1.05 | **3 900** | 99.9 % | 1.04 | **6.93** | **36.4 %** |
| **Option 2 — accumulate W = 8** | **9.07** | **1.00** | 480 | 98.5 % | **1.01** | **0.98** | **0.0 %** |
| Option 4 — hybrid w = 3, W = 3 | 9.41 | 1.07 | 1 260 | **94.8 %** | **1.06** | 2.92 | 36.1 % |

### D6 — the diagnostic that settled it: *information* vs *repetition*

A sampler weight does not create labels. **There are 75 and there will be 75.** Weight `w` runs `w`
age-epochs inside every fate-epoch: across the run the status quo and Option 2 make **59** passes
over those 75 labels (one per epoch, i.e. what "60 epochs" means), Option 4 makes 175, and **Option 1
makes 416**.

Option 1's extra gradient updates are bought entirely by re-showing the same 75 labels 7× per epoch
— 416 effective epochs over 75 examples, a memorisation regime. Worse for step 6 specifically: it
changes *three* things at once (delivery, exposure, and the fate head's sampling), so a
`dage_mae_model` move could not be attributed to the fix. **Step 6 is a diagnostic retrain whose
entire purpose is attribution**, and the one-change rule applies.

Option 2 changes **delivery only** — same labels, same one pass per epoch, same fate training set,
regrouped so no update is empty. Asserted as a test
(`test_accumulation_changes_delivery_and_not_exposure`): if it stops being true, reopen the decision.

### D7 — step 5's cost comparison would have been overstated by 47 %

The status quo does **not** get 3 900 age updates. 32 % of its batches hit the hard zero at
`models/losses.py:55-57`, so it gets **2 660**. Comparing Option 2's 480 against 3 900 inflates the
apparent cost of accumulation by nearly half.

### 🟡 D5 came out WEAKER than step 5 implied — corrected

Step 5 quantified Option 1's fate cost as a 1.34 % batch-composition shift and I expected the
bootstrap to be the larger, unmeasured cost. Per epoch it is — **36.4 %** of fate cells are missed.
But a bootstrap **re-rolls its misses every epoch**: over 60 epochs `P(a cell is never seen)` is
**1.8 × 10⁻²⁶**. Nothing is deleted. The real cost is variance — **59.3 ± 7.7 visits, CV 13 %** —
against a permutation's exact 60. A genuine cost, but *sampling variance*, not *data loss*; the 36 %
alone would have been an overclaim.

### Option 4 (the hybrid) loses on its own merits, not by exclusion

Added to break a forced choice, it is **dominated**: it still pays Option 1's full bootstrap cost
(36.1 %), still repeats labels 3×, and posts the **worst label coverage** (94.8 % — it misses 4 of
the 75 labels in an average epoch) and **worst donor balance** (1.06) of any candidate, for 1 260
updates. No axis makes it the best choice.

### 🔵 W = 8, and W = 7 rejected for a measured reason

W = 8 was chosen for comfort, not derived, so the whole range was swept. **W = 7 is the smallest that
clears B2** (95.6 % vs the 95 % bar) and buys 540 updates instead of 480 — and is still wrong, because
75 is not a constant, it is what survives C-1 masking *on this fold*:

| n_age | W = 7 | W = 8 |
|---:|---:|---:|
| 75 | 95.6 % ✅ | 98.1 % ✅ |
| 70 | **93.7 % ❌** | 97.2 % ✅ |
| 65 | 91.3 % ❌ | 95.6 % ✅ |
| 60 | 88.0 % ❌ | 93.4 % ❌ |

W = 7 falls below its own bar as soon as the label count moves at all. 12.5 % more updates is not
worth sitting 0.6 pp above the bar. W = 8 holds to n_age ≥ 65; below that C-5 needs revisiting
regardless of W — recorded as a known boundary rather than a step-6 surprise.

### ✅ Decision: **Option 2, W = 8** — confirmed on seven axes rather than one

### 🟠 Residual risk + the implementation trap, both pre-registered

**480 age updates may be too few to converge**, and no simulation can tell — that is
`dage_mae_model` at step 6. Contingency fixed in advance so it is not decided after seeing the
answer: if the age head is still underfit at the final epoch, the remedy is a higher **age learning
rate**, *not* a smaller W.

`huber_age_loss` (`src/cellfate/models/losses.py:48-58`) ends in `F.huber_loss(...)` —
**`reduction='mean'` by default**, over the valid cells *in that batch*. So averaging the per-batch
age losses over the window is **wrong**: it weights a 1-cell batch as heavily as a 9-cell batch,
which is the very defect C-5 exists to remove, moved up one level. The window's loss must be
`Σ(per-cell losses) / Σ(valid cells)`. And **the fate term must keep stepping every batch** — if both
accumulate, Option 2 has silently become "train 8× less" and its claim to cost the fate task nothing
is void.

18 unit tests in `tests/test_c5_deeper_tests.py`, graded against closed forms and constructions with
known answers — including the bootstrap spread checked against a direct 3 000-run simulation.

---

## 2026-08-02 — 🛑 Readiness audit for step 6: **NOT ready.** Two problems, both found by checking

Asked whether we were ready to run G-c step 2, I audited instead of answering. We are not. Neither
problem is in the code that shipped at steps 1–5b — both are in what step 6 would have done next.

### Problem 1 — no step actually implements C-5

The step table ran 1, 2, 3, 4, 5, 5b, 6. Step 5 is *"C-5 **design** + its bar"*; 5b chose the option.
**PART D's manifest lists `training/train.py` as a file this stage changes, but no step scheduled
that change.** `src/cellfate/training/train.py:117` is still
`train_dl = loader(train_ds, cfg.batch_size, shuffle=True)` — plain shuffling, exactly as E26
recorded it. C-5 is graded and unbuilt.

Not bookkeeping: **step 6's arm B *is* the starved regime C-5 exists to fix** — 75 labels, 1.14 per
batch, 32 % of updates a hard zero. Running step 6 as it stands would measure "do HFF's labels help?"
**confounded with** "is the age head trainable at 75 labels with the current loader?", and the
pre-registered reading *"A better ⇒ HFF's labels help, keep them"* would be wrong for a reason the
outcome table cannot express.

### Problem 2 — 🔴 a fixed W = 8 biases step 6 toward its own treatment

This one I got wrong in 5b, and it is the more serious. I pinned W = 8 by asking what the **masked**
regime needs. Step 6 runs **two** arms, and arm A is not masked:

| | age-valid cells | age cells/batch | age updates/epoch at fixed W = 8 | vs today |
|---|---:|---:|---:|---|
| **arm A** (control) | **33 688 of 33 688** | ~512 | 8 | **65 → 8, an 8× cut for no reason** |
| **arm B** (treatment) | 75 of 33 688 | 1.14 | 8 | 44 → 8, but each is usable |

Arm A has **no occupancy problem** — every batch is full. Fixed W = 8 buys it nothing and costs it 8×
its age optimisation. **The mechanism would handicap the control and help the treatment**, pushing
`dage_mae_model` toward *"B better, CI excludes 0"* — one of the three pre-registered outcomes, and
the one concluding *"99.7 % of the labels were net-negative."* A mechanism that tilts the result
toward the treatment conclusion is a validity threat, not a detail.

### The fix — one rule, not one constant

Trigger on the **accumulated age-cell count**, not a batch count: *step the age term once the window
holds ≥ k age cells, or after W_max batches, whichever comes first.*

* **arm A** — the first batch already holds ~512 ≥ k, so W = 1: **identical to today**, the control is
  left alone and `scorecard/baseline.json` stays meaningful.
* **arm B** — ~7–8 batches to reach k, so W ≈ 8: exactly the regime 5b validated.

One policy applied identically to both arms; it only *behaves* differently because the data differ,
which is what a controlled comparison is. It also satisfies B2 **by construction** rather than at
98.1 % probability, and `W_max = 8` from 5b's sensitivity table becomes the cap.

**5b's W = 8 analysis is not withdrawn** — it still fixes `W_max`, and the n_age ≥ 65 boundary still
holds. W = 8 becomes a **ceiling**, not a constant.

### New step 5c, blocking step 6

Added to the step table: implement C-5 Option 2 in `training/train.py`. Gates — `k` registered via
`bar_verdict`; **arm-A behaviour bit-identical to today**; the window loss is `Σloss/Σcells`, not a
mean of means; the fate term still steps every batch; the data-dependent stop asserted deterministic
under a fixed shuffle seed; and a test that every label is still used exactly once per epoch, so the
rule selects *windows*, not *labels*.

No `src/` file touched by this entry — it is a plan correction. Step 6 stays blocked until 5c ships.

---

## 2026-08-02 — Stage 1.5.3 **step 5c**: C-5 Option 2 implemented, and it ships **inert**

`python plan_tests/register_c5c_bar.py` -> `results/register_c5c_bar_results.json`, then the code.
**743 tests pass** (+22: 18 new, +4 auto-discovered by the `test_results_paths.py` write-path guard).
Ruff clean. **No label moved, no retrain, nothing rebuilt.**

### The bar went first — and it failed, which is why it goes first

`REF_GROUND_RULES.md` §5b: the bar is registered before the change it grades. Attempt 1 forced the
age window to close at each epoch's last batch, so every label would be consumed inside its own
epoch. It scored **93.9 %** against A2's 95 % bar and **failed**.

The bar was not lowered. Attributing the shortfall: the epoch-end window accounted for **4.44 pp** of
the 6.12 pp gap and the irreducible `W_max` limit for only **1.67 pp** — the *mechanism* was wrong,
not the bar. Letting the window **carry across the epoch boundary** removes the artificial partial
window entirely, and is *simpler code* (one fewer special case). Re-run: **98.2 %, PASS.**

| bar | what it grades | result |
|---|---|---|
| **A1** | control arm closes every window at W = 1 — an **equality**, not a rate | **1.0000** ✅ |
| **A2** | P(window holds ≥ 4 age cells) in the masked arm | **98.2 %** ✅ (bar 0.95) |
| **A3** | masked arm gets *more* age updates than the fixed W = 8 it replaced | **980 vs 480** ✅ |

### A3 was an unexpected bonus: the bias fix also doubles the age optimisation

Triggering on accumulated **cells** rather than **batches** closes a window as soon as it is worth
stepping on, so the masked arm gets **16.3 updates/epoch (980 over the run)** instead of fixed-W's
8/epoch (480) — at the same per-update quality. **That directly reduces the "480 updates may be too
few to converge" risk that 5b had to leave open**, without touching the learning rate.

### What shipped

| file | change |
|---|---|
| `models/losses.py` | `+ huber_age_window()` — one Huber over the window's cells, `Σloss/Σcells` |
| `models/__init__.py` | export it |
| `training/train.py` | `+ _AgeWindow`, and 6 lines in the batch loop |
| `training/train_model.py` | `+ age_window_k: int = 1`, `+ age_window_max_batches: int = 8` |

**`age_window_k = 1` is the default, and 1 means OFF — the pre-1.5.3 path, bit for bit.** It ships
inert on purpose: this stage's guard is that nothing moves until step 6 turns it on deliberately in
**both** arms. It also makes the rollback a one-value edit rather than a revert.

### The gate, proved rather than asserted

`test_arm_a_is_bit_identical_when_every_cell_is_age_valid` runs `train_member` twice — mechanism off,
then on — and compares **every parameter tensor** with `torch.equal`. It passes, and holds for
k ∈ {2, 4, 8, 16}.

A test that only asserts invariance can pass on a no-op, so two more sit beside it: one confirming
the mechanism **does** move a sparsely-labelled run, and — the real check — **the exact bug the
readiness audit found was re-injected** (a fixed-W window ignoring the cell count) and confirmed to
fail **both** arm-A identity tests plus the drift check, then restored.

18 tests in `tests/test_c5c_age_accumulation.py` + 5 rows in `tests/test_bars_resolvable.py`, covering
all five gates: arm-A identity, `Σloss/Σcells` (constructed so a mean-of-means gives a visibly
different answer), the fate head still stepping on a held-back batch, determinism under a fixed seed,
and windows-not-labels. One test drives the **shipped** `_AgeWindow` against the bar script's
`close_windows` over 30 random sequences, so the simulation the decision rests on cannot drift from
the code that ships.

### Still open, unchanged by this step

Whether the age head **learns** from 75 labels is `dage_mae_model` at step 6. 5c improves the odds
(980 updates, not 480) and removes a bias; it settles nothing about the outcome. Step 6 remains the
first thing that moves a label, and needs ~2× a full LOOCV run.

---

## 2026-08-02 — Review of the 5 incoming commits: 4 verified, 1 correction

Pulled `9592db4..28565b7` and checked each independently rather than accepting the claims.

### Verified by re-derivation, not taken on trust

* **Both changed results files moved only their timestamp.** Compared every numeric leaf:
  `diag_m2a_calibratability_results.json` (851 leaves) drifted at most **2.44e-15** relative, and
  `verify_rev_final_4_4_results.json` (389 leaves) by **exactly 0.00e+00**. No conclusion moved; the
  M-2a SPLIT verdict now has an independent cross-machine reproduction.
* **The `verify_rev_final_4_4.py` path bug was real, and it was mine.** The repo tidy-up moved
  `diag_methylation_anchor_results.json` into `results/` but the script still read it from the root,
  so it would have hard-errored *and* dropped a stray JSON in the root. **A fifth move-induced
  breakage from that reorg** — I had found four.
* **The `test_results_paths.py` hole was real.** My `pytest.skip("reads results but does not write
  any")` asserted something it never checked. Confirmed against the pre-fix file: it mentioned a
  results JSON, defined no `_RESULTS`, and *did* write — skipped silently. I also scanned all 50
  scripts for write idioms the new regex misses (`open(...,'w')`, `to_csv`, `np.save`, `savefig`,
  `to_parquet`): **none slip through** today.
* **Step 6's power arithmetic reproduces exactly.** Independently re-simulated: MDE multiplier
  1.0494 (vs 1.05), power 0.9338 / 0.6476 / 0.0752 at SD 2.0 / 3.0 / 13.7, FPR 0.0505 (vs 0.0508).
  Δ\* = 3.57 is 25 % of the 14.29 yr baseline mean, and the independent-arms bound 13.68 ≈ 13.7 —
  both check out. Their 7.5 % figure is the **correct-sign** definition, stricter than the plain
  "CI excludes 0" I first computed, and the right one.
* **`experiments/` lint went 11 → 2 errors**, and CI does not lint that directory anyway
  (`ruff check src/ tests/ scripts/ plan_tests/`). No regression.

### GAP 2 was a real hole in my own step 5c

5c ships inert at `age_window_k = 1`, and 1 means OFF. Confirmed directly: the step-6 command block
sets `AGE_MASKED_DATASETS` and **never sets `age_window_k`**. Run as written, both arms would have
used k = 1, arm B would be starved, and problem #1 from my own readiness audit would have returned
silently. Shipping inert was right; failing to schedule turning it on was not. Now pinned in the
step-6 gate.

### 🔵 The one correction: "SD ≤ ~1.0 yr" understated the usable SD by ~2×

`register_gc_step2_bar.py` computed its headline as `max(passing gridpoint)` over
`CANDIDATE_SDS = (0.5, 1.0, 2.0, ...)`. The grid **jumps 1.0 → 2.0 and never samples between**, so
it reported 1.0. Solving for the crossover by bisection: **1.91 yr** — independently cross-checked at
power 0.9609 / 0.9523 / 0.9428 for SD 1.85 / 1.90 / 1.95.

Every number in the sweep table is correct and I reproduced all of them. Only the *conclusion drawn
from it* was wrong — and only conservatively. **But it was decision-relevant:** an observed SD in
**(1.0, 1.91]** would have been declared INCONCLUSIVE while the run was in fact ≥95 % powered,
discarding a real result on a reporting artefact — in the step that decides whether 99.7 % of the
age labels are thrown away.

Fixed: `max_resolvable_sd()` bisects instead of reading a gridpoint, both figures go into the results
JSON, and two tests pin it — one that the solved value exceeds the gridpoint, one that independently
re-simulates the power actually delivered at the reported crossover. The plan's original sentence is
left as written with a correction box beside it, per the annotate-never-rewrite rule.

758 tests pass, ruff clean on the CI scope.

---

## 2026-08-02 — 🛑 Step 6 pre-flight: **STOPPED before running.** The run as documented would fabricate a null

Cleared to run step 6 with `age_window_k = 4` in both arms. I ran the pre-flight first and did **not**
start the retrain. Three blockers; the first is the dangerous kind — it does not fail, it returns a
plausible answer.

### B-1 🔴 The retrain path cannot see the arm change — proved, not inferred

`retrain_stage1.py` **reuses the existing shards** (its own docstring says so) and redoes only
train → calibrate → bundle. Its `retrain()` imports exactly `TrainConfig` and `train_model.run`, and
`run()` only calls `load_split_tensors`. **There is no build step on that path.** But `age_mask` is
computed at *build* time in `build_dataset.py` and written into the shards, then read back at
`training/dataset.py:57`.

Measured on `runs/cellfate_loocv_N2` (103 shards), reading `age_mask` off disk with the constant set
both ways: **127 815 / 127 815 age-valid either way — identical.**

Both arms would train on the same data. With `base_seed = 0` and deterministic algorithms the two
snapshots would differ by ~nothing, the paired CI would include 0, and the pre-registered outcome
table reads that as *"HFF's labels are not contributing → mask them anyway"* — **licensing the
discard of 99.7 % of the project's age labels on a run where the treatment was never applied.**

The step-6 bar cannot catch this: it grades the comparison, not whether the arms differ at all.

### B-2 `age_window_k` cannot be set through the documented path

`retrain_stage1.py:147` builds `TrainConfig(...)` from a fixed kwarg list with **no `age_window_k`**,
so it takes the default `1` = OFF. GAP 2 is closed in the plan and still open in the code.

### B-3 A faithful arm B is not a mask flip, and the material to build it is gone

`_deconfound_train_only` (`build_dataset.py:448`) fits the deconfounder on age-valid TRAIN cells and
rewrites `y_age` on every shard, so masking HFF moves `y_age` itself — C-5's "second consequence".
Redoing that needs the `_cc_cache` sidecars, and **all six folds hold 0 of them** (deleted at the end
of a build). There is also no `data/` directory — **no raw GEO input on this machine.** So neither
the cheap route (re-mask in place) nor the full rebuild can run here today.

### What step 6 requires

1. Restore the raw GEO inputs and rebuild **per arm** (6 folds × 2 arms). `retrain_stage1.py` is the
   wrong driver; PART E must name the rebuild path.
2. Plumb `age_window_k` through that driver and **assert** it is 4 in both arms.
3. Stop deleting `_cc_cache`, so a future arm-B build is cheap and needs no re-download.
4. **Guard B-1 directly:** step 6 must assert the two arms' age-valid label counts differ
   (≈33 688 vs ≈75 on train) *before* training. Two identical arms must fail loudly, not return null.

Nothing was run, nothing was rebuilt, no snapshot taken, no bundle touched. 758 tests still pass.

---

## 2026-08-02 — Correction: **B-3 was wrong. Step 6 is runnable; the GEO data was there all along**

I wrote *"no raw GEO input on this machine."* **That was my error** — I searched `find . -maxdepth 2`,
inside the repo, when the pipeline's own defaults are `D:\GSE242423` and `D:\Gill`, outside it.

Verified present: the GSE242423 genes file, 9 matrix + 9 barcode files, the Gill series matrix, the
Gill expression matrix, and the Fleischer clock. Nothing was missing and nothing had changed.

Also verified rather than assumed — having just proved the *opposite* for the retrain path:
**the arm switch does reach the build.** `aging.py:304` reads `C.AGE_MASKED_DATASETS` at call time
inside `delta_age`, which `build_dataset.py` calls during the build, and `sources.py:730` emits
`dataset_id="hff_sc"`, so the filter string matches. The N2 fold's manifest holds **42 481 HFF cells
vs 124 Gill donor samples**, consistent with the ~75 training labels expected after masking.

**Unchanged from the pre-flight:**
* **B-1 stands** — `retrain_stage1.py` is the wrong driver; PART E must call the rebuild path
  (`local_runners/run_loocv.py`). The proof that the retrain path cannot see the arm change holds.
* **B-2 stands and is broader** — `run_multi_local.py:189` also builds `TrainConfig` without
  `age_window_k`, so the rebuild driver would silently run at `k = 1` = OFF too.
* **B-3's substance stands** — arm B is not a mask flip, because `_deconfound_train_only` moves
  `y_age` itself. Only my "cannot be done here" conclusion was wrong; a full rebuild regenerates
  `_cc_cache`.

Remaining before the run: plumb `age_window_k = 4` and the arm switch into the rebuild driver, and
add the guard asserting the arms' age-valid counts differ before training. Cost, from
`run_loocv.py`'s own docstring: *"~6 full builds. Expect a few hours; run it overnight"* — twice,
once per arm.

---

## 2026-08-02 — Step 6 plumbing + the arm-contrast proof, and arm A launched

### The three fixes

1. **The arm switch now reaches the data.** `run_multi_local.py` gains `AGE_MASKED` /
   `AGE_WINDOW_K` / `AGE_WINDOW_MAX_BATCHES` and sets `constants.AGE_MASKED_DATASETS` **before the
   build** — where `age_label_policy` reads it, via `delta_age` at `aging.py:304`. This is precisely
   why step 6 must run through this driver and not `retrain_stage1.py`, which reuses shards.
2. **`age_window_k` is plumbed into `TrainConfig`.** The rebuild driver omitted it too, so step 6
   would have run at `k = 1` = OFF whatever the plan said.
3. **The B-1 guard**, placed *before* training so it fails before the compute is spent. Arm B must
   leave under 5 % of train cells age-valid and more than zero; arm A must leave over 50 %. Either
   way it writes `step6_arm_census.json`.

`run_loocv.py` takes `--arm A|B --age-window-k 4`. `run_step6_arm.sh` chains the snapshot onto the
run, because arm B overwrites arm A's builds (`scorecard.py:132` resolves `cellfate_loocv_<donor>`
exactly) and a forgotten snapshot costs hours of recompute.

### 🔬 The contrast proof — two scratch single-fold builds, identical geometry

Not inferred from older builds: both arms were actually built, holdout donor O1, same scratch
config (800 cells/timepoint, 2 epochs, 1 member), differing **only** in the mask.

| arm | age-valid / train cells | | `age_window_k` |
|---|---:|---:|---:|
| **A** (control) | **5 718 / 5 718** | 100.00 % | 4 |
| **B** (treatment) | **78 / 5 718** | 1.36 % | 4 |

**Identical train-cell count, 98.64 % of labels removed.** One change, and it lands. Both branches of
the B-1 guard executed and passed — the arm-A branch had never run before. Both smoke runs went
end-to-end (build → train → evaluate → bundle, exit 0).

This is the direct refutation of the failure mode the pre-flight found: through `retrain_stage1.py`
the two arms were provably identical (127 815 / 127 815 either way); through the rebuild driver they
differ by construction.

### Launched

Arm A (control) is running: 6 folds, `--age-window-k 4`, auto-snapshotting to `gc2_A_keep_hff`.
`xdonor_calibration` defaults to `True`, so each fold trains **6 ensembles** (5 inner + 1 deployed).
Expect hours per arm, twice. Scratch dirs and the root `cellfate_multi_bundle.zip` cleaned up.

Still to report when both arms land: the **observed SD and MDE alongside the effect**, per the
registered bar — and if `|effect| <= MDE`, the pre-registered reading is INCONCLUSIVE, not "the
labels make no difference".

---

## 2026-08-02 — Step 6 **arm A complete**, and bar A1 is confirmed in production

`./run_step6_arm.sh A` → `scorecard/gc2_A_keep_hff.json`. **All 6 folds rebuilt from raw GEO and
retrained**, 11:23 → 16:27 (~1 h/fold: ~25 min build + ~35 min train, 6 ensembles/fold on a GTX 1080).
No fold SKIPPED or FAILED.

### The guard fired on every fold

```
[arm ] AGE_MASKED_DATASETS = (empty -> arm A, control) | age_window_k = 4
[guard] train split: 33,688 age-valid of 33,688 cells (100.00%)
[guard] OK: arm A (control) carries its full label set
```

**33 688** — exactly the plan's E11 training-split figure, arrived at independently by a fresh build.

### 🔵 Arm A reproduces `baseline.json` to 3 decimals, on all six folds

| donor | arm A (k=4) | baseline |
|---|---:|---:|
| N2 | 21.794 | 21.79 |
| N3 | 29.695 | 29.69 |
| O1 | 5.388 | 5.39 |
| O2 | 7.535 | 7.54 |
| Y1 | 7.279 | 7.28 |
| Y2 | 14.057 | 14.06 |

Two things follow, and I checked they were not the same thing as "the rebuild silently didn't happen":

1. **Bar A1 holds in production.** 5c predicted that with every cell age-valid the window closes at
   W = 1 and `huber_age_window` reduces to `huber_age_loss`, leaving the control arm bit-identical.
   That was a simulation plus a unit test on a toy model; this is a full 6-fold LOOCV rebuild with
   the mechanism ON. **`age_window_k = 4` costs the control arm nothing**, so `baseline.json` stays a
   valid reference and the step-6 comparison is not confounded by the mechanism.
2. **Determinism holds across ~3 weeks and a full rebuild** from raw GEO.

**Verified it was a real rebuild, not reuse:** the run writes to `./cellfate_loocv_<donor>` (ROOT is
relative), all six built today 11:23–16:27, and `scorecard.resolve_root` searches `"."` first so the
snapshot measured *those* builds — not the July copies still sitting under `runs/`.

⚠️ **Housekeeping to settle after step 6:** there are now two fold sets — stale July builds in
`runs/` and today's in the repo root, with the root ones shadowing. The root is also the tidy-up
target from earlier. To be resolved once both arms are snapshotted, not mid-run.

Arm B is now running; it overwrites the root fold set, which is why arm A's snapshot was chained
onto its own run.

---

## 2026-08-02 — 🔬 **STEP 6 / G-c step 2 RAN. Primary result: INCONCLUSIVE — and the design is confounded.**

Both arms rebuilt from raw GEO and retrained, 6 folds each, `age_window_k = 4` in both.
`scorecard/gc2_A_keep_hff.json`, `scorecard/gc2_B_mask_hff.json`. The B-1 guard passed on all 12
folds: arm A **33 688 / 33 688** age-valid, arm B **75 / 33 688** — 75 exactly, the plan's E11
prediction, reached independently by a fresh build.

### The primary metric, with the SD and MDE reported alongside as required

| | |
|---|---|
| per-fold deltas (B − A) | `+21.37, −7.22, +3.63, +1.84, +5.35, −1.15` |
| **observed effect** | **+3.971 yr** (arm B worse) |
| **observed SD** | **9.599 yr** — *never measured before this run* |
| **MDE** (1.0494 × SD) | **10.074 yr** |
| 95 % CI | `[−6.102, +14.045]` — includes 0 |
| power for Δ\* = 3.57 at this SD | **11.3 %** |

**|effect| 3.97 ≤ MDE 10.07.** The pre-registered rule is explicit and was fixed before any number
existed: *a null with MDE > Δ\* is **INCONCLUSIVE** and licenses nothing.* This is **not** evidence
that HFF's labels contribute nothing, and it does **not** license discarding 99.7 % of the age labels.

The registered crossover for ≥95 % power was SD ≤ 1.91 yr. The observed SD is **5× that**. The arms
do not track each other, which is exactly the risk the bar was registered to expose — and it took
the run to measure it. N2 alone (+21.37 against a −7…+5 spread elsewhere) carries most of the SD.

### 🔴 Worse than underpowered: the comparison is CONFOUNDED. Two changes, not one.

**`ridge` regressed almost as much as the model** — and ridge never touches the trained age head:

| | arm A | arm B | diff |
|---|---:|---:|---:|
| `rank_model_dage` | 0.948 | 0.761 | **−0.186** |
| `rank_ridge_dage` | 0.955 | 0.808 | **−0.146** |
| `dage_mae_ridge` | 14.05 | 23.27 | **+9.21** |

The cause is C-5's **second consequence**, landing in full: `_deconfound_train_only` refits on
age-valid TRAIN cells, so masking HFF drops it from **33 613 single cells to 75 bulk samples** — and
the refitted transform is applied to *every* cell's `y_age`, including the held-out evaluation
targets. The coefficient does not shift, it changes character:

| fold | arm A (slope, intercept) | arm B (slope, intercept) |
|---|---|---|
| N2 | −3.93, −3.42 | **−24.20, +10.12** |
| O1 | −9.27, −10.09 | **−24.80, +6.62** |
| Y2 | −9.66, −9.92 | **−24.88, +3.06** |

Slope ~2.5–6× steeper, intercept sign-flipped in all three.

> **So arm B is not "arm A minus HFF's labels." It is a different target variable.** The run measures
> the age head having 75 labels **and** a different ΔAge definition for everything, simultaneously —
> a one-change-rule violation baked into the design, not into its execution. **The pre-registered
> outcome table cannot be applied to this result**, in any of its three branches.

### Fate guards: three of four hold

`fate_prauc` (0.992→0.981), `fate_roc` (0.983→0.961) and `fate_ece` (0.249→0.326) all read **noise**,
as required. But **`fate ECE (Platt)` REGRESSED** — 0.140 → 0.288, CI `[+0.008, +0.288]`. The plan
says a move there "is a finding to explain, not a trade-off." It is consistent with the same cause:
`y_age` moving changes which cells the calibration path sees. Recorded, not waved through.

Also regressed: `interval width` 65.9 → 91.9 `[+3.59, +48.53]`.

### What this licenses

**Nothing about discarding HFF's labels.** The honest reading is that step 6 as designed cannot
answer its question: the treatment is entangled with a refit of the ΔAge target itself. Making it
answerable needs the deconfounder held FIXED across arms — fit it once on arm A and reuse the
coefficient in arm B — so the only difference is which labels the age head sees. That is a new
change with its own bar, not a re-read of this run.

⚠️ **Gap in this record:** arm B overwrote arm A's builds, so arm A's `proliferation_coef` above is
the July `runs/` build used as a proxy (arm A reproduced baseline to 3 decimals, so it is a close
one, but it is a proxy and is labelled as such). Future arm runs should copy `scalers.json` out
before the next arm starts.

---

## 2026-08-02 — Step 6 deep analysis: **three confounds**, and the report

Full report: `results/STEP6_REPORT.md`.

### The headline, restated with everything alongside it

`dage_mae_model`, paired B − A, 6 folds: **effect +3.971 yr | SD 9.599 | MDE 10.074 |
CI [−6.102, +14.045] | power for Δ\* at this SD 11.3 %**. |effect| ≤ MDE ⇒ **INCONCLUSIVE**, per the
rule registered before any number existed. Dropping the N2 outlier does **not** rescue it
(n=5: effect +0.491, SD 4.933, MDE 6.125 — still inconclusive).

### Three confounds, not one

* **C-I — the ΔAge target moved.** The deconfounder refits from 33 613 single cells to 75 bulk
  samples; slope goes −3.93 → −24.20 (N2) with the intercept sign-flipping, and the new transform is
  applied to every cell including held-out evaluation targets. Arm B is a different target variable.
* **C-II — 🆕 the label pool's composition shifts into the clock's extrapolation zone.** N2 and N3 are
  donor age **0**, outside `fleischer_clock.json`'s `age_range = [1.0, 96.0]`. Their share of the age
  labels goes from **0.09 % in arm A to 40 % in arm B** (30 of 75). Masking HFF does not merely
  reduce labels — it up-weights out-of-clock-range donors by ~400×. Nothing registered this.
* **C-III — ridge is not a control.** `scorecard.py:95` fits it on `tr.y_age[tr.mask]`, the same
  masked labels, so it suffers both changes too. **I had framed it as isolating the target change;
  that was wrong.** What it does show: ridge degraded on all six folds, the model on four, mean
  excess −5.24 yr — the model degrades *less* than the linear baseline under identical damage, but
  the DiD CI [−12.58, +2.10] includes 0, so it is suggestive, not established.

### Also found

* **`scorecard.py`'s `level shift` row prints the mean without its sign.** A reads `5.713` but is
  **−5.713**; signed means move −5.713 → +2.267 (looks better) while magnitudes move 13.12 → 18.66
  (actually worse). A reader trusting that row would draw the opposite conclusion.
* **N2 breaks in arm B**: MAE 21.79 → 43.17, level shift 15.0 → 42.2, **conformal coverage
  1.00 → 0.095**. Stated as observation; that C-II explains it is a hypothesis this run cannot test.
* **Fate guards: 3 of 4 hold.** `fate_ece` (Platt) regressed 0.140 → 0.288 [+0.008, +0.288].

### Conclusion

The primary result licenses nothing. **The experiment as designed cannot answer its question** — the
treatment is entangled with a refit of the target and a 400× reweighting of extrapolated donors. The
*machinery* is sound: guards fired, both arms hit their predicted label counts, arm A reproduced
baseline exactly, C-5 Option 2 cost the control nothing. What failed is the comparison's validity,
and the plan predicted the dominant cause in writing before the run.

---

## 2026-08-03 — 🔬 **STEP 6 RERUN (post-C-I). C-I verified. Primary still INCONCLUSIVE — honestly this time.**

Both arms rebuilt and retrained, 6 folds each, `age_window_k = 4`, arm-suffixed roots. Full report:
`results/STEP6_REPORT.md` (rerun section).

### C-I is fixed — verified three ways, not asserted

| check | first run | rerun |
|---|---|---|
| deconfounder coef, N2 | A `−3.93,−3.42` vs B **`−24.20,+10.12`** | **identical in both arms** |
| …all six folds | wildly different | **identical in all 6** |
| `y_age` across arms (row-exact) | — | **`max|Δ| = 0.000e+00` over 7 062 rows** |
| arm A vs pre-C-I arm A | — | **`max|Δ| = 0.000e+00` on all 6 folds** |

C-I was a no-op on the control, as predicted, and removed the target drift in the treatment. The
comparison is now genuinely one-change.

### Primary — all 6 folds

| | rerun | first run |
|---|---|---|
| **effect** | **+0.661 yr** | +3.971 |
| **observed SD** | **4.808 yr** | 9.599 |
| **MDE** | **5.045 yr** | 10.074 |
| 95 % CI | **[−4.384, +5.707]** | [−6.102, +14.045] |
| power for Δ\* = 3.57 | 31.5 % | 11.3 % |

**The SD halved; the effect shrank six-fold to 4.6 % of baseline.** Most of the confounded run's
+3.97 was the confound, not the labels. **CI includes 0 and MDE > Δ\* → INCONCLUSIVE, licenses
nothing.**

### 🔬 C-II confirmed, and it is now the binding constraint

The pre-registered 4-in-range-fold secondary: effect +0.843, **SD 1.130**, MDE 1.799,
CI [−0.956, +2.642].

**Dropping N2 and N3 collapses the SD 4.808 → 1.130 — a factor of 4.3 on 2 of 6 folds.** The two
donors outside the clock's validated `age_range = [1.0, 96.0]` carry almost all the fold variance.
That was a hypothesis when the first run ended; it is now measured.

### ⚠️ …but the "we were powered" reading does NOT survive, and I am not taking it

At face value the secondary's MDE 1.799 ≤ Δ\* 3.572, which under the registered table would read
*"the labels are genuinely not contributing → mask them."* **That is not safe.** With n = 4 the SD is
itself a noisy estimate; the χ² 95 % interval on σ is **[0.640, 4.213]**, and at the upper end the
MDE is **6.704 ≫ Δ\***. Same for the primary (σ ∈ [3.001, 11.792], MDE up to 12.375).

**Neither analysis is robustly powered once σ is admitted to be an estimate.** The secondary was
pre-registered as "underpowered by construction"; it turned out *better* powered than the primary and
still not enough. It licenses nothing either.

### Guards

`fate_prauc`, `fate_roc`, `fate_ece` all noise ✅. **`fate_ece` (Platt) regressed again**
(+0.096 [+0.011, +0.182]) — it survived C-I, so it is now a **standing anomaly** needing explanation
rather than another run. `ood_flag_rate` nearly doubled (+0.243 [+0.068, +0.419]).
`dage_mae_ridge` and `interval_width` both fell back to noise (were REGRESSIONs).

**Ranking is the one consistent cost:** `rank_model_dage` −0.069 [−0.100, −0.037], with ridge at
−0.064 — two learners both ranking worse on 75 labels than on 33 688, no longer confounded.

### Where this leaves step 6

Not another 10 h of the same design. The binding limit is 6 paired folds with 2 outside the clock's
validated range, and no re-run fixes that. The options are: accept and report a bounded estimate;
fix C-II at source (a clock valid at age 0, or a 4-donor design with the power that implies); or
change the estimand to ranking, where a consistent effect does show — as a **new** pre-registration,
not a re-read of this one.

---

## 2026-08-04 — ARM D built and pre-registered: the stratified shuffle (step-6 follow-up)

The experiment arm C proposed and that Stage 6 (REV) gates its acquisition target on. It was
referenced in two plans as the pending gate but had **never been built or run** — confirmed by
search before implementing.

**What it does.** Permutes HFF's ΔAge labels **within each `(cell_line, time_h)` stratum** instead of
globally (arm C). The between-timepoint trajectory (ρ(day, ΔAge) = −0.905) survives exactly; only the
within-stratum cell-level pairing is destroyed. That separates a *day-level* effect (real
rejuvenation **or** systematic artefact — both produce the trajectory) from *within-timepoint*
cell-level signal (only real per-cell signal produces it). Stage 1.5.5 already removed identity and
sequencing depth as the within-timepoint candidates.

**Implementation** — one code path for both arms, so they cannot drift:
- `DataConfig.age_shuffle_strata: bool = False`; `ChunkAux.stratum` (`f"{cell_line}|{time_h}"`,
  defaulting to one global group so pre-arm-D callers and old sidecars get arm C's behaviour).
- `_shuffle_age_labels` now groups target slots by stratum and permutes within each group with one
  seeded generator consumed in sorted-key order. Unstratified = every slot in one group = arm C,
  byte-for-byte.
- Runners: `run_loocv.py --arm D`, `run_step6_arm.sh D <seed>`, tag `gc2_D_stratshuffle_hff_s<seed>`.

**Bar registered BEFORE the run** (`plan_tests/register_arm_d_bar.py`), with both lessons from arm C
applied: a **pre-registered outcome for D landing outside [A, C]** (arm C's table had none and C
landed 540 % of the way to B), and **"D is like A" treated as an equivalence claim** — TOST, margin
Δ_eq = 0.1856 fixed in advance, not a CI containing zero. Difference branch ("D like C") is fully
powered (P = 100 % at the A→C SD 0.1213); equivalence branch resolves for any SD(D−A) ≤ ~0.107
(solved by bisection, not read off the sweep — the gridpoint bug from the Δ\* bar); false-equivalence
0.0 %. The achieved SD(D−A) will be reported alongside the verdict.

**Validity guards** — 10 tests (`tests/test_arm_d_stratified_shuffle.py`): stratum multisets
preserved (so each day's mean ΔAge is untouched), no label crosses a stratum boundary, strata
respected across shards, singleton strata left alone, deterministic under seed, arm C path unchanged.
A contrast test confirms the global shuffle does **not** preserve stratum means.

857 tests pass, ruff clean, no label moved by this commit (arm D is inert until a run sets it).

---

## 2026-08-07 — Arm D pre-flight: a false alarm, traced to ground, and a better gate

Before launching arm D's ~5 h run I built it at smoke scale and compared against an arm-A build.
**It failed** — the per-stratum ΔAge multiset was not preserved (4.58 yr mean shift). The shuffle is
the whole experiment, so I stopped and traced it rather than launch.

**It was not arm D.** Three isolated measurements:

1. **Determinism** — two in-process arm-A builds are bit-identical (`max |Δ| = 0.0`). So the ETL is
   reproducible; the discrepancy was not run-to-run noise.
2. **The shuffle is a pure permutation** — instrumenting `_shuffle_age_labels` inside one real build:
   6938 target labels in, 6938 out, global ΔAge multiset preserved **exactly** (`max sorted diff
   0.0`), pre mean/sd == post mean/sd.
3. **Per-stratum, intrinsic** — within one build, all **9 timepoint strata** preserve their ΔAge
   multiset, so **each timepoint's mean ΔAge is preserved to machine precision** (worst
   `|Δmean| = 1.07e-14`). The trajectory is intact: +4.1 yr at day 4 → −33.2 at day 21.

**The failure was in the harness, not the code.** Comparing HFF *training-label values* across two
separate smoke builds is unreliable (the builds even differed in HFF cell count, 900 vs 6938, at the
same `MAX_CELLS`), and it is the wrong thing to check anyway: arm A and arm D differ only in HFF
*training* labels, while `rank_model_dage` — the metric — is scored on the **never-shuffled held-out
Gill donor**, which is deconfounded identically in both arms. So the cross-build comparison is both
flaky and irrelevant.

**Fix:** replaced the cross-build check with `plan_tests/armd_intrinsic_preflight.py` — a
single-build gate that captures the ΔAge pool immediately before and after the shuffle and asserts
per-stratum multiset preservation. No second build, no confound. Result recorded in
`results/armd_intrinsic_preflight_results.json`.

The unit tests (`tests/test_arm_d_stratified_shuffle.py`, 10/10) already proved the permutation logic
on synthetic data; the intrinsic pre-flight confirms the real stratum keys are built correctly from
`raw.obs` end to end. **Arm D is validated. The full run is clear to launch.**

Lesson, recorded because it will recur: **validate a data transform intrinsically within one build**,
not by diffing two builds — build-to-build value comparison at reduced scale carries confounds that
have nothing to do with the transform under test.

---

## 2026-08-07 — ARM D RESULT: HFF structure is WITHIN-timepoint, not the day trajectory

Stratified shuffle ran, 6 folds, snapshot gc2_D_stratshuffle_hff_s0. Fate guards held
(fate_prauc 0.992->0.993), so the age-ranking collapse is label-specific.

rank_model_dage means: A 0.948  C 0.577  D 0.610.
  D - A  -0.338  95%% CI [-0.551,-0.125]  EXCLUDES 0  -> D differs from A
  D - C  +0.033  95%% CI [-0.221,+0.287]  INCLUDES 0  -> D not distinguishable from C
D sits 91%% of the way from A to C.

PRE-REGISTERED OUTCOME #2 fires (D-A detectable AND D-C not detectable): the exploitable structure
is WITHIN-timepoint and cell-level, not the between-timepoint trajectory. Arm D preserved the
trajectory (rho -0.905, intact to 1e-14) and destroyed only within-timepoint pairing; ranking still
collapsed 91%% of the way to the fully-shuffled control. A day-level systematic artefact cannot
produce that and is REJECTED.

The mundane-explanation space is now largely closed: volume (arm C), identity and depth (1.5.5),
and now the day-level artefact (arm D) are all rejected. NOT established -- same limit as arm C --
is that the labels are CORRECT age: real per-cell signal vs clock noise vs an untested
within-timepoint artefact remain. Arm D narrows hard but does not separate real signal from clock
noise; that is the remaining open question (1.5.6 clock density bears on it).

Caveat recorded: 'D like C' is a wide-CI non-detection (SD 0.24), and the equivalence branch was
NOT resolvable (achieved SD(D-A) 0.203 vs the <=0.107 needed). The load-bearing facts are the two
that do not need equivalence: D differs from A decisively, and D lands 91%% of the way to C.
