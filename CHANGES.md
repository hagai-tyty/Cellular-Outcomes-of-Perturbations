# CHANGES

Running log of every modification to this repository, newest first.

**Convention.** One entry per stage or task. Every entry states **what** changed, **why**, and
**whether it has been executed**. Nothing is marked verified until it has actually run on the data
machine — "the code looks right" is not verification and is recorded as such.

Files added by the user (`scorecard.py`, `test18_forward_gate.py`, `plans/*` except the deviation
log, `experiments/score + test 18.docx`) are noted where relevant but are not entries here.

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
