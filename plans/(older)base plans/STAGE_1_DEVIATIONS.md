# STAGE 1 — Deviation log

**What this is.** Every place the implementation departs from `STAGE_1_CALIBRATION.md`, with the
reason and how to check it. The stage documents are **unmodified**; this file is additive, so the
original plan stays auditable next to what was actually built.

**Status.** Code complete, **nothing executed** — the machine this was written on has no Python,
no `D:` drive and no shards. Every claim below is from reading the code, not from running it.
The first run on the data machine is the real test.

---

## A. Defects found in the plan, and what was done instead

### A1. Inner-LODO leaked the held-out donor into early stopping  — **FIXED**

`S1b.1` specifies:

```python
members, _ = train_fn(make_model, inner_tr, inner_te, cfg, device)
```

The second dataset argument of `train_ensemble` is the **monitoring split**: `train_member` sets
`monitor_dl` from it and uses it for early stopping and best-checkpoint selection. Passing
`inner_te` means each inner model **selects its checkpoint on the very donor whose residuals then
fit `q` and `sigma_scale`.**

Those residuals would be best-case, not honest out-of-donor — under-stating exactly the
quantity Stage 1 exists to widen. Deployment gets no such privilege.

**Implemented instead:** the outer `val_ds` with donor `d` removed. Legitimate (val comes from
training donors, never the outer held-out one), preserves early stopping, and carries no
information about the inner held-out donor.

```python
inner_val = _subset(val_ds, val_ds.tensors[DONOR_I] != d) if len(val_ds) else val_ds
```

If `val_ds` is empty this passes it straight through, and `train_member` falls back to monitoring
on train — no leakage either way.

> **Watch for:** if coverage after Stage 1 lands well *above* 0.95, suspect this. Over-covering is
> the signature of residuals that are too large; under-covering was the signature of the leak.

### A2. Fitting OOD on cross-donor features is not implementable — **NOT DONE, deliberately**

`S1b.2` Edit 4 says to fit the OOD reference on `xstats.feats`. Those features come from **five
independently-seeded inner models**. Independently trained networks have latent bases related only
by arbitrary rotation and permutation, so their coordinates are not comparable.

`OODDetector` compares the **deployed** `members[0]`'s latent `z` against the stored mean and
precision matrix. Fitting that Gaussian on pooled, incomparable coordinates makes the Mahalanobis
distance meaningless — and silently so. That is worse than the defect it was meant to fix.

**Implemented instead:** OOD still fitted on the deployed model's train features (original
behaviour), with the reasoning in a comment at the call site.

This is sanctioned by the plan's own structure: `S3` states the refits are **independent** ("if one
fails, adopt the other two"), and `S1b.4` already predicts the outcome — *"the detector is
uninformative regardless of fitting split → **disable the gate** (Stage 3d) rather than chase
it."* T15's AUC of 0.47 is a property of the representation, not of the fitting split.

`xstats.feats` is still collected, marked diagnostic-only, in case a future redesign wants it.

**Consequence for acceptance:** `ood_rate` is expected to be **unchanged** (~0.273). Three of the
four refits ship; the OOD one is deferred to the Stage 3d disable decision.

### A3. `cfg.inference_mode` did not exist — **ADDED, plus a runtime guard**

`S3` says to `assert cfg.inference_mode == "ensemble"` at bundle-write time. `TrainConfig` had no
such field; the assert would have raised `AttributeError`.

Added `inference_mode: str = "ensemble"` to `TrainConfig`. But a train-time assert only catches
train-time misuse — the real risk is someone constructing `Predictor(..., mode="mc_dropout")`
later, against a bundle whose `sigma_scale` was calibrated on ensemble spread. So the mode is
**recorded in the bundle** (`ConformalParams.sigma_scale_mode`) and asserted at load:

```python
if self.sigma_scale != 1.0 and scale_mode != self.mode:
    raise ConfigError(...)
```

A unit factor is mode-agnostic and never trips the guard, so pre-Stage-1b bundles are unaffected.

> **The plan's assert was still needed — dropping it opened a hole.** First implementation wrote
> `sigma_scale_mode=cfg.inference_mode`, i.e. the label the *caller declared*. But `sigma_pred` is
> always collected as the spread across ENSEMBLE MEMBERS, so setting
> `inference_mode="mc_dropout"` would have stamped an ensemble-derived factor with an
> `mc_dropout` label — and the load-time guard, finding a matching label, would have waved it
> through. The runtime guard would have been defeated by its own input.
>
> **Fixed:** `SIGMA_SCALE_MODE = "ensemble"` is a module constant in `xdonor_calib`, not a config
> value; the label written always describes what was computed. `assert_mode_matches()` implements
> §3's check as a real `ConfigError`, and is unit-tested. Both halves are now needed: the
> write-time check keeps the label honest, the load-time check keeps the consumer honest.

### A4. `ResParams` was constructed twice — **FIXED**

Edit 4 builds `ResParams(**cfg.res)` for `z_conf` at ~line 134; `train_model.py:136` already built
it. Moved the original construction above the calibration block and reused it.

### A5. Stage 1b breaks an invariant an existing test asserts — **FOUND IN AUDIT, FIXED**

`fit_temperature` guarantees "never worse than T=1" **on the data it was fitted on**.
`tests/test_training.py` asserted that as `nll_after_temp <= nll_before_temp`, where both are
computed on the **calib/val** split — which used to be the fitting split and no longer is.

The invariant is now legitimately false there, and for a reason worth stating: baseline
`temperature` is **0.542**, i.e. below 1, which *sharpens*. The model is **under-confident
in-distribution** and **over-confident out-of-donor**, so the cross-donor fit will push T above 1
to soften. One scalar cannot serve both regimes — in-distribution NLL rising is the price of the
trade, not a bug.

**Fixed** by reporting the guarantee where it now holds: `_report` gains
`xdonor_nll_before_temp` / `xdonor_nll_after_temp`, and the test asserts on whichever split the
temperature was actually fitted on. The in-distribution numbers are kept as the *contrast*, which
is what Stage 1 is measuring.

> **If `nll_after_temp > nll_before_temp` in the run, that is expected.** The number to watch is
> `xdonor_ece_after_temp` — that is the one Stage 1 is trying to move, from 0.281 toward ~0.13.

---

## B. Deliberate design choices

| # | Plan | Implemented | Why |
|---|---|---|---|
| B1 | `crossdonor_stats(..., train_fn, ensemble_logits, ensemble_age, member_outputs)` — functions injected | direct imports from `.train` | No circular import exists (`xdonor_calib → train → dataset`), and direct imports match how `train_model` gets every other helper. Four injected callables were guarding against a cycle that isn't there. |
| B2 | `crossdonor_stats` raises on <2 donors, full stop | still raises, but `run()` checks `n_train_donors` **first** and takes the documented fallback | A single-cell-line dataset is legitimate for the library (the synthetic test sources are). Crashing `run()` on it would be wrong. The guard stays for direct callers — it is a programming error there. |
| B3 | rollback = "force `xstats` empty" by editing code | `TrainConfig.xdonor_calibration: bool = True` | `S1b.5` requires the fallbacks to *be* the rollback mechanism. A flag makes that switch explicit and testable instead of requiring a code edit under pressure. |
| B4 | `sigma_scale` added to the conformal schema | same, **plus** `sigma_scale_mode`, threaded through `fit_conformal(...)` | Keeps `ConformalParams` construction in one place rather than mutating it afterwards. |
| B5 | — | `metrics.json` records `xdonor_calibrated`, `xdonor_n_donors`, `xdonor_n_residuals`, `sigma_scale`, and cross-donor ECE before/after temperature | A silent fallback to in-distribution calibration must be **auditable after the fact**, not merely visible in a log line that scrolled past. |

### B6. `SCHEMA_VERSION` deliberately NOT bumped

`schemas.py` says changing a contract is breaking and requires a bump. But `Predictor.__init__`
**raises** on a version mismatch, so bumping would make **every bundle in `runs/` fail to load** —
directly contradicting `S1b.3`'s stated intent that "every existing bundle keeps working
unchanged."

Both new fields are additive and defaulted, so old artefacts validate unchanged. Locked in by
`test_conformal_schema_still_loads_pre_stage1b_bundles`, which loads a `conformal.json` written
without either field.

### B9. Two small differences found in a line-by-line re-read of §1b

| Plan | Implemented | Note |
|---|---|---|
| `n_donors=len(uniq)` | `n_donors=used` (folds actually run) | §1b.4's diagnostic is *"did we really pool from 5 donors?"*, which `used` answers and `len(uniq)` does not. **They can never differ in practice** — the `continue` guard is unreachable once the ≥2 check has passed — so this is intent, not behaviour |
| Edit 2's `if / elif` has **no `else`** | added the original `TemperatureParams(temperature=1.0)` branch | As written the plan leaves `temperature` unassigned when xdonor logits *and* both in-distribution splits are empty → `NameError`. The pre-Stage-1 code had that else; dropping it would have been a regression |

`sigma_scale_factor` matches the plan exactly, plus `isfinite` guards on both operands — a strict
superset of the plan's `if have > 0` check.

### B8. The donor-source fallback chain was NOT implemented

§1a.2 offers a chain: `cell_line`, else `donor`, else parse a prefix off `cell_id`. Only the first
is implemented, because:

- `shard_to_numpy` **always** returns `cell_line` (it is in `SHARD_SCHEMA` and set
  unconditionally), and `Sample` requires it non-null. The other two branches are unreachable.
- The `cell_id`-parsing branch contradicts §1a.2's own instruction two paragraphs later — *"If
  `cell_line` is absent and `cell_id` carries no donor prefix, **STOP. Do not guess.**"* A silent
  prefix-parse is exactly the guess that produces in-distribution calibration wearing a
  cross-donor label.

`verify_1a.py` surfaces the real question instead: what `cell_line` actually *contains*.

### B7. No knob to shrink the inner ensembles

Cost is ~6× training time (5 inner ensembles plus the final one, per fold). Tempting to run the
inner passes with fewer members — but `sigma_scale` calibrates the **ensemble spread**, and the
spread over *k* members is a different quantity for a different *k*. The inner ensembles must use
the deployed `ensemble_size`. Stated in the module docstring so nobody adds the knob later.

---

## C. Bookkeeping errors found in the plan docs while reading the baseline

**Not fixed in the source documents** — flagged here for a decision, since two of them would
otherwise reach the manuscript.

### C1. "±12.7 yr per-donor shift" is the **ridge** baseline's shift, not the model's

`STAGE_2_LEVEL_CORRECTION.md` §1's table matches the scorecard's `level shift (ridge)` row to
three decimals on all six donors, including the mean of +0.230. The model's own shifts are:

| | N2 | N3 | O1 | O2 | Y1 | Y2 | mean | mean abs |
|---|---|---|---|---|---|---|---|---|
| **model** | +15.03 | −28.35 | +0.64 | +6.56 | −8.13 | −20.02 | **−5.71** | **13.12** |
| ridge | +20.11 | −24.40 | +5.72 | +13.04 | −4.28 | −8.81 | +0.23 | 12.72 |

§4 of the *same document* uses the **model** numbers (O1's 0.64, N2's 15.03). So §1 and §4 quote
different models.

The ±12.7 figure propagates into `REF_ARCHITECTURE.md` §1, `STAGE_3_TOOL.md` §0.2,
`MASTER_PLAN.md`, and the Stage 5 claim list. The mislabelling is real and worth fixing.

> **CORRECTION (same day).** An earlier version of this section went further and claimed *"part of
> the model's shift IS global, so 'no global correction can fix it' is false."* **That was
> overstated and is withdrawn.** The point estimate is −5.71, but with n=6 and sd = 16.39 the
> standard error is 6.69, so the 95% CI is **[−22.9, +11.5] — it includes zero.** The model's mean
> shift is *not* distinguishable from no global bias at this sample size.
>
> What survives: §1 and §4 of `STAGE_2` quote **different models**, and the headline ±12.7 yr is
> ridge's. What does not survive: any inference that a free global correction is available. Per
> ground rule §1, that would have entered the plan as a finding when it was only a hypothesis with
> n=6 behind it — the exact failure mode this project has caught five times before.

### C2. `conformal_width` is 2q, not q

Baseline width 17.717 = 2 × 8.86. `REF_GROUND_RULES.md` §4 reads "conformal_width … 8.9 → ~35–43
yr", conflating the half-width with the width. If `q` reaches 35–43 the **width** row will read
70–86. Read Stage 1's result against `q`, not the width row, or a correct outcome will look
catastrophic.

### C3. RES over-approval is 3 vs 0 at baseline, not the 14 vs 11 quoted from T7.4.3

Direction confirmed (model over-approves), magnitude differs — a different test configuration.
`REF_GROUND_RULES.md` §4 and `STAGE_2` §3 both cite 14 vs 11.

---

## D. What was not attempted

- **Stage 2's A′ refit** (`STAGE_2` §9) — belongs to Stage 2, not here.
- **Caching `xstats` per fold** (`S1b.4`) — accepted the ~6× cost rather than adding a cache whose
  invalidation rules would need their own tests.
- **Anything to fix OOD** — see A2; the decision is Stage 3d's.

---

## E. Verification checklist

Run in order, from the repo root with the venv active:

```powershell
python verify_1a.py                                   # gates 1b: >=2 training donors?
python -m pytest tests/ -q                            # 198 + 8 new
python scorecard.py snapshot --tag 1a_donorlabels     # after 1a only, if taken separately
python scorecard.py snapshot --tag A_xdonor
python scorecard.py compare baseline A_xdonor
```

**Pre-registered bars** (`STAGE_1` §3), unchanged:

| Role | Metric | Bar |
|---|---|---|
| TARGET | `conformal_coverage` | reach 0.85–0.95 (from 0.401) |
| TARGET | `fate_ece` | ACCEPT + ≥40% drop (from 0.281, so ≲0.17) |
| GUARD | `fate_prauc`, `fate_roc` | must be **noise** |
| GUARD | `rank_model_dage`, `dage_mae_model` | must be **noise** |

**Read against this log, not the scorecard arrows:**

- `conformal_width` **rising is correct** — expect roughly 70–86 (see C2). The narrow interval is
  the defect.
- `ood_rate` **should not move** — the OOD refit was not implemented (see A2).
- `res_approvals` — more is not better; judge the over-approval gap.

Expect `sigma_scale` around **5–6** (≈14 yr of true error over ≈2.4 yr of ensemble spread at
`z_conf=1.0`). A factor of exactly 1.0 means the clamp fired and something upstream is wrong.

### The RES rows will move, and that is the predicted correct result

`sigma_scale` multiplies `sigma_age`, which `compute_res_batch` consumes via
`R_eff = max(0, −(mu + z·sigma))`. A ~5× wider sigma drives `R_eff` to zero, so expect:

| Metric | Baseline | After Stage 1 | Reading |
|---|---|---|---|
| `res_approvals` | 3 (N3: 2, Y2: 1) | **0** | over-approval gap closes to 0 — an improvement |
| `res_median` / `res_max` | 0.000 / 0.031 | **0.000 / 0.000** | expected |
| `rank_res` | 0.686 | may fall | **not a guard** |

This is `MASTER_PLAN` §5b-ter playing out exactly as written: with honest per-cell uncertainty
(~19 yr) against a real effect (~11 yr), per-cell `R_eff = 0`. **Not a regression** — it is the
finding that per-cell confident rejuvenation is arithmetically unreachable at this data scale.

`rank_res` is deliberately **not** in §3's guard list, and the RES verdict is deferred to Change C
(Stage 4) precisely because it cannot be judged until its inputs are fixed. Do not treat a fall
in `rank_res` as grounds to revert.

**The four guards that do matter** — `fate_prauc`, `fate_roc`, `rank_model_dage`,
`dage_mae_model` — must all read `noise`. Stage 1 touches calibration only: never discrimination,
never point predictions. Any of those four moving means the change reached something it shouldn't.
