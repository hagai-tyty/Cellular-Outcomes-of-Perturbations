# STAGE 1 — Cross-donor calibration (Change A)

**Implements:** `MASTER_PLAN.md` §5a, §5b, §5b-bis.
**Blocking for:** every other stage, including the tool.
**Scope:** 2 new files, 4 modified files, ~250 lines of new code, 2 sub-stages.

---

## 0. READ THIS FIRST — a prerequisite the earlier draft missed

Stage 1 requires an **inner leave-one-donor-out** over the training donors. Writing that as a loop
assumes the training data carries donor identity. **It does not:**

```python
# src/cellfate/training/dataset.py:21
X_I, FP_I, DT_I, YC_I, YA_I, AM_I = range(6)
#  X, fingerprint, dose_time, y_cls, y_age, age_mask   <-- NO donor column
```

`load_split_tensors` reads `arr["cell_id"]` only to filter, then discards it. **There is nothing to
group by.**

**So Stage 1 has two sub-stages, and 1a must ship first:**

| Sub-stage | What | Why |
|---|---|---|
| **1a** | add a donor label column to the training tensors | inner-LODO is impossible without it |
| **1b** | cross-donor calibration proper | the actual fix |

---

## 1. The defect

**Every calibration parameter is fitted on data from donors the model already trained alongside,
then applied to a held-out donor whose error regime is completely different.**

| Parameter | Fitted on | Line | Symptom out-of-donor | Test |
|---|---|---|---|---|
| `temperature` | `val` (fallback `calib`) | `train_model.py:116` | fate ECE **0.28** | T8.2 |
| `conformal.q` | `calib` residuals (MAE≈4 yr) | `:130` | coverage **0.40 vs 0.90**; **0.00** on N2/N3 | T14 |
| `ood` reference | `train` trunk features | `:134` | AUC **0.47** ≈ chance | T15 |
| `sigma_age` | ensemble spread — never calibrated | — | **2.4 yr** vs true error **14 yr** | T7.4.2 |

**One architectural mistake, four manifestations.**

---

# SUB-STAGE 1a — Donor labels in the training tensors

## 1a.1 The edit

`src/cellfate/training/dataset.py`

**Before (line 21):**
```python
X_I, FP_I, DT_I, YC_I, YA_I, AM_I = range(6)
```

**After:**
```python
X_I, FP_I, DT_I, YC_I, YA_I, AM_I, DONOR_I = range(7)
```

## 1a.2 Where the donor value comes from — verify before coding

`load_split_tensors` already reads `arr["cell_id"]`. Try these in order:

```python
donor_src = arr.get("cell_line")                 # SplitData exposes cell_line
if donor_src is None:
    donor_src = arr.get("donor")
if donor_src is None:                            # last resort: parse the id
    donor_src = np.array([str(c).split("_")[0] for c in ids])
```

**Run this first and confirm you get exactly the six expected donors:**

```python
from cellfate.common import io
from cellfate.common.io import ArtifactPaths
p = ArtifactPaths.of("cellfate_loocv_N2")
shard = sorted(p.shards_dir.glob("*.parquet"))[0]
arr = io.shard_to_numpy(io.read_shard(shard))
print("keys:", sorted(arr.keys()))
print("sample cell_id:", arr["cell_id"][:3])
```

> **If `cell_line` is absent and `cell_id` carries no donor prefix, STOP.** Do not guess. A wrong
> grouping silently produces in-distribution calibration wearing a cross-donor label — the exact
> defect this stage exists to fix, now invisible.

## 1a.3 Encoding

Donor is a string; `TensorDataset` needs numeric. Use a stable integer code:

```python
DONOR_VOCAB: dict[str, int] = {}          # module level, built once per run

def _donor_code(name: str) -> int:
    return DONOR_VOCAB.setdefault(str(name), len(DONOR_VOCAB))
```

Append as the 7th tensor:
```python
torch.from_numpy(np.concatenate(donors)).long(),
```

**And in the empty-split branch** — it builds its own tensors and is easy to miss:
```python
return TensorDataset(torch.empty(0, g), torch.empty(0, w),
                     torch.empty(0, N_DOSE_TIME), torch.empty(0, 3), z, z,
                     torch.empty(0, dtype=torch.long))     # <-- add
```

## 1a.4 Blast radius

Positional indexing (`ds.tensors[X_I]`) is safe. **Anything that unpacks all columns breaks.**

**The exact sites, already located:**

| File:line | Code | Safe? |
|---|---|---|
| `training/train.py:47` | `for x, fp, dt, *_ in DataLoader(...)` | ✅ **safe** — `*_` absorbs extra columns |
| `training/train.py:82` | `for x, fp, dt, yc, ya, am in dl:` | ❌ **BREAKS** — exactly 6 names |
| `training/train.py:107` | `for x, fp, dt, yc, ya, am in train_dl:` | ❌ **BREAKS** — exactly 6 names |

**Two edits required.** Re-run this after your change to confirm nothing new appeared:

```powershell
Select-String -Path src\cellfate\**\*.py -Pattern "for x, fp, dt, yc, ya, am in"
```

**Fix pattern — never positional unpack:**
```python
# BAD  — breaks when a column is added
x, fp, dt, yc, ya, am = batch
# GOOD — survives schema growth
x, fp, dt = batch[X_I], batch[FP_I], batch[DT_I]
yc, ya, am = batch[YC_I], batch[YA_I], batch[AM_I]
```

## 1a.5 Verification

```python
from cellfate.training.dataset import load_split_tensors, DONOR_I
ds = load_split_tensors(paths, scalers, "holdout", "train")
d = ds.tensors[DONOR_I]
assert len(ds.tensors) == 7,          "donor column not added"
assert d.dtype == torch.long,         "donor codes must be integer"
assert len(d) == len(ds.tensors[0]),  "length mismatch"
assert len(set(d.tolist())) >= 2,     "need >=2 training donors for inner-LODO"
print("donors:", sorted(set(d.tolist())), "counts:", torch.bincount(d).tolist())
```

**Expected:** 5 distinct codes in a LOOCV training split (6 donors minus the held-out one).
**If only 1 appears, inner-LODO is impossible — report it rather than working around it.**

## 1a.6 Acceptance for 1a

- [ ] `len(ds.tensors) == 7` on all four splits (train/val/calib/test)
- [ ] the empty-split branch also returns 7 tensors
- [ ] every positional unpack replaced with indexed access
- [ ] `python -m pytest tests/ -q` still passes (198 tests)
- [ ] `python scorecard.py snapshot --tag 1a_donorlabels` → compare vs baseline shows **every
      metric as `noise`**

> **1a must change no metric.** It only adds a column. If anything moves, something else broke.

---

# SUB-STAGE 1b — Cross-donor calibration

## 1b.1 New file: `src/cellfate/training/xdonor_calib.py`

```python
"""Cross-donor calibration statistics via inner leave-one-donor-out.

The bundle's calibrators are fitted on val/calib/train — all from donors the model trained
alongside. Out-of-donor they fail (coverage 0.40 vs 0.90; fate ECE 0.28; OOD AUC 0.47).
This module produces statistics from the regime deployment actually faces.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import torch
from torch.utils.data import TensorDataset

from .dataset import YC_I, YA_I, AM_I, DONOR_I


@dataclass
class XDonorStats:
    abs_residuals: np.ndarray   # (M,)  |ΔAge error| pooled over held-out donors
    logits:        np.ndarray   # (M,3) fate logits, out-of-donor
    targets:       np.ndarray   # (M,3) matching labels
    feats:         np.ndarray   # (M,D) trunk features, out-of-donor
    sigma_pred:    np.ndarray   # (M,)  the model's own sigma_age on those rows
    n_donors:      int


def _subset(ds: TensorDataset, mask: torch.Tensor) -> TensorDataset:
    return TensorDataset(*[t[mask] for t in ds.tensors])


def crossdonor_stats(train_ds, make_model, cfg, device,
                     train_fn, ensemble_logits, ensemble_age, member_outputs) -> XDonorStats:
    """For each donor d in the TRAINING set: train on the others, predict on d, keep the
    residuals / logits / features. Pool them.

    Cost: one extra training pass per training donor (5 for a 6-donor LOOCV fold)."""
    donors = train_ds.tensors[DONOR_I]
    uniq = sorted(set(donors.tolist()))
    if len(uniq) < 2:
        raise ValueError(f"inner-LODO needs >=2 training donors, found {len(uniq)}")

    res, log, tgt, fts, sig = [], [], [], [], []
    for d in uniq:
        hold = donors == d
        inner_tr, inner_te = _subset(train_ds, ~hold), _subset(train_ds, hold)
        if len(inner_te) == 0 or len(inner_tr) == 0:
            continue

        members, _ = train_fn(make_model, inner_tr, inner_te, cfg, device)

        # --- ΔAge residuals on the held-out donor (age-valid rows only) ---
        age = ensemble_age(members, inner_te, device).numpy()
        ya  = inner_te.tensors[YA_I].numpy()
        am  = inner_te.tensors[AM_I].numpy().astype(bool)
        if am.any():
            res.append(np.abs(age[am] - ya[am]))
            per = np.stack([member_outputs(m, inner_te, device)[1].numpy() for m in members])
            sig.append(per.std(axis=0)[am])          # ensemble spread, for the scale factor

        log.append(ensemble_logits(members, inner_te, device).numpy())
        tgt.append(inner_te.tensors[YC_I].numpy())
        fts.append(member_outputs(members[0], inner_te, device)[2].numpy())

    return XDonorStats(
        abs_residuals=np.concatenate(res) if res else np.array([]),
        logits=np.vstack(log) if log else np.zeros((0, 3)),
        targets=np.vstack(tgt) if tgt else np.zeros((0, 3)),
        feats=np.vstack(fts) if fts else np.zeros((0, 1)),
        sigma_pred=np.concatenate(sig) if sig else np.array([]),
        n_donors=len(uniq),
    )


def sigma_scale_factor(stats: XDonorStats, z_conf: float, level: float = 0.90) -> float:
    """Multiplier s such that  mu ± z_conf·(s·sigma_age)  attains nominal coverage.

    ⚠ MODE DEPENDENCY. `Predictor._raw_batch` produces sigma_age two different ways:
        mode="ensemble"   (DEFAULT) -> spread across ensemble members
        mode="mc_dropout"           -> spread across T dropout passes of member[0]
    `stats.sigma_pred` below is the ENSEMBLE spread, so this factor is only valid for
    mode="ensemble". If the tool is ever run with mc_dropout, the factor must be
    recomputed against MC samples or it calibrates the wrong quantity.

    Why this exists: `sigma_age` is the ENSEMBLE SPREAD (~2.4 yr) while true out-of-donor
    error is ~14 yr. R_eff = max(0, -(mu + z·sigma)) consumes sigma, NOT the conformal q —
    so refitting q alone leaves RES exactly as broken (MASTER_PLAN §5b-bis)."""
    if stats.abs_residuals.size == 0 or stats.sigma_pred.size == 0:
        return 1.0
    need = float(np.quantile(stats.abs_residuals, level))    # half-width for `level` coverage
    have = float(np.median(stats.sigma_pred)) * z_conf
    return max(1.0, need / have) if have > 0 else 1.0
```

## 1b.2 The four edits to `train_model.py`

**Edit 1 — compute the stats once, before any calibrator is fitted.**
Insert right after `members, val_losses = train_ensemble(...)` (~line 109):

```python
    # -- cross-donor calibration statistics (inner LODO over training donors) --
    from .xdonor_calib import crossdonor_stats, sigma_scale_factor
    xstats = crossdonor_stats(train_ds, make_model, cfg, device,
                              train_ensemble, ensemble_logits, ensemble_age, member_outputs)
```

**Edit 2 — temperature** (replaces ~111-120)

*Before:*
```python
    cal_ds = val_ds if len(val_ds) else calib_ds
    if len(cal_ds):
        cal_logits = ensemble_logits(members, cal_ds, device).numpy()
        cal_target = cal_ds.tensors[YC_I].numpy()
        temperature = fit_temperature(cal_logits, cal_target)
```
*After:*
```python
    if len(xstats.logits):
        temperature = fit_temperature(xstats.logits, xstats.targets)      # CROSS-DONOR
    elif len(val_ds) or len(calib_ds):                                     # fallback, logged
        cal_ds = val_ds if len(val_ds) else calib_ds
        temperature = fit_temperature(ensemble_logits(members, cal_ds, device).numpy(),
                                      cal_ds.tensors[YC_I].numpy())
        log.warning("xdonor logits empty; fell back to in-distribution temperature")
```

**Edit 3 — conformal** (replaces ~122-130)

*Before:*
```python
    if len(calib_ds):
        age_pred = ensemble_age(members, calib_ds, device).numpy()
        ...
        abs_res = np.abs(age_pred[am] - ya[am])
    conformal = fit_conformal(abs_res, cfg.conformal_levels)
```
*After:*
```python
    abs_res = xstats.abs_residuals                                        # CROSS-DONOR
    if abs_res.size == 0:
        log.warning("xdonor residuals empty; falling back to calib residuals")
        # ...existing calib block, unchanged, as the fallback...
    conformal = fit_conformal(abs_res, cfg.conformal_levels)
```

**Edit 4 — OOD + the new sigma scale** (replaces ~132-134)

*Before:*
```python
    train_feats = member_outputs(members[0], train_ds, device)[2].numpy()
    ood = fit_ood(train_feats)
```
*After:*
```python
    ood = fit_ood(xstats.feats if len(xstats.feats) else
                  member_outputs(members[0], train_ds, device)[2].numpy())   # CROSS-DONOR
    sigma_scale = sigma_scale_factor(xstats, ResParams(**cfg.res).z_conf)
```

## 1b.3 Persisting `sigma_scale` to inference

Least-invasive route — store it in the conformal artifact rather than inventing a new one:

1. add `sigma_scale: float = 1.0` to the conformal schema (**defaulted, so old bundles still load**)
2. pass it through `assemble_bundle`
3. in `Predictor.__init__`, after `self.q = conf.q[key]`:
   ```python
   self.sigma_scale = float(getattr(conf, "sigma_scale", 1.0))
   ```
4. apply where `sigma_age` is produced (`predictor.py:133`):
   ```python
   "sigma_age": age.std(0, unbiased=False).cpu().numpy() * self.sigma_scale,
   ```

**The 1.0 default is deliberate:** every existing bundle keeps working unchanged.

## 1b.4 Failure modes

| Symptom | Cause | Action |
|---|---|---|
| `ValueError: inner-LODO needs >=2 training donors` | donor label wrong, or degenerate split | fix 1a; do **not** bypass |
| coverage overshoots (>0.95) | `q` inflated by the two outlier donors (N2/N3) | expected; note it. **Do not tune `q` down** — that is fitting the test |
| fate ECE unchanged | logits pooled from too few donors | assert `xstats.n_donors == 5` |
| OOD AUC still ≈0.5 | the detector is uninformative regardless of fitting split | **disable the gate** (Stage 3d) rather than chase it |
| training time ×6 | expected — 5 inner passes plus the final | cache `xstats` per fold |
| a guard metric REGRESSES | the change touched something beyond calibration | **revert.** A bug, not a trade-off |

## 1b.5 Rollback

Every edit is additive with a fallback branch. To revert: force `xstats` empty and all four call
sites take their original in-distribution path. **Keep the fallbacks in the code — they are the
rollback mechanism.**

---

## 2. Expected effects — record before running

| Metric | Now | After | Class |
|---|---|---|---|
| conformal coverage | **0.40** | 0.85–0.95 | derived |
| conformal `q` | 8.9 yr | ~30–40 yr (uncorrected model) | derived |
| fate ECE | **0.28** | ~0.13 (bar: **≲0.17**) | **measured** (T8.2) |
| OOD AUC | **0.47** | >0.6, or the gate is dropped | measured |
| effective `sigma_age` | 2.4 yr | ~10–15 yr | derived |

**How the q range was derived** (so this is falsifiable, not a guess). For a 90% interval,
`q` ≈ the 90th percentile of |residual|. P90/mean was verified numerically:

| residual distribution | P90 / mean |
|---|---|
| normal \|error\| | **2.07** |
| exponential | **2.31** |
| **heavy-tailed mixture (our case: O1 MAE 5.4 vs N3 29.7)** | **2.67** |

Hence q ≈ 2.5–3.0 × mean error — a **range, not a point estimate**.

**Intervals getting wider is CORRECT, not a regression.** The current narrow interval is the
defect. `scorecard.py` marks `conformal_width` "lower is better" — read that row against this
plan, not the arrow.

## 3. Acceptance

```powershell
python scorecard.py snapshot --tag A_xdonor
python scorecard.py compare baseline A_xdonor
```

| Role | Metric | Bar |
|---|---|---|
| **TARGET** | `conformal_coverage` | reach 0.85–0.95 |
| **TARGET** | `fate_ece` | ACCEPT + ≥40% drop |
| **GUARD** | `fate_prauc`, `fate_roc` | must be **noise** |
| **GUARD** | `rank_model_dage`, `dage_mae_model` | must be **noise** |

**Assumption on record:** `sigma_scale` is calibrated for `mode="ensemble"` (the Predictor default). Assert it at bundle-write time: `assert cfg.inference_mode == "ensemble"`.

**Any guard in REGRESSION means the change touched something it shouldn't.** Stage 1 alters only
calibration — never discrimination, never point predictions.

The three refits are **independent**: if one fails, adopt the other two and diagnose separately.
Most likely culprit is too few donors for a stable inner-LODO pool.

## 4. Why the whole program waits on this

| Tool feature (Stage 3) | Requires |
|---|---|
| risk threshold on `P(unsafe)` | calibrated `temperature` |
| honest error bars | calibrated `q` |
| a meaningful `R_eff` | corrected `sigma_age` |
| "unlike training" warning | working OOD |

**Building Stage 3 first would ship a recommender whose stated 90% intervals cover 40%** — a
correctness failure, not a polish issue.

## 5. Done when

1a acceptance passes, 1b acceptance passes, the five expected effects are measured, and the result
is in the lab notebook **including if it fails**.
