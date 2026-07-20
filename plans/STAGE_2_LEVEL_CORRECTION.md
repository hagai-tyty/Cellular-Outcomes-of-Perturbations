# STAGE 2 — Per-donor level correction (Change B + A′)

**Implements:** `MASTER_PLAN.md` §5c, §5d.
**Depends on:** Stage 1 — a **coupling**, not a preference (§3).
**Blocking for:** absolute ΔAge claims only. **Not required for the tool.**
**Scope:** 1 new file, 2 modified, ~150 lines, plus a wet-lab protocol decision.

---

## 0. READ THIS FIRST — this stage has a non-code prerequisite

**Fix B requires k ≈ 3 cells per new donor with known true ΔAge.** True ΔAge means the clock was
run on **both a control and a perturbed sample** for those cells.

> **Decide whether that is experimentally acceptable BEFORE writing any code.**
> If it is not, **skip Stage 2 entirely.** The tool still works — you simply cannot report
> absolute ΔAge, only within-donor comparisons. That is a supported, honest product.

Everything below assumes the answer is yes.

---

## 1. The defect

The model is **unbiased in-distribution** (calib offset **−0.03 yr**) and shifted only
out-of-donor, by a **donor-specific amount that cancels on average**:

| donor | med PRED − med TRUE |
|---|---|
| N2 | **+20.11** |
| N3 | **−24.40** |
| O1 | +5.71 |
| O2 | +13.04 |
| Y1 | −4.27 |
| Y2 | −8.80 |

mean **+0.23** (cancels) · mean absolute **12.72** · std **14.71**

**No global correction can fix this** — there is no global bias. It is invisible from calib
because calib is in-distribution. **Only per-donor reference cells reveal it.**

## 2. The fix, and its measured effect

From k reference cells: `d = median(pred − true)`, then subtract.

**Measured (T16, scalar, k=5):** aggregate ΔAge MAE **14.3 → 6.9 (−52%)**.
At k=3: **14.3 → 7.1 (−50%)**. Both far exceed the ≥25% bar.

## 3. Why Stage 1 must come first — the coupling

`R_eff = max(0, −(mu + z·sigma))`. A **smaller** sigma makes `R_eff` **easier** to be positive:

| state | mu | sigma | R_eff | g |
|---|---|---|---|---|
| now (shifted mu, small sigma) | +8.0 | 2.4 | 0.0 | 0.00 |
| **level fixed, sigma still small** | −11.0 | 2.4 | **8.6** | **0.63** |
| level fixed, sigma correct | −11.0 | 9.0 | 2.0 | **0.29** |

**Applying Stage 2 without Stage 1 makes the safety score ~2× more permissive** — worsening the
over-approval already measured in T7.4.3 (14 approvals vs 11 oracle). The shifted `mu` is
currently *masking* the sigma defect.

**And A′:** Stage 2 changes the error scale, so calibration must be **refitted afterwards**
(q ~30–40 → ~17–21 yr). Hence **Stage 2 = B bundled with A′**, snapshotted together — the one
sanctioned exception to one-change-per-snapshot, because splitting them measures a state that
will never be deployed.

## 4. The conditional rule — do not apply blindly

T16 shows the correction **helps 4 donors and hurts 2**:

| donor | \|shift\| before | MAE before | MAE after (k=5) | effect |
|---|---|---|---|---|
| N2 | 15.0 | 21.8 | **7.1** | huge gain |
| N3 | 28.3 | 29.7 | **10.0** | huge gain |
| O2 | 6.6 | 7.5 | **4.3** | gain |
| Y2 | 20.0 | 14.1 | **5.4** | gain |
| **O1** | **0.6** | 5.4 | 6.1 | **HURT** |
| **Y1** | 8.1 | 7.3 | 8.5 | **HURT** |

O1's shift was already 0.6 yr — essentially calibrated. Estimating an offset from k cells
**injects noise where there was none**.

> **Rule:** compute `d = median(pred − true)` **and** its standard error
> `SE ≈ 1.253·sd/√k`. **Apply only if `|d| > 2·SE`.**

*Numerically verified: at k=3 the asymptotic formula **overstates** the true SE of a median by ~8%
(empirical/asymptotic = 0.92 over 20k draws) — it errs toward **not** correcting, the safe
direction. `SE = sd/√k` is exact if you use the mean; the median matches what T16 measured.*

---

## 5. New file: `src/cellfate/inference/donor_calib.py`

```python
"""Per-donor ΔAge level correction from a handful of labelled reference cells.

The model is unbiased in-distribution but carries a donor-specific level shift out-of-donor
(±12.7 yr, cancelling on average — MASTER_PLAN §5c). The shift is invisible from calib, so the
only way to estimate it for a new donor is from cells of THAT donor with known true ΔAge.

Correction is CONDITIONAL: applying it to an already-calibrated donor injects noise (T16: O1
went 5.4 -> 6.1 MAE). Only correct when the estimate is distinguishable from zero.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

SE_MEDIAN_FACTOR = 1.253      # asymptotic SE(median)/SE(mean) for a normal sample
MIN_K = 2                     # below this, no SE can be formed


@dataclass(frozen=True)
class DonorOffset:
    d: float          # median(pred - true) on the reference cells
    se: float         # standard error of that estimate
    k: int            # number of reference cells used
    applied: bool     # |d| > 2*se
    reason: str       # human-readable explanation, for the report

    def correct(self, mu: np.ndarray) -> np.ndarray:
        return np.asarray(mu, float) - self.d if self.applied else np.asarray(mu, float)


def estimate_offset(pred_ref, true_ref, z: float = 2.0) -> DonorOffset:
    """Estimate the donor's level shift from k reference cells.

    pred_ref : model ΔAge predictions for the reference cells
    true_ref : their measured true ΔAge (requires control + perturbed clock readings)
    z        : how many standard errors the estimate must clear (default 2)
    """
    p = np.asarray(pred_ref, float)
    t = np.asarray(true_ref, float)
    if p.shape != t.shape:
        raise ValueError(f"shape mismatch: pred {p.shape} vs true {t.shape}")
    k = int(p.size)
    if k < MIN_K:
        return DonorOffset(0.0, float("inf"), k, False,
                           f"k={k} below minimum {MIN_K}; no correction")

    resid = p - t
    d = float(np.median(resid))
    sd = float(np.std(resid, ddof=1))
    se = SE_MEDIAN_FACTOR * sd / np.sqrt(k)

    if not np.isfinite(se) or se == 0.0:
        return DonorOffset(d, se, k, False, "degenerate SE; no correction")

    applied = abs(d) > z * se
    reason = (f"|d|={abs(d):.2f} > {z}*SE={z*se:.2f} -> correcting" if applied
              else f"|d|={abs(d):.2f} <= {z}*SE={z*se:.2f} -> donor already calibrated, "
                   f"correction would inject noise")
    return DonorOffset(d, se, k, applied, reason)
```

## 6. Where it hooks into inference

`Predictor.predict_encoded` returns rows containing `mu_age` (`predictor.py:156`). The offset is
applied **after** prediction — the bundle is never modified, so the correction is per-session and
per-donor.

**Recommended: a thin wrapper, not an edit to `Predictor`.**

```python
def predict_with_donor_offset(pred, X, fp, dose_time, offset: DonorOffset) -> list[dict]:
    rows = pred.predict_encoded(X, fp, dose_time)
    if not offset.applied:
        return rows
    for r in rows:
        r["mu_age"] = float(r["mu_age"] - offset.d)
        r["donor_offset_applied"] = True
    return rows
```

**Why a wrapper rather than editing `Predictor`:**
- `Predictor` is used by `scorecard.py`, every diagnostic test, and Stage 3. Changing its output
  silently would alter every one of them.
- The offset is **session state**, not model state — it belongs at the call site.
- Rollback is deleting one function.

> **Do not bake the offset into the bundle.** It is donor-specific and would be wrong for the next
> donor.

## 7. The reference-cell collection protocol

The wet-lab side, stated precisely so it can be handed over:

| Requirement | Detail |
|---|---|
| **How many** | k = 3 minimum (T16: PASS at k=3). k = 5 marginally better; k = 10 no further gain |
| **What is measured** | scRNA-seq of **k control cells** and **k perturbed cells** from the same donor |
| **Why both** | true ΔAge = clock(perturbed) − clock(control). Without the control, no ground truth |
| **When** | any timepoint, but the same protocol as the cells being corrected |
| **Cost** | one extra small sequencing run per donor |

**Sanity check on the reference cells before trusting the offset:**
```python
assert np.std(true_ref) > 0, "reference cells have identical true ΔAge - suspicious"
assert np.all(np.isfinite(pred_ref)) and np.all(np.isfinite(true_ref))
```

## 8. Verification

```python
import numpy as np
from cellfate.inference.donor_calib import estimate_offset

rng = np.random.default_rng(0)

# CASE 1 - a genuinely shifted donor (like N2, +20): should apply
true = rng.normal(-8, 4, 3)
off = estimate_offset(true + 20.0 + rng.normal(0, 1, 3), true)
assert off.applied, off.reason
print("N2-like:", off)

# CASE 2 - an already-calibrated donor (like O1, +0.6): should DECLINE
true = rng.normal(-8, 4, 3)
off = estimate_offset(true + 0.6 + rng.normal(0, 4, 3), true)
assert not off.applied, "must decline a donor that needs no correction"
print("O1-like:", off.reason)

# CASE 3 - too few cells
off = estimate_offset([1.0], [0.0])
assert not off.applied and off.k == 1
```

**All three must pass before the correction is used on real data.** Case 2 is the important one —
it is the failure mode that damaged O1 and Y1 in T16.

## 9. The A′ half — refit calibration afterwards

Stage 2 halves the error, so `q` fitted on the *uncorrected* residuals is now ~2× too wide.

**Rerun Stage 1's `crossdonor_stats` on the corrected model**, i.e. inside the inner-LODO loop
apply the per-donor offset before collecting residuals:

```python
# inside crossdonor_stats, after computing `age` for the held-out inner donor
off = estimate_offset(age[am][:k_ref], ya[am][:k_ref])      # k_ref reference cells
age_corrected = off.correct(age)
res.append(np.abs(age_corrected[am] - ya[am]))
```

> ⚠️ **Leakage warning.** The reference cells used to estimate the offset **must be excluded from
> the residuals** that fit `q`. Otherwise `q` is fitted on cells the correction already saw, and
> coverage will look better than it is.
>
> ```python
> ref_idx = np.arange(k_ref)
> eval_idx = np.setdiff1d(np.arange(am.sum()), ref_idx)
> res.append(np.abs(age_corrected[am][eval_idx] - ya[am][eval_idx]))
> ```

## 10. Failure modes

| Symptom | Cause | Action |
|---|---|---|
| MAE improves but `rank_model_dage` **REGRESSES** | the implementation is doing more than shifting | **revert.** A level shift is rank-invariant by construction |
| Correction applied to every donor including O1 | `z` too low, or SE mis-computed | check case 2 of §8 passes |
| Correction never applied | `z` too high, or reference cells too noisy | inspect `off.reason`; consider k=5 |
| Coverage now overshoots after A′ | `q` refit on corrected residuals, but reference cells leaked in | apply the §9 exclusion |
| MAE improves, `level_shift` does not | the gain came from elsewhere | **investigate before keeping** — the diagnosis predicts they move together. An unexplained win is a red flag |

## 11. Rollback

Delete the wrapper call. `Predictor` is untouched, the bundle is untouched, and every diagnostic
reverts automatically because none of them use the wrapper.

## 12. Acceptance

```powershell
python scorecard.py snapshot --tag B_percalib
python scorecard.py compare A_xdonor B_percalib
```

| Role | Metric | Bar |
|---|---|---|
| **TARGET** | `dage_mae_model` | ACCEPT + ≥25% drop (T16 predicts ~50%) |
| **GUARD** | `rank_model_dage` | **noise or ACCEPT — never REGRESSION** |
| **GUARD** | per-fold | no fold may worsen by >20% |
| **WATCH** | `conformal_width` | should **fall** (~35–43 → ~17–21 yr) |

**The ranking guard is the real test.** A pure level shift is rank-invariant, so ranking *must not
move*. If it does, the implementation is wrong.

**If it fails:** revert, keep Stage 1, report the model as a within-donor ranker with honest wide
intervals. That remains a valid product.

## 13. Interaction with the tool

**Stage 3 does not need this stage.** Comparing "day 15 vs day 21 **for this donor**" cancels a
constant offset. Stage 2 is needed only to say *"you gained 9 years"* rather than *"day 15 beats
day 21."*

Stage 3d must therefore **warn** when absolute ΔAge is reported without a donor offset applied.

## 14. Done when

Acceptance passes, all three verification cases in §8 pass, the leakage exclusion in §9 is
implemented, and the reference-cell protocol is documented wherever the tool is described.
