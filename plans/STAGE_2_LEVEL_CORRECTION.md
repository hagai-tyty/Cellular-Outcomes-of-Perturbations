# STAGE 2 — Per-donor level correction (Change B + A′)

**Implements:** `MASTER_PLAN.md` §5c, §5d.
**Depends on:** Stage 1 — a **coupling**, not a preference (§3).
**Blocking for:** absolute ΔAge claims only. **Not required for the tool.**
**Scope:** 1 new file, 2 modified, ~150 lines, plus a wet-lab protocol decision.

---

> ## ⚠️ UPDATED 2026-07-26 — the control-swap threat is gone; the premise is still unverified
>
> *Additive note. Nothing below this box is modified. Supersedes the earlier warning that cited
> `STAGE_1_5_1_REVISED.md`, which proposed a fix that has since been rejected.*
>
> **What is settled.** The earlier warning demanded this stage wait for a **non-responder control
> swap**. That swap is **REJECTED as unsafe** — `STAGE_1_5_1_REVISED_REVIEW.md` showed the arms are
> not identity-matched, and `STAGE_1_5_1_REV_FINAL.md` then confirmed on DNA methylation
> (`GSE165179`) that the **+36.5 yr non-responder reading is a transcriptomic artefact**: measured
> against a real untreated control, failed-to-reprogram cells are **inert (+0.5 / −2.4 yr)**.
> Rejuvenation itself is **real and large** (−24 to −28 yr). **So this stage is no longer blocked by
> a pending label change.** `src/` remains **100% untouched**.
>
> **What is NOT settled — do not read the above as clearance.** The original blocker asked for one
> specific thing: *re-measure the per-donor level shift on corrected labels before spending.*
> **That was never done, because the labels were never corrected.** `REV FINAL` §8 states the
> position plainly: methylation anchors ΔAge only *where methylation exists* (Gill's 3 donors), and
> **HFF's age labels — ~99.8% of all age labels — remain unanchored.** The shift in §1 is still
> measured on RNA labels produced by a clock now **proven out of domain on reprogramming cells**.
>
> **Two live findings that still apply to §1's numbers:**
>
> 1. **The artefact varies enormously between donors**, which is the exact shape this stage exists to
>    correct — so the two remain confounded in the RNA labels:
>
>    | donor | N2 | N3 | O1 | O2 | Y1 | Y2 |
>    |---|---|---|---|---|---|---|
>    | non-responder day0→peak (RNA, now known artefactual) | −15.7 | +36.5 | +33.3 | +44.8 | +50.3 | +69.8 |
>
> 2. **§1's ±12.72 yr is mis-attributed** (`STAGE_1_DEVIATIONS.md` §C1, independent of all the above):
>    it is the **ridge** baseline's shift. The **model's** is mean **−5.71**, 95% CI
>    **[−22.9, +11.5] — which includes zero.** Use the model numbers when judging whether this stage
>    is worth its cost.
>
> ### Why this stage may nonetheless proceed
>
> Not because the premise was verified — it wasn't — but because **k ≈ 3 reference cells per donor
> helps under both hypotheses**:
>
> * if the per-donor offset is **real donor biology**, correcting it is exactly the point;
> * if it is **n = 1 baseline noise** (every donor's zero-point is a single unreplicated day-0
>   sample), then replicating the baseline fixes it directly.
>
> **The justification changes; the action does not.** State it that way in any write-up — do not
> claim the ±12.7 yr offset has been established as biology.
>
> **Worth weighing first:** now that methylation is in hand, budget spent on **methylation for more
> donors** would settle more than reference cells will — it would resolve the open retention question
> (`REV FINAL` §5, needs ≈ 16 pairs) *and* begin anchoring the labels. Reference cells address only
> the baseline-noise half.

> ## 🆕 ADDED 2026-07-30 — three things a reviewer will find in ten minutes
>
> *Additive note; nothing below is modified. All three bear on the **§0 spend decision**, which is
> why they sit above it.*
>
> ### A. §1 and §4 quote DIFFERENT MODELS — this document never says so
>
> Every donor number disagrees between the two tables, with no explanation in the file:
>
> | donor | N2 | N3 | O1 | O2 | Y1 | Y2 | mean | mean abs |
> |---|---|---|---|---|---|---|---|---|
> | **§1's table** = **ridge** baseline | +20.11 | −24.40 | +5.72 | +13.04 | −4.28 | −8.81 | +0.23 | **12.72** |
> | **§4's table** = the **model** | +15.03 | −28.35 | **+0.64** | +6.56 | −8.13 | −20.02 | **−5.71** | 13.12 |
>
> This is `STAGE_1_DEVIATIONS.md` **§C1**, and it is already recorded there — but a reader of *this*
> file sees only two contradictory tables. **§4's are the deployed model's; §1's are not.** Use §4's
> for any decision. (§C1's own correction also applies: the model's mean −5.71 has 95% CI
> **[−22.9, +11.5]**, which includes zero, so it is *not* evidence of a global bias either.)
>
> ### B. §2's benefit was measured with the correction applied to EVERY donor — before §4's rule existed
>
> §2 reports **14.3 → 6.9 (−52%)** at k=5 and **14.3 → 7.1 (−50%)** at k=3, and §4 then says the same
> T16 run *"helps 4 donors and hurts 2"*. **So §2's headline is the unconditional number.** §4's
> `|d| > 2·SE` rule exists precisely to suppress some of those corrections — and **nobody re-measured
> the benefit with the rule active.** The true figure is **between** §2's −50% and zero, and is
> currently unknown. **Do not quote −50% as the expected gain of the stage as specified.**
>
> ### C. Whether the rule fires at k = 3 is never computed — and it decides the spend
>
> §4's rule is `|d| > 2·SE` with `SE ≈ 1.253·s/√k`, where **`s` is the within-donor sd of
> `pred − true` across the reference cells**. Substituting k = 3:
>
> ```
> fires  <=>  |d|  >  2 · 1.253 · s / sqrt(3)  =  1.447 · s
> ```
>
> So a donor is corrected only if its shift exceeds **~1.45× the within-donor spread**. Against §4's
> model shifts, the break-even `s` per donor is `|d| / 1.447`:
>
> | donor | model \|shift\| | needs `s` below | T16 verdict |
> |---|---|---|---|
> | N3 | 28.35 | **19.6** | huge gain |
> | Y2 | 20.02 | **13.8** | gain |
> | N2 | 15.03 | **10.4** | huge gain |
> | Y1 | 8.13 | 5.6 | **HURT** |
> | O2 | 6.56 | 4.5 | gain |
> | O1 | 0.64 | 0.44 | **HURT** |
>
> **Best available estimate of `s`:** §2's corrected MAE is 6.9 yr, so for roughly normal residuals
> `s ≈ 1.253 × 6.9 ≈ 8.6 yr`. **This is inferred, not measured** — see the pre-registration below.
>
> **At `s ≈ 8.65` the threshold is 12.5 yr, and the rule fires for 3 of 6 donors: N3, Y2, N2.**
> That is the rule working as designed — **it correctly declines both donors T16 damaged** (O1, Y1).
> But it **also declines O2**, which T16 *helped* (7.5 → 4.3). Capturing O2 as well needs
> `k > 6.28·s²/d²` ≈ **11 cells**, not 3; Y1 would need ≈ 7 — and Y1 is a donor the rule *should*
> decline.
>
> **Consequence for §0 and §7:** "k = 3 minimum" is the number at which the *unconditional*
> correction passed T16. It is **not** established as sufficient for the *conditional* correction
> this document actually specifies. k = 3 buys the three large-shift donors and nothing else.
>
> ### Pre-registered before any wet-lab spend (ground rule §5)
>
> 1. **Measure `s` directly** — the within-donor sd of `pred − true`, per donor, on existing held-out
>    data. It needs no new cells and it is the single number this stage's cost-effectiveness turns on.
> 2. **Re-measure the benefit with §4's rule active**, at k = 3 and k = 5, and report it beside §2's
>    unconditional figure. The **≥25% TARGET bar in §12 must be graded on the conditional number**,
>    since that is what would ship.
> 3. **If the conditional benefit misses the bar**, the honest options are k > 3 (cost rises), or
>    §0's stated fallback — ship as a within-donor ranker. **Not** relaxing the bar (§5), and **not**
>    quietly reverting to the unconditional correction that T16 showed damages O1 and Y1.

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
