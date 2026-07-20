# STAGE 4 — Validation

**Implements:** `MASTER_PLAN.md` §9 criteria, §4 (Change C — the deferred RES verdict).
**Depends on:** Stage 3 producing recommendations.
**Feeds:** the decision to ship, narrow, or abandon.
**Scope:** 1 new script (~200 lines), a scorecard extension, and one honest verdict.

---

## 0. What this stage is for

Everything up to here has been *building*. This stage asks whether it works — and is written so
that **"it doesn't" is a reportable outcome, not a failure of the plan.**

---

## 1. Three questions, in increasing order of what they're worth

| # | Question | Metric | Pre-registered bar |
|---|---|---|---|
| **V1** | Are the uncertainties honest? | conformal coverage, held-out donors | **0.85–0.95** (from 0.40) |
| **V2** | Are the safety probabilities honest? | fate ECE, held-out donors | **≤ 0.15** (from 0.28) |
| **V3** | **Does following the recommendation beat a fixed protocol?** | ΔAge and P(unsafe) at the recommended day vs a fixed day | **wins on ≥4/6 donors** |

**V3 decides whether this is a product.** V1 and V2 are prerequisites — a recommendation with
dishonest error bars is worse than none, because it invites misplaced confidence.

---

## 2. New file: `validate_stopping.py`

```python
"""Stage 4 validation: does the stopping recommendation beat a fixed protocol?

Leave-one-donor-out. For each held-out donor:
  1. take its EARLIEST timepoint as "the sample the researcher submits"
  2. ask the tool for a recommendation
  3. look up what ACTUALLY happened at the recommended day in that donor's real trajectory
  4. compare against baselines

The comparison is on OBSERVED outcomes, never on predictions — otherwise the tool is being
graded by its own homework.
"""
from __future__ import annotations
import numpy as np
from cellfate.inference.stopping import recommend_stopping

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
RISK = 0.20


def observed_at(rows, day):
    """The donor's TRUE (delta_age, unsafe_fraction) at the timepoint nearest `day`."""
    i = int(np.argmin([abs(r["t"] - day) for r in rows]))
    return rows[i]["y_age_true"], rows[i]["frac_unsafe_true"], rows[i]["t"]


def wins(rec, base):
    """Recommendation beats baseline if it achieves MORE rejuvenation at equal-or-lower true
    unsafe fraction, OR lower unsafe fraction at equal-or-better delta_age.
    TIES COUNT AS LOSSES - the tool must earn its complexity."""
    r_age, r_unsafe = rec
    b_age, b_unsafe = base
    return (r_age < b_age and r_unsafe <= b_unsafe) or (r_unsafe < b_unsafe and r_age <= b_age)


def main():
    fixed_day = None      # median withdrawal day across TRAINING donors, computed per fold
    results = []
    for held in DONORS:
        rows = load_true_trajectory(held)              # from the held-out donor's test split
        train_days = [median_day(d) for d in DONORS if d != held]
        fixed_day = float(np.median(train_days))

        report = recommend_stopping(X_now=rows[0]["x"], bundle_root=f"forward_{held}",
                                    risk_threshold=RISK)
        if report.recommendation is None:
            results.append({"donor": held, "outcome": "no safe window", "win": False})
            continue

        rec_obs   = observed_at(rows, report.recommendation.withdraw_day)[:2]
        fixed_obs = observed_at(rows, fixed_day)[:2]
        late_obs  = observed_at(rows, max(r["t"] for r in rows))[:2]
        early_obs = observed_at(rows, min(r["t"] for r in rows))[:2]

        results.append({
            "donor": held,
            "rec_day": report.recommendation.withdraw_day,
            "fixed_day": fixed_day,
            "rec": rec_obs, "fixed": fixed_obs, "late": late_obs, "early": early_obs,
            "win": wins(rec_obs, fixed_obs),
        })
    report_table(results)


if __name__ == "__main__":
    main()
```

## 3. The three baselines it must beat

| Baseline | Why it belongs |
|---|---|
| **Fixed day** (median training withdrawal day) | current practice — **the real alternative** |
| **Always latest** | maximum rejuvenation, ignoring safety |
| **Always earliest** | maximum safety, ignoring rejuvenation |

Beating "fixed day" is the bar that matters. **If the tool cannot beat a constant, it has no reason
to exist.** This project has repeatedly found simple baselines winning (T5, T6, T7, T9) — expect
that here too, and be pleased if it does not happen.

## 4. The comparison rule, stated precisely

The recommendation **wins on a donor** if it achieves:
- **more rejuvenation at equal-or-lower true unsafe fraction**, **or**
- **lower unsafe fraction at equal-or-better ΔAge**

**Ties count as losses.** The tool must earn its complexity, not draw with a constant.

**Everything is scored on the donor's OBSERVED trajectory**, never on predictions. Grading the
recommendation with the model that produced it is circular.

## 5. Scorecard extension

Add to `scorecard.py` so forward metrics live alongside the existing ones:

```python
"forward_coverage":   ("higher", "forward conformal coverage"),
"forward_dage_mae":   ("lower",  "forward ΔAge MAE"),
"forward_unsafe_ece": ("lower",  "ECE of predicted unsafe fraction"),
"dt_response":        ("neutral","ΔAge swing across the Δt sweep"),   # must exceed ~2 yr
"rec_beats_fixed":    ("higher", "donors where the recommendation wins"),
```

Same discipline: paired 95% CI across donors, **ACCEPT** only when the CI excludes zero in the
improving direction.

## 6. Honest failure modes

| Result | Meaning | Action |
|---|---|---|
| V1 fails (coverage still ~0.4) | cross-donor calibration did not transfer | **do not ship** — the uncertainty is unusable |
| V2 fails | safety probabilities uncalibrated | ship **relative ranking only**, no risk threshold |
| V3 fails, V1+V2 pass | honest but not useful | ship as a **calibrated readout, not a recommender** |
| `dt_response` still ≈ 0 | the model never learned the forward signal | Stage 3c failed; revisit the pairing |
| **Recommendation wins 6/6** | **suspect leakage first** | re-audit the donor split before celebrating |

That last row is deliberate. Given how consistently simple baselines have won here, a clean sweep
is more likely a bug than a breakthrough. **Specific things to re-check:** did any pair cross
donors? Did the reference cells leak into the residuals (Stage 2 §9)? Was `fixed_day` computed
from training donors only?

## 7. Change C — the deferred RES verdict (`MASTER_PLAN.md` §4)

RES was never refuted, only **untested on defective inputs**.
`RES = Φ(S)·S^k·g(R_eff)·exp(−λ·P_loss)` — **every input it consumes is defective:**

| RES input | Defect | Fixed by |
|---|---|---|
| `mu_age` (via `R_eff`) | per-donor level shift ±12.7 yr → `R_eff = 0` everywhere (T7.4.2) | Stage 2 |
| `sigma_age` (via `R_eff`) | intervals cover 0.40 vs 0.90 → sigma understated (T14) | Stage 1 |
| `S`, `P_loss` | ECE 0.28 out-of-donor (T8.2) | Stage 1 |
| `in_dist` | AUC 0.47 ≈ chance, zeroes RES arbitrarily (T15) | Stage 1 |
| `lam` | `lam = 0` makes the `P_loss` term **inert** (T7.4.2) | config, trivially |

**A score built on four defective inputs cannot be judged on its formula.** Stages 1–3 fix all of
them. **Only now can RES be judged — this is Test 7.6.**

Rerun Tests 7, 7.1 and 7.4.1 on the corrected model:

| Outcome | Decision |
|---|---|
| `rank_res` within **0.03** of the ΔAge sort, over-approval gap ≤ 0.5 | **RES vindicated** — keep as the headline score |
| within **0.05**, gap ≤ 1.0 | keep as a **secondary safety filter**; rank by ΔAge |
| still ≥0.10 below the ΔAge sort | **now** demotion is justified, with evidence |

**Expect per-cell RES to still approve nothing — and for that to be correct.** Honest per-cell
uncertainty (~19 yr) exceeds the effect (~11 yr). **The real test is whether the condition-level
score from Stage 3d beats a condition-level ΔAge sort.**

## 8. The sample-size caveat that must travel with any claim

**n = 6 donors. V3 is 6 comparisons.** Even 5/6 is weak — a fair coin lands 5/6 about **11%** of
the time. **This validation can rule the tool out; it cannot establish that it works.**

**What would establish it:** prospective use on donors never seen during development, ideally at a
different site. Until then the honest phrasing is *"recommendations were consistent with observed
outcomes on 6 retrospective donors."*

## 9. Interaction with other stages

- **← Stage 3d** must expose the recommendation **programmatically**, not only as printed text.
- **← Stage 3b** owns the donor-split guarantee. **Stage 4 re-verifies it independently** — this is
  the last point at which leakage can be caught, and the most likely explanation for a
  suspiciously good result.
- **→ Stage 5** determines which claims survive.

## 10. Done when

All three questions are answered with recorded numbers and verdicts against §1's bars, Change C
has a verdict, and the result is written into the lab notebook **whether it passes or fails**.
