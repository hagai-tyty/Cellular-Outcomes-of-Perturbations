"""STAGE 14 -- the ΔAge scale correction, applied at the REPORTING boundary only.

Stage 11 measured that the dense transcriptomic clock was never broken, only mis-scaled: raw ΔAge
carries the right ordering (Spearman 0.770 vs methylation) with a 66 % magnitude inflation
(SD ratio 1.658). A single multiplicative factor, fitted leave-one-donor-out, corrects it.

WHY THE REPORTING BOUNDARY AND NOT THE TARGET
---------------------------------------------
Rescaling `y_age` would need a rebuild + retrain, and Stage 14's pre-flight measured that it also
changes the LOSS REGIME as a side effect: `huber_delta = 2.0` is a knee fixed in absolute years,
and training residuals are 1.36-2.40 yr, so shrinking the target moves the fraction of the loss
inside the quadratic region from 43-67 % to 85-97 %. That is a second change wearing the costume
of the first.

Applied here instead, the correction touches only what a human reads. It cannot change model
skill, cannot change ranking (a positive rescale is rank-invariant), and cannot be mistaken for
either. Crucially it also leaves `res.py` untouched: `kappa = 5.0` is a rejuvenation
half-saturation **in years**, so rescaling ΔAge upstream of RES would silently reinterpret it --
the same class of defect as the Huber knee.

WHICH FACTOR
------------
Two were available. `k_LS` minimises MAE against methylation and reaches 6.78 yr, beating the
7.30 yr instrument floor -- but its SD ratio is 0.597, i.e. it wins partly by UNDER-REPORTING
magnitude by 40 %. For a reporting transform the objective is an unbiased magnitude, not a minimal
error, so the variance-matched factor ships and its worse MAE (10.68) is the honest number rather
than a defeat. Choosing `k_LS` would be selecting the estimator that flatters the headline -- the
shrinkage trap this project has already hit once.

WHAT IS NOT ESTABLISHED
-----------------------
`k` was fitted on donors **O1/O2/O3 of the transient arm** -- the only rows in
`results/dage_ledger.csv` carrying `TRUTH_meth_dage_mt` (68 of 90; the Sendai cohort carries
none). The cohort with methylation truth and the cohort the model trains on are **disjoint**, and
Stage 11 §11.4 explicitly forbade claiming transfer. Calibrated values are therefore reported
ALONGSIDE the raw ones, never in place of them, and always carrying `CAVEAT`.
"""

from __future__ import annotations

# Stage 11, `results/diag_stage11_scale_results.json`, variant `raw`, mean of the LODO
# variance-matched factors (O1 0.5836, O2 0.5406, O3 0.6730). Spread 1.19x across donors,
# against a pre-registered stability bar of 2x.
K_VAR = 0.5991

# The least-squares alternative, recorded so the choice stays visible and reversible. NOT used:
# SD ratio 0.597 means it under-reports magnitude by 40 %.
K_LS = 0.3637

CAVEAT = (
    f"delta_age_calibrated applies a single scale factor k={K_VAR:.4f} fitted leave-one-donor-out "
    "on 3 donors of the transient arm against methylation truth. That cohort is DISJOINT from the "
    "one this model trains on, so the transfer is UNTESTED. Ranking is unaffected (a positive "
    "rescale cannot reorder). The raw delta_age_mean is reported alongside and is unchanged."
)


def calibrate(value: float, k: float = K_VAR) -> float:
    """Scale one ΔAge-like quantity (a mean, a half-width, a bound) into calibrated years."""
    return float(value) * float(k)


def calibrate_interval(lo: float, hi: float, k: float = K_VAR) -> list[float]:
    """Scale a conformal interval.

    NOTE ON WHAT THE 90 % THEN MEANS. `q` was calibrated so that `|mu - y| <= q` covers the
    nominal fraction of the **RNA-clock label** `y`. Scaling both endpoints keeps the interval
    coherent with the scaled point estimate and preserves that coverage exactly, because the
    label scales with them. It does **not** establish coverage against methylation truth, which
    is a different quantity and has never been measured here. Stated in `CAVEAT` rather than
    implied by the units.
    """
    return [calibrate(lo, k), calibrate(hi, k)]
