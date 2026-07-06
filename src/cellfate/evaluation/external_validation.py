"""External ΔAge validation (Document 5, S4).

The model must track *real* rejuvenation, not merely reproduce the training clock.
Two independent, data-dependent checks: correlation of predicted ΔAge against a
DNA-methylation clock, and recovery of the known trajectory on a held-out OSKM
partial-reprogramming time course (ΔAge should fall while identity-loss probability
rises only late -- the project's core age/identity-decoupling thesis). These require
external data; when none is configured the evaluator records them as not-available
rather than fabricating a number.
"""

from __future__ import annotations

import numpy as np

NAN = float("nan")


def validate_against_methylation(pred_dage, meth_dage) -> dict[str, float]:
    a = np.asarray(pred_dage, dtype=np.float64)
    b = np.asarray(meth_dage, dtype=np.float64)
    pear = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else NAN
    return {"pearson": pear, "bias": float((a - b).mean()), "n": int(len(a))}


def validate_oskm_holdout(predict_fn, timecourse) -> dict:
    """``timecourse``: ordered list of (X, fp, dose_time) at increasing reprogramming
    time. ``predict_fn(X, fp, dose_time) -> (probs (N,3), age (N,))``. Checks that
    mean ΔAge falls along the course and that identity-loss probability rises later
    than rejuvenation begins (decoupling)."""
    ages, ploss = [], []
    for X, fp, dt in timecourse:
        probs, age = predict_fn(X, fp, dt)
        ages.append(float(np.mean(age)))
        ploss.append(float(np.mean(np.asarray(probs)[:, 1])))
    ages = np.asarray(ages)
    ploss = np.asarray(ploss)
    third = max(1, len(ploss) // 3)
    early_rise = float(ploss[third] - ploss[0]) if len(ploss) > third else 0.0
    late_rise = float(ploss[-1] - ploss[-1 - third]) if len(ploss) > third else 0.0
    return {
        "delta_age_trajectory": ages.tolist(),
        "p_loss_trajectory": ploss.tolist(),
        "rejuvenates": bool(ages[-1] < ages[0]) if len(ages) > 1 else False,
        "age_identity_decoupled": bool(late_rise >= early_rise),
    }
