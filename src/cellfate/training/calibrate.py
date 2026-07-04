"""Post-hoc confidence calibration (Document 3, S5).

Temperature scaling fits a single scalar T on the ensemble-averaged held-out
logits by minimising NLL with a bounded L-BFGS-B optimiser (L-BFGS with box
constraints on T). Because T=1 is feasible and we keep an explicit guard, the
calibrated NLL can never be worse than the uncalibrated one. T is shipped in the
bundle and divides the logits before softmax at inference.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_softmax

from cellfate.common.schemas import TemperatureParams

_BOUNDS = (1e-2, 1e2)


def _nll(logits: np.ndarray, target: np.ndarray, t: float) -> float:
    lp = log_softmax(logits / t, axis=1)
    return float(-(target * lp).sum(axis=1).mean())


def fit_temperature(logits, target, bounds: tuple[float, float] = _BOUNDS) -> TemperatureParams:
    """Fit the calibration temperature on (logits, soft target) with L-BFGS-B.

    A single scalar T is optimised by box-constrained L-BFGS (L-BFGS-B) to
    minimise NLL, starting from T=1. Returns ``TemperatureParams``; guaranteed
    never worse than T=1.
    """
    logits = np.asarray(logits, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if logits.size == 0:
        return TemperatureParams(temperature=1.0)

    res = minimize(
        lambda x: _nll(logits, target, float(x[0])),
        x0=np.array([1.0]),
        method="L-BFGS-B",
        bounds=[bounds],
    )
    t = float(np.clip(res.x[0], bounds[0], bounds[1])) if res.success else 1.0
    # never ship a temperature that is worse than no calibration
    if not np.isfinite(t) or _nll(logits, target, t) > _nll(logits, target, 1.0) + 1e-12:
        t = 1.0
    return TemperatureParams(temperature=t)
