"""Post-hoc confidence calibration (Document 3, S5).

Temperature scaling fits a single scalar T on the ensemble-averaged held-out
logits by minimising NLL with a bounded L-BFGS-B optimiser (L-BFGS with box
constraints on T). Because T=1 is feasible and we keep an explicit guard, the
calibrated NLL can never be worse than the uncalibrated one. T is shipped in the
bundle and divides the logits before softmax at inference.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_softmax

from cellfate.common.schemas import TemperatureParams

_BOUNDS = (1e-2, 1e2)


def _nll(logits: np.ndarray, target: np.ndarray, t: float) -> float:
    lp = log_softmax(logits / t, axis=1)
    return float(-(target * lp).sum(axis=1).mean())


def has_class_variation(target, min_mass_frac: float = 0.01) -> bool:
    """True iff at least two classes carry real mass in ``target``.

    Temperature scaling is UNIDENTIFIABLE without it: on effectively single-class data the NLL
    falls monotonically as T -> 0, because "always this class, with certainty" is optimal. The
    optimiser then runs to the lower bound and the "never worse than T=1" guard *passes* -- the
    fit really is better on that data -- so a maximally overconfident temperature ships silently.
    """
    mass = np.asarray(target, dtype=np.float64).sum(axis=0)
    total = float(mass.sum())
    if not np.isfinite(total) or total <= 0:
        return False
    return int((mass / total >= min_mass_frac).sum()) >= 2


def fit_temperature(logits, target, bounds: tuple[float, float] = _BOUNDS) -> TemperatureParams:
    """Fit the calibration temperature on (logits, soft target) with L-BFGS-B.

    A single scalar T is optimised by box-constrained L-BFGS (L-BFGS-B) to
    minimise NLL, starting from T=1. Returns ``TemperatureParams``; guaranteed
    never worse than T=1.

    Returns T=1 (no calibration) when the data cannot identify a temperature -- empty, or
    effectively single-class. Stage 1b made that reachable: the cross-donor pool is small
    (~75 cells over 5 donors here, against ~4,400 for the old in-distribution split), so a
    fold whose held-out donors happen to be nearly all one class is a real possibility rather
    than a theoretical one.
    """
    logits = np.asarray(logits, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if logits.size == 0:
        return TemperatureParams(temperature=1.0)
    if not has_class_variation(target):
        warnings.warn(
            "temperature scaling skipped: the calibration targets carry mass in fewer than two "
            "classes, which makes T unidentifiable (NLL falls monotonically toward T=0, i.e. "
            "maximal confidence). Shipping T=1.0 uncalibrated instead of an overconfident fit.",
            stacklevel=2,
        )
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
