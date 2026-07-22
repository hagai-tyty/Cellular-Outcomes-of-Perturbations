"""Post-hoc confidence calibration (Document 3, S5).

Two calibrators live here.

``fit_temperature`` -- a single scalar T fitted on the ensemble-averaged held-out logits by
minimising MULTI-CLASS NLL with a bounded L-BFGS-B optimiser. T=1 is feasible and guarded, so
the calibrated NLL can never be worse than the uncalibrated one.

``fit_platt_binary`` -- a 2-parameter Platt fit on the SAFE-vs-rest boundary, which is the
quantity the product actually ships: ``res.py`` consumes ``S = P(safe)`` and ``P_loss``,
``STAGE_3`` S0.1 requires a risk threshold on ``P(unsafe)``, and ``scorecard.py`` grades binary
ECE on ``P(safe)``. Stage 1 run 2 calibrated multi-class NLL instead and REGRESSED that metric
(0.281 -> 0.364 on every fold), because the two objectives disagree about how much to sharpen.

Platt's slope IS a temperature on the binary logit, so it subsumes the scalar and adds the
intercept a scalar cannot express -- the reason in-distribution Platt reaches 0.153 where
temperature reaches 0.281 on this data. It is monotone in ``P(safe)``, so ranking metrics
(fate PR-AUC / ROC-AUC) are preserved exactly.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_softmax

from cellfate.common.calibration import EPS as _EPS
from cellfate.common.calibration import apply_platt, platt_safe
from cellfate.common.schemas import TemperatureParams

__all__ = ["fit_temperature", "fit_platt_binary", "has_class_variation", "apply_platt"]

_BOUNDS = (1e-2, 1e2)
_PLATT_BOUNDS = ((1e-2, 1e2), (-1e2, 1e2))   # slope must stay positive: rank preservation


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


def _binary_logloss(p_safe: np.ndarray, y_safe: np.ndarray, a: float, b: float) -> float:
    p = np.clip(platt_safe(p_safe, a, b), _EPS, 1.0 - _EPS)
    return float(-(y_safe * np.log(p) + (1.0 - y_safe) * np.log(1.0 - p)).mean())


def fit_platt_binary(p_safe, y_safe) -> tuple[float, float]:
    """Fit ``P(safe) -> sigmoid(a*logit(P(safe)) + b)`` on cross-donor probabilities.

    WHY THIS AND NOT A TEMPERATURE. ``P(safe)`` is what the product ships and what the bar
    grades: ``res.py`` consumes ``S`` and ``P_loss``, ``STAGE_3`` S0.1 needs a risk threshold on
    ``P(unsafe)``, and ``scorecard.py`` scores binary ECE on ``P(safe)``. MASTER_PLAN S5a names
    the defective quantity as "``S``, ``P_loss``" and records "YES -- Platt halves it" (T8.2);
    the ~0.13 that STAGE_1's <=0.17 bar is derived from was MEASURED with Platt. Stage 1 run 2
    fitted a scalar on multi-class NLL instead and regressed the metric on every fold
    (0.281 -> 0.364), because those two objectives disagree about how much to sharpen.

    RANK PRESERVATION. The slope is constrained positive, so the map is strictly increasing in
    ``P(safe)``. Every ranking metric (fate PR-AUC, ROC-AUC) is therefore mathematically
    unchanged -- which is what lets this ship without disturbing Stage 1's guards.

    Returns the identity ``(1.0, 0.0)`` when the data cannot identify a fit (empty, or a single
    class present), and never returns a fit worse than the identity on its own objective.
    """
    p = np.asarray(p_safe, dtype=np.float64).ravel()
    y = np.asarray(y_safe, dtype=np.float64).ravel()
    if p.size == 0 or p.size != y.size:
        return 1.0, 0.0
    # a single class carries no information about where the boundary should sit
    pos = float(y.sum())
    if not (0.0 < pos < float(y.size)):
        warnings.warn(
            "Platt calibration skipped: the safe/unsafe target carries a single class, so the "
            "boundary is unidentifiable. Shipping the identity (1.0, 0.0) uncalibrated.",
            stacklevel=2,
        )
        return 1.0, 0.0

    res = minimize(lambda w: _binary_logloss(p, y, float(w[0]), float(w[1])),
                   x0=np.array([1.0, 0.0]), method="L-BFGS-B", bounds=_PLATT_BOUNDS)
    a, b = (float(res.x[0]), float(res.x[1])) if res.success else (1.0, 0.0)
    if not (np.isfinite(a) and np.isfinite(b)):
        return 1.0, 0.0
    # never ship a calibration worse than none, mirroring fit_temperature's guard
    if _binary_logloss(p, y, a, b) > _binary_logloss(p, y, 1.0, 0.0) + 1e-12:
        return 1.0, 0.0
    return a, b


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
