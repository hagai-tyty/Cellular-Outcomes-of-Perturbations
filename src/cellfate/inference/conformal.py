"""Guaranteed-coverage ΔAge interval from the stored conformal quantile (Doc 4, S4).

Training calibrates ``q`` as a quantile of the absolute residual
``|ensemble-mean ΔAge − true ΔAge|`` on a held-out split, so the marginal-coverage
interval for the ensemble-mean prediction is ``[µ − q, µ + q]`` (``q`` is already
in age units). The epistemic spread ``σ`` is reported separately and drives the
RES score; it is not multiplied into the interval, which would break the coverage
guarantee that the raw-residual calibration established.
"""

from __future__ import annotations

import numpy as np


def interval(mu: float, q: float) -> list[float]:
    return [float(mu) - float(q), float(mu) + float(q)]


def intervals(mu, q) -> np.ndarray:
    mu = np.asarray(mu, dtype=np.float64)
    return np.stack([mu - float(q), mu + float(q)], axis=1)
