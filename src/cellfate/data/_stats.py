"""Small, dependency-light numerical helpers reused across the ETL."""

from __future__ import annotations

import numpy as np

_EPS = 1e-8


def standardize(expr: np.ndarray) -> np.ndarray:
    """Per-column (per-gene) z-score across cells. ``expr`` is (N, G)."""
    expr = np.asarray(expr, dtype=np.float64)
    mu = expr.mean(axis=0, keepdims=True)
    sd = expr.std(axis=0, keepdims=True)
    return (expr - mu) / (sd + _EPS)


def softmax(z: np.ndarray, axis: int = -1) -> np.ndarray:
    """Numerically stable softmax."""
    z = np.asarray(z, dtype=np.float64)
    z = z - z.max(axis=axis, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=axis, keepdims=True)


def ols_slope_intercept(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Ordinary least squares for ``y ~ a*x + b``; returns (a, b).

    Falls back to (0, mean(y)) when ``x`` has no variance (cannot fit a slope).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2 or np.allclose(x.std(), 0.0):
        return 0.0, float(y.mean()) if y.size else 0.0
    a, b = np.polyfit(x, y, 1)
    return float(a), float(b)
