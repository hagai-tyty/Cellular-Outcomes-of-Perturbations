"""Proliferation deconfounding (Document 2, S11).

Transcriptomic age clocks are confounded by proliferation: rapidly dividing
cells read as "younger". If left in, the model could learn to score a
perturbation as rejuvenating simply because it pushes cells to divide -- exactly
the cancer-risk failure mode this project tries to avoid. We therefore remove the
linear cell-cycle component from the ΔAge label.

The cell-cycle score is the mean standardised expression of the S + G2M gene
sets (Tirosh et al.); the deconfounder is ``ΔAge - (a*cc + b)`` with (a, b) fit by
OLS on the aged cells (the slope ``a`` is the artefact to remove).
"""

from __future__ import annotations

import numpy as np

from cellfate.common import constants as C

from ._stats import ols_slope_intercept, standardize


def cell_cycle_score(
    norm_expr: np.ndarray,
    genes: list[str],
    s_genes: tuple[str, ...] = C.S_GENES,
    g2m_genes: tuple[str, ...] = C.G2M_GENES,
) -> np.ndarray:
    """Mean standardised expression over S + G2M genes present (N,)."""
    cc = list(s_genes) + list(g2m_genes)
    idx = [genes.index(g) for g in cc if g in genes]
    if not idx:
        return np.zeros(norm_expr.shape[0], dtype=np.float64)
    return standardize(norm_expr[:, idx]).mean(axis=1)


def fit_deconfounder(delta_age: np.ndarray, cc_score: np.ndarray) -> tuple[float, float]:
    """Fit ``delta_age ~ a*cc + b`` by OLS; returns (a, b)."""
    return ols_slope_intercept(cc_score, delta_age)


def deconfound_age(
    delta_age: np.ndarray, cc_score: np.ndarray, coef: tuple[float, float]
) -> np.ndarray:
    """Remove the linear cell-cycle component: ``delta_age - (a*cc + b)``."""
    a, b = coef
    return np.asarray(delta_age, dtype=np.float64) - (a * np.asarray(cc_score, dtype=np.float64) + b)
