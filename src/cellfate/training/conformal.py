"""Conformal calibration of the ΔAge interval (Document 3, S5).

Split-conformal regression on the calibration set: for coverage level ``lvl`` the
half-width ``q[lvl]`` is the ``ceil((n+1)*lvl)``-th smallest absolute residual, so
the interval ``age ± q`` has finite-sample coverage >= ``lvl``. The quantiles feed
the RES score's confidence discount at inference.
"""

from __future__ import annotations

import warnings

import numpy as np

from cellfate.common.schemas import ConformalParams


def fit_conformal(abs_residuals, levels, default_q: float = 1e3,
                  sigma_scale: float = 1.0,
                  sigma_scale_mode: str = "ensemble") -> ConformalParams:
    """Build ConformalParams from absolute age residuals.

    From Stage 1b these residuals are CROSS-DONOR (inner leave-one-donor-out), not the
    in-distribution calib residuals -- see ``xdonor_calib``.

    If there are no age-valid calibration cells (e.g. an all-cancer dataset), a
    large ``default_q`` is stored so the bundle stays valid and intervals are
    conservatively wide. Split conformal can only guarantee coverage up to
    ``n/(n+1)``; a requested level above that cannot be met with this many
    calibration points, so we warn and return the widest finite interval (the
    largest residual) rather than silently under-covering.

    ``sigma_scale`` rides along in this artefact rather than in one of its own: it is
    calibration, it is fitted at the same moment, and adding a defaulted field keeps every
    pre-Stage-1b bundle loadable.
    """
    r = np.asarray(abs_residuals, dtype=np.float64)
    n = r.size
    rs = np.sort(r) if n else None
    levels = [float(lvl) for lvl in levels]
    q: dict[str, float] = {}
    for lvl in levels:
        if n == 0:
            q[str(lvl)] = float(default_q)
            continue
        rank = int(np.ceil((n + 1) * lvl))
        if rank > n:
            warnings.warn(
                f"conformal level {lvl} needs > {n} calibration points to guarantee "
                f"coverage (max attainable ~{n/(n+1):.3f}); using the widest finite interval",
                stacklevel=2,
            )
        k = min(max(rank, 1), n)
        q[str(lvl)] = float(rs[k - 1])
    return ConformalParams(levels=levels, q=q, sigma_scale=float(sigma_scale),
                           sigma_scale_mode=str(sigma_scale_mode))


def coverage(abs_residuals, q: float) -> float:
    """Empirical fraction of residuals within ``q`` (sanity metric)."""
    r = np.asarray(abs_residuals, dtype=np.float64)
    return float((r <= q).mean()) if r.size else float("nan")
