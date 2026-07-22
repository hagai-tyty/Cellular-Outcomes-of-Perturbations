"""Calibration transforms shared by training (which FITS them) and inference (which APPLIES).

Lives in ``common`` because both layers need it and the dependency only runs one way:
``cellfate.inference`` imports nothing from ``cellfate.training`` (see the ``inference`` package
docstring). Fitting stays in ``training.calibrate``; only the transform is shared, so the
function that produces a calibration and the function that consumes it can never drift apart.
"""

from __future__ import annotations

import numpy as np

# Clamp for logit(). LOAD-BEARING: this model's P(safe) saturates (fate PR-AUC ~1.0), and values
# that reach exactly 1.0 or 0.0 would give an infinite logit and a NaN probability.
#
# SIZED so it does not collapse resolution the input actually carries. Probabilities arrive from
# a float32 softmax, whose ulp just below 1.0 is ~6e-08; a clamp coarser than that maps distinct
# inputs onto one value and CREATES TIES. Measured at the previous 1e-6: four of eight
# float32-representable values near 1.0 collapsed onto a single number, and ties in a 400-cell
# sample rose from 72 to ~220. Those ties happened to fall within one class, so PR-AUC/ROC-AUC
# were unaffected -- but that is a property of the data, not a guarantee. 1e-9 sits two orders
# below the float32 ulp, so every representable input except exact 0/1 survives distinct.
EPS = 1e-9


def _stable_sigmoid(x: np.ndarray) -> np.ndarray:
    """1/(1+exp(-x)) without overflow. A tighter EPS admits |x| of a few hundred."""
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def platt_safe(p_safe: np.ndarray, a: float, b: float) -> np.ndarray:
    """``sigmoid(a*logit(P(safe)) + b)`` -- the calibrated safe probability.

    ⚠ NOISE AMPLIFICATION. Working in logit space multiplies input perturbations by roughly the
    slope ``a``. Measured on a trained bundle: torch's batch-size-dependent CPU kernels move the
    raw ensemble probability by ~9e-08, and this map takes that to ~5e-07 at a ~ 8. Harmless at
    these magnitudes, and the reason batch-agreement is asserted to a tolerance rather than
    bit-exactly (tests/test_inference.py), but it scales with ``a`` -- worth remembering if the
    slope bound in ``fit_platt_binary`` is ever raised.

    ⚠ MONOTONE, NOT INJECTIVE. With ``a > 0`` the map is strictly increasing in ``P(safe)``, so
    it never REORDERS cells. It can still merge them: the clamp above ties exact 0/1 to their
    neighbours, and the float32 cast on the way out can tie values the map left distinct. Rank
    metrics are therefore stable in practice but not guaranteed bit-identical -- do not claim
    they are.
    """
    p = np.clip(np.asarray(p_safe, dtype=np.float64), EPS, 1.0 - EPS)
    return _stable_sigmoid(a * np.log(p / (1.0 - p)) + b)


def apply_platt(probs, a: float, b: float, safe_idx: int) -> np.ndarray:
    """Recalibrate ``P(safe)`` and rescale the other classes so rows still sum to 1.

    The loss/death RATIO is preserved -- only the safe-vs-unsafe boundary moves -- so
    ``P_loss`` stays a meaningful input to RES instead of being redistributed arbitrarily.

    With ``a > 0`` the map is strictly increasing in ``P(safe)``, so every ranking metric
    (fate PR-AUC, ROC-AUC) is stable -- though see the injectivity note on `platt_safe`.
    """
    p = np.array(probs, dtype=np.float64, copy=True)
    flat = p.ndim == 1
    if flat:
        p = p[None, :]
    s_cal = platt_safe(p[:, safe_idx], a, b)

    other = [i for i in range(p.shape[1]) if i != safe_idx]
    rest = p[:, other]
    denom = rest.sum(axis=1, keepdims=True)
    # all mass already on `safe` -> nothing to redistribute; split the remainder evenly
    even = np.full_like(rest, 1.0 / max(len(other), 1))
    share = np.where(denom > 0, rest / np.where(denom > 0, denom, 1.0), even)

    out = np.empty_like(p)
    out[:, safe_idx] = s_cal
    out[:, other] = share * (1.0 - s_cal)[:, None]
    return out[0] if flat else out
