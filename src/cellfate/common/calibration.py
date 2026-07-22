"""Calibration transforms shared by training (which FITS them) and inference (which APPLIES).

Lives in ``common`` because both layers need it and the dependency only runs one way:
``cellfate.inference`` imports nothing from ``cellfate.training`` (see the ``inference`` package
docstring). Fitting stays in ``training.calibrate``; only the transform is shared, so the
function that produces a calibration and the function that consumes it can never drift apart.
"""

from __future__ import annotations

import numpy as np

# Clamp for logit(). LOAD-BEARING, not cosmetic: this model's P(safe) saturates (fate PR-AUC
# ~1.0), and float32 values that round to exactly 1.0 would give an infinite logit and a NaN
# probability. The clamp also bounds the sensitivity below -- d(logit)/dp = 1/(p(1-p)) reaches
# ~1e7 at p = 1 - 1e-7, where a single float32 ulp moves the logit by ~0.9.
EPS = 1e-6


def platt_safe(p_safe: np.ndarray, a: float, b: float) -> np.ndarray:
    """``sigmoid(a*logit(P(safe)) + b)`` -- the calibrated safe probability.

    ⚠ NOISE AMPLIFICATION. Working in logit space multiplies input perturbations by roughly the
    slope ``a``. Measured on a trained bundle: torch's batch-size-dependent CPU kernels move the
    raw ensemble probability by ~9e-08, and this map takes that to ~5e-07 at a ~ 8. Harmless at
    these magnitudes, and the reason batch-agreement is asserted to a tolerance rather than
    bit-exactly (tests/test_inference.py), but it scales with ``a`` -- worth remembering if the
    slope bound in ``fit_platt_binary`` is ever raised.
    """
    p = np.clip(np.asarray(p_safe, dtype=np.float64), EPS, 1.0 - EPS)
    z = np.log(p / (1.0 - p))
    return 1.0 / (1.0 + np.exp(-(a * z + b)))


def apply_platt(probs, a: float, b: float, safe_idx: int) -> np.ndarray:
    """Recalibrate ``P(safe)`` and rescale the other classes so rows still sum to 1.

    The loss/death RATIO is preserved -- only the safe-vs-unsafe boundary moves -- so
    ``P_loss`` stays a meaningful input to RES instead of being redistributed arbitrarily.

    With ``a > 0`` the map is strictly increasing in ``P(safe)``, so every ranking metric
    (fate PR-AUC, ROC-AUC) is mathematically unchanged.
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
