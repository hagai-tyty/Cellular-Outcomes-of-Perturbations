"""The Rejuvenation Efficacy Score (Document 4, S5; Technical Architecture 4.3).

A safety-dominant, uncertainty-aware score in [0, 1):

    RES = Φ(S) · S**k · g(R_eff) · exp(−λ·P_loss)

with Φ(S) = sigmoid((S − τ_safe)/w) a smooth safety floor (no cliff),
R_eff = max(0, −(µ_age + z_conf·σ_age)) crediting only *confident* rejuvenation
(the upper age bound must be negative), and g(R_eff) = R_eff/(R_eff + κ) a concave,
bounded rejuvenation term. Because uncertainty is subtracted before credit is
given, a confident modest rejuvenator can outrank an uncertain large one -- the
opposite of a naive S·|Δ| score.
"""

from __future__ import annotations

import numpy as np

from cellfate.common.schemas import ResParams

APPROVED = "APPROVED"
REJECTED_OOD = "REJECTED_OOD"
REJECTED_UNSAFE = "REJECTED_UNSAFE"
REJECTED_NO_REJUVENATION = "REJECTED_NO_REJUVENATION"


def _sigmoid(x):
    # clip to the saturated region (|x|>500 is 0/1 to full float precision); this
    # is exact for real RES inputs (|x| < ~30) and avoids overflow warnings.
    x = np.clip(np.asarray(x, dtype=np.float64), -500.0, 500.0)
    return 1.0 / (1.0 + np.exp(-x))


def compute_res(S, P_loss, mu_age, sigma_age, in_dist, p: ResParams):
    """Scalar RES + status for one query."""
    if not in_dist:
        return 0.0, REJECTED_OOD
    R_eff = max(0.0, -(float(mu_age) + p.z_conf * float(sigma_age)))
    g = R_eff / (R_eff + p.kappa)
    phi = float(_sigmoid((float(S) - p.tau_safe) / p.w))
    res = float(phi * (float(S) ** p.k) * g * np.exp(-p.lam * float(P_loss)))
    if S < p.tau_safe - 3.0 * p.w:
        status = REJECTED_UNSAFE
    elif R_eff == 0.0:
        status = REJECTED_NO_REJUVENATION
    else:
        status = APPROVED
    return res, status


def compute_res_batch(S, P_loss, mu_age, sigma_age, in_dist, p: ResParams):
    """Vectorised RES + status over arrays (for batch ranking)."""
    S = np.asarray(S, dtype=np.float64)
    P_loss = np.asarray(P_loss, dtype=np.float64)
    mu = np.asarray(mu_age, dtype=np.float64)
    sig = np.asarray(sigma_age, dtype=np.float64)
    ind = np.asarray(in_dist, dtype=bool)

    R_eff = np.maximum(0.0, -(mu + p.z_conf * sig))
    g = R_eff / (R_eff + p.kappa)
    phi = _sigmoid((S - p.tau_safe) / p.w)
    res = phi * (S ** p.k) * g * np.exp(-p.lam * P_loss)
    res = np.where(ind, res, 0.0)

    status = np.where(
        ~ind, REJECTED_OOD,
        np.where(S < p.tau_safe - 3.0 * p.w, REJECTED_UNSAFE,
                 np.where(R_eff == 0.0, REJECTED_NO_REJUVENATION, APPROVED)),
    )
    return res, status
