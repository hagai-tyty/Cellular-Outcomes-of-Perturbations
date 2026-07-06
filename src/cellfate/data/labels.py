"""Fate labels (Document 2, S9).

Two ways to turn expression into an (N, 3) probability vector over
(safe, loss, death):

* ``soft_labels`` -- a plain temperature softmax over precomputed signature
  scores (kept for tests / simple use).
* ``fate_labels`` -- the **control-relative** labeller used by the ETL. It scores
  each cell on a pluripotency, somatic-identity, and apoptosis program using the
  ``LABEL_HOLDOUT`` marker genes, then thresholds each cell **relative to its own
  cell line's control distribution** (a distributional threshold, not an absolute
  one) -- so it is robust to cross-dataset batch effects, mirroring the
  control-relative ΔAge strategy. Identity loss = gaining pluripotency AND losing
  somatic identity; death = elevated apoptosis; safe = neither.
"""

from __future__ import annotations

import numpy as np

from cellfate.common import constants as C

from ._stats import softmax
from .signatures import score_one


def soft_labels(sig_scores: np.ndarray, tau: float = 1.0) -> np.ndarray:
    """Convert (N, 3) signature scores into (N, 3) probabilities summing to 1."""
    if tau <= 0:
        raise ValueError("tau must be > 0")
    return softmax(np.asarray(sig_scores, dtype=np.float64) / tau, axis=1)


def _control_relative_z(score: np.ndarray, lines: np.ndarray, is_ctrl: np.ndarray) -> np.ndarray:
    """z-score each cell's program score against its cell line's CONTROL cells.

    Lines with too few controls fall back to the global control distribution
    (or, absent any controls, the global distribution). This is what makes the
    threshold distributional rather than an absolute count.
    """
    score = np.asarray(score, dtype=np.float64)
    z = np.zeros_like(score)
    ctrl = score[is_ctrl]
    g_mu, g_sd = (float(ctrl.mean()), float(ctrl.std())) if ctrl.size else (
        float(score.mean()), float(score.std()))
    g_sd = g_sd if g_sd > 1e-8 else 1.0
    for line in np.unique(lines):
        m = lines == line
        cm = m & is_ctrl
        if int(cm.sum()) >= 2:
            mu = float(score[cm].mean())
            sd = float(score[cm].std())
            sd = sd if sd > 1e-8 else g_sd
        else:
            mu, sd = g_mu, g_sd
        z[m] = (score[m] - mu) / sd
    return z


def fate_labels(
    norm_expr: np.ndarray,
    genes: list[str],
    obs,
    tau: float = 1.0,
    baseline_safe: float = 2.0,
) -> np.ndarray:
    """Control-relative soft labels over (safe, loss, death).

    ``obs`` must have ``cell_line`` and ``is_control`` columns. Uses the
    ``LABEL_HOLDOUT`` markers (held out of the model panel) on the full profile.
    """
    if tau <= 0:
        raise ValueError("tau must be > 0")
    plu = score_one(norm_expr, genes, C.DEFAULT_SIGNATURES["loss"])   # pluripotency
    som = score_one(norm_expr, genes, C.DEFAULT_SIGNATURES["safe"])   # somatic identity
    apo = score_one(norm_expr, genes, C.DEFAULT_SIGNATURES["death"])  # apoptosis
    lines = obs["cell_line"].to_numpy()
    is_ctrl = obs["is_control"].to_numpy().astype(bool)
    z_plu = _control_relative_z(plu, lines, is_ctrl)
    z_som = _control_relative_z(som, lines, is_ctrl)
    z_apo = _control_relative_z(apo, lines, is_ctrl)
    n = norm_expr.shape[0]
    # columns MUST follow CLASSES = (safe, loss, death)
    logit_safe = np.full(n, float(baseline_safe))
    logit_loss = z_plu - z_som           # pluripotency up AND identity down
    logit_death = z_apo
    logits = np.stack([logit_safe, logit_loss, logit_death], axis=1)
    return softmax(logits / tau, axis=1)
