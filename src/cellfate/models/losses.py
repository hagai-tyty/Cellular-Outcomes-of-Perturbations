"""Loss functions for the multi-task head (Document 3, S2).

* class-balanced **focal** loss for the safe/loss/death head (soft-label aware),
* masked **Huber** loss for the ΔAge head (only on cells with a valid age),
* **Kendall & Gal** uncertainty weighting that learns to balance the two tasks.
"""

from __future__ import annotations

import numpy as np
import torch
import torch.nn.functional as F
from torch import nn


def class_balanced_weights(class_mass: np.ndarray, beta: float = 0.999) -> np.ndarray:
    """Effective-number class weights (Cui et al. 2019), normalised to mean 1.

    ``class_mass`` is the per-class total soft-label mass over the train split
    (sum of the y_cls columns), so rare outcomes (death) are up-weighted.
    """
    mass = np.asarray(class_mass, dtype=np.float64)
    eff = 1.0 - np.power(beta, np.maximum(mass, 1.0))
    w = (1.0 - beta) / np.maximum(eff, 1e-12)
    return (w / w.sum() * len(mass)).astype(np.float32)


def focal_loss(logits: torch.Tensor, target: torch.Tensor,
               class_weights: torch.Tensor, gamma: float = 2.0) -> torch.Tensor:
    """Soft-label focal loss. ``target`` is a probability vector per row (sums to 1).

    Reduces to standard focal loss when ``target`` is one-hot; with soft labels it
    is the focal-modulated cross-entropy summed over the label distribution.

    Note: focal loss is fundamentally a hard-label heuristic — with soft targets its
    minimum is not exactly at ``p == target`` (the modulating ``(1-p)^gamma`` term
    shifts it), so it trades a small calibration bias for its intended focus on the
    rare, hard death class. The post-hoc temperature scaling step corrects the
    residual miscalibration, and ``gamma=0`` recovers exact (calibrated) soft CE.
    """
    logp = F.log_softmax(logits, dim=1)
    p = logp.exp()
    focal = (1.0 - p).clamp(min=0.0, max=1.0) ** gamma
    w = class_weights.unsqueeze(0)
    return -(w * focal * target * logp).sum(dim=1).mean()


def huber_age_loss(age_pred: torch.Tensor, age_true: torch.Tensor,
                   mask: torch.Tensor, delta: float = 2.0) -> torch.Tensor:
    """Huber (smooth-L1) loss on ΔAge, averaged over masked (age-valid) cells only.

    Returns a differentiable zero (tied to the graph) when no cell in the batch
    has a valid age, so masked-only batches contribute nothing to the age task.
    """
    m = mask.bool()
    if not torch.any(m):
        return age_pred.sum() * 0.0
    return F.huber_loss(age_pred[m], age_true[m], delta=delta)


class MultiTaskLoss(nn.Module):
    """Kendall & Gal (2018) homoscedastic uncertainty weighting of two tasks.

    total = exp(-s_cls)*L_cls + 0.5*exp(-s_age)*L_age + 0.5*(s_cls + s_age),
    where s_* = log(variance) are learned, so the network balances the tasks
    instead of using a hand-tuned weight.

    The asymmetry is intentional and follows the Gaussian/softmax derivation:
    the regression (Gaussian) likelihood contributes ``1/(2*sigma^2)`` -> the
    ``0.5`` factor on ``exp(-s_age)``, whereas the classification (softmax)
    likelihood contributes ``1/sigma^2`` -> a full ``exp(-s_cls)`` with no 0.5.
    """

    def __init__(self) -> None:
        super().__init__()
        self.log_var_cls = nn.Parameter(torch.zeros(()))
        self.log_var_age = nn.Parameter(torch.zeros(()))

    def forward(self, l_cls: torch.Tensor, l_age: torch.Tensor) -> torch.Tensor:
        c = torch.exp(-self.log_var_cls) * l_cls + 0.5 * self.log_var_cls
        a = 0.5 * torch.exp(-self.log_var_age) * l_age + 0.5 * self.log_var_age
        return c + a
