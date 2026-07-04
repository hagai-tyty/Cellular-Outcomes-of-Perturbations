"""Output heads (Document 3, S1).

Two heads sit on the shared trunk latent: a classification head over the fixed
safe/loss/death classes (returns logits), and an age head returning a single
ΔAge point estimate. Epistemic uncertainty comes from the deep ensemble +
MC-dropout, not from a predicted variance, so the age head is a plain scalar.
"""

from __future__ import annotations

import torch
from torch import nn


class ClassificationHead(nn.Module):
    """Trunk latent (B, d) -> class logits (B, n_classes)."""

    def __init__(self, d: int, n_classes: int) -> None:
        super().__init__()
        self.fc = nn.Linear(d, n_classes)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.fc(z)


class AgeHead(nn.Module):
    """Trunk latent (B, d) -> ΔAge (B,)."""

    def __init__(self, d: int) -> None:
        super().__init__()
        self.fc = nn.Linear(d, 1)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.fc(z).squeeze(-1)
