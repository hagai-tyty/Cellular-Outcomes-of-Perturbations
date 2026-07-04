"""``cellfate.models`` -- the network + losses (Document 3).

Defines the multi-task CellFate-Rx network (shared encoders + trunk, a
safe/loss/death classification head, and a ΔAge head) and its training losses
(class-balanced focal, masked Huber, Kendall multi-task weighting). Requires
PyTorch; ``cellfate.common`` stays torch-free, so importing this package is what
pulls torch in.
"""

from __future__ import annotations

from .encoders import CellEncoder, ChemEncoder, TFEncoder
from .heads import AgeHead, ClassificationHead
from .losses import (
    MultiTaskLoss,
    class_balanced_weights,
    focal_loss,
    huber_age_loss,
)
from .network import CellFateNet, mc_dropout_predict

__all__ = [
    "CellFateNet", "mc_dropout_predict",
    "CellEncoder", "ChemEncoder", "TFEncoder", "ClassificationHead", "AgeHead",
    "class_balanced_weights", "focal_loss", "huber_age_loss", "MultiTaskLoss",
]
