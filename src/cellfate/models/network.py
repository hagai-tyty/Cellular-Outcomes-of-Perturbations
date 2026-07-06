"""The CellFate-Rx network (Document 3, S1).

Shared encoders feed a shared trunk; two heads produce the calibrated
safe/loss/death distribution and the ΔAge estimate. The trunk latent is also
returned as the feature used for Mahalanobis OOD detection. ``mc_dropout_predict``
runs the stochastic forward passes used for epistemic uncertainty at inference.
"""

from __future__ import annotations

from pathlib import Path

import torch
import torch.nn.functional as F
from torch import nn

from cellfate.common import constants as C
from cellfate.common.constants import N_DOSE_TIME, N_FINGERPRINT_BITS

from .encoders import ChemEncoder, _mlp_block
from .heads import AgeHead, ClassificationHead


class CellFateNet(nn.Module):
    """(expression, fingerprint, dose/time) -> (class logits, ΔAge, trunk feature)."""

    def __init__(
        self,
        g: int,
        n_fp: int = N_FINGERPRINT_BITS,
        n_dt: int = N_DOSE_TIME,
        d_cell: int = 256,
        d_u: int = 256,
        latent_dim: int = 256,
        p_drop: float = 0.2,
        n_classes: int = C.N_CLASSES,
        pert_kind: str = "chem",
        n_pert: int | None = None,
    ) -> None:
        super().__init__()
        # perturbation input width: fingerprint bits (chem) or TF vocab (tf)
        if n_pert is None:
            n_pert = len(C.TF_VOCAB) if pert_kind == "tf" else n_fp
        self.arch = dict(g=g, n_fp=n_fp, n_dt=n_dt, d_cell=d_cell, d_u=d_u,
                         latent_dim=latent_dim, p_drop=p_drop, n_classes=n_classes,
                         pert_kind=pert_kind, n_pert=n_pert)
        from .encoders import CellEncoder, TFEncoder  # local import keeps file order clean
        self.cell = CellEncoder(g, d_cell, p_drop)
        self.pert_kind = pert_kind
        self.pert = (TFEncoder(n_pert, n_dt, d_u, p_drop) if pert_kind == "tf"
                     else ChemEncoder(n_pert, n_dt, d_u, p_drop))
        self.trunk = nn.Sequential(
            _mlp_block(d_cell + d_u, latent_dim, p_drop),
            _mlp_block(latent_dim, latent_dim, p_drop),
        )
        self.cls_head = ClassificationHead(latent_dim, n_classes)
        self.age_head = AgeHead(latent_dim)

    def forward(self, x: torch.Tensor, u: torch.Tensor, dose_time: torch.Tensor):
        # u is the perturbation vector: a fingerprint (chem) or a TF multi-hot (tf)
        z = self.trunk(torch.cat([self.cell(x), self.pert(u, dose_time)], dim=1))
        return self.cls_head(z), self.age_head(z), z

    # -- (de)serialisation -------------------------------------------------- #
    def save_member(self, path: str | Path) -> None:
        torch.save({"state_dict": self.state_dict(), "arch": self.arch}, path)

    @classmethod
    def load_member(cls, path: str | Path, map_location: str = "cpu") -> CellFateNet:
        blob = torch.load(path, map_location=map_location, weights_only=False)
        model = cls(**blob["arch"])
        model.load_state_dict(blob["state_dict"])
        model.eval()
        return model


@torch.no_grad()
def mc_dropout_predict(model: CellFateNet, x, fp, dose_time, n_samples: int = 20):
    """Run ``n_samples`` stochastic forward passes (dropout ON) for epistemic uncertainty.

    Returns stacked class probabilities (T, B, C) and ages (T, B). The model is
    left in eval mode on exit.
    """
    was_training = model.training
    model.train()  # enable dropout; LayerNorm has no running stats so this is safe
    probs, ages = [], []
    for _ in range(n_samples):
        logits, age, _ = model(x, fp, dose_time)
        probs.append(F.softmax(logits, dim=1))
        ages.append(age)
    if not was_training:
        model.eval()
    return torch.stack(probs), torch.stack(ages)
