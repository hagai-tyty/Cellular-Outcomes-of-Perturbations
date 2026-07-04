"""Input encoders (Document 3, S1).

A cell encoder maps the HVG expression vector to a latent, and a modality-specific
perturbation encoder maps the chemical fingerprint + dose/time to a latent. v1
ships the chemical path; genetic / TF encoders plug in alongside it for v2.
"""

from __future__ import annotations

import torch
from torch import nn


def _mlp_block(d_in: int, d_out: int, p_drop: float) -> nn.Sequential:
    return nn.Sequential(
        nn.Linear(d_in, d_out),
        nn.LayerNorm(d_out),
        nn.GELU(),
        nn.Dropout(p_drop),
    )


class CellEncoder(nn.Module):
    """HVG expression (B, G) -> cell latent (B, d_cell)."""

    def __init__(self, g: int, d_cell: int, p_drop: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            _mlp_block(g, d_cell, p_drop),
            _mlp_block(d_cell, d_cell, p_drop),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ChemEncoder(nn.Module):
    """Morgan fingerprint (B, n_fp) + dose/time (B, n_dt) -> pert latent (B, d_u)."""

    def __init__(self, n_fp: int, n_dt: int, d_u: int, p_drop: float) -> None:
        super().__init__()
        self.fp_net = _mlp_block(n_fp, d_u, p_drop)
        self.head = _mlp_block(d_u + n_dt, d_u, p_drop)

    def forward(self, fp: torch.Tensor, dose_time: torch.Tensor) -> torch.Tensor:
        h = self.fp_net(fp)
        return self.head(torch.cat([h, dose_time], dim=1))


class TFEncoder(nn.Module):
    """TF-cocktail (B, vocab) + dose/time (B, n_dt) -> pert latent (B, d_u).

    A learned **vocabulary embedding**: the first linear layer's weight is one
    embedding column per factor, and the multi-hot (dose-scaled) input selects and
    sums the embeddings of the factors present. So OSKM and OSK share factor
    embeddings -- the encoder learns factor identity and generalises across
    cocktails, unlike an opaque hashed token. Same call signature as ChemEncoder.
    """

    def __init__(self, vocab: int, n_dt: int, d_u: int, p_drop: float) -> None:
        super().__init__()
        self.embed = nn.Linear(vocab, d_u, bias=False)   # per-factor embedding columns
        self.head = _mlp_block(d_u + n_dt, d_u, p_drop)

    def forward(self, factor_vec: torch.Tensor, dose_time: torch.Tensor) -> torch.Tensor:
        h = self.embed(factor_vec)                       # sum of present factors' embeddings x dose
        return self.head(torch.cat([h, dose_time], dim=1))
