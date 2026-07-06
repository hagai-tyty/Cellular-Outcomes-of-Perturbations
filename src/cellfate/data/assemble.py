"""Assemble validated :class:`Sample` objects from per-chunk arrays (S13)."""

from __future__ import annotations

import numpy as np

from cellfate.common.constants import Modality
from cellfate.common.schemas import Sample


def assemble_samples(
    *,
    cell_ids: list[str],
    x_panel: np.ndarray,        # (N, G) normalised, panel order
    fingerprints: np.ndarray,   # (N, n_bits) uint8 (chem)
    dose_time: np.ndarray,      # (N, 2)
    y_cls: np.ndarray,          # (N, 3) sums to 1
    y_age: np.ndarray,          # (N,) ΔAge (used only where age_mask)
    age_mask: np.ndarray,       # (N,) bool
    sig_scores: np.ndarray,     # (N, 3)
    cell_line: list[str],
    pert_id: list[str],
    scaffold_id: list[str],
    source: str,
    modality: Modality = Modality.CHEM,
    tf_emb: np.ndarray | None = None,   # (N, |TF_VOCAB|) multi-hot (tf modality)
) -> list[Sample]:
    """Build one validated Sample per cell. Raises if any row violates the schema."""
    n = len(cell_ids)
    samples: list[Sample] = []
    for i in range(n):
        masked = bool(age_mask[i])
        samples.append(
            Sample(
                cell_id=cell_ids[i],
                X=x_panel[i].astype(float).tolist(),
                u_modality=modality,
                u_chem_fp=fingerprints[i].astype(int).tolist() if modality is Modality.CHEM else None,
                u_tf_emb=(tf_emb[i].astype(float).tolist()
                          if (modality is Modality.TF and tf_emb is not None) else None),
                dose_time=dose_time[i].astype(float).tolist(),
                y_cls=y_cls[i].astype(float).tolist(),
                y_age=(float(y_age[i]) if masked else None),
                age_mask=masked,
                sig_scores=sig_scores[i].astype(float).tolist(),
                cell_line=cell_line[i],
                pert_id=pert_id[i],
                scaffold_id=scaffold_id[i],
                source=source,
            )
        )
    return samples
