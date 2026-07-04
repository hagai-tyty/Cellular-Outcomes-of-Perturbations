"""Quality control (Document 2, S5). Pure numpy; operates on a RawChunk."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sources import RawChunk


@dataclass
class QCConfig:
    min_genes: int = 200
    max_mito_frac: float = 0.10
    max_counts: float | None = None  # None disables the upper count cap


def compute_qc_metrics(counts: np.ndarray, genes: list[str]) -> dict[str, np.ndarray]:
    """Per-cell total counts, number of expressed genes, and mito fraction."""
    counts = np.asarray(counts, dtype=np.float64)
    total = counts.sum(axis=1)
    n_genes = (counts > 0).sum(axis=1)
    mt_idx = [i for i, g in enumerate(genes) if g.upper().startswith("MT-")]
    mito = counts[:, mt_idx].sum(axis=1) / np.maximum(total, 1.0) if mt_idx else np.zeros(len(total))
    return {"total": total, "n_genes": n_genes, "mito_frac": mito}


def qc_mask(counts: np.ndarray, genes: list[str], cfg: QCConfig) -> np.ndarray:
    """Boolean keep-mask of cells passing all QC thresholds."""
    m = compute_qc_metrics(counts, genes)
    keep = (m["n_genes"] >= cfg.min_genes) & (m["mito_frac"] <= cfg.max_mito_frac)
    if cfg.max_counts is not None:
        keep &= m["total"] <= cfg.max_counts
    return keep


def apply_qc(chunk: RawChunk, cfg: QCConfig) -> RawChunk:
    """Return a new RawChunk with only the cells that pass QC."""
    keep = qc_mask(chunk.counts, chunk.genes, cfg)
    return RawChunk(
        chunk_id=chunk.chunk_id,
        source=chunk.source,
        counts=chunk.counts[keep],
        genes=chunk.genes,
        obs=chunk.obs.loc[keep].reset_index(drop=True),
    )
