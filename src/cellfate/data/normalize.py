"""Normalisation + the frozen gene panel (Document 2, S6).

Per-cell normalisation (library-size + log1p) and projection onto the frozen
:class:`GenePanel` are implemented in numpy so the pipeline does not hard-depend
on scanpy. ``GenePanel`` itself is a shared contract and is re-exported from
``cellfate.common``. Reading real ``.h5ad`` atlases uses anndata/scanpy (lazy,
in the source connectors).
"""

from __future__ import annotations

import numpy as np

from cellfate.common import constants as C
from cellfate.common.panel import GenePanel  # re-export

__all__ = ["GenePanel", "normalize_counts", "fit_gene_panel", "to_panel_matrix"]


def normalize_counts(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """Library-size normalise to ``target_sum`` counts per cell, then log1p.

    This is the standard CP10k + log1p transform. The result is what is stored
    as ``Sample.X`` (per-gene z-scoring is applied later by ``Scalers``).
    """
    counts = np.asarray(counts, dtype=np.float64)
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(counts / lib * target_sum).astype(np.float32)


def fit_gene_panel(
    norm_expr: np.ndarray,
    genes: list[str],
    n_top: int = C.DEFAULT_N_GENES,
    must_include: tuple[str, ...] = (),
    must_exclude: tuple[str, ...] = (),
) -> GenePanel:
    """Select the top-``n_top`` most variable genes from normalised data, ONCE.

    ``must_include`` genes are guaranteed to be in the panel. ``must_exclude``
    genes (e.g. the fate-label markers, ``LABEL_HOLDOUT``) are guaranteed to be
    kept OUT -- so the model cannot read its own label off its input. The
    resulting panel is frozen and committed; its hash is written into every
    downstream artefact.
    """
    norm_expr = np.asarray(norm_expr, dtype=np.float64)
    if norm_expr.shape[1] != len(genes):
        raise ValueError("norm_expr columns must match len(genes)")
    exclude = set(must_exclude)
    var = norm_expr.var(axis=0)
    order = list(np.argsort(-var))
    chosen: list[str] = []
    seen: set[str] = set()
    for g in must_include:
        if g in genes and g not in seen and g not in exclude:
            chosen.append(g)
            seen.add(g)
    for i in order:
        if len(chosen) >= n_top:
            break
        g = genes[i]
        if g not in seen and g not in exclude:
            chosen.append(g)
            seen.add(g)
    return GenePanel(chosen)


def to_panel_matrix(norm_expr: np.ndarray, genes: list[str], panel: GenePanel) -> np.ndarray:
    """Project normalised expression onto the panel's gene order.

    Genes present in the panel but absent from ``genes`` are filled with 0
    (a gene the assay did not measure contributes nothing). Output is
    (N, len(panel)) float32 in panel order -- the model input ``X``.
    """
    norm_expr = np.asarray(norm_expr, dtype=np.float32)
    idx = {g: i for i, g in enumerate(genes)}
    out = np.zeros((norm_expr.shape[0], len(panel)), dtype=np.float32)
    for j, g in enumerate(panel.genes):
        col = idx.get(g)
        if col is not None:
            out[:, j] = norm_expr[:, col]
    return out
