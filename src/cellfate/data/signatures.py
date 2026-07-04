"""Signature scoring (Document 2, S8).

Each outcome class (``safe`` / ``loss`` / ``death``) has a marker gene set.
A cell's signature score is the mean standardised expression of that set --
high ``loss`` means pluripotency/dedifferentiation markers are up, high ``death``
means apoptosis markers are up, high ``safe`` means somatic-identity markers are
retained. Scores are computed on the FULL normalised matrix (not the HVG panel)
so signature genes are never lost.
"""

from __future__ import annotations

import numpy as np

from cellfate.common import constants as C

from ._stats import standardize


def score_one(norm_expr: np.ndarray, genes: list[str], gene_set: tuple[str, ...]) -> np.ndarray:
    """Mean standardised expression over the genes of ``gene_set`` present."""
    idx = [genes.index(g) for g in gene_set if g in genes]
    if not idx:
        return np.zeros(norm_expr.shape[0], dtype=np.float64)
    return standardize(norm_expr[:, idx]).mean(axis=1)


def signature_scores(
    norm_expr: np.ndarray,
    genes: list[str],
    signatures: dict[str, tuple[str, ...]] | None = None,
) -> np.ndarray:
    """Return an (N, 3) score matrix in the fixed ``CLASSES`` order.

    Column order is guaranteed to be (safe, loss, death) -- the same order as
    ``Sample.sig_scores`` and ``Sample.y_cls``.
    """
    sigs = signatures or C.DEFAULT_SIGNATURES
    cols = [score_one(norm_expr, genes, sigs[c]) for c in C.CLASSES]
    return np.stack(cols, axis=1)
