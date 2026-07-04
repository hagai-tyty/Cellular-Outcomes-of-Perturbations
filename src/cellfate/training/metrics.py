"""Small evaluation helpers for the training report (Document 3, S5)."""

from __future__ import annotations

import numpy as np


def soft_nll(probs: np.ndarray, target: np.ndarray) -> float:
    """Mean soft-label negative log-likelihood: -sum_c target_c log p_c."""
    p = np.clip(np.asarray(probs, dtype=np.float64), 1e-12, 1.0)
    return float(-(np.asarray(target, dtype=np.float64) * np.log(p)).sum(axis=1).mean())


def ece(probs: np.ndarray, labels: np.ndarray, n_bins: int = 10) -> float:
    """Expected calibration error of the top-1 prediction (``labels`` = int classes)."""
    probs = np.asarray(probs, dtype=np.float64)
    labels = np.asarray(labels)
    conf = probs.max(axis=1)
    correct = (probs.argmax(axis=1) == labels).astype(np.float64)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    n = len(labels)
    err = 0.0
    for i in range(n_bins):
        m = (conf > edges[i]) & (conf <= edges[i + 1])
        if m.any():
            err += abs(correct[m].mean() - conf[m].mean()) * m.sum() / n
    return float(err)
