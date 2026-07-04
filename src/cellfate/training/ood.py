"""Out-of-distribution detection (Document 3, S6).

Fits a single Gaussian to the trunk features of the training cells and uses the
Mahalanobis distance as the OOD score (Lee et al. 2018). The mean, precision
matrix, and a high-quantile distance threshold are stored in the bundle so
inference can flag inputs that lie outside the training manifold (e.g. a wrong
gene panel or an unseen assay) instead of silently extrapolating.
"""

from __future__ import annotations

import numpy as np

from cellfate.common.io import ArtifactPaths

_OOD_FILENAME = "mahalanobis.npz"


def fit_ood(features, ridge: float = 1e-3, quantile: float = 0.99) -> dict:
    """Fit mean + (ridge-regularised) precision; threshold at a train-distance quantile."""
    f = np.asarray(features, dtype=np.float64)
    mu = f.mean(axis=0)
    centered = f - mu
    cov = centered.T @ centered / max(len(f) - 1, 1)
    cov += ridge * np.eye(f.shape[1])
    precision = np.linalg.inv(cov)
    d2 = np.einsum("ni,ij,nj->n", centered, precision, centered)
    return {
        "mean": mu.astype(np.float32),
        "precision": precision.astype(np.float32),
        "threshold": float(np.quantile(d2, quantile)),
        "dim": int(f.shape[1]),
    }


def mahalanobis(features, ood: dict) -> np.ndarray:
    """Squared Mahalanobis distance of each row to the training Gaussian."""
    centered = np.asarray(features, dtype=np.float64) - ood["mean"]
    return np.einsum("ni,ij,nj->n", centered, ood["precision"], centered)


def save_ood(paths: ArtifactPaths, ood: dict) -> None:
    paths.bundle_ood_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        paths.bundle_ood_dir / _OOD_FILENAME,
        mean=ood["mean"], precision=ood["precision"],
        threshold=np.array(ood["threshold"]), dim=np.array(ood["dim"]),
    )


def load_ood(paths: ArtifactPaths) -> dict:
    d = np.load(paths.bundle_ood_dir / _OOD_FILENAME)
    return {"mean": d["mean"], "precision": d["precision"],
            "threshold": float(d["threshold"]), "dim": int(d["dim"])}
