"""Out-of-distribution gate for inference (Document 4, S3).

A real detector, not dropout variance: the (squared) Mahalanobis distance of the
query's latent trunk feature to the training-set Gaussian stored in the bundle.
A query is OOD when the distance exceeds the calibrated threshold; the model then
declines to extrapolate (status ``REJECTED_OOD``) instead of guessing.

The bundle stores a single-Gaussian reference (mean, precision, squared-distance
threshold); this loader reads it directly, so ``inference`` does not import from
``training``. Chemical-space novelty (max-Tanimoto to the training fingerprints)
is a documented extension point: it activates only if the bundle carries a
fingerprint index, and otherwise contributes nothing, leaving the latent
Mahalanobis gate -- which already flags a novel compound whose *effect* is unlike
anything seen in training.
"""

from __future__ import annotations

import numpy as np

from cellfate.common.io import ArtifactPaths


class OODDetector:
    def __init__(self, paths: ArtifactPaths):
        npz_files = sorted(paths.bundle_ood_dir.glob("*.npz"))
        if not npz_files:
            raise FileNotFoundError(f"no OOD reference (*.npz) in {paths.bundle_ood_dir}")
        d = np.load(npz_files[0])
        self.mean = np.asarray(d["mean"], dtype=np.float64)
        self.precision = np.asarray(d["precision"], dtype=np.float64)
        self.threshold = float(d["threshold"])
        self.dim = int(d["dim"])

    def distances(self, Z) -> np.ndarray:
        """Squared Mahalanobis distance of each row of ``Z`` (N, dim) to the mean."""
        c = np.asarray(Z, dtype=np.float64).reshape(-1, self.dim) - self.mean
        return np.einsum("ni,ij,nj->n", c, self.precision, c)

    def is_in_distribution(self, z) -> bool:
        return bool(self.distances(z)[0] <= self.threshold)

    def in_distribution_mask(self, Z) -> np.ndarray:
        return self.distances(Z) <= self.threshold
