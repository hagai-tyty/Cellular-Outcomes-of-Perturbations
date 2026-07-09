"""Cross-modality harmonization (architecture spec, committed build).

Phase 1 -- Dataset-Aware Control Anchoring: per-dataset control statistics
(pseudobulk mu/sigma for single-cell, real replicates for bulk), a variance
floor, an admissible-gene mask, and the per-cell Z-transform.

Phase 2 -- The Gill Projection: reverse the Z-score into the reference (bulk)
biological scale before the frozen clock. Requires NO fitted parameters and NO
donor-age anchors, so control-relative Delta-Age is batch-immune by construction:

    Delta-Age = (X_pert,scaled - X_ctrl,scaled) . (sigma_ref . w)

All statistics are estimated from TRAINING control cells only (the caller must
exclude the held-out donor); the transform is then applied to every cell.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

EPS = 1e-6
DEFAULT_EXPR_FLOOR = 0.1  # min mean control expression (log1p-CP10k) to be admissible
MIN_REPLICATES = 3        # fewest control observations that still define a sigma


@dataclass
class _DatasetStats:
    genes: list[str]      # aligned to the common gene space G
    mu: np.ndarray        # (|G|,)
    sigma: np.ndarray     # (|G|,) already variance-floored


class Harmonizer:
    """Holds per-dataset control mu/sigma on a shared gene space and applies the
    Z-transform + Gill Projection. Fit on training controls only."""

    def __init__(self, genes: list[str], stats: dict[str, _DatasetStats], ref_dataset: str):
        self.genes = genes
        self._stats = stats
        self.ref_dataset = ref_dataset
        self._idx = {g: i for i, g in enumerate(genes)}
        if ref_dataset not in stats:
            raise ValueError(f"reference dataset {ref_dataset!r} not among {list(stats)}")

    # ---------------------------------------------------------------- fit ---- #
    @classmethod
    def fit(
        cls,
        controls: dict[str, list[tuple[np.ndarray, list[str]]]],
        ref_dataset: str = "gill_bulk",
        expr_floor: float = DEFAULT_EXPR_FLOOR,
    ) -> Harmonizer:
        """Estimate control statistics per dataset.

        ``controls[dataset_id]`` is a list of (norm, genes) blocks of **training
        control cells** for that dataset (already excludes the held-out donor).
        """
        if not controls:
            raise ValueError("no control data provided to Harmonizer.fit")

        # 1. pool each dataset's controls into one (cells x genes) matrix on a per-
        #    dataset gene order, and record which genes clear the expression floor.
        pooled: dict[str, tuple[np.ndarray, list[str]]] = {}
        admissible: dict[str, set[str]] = {}
        for ds, blocks in controls.items():
            if not blocks:
                raise ValueError(f"dataset {ds!r} has no control blocks")
            genes0 = blocks[0][1]
            mats = []
            for arr, genes_b in blocks:
                if genes_b == genes0:
                    mats.append(np.asarray(arr, dtype=np.float64))
                else:  # align a block with a different gene order
                    a = np.zeros((arr.shape[0], len(genes0)), dtype=np.float64)
                    src = {g: j for j, g in enumerate(genes_b)}
                    for j, g in enumerate(genes0):
                        s = src.get(g)
                        if s is not None:
                            a[:, j] = arr[:, s]
                    mats.append(a)
            M = np.vstack(mats)
            pooled[ds] = (M, genes0)
            mean_expr = M.mean(axis=0)
            admissible[ds] = {genes0[i] for i in range(len(genes0)) if mean_expr[i] >= expr_floor}

        # 2. common gene space G = genes admissible in EVERY dataset (sorted, stable)
        common = set.intersection(*admissible.values())
        genes_G = sorted(common)
        if not genes_G:
            raise ValueError("no genes are admissible in all datasets (expr_floor too high?)")

        # 3. per-dataset mu/sigma on G, computed over INDIVIDUAL control observations
        #    (single cells for scRNA-seq, real samples for bulk). NOTE: sigma must be
        #    the individual-observation spread, NOT a pseudobulk spread -- pseudobulk
        #    sigma is ~sqrt(bucket_size) too small and inflates per-cell Z-scores. The
        #    zero-inflation explosion (near-constant genes -> sigma~0) is instead handled
        #    by the variance floor below.
        stats: dict[str, _DatasetStats] = {}
        for ds, (M, genes0) in pooled.items():
            col = {g: i for i, g in enumerate(genes0)}
            MG = M[:, [col[g] for g in genes_G]]              # (obs, |G|) on G
            if MG.shape[0] < MIN_REPLICATES:
                raise ValueError(
                    f"dataset {ds!r} has {MG.shape[0]} control observations "
                    f"(< {MIN_REPLICATES}); sigma undefined")
            mu = MG.mean(axis=0)
            sigma = MG.std(axis=0)
            floor = float(np.median(sigma))                   # median over admissible genes
            sigma = np.maximum(sigma, floor)
            stats[ds] = _DatasetStats(genes=genes_G, mu=mu, sigma=sigma)
        return cls(genes_G, stats, ref_dataset)

    # ---------------------------------------------------------- transform ---- #
    def transform(self, norm: np.ndarray, genes: list[str], dataset_id: str) -> np.ndarray:
        """Z-score a chunk against its own dataset's control stats. Returns
        (N, |G|) on the harmonizer's gene space ``self.genes``."""
        st = self._stats.get(dataset_id)
        if st is None:
            raise KeyError(f"no harmonization stats for dataset {dataset_id!r}")
        aligned = self._align(norm, genes)                    # (N, |G|)
        return (aligned - st.mu) / (st.sigma + EPS)

    def project_to_clock(self, x_scaled: np.ndarray) -> np.ndarray:
        """Gill Projection: reverse the Z-score into the reference (bulk) scale so
        the frozen clock reads biologically-weighted features. Uses the reference
        dataset's mu/sigma for ALL cells."""
        ref = self._stats[self.ref_dataset]
        return x_scaled * ref.sigma + ref.mu

    def _align(self, norm: np.ndarray, genes: list[str]) -> np.ndarray:
        norm = np.asarray(norm, dtype=np.float64)
        src = {g: i for i, g in enumerate(genes)}
        out = np.zeros((norm.shape[0], len(self.genes)), dtype=np.float64)
        for j, g in enumerate(self.genes):
            s = src.get(g)
            if s is not None:
                out[:, j] = norm[:, s]
        return out

    # --------------------------------------------------------- persistence -- #
    def to_json(self, path: str | Path) -> None:
        payload = {
            "genes": self.genes,
            "ref_dataset": self.ref_dataset,
            "stats": {ds: {"mu": st.mu.tolist(), "sigma": st.sigma.tolist()}
                      for ds, st in self._stats.items()},
        }
        Path(path).write_text(json.dumps(payload), encoding="utf-8")

    @classmethod
    def from_json(cls, path: str | Path) -> Harmonizer:
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        genes = d["genes"]
        stats = {ds: _DatasetStats(genes=genes, mu=np.array(s["mu"]), sigma=np.array(s["sigma"]))
                 for ds, s in d["stats"].items()}
        return cls(genes, stats, d["ref_dataset"])

