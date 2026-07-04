"""``cellfate.data`` -- the ETL package (Document 2).

Turns raw single-cell perturbation atlases into validated, sharded training
artefacts that satisfy the ``cellfate.common`` contracts (shards of ``Sample``
rows + manifest + per-regime splits + fitted scalers). Heavy dependencies
(scanpy, rdkit, Hugging Face) are imported lazily inside the functions that need
them, so importing this package is cheap and the scientific core is testable
without them.
"""

from __future__ import annotations

from .aging import AgingClock, LinearClock, delta_age, recenter_on_controls
from .assemble import assemble_samples
from .build_dataset import DataConfig, process_chunk, run
from .chunking import CellChunk, ProgressTracker, plan_all
from .labels import soft_labels
from .normalize import GenePanel, fit_gene_panel, normalize_counts, to_panel_matrix
from .perturbation import (
    bemis_murcko_scaffold,
    encode_dose_time,
    encode_fingerprints,
    hashed_fingerprint,
    morgan_fingerprint,
    resolve_scaffolds,
)
from .proliferation import cell_cycle_score, deconfound_age, fit_deconfounder
from .qc import QCConfig, apply_qc, compute_qc_metrics, qc_mask
from .signatures import score_one, signature_scores
from .sources import (
    SOURCE_REGISTRY,
    DataSource,
    GillReprogrammingSource,
    GSE242423SingleCellSource,
    RawChunk,
    ReprogrammingSource,
    SciplexSource,
    SyntheticSource,
    TahoeSource,
)
from .splits import cell_line_split, make_splits, scaffold_split

__all__ = [
    # orchestration
    "DataConfig", "run", "process_chunk",
    # sources / chunking
    "DataSource", "RawChunk", "SyntheticSource", "TahoeSource", "SciplexSource",
    "ReprogrammingSource", "GillReprogrammingSource", "GSE242423SingleCellSource",
    "SOURCE_REGISTRY", "CellChunk", "ProgressTracker", "plan_all",
    # qc / normalise / panel
    "QCConfig", "apply_qc", "compute_qc_metrics", "qc_mask",
    "GenePanel", "fit_gene_panel", "normalize_counts", "to_panel_matrix",
    # perturbation
    "morgan_fingerprint", "hashed_fingerprint", "bemis_murcko_scaffold",
    "encode_fingerprints", "encode_dose_time", "resolve_scaffolds",
    # labels / signatures / aging / proliferation
    "score_one", "signature_scores", "soft_labels",
    "AgingClock", "LinearClock", "delta_age", "recenter_on_controls",
    "cell_cycle_score", "fit_deconfounder", "deconfound_age",
    # splits / assemble
    "scaffold_split", "cell_line_split", "make_splits", "assemble_samples",
]
