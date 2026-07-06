"""``cellfate.common`` -- the foundation / contract layer.

Everything every other package depends on lives here: constants, schemas,
artefact IO, normalisation, determinism, logging, and the progress ledger.
No module in this package imports from ``cellfate.data/models/training/
inference/evaluation`` -- dependencies flow one way, into ``common``.
"""

from __future__ import annotations

from . import constants, errors, io, logging, scalers, schemas, seeding
from .constants import (
    CLASS_TO_IDX,
    CLASSES,
    IDX_TO_CLASS,
    N_CLASSES,
    SCHEMA_VERSION,
    Modality,
    Regime,
    Split,
)
from .errors import (
    BundleError,
    CellFateError,
    ConfigError,
    ContractViolation,
    DataSourceError,
    GenePanelMismatch,
    NotImplementedInFoundation,
    SchemaError,
    ShardIOError,
)
from .io import ArtifactPaths
from .panel import GenePanel
from .progress import ProgressTracker
from .scalers import Scalers
from .schemas import (
    BundleMeta,
    ConformalParams,
    ManifestRow,
    ProgressState,
    ResParams,
    Sample,
    ScalerParams,
    TemperatureParams,
)
from .seeding import set_global_seed

__all__ = [
    # submodules
    "constants", "errors", "io", "logging", "schemas", "scalers", "seeding",
    # constants / enums
    "SCHEMA_VERSION", "CLASSES", "N_CLASSES", "CLASS_TO_IDX", "IDX_TO_CLASS",
    "Modality", "Split", "Regime",
    # schemas
    "Sample", "ManifestRow", "ScalerParams", "ResParams", "ConformalParams",
    "TemperatureParams", "ProgressState", "BundleMeta",
    # helpers
    "ArtifactPaths", "GenePanel", "Scalers", "ProgressTracker", "set_global_seed",
    # errors
    "CellFateError", "SchemaError", "ContractViolation", "GenePanelMismatch",
    "ShardIOError", "BundleError", "ConfigError", "DataSourceError",
    "NotImplementedInFoundation",
]
