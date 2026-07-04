"""Pydantic data contracts. These are the *single source of truth* for the
shape of every object that crosses a package boundary or is written to disk.

Changing any of these is a breaking change: bump ``SCHEMA_VERSION`` in
``constants`` and migrate existing artefacts.
"""

from __future__ import annotations

import math
import time

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from . import constants as C
from .constants import Modality

# Re-export so callers can do ``from cellfate.common.schemas import Modality``.
__all__ = [
    "Modality",
    "Sample",
    "ManifestRow",
    "ScalerParams",
    "ResParams",
    "ConformalParams",
    "TemperatureParams",
    "ChunkDone",
    "ChunkFailed",
    "ProgressState",
    "BundleMeta",
]


def _all_finite(xs: list[float]) -> bool:
    return all(math.isfinite(x) for x in xs)


# --------------------------------------------------------------------------- #
# The atomic unit: one perturbed cell                                         #
# --------------------------------------------------------------------------- #
class Sample(BaseModel):
    """One perturbed cell -- the row written into a shard (Document 1, S6.2)."""

    model_config = ConfigDict(extra="forbid")

    cell_id: str
    X: list[float]                      # normalised HVG expression (model input)
    u_modality: Modality
    u_chem_fp: list[int] | None = None  # Morgan fingerprint bits (chem)
    u_gene_emb: list[float] | None = None
    u_tf_emb: list[float] | None = None
    dose_time: list[float]              # [log10(dose_uM), log(time_h)]
    y_cls: list[float]                  # soft label over CLASSES, sums to 1
    y_age: float | None                 # dAge in years; None when masked
    age_mask: bool                      # True iff y_age is a valid label
    sig_scores: list[float]             # raw Safe/Loss/Death signature scores
    cell_line: str
    pert_id: str
    scaffold_id: str | None = None
    source: str

    # -- field-level checks ------------------------------------------------- #
    @field_validator("X")
    @classmethod
    def _x_ok(cls, v: list[float]) -> list[float]:
        if not v:
            raise ValueError("X must be non-empty")
        if not _all_finite(v):
            raise ValueError("X contains non-finite values")
        return v

    @field_validator("dose_time")
    @classmethod
    def _dt_ok(cls, v: list[float]) -> list[float]:
        if len(v) != C.N_DOSE_TIME:
            raise ValueError(f"dose_time must have length {C.N_DOSE_TIME}, got {len(v)}")
        if not _all_finite(v):
            raise ValueError("dose_time contains non-finite values")
        return v

    @field_validator("y_cls")
    @classmethod
    def _ycls_ok(cls, v: list[float]) -> list[float]:
        if len(v) != C.N_CLASSES:
            raise ValueError(f"y_cls must have length {C.N_CLASSES}, got {len(v)}")
        if any(p < -1e-6 for p in v):
            raise ValueError("y_cls has negative probabilities")
        if abs(sum(v) - 1.0) > 1e-3:
            raise ValueError(f"y_cls must sum to 1, got {sum(v):.4f}")
        return v

    @field_validator("sig_scores")
    @classmethod
    def _sig_ok(cls, v: list[float]) -> list[float]:
        if len(v) != C.N_CLASSES:
            raise ValueError(f"sig_scores must have length {C.N_CLASSES}, got {len(v)}")
        return v

    # -- cross-field checks ------------------------------------------------- #
    @model_validator(mode="after")
    def _consistency(self) -> Sample:
        # modality <-> descriptor: the selected modality's field is required, and
        # the chemical fingerprint must be absent for non-chemical modalities.
        if self.u_modality is Modality.CHEM:
            if self.u_chem_fp is None:
                raise ValueError("chem modality requires u_chem_fp")
            if len(self.u_chem_fp) != C.N_FINGERPRINT_BITS:
                raise ValueError(
                    f"u_chem_fp must have length {C.N_FINGERPRINT_BITS}, got {len(self.u_chem_fp)}"
                )
            if any(b not in (0, 1) for b in self.u_chem_fp):
                raise ValueError("u_chem_fp must be a 0/1 bit vector")
        elif self.u_modality is Modality.GENETIC:
            if self.u_gene_emb is None:
                raise ValueError("genetic modality requires u_gene_emb")
            if self.u_chem_fp is not None:
                raise ValueError("genetic modality must not carry u_chem_fp")
        elif self.u_modality is Modality.TF:
            if self.u_tf_emb is None:
                raise ValueError("tf modality requires u_tf_emb")
            if self.u_chem_fp is not None:
                raise ValueError("tf modality must not carry u_chem_fp")

        # age label <-> mask (both directions): a valid mask requires a finite
        # label, and a masked-out row must not carry a spurious finite label.
        if self.age_mask:
            if self.y_age is None or not math.isfinite(self.y_age):
                raise ValueError("age_mask=True requires a finite y_age")
        elif self.y_age is not None and math.isfinite(self.y_age):
            raise ValueError("age_mask=False requires y_age to be None or NaN")
        # scaffold present for chem (used by the leave-drug-out split)
        if self.u_modality is Modality.CHEM and not self.scaffold_id:
            raise ValueError("chem samples require a scaffold_id for the scaffold split")
        return self


# --------------------------------------------------------------------------- #
# Manifest: grouping keys only (so splitting never loads X)                    #
# --------------------------------------------------------------------------- #
class ManifestRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cell_id: str
    cell_line: str
    pert_id: str
    scaffold_id: str | None
    source: str
    age_mask: bool
    shard_id: str        # sanitised file stem of the shard holding this sample
    row_idx: int         # row position within that shard

    @classmethod
    def from_sample(cls, s: Sample, shard_id: str, row_idx: int) -> ManifestRow:
        return cls(
            cell_id=s.cell_id, cell_line=s.cell_line, pert_id=s.pert_id,
            scaffold_id=s.scaffold_id, source=s.source, age_mask=s.age_mask,
            shard_id=shard_id, row_idx=row_idx,
        )


# --------------------------------------------------------------------------- #
# Normalisation parameters (fit on TRAIN only; shipped to inference)          #
# --------------------------------------------------------------------------- #
class ScalerParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    x_mean: list[float]
    x_std: list[float]
    dt_mean: list[float]
    dt_std: list[float]
    proliferation_coef: list[float] = Field(..., description="(a, b): dAge ~ a*cc + b")
    gene_panel_hash: str

    @model_validator(mode="after")
    def _ok(self) -> ScalerParams:
        if len(self.x_mean) != len(self.x_std):
            raise ValueError("x_mean and x_std must have equal length")
        if len(self.dt_mean) != C.N_DOSE_TIME or len(self.dt_std) != C.N_DOSE_TIME:
            raise ValueError(f"dt_mean/dt_std must have length {C.N_DOSE_TIME}")
        if len(self.proliferation_coef) != 2:
            raise ValueError("proliferation_coef must be (a, b)")
        if any(s < 0 for s in self.x_std) or any(s < 0 for s in self.dt_std):
            raise ValueError("standard deviations must be non-negative")
        return self


# --------------------------------------------------------------------------- #
# RES parameters (serialised into the bundle so scoring is reproducible)       #
# --------------------------------------------------------------------------- #
class ResParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    tau_safe: float = 0.85   # safety floor
    w: float = 0.03          # floor width
    k: float = 2.0           # safety-dominant exponent (>=1)
    kappa: float = 5.0       # concavity scale for rejuvenation (years)
    z_conf: float = 1.0      # confidence multiplier on sigma_age (>=0)
    lam: float = 0.0         # optional cancer-risk penalty weight (>=0)

    @model_validator(mode="after")
    def _ok(self) -> ResParams:
        if not 0.0 < self.tau_safe < 1.0:
            raise ValueError("tau_safe must be in (0, 1)")
        if self.w <= 0:
            raise ValueError("w must be > 0")
        if self.k < 1:
            raise ValueError("k must be >= 1 (safety must dominate)")
        if self.kappa <= 0:
            raise ValueError("kappa must be > 0")
        if self.z_conf < 0 or self.lam < 0:
            raise ValueError("z_conf and lam must be >= 0")
        return self


class ConformalParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    levels: list[float]
    q: dict[str, float]      # keyed by str(level), e.g. {"0.9": 1.83}

    @model_validator(mode="after")
    def _ok(self) -> ConformalParams:
        for lvl in self.levels:
            if not 0.0 < lvl < 1.0:
                raise ValueError(f"conformal level {lvl} must be in (0, 1)")
            if str(lvl) not in self.q:
                raise ValueError(f"missing quantile for level {lvl} (key {str(lvl)!r})")
        return self


class TemperatureParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    temperature: float = 1.0

    @field_validator("temperature")
    @classmethod
    def _pos(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("temperature must be > 0")
        return v


# --------------------------------------------------------------------------- #
# Fault-tolerance ledger (used by the ETL ProgressTracker)                     #
# --------------------------------------------------------------------------- #
class ChunkDone(BaseModel):
    n: int
    ts: float = Field(default_factory=time.time)


class ChunkFailed(BaseModel):
    err: str
    ts: float = Field(default_factory=time.time)


class ProgressState(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: str = C.SCHEMA_VERSION
    done: dict[str, ChunkDone] = Field(default_factory=dict)
    failed: dict[str, ChunkFailed] = Field(default_factory=dict)


# --------------------------------------------------------------------------- #
# Deployment-bundle metadata (the sole input to inference)                    #
# --------------------------------------------------------------------------- #
class BundleMeta(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: str = C.SCHEMA_VERSION
    n_members: int
    gene_panel_hash: str
    classes: list[str] = Field(default_factory=lambda: list(C.CLASSES))
    conformal_levels: list[float] = Field(default_factory=list)
    created_at: float = Field(default_factory=time.time)
    deps_hash: str | None = None
    config_hash: str | None = None

    @model_validator(mode="after")
    def _ok(self) -> BundleMeta:
        if self.classes != list(C.CLASSES):
            raise ValueError(f"classes must equal {list(C.CLASSES)}, got {self.classes}")
        if self.n_members < 1:
            raise ValueError("n_members must be >= 1")
        return self
