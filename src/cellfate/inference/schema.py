"""Request / Response contract for the inference service (Document 4).

A ``Request`` describes one (cell, perturbation) query. ``X_raw`` is the
panel-restricted, library-size-normalised ``log1p`` expression -- i.e. exactly the
space of the training ``Sample.X`` -- and the bundled scaler applies the final
standardisation. (Library normalisation needs the full transcriptome, which a
single panel vector cannot supply, so it is done upstream; the name ``X_raw``
means "raw input to the model's scaler", not raw counts.)
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from cellfate.common.constants import Modality


class Request(BaseModel):
    model_config = ConfigDict(extra="forbid")

    X_raw: list[float]                 # length G: log-normalised panel expression (Sample.X space)
    u_modality: Modality
    u_descriptor: str | list[float]    # SMILES (chem) or a precomputed 2048-bit fingerprint
    dose_uM: float
    time_h: float


class Response(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: str                        # APPROVED | REJECTED_OOD | REJECTED_UNSAFE | REJECTED_NO_REJUVENATION
    rejuvenation_efficacy_score: float # 10 * RES, in [0, 10)
    p_identity_preserved: float
    p_identity_loss: float
    p_apoptosis: float
    delta_age_mean: float
    delta_age_interval: list[float]    # conformal, guaranteed marginal coverage
    in_distribution: bool
    epistemic_std: float
    predictive_entropy: float
    warning: str | None = None
