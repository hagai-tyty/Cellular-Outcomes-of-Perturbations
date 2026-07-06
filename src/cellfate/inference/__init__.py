"""``cellfate.inference`` -- serving, OOD gating, conformal intervals and the
Rejuvenation Efficacy Score (Document 4).

Loads a deployment ``bundle/`` and turns a (cell, perturbation) query into a
calibrated, safety-gated decision. Purely a consumer of the training bundle: it
imports nothing from ``cellfate.training``.
"""

from __future__ import annotations

from .conformal import interval, intervals
from .encode import encode_batch
from .ood import OODDetector
from .predictor import Predictor, enable_mc_dropout
from .res import (
    APPROVED,
    REJECTED_NO_REJUVENATION,
    REJECTED_OOD,
    REJECTED_UNSAFE,
    compute_res,
    compute_res_batch,
)
from .schema import Request, Response
from .service import (
    build_response,
    create_app,
    predict_one,
    score_requests,
    score_shard,
)

__all__ = [
    "Request", "Response", "Predictor", "enable_mc_dropout", "encode_batch",
    "OODDetector", "interval", "intervals",
    "compute_res", "compute_res_batch",
    "APPROVED", "REJECTED_OOD", "REJECTED_UNSAFE", "REJECTED_NO_REJUVENATION",
    "build_response", "predict_one", "score_requests", "score_shard", "create_app",
]
