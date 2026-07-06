"""``cellfate.evaluation`` -- baselines, metrics, external validation, and the
automated acceptance gates that turn the Document-1 success criteria into a
falsifiable pass/fail verdict (Document 5).

Consumes ``bundle/`` + ``splits/`` + ``shards/``; produces ``reports/``.
"""

from __future__ import annotations

from .baselines import (
    BASELINE_NAMES,
    Estimator,
    KNNFingerprint,
    MeanBaseline,
    ModelEstimator,
    PredictControl,
    RidgeLinear,
    UOnly,
    XOnly,
    make_baselines,
)
from .data import SplitData, gather_split
from .evaluate_cli import EvalConfig, check_gates, evaluate, run_external
from .external_validation import validate_against_methylation, validate_oskm_holdout
from .metrics import (
    brier,
    coverage,
    des_pds,
    ece,
    edistance_to_control,
    energy_distance,
    per_class_auroc,
    per_class_prauc,
    precision_at_k,
    ranking_metrics,
    regression_metrics,
)
from .regimes import REGIMES, describe, iter_regimes
from .report import leaderboard, reliability_diagram, write_report, write_summary

__all__ = [
    "EvalConfig", "evaluate", "check_gates", "run_external",
    "SplitData", "gather_split",
    "Estimator", "ModelEstimator", "make_baselines", "BASELINE_NAMES",
    "MeanBaseline", "RidgeLinear", "XOnly", "UOnly", "KNNFingerprint", "PredictControl",
    "per_class_auroc", "per_class_prauc", "brier", "ece", "regression_metrics",
    "coverage", "ranking_metrics", "precision_at_k", "des_pds",
    "energy_distance", "edistance_to_control",
    "validate_against_methylation", "validate_oskm_holdout",
    "REGIMES", "describe", "iter_regimes",
    "leaderboard", "reliability_diagram", "write_report", "write_summary",
]
