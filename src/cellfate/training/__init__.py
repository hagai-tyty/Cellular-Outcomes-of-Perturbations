"""``cellfate.training`` -- training, ensembling, calibration, bundling (Document 3).

Trains the deep ensemble, fits temperature + conformal calibration and the OOD
reference, and writes the contract-valid deployment bundle that inference
(Document 4) consumes. Requires PyTorch.
"""

from __future__ import annotations

from .bundle import assemble_bundle
from .calibrate import fit_temperature, has_class_variation
from .conformal import coverage, fit_conformal
from .dataset import load_split_tensors, loader
from .metrics import ece, soft_nll
from .ood import fit_ood, load_ood, mahalanobis, save_ood
from .train import (
    class_mass,
    ensemble_age,
    ensemble_logits,
    member_outputs,
    train_ensemble,
    train_member,
)
from .train_model import TrainConfig, run
from .xdonor_calib import XDonorStats, crossdonor_stats, n_train_donors, sigma_scale_factor

__all__ = [
    "TrainConfig", "run",
    "load_split_tensors", "loader",
    "train_member", "train_ensemble", "member_outputs",
    "ensemble_logits", "ensemble_age", "class_mass",
    "fit_temperature", "has_class_variation", "fit_conformal", "coverage",
    "fit_ood", "mahalanobis", "save_ood", "load_ood",
    "assemble_bundle", "ece", "soft_nll",
    "crossdonor_stats", "sigma_scale_factor", "n_train_donors", "XDonorStats",
]
