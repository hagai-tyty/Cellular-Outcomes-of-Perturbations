"""Training orchestrator (Document 3, S7).

Loads the Document-2 dataset, trains the ensemble, calibrates (temperature +
conformal), fits the OOD reference, and writes a complete, contract-valid
deployment bundle plus a metrics report. ``run`` is plain-Python and testable;
``cli`` is the Hydra entry point.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch

from cellfate.common import constants as C
from cellfate.common import io
from cellfate.common.errors import ConfigError
from cellfate.common.io import ArtifactPaths
from cellfate.common.logging import get_logger, log_event
from cellfate.common.scalers import Scalers
from cellfate.common.schemas import ResParams
from cellfate.models import CellFateNet

from .bundle import assemble_bundle
from .calibrate import fit_temperature
from .conformal import coverage, fit_conformal
from .dataset import AM_I, FP_I, YA_I, YC_I, load_split_tensors
from .metrics import ece, soft_nll
from .ood import fit_ood
from .train import ensemble_age, ensemble_logits, member_outputs, train_ensemble

log = get_logger("cellfate.training")


@dataclass
class TrainConfig:
    dataset_dir: str
    out: str | None = None              # bundle root; defaults to dataset_dir
    regime: str = "scaffold"
    # model
    d_cell: int = 256
    d_u: int = 256
    latent_dim: int = 256
    p_drop: float = 0.2
    # optimisation
    lr: float = 3e-4
    wd: float = 1e-5
    epochs: int = 60
    patience: int = 8
    min_delta: float = 1e-4
    clip: float = 1.0
    batch_size: int = 512
    ensemble_size: int = 5
    base_seed: int = 0
    # losses
    focal_gamma: float = 2.0
    class_weight_beta: float = 0.999
    huber_delta: float = 2.0
    # calibration / device / RES
    conformal_levels: tuple[float, ...] = (0.90,)
    device: str = "cpu"
    res: dict = field(default_factory=dict)


def _resolve_device(requested: str) -> str:
    if requested == "cuda" and not torch.cuda.is_available():
        warnings.warn("cuda requested but unavailable; using cpu", stacklevel=2)
        return "cpu"
    return requested


def run(cfg: TrainConfig) -> dict:
    """Train + calibrate + bundle; returns a summary/metrics dict."""
    device = _resolve_device(cfg.device)
    data_paths = ArtifactPaths.of(cfg.dataset_dir)
    out_root = cfg.out or cfg.dataset_dir
    scalers = Scalers.load(data_paths.scalers_file)
    g = len(scalers.params.x_mean)

    train_ds = load_split_tensors(data_paths, scalers, cfg.regime, "train")
    val_ds = load_split_tensors(data_paths, scalers, cfg.regime, "val")
    calib_ds = load_split_tensors(data_paths, scalers, cfg.regime, "calib")
    if len(train_ds) == 0:
        raise ConfigError(f"empty train split for regime {cfg.regime!r}")
    log_event(log, "data.loaded", g=g, n_train=len(train_ds),
              n_val=len(val_ds), n_calib=len(calib_ds))

    # perturbation modality: TF-cocktail vs chemical fingerprint. The width comes
    # from the loaded tensor; the kind from the shards' u_modality. Both go into
    # the model config, so the member arch (and inference) get the right encoder.
    def _detect_pert_kind() -> str:
        for shard in sorted(data_paths.shards_dir.glob("*.parquet")):
            mod = io.shard_to_numpy(io.read_shard(shard))["u_modality"]
            return "tf" if (len(mod) and str(mod[0]) == "tf") else "chem"
        return "chem"

    pert_kind = _detect_pert_kind()
    n_pert = int(train_ds.tensors[FP_I].shape[1])
    model_cfg = dict(g=g, d_cell=cfg.d_cell, d_u=cfg.d_u,
                     latent_dim=cfg.latent_dim, p_drop=cfg.p_drop,
                     pert_kind=pert_kind, n_pert=n_pert)

    def make_model() -> CellFateNet:
        return CellFateNet(**model_cfg)

    members, val_losses = train_ensemble(make_model, train_ds, val_ds, cfg, device)

    # -- temperature scaling on val (fallback to calib) --------------------- #
    cal_ds = val_ds if len(val_ds) else calib_ds
    if len(cal_ds):
        cal_logits = ensemble_logits(members, cal_ds, device).numpy()
        cal_target = cal_ds.tensors[YC_I].numpy()
        temperature = fit_temperature(cal_logits, cal_target)
    else:
        from cellfate.common.schemas import TemperatureParams
        cal_logits = cal_target = None
        temperature = TemperatureParams(temperature=1.0)

    # -- conformal on calib age residuals (ensemble age) -------------------- #
    if len(calib_ds):
        age_pred = ensemble_age(members, calib_ds, device).numpy()
        ya = calib_ds.tensors[YA_I].numpy()
        am = calib_ds.tensors[AM_I].numpy().astype(bool)
        abs_res = np.abs(age_pred[am] - ya[am])
    else:
        abs_res = np.array([])
    conformal = fit_conformal(abs_res, cfg.conformal_levels)

    # -- OOD reference from train trunk features ---------------------------- #
    train_feats = member_outputs(members[0], train_ds, device)[2].numpy()
    ood = fit_ood(train_feats)

    res_params = ResParams(**cfg.res)
    config_hash = io.hash_config({"model": model_cfg, "train": vars(cfg)})

    paths = assemble_bundle(
        out_root, members, temperature, conformal, res_params, ood,
        gene_panel_hash=scalers.params.gene_panel_hash,
        scalers_src=str(data_paths.scalers_file),
        model_cfg=model_cfg, train_cfg=_jsonable(vars(cfg)),
        deps_hash=io.deps_hash(), config_hash=config_hash,
    )

    metrics = _report(members, val_losses, train_ds, val_ds, calib_ds,
                      cal_logits, cal_target, temperature, conformal, abs_res)
    io.write_json(paths.bundle_dir / C.BUNDLE_METRICS_FILENAME, metrics)
    log_event(log, "bundle.done", bundle=str(paths.bundle_dir),
              n_members=len(members), temperature=round(temperature.temperature, 4))
    return {"bundle": str(paths.bundle_dir), **metrics}


def _jsonable(d: dict) -> dict:
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in d.items()}


def _report(members, val_losses, train_ds, val_ds, calib_ds,
            cal_logits, cal_target, temperature, conformal, abs_res) -> dict:
    out = {
        "n_members": len(members),
        "n_train": len(train_ds),
        "n_val": len(val_ds),
        "n_calib": len(calib_ds),
        "val_loss_mean": float(np.mean(val_losses)),
        "val_loss_per_member": [float(v) for v in val_losses],
        "temperature": temperature.temperature,
        "conformal_q": dict(conformal.q),
    }
    if cal_logits is not None:
        from scipy.special import softmax  # available via the data-stage deps
        labels = cal_target.argmax(axis=1)
        p_before = softmax(cal_logits, axis=1)
        p_after = softmax(cal_logits / temperature.temperature, axis=1)
        out["nll_before_temp"] = soft_nll(p_before, cal_target)
        out["nll_after_temp"] = soft_nll(p_after, cal_target)
        out["ece_before_temp"] = ece(p_before, labels)
        out["ece_after_temp"] = ece(p_after, labels)
    if abs_res.size:
        lvl = conformal.levels[0]
        out["conformal_coverage_calib"] = coverage(abs_res, conformal.q[str(lvl)])
    return out


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def _config_from_omegaconf(cfg) -> TrainConfig:
    m, t = cfg.model, cfg.train
    return TrainConfig(
        dataset_dir=cfg.data.out,
        out=cfg.get("bundle_out", None) or cfg.data.out,
        regime=str(cfg.data.splits.get("primary", "scaffold")),
        d_cell=int(m.d_cell), d_u=int(m.d_u), latent_dim=int(m.latent_dim),
        p_drop=float(m.p_drop),
        lr=float(t.lr), wd=float(t.wd), epochs=int(t.epochs), patience=int(t.patience),
        min_delta=float(t.min_delta), clip=float(t.clip), batch_size=int(t.batch_size),
        ensemble_size=int(t.ensemble_size), base_seed=int(t.base_seed),
        focal_gamma=float(t.focal.gamma), class_weight_beta=float(t.focal.class_weight_beta),
        huber_delta=float(t.huber_delta),
        conformal_levels=(float(t.conformal_level),),
        device=str(cfg.get("device", "cpu")),
        res=dict(cfg.get("res", {})),
    )


def cli() -> None:  # pragma: no cover - exercised in production via Hydra
    try:
        import hydra
        from omegaconf import DictConfig
    except ImportError as exc:
        raise ConfigError("hydra-core/omegaconf required for the CLI") from exc

    config_dir = str(Path(__file__).resolve().parents[3] / "configs")

    @hydra.main(version_base=None, config_path=config_dir, config_name="config")
    def _main(cfg: DictConfig) -> None:
        summary = run(_config_from_omegaconf(cfg))
        log_event(log, "cli.done", bundle=summary["bundle"], n_members=summary["n_members"])

    _main()


if __name__ == "__main__":
    cli()
