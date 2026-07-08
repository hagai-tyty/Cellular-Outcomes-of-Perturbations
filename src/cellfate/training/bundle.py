"""Deployment-bundle assembly (Document 3, S7).

Writes everything inference (Document 4) needs into ``out/bundle/`` and then calls
the foundation's ``assert_bundle_complete`` so a half-written bundle can never be
shipped: the K ensemble members, the calibration temperature, the conformal
quantiles, the RES parameters, the OOD reference, the (unchanged) scalers, a YAML
copy of the config, and the metadata.
"""

from __future__ import annotations

from omegaconf import OmegaConf

from cellfate.common import constants as C
from cellfate.common import io
from cellfate.common.io import ArtifactPaths
from cellfate.common.scalers import Scalers
from cellfate.common.schemas import (
    BundleMeta,
    ConformalParams,
    ResParams,
    TemperatureParams,
)

from .ood import save_ood


def assemble_bundle(
    out_root: str,
    members,
    temperature: TemperatureParams,
    conformal: ConformalParams,
    res: ResParams,
    ood: dict,
    *,
    gene_panel_hash: str,
    scalers_src: str,
    model_cfg: dict,
    train_cfg: dict,
    deps_hash: str | None = None,
    config_hash: str | None = None,
) -> ArtifactPaths:
    """Write the bundle and assert it is complete; returns the ArtifactPaths."""
    paths = ArtifactPaths.of(out_root)
    for d in (paths.bundle_dir, paths.bundle_members_dir, paths.bundle_ood_dir):
        d.mkdir(parents=True, exist_ok=True)

    for i, model in enumerate(members):
        model.save_member(paths.bundle_members_dir / f"member_{i}.pt")

    io.save_temperature(paths, temperature)
    io.save_conformal(paths, conformal)
    io.save_res_params(paths, res)
    Scalers.load(scalers_src).save(paths.bundle_scalers_file)   # shipped unchanged
    save_ood(paths, ood)

    (paths.bundle_dir / C.BUNDLE_CONFIG_FILENAME).write_text(
        OmegaConf.to_yaml(OmegaConf.create({"model": model_cfg, "train": train_cfg})),
        encoding="utf-8",
    )

    meta = BundleMeta(
        n_members=len(members),
        gene_panel_hash=gene_panel_hash,
        conformal_levels=list(conformal.levels),
        deps_hash=deps_hash,
        config_hash=config_hash,
    )
    io.save_bundle_meta(paths, meta)
    io.assert_bundle_complete(paths)   # raises BundleError if anything is missing
    return paths
