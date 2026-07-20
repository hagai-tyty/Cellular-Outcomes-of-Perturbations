"""Load a deployment bundle and produce prediction summaries (Document 4, S2).

The Predictor loads the bundle once (members, scalers, temperature, OOD, conformal
quantile, RES params) and answers "what happens to this cell under this
perturbation?" with batched stochastic passes. Two modes:

* ``ensemble`` (primary, coverage-guaranteed): one forward per deep-ensemble member;
  epistemic σ is the spread across members.
* ``mc_dropout`` (fallback, single member): T dropout passes, run as ONE batched
  forward over the tiled input -- never a Python per-sample loop.

``enable_mc_dropout`` puts *only* dropout sub-modules in train mode; everything
else (here LayerNorm, which is stateless) stays in eval, so no running statistics
are ever corrupted.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import torch

from cellfate.common.constants import SCHEMA_VERSION
from cellfate.common.errors import ConfigError, GenePanelMismatch, SchemaError
from cellfate.common.io import (
    ArtifactPaths,
    assert_bundle_complete,
    load_bundle_meta,
    load_conformal,
    load_res_params,
    load_temperature,
)
from cellfate.common.scalers import Scalers
from cellfate.models import CellFateNet

from .encode import encode_batch
from .ood import OODDetector


def enable_mc_dropout(model: torch.nn.Module) -> None:
    """Put ONLY dropout sub-modules in train mode (everything else stays eval)."""
    model.eval()
    for m in model.modules():
        if isinstance(m, torch.nn.Dropout):
            m.train()


def _resolve_root(bundle_dir) -> Path:
    """Accept the artefact root, the ``bundle/`` dir, or a path containing it."""
    p = Path(bundle_dir)
    if (p / "bundle").is_dir():
        return p
    if p.name == "bundle":
        return p.parent
    return p


class Predictor:
    def __init__(self, bundle_dir, mode: str = "ensemble", T: int = 50,
                 conformal_level: float = 0.90):
        if mode not in ("ensemble", "mc_dropout"):
            raise ConfigError(f"unknown inference mode '{mode}' (use 'ensemble' or 'mc_dropout')")
        self.paths = ArtifactPaths.of(_resolve_root(bundle_dir))
        assert_bundle_complete(self.paths)                      # fail fast, clear error

        self.meta = load_bundle_meta(self.paths)
        if self.meta.schema_version != SCHEMA_VERSION:
            raise SchemaError(
                f"bundle schema_version {self.meta.schema_version!r} != runtime "
                f"{SCHEMA_VERSION!r}; rebuild the bundle with this version of cellfate."
            )
        self.scalers = Scalers.load(self.paths.bundle_scalers_file)
        if self.meta.gene_panel_hash != self.scalers.params.gene_panel_hash:
            raise GenePanelMismatch(
                self.meta.gene_panel_hash, self.scalers.params.gene_panel_hash,
                "bundle meta and its shipped scalers disagree.",
            )
        self.members = [CellFateNet.load_member(p)
                        for p in sorted(self.paths.bundle_members_dir.glob("*.pt"))]
        for m in self.members:
            m.eval()
        if not self.members:
            raise ConfigError("bundle contains no ensemble members")

        self.temperature = load_temperature(self.paths).temperature
        self.ood = OODDetector(self.paths)
        self.res_params = load_res_params(self.paths)

        conf = load_conformal(self.paths)
        key = str(float(conformal_level))
        if key not in conf.q:                                   # fall back to a stored level
            key = str(conf.levels[0])
        self.q = conf.q[key]
        self.conformal_level = float(key)
        self.mode, self.T = mode, int(T)

        # Stage 1b: sigma_age is the ENSEMBLE SPREAD (~2.4 yr) while true out-of-donor error is
        # ~14 yr, and RES consumes sigma -- not q -- so it needs its own rescaling. The factor is
        # calibrated against one specific spread, so applying it under the other mode would
        # silently calibrate the wrong quantity. getattr keeps pre-Stage-1b bundles loading.
        self.sigma_scale = float(getattr(conf, "sigma_scale", 1.0))
        scale_mode = str(getattr(conf, "sigma_scale_mode", "ensemble"))
        if self.sigma_scale != 1.0 and scale_mode != self.mode:
            raise ConfigError(
                f"bundle's sigma_scale={self.sigma_scale:.3f} was calibrated for "
                f"mode={scale_mode!r} but this Predictor runs mode={self.mode!r}. "
                f"Re-run with mode={scale_mode!r}, or recompute the factor against "
                f"{self.mode!r} samples -- it calibrates a different spread in each mode."
            )

    # -- core stochastic passes -------------------------------------------- #
    @torch.no_grad()
    def _raw_batch(self, X: torch.Tensor, fp: torch.Tensor, dt: torch.Tensor):
        """Return (probs, age, latent): probs (K, N, 3), age (K, N), latent (N, d)."""
        if self.mode == "ensemble":
            probs, ages = [], []
            for m in self.members:
                lg, ag, _ = m(X, fp, dt)
                probs.append(torch.softmax(lg / self.temperature, dim=-1))
                ages.append(ag)
            _, _, z = self.members[0](X, fp, dt)
            return torch.stack(probs), torch.stack(ages), z

        # mc_dropout: latent from the deterministic (eval) pass, then T tiled passes
        m = self.members[0]
        _, _, z = m(X, fp, dt)
        enable_mc_dropout(m)
        try:
            n = X.shape[0]
            lg, ag, _ = m(X.repeat(self.T, 1), fp.repeat(self.T, 1), dt.repeat(self.T, 1))
            probs = torch.softmax(lg / self.temperature, dim=-1).view(self.T, n, 3)
            age = ag.view(self.T, n)
        finally:
            m.eval()
        return probs, age, z

    def _summaries(self, probs, age, z):
        pbar = probs.mean(0)                                    # (N, 3)
        ent = -(pbar * (pbar + 1e-12).log()).sum(1)
        Z = z.detach().cpu().numpy()
        return {
            "S": pbar[:, 0].cpu().numpy(),
            "P_loss": pbar[:, 1].cpu().numpy(),
            "P_death": pbar[:, 2].cpu().numpy(),
            "mu_age": age.mean(0).cpu().numpy(),
            "sigma_age": age.std(0, unbiased=False).cpu().numpy() * self.sigma_scale,
            "entropy": ent.cpu().numpy(),
            "in_dist": self.ood.in_distribution_mask(Z),
            "latent": Z,
        }

    @staticmethod
    def _rows(s, n):
        return [{
            "S": float(s["S"][i]), "P_loss": float(s["P_loss"][i]), "P_death": float(s["P_death"][i]),
            "mu_age": float(s["mu_age"][i]), "sigma_age": float(s["sigma_age"][i]),
            "entropy": float(s["entropy"][i]), "in_dist": bool(s["in_dist"][i]),
        } for i in range(n)]

    # -- public API -------------------------------------------------------- #
    def predict_batch(self, reqs) -> list[dict]:
        X, fp, dt = encode_batch(reqs, self.scalers)
        probs, age, z = self._raw_batch(X, fp, dt)
        return self._rows(self._summaries(probs, age, z), len(reqs))

    def predict(self, req) -> dict:
        return self.predict_batch([req])[0]

    def predict_encoded(self, X, fp, dose_time) -> list[dict]:
        """Batch summaries directly from Sample.X-space arrays (as stored in shards):
        X log-normalised expression, fp the 2048-bit fingerprint, dose_time the
        [log10 dose, log time] encoding. The scaler standardises X and dose_time."""
        Xz = torch.from_numpy(self.scalers.transform_x(np.asarray(X, np.float32))).float()
        dtz = torch.from_numpy(self.scalers.transform_dose_time(np.asarray(dose_time, np.float32))).float()
        fpt = torch.from_numpy(np.asarray(fp, np.float32)).float()
        probs, age, z = self._raw_batch(Xz, fpt, dtz)
        return self._rows(self._summaries(probs, age, z), Xz.shape[0])
