"""Runtime normalisation. ``Scalers`` wraps the serialisable :class:`ScalerParams`
and provides the actual transforms used identically at training and inference
time. Fit on the TRAIN split only (Document 2, S12).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .io import read_json, write_json
from .panel import GenePanel
from .schemas import ScalerParams

_EPS = 1e-6


class Scalers:
    """Per-gene and dose/time standardisation, plus the proliferation deconfounder.

    The deconfounder coefficients (a, b) are carried for provenance/inverse use;
    deconfounding itself is applied at *label* time in the data package.
    """

    def __init__(self, params: ScalerParams) -> None:
        self.params = params
        self._x_mean = np.asarray(params.x_mean, dtype=np.float32)
        self._x_std = np.asarray(params.x_std, dtype=np.float32)
        self._dt_mean = np.asarray(params.dt_mean, dtype=np.float32)
        self._dt_std = np.asarray(params.dt_std, dtype=np.float32)

    # -- transforms --------------------------------------------------------- #
    def transform_x(self, x: np.ndarray) -> np.ndarray:
        """z-score expression. Accepts (G,) or (N, G)."""
        x = np.asarray(x, dtype=np.float32)
        if x.shape[-1] != self._x_mean.shape[0]:
            raise ValueError(
                f"X has {x.shape[-1]} genes but scalers expect {self._x_mean.shape[0]}"
            )
        return (x - self._x_mean) / (self._x_std + _EPS)

    def transform_dose_time(self, dt: np.ndarray) -> np.ndarray:
        """z-score the (already log-transformed) dose/time vector. (2,) or (N, 2)."""
        dt = np.asarray(dt, dtype=np.float32)
        if dt.shape[-1] != self._dt_mean.shape[0]:
            raise ValueError(f"dose_time must have {self._dt_mean.shape[0]} entries")
        return (dt - self._dt_mean) / (self._dt_std + _EPS)

    @property
    def proliferation_coef(self) -> tuple[float, float]:
        a, b = self.params.proliferation_coef
        return float(a), float(b)

    # -- fitting (train split only) ----------------------------------------- #
    @classmethod
    def fit(
        cls,
        x_train: np.ndarray,
        dose_time_train: np.ndarray,
        gene_panel: GenePanel,
        proliferation_coef: tuple[float, float] = (0.0, 0.0),
    ) -> Scalers:
        x_train = np.asarray(x_train, dtype=np.float64)
        dt_train = np.asarray(dose_time_train, dtype=np.float64)
        params = ScalerParams(
            x_mean=x_train.mean(0).astype(np.float32).tolist(),
            x_std=x_train.std(0).astype(np.float32).tolist(),
            dt_mean=dt_train.mean(0).astype(np.float32).tolist(),
            dt_std=dt_train.std(0).astype(np.float32).tolist(),
            proliferation_coef=[float(proliferation_coef[0]), float(proliferation_coef[1])],
            gene_panel_hash=gene_panel.hash(),
        )
        return cls(params)

    # -- persistence -------------------------------------------------------- #
    def save(self, path: str | Path) -> None:
        write_json(path, self.params.model_dump())

    @classmethod
    def load(cls, path: str | Path) -> Scalers:
        return cls(ScalerParams.model_validate(read_json(path)))
