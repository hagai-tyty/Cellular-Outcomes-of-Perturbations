"""Transcriptomic aging clock + ΔAge (Document 2, S10).

The clock predicts a transcriptomic age from expression; ΔAge is the predicted
age relative to the matched vehicle-control baseline of the same cell line
(negative = rejuvenated). On cancer / transformed lines the clock is
out-of-distribution, so ΔAge is masked (``age_mask = False``) -- the safety head
still trains on those cells, the age head does not.

``LinearClock`` (age = w.x + b over panel genes) is the dependency-free default
and the interface real clocks (Buckley et al.; scAgeClock) plug into.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd

from cellfate.common import constants as C
from cellfate.common.panel import GenePanel


class AgingClock(ABC):
    """Predicts a transcriptomic age (years) from expression.

    The clock consumes the **full normalised profile** and its own gene panel --
    it aligns its weights against the gene symbols it is handed, NOT the model's
    2000-HVG input. So it is decoupled from the model input: the model sees the
    HVG panel; the clock sees every gene it can match.
    """

    @abstractmethod
    def predict_age(self, expr: np.ndarray, genes: list[str]) -> np.ndarray:
        """Return (N,) predicted ages for an (N, len(genes)) matrix in ``genes`` order."""


class LinearClock(AgingClock):
    """age = sum_g w_g * x_g + b, with weights keyed by gene symbol."""

    def __init__(self, weights: dict[str, float], intercept: float = 0.0) -> None:
        self.weights = weights
        self.intercept = float(intercept)

    def predict_age(self, expr: np.ndarray, genes: list[str]) -> np.ndarray:
        w = np.array([self.weights.get(g, 0.0) for g in genes], dtype=np.float64)
        return np.asarray(expr, dtype=np.float64) @ w + self.intercept

    @classmethod
    def random(cls, panel: GenePanel, seed: int = 0, scale: float = 1.0) -> LinearClock:
        """A deterministic **random** clock over a panel.

        For synthetic/smoke runs and tests ONLY -- its ages are meaningless. Real
        runs must load fitted weights via :meth:`from_json` (see scripts/fit_clock.py).
        """
        rng = np.random.default_rng(seed)
        w = rng.normal(0.0, scale, size=len(panel)) / np.sqrt(len(panel))
        return cls({g: float(wi) for g, wi in zip(panel.genes, w, strict=True)}, intercept=40.0)

    @classmethod
    def from_json(cls, path: str | Path) -> LinearClock:
        """Load a fitted clock: ``{"weights": {gene: w}, "intercept": b, ...}``."""
        d = json.loads(Path(path).read_text())
        if "weights" not in d:
            raise ValueError(f"clock file {path} has no 'weights' key")
        weights = {str(k): float(v) for k, v in d["weights"].items()}
        if not weights:
            raise ValueError(f"clock file {path} has empty weights")
        return cls(weights, intercept=float(d.get("intercept", 0.0)))

    def to_json(self, path: str | Path, meta: dict | None = None) -> None:
        """Serialise fitted weights (+ optional provenance ``meta``) to JSON."""
        payload = {"weights": self.weights, "intercept": self.intercept}
        if meta:
            payload["meta"] = meta
        Path(path).write_text(json.dumps(payload, indent=2))


def _control_baseline(values: np.ndarray, lines: np.ndarray, is_ctrl: np.ndarray) -> np.ndarray:
    """Per-line mean over vehicle controls. Falls back to the line's own mean when
    a line has no controls in this chunk (values then centred within the line)."""
    baseline = np.empty_like(values, dtype=np.float64)
    for line in np.unique(lines):
        in_line = lines == line
        ctrl = in_line & is_ctrl
        ref = values[ctrl] if ctrl.any() else values[in_line]
        baseline[in_line] = ref.mean()
    return baseline


def recenter_on_control_arrays(
    values: np.ndarray, lines: np.ndarray, is_ctrl: np.ndarray
) -> np.ndarray:
    """Array form of :func:`recenter_on_controls` (no DataFrame needed).

    Subtracts the per-line vehicle-control baseline from ``values`` given the
    per-cell ``lines`` and boolean ``is_ctrl`` arrays.
    """
    return np.asarray(values, dtype=np.float64) - _control_baseline(values, lines, is_ctrl)


def recenter_on_controls(values: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    """Subtract the per-line vehicle-control baseline from ``values``.

    Used to re-anchor ΔAge after cell-cycle deconfounding: ``deconfound_age``
    removes the regression intercept and so re-centres the whole population to
    mean 0, which shifts the controls off zero. Re-applying the control baseline
    restores the invariant that ΔAge is *control-relative* (controls ~ 0), which
    is the zero-point the rejuvenation score depends on -- without reintroducing
    the cell-cycle slope that was just removed.
    """
    lines = obs["cell_line"].to_numpy()
    is_ctrl = obs["is_control"].to_numpy().astype(bool)
    return recenter_on_control_arrays(values, lines, is_ctrl)


def delta_age(
    clock: AgingClock,
    expr: np.ndarray,
    genes: list[str],
    obs: pd.DataFrame,
    source: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute (ΔAge, age_mask) from the **full normalised profile**.

    The clock is handed ``expr`` (N, len(genes)) with its gene symbols ``genes``
    -- the full profile, not the 2000-HVG model input -- so it can dot its
    weights against every gene it matches.

    ΔAge[i] = age[i] - mean(age over vehicle controls of the same cell line).
    If a cell line has no controls in this chunk, its own mean age is used as the
    baseline (ΔAge centred within the line). ``age_mask`` is all-False for
    sources in ``CANCER_SOURCES``.

    Note: if the caller subsequently deconfounds ΔAge for cell cycle, it must
    re-anchor with ``recenter_on_controls`` to preserve the control-relative
    zero-point (deconfounding otherwise re-centres the whole population).
    """
    age = clock.predict_age(expr, genes)
    lines = obs["cell_line"].to_numpy()
    is_ctrl = obs["is_control"].to_numpy().astype(bool)
    d = age - _control_baseline(age, lines, is_ctrl)
    age_mask = np.full(age.shape[0], source not in C.CANCER_SOURCES, dtype=bool)
    return d, age_mask
