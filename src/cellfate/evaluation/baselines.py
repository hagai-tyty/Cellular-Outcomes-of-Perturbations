"""Mandatory baselines + the model wrapper, behind one ``Estimator`` interface
(Document 5, S2).

Nature Methods benchmarks and the Arc Virtual Cell Challenge repeatedly found
sophisticated models failing to beat trivial baselines, so "beats every baseline"
is an automated acceptance gate. Diagnostic reads: if X-only ties the model the
perturbation is ignored; if U-only ties it the cell state is ignored; if mean or
predict-control ties it the model learned nothing beyond central tendency.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from cellfate.data.splits import CONTROL_SCAFFOLD

from .data import SplitData


@runtime_checkable
class Estimator(Protocol):
    def fit(self, train: SplitData) -> Estimator: ...
    def predict(self, X, fp, dose_time): ...  # -> (probs (N,3), age (N,))


def _full_proba(clf: LogisticRegression, Xf: np.ndarray) -> np.ndarray:
    """LogisticRegression.predict_proba widened to a full (N, 3) matrix, filling
    any class absent from the training split with probability 0."""
    p = clf.predict_proba(Xf)
    full = np.zeros((Xf.shape[0], 3), dtype=np.float64)
    for j, c in enumerate(clf.classes_.astype(int)):
        full[:, c] = p[:, j]
    return full


class MeanBaseline:
    """Train marginal class frequency + mean ΔAge; ignores all inputs."""

    def fit(self, tr: SplitData) -> MeanBaseline:
        self._freq = (np.bincount(tr.y_cls.astype(int), minlength=3) / max(tr.n, 1)).astype(np.float64)
        self._age = float(tr.y_age[tr.mask].mean()) if tr.mask.any() else 0.0
        return self

    def predict(self, X, fp, dose_time):
        n = len(X)
        return np.tile(self._freq, (n, 1)), np.full(n, self._age)


class PredictControl:
    """'No change': ΔAge = 0, class = matched vehicle-control distribution."""

    def fit(self, tr: SplitData) -> PredictControl:
        ctrl = tr.scaffold_id == CONTROL_SCAFFOLD
        y = tr.y_cls[ctrl] if ctrl.any() else tr.y_cls
        self._freq = (np.bincount(y.astype(int), minlength=3) / max(len(y), 1)).astype(np.float64)
        return self

    def predict(self, X, fp, dose_time):
        n = len(X)
        return np.tile(self._freq, (n, 1)), np.zeros(n)


class _LinearBase:
    """Logistic (class) + Ridge (age) on a configurable feature view."""

    view = "all"  # "all" | "x" | "u"

    def fit(self, tr: SplitData) -> _LinearBase:
        self._sx = StandardScaler().fit(tr.X) if tr.n else None
        self._sdt = StandardScaler().fit(tr.dose_time) if tr.n else None
        Xf = self._features(tr.X, tr.fp, tr.dose_time)
        self._clf = LogisticRegression(max_iter=1000, C=1.0).fit(Xf, tr.y_cls.astype(int))
        if tr.mask.sum() >= 2:
            self._reg = Ridge(alpha=1.0).fit(Xf[tr.mask], tr.y_age[tr.mask])
            self._agebar = float(tr.y_age[tr.mask].mean())
        else:
            self._reg, self._agebar = None, 0.0
        return self

    def _features(self, X, fp, dose_time) -> np.ndarray:
        parts = []
        if self.view in ("all", "x"):
            parts.append(self._sx.transform(X))
        if self.view in ("all", "u"):
            parts.append(np.asarray(fp, dtype=np.float64))
            parts.append(self._sdt.transform(dose_time))
        return np.hstack(parts)

    def predict(self, X, fp, dose_time):
        Xf = self._features(X, fp, dose_time)
        probs = _full_proba(self._clf, Xf)
        age = self._reg.predict(Xf) if self._reg is not None else np.full(len(X), self._agebar)
        return probs, np.asarray(age, dtype=np.float64)


class RidgeLinear(_LinearBase):
    view = "all"


class XOnly(_LinearBase):
    view = "x"


class UOnly(_LinearBase):
    view = "u"


class KNNFingerprint:
    """k-NN over Morgan-fingerprint (Jaccard) space; average neighbours' outcomes."""

    def __init__(self, k: int = 15):
        self.k = k

    def fit(self, tr: SplitData) -> KNNFingerprint:
        k = max(1, min(self.k, tr.n))
        self._nn = NearestNeighbors(n_neighbors=k, metric="jaccard").fit(tr.fp.astype(bool))
        self._ycls = tr.y_cls.astype(int)
        self._yage = tr.y_age
        self._mask = tr.mask
        return self

    def predict(self, X, fp, dose_time):
        _, idx = self._nn.kneighbors(np.asarray(fp).astype(bool))
        probs = np.zeros((len(fp), 3), dtype=np.float64)
        age = np.zeros(len(fp), dtype=np.float64)
        for i, nb in enumerate(idx):
            probs[i] = np.bincount(self._ycls[nb], minlength=3) / len(nb)
            valid = self._mask[nb]
            age[i] = float(self._yage[nb][valid].mean()) if valid.any() else 0.0
        return probs, age


class ModelEstimator:
    """Wraps a trained ``Predictor`` behind the Estimator interface."""

    name = "model"

    def __init__(self, predictor):
        self.pred = predictor

    def fit(self, tr: SplitData) -> ModelEstimator:
        return self  # already trained

    def rows(self, X, fp, dose_time) -> list[dict]:
        return self.pred.predict_encoded(
            np.asarray(X, np.float32), np.asarray(fp, np.float32), np.asarray(dose_time, np.float32)
        )

    def predict(self, X, fp, dose_time):
        rows = self.rows(X, fp, dose_time)
        probs = np.array([[r["S"], r["P_loss"], r["P_death"]] for r in rows], dtype=np.float64)
        age = np.array([r["mu_age"] for r in rows], dtype=np.float64)
        return probs, age


_REGISTRY = {
    "mean": MeanBaseline,
    "ridge": RidgeLinear,
    "x_only": XOnly,
    "u_only": UOnly,
    "knn": KNNFingerprint,
    "predict_control": PredictControl,
}
BASELINE_NAMES = tuple(_REGISTRY)


def make_baselines(names) -> dict[str, Estimator]:
    return {n: _REGISTRY[n]() for n in names}
