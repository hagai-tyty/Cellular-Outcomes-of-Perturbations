"""CellFate-Rx Gen-1 predictor — the Stage-23.5 §6 tool contract.

Scores one starting clone against the six observed WM989 conditions:

    f(X, B, U) -> future_detection_score,  for U in the six frozen conditions

The model is W5 (`X + B + U + X*U`), serialized by Stage 24C from the frozen Stage-23
implementation and verified to regenerate the frozen out-of-fold predictions to 5e-16.

Three contract points do real work here, and all three are refusals rather than conveniences:

  * `B` is required. Expression alone is NOT equivalent to the evaluated model, so a missing or
    incomplete nuisance vector returns `MISSING_REQUIRED_NUISANCE` and no score. It is never
    imputed to a plausible default (§6.2).
  * unknown conditions return `UNSUPPORTED_TREATMENT`. They are never embedded, nearest-neighboured
    or mapped onto a known condition (§6.3).
  * `validated_condition_order` is withheld unless Stage 25 has recorded `RANKING_SUPPORTED`.
    Until then the six scores are returned and their ORDER carries no validated meaning (§6.1).

The score is not a calibrated probability, not a measure of death or sensitivity, and not a
clinical recommendation. See `known_limitations` on every response.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["Gen1Predictor", "PredictionError", "clone_pseudobulk_from_counts",
           "clone_input_from_cells", "NAIVE_SAMPLES"]

# The three WM989 pretreatment libraries the frozen nuisance block counts over.
# These are a property of THIS experiment, not of the biology -- see the note in
# clone_input_from_cells about what that means for other datasets.
NAIVE_SAMPLES = ("Naive1", "Naive2", "Naive3")

CP10K = 10_000.0


class PredictionError(RuntimeError):
    """Raised only for contract violations that cannot be expressed as a support flag."""


def clone_pseudobulk_from_counts(counts: np.ndarray) -> np.ndarray:
    """§6.2 form A -> form B: sum raw pretreatment counts over a clone's cells, then CP10K + log1p.

    Applied EXACTLY ONCE, matching the frozen Stage-23 representation. Summing already-normalised
    cells, or log1p-ing a second time, produces a different feature space and therefore a model
    input the benchmark never evaluated.
    """
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim == 1:
        counts = counts[None, :]
    total = counts.sum()
    if total <= 0:
        raise PredictionError("clone has zero total counts; the frozen pipeline treats this as a "
                              "blocking condition rather than an imputation case")
    return np.log1p(counts.sum(axis=0) * (CP10K / total))


def clone_input_from_cells(counts: np.ndarray, samples: list[str]) -> tuple[np.ndarray,
                                                                                  np.ndarray]:
    """§6.2 form A -- the PREFERRED input. Build both `X` and `B` from raw pretreatment cells.

    `counts`   (n_cells, 36601) raw pretreatment Gene-Expression counts for ONE clone
    `samples`  the naive library each of those cells came from, e.g. ["Naive1", "Naive3", ...]

    Returns `(X, B)` ready for `Gen1Predictor.predict`. This is why the caller does not have to
    hand-compute the nuisance block: `B` is just cell counts, and the tool can count.

    IMPORTANT, and not solved by this function. `B` counts cells per WM989 naive library. Those
    three libraries are the structure of one experiment, not a property of melanoma. Supplying
    sample labels that merely LOOK like Naive1/2/3 -- from a different lab, a different depth, a
    different number of libraries -- produces a `B` the model never saw and a score the benchmark
    never evaluated. This function removes a chore for someone working with WM989-structured data;
    it does not make the model transferable.
    """
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim == 1:
        counts = counts[None, :]
    if len(samples) != counts.shape[0]:
        raise PredictionError(
            f"{counts.shape[0]} cells but {len(samples)} sample labels; every cell needs one")
    unknown = sorted(set(samples) - set(NAIVE_SAMPLES))
    if unknown:
        raise PredictionError(
            f"unknown naive libraries {unknown}; the frozen nuisance block is defined over "
            f"{list(NAIVE_SAMPLES)} and cannot be computed for anything else")

    x = clone_pseudobulk_from_counts(counts)
    per = [float(sum(1 for s in samples if s == lib)) for lib in NAIVE_SAMPLES]
    b = np.log1p(np.array([float(len(samples)), *per]))
    return x, b


@dataclass
class _Component:
    """One serialized W5 component: a gene filter, a PCA, three scalers and a linear model."""

    keep: np.ndarray
    gene_mu: np.ndarray
    gene_sd: np.ndarray
    pca_mean: np.ndarray
    pca_components: np.ndarray
    pc_mu: np.ndarray
    pc_sd: np.ndarray
    nuis_mu: np.ndarray
    nuis_sd: np.ndarray
    K: int
    C: float
    coef: np.ndarray
    intercept: float
    train_clones_sha256: str = ""

    def score(self, x: np.ndarray, b: np.ndarray, dummies: np.ndarray) -> np.ndarray:
        """Apply the frozen transform chain, then the linear model, for one clone x N conditions."""
        z = (x[self.keep] - self.gene_mu) / self.gene_sd
        pcs = ((z - self.pca_mean) @ self.pca_components.T - self.pc_mu) / self.pc_sd
        p = pcs[: self.K]
        bs = (b - self.nuis_mu) / self.nuis_sd
        n = dummies.shape[0]
        P = np.repeat(p[None, :], n, axis=0)
        B = np.repeat(bs[None, :], n, axis=0)
        inter = np.hstack([P * dummies[:, [t]] for t in range(dummies.shape[1])])
        A = np.hstack([P, B, dummies, inter])
        return 1.0 / (1.0 + np.exp(-(A @ self.coef + self.intercept)))


@dataclass
class Gen1Predictor:
    """Loaded from the Stage-24C artifact. See `Gen1Predictor.load`."""

    meta: dict
    components: dict[str, _Component]
    ranking_status: str = "NOT_SUPPORTED"
    _limitations: list[str] = field(default_factory=list)

    # ---- construction ------------------------------------------------------------------------ #
    @classmethod
    def load(cls, artifact_npz: Path | str, artifact_meta: Path | str,
             stage25_verdict: Path | str | None = None) -> Gen1Predictor:
        meta = json.loads(Path(artifact_meta).read_text(encoding="utf-8"))
        z = np.load(artifact_npz, allow_pickle=False)

        def build(prefix: str) -> _Component:
            g = {k.split("__", 1)[1]: z[k] for k in z.files if k.startswith(prefix + "__")}
            return _Component(
                keep=g["keep"].astype(int), gene_mu=g["gene_mu"], gene_sd=g["gene_sd"],
                pca_mean=g["pca_mean"], pca_components=g["pca_components"],
                pc_mu=g["pc_mu"], pc_sd=g["pc_sd"],
                nuis_mu=g["nuis_mu"], nuis_sd=g["nuis_sd"],
                K=int(g["K"]), C=float(g["C"]),
                coef=g["coef"], intercept=float(g["intercept"]),
                train_clones_sha256=str(g["train_clones_sha256"])
                if "train_clones_sha256" in g else "")

        prefixes = sorted({k.split("__", 1)[0] for k in z.files})
        components = {p: build(p) for p in prefixes}

        # Ranking stays unvalidated unless Stage 25 has actually said otherwise.
        status = "NOT_SUPPORTED"
        if stage25_verdict is not None and Path(stage25_verdict).exists():
            v = json.loads(Path(stage25_verdict).read_text(encoding="utf-8"))
            if v.get("verdict") == "STAGE_25_RANKING_SUPPORTED":
                status = "SUPPORTED"
        return cls(meta=meta, components=components, ranking_status=status,
                   _limitations=list(meta.get("known_limitations", [])))

    # ---- contract helpers -------------------------------------------------------------------- #
    @property
    def conditions(self) -> list[str]:
        return list(self.meta["treatment_vocabulary"])

    @property
    def reference_condition(self) -> str:
        return self.meta["reference_treatment"]

    def _dummies(self, treatments: list[str]) -> np.ndarray:
        non_ref = [t for t in self.conditions if t != self.reference_condition]
        arr = np.asarray(treatments)
        return np.column_stack([(arr == t).astype(float) for t in non_ref])

    # ---- the API ----------------------------------------------------------------------------- #
    def predict(self, expression: np.ndarray, nuisance: np.ndarray | None = None, *,
                treatments: list[str] | None = None,
                component: str = "deployment") -> list[dict]:
        """Score one clone against the requested conditions.

        `expression`  clone-level CP10K/log1p vector of length `n_expression_features_expected`
        `nuisance`    the complete frozen nuisance block, in `meta["nuisance_columns"]` order
        `treatments`  defaults to all six supported conditions
        `component`   "deployment" for a new clone; "fold{0..4}" reproduces a benchmark clone's
                      out-of-fold score using the one model that did not train on it
        """
        base = {"model_version": self.meta["model_version"],
                "feature_contract_version": self.meta["feature_contract_version"],
                "ranking_status": self.ranking_status,
                "known_limitations": self._limitations}
        requested = list(treatments) if treatments is not None else self.conditions

        if component not in self.components:
            raise PredictionError(
                f"unknown component {component!r}; available: {sorted(self.components)}")

        x = np.asarray(expression, dtype=np.float64).ravel()
        expected = int(self.meta["n_expression_features_expected"])
        if x.shape[0] != expected:
            return [{**base, "condition": t, "future_detection_score": None,
                     "support_status": "UNSUPPORTED_FEATURE_SCHEMA",
                     "detail": f"expected {expected} expression features, got {x.shape[0]}"}
                    for t in requested]

        n_nuis = len(self.meta["nuisance_columns"])
        b = None if nuisance is None else np.asarray(nuisance, dtype=np.float64).ravel()
        if b is None or b.shape[0] != n_nuis or not np.isfinite(b).all():
            # §6.2: fail closed. B is part of the evaluated model and is never imputed.
            return [{**base, "condition": t, "future_detection_score": None,
                     "support_status": "MISSING_REQUIRED_NUISANCE",
                     "detail": f"the complete nuisance block is required "
                               f"({self.meta['nuisance_columns']}); it is never imputed"}
                    for t in requested]

        known = [t for t in requested if t in self.conditions]
        out: list[dict] = []
        scores: dict[str, float] = {}
        if known:
            s = self.components[component].score(x, b, self._dummies(known))
            scores = dict(zip(known, (float(v) for v in s), strict=True))

        for t in requested:
            if t not in self.conditions:
                out.append({**base, "condition": t, "future_detection_score": None,
                            "support_status": "UNSUPPORTED_TREATMENT",
                            "detail": "not one of the six observed conditions; unknown conditions "
                                      "are never embedded or mapped to a known condition"})
            else:
                out.append({**base, "condition": t,
                            "future_detection_score": scores[t],
                            "support_status": "SUPPORTED_KNOWN_CONDITION",
                            "component": component})
        return out

    def rank_conditions(self, expression: np.ndarray, nuisance: np.ndarray | None = None, *,
                        component: str = "deployment") -> dict:
        """Return the six scores, and an ORDER only if Stage 25 validated ordering.

        Withholding the order is the point. Until Stage 25 records RANKING_SUPPORTED, sorting the
        six scores is an unvalidated operation and the tool must not present it as condition
        selection (§6.1, §6.4).
        """
        rows = self.predict(expression, nuisance, component=component)
        scored = [r for r in rows if r["future_detection_score"] is not None]
        res = {"scores": {r["condition"]: r["future_detection_score"] for r in scored},
               "ranking_status": self.ranking_status,
               "model_version": self.meta["model_version"],
               "known_limitations": self._limitations}
        if not scored:
            res["support_status"] = rows[0]["support_status"] if rows else "OUT_OF_CONTRACT_INPUT"
            return res
        if self.ranking_status == "SUPPORTED":
            # lower predicted detection = the condition the model expects to leave the clone
            # undetected. Experimental-condition selection, NOT clinical recommendation.
            res["validated_condition_order"] = [
                c for c, _ in sorted(res["scores"].items(), key=lambda kv: (kv[1], kv[0]))]
        else:
            res["validated_condition_order"] = None
            res["detail"] = ("Stage 25 has not recorded RANKING_SUPPORTED, so the order of these "
                             "scores is not a validated condition ranking.")
        return res
