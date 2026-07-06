"""Orchestrator + acceptance gates (Document 5, S5).

Runs the model and every baseline on each held-out regime, computes all metrics,
writes reports/<regime>.{json,md} + a summary, and turns the Document-1 success
criteria into pass/fail gates. The test split is touched only here, once per regime;
thresholds come from val and the conformal quantile from calib -- that separation is
what makes the gates trustworthy.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field

import numpy as np

from cellfate.common.constants import Split
from cellfate.common.io import ArtifactPaths
from cellfate.inference import Predictor, compute_res_batch

from .baselines import BASELINE_NAMES, ModelEstimator, make_baselines
from .data import gather_split
from .external_validation import validate_against_methylation, validate_oskm_holdout
from .metrics import (
    brier,
    coverage,
    ece,
    mean_finite,
    per_class_auroc,
    per_class_prauc,
    ranking_metrics,
    regression_metrics,
)
from .regimes import REGIMES, iter_regimes
from .report import write_report, write_summary

NAN = float("nan")


@dataclass
class EvalConfig:
    bundle: str
    dataset: str
    regimes: tuple = tuple(REGIMES)
    baselines: tuple = BASELINE_NAMES
    out: str = "reports"
    max_ece: float = 0.05
    level: float = 0.90
    cov_tol: float = 0.03
    min_spearman: float = 0.3
    external: dict | None = field(default=None)


def _estimator_metrics(test, p, age) -> dict:
    m: dict[str, float] = {}
    for c, v in per_class_prauc(test.y_cls, p).items():
        m[f"prauc_{c}"] = v
    for c, v in per_class_auroc(test.y_cls, p).items():
        m[f"auroc_{c}"] = v
    m["brier"] = brier(test.y1h, p)
    m["ece"] = ece(test.y_cls, p)
    for k, v in regression_metrics(test.y_age, age, test.mask).items():
        m[f"reg_{k}"] = v
    return m


def _mean_prauc(metrics: dict) -> float:
    return mean_finite(metrics.get(f"prauc_{c}") for c in range(3))


def run_external(model_est: ModelEstimator, cfg: EvalConfig) -> dict:
    ext = cfg.external or {}
    if not ext:
        return {"status": "not_available",
                "note": "no methylation/OSKM data configured; external anchors require real data"}
    out: dict = {"status": "evaluated"}
    if "methylation" in ext:
        out["methylation"] = validate_against_methylation(
            ext["methylation"]["pred"], ext["methylation"]["meth"])
    if "oskm" in ext:
        out["oskm"] = validate_oskm_holdout(model_est.predict, ext["oskm"])
    return out


def evaluate(cfg: EvalConfig) -> dict:
    pred = Predictor(cfg.bundle)
    paths = ArtifactPaths.of(cfg.dataset)
    model_est = ModelEstimator(pred)
    baselines = make_baselines(cfg.baselines)
    results: dict = {}

    for regime in iter_regimes(cfg.regimes):
        train = gather_split(paths, regime, Split.TRAIN.value)
        test = gather_split(paths, regime, Split.TEST.value)
        if test.n == 0:
            results[regime] = {"_empty": True}
            continue

        R: dict = {}
        rows = model_est.rows(test.X, test.fp, test.dose_time)
        p_model = np.array([[r["S"], r["P_loss"], r["P_death"]] for r in rows], dtype=np.float64)
        age_model = np.array([r["mu_age"] for r in rows], dtype=np.float64)
        R["model"] = _estimator_metrics(test, p_model, age_model)
        R["_y_true"] = test.y_cls.tolist()
        R["_p_model"] = p_model.tolist()

        for name, est in baselines.items():
            est.fit(train)
            p, age = est.predict(test.X, test.fp, test.dose_time)
            R[name] = _estimator_metrics(test, p, age)

        lo, hi = age_model - pred.q, age_model + pred.q
        R["coverage"] = coverage(test.y_age, lo, hi, test.mask)
        res, _ = compute_res_batch(
            [r["S"] for r in rows], [r["P_loss"] for r in rows], age_model,
            [r["sigma_age"] for r in rows], [r["in_dist"] for r in rows], pred.res_params)
        R["ranking"] = (ranking_metrics(res[test.mask], test.y_age[test.mask])
                        if test.mask.any() else {"spearman": NAN, "precision_at_k": NAN})

        results[regime] = R
        write_report(regime, R, cfg.out)

    results["external"] = run_external(model_est, cfg)
    gates = check_gates(results, cfg)
    write_summary(results, gates, cfg.out)
    return gates


def check_gates(results: dict, cfg: EvalConfig) -> dict:
    gates: dict = {}
    for regime, R in results.items():
        if regime == "external" or "model" not in R:
            continue
        model_prauc = _mean_prauc(R["model"])
        model_mae = R["model"].get("reg_mae", NAN)

        beats = bool(np.isfinite(model_prauc))
        for b in cfg.baselines:
            bp = _mean_prauc(R[b])
            bmae = R[b].get("reg_mae", NAN)
            if np.isfinite(bp) and not (model_prauc > bp):
                beats = False
            if np.isfinite(bmae) and np.isfinite(model_mae) and not (model_mae <= bmae):
                beats = False

        cov = R.get("coverage", NAN)
        rank = R.get("ranking", {}).get("spearman", NAN)
        model_ece = R["model"].get("ece", NAN)
        gates[regime] = {
            "beats_all_baselines": beats,                                             # crit. 4
            "ece_ok": bool(np.isfinite(model_ece) and model_ece < cfg.max_ece),       # crit. 1
            "coverage_ok": bool(np.isfinite(cov) and abs(cov - cfg.level) < cfg.cov_tol),  # crit. 6
            "ranking_ok": bool(np.isfinite(rank) and rank > cfg.min_spearman),        # crit. 5
        }
    return gates


def cli() -> None:
    ap = argparse.ArgumentParser(
        prog="cellfate-evaluate",
        description="Run the model + baselines across regimes and check acceptance gates.",
    )
    ap.add_argument("--bundle", required=True, help="artefact root (contains bundle/)")
    ap.add_argument("--dataset", required=True, help="dataset root (shards/ + splits/)")
    ap.add_argument("--regimes", nargs="+", default=list(REGIMES))
    ap.add_argument("--baselines", nargs="+", default=list(BASELINE_NAMES))
    ap.add_argument("--out", default="reports")
    ap.add_argument("--max-ece", type=float, default=0.05)
    ap.add_argument("--level", type=float, default=0.90)
    ap.add_argument("--cov-tol", type=float, default=0.03)
    ap.add_argument("--min-spearman", type=float, default=0.3)
    args = ap.parse_args()

    cfg = EvalConfig(
        bundle=args.bundle, dataset=args.dataset, regimes=tuple(args.regimes),
        baselines=tuple(args.baselines), out=args.out, max_ece=args.max_ece,
        level=args.level, cov_tol=args.cov_tol, min_spearman=args.min_spearman,
    )
    gates = evaluate(cfg)
    for regime, g in gates.items():
        status = "PASS" if all(g.values()) else "FAIL"
        detail = ", ".join(f"{k}={'ok' if v else 'no'}" for k, v in g.items())
        print(f"[{status}] {regime}: {detail}")
    print(f"reports written to {args.out}/")


if __name__ == "__main__":
    cli()
