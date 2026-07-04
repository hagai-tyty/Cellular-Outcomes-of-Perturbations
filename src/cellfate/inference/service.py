"""Service layer (Document 4, S6): one Predictor backs an offline batch scorer /
CLI and (optionally) an HTTP endpoint. FastAPI is imported lazily so the package
is fully usable -- and testable -- without the web dependency installed.

An out-of-distribution query is *not* an error: it returns a valid Response with
``status=REJECTED_OOD`` and a warning. Malformed requests raise typed CellFate
errors (surfaced as HTTP 422 when served).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cellfate.common import io
from cellfate.common.errors import ContractViolation

from .conformal import interval
from .predictor import Predictor
from .res import compute_res
from .schema import Request, Response


def build_response(pred: Predictor, s: dict) -> Response:
    """Assemble a Response from a Predictor summary (RES + status + interval)."""
    res, status = compute_res(s["S"], s["P_loss"], s["mu_age"], s["sigma_age"],
                              s["in_dist"], pred.res_params)
    lo, hi = interval(s["mu_age"], pred.q)
    return Response(
        status=status,
        rejuvenation_efficacy_score=round(10.0 * res, 2),
        p_identity_preserved=round(s["S"], 3),
        p_identity_loss=round(s["P_loss"], 3),
        p_apoptosis=round(s["P_death"], 3),
        delta_age_mean=round(s["mu_age"], 2),
        delta_age_interval=[round(lo, 2), round(hi, 2)],
        in_distribution=s["in_dist"],
        epistemic_std=round(s["sigma_age"], 3),
        predictive_entropy=round(s["entropy"], 3),
        warning=None if s["in_dist"] else "Out-of-distribution: prediction not trustworthy.",
    )


def predict_one(pred: Predictor, req: Request) -> Response:
    return build_response(pred, pred.predict(req))


def score_requests(pred: Predictor, reqs: list[Request]) -> list[Response]:
    return [build_response(pred, s) for s in pred.predict_batch(reqs)]


def score_shard(pred: Predictor, shard_path):
    """Rank every cell in a data shard (chem or TF). Returns (responses, cell_ids)."""
    arr = io.shard_to_numpy(io.read_shard(shard_path))
    # perturbation input: chem fingerprint OR TF-cocktail multi-hot
    pert = arr["u_chem_fp"] if arr["u_chem_fp"] is not None else arr["u_tf_emb"]
    if pert is None:
        raise ContractViolation("shard has no perturbation features (neither u_chem_fp nor u_tf_emb)")
    summaries = pred.predict_encoded(arr["X"], pert, arr["dose_time"])
    return [build_response(pred, s) for s in summaries], list(arr["cell_id"])


def create_app(bundle_dir, mode: str = "ensemble", T: int = 50):
    """Build a FastAPI app exposing POST /predict. Requires `fastapi` (+ a server)."""
    try:
        from fastapi import FastAPI, HTTPException
    except ImportError as e:  # pragma: no cover - optional dependency
        raise RuntimeError("FastAPI is not installed; `pip install fastapi uvicorn` to serve HTTP") from e
    from cellfate.common.errors import CellFateError

    app = FastAPI(title="CellFate-Rx Inference", version="1.0")
    pred = Predictor(bundle_dir, mode=mode, T=T)

    @app.post("/predict", response_model=Response)
    def predict(req: Request):  # pragma: no cover - exercised via TestClient if fastapi present
        try:
            return predict_one(pred, req)
        except CellFateError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    return app


# --------------------------------------------------------------------------- #
# CLI: rank candidate perturbations in a shard by rejuvenation efficacy         #
# --------------------------------------------------------------------------- #
def _rank_shard(pred: Predictor, shard_path: Path, top: int) -> list[dict]:
    responses, cell_ids = score_shard(pred, shard_path)
    rows = [{"cell_id": cid, **resp.model_dump()} for cid, resp in zip(cell_ids, responses, strict=True)]
    # APPROVED first, then by score descending
    rows.sort(key=lambda r: (r["status"] != "APPROVED", -r["rejuvenation_efficacy_score"]))
    return rows[:top] if top else rows


def cli() -> None:
    ap = argparse.ArgumentParser(
        prog="cellfate-serve",
        description="Score/rank candidate perturbations from a CellFate-Rx bundle.",
    )
    ap.add_argument("--bundle", required=True, help="artefact root (contains bundle/)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--shard", help="a single .parquet shard to score")
    src.add_argument("--dataset", help="a dataset root; scores every shard under it")
    ap.add_argument("--mode", default="ensemble", choices=["ensemble", "mc_dropout"])
    ap.add_argument("--mc-samples", type=int, default=50, dest="T")
    ap.add_argument("--top", type=int, default=20, help="keep the top-N (0 = all)")
    ap.add_argument("--out", help="write ranked results as JSON to this path")
    args = ap.parse_args()

    pred = Predictor(args.bundle, mode=args.mode, T=args.T)
    if args.shard:
        shards = [Path(args.shard)]
    else:
        shards = sorted(io.ArtifactPaths.of(args.dataset).shards_dir.glob("*.parquet"))
        if not shards:
            raise SystemExit(f"no shards found under {args.dataset}")

    rows: list[dict] = []
    for sh in shards:
        rows.extend(_rank_shard(pred, sh, top=0))
    rows.sort(key=lambda r: (r["status"] != "APPROVED", -r["rejuvenation_efficacy_score"]))
    if args.top:
        rows = rows[: args.top]

    n_appr = sum(r["status"] == "APPROVED" for r in rows)
    print(f"scored {len(shards)} shard(s); {n_appr}/{len(rows)} shown are APPROVED "
          f"(mode={args.mode}, conformal={pred.conformal_level})")
    for r in rows[:20]:
        print(f"  {r['status']:<26} RES={r['rejuvenation_efficacy_score']:<5} "
              f"P_safe={r['p_identity_preserved']:<5} ΔAge={r['delta_age_mean']:<7} "
              f"{r['delta_age_interval']}  {r['cell_id']}")
    if args.out:
        Path(args.out).write_text(json.dumps(rows, indent=2))
        print(f"wrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    cli()
