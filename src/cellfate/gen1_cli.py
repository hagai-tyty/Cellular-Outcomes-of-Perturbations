"""Command-line interface for the CellFate-Rx Gen-1 predictor (Stage-23.5 §6.5).

    python -m cellfate.gen1_cli --artifact <npz> --meta <json> \\
        --expression clone.npy --nuisance 12.0,4.0,3.0,2.0

Prints one JSON object per requested condition, carrying its score, its support flag, its
provenance and its limitations. Refusals are printed the same way successes are — a caller that
pipes this into a table gets the `support_status` column whether or not a score came back, and a
missing score is never silently a zero.

Exit codes:
    0  every requested condition was scored
    2  at least one condition was refused (unsupported treatment, missing nuisance, bad schema)
    3  the input could not be read at all
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from cellfate.gen1_predictor import Gen1Predictor, PredictionError


def _load_vector(spec: str) -> np.ndarray:
    """Accept a .npy path, a .csv/.txt path, or a bare comma-separated list."""
    p = Path(spec)
    if p.exists():
        if p.suffix == ".npy":
            return np.load(p).ravel()
        return np.loadtxt(p, delimiter=",").ravel()
    return np.array([float(x) for x in spec.split(",")])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="cellfate.gen1_cli",
        description="Score one clone against the six observed WM989 conditions.")
    ap.add_argument("--artifact", required=True, help="stage24_w5_artifact.npz")
    ap.add_argument("--meta", required=True, help="stage24_w5_artifact.json")
    ap.add_argument("--expression", required=True,
                    help=".npy/.csv path, or comma-separated clone-level CP10K/log1p values")
    ap.add_argument("--nuisance", default=None,
                    help="the complete frozen nuisance block; REQUIRED for a score")
    ap.add_argument("--treatments", default=None,
                    help="comma-separated; defaults to all six supported conditions")
    ap.add_argument("--component", default="deployment",
                    help="'deployment' for a new clone, or fold0..fold4 to reproduce a "
                         "benchmark clone's out-of-fold score")
    ap.add_argument("--stage25-verdict", default=None,
                    help="Stage-25 verdict JSON; without it, ranking stays NOT_SUPPORTED")
    ap.add_argument("--rank", action="store_true",
                    help="return the six scores plus a condition order IF Stage 25 validated it")
    a = ap.parse_args(argv)

    try:
        x = _load_vector(a.expression)
        b = _load_vector(a.nuisance) if a.nuisance else None
    except (OSError, ValueError) as exc:
        print(json.dumps({"error": f"could not read input: {exc}",
                          "support_status": "OUT_OF_CONTRACT_INPUT"}), file=sys.stderr)
        return 3

    p = Gen1Predictor.load(a.artifact, a.meta, stage25_verdict=a.stage25_verdict)
    treatments = a.treatments.split(",") if a.treatments else None

    try:
        if a.rank:
            out = p.rank_conditions(x, b, component=a.component)
            print(json.dumps(out, indent=2))
            return 0 if out.get("scores") else 2
        rows = p.predict(x, b, treatments=treatments, component=a.component)
    except PredictionError as exc:
        print(json.dumps({"error": str(exc), "support_status": "OUT_OF_CONTRACT_INPUT"}),
              file=sys.stderr)
        return 3

    for r in rows:
        print(json.dumps(r))
    return 0 if all(r["support_status"] == "SUPPORTED_KNOWN_CONDITION" for r in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
