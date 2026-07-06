"""
Re-run ONLY the evaluation on an already-trained bundle (no rebuild, no retrain).

Your model is already trained in cellfate_run/bundle. This loads it + the existing
dataset and runs the evaluation that previously returned {} (fixed now), printing the
full gates + per-metric comparison against the baselines.

USAGE (from the repo root, env active):
    python evaluate_only.py
    python evaluate_only.py cellfate_run random     # (folder, regime) if different
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = sys.argv[1] if len(sys.argv) > 1 else "cellfate_run"
REGIME = sys.argv[2] if len(sys.argv) > 2 else "random"


def main() -> None:
    from cellfate.evaluation.evaluate_cli import EvalConfig, evaluate

    if not Path(ROOT, "bundle").is_dir():
        raise SystemExit(f"no trained bundle at {ROOT}/bundle — run training first")

    print(f"[eval] bundle={ROOT}/bundle  regime={REGIME}\n")
    gates = evaluate(EvalConfig(bundle=ROOT, dataset=ROOT, regimes=(REGIME,), out=f"{ROOT}/reports"))

    print("===== GATES (Document-1 success criteria) =====")
    print(json.dumps(gates, indent=2, default=str))
    if not gates:
        raise SystemExit("gates still empty — check that the dataset shards + splits exist under "
                         f"{ROOT} and the regime name matches the one used at build time")

    # detailed per-metric report (model vs baselines), written by evaluate()
    rep = Path(ROOT, "reports", f"{REGIME}.json")
    if rep.exists():
        R = json.loads(rep.read_text())
        print("\n===== MODEL vs BASELINES (held-out test split) =====")

        def prauc(d: dict) -> float:
            vals = [v for k, v in d.items() if k.startswith("prauc_")]
            return sum(vals) / len(vals) if vals else float("nan")

        estimators = [k for k in R if not k.startswith("_") and isinstance(R[k], dict)
                      and any(kk.startswith("prauc_") for kk in R[k])]
        print(f"{'estimator':<18}{'mean PR-AUC':>12}{'ECE':>8}{'reg_MAE':>10}")
        for name in ["model", *[e for e in estimators if e != "model"]]:
            if name not in R:
                continue
            d = R[name]
            print(f"{name:<18}{prauc(d):>12.3f}{d.get('ece', float('nan')):>8.3f}"
                  f"{d.get('reg_mae', float('nan')):>10.2f}")
        print(f"\ncoverage@0.90: {R.get('coverage', 'n/a')}"
              f"   ranking spearman: {R.get('ranking', {}).get('spearman', 'n/a')}")

    print("\nDONE. Paste the GATES + the MODEL vs BASELINES table back for interpretation.")


if __name__ == "__main__":
    main()
