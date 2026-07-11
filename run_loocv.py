"""
Leave-one-donor-out cross-validation for CellFate-Rx.

Rotates the held-out Gill donor through all six (N2, N3, O1, O2, Y1, Y2), running
the full harmonized leave-cell-line-out pipeline once per fold, and reports the
DISTRIBUTION of results (mean +/- std) instead of a single held-out number. This
hardens the generalization claim against the "n=1 held-out donor" objection.

Each fold rebuilds the dataset from scratch, refitting harmonization statistics on
that fold's training donors only (the held-out donor never contributes) -- so this
is compute-heavy: ~6 full builds. Expect a few hours; run it overnight.

USAGE (repo root, venv active):
    python run_loocv.py "D:\\GSE242423" "D:\\Gill"
"""
from __future__ import annotations

import json
import statistics as stats
import sys
import traceback
from pathlib import Path

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]


def _prauc(est: dict) -> float:
    vals = [v for k, v in est.items() if k.startswith("prauc_")]
    return sum(vals) / len(vals) if vals else float("nan")


def main() -> None:
    gse_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
    gill_dir = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"

    # import the existing pipeline and drive it fold by fold
    import importlib.util
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("run_multi_local", here / "run_multi_local.py")
    rml = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rml)
    rml.GSE_DIR, rml.GILL_DIR = gse_dir, gill_dir
    rml.HARMONIZE = True

    Path("loocv_results").mkdir(exist_ok=True)
    folds: list[dict] = []

    for i, donor in enumerate(DONORS, 1):
        print("\n" + "=" * 78)
        print(f"[LOOCV fold {i}/{len(DONORS)}]  held-out donor = {donor}")
        print("=" * 78)
        rml.HOLDOUT_DONOR = donor
        rml.ROOT = f"cellfate_loocv_{donor}"
        try:
            rml.main()
        except SystemExit as e:
            print(f"[fold {donor}] SKIPPED: {e}")
            folds.append({"donor": donor, "ok": False, "error": str(e)})
            continue
        except Exception:
            print(f"[fold {donor}] FAILED:\n{traceback.format_exc()}")
            folds.append({"donor": donor, "ok": False, "error": "exception"})
            continue

        # collect this fold's metrics from the report + bundle
        rep = Path(rml.ROOT, "reports", "holdout.json")
        row: dict = {"donor": donor, "ok": True}
        if rep.exists():
            R = json.loads(rep.read_text())
            model = R.get("model", {})
            row.update({
                "spearman": R.get("ranking", {}).get("spearman", float("nan")),
                "reg_mae": model.get("reg_mae", float("nan")),
                "prauc": _prauc(model),
                "ece": model.get("ece", float("nan")),
                "coverage": R.get("coverage", float("nan")),
                "ridge_mae": R.get("ridge", {}).get("reg_mae", float("nan")),
            })
        bundle_metrics = Path(rml.ROOT, "bundle", "metrics.json")
        if bundle_metrics.exists():
            m = json.loads(bundle_metrics.read_text())
            row.update({"n_train": m.get("n_train"), "n_val": m.get("n_val"),
                        "n_calib": m.get("n_calib"), "temp": m.get("temperature")})
        folds.append(row)
        Path("loocv_results", "folds.json").write_text(json.dumps(folds, indent=2, default=str))

    # ---- aggregate ---------------------------------------------------------- #
    ok = [f for f in folds if f.get("ok")]
    print("\n" + "=" * 78)
    print("LEAVE-ONE-DONOR-OUT CROSS-VALIDATION  —  per-fold results")
    print("=" * 78)
    print(f"{'held-out':<10}{'Spearman':>10}{'reg_MAE':>10}{'ridge_MAE':>11}"
          f"{'PR-AUC':>9}{'ECE':>8}{'cover':>8}")
    for f in folds:
        if not f.get("ok"):
            print(f"{f['donor']:<10}{'  (fold skipped/failed)':<50}")
            continue
        print(f"{f['donor']:<10}{_g(f,'spearman'):>10}{_g(f,'reg_mae'):>10}"
              f"{_g(f,'ridge_mae'):>11}{_g(f,'prauc'):>9}{_g(f,'ece'):>8}{_g(f,'coverage'):>8}")

    def agg(key: str):
        xs = [f[key] for f in ok if isinstance(f.get(key), (int, float))
              and f[key] == f[key]]   # drop nan
        if not xs:
            return None
        m = stats.mean(xs)
        s = stats.stdev(xs) if len(xs) > 1 else 0.0
        return m, s, min(xs), max(xs), len(xs)

    print("\n" + "-" * 78)
    print("AGGREGATE across folds (mean +/- std [min, max], n valid folds)")
    for label, key in [("Spearman (ranking)", "spearman"), ("reg_MAE (model)", "reg_mae"),
                       ("reg_MAE (ridge)", "ridge_mae"), ("ECE", "ece"), ("coverage@0.90", "coverage")]:
        a = agg(key)
        if a:
            m, s, lo, hi, n = a
            print(f"  {label:<22} {m:+.3f} +/- {s:.3f}   [{lo:+.3f}, {hi:+.3f}]   (n={n})")
        else:
            print(f"  {label:<22} (no valid folds)")

    # headline line for the paper
    sp = agg("spearman")
    if sp:
        print(f"\n>>> Ranking generalizes across held-out donors: "
              f"Spearman {sp[0]:.2f} +/- {sp[1]:.2f} (n={sp[4]} donors) <<<")
    Path("loocv_results", "summary.json").write_text(
        json.dumps({"folds": folds,
                    "aggregate": {k: agg(k) for k in
                                  ["spearman", "reg_mae", "ridge_mae", "ece", "coverage"]}},
                   indent=2, default=str))
    print("\nSaved: loocv_results/folds.json, loocv_results/summary.json")


def _g(f: dict, key: str) -> str:
    v = f.get(key)
    return f"{v:+.3f}" if isinstance(v, (int, float)) and v == v else "   nan"


if __name__ == "__main__":
    main()
