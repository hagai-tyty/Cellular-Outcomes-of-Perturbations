"""
Retrain the six LOOCV bundles IN PLACE, so Stage 1's changes are actually measurable.

    python retrain_stage1.py --donors N2          # smoke test ONE fold first
    python retrain_stage1.py                      # all six  (Stage 1: A_xdonor)
    python retrain_stage1.py --no-xdonor          # rollback path (Stage 1a only)

WHY THIS EXISTS. `scorecard.py` does NOT train -- `measure_fold` calls `Predictor(root)`, which
loads the bundle already on disk. Every Stage 1 change is in the TRAINING path, so taking a
snapshot without retraining measures the OLD bundles and shows no change at all. That would look
like "Stage 1 did nothing" when in fact Stage 1 never ran.

WHY NOT `run_loocv.py`. That calls `run_multi_local.main()`, which does `shutil.rmtree(ROOT)` and
rebuilds each fold's dataset from raw GEO files -- harmonization, QC, clock, shards. Stage 1
changed none of that. This script reuses the existing shards/scalers/splits and redoes only
train -> calibrate -> bundle, which is exactly what changed.

THE CONFIG BELOW MIRRORS `run_multi_local.py` EXACTLY (d_cell=256, lr=1e-3, epochs=80,
ensemble_size=5, base_seed=0). Do not "improve" it: changing a hyperparameter here would mean the
comparison against `baseline` measures two changes at once, which violates the one-change rule in
REF_GROUND_RULES §2 and makes the result uninterpretable.

⚠ COST. Stage 1b trains one extra ensemble per training donor -- 5 inner + 1 deployed = 6x the
training work per fold, and the inner ensembles cannot be shrunk (their spread IS what
`sigma_scale` calibrates). Budget accordingly and run ONE fold first.

⚠ THE OLD BUNDLE IS OVERWRITTEN. Each fold's `bundle/` is copied to `bundle_pre_stage1/` first,
once, so the baseline stays reproducible. `scorecard/baseline.json` already holds the numbers.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

# Must be set BEFORE torch initialises CUDA, hence before any torch import. `set_global_seed`
# already requests deterministic algorithms, but on CUDA >= 10.2 cuBLAS GEMMs remain
# nondeterministic without this workspace setting -- torch warns about exactly that, and run 1
# printed the warning. Determinism is not cosmetic here: Stage 1's guards are supposed to come
# back BIT-IDENTICAL, a far sharper test than "the CI includes zero", and that only holds if the
# deployed ensemble is reproducible.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
BACKUP_DIRNAME = "bundle_pre_stage1"


def resolve_root(name: str) -> Path | None:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return p
    return None


def backup_bundle(root: Path) -> str:
    """Copy bundle/ -> bundle_pre_stage1/ once, so the pre-Stage-1 state stays reproducible."""
    src, dst = root / "bundle", root / BACKUP_DIRNAME
    if not src.is_dir():
        return "no bundle to back up"
    if dst.exists():
        return "backup already exists (kept)"
    shutil.copytree(src, dst)
    return f"backed up -> {dst.name}"


def retrain(root: Path, xdonor: bool, device: str) -> dict:
    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run

    # EXACTLY run_multi_local.py's config -- see the module docstring before touching this.
    cfg = TrainConfig(
        dataset_dir=str(root), out=str(root), regime=REGIME,
        d_cell=256, d_u=256, latent_dim=256, p_drop=0.2, lr=1e-3,
        epochs=80, patience=10, batch_size=256, ensemble_size=5,
        base_seed=0, conformal_levels=(0.90,), device=device,
        xdonor_calibration=xdonor, inference_mode="ensemble",
    )
    return train_run(cfg)


def main() -> None:
    ap = argparse.ArgumentParser(description="Retrain LOOCV bundles in place for Stage 1")
    ap.add_argument("--donors", default=",".join(DONORS),
                    help="comma-separated subset, e.g. N2 for a smoke test")
    ap.add_argument("--no-xdonor", action="store_true",
                    help="disable cross-donor calibration (the Stage 1a-only rollback path)")
    args = ap.parse_args()

    import torch

    from cellfate.common.console import install_pretty_console
    install_pretty_console()
    # cp1255 (Hebrew) console can't encode some glyphs; emit UTF-8 so a stray print
    # can never kill a multi-hour training run. (JSON results are written per fold.)
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    device = "cuda" if torch.cuda.is_available() else "cpu"
    xdonor = not args.no_xdonor
    donors = [d.strip() for d in args.donors.split(",") if d.strip()]

    print("\nRETRAIN for Stage 1 — reusing existing shards, redoing train -> calibrate -> bundle")
    print(f"  device               : {device}")
    print(f"  cross-donor calib    : {'ON  (Stage 1b -> tag A_xdonor)' if xdonor else 'OFF (Stage 1a only -> tag 1a_donorlabels)'}")
    print(f"  folds                : {', '.join(donors)}")
    if xdonor:
        print("  NOTE: ~6x the usual training time per fold (5 inner ensembles + 1 deployed).")
        print("        Run a single fold first (--donors N2) before committing to all six.")

    results, t0 = [], time.time()
    for i, d in enumerate(donors, 1):
        root = resolve_root(f"cellfate_loocv_{d}")
        if root is None:
            print(f"\n[{i}/{len(donors)}] {d}: SKIPPED — fold directory not found")
            results.append({"donor": d, "ok": False, "error": "directory not found"})
            continue
        if not any((root / "shards").glob("*.parquet")):
            print(f"\n[{i}/{len(donors)}] {d}: SKIPPED — no shards in {root/'shards'}")
            print("      This script reuses the dataset; it cannot rebuild it. Use run_loocv.py.")
            results.append({"donor": d, "ok": False, "error": "no shards"})
            continue

        print(f"\n{'=' * 78}\n[{i}/{len(donors)}] fold {d}  ({root})\n{'=' * 78}")
        print(f"  {backup_bundle(root)}")
        t = time.time()
        try:
            summary = retrain(root, xdonor, device)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED: {exc!r}")
            results.append({"donor": d, "ok": False, "error": repr(exc)[:200]})
            continue

        q = list(summary.get("conformal_q", {}).values())
        row = {
            "donor": d, "ok": True, "minutes": round((time.time() - t) / 60, 1),
            "temperature": summary.get("temperature"),
            "conformal_q": q[0] if q else None,
            "sigma_scale": summary.get("sigma_scale"),
            "xdonor_calibrated": summary.get("xdonor_calibrated"),
            "xdonor_n_donors": summary.get("xdonor_n_donors"),
            "xdonor_ece_before_temp": summary.get("xdonor_ece_before_temp"),
            "xdonor_ece_after_temp": summary.get("xdonor_ece_after_temp"),
        }
        results.append(row)
        print(f"  done in {row['minutes']} min | temperature {row['temperature']} | "
              f"q {row['conformal_q']} | sigma_scale {row['sigma_scale']} | "
              f"xdonor donors {row['xdonor_n_donors']}")
        Path("retrain_stage1_results.json").write_text(
            json.dumps(results, indent=2, default=str), encoding="utf-8")

    # ---- summary ----
    ok = [r for r in results if r.get("ok")]
    print(f"\n{'=' * 78}\nRETRAIN COMPLETE — {len(ok)}/{len(results)} folds in "
          f"{(time.time() - t0) / 60:.1f} min\n{'=' * 78}")
    if ok:
        print(f"{'fold':<8}{'temp':>8}{'q':>10}{'sigma_x':>10}{'donors':>8}"
              f"{'ECE pre':>10}{'ECE post':>10}")
        for r in ok:
            def _f(v, w=10, p=3):
                return f"{v:>{w}.{p}f}" if isinstance(v, (int, float)) else f"{'n/a':>{w}}"
            print(f"{r['donor']:<8}{_f(r['temperature'], 8)}{_f(r['conformal_q'])}"
                  f"{_f(r['sigma_scale'])}{str(r['xdonor_n_donors']):>8}"
                  f"{_f(r['xdonor_ece_before_temp'])}{_f(r['xdonor_ece_after_temp'])}")

        n_bad = [r["donor"] for r in ok if (r.get("xdonor_n_donors") or 0) != 5]
        if xdonor and n_bad:
            print(f"\n  !! folds with != 5 inner donors: {n_bad}")
            print("     Expected 5 (six donors minus the held-out one). Fewer means a thin")
            print("     inner-LODO pool; more means `cell_line` is finer-grained than donor and")
            print("     the residuals UNDERSTATE cross-donor error. Check verify_1a.py output.")

    print("\n  NEXT:")
    if xdonor:
        print("    python scorecard.py snapshot --tag A_xdonor")
        print("    python scorecard.py compare baseline A_xdonor")
    else:
        print("    python scorecard.py snapshot --tag 1a_donorlabels")
        print("    python scorecard.py compare baseline 1a_donorlabels    # expect ALL noise")
    print(f"\n  Pre-Stage-1 bundles are preserved in each fold's {BACKUP_DIRNAME}/ if you need")
    print("  to restore: remove bundle/ and rename that directory back.")


if __name__ == "__main__":
    main()
