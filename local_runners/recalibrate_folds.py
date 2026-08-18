"""STAGE 16 -- re-fit the shipped fate calibrator on the HARD class, without retraining.

`train_model.py` now fits Platt on `argmax(y_cls)` instead of the soft probability, but Platt is
fitted DURING TRAINING, so every bundle already on disk still carries soft-fitted coefficients.
This re-runs *only* that step against the corrected target, so the change can be evaluated on real
artefacts instead of asserted from unit tests.

WHAT IT REPRODUCES. Exactly the shipped calibration block:

    p_all = concat( ensemble P(safe) on calib , xdonor pool P(safe) )
    y_all = concat( hard(calib target)        , hard(xdonor target) )
    a, b  = fit_platt_binary(p_all, y_all)

The cross-donor pool comes from the bundle's own `xdonor_stats.npz` (`probs_mean`, `targets`), so
it is the same 100 cells the training run used -- not a re-derivation.

WHY THE CALIB PROBABILITIES ARE OBTAINED BY INVERSION. `platt_safe(p,a,b) = sigmoid(a*logit(p)+b)`
is exactly invertible, so the raw ensemble probability is recovered from the served `S` as
`sigmoid((logit(S) - b)/a)`. That avoids rebuilding the training Dataset just to redo a forward
pass, and the round-trip is asserted to 1e-9 before anything is fitted.

SAFETY. Writes into `<root>_s16`, created by HARDLINKING the `_s12` fold (so the 242 MB of shards
cost nothing) and then REPLACING `bundle/temperature.json` -- the old file is unlinked first, so
the original fold is untouched. `_s12` and the `c7t_stage12` snapshot taken from it stay valid.

    python local_runners/recalibrate_folds.py
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
SRC_SUFFIX = "_s12"
DST_SUFFIX = "_s16"


def _logit(p):
    p = np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    return np.log(p / (1 - p))


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(x, float), -500, 500)))


def clone_fold(src: Path, dst: Path) -> None:
    """Hardlink the fold, then unlink the one file we are going to replace."""
    if dst.exists():
        shutil.rmtree(dst)
    # `cp -al` keeps 242 MB of shards at zero extra cost. Fall back to a real copy if the
    # filesystem refuses hardlinks.
    r = subprocess.run(["cp", "-al", str(src), str(dst)], capture_output=True, text=True)
    if r.returncode != 0:
        shutil.copytree(src, dst)
    # unlink (never truncate) so the ORIGINAL bundle's file is not modified through the link
    (dst / "bundle" / "temperature.json").unlink()


def recalibrate(donor: str) -> dict:
    from cellfate.common.calibration import platt_safe
    from cellfate.common.constants import SAFE_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.baselines import ModelEstimator
    from cellfate.evaluation.data import gather_split
    from cellfate.inference import Predictor
    from cellfate.training.calibrate import fit_platt_binary

    src = ROOT / f"cellfate_loocv_{donor}{SRC_SUFFIX}"
    paths = ArtifactPaths.of(str(src))
    pred = Predictor(str(src))
    if pred.platt is None:
        raise SystemExit(f"{donor}: bundle carries no Platt coefficients; nothing to recalibrate")
    a1, b1 = pred.platt

    cal = gather_split(paths, "holdout", "calib")
    S = np.array([r["S"] for r in ModelEstimator(pred).rows(cal.X, cal.fp, cal.dose_time)])
    raw = _sigmoid((_logit(S) - b1) / a1)
    # the inversion must be exact, or the refit is being done on the wrong quantity
    err = float(np.max(np.abs(platt_safe(raw, a1, b1) - S)))
    if err > 1e-9:
        raise SystemExit(f"{donor}: Platt inversion round-trip failed ({err:.2e})")

    y_cal_hard = (cal.y_cls.astype(int) == SAFE_IDX).astype(np.float64)

    z = np.load(src / "bundle" / "xdonor_stats.npz", allow_pickle=True)
    p_xd = np.asarray(z["probs_mean"], float)[:, SAFE_IDX]
    y_xd = (np.argmax(np.asarray(z["targets"], float), axis=1) == SAFE_IDX).astype(np.float64)

    p_all = np.concatenate([raw, p_xd])
    y_all = np.concatenate([y_cal_hard, y_xd])
    a2, b2 = fit_platt_binary(p_all, y_all)

    dst = ROOT / f"cellfate_loocv_{donor}{DST_SUFFIX}"
    clone_fold(src, dst)
    (dst / "bundle" / "temperature.json").write_text(
        json.dumps({"temperature": 1.0, "platt_a": a2, "platt_b": b2}, indent=2), encoding="utf-8")
    (dst / "bundle" / "recalibration.json").write_text(json.dumps({
        "stage": 16, "target": "hard", "source_fold": src.name,
        "shipped_platt_a": a1, "shipped_platt_b": b1,
        "recalibrated_platt_a": a2, "recalibrated_platt_b": b2,
        "n_calib": int(len(raw)), "n_xdonor": int(len(p_xd)),
        "calib_hard_safe_frac": float(y_cal_hard.mean()),
        "xdonor_hard_safe_frac": float(y_xd.mean()),
    }, indent=2), encoding="utf-8")
    return {"donor": donor, "a_old": a1, "b_old": b1, "a_new": a2, "b_new": b2,
            "n_calib": int(len(raw)), "n_xdonor": int(len(p_xd)),
            "calib_hard_safe": float(y_cal_hard.mean())}


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    print("\nSTAGE 16 -- re-fitting the fate calibrator on the HARD class (no retraining)")
    print(f"  {SRC_SUFFIX} -> {DST_SUFFIX}; shards hardlinked, only bundle/temperature.json replaced")
    print(f"\n     {'fold':<6}{'a (shipped)':>13}{'b':>8}{'a (hard)':>11}{'b':>8}"
          f"{'n_calib':>9}{'n_xd':>6}{'calib safe':>12}")
    rows = []
    for d in DONORS:
        r = recalibrate(d)
        rows.append(r)
        print(f"     {d:<6}{r['a_old']:>13.4f}{r['b_old']:>8.4f}{r['a_new']:>11.4f}"
              f"{r['b_new']:>8.4f}{r['n_calib']:>9}{r['n_xdonor']:>6}{r['calib_hard_safe']:>12.3f}")
    out = ROOT / "results" / "stage16_recalibration_results.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\n  saved -> {out.relative_to(ROOT)}")
    print(f"  next: evaluate against {DST_SUFFIX} and compare to the {SRC_SUFFIX} baseline")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
