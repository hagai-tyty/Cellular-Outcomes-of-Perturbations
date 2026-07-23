"""Collect everything needed to choose the next calibrator OFFLINE, into a small sendable dump.

Run 3 missed `fate_ece` (0.249 vs the <=0.169 bar). Deciding what to change next needs more than
one in-sample number, and every further guess costs 3.8 h of GPU if it has to be answered by a
retrain. It does not have to be: `save_xstats` already persists the cross-donor pool for exactly
this reason, and the calib/test probabilities are one inference pass away.

This writes `diag_dump/<DONOR>.npz` + `diag_dump/manifest.json` (~2 MB total, zip and send).
It is READ-ONLY with respect to the run: no training, no bundle is modified, nothing is refitted.

    python dump_diag_bundle.py            # pool + calib + test  (a few minutes, needs torch)
    python dump_diag_bundle.py --pool-only  # pool only          (seconds, no torch, no splits)

WHAT IT ENABLES, and the discipline that goes with it
-----------------------------------------------------
The pool (`probs_mean`, ~103 rows, 5 donors) is the FITTING set: calibrator families and
hyper-parameters may be compared on it, honestly, by leave-one-donor-out WITHIN the pool.

The test arrays are the GRADED set. They are dumped so a chosen calibrator's effect on
`fate_ece` can be CONFIRMED once -- not so that several can be tried until one passes. Choosing on
them is choosing on the test set, which is the one forbidden move (`save_xstats`'s own warning,
and MASTER_PLAN's ground rules). The manifest records which arrays are which so that a later
reader cannot mistake one for the other.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
OUT = Path("diag_dump")

# Scalars worth carrying alongside the arrays so the dump is self-describing.
SCALAR_KEYS = (
    "temperature", "platt_a", "platt_b", "conformal_q", "sigma_scale", "sigma_scale_mc",
    "sigma_scale_mode", "xdonor_calibrated", "xdonor_n_donors", "fate_calib_n",
    "xdonor_only_platt_a", "xdonor_only_platt_b", "xdonor_only_n", "xdonor_only_n_donors",
    "xdonor_only_safe_ece_insample", "shipped_safe_ece_on_pool",
    "xdonor_safe_ece_before", "xdonor_safe_ece_after",
    "xdonor_ece_before_temp", "xdonor_ece_after_temp",
    "val_loss_mean", "n_train", "n_val", "n_calib",
)


def donor_ids_from_counts(counts: dict[int, int], n_rows: int) -> np.ndarray | None:
    """Rebuild a per-row donor label for the pool, or None if it cannot be trusted.

    The pool has no per-row donor column: `crossdonor_stats` appends one block per donor, in
    sorted donor order, and records only the COUNTS in `residuals_per_donor`. Leave-one-donor-out
    within the pool needs the labels, so they are reconstructed from those counts.

    This is only valid when every held-out cell was age-valid. Residual rows are masked by
    `am` while `logits`/`probs_mean` are NOT (xdonor_calib.py:326), so the two blocks can differ
    in length; when they do, the counts describe the residual rows and cannot index the fate
    rows. The totals matching is exactly the condition that rules that out -- per-donor
    `am.sum() <= len(inner_te)`, so equal totals force equality donor by donor.
    """
    if not counts:
        return None
    if sum(counts.values()) != n_rows:
        return None
    return np.concatenate([np.full(n, d, dtype=np.int32) for d, n in counts.items()])


def _pool(root: Path) -> dict:
    """Read the persisted cross-donor pool. Nothing here needs torch."""
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
    from cellfate.training.xdonor_calib import load_xstats

    xs = load_xstats(root / "bundle")
    n = len(xs.probs_mean)
    ids = donor_ids_from_counts(xs.residuals_per_donor, n)
    out = {
        "pool_probs_mean": xs.probs_mean,          # raw ensemble probs, == Predictor pbar at T=1
        "pool_targets": xs.targets,
        "pool_logits": xs.logits,
        "pool_abs_residuals": xs.abs_residuals,
        "pool_sigma_pred": xs.sigma_pred,
        "pool_sigma_pred_mc": xs.sigma_pred_mc,
    }
    if ids is not None:
        out["pool_donor_id"] = ids
    return out, {
        "pool_n": int(n),
        "pool_n_donors": int(xs.n_donors),
        "pool_residuals_per_donor": xs.residuals_per_donor,
        "pool_donor_scales": xs.donor_scales,
        "pool_donor_ids_reconstructed": ids is not None,
        "pool_donor_ids_note": (
            "per-row donor labels rebuilt from residuals_per_donor; totals matched"
            if ids is not None else
            "NOT reconstructible -- residual rows and fate rows differ in length, so "
            "leave-one-donor-out within the pool is unavailable for this fold"),
    }


def _splits(donor: str) -> tuple[dict, dict]:
    """Calib + test arrays, RAW and CALIBRATED, via the same loaders scorecard.py uses."""
    # Imported THROUGH scorecard rather than from cellfate directly: these must be the exact
    # objects the graded metric is computed with, and scorecard is where that set is defined.
    from scorecard import (  # noqa: PLC0415
        REGIME, ArtifactPaths, ModelEstimator, Predictor, gather_split, resolve_root,
    )

    root = resolve_root(f"cellfate_loocv_{donor}")
    paths = ArtifactPaths.of(root)
    te = gather_split(paths, REGIME, "test")
    cal = gather_split(paths, REGIME, "calib")
    pred = Predictor(root)

    def rows(split):
        r = ModelEstimator(pred).rows(split.X, split.fp, split.dose_time)
        return (np.array([x["S"] for x in r]), np.array([x["P_loss"] for x in r]),
                np.array([x["mu_age"] for x in r]), np.array([x["sigma_age"] for x in r]),
                np.array([x["in_dist"] for x in r]))

    s_te, pl_te, mu_te, sg_te, ind_te = rows(te)
    s_cal, pl_cal, _, _, _ = rows(cal)

    # ...and again with the bundle's calibration switched off IN MEMORY, so any candidate
    # calibrator can be fitted on the same raw probabilities the shipped one saw. Nothing is
    # written back; the object is discarded at the end of this function.
    shipped_platt = pred.platt
    pred.platt = None
    s_te_raw, pl_te_raw, _, _, _ = rows(te)
    s_cal_raw, pl_cal_raw, _, _, _ = rows(cal)
    pred.platt = shipped_platt

    arrays = {
        "test_S": s_te, "test_P_loss": pl_te, "test_mu_age": mu_te,
        "test_sigma_age": sg_te, "test_in_dist": ind_te,
        "test_S_raw": s_te_raw, "test_P_loss_raw": pl_te_raw,
        "test_y_cls": te.y_cls, "test_y_age": te.y_age, "test_mask": te.mask,
        "calib_S": s_cal, "calib_P_loss": pl_cal,
        "calib_S_raw": s_cal_raw, "calib_P_loss_raw": pl_cal_raw,
        "calib_y_cls": cal.y_cls,
    }
    meta = {
        "test_n": int(len(s_te)), "calib_n": int(len(s_cal)),
        "test_n_age_valid": int(np.asarray(te.mask).sum()),
        "bundle_had_platt": shipped_platt is not None,
        "conformal_q": float(pred.q),
    }
    return arrays, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool-only", action="store_true",
                    help="skip the inference pass; dump only the persisted cross-donor pool")
    args = ap.parse_args()

    from scorecard import resolve_root                  # noqa: PLC0415  (optional dependency)

    OUT.mkdir(exist_ok=True)
    manifest: dict[str, dict] = {}
    for d in DONORS:
        root = Path(resolve_root(f"cellfate_loocv_{d}"))
        arrays: dict[str, np.ndarray] = {}
        meta: dict = {"donor": d, "root": str(root)}

        # Per-fold error containment: one missing bundle must not discard the folds that worked.
        try:
            a, m = _pool(root)
            arrays.update(a)
            meta.update(m)
        except Exception as exc:  # noqa: BLE001
            meta["pool_error"] = repr(exc)[:160]

        if not args.pool_only:
            try:
                a, m = _splits(d)
                arrays.update(a)
                meta.update(m)
            except Exception as exc:  # noqa: BLE001
                meta["splits_error"] = repr(exc)[:160]

        try:
            mj = json.loads((root / "bundle" / "metrics.json").read_text(encoding="utf-8"))
            meta["metrics"] = {k: mj.get(k) for k in SCALAR_KEYS}
        except Exception as exc:  # noqa: BLE001
            meta["metrics_error"] = repr(exc)[:160]

        if arrays:
            np.savez_compressed(OUT / f"{d}.npz", **arrays)
        manifest[d] = meta
        bits = [k for k in ("pool_error", "splits_error", "metrics_error") if k in meta]
        print(f"  {d}: {len(arrays)} arrays"
              + (f"   [!] {', '.join(bits)}" if bits else "   ok"))

    manifest["_README"] = {
        "FITTING set (choose on these)": "pool_* -- the cross-donor pool, ~103 rows over 5 donors. "
            "Use leave-one-donor-out WITHIN the pool via pool_donor_id.",
        "GRADED set (confirm once, never choose)": "test_* -- the held-out donor. Choosing a "
            "calibrator by reading these is choosing on the test set.",
        "IN-DISTRIBUTION": "calib_* -- the calib split, same regime as training.",
        "raw vs calibrated": "*_raw are with the bundle's Platt switched off; the others are "
            "what the bundle actually ships. pool_probs_mean is always raw.",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str),
                                       encoding="utf-8")
    total = sum(f.stat().st_size for f in OUT.glob("*")) / 1e6
    print(f"\n  wrote {OUT}/ ({total:.1f} MB) -- zip this folder and send it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
