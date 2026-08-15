"""Target-path audit + the SCALE-MISMATCH test.  (read-only; inference on ~20 cells/fold)

    python experiments/diag_target_path.py

TWO QUESTIONS
-------------
1. TARGET PATH. What does the network actually receive as `y_age`, per source, and is anything
   in the pipeline reconciling HFF's and Gill's scales?

2. SCALE MISMATCH (the hypothesis). `diag_target_shift` measured that C-7 halved HFF's target
   spread (SD 19.13 -> 8.77) while Gill's did not move (25.78 -> 27.65). HFF is ~99.8% of the
   age-valid training cells and the age loss is an unweighted mean over them, so the age head is
   fitted to HFF's spread and then scored on Gill's. If that is what hurt the ΔAge MAE, the
   model's predictions on the held-out donor should be COMPRESSED relative to the truth, and
   MORE compressed under C-7 than before it.

WHAT THE CODE PATH ALREADY TELLS US (read, not inferred -- reported by this script for the record)
  * `scalers.json` carries x_mean/x_std, dt_mean/dt_std, proliferation_coef -- and NO age scaler.
  * `training/dataset.py:58`  ya = np.where(am, arr["y_age"], 0.0)   <- raw YEARS, unscaled,
    while X and dose_time both go through `scalers.transform_*`.
  * `models/losses.py:58`     F.huber_loss(age_pred[m], age_true[m], delta=2.0)
    -- a mean over masked cells, with NO per-source weighting, and a delta fixed in YEARS.
  So nothing normalises the target and nothing reweights the sources: whatever scale HFF's labels
  have IS the scale the age head learns.

PRE-REGISTERED READING (constants below, fixed before running)
  compression = SD(pred) / SD(true) on the held-out donor.
  H-SUPPORTED   if C-7 compression < pre-C-7 compression in a MAJORITY of folds
                AND median C-7 compression <= COMPRESSION_CEILING.
  H-REFUTED     if C-7 compression is not below pre-C-7 in a majority of folds.
Reported alongside slope(pred ~ true), which is the same claim in regression form.

NOT CLAIMED: that fixing the scale would make the model correct, or that the pre-C-7 target was
the right one. C-7's justification is the degenerate control and is untouched by any of this.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))

from cellfate.common.io import ArtifactPaths  # noqa: E402
from cellfate.evaluation.baselines import ModelEstimator  # noqa: E402
from cellfate.evaluation.data import gather_split  # noqa: E402
from cellfate.inference import Predictor  # noqa: E402

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
ARMS = {"pre-C-7": "_armA", "C-7": "_c7t"}
REGIME = "holdout"

COMPRESSION_CEILING = 0.80    # SD(pred)/SD(true) at or below this is materially compressed


def sd(a: np.ndarray) -> float:
    return float(np.std(a, ddof=1)) if len(a) > 1 else float("nan")


def slope_of(true: np.ndarray, pred: np.ndarray) -> float:
    if len(true) < 3 or np.std(true) == 0:
        return float("nan")
    return float(np.polyfit(true, pred, 1)[0])


def hypothesis_verdict(comp_old: dict, comp_new: dict) -> dict:
    """The pre-registered rule, as a pure function so every branch is testable.

    A caveat recorded with the rule rather than after seeing the result: with 5 comparable folds
    the MAJORITY clause is nearly powerless -- P(>=3 of 5 in one direction | no effect) = 0.5, a
    coin flip. The load-bearing evidence is therefore the MAGNITUDE (the median), and the count is
    reported so a heterogeneous effect cannot hide inside it.
    """
    shared = sorted(set(comp_old) & set(comp_new))
    if not shared:
        return {"verdict": "UNDETERMINED", "folds": []}
    worse = [d for d in shared if comp_new[d] < comp_old[d]]
    med_old = float(np.median([comp_old[d] for d in shared]))
    med_new = float(np.median([comp_new[d] for d in shared]))
    supported = len(worse) > len(shared) / 2 and med_new <= COMPRESSION_CEILING
    return {"verdict": "H-SUPPORTED" if supported else "H-REFUTED", "folds": shared,
            "more_compressed": worse, "n_more_compressed": len(worse), "n_folds": len(shared),
            "median_compression_pre": med_old, "median_compression_c7": med_new}


def train_target_by_source(root: Path) -> dict:
    """The audit table: what the age head is actually fitted on, per source."""
    tr = gather_split(ArtifactPaths.of(str(root)), REGIME, "train")
    m = tr.mask.astype(bool)
    y, line = tr.y_age[m], np.asarray(tr.cell_line)[m]
    is_hff = line == "HFF"
    total = int(m.sum())
    out = {}
    for name, sel in (("HFF", is_hff), ("Gill", ~is_hff), ("combined", np.ones_like(is_hff))):
        v = y[sel]
        out[name] = {"n_age_valid": int(len(v)),
                     "weight": float(len(v) / total) if total else float("nan"),
                     "mean": float(v.mean()) if len(v) else float("nan"),
                     "sd": sd(v)}
    return out


def test_predictions(root: Path) -> dict:
    """ΔAge truth vs prediction on the held-out donor. ~20 cells, so this is seconds."""
    paths = ArtifactPaths.of(str(root))
    te = gather_split(paths, REGIME, "test")
    m = te.mask.astype(bool)
    if te.n == 0 or not m.any():
        return {"n": 0}
    est = ModelEstimator(Predictor(str(root / "bundle")))
    rows = est.rows(te.X, te.fp, te.dose_time)
    pred = np.array([r["mu_age"] for r in rows], dtype=np.float64)[m]
    true = np.asarray(te.y_age, dtype=np.float64)[m]
    s_t, s_p = sd(true), sd(pred)
    return {"n": int(m.sum()), "mean_true": float(true.mean()), "mean_pred": float(pred.mean()),
            "sd_true": s_t, "sd_pred": s_p,
            "compression": float(s_p / s_t) if s_t else float("nan"),
            "slope": slope_of(true, pred),
            "mae": float(np.abs(pred - true).mean())}


def main() -> None:
    res: dict = {"arms": {}, "folds": {}}
    print("=" * 104)
    print("TARGET-PATH AUDIT  --  what the age head is fitted on (TRAIN split, age-valid cells)")
    print("  target normalization: NONE (no age scaler in scalers.json; dataset.py:58 passes raw years)")
    print("  loss: Huber(delta=2.0 YEARS), unweighted mean over masked cells (models/losses.py:58)")
    print("=" * 104)
    for arm, sfx in ARMS.items():
        print(f"\n[{arm}]")
        print(f"  {'fold':<6}{'source':<10}{'n_age_valid':>13}{'weight':>9}{'mean':>10}{'sd':>10}")
        res["arms"][arm] = {}
        for d in DONORS:
            root = ROOT / f"cellfate_loocv_{d}{sfx}"
            if not (root / "shards").is_dir():
                continue
            t = train_target_by_source(root)
            res["arms"][arm][d] = t
            for src in ("HFF", "Gill", "combined"):
                s = t[src]
                print(f"  {d if src == 'HFF' else '':<6}{src:<10}{s['n_age_valid']:>13}"
                      f"{s['weight']:>9.4f}{s['mean']:>10.3f}{s['sd']:>10.3f}")

    print("\n" + "=" * 104)
    print("SCALE-MISMATCH TEST  --  held-out donor: is the prediction compressed vs the truth?")
    print("  pre-registered: H supported if C-7 compression < pre-C-7 in a majority of folds")
    print(f"                  AND median C-7 compression <= {COMPRESSION_CEILING}")
    print("=" * 104)
    print(f"  {'fold':<6}{'arm':<9}{'n':>4}{'sd_true':>10}{'sd_pred':>10}{'compress':>10}"
          f"{'slope':>9}{'mean_true':>11}{'mean_pred':>11}{'MAE':>9}")
    comp: dict[str, dict[str, float]] = {a: {} for a in ARMS}
    for d in DONORS:
        for arm, sfx in ARMS.items():
            root = ROOT / f"cellfate_loocv_{d}{sfx}"
            if not (root / "bundle").is_dir():
                continue
            try:
                r = test_predictions(root)
            except Exception as exc:                                    # noqa: BLE001
                print(f"  {d:<6}{arm:<9}  FAILED: {type(exc).__name__}: {exc}")
                continue
            res["folds"].setdefault(d, {})[arm] = r
            if r["n"] == 0:
                print(f"  {d:<6}{arm:<9}{0:>4}   (no age-valid test cell -- C-7 masks N2)")
                continue
            comp[arm][d] = r["compression"]
            print(f"  {d:<6}{arm:<9}{r['n']:>4}{r['sd_true']:>10.3f}{r['sd_pred']:>10.3f}"
                  f"{r['compression']:>10.3f}{r['slope']:>9.3f}{r['mean_true']:>11.3f}"
                  f"{r['mean_pred']:>11.3f}{r['mae']:>9.3f}")

    v = hypothesis_verdict(comp["pre-C-7"], comp["C-7"])
    if v["folds"]:
        print(f"\n  folds comparable: {v['folds']}")
        print(f"  median compression   pre-C-7 {v['median_compression_pre']:.3f}"
              f"  ->  C-7 {v['median_compression_c7']:.3f}")
        print(f"  more compressed under C-7 in {v['n_more_compressed']} of {v['n_folds']} "
              f"folds: {v['more_compressed']}")
        print(f"\n  VERDICT: {v['verdict']}")
        print("  CAVEAT (stated with the rule, not after the result): 3-of-5 is chance level as a "
              "sign\n  test, P=0.5 under the null. The magnitude carries this, not the count.")
        res["verdict"] = v

    _RESULTS.mkdir(exist_ok=True)
    out = _RESULTS / "diag_target_path_results.json"
    out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
