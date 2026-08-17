"""STAGE 14 PRE-FLIGHT — what would actually happen if we adopted a calibrated ΔAge target?

Stage 11 found the dense clock was never broken, only mis-scaled: LODO-calibrated raw ΔAge lands
at MAE 6.78 against a 7.30 yr methylation floor, with `k` stable across donors (spread 1.19x).
The obvious next move is to adopt `y_age -> k * y_age`. This script exists because that move has
a trap in it, and the trap is the kind this project has already fallen into twice.

**THE TRAP.** ΔAge is the regression target. Multiply it by a constant and every ΔAge metric on
the scorecard moves by that constant BY ARITHMETIC: MAE, level shift and conformal width all
scale by k, while Spearman is exactly unchanged. With k ~ 0.37 that is a 63 % "improvement" in
ΔAge MAE representing **no gain in model skill whatsoever**. A scorecard compare run after such a
rebuild would print a large ACCEPT.

**WHAT MAKES IT NOT PURELY A UNITS CHANGE.** `huber_age_loss` uses `huber_delta = 2.0`
(`configs/train/default.yaml:17`) — a knee fixed in ABSOLUTE YEARS. Shrinking the target moves the
residual distribution relative to that knee, so the loss changes character (more of it quadratic),
and `MultiTaskLoss`'s learned `log_var_age` re-weights age against fate. Those effects are real
and their sign is not knowable in advance.

So this script measures the two things that decide how to pre-register the Change:

  E1  EXACT EQUIVARIANCE, demonstrated on the linear path. Ridge with an L2 loss must satisfy
      pred(k*y) == k*pred(y) exactly. If it does, the "units change" claim is established on real
      data rather than argued from algebra.
  E2  HOW FAR the neural path is from that, quantified as the fraction of residuals inside the
      Huber knee before and after rescaling. Near 0 both times -> the loss is effectively L1 and
      the change is close to pure units. Materially larger after -> the optimisation genuinely
      changes and the rebuild can move for reasons unrelated to scale.

Read-only. Reads built folds, refits only a cheap ridge, writes one results file. Trains nothing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage14_calibration_equivariance_results.json"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

# Stage 11, `results/diag_stage11_scale_results.json`, variant `raw`, mean of the LODO k's.
K_LS = 0.3637          # least squares -- best MAE, but SD ratio 0.597 (it SHRINKS)
K_VAR = 0.5991         # variance-matched -- SD ratio 0.990, preserves the spread
HUBER_DELTA = 2.0      # configs/train/default.yaml:17, in YEARS
DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
SUFFIX = "_c7t"


def _sc():
    import importlib.util
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def ridge_equivariance(tr, te, k: float) -> dict:
    """E1: refit the SAME ridge on k*y and check it is exactly k times the original.

    This is the load-bearing demonstration. If a linear model's predictions scale exactly, then
    for the linear part of the system a calibrated target buys precisely nothing beyond units.
    """
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)

    def feats(s):
        return np.hstack([sx.transform(s.X), np.asarray(s.fp, float),
                          sdt.transform(s.dose_time)])

    ftr, fte = feats(tr), feats(te)
    y = tr.y_age[tr.mask]
    p1 = Ridge(alpha=1.0).fit(ftr[tr.mask], y).predict(fte)
    pk = Ridge(alpha=1.0).fit(ftr[tr.mask], k * y).predict(fte)

    m = te.mask
    yt = te.y_age[m]
    mae1 = float(np.abs(p1[m] - yt).mean())
    maek = float(np.abs(pk[m] - k * yt).mean())
    sc = _sc()
    return {
        "max_abs_dev_from_k_times_pred": float(np.max(np.abs(pk - k * p1))),
        "rel_dev": float(np.max(np.abs(pk - k * p1)) / max(np.max(np.abs(k * p1)), 1e-12)),
        "mae_unscaled": mae1, "mae_scaled": maek,
        "mae_ratio": maek / mae1 if mae1 else float("nan"),
        "k": k,
        "spearman_unscaled": sc._sp(-p1[m], -yt),
        "spearman_scaled": sc._sp(-pk[m], -(k * yt)),
    }


def huber_region(residuals: np.ndarray, k: float, delta: float = HUBER_DELTA) -> dict:
    """E2: what fraction of the loss sits inside the quadratic knee, before and after rescaling.

    The knee is fixed in years, so this is the ONLY channel through which a pure rescale can
    change what the neural model optimises.
    """
    r = np.abs(np.asarray(residuals, float))
    return {"frac_inside_before": float((r < delta).mean()),
            "frac_inside_after": float((r * k < delta).mean()),
            "median_abs_resid_before": float(np.median(r)),
            "median_abs_resid_after": float(np.median(r) * k),
            "delta": delta}


def predicted_scorecard_effect(fold_metrics: dict, k: float) -> dict:
    """What a rebuild on k*y MUST produce if it is a pure rescale. This is the guard: any
    deviation from these numbers means the rebuild did something other than change units."""
    return {
        "dage_mae_model": fold_metrics["dage_mae_model"] * k,
        "level_shift_model": fold_metrics["level_shift_model"] * k,
        "conformal_width": fold_metrics["conformal_width"] * k,
        "rank_model_dage": fold_metrics["rank_model_dage"],      # EXACTLY unchanged
    }


def run(donors=None, k: float = K_LS) -> dict:
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split
    sc = _sc()
    donors = donors or DONORS
    per_fold, errs = {}, {}
    for d in donors:
        root = sc.resolve_root(f"cellfate_loocv_{d}{SUFFIX}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, "holdout", "train")
            te = gather_split(paths, "holdout", "test")
        except Exception as exc:                              # noqa: BLE001
            errs[d] = repr(exc)[:120]
            continue
        if te.mask.sum() < 3:
            errs[d] = "too few age-valid cells"
            continue
        eq_ls = ridge_equivariance(tr, te, k)
        eq_var = ridge_equivariance(tr, te, K_VAR)
        # Residuals of the ridge stand in for "a model's residuals on this fold"; the Huber
        # question is about the SIZE of residuals relative to a fixed knee, not about which
        # estimator produced them.
        from sklearn.linear_model import Ridge
        from sklearn.preprocessing import StandardScaler
        sx = StandardScaler().fit(tr.X)
        sdt = StandardScaler().fit(tr.dose_time)
        f = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float),
                       sdt.transform(tr.dose_time)])
        yy = tr.y_age[tr.mask]
        resid = Ridge(alpha=1.0).fit(f[tr.mask], yy).predict(f[tr.mask]) - yy
        per_fold[d] = {"equivariance_k_ls": eq_ls, "equivariance_k_var": eq_var,
                       "huber": huber_region(resid, k),
                       "y_sd": float(np.std(yy, ddof=1)), "n_train_age": int(len(yy))}

    # The units-effect guard, computed off the committed snapshot rather than a fresh run.
    snap_path = ROOT / "scorecard" / "c7_A_keep_hff.json"
    predicted = {}
    if snap_path.exists():
        folds = json.loads(snap_path.read_text(encoding="utf-8"))["folds"]
        for d, f in folds.items():
            if isinstance(f, dict) and "_error" not in f and f.get("dage_mae_model") is not None:
                predicted[d] = predicted_scorecard_effect(f, k)

    return {"k_ls": K_LS, "k_var": K_VAR, "huber_delta": HUBER_DELTA,
            "folds": per_fold, "errors": errs,
            "predicted_scorecard_under_pure_rescale": predicted}


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                     # noqa: BLE001
            pass
    r = run()
    print("\nSTAGE 14 PRE-FLIGHT — what adopting a calibrated ΔAge would actually do")
    print(f"  k_ls {r['k_ls']}   k_var {r['k_var']}   huber_delta {r['huber_delta']} yr")
    if r["errors"]:
        print(f"  folds unavailable: {r['errors']}")

    print("\n  E1 — EXACT EQUIVARIANCE of the linear path (ridge, L2)")
    print(f"     {'fold':<6}{'max|pred_k - k*pred|':>22}{'MAE ratio':>12}{'k':>8}"
          f"{'Δrho':>10}")
    for d, f in r["folds"].items():
        e = f["equivariance_k_ls"]
        drho = (e["spearman_scaled"] - e["spearman_unscaled"]
                if None not in (e["spearman_scaled"], e["spearman_unscaled"]) else float("nan"))
        print(f"     {d:<6}{e['max_abs_dev_from_k_times_pred']:>22.2e}{e['mae_ratio']:>12.4f}"
              f"{e['k']:>8.4f}{drho:>10.2e}")
    print("     MAE ratio == k and Δrho == 0 means the rescale bought UNITS, nothing else.")

    print("\n  E2 — how far the neural path can deviate: residuals vs the fixed Huber knee")
    print(f"     {'fold':<6}{'med|resid|':>12}{'-> after':>10}{'% inside knee':>15}"
          f"{'-> after':>10}{'y SD':>9}")
    for d, f in r["folds"].items():
        h = f["huber"]
        print(f"     {d:<6}{h['median_abs_resid_before']:>12.2f}"
              f"{h['median_abs_resid_after']:>10.2f}{h['frac_inside_before']:>14.1%}"
              f"{h['frac_inside_after']:>10.1%}{f['y_sd']:>9.2f}")

    print("\n  GUARD for the eventual rebuild — a PURE rescale must produce exactly these:")
    print(f"     {'fold':<6}{'ΔAge MAE':>11}{'level shift':>13}{'conformal width':>17}")
    for d, p in r["predicted_scorecard_under_pure_rescale"].items():
        print(f"     {d:<6}{p['dage_mae_model']:>11.2f}{p['level_shift_model']:>13.2f}"
              f"{p['conformal_width']:>17.2f}")
    print("     rank_model_dage must be EXACTLY unchanged. Any deviation from this table means")
    print("     the rebuild did something other than change units -- investigate before keeping.")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
