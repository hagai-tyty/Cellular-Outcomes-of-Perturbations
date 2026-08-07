"""
DIAGNOSTIC — why is arm A's RES score CONSTANT (Spearman = nan) on 3 of 6 folds?

The reproduction test (experiments/repro_test7_res_armA.py, outcome O3) found that
model_dAge and ridge_dAge reproduce the July numbers EXACTLY on all six folds, while
model_RES collapsed: N3, O2 and Y1 return Spearman = nan, which `ranking_metrics`
emits only when the score has zero variance. So RES is not merely worse on those
folds -- it is a CONSTANT, and a constant carries no ranking at all.

RES (src/cellfate/inference/res.py:59-63):

    R_eff = max(0, -(mu_age + z_conf * sigma_age))
    g     = R_eff / (R_eff + kappa)
    phi   = sigmoid((S - tau_safe) / w)
    res   = phi * S**k * g * exp(-lam * P_loss)
    res   = where(in_dist, res, 0)                     <-- OOD gate

There are exactly two ways every cell in a fold can land on the same value, and they
have completely different meanings:

  M1  OOD GATE. `in_dist` is False for every test cell -> res = 0 everywhere. Means the
      held-out donor looks out-of-distribution to the model; RES is switched off by
      design and the constant is the gate working, not RES failing.
  M2  NO-REJUVENATION FLOOR. `mu + z*sigma >= 0` for every cell -> R_eff = 0 -> g = 0 ->
      res = 0 everywhere. Means the model's uncertainty sigma_age has grown (or mu_age
      shifted up) enough to swallow the whole predicted rejuvenation. This is a real
      degradation and would be caused by something in the 1.5.3-1.5.6 changes.

They can also co-occur. This script separates them per fold, and reports the same
quantities for the folds that still produce a finite RES so the contrast is visible.

It also reports, per fold, what fraction of cells each `status` code takes, since
status names the mechanism directly (REJECTED_OOD vs REJECTED_NO_REJUVENATION).

READ:
  - M2 dominant  -> sigma_age inflation is the cause. RES's degeneracy is downstream of
                    the uncertainty head, not of the OOD gate. This is the escalation
                    path: something between July and arm A widened sigma_age.
  - M1 dominant  -> the OOD gate is firing on whole held-out donors. RES is off by
                    design; the July finite values mean the gate did NOT fire in July,
                    so the OOD calibration changed.
  - neither      -> the constant comes from phi/S**k/exp(-lam*P_loss) saturating; look
                    at the fate head.

READ-ONLY. Touches no build, writes results/diag_res_degenerate_armA_results.json.

USAGE (repo root, venv active):
    python experiments/diag_res_degenerate_armA.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(REPO / "experiments"))

from cellfate.common.io import ArtifactPaths  # noqa: E402
from cellfate.evaluation.baselines import ModelEstimator  # noqa: E402
from cellfate.evaluation.data import gather_split  # noqa: E402
from cellfate.inference import Predictor, compute_res_batch  # noqa: E402

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
ARM = "armA"
# folds whose RES came back nan in the reproduction test
DEGENERATE = {"N3", "O2", "Y1"}
OUT = _RESULTS / "diag_res_degenerate_armA_results.json"


def one_fold(donor: str) -> dict | None:
    root = str(REPO / f"cellfate_loocv_{donor}_{ARM}")
    paths = ArtifactPaths.of(root)
    te = gather_split(paths, REGIME, "test")
    m = te.mask
    if m.sum() < 3:
        return None

    pred = Predictor(root)
    rows = ModelEstimator(pred).rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    mu = np.array([r["mu_age"] for r in rows])
    sig = np.array([r["sigma_age"] for r in rows])
    ind = np.array([r["in_dist"] for r in rows])
    res, status = compute_res_batch(S, P_loss, mu, sig, ind, pred.res_params)

    p = pred.res_params
    R_eff = np.maximum(0.0, -(mu + p.z_conf * sig))
    g = R_eff / (R_eff + p.kappa)
    phi = 1.0 / (1.0 + np.exp(-(S - p.tau_safe) / p.w))

    # restrict to the age-valid cells the ranking actually scored
    sl = m
    resm = res[sl]
    n = int(sl.sum())

    m1 = float((~ind[sl]).mean())                 # OOD-gated fraction
    m2 = float((R_eff[sl] == 0.0).mean())         # no-rejuvenation fraction
    both = float(((~ind[sl]) & (R_eff[sl] == 0.0)).mean())

    st, cnt = np.unique(status[sl], return_counts=True)
    status_frac = {str(k): float(v) / n for k, v in zip(st, cnt, strict=True)}

    return {
        "donor": donor,
        "n_age_valid": n,
        "res_constant": bool(np.std(resm) == 0.0),
        "res_n_unique": int(np.unique(np.round(resm, 12)).size),
        "res_frac_zero": float((resm == 0.0).mean()),
        "res_min": float(resm.min()), "res_max": float(resm.max()),
        "M1_frac_ood": m1,
        "M2_frac_no_rejuv": m2,
        "frac_both": both,
        "status_frac": status_frac,
        "mu_age": {"mean": float(mu[sl].mean()), "min": float(mu[sl].min()),
                   "max": float(mu[sl].max())},
        "sigma_age": {"mean": float(sig[sl].mean()), "min": float(sig[sl].min()),
                      "max": float(sig[sl].max())},
        "mu_plus_z_sigma": {"mean": float((mu + p.z_conf * sig)[sl].mean()),
                            "frac_ge_0": float(((mu + p.z_conf * sig)[sl] >= 0).mean())},
        "R_eff": {"mean": float(R_eff[sl].mean()), "max": float(R_eff[sl].max())},
        "g": {"mean": float(g[sl].mean()), "max": float(g[sl].max())},
        "phi": {"mean": float(phi[sl].mean()), "min": float(phi[sl].min())},
        "S": {"mean": float(S[sl].mean()), "min": float(S[sl].min()),
              "max": float(S[sl].max())},
        "P_loss": {"mean": float(P_loss[sl].mean())},
        "res_params": {"z_conf": p.z_conf, "kappa": p.kappa, "tau_safe": p.tau_safe,
                       "w": p.w, "k": p.k, "lam": p.lam},
    }


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\n" + "=" * 78)
    print("DIAGNOSTIC — why RES is constant on 3/6 arm-A folds")
    print("=" * 78)

    per = {}
    for d in DONORS:
        r = one_fold(d)
        if r is not None:
            per[d] = r
    if not per:
        print("   no folds found")
        return 2

    rows = []
    for d in DONORS:
        if d not in per:
            continue
        r = per[d]
        rows.append([
            d,
            "nan(const)" if r["res_constant"] else "finite",
            str(r["res_n_unique"]),
            f"{r['res_frac_zero']:.0%}",
            f"{r['M1_frac_ood']:.0%}",
            f"{r['M2_frac_no_rejuv']:.0%}",
        ])
    print("\n" + render_table(
        ["fold", "RES", "uniq", "%zero", "M1 %OOD", "M2 %no-rejuv"],
        rows, aligns=["l", "l", "r", "r", "r", "r"]))

    rows = []
    for d in DONORS:
        if d not in per:
            continue
        r = per[d]
        rows.append([
            d,
            f"{r['mu_age']['mean']:+.2f}",
            f"{r['sigma_age']['mean']:.2f}",
            f"{r['mu_plus_z_sigma']['mean']:+.2f}",
            f"{r['mu_plus_z_sigma']['frac_ge_0']:.0%}",
            f"{r['R_eff']['max']:.2f}",
            f"{r['g']['mean']:.3f}",
        ])
    print("\n  the R_eff chain   (R_eff = max(0, -(mu + z*sigma)),  z = "
          f"{per[next(iter(per))]['res_params']['z_conf']})")
    print(render_table(
        ["fold", "mean mu", "mean sigma", "mean mu+z·sig", "%>=0", "max R_eff", "mean g"],
        rows, aligns=["l", "r", "r", "r", "r", "r", "r"]))

    print("\n  status composition (the gate names its own mechanism)")
    keys = sorted({k for r in per.values() for k in r["status_frac"]})
    rows = [[d] + [f"{per[d]['status_frac'].get(k, 0.0):.0%}" for k in keys]
            for d in DONORS if d in per]
    print(render_table(["fold"] + keys, rows, aligns=["l"] + ["r"] * len(keys)))

    deg = [d for d in per if per[d]["res_constant"]]
    fin = [d for d in per if not per[d]["res_constant"]]
    print(f"\n   constant-RES folds: {deg or 'none'}")
    print(f"   finite-RES   folds: {fin or 'none'}")
    print(f"   (reproduction test reported nan on {sorted(DEGENERATE)})")

    print("\n   MECHANISM:")
    if deg:
        m1 = float(np.mean([per[d]["M1_frac_ood"] for d in deg]))
        m2 = float(np.mean([per[d]["M2_frac_no_rejuv"] for d in deg]))
        print(f"     on the constant folds: OOD-gated {m1:.0%}   no-rejuvenation {m2:.0%}")
        if m2 > 0.99 and m1 < 0.99:
            print("     -> M2. R_eff = 0 for EVERY cell: mu + z*sigma >= 0 throughout, so the")
            print("        predicted rejuvenation is entirely swallowed by the uncertainty")
            print("        margin. RES degeneracy is downstream of sigma_age, not the OOD gate.")
        elif m1 > 0.99 and m2 < 0.99:
            print("     -> M1. Every test cell is flagged out-of-distribution; RES is switched")
            print("        off by the gate. The OOD calibration changed since July.")
        elif m1 > 0.99 and m2 > 0.99:
            print("     -> BOTH fire on every cell; they are not separable from this fold alone.")
        else:
            print("     -> neither M1 nor M2 is universal; the constant must come from phi,")
            print("        S**k or exp(-lam*P_loss) saturating. Inspect the fate head.")
    else:
        print("     no constant folds reproduced here — investigate the harness, not the model.")

    payload = {"arm": ARM, "per_fold": per,
               "constant_folds": deg, "finite_folds": fin}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
