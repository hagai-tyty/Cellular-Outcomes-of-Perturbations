"""STAGE 15 — why RES is identically zero. Closing the last open item from the C-7 arc.

`res_median`, `res_max` and `res_approvals` have read 0.000 in every snapshot ever taken, and the
cause was never established. (An earlier "Spearman 0.40" headline for RES turned out to be a
correlation over floating-point residue, which is what made the zero worth chasing rather than
shrugging at.)

RES = phi(S) * S**k * g(R_eff) * exp(-lam * P_loss)

Four factors. Three of them CANNOT be zero on real inputs:

  phi(S) = sigmoid((S - tau_safe)/w)   a sigmoid is never exactly 0
  S**k                                 zero only if S is exactly 0
  exp(-lam * P_loss)                   lam ships at 0.0, so this is identically 1

so `g(R_eff) = R_eff/(R_eff + kappa)` is the ONLY factor that can produce a zero, and
`R_eff = max(0, -(mu_age + z_conf * sigma_age))`.

This script decomposes RES factor by factor on the real folds, so the answer is attributed rather
than inferred, and reports how far the system is from a non-zero RES in units of its own design
parameter.

Read-only. Runs inference on ~20 held-out cells per fold. Trains nothing, writes one results file.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage15_res_zero_results.json"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
SUFFIX = "_c7t"


def _sc():
    spec = importlib.util.spec_from_file_location("scorecard_mod", ROOT / "scorecard.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def factor_decomposition(S, P_loss, mu, sig, ind, p) -> dict:
    """Every factor of the RES product, so the zero can be ATTRIBUTED to one of them."""
    from cellfate.inference.res import _sigmoid
    phi = _sigmoid((S - p.tau_safe) / p.w)
    s_k = S ** p.k
    R_eff = np.maximum(0.0, -(mu + p.z_conf * sig))
    g = R_eff / (R_eff + p.kappa)
    loss_term = np.exp(-p.lam * P_loss)
    return {
        "phi_min": float(phi.min()), "phi_median": float(np.median(phi)),
        "s_k_min": float(s_k.min()), "s_k_median": float(np.median(s_k)),
        "loss_term_min": float(loss_term.min()), "loss_term_max": float(loss_term.max()),
        "g_min": float(g.min()), "g_max": float(g.max()),
        "n_phi_zero": int((phi == 0).sum()),
        "n_s_k_zero": int((s_k == 0).sum()),
        "n_g_zero": int((g == 0).sum()),
        "n_cells": int(len(S)),
    }


def headroom(mu: np.ndarray, sig: np.ndarray, z_conf: float) -> dict:
    """How far from a non-zero RES, stated three ways.

    `min_upper_bound` is the closest any cell comes to the credit threshold: R_eff > 0 needs
    `mu + z*sigma < 0`, so this number is the shortfall in years.

    `z_required` is the same distance expressed in the design parameter: the largest `z` at which
    SOME cell would still qualify. If it is below the shipped `z_conf`, the gate is the binding
    constraint; how far below says by how much.
    """
    ub = mu + z_conf * sig
    rej = mu < 0
    z_req = float(np.max(-mu[rej] / sig[rej])) if rej.any() and (sig[rej] > 0).all() else None
    return {
        "min_upper_bound": float(ub.min()),
        "n_upper_bound_negative": int((ub < 0).sum()),
        "n_mu_negative": int(rej.sum()),
        "z_required_for_any_cell": z_req,
        "z_conf_shipped": float(z_conf),
        "sigma_median": float(np.median(sig)),
        "sigma_max": float(sig.max()),
        "mu_min": float(mu.min()),
        "abs_mu_median": float(np.median(np.abs(mu))),
        "sigma_over_abs_mu_median": (float(np.median(sig) / np.median(np.abs(mu)))
                                     if np.median(np.abs(mu)) > 0 else None),
    }


def run(donors=None) -> dict:
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.baselines import ModelEstimator
    from cellfate.evaluation.data import gather_split
    from cellfate.inference import Predictor, compute_res_batch
    sc = _sc()
    per_fold, errs = {}, {}
    params = None
    for d in donors or DONORS:
        root = sc.resolve_root(f"cellfate_loocv_{d}{SUFFIX}")
        try:
            paths = ArtifactPaths.of(root)
            te = gather_split(paths, "holdout", "test")
            pred = Predictor(root)
        except Exception as exc:                              # noqa: BLE001
            errs[d] = repr(exc)[:120]
            continue
        p = pred.res_params
        params = {"tau_safe": p.tau_safe, "w": p.w, "k": p.k,
                  "kappa": p.kappa, "z_conf": p.z_conf, "lam": p.lam}
        rows = ModelEstimator(pred).rows(te.X, te.fp, te.dose_time)
        S = np.array([r["S"] for r in rows])
        P_loss = np.array([r["P_loss"] for r in rows])
        mu = np.array([r["mu_age"] for r in rows])
        sig = np.array([r["sigma_age"] for r in rows])
        ind = np.array([r["in_dist"] for r in rows])
        res, status = compute_res_batch(S, P_loss, mu, sig, ind, p)
        uniq, cnt = np.unique(status, return_counts=True)
        per_fold[d] = {
            "res_max": float(np.max(res)), "res_median": float(np.median(res)),
            "res_all_zero": bool(np.all(res == 0.0)),
            "factors": factor_decomposition(S, P_loss, mu, sig, ind, p),
            "headroom": headroom(mu, sig, p.z_conf),
            "status_counts": dict(zip([str(u) for u in uniq], [int(c) for c in cnt],
                                      strict=True)),
            "n_out_of_distribution": int((~ind).sum()),
        }
    return {"params": params, "folds": per_fold, "errors": errs}


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                     # noqa: BLE001
            pass
    r = run()
    print("\nSTAGE 15 — why RES is identically zero")
    print(f"  RES params: {r['params']}")
    if r["errors"]:
        print(f"  folds unavailable: {r['errors']}")

    print("\n  WHICH FACTOR IS ZERO (RES = phi * S^k * g(R_eff) * exp(-lam*P_loss))")
    print(f"     {'fold':<6}{'n':>4}{'phi>0':>8}{'S^k>0':>8}{'exp()':>8}{'g=0 for':>10}"
          f"{'  verdict'}")
    for d, f in r["folds"].items():
        x = f["factors"]
        print(f"     {d:<6}{x['n_cells']:>4}"
              f"{x['n_cells'] - x['n_phi_zero']:>8}{x['n_cells'] - x['n_s_k_zero']:>8}"
              f"{x['loss_term_min']:>8.2f}{x['n_g_zero']:>7}/{x['n_cells']:<3}"
              f"  {'g(R_eff) alone' if x['n_g_zero'] == x['n_cells'] else 'MIXED'}")

    print("\n  WHY R_eff = 0:  it needs mu + z*sigma < 0, and sigma dwarfs mu")
    print(f"     {'fold':<6}{'mu min':>9}{'|mu| med':>10}{'sigma med':>11}{'sigma/|mu|':>12}"
          f"{'min(mu+z*s)':>13}{'z needed':>10}")
    for d, f in r["folds"].items():
        h = f["headroom"]
        zr = h["z_required_for_any_cell"]
        print(f"     {d:<6}{h['mu_min']:>9.2f}{h['abs_mu_median']:>10.2f}"
              f"{h['sigma_median']:>11.2f}{h['sigma_over_abs_mu_median']:>12.2f}"
              f"{h['min_upper_bound']:>13.2f}{(zr if zr is not None else float('nan')):>10.3f}")
    z_ship = r["params"]["z_conf"]
    print(f"     'z needed' is the largest z_conf at which ANY cell would qualify. Shipped: "
          f"{z_ship}.")

    print("\n  STATUS BREAKDOWN")
    for d, f in r["folds"].items():
        print(f"     {d:<6}{f['status_counts']}")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
