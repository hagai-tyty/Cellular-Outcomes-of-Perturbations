"""
Test 11.1 (ΔAge lab notebook) — ISOLATE THE TIME DIMENSION.  [revised after code review]

Question (Test 11 correction): does the varying TIME carry ΔAge signal beyond cell state, and does
the TRAINED NET actually use it?  (A) time redundant with state (benign) vs (B) net under-uses a
real signal (a fixable bug).

CAVEAT THAT SHAPES THE DESIGN: each held-out donor's data is its reprogramming TIME-COURSE (~21
samples), so within a donor, expression state ≈ a readout of progress ≈ time — state and time are
COLLINEAR. The data-side statistics (Parts 0-2) therefore can barely separate "time redundant" from
"time under-used": a near-zero partial-R² is expected whenever state proxies time, regardless of
whether the net wastes signal. So the DECISIVE test is Part 3 (intervene on the model's time input
directly); Parts 0-2 are diagnostics that report how well the data can even answer.

Data side — pooled over all 6 donors' held-out time-courses (~126 samples), leave-one-donor-out:
  Part 0  collinearity: R²(time ~ state). High => state and time are near-interchangeable.
  Part 1  variance decomposition: R²(ΔAge~state), R²(ΔAge~time), partial R²(time | state).
  Part 2  state-controlled kNN: within expression-neighborhoods, the TIME SPREAD (can we hold state
          and still move time?) and the local ΔAge~time corr where spread allows.
Model side — per fold (needs each fold's trained net):
  Part 3  DECISIVE — hold X + fingerprint fixed, sweep the time input 10th->90th pct, measure the
          net's ΔAge shift. (Deliberately moves time off the state-time manifold; fine for a pure
          sensitivity test of the model's use of the time channel.)

VERDICT (Part 3 leads; Parts 0-1 disambiguate the "net ignores time" case):
  - net sweep RESPONDS                            -> the net USES the time input; Test 11's u_only
                                                     failure was the linear ablation, not the net.
  - net sweep ~ 0 AND high collinearity/partial~0 -> (A) time redundant with state; net rightly rides
                                                     state (which already encodes time). Benign.
  - net sweep ~ 0 BUT partial-R²(time|state) > 0  -> (B) net wastes a real, state-independent time
                                                     signal. Fixable: strengthen the time encoder.

USAGE (repo root, venv active; place in experiments/. needs cellfate_loocv_* bundles + splits):
    python test11_1_time_isolation.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
K_NN = 10
SWEEP_LO, SWEEP_HI = 0.10, 0.90


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def load_donor_timecourses():
    """Each donor's own held-out 'test' split = its reprogramming time-course (age-valid cells)."""
    per = {}
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            te = gather_split(ArtifactPaths.of(root), REGIME, "test")
        except Exception:  # noqa: BLE001
            continue
        m = te.mask
        if m.sum() < 3 or np.std(te.dose_time[m, 1]) == 0:
            continue
        per[d] = {"X": np.asarray(te.X[m], float),
                  "time": np.asarray(te.dose_time[m, 1], float),
                  "y": np.asarray(te.y_age[m], float)}
    return per


def _lodo_pred(feats, targets, donors):
    """Leave-one-donor-out ridge: out-of-donor predictions per donor."""
    out = {}
    for hd in donors:
        tr = [d for d in donors if d != hd]
        Xtr = np.vstack([feats[d] for d in tr])
        ytr = np.concatenate([targets[d] for d in tr])
        if np.std(ytr) == 0:
            out[hd] = np.full(len(targets[hd]), float(np.mean(ytr)))
            continue
        sc = StandardScaler().fit(Xtr)
        out[hd] = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr).predict(sc.transform(feats[hd]))
    return out


def _agg_r2(pred, targets, donors):
    r2s = [r2_score(targets[d], pred[d]) for d in donors
           if len(targets[d]) >= 3 and np.std(targets[d]) > 0]
    return float(np.mean(r2s)) if r2s else float("nan")


def data_side(per):
    donors = list(per)
    Xf = {d: per[d]["X"] for d in donors}
    Tf = {d: per[d]["time"].reshape(-1, 1) for d in donors}
    y = {d: per[d]["y"] for d in donors}
    time = {d: per[d]["time"] for d in donors}

    r2_collin = _agg_r2(_lodo_pred(Xf, time, donors), time, donors)          # Part 0
    pred_state = _lodo_pred(Xf, y, donors)                                   # Part 1
    r2_state = _agg_r2(pred_state, y, donors)
    r2_time = _agg_r2(_lodo_pred(Tf, y, donors), y, donors)
    resid = {d: y[d] - pred_state[d] for d in donors}
    r2_partial = _agg_r2(_lodo_pred(Tf, resid, donors), resid, donors)

    Xall = np.vstack([per[d]["X"] for d in donors])                          # Part 2
    tall = np.concatenate([per[d]["time"] for d in donors])
    yall = np.concatenate([per[d]["y"] for d in donors])
    Xs = StandardScaler().fit_transform(Xall)
    k = min(K_NN, len(Xall) - 1)
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(Xs).kneighbors(Xs)
    spreads, locals_ = [], []
    for i in range(len(Xall)):
        nbr = idx[i, 1:]
        spreads.append(float(np.std(tall[nbr])))
        if np.std(tall[nbr]) > 0 and np.std(yall[nbr]) > 0:
            locals_.append(float(np.corrcoef(tall[nbr], yall[nbr])[0, 1]))

    return {"r2_collin": r2_collin, "r2_state": r2_state, "r2_time": r2_time,
            "r2_partial": r2_partial, "nbr_spread": float(np.mean(spreads)),
            "overall_spread": float(np.std(tall)),
            "local_corr": float(np.mean(locals_)) if locals_ else float("nan"), "n": int(len(Xall))}


def model_side():
    """Per fold: hold X + fingerprint fixed, sweep time input, measure the net's ΔAge shift."""
    per = {}
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            te = gather_split(ArtifactPaths.of(root), REGIME, "test")
        except Exception:  # noqa: BLE001
            continue
        m = te.mask
        tt = te.dose_time[m, 1]
        if m.sum() < 3 or np.std(tt) == 0:
            continue
        est = ModelEstimator(Predictor(root))
        X, fp, dt = te.X[m], te.fp[m], te.dose_time[m]
        t_lo, t_hi = np.quantile(tt, SWEEP_LO), np.quantile(tt, SWEEP_HI)
        dt_lo, dt_hi = dt.copy(), dt.copy()
        dt_lo[:, 1] = t_lo
        dt_hi[:, 1] = t_hi
        age_lo = np.array([r["mu_age"] for r in est.rows(X, fp, dt_lo)], float)
        age_hi = np.array([r["mu_age"] for r in est.rows(X, fp, dt_hi)], float)
        per[d] = float(np.mean(np.abs(age_hi - age_lo)))
    return per


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 11.1 — ISOLATE THE TIME DIMENSION (does the net use the time input?)")
    print("time = dose_time[:,1] = log(time_h). Part 3 is decisive; Parts 0-2 say if the data can separate.")

    per = load_donor_timecourses()
    if len(per) < 2:
        print("\n   Not enough donor time-courses found (need cellfate_loocv_*/test splits).")
        return
    ds = data_side(per)
    ms = model_side()

    ratio = ds["nbr_spread"] / ds["overall_spread"] if ds["overall_spread"] > 0 else float("nan")
    print(f"\n  DATA SIDE (pooled {ds['n']} donor samples, leave-one-donor-out):")
    print(f"    Part 0  collinearity  R2(time ~ state) = {ds['r2_collin']:+.3f}   "
          f"(high => state and time interchangeable)")
    print(f"    Part 1  R2(dAge~state)={ds['r2_state']:+.3f}  R2(dAge~time)={ds['r2_time']:+.3f}  "
          f"partial R2(time|state)={ds['r2_partial']:+.3f}")
    print(f"    Part 2  kNN nbhd time spread={ds['nbr_spread']:.3f} vs overall {ds['overall_spread']:.3f} "
          f"(ratio {ratio:.2f})  local dAge~time r={ds['local_corr']:+.3f}")

    print("\n  MODEL SIDE (per fold: net dAge shift over a 10th->90th pct time sweep, X held fixed):")
    if ms:
        print(render_table(["fold", "net dAge sweep (yr)"], [[d, f"{ms[d]:.3f}"] for d in ms],
                           aligns=["l", "r"]))
        net_sweep = float(np.mean(list(ms.values())))
        print(f"    aggregate net dAge sweep = {net_sweep:.3f} yr")
    else:
        net_sweep = float("nan")
        print("    (no trained folds loaded — model side unavailable)")

    print("\n   VERDICT:")
    responds = np.isfinite(net_sweep) and net_sweep > 0.5
    high_collin = ds["r2_collin"] > 0.5
    partial_signal = np.isfinite(ds["r2_partial"]) and ds["r2_partial"] > 0.05
    if responds:
        print("     NET USES TIME — its dAge shifts when the time input is swept, so Test 11's u_only")
        print("     failure reflects the LINEAR ablation, not the net ignoring time. (Not a redundancy proof.)")
    elif np.isfinite(net_sweep):
        if partial_signal and not high_collin:
            print("     (B) NET UNDER-USES TIME — time carries state-independent dAge signal (partial-R2>0,")
            print("     state does NOT fully predict time) yet the net's dAge barely moves. Fixable bug.")
        else:
            print("     (A) TIME REDUNDANT WITH STATE — state ~ time (collinear) and little residual time")
            print("     signal, so the net riding state loses nothing. Benign; nothing to fix.")
        print("     NOTE: with state/time collinear here, (A) is hard to distinguish from 'data cannot")
        print("     separate' — treat as suggestive; confirm with varied-time-at-fixed-state data.")
    else:
        print("     Model side unavailable; the data side alone cannot decide (collinearity). Re-run with")
        print("     the trained folds present to get the decisive Part 3.")


if __name__ == "__main__":
    main()
