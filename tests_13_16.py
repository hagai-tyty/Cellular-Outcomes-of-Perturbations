"""
TESTS 13–16 — the complete pre-change battery, in ONE run.

    python tests_13_16.py

Everything that must be known BEFORE changing any model code. Four tests, one command:

  TEST 13  ΔAge TRAJECTORY SHAPE + the Test-16 gate
           Is the trajectory biphasic (stress then rejuvenation) or monotonic? And — the part
           that actually gates the plan — is the per-donor error a CONSTANT offset or does it
           vary with time? Test 16 assumes constant; if it is not, a scalar correction is
           MISSPECIFIED and would distort trajectories instead of fixing them.

  TEST 14  CONFORMAL INTERVAL VALIDATION  (never measured before)
           Do the uncertainty intervals cover at their nominal level? Are they narrow enough to
           be informative? RES consumes sigma_age inside R_eff, so broken uncertainty is an
           independent upstream cause of RES's collapse.

  TEST 15  OOD DETECTOR VALIDATION  (never measured before)
           It flags ~27% of held-out cells and OOD zeroes RES outright. Does flagging actually
           track higher prediction error, or is it firing indiscriminately on anything foreign?

  TEST 16  PER-DONOR CALIBRATION FEASIBILITY  ***the gate for the main fix***
           With k labelled reference cells from a held-out donor, how well can the level shift be
           estimated? Sweep k = 1, 3, 5, 10, with repeated random draws. Two correction variants:
             (a) SCALAR   — one offset from k random cells
             (b) MATCHED  — offset from the k cells nearest in TIME to each target cell
           If (b) clearly beats (a), the error is time-varying and the fix must be time-aware.
           Correction cells are ALWAYS excluded from evaluation — no leakage.

PRE-REGISTERED CRITERIA (from MASTER_PLAN.md §7b, set before running):
  Test 16 PASS  : |level shift| reduced >=50% on >=4/6 folds at k <= 5  -> implement the fix
  Test 16 BORDER: only achieved at k >= 10                              -> document cost, decide
  Test 16 FAIL  : not achieved even at k = 10                           -> STOP; within-donor
                                                                           ranker only

HONEST POWER CAVEAT: ~21 cells per donor. Treat "not significant" as NOT DEMONSTRATED, never as
"disproven".

USAGE (repo root, venv active; needs cellfate_loocv_* bundles with calib/test splits).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
K_SWEEP = [1, 3, 5, 10]
N_REPEATS = 40
SEED = 0
T_CRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 10: 2.228, 15: 2.131,
          18: 2.101, 20: 2.086}


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def _tcrit(df: int) -> float:
    if df <= 0:
        return float("inf")
    for k in sorted(T_CRIT):
        if df <= k:
            return T_CRIT[k]
    return 1.96


def ols(X: np.ndarray, y: np.ndarray):
    n, p = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = n - p
    if df <= 0:
        return beta, np.full(p, np.nan), float("nan")
    s2 = float(resid @ resid) / df
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return beta, np.full(p, np.nan), float("nan")
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        tvals = beta / se
    ss = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss if ss > 0 else float("nan")
    return beta, tvals, r2


def fit_shape(t: np.ndarray, y: np.ndarray):
    """Linear + quadratic fit of y on standardised time; classify the shape."""
    ts = (t - t.mean()) / (t.std() + 1e-9)
    _, t_lin, _ = ols(np.column_stack([np.ones_like(ts), ts]), y)
    b_q, t_q, r2q = ols(np.column_stack([np.ones_like(ts), ts, ts ** 2]), y)
    crit_q, crit_l = _tcrit(len(ts) - 3), _tcrit(len(ts) - 2)
    quad_sig = bool(np.isfinite(t_q[2]) and abs(t_q[2]) > crit_q)
    vertex = -b_q[1] / (2 * b_q[2]) if abs(b_q[2]) > 1e-12 else np.inf
    interior = bool(np.isfinite(vertex) and ts.min() < vertex < ts.max())
    if quad_sig and b_q[2] < 0 and interior:
        shape = "HUMP (rise then fall)"
    elif quad_sig and b_q[2] > 0 and interior:
        shape = "U (fall then rise)"
    elif np.isfinite(t_lin[1]) and abs(t_lin[1]) > crit_l:
        shape = "monotonic " + ("up" if t_lin[1] > 0 else "down")
    else:
        shape = "no clear trend"
    return {"slope_t": float(t_lin[1]), "quad_t": float(t_q[2]), "r2_quad": r2q,
            "shape": shape}


def load_fold(donor: str):
    """Everything the four tests need from one fold."""
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        pred = Predictor(root)
    except Exception:  # noqa: BLE001
        return None
    m = te.mask
    if m.sum() < 6:
        return None
    rows = ModelEstimator(pred).rows(te.X, te.fp, te.dose_time)
    mu = np.array([r["mu_age"] for r in rows], float)
    ind = np.array([r["in_dist"] for r in rows], bool)
    # ridge ΔAge on the pipeline's feature view
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(te.X), np.asarray(te.fp, float), sdt.transform(te.dose_time)])
    r_te = Ridge(alpha=1.0).fit(ftr[tr.mask], tr.y_age[tr.mask]).predict(fte)
    return {"mu": mu[m], "ridge": r_te[m], "true": np.asarray(te.y_age, float)[m],
            "time": np.asarray(te.dose_time[:, 1], float)[m], "in_dist": ind[m],
            "mu_all": mu, "ind_all": ind, "true_all": np.asarray(te.y_age, float),
            "mask": m, "q": float(pred.q), "level": float(pred.conformal_level),
            "n": int(m.sum())}


# --------------------------------------------------------------------------- #
# TEST 13 — trajectory shape + the Test-16 gate
# --------------------------------------------------------------------------- #
def test13(F, render_table):
    print("\n" + "=" * 78)
    print("TEST 13 — ΔAge TRAJECTORY SHAPE, and is the per-donor error time-varying?")
    print("=" * 78)
    print("Part 3 is the gate: a CURVED residual misspecifies a scalar offset as badly as a")
    print("SLOPED one, so the gate fires on EITHER a significant slope OR a significant curve.")

    res = {}
    for d, f in F.items():
        if np.std(f["time"]) == 0:
            continue
        res[d] = {"true": fit_shape(f["time"], f["true"]),
                  "pred": fit_shape(f["time"], f["mu"]),
                  "resid": fit_shape(f["time"], f["mu"] - f["true"]),
                  "rng": (float(f["true"].min()), float(f["true"].max()))}
    if not res:
        print("   no folds with varying time")
        return {}

    print("\n  PART 1 — TRUE ΔAge vs time (is the biphasic claim supported on OUR data?)")
    print(render_table(["fold", "true ΔAge range", "curve t", "R2 quad", "shape"],
                       [[d, f"[{res[d]['rng'][0]:+.1f},{res[d]['rng'][1]:+.1f}]",
                         f"{res[d]['true']['quad_t']:+.2f}", f"{res[d]['true']['r2_quad']:.2f}",
                         res[d]["true"]["shape"]] for d in res],
                       aligns=["l", "r", "r", "r", "l"]))

    print("\n  PART 2 — PREDICTED ΔAge vs time (does the model reproduce the shape?)")
    print(render_table(["fold", "curve t", "R2 quad", "shape"],
                       [[d, f"{res[d]['pred']['quad_t']:+.2f}", f"{res[d]['pred']['r2_quad']:.2f}",
                         res[d]["pred"]["shape"]] for d in res], aligns=["l", "r", "r", "l"]))

    print("\n  PART 3 — RESIDUAL (pred - true) vs time   ***THE TEST-16 GATE***")
    rows = []
    for d in res:
        r = res[d]["resid"]
        sl, qd = abs(r["slope_t"]) > 2.1, abs(r["quad_t"]) > 2.1
        why = ",".join([x for x, c in (("slope", sl), ("curve", qd)) if c])
        rows.append([d, f"{r['slope_t']:+.2f}", f"{r['quad_t']:+.2f}",
                     f"TIME-VARYING ({why})" if (sl or qd) else "~constant"])
    print(render_table(["fold", "slope t", "curve t", "error is..."], rows,
                       aligns=["l", "r", "r", "l"]))

    n_hump = sum(1 for d in res if "HUMP" in res[d]["true"]["shape"])
    n_tv = sum(1 for d in res if abs(res[d]["resid"]["slope_t"]) > 2.1
               or abs(res[d]["resid"]["quad_t"]) > 2.1)
    print(f"\n   biphasic (hump) TRUE trajectories: {n_hump}/{len(res)}")
    if n_hump == 0:
        print("   => BIPHASIC CLAIM NOT SUPPORTED here. Do NOT build the plan on it.")
    elif n_hump >= len(res) / 2:
        print("   => BIPHASIC SUPPORTED. RES's absolute per-cell gating is biologically")
        print("      misspecified (it scores trajectory points in isolation) — say so in the writeup.")
    else:
        print("   => MIXED: donor-specific, not a universal protocol effect.")
    print(f"\n   time-varying per-donor error: {n_tv}/{len(res)}")
    if n_tv >= 2:
        print("   => TEST 16 MUST BE TIME-AWARE. A single scalar offset is misspecified.")
        print("      (Test 16 below reports a MATCHED-time variant precisely for this case.)")
    else:
        print("   => TEST 16's scalar design is well-specified.")
    return {"n_hump": n_hump, "n_tv": n_tv, "n": len(res)}


# --------------------------------------------------------------------------- #
# TEST 14 — conformal interval validation
# --------------------------------------------------------------------------- #
def test14(F, render_table):
    print("\n" + "=" * 78)
    print("TEST 14 — CONFORMAL INTERVAL VALIDATION (never measured before)")
    print("=" * 78)
    print("coverage = fraction of TRUE ΔAge inside [mu-q, mu+q]. Should match the nominal level.")
    rows, covs = [], []
    for d, f in F.items():
        lo, hi = f["mu"] - f["q"], f["mu"] + f["q"]
        cov = float(((f["true"] >= lo) & (f["true"] <= hi)).mean())
        covs.append(cov)
        width = 2.0 * f["q"]
        spread = float(np.std(f["true"]))
        ratio = width / (2 * spread) if spread > 0 else float("nan")
        rows.append([d, f"{f['level']:.2f}", f"{cov:.2f}", f"{cov - f['level']:+.2f}",
                     f"{width:.1f}", f"{ratio:.2f}",
                     "UNDER-covers" if cov < f["level"] - 0.05
                     else ("over-covers" if cov > f["level"] + 0.05 else "calibrated")])
    print(render_table(["fold", "nominal", "coverage", "gap", "width(yr)", "width/2sd", "verdict"],
                       rows, aligns=["l", "r", "r", "r", "r", "r", "l"]))
    mean_cov = float(np.mean(covs)) if covs else float("nan")
    nominal = float(np.mean([f["level"] for f in F.values()]))
    print(f"\n   mean coverage = {mean_cov:.2f} vs nominal {nominal:.2f}")
    if mean_cov < nominal - 0.05:
        print("   => INTERVALS UNDER-COVER (overconfident). sigma_age is too small out-of-donor.")
        print("      This is an INDEPENDENT upstream cause of RES's R_eff collapse — R_eff uses")
        print("      mu + z*sigma, so understated sigma distorts the rejuvenation credit.")
    elif mean_cov > nominal + 0.05:
        print("   => INTERVALS OVER-COVER (too wide). Honest but uninformative; check width/2sd —")
        print("      a ratio >> 1 means the interval spans more than the data's own spread.")
    else:
        print("   => INTERVALS ARE CALIBRATED. Uncertainty is not the problem; rule it out.")
    return {"mean_cov": mean_cov, "nominal": nominal}


# --------------------------------------------------------------------------- #
# TEST 15 — OOD detector validation
# --------------------------------------------------------------------------- #
def test15(F, render_table):
    print("\n" + "=" * 78)
    print("TEST 15 — OOD DETECTOR VALIDATION (never measured before)")
    print("=" * 78)
    print("Does flagging track HIGHER error, or fire indiscriminately? OOD zeroes RES outright.")
    rows, aucs = [], []
    for d, f in F.items():
        ind, err = f["in_dist"], np.abs(f["mu"] - f["true"])
        rate = float((~ind).mean())
        if ind.all() or (~ind).all():
            rows.append([d, f"{rate:.2f}", "n/a", "n/a", "n/a", "single class"])
            continue
        e_flag, e_ok = float(err[~ind].mean()), float(err[ind].mean())
        auc = float(roc_auc_score((~ind).astype(int), err))   # does error predict flagging?
        aucs.append(auc)
        rows.append([d, f"{rate:.2f}", f"{e_flag:.1f}", f"{e_ok:.1f}",
                     f"{auc:.2f}", "informative" if auc > 0.6 else
                     ("anti-informative" if auc < 0.4 else "uninformative")])
    print(render_table(["fold", "OOD rate", "MAE flagged", "MAE kept", "AUC", "verdict"],
                       rows, aligns=["l", "r", "r", "r", "r", "l"]))
    mauc = float(np.mean(aucs)) if aucs else float("nan")
    print(f"\n   mean AUC(error -> flagged) = {mauc:.2f}   (0.5 = flags are random wrt error)")
    if np.isfinite(mauc) and mauc > 0.6:
        print("   => OOD DETECTOR WORKS: flagged cells genuinely have higher error. Zeroing RES")
        print("      on them is defensible.")
    elif np.isfinite(mauc) and mauc < 0.45:
        print("   => OOD DETECTOR IS MISLEADING: flagged cells are NOT the erroneous ones, yet it")
        print("      zeroes their RES. This is an independent cause of RES's collapse — FIX OR")
        print("      DISABLE the OOD gate.")
    else:
        print("   => OOD DETECTOR IS UNINFORMATIVE (~chance wrt error). It discards ~27% of cells")
        print("      for no measurable benefit. Consider disabling it in the RES path.")
    return {"mean_auc": mauc}


# --------------------------------------------------------------------------- #
# TEST 16 — per-donor calibration feasibility  (THE GATE FOR THE MAIN FIX)
# --------------------------------------------------------------------------- #
def test16(F, render_table):
    print("\n" + "=" * 78)
    print("TEST 16 — PER-DONOR CALIBRATION FEASIBILITY  ***THE GATE FOR THE MAIN FIX***")
    print("=" * 78)
    print("With k labelled reference cells from the held-out donor, estimate the level shift and")
    print("correct the REMAINING cells (reference cells excluded from evaluation - no leakage).")
    print(f"{N_REPEATS} random draws per k. SCALAR = one global offset; MATCHED = offset from the")
    print("k reference cells nearest in TIME (handles a time-varying error).")

    rng = np.random.default_rng(SEED)
    out = {}
    for d, f in F.items():
        mu, y, tt = f["mu"], f["true"], f["time"]
        n = len(y)
        base_shift = abs(float(np.median(mu) - np.median(y)))
        base_mae = float(np.abs(mu - y).mean())
        per_k = {}
        for k in K_SWEEP:
            if n - k < 3:
                continue
            sc_shift, sc_mae, mt_shift, mt_mae = [], [], [], []
            for _ in range(N_REPEATS):
                idx = rng.choice(n, k, replace=False)
                ev = np.setdiff1d(np.arange(n), idx)
                # (a) SCALAR: one offset from the k reference cells
                off = float(np.median(mu[idx] - y[idx]))
                c = mu[ev] - off
                sc_shift.append(abs(float(np.median(c) - np.median(y[ev]))))
                sc_mae.append(float(np.abs(c - y[ev]).mean()))
                # (b) MATCHED: per-cell offset from the nearest reference cell in time
                near = idx[np.argmin(np.abs(tt[ev][:, None] - tt[idx][None, :]), axis=1)]
                cm = mu[ev] - (mu[near] - y[near])
                mt_shift.append(abs(float(np.median(cm) - np.median(y[ev]))))
                mt_mae.append(float(np.abs(cm - y[ev]).mean()))
            per_k[k] = {"sc_shift": float(np.mean(sc_shift)), "sc_mae": float(np.mean(sc_mae)),
                        "mt_shift": float(np.mean(mt_shift)), "mt_mae": float(np.mean(mt_mae))}
        out[d] = {"base_shift": base_shift, "base_mae": base_mae, "k": per_k}

    for k in K_SWEEP:
        rows = []
        for d in out:
            if k not in out[d]["k"]:
                continue
            b, r = out[d]["base_shift"], out[d]["k"][k]
            red = 100 * (1 - r["sc_shift"] / b) if b > 1e-9 else float("nan")
            rows.append([d, f"{b:.1f}", f"{r['sc_shift']:.1f}", f"{red:+.0f}%",
                         f"{out[d]['base_mae']:.1f}", f"{r['sc_mae']:.1f}",
                         f"{r['mt_shift']:.1f}", f"{r['mt_mae']:.1f}"])
        if not rows:
            continue
        print(f"\n  k = {k} reference cells")
        print(render_table(["fold", "|shift| before", "after (scalar)", "reduction",
                            "MAE before", "MAE after", "|shift| matched", "MAE matched"],
                           rows, aligns=["l"] + ["r"] * 7))

    print("\n  PRE-REGISTERED CRITERION: |shift| reduced >=50% on >=4/6 folds at k <= 5")
    verdict_rows, passed_k = [], None
    for k in K_SWEEP:
        good = 0
        for d in out:
            if k not in out[d]["k"]:
                continue
            b = out[d]["base_shift"]
            if b > 1e-9 and (1 - out[d]["k"][k]["sc_shift"] / b) >= 0.50:
                good += 1
        verdict_rows.append([f"k={k}", f"{good}/{len(out)}",
                             "PASS" if good >= 4 else "fail"])
        if good >= 4 and passed_k is None:
            passed_k = k
    print(render_table(["k", "folds >=50% reduced", "criterion"], verdict_rows,
                       aligns=["l", "r", "l"]))

    # scalar vs matched, to detect a time-varying error
    diffs = [out[d]["k"][k]["sc_mae"] - out[d]["k"][k]["mt_mae"]
             for d in out for k in out[d]["k"]]
    md = float(np.mean(diffs)) if diffs else float("nan")
    print(f"\n   mean(MAE scalar - MAE matched) = {md:+.2f} yr "
          f"({'MATCHED is better -> error is TIME-VARYING' if md > 0.5 else 'scalar is adequate'})")

    print("\n   VERDICT:")
    if passed_k is not None and passed_k <= 5:
        print(f"     => PASS at k={passed_k}. Per-donor calibration is PRACTICAL with a handful of")
        print("        reference cells. IMPLEMENT the fix (MASTER_PLAN.md C2), then snapshot and")
        print("        compare. Use the MATCHED variant if the line above says time-varying.")
    elif passed_k is not None:
        print(f"     => BORDERLINE: only reached at k={passed_k}. Document the cost per donor and")
        print("        decide whether that is experimentally practical before writing code.")
    else:
        print("     => FAIL even at k=10. STOP: do not implement per-donor calibration. Report the")
        print("        model as a WITHIN-DONOR RANKER and state the cross-donor limit openly.")
    return {"passed_k": passed_k, "matched_gain": md}


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTESTS 13-16 — complete pre-change battery (one run)")
    F = {}
    for d in DONORS:
        f = load_fold(d)
        if f is not None:
            F[d] = f
    if not F:
        print("\n   No usable folds found (need cellfate_loocv_* or runs/cellfate_loocv_*).")
        return
    print(f"   loaded {len(F)} folds: {', '.join(F)}")

    r13 = test13(F, render_table)
    r14 = test14(F, render_table)
    r15 = test15(F, render_table)
    r16 = test16(F, render_table)

    print("\n" + "=" * 78)
    print("FINAL SUMMARY — what to do next")
    print("=" * 78)
    if r16.get("passed_k") is not None and r16["passed_k"] <= 5:
        print(f"  1. MAIN FIX IS GO (Test 16 passed at k={r16['passed_k']}). Implement per-donor")
        print("     calibration" + (" using the MATCHED-time variant."
              if r16.get("matched_gain", 0) > 0.5 or r13.get("n_tv", 0) >= 2
              else " (scalar offset is adequate)."))
    elif r16.get("passed_k") is not None:
        print(f"  1. MAIN FIX IS BORDERLINE (needs k={r16['passed_k']}). Decide on practicality.")
    else:
        print("  1. MAIN FIX IS OFF. Report as a within-donor ranker; do not implement C2.")
    if np.isfinite(r14.get("mean_cov", np.nan)) and \
            abs(r14["mean_cov"] - r14["nominal"]) > 0.05:
        print("  2. UNCERTAINTY IS MISCALIBRATED -> an independent upstream cause of RES's")
        print("     collapse. Fix before judging RES further.")
    else:
        print("  2. Uncertainty is calibrated -> ruled out as a cause. ")
    if np.isfinite(r15.get("mean_auc", np.nan)) and r15["mean_auc"] < 0.45:
        print("  3. OOD GATE IS MISLEADING -> disable it in the RES path or fix it.")
    else:
        print("  3. OOD gate behaves acceptably wrt error.")
    print("  4. Then: scorecard snapshot --tag baseline (if not already), make ONE change,")
    print("     snapshot again, and compare. Criteria are in MASTER_PLAN.md §7b.")


if __name__ == "__main__":
    main()
