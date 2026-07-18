"""
Test 7.4.2 (ΔAge lab notebook) — WHY is raw RES collapsed to ~0 out-of-donor?

Test 7.4.1 found raw RES is degenerate: it approves ZERO cells on every held-out donor, while
Platt-recalibrated RES approves a few (mostly correct) cells. From that we DEDUCED the cause:
recalibration changes only S and P_loss (never mu_age/sigma_age), so since recalibrated RES is
non-zero, g(R_eff) must be > 0 and the collapse must come from the SAFETY terms.

This test MEASURES that instead of deducing it, by decomposing RES into its four multiplicative
factors and reporting each one's magnitude, raw vs recalibrated:

    RES = Phi(S) * S**k * g(R_eff) * exp(-lam * P_loss)          [zeroed if not in_dist]

    Phi(S)   = sigmoid((S - tau_safe)/w)   safety floor       <- changes under recalibration
    S**k     = safety-dominant exponent                        <- changes under recalibration
    g(R_eff) = R_eff/(R_eff + kappa),                          <- IDENTICAL in both arms
               R_eff = max(0, -(mu_age + z_conf*sigma_age))       (depends only on age terms)
    exp(-lam*P_loss)                                           <- INERT if lam == 0

Whichever factor is ~0 in the raw arm is the killer. Four candidate diagnoses:
  (1) Phi(S) ~ 0            -> out-of-donor S sits below tau_safe; the narrow safety floor
                               annihilates RES. Fix = recalibrate fate before RES (or widen w).
  (2) g(R_eff) ~ 0          -> sigma_age is so large that the upper age bound is never negative;
                               the model is never CONFIDENTLY rejuvenating. Fix = uncertainty,
                               not calibration (and recalibration could NOT have helped).
  (3) in_dist mostly False  -> the OOD gate is rejecting everything. Fix = the OOD detector.
  (4) exp(-lam*P_loss) ~ 0  -> the cancer-risk penalty dominates (only possible if lam > 0).

Also reports the STATUS breakdown (APPROVED / REJECTED_OOD / REJECTED_UNSAFE /
REJECTED_NO_REJUVENATION), which names the rejection reason directly, and checks whether
lam == 0 (in which case recalibrating P_loss in 7.3/7.4/7.4.1 was a no-op by construction).

Same isolation as 7.4.1: ΔAge held identical (ridge) in both arms; only fate calibration differs.

USAGE (repo root, venv active; needs cellfate_loocv_* bundles + calib/test splits):
    python test7_4_2_res_decompose.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import LOSS_IDX, SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
STATUSES = ["APPROVED", "REJECTED_OOD", "REJECTED_UNSAFE", "REJECTED_NO_REJUVENATION"]


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def ridge_pred(tr, target) -> np.ndarray:
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(target.X), np.asarray(target.fp, float),
                     sdt.transform(target.dose_time)])
    m = tr.mask
    return Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m]).predict(fte)


def _platt(p_cal, is_pos, p_te) -> np.ndarray:
    p_cal = np.asarray(p_cal, float).reshape(-1, 1)
    is_pos = np.asarray(is_pos, int)
    if not (0 < is_pos.sum() < len(is_pos)):
        return np.asarray(p_te, float)
    lr = LogisticRegression(max_iter=1000).fit(p_cal, is_pos)
    return lr.predict_proba(np.asarray(p_te, float).reshape(-1, 1))[:, 1]


def _sigmoid(x):
    x = np.clip(np.asarray(x, dtype=np.float64), -500.0, 500.0)
    return 1.0 / (1.0 + np.exp(-x))


def factors(S, P_loss, mu, sig, p):
    """The four multiplicative RES terms, computed exactly as res.py does."""
    R_eff = np.maximum(0.0, -(mu + p.z_conf * sig))
    return {
        "phi": _sigmoid((S - p.tau_safe) / p.w),
        "s_k": np.asarray(S, float) ** p.k,
        "g": R_eff / (R_eff + p.kappa),
        "ploss": np.exp(-p.lam * np.asarray(P_loss, float)),
        "R_eff": R_eff,
    }


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        cal = gather_split(paths, REGIME, "calib")
    except Exception:  # noqa: BLE001
        return None
    m = te.mask
    if m.sum() < 3:
        return None

    pred = Predictor(root)
    est = ModelEstimator(pred)
    p = pred.res_params

    rows = est.rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    sig = np.array([r["sigma_age"] for r in rows])
    ind = np.array([r["in_dist"] for r in rows])

    cal_rows = est.rows(cal.X, cal.fp, cal.dose_time)
    S_cal = np.array([r["S"] for r in cal_rows])
    Pl_cal = np.array([r["P_loss"] for r in cal_rows])
    cls_cal = cal.y_cls.astype(int)
    S_re = _platt(S_cal, cls_cal == SAFE_IDX, S)
    Pl_re = _platt(Pl_cal, cls_cal == LOSS_IDX, P_loss)

    mu = ridge_pred(tr, te)

    fB = factors(S, P_loss, mu, sig, p)
    fC = factors(S_re, Pl_re, mu, sig, p)
    resB, statB = compute_res_batch(S, P_loss, mu, sig, ind, p)
    resC, statC = compute_res_batch(S_re, Pl_re, mu, sig, ind, p)

    unsafe_thr = p.tau_safe - 3.0 * p.w
    med = lambda a: float(np.median(np.asarray(a)[m]))  # noqa: E731

    out = {
        "params": p, "n": int(m.sum()), "n_ood": int((~ind[m]).sum()),
        "S_raw_med": med(S), "S_re_med": med(S_re),
        "S_raw_max": float(np.max(S[m])), "S_re_max": float(np.max(S_re[m])),
        "unsafe_thr": unsafe_thr,
        "n_below_unsafe_raw": int((S[m] < unsafe_thr).sum()),
        "n_below_unsafe_re": int((S_re[m] < unsafe_thr).sum()),
        "sigma_med": med(sig), "mu_med": med(mu), "R_eff_med": med(fB["R_eff"]),
        "resB_med": med(resB), "resC_med": med(resC),
        "resB_max": float(np.max(resB[m])), "resC_max": float(np.max(resC[m])),
    }
    for arm, f in (("B", fB), ("C", fC)):
        for key in ("phi", "s_k", "g", "ploss"):
            out[f"{arm}_{key}"] = med(f[key])
    for arm, st in (("B", statB), ("C", statC)):
        for s in STATUSES:
            out[f"{arm}_{s}"] = int((np.asarray(st)[m] == s).sum())
    return out


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7.4.2 — WHY is raw RES ~0 out-of-donor? (decompose the four RES factors)")
    print("RES = Phi(S) * S^k * g(R_eff) * exp(-lam*P_loss), zeroed if not in_dist.")
    print("Phi and S^k change under recalibration; g is IDENTICAL in both arms (age terms only).")

    per = {d: r for d in DONORS if (r := one_fold(d)) is not None}
    if not per:
        print("\n   No usable folds (need cellfate_loocv_* with calib+test splits).")
        return

    p = per[next(iter(per))]["params"]
    print(f"\n  RES PARAMS: tau_safe={p.tau_safe}  w={p.w}  k={p.k}  kappa={p.kappa}  "
          f"z_conf={p.z_conf}  lam={p.lam}")
    print(f"  => REJECTED_UNSAFE fires when S < tau_safe - 3w = {p.tau_safe - 3*p.w:.3f}")
    if p.lam == 0.0:
        print("  !! lam = 0 -> exp(-lam*P_loss) = 1 ALWAYS. The P_loss term is INERT, so")
        print("     recalibrating P_loss in 7.3/7.4/7.4.1 was a NO-OP by construction.")

    print("\n  FACTOR MEDIANS (the ~0 factor is the killer; g is the same in both arms)")
    rows = [[d, f"{per[d]['B_phi']:.2e}", f"{per[d]['C_phi']:.2e}",
             f"{per[d]['B_s_k']:.3f}", f"{per[d]['C_s_k']:.3f}",
             f"{per[d]['B_g']:.3f}", f"{per[d]['B_ploss']:.3f}",
             f"{per[d]['resB_med']:.2e}", f"{per[d]['resC_med']:.2e}"] for d in per]
    print(render_table(["fold", "Phi B", "Phi C", "S^k B", "S^k C", "g (both)", "Ploss",
                        "RES B", "RES C"], rows,
                       aligns=["l", "r", "r", "r", "r", "r", "r", "r", "r"]))

    print("\n  SAFETY GATE: where does S sit vs the thresholds?")
    rows = [[d, f"{per[d]['S_raw_med']:.3f}", f"{per[d]['S_re_med']:.3f}",
             f"{per[d]['S_raw_max']:.3f}", f"{per[d]['S_re_max']:.3f}",
             f"{per[d]['n_below_unsafe_raw']}/{per[d]['n']}",
             f"{per[d]['n_below_unsafe_re']}/{per[d]['n']}"] for d in per]
    print(render_table(["fold", "med S raw", "med S recal", "max S raw", "max S recal",
                        "n unsafe raw", "n unsafe recal"], rows,
                       aligns=["l", "r", "r", "r", "r", "r", "r"]))

    print("\n  REJUVENATION TERM: is the model ever CONFIDENTLY rejuvenating?")
    rows = [[d, f"{per[d]['mu_med']:+.2f}", f"{per[d]['sigma_med']:.2f}",
             f"{per[d]['R_eff_med']:.2f}", f"{per[d]['B_g']:.3f}",
             f"{per[d]['n_ood']}/{per[d]['n']}"] for d in per]
    print(render_table(["fold", "med mu_age", "med sigma", "med R_eff", "med g", "n OOD"],
                       rows, aligns=["l", "r", "r", "r", "r", "r"]))

    print("\n  STATUS BREAKDOWN (the rejection reason, named directly)")
    rows = []
    for d in per:
        rows.append([d] + [f"{per[d][f'B_{s}']}" for s in STATUSES]
                    + [f"{per[d][f'C_{s}']}" for s in STATUSES])
    print(render_table(["fold", "B appr", "B ood", "B unsafe", "B no-rejuv",
                        "C appr", "C ood", "C unsafe", "C no-rejuv"], rows,
                       aligns=["l"] + ["r"] * 8))

    # ---- verdict ----
    mphi = float(np.median([per[d]["B_phi"] for d in per]))
    mg = float(np.median([per[d]["B_g"] for d in per]))
    mood = sum(per[d]["n_ood"] for d in per)
    mploss = float(np.median([per[d]["B_ploss"] for d in per]))
    print("\n   VERDICT — which factor kills raw RES?")
    killers = []
    if mphi < 0.01:
        killers.append(("Phi(S) safety floor", mphi))
    if mg < 0.01:
        killers.append(("g(R_eff) rejuvenation term", mg))
    if mploss < 0.01:
        killers.append(("exp(-lam*P_loss) risk penalty", mploss))
    if mood > 0:
        print(f"     - OOD gate rejects {mood} cell(s) across folds "
              f"({'a factor' if mood else 'not a factor'}).")
    if not killers:
        print("     - No single factor is ~0: the collapse is the PRODUCT of several small terms.")
    for name, val in killers:
        print(f"     - **{name}** has median {val:.2e} -> this is a primary killer.")
    if mphi < 0.01 <= mg:
        print("\n     => MECHANISM CONFIRMED (diagnosis 1): out-of-donor S sits below tau_safe, so")
        print("        the narrow safety floor annihilates RES while the rejuvenation term is")
        print("        healthy. Recalibration raises S into the usable range -> RES functions.")
        print("        This is a CALIBRATION defect, not an uncertainty or OOD defect.")
    elif mg < 0.01:
        print("\n     => DIAGNOSIS 2: g(R_eff) ~ 0 -> sigma_age is too large for the upper age bound")
        print("        ever to be negative; the model is never CONFIDENTLY rejuvenating. This is an")
        print("        UNCERTAINTY defect; recalibrating fate cannot fix it.")
    print("\n   (Compare 'RES B' vs 'RES C' in the first table: the size of the gap is how much")
    print("    of the collapse recalibration actually undoes.)")


if __name__ == "__main__":
    main()
