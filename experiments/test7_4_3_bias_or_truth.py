"""
Test 7.4.3 (ΔAge lab notebook) — is the positive predicted ΔAge model BIAS or biological TRUTH?

Test 7.4.2 found the primary reason raw RES approves nothing: predicted ΔAge (mu_age) is POSITIVE
on every fold (median +4.6 to +21 yr), so R_eff = max(0, -(mu + z*sigma)) = 0, so g = 0, so
RES = 0 before the safety floor is even consulted. Recalibration cannot touch this.

But 7.4.2 could not tell WHY mu_age is positive, and the two possibilities need opposite fixes:

  (BIAS)  the model systematically under-predicts rejuvenation out-of-donor (it must infer each
          donor's control baseline and gets it wrong). Evidence for: N2 has 16/21 truly
          safe-and-rejuvenating cells (true ΔAge < 0) yet ridge predicts +8.76.
          -> Fix the offset and RES revives. Also explains MAE ~14 with ranking ~0.95 (an offset
             inflates MAE but is invisible to rank correlation).

  (TRUTH) the cells genuinely are not rejuvenating by this clock. Evidence for: O1 has only 2/21
          truly safe-and-rejuvenating cells, so a positive prediction may be CORRECT there.
          -> RES is correctly approving nothing; the defect is in the DATA/clock (aging clocks are
             trained on natural aging and extrapolate badly to reprogramming), not in RES.

THREE PARTS:
  Part 1  TRUE vs PREDICTED ΔAge distributions per fold, and the offset (median pred - median
          true). A large positive offset with good rank correlation = BIAS.
  Part 2  ORACLE RES — feed RES the TRUE ΔAge (with the model's sigma, recalibrated fate, real
          OOD flags). This is the decisive part:
            oracle approves a healthy number -> RES's logic is SOUND; ΔAge prediction is the
              defect -> fixing the offset revives the product.
            oracle ALSO approves ~nothing    -> either the cells genuinely do not rejuvenate by
              this clock, or the RES thresholds are mis-tuned for this biology. Either way RES is
              not the thing to fix first.
  Part 3  OFFSET-CORRECTED RES — subtract each fold's offset (estimated on the CALIB split, never
          on test) from the predicted ΔAge and re-run RES. Shows how much of the gap a legitimate,
          leak-free offset correction would actually recover.

NOTE ON PART 3 HONESTY: the offset is estimated on calib (in-distribution cells the model trained
alongside), NOT on the held-out donor. That is the strongest correction available without leaking
test labels. If the true offset is donor-specific, this will under-correct — report, do not tune.

USAGE (repo root, venv active; needs cellfate_loocv_* bundles + calib/test splits):
    python test7_4_3_bias_or_truth.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor, compute_res_batch

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
APPROVED = "APPROVED"


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def ridge_fit_pred(tr, targets):
    """Fit ridge ΔAge once on train; predict for each target split."""
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    m = tr.mask
    reg = Ridge(alpha=1.0).fit(ftr[m], tr.y_age[m])
    out = []
    for s in targets:
        f = np.hstack([sx.transform(s.X), np.asarray(s.fp, float), sdt.transform(s.dose_time)])
        out.append(reg.predict(f))
    return out


def _platt(p_cal, is_pos, p_te) -> np.ndarray:
    p_cal = np.asarray(p_cal, float).reshape(-1, 1)
    is_pos = np.asarray(is_pos, int)
    if not (0 < is_pos.sum() < len(is_pos)):
        return np.asarray(p_te, float)
    lr = LogisticRegression(max_iter=1000).fit(p_cal, is_pos)
    return lr.predict_proba(np.asarray(p_te, float).reshape(-1, 1))[:, 1]


def one_fold(donor: str):
    root = resolve_root(f"cellfate_loocv_{donor}")
    try:
        paths = ArtifactPaths.of(root)
        tr = gather_split(paths, REGIME, "train")
        te = gather_split(paths, REGIME, "test")
        cal = gather_split(paths, REGIME, "calib")
    except Exception:  # noqa: BLE001
        return None
    m, mc = te.mask, cal.mask
    if m.sum() < 3 or mc.sum() < 3:
        return None

    pred = Predictor(root)
    est = ModelEstimator(pred)
    p = pred.res_params

    rows = est.rows(te.X, te.fp, te.dose_time)
    S = np.array([r["S"] for r in rows])
    P_loss = np.array([r["P_loss"] for r in rows])
    sig = np.array([r["sigma_age"] for r in rows])
    ind = np.array([r["in_dist"] for r in rows])
    mu_model = np.array([r["mu_age"] for r in rows])

    cal_rows = est.rows(cal.X, cal.fp, cal.dose_time)
    S_cal = np.array([r["S"] for r in cal_rows])
    S_re = _platt(S_cal, cal.y_cls.astype(int) == SAFE_IDX, S)

    mu_te, mu_cal = ridge_fit_pred(tr, [te, cal])

    y_true = te.y_age[m]
    # offset estimated on CALIB only (leak-free): how much does the predictor overshoot there?
    offset = float(np.median(mu_cal[mc] - cal.y_age[mc]))
    mu_corr = mu_te - offset

    thr = -p.z_conf * float(np.median(sig[m]))   # mu must be below this for any RES credit

    def approvals(mu_vec):
        res, stat = compute_res_batch(S_re, P_loss, mu_vec, sig, ind, p)
        return int((np.asarray(stat)[m] == APPROVED).sum()), res

    # ORACLE: use the TRUE ΔAge where it is valid (cells outside the mask are not counted anyway)
    mu_oracle = np.where(m, te.y_age, mu_te)
    n_oracle, _ = approvals(mu_oracle)
    n_pred, _ = approvals(mu_te)
    n_corr, _ = approvals(mu_corr)

    sp = (float(spearmanr(mu_te[m], y_true).correlation)
          if np.std(mu_te[m]) > 0 and np.std(y_true) > 0 else float("nan"))

    return {
        "n": int(m.sum()),
        "true_med": float(np.median(y_true)), "pred_med": float(np.median(mu_te[m])),
        "true_frac_neg": float((y_true < 0).mean()),
        "pred_frac_neg": float((mu_te[m] < 0).mean()),
        "true_frac_credit": float((y_true < thr).mean()),
        "pred_frac_credit": float((mu_te[m] < thr).mean()),
        "corr_frac_credit": float((mu_corr[m] < thr).mean()),
        "offset": offset, "spearman": sp, "thr": thr,
        "mae_raw": float(np.abs(mu_te[m] - y_true).mean()),
        "mae_corr": float(np.abs(mu_corr[m] - y_true).mean()),
        "n_appr_pred": n_pred, "n_appr_corr": n_corr, "n_appr_oracle": n_oracle,
        "model_med": float(np.median(mu_model[m])),
    }


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 7.4.3 — is the positive predicted ΔAge model BIAS or biological TRUTH?")
    print("'credit' = ΔAge below -z_conf*sigma, the threshold for ANY RES rejuvenation credit.")
    per = {d: r for d in DONORS if (r := one_fold(d)) is not None}
    if not per:
        print("\n   No usable folds.")
        return

    print("\n  PART 1 — TRUE vs PREDICTED ΔAge (is the prediction offset from truth?)")
    rows = [[d, f"{per[d]['true_med']:+.2f}", f"{per[d]['pred_med']:+.2f}",
             f"{per[d]['offset']:+.2f}", f"{per[d]['true_frac_neg']:.2f}",
             f"{per[d]['pred_frac_neg']:.2f}", f"{per[d]['spearman']:+.3f}"] for d in per]
    print(render_table(["fold", "med TRUE", "med PRED", "offset(calib)", "frac true<0",
                        "frac pred<0", "spearman"], rows,
                       aligns=["l", "r", "r", "r", "r", "r", "r"]))

    print("\n  Fraction of cells eligible for ANY RES credit (ΔAge < -z*sigma)")
    rows = [[d, f"{per[d]['thr']:+.2f}", f"{per[d]['true_frac_credit']:.2f}",
             f"{per[d]['pred_frac_credit']:.2f}", f"{per[d]['corr_frac_credit']:.2f}"]
            for d in per]
    print(render_table(["fold", "credit thr", "TRUE eligible", "PRED eligible",
                        "offset-corrected"], rows, aligns=["l", "r", "r", "r", "r"]))

    print("\n  PART 2/3 — APPROVALS: predicted vs offset-corrected vs ORACLE (true ΔAge)")
    rows = [[d, f"{per[d]['n_appr_pred']}/{per[d]['n']}", f"{per[d]['n_appr_corr']}/{per[d]['n']}",
             f"{per[d]['n_appr_oracle']}/{per[d]['n']}", f"{per[d]['mae_raw']:.2f}",
             f"{per[d]['mae_corr']:.2f}"] for d in per]
    print(render_table(["fold", "appr PRED", "appr CORRECTED", "appr ORACLE",
                        "MAE raw", "MAE corrected"], rows,
                       aligns=["l", "r", "r", "r", "r", "r"]))

    tot_o = sum(per[d]["n_appr_oracle"] for d in per)
    tot_p = sum(per[d]["n_appr_pred"] for d in per)
    tot_c = sum(per[d]["n_appr_corr"] for d in per)
    med_off = float(np.median([per[d]["offset"] for d in per]))
    true_elig = float(np.mean([per[d]["true_frac_credit"] for d in per]))

    print(f"\n   totals: approvals PRED={tot_p}  CORRECTED={tot_c}  ORACLE={tot_o}")
    print(f"   median calib offset = {med_off:+.2f} yr   mean TRUE-eligible fraction = "
          f"{true_elig:.2f}")

    print("\n   VERDICT:")
    if tot_o > max(3, 2 * tot_p):
        print("     => BIAS-DOMINATED. With the TRUE ΔAge, RES approves substantially more cells,")
        print("        so RES's logic is SOUND and the bottleneck is ΔAge PREDICTION. The offset")
        print("        (positive mu_age) is a model error, and fixing it revives the product.")
        print("        This also explains MAE ~14 alongside ranking ~0.95 (offset is rank-invisible).")
    elif true_elig < 0.15:
        print("     => TRUTH-DOMINATED. Even the TRUE ΔAge rarely clears the credit threshold:")
        print("        these cells largely do NOT rejuvenate by this clock. RES is CORRECTLY")
        print("        approving nothing. The defect is the DATA/CLOCK (natural-aging clock")
        print("        extrapolated to reprogramming), not the RES formula. Fixing RES is the")
        print("        wrong target; the product claim needs data where rejuvenation occurs.")
    else:
        print("     => MIXED. Oracle approvals are only modestly higher than predicted, and a")
        print("        meaningful fraction of TRUE ΔAge clears the threshold. Both the prediction")
        print("        offset AND genuine lack of rejuvenation contribute. Read the per-fold rows:")
        print("        folds with high 'frac true<0' but low 'frac pred<0' are BIAS; folds where")
        print("        both are low are TRUTH.")
    print("\n   (offset is estimated on CALIB only — never on the held-out donor — so Part 3 is")
    print("    leak-free. If the true offset is donor-specific it will UNDER-correct.)")


if __name__ == "__main__":
    main()
