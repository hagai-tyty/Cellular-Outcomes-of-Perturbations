"""STAGE 3a — DIAGNOSIS of the withdrawn STOP, and the resolvability check that was never run.

    python experiments/stage3a_diagnose.py            # operative arm (_c7)
    python experiments/stage3a_diagnose.py _armA      # contaminated arm, for contrast

READ-ONLY. Writes `results/stage3a_diagnose_results.json`. `src/` untouched, no build touched,
no retrain. `experiments/test18_forward_gate.py` is imported UNMODIFIED and its own
`timepoint_table` / `build_pairs` / `feats` / `paired_ci` are the primitives used here, so what
is diagnosed is the estimator 3a actually ran -- not a re-implementation of it.

WHY THIS EXISTS
---------------
`STAGE_3_TOOL.md` §3a returned STOP, which is the one TERMINAL verdict in this project ("do not
write tool code; ship the scoring model; go to Stage 5"). That verdict was WITHDRAWN on
2026-08-08 because it is produced by a fit that diverges, not by an absence of signal. The
decisive fact: Part C's target is `unsafe.mean()`, a FRACTION bounded `[0,1]` by construction
(`test18_forward_gate.py:85`), yet held-out Y1 reports MAE **7.589** -- the model is emitting
values around 8 on a quantity that cannot exceed 1.

Four items were accepted when the STOP was withdrawn. This script is all four, plus one section
(D0) added on 2026-08-12 after the shape of the target became visible -- it is NOT pre-registered
and is labelled as such wherever it appears:

  D0  What Part C's target actually is (the mean of ~1.7 binary cells per timepoint, ~90 %
      saturated at 0 or 1, nearly a monotone step in time), and the MODEL-FREE CEILING on the
      forward question: predict the held-out donor's value at t_j from the OTHER donors at the
      same t_j. That uses only t_j = t_i + Δt -- no genes, no fitting, no leakage.
  D1  Per fold, the held-out Δt range against the training range, plus the fitted Δt
      coefficient -- and a decomposition of each prediction into its GENE block and its Δt
      block, so the divergence is localised rather than described.
  D2  Bound Part C's prediction to [0,1]. The target is a fraction by construction, so a
      predictor that leaves the range is misspecified. `clip` is the minimum; a logit link is
      the honest fix. Both are measured, against the raw estimator and against a mean-only
      baseline.
  D3  Part B's identical −269.13 across all five folds, explained ARITHMETICALLY and then
      re-run correctly (per-fold LODO fits, each fold swept over its own Δt range).
  D4  `bar_verdict` at the real geometry -- 5 folds, 55-66 pairs, 11-12 timepoints -- with a
      SIMULATED TRUE Δt effect. `REF_GROUND_RULES.md` §5b requires this BEFORE a bar is graded
      and it was never run for 3a. If a correct system cannot clear the bar at this scale then
      the dataset cannot answer the question, and THAT is the finding, not STOP.

PRE-REGISTERED EXPECTATIONS (written before the script was run; graded in §D5 of the output)
-------------------------------------------------------------------------------------------
  P1  Y1's Δt range is NESTED INSIDE the training range, not outside it. Y1 is missing d29, so
      it has FEWER pairs spanning a NARROWER span -- the opposite of extrapolation in Δt. If
      P1 holds, the "Y1's Δt sits outside the training support" story is WRONG and the
      divergence must be located somewhere else. Stated so it can fail.
  P2  The divergence is carried by the Δt BLOCK, not the gene block: the mean absolute Δt-block
      contribution for Y1 exceeds every other fold's by more than an order of magnitude.
  P3  Part B's swing is IDENTICAL across folds to machine precision, and equals the analytic
      value w_dt·Δz(dt) + w_dt²·Δz(dt²), which contains no x0 term at all. SD across folds = 0.
  P4  Part C's "state only" arm is ALSO broken: its MAE exceeds the mean-only baseline. If so,
      3a had no working arm to compare against, and "state+Δt ties state" was never a
      comparison between two functioning predictors.
  P5  At this geometry the raw estimator FAILS `bar_verdict` even when the target is a pure
      function of Δt with near-zero noise -> the bar is UNRESOLVABLE and STOP was measuring the
      instrument.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
It does not re-take the 3a verdict. It establishes whether 3a's instrument can take one at all.
A repaired 3a run is a separate, pre-registered act.
"""
from __future__ import annotations

import contextlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "stage3a_diagnose_results.json"

ALPHA = 1.0             # test18's ridge alpha, unchanged
LOGIT_EPS = 0.01        # clamp before logit; the target is a fraction over <=2 cells/timepoint
SIM_TRIALS = 2000       # D4 trials per (rho, sigma) cell
SIM_SEED = 0
RHOS = (0.0, 0.25, 0.50, 0.75, 1.0)      # share of the simulated signal carried by Δt
SIGMAS = (0.05, 0.001)                    # observation noise on the simulated fraction


# ---------------------------------------------------------------------------------------------
# arm loading -- test18's own primitives, redirected exactly as stage3a_forward_gate.py does
# ---------------------------------------------------------------------------------------------
def load_t18(suffix: str):
    if "test18_forward_gate" in sys.modules:
        del sys.modules["test18_forward_gate"]
    t18 = importlib.import_module("test18_forward_gate")
    t18.resolve_root = lambda name, _s=suffix: str(REPO / f"{name}{_s}")
    return t18


def rows_ignoring_age_mask(t18, donor: str):
    """Part C's target needs CLASS labels only -- rebuild rows without the AGE mask.

    test18's `timepoint_table` filters on `te.mask` (age validity) before computing the unsafe
    fraction. Under C-7 rule 4 donor N2's ΔAge is masked in every fold, so N2 is dropped from
    the run entirely -- including from Part C, whose target does not depend on ΔAge at all.
    That costs a whole fold of SAFETY data for an AGE reason. This rebuild keeps it.
    """
    from cellfate.common.constants import DEATH_IDX, LOSS_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split

    try:
        te = gather_split(ArtifactPaths.of(t18.resolve_root(f"cellfate_loocv_{donor}")),
                          t18.REGIME, "test")
    except Exception:                                                   # noqa: BLE001
        return None
    if te.n < 4:
        return None
    t = np.asarray(te.dose_time[:, 1], float)
    X = np.asarray(te.X, float)
    cls = te.y_cls.astype(int)
    unsafe = ((cls == LOSS_IDX) | (cls == DEATH_IDX)).astype(float)
    tps = np.unique(np.round(t, 6))
    if len(tps) < 3:
        return None
    rows = []
    for tp in tps:
        sel = np.isclose(t, tp)
        if sel.sum() < t18.MIN_CELLS_PER_TP:
            continue
        rows.append({"t": float(tp), "x": X[sel].mean(0), "y": float("nan"),
                     "u": float(unsafe[sel].mean()), "n": int(sel.sum())})
    return rows if len(rows) >= 3 else None


def build_per(t18, *, ignore_age_mask: bool) -> dict:
    per = {}
    for d in t18.DONORS:
        rows = rows_ignoring_age_mask(t18, d) if ignore_age_mask else t18.timepoint_table(d)
        if rows is None:
            continue
        per[d] = {"rows": rows, "pairs": t18.build_pairs(rows)}
    return per


def lodo_folds(per: dict) -> list[str]:
    """The folds test18 actually grades: >=8 training pairs and >=3 held-out pairs."""
    out = []
    for hd in per:
        tr = sum(len(per[d]["pairs"]) for d in per if d != hd)
        if tr >= 8 and len(per[hd]["pairs"]) >= 3:
            out.append(hd)
    return out


# ---------------------------------------------------------------------------------------------
# a ridge whose PREDICTIONS are precomputable in y -- validated against sklearn before use
# ---------------------------------------------------------------------------------------------
class FrozenRidge:
    """Ridge with fixed features: pred = H @ (y_tr - ȳ_tr) + ȳ_tr, H independent of y.

    Ridge with `fit_intercept=True` centres X and y on the training mean, so with the design
    matrix held fixed the whole map from `y_tr` to `pred_te` is the single matrix
    `H = Zc_te · Zc_tr' · (Zc_tr · Zc_tr' + αI)⁻¹`. D4 refits the same folds thousands of times
    with only the TARGET changing, so H is built once per (fold, feature set).

    The dual form is used because n_pairs (~264) << n_features (~5k). `max_abs_err_vs_sklearn`
    records the agreement with `sklearn.linear_model.Ridge` on the REAL target, so the fast path
    is audited rather than assumed.
    """

    def __init__(self, Xtr: np.ndarray, Xte: np.ndarray, alpha: float = ALPHA) -> None:
        sc = StandardScaler().fit(Xtr)
        Ztr, Zte = sc.transform(Xtr), sc.transform(Xte)
        self.mu = Ztr.mean(0)
        Zc_tr, Zc_te = Ztr - self.mu, Zte - self.mu
        K = Zc_tr @ Zc_tr.T + alpha * np.eye(len(Zc_tr))
        self.H = Zc_te @ (Zc_tr.T @ np.linalg.inv(K))
        self.scaler, self.Ztr, self.Zte, self.Zc_te, self.alpha = sc, Ztr, Zte, Zc_te, alpha

    def predict(self, ytr: np.ndarray) -> np.ndarray:
        ybar = float(np.mean(ytr))
        return self.H @ (np.asarray(ytr, float) - ybar) + ybar

    def sklearn_predict(self, ytr: np.ndarray) -> np.ndarray:
        m = Ridge(alpha=self.alpha).fit(self.Ztr, ytr)
        return m.predict(self.Zte)

    def coef(self, ytr: np.ndarray) -> np.ndarray:
        return Ridge(alpha=self.alpha).fit(self.Ztr, ytr).coef_


def logit(p: np.ndarray, eps: float = LOGIT_EPS) -> np.ndarray:
    q = np.clip(np.asarray(p, float), eps, 1.0 - eps)
    return np.log(q / (1.0 - q))


def expit(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(z, float), -60.0, 60.0)))


def mae(a, b) -> float:
    return float(np.abs(np.asarray(a, float) - np.asarray(b, float)).mean())


def oracle_on_tj(per_map: dict, fold_list: list[str], key: str, paired_ci) -> dict:
    """The model-free ceiling: predict the held-out donor at t_j from the OTHERS at t_j.

    Uses only `t_j = t_i + Δt` -- no genes, no fitting, no held-out information -- so it is the
    most forward-in-time information Δt can carry at this geometry. Compared against the pooled
    training mean, which is what a predictor that ignores time would give. If the oracle does not
    beat that, the forward question has no answer in this data; if it does, any estimator that
    reads "tied" is failing to find a signal that is present.
    """
    rows_o, o_mae, m_mae, kept = [], [], [], {}
    for hd in fold_list:
        tr_d = [d for d in per_map if d != hd]
        R = per_map[hd]["rows"]
        tr_vals = [np.array([r[key] for r in per_map[d]["rows"]], float) for d in tr_d]
        tr_ts = [np.array([r["t"] for r in per_map[d]["rows"]], float) for d in tr_d]
        preds, truth = [], []
        for i in range(len(R)):
            for j in range(len(R)):
                if R[j]["t"] <= R[i]["t"]:
                    continue
                tj = R[j]["t"]
                preds.append(float(np.mean([v[int(np.argmin(np.abs(ts - tj)))]
                                            for v, ts in zip(tr_vals, tr_ts, strict=True)])))
                truth.append(float(R[j][key]))
        if not preds or not np.isfinite(truth).all():
            continue
        pooled = float(np.mean(np.concatenate(tr_vals)))
        o_mae.append(mae(preds, truth))
        m_mae.append(mae(np.full(len(truth), pooled), truth))
        kept[hd] = preds
        rows_o.append([hd, str(len(truth)), f"{m_mae[-1]:.3f}", f"{o_mae[-1]:.3f}",
                       f"{m_mae[-1] - o_mae[-1]:+.3f}"])
    md, (lo, hi), n = paired_ci([b - a for a, b in zip(m_mae, o_mae, strict=True)])
    return {"rows": rows_o, "preds": kept,
            "oracle_mean": float(np.mean(o_mae)) if o_mae else float("nan"),
            "pooled_mean": float(np.mean(m_mae)) if m_mae else float("nan"),
            "paired_mean": md, "ci": [lo, hi], "n": n,
            "verdict": ("t_j HELPS" if hi < 0 else "t_j HURTS" if lo > 0 else "tied")}


# ---------------------------------------------------------------------------------------------
def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    suffix = sys.argv[1] if len(sys.argv) > 1 else "_c7"
    t18 = load_t18(suffix)
    missing = [d for d in t18.DONORS if not (REPO / f"cellfate_loocv_{d}{suffix}").exists()]
    if missing:
        print(f"missing fold roots for arm {suffix}: {missing}")
        return 1

    print("\n" + "=" * 90)
    print(f"STAGE 3a — DIAGNOSIS of the withdrawn STOP     arm: {suffix}")
    print("=" * 90)
    print("READ-ONLY. Does NOT re-take the verdict; establishes whether 3a's instrument can")
    print("take one. Pre-registered expectations P1-P5 are in the docstring and graded in D5.")

    out: dict = {"script": "stage3a_diagnose", "arm": suffix, "alpha": ALPHA}

    per = build_per(t18, ignore_age_mask=False)      # exactly what 3a graded
    per_all = build_per(t18, ignore_age_mask=True)   # Part C without the AGE mask
    folds = lodo_folds(per)
    folds_all = lodo_folds(per_all)
    G = per[folds[0]]["rows"][0]["x"].shape[0]
    out["geometry"] = {
        "genes": int(G),
        "folds_graded_by_3a": folds,
        "folds_available_for_safety_without_age_mask": folds_all,
        "pairs": {d: len(per[d]["pairs"]) for d in per},
        "timepoints": {d: len(per[d]["rows"]) for d in per},
        "pairs_no_age_mask": {d: len(per_all[d]["pairs"]) for d in per_all},
        "timepoints_no_age_mask": {d: len(per_all[d]["rows"]) for d in per_all},
    }
    print(f"\n  genes={G}  folds 3a graded={folds}  folds available without the age mask="
          f"{folds_all}")

    # =========================================================================================
    # D0 — what the Part C target actually IS, and the ceiling Δt could reach at best
    # =========================================================================================
    # NOT pre-registered. Added 2026-08-12 after the target's shape became visible: it is the
    # mean of ~1.7 binary cells per timepoint, ~90% saturated at 0 or 1, and per donor it is
    # very close to a monotone STEP in time. If that is what it is, then "can Δt predict it
    # forward" has a model-free ceiling: predict the held-out donor's u at t_j from the OTHER
    # donors' u at the same timepoint. That oracle uses ONLY t_j = t_i + Δt -- no genes, no
    # fitting, no leakage -- so it is the most forward information Δt can carry at this
    # geometry. If the oracle beats state-only, a forward signal EXISTS and 3a's estimator
    # merely failed to find it. If it does not, the data genuinely lacks one.
    print("\n" + "-" * 90)
    print("D0 — WHAT THE PART C TARGET IS, AND THE MODEL-FREE CEILING ON Δt (not pre-registered)")
    print("-" * 90)

    prof_rows, n_sat, n_tot = [], 0, 0
    for d in per_all:
        u = np.array([r["u"] for r in per_all[d]["rows"]], float)
        cells = sum(r["n"] for r in per_all[d]["rows"])
        sat = int(((u == 0.0) | (u == 1.0)).sum())
        n_sat, n_tot = n_sat + sat, n_tot + len(u)
        mono = bool(np.all(np.diff(u) >= -1e-12))
        prof_rows.append([d, str(len(u)), str(cells), f"{cells / len(u):.1f}",
                          " ".join(f"{v:g}" for v in u), f"{sat}/{len(u)}",
                          f"{u.std():.3f}", "yes" if mono else "no"])
    print(render_table(
        ["donor", "timepoints", "cells", "cells/tp", "unsafe fraction by timepoint",
         "at 0 or 1", "std", "monotone?"], prof_rows,
        aligns=["l", "r", "r", "r", "l", "r", "r", "l"]))
    print(f"   {n_sat}/{n_tot} = {100 * n_sat / n_tot:.0f}% of target values sit exactly at 0 or 1,")
    print("   each estimated from 1-2 binary cells. Least-squares ridge on that is misspecified")
    print("   twice over: unbounded output, and a target that is nearly a binary step in time.")

    out["D0"] = {"saturated_at_bounds": [n_sat, n_tot],
                 "profiles": {d: [float(r["u"]) for r in per_all[d]["rows"]] for d in per_all},
                 "cells_per_timepoint": {d: sum(r["n"] for r in per_all[d]["rows"])
                                         / len(per_all[d]["rows"]) for d in per_all}}
    for lbl, pm, fl, key in (
            ("UNSAFE FRACTION, as 3a graded it (age mask on)", per, folds, "u"),
            ("UNSAFE FRACTION, without the age mask", per_all, folds_all, "u"),
            ("ΔAge (Part A's target)", per, folds, "y")):
        r = oracle_on_tj(pm, fl, key, t18.paired_ci)
        print(f"\n  ORACLE ON t_j ALONE — {lbl}")
        print(render_table(["held-out", "pairs", "pooled mean", "oracle on t_j", "gain"],
                           r["rows"], aligns=["l", "r", "r", "r", "r"]))
        print(f"   pooled-mean baseline={r['pooled_mean']:.3f}  oracle={r['oracle_mean']:.3f}  "
              f"paired={r['paired_mean']:+.3f} CI=[{r['ci'][0]:+.3f},{r['ci'][1]:+.3f}] "
              f"(n={r['n']}) -> {r['verdict']}")
        out["D0"].setdefault("oracle", {})[lbl] = r

    # =========================================================================================
    # D1 — where the divergence lives
    # =========================================================================================
    print("\n" + "-" * 90)
    print("D1 — Δt SUPPORT, FITTED Δt COEFFICIENTS, AND WHICH BLOCK CARRIES THE DIVERGENCE")
    print("-" * 90)

    d1_rows, d1 = [], {}
    repro = {"partA": {}, "partC": {}}
    fast_err = 0.0
    for hd in folds:
        tr = [p for d in per if d != hd for p in per[d]["pairs"]]
        te = per[hd]["pairs"]
        dt_tr = np.array([p["dt"] for p in tr], float)
        dt_te = np.array([p["dt"] for p in te], float)
        outside = int(((dt_te < dt_tr.min()) | (dt_te > dt_tr.max())).sum())

        Xtr, Xte = t18.feats(tr, True), t18.feats(te, True)
        fr = FrozenRidge(Xtr, Xte)
        ytr = np.array([p["y_j"] for p in tr], float)
        yte = np.array([p["y_j"] for p in te], float)
        utr = np.array([p["u_j"] for p in tr], float)
        ute = np.array([p["u_j"] for p in te], float)

        # the fast path is only usable if it reproduces sklearn on the real target
        fast_err = max(fast_err, float(np.abs(fr.predict(ytr) - fr.sklearn_predict(ytr)).max()))
        fast_err = max(fast_err, float(np.abs(fr.predict(utr) - fr.sklearn_predict(utr)).max()))

        # reproduce 3a's own two numbers for this fold, so the diagnosis is anchored
        Xtr_s, Xte_s = t18.feats(tr, False), t18.feats(te, False)
        fr_s = FrozenRidge(Xtr_s, Xte_s)
        repro["partA"][hd] = {"state": mae(fr_s.predict(ytr), yte),
                              "state_dt": mae(fr.predict(ytr), yte)}
        repro["partC"][hd] = {"state": mae(fr_s.predict(utr), ute),
                              "state_dt": mae(fr.predict(utr), ute)}

        # decomposition of the Part C prediction: gene block vs Δt block
        w = fr.coef(utr)
        gene_c = fr.Zc_te[:, :G] @ w[:G]
        dt_c = fr.Zc_te[:, G:] @ w[G:]
        pred_u = fr.predict(utr)

        # NOT pre-registered — added after P2 failed, to locate the divergence properly.
        # The same gene block WITHOUT the Δt columns, so the two can be compared directly,
        # plus how far the held-out STATE sits outside the training scaler's support.
        w_s = fr_s.coef(utr)
        gene_c_s = fr_s.Zc_te @ w_s
        cos = float(w[:G] @ w_s / max(np.linalg.norm(w[:G]) * np.linalg.norm(w_s), 1e-30))
        zg_te, zg_tr = fr.Zc_te[:, :G], fr.Ztr[:, :G] - fr.mu[:G]

        # the same decomposition for Part A, whose Y1 MAE is 311.47
        wA = fr.coef(ytr)
        gene_cA = fr.Zc_te[:, :G] @ wA[:G]
        dt_cA = fr.Zc_te[:, G:] @ wA[G:]

        z_dt = fr.Zc_te[:, G]
        d1[hd] = {
            "dt_train_range": [float(dt_tr.min()), float(dt_tr.max())],
            "dt_heldout_range": [float(dt_te.min()), float(dt_te.max())],
            "heldout_pairs_outside_train_dt_support": outside,
            "z_dt_heldout_absmax": float(np.abs(z_dt).max()),
            "z_dt_train_absmax": float(np.abs(fr.Ztr[:, G] - fr.mu[G]).max()),
            "coef_dt_partC": float(w[G]), "coef_dt2_partC": float(w[G + 1]),
            "coef_dt_partA": float(wA[G]), "coef_dt2_partA": float(wA[G + 1]),
            "gene_block_norm_partC": float(np.linalg.norm(w[:G])),
            "dt_block_norm_partC": float(np.linalg.norm(w[G:])),
            "mean_abs_gene_contrib_partC": float(np.abs(gene_c).mean()),
            "mean_abs_dt_contrib_partC": float(np.abs(dt_c).mean()),
            "mean_abs_gene_contrib_partA": float(np.abs(gene_cA).mean()),
            "mean_abs_dt_contrib_partA": float(np.abs(dt_cA).mean()),
            "partC_pred_range": [float(pred_u.min()), float(pred_u.max())],
            "partC_frac_pred_outside_unit": float(((pred_u < 0) | (pred_u > 1)).mean()),
            # not pre-registered; added after P2 failed
            "mean_abs_gene_contrib_partC_state_only": float(np.abs(gene_c_s).mean()),
            "gene_block_norm_partC_state_only": float(np.linalg.norm(w_s)),
            "cos_gene_weights_state_vs_state_dt": cos,
            "z_gene_heldout_absmax": float(np.abs(zg_te).max()),
            "z_gene_train_absmax": float(np.abs(zg_tr).max()),
            "z_gene_heldout_mean_abs": float(np.abs(zg_te).mean()),
            "z_gene_train_mean_abs": float(np.abs(zg_tr).mean()),
        }
        d1_rows.append([
            hd, f"[{dt_tr.min():.2f},{dt_tr.max():.2f}]", f"[{dt_te.min():.2f},{dt_te.max():.2f}]",
            str(outside), f"{np.abs(z_dt).max():.2f}",
            f"{w[G]:+.3f}", f"{np.abs(gene_c).mean():.3f}", f"{np.abs(dt_c).mean():.3f}",
            f"[{pred_u.min():.2f},{pred_u.max():.2f}]",
        ])

    print(f"\n   fast-path agreement with sklearn.Ridge on the real targets: "
          f"max|Δ| = {fast_err:.3e}")
    out["frozen_ridge_max_abs_err_vs_sklearn"] = fast_err

    print("\n  Δt SUPPORT AND PART-C PREDICTION DECOMPOSITION (target is a fraction in [0,1])")
    print(render_table(
        ["held-out", "train Δt", "held-out Δt", "outside", "max|z_Δt|", "coef Δt",
         "|gene contrib|", "|Δt contrib|", "pred range"],
        d1_rows, aligns=["l", "r", "r", "r", "r", "r", "r", "r", "r"]))
    out["D1"] = d1

    print("\n  D1b — WHICH BLOCK ACTUALLY DIVERGES (added after P2 failed; NOT pre-registered)")
    print("   the same gene block with and without the Δt columns, and how far the held-out")
    print("   STATE sits outside the training scaler's support.")
    print(render_table(
        ["held-out", "|gene| state-only", "|gene| state+Δt", "ratio", "cos(w_gene)",
         "mean|z_gene| tr", "mean|z_gene| te", "pred outside [0,1]"],
        [[hd, f"{d1[hd]['mean_abs_gene_contrib_partC_state_only']:.3f}",
          f"{d1[hd]['mean_abs_gene_contrib_partC']:.3f}",
          f"{d1[hd]['mean_abs_gene_contrib_partC'] / max(d1[hd]['mean_abs_gene_contrib_partC_state_only'], 1e-12):.1f}x",
          f"{d1[hd]['cos_gene_weights_state_vs_state_dt']:+.3f}",
          f"{d1[hd]['z_gene_train_mean_abs']:.3f}", f"{d1[hd]['z_gene_heldout_mean_abs']:.3f}",
          f"{100 * d1[hd]['partC_frac_pred_outside_unit']:.0f}%"] for hd in folds],
        aligns=["l", "r", "r", "r", "r", "r", "r", "r"]))

    print("\n  3a REPRODUCTION CHECK — these must match the recorded run")
    print(render_table(
        ["held-out", "A state", "A state+Δt", "C state", "C state+Δt"],
        [[hd, f"{repro['partA'][hd]['state']:.2f}", f"{repro['partA'][hd]['state_dt']:.2f}",
          f"{repro['partC'][hd]['state']:.3f}", f"{repro['partC'][hd]['state_dt']:.3f}"]
         for hd in folds], aligns=["l", "r", "r", "r", "r"]))
    out["reproduction"] = repro

    # =========================================================================================
    # D2 — bound the prediction; and the baseline 3a never printed
    # =========================================================================================
    print("\n" + "-" * 90)
    print("D2 — PART C WITH A BOUNDED PREDICTOR, AND AGAINST A MEAN-ONLY BASELINE")
    print("-" * 90)
    print("   the target is `unsafe.mean()`, a FRACTION. raw ridge is unbounded; clip is the")
    print("   minimum fix; a logit link is the honest one. `mean only` predicts the training")
    print("   mean and is the floor any usable predictor must beat.")

    def part_c(per_map: dict, fold_list: list[str]) -> dict:
        acc = {k: [] for k in ("mean_only", "state_raw", "dt_raw", "state_clip", "dt_clip",
                               "state_logit", "dt_logit")}
        rows = []
        for hd in fold_list:
            tr = [p for d in per_map if d != hd for p in per_map[d]["pairs"]]
            te = per_map[hd]["pairs"]
            utr = np.array([p["u_j"] for p in tr], float)
            ute = np.array([p["u_j"] for p in te], float)
            if np.std(utr) == 0 or np.std(ute) == 0:
                continue
            fr_s = FrozenRidge(t18.feats(tr, False), t18.feats(te, False))
            fr_d = FrozenRidge(t18.feats(tr, True), t18.feats(te, True))
            vals = {
                "mean_only": mae(np.full(len(ute), utr.mean()), ute),
                "state_raw": mae(fr_s.predict(utr), ute),
                "dt_raw": mae(fr_d.predict(utr), ute),
                "state_clip": mae(np.clip(fr_s.predict(utr), 0, 1), ute),
                "dt_clip": mae(np.clip(fr_d.predict(utr), 0, 1), ute),
                "state_logit": mae(expit(fr_s.predict(logit(utr))), ute),
                "dt_logit": mae(expit(fr_d.predict(logit(utr))), ute),
            }
            for k, v in vals.items():
                acc[k].append(v)
            rows.append([hd, str(len(te))] + [f"{vals[k]:.3f}" for k in
                         ("mean_only", "state_raw", "dt_raw", "state_clip", "dt_clip",
                          "state_logit", "dt_logit")])
        res = {"per_fold": rows, "means": {k: float(np.mean(v)) for k, v in acc.items() if v}}
        for lbl, s, d in (("raw", "state_raw", "dt_raw"), ("clip", "state_clip", "dt_clip"),
                          ("logit", "state_logit", "dt_logit")):
            diffs = [b - a for a, b in zip(acc[s], acc[d], strict=True)]
            md, (lo, hi), n = t18.paired_ci(diffs)
            res[lbl] = {"mean": md, "ci": [lo, hi], "n": n,
                        "verdict": ("Δt HELPS" if hi < 0 else "Δt HURTS" if lo > 0
                                    else "tied (no Δt signal)")}
        res["fold_gains_raw"] = [float(b - a) for a, b in
                                 zip(acc["state_raw"], acc["dt_raw"], strict=True)]
        res["folds"] = [r[0] for r in rows]
        return res

    hdrs = ["held-out", "pairs", "mean only", "state raw", "+Δt raw", "state clip", "+Δt clip",
            "state logit", "+Δt logit"]
    for label, pm, fl in (("as 3a graded it (age mask on)", per, folds),
                          ("safety target WITHOUT the age mask", per_all, folds_all)):
        r = part_c(pm, fl)
        print(f"\n  PART C — {label}   [{len(r['folds'])} folds]")
        print(render_table(hdrs, r["per_fold"], aligns=["l"] + ["r"] * 8))
        m = r["means"]
        print(f"   means: mean-only={m['mean_only']:.3f}  state raw={m['state_raw']:.3f}  "
              f"state clip={m['state_clip']:.3f}  state logit={m['state_logit']:.3f}")
        for lbl in ("raw", "clip", "logit"):
            v = r[lbl]
            print(f"   paired ({lbl:5s}) mean={v['mean']:+.3f} "
                  f"95% CI=[{v['ci'][0]:+.3f},{v['ci'][1]:+.3f}] (n={v['n']}) -> {v['verdict']}")
        out.setdefault("D2", {})["age_masked" if pm is per else "no_age_mask"] = r

    # =========================================================================================
    # D3 — Part B's identical swing
    # =========================================================================================
    print("\n" + "-" * 90)
    print("D3 — PART B's IDENTICAL SWING, EXPLAINED AND THEN RE-RUN CORRECTLY")
    print("-" * 90)

    all_pairs = [p for d in per for p in per[d]["pairs"]]
    ytr_all = np.array([p["y_j"] for p in all_pairs], float)
    Xtr_all = t18.feats(all_pairs, True)
    sc = StandardScaler().fit(Xtr_all)
    mdl = Ridge(alpha=ALPHA).fit(sc.transform(Xtr_all), ytr_all)
    dts = np.array([p["dt"] for p in all_pairs], float)
    lo_dt, hi_dt = float(np.quantile(dts, 0.10)), float(np.quantile(dts, 0.90))

    swings, b_rows = [], []
    for d in per:
        x0 = per[d]["rows"][0]["x"]
        p_lo, p_hi = (float(mdl.predict(sc.transform(
            np.hstack([x0, [q, q ** 2]]).reshape(1, -1)))[0]) for q in (lo_dt, hi_dt))
        swings.append(p_hi - p_lo)
        b_rows.append([d, f"{p_lo:+.2f}", f"{p_hi:+.2f}", f"{p_hi - p_lo:+.2f}"])
    # the analytic swing: the x0 term cancels because the model is linear in the z-features
    zl = (np.array([lo_dt, lo_dt ** 2]) - sc.mean_[G:]) / sc.scale_[G:]
    zh = (np.array([hi_dt, hi_dt ** 2]) - sc.mean_[G:]) / sc.scale_[G:]
    analytic = float(mdl.coef_[G:] @ (zh - zl))

    print(f"\n  as 3a ran it: ONE global fit on all {len(all_pairs)} pairs, only x0 varying")
    print(render_table(["fold", "ΔAge @ short Δt", "ΔAge @ long Δt", "swing"], b_rows,
                       aligns=["l", "r", "r", "r"]))
    print(f"   SD of the swing across folds = {np.std(swings, ddof=0):.3e}")
    print(f"   analytic swing  w_Δt·Δz(Δt) + w_Δt²·Δz(Δt²) = {analytic:+.5f}   "
          f"(no x0 term exists)")
    print(f"   max |observed − analytic| = {max(abs(s - analytic) for s in swings):.3e}")
    print("   => the five numbers are ONE fit printed five times. They cannot disagree.")

    # the corrected Part B: a fit per held-out fold, swept over that fold's own Δt range
    c_rows, c_swings = [], []
    for hd in folds:
        tr = [p for d in per if d != hd for p in per[d]["pairs"]]
        ytr = np.array([p["y_j"] for p in tr], float)
        Xtr = t18.feats(tr, True)
        sc_h = StandardScaler().fit(Xtr)
        m_h = Ridge(alpha=ALPHA).fit(sc_h.transform(Xtr), ytr)
        dt_h = np.array([p["dt"] for p in per[hd]["pairs"]], float)
        l_h, h_h = float(np.quantile(dt_h, 0.10)), float(np.quantile(dt_h, 0.90))
        x0 = per[hd]["rows"][0]["x"]
        p_lo, p_hi = (float(m_h.predict(sc_h.transform(
            np.hstack([x0, [q, q ** 2]]).reshape(1, -1)))[0]) for q in (l_h, h_h))
        c_swings.append(p_hi - p_lo)
        c_rows.append([hd, f"{l_h:.2f}-{h_h:.2f}", f"{p_lo:+.2f}", f"{p_hi:+.2f}",
                       f"{p_hi - p_lo:+.2f}"])
    print("\n  CORRECTED — leave-one-donor-out fit per fold, each swept over its OWN Δt range")
    print(render_table(["held-out", "Δt swept", "ΔAge @ short", "ΔAge @ long", "swing"],
                       c_rows, aligns=["l", "r", "r", "r", "r"]))
    print(f"   SD of the swing across folds = {np.std(c_swings, ddof=1):.2f} yr   "
          f"mean |swing| = {np.mean(np.abs(c_swings)):.2f} yr")
    print("   NOTE: a swing of hundreds of years is NONPHYSICAL. It does not 'clear' the >2 yr")
    print("   threshold — it is out of range, exactly as Part C's MAE 7.589 is.")
    out["D3"] = {
        "global_fit_swings": {d: float(s) for d, s in zip(per, swings, strict=True)},
        "swing_sd_global_fit": float(np.std(swings, ddof=0)),
        "analytic_swing": analytic,
        "max_abs_observed_minus_analytic": float(max(abs(s - analytic) for s in swings)),
        "dt_sweep": [lo_dt, hi_dt],
        "corrected_lodo_swings": {hd: float(s) for hd, s in zip(folds, c_swings, strict=True)},
        "corrected_swing_sd": float(np.std(c_swings, ddof=1)),
        "corrected_mean_abs_swing": float(np.mean(np.abs(c_swings))),
    }

    # =========================================================================================
    # D4 — bar_verdict at the real geometry, with a SIMULATED TRUE Δt effect
    # =========================================================================================
    print("\n" + "-" * 90)
    print("D4 — bar_verdict AT THE REAL GEOMETRY (REF_GROUND_RULES §5b), NEVER RUN FOR 3a")
    print("-" * 90)
    print("   the real features, the real Δt, the real fold/pair structure. Only the TARGET is")
    print("   simulated, so that a Δt effect is present BY CONSTRUCTION. ρ is the share of the")
    print("   simulated signal carried by Δt; ρ=1 is a target that is a pure function of Δt.")
    print("   3a's rule is graded verbatim: PASS iff the paired 95% CI upper end < 0.")

    sys.path.insert(0, str(REPO))
    from audit_metrics import MIN_PASS_RATE, bar_verdict

    def sim_bar(per_map: dict, fold_list: list[str]) -> dict:
        """Precompute H per (fold, feature set); then only the target changes per trial.

        The simulated truth is ONE GLOBAL function of (state, Δt), standardised on statistics
        pooled over every pair in the arm and then applied identically to the training and the
        held-out side of every fold. Standardising each side separately would make the same Δt
        mean a different target value in train and test -- which is a broken generative model,
        not a hard one, and would charge the estimator for the simulation's own defect.
        """
        donors = list(per_map)
        prep = []
        for hd in fold_list:
            tr = [p for d in per_map if d != hd for p in per_map[d]["pairs"]]
            te = per_map[hd]["pairs"]
            prep.append({
                "s": FrozenRidge(t18.feats(tr, False), t18.feats(te, False)),
                "d": FrozenRidge(t18.feats(tr, True), t18.feats(te, True)),
                "hd": hd,
            })

        u_by = {d: np.array([p["u_j"] for p in per_map[d]["pairs"]], float) for d in donors}
        dt_by = {d: np.array([p["dt"] for p in per_map[d]["pairs"]], float) for d in donors}
        u_all = np.concatenate([u_by[d] for d in donors])
        dt_all = np.concatenate([dt_by[d] for d in donors])

        def gz(a, ref):
            s = float(ref.std())
            return (np.asarray(a, float) - float(ref.mean())) / (s if s > 0 else 1.0)

        def latents(rho: float) -> dict:
            mix_all = (1 - rho) * gz(u_all, u_all) + rho * gz(dt_all, dt_all)
            return {d: np.clip(0.5 + 0.25 * gz((1 - rho) * gz(u_by[d], u_all)
                                               + rho * gz(dt_by[d], dt_all), mix_all), 0, 1)
                    for d in donors}

        rng = np.random.default_rng(SIM_SEED)
        cells = {}
        for rho in RHOS:
            lb = latents(rho)
            lat = [{"tr": np.concatenate([lb[d] for d in donors if d != p["hd"]]),
                    "te": lb[p["hd"]]} for p in prep]
            for sigma in SIGMAS:
                hi_raw, hi_log = [], []
                for _ in range(SIM_TRIALS):
                    ds_raw, dd_raw, ds_log, dd_log = [], [], [], []
                    for p, L in zip(prep, lat, strict=True):
                        ytr = np.clip(L["tr"] + sigma * rng.standard_normal(len(L["tr"])), 0, 1)
                        yte = np.clip(L["te"] + sigma * rng.standard_normal(len(L["te"])), 0, 1)
                        ds_raw.append(mae(p["s"].predict(ytr), yte))
                        dd_raw.append(mae(p["d"].predict(ytr), yte))
                        ds_log.append(mae(expit(p["s"].predict(logit(ytr))), yte))
                        dd_log.append(mae(expit(p["d"].predict(logit(ytr))), yte))
                    hi_raw.append(t18.paired_ci([b - a for a, b in
                                                 zip(ds_raw, dd_raw, strict=True)])[1][1])
                    hi_log.append(t18.paired_ci([b - a for a, b in
                                                 zip(ds_log, dd_log, strict=True)])[1][1])
                cells[f"rho={rho}|sigma={sigma}"] = {
                    "raw": bar_verdict(np.array(hi_raw, float), 0.0, lower_is_better=True),
                    "logit": bar_verdict(np.array(hi_log, float), 0.0, lower_is_better=True),
                }
        return cells

    for label, pm, fl in (("5 folds — the geometry 3a graded", per, folds),
                          ("6 folds — safety target without the age mask", per_all, folds_all)):
        cells = sim_bar(pm, fl)
        rows = []
        for key, v in cells.items():
            rho, sig = (s.split("=")[1] for s in key.split("|"))
            rows.append([rho, sig,
                         f"{v['raw']['pass_rate']:.3f}", v["raw"]["verdict"],
                         f"{v['logit']['pass_rate']:.3f}", v["logit"]["verdict"]])
        print(f"\n  {label}   [{len(fl)} folds, {SIM_TRIALS} trials/cell, "
              f"MIN_PASS_RATE={MIN_PASS_RATE}]")
        print(render_table(["ρ (Δt share)", "noise σ", "pass rate (raw)", "verdict (raw)",
                            "pass rate (logit)", "verdict (logit)"], rows,
                           aligns=["r", "r", "r", "l", "r", "l"]))
        out.setdefault("D4", {})["age_masked" if pm is per else "no_age_mask"] = cells

    # =========================================================================================
    # D5 — the pre-registered expectations, graded
    # =========================================================================================
    print("\n" + "-" * 90)
    print("D5 — PRE-REGISTERED EXPECTATIONS, GRADED")
    print("-" * 90)
    grades = {}

    y1 = d1.get("Y1")
    if y1:
        nested = (y1["dt_heldout_range"][0] >= y1["dt_train_range"][0]
                  and y1["dt_heldout_range"][1] <= y1["dt_train_range"][1])
        grades["P1_Y1_dt_nested_inside_training"] = {
            "held": bool(nested), "outside_pairs": y1["heldout_pairs_outside_train_dt_support"]}
        others = [d1[h]["mean_abs_dt_contrib_partC"] for h in folds if h != "Y1"] or [0.0]
        ratio = y1["mean_abs_dt_contrib_partC"] / max(max(others), 1e-12)
        grades["P2_dt_block_carries_divergence"] = {
            "held": bool(ratio > 10.0), "Y1_over_next_worst": float(ratio)}
    grades["P3_partB_swing_identical_and_analytic"] = {
        "held": bool(out["D3"]["swing_sd_global_fit"] < 1e-9
                     and out["D3"]["max_abs_observed_minus_analytic"] < 1e-6),
        "sd": out["D3"]["swing_sd_global_fit"],
        "max_abs_dev": out["D3"]["max_abs_observed_minus_analytic"]}
    m5 = out["D2"]["age_masked"]["means"]
    grades["P4_state_only_arm_also_broken"] = {
        "held": bool(m5["state_raw"] > m5["mean_only"]),
        "state_raw": m5["state_raw"], "mean_only": m5["mean_only"]}
    best = out["D4"]["age_masked"][f"rho={RHOS[-1]}|sigma={SIGMAS[-1]}"]
    grades["P5_bar_unresolvable_even_at_pure_dt"] = {
        "held": bool(best["raw"]["verdict"] == "UNRESOLVABLE"),
        "pass_rate_raw": best["raw"]["pass_rate"],
        "pass_rate_logit": best["logit"]["pass_rate"]}

    print(render_table(["expectation", "held?", "evidence"],
                       [[k, "YES" if v["held"] else "NO",
                         ", ".join(f"{a}={b:.4g}" if isinstance(b, float) else f"{a}={b}"
                                   for a, b in v.items() if a != "held")]
                        for k, v in grades.items()], aligns=["l", "l", "l"]))
    out["D5_pre_registered"] = grades

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   This diagnoses the instrument. It does NOT re-take 3a's verdict; a repaired 3a")
    print("   run is a separate, pre-registered act.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
