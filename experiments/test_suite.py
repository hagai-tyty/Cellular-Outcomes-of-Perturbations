"""
Test suite 8-9 (ΔAge lab notebook) — does the model earn its keep, and do prior-knowledge
gene embeddings help? Batched into one file with subcommands to save time.

  python test_suite.py fate_baseline    # Test 8    model fate vs logistic regression (no embeddings)
  python test_suite.py indist_vs_donor  # Test 8.1  in-dist cells vs out-of-donor (fitting vs generalization)
  python test_suite.py fate_cal_disc    # Test 8.2  fate discrimination vs calibration (which failure?)
  python test_suite.py string_dage      # Test 9    gene embeddings vs raw, on ΔAge
  python test_suite.py string_fate      # Test 9.1  gene embeddings vs raw, on fate
  python test_suite.py input_ablation   # Test 11   cell-only vs perturbation-only vs both
  python test_suite.py per_donor        # Test 12   per-donor jackknife (are aggregates donor-robust?)
  # Test 10 (can ridge predict full transcriptome response) — PENDING target definition

HONEST GUARDRAILS (see notebook):
- Data-scale law: flexible models need MORE data to beat linear; a win at scale does NOT imply
  a win on 6 donors. None of these are assumed wins.
- STRING is unreachable from the build env; we use gene2vec (200-dim, co-expression-based,
  external/frozen, ~94% panel coverage) as a PROXY for the prior-knowledge-embedding principle.
  Not identical to STRING (protein-interaction). If a signal appears, confirm with STRING later.
- The clock was trained on NATURAL aging; reprogramming is EXTREME rejuvenation outside its
  training domain — read clock-based results with that caveat.

All tests use leave-one-donor-out, the same held-out cells, deterministic. Reads folds from
repo root or runs/. gene2vec downloads once to gene2vec_cache.txt.
"""
from __future__ import annotations

import math
import sys
import urllib.request
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import SAFE_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.common.panel import GenePanel
from cellfate.evaluation.baselines import ModelEstimator
from cellfate.evaluation.data import gather_split
from cellfate.inference import Predictor

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
GENE2VEC_URL = ("https://raw.githubusercontent.com/jingcheng-du/Gene2vec/master/"
                "pre_trained_emb/gene2vec_dim_200_iter_9.txt")
CACHE = "gene2vec_cache.txt"


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def load_gene2vec() -> dict:
    if not Path(CACHE).exists():
        print("   downloading gene2vec embeddings (once)...")
        req = urllib.request.Request(GENE2VEC_URL, headers={"User-Agent": "Mozilla/5.0"})
        Path(CACHE).write_bytes(urllib.request.urlopen(req, timeout=120).read())
    emb = {}
    for ln in Path(CACHE).read_text().splitlines():
        parts = ln.split()
        if len(parts) > 2:
            emb[parts[0]] = np.array(parts[1:], dtype=np.float32)
    return emb


def embed_matrix(genes, emb) -> np.ndarray:
    """(n_genes x 200) rows = each panel gene's embedding (0 if missing)."""
    dim = len(next(iter(emb.values())))
    E = np.zeros((len(genes), dim), dtype=np.float32)
    for i, g in enumerate(genes):
        if g in emb:
            E[i] = emb[g]
    return E


def paired_ci(diffs):
    diffs = [d for d in diffs if np.isfinite(d)]
    n = len(diffs)
    if n < 2:
        return float("nan"), (float("nan"), float("nan")), n
    md = sum(diffs) / n
    sd = math.sqrt(sum((x - md) ** 2 for x in diffs) / (n - 1))
    se = sd / math.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n


def feats(tr, te):
    sx = StandardScaler().fit(tr.X)
    sdt = StandardScaler().fit(tr.dose_time)
    ftr = np.hstack([sx.transform(tr.X), np.asarray(tr.fp, float), sdt.transform(tr.dose_time)])
    fte = np.hstack([sx.transform(te.X), np.asarray(te.fp, float), sdt.transform(te.dose_time)])
    return ftr, fte


# ---------------------------------------------------------------------------- #
# 8.1 — fate: model vs logistic regression
# ---------------------------------------------------------------------------- #
def fate_baseline():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 8 — does the model's FATE beat logistic regression? (leave-one-donor-out)")
    print("safe-class PR-AUC (higher better); model = neural net, logreg = linear baseline.")

    rows, d_model, d_lr = [], [], []
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, REGIME, "train")
            te = gather_split(paths, REGIME, "test")
        except Exception:  # noqa: BLE001
            rows.append([d, "n/a", "n/a", "n/a"])
            continue
        y = te.y_cls.astype(int)
        safe_true = (y == SAFE_IDX).astype(int)
        if safe_true.sum() in (0, len(safe_true)):
            rows.append([d, "n/a (1 class)", "n/a", "n/a"])
            continue
        ftr, fte = feats(tr, te)
        lr = LogisticRegression(max_iter=2000, class_weight="balanced").fit(ftr, tr.y_cls.astype(int))
        p_lr_safe = lr.predict_proba(fte)[:, list(lr.classes_).index(SAFE_IDX)]
        rows_m = ModelEstimator(Predictor(root)).rows(te.X, te.fp, te.dose_time)
        p_model_safe = np.array([r["S"] for r in rows_m])
        pr_m = average_precision_score(safe_true, p_model_safe)
        pr_lr = average_precision_score(safe_true, p_lr_safe)
        d_model.append(pr_m)
        d_lr.append(pr_lr)
        rows.append([d, f"{pr_m:.3f}", f"{pr_lr:.3f}", "model" if pr_m > pr_lr else "logreg"])

    print("\n" + render_table(["fold", "model PR-AUC", "logreg PR-AUC", "winner"],
                              rows, aligns=["l", "r", "r", "l"]))
    if len(d_model) >= 2:
        diffs = [m - lg for m, lg in zip(d_model, d_lr, strict=True)]
        md, (lo, hi), n = paired_ci(diffs)
        print(f"\n   aggregate: model={np.mean(d_model):.3f}  logreg={np.mean(d_lr):.3f}")
        print(f"   paired (model−logreg): mean={md:+.3f} 95% CI=[{lo:+.3f},{hi:+.3f}] (n={n})")
        v = ("model BEATS logreg" if lo > 0 else "model WORSE than logreg" if hi < 0
             else "tied (noise)")
        print(f"   -> {v}")
        print("\n   READ: model beats logreg beyond noise -> fate is a real contribution.")
        print("         tied/worse -> the fate head doesn't earn its keep vs a linear classifier.")


# ---------------------------------------------------------------------------- #
# shared: per-fold panel gene names (X-column order)
# ---------------------------------------------------------------------------- #
def _panel_genes(root):
    return list(GenePanel.load(f"{root}/panel.json").genes)


# ---------------------------------------------------------------------------- #
# 9.1 — ΔAge: ridge-raw vs ridge-emb vs mlp-emb
# ---------------------------------------------------------------------------- #
def string_dage():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 9 — do gene embeddings help on ΔAge? (ridge-raw vs ridge-emb vs MLP-emb)")
    print("gene2vec embeddings (proxy for STRING). MAE in years, lower better.")
    emb = load_gene2vec()

    names = ["ridge_raw", "ridge_emb", "mlp_emb"]
    per = {}
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, REGIME, "train")
            te = gather_split(paths, REGIME, "test")
            genes = _panel_genes(root)
        except Exception:  # noqa: BLE001
            continue
        mtr, mte = tr.mask, te.mask
        if mtr.sum() < 10 or mte.sum() < 3:
            continue
        E = embed_matrix(genes, emb)                     # (n_genes x 200)
        Xtr_raw, Xte_raw = tr.X[mtr], te.X[mte]
        Xtr_emb, Xte_emb = tr.X[mtr] @ E, te.X[mte] @ E  # expression-weighted gene embeddings
        ytr, yte = tr.y_age[mtr], te.y_age[mte]
        sr = StandardScaler().fit(Xtr_raw)
        se = StandardScaler().fit(Xtr_emb)
        r_raw = Ridge(alpha=1.0).fit(sr.transform(Xtr_raw), ytr).predict(sr.transform(Xte_raw))
        r_emb = Ridge(alpha=1.0).fit(se.transform(Xtr_emb), ytr).predict(se.transform(Xte_emb))
        m_emb = MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=400, random_state=0,
                             early_stopping=True).fit(se.transform(Xtr_emb), ytr).predict(
            se.transform(Xte_emb))
        per[d] = {"ridge_raw": float(np.abs(r_raw - yte).mean()),
                  "ridge_emb": float(np.abs(r_emb - yte).mean()),
                  "mlp_emb": float(np.abs(m_emb - yte).mean())}

    if not per:
        print("\n   No folds found.")
        return
    rows = [[d] + [f"{per[d][n]:.2f}" for n in names] for d in per]
    print("\n" + render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
    agg = {n: np.mean([per[d][n] for d in per]) for n in names}
    print("   aggregate MAE: " + "   ".join(f"{n}={agg[n]:.2f}" for n in names))
    for n in ["ridge_emb", "mlp_emb"]:
        diffs = [per[d][n] - per[d]["ridge_raw"] for d in per]
        md, (lo, hi), k = paired_ci(diffs)
        v = ("EMB BEATS raw" if hi < 0 else "EMB WORSE than raw" if lo > 0 else "tied")
        print(f"   {n} vs ridge_raw: mean={md:+.2f} CI=[{lo:+.2f},{hi:+.2f}] -> {v}")
    print("\n   READ: for a LINEAR ΔAge target, embeddings (which encode interactions) are")
    print("         expected to NOT help. If they don't -> confirms ΔAge is purely linear.")


# ---------------------------------------------------------------------------- #
# 9.2 — fate: logreg-raw vs logreg-emb vs mlp-emb
# ---------------------------------------------------------------------------- #
def string_fate():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 9.1 — do gene embeddings help FATE? (logreg-raw vs logreg-emb vs MLP-emb)")
    print("safe-class PR-AUC, higher better. gene2vec embeddings (proxy for STRING).")
    emb = load_gene2vec()

    names = ["logreg_raw", "logreg_emb", "mlp_emb"]
    per = {}
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, REGIME, "train")
            te = gather_split(paths, REGIME, "test")
            genes = _panel_genes(root)
        except Exception:  # noqa: BLE001
            continue
        y_te = te.y_cls.astype(int)
        safe_true = (y_te == SAFE_IDX).astype(int)
        if safe_true.sum() in (0, len(safe_true)):
            continue
        E = embed_matrix(genes, emb)
        Xtr_raw, Xte_raw = tr.X, te.X
        Xtr_emb, Xte_emb = tr.X @ E, te.X @ E
        ytr = tr.y_cls.astype(int)
        sr = StandardScaler().fit(Xtr_raw)
        se = StandardScaler().fit(Xtr_emb)

        def prauc(model, Xte_s, st=safe_true):
            idx = list(model.classes_).index(SAFE_IDX)
            return average_precision_score(st, model.predict_proba(Xte_s)[:, idx])
        lr_raw = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
            sr.transform(Xtr_raw), ytr)
        lr_emb = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
            se.transform(Xtr_emb), ytr)
        mlp = MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=400, random_state=0,
                            early_stopping=True).fit(se.transform(Xtr_emb), ytr)
        per[d] = {"logreg_raw": prauc(lr_raw, sr.transform(Xte_raw)),
                  "logreg_emb": prauc(lr_emb, se.transform(Xte_emb)),
                  "mlp_emb": prauc(mlp, se.transform(Xte_emb))}

    if not per:
        print("\n   No folds with fate variation found.")
        return
    rows = [[d] + [f"{per[d][n]:.3f}" for n in names] for d in per]
    print("\n" + render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
    agg = {n: np.mean([per[d][n] for d in per]) for n in names}
    print("   aggregate PR-AUC: " + "   ".join(f"{n}={agg[n]:.3f}" for n in names))
    for n in ["logreg_emb", "mlp_emb"]:
        diffs = [per[d][n] - per[d]["logreg_raw"] for d in per]
        md, (lo, hi), k = paired_ci(diffs)
        v = ("EMB BEATS raw" if lo > 0 else "EMB WORSE than raw" if hi < 0 else "tied")
        print(f"   {n} vs logreg_raw: mean={md:+.3f} CI=[{lo:+.3f},{hi:+.3f}] -> {v}")
    print("\n   READ: fate is classification, plausibly interaction-dependent — the place")
    print("         embeddings COULD help. If they help beyond noise -> pursue STRING; else")
    print("         fate is limited by data/donors, not representation.")


def _ece(p_safe, is_safe, bins=10):
    """Expected calibration error on the safe-class probability."""
    p_safe = np.asarray(p_safe, float)
    is_safe = np.asarray(is_safe, float)
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for i in range(bins):
        m = (p_safe >= edges[i]) & (p_safe < edges[i + 1] if i < bins - 1 else p_safe <= 1.0)
        if m.sum() == 0:
            continue
        ece += (m.sum() / len(p_safe)) * abs(p_safe[m].mean() - is_safe[m].mean())
    return float(ece)


def _model_safe_probs(root, sd):
    rows = ModelEstimator(Predictor(root)).rows(sd.X, sd.fp, sd.dose_time)
    return np.array([r["S"] for r in rows]), np.array([r["mu_age"] for r in rows])


# ---------------------------------------------------------------------------- #
# 10 — in-distribution (held-out CELLS) vs out-of-donor (held-out DONOR)
# ---------------------------------------------------------------------------- #
def indist_vs_donor():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 8.1 — IN-DIST (held-out cells) vs OUT-OF-DONOR (held-out donor)")
    print("isolates FITTING failure from GENERALIZATION failure. fate=PR-AUC(safe), ΔAge=MAE.")
    rows = []
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            val = gather_split(paths, REGIME, "val")     # in-dist held-out cells
            te = gather_split(paths, REGIME, "test")     # out-of-donor
        except Exception:  # noqa: BLE001
            rows.append([d, "n/a", "n/a", "n/a", "n/a"])
            continue
        out = []
        for sd in (val, te):
            ps, mu = _model_safe_probs(root, sd)
            m = sd.mask
            mae = float(np.abs(mu[m] - sd.y_age[m]).mean()) if m.sum() else float("nan")
            st = (sd.y_cls.astype(int) == SAFE_IDX).astype(int)
            pr = (average_precision_score(st, ps) if 0 < st.sum() < len(st) else float("nan"))
            out += [pr, mae]
        rows.append([d, f"{out[0]:.3f}" if np.isfinite(out[0]) else "n/a", f"{out[1]:.2f}",
                     f"{out[2]:.3f}" if np.isfinite(out[2]) else "n/a", f"{out[3]:.2f}"])
    print("\n" + render_table(
        ["fold", "fate PR in-dist", "ΔAge MAE in-dist", "fate PR out-donor", "ΔAge MAE out-donor"],
        rows, aligns=["l", "r", "r", "r", "r"]))
    print("\n   WHAT THIS PINPOINTS: if fate PR is high in-dist but low out-of-donor while ΔAge")
    print("   MAE is similar both -> the failure is DONOR SHIFT (fix=more donors), not the head.")


# ---------------------------------------------------------------------------- #
# 11 — input ablation: cell-only vs perturbation-only vs both
# ---------------------------------------------------------------------------- #
def input_ablation():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 11 — INPUT ABLATION: does using cell+perturbation beat either alone?")
    print("ΔAge MAE (lower better) and fate PR-AUC (higher better), per fold. x=cell, u=perturbation.")
    rows_age, rows_fate = [], []
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            tr = gather_split(paths, REGIME, "train")
            te = gather_split(paths, REGIME, "test")
        except Exception:  # noqa: BLE001
            continue
        U_tr = np.hstack([np.asarray(tr.fp, float), tr.dose_time])
        U_te = np.hstack([np.asarray(te.fp, float), te.dose_time])
        views = {"x_only": (tr.X, te.X), "u_only": (U_tr, U_te),
                 "x+u": (np.hstack([tr.X, U_tr]), np.hstack([te.X, U_te]))}
        # ΔAge (ridge)
        ma, mte = tr.mask, te.mask
        age_row = [d]
        for _, (Xt, Xe) in views.items():
            s = StandardScaler().fit(Xt[ma])
            pred = Ridge(alpha=1.0).fit(s.transform(Xt[ma]), tr.y_age[ma]).predict(s.transform(Xe[mte]))
            age_row.append(f"{np.abs(pred - te.y_age[mte]).mean():.2f}")
        rows_age.append(age_row)
        # fate (logreg)
        st = (te.y_cls.astype(int) == SAFE_IDX).astype(int)
        if 0 < st.sum() < len(st):
            fate_row = [d]
            for _, (Xt, Xe) in views.items():
                s = StandardScaler().fit(Xt)
                lr = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
                    s.transform(Xt), tr.y_cls.astype(int))
                idx = list(lr.classes_).index(SAFE_IDX)
                fate_row.append(f"{average_precision_score(st, lr.predict_proba(s.transform(Xe))[:, idx]):.3f}")
            rows_fate.append(fate_row)
    print("\n  ΔAge MAE (lower better)")
    print(render_table(["fold", "x_only", "u_only", "x+u"], rows_age, aligns=["l", "r", "r", "r"]))
    if rows_fate:
        print("\n  fate PR-AUC (higher better)")
        print(render_table(["fold", "x_only", "u_only", "x+u"], rows_fate, aligns=["l", "r", "r", "r"]))
    print("\n   WHAT THIS PINPOINTS: if x+u ~ x_only, the perturbation adds nothing (OSKM ~fixed)")
    print("   and the state\u00d7perturbation interaction is NOT load-bearing on this data.")


# ---------------------------------------------------------------------------- #
# 12 — fate: discrimination vs calibration (+ does recalibration fix it?)
# ---------------------------------------------------------------------------- #
def fate_cal_disc():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 8.2 — fate DISCRIMINATION vs CALIBRATION (which kind of failure?)")
    print("PR-AUC/ROC = ranking; ECE = probability accuracy. + Platt recalibration on calib split.")
    from sklearn.linear_model import LogisticRegression as PlattLR
    from sklearn.metrics import roc_auc_score
    rows = []
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            cal = gather_split(paths, REGIME, "calib")
            te = gather_split(paths, REGIME, "test")
        except Exception:  # noqa: BLE001
            rows.append([d] + ["n/a"] * 5)
            continue
        ps_te, _ = _model_safe_probs(root, te)
        st_te = (te.y_cls.astype(int) == SAFE_IDX).astype(int)
        if not (0 < st_te.sum() < len(st_te)):
            rows.append([d] + ["n/a"] * 5)
            continue
        pr = average_precision_score(st_te, ps_te)
        roc = roc_auc_score(st_te, ps_te)
        ece0 = _ece(ps_te, st_te)
        # Platt recalibration: fit logistic on calib safe-probs -> safe outcome, apply to test
        ps_cal, _ = _model_safe_probs(root, cal)
        st_cal = (cal.y_cls.astype(int) == SAFE_IDX).astype(int)
        if 0 < st_cal.sum() < len(st_cal):
            platt = PlattLR(max_iter=1000).fit(ps_cal.reshape(-1, 1), st_cal)
            ps_re = platt.predict_proba(ps_te.reshape(-1, 1))[:, 1]
            ece1 = _ece(ps_re, st_te)
        else:
            ece1 = float("nan")
        rows.append([d, f"{pr:.3f}", f"{roc:.3f}", f"{ece0:.3f}",
                     f"{ece1:.3f}" if np.isfinite(ece1) else "n/a",
                     "recal helps" if np.isfinite(ece1) and ece1 < ece0 else "no"])
    print("\n" + render_table(
        ["fold", "PR-AUC", "ROC-AUC", "ECE raw", "ECE recal", "recal?"],
        rows, aligns=["l", "r", "r", "r", "r", "l"]))
    print("\n   WHAT THIS PINPOINTS: high ROC/PR but high ECE -> fate RANKS ok but is MISCALIBRATED")
    print("   -> Platt/temperature recalibration is the cheap fix; RES could use recalibrated fate.")


# ---------------------------------------------------------------------------- #
# 13 — per-donor jackknife: which donors drive the numbers?
# ---------------------------------------------------------------------------- #
def per_donor():
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    print("\nTEST 12 — PER-DONOR JACKKNIFE: are aggregates robust or hostage to 1-2 donors?")
    print("model ΔAge MAE + fate PR-AUC per fold, then aggregate with each donor DROPPED.")
    per = {}
    for d in DONORS:
        root = resolve_root(f"cellfate_loocv_{d}")
        try:
            paths = ArtifactPaths.of(root)
            te = gather_split(paths, REGIME, "test")
        except Exception:  # noqa: BLE001
            continue
        ps, mu = _model_safe_probs(root, te)
        m = te.mask
        mae = float(np.abs(mu[m] - te.y_age[m]).mean()) if m.sum() else float("nan")
        st = (te.y_cls.astype(int) == SAFE_IDX).astype(int)
        pr = average_precision_score(st, ps) if 0 < st.sum() < len(st) else float("nan")
        per[d] = (mae, pr)
    rows = [[d, f"{per[d][0]:.2f}", f"{per[d][1]:.3f}" if np.isfinite(per[d][1]) else "n/a"]
            for d in per]
    print("\n" + render_table(["fold", "ΔAge MAE", "fate PR-AUC"], rows, aligns=["l", "r", "r"]))
    print("\n  aggregate with each donor DROPPED (jackknife):")
    jk = []
    ds = list(per)
    for drop in ds:
        keep = [x for x in ds if x != drop]
        maes = [per[x][0] for x in keep if np.isfinite(per[x][0])]
        prs = [per[x][1] for x in keep if np.isfinite(per[x][1])]
        jk.append([f"drop {drop}", f"{np.mean(maes):.2f}", f"{np.mean(prs):.3f}" if prs else "n/a"])
    jk.append(["ALL", f"{np.nanmean([per[x][0] for x in ds]):.2f}",
               f"{np.nanmean([per[x][1] for x in ds]):.3f}"])
    print(render_table(["excluded", "ΔAge MAE", "fate PR-AUC"], jk, aligns=["l", "r", "r"]))
    print("\n   WHAT THIS PINPOINTS: if dropping one donor (e.g. N3) swings the aggregate a lot,")
    print("   the conclusions are hostage to that donor -> flag in the writeup; don't over-claim.")


CMDS = {"fate_baseline": fate_baseline, "string_dage": string_dage, "string_fate": string_fate,
        "indist_vs_donor": indist_vs_donor, "input_ablation": input_ablation,
        "fate_cal_disc": fate_cal_disc, "per_donor": per_donor}

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd in CMDS:
        CMDS[cmd]()
    else:
        print("usage: python test_suite.py <cmd>")
        print("  Test 8   fate_baseline    (+ 8.1 indist_vs_donor, 8.2 fate_cal_disc)")
        print("  Test 9   string_dage      (+ 9.1 string_fate)")
        print("  Test 11  input_ablation")
        print("  Test 12  per_donor")
        print("  Test 10 (full transcriptome) — pending target definition; see notebook")
