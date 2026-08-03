"""STAGE 1.5.4 — can a model LEARN methylation age from RNA? (GSE165177 x GSE165179)

    python experiments/diag_learned_clock.py "D:\\GSE165177" "D:\\GSE165179"

READ-ONLY. Writes `results/diag_learned_clock_results.json`. `src/` untouched, no label moves.

THE QUESTION, AND HOW IT DIFFERS FROM M-2a
------------------------------------------
M-2a asked whether the **existing Fleischer clock's scalar output** tracks methylation age here.
It does not (rho_partial 0.267 / 0.516 -> SPLIT). That is a fact about one clock: it was fitted on
quiescent fibroblasts, against CHRONOLOGICAL age, over 33,155 genes from 133 samples.

This asks the second question, which nobody has: **can a model TRAINED on the transcriptome predict
methylation age on these cells?** If yes, ΔAge labels become obtainable at scale. If no, the RNA
route is closed on evidence rather than inferred from one clock.

DESIGN (pre-registered in plans/STAGE_1_5_4_LEARNED_CLOCK.md before this was written)
-------------------------------------------------------------------------------------
  target     the clock's LINEAR PREDICTOR (not years), both Horvath clocks, both always reported
  split      leave-one-donor-out over O1/O2/O3 -- everything fitted INSIDE the fold
  estimand   pooled rho_partial over all 68 out-of-fold predictions
  bar        >= 0.50 on BOTH clocks. SPLIT counts as failure.

Why the estimand is pooled while the split is per-donor: they are different choices. LODO prevents
leakage; it is not what the statistic is computed on. sec 5b's record already settles which geometry
works here -- pooled n=68 is RESOLVABLE at 0.9940, every smaller geometry tried was not -- and it
makes the number directly comparable to the Fleischer clock's 0.267 / 0.516.

THE TWO GUARDS
--------------
G1  rho_partial is the ONLY pass criterion. Both modalities move with reprogramming progress, so a
    model that learns only "how far along is this cell" would score high on a raw correlation while
    carrying no age information. rho_all and rho_within are reported, never graded.

G2  LABEL-SHUFFLE NULL. With ~20k features and ~45 training samples a pipeline can manufacture
    correlation from nothing. The whole procedure is re-run with TRAINING targets shuffled inside
    each fold. If that does not collapse to ~0, the pipeline is broken and NO positive result from
    it may be believed. This is checked BEFORE the real result is read.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

# Frozen in plans/STAGE_1_5_4_LEARNED_CLOCK.md sec 5. Re-used from M-2a deliberately: the bar, the
# metric and the geometry are identical, which removes any freedom to pick a friendlier one.
RHO_BAR = 0.50
N_PERM = 20               # permutation-null draws per family x clock (see G2 in the plan)
PERM_Q = 0.95             # a real rho must beat this quantile of its OWN null
ALPHAS = (1e1, 1e2, 1e3, 1e4, 1e5, 1e6)
PCA_K = 10
SEED = 0


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# Pure logic -- unit-tested with no repo data present                          #
# --------------------------------------------------------------------------- #
def ridge_fit(X: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    """Closed-form ridge on centred data. Returns coefficients for CENTRED X."""
    n, p = X.shape
    if p <= n:
        A = X.T @ X + alpha * np.eye(p)
        return np.linalg.solve(A, X.T @ y)
    # p >> n: solve in sample space (Woodbury), which is the regime this stage lives in
    K = X @ X.T + alpha * np.eye(n)
    return X.T @ np.linalg.solve(K, y)


def inner_cv_alpha(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                   alphas=ALPHAS) -> float:
    """Pick alpha by leave-one-GROUP-out inside the training set only.

    Grouped, not random: the training folds are whole donors, so a random inner split would let a
    donor appear on both sides and choose an alpha that does not generalise across donors -- the
    same leak the outer split exists to prevent, one level down.
    """
    uniq = np.unique(groups)
    if len(uniq) < 2:
        return float(alphas[len(alphas) // 2])
    best, best_err = alphas[0], np.inf
    for a in alphas:
        err = []
        for g in uniq:
            tr, te = groups != g, groups == g
            mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-8
            ym = y[tr].mean()
            w = ridge_fit((X[tr] - mu) / sd, y[tr] - ym, a)
            err.append(np.mean(((X[te] - mu) / sd @ w + ym - y[te]) ** 2))
        m = float(np.mean(err))
        if m < best_err:
            best, best_err = a, m
    return float(best)


def _ranks(v: np.ndarray) -> np.ndarray:
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v), float)
    r[order] = np.arange(1, len(v) + 1)
    # average ties
    for val in np.unique(v):
        m = v == val
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    rx, ry = _ranks(x), _ranks(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    den = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float(rx @ ry / den) if den else float("nan")


def _resid_on(a: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """Residual of `a` after least-squares removal of every column of `Z` (intercept included)."""
    Z = np.column_stack([np.ones(len(a)), Z])
    beta, *_ = np.linalg.lstsq(Z, a, rcond=None)
    return a - Z @ beta


def partial_spearman(x, y, z, donor=None) -> float:
    """Spearman between x and y, partialling out z -- and DONOR when supplied.

    Donor must be partialled out whenever `x` is a LEAVE-ONE-DONOR-OUT prediction, and this is a
    correctness requirement rather than a refinement. A fold trained without donor d predicts d
    with something close to the mean of the OTHER donors, so if the donor means are y1<y2<y3 then
    d=1 receives a high prediction and d=3 a low one: **the prediction is anti-correlated with the
    donor mean by construction.** Measured on this data with the training labels shuffled, that
    artefact alone produces rho_partial as low as -0.45 (G2, first run) -- a large correlation from
    a model that has learned nothing. Removing donor removes the artefact and leaves the question
    that matters here anyway: does RNA track methylation age WITHIN a donor, beyond reprogramming
    progress? Between-donor variance is only 3 points at ages 53/53/38 and can settle nothing.
    """
    rx, ry = _ranks(np.asarray(x, float)), _ranks(np.asarray(y, float))
    Z = _ranks(np.asarray(z, float)).reshape(-1, 1)
    if donor is not None:
        d = np.asarray(donor)
        # one-hot minus a reference level; the intercept in `_resid_on` carries the dropped one
        levels = sorted(set(d.tolist()))[1:]
        if levels:
            Z = np.column_stack([Z] + [(d == lv).astype(float) for lv in levels])
    return spearman(_resid_on(rx, Z), _resid_on(ry, Z))


def verdict(rhos: dict[str, float], bar: float = RHO_BAR) -> str:
    """PASS only if every clock clears the bar. One of two clearing it is SPLIT, which sec 6
    counts as a failure -- the same rule M-2a operated under."""
    ok = [v >= bar for v in rhos.values()]
    return "PASS" if all(ok) else ("SPLIT" if any(ok) else "FAIL")


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def build_matrix(rna_dir: Path, meth_dir: Path):
    """Join RNA to methylation, returning (X, targets per clock, pluripotency, donor, keys)."""
    m2a = _load("diag_m2a", ROOT / "experiments" / "diag_m2a_calibratability.py")
    from cellfate.data.normalize import normalize_counts

    samples, genes, lin = m2a.load_rna(rna_dir)
    norm = normalize_counts(lin)                       # log1p CP10k, per sample
    meth = m2a.meth_ages(meth_dir)                     # {clock: {key: lp}}
    dcv = _load("dcv", ROOT / "experiments" / "diag_clock_validity.py")
    idx = [i for i, g in enumerate(genes) if g in dcv.OSKM_PLURIPOTENCY]
    plu_all = norm[:, idx].mean(axis=1) if idx else np.zeros(len(samples))

    # M-2a's join, exactly: RNA and methylation share the SAME sample title, so pair on the raw
    # title and only then average exp1/exp2 into (donor, arm, day) conditions. Averaging the
    # EXPRESSION as well as the target keeps one row per condition, which is the unit the 68-point
    # geometry -- and its registered resolvability -- is defined on.
    any_clock = next(iter(meth))
    groups: dict[tuple, list[int]] = {}
    for i, s in enumerate(samples):
        p = m2a.parse_title(s)
        if p and s in meth[any_clock]:
            groups.setdefault((p["donor"], p["arm"], p["day"]), []).append(i)

    keys = sorted(groups)
    X = np.vstack([norm[groups[k], :].mean(axis=0) for k in keys])
    plu = np.array([float(plu_all[groups[k]].mean()) for k in keys])
    donor = np.array([k[0] for k in keys])
    title_of = {i: s for i, s in enumerate(samples)}
    y = {c: np.array([float(np.mean([meth[c][title_of[i]] for i in groups[k]])) for k in keys])
         for c in meth}
    return X, y, plu, donor, [f"{a}_{b}_{d}" for a, b, d in keys], genes


def lodo_predict(X, y, donor, family: str, genes, clock_genes, shuffle=False, seed=SEED):
    """Out-of-fold predictions. Everything -- scaling, feature choice, alpha -- fitted in-fold."""
    rng = np.random.default_rng(seed)
    pred = np.full(len(y), np.nan)
    used_alpha = {}
    if family == "clockgenes":
        sel = np.array([i for i, g in enumerate(genes) if g in clock_genes])
        Xf = X[:, sel] if len(sel) else X
    else:
        Xf = X
    for d in np.unique(donor):
        tr, te = donor != d, donor == d
        ytr = y[tr].copy()
        if shuffle:
            ytr = rng.permutation(ytr)          # destroy the pairing INSIDE the training set only
        mu, sd = Xf[tr].mean(0), Xf[tr].std(0) + 1e-8
        Ztr, Zte = (Xf[tr] - mu) / sd, (Xf[te] - mu) / sd
        if family == "pca":
            # Components come from CENTRED training data, so both matrices must be centred by the
            # SAME training mean before projection. Projecting uncentred data onto centred
            # components mixes the mean into every score -- an implementation bug that left a
            # residual artefact visible in G2 (pca -0.171/-0.261 while the other families sat
            # near zero).
            m = Ztr.mean(0)
            _, _, vt = np.linalg.svd(Ztr - m, full_matrices=False)
            P = vt[:PCA_K].T
            Ztr, Zte = (Ztr - m) @ P, (Zte - m) @ P
        a = inner_cv_alpha(Ztr, ytr, donor[tr])
        used_alpha[str(d)] = a
        ym = ytr.mean()
        w = ridge_fit(Ztr - Ztr.mean(0), ytr - ym, a)
        pred[te] = (Zte - Ztr.mean(0)) @ w + ym
    return pred, used_alpha


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    rna_dir, meth_dir = Path(sys.argv[1]), Path(sys.argv[2])

    X, Y, plu, donor, keys, genes = build_matrix(rna_dir, meth_dir)
    print(f"\n[shape before statistic]  {X.shape[0]} conditions x {X.shape[1]} genes")
    print(f"  donors: {dict(zip(*np.unique(donor, return_counts=True), strict=True))}")
    print(f"  clocks: {list(Y)}")
    if X.shape[0] != 68:
        print(f"  WARNING: expected 68 joined conditions, got {X.shape[0]}")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    clock_genes = set(clock["weights"])

    out = {"script": "diag_learned_clock", "utc": datetime.now(UTC).isoformat(),
           "n_conditions": int(X.shape[0]), "n_genes": int(X.shape[1]),
           "bar": RHO_BAR, "n_perm": N_PERM, "perm_q": PERM_Q, "families": {}}

    # ---- G2: a PERMUTATION NULL, not a fixed threshold -------------------------------------- #
    # The pre-registered form was `|rho| <= 0.20` on the worst of 6 comparisons. That bar is
    # UNRESOLVABLE: at n=68 a correct pipeline has SD(rho) ~ 0.122, so P(|rho| > 0.20) ~ 0.102
    # per comparison and ~0.474 across six -- it would fire on a perfectly sound pipeline almost
    # half the time. Replaced, while still blind to every real rho, by the null the design implies:
    # re-run the whole procedure with training labels permuted and require the real value to beat
    # its OWN null's 95th percentile. This self-calibrates, needs no chosen constant, and prices in
    # both the LODO mean-reversion artefact and the p >> n regime.
    print(f"\n[G2] permutation null: {N_PERM} draws per family x clock, blind to the real result")
    null = {}
    for fam in ("full", "clockgenes", "pca"):
        null[fam] = {}
        for cname, yv in Y.items():
            draws = []
            for s_i in range(N_PERM):
                pr, _ = lodo_predict(X, yv, donor, fam, genes, clock_genes,
                                     shuffle=True, seed=SEED + s_i)
                draws.append(partial_spearman(pr, yv, plu, donor))
            draws = np.asarray(draws, float)
            null[fam][cname] = {"median": float(np.median(draws)),
                                "q95": float(np.quantile(draws, PERM_Q)),
                                "max_abs": float(np.max(np.abs(draws)))}
            print(f"     {fam:10s} {cname[:26]:26s} null median {np.median(draws):+.3f}  "
                  f"q95 {np.quantile(draws, PERM_Q):+.3f}")
    out["permutation_null"] = null

    # ---- the real run ----------------------------------------------------------------------- #
    for fam in ("full", "clockgenes", "pca"):
        print(f"\n[{fam}]")
        blk = {}
        for cname, yv in Y.items():
            pr, alphas = lodo_predict(X, yv, donor, fam, genes, clock_genes)
            r_all = spearman(pr, yv)
            r_par = partial_spearman(pr, yv, plu, donor)
            q95 = null[fam][cname]["q95"]
            blk[cname] = {"rho_all": r_all, "rho_partial": r_par,
                          "null_q95": q95, "beats_null": bool(r_par > q95),
                          "clears_bar": bool(r_par >= RHO_BAR), "alpha_per_fold": alphas}
            print(f"   {cname[:30]:30s} rho_all {r_all:+.3f}   rho_partial {r_par:+.3f}   "
                  f"null q95 {q95:+.3f}   beats-null {r_par > q95}   bar {r_par >= RHO_BAR}")
        v = verdict({c: b["rho_partial"] for c, b in blk.items()})
        out["families"][fam] = {"per_clock": blk, "verdict": v}
        print(f"   => {fam}: {v}")

    best = max(out["families"], key=lambda f: min(
        b["rho_partial"] for b in out["families"][f]["per_clock"].values()))
    out["verdict"] = out["families"][best]["verdict"]
    out["best_family"] = best
    out["fleischer_baseline"] = {"skin_blood": 0.267, "multi_tissue": 0.516,
                                 "note": "M-2a, same metric and geometry"}
    (_RESULTS / "diag_learned_clock_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print(f"\n  OVERALL: {out['verdict']}  (best family: {best})")
    print("  Fleischer baseline for comparison: 0.267 / 0.516")
    print("  wrote results/diag_learned_clock_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
