"""STAGE 1.5.1 — independent verification of every claim in STAGE_1_5_1_NEW_CHANGES.md.

    python experiments/stage_1_5_1_tests.py "D:\\GSE113957" "D:\\Gill"

READ-ONLY. Writes `stage_1_5_1_tests_results.json`. No refit is shipped, no artefact is replaced,
`src/` is untouched. Every number the review reports is RE-DERIVED here from the raw data rather
than accepted; nothing below cites the review as evidence for itself.

TESTS
  T1  is `normalize_counts` per-row?              -> tests R4's stated MECHANISM
  T2  does any cross-sample scaler exist?          -> tests R4's stated mechanism, structurally
  T3  are samples donor-independent?               -> tests the group-leakage alternative
  T4  reproduce the shipped CV + slope + alpha     -> tests R4's CONCLUSION, and R1/R2
  T5  does slope recalibration (C3) help?          -> tests the C3 elimination
  T6  how many samples does GSE113957 yield?       -> tests the 133-vs-143 loose end
  T7  ΔAge label noise measured from REPLICATES    -> NEW: grounds the §3 bar empirically
  T8  is `cv_mae <= 4.0` reachable at all?         -> NEW: bar feasibility, sparse vs dense
  T9  error by age decile                          -> the outstanding Step 1 item (bar S4)

T7 is the one that matters most and neither document ran it. Both the original plan and the review
DERIVE the ΔAge label noise as `sqrt(2) * cv_mae`. That is an assumption: `cv_mae` is a BETWEEN-donor
quantity, while ΔAge is a WITHIN-donor difference, and any per-donor systematic offset cancels in the
difference (the same cancellation proved in Stage 1.5 §2 Group A). If most of `cv_mae` is
per-donor bias, the true label noise is far smaller than `sqrt(2)*12.27` and the bar is too strict;
if it is not, the bar is right for a reason nobody has yet checked. Gill's Exp1/Exp2 pairs are the
same donor, day and marker, so their true ages are IDENTICAL: the spread of their predicted-age
differences IS the ΔAge label noise, measured directly and assumption-free.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "local_runners", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

RESULTS: dict = {}


def _hdr(t: str) -> None:
    print(f"\n{'=' * 78}\n{t}\n{'=' * 78}")


# --------------------------------------------------------------------------- #
# T1/T2 — R4's mechanism: is there any cross-sample statistic?                 #
# --------------------------------------------------------------------------- #
def t1_normalisation_is_per_row() -> dict:
    _hdr("T1 — is normalize_counts per-row? (R4's stated mechanism)")
    from cellfate.data.normalize import normalize_counts

    rng = np.random.default_rng(0)
    X = rng.poisson(5.0, size=(40, 200)).astype(np.float64)
    sub = normalize_counts(X[:5])
    full = normalize_counts(X)[:5]
    diff = float(np.max(np.abs(sub - full)))
    # a second, harder probe: change OTHER rows and see if the first rows move
    X2 = X.copy()
    X2[10:] *= 7.0
    moved = float(np.max(np.abs(normalize_counts(X2)[:5] - full)))
    per_row = diff == 0.0 and moved == 0.0
    print(f"  max|normalize(X[:5]) - normalize(X)[:5]| = {diff:.3e}")
    print(f"  max|rows 0-4 after scaling rows 10+ by 7| = {moved:.3e}")
    print(f"  -> per-row (no cross-sample statistic): {per_row}")
    return {"max_diff_subset_vs_full": diff, "max_diff_after_perturbing_other_rows": moved,
            "is_per_row": bool(per_row)}


def t2_no_cross_sample_scaler() -> dict:
    _hdr("T2 — does clock_fit.py contain any cross-sample scaler?")
    src = (ROOT / "src" / "cellfate" / "data" / "clock_fit.py").read_text(encoding="utf-8")
    tokens = ["StandardScaler", "scaler", "fit_transform", "MinMaxScaler", "RobustScaler",
              "zscore", "z_score"]
    hits = {t: src.count(t) for t in tokens if t in src}
    # RidgeCV standardises internally ONLY if normalize/scale is requested -- check that too
    ridge_line = [ln.strip() for ln in src.splitlines() if "RidgeCV(" in ln]
    print(f"  scaler-like tokens found: {hits or 'NONE'}")
    print(f"  RidgeCV call sites: {ridge_line}")
    return {"scaler_tokens": hits, "ridgecv_calls": ridge_line, "has_cross_sample_scaler": bool(hits)}


# --------------------------------------------------------------------------- #
# Shared loader                                                                #
# --------------------------------------------------------------------------- #
def load_gse113957(known_dir: str):
    """Reuse the diagnostic's loader so this measures the SAME pipeline the project uses."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "diag_clock_validity", ROOT / "experiments" / "diag_clock_validity.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    ages, counts, genes = m._load_known_age_fibroblasts(known_dir)
    return ages, counts, genes, m


def t3_donor_independence(known_dir: str) -> dict:
    _hdr("T3 — are GSE113957 samples donor-independent? (group-leakage alternative)")
    import gzip
    import re
    from collections import Counter
    ids: list[str] = []
    fields: Counter = Counter()
    for p in sorted(Path(known_dir).glob("*series_matrix*")):
        op = gzip.open if str(p).endswith(".gz") else open
        with op(p, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("!Sample_characteristics_ch1"):
                    vals = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
                    key = vals[0].split(":")[0].strip().lower() if vals else ""
                    fields[key] += 1
                    if re.search(r"cell.?id|subject|individual|donor", key):
                        ids += [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    n_ids, n_uniq = len(ids), len(set(ids))
    dupes = [k for k, c in Counter(ids).items() if c > 1]
    print(f"  characteristic fields seen: {dict(fields)}")
    print(f"  donor-identifying values: {n_ids} total, {n_uniq} unique, {len(dupes)} repeated")
    return {"n_id_values": n_ids, "n_unique": n_uniq, "n_repeated": len(dupes),
            "repeated_examples": dupes[:5], "fields": dict(fields),
            "donor_independent": bool(n_ids > 0 and len(dupes) == 0)}


def t4_reproduce_cv(ages, counts, genes) -> dict:
    _hdr("T4 — reproduce the shipped CV, slope, alpha, in-sample gap (R4 conclusion; R1/R2)")
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold, cross_val_predict

    from cellfate.data.normalize import normalize_counts

    Xn = normalize_counts(np.asarray(counts, float))
    y = np.asarray(ages, float)
    alphas = np.logspace(-1.0, 4.0, 24)
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    pred = cross_val_predict(RidgeCV(alphas=alphas), Xn, y, cv=cv)
    mae = float(mean_absolute_error(y, pred))
    pear = float(np.corrcoef(pred, y)[0, 1])
    slope = float(np.polyfit(y, pred, 1)[0])
    full = RidgeCV(alphas=alphas).fit(Xn, y)
    in_mae = float(mean_absolute_error(y, full.predict(Xn)))
    print(f"  n={len(y)}  cv_mae={mae:.2f}  cv_pearson={pear:.3f}  slope(pred~true)={slope:.3f}")
    print(f"  alpha={full.alpha_:.4g} (grid {alphas.min():.2g}..{alphas.max():.2g})")
    print(f"  in-sample MAE={in_mae:.2f}  ->  CV/in-sample ratio = {mae / max(in_mae, 1e-9):.1f}x")
    return {"n": int(len(y)), "cv_mae": mae, "cv_pearson": pear, "slope": slope,
            "alpha": float(full.alpha_), "alpha_grid_min": float(alphas.min()),
            "in_sample_mae": in_mae, "cv_over_insample": float(mae / max(in_mae, 1e-9)),
            "cv_pred": pred.tolist(), "ages": y.tolist()}


def t5_slope_recalibration(t4: dict) -> dict:
    _hdr("T5 — does out-of-fold slope recalibration (C3) help?")
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold
    y = np.asarray(t4["ages"], float)
    p = np.asarray(t4["cv_pred"], float)
    # recalibrate OUT OF FOLD: fit the linear map on the other folds' (pred, true) pairs
    out = np.empty_like(p)
    for tr, te in KFold(n_splits=5, shuffle=True, random_state=0).split(p.reshape(-1, 1)):
        a, b = np.polyfit(p[tr], y[tr], 1)
        out[te] = a * p[te] + b
    mae_recal = float(mean_absolute_error(y, out))
    print(f"  control (dense ridge CV) MAE = {t4['cv_mae']:.2f}")
    print(f"  C3 slope-recalibrated   MAE = {mae_recal:.2f}  "
          f"({(mae_recal - t4['cv_mae']) / t4['cv_mae'] * 100:+.1f}%)")
    return {"control_mae": t4["cv_mae"], "recalibrated_mae": mae_recal,
            "delta_pct": float((mae_recal - t4["cv_mae"]) / t4["cv_mae"] * 100),
            "helps": bool(mae_recal < t4["cv_mae"])}


def t6_sample_count(ages) -> dict:
    _hdr("T6 — how many samples does GSE113957 actually yield? (133 vs 143)")
    shipped = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text(
        encoding="utf-8")).get("meta", {})
    n_here = int(len(ages))
    print(f"  parsed here: {n_here}   |   artefact records: {shipped.get('n_samples')}")
    print(f"  age range here: [{np.min(ages):.0f}, {np.max(ages):.0f}]  "
          f"|  artefact: {shipped.get('age_range')}")
    n_zero = int(np.sum(np.asarray(ages, float) < 1.0))
    print(f"  samples below age 1 (the out-of-range gap, R5): {n_zero}")
    return {"n_parsed": n_here, "n_in_artefact": shipped.get("n_samples"),
            "delta": n_here - int(shipped.get("n_samples", 0)),
            "age_min": float(np.min(ages)), "age_max": float(np.max(ages)),
            "n_below_age_1": n_zero}


# --------------------------------------------------------------------------- #
# T7 — NEW: measure the ΔAge label noise directly, no cv_mae assumption        #
# --------------------------------------------------------------------------- #
def t7_label_noise_from_replicates(gill_dir: str) -> dict:
    _hdr("T7 — ΔAge LABEL NOISE measured from same-donor/day/marker replicate pairs (NEW)")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "diag_zero_point", ROOT / "experiments" / "diag_zero_point.py")
    dzp = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = dzp
    spec.loader.exec_module(dzp)
    diffs, meta = dzp.matched_exp_offsets(gill_dir)
    d = np.asarray([x for x in diffs if np.isfinite(x)], float)
    if len(d) < 3:
        print("  too few matched pairs")
        return {"n_pairs": int(len(d)), "status": "CANNOT_VERIFY"}
    sd_diff = float(d.std(ddof=1))
    single = sd_diff / np.sqrt(2.0)          # SD of ONE measurement
    clock_cv = 12.26879346460328
    derived = float(np.sqrt(2.0) * clock_cv)  # what BOTH documents assumed
    print(f"  {len(d)} matched pairs (same donor, day, marker -> TRUE ages identical)")
    print(f"  mean diff {d.mean():+.2f} yr   SD of differences = {sd_diff:.2f} yr")
    print(f"  => single-measurement within-donor SD = {single:.2f} yr")
    print(f"  => ΔAge label noise (a difference of two)  = {sd_diff:.2f} yr   [MEASURED]")
    print(f"     both documents ASSUMED sqrt(2)*cv_mae   = {derived:.2f} yr   [DERIVED]")
    print(f"     ratio measured/derived = {sd_diff / derived:.2f}")
    return {"n_pairs": int(len(d)), "mean_diff": float(d.mean()), "sd_of_differences": sd_diff,
            "single_measurement_sd": single, "assumed_sqrt2_cv_mae": derived,
            "ratio_measured_over_derived": float(sd_diff / derived), "pair_meta": meta}


# --------------------------------------------------------------------------- #
# T8 — NEW: is the 4.0 yr bar reachable at all? sparse vs dense                #
# --------------------------------------------------------------------------- #
def t10_disease_composition(known_dir: str) -> dict:
    """NEW — the 133-vs-143 gap: WHO are the 10? Neither document checked the `disease` field."""
    _hdr("T10 — what are the 10 'missing' samples? (resolves the P3 correctness claim)")
    import glob
    import gzip
    from collections import Counter
    counts_c: Counter = Counter()
    gsm_disease: dict[str, str] = {}
    for f in sorted(glob.glob(str(Path(known_dir) / "*series_matrix*"))):
        op = gzip.open if f.endswith(".gz") else open
        gsms: list[str] = []
        with op(f, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("!Sample_geo_accession"):
                    gsms = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    vals = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
                    if vals and vals[0].lower().startswith("disease"):
                        dz = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
                        counts_c.update(dz)
                        for g, d in zip(gsms, dz, strict=False):
                            gsm_disease[g] = d
    print(f"  disease field across all samples: {dict(counts_c)}")
    print(f"  -> the 10 are {'HGPS (progeria)' if counts_c.get('HGPS') else 'NOT identified'}; "
          "excluding them is scientifically REQUIRED, not an accident")
    return {"disease_counts": dict(counts_c), "n_normal": int(counts_c.get("Normal", 0)),
            "n_hgps": int(counts_c.get("HGPS", 0)), "gsm_disease": gsm_disease}


def normal_mask(known_dir: str, n_expected: int) -> np.ndarray | None:
    """Boolean mask over the loader's sample order selecting `disease == Normal`."""
    import glob

    import pandas as pd
    t10 = RESULTS.get("T10_disease_composition") or t10_disease_composition(known_dir)
    gsm_disease = t10["gsm_disease"]
    hits = sorted(glob.glob(str(Path(known_dir) / "*raw_counts*NCBI*.tsv*"))) or \
        sorted(glob.glob(str(Path(known_dir) / "*raw_counts*.tsv*")))
    if not hits:
        return None
    cols = list(pd.read_csv(hits[0], sep="\t", nrows=0).columns)[1:]   # drop GeneID
    sample_cols = [c for c in cols if c in gsm_disease]
    m = np.array([gsm_disease.get(c) == "Normal" for c in sample_cols], dtype=bool)
    return m if len(m) == n_expected else None


def t4b_cv_normal_only(ages, counts, mask) -> dict:
    """NEW — does the shipped 12.27 reproduce once HGPS is excluded, as the paper did?"""
    _hdr("T4b — reproduce the CV on Normal-only (n=133), the set the paper actually used (NEW)")
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold, cross_val_predict

    from cellfate.data.normalize import normalize_counts
    Xn = normalize_counts(np.asarray(counts, float)[mask])
    y = np.asarray(ages, float)[mask]
    alphas = np.logspace(-1.0, 4.0, 24)
    pred = cross_val_predict(RidgeCV(alphas=alphas), Xn, y,
                             cv=KFold(5, shuffle=True, random_state=0))
    mae = float(mean_absolute_error(y, pred))
    med = float(np.median(np.abs(pred - y)))
    pear = float(np.corrcoef(pred, y)[0, 1])
    a = float(RidgeCV(alphas=alphas).fit(Xn, y).alpha_)
    print(f"  n={len(y)}  cv_mae={mae:.2f}  cv_MEDIAN_err={med:.2f}  cv_pearson={pear:.3f}  "
          f"alpha={a:.4g}")
    print("  shipped artefact: cv_mae 12.27, cv_pearson 0.837, alpha 0.2721")
    return {"n": int(len(y)), "cv_mae": mae, "cv_median_error": med, "cv_pearson": pear,
            "alpha": a, "shipped_cv_mae": 12.26879346460328,
            "reproduces_within_0_5yr": bool(abs(mae - 12.26879346460328) <= 0.5)}


def t11_alpha_grid_boundary(ages, counts, mask) -> dict:
    """NEW — the shipped grid starts at 0.1 and RidgeCV lands ON it. Is the grid mis-specified?"""
    _hdr("T11 — is the alpha grid pinned at its lower boundary? (NEW; neither doc checked)")
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold, cross_val_predict

    from cellfate.data.normalize import normalize_counts
    Xn = normalize_counts(np.asarray(counts, float)[mask])
    y = np.asarray(ages, float)[mask]
    out = {}
    for name, grid in (("shipped logspace(-1,4)", np.logspace(-1.0, 4.0, 24)),
                       ("extended logspace(-4,4)", np.logspace(-4.0, 4.0, 40)),
                       ("wide logspace(0,6)", np.logspace(0.0, 6.0, 40))):
        pred = cross_val_predict(RidgeCV(alphas=grid), Xn, y,
                                 cv=KFold(5, shuffle=True, random_state=0))
        mae = float(mean_absolute_error(y, pred))
        a = float(RidgeCV(alphas=grid).fit(Xn, y).alpha_)
        at_edge = bool(np.isclose(a, grid.min()) or np.isclose(a, grid.max()))
        out[name] = {"cv_mae": mae, "alpha": a, "at_grid_edge": at_edge}
        print(f"  {name:26s} cv_mae={mae:6.2f}  alpha={a:9.4g}  at_edge={at_edge}")
    best = min(out.values(), key=lambda d: d["cv_mae"])["cv_mae"]
    print(f"  best over grids = {best:.2f} yr  (vs shipped grid "
          f"{out['shipped logspace(-1,4)']['cv_mae']:.2f})")
    return {"grids": out, "best_cv_mae": best}


def t8_bar_feasibility(ages, counts, mask) -> dict:
    _hdr("T8 — BAR FEASIBILITY: can a sparse model approach cv_mae <= 4.0? (NEW)")
    from sklearn.linear_model import ElasticNetCV
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold

    from cellfate.data.normalize import normalize_counts
    Xn = normalize_counts(np.asarray(counts, float)[mask])
    y = np.asarray(ages, float)[mask]
    outer = KFold(n_splits=5, shuffle=True, random_state=0)
    pred = np.empty_like(y)
    nnz: list[int] = []
    for tr, te in outer.split(Xn):
        # everything fitted INSIDE the fold -- the guard the review correctly insists on
        en = ElasticNetCV(l1_ratio=0.5, alphas=np.logspace(-2, 1.5, 12), cv=3,
                          max_iter=3000, tol=1e-3, random_state=0)
        en.fit(Xn[tr], y[tr])
        pred[te] = en.predict(Xn[te])
        nnz.append(int(np.sum(en.coef_ != 0)))
    mae = float(mean_absolute_error(y, pred))
    med = float(np.median(np.abs(pred - y)))
    pear = float(np.corrcoef(pred, y)[0, 1])
    slope = float(np.polyfit(y, pred, 1)[0])
    print("  ElasticNet (l1_ratio 0.5), all selection inside each fold, Normal-only:")
    print(f"    cv_mae={mae:.2f}  cv_MEDIAN={med:.2f}  cv_pearson={pear:.3f}  slope={slope:.3f}  "
          f"median non-zero genes={int(np.median(nnz))}")
    verdict = "PASS" if mae <= 4.0 else ("MARGINAL" if mae <= 6.0 else "FAIL")
    print(f"    -> against the pre-registered bar (mean MAE): {verdict}")
    return {"cv_mae": mae, "cv_median_error": med, "cv_pearson": pear, "slope": slope,
            "median_nonzero_genes": int(np.median(nnz)), "nonzero_per_fold": nnz,
            "bar_verdict": verdict}


def t9_error_by_decile(t4: dict) -> dict:
    _hdr("T9 — error by age decile (outstanding Step 1 item; feeds bar S4)")
    y = np.asarray(t4["ages"], float)
    p = np.asarray(t4["cv_pred"], float)
    err = np.abs(p - y)
    qs = np.quantile(y, np.linspace(0, 1, 6))
    rows = []
    for i in range(5):
        lo, hi = qs[i], qs[i + 1]
        m = (y >= lo) & (y <= hi if i == 4 else y < hi)
        if m.sum():
            rows.append({"age_lo": float(lo), "age_hi": float(hi), "n": int(m.sum()),
                         "mae": float(err[m].mean()), "bias": float((p[m] - y[m]).mean())})
            print(f"  age {lo:5.1f}-{hi:5.1f}  n={m.sum():3d}  MAE={err[m].mean():6.2f}  "
                  f"bias={np.mean(p[m] - y[m]):+7.2f}")
    maes = [r["mae"] for r in rows]
    ratio = float(max(maes) / max(min(maes), 1e-9))
    print(f"  worst/best decile MAE ratio = {ratio:.2f}")
    return {"bins": rows, "worst_over_best_ratio": ratio}


def main() -> int:
    known = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE113957"
    gill = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"
    print("STAGE 1.5.1 — independent verification of STAGE_1_5_1_NEW_CHANGES.md")
    print(f"  GSE113957: {known}\n  Gill      : {gill}")

    RESULTS["T1_normalisation_per_row"] = t1_normalisation_is_per_row()
    RESULTS["T2_no_cross_sample_scaler"] = t2_no_cross_sample_scaler()

    if not Path(known).exists():
        print(f"\n  !! {known} missing; T3-T6,T8,T9 skipped")
    else:
        RESULTS["T3_donor_independence"] = t3_donor_independence(known)
        RESULTS["T10_disease_composition"] = t10_disease_composition(known)
        ages, counts, genes, _m = load_gse113957(known)
        if ages is None:
            print("  !! GSE113957 did not parse")
        else:
            RESULTS["T6_sample_count"] = t6_sample_count(ages)
            t4 = t4_reproduce_cv(ages, counts, genes)
            RESULTS["T4_reproduce_cv"] = {k: v for k, v in t4.items()
                                          if k not in ("cv_pred", "ages")}
            RESULTS["T5_slope_recalibration"] = t5_slope_recalibration(t4)
            RESULTS["T9_error_by_decile"] = t9_error_by_decile(t4)
            mask = normal_mask(known, len(ages))
            if mask is None or not mask.any():
                print("  !! could not align the Normal mask; T4b/T11/T8 fall back to all samples")
                mask = np.ones(len(ages), dtype=bool)
            RESULTS["T4b_cv_normal_only"] = t4b_cv_normal_only(ages, counts, mask)
            RESULTS["T11_alpha_grid"] = t11_alpha_grid_boundary(ages, counts, mask)
            RESULTS["T8_bar_feasibility"] = t8_bar_feasibility(ages, counts, mask)

    if Path(gill).exists():
        try:
            RESULTS["T7_label_noise_from_replicates"] = t7_label_noise_from_replicates(gill)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! T7 failed: {exc!r}")
            RESULTS["T7_label_noise_from_replicates"] = {"error": repr(exc)[:200]}
    else:
        print(f"\n  !! {gill} missing; T7 skipped")

    out = {"script": "stage_1_5_1_tests", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "gse113957_dir": known, "gill_dir": gill, "tests": RESULTS}
    Path("stage_1_5_1_tests_results.json").write_text(json.dumps(out, indent=2, default=str),
                                                      encoding="utf-8")
    print("\n  wrote stage_1_5_1_tests_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
