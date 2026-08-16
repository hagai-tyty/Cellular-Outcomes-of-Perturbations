"""Does the age representation TRANSFER out of cohort?  (read-only)

    python experiments/diag_age_transfer.py

WHY THIS EXISTS
---------------
`diag_age_capacity` showed the pipeline's panel + ridge reads chronological age at MAE 11.95 yr
(r 0.846, n=133, donors held out) -- matching the published cv_mae of 12.27. But that was measured
ON GSE113957, the cohort the Fleischer clock was FITTED on, so age signal is guaranteed present
there. It establishes CAPACITY, not TRANSFER, and the difference decides what "we just need more
donors" is worth:

  TRANSFERS      -> more donors is the correct and complete answer.
  DOES NOT       -> the representation is cohort-specific, and more donors OF THE SAME KIND
                    will not help. The requirement changes.

THE TEST
--------
Train on GSE113957 (133 normal donors, ages 1-96). Predict the chronological age of DAY-0
FIBROBLASTS from cohorts the model never saw:

  GSE165176 (Gill Sendai)  day-0 `_Fib_` samples, donor ages from the series matrix
  GSE165177 (transient)    `O1/O2/O3 Fib` day-0 samples, ages 38 / 53

Day-0 only. Reprogramming samples are excluded by construction: their expression has been
perturbed, so a miss there would be uninterpretable -- it could mean the representation failed OR
that reprogramming genuinely moved the transcriptome. A resting fibroblast's chronological age is
the only clean out-of-cohort target available.

C-7 APPLIES. `N2_Fib_Sendai_Exp2` is not a transcriptome (library 1.03e8 against the 1e6 that RPM
means; dynamic range 1.74 log2). It is rejected here by the same gate the pipeline uses, not
dropped by hand.

CROSS-COHORT SCALE. Training counts and Gill's Log2-RPM reach the clock through different
pipelines, so two feature treatments are reported and neither is selected after the fact:
  raw      both cohorts CP10k + log1p, nothing further -- transfer including any batch offset
  zscore   features standardised WITHIN each cohort -- removes a constant per-gene batch shift,
           which is the crudest correction that could rescue an otherwise-working representation

PRE-REGISTERED READING (constants below, fixed before running)
  TRANSFERS      MAE <= TRANSFER_MAE_BAR yr AND Spearman >= TRANSFER_RHO_BAR, for a majority of
                 the alpha grid, in at least one feature treatment.
  DOES NOT       otherwise.
The comparison that matters is against predicting the TRAINING cohort's median age, which is what
a model that transfers nothing degrades to. Reported alongside.

n is 8-9 donors, so this can detect a gross failure or a clear success and little in between.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


dre = _load("dre", "experiments/diag_residual_expression.py")
dac = _load("dac", "experiments/diag_age_capacity.py")

from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.integrity import bulk_column_verdict  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402

ALPHAS = (1.0, 10.0, 100.0, 1_000.0, 10_000.0)
TRANSFER_MAE_BAR = 20.0
TRANSFER_RHO_BAR = 0.6

GILL = (r"D:\Gill\GSE165176_Log2_RPM_Sendai_reprogramming (1).txt.gz",
        r"D:\Gill\GSE165176_series_matrix.txt.gz")
# GSE165177 ships its samples across TWO matrices. The day-0 `O1/O2/O3 Fib` columns are in
# part2 only -- the main file has 24 sample columns and no Fib at all, so pointing at it silently
# yields an empty cohort rather than an error.
TRANS = (r"D:\GSE165177\GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz",
         r"D:\GSE165177\GSE165177_series_matrix.txt.gz")


def series_map(path: str) -> tuple[list[str], dict]:
    rows: dict[str, list] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                rows["title"] = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                k = v[0].split(":")[0].strip()
                rows[k] = [x.split(":", 1)[1].strip() if ":" in x else x for x in v]
    age_key = next(k for k in rows if k.lower().startswith("donor age"))
    # GSE165177 titles its day-0 samples `O1_Fib` in the series matrix and `O1 Fib` in the
    # expression header. Keying on the raw title makes the age lookup miss silently and the whole
    # cohort vanish with no error, so both forms are normalised to one.
    norm = {_norm(t): a for t, a in zip(rows["title"], rows[age_key], strict=False)}
    return rows["title"], norm


def _norm(s: str) -> str:
    return re.sub(r"[\s_]+", "_", s.strip()).lower()


def load_bulk_day0(expr_path: str, series_path: str, annot_cols: int) -> pd.DataFrame:
    """Day-0 fibroblast columns only, gated by C-7, as CP10k+log1p over gene symbols."""
    titles, age_of = series_map(series_path)
    df = pd.read_csv(expr_path, sep="\t", low_memory=False)
    sample_cols = list(df.columns[annot_cols:])
    log2 = df[sample_cols].to_numpy(dtype=np.float64)
    keep_cols, ages = [], []
    for j, c in enumerate(sample_cols):
        if not re.search(r"(^|_| )Fib($|_| )", c):
            continue
        ok, _ = bulk_column_verdict(log2[:, j])          # C-7, same gate as the pipeline
        if not ok:
            continue
        a = age_of.get(_norm(c))
        if a is None:
            continue
        keep_cols.append(c)
        ages.append(float(a))
    if not keep_cols:
        return pd.DataFrame()
    idx = [sample_cols.index(c) for c in keep_cols]
    lin = np.power(2.0, log2[:, idx]) - 1.0
    lin[lin < 0] = 0.0
    out = pd.DataFrame(lin, columns=keep_cols)
    out["__sym__"] = df[df.columns[0]].astype(str).to_numpy()
    out["__tot__"] = lin.sum(1)
    out = out.sort_values("__tot__", ascending=False).drop_duplicates("__sym__", keep="first")
    expr = normalize_counts(out[keep_cols].to_numpy(dtype=np.float64).T, target_sum=1e4)
    return pd.DataFrame(expr, index=keep_cols, columns=out.__sym__.tolist()).assign(__age__=ages)


def zscore(a: np.ndarray) -> np.ndarray:
    mu, sd = a.mean(0), a.std(0)
    return (a - mu) / np.where(sd < 1e-12, 1.0, sd)


def main() -> None:
    meta = dac.load_meta()
    tr_expr, tr_genes = dac.load_expression(meta.gsm.tolist())
    keep = meta.disease.eq("Normal").to_numpy() & meta.age.notna().to_numpy()
    tr_expr, tr_age = tr_expr[keep], meta.age.to_numpy(float)[keep]

    cohorts = {}
    for name, (e, s), nannot in (("GSE165176 Gill", GILL, 12), ("GSE165177 transient", TRANS, 12)):
        d = load_bulk_day0(e, s, nannot)
        if not d.empty:
            cohorts[name] = d

    print("=" * 96)
    print(f"OUT-OF-COHORT AGE TRANSFER -- train GSE113957 (n={len(tr_age)}), predict day-0 "
          "fibroblasts elsewhere")
    print(f"  pre-registered: TRANSFERS iff MAE <= {TRANSFER_MAE_BAR} AND spearman >= "
          f"{TRANSFER_RHO_BAR}, majority of alphas")
    print("=" * 96)

    panel = set(GenePanel.load(str(dre.PANEL_PATH)).genes)
    tr_idx = {g: i for i, g in enumerate(tr_genes)}
    res: dict = {"n_train": int(len(tr_age)), "cohorts": {}}

    for name, d in cohorts.items():
        te_age = d["__age__"].to_numpy(float)
        te = d.drop(columns="__age__")
        shared = [g for g in te.columns if g in tr_idx and g in panel]
        Xtr_all = tr_expr[:, [tr_idx[g] for g in shared]]
        Xte_all = te[shared].to_numpy(dtype=np.float64)
        base = float(np.abs(te_age - np.median(tr_age)).mean())
        print(f"\n[{name}]  n={len(te_age)} donors  ages {sorted(te_age.astype(int))}  "
              f"shared panel genes {len(shared)}")
        print(f"  predicting the TRAINING median ({np.median(tr_age):.0f} yr) gives MAE {base:.2f}")
        print(f"  {'features':<9}{'alpha':>8}{'MAE':>9}{'spearman':>10}{'pearson':>9}   pass")
        res["cohorts"][name] = {"n": len(te_age), "ages": te_age.tolist(),
                                "baseline_mae": base, "treatments": {}}
        for treat in ("raw", "zscore"):
            Xtr = zscore(Xtr_all) if treat == "zscore" else Xtr_all
            Xte = zscore(Xte_all) if treat == "zscore" else Xte_all
            n_pass, per = 0, {}
            for a in ALPHAS:
                p = dre.ridge_fit_predict(Xtr, tr_age, Xte, a)
                mae = float(np.abs(p - te_age).mean())
                rho = (float(pd.Series(p).corr(pd.Series(te_age), method="spearman"))
                       if np.std(p) > 0 and np.std(te_age) > 0 else float("nan"))
                pe = (float(np.corrcoef(p, te_age)[0, 1])
                      if np.std(p) > 0 and np.std(te_age) > 0 else float("nan"))
                ok = bool(mae <= TRANSFER_MAE_BAR and np.isfinite(rho) and rho >= TRANSFER_RHO_BAR)
                n_pass += ok
                per[str(a)] = {"mae": mae, "spearman": rho, "pearson": pe, "pass": ok,
                               "pred": p.tolist()}
                print(f"  {treat:<9}{a:>8.0f}{mae:>9.2f}{rho:>10.3f}{pe:>9.3f}   "
                      f"{'YES' if ok else 'no'}")
            verdict = "TRANSFERS" if n_pass > len(ALPHAS) / 2 else "does not"
            print(f"  {'':<9}{'-> ':>8}{verdict}  ({n_pass}/{len(ALPHAS)} alphas)")
            res["cohorts"][name]["treatments"][treat] = {"alphas": per, "n_pass": n_pass,
                                                         "verdict": verdict}

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_age_transfer_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
