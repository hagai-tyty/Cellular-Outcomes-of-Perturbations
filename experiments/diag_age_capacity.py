"""CAN THE REPRESENTATION LEARN AGE AT ALL?  n=143 donors, held out.  (read-only)

    python experiments/diag_age_capacity.py

THE QUESTION THIS SETTLES
-------------------------
Every negative result in this project so far has been reported against 6 reprogramming donors, and
the standing conclusion was "the constraint is n, go get more donors". That conclusion was never
tested, because the reprogramming timecourse and the ABILITY TO READ AGE FROM EXPRESSION are two
different questions, and only the first one is limited to 6.

GSE113957 (Fleischer 2018) has **143 dermal fibroblast samples with declared donor ages 1-96**.
That is 24x the n, and it answers the question directly:

  MAE near the published cv_mae (12.27 yr)  ->  the representation CARRIES age. n was the problem,
                                                and the reprogramming work is a data-collection
                                                question, not a method question.
  MAE near the mean baseline               ->  the representation CANNOT carry age, and no number
                                                of donors fixes it.

NOT CIRCULAR. The target is the **GEO-declared chronological age of the donor**, not a clock
output. `diag_clock_circularity` found the ΔAge regression predicts a linear readout of its own
input; that cannot happen here, because donor age is metadata that no transform of the expression
produced.

ONE HONEST CAVEAT, STATED UP FRONT
----------------------------------
GSE113957 is the dataset the Fleischer clock was FITTED on. That does not make this circular --
the target is the donor's real age -- but it does mean age signal is guaranteed to be present in
this cohort. So a PASS says "the representation can carry age when age signal is there", which is
exactly the capacity question. It does NOT establish out-of-cohort generalisation.

CONFOUNDS HANDLED
  * HGPS (progeria) samples age abnormally fast. Excluded from the primary and reported separately;
    including them would inflate apparent performance.
  * TWO PLATFORMS (GPL18573 n=130, GPL16791 n=13). A batch effect confounded with age would be
    indistinguishable from signal, so per-platform results are reported.
  * Sex casing is inconsistent in GEO ('Male'/'male'); normalised before use.

PRE-REGISTERED READING (constants below, fixed before running)
  CARRIES AGE   cross-validated MAE <= MAE_RATIO_BAR x the mean-baseline MAE, for a majority of the
                alpha grid. Every alpha is reported; none is selected.
  DOES NOT      otherwise.
The published 12.27 yr is a REFERENCE POINT, not a bar -- it came from a different method on the
same data, so beating or missing it is informative rather than decisive.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "dre", ROOT / "experiments" / "diag_residual_expression.py")
dre = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dre)

from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402

DDIR = Path(r"D:\GSE113957")
COUNTS = DDIR / "GSE113957_raw_counts_GRCh38.p13_NCBI.tsv.gz"
ANNOT = DDIR / "Human.GRCh38.p13.annot.tsv.gz"
SERIES = ["GSE113957-GPL18573_series_matrix.txt.gz", "GSE113957-GPL16791_series_matrix.txt.gz"]

ALPHAS = (1.0, 10.0, 100.0, 1_000.0, 10_000.0)
N_SPLITS = 10
SEED = 0
MAE_RATIO_BAR = 0.75
PUBLISHED_CV_MAE = 12.27       # reference only: Fleischer's own CV on this cohort
N_HVG = 2000


def load_meta() -> pd.DataFrame:
    frames = []
    for s in SERIES:
        rows: dict[str, list] = {}
        with gzip.open(DDIR / s, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_geo_accession"):
                    rows["gsm"] = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    v = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                    key = v[0].split(":")[0].strip()
                    rows[key] = [x.split(":", 1)[1].strip() if ":" in x else x for x in v]
        d = pd.DataFrame(rows)
        d["platform"] = s.split("-")[1].split("_")[0]
        frames.append(d)
    m = pd.concat(frames, ignore_index=True)
    m["age"] = pd.to_numeric(m["age"], errors="coerce")
    m["Sex"] = m["Sex"].str.strip().str.lower()
    m["disease"] = m["disease"].str.strip()
    return m


def load_expression(gsms: list[str]) -> tuple[np.ndarray, list[str]]:
    """Samples x genes, CP10k + log1p, gene SYMBOLS, deduped by highest total expression.

    Dedup rule matches the clock's own metadata (`dedup: highest_expressed`), so the gene space is
    the one the rest of this project uses rather than a differently-derived one.
    """
    ann = pd.read_csv(ANNOT, sep="\t", usecols=["GeneID", "Symbol"], low_memory=False)
    sym = dict(zip(ann.GeneID.astype(str), ann.Symbol.astype(str), strict=False))
    df = pd.read_csv(COUNTS, sep="\t", low_memory=False)
    df["__sym__"] = df.GeneID.astype(str).map(sym)
    df = df[df.__sym__.notna() & (df.__sym__ != "nan")]
    counts = df[gsms].to_numpy(dtype=np.float64)
    df["__tot__"] = counts.sum(1)
    keep = df.sort_values("__tot__", ascending=False).drop_duplicates("__sym__", keep="first")
    genes = keep.__sym__.tolist()
    expr = normalize_counts(keep[gsms].to_numpy(dtype=np.float64).T, target_sum=1e4)
    return expr, genes


def cv_predict(X: np.ndarray, y: np.ndarray, alpha: float,
               n_splits: int = N_SPLITS, seed: int = SEED) -> np.ndarray:
    """K-fold with every sample a distinct donor (cell ids are unique), so a plain shuffle split
    already holds donors out."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(y))
    folds = np.array_split(order, n_splits)
    pred = np.empty(len(y))
    for te in folds:
        tr = np.setdiff1d(order, te)
        pred[te] = dre.ridge_fit_predict(X[tr], y[tr], X[te], alpha)
    return pred


def scores(y: np.ndarray, p: np.ndarray) -> dict:
    return {"mae": float(np.abs(p - y).mean()),
            "pearson": float(np.corrcoef(y, p)[0, 1]),
            "spearman": float(pd.Series(y).corr(pd.Series(p), method="spearman"))}


def main() -> None:
    meta = load_meta()
    expr, genes = load_expression(meta.gsm.tolist())
    print("=" * 100)
    print("AGE CAPACITY TEST -- can the representation learn chronological age?  (n donors held out)")
    print(f"  target: GEO-declared donor age (NOT a clock output).  reference cv_mae "
          f"{PUBLISHED_CV_MAE} yr")
    print(f"  pre-registered: CARRIES AGE iff CV MAE <= {MAE_RATIO_BAR} x mean-baseline MAE for a "
          "majority of alphas")
    print("=" * 100)
    print(f"  loaded {expr.shape[0]} samples x {expr.shape[1]} genes")
    print(f"  disease: {meta.disease.value_counts().to_dict()}")
    print(f"  platform: {meta.platform.value_counts().to_dict()}")
    print(f"  age: min {meta.age.min():.0f} max {meta.age.max():.0f} "
          f"mean {meta.age.mean():.1f} sd {meta.age.std():.1f}")

    idx_of = {g: i for i, g in enumerate(genes)}
    panel = [g for g in GenePanel.load(str(dre.PANEL_PATH)).genes if g in idx_of]
    feats = {"pipeline panel": np.array([idx_of[g] for g in panel])}

    res: dict = {"n_total": int(expr.shape[0]), "runs": {}}
    for label, sel in (("NORMAL only (primary)", meta.disease.eq("Normal").to_numpy()),
                       ("HGPS only", meta.disease.eq("HGPS").to_numpy()),
                       ("GPL18573 only", meta.platform.eq("GPL18573").to_numpy()
                        & meta.disease.eq("Normal").to_numpy())):
        m = sel & meta.age.notna().to_numpy()
        if m.sum() < 30:
            print(f"\n[{label}] n={int(m.sum())} -- too few for {N_SPLITS}-fold, skipped")
            continue
        y = meta.age.to_numpy(float)[m]
        base = float(np.abs(y - np.median(y)).mean())
        # HVGs are refit within each cohort so the feature set is not chosen using other cohorts
        sub = expr[m]
        hv = np.argsort(-sub.std(0))[:N_HVG]
        runs = dict(feats, **{f"top-{N_HVG} HVG": hv})
        print(f"\n[{label}]  n={int(m.sum())}   mean-baseline MAE {base:.2f} yr "
              f"(predict the median age)")
        print(f"  {'features':<18}{'alpha':>8}{'MAE':>9}{'ratio':>8}{'pearson':>9}"
              f"{'spearman':>10}   pass")
        res["runs"][label] = {"n": int(m.sum()), "baseline_mae": base, "features": {}}
        for fname, cols in runs.items():
            X = sub[:, cols]
            n_pass, per = 0, {}
            for a in ALPHAS:
                s = scores(y, cv_predict(X, y, a))
                ratio = s["mae"] / base
                ok = ratio <= MAE_RATIO_BAR
                n_pass += ok
                per[str(a)] = {**s, "ratio": ratio, "pass": ok}
                print(f"  {fname:<18}{a:>8.0f}{s['mae']:>9.2f}{ratio:>8.2f}"
                      f"{s['pearson']:>9.3f}{s['spearman']:>10.3f}   {'YES' if ok else 'no'}")
            verdict = "CARRIES AGE" if n_pass > len(ALPHAS) / 2 else "does not"
            print(f"  {'':<18}{'-> ':>8}{verdict}  ({n_pass}/{len(ALPHAS)} alphas)")
            res["runs"][label]["features"][fname] = {"alphas": per, "n_pass": n_pass,
                                                     "verdict": verdict}

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_age_capacity_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
