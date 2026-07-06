"""
CellFate-Rx — Multi-donor training run (Gill + GSE242423), leave-cell-line-out.

This is the publishable generalization test. It combines:
  * GSE242423 (Kundaje) : 1 single-cell line "HFF"  (D0->iPSC, ~volume)
  * Gill 2022 (GSE165176): 6 bulk donors             (age-valid rejuvenation)
=> 7 cell lines. The `cell_line` regime leaves WHOLE donors out of training, so
the test donor is biologically unseen -- no leakage. This is where a linear clock
(Ridge) should collapse on cross-donor shift and the network should pull ahead.

USAGE (repo root, venv active):
    python run_multi_local.py "D:\\GSE242423" "D:\\Gill"

- arg1 = GSE242423 folder (matrix/barcodes/genes, as before)
- arg2 = Gill folder (must contain the expression matrix + the *series_matrix* file)
  If auto-discovery picks the wrong Gill expr file, set GILL_EXPR/GILL_SERIES below.
"""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Config                                                                       #
# --------------------------------------------------------------------------- #
GSE_DIR = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
GILL_DIR = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"
GILL_EXPR: str | None = None      # set explicitly if auto-discovery is wrong
GILL_SERIES: str | None = None    # set explicitly if auto-discovery is wrong

MAX_CELLS = 5000          # GSE242423 cells per timepoint (volume)
CELLS_PER_RUN = 1000      # densify only this many at a time (bounds RAM)
REGIME = "cell_line"      # leave-cell-line-out -> real cross-donor generalization
N_GENES = 2000
EPOCHS = 80
ENSEMBLE = 5
BATCH = 256
ROOT = "cellfate_multi"

HERE = Path(__file__).resolve().parent
CLOCK = str(HERE / "configs" / "clocks" / "fleischer_clock.json")


def discover_gse(data_dir: str):
    genes = glob.glob(os.path.join(data_dir, "*genes*.tsv.gz"))
    if not genes:
        raise SystemExit(f"no *genes*.tsv.gz in {data_dir}")
    samples = []
    for mtx in sorted(glob.glob(os.path.join(data_dir, "*matrix.mtx.gz"))):
        base = mtx[: -len(".matrix.mtx.gz")]
        bc = base + ".barcodes.tsv.gz"
        label = re.search(r"_([^_]+)$", os.path.basename(base)).group(1)
        if not os.path.exists(bc):
            raise SystemExit(f"missing barcodes for {mtx}")
        samples.append({"matrix": mtx, "barcodes": bc, "label": label})
    if not samples:
        raise SystemExit(f"no *.matrix.mtx.gz in {data_dir}")
    return samples, genes[0]


def discover_gill(data_dir: str):
    series = GILL_SERIES or next(iter(glob.glob(os.path.join(data_dir, "*series_matrix*"))), None)
    if not series:
        raise SystemExit(f"no *series_matrix* file in {data_dir} (need the GEO series matrix)")
    if GILL_EXPR:
        return GILL_EXPR, series
    cands = [p for ext in ("*.txt.gz", "*.tsv.gz", "*.txt", "*.tsv")
             for p in glob.glob(os.path.join(data_dir, ext)) if "series_matrix" not in p]
    if len(cands) != 1:
        raise SystemExit(
            f"could not uniquely identify the Gill expression matrix in {data_dir}.\n"
            f"candidates: {[os.path.basename(c) for c in cands]}\n"
            "-> set GILL_EXPR at the top of this script to the right file.")
    return cands[0], series


def main() -> None:
    import pandas as pd
    import torch

    print(f"[env] torch {torch.__version__} | CUDA: {torch.cuda.is_available()}"
          + (f" | GPU: {torch.cuda.get_device_name(0)}" if torch.cuda.is_available() else ""))
    device = "cuda" if torch.cuda.is_available() else "cpu"

    gse_samples, gse_genes = discover_gse(GSE_DIR)
    gill_expr, gill_series = discover_gill(GILL_DIR)
    print(f"[data] GSE242423 timepoints: {[s['label'] for s in gse_samples]}")
    print(f"[data] Gill expr  : {os.path.basename(gill_expr)}")
    print(f"[data] Gill series: {os.path.basename(gill_series)}")

    from cellfate.common import io
    from cellfate.common.io import ArtifactPaths
    from cellfate.data import DataConfig, QCConfig
    from cellfate.data import run as build_run
    from cellfate.data.sources import GillReprogrammingSource, GSE242423SingleCellSource

    gse = GSE242423SingleCellSource(gse_samples, gse_genes, cell_line="HFF", min_genes=500,
                                    max_cells_per_sample=MAX_CELLS, cells_per_run=CELLS_PER_RUN, seed=0)
    gill = GillReprogrammingSource(gill_expr, gill_series)

    if os.path.isdir(ROOT):
        shutil.rmtree(ROOT)
    print(f"\n[1/5] BUILD combined dataset (regime={REGIME}) — streaming, ~15-25 min ...")
    # GSE242423 first so the gene panel is fit on rich single-cell data.
    build_run(DataConfig(
        out=ROOT, gene_panel=f"{ROOT}/panel.json", n_genes=N_GENES, clock=CLOCK, modality="tf",
        qc=QCConfig(min_genes=500, max_mito_frac=0.20), label_tau=0.7,
        split_fracs=(0.7, 0.1, 0.1, 0.1), split_regimes=(REGIME,), primary_regime=REGIME,
        deconfound=True, seed=0), sources=[gse, gill])

    # ---- composition: which cell line (donor) landed in which split ----
    paths = ArtifactPaths.of(ROOT)
    man = io.manifest_rows(io.load_manifest(paths))
    splits = io.load_splits(paths, REGIME)
    cl_of = {r.cell_id: r.cell_line for r in man}
    n_by = defaultdict(int)
    split_by = {}
    for r in man:
        n_by[r.cell_line] += 1
        split_by[r.cell_line] = splits.get(r.cell_id, "?")

    rows = []
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i in range(len(a["cell_id"])):
            rows.append({"cell_id": a["cell_id"][i], "P_loss": a["y_cls"][i][1],
                         "y_age": a["y_age"][i], "age_mask": bool(a["age_mask"][i])})
    df = pd.DataFrame(rows)
    df["cell_line"] = df.cell_id.map(cl_of)

    print("\n[sanity] per cell line (donor)  ->  split | n | mean P_loss | mean ΔAge")
    for cl in sorted(n_by):
        sub = df[df.cell_line == cl]
        age = sub[sub.age_mask]["y_age"]
        print(f"   {cl:<18} {split_by[cl]:<6} n={n_by[cl]:<6} "
              f"P_loss={sub.P_loss.mean():.3f}  ΔAge={age.mean():+.2f}"
              if len(age) else
              f"   {cl:<18} {split_by[cl]:<6} n={n_by[cl]:<6} P_loss={sub.P_loss.mean():.3f}  ΔAge=n/a")
    heldout = [cl for cl, sp in split_by.items() if sp == "test"]
    print(f"\n   >>> HELD-OUT (test) donor(s): {heldout}  <<<  (model never trains on these)")

    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run
    print(f"\n[2/5] TRAIN {ENSEMBLE}-member ensemble on {device} (leave-cell-line-out) ...")
    train_run(TrainConfig(dataset_dir=ROOT, out=ROOT, regime=REGIME,
                          d_cell=256, d_u=256, latent_dim=256, p_drop=0.2, lr=1e-3,
                          epochs=EPOCHS, patience=10, batch_size=BATCH, ensemble_size=ENSEMBLE,
                          base_seed=0, conformal_levels=(0.90,), device=device))
    m = json.load(open(f"{ROOT}/bundle/metrics.json"))
    print(f"\n[check] n_train={m['n_train']} n_val={m['n_val']} n_calib={m['n_calib']} "
          f"temp={m['temperature']:.3f} conformal_q={m['conformal_q']}")

    from cellfate.evaluation.evaluate_cli import EvalConfig, evaluate
    print("\n[3/5] EVALUATE — leave-cell-line-out gates + baselines ...")
    gates = evaluate(EvalConfig(bundle=ROOT, dataset=ROOT, regimes=(REGIME,), out=f"{ROOT}/reports"))
    print(json.dumps(gates, indent=2, default=str))

    rep = Path(ROOT, "reports", f"{REGIME}.json")
    if rep.exists():
        R = json.loads(rep.read_text())

        def prauc(d):
            v = [x for k, x in d.items() if k.startswith("prauc_")]
            return sum(v) / len(v) if v else float("nan")
        ests = [k for k in R if isinstance(R[k], dict) and any(kk.startswith("prauc_") for kk in R[k])]
        print("\n===== MODEL vs BASELINES (held-out donor) =====")
        print(f"{'estimator':<18}{'PR-AUC':>10}{'ECE':>8}{'reg_MAE':>10}")
        for name in ["model", *[e for e in ests if e != "model"]]:
            if name in R:
                d = R[name]
                print(f"{name:<18}{prauc(d):>10.3f}{d.get('ece', float('nan')):>8.3f}"
                      f"{d.get('reg_mae', float('nan')):>10.2f}")
        print(f"\ncoverage@0.90: {R.get('coverage', 'n/a')}   "
              f"ranking spearman: {R.get('ranking', {}).get('spearman', 'n/a')}")

    print("\n[4/5] SAVE checkpoint ...")
    shutil.make_archive("cellfate_multi_bundle", "zip", f"{ROOT}/bundle")
    print(f"   bundle -> {os.path.abspath('cellfate_multi_bundle.zip')}")
    print("\n[5/5] DONE. This is the leave-cell-line-out result — send back the [sanity]")
    print("composition, the [check] line, the GATES, and the MODEL vs BASELINES table.")


if __name__ == "__main__":
    main()
