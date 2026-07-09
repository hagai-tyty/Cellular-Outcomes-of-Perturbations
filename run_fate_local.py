"""
CellFate-Rx — Fate-classification test on a held-out HFF slice.

The donor-holdout run left fate PR-AUC undefined because the held-out Gill donor
(21 bulk samples) has no fate diversity. THIS run instead holds out a random slice
of the HFF single-cell line (which spans D0->iPSC, so it contains safe / loss /
death), giving a MEASURABLE PR-AUC. It answers: "does the fate head classify?"

It is a held-out-CELLS test on one line, NOT a cross-donor test. If fate PR-AUC is
strong here, the fate head works and we move on to the ΔAge harmonisation; if it is
weak even here, the fate head itself needs debugging.

USAGE (repo root, venv active):
    python run_fate_local.py "D:\\GSE242423" "D:\\Gill"
"""

from __future__ import annotations

import glob
import json
import math
import os
import re
import shutil
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

GSE_DIR = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
GILL_DIR = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"
GILL_EXPR: str | None = None
GILL_SERIES: str | None = None

MAX_CELLS = 5000
CELLS_PER_RUN = 1000
FATE_TEST_LINE = "HFF"     # hold out a slice of this (fate-diverse) single-cell line
FATE_TEST_FRAC = 0.15      # 15% of HFF -> test
N_GENES = 2000
EPOCHS = 80
ENSEMBLE = 5
BATCH = 256
ROOT = "cellfate_fate"

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
        raise SystemExit(f"no *series_matrix* file in {data_dir}")
    if GILL_EXPR:
        return GILL_EXPR, series
    cands = [p for ext in ("*.txt.gz", "*.tsv.gz", "*.txt", "*.tsv")
             for p in glob.glob(os.path.join(data_dir, ext)) if "series_matrix" not in p]
    if len(cands) != 1:
        raise SystemExit(f"could not identify the Gill expression matrix; candidates: "
                         f"{[os.path.basename(c) for c in cands]} -- set GILL_EXPR.")
    return cands[0], series


def main() -> None:
    import numpy as np
    import pandas as pd
    import torch

    print(f"[env] torch {torch.__version__} | CUDA: {torch.cuda.is_available()}"
          + (f" | GPU: {torch.cuda.get_device_name(0)}" if torch.cuda.is_available() else ""))
    device = "cuda" if torch.cuda.is_available() else "cpu"

    gse_samples, gse_genes = discover_gse(GSE_DIR)
    gill_expr, gill_series = discover_gill(GILL_DIR)
    print(f"[data] GSE242423 timepoints: {[s['label'] for s in gse_samples]}")
    print(f"[data] Gill expr: {os.path.basename(gill_expr)}")
    print(f"[data] fate test: hold out {FATE_TEST_FRAC:.0%} of '{FATE_TEST_LINE}' (fate-diverse) as TEST")

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
    print(f"\n[1/5] BUILD combined dataset (regime=line_holdout, test={FATE_TEST_LINE}) — streaming ...")
    build_run(DataConfig(
        out=ROOT, gene_panel=f"{ROOT}/panel.json", n_genes=N_GENES, clock=CLOCK, modality="tf",
        qc=QCConfig(min_genes=500, max_mito_frac=0.20), label_tau=0.7,
        split_fracs=(0.8, 0.1, 0.1, 0.0), split_regimes=("line_holdout",), primary_regime="line_holdout",
        fate_test_line=FATE_TEST_LINE, fate_test_frac=FATE_TEST_FRAC,
        deconfound=True, seed=0), sources=[gse, gill])

    paths = ArtifactPaths.of(ROOT)
    man = io.manifest_rows(io.load_manifest(paths))
    splits = io.load_splits(paths, "line_holdout")
    per_cl_split: dict[str, Counter] = defaultdict(Counter)
    for r in man:
        per_cl_split[r.cell_line][splits.get(r.cell_id, "?")] += 1

    # fate-class composition of the TEST set (this is the whole point)
    rows = []
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i in range(len(a["cell_id"])):
            rows.append({"cell_id": a["cell_id"][i], "cls": int(np.argmax(a["y_cls"][i])),
                         "split": splits.get(a["cell_id"][i], "?")})
    df = pd.DataFrame(rows)
    test_cls = Counter(df[df.split == "test"]["cls"])
    names = {0: "safe", 1: "loss", 2: "death"}
    print("\n[sanity] per cell line -> split distribution")
    for cl in sorted(per_cl_split):
        dist = " ".join(f"{k}={v}" for k, v in sorted(per_cl_split[cl].items()))
        print(f"   {cl:<8} [{dist}]")
    print(f"\n[sanity] TEST-set fate classes: "
          f"{ {names[k]: int(v) for k, v in sorted(test_cls.items())} }")
    if len([c for c in test_cls if test_cls[c] > 0]) < 2:
        print("   [WARN] test set has <2 fate classes -> PR-AUC may still be undefined.")

    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run
    print(f"\n[2/5] TRAIN {ENSEMBLE}-member ensemble on {device} ...")
    train_run(TrainConfig(dataset_dir=ROOT, out=ROOT, regime="line_holdout",
                          d_cell=256, d_u=256, latent_dim=256, p_drop=0.2, lr=1e-3,
                          epochs=EPOCHS, patience=10, batch_size=BATCH, ensemble_size=ENSEMBLE,
                          base_seed=0, conformal_levels=(0.90,), device=device))
    m = json.load(open(f"{ROOT}/bundle/metrics.json"))
    print(f"\n[check] n_train={m['n_train']} n_val={m['n_val']} n_calib={m['n_calib']} "
          f"temp={m['temperature']:.3f}")

    from cellfate.evaluation.evaluate_cli import EvalConfig, evaluate
    print("\n[3/5] EVALUATE — fate classification on held-out HFF ...")
    gates = evaluate(EvalConfig(bundle=ROOT, dataset=ROOT, regimes=("line_holdout",), out=f"{ROOT}/reports"))
    print(json.dumps(gates, indent=2, default=str))

    rep = Path(ROOT, "reports", "line_holdout.json")
    if rep.exists():
        R = json.loads(rep.read_text())

        def prauc(d):
            v = [x for k, x in d.items() if k.startswith("prauc_")]
            return sum(v) / len(v) if v else float("nan")

        def perclass(d, kind):
            return {names[c]: d.get(f"{kind}_{c}", float("nan")) for c in range(3)}
        ests = [k for k in R if isinstance(R[k], dict) and any(kk.startswith("prauc_") for kk in R[k])]
        mpr = prauc(R.get("model", {}))
        print("\n===== FATE: per-class PR-AUC (held-out HFF) =====")
        if "model" in R:
            print(f"   model PR-AUC by class: "
                  f"{ {k: round(v, 3) for k, v in perclass(R['model'], 'prauc').items()} }")

        def verdict(mv, bv, lower_better):
            if not (math.isfinite(mv) and math.isfinite(bv)):
                return "  -  "
            return " WIN " if ((mv < bv) if lower_better else (mv > bv)) else " loss"
        print("\n===== MODEL vs BASELINES (fate) — per-metric win/loss =====")
        print(f"{'estimator':<16}{'PR-AUC':>9}{'beatsPR':>9}{'ECE':>8}{'beatsECE':>10}")
        me = R.get("model", {}).get("ece", float("nan"))
        print(f"{'model':<16}{mpr:>9.3f}{'  -  ':>9}{me:>8.3f}{'  -  ':>10}")
        beats_pr = math.isfinite(mpr)
        for name in [e for e in ests if e != "model"]:
            d = R[name]
            bp, be = prauc(d), d.get("ece", float("nan"))
            vp, ve = verdict(mpr, bp, False), verdict(me, be, True)
            if vp.strip() == "loss":
                beats_pr = False
            print(f"{name:<16}{bp:>9.3f}{vp:>9}{be:>8.3f}{ve:>10}")
        print(f"\nmodel beats ALL baselines on fate PR-AUC: {beats_pr}")

    print("\n[4/5] SAVE ...")
    shutil.make_archive("cellfate_fate_bundle", "zip", f"{ROOT}/bundle")
    print(f"   bundle -> {os.path.abspath('cellfate_fate_bundle.zip')}")
    print("\n[5/5] DONE. Send back the [sanity] TEST-set fate classes, the [check] line,")
    print("the per-class PR-AUC, and the win/loss table.")


if __name__ == "__main__":
    main()
