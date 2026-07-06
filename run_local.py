"""
CellFate-Rx — Stage 3 local training run (GPU).

Runs the whole pipeline on your machine: build the GSE242423 dataset -> sanity-check ->
train the calibrated ensemble -> evaluate -> save the checkpoint.

USAGE (from the repo root, with your env active):
    python run_local.py "D:\\GSE242423"

- The single argument is the folder holding the GSE242423 10x files
  (the *.matrix.mtx.gz + *.barcodes.tsv.gz per timepoint, plus one *genes*.tsv.gz).
- Edit MAX_CELLS below if you hit GPU/RAM limits (lower it) or want more volume (raise it).
"""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Config — the two values that broke the last run are set correctly here.      #
# --------------------------------------------------------------------------- #
DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
MAX_CELLS = 1500          # per timepoint. NOT 5. Lower to 3000/2000 if GPU/RAM is tight.
REGIME = "random"         # cell-level split -> fills train/val/calib/test on 1 donor
N_GENES = 2000
EPOCHS = 80
ENSEMBLE = 5
BATCH = 256               # lower to 128 if the GPU runs out of memory
ROOT = "cellfate_run"     # output folder (created next to this script)

HERE = Path(__file__).resolve().parent
CLOCK = str(HERE / "configs" / "clocks" / "fleischer_clock.json")


def discover_samples(data_dir: str) -> list[dict]:
    samples = []
    for mtx in sorted(glob.glob(os.path.join(data_dir, "*matrix.mtx.gz"))):
        base = mtx[: -len(".matrix.mtx.gz")]
        bc = base + ".barcodes.tsv.gz"
        label = re.search(r"_([^_]+)$", os.path.basename(base)).group(1)
        if not os.path.exists(bc):
            raise SystemExit(f"missing barcodes file for {mtx}")
        samples.append({"matrix": mtx, "barcodes": bc, "label": label})
    return samples


def main() -> None:
    import torch

    print(f"[env] torch {torch.__version__} | CUDA available: {torch.cuda.is_available()}"
          + (f" | GPU: {torch.cuda.get_device_name(0)}" if torch.cuda.is_available() else ""))
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cpu":
        print("[warn] no CUDA GPU detected -> training on CPU will be SLOW.")

    if not os.path.isdir(DATA_DIR):
        raise SystemExit(f"DATA_DIR not found: {DATA_DIR}\nPass your folder: python run_local.py <folder>")
    genes_hits = glob.glob(os.path.join(DATA_DIR, "*genes*.tsv.gz"))
    if not genes_hits:
        raise SystemExit(f"no *genes*.tsv.gz in {DATA_DIR} (must end in .tsv.gz)")
    genes = genes_hits[0]
    samples = discover_samples(DATA_DIR)
    print(f"[data] genes: {os.path.basename(genes)} | timepoints: {[s['label'] for s in samples]}")
    if not samples:
        raise SystemExit("no *.matrix.mtx.gz files found — check DATA_DIR")

    from cellfate.data import DataConfig, QCConfig
    from cellfate.data import run as build_run
    from cellfate.data.sources import GSE242423SingleCellSource

    if os.path.isdir(ROOT):
        shutil.rmtree(ROOT)

    print(f"\n[1/5] BUILD dataset (MAX_CELLS={MAX_CELLS}/timepoint) — streams droplets, ~10-20 min ...")
    src = GSE242423SingleCellSource(samples, genes, cell_line="HFF",
                                    min_genes=500, max_cells_per_sample=MAX_CELLS, seed=0)
    build_run(DataConfig(
        out=ROOT, gene_panel=f"{ROOT}/panel.json", n_genes=N_GENES,
        clock=CLOCK, modality="tf",
        qc=QCConfig(min_genes=500, max_mito_frac=0.20), label_tau=0.7,
        split_fracs=(0.7, 0.1, 0.1, 0.1), split_regimes=(REGIME,), primary_regime=REGIME,
        deconfound=True, seed=0), sources=[src])

    # sanity: ΔAge + fate by timepoint
    import numpy as np
    import pandas as pd

    from cellfate.common import io
    rows = []
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i in range(len(a["cell_id"])):
            rows.append({"tp": round(float(a["dose_time"][i][1]), 3), "y_age": a["y_age"][i],
                         "age_mask": bool(a["age_mask"][i]), "P_safe": a["y_cls"][i][0],
                         "P_loss": a["y_cls"][i][1], "P_death": a["y_cls"][i][2]})
    df = pd.DataFrame(rows)
    tps = sorted(df.tp.unique())
    df["timepoint"] = df.tp.map(dict(zip(tps, [s["label"] for s in samples], strict=False)))
    print(f"\n[sanity] cells: {len(df)} | age-valid: {int(df.age_mask.sum())}")
    print(df.groupby("timepoint", sort=False)[["P_safe", "P_loss", "P_death"]].mean().round(3).to_string())
    print(df[df.age_mask].groupby("timepoint", sort=False)["y_age"].mean().round(2).to_string())

    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run
    print(f"\n[2/5] TRAIN {ENSEMBLE}-member ensemble on {device} — this is the long step ...")
    train_run(TrainConfig(
        dataset_dir=ROOT, out=ROOT, regime=REGIME,
        d_cell=256, d_u=256, latent_dim=256, p_drop=0.2, lr=1e-3,
        epochs=EPOCHS, patience=10, batch_size=BATCH, ensemble_size=ENSEMBLE,
        base_seed=0, conformal_levels=(0.90,), device=device))

    # show the split sizes + calibration (the numbers that were broken last time)
    m = json.load(open(f"{ROOT}/bundle/metrics.json"))
    print("\n[check] split sizes + calibration (these must be non-degenerate now):")
    print(f"   n_train={m['n_train']}  n_val={m['n_val']}  n_calib={m['n_calib']}"
          f"  temperature={m['temperature']}  conformal_q={m['conformal_q']}")
    if m["n_calib"] == 0 or m["temperature"] == 1.0:
        print("   [WARN] calibration split still empty — check REGIME/MAX_CELLS.")

    from cellfate.evaluation.evaluate_cli import EvalConfig, evaluate
    print("\n[3/5] EVALUATE — baselines, gates, calibration, coverage ...")
    gates = evaluate(EvalConfig(bundle=ROOT, dataset=ROOT, regimes=(REGIME,), out=f"{ROOT}/reports"))
    print(json.dumps(gates, indent=2, default=str))

    # E-distance ranking (measured effect size)
    from cellfate.evaluation.metrics import edistance_to_control
    Xs, tp = [], []
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        Xs.append(a["X"])
        tp += [round(float(t[1]), 3) for t in a["dose_time"]]
    X = np.vstack(Xs)
    tp = np.array(tp)
    groups = np.array([dict(zip(sorted(set(tp)), [s["label"] for s in samples], strict=False))[t] for t in tp])
    try:
        ed = edistance_to_control(X, groups, "D0", n_pcs=50, max_cells=1000)
        print("\n[4/5] E-distance from D0:", {k: round(v, 2) for k, v in sorted(ed.items(), key=lambda x: x[1])})
    except Exception as e:
        print(f"\n[4/5] E-distance skipped ({e})")

    print("\n[5/5] SAVE checkpoint ...")
    shutil.make_archive("cellfate_bundle", "zip", f"{ROOT}/bundle")
    print(f"   bundle -> {os.path.abspath('cellfate_bundle.zip')}")
    print("\nDONE. Send the [sanity] table + the EVALUATE JSON back for interpretation.")


if __name__ == "__main__":
    main()
