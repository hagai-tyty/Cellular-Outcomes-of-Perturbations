"""
CellFate-Rx - Multi-donor training run (Gill + GSE242423), leave-cell-line-out.

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
from collections import Counter, defaultdict
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
REGIME = "holdout"        # hold out ONE Gill donor as test; HFF + rest -> train/val/calib
HOLDOUT_DONOR: str | None = None   # None -> auto-pick (first 'old' donor, else first). Set e.g. "O1".
HARMONIZE = True          # cross-modality control-anchoring + Gill Projection (fixes ΔAge scale)
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

    from cellfate.common.console import install_pretty_console
    install_pretty_console()   # cosmetic: progress bars + clean lines instead of JSON spam

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
    gill_donors = [c["cell_line"] for c in gill.plan()]        # e.g. ['N2','N3','O1','O2','Y1','Y2']
    # pick the held-out test donor: honour HOLDOUT_DONOR, else prefer an 'old' (O*) donor
    # (rejuvenation signal is largest there), else the first donor.
    if HOLDOUT_DONOR and HOLDOUT_DONOR in gill_donors:
        test_donor = HOLDOUT_DONOR
    else:
        test_donor = next((d for d in gill_donors if d.upper().startswith("O")), gill_donors[0])
    print(f"[data] Gill donors: {gill_donors}  ->  holding out '{test_donor}' as TEST")

    if os.path.isdir(ROOT):
        shutil.rmtree(ROOT)
    print(f"\n[1/5] BUILD combined dataset (regime={REGIME}, holdout={test_donor}, "
          f"harmonize={'ON' if HARMONIZE else 'off'}) - streaming ...")
    # GSE242423 first so the gene panel is fit on rich single-cell data. HFF + the
    # non-test donors split at the cell level (train/val/calib); test_donor -> test.
    build_run(DataConfig(
        out=ROOT, gene_panel=f"{ROOT}/panel.json", n_genes=N_GENES, clock=CLOCK, modality="tf",
        qc=QCConfig(min_genes=500, max_mito_frac=0.20), label_tau=0.7,
        split_fracs=(0.8, 0.1, 0.1, 0.0), split_regimes=(REGIME,), primary_regime=REGIME,
        holdout_cell_lines=(test_donor,), harmonize=HARMONIZE, harmonize_ref_dataset="gill_bulk",
        deconfound=True, seed=0), sources=[gse, gill])

    # ---- composition: how each cell line is distributed across splits ----
    paths = ArtifactPaths.of(ROOT)
    man = io.manifest_rows(io.load_manifest(paths))
    cell_lines = sorted({r.cell_line for r in man})
    if len(cell_lines) < 4:
        raise SystemExit(
            f"\n[FATAL] only {len(cell_lines)} cell line(s) built: {cell_lines}\n"
            "Leave-cell-line-out needs >=4 (HFF + several Gill donors). The Gill donors did not "
            "make it into the dataset -- check the [data] lines above and confirm the Gill "
            "expression file is the SeqMonk Log2-RPM report (sample columns like "
            "'N2_d11_CD13_Sendai_Exp1'), not a raw-count / Entrez-ID matrix.")
    splits = io.load_splits(paths, REGIME)
    cl_of = {r.cell_id: r.cell_line for r in man}
    per_cl_split: dict[str, Counter] = defaultdict(Counter)
    for r in man:
        per_cl_split[r.cell_line][splits.get(r.cell_id, "?")] += 1

    rows = []
    for sh in sorted(glob.glob(f"{ROOT}/shards/*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i in range(len(a["cell_id"])):
            rows.append({"cell_id": a["cell_id"][i], "P_loss": a["y_cls"][i][1],
                         "y_age": a["y_age"][i], "age_mask": bool(a["age_mask"][i])})
    df = pd.DataFrame(rows)
    df["cell_line"] = df.cell_id.map(cl_of)

    from cellfate.common.console import render_table
    print("\n[sanity] where each cell line's cells landed, and their average predictions")
    srows = []
    for cl in cell_lines:
        sub = df[df.cell_line == cl]
        age = sub[sub.age_mask]["y_age"]
        dist = " ".join(f"{k}={v}" for k, v in sorted(per_cl_split[cl].items()))
        age_s = f"{age.mean():+.2f}" if len(age) else "n/a"
        tag = "  <- TEST" if cl == test_donor else ""
        srows.append([cl, f"{len(sub):,}", dist, f"{sub.P_loss.mean():.3f}", age_s + tag])
    print(render_table(
        ["cell line", "cells", "where they went (train/val/calib/test)", "P(loss)", "ΔAge (yrs)"],
        srows, aligns=["l", "r", "l", "r", "r"]))
    print("   P(loss) = avg predicted chance the cell loses its identity (0-1).  "
          "ΔAge = avg years younger(-)/older(+).")
    print(f"\n   >>> HELD-OUT (test) donor: '{test_donor}'  <<<  (model never trains on it)")
    if HARMONIZE:
        print("   [acceptance] with harmonization ON, the Gill donors' mean ΔAge above should now be")
        print("   on a comparable scale to HFF (no +16..+64 artifact) -- that is A1-A3 from the spec.")

    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run
    print(f"\n[2/5] TRAIN {ENSEMBLE}-member ensemble on {device} (leave-cell-line-out) ...")
    train_run(TrainConfig(dataset_dir=ROOT, out=ROOT, regime=REGIME,
                          d_cell=256, d_u=256, latent_dim=256, p_drop=0.2, lr=1e-3,
                          epochs=EPOCHS, patience=10, batch_size=BATCH, ensemble_size=ENSEMBLE,
                          base_seed=0, conformal_levels=(0.90,), device=device))
    m = json.load(open(f"{ROOT}/bundle/metrics.json"))
    print(f"\n[check] trained on {m['n_train']:,} cells | validated on {m['n_val']:,} | "
          f"calibrated on {m['n_calib']:,}")
    print(f"        (temperature {m['temperature']:.2f} = confidence rescaling; "
          f"conformal_q {list(m['conformal_q'].values())[0]:.2f} = +/-interval width in years)")

    from cellfate.evaluation.evaluate_cli import EvalConfig, evaluate
    print("\n[3/5] EVALUATE  (held-out donor - model never saw it) ...")
    gates = evaluate(EvalConfig(bundle=ROOT, dataset=ROOT, regimes=(REGIME,), out=f"{ROOT}/reports"))
    g = gates.get(REGIME, gates)
    _tick = lambda ok: "PASS" if ok else "----"   # noqa: E731
    print("\n   gates:  "
          f"ranking {_tick(g.get('ranking_ok'))}   "
          f"calibration {_tick(g.get('ece_ok'))}   "
          f"coverage {_tick(g.get('coverage_ok'))}   "
          f"beats-all {_tick(g.get('beats_all_baselines'))}")

    rep = Path(ROOT, "reports", f"{REGIME}.json")
    if rep.exists():
        R = json.loads(rep.read_text())

        def prauc(d):
            v = [x for k, x in d.items() if k.startswith("prauc_")]
            return sum(v) / len(v) if v else float("nan")
        ests = [k for k in R if isinstance(R[k], dict) and any(kk.startswith("prauc_") for kk in R[k])]
        m_pr = prauc(R["model"]) if "model" in R else float("nan")
        m_ece = R.get("model", {}).get("ece", float("nan"))
        m_mae = R.get("model", {}).get("reg_mae", float("nan"))
        spearman = R.get("ranking", {}).get("spearman", float("nan"))
        coverage = R.get("coverage", float("nan"))

        import math

        def cell(v, width=8, prec=3):
            return (f"{v:>{width}.{prec}f}" if isinstance(v, (int, float)) and math.isfinite(v)
                    else f"{'n/a':>{width}}")

        def mark(model_v, base_v, lower_is_better):
            if not (math.isfinite(model_v) and math.isfinite(base_v)):
                return "  - "
            win = (model_v < base_v) if lower_is_better else (model_v > base_v)
            return " win" if win else "  - "

        # ---- readable table: model row on top, baselines below, aligned --------- #
        # ---- readable box table with a note column ----------------------------- #
        beats_all = {"pr": True, "mae": True}
        trows = [["model", cell(m_pr, 6), cell(m_ece, 6), cell(m_mae, 8, 2), "<- our model"]]
        for name in [e for e in ests if e != "model"]:
            d = R[name]
            bp, be, bm = prauc(d), d.get("ece", float("nan")), d.get("reg_mae", float("nan"))
            vp, vm = mark(m_pr, bp, False), mark(m_mae, bm, True)
            if math.isfinite(m_pr) and math.isfinite(bp) and not (m_pr > bp):
                beats_all["pr"] = False
            if math.isfinite(m_mae) and math.isfinite(bm) and not (m_mae < bm):
                beats_all["mae"] = False
            note = " ".join(x for x in [("model wins PR-AUC" if vp.strip() == "win" else ""),
                                        ("model wins MAE" if vm.strip() == "win" else "")] if x)
            trows.append([name, cell(bp, 6), cell(be, 6), cell(bm, 8, 2), note])
        print("\n" + render_table(
            ["estimator", "PR-AUC", "ECE", "ΔAge MAE", "vs model"], trows,
            aligns=["l", "r", "r", "r", "l"],
            title="   MODEL vs BASELINES on the held-out donor:"))
        print("   PR-AUC: fate accuracy, higher better (0-1).   "
              "ECE: calibration error, lower better (<0.05 good).")
        print("   ΔAge MAE: avg error in years, lower better.   "
              "'model' is our tool; the rest are simple baselines to beat.")
        from cellfate.common.console import baselines_legend
        print(baselines_legend())

        # ---- plain-language verdict (so the numbers interpret themselves) ------- #
        def rank_word(s):
            return "STRONG" if s >= 0.6 else "moderate" if s >= 0.4 else "weak"
        print("\n   WHAT THIS MEANS:")
        print(f"     - Ranking (the tool's main job): {rank_word(spearman)}  "
              f"(Spearman {cell(spearman,0,2).strip()}) - can it order perturbations correctly?")
        print(f"     - ΔAge magnitude: MAE {cell(m_mae,0,1).strip()}  "
              f"- {'beats' if beats_all['mae'] else 'ties/loses to'} the linear baseline (ridge)")
        cov_word = ("good" if isinstance(coverage, float) and coverage >= 0.85
                    else "off-target (in-distribution only)")
        print(f"     - Calibration on this unseen donor: {cov_word}  "
              f"(coverage {cell(coverage,0,2).strip()}, ECE {cell(m_ece,0,2).strip()})")
        print("     Note: one held-out donor is noisy; the leave-one-donor-out run "
              "(run_loocv.py) gives the honest mean+/-std.")

    print("\n[4/5] SAVE checkpoint ...")
    shutil.make_archive("cellfate_multi_bundle", "zip", f"{ROOT}/bundle")
    print(f"   bundle -> {os.path.abspath('cellfate_multi_bundle.zip')}")
    print("\n[5/5] DONE.  (one held-out donor - run_loocv.py rotates all 6 for the honest mean+/-std)")


if __name__ == "__main__":
    main()
