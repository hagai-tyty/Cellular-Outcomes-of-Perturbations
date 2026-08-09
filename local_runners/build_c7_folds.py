"""CHANGE C-7 ADOPTION — rebuild the six LOOCV folds with the integrity gate ON.

    python local_runners/build_c7_folds.py "D:\\GSE242423" "D:\\Gill"

Writes `cellfate_loocv_<donor>_c7/`. **DATASET ONLY -- no training, no bundle.**

WHY DATASET-ONLY IS ENOUGH
--------------------------
The next thing that needs these folds is Stage 3a's gate, `experiments/test18_forward_gate.py`.
It imports `ArtifactPaths` and `gather_split` and fits its own **ridge** -- no `Predictor`, no
`ModelEstimator`, no bundle. Its docstring says "needs cellfate_loocv_* bundles", but that is
loose wording: what it reads is the built dataset's shards and splits.

So adopting C-7 far enough to unblock 3a costs a **build**, not a retrain. The full retrain (and
the Stage 1 guard re-report over SIX folds, per option (c)) is still required before C-7 is
adopted for anything that consumes a trained model -- that is a separate, pre-registered run.

WHAT CHANGES vs ARM A
---------------------
Exactly one config field: `bulk_integrity_gate=True`. Everything else is copied from
`run_multi_local.main` verbatim so the comparison is one change, not two.

Expected effect, from the verified end-to-end check:
  * 5 of Gill's 124 columns rejected (G1 library band and/or G2 dynamic range);
  * `N2_Fib_Sendai_Exp2` among them, so donor N2 is left with **zero** admissible controls;
  * rule 4 masks N2's ΔAge labels -- **the donor and its fold survive** (option (c));
  * B2' asserts no line both falls back and keeps a label.

`AGE_MASKED_DATASETS` is left EMPTY, i.e. arm A's setting, so `_c7` differs from `_armA` by the
gate alone.
"""
from __future__ import annotations

import os
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
SUFFIX = "_c7"


def main() -> int:
    gse_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
    gill_dir = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"

    import importlib.util
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("run_multi_local", here / "run_multi_local.py")
    rml = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rml)

    from cellfate.common import constants as _C
    from cellfate.common.console import install_pretty_console
    from cellfate.data import DataConfig, QCConfig
    from cellfate.data import run as build_run
    from cellfate.data.sources import GillReprogrammingSource, GSE242423SingleCellSource
    install_pretty_console()

    # Arm A's setting: the gate is the ONLY difference between `_c7` and `_armA`.
    _C.AGE_MASKED_DATASETS = frozenset()

    gse_samples, gse_genes = rml.discover_gse(gse_dir)
    gill_expr, gill_series = rml.discover_gill(gill_dir)
    print(f"[data] Gill expr  : {os.path.basename(gill_expr)}")

    t0 = time.time()
    for donor in DONORS:
        root = f"cellfate_loocv_{donor}{SUFFIX}"
        if os.path.isdir(root):
            shutil.rmtree(root)
        gse = GSE242423SingleCellSource(gse_samples, gse_genes, cell_line="HFF", min_genes=500,
                                        max_cells_per_sample=rml.MAX_CELLS,
                                        cells_per_run=rml.CELLS_PER_RUN, seed=0)
        gill = GillReprogrammingSource(gill_expr, gill_series)
        print(f"\n=== C-7 build: holdout {donor} -> {root} ===")
        build_run(DataConfig(
            out=root, gene_panel=f"{root}/panel.json", n_genes=rml.N_GENES, clock=rml.CLOCK,
            modality="tf", qc=QCConfig(min_genes=500, max_mito_frac=0.20), label_tau=0.7,
            split_fracs=(0.8, 0.1, 0.1, 0.0), split_regimes=(rml.REGIME,),
            primary_regime=rml.REGIME, holdout_cell_lines=(donor,), harmonize=True,
            harmonize_ref_dataset="gill_bulk", deconfound=True, seed=0,
            bulk_integrity_gate=True),
            sources=[gse, gill])
        rej = sorted(getattr(gill, "rejected_samples", {}))
        print(f"    rejected {len(rej)} Gill samples: {rej}")
        print(f"    lines_without_controls: {sorted(gill.lines_without_controls())}")
        print(f"    elapsed {time.time() - t0:.0f}s")

    print(f"\n=== all six C-7 folds built in {time.time() - t0:.0f}s ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
