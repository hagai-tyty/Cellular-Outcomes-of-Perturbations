#!/usr/bin/env python
"""Fit a real linear aging clock and write weights JSON for LinearClock.from_json.

Turnkey for GSE113957 (Fleischer 2018, human dermal fibroblasts, ages 1-96).

HOW TO GET GSE113957:
  1. Download the processed counts from NCBI GEO accession GSE113957
     (genes x samples table; gene SYMBOLS in the first column).
  2. Get per-sample ages from the GEO series matrix ("characteristics: age").
     Save a CSV with columns: sample_id,age   (one row per sample/column).

USAGE:
  python scripts/fit_clock.py \
      --counts GSE113957_counts.tsv --genes-axis rows \
      --metadata GSE113957_ages.csv --sample-col sample_id --age-col age \
      --panel artifacts/<run>/panel.json \
      --out configs/clocks/fleischer_clock.json

Then set  clock: configs/clocks/fleischer_clock.json  in your data config.

Works for ANY age-labelled matrix (e.g. Tabula Muris Senis pseudobulk) -- only the
column/orientation flags change.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cellfate.data.clock_fit import fit_linear_clock  # noqa: E402


def _read_matrix(path: str, genes_axis: str) -> tuple[np.ndarray, list[str], list[str]]:
    """Return (counts[n_samples, n_genes], genes, sample_ids)."""
    sep = "\t" if path.endswith((".tsv", ".txt", ".gz")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    if genes_axis == "rows":          # rows=genes, cols=samples -> transpose
        df = df.T
    # now rows=samples, cols=genes
    genes = [str(c) for c in df.columns]
    samples = [str(i) for i in df.index]
    return df.to_numpy(dtype=np.float64), genes, samples


def main() -> None:
    ap = argparse.ArgumentParser(description="Fit a linear transcriptomic aging clock.")
    ap.add_argument("--counts", required=True, help="counts matrix (CSV/TSV)")
    ap.add_argument("--genes-axis", choices=["rows", "cols"], default="rows",
                    help="are genes on rows (GEO default) or columns?")
    ap.add_argument("--metadata", required=True, help="CSV with sample id + age columns")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--age-col", default="age")
    ap.add_argument("--panel", default=None, help="panel.json to align the clock to (optional)")
    ap.add_argument("--normalized", action="store_true", help="counts already log1p CP10k")
    ap.add_argument("--out", required=True, help="output clock JSON")
    args = ap.parse_args()

    counts, genes, samples = _read_matrix(args.counts, args.genes_axis)
    meta = pd.read_csv(args.metadata)
    age_by_sample = dict(zip(meta[args.sample_col].astype(str),
                             meta[args.age_col].astype(float), strict=False))
    missing = [s for s in samples if s not in age_by_sample]
    if missing:
        raise SystemExit(f"{len(missing)} samples in the matrix have no age in metadata "
                         f"(first few: {missing[:5]})")
    ages = np.array([age_by_sample[s] for s in samples], dtype=np.float64)

    panel_genes = None
    if args.panel:
        pj = json.loads(Path(args.panel).read_text())
        panel_genes = pj.get("genes") or pj.get("panel") or pj
        if isinstance(panel_genes, dict):
            panel_genes = list(panel_genes)

    clock, metrics = fit_linear_clock(counts, genes, ages,
                                      panel_genes=panel_genes, normalized=args.normalized)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    clock.to_json(args.out, meta={"source": Path(args.counts).name, **metrics})

    print(f"fit clock on {metrics['n_samples']} samples x {metrics['n_genes']} genes")
    print(f"  CV MAE   = {metrics['cv_mae_years']:.2f} years")
    print(f"  CV Pearson = {metrics['cv_pearson']:.3f}   (ridge alpha={metrics['alpha']:.3g})")
    print(f"  wrote {args.out}  ->  set  clock: {args.out}  in your data config")


if __name__ == "__main__":
    main()
