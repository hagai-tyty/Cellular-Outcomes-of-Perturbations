# CellFate-Rx — Document 2: the `cellfate.data` ETL package

This is the **data pipeline** that turns raw single-cell perturbation atlases into
validated, sharded training artefacts. It is built on top of the foundation
(`cellfate.common`, Document 1) and writes **only** the on-disk artefacts the other
packages read — it never imports models / training / inference / evaluation.

## What it produces

Running the pipeline writes a complete dataset under `out/`:

```
out/
  shards/<chunk>.parquet          # rows of Sample (the atomic unit: one perturbed cell)
  manifest_parts/<chunk>.parquet  # per-chunk manifest fragments
  manifest.parquet                # consolidated manifest (one row per cell)
  splits/{scaffold,cell_line,both}.json   # cell_id -> train/val/calib/test/drop
  scalers.json                    # per-gene + dose/time standardisation, prolif. coef
  gene_panel.json                 # the frozen HVG feature order (its hash is in every artefact)
  progress_tracker.json           # resumable chunk state
  dataset_summary.json            # counts, label distribution, panel hash, split sizes
```

Every shard row satisfies the foundation `Sample` contract: `X` (normalised HVG
expression), `u_chem_fp` (2048-bit Morgan fingerprint), `dose_time` =
`[log10(dose_uM), log(time_h)]`, `y_cls` (soft label over **safe/loss/death**, sums to 1),
`sig_scores` (same class order), `y_age` (ΔBiological age, **None** when masked),
`age_mask`, and the `cell_line`/`pert_id`/`scaffold_id`/`source` keys.

## The per-chunk pipeline

`fetch → QC → CP10k+log1p normalise → signature scores → soft labels →
cell-cycle score → ΔAge (vs matched controls; masked for cancer lines) →
proliferation deconfounding → project onto frozen panel → encode perturbation →
assemble Sample rows → write shard + manifest part`. After all chunks:
`consolidate manifest → group-aware splits (3 regimes) → fit scalers on the
primary regime's TRAIN rows only → summary`.

Two scientific safeguards are baked in: **group-aware splits** (a Bemis-Murcko
scaffold / cell line is entirely train *or* test — no leakage), and
**proliferation deconfounding** (the linear cell-cycle component is removed from
the age label, so "rejuvenation" can't be faked by making cells divide).

## How to run

Zero-dependency synthetic run (no downloads, exercises the whole pipeline + CLI):

```bash
python scripts/build_dataset.py data=synthetic device=cpu \
    data.out=/tmp/synth_ds paths.gene_panel=/tmp/synth_ds/panel.json
```

Programmatic / testable entry point (inject sources, no Hydra needed):

```python
from cellfate.data import DataConfig, run, SyntheticSource
summary = run(DataConfig(out="ds", gene_panel="ds/panel.json"),
              sources=[SyntheticSource(name="synth"), SyntheticSource(name="tahoe")])
```

Real chem dataset (Tahoe-100M + sci-Plex): `data=chem_v1`. This needs the heavy
extras and data access — see "Wiring real sources" below.

## Dependencies

Light deps (numpy / pandas / scipy / scikit-learn / pydantic / pyarrow) are
enough to import the package, run the full synthetic pipeline, and pass the whole
test suite. Heavy deps are **lazy-imported** inside the functions that need them:

* `rdkit` — real Morgan fingerprints + Bemis-Murcko scaffolds. Without it the
  pipeline falls back to a deterministic hashed fingerprint (and warns once), so
  CI still runs; install it for production.
* `datasets` / `anndata` / `scanpy` — used only by the real source connectors to
  stream Tahoe-100M and read `.h5ad` files.

## Wiring real sources (the one open task for the data engineer)

`SyntheticSource` is complete and used everywhere for tests. `TahoeSource` and
`SciplexSource` are structured skeletons: implement their `plan()` (enumerate
cell-line × drug chunks) and `fetch()` (return a `RawChunk` of raw counts + gene
symbols + the `obs` columns). Everything downstream — QC, normalisation,
labelling, splitting, scaler fitting — is source-agnostic and already done.
`tahoe` is in `CANCER_SOURCES`, so its age labels mask automatically.

## Tests

`pytest` → **69 passed** (41 foundation + 28 data). The data tests cover every
scientific function (soft-label normalisation, signature ordering, deconfounding
actually removing the cell-cycle correlation, split group-disjointness, QC,
normalisation, fingerprint determinism, dose/time encoding, age masking) plus a
full **end-to-end build** that asserts every produced artefact validates against
the foundation contracts. A printed 17-point "contract match" report is
reproducible from the snippet in the delivery notes.
