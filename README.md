# CellFate-Rx

[![CI](https://github.com/hagai-tyty/Cellular-Outcomes-of-Perturbations/actions/workflows/ci.yml/badge.svg)](https://github.com/hagai-tyty/Cellular-Outcomes-of-Perturbations/actions/workflows/ci.yml)

**Probabilistic single-cell triage for safe rejuvenation.**

CellFate-Rx takes a single-cell transcriptome and a proposed perturbation (a drug,
gene, or transcription factor, with dose and time) and returns a **calibrated**
prediction: the probabilities that the perturbation leaves the cell's identity intact,
causes identity loss, or kills it — plus a **ΔBiological-age with honest uncertainty**.
On top of these it computes a **Rejuvenation Efficacy Score (RES)** that flags a
perturbation only when it is simultaneously *safe*, *rejuvenating*, and *in-distribution*.

The goal is in-silico triage: score many candidate perturbations first, and spend scarce
wet-lab runs only on the ones that clear the bar. The scientific heart is **partial
reprogramming without losing cell identity or tipping into cancer** — rewinding
biological age while keeping the cell what it is.

> **Scope, stated honestly.** This repository is a complete, tested ML *system*: every
> stage runs end-to-end and reports honestly on itself. It is validated on a synthetic
> data source, **not** yet on real biology. The acceptance gates (below) are built so
> that when you connect real perturbation data, they tell you truthfully whether the
> model is trustworthy — rather than assuming it is.

---

## How it works: five stages, one artefact chain

The system is one Python package (`cellfate`) of five sub-packages that communicate
**only through on-disk artefacts** — each stage reads the previous stage's output from
disk and writes its own. This keeps the stages independently testable and reproducible.

```
raw perturbation data
        │
        ▼  cellfate.data ............ QC → normalise → aging clock (ΔAge vs vehicle
        │                             control, cell-cycle deconfounded) → soft fate
        │                             labels → scaffold / cell-line / both splits
   shards/ + splits/ + panel
        │
        ▼  cellfate.training ........ shared encoder + multi-task heads (3-way fate +
        │                             heteroscedastic ΔAge), deep ensemble, temperature
        │                             scaling, split-conformal calibration
     bundle/
        │
        ▼  cellfate.inference ....... encode → ensemble predict → OOD flag → conformal
        │                             ΔAge interval → Rejuvenation Efficacy Score + gate
   Response / RES
        │
        ▼  cellfate.evaluation ...... baselines + metrics + external ΔAge anchors +
        │                             automated acceptance gates
     reports/
```

`cellfate.common` underpins all of them: artefact schemas, the frozen gene panel and its
hash, scalers, and shared constants.

### The acceptance gates

Evaluation turns "does it work?" into a falsifiable, per-regime verdict:

| gate | meaning |
|---|---|
| `beats_all_baselines` | beats **every** baseline (mean, ridge, X-only, U-only, kNN-fp, predict-control) on PR-AUC **and** ΔAge MAE |
| `ece_ok` | classifier is calibrated (ECE < 0.05) |
| `coverage_ok` | conformal ΔAge interval covers at its nominal 90% (±3%) |
| `ranking_ok` | RES ranks rejuvenation (Spearman vs measured ΔAge > 0.3) |

Reported separately for three held-out regimes — **leave-drug-out**,
**leave-cell-line-out**, and **both-unseen** — never as a single pooled number.
"Beats every baseline" is an automated gate on purpose: benchmarks repeatedly find
sophisticated models that fail to beat trivial ones, so the harness checks.

---

## Install

Requires Python ≥ 3.11.

```bash
git clone https://github.com/<you>/cellfate-rx.git
cd cellfate-rx
pip install -e ".[all]"      # or: pip install -e .   (core only)
```

Optional extras are grouped by stage: `data` (anndata/scanpy/rdkit/hf),
`model` (torch), `serve` (fastapi/uvicorn), `dev` (pytest/ruff/mypy). `all` pulls
everything. rdkit is optional — the fingerprint path falls back to a hashed descriptor.

---

## Quickstart

Run the whole chain on the built-in synthetic source. **Python:**

```python
from cellfate.data import DataConfig, QCConfig, SyntheticSource, run as build
from cellfate.training.train_model import TrainConfig, run as train
from cellfate.inference import Predictor, score_shard
from cellfate.evaluation import evaluate, EvalConfig

# 1. build a dataset (shards/ + splits/ + panel)
build(DataConfig(out="run", gene_panel="run/panel.json", n_genes=96,
                 qc=QCConfig(min_genes=5, max_mito_frac=0.5),
                 split_regimes=("scaffold", "cell_line", "both"), primary_regime="cell_line"),
      sources=[SyntheticSource(name="synth", n_scaffold_families=7)])

# 2. train a calibrated bundle
train(TrainConfig(dataset_dir="run", out="run", regime="cell_line", ensemble_size=3))

# 3. serve: score cells -> safe/loss/death + ΔAge interval + RES
pred = Predictor("run")
responses, cell_ids = score_shard(pred, "run/shards/0000.parquet")

# 4. evaluate: baselines + metrics + acceptance gates
gates = evaluate(EvalConfig(bundle="run", dataset="run", out="run/reports"))
```

**Command line** (installed as console scripts):

```bash
cellfate-build-dataset          # Hydra-configured; see configs/data/
cellfate-train                  # Hydra-configured; see configs/train/
cellfate-serve   --bundle run --shard run/shards/0000.parquet --top 10
cellfate-evaluate --bundle run --dataset run --regimes scaffold cell_line both
```

---

## The aging clock (real vs. random)

ΔAge is measured by a transcriptomic aging clock. For synthetic/smoke runs the clock is
an explicit **random** placeholder (`clock: random`) whose ages are not meaningful. For
real runs, fit a real clock and point `clock:` at the weights file — anything else fails
loud (no silent fallback):

```bash
# fit a real clock from an age-labelled matrix (e.g. GSE113957 human fibroblasts, ages 1-96)
python scripts/fit_clock.py --counts GSE113957_counts.tsv --genes-axis rows \
    --metadata ages.csv --sample-col sample_id --age-col age \
    --panel artifacts/<run>/panel.json --out configs/clocks/fleischer_clock.json
# then set  clock: configs/clocks/fleischer_clock.json  in your data config
```

The **reprogramming connector** (`ReprogrammingSource`) is the age-valid arm (OSKM / partial
reprogramming time-courses): unlike cancer sources, its ΔAge is kept, and OSKM is encoded as
a transcription-factor cocktail (timepoint → `time_h`).

## Testing

```bash
pip install -e ".[dev]"
pytest            # 170 tests
ruff check .      # lint
```

The suite covers unit correctness (metrics cross-checked against scikit-learn/scipy,
calibration, conformal coverage, RES gating logic), reproducibility, and full
build→train→serve→evaluate integration on synthetic data.

---

## A note on the synthetic demo

Because the bundled synthetic data is deliberately simple, two honest things happen, and
the harness is designed to surface both: (1) the synthetic expression is trivially
linearly separable, so linear baselines are near-perfect and the model does **not** beat
them — `beats_all_baselines` correctly reports FALSE (the gate doing its job); and
(2) the two external ΔAge anchors (a DNA-methylation clock and an OSKM
partial-reprogramming holdout) require real data that isn't shipped, so they report
`not_available` rather than a fabricated number. On real data with genuine signal, these
gates are what separate a trustworthy model from a plausible-looking one.

---

## Repository layout

```
src/cellfate/
├── common/       # artefact schemas, gene panel + hash, scalers, constants
├── data/         # QC, aging clock, fate labels, splits, dataset builder
├── models/       # encoder + multi-task heads
├── training/     # training loop, deep ensemble, temperature + conformal calibration
├── inference/    # Predictor, OOD, conformal intervals, RES + serving
└── evaluation/   # baselines, metrics, external validation, acceptance gates
configs/          # Hydra configs (data / model / train / infer / eval)
scripts/          # thin CLI entry points
tests/            # 170 tests
DATA_README.md    # data-layer design notes
```

## License

See `pyproject.toml` (`license`). Set this to whatever suits you before making the repo
public — e.g. MIT for permissive reuse, or keep it private/all-rights-reserved.
