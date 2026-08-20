# CellFate-Rx

Predicting the cellular outcome of OSKM reprogramming from single-cell and bulk transcriptomes:
a **fate** call (identity preserved / lost / apoptotic) and a **ΔAge** estimate, with calibrated
uncertainty and an explicit safety gate.

> **Status: a research instrument, not a validated predictor.**
> The pipeline is heavily tested (1583 tests) and its measurements are reproducible bit-for-bit.
> Most of the project's original claims did **not** survive being measured properly.
> **[`ARCHITECTURE.md` §11](ARCHITECTURE.md#11-where-the-project-actually-stands--2026-08-18) is
> the honest status; read it before trusting any number here.**

---

## What actually holds, as of 2026-08-18

| | |
|---|---|
| ✅ **Fate, within a timepoint** | stratified AUC **0.917**, permutation p = 0.0091 — on **12 pairs**. Real, thin |
| ✅ **Safety gate** | sensitivity **0.670**, specificity 0.929 (Stage 16) |
| ✅ **Within-donor ΔAge ranking** | Spearman **0.942** |
| ✅ **`top100` clock vs methylation** | MAE **7.15** against a **7.30** instrument floor — CI spans zero |
| ❌ **Same-timepoint ΔAge prediction** | **circular** (ρ 0.96–0.99): the label is a linear functional of the input |
| ❌ **Predicting future ΔAge** | partial **−0.064** after donor age. Donor age is free |
| ❌ **RES** (the headline output) | **identically zero** — σ_age is 2–4.5× \|µ_age\|, so no cell earns confident-rejuvenation credit |

**Headline `fate_prauc` ≈ 0.96 is mostly a clock.** `dose_time` is a model input, and on 4 of 6
held-out donors *zero* timepoints carry more than one class — a lookup table on the hour scores
1.000 there. Use the stratified number instead.

---

## Install

```bash
pip install -e ".[dev]"
```

Python 3.11/3.12. PyTorch, and `[data]` extra for scanpy/rdkit.

## Quickstart

```bash
python scripts/build_dataset.py    # shards/ + splits/ + panel
python scripts/train.py            # calibrated 5-member bundle
python scripts/evaluate.py         # baselines, metrics, gates
python scripts/serve.py            # score cells
```

Real research runs go through `local_runners/`, not these CLIs:

```bash
CELLFATE_FOLD_SUFFIX=_s16 CELLFATE_BULK_GATE=1 python local_runners/run_loocv.py "D:\GSE242423" "D:\Gill" --arm A
```

## Grading a change

Never by eye. Snapshot, then diff:

```bash
python scorecard.py snapshot --tag my_change
python scorecard.py compare c7t_stage16 my_change
```

A change is real only if the paired 95 % CI across folds excludes zero.

## Tests

```bash
pytest -q
```

---

## Layout

| | |
|---|---|
| `src/cellfate/` | the package — `data/` `models/` `training/` `evaluation/` `inference/` |
| `local_runners/` | the drivers actually used day to day |
| `experiments/` | 86 read-only diagnostics, one per scientific question |
| `plans/` | the decision record, one file per stage |
| `results/`, `scorecard/` | every measurement ever taken |
| `CHANGES.md` | append-only log; wrong claims stay, with corrections beside them |
| **`ARCHITECTURE.md`** | **the real documentation** |

Current fold set: **`_s16`**. `_c7t` and `_s12` are retained as comparison baselines only.

## Ground rules

Bars are registered before runs. Comparisons are paired across folds. Records are appended, never
rewritten. One change per experiment. See `plans/REF_GROUND_RULES.md`.

## License

See [LICENSE](LICENSE).
