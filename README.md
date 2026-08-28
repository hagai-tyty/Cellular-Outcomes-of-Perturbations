# CellFate-Rx

Two lines of work live in this repository, and they are kept separate on purpose.

**Generation 1 — shipped.** Clone-level prospective prediction in a lineage-traced melanoma line:
does a clone's pretreatment transcriptional state say anything about *which* of six observed
experimental conditions it survives? Complete, locked, and summarised immediately below.

**The earlier reprogramming line.** OSKM fate and ΔAge, with a safety gate. Most of its original
claims did **not** survive being measured properly; its honest status is kept in full further down
and in [`ARCHITECTURE.md` §11](ARCHITECTURE.md#11-where-the-project-actually-stands--2026-08-18).
Nothing there was deleted when Generation 1 shipped.

The whole repository is heavily tested (2,259 tests) and its measurements are reproducible
bit-for-bit.

---

## Generation 1 — shipped

```text
  GEN1_MANUSCRIPT_READY

  evidence lock   455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9
  claim lock      23ea00b808d1ae6a5f3b19e186a9fd0b327ba4b500c5c2a65d4c254142b431ab
  package         68a1fca2e260fb8e834aaca5949fec4fb05668d823cb0792564b10df1e99bf76
```

**What was measured.** In WM989 (GSE279162), 1,401 barcoded clones were split across six observed
experimental conditions — Acid, Cisplatin, CoCl2, Dabrafenib, Doxorubicin, Trametinib. Under a test
preregistered in full, a frozen state-by-condition interaction model improves clone-specific
ordering of those conditions over a non-interactive additive model:

```text
  delta_RANK   +0.051605   CI95 [+0.037197, +0.065571]
  null         0 of 1000 full-refit permutation draws reached the observed value, p < 0.001
  population   892 of 1,401 clones, clone-held-out, frozen before any result existed
```

**What it is not.** The outcome is an observed post-treatment clone-detection proxy and is **not
death**, sensitivity, resistance or clinical response. The six conditions are the entire supported
vocabulary; there is no claim about unseen conditions, other cell lines, or patients, and the model
emits no calibrated probability. Captured pretreatment clone abundance remains roughly 3.45× the
whole state contribution — state adds something specific on top of abundance, it does not dominate
it. Independent biological replication has not been performed and is Generation 2 work.

**Read it.**

| | |
|---|---|
| [`results/manuscript/MANUSCRIPT.md`](results/manuscript/MANUSCRIPT.md) | the write-up |
| [`results/manuscript/REPRODUCIBILITY.md`](results/manuscript/REPRODUCIBILITY.md) | how to re-run it |
| [`results/claim_lock/GEN1_CLAIMS.md`](results/claim_lock/GEN1_CLAIMS.md) | what may and may not be said |
| [`results/evidence_lock/GEN1_EVIDENCE_LOCK.md`](results/evidence_lock/GEN1_EVIDENCE_LOCK.md) | the 54 locked artifacts |
| `plans/(newer)practical plans/RECORDs/stage_25_RECORD.md` | the result itself |

**Verify before trusting any of it.** Each command re-hashes its layer and refuses if anything
moved; each was shown to catch a one-bit change before it was issued.

```bash
python experiments/run_stage24_gen1_tool.py --stage 24c   # rebuild the gitignored model artifact
python experiments/run_gen1_evidence_lock.py --verify
python experiments/run_gen1_claim_lock.py --verify
python experiments/run_gen1_manuscript.py --verify
```

---

## The earlier reprogramming line — status as of 2026-08-18

Kept as recorded. Where a claim turned out to be wrong it stays here with the correction beside it,
rather than being edited out.

### What actually holds

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
| `experiments/` | 99 read-only diagnostics and stage executors, one per scientific question |
| `plans/` | the decision record, one file per stage; `RECORDs/` is what each stage found |
| `results/`, `scorecard/` | every measurement ever taken |
| `CHANGES.md` | append-only log; wrong claims stay, with corrections beside them |
| **`ARCHITECTURE.md`** | **the real documentation** |

Current fold set: **`_s16`**. `_c7t` and `_s12` are retained as comparison baselines only.

## Ground rules

Bars are registered before runs. Comparisons are paired across folds. Records are appended, never
rewritten. One change per experiment. See `plans/REF_GROUND_RULES.md`.

## License

See [LICENSE](LICENSE).
