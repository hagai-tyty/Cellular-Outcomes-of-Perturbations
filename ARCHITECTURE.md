# CellFate-Rx — Architecture

**A map of the system as it actually is**, read from the code on 2026-08-03. Where this document
and `README.md` disagree, this one is right — the README predates most of the pipeline.

---

## Contents

- [1. What this system does](#1-what-this-system-does)
- [2. The 30-second picture](#2-the-30-second-picture)
- [3. Directory map](#3-directory-map)
- [4. The core data flow](#4-the-core-data-flow)
- [5. Key components, module by module](#5-key-components-module-by-module)
- [6. Every loop in the system](#6-every-loop-in-the-system)
- [7. On-disk artefacts](#7-on-disk-artefacts)
- [8. The scientific-method layer](#8-the-scientific-method-layer)
- [9. Entry points](#9-entry-points)
- [10. Cross-cutting invariants](#10-cross-cutting-invariants)

---

## 1. What this system does

CellFate-Rx predicts the **cellular outcome of a perturbation** on two axes at once:

| head | question | output |
|---|---|---|
| **fate** | does the cell keep its identity, lose it, or die? | 3-class probability |
| **age** | does the perturbation make the cell biologically younger? | **ΔAge** in years, with an interval |

The two are combined into a single decision score, **RES**, which ranks candidate perturbations by
*confident rejuvenation without loss of identity*:

```
RES = Φ(S) · S^k · g(R_eff) · exp(−λ · P_loss)

  Φ(S)   = sigmoid((S − τ_safe)/w)      smooth safety floor, no cliff
  R_eff  = max(0, −(µ_age + z_conf·σ_age))   credits only CONFIDENT rejuvenation
  g(R)   = R / (R + κ)                  concave, saturating
```

Two properties make the system unusual and drive most of its design:

1. **ΔAge is control-relative.** It is `clock(cell) − clock(that line's control)`, so a per-line
   baseline is part of the label's definition, not a post-hoc adjustment.
2. **Deployment is out-of-distribution by construction.** The tool is asked about donors it has
   never seen, so every calibrator is fitted on **cross-donor** statistics rather than on a random
   held-out split (see [inner LODO](#calibration-loops--trainingxdonor_calibpy)).

**Data.** Two real corpora, harmonised into one feature space:

| corpus | modality | size | role |
|---|---|---|---|
| **GSE242423** (`hff_sc`) | single-cell RNA, HFF line, D0→iPSC | ~42 500 cells | volume; 99.7 % of all ΔAge labels |
| **Gill 2022 / GSE165176** (`gill_bulk`) | bulk RNA, 6 donors (N2, N3, O1, O2, Y1, Y2) | ~124 samples | the age-valid rejuvenation signal |

The clock is **Fleischer 2018** (`configs/clocks/fleischer_clock.json`), a frozen external artefact
with `cv_mae_years = 12.27` and a validated `age_range = [1.0, 96.0]`. It is **validated, never
refitted** — refitting it to improve our own numbers would be fitting the test.

> ⚠️ **The clock is under active scrutiny, and two findings bear on every ΔAge number in this repo.**
> It is a *dense* ridge over 33 155 genes fitted from 133 samples, and (Stage 1.5.6) restricting it
> to its ~100 largest-weight genes cuts MAE against methylation from **16.61 → 5.36 yr** on the Gill
> transient arm — thousands of near-zero weights each drift slightly and sum to a −14 yr offset.
> Separately (Stage 1.5.6 step 1c) **harmonization applies a per-*gene* reweighting, not a rescaling**,
> worth a ×2.15 gain on HFF's ΔAge. The two interact: sparsifying concentrates the clock onto exactly
> the genes harmonization amplifies most, so the sparse clock **improves Gill and degrades HFF**.
> Nothing has been adopted. See `plans/STAGE_1_5_6_SPARSE_CLOCK.md`.

---

## 2. The 30-second picture

```
   D:\GSE242423   D:\Gill                    ← raw GEO, OUTSIDE the repo
        │             │
        ▼             ▼
  ┌──────────────────────────┐
  │  ETL  (data/)            │  chunked, resumable, two-pass
  │  QC → normalise → clock  │
  │  → ΔAge → deconfound     │
  └────────────┬─────────────┘
               ▼
        shards/*.parquet  +  manifest  +  splits  +  scalers
               │
               ▼
  ┌──────────────────────────┐
  │  TRAIN  (training/)      │  5-member deep ensemble
  │  shared trunk → 2 heads  │
  └────────────┬─────────────┘
               ▼
  ┌──────────────────────────┐
  │  CALIBRATE               │  inner LODO  ← ~6× the training cost
  │  temperature · Platt     │
  │  conformal · σ-scale     │
  └────────────┬─────────────┘
               ▼
          bundle/          ← the deployable unit
               │
      ┌────────┴────────┐
      ▼                 ▼
  EVALUATE          INFERENCE
 (evaluation/)      (inference/)  →  RES ranking, HTTP or CLI
      │
      ▼
  scorecard.py  →  snapshots  →  paired 6-fold comparisons
```

---

## 3. Directory map

### Source

| path | what lives here |
|---|---|
| **`src/cellfate/common/`** | Foundation shared by every stage: `io.py` (Parquet shard/manifest schemas, atomic writes), `schemas.py` (Pydantic `Sample`, `ManifestRow`, `AgeProvenance` + validators), `constants.py`, `scalers.py`, `calibration.py`, `seeding.py` (global determinism), `logging.py`, `console.py`, `errors.py`, `panel.py`, `progress.py` |
| **`src/cellfate/data/`** | The ETL. `build_dataset.py` is the orchestrator; `sources.py` the corpus adapters; `aging.py` the ΔAge definition and label policy; `harmonize.py` the cross-modality alignment; `proliferation.py` the cell-cycle deconfounder; plus `qc`, `normalize`, `labels`, `signatures`, `perturbation`, `splits`, `chunking`, `clock_fit`, `assemble` |
| **`src/cellfate/models/`** | The network. `network.py` (`CellFateNet`: shared trunk → two heads), `encoders.py` (cell / chem / TF), `heads.py`, `losses.py` (focal, masked-Huber, Kendall–Gal `MultiTaskLoss`) |
| **`src/cellfate/training/`** | `train_model.py` (orchestrator), `train.py` (loops + `_AgeWindow`), `dataset.py` (shard→tensor), `xdonor_calib.py` (**inner LODO**), `calibrate.py`, `conformal.py`, `ood.py`, `bundle.py`, `metrics.py` |
| **`src/cellfate/evaluation/`** | `evaluate_cli.py` (gates), `baselines.py` (mean / ridge / x-only / u-only / kNN / predict-control), `metrics.py`, `regimes.py`, `report.py`, `external_validation.py`, `data.py` |
| **`src/cellfate/inference/`** | `predictor.py` (loads a bundle), `res.py` (**RES**), `conformal.py`, `ood.py`, `encode.py`, `schema.py` (request/response), `service.py` (CLI + optional FastAPI) |

### Configuration, entry points, and governance

| path | what lives here |
|---|---|
| **`configs/`** | Hydra config tree — `config.yaml` composes `data/`, `model/`, `train/`, `infer/`, `eval/`. **`configs/clocks/`** holds the frozen clocks (Fleischer RNA; Horvath skin&blood 2018 and multi-tissue 2013 for methylation cross-checks). `configs/panels/` holds the frozen gene-panel order |
| **`scripts/`** | Thin Hydra CLIs: `build_dataset.py`, `train.py`, `evaluate.py`, `fit_clock.py`, `serve.py` |
| **`local_runners/`** | The drivers actually used day to day. **`run_multi_local.py`** = one full fold (build → train → evaluate → bundle). **`run_loocv.py`** = rotates all 6 donors. Plus `run_local.py`, `run_fate_local.py`, `evaluate_only.py`, `diag_harmonize.py`, `show_ui.py` |
| **`plans/`** | The project's decision record. `00_START_HERE.md`, `MASTER_PLAN.md`, `REF_GROUND_RULES.md` (the rules everything is graded against), `REF_ARCHITECTURE.md`, `REF_DATA_STRATEGY.md`, and one file per stage. `plans/archive/` holds superseded drafts — kept, never rewritten |
| **`tests/`** | 47 files, **843 tests**. Unit tests plus **registered-bar tests** (`test_bars_resolvable.py`) and invariance guards (`test_ci_deconfounder_arm_invariance.py`, `test_c5c_age_accumulation.py`, `test_arm_c_label_shuffle.py`, `test_results_paths.py`) |
| **`plan_tests/`** | Scripts a *plan* requires: pre-registered bars (`register_*_bar.py`), pre-flight gates (`step6_preflight.py`), bit-identity verifiers (`verify_age_mask_identical.py`, `verify_stage1_5.py`, `verify_1a.py`), smoke tests |
| **`experiments/`** | **The active research frontier** — 49 read-only diagnostics, one per scientific question (`diag_*.py`), plus `DELTAAGE_LAB_NOTEBOOK.md`. Nothing here touches `src/`; each answers a question and writes JSON to `results/`. *(Corrected 2026-08-03: an earlier revision called this "historical, not on the forward path." That was wrong — the clock-density, label-provenance and harmonization-gain work all lives here, and it is where the open questions are currently being resolved.)* |
| **`results/`** | Every diagnostic's JSON output, plus the written reports (`STEP6_FULL_REPORT.md`, `STEP6_REPORT.md`, `DAGE_LEDGER.md` + `dage_ledger.csv`) |
| **`scorecard/`** | Metric snapshots — `baseline.json`, `A_xdonor.json`, `gc2_A_keep_hff.json`, … Each is one 6-fold measurement; comparisons are always snapshot-vs-snapshot |
| repo root | `scorecard.py` (the grading tool), `audit_metrics.py` (`MIN_PASS_RATE`, `bar_verdict`, `sensitivity_multiplier`), `retrain_stage1.py` (retrain-only path — **reuses shards, cannot see a data-config change**), `run_step6_arm.sh`, `CHANGES.md` (the append-only log) |

### Generated, not tracked

`cellfate_loocv_<donor>[_arm{A,B,C}]/` — fold builds (~260 MB each). `runs/` — older builds.
`loocv_results/`, `diag_dump/`, `diag_logs/`, `*.zip`. All gitignored.

---

## 4. The core data flow

### Phase 1 — ETL (`data/build_dataset.py::run`)

```
sources.plan()  →  work list of CellChunks
      │
      ├─ load_or_fit_panel()     fit the 2000-HVG panel on the first N chunks, then FREEZE it
      ├─ build_clock()           load the frozen clock JSON
      ├─ fit_harmonizer()        control-anchor each dataset + fit the Gill projection
      │
      ▼  for each chunk  (resumable — ProgressTracker skips completed ids)
   ┌─────────────────────────────────────────────────────────────┐
   │ fetch → QC → normalise → signatures → fate labels           │
   │ → cell-cycle score → CLOCK → ΔAge (control-relative)        │
   │ → age_label_policy → panel matrix → encode perturbation     │
   │ → assemble Samples                                          │
   └───────────────────┬─────────────────────────────────────────┘
                       ├─→ shards/<chunk>.parquet
                       ├─→ manifest_parts/<chunk>.parquet
                       └─→ _cc_cache/<chunk>.npz   (cell-cycle sidecar)
      │
      ▼  after every chunk
   consolidate manifest  →  make_splits()  →  splits/<regime>.json
      │
      ▼  TWO-PASS deconfounder (needs splits, so it cannot run per-chunk)
   pass 1: fit ΔAge ~ a·cc + b on the primary regime's TRAIN cells
   pass 2: apply that one transform to EVERY shard, re-centre on controls,
           [arm C only: permute labels], rewrite y_age
      │
      ▼
   fit scalers on TRAIN only  →  scalers.json  →  dataset_summary.json
```

**Two masks, deliberately distinct** — conflating them invalidated a 10-hour experiment:

| mask | question | governs |
|---|---|---|
| `age_mask` | may the age head **train** on this cell? | which cells the loss uses |
| `deconfound_mask` | is this cell's ΔAge **computable**? | the **value** of `y_age` |

They differ by exactly the dataset-policy exclusions. Because the deconfounder and the control
re-centring use `deconfound_mask`, **withholding labels no longer changes the target variable.**

### Phase 2 — Train (`training/train_model.py::run`)

```
load_split_tensors(train/val/calib)   →  detect perturbation modality (chem | tf)
      ▼
train_ensemble()  ──  5 independently-seeded members
      ▼
inner LODO (xdonor_calib)  ──  one extra ensemble PER TRAINING DONOR
      ▼
calibrators, all fitted on cross-donor statistics:
   · temperature      · Platt on P(safe)      · conformal q      · σ-scale
      ▼
OOD: Mahalanobis on trunk features
      ▼
bundle/  =  members/*.pt + scalers + temperature + conformal + res_params + ood + xdonor_stats + meta
```

### Phase 3 — Evaluate & Infer

```
evaluate_cli   →  model vs 6 baselines on the held-out donor  →  gates
scorecard.py   →  snapshot (6 folds)  →  compare(A, B) → paired 95 % CI per metric
predictor      →  encode → ensemble → calibrate → conformal interval → RES → rank
```

---

## 5. Key components, module by module

### `models/network.py` — the shared trunk **(architecturally load-bearing)**

```python
self.trunk = nn.Sequential(_mlp_block(d_cell + d_u, latent_dim, p_drop),
                           _mlp_block(latent_dim, latent_dim, p_drop))
...
z = self.trunk(torch.cat([self.cell(x), self.pert(u, dose_time)], dim=1))
return self.cls_head(z), self.age_head(z), z
```

**One trunk feeds both heads**, and `MultiTaskLoss` sums their losses. So a change to the *age*
labels propagates into the *fate* head's representation. This single fact explains why fate ECE,
ranking, and OOD-flag rate all move together when age labels change — they are one mechanism, not
independent findings.

### `models/losses.py`

- `focal_loss` — class-balanced focal for fate.
- `huber_age_loss` — Huber over **masked** cells; returns a *differentiable* zero when a batch has
  no age-valid cell (never a detached constant, or `.backward()` breaks).
- `huber_age_window` — one Huber over an accumulation window: `Σloss / Σcells`, **never a mean of
  means** (that would weight a 1-cell batch as heavily as a 9-cell one).
- `MultiTaskLoss` — Kendall & Gal homoscedastic weighting, `total = exp(−s_cls)·L_cls +
  0.5·exp(−s_age)·L_age + 0.5·(s_cls + s_age)`. The asymmetric `0.5` follows the
  Gaussian-vs-softmax likelihood derivation.

### `data/aging.py` — the ΔAge definition

`delta_age()` returns `(ΔAge, age_mask, age_mask_reason)`. `age_label_policy()` applies three
exclusion rules **in order**, recording only the *first* reason that fires:

1. `cancer_source` 2. `dataset_policy` (`AGE_MASKED_DATASETS`) 3. `donor_out_of_clock_range`

> ⚠️ Because only the first reason is recorded, you **cannot** invert `age_mask_reason` to recover
> "excluded only by policy" — a cell that is both out-of-range and policy-masked reads
> `"dataset_policy"`. Recompute with a different `masked_datasets` argument instead.

### `training/train.py::_AgeWindow` — the age-accumulation window

With HFF's labels masked, 75 age labels sit among 33 688 cells: ~1.14 per 512-cell batch, so ~32 %
of updates carry a hard zero. `_AgeWindow` holds age cells back until a window carries **≥ k** of
them (or `w_max` batches pass), then takes one Huber over the whole window.

Three deliberate properties:

- **Triggers on cell count, not batch count.** In an unmasked arm the first batch already clears
  `k`, so the window closes at W=1 and training is **bit-identical to before** — the control arm is
  never handicapped by a mechanism it does not need.
- **Buffered cells are re-run, not replayed from a stale graph.** The optimiser steps every batch,
  so a held gradient would be computed against parameters that have since moved.
- **The window carries across the epoch boundary.** Forcing a close at each epoch's end manufactured
  one deliberately-partial window per epoch and failed its own bar.

Default `age_window_k = 1` — **off**, i.e. the pre-existing path exactly.

### `training/xdonor_calib.py` — calibrating for the regime deployment faces

For each training donor `d`: train an inner ensemble on `train − d`, early-stop on `val − d`,
collect residuals on **all** of `d`'s cells. Donors whose removal would leave under
`MIN_INNER_TRAIN_FRAC = 0.5` of the training set are **skipped** — holding out a bulk corpus
measures data starvation, not donor shift.

### `data/harmonize.py`

Control-anchors each dataset onto a shared scale, then projects into the clock's space
(`project_to_clock`). Without it the Gill donors' ΔAge sits +16…+64 years off HFF's.

### `inference/res.py`

Pure, vectorised, and the only place decision policy lives. Returns a status, not just a score:
`REJECTED_OOD`, `REJECTED_UNSAFE`, `REJECTED_NO_REJUVENATION`.

---

## 6. Every loop in the system

Grouped by phase. **Bold** loops dominate runtime.

### ETL loops — `data/build_dataset.py`

| # | loop | over | notes |
|---|---|---|---|
| L1 | source planning | every source's `plan()` | flattened into one ordered work list; chunk ids must be globally unique |
| L2 | panel fit | first `panel_ref_chunks` chunks | fits the 2000-HVG panel, then it is **frozen** for the whole project |
| L3 | harmonizer fit | all chunks (+ inner loop over `dataset_id`) | control-anchoring statistics |
| **L4** | **main chunk loop** | every chunk | the ETL body. **Resumable** — `ProgressTracker` skips done ids; failures are recorded and can continue |
| L5 | baseline census | lines within a chunk | keyed `chunk::line`, never line alone (HFF spans 45 chunks) |
| L6 | label tally | samples in a chunk | fate class counts + `n_age_labeled` |
| L7 | split write | each regime | `scaffold`, `cell_line`, `both`, `holdout`, `line_holdout` |
| **L8** | **deconfounder pass 1** | sidecars × cells | fits `ΔAge ~ a·cc + b` on TRAIN cells with `deconfound_mask` |
| **L9** | **deconfounder pass 2** | every sidecar | applies the single transform + control re-centring |
| L10 | arm-C shuffle | sorted sidecars → target cells → scatter | **global across shards**, never per chunk (chunks are timepoint-homogeneous, so a within-chunk shuffle would leave between-timepoint structure intact) |
| L11 | shard write-back | each `y_age` array | `rewrite_shard_yage` |
| L12 | scaler fit | every shard | TRAIN rows only |
| L13 | sidecar load | `_cc_cache/*.npz` | sorted, so order never depends on dict ordering |

### Training loops — `training/train.py`, `train_model.py`

| # | loop | over | notes |
|---|---|---|---|
| **L14** | **ensemble** | `ensemble_size` (5) | independently seeded; the backbone of the epistemic estimate |
| **L15** | **epoch** | `epochs` (80) | with early stopping |
| **L16** | **batch** | `train_dl` | forward → focal + masked-Huber → `MultiTaskLoss` → clip → step |
| L17 | **age window** | spans batches (stateful) | `_AgeWindow`; not a `for` loop but a genuine cross-iteration cycle |
| L18 | eval loss | validation batches | Kendall-free fixed-weight objective, so early stopping cannot be gamed by inflating σ |
| L19 | early stop | patience counter | terminates L15 |
| L20 | member outputs | inference batches (2048) | eval mode, dropout off |
| L21 | ensemble aggregation | members | `ensemble_logits` / `ensemble_probs` / `ensemble_age`. **`ensemble_probs` averages softmaxes, not logits** — averaging logits then softmaxing is a different quantity by Jensen, and mixing them was a real shipped bug |

### Calibration loops — `training/xdonor_calib.py`


| # | loop | over | notes |
|---|---|---|---|
| **L22** | **inner leave-one-donor-out** | each training donor | ⏱ **the cost centre — one extra ensemble per donor, ~6× training time per fold.** Deliberately has no shrink knob: the inner spread *is* what σ-scale calibrates |
| L23 | inner ensemble | nested inside L22 | same `ensemble_size` as deployed, or it measures a different quantity |
| L24 | conformal quantile | each level | `rank = ceil((n+1)·lvl)` |

### Evaluation & experiment loops

| # | loop | over | notes |
|---|---|---|---|
| **L25** | **LOOCV outer** | 6 donors — `run_loocv.py` | **each fold rebuilds from raw GEO**, refitting harmonisation on that fold's training donors only. ~1 h/fold ⇒ ~5 h/arm |
| L26 | scorecard folds | 6 donors | loads each bundle, measures, writes one snapshot |
| L27 | baselines | 6 estimators | mean, ridge, x-only, u-only, kNN, predict-control |
| L28 | pooled-ECE bootstrap | 4000 trials | |
| L29 | **arm loop** | A → B → C | sequential; each arm is a full LOOCV. Arms write suffixed roots (`CELLFATE_FOLD_SUFFIX`) so they do not overwrite each other |
| L30 | bar resolvability | 20 000–40 000 simulations | Monte-Carlo under the null, *before* the real run |

### Inference loops — `inference/`

| # | loop | over | notes |
|---|---|---|---|
| L31 | MC-dropout | `T` stochastic passes (default 20–50) | dropout ON at inference; epistemic uncertainty on top of the ensemble |
| L32 | batch predict | requests | |
| L33 | shard ranking | cells in a shard | top-N by RES |

### The governance loop

| # | loop | notes |
|---|---|---|
| **L34** | **pre-register → gate → run → grade → record → push** | The outermost loop, and the one that catches the others' failures. A bar without a resolvability test is **not** considered pre-registered (`REF_GROUND_RULES.md` §5b). Records are **annotated, never rewritten** |

---

## 7. On-disk artefacts

### A fold build — `cellfate_loocv_<donor>[_arm{A,B,C}]/`

```
shards/*.parquet        the cells: X, perturbation, y_cls, y_age, age_mask, age_mask_reason, …
manifest_parts/ →       consolidated manifest.parquet
splits/<regime>.json    cell_id → train | val | calib | test
scalers.json            TRAIN-fit normalisation + the deconfounder coefficient
panel.json              the frozen 2000-gene order
harmonization.json      per-dataset anchoring + clock projection
progress_tracker.json   resume state
_cc_cache/*.npz         cell-cycle sidecars (deleted after the build)
dataset_summary.json    n_samples, n_shards, n_age_labeled, split sizes, panel hash
step6_arm_census.json   which arm, k, label counts, shuffle seed
reports/<regime>.json   evaluation output
bundle/                 ← the deployable unit
   members/member_{0..4}.pt   config.yaml   meta.json   metrics.json
   scalers.json  temperature.json  conformal.json  res_params.json
   ood/mahalanobis.npz         xdonor_stats.npz
```

### `scorecard/<tag>.json`

One 6-fold measurement: per-fold `dage_mae_model`, `dage_mae_ridge`, `rank_model_dage`,
`fate_prauc`, `fate_roc`, `fate_ece`, `fate_ece_platt`, `conformal_coverage`, `level_shift_model`,
`ood_flag_rate` — plus environment provenance. **Comparisons are always snapshot-vs-snapshot**, so a
result survives the builds being overwritten.

---

## 8. The scientific-method layer

This is as much a part of the architecture as the network, and most of the recent code exists to
serve it.

### Bars are registered before runs

`audit_metrics.py` provides `bar_verdict()` and `MIN_PASS_RATE = 0.95`. A bar must come with a
**resolvability simulation**: *would a system that meets the intent exactly actually clear this?*
A bar that a correct system fails is measuring sample size, not quality — which happened twice on
`fate_ece` and is why the rule exists. Registered bars live in `tests/test_bars_resolvable.py`.

### Comparisons are paired across folds

`scorecard.py compare` reports, per metric, the paired 95 % CI across 6 donor folds. Sensitivity is
set by **consistency**, not size: `sensitivity_multiplier(6) = t(.975,5)/√6 ≈ 1.049`, so
**MDE ≈ 1.05 × SD** of the per-fold difference. A change that helps some folds and hurts others can
be large in the mean and still read as noise — the per-fold column must be checked before trusting
any aggregate.

### Nulls are only informative when powered

A CI containing zero is *absence of evidence*. It licenses a conclusion only when
**MDE ≤ Δ\***, the smallest effect worth acting on. And because SD is itself an estimate at n=6, the
χ² interval on σ must be checked before claiming a study was powered.

### Equivalence claims need a margin

"A is like B" cannot be established by a CI containing zero. It requires a pre-registered margin and
a TOST (the 90 % CI inside ±Δ_eq) — with the margin **derived** from a measured quantity, never
widened until it passes.

### Guards that fail loudly

- **Bit-identity gates** — a change that should record and not compute must leave `y_age`
  bit-identical (`max|Δ| = 0.00e+00`), asserted row-exact.
- **Pre-flight gates** — `plan_tests/step6_preflight.py` builds both arms at reduced scale and
  refuses the real run unless the arms differ in the intended way **and only** in that way.
- **Arm census** — every fold writes what arm it ran, its label counts, and any shuffle seed.
- **Mutation testing** — guards are validated by re-injecting the exact bug they exist to catch.

---

## 9. Entry points

```bash
# Full LOOCV, one arm, build + train + snapshot           (~5 h)
./run_step6_arm.sh A                 # control
./run_step6_arm.sh B                 # HFF age labels masked
./run_step6_arm.sh C 0               # HFF age labels shuffled, seed 0

# One fold only
python local_runners/run_multi_local.py "D:\GSE242423" "D:\Gill"

# Grade and compare
python scorecard.py snapshot --tag <tag>
python scorecard.py compare <tag_a> <tag_b>

# Gates and bars
python plan_tests/step6_preflight.py         # refuses a run that would waste hours
python plan_tests/register_*_bar.py          # register a bar BEFORE its run

# Hydra CLIs
python scripts/build_dataset.py  |  scripts/train.py  |  scripts/evaluate.py  |  scripts/serve.py

# Tests
python -m pytest -q                          # ~780 tests
ruff check src/ tests/ scripts/ plan_tests/  # the CI lint scope
```

> ⚠️ **`retrain_stage1.py` reuses existing shards and has no build step.** It is correct for
> training-path changes only. A change to any *data* config (masks, policies, the clock) will not
> reach it, and the run will silently produce two identical arms.

---

## 10. Cross-cutting invariants

| invariant | why it exists |
|---|---|
| **The gene panel is fit once, then frozen** | feature order must be stable across every fold and every stage |
| **The clock is validated, never refitted** | refitting to improve our numbers is fitting the test |
| **Scalers and the deconfounder are fit on TRAIN only** | applying the identical transform to val/test is what keeps them leak-free |
| **Calibrators are fit cross-donor** | deployment is out-of-donor; an in-distribution calibrator is measuring the wrong regime — decisive for ΔAge (~4 yr in-distribution vs ~14 yr out-of-donor) |
| **Determinism is engineered** | `set_global_seed` + `CUBLAS_WORKSPACE_CONFIG`; guards assert *bit-identical*, a far sharper test than "the CI includes zero". Demonstrated: a full rebuild three weeks later reproduced all six folds exactly |
| **One change per experiment** | violated *by the design* once, which cost a 10-hour run |
| **Every `y_age` consumer gates on `age_mask`** | withheld cells carry a value nobody may read |
| **Records are appended, never rewritten** | `CHANGES.md` and `plans/` keep wrong claims visible with corrections beside them |
| **Chunk ids are globally unique; the manifest is keyed `chunk::line`** | `cell_line` is *not* unique across chunks — HFF spans 45 |
| **`cell_id` is NOT unique across the dataset** | the same cell recurs under different perturbations. Keying a dict on it silently deduplicates — this turned a 7 062-row check into a 1 024-row one |
