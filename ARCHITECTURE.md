# CellFate-Rx — Architecture

**A map of the system as it actually is**, re-read from the code on **2026-08-18** (Stages 12–18).
Where this document and `README.md` disagree, this one is right — the README predates most of the
pipeline.

> **Read §11 first if you want the scientific status rather than the machinery.** The code below is
> in good order; what it has established is a much shorter list, and §11 is the honest version.

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
- [8a. Step 6 — the arm experiment](#8a-step-6-the-arm-experiment)
- [9. Entry points](#9-entry-points)
- [10. Cross-cutting invariants](#10-cross-cutting-invariants)
- [**11. Where the project actually stands — 2026-08-18**](#11-where-the-project-actually-stands--2026-08-18)
- [12. The data — every corpus, and what it can and cannot answer](#12-the-data--every-corpus-and-what-it-can-and-cannot-answer)
- [13. The result ledger — Stages 10–18](#13-the-result-ledger--stages-1018)
- [14. Failure modes this project has actually hit](#14-failure-modes-this-project-has-actually-hit)
- [15. Reading the experiments directory](#15-reading-the-experiments-directory)
- [16. If you are picking this up cold](#16-if-you-are-picking-this-up-cold)

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
| **`src/cellfate/data/`** | The ETL. `build_dataset.py` is the orchestrator; `sources.py` the corpus adapters; `aging.py` the ΔAge definition and label policy; `harmonize.py` the cross-modality alignment; `proliferation.py` the cell-cycle deconfounder; **`integrity.py` the bulk-sample gate (Change C-7)**; plus `qc`, `normalize`, `labels`, `signatures`, `perturbation`, `splits`, `chunking`, `clock_fit`, `assemble` |
| **`src/cellfate/models/`** | The network. `network.py` (`CellFateNet`: shared trunk → two heads), `encoders.py` (cell / chem / TF), `heads.py`, `losses.py` (focal, masked-Huber, Kendall–Gal `MultiTaskLoss`) |
| **`src/cellfate/training/`** | `train_model.py` (orchestrator), `train.py` (loops + `_AgeWindow`), `dataset.py` (shard→tensor), `xdonor_calib.py` (**inner LODO**), `calibrate.py`, `conformal.py`, `ood.py`, `bundle.py`, `metrics.py` |
| **`src/cellfate/evaluation/`** | `evaluate_cli.py` (gates), `baselines.py` (mean / ridge / x-only / u-only / kNN / predict-control), `metrics.py`, `regimes.py`, `report.py`, `external_validation.py`, `data.py` |
| **`src/cellfate/inference/`** | `predictor.py` (loads a bundle; **applies the fate Platt calibrator** at line ~176), `res.py` (**RES**), `conformal.py`, `ood.py`, `encode.py`, `schema.py` (request/response), `service.py` (CLI + optional FastAPI), **`dage_calibration.py` (Stage 14 — the ΔAge scale factor, applied at the reporting boundary and nowhere else)** |

### Configuration, entry points, and governance

| path | what lives here |
|---|---|
| **`configs/`** | Hydra config tree — `config.yaml` composes `data/`, `model/`, `train/`, `infer/`, `eval/`. **`configs/clocks/`** holds the frozen clocks (Fleischer RNA; Horvath skin&blood 2018 and multi-tissue 2013 for methylation cross-checks). `configs/panels/` holds the frozen gene-panel order |
| **`scripts/`** | Thin Hydra CLIs: `build_dataset.py`, `train.py`, `evaluate.py`, `fit_clock.py`, `serve.py` |
| **`local_runners/`** | The drivers actually used day to day. **`run_multi_local.py`** = one full fold (build → train → evaluate → bundle). **`run_loocv.py`** = rotates all 6 donors. Plus `build_c7_folds.py` (dataset-only C-7 builds), **`recalibrate_folds.py` (Stage 16 — re-fits the shipped fate calibrator on the HARD class without retraining)**, `run_local.py`, `run_fate_local.py`, `evaluate_only.py`, `diag_harmonize.py`, `show_ui.py` |
| **`plans/`** | The project's decision record. `00_START_HERE.md`, `MASTER_PLAN.md`, `REF_GROUND_RULES.md` (the rules everything is graded against), `REF_ARCHITECTURE.md`, `REF_DATA_STRATEGY.md`, and one file per stage. `plans/archive/` holds superseded drafts — kept, never rewritten |
| **`tests/`** | **80 files, 1583 tests**. Unit tests plus **registered-bar tests** (`test_bars_resolvable.py`) and invariance guards (`test_ci_deconfounder_arm_invariance.py`, `test_c5c_age_accumulation.py`, `test_arm_c_label_shuffle.py`, `test_arm_d_stratified_shuffle.py`, `test_results_paths.py`) |
| **`plan_tests/`** | Scripts a *plan* requires: pre-registered bars (`register_*_bar.py`, incl. `register_arm_c_bar.py` / `register_arm_d_bar.py` / `register_gc_step2_bar.py`), pre-flight gates (`step6_preflight.py`, `armd_intrinsic_preflight.py`), bit-identity verifiers (`verify_age_mask_identical.py`, `verify_stage1_5.py`, `verify_1a.py`), smoke tests |
| **`experiments/`** | **The active research frontier** — **86** read-only scripts, one per scientific question (`diag_*.py`, `test*.py`, and since 2026-08-07 the `repro_*.py` reproduction checks), plus `DELTAAGE_LAB_NOTEBOOK.md`. Nothing here touches `src/`; each answers a question and writes JSON to `results/`. *(Corrected 2026-08-03: an earlier revision called this "historical, not on the forward path." That was wrong — the clock-density, label-provenance and harmonization-gain work all lives here, and it is where the open questions are currently being resolved.)* *(Count corrected 2026-08-07: the figure read 49 when there were 45 files; it is now an actual `ls` count, not an estimate.)* **`repro_*.py` re-run a recorded result against a current build and grade it against a pre-registered bar — `repro_test7_res_armA.py` (RES/ranking) and `repro_hff_signature_armA.py` (HFF ΔAge trajectory). They import the original scripts UNMODIFIED and redirect only the run directory, so the arithmetic is identical on both sides.** |
| **`results/`** | Every diagnostic's JSON output, plus the written reports (`STEP6_FULL_REPORT.md`, `STEP6_REPORT.md`, `DAGE_LEDGER.md` + `dage_ledger.csv`) |
| **`scorecard/`** | Metric snapshots — `baseline.json`, `A_xdonor.json`, `B_fatecal*.json`, the step-6 arms `gc2_{A,B,C,D}_*`, and the current line: **`c7_A_keep_hff`** (C-7 gate on, `_c7t` folds) → **`c7t_stage12`** (Stage 12 key fix, `_s12`) → **`c7t_stage16`** (hard-label calibrator, `_s16`). Each is one 6-fold measurement; comparisons are always snapshot-vs-snapshot |
| repo root | `scorecard.py` (the grading tool), `audit_metrics.py` (`MIN_PASS_RATE`, `bar_verdict`, `sensitivity_multiplier`), `retrain_stage1.py` (retrain-only path — **reuses shards, cannot see a data-config change**), `run_step6_arm.sh`, `CHANGES.md` (the append-only log) |

### Generated, not tracked

`cellfate_loocv_<donor><suffix>/` — fold builds (~260 MB each), selected by `CELLFATE_FOLD_SUFFIX`.
Suffixes in use, oldest first: `_armA`/`_armB` (step 6), `_c7` (dataset-only), **`_c7t`** (C-7
trained), **`_s12`** (Stage 12 `cell_id` fix), **`_s16`** (Stage 16 hard-label calibrator —
hardlinked from `_s12` with only `bundle/temperature.json` replaced). `runs/` — older builds.
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
           [arm C: global label permute | arm D: within-(line,time) permute], rewrite y_age
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

### `data/integrity.py` — the bulk-sample gate (Change C-7, ships OFF)

**Why it exists.** `GSE165176` contains columns that are **not transcriptomes**.
`N2_Fib_Sendai_Exp2` — donor N2's day-0 control, and therefore N2's entire ΔAge zero-point — is
nearly a constant vector: `min = median = mean = 11.490`, a log2 range of **1.74** where sound
controls span 13–15, and a linear library **68×** the cohort. The clock reads it as **98.65 yr**
for a donor of age **0**.

**Why one bad column mattered so much.** The day-0 `_Fib_` sample is `is_control`, so it is *two
things at once*: that donor's zero-point, **and** one of the five or six controls `sigma_gill` is
fitted on. `sigma_gill / sigma_hff` is the gain applied to **HFF's** labels — 99.7 % of the
age-labelled corpus — in every fold that does not hold N2 out.

**Two conditions, justified by units rather than by this cohort's quantiles:**

| | |
|---|---|
| **G1** library | the matrix is RPM, so a sound column's linear values sum to ≈ 1e6 *by definition*. Accept `[1e5, 1e7]` |
| **G2** dynamic range | a real transcriptome spans orders of magnitude; require `log2(max) − log2(min) ≥ 8` |

On the 124 Gill columns these reject **exactly 5** with **0** false positives, each condition
independently rejecting all five.

**Where it acts.** In `GillReprogrammingSource._load` — the single place the matrix is read — so
`plan()` and all **three** `src.fetch` call sites are covered by one edit. Missing the third
(`fit_harmonizer`) would leave the degenerate column inside `sigma_ref`, which is the entire
defect.

**What it forces downstream.** Rejecting a control can leave a line with **no zero-point**, so
C-7 ships with two companions that must land with it:

* **`age_label_policy` rule 4 `no_control_baseline`** — a `cell_line` with zero admissible
  controls **anywhere in the corpus** has undefined ΔAge and is masked. GLOBAL, never
  chunk-local: `_control_baseline` *also* falls back when a line has no controls *in this chunk*
  (Stage 1.5's Group E), and conflating the two would block C-7 for the wrong reason.
* **`assert_no_unmasked_fallback` (bar B2′)** — no line may reach the fallback *and* keep its
  label. Guards **both** of `_control_baseline`'s call sites, including the S4 re-centring in
  `recenter_on_control_arrays`, which previously passed no census and so had an invisible
  fallback.

**One flag** (`DataConfig.bulk_integrity_gate`), applied by `apply_source_flags` from **both**
`build_sources` and `run` — because callers inject sources, and a flag set only on the
construction path silently never reaches a production build.

**Measured effect** (six dataset-only folds, `_c7`): HFF's day-14 ΔAge spread across folds falls
from **16.671 → 3.686 yr**, and the five contaminated folds move from −22…−24 to −4.6…−8.3, i.e.
**toward the N2 fold, which barely moves** because it already excluded the bad control. The
residual 3.686 yr is what an `n = 1` unreplicated control per donor predicts — Stage 1.5's D2,
owned by Stage 6.

---

### `data/harmonize.py`

Control-anchors each dataset onto a shared scale, then projects into the clock's space
(`project_to_clock`). Without it the Gill donors' ΔAge sits +16…+64 years off HFF's.

### `inference/res.py`

Pure, vectorised, and the only place decision policy lives. Returns a status, not just a score:
`REJECTED_OOD`, `REJECTED_UNSAFE`, `REJECTED_NO_REJUVENATION`.

**RES is identically zero on every fold, and Stage 15 attributed why.** Of the four factors in
`RES = φ(S)·S^k·g(R_eff)·exp(−λ·P_loss)`, three cannot be zero on real inputs — φ is a sigmoid,
`S^k > 0`, and `λ` ships at **0.0** so the loss term is exactly 1. The zero is **`g(R_eff)` alone,
for 119 of 119 cells**: `R_eff = max(0, −(µ_age + z·σ_age))` and **σ_age is 2.0–4.5× larger than
|µ_age|**, so no cell ever earns confident-rejuvenation credit. The closest miss across all folds
is **+2.0 yr**. It is not a bug — it is the formula correctly reporting that there is no confident
rejuvenation to credit, and it is **over-determined**: the OOD and safety gates independently zero
most cells too. See `experiments/diag_stage15_res_zero.py`.

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
| L10 | age-label shuffle (arm C/D) | sorted sidecars → target cells grouped by stratum → permute within each | **global across shards**, never per chunk. Arm C = one global stratum; arm D = one stratum per `(cell_line, time_h)`, so the between-timepoint trajectory survives and only within-timepoint pairing dies |
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
| L29 | **step-6 arm loop** | A → B → C → D | sequential; each arm is a full LOOCV (~5 h). Arms write suffixed roots (`CELLFATE_FOLD_SUFFIX`) so they do not overwrite each other. See [§8a](#8a-step-6-the-arm-experiment) |
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

### A fold build — `cellfate_loocv_<donor><suffix>/`

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
   recalibration.json          ← _s16 only: Stage 16 provenance (old/new Platt coefficients)
```

### `scorecard/<tag>.json`

One 6-fold measurement: per-fold `dage_mae_model`, `dage_mae_ridge`, `rank_model_dage`,
`fate_prauc`, `fate_roc`, `fate_ece`, `fate_ece_platt`, `conformal_coverage`, `level_shift_model`,
`ood_flag_rate` — plus environment provenance. **Comparisons are always snapshot-vs-snapshot**, so a
result survives the builds being overwritten.

**`scorecard.py` was materially wrong until Stages 13 and 17, and both defects INVERTED decisions
rather than degrading numbers.** What it does now:

| direction | meaning | metrics |
|---|---|---|
| `lower` / `higher` | monotone | MAE, ranking, PR-AUC, ECE, width |
| **`abs`** | judged on the **per-fold magnitude**, with the signed mean kept as a never-judged context row | `level_shift_*` |
| **`target`** | judged on **distance to a per-fold target** read from the fold itself | `conformal_coverage` → `conformal_level` |
| `neutral` | reported, never judged | RES rows, `ood_rate`, `n_cells` |

* **Stage 13** — `abs()` was applied to the *aggregate*, so a per-donor bias whose sign varies by
  donor measured **cancellation**, not error: it printed **0.230** for a shift whose true mean
  magnitude is **12.72 yr**. And the paired CI was built on signed values, so `−28 → −22` read as an
  increase. **12 of 20 past verdicts changed; 8 of those were comparisons in which a SHUFFLE CONTROL
  had scored `ACCEPT (better)`.**
* **Stage 17** — `conformal_coverage` was `("higher", …)`, so **widening every interval until
  nothing escaped would have scored ACCEPT**. 4–5 of 6 folds are over-covering, so this was live.
  Also added the **`b/w` fold-direction tally** (`*` = unanimous across ≥4 folds), because the
  header told readers to "check the per-fold column" and only one metric ever printed one. It is
  **context, never a verdict** — with 5 folds the best possible sign-test p is 0.0625, so no
  p-value is printed.
* **Aggregate columns are averaged over the paired fold set**, so `col_B − col_A == mean diff`
  exactly, on every row. That identity was false for 13 of 18 rows before Stage 13 and is the
  single assertion that pins both fixes.

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
  refuses the real run unless the arms differ in the intended way **and only** in that way. Arm D
  adds `plan_tests/armd_intrinsic_preflight.py`, an *intrinsic* single-build gate (see §8a).
- **Arm census** — every fold writes what arm it ran, its label counts, and any shuffle seed.
- **Mutation testing** — guards are validated by re-injecting the exact bug they exist to catch.

---

## 8a. Step 6 — the arm experiment

Step 6 (a.k.a. Stage 1.5.2 "G-c step 2") is the largest experiment the system runs, and much of the
recent `data/` and `training/` machinery exists to serve it. **The question:** HFF supplies 33 613 of
the project's 33 688 ΔAge training labels (99.7 %) and they carry an artefact signature — do they
*help* the age head, or is the model learning artefact from them? It is answered not by one run but by
a family of **arms**, each a full 6-fold LOOCV that changes exactly one thing about HFF's labels and
is graded by `scorecard.py compare` on the held-out donor.

| arm | what it does to HFF's labels | isolates | driver |
|---|---|---|---|
| **A** — control | nothing (all 33 688 age-valid) | the baseline | `run_step6_arm.sh A` |
| **B** — mask | withhold all HFF labels (75 remain) | *does removing them change ΔAge MAE?* | `run_step6_arm.sh B` |
| **C** — global shuffle | permute HFF's labels across all cells | *volume/regularisation vs information* | `run_step6_arm.sh C <seed>` |
| **D** — stratified shuffle | permute HFF's labels **within `(cell_line, time_h)`** | *between-timepoint trajectory vs within-timepoint cell signal* | `run_step6_arm.sh D <seed>` |

Each arm is one line in `DataConfig` — `AGE_MASKED_DATASETS` (B), `age_shuffle_datasets` (C/D),
`age_shuffle_strata` (D). Arms A/C/D all keep the full label **count** (`age_window_k = 4` in every
arm so the mechanism is constant); only what the labels *mean* changes. The metric that discriminates
is **`rank_model_dage`** — can the age head still *order* perturbations?

**Why the arms are a ladder.** Each rules out one mundane explanation that the previous could not:

```
A 0.95  ── true labels
   │  B: mask 99.7% → MAE inconclusive (underpowered), but ranking drops
C 0.58  ── global shuffle: ranking collapses ⇒ labels are INFORMATIVE, not just volume
D 0.61  ── stratified shuffle: still collapses (91% of the way to C)
   ⇒ the exploitable structure is WITHIN-timepoint & cell-level, NOT the day trajectory
      ⇒ a day-level systematic artefact is rejected
```

Combined with Stage 1.5.5 (the within-timepoint signal is **not** identity or sequencing depth), the
mundane-explanation space is largely closed. **What remains open** — and is *not* settled by any arm —
is whether the surviving cell-level signal is real rejuvenation or clock noise; that is the
clock-density thread (Stage 1.5.6), not another shuffle. Full write-up: `results/STEP6_FULL_REPORT.md`.

**Two design rules the arm experiment forced into the code**, both now load-bearing invariants:

1. **`age_mask` vs `deconfound_mask`** (§4). Arm B's first run was void because withholding labels
   also moved `y_age` itself — the deconfounder fit on `age_mask`. Now it fits on `deconfound_mask`,
   so the arms differ in exactly one thing.
2. **Intrinsic validation, not cross-build.** Arm D's shuffle must permute *within* each timepoint so
   the trajectory survives. This is checked inside **one** build (`armd_intrinsic_preflight.py`:
   per-stratum ΔAge multiset preserved to ~1e-14) rather than by diffing two builds — build-to-build
   value comparison at reduced scale carries confounds unrelated to the transform under test.

---

## 9. Entry points

```bash
# Full LOOCV, one arm, build + train + snapshot           (~5 h)
./run_step6_arm.sh A                 # control
./run_step6_arm.sh B                 # HFF age labels masked
./run_step6_arm.sh C 0               # HFF age labels globally shuffled, seed 0
./run_step6_arm.sh D 0               # HFF age labels shuffled WITHIN (line,time), seed 0

# The current line (C-7 gate ON, fold set selected by suffix)
CELLFATE_FOLD_SUFFIX=_s16 CELLFATE_BULK_GATE=1   python local_runners/run_loocv.py "D:\GSE242423" "D:\Gill" --arm A

# Re-fit the fate calibrator on existing bundles, no retraining   (~2 min/fold)
python local_runners/recalibrate_folds.py          # _s12 -> _s16

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
python -m pytest -q                          # 1583 tests
ruff check src/ tests/ scripts/ plan_tests/  # the CI lint scope
```

> ⚠️ **Pick the snapshot tag by hand.** `run_loocv.py` derives the tag it advertises from
> arm + gate only, so re-running an arm prints the tag of the *previous* run — which would
> overwrite the only comparator the new run has. See §14 trap 10.

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
| **`cell_id` IS unique — as of Stage 12 (2026-08-17)** | ⚠️ *This row previously read "`cell_id` is NOT unique across the dataset", describing the defect as if it were a design property. It was a bug.* The key omitted the chunk, so **42,600 cells carried 1,100 distinct ids** and `splits/holdout.json` held 1,100 entries for 42,600 cells. It is now `{chunk_id}:{row}`, with a **build-time guard** before `make_splits`. Builds before `_s12` still carry colliding ids and remain readable — the fix is forward-only |
| **A calibrator must be fitted against the target its consumers read** | Stage 16: the fate Platt was fitted on the **soft** `y_cls` probability while the safety gate, `fate_ece` and the served `p_identity_preserved` all read `S` as P(**hard** class = safe). It was near-perfect against soft (ECE 0.009–0.013) and badly wrong against hard (0.106–0.113) |
| **A units correction never touches a decision path** | Stage 14: ΔAge calibration is applied in `build_response` only. `res.py`'s `kappa` is a half-saturation **in years**, and the training loss's `huber_delta = 2.0` is a knee in years — rescaling upstream of either silently reinterprets it |


---

## 11. Where the project actually stands — 2026-08-18

*The machinery above is sound and heavily tested. This section is what it has and has not
established. It is deliberately blunt; `CHANGES.md` carries the full evidence for every line.*

### 11.1 What is WORKING

| | evidence |
|---|---|
| **Fate classification, within a timepoint** | stratified AUC **0.917** over **12 (safe, unsafe) pairs** from the same timepoint, permutation **p = 0.0091** (`diag_stage18_fate_beyond_day`). Real signal, very thin base |
| **The fate safety gate, after Stage 16** | sensitivity **0.275 → 0.670** with specificity **unchanged at 0.929** and false approvals unchanged at 2; `fate_ece` 0.276 → 0.182 (ACCEPT, 5/0 unanimous) |
| **Within-donor ΔAge RANKING** | `rank_model_dage` **0.942**, `rank_ridge_dage` 0.981 — stable across every change in this arc |
| **`top100` clock variant vs methylation** | MAE **7.15** against an instrument floor of **7.30**, CI spanning zero — statistically indistinguishable from the disagreement between the two methylation clocks **on these 44 samples**, and better than the predict-zero control (MAE 9.89 vs skin & blood, 11.71 vs pan-tissue). **The 7.30 yr figure is our measurement, and calling it a floor overstates it.** The two clocks are not equally suited to this tissue: Horvath's 2018 skin & blood clock reports MAE **2.6 yr** and r 0.91 in fibroblasts and outperforms the 2013 pan-tissue clock there, so a large share of the 7.30 is the pan-tissue clock being off-tissue rather than an irreducible limit on measuring ΔAge |
| **The instrument itself** | determinism is bit-exact; C-7 gate reproduces exactly; the scorecard's two inverted decision rules are fixed and its verdicts are re-derivable from committed snapshots |

### 11.2 What we GAVE UP ON — scrapped, with the reason

| scrapped | why | where |
|---|---|---|
| **Same-timepoint ΔAge *prediction* as a headline** | **circular.** 1,956 of the 2,000 panel genes carry clock weights; the clock's own weights reconstruct the label at ρ **0.96–0.97**, ridge reproduces it at ρ 0.96–0.99. Predicting ΔAge from expression is reading back a linear functional of the input. Per arm: under **C-7 all five arms verdict CIRCULAR** (ridge-vs-label ρ 0.965–0.995). In the earlier **pre-C-7** set the split is N2 and Y2 CIRCULAR, O1/O2/Y1 LABEL-RECOVERABLE, and N3 NOT CIRCULAR — 5 of 6 recover the label. The single dissenting arm, N3, is itself CIRCULAR under C-7, so the headline word "circular" is supported rather than compressed | `diag_clock_circularity` |
| **"early ΔAge → late ΔAge"** | partial correlation **−0.064** after controlling for donor chronological age. Donor age does all the work, is known at t=0, and needs no model | `diag_early_late_forward` |
| **The 10 % ΔAge accuracy target** | unverifiable, not merely hard. 10 % of a truth with SD 12.66 yr (pan-tissue) is MAE ≤ **1.27** yr; the 1.36 yr figure quoted earlier is 10 % of the skin & blood SD of 13.55 yr — the two clocks' SDs were mixed in one sentence, and both bounds are far below anything achievable here; the two reference clocks disagree with each other by **7.30 yr** | `diag_instrument_floor` |
| **Removing the pluripotency signature from ΔAge** | recommendation **WITHDRAWN**: in this dataset pluripotency behaves as **mediation, not contamination** (3/3 tests) — removing it deletes signal (ρ 0.770 → 0.354). **Scope this carefully.** It is a statement about a full OSKM timecourse and about our measurement path, NOT a claim that pluripotency is necessary for rejuvenation. The partial-reprogramming literature shows the opposite is achievable: Ocampo 2016 and Gill 2022 (eLife, ~30 yr rejuvenation retaining cell identity) report substantial age reversal without stable pluripotency, so the rejuvenation and pluripotency programs are separable in general even though they co-vary here | Stage 10 |
| **A learned replacement clock** | **NOT LEARNABLE** — split on all three model families | Stage 1.5.4 |
| **The sparse clock as a shipped default** | validated leave-one-donor-out on Gill; **does not transfer to HFF**, and interacts adversarially with harmonization | Stage 1.5.6 |
| **`p_unsafe` regression on GSE165177** | structurally impossible: `p_unsafe` is a fraction of *cells*, a bulk sample is already a population average, so a per-sample hard label collapses it to 0/1. `unsafe_sd_by_donor` = 0.10 / 0.00 / 0.00 | Stage 3a, `P0_void` |
| **Rescaling the ΔAge training target** | would change the loss regime as a side effect — `huber_delta = 2.0` is a knee in years, and residuals are 1.36–2.40 yr, so the quadratic share moves 43–67 % → 85–97 %. Two changes in one costume | Stage 14 pre-flight |
| **Stacking a second Platt at inference** | two Platts compose exactly into one, so it would have worked numerically while hiding that the first was fitted against the wrong target | Stage 16 |
| **The "Ranking generalizes: Spearman 0.40" headline** | **RETRACTED.** The correlation was computed over floating-point residue — RES was 0 everywhere except dust at 1e-11 | `run_loocv.py` |

### 11.3 What we are STILL TRYING to make work

| open | status |
|---|---|
| **Fate signal beyond the clock** | **the live question.** Only **7 of 70** held-out timepoints carry more than one class; on O1 and O2 the timepoint ALONE reaches PR-AUC **1.000**. The 12 within-timepoint pairs are all the evidence there is, and more mixed-timepoint data is the only way to grow it |
| **Whether `k_var = 0.5991` transfers** | ΔAge calibration is fitted on O1/O2/O3 of the transient arm — the only rows with methylation truth — and that cohort is **disjoint** from the training one. Reported alongside raw with an UNTESTED caveat |
| **RES becoming non-zero** | requires σ_age < |µ_age| for some cell. That is the same signal-vs-noise wall as ΔAge, not a formula problem |
| **`top100` vs skin & blood** | passes multi-tissue (−0.16, CI spans 0), **fails** skin & blood (+3.97, CI excludes 0). The split is quantified but unexplained |
| **A second timecourse** | `DATA_REQUIREMENT_SECOND_TIMECOURSE.md` — an open data-acquisition requirement, not an analysis one |

### 11.4 What we NO LONGER USE

| | |
|---|---|
| **`_c7t`, `_s12` fold sets** | superseded by **`_s16`** for anything fate-related. `_s12` is retained as the Stage 12 comparison baseline; `_c7t` is the pre-Stage-12 baseline. Neither should be used for new work |
| **`retrain_stage1.py`** | reuses shards and **cannot see a data-config change** — unsafe for anything touching ETL |
| **`README.md`** | predates most of the pipeline. This document supersedes it |
| **`STAGE_6_NEW_DATA.md`** | superseded by `STAGE_6_NEW_DATA_REV.md` |
| **`res_approvals` as a quality metric** | approvals are 0 everywhere; the meaningful quantity was always approvals *relative to oracle* |
| **Marginal `fate_prauc` as a headline** | inflated by an input the model was handed. Use the stratified number |

### 11.4b Checked against the published record — 2026-08-28

Every §11 claim that could be checked against the literature rather than against our own artifacts
was checked. Four held, three needed scoping, one was verified as exactly right.

| claim | verdict |
|---|---|
| **Same-timepoint ΔAge is circular** | **HOLDS, and needs no literature support.** The clock is an elastic-net linear model on log1p-CP10K expression and 1,956 of the 2,000 panel genes carry weights, so predicting ΔAge from that expression recovers a linear functional of the input. That is arithmetic. Per-arm: all five C-7 arms verdict CIRCULAR |
| **Two Platts compose exactly into one** | **VERIFIED.** Our Platt is `sigmoid(a·logit(P)+b)`, and composition on the logit is exactly a single Platt — reproduced numerically to 8e-17. It would NOT hold for a Platt applied to a raw probability |
| **Donor age does the work; partial −0.064** | **HOLDS.** Chronological age dominating a clock-derived measure is the expected behaviour, not an anomaly |
| **`p_unsafe` regression structurally impossible on bulk** | **HOLDS as stated** — the point is that a per-sample *hard label* collapses a cell fraction to 0/1, which is a property of our label construction |
| **The 7.30 yr "instrument floor"** | **OVERSTATED, now scoped.** It is our measurement on 44 samples, and the two clocks are not equally suited to fibroblasts — the skin & blood clock reports 2.6 yr MAE and beats the pan-tissue clock in that tissue |
| **Pluripotency is mediation, not contamination** | **TRUE HERE, but must not be read as general.** Partial-reprogramming work (Ocampo 2016; Gill 2022) shows substantial rejuvenation without stable pluripotency, so the two programs are separable even though they co-vary in a full OSKM timecourse |
| **"10 % target = MAE ≤ 1.36 yr"** | **ARITHMETIC MIXED.** 1.36 is 10 % of the skin & blood SD (13.55), not of the quoted pan-tissue SD (12.66, which gives 1.27). The conclusion is unaffected |

**One unstated weakness, now stated.** The clock underneath all ΔAge work is our own
reimplementation from GSE113957, and it is materially weaker than the published one: our CV MAE is
**12.27 yr** (Pearson 0.837) against Fleischer 2018's ensemble at **7.7 yr mean / 4.0 yr median**
error and r² 0.81 on the same 133 samples. They used an LDA ensemble; we used a single elastic net.
Nothing above depends on the clock being *good* — the circularity result depends only on its being
*linear* — but any claim of the form "our estimate sits on the measurement floor" is easier to
satisfy with a noisy instrument, and should be read with that in mind.

### 11.5 The one-sentence version

**The instrument is now trustworthy and most of the original claims are not.** ΔAge survives as a
*measurement* (top100 sits on the methylation floor) and as a *within-donor ranking*, but not as a
prediction; fate classification has real within-timepoint signal resting on twelve pairs; and RES,
the product's headline output, is structurally zero until the model's uncertainty falls below its
signal.

---

## 12. The data — every corpus, and what it can and cannot answer

The `DataSource` adapters live in `data/sources.py`. Only two are load-bearing; the rest are
scaffolding or unused.

| adapter | corpus | shape | role |
|---|---|---|---|
| **`GSE242423SingleCellSource`** | GSE242423 | HFF single cells, 9 timepoints, chunked into **45** blocks (~42.5 k cells at `MAX_CELLS = 5000`) | **the training mass.** Supplies train/val/calib. Every fold shares it |
| **`GillReprogrammingSource`** | GSE165176 (Gill 2022, Sendai) | 6 donors × ~20 **bulk** samples | **the held-out donors.** One donor per LOOCV fold is the test set |
| `ReprogrammingSource` | base class | — | shared parsing: `time_h`, `cell_type`, control detection |
| `SyntheticSource` | generated | arbitrary | tests and smoke builds |
| `TahoeSource`, `SciplexSource` | perturbation atlases | — | **not used in any recorded result.** Present for the chem modality, which the reprogramming line never exercises |

### Corpora used in analysis but NOT in any build

| corpus | used for | why it never enters a build |
|---|---|---|
| **GSE165177** (Gill transient arm, donors O1/O2/O3) | the **ΔAge ledger** — the only rows carrying methylation truth (68 of 90). Source of `k_var = 0.5991` | different protocol; and per-sample labels cannot supply a per-cell `p_unsafe` |
| **GSE165178 / GSE165179** | the two **methylation twins** — Horvath skin&blood and multi-tissue. They define the **instrument floor** | reference measurements, not model inputs |
| **GSE113957** (143 donors, ages 1–96) | the Fleischer clock's own training cohort | scoring the clock here is circular — measured at MAE **0.13** against a published `cv_mae` of **12.27**, a 94× gap |
| **GSE297234** (2 donors, D0 only) | a clock sanity check | no timecourse, so no fate labels |

### The five structural limits this data imposes

1. **Six donors.** Every paired comparison is n ≤ 6, so the minimum detectable mean is ≈ 1.05 × SD
   of the effect. A change that helps some folds and hurts others is invisible by construction.
2. **~20 held-out cells per fold.** `conformal_coverage` is therefore quantised in steps of ~0.05,
   and a change smaller than one cell cannot be seen. Stage 12's null is a null *at that
   resolution* — not proof that nothing changed.
3. **The held-out sets are bulk.** A bulk sample is already a population average, so any per-cell
   *fraction* (`p_unsafe`) collapses to 0/1 before it can be estimated. This is why Stage 3a's
   `P0` is void, and why sample replication cannot fix it.
4. **Fate is nearly a function of timepoint.** Only **7 of 70** held-out timepoints carry more than
   one class. This is the binding constraint on the entire fate claim (§11.3), and it is a property
   of the experimental design, not of the model.
5. **The methylation cohort is disjoint from the training cohort.** No Sendai condition carries
   methylation truth, so every cross-instrument number is transient-arm only, n = 44 — and the
   ΔAge calibration factor inherits that as an untested transfer.

---

## 13. The result ledger — Stages 10–18

*One line per stage: the question, the answer, and whether it moved `src/`. Full evidence in
`CHANGES.md`; each stage has a plan in `plans/`.*

| stage | question | answer | `src/`? |
|---|---|---|---|
| **10** | Is pluripotency contaminating ΔAge? | **No — MEDIATION** (3/3 tests). Removing it deletes signal, ρ 0.770 → 0.354. Recommendation **withdrawn** | no |
| **11** | Is the ΔAge scale error calibratable? | **Scale is a major component of dense-clock error — not the whole of it.** Variance-preserving calibration reduces raw MAE **22.69 → 10.68** but **does not reach the 7.30 yr methylation floor**. LODO `k` is stable across donors (spread 1.19×, bar was 2×), but **transfer of the reporting calibration outside its methylation cohort remains untested** | no |
| **12** | Is `cell_id` unique? | **No.** 42,600 cells carried **1,100** ids; the split map held 1,100 entries. Fixed, plus a build-time guard. The rebuild changed the split **exactly as predicted** and moved **no model metric** | **yes** |
| **13** | Does the scorecard judge `level_shift` correctly? | **No.** `abs()` applied to the aggregate measured *cancellation* — it printed **0.230** for a shift whose magnitude is **12.72 yr** — and the CI was built on signed values. **12 of 20 past verdicts changed; 8 were shuffle controls scored as improvements** | scorecard |
| **14** | Should we adopt a calibrated ΔAge? | **Yes — at the reporting boundary only.** Rescaling the *target* would also change the loss regime (Huber quadratic share 43–67 % → 85–97 %). Ships `k_var = 0.5991` **alongside** raw | **yes** |
| **15** | Why is RES zero? | **`g(R_eff)` alone, 119 of 119 cells.** σ_age is 2.0–4.5× \|µ_age\|. Not a bug, and **over-determined** by three independent gates | no |
| **16** | Why does the safety gate reject safe cells? | **Target mismatch.** Platt was fitted on the **soft** label while every consumer reads `S` as P(**hard** = safe). Fixed and verified on recalibrated folds: sensitivity 0.275 → **0.670**, specificity **unchanged** | **yes** |
| **17** | Does the scorecard judge coverage correctly? | **No.** `conformal_coverage` was `("higher", …)`, so widening every interval to 1.0 would have scored ACCEPT. Added the `target` direction and the `b/w` fold tally | scorecard |
| **18** | Is the fate head biology or a clock? | **Both — the clock is most of it.** Marginal 0.93–0.96 is largely `dose_time`; within-timepoint AUC **0.917**, p = 0.0091, on **12 pairs** | no |

### Predictions this arc got wrong, and what corrected them

Kept, because a method that records only its successes is not a method.

| predicted | actual | corrected by |
|---|---|---|
| Rescaling ΔAge is a pure units change | residuals are *comparable* to the Huber knee; the quadratic share more than doubles | Stage 14 pre-flight |
| The fate gate is plainly miscalibrated | Platt was **already applied** — the real defect was the *target* it was fitted against | Stage 16 diagnosis correction |
| Shipping Stage 16 costs specificity 0.929 → 0.821 | specificity **unchanged**; false approvals unchanged at 2 | Stage 16 verification |
| The pluripotency signature should be removed | mediation, not contamination | Stage 10 |
| D0 occupies 112 index slots | 117 across the union of chunks; 112 was one shard | Stage 12 effect |
| Only `fate_ece` would move under recalibration | `fate_ece_platt` moves per fold too (aggregate unchanged) | Stage 16 verification |

---

## 14. Failure modes this project has actually hit

*Each of these cost real time and is now guarded. They are listed because they recur.*

| # | trap | how it showed up | guard now in place |
|---|---|---|---|
| 1 | **A metric that measures sample size** | per-fold `fate_ece` on ~21 cells in 10 bins: a *perfectly* calibrated model scores 0.183 and clears a 0.169 bar only **26.9 %** of the time | `pooled_fate_ece` + `audit_metrics.bar_verdict`; a bar must be shown *resolvable* before it is registered |
| 2 | **A correlation over floating-point dust** | "Ranking generalizes: Spearman 0.40" — computed where RES was 0 everywhere except residue at 1e-11 | headline retracted in `run_loocv.py`; tests assert dust as `< 1e-9`, never `== 0.0` |
| 3 | **A guard that cannot fail** | `verify_1a` graded PASS on a warning it had itself printed | every gate now has a test that *constructs a failure* and requires it to raise |
| 4 | **A flag that silently does nothing** | the C-7 gate printed ON and did nothing **three times** — the source had cached its matrix before the flag arrived | cache-invalidating property setter; exact-match constants (`C7_EXPECT_REJECTED`, counts) asserted at run time |
| 5 | **A many-to-many join** | `diag_target_shift` joined on `cell_id`, giving a 43× explosion and entirely *plausible* fake numbers | positional pairing, plus a canary comparing a known-identical pair — it read 73.77 instead of 0 |
| 6 | **A near-zero std passing an `== 0` check** | residuals of ~1e-14 sailed through and produced a Spearman over dust | scale-relative guards, never exact-zero |
| 7 | **A shrinkage trap** | MAE improves while ordering degrades, because `k_LS = ρ·SD(y)/SD(p) < 1` | variance-matched `k` preferred for reporting; SD ratio reported beside every MAE |
| 8 | **A circular test** | scoring the clock on its own training cohort: MAE 0.13 against a published `cv_mae` of 12.27 | cohort provenance stated for every clock number |
| 9 | **A comparison spanning two changes** | arm B changed the target *and* reweighted donors 400× | "one change per experiment"; the validity precondition is checked from `git log -- src/` before a rebuild |
| 10 | **A tag that destroys its own baseline** | `run_loocv.py` derives its snapshot tag from arm + gate only, so re-running an arm overwrites the comparator | recorded in `plans/STAGE_12_CELL_ID_UNIQUENESS.md` §12.9; a fresh tag is chosen deliberately |
| 11 | **A directory glob as an audit scope** | writing a new snapshot silently enlarged Stage 13's retro-audit and broke three pinned counts hours after it shipped | scope frozen to the nine snapshots the broken rule actually judged |
| 12 | **Validating a transform by diffing two builds** | build-to-build value comparison at reduced scale carries confounds unrelated to the transform under test | validate **intrinsically within one build** — the arm-D lesson |
| 13 | **Judging a target-seeking metric directionally** | `conformal_coverage` as "higher is better": widening every interval to 1.0 would have scored ACCEPT | the `target` direction, with the target read **per fold** from the data |
| 14 | **Believing unit tests are deployment evidence** | Stage 16's fix passed 13 tests while every bundle on disk still carried the old coefficients | a stage is "implemented and tested" until artefacts are rebuilt; only then "empirically validated" |

---

## 15. Reading the experiments directory

86 read-only scripts. They are not a pile — they follow four conventions, and the conventions are
what make them auditable.

1. **Pure functions, then a printed table, then one JSON dump.** `.write_text(` appears exactly
   once per script, and a test asserts that.
2. **Bars in the docstring, before the numbers.** A `verdict_from` function applies the
   pre-registered rule *mechanically*, so nobody picks a branch by hand after seeing the result.
3. **`_RESULTS = Path(__file__).resolve().parents[N] / "results"`** — enforced by
   `tests/test_results_paths.py`, which exists because a regex rewrite once produced
   `_RESULTS / "x.json".write_text(...)` in 20 places (`.` binds tighter than `/`, so the method
   bound to the *string*).
4. **Every script has a `tests/test_<name>.py`** pinning its decision branches on constructed
   input, so a branch that never fires in production is still exercised.

### The scripts carrying current conclusions

| script | conclusion |
|---|---|
| `diag_clock_circularity.py` | ΔAge prediction is circular — ρ 0.96–0.99 |
| `diag_early_late_forward.py` | early→late ΔAge is donor age — partial −0.064 |
| `diag_instrument_floor.py` | the floor is 7.30 yr; `top100` sits on it, the dense clock is 3× outside |
| `diag_stage10_pluri.py` | pluripotency is mediation, not contamination |
| `diag_stage11_scale.py` | scale is a major component of dense-clock error; calibration reaches 10.68 against a 7.30 floor, so a residual gap remains |
| `diag_stage12_split_effect.py` | the split-composition harm, measured with no rebuild |
| `diag_stage12_rebuild_verdict.py` | §12.9's pre-registered null, applied mechanically |
| `diag_stage13_retro_verdicts.py` | 12 of 20 past verdicts changed; 8 flattered shuffle controls |
| `diag_stage14_calibration_equivariance.py` | rescaling the target is not a units change |
| `diag_stage15_res_zero.py` | RES = 0 attributed to `g(R_eff)` alone |
| `diag_stage16_safety_floor.py` | the safety gate rejects 70 % of demonstrably safe cells |
| `diag_stage18_fate_beyond_day.py` | the fate head is partly a clock; 0.917 on 12 within-timepoint pairs |

### Naming

`diag_*` — a question. `repro_*` — re-runs a recorded result against a current build and grades it
against a pre-registered bar, importing the original script **unmodified** and redirecting only the
run directory, so the arithmetic is identical on both sides. `test*.py` inside `experiments/` are
historical numbered tests (`test7_4_*`), **not** pytest files.

---

## 16. If you are picking this up cold

1. Read **§11**. It is the only section that says what is true.
2. Run `pytest -q`. 1583 tests, and they encode most of what was learned the hard way.
3. Use **`_s16`** folds. `_c7t` and `_s12` are baselines, not working artefacts.
4. Grade every change with `scorecard.py compare`, never by eye — and read the `b/w` column before
   trusting an aggregate verdict.
5. Before spending hours of compute, check the validity precondition: `git log -- src/` between the
   baseline build and now must show **exactly one** change.
6. **Stages 21–23 are fate-only** — ΔAge and RES are not prerequisites and not prospective
   targets; nothing in them depends on the aging clock (`plans/THE_PATH.md`).
7. The binding constraint is **data, not code**: 7 of 70 timepoints carry mixed fate, and 12 pairs
   is the entire evidence base for the one live claim.
