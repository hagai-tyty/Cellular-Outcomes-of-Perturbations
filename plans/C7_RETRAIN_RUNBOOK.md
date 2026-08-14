# C-7 RETRAIN RUNBOOK — retrain → snapshot → compare

**Why:** every performance number this project would ship on is **pre-C-7**. The newest snapshot,
`scorecard/gc2_D_stratshuffle_hff_s0.json`, is dated **2026-08-07**; C-7 was adopted after it. So
the fate head's ROC 0.983 / PR-AUC 0.992, RES, and all four Stage 1 guards were measured against
labels that C-7 subsequently changed — N2's ΔAge masked, five degenerate columns rejected, HFF
day-14 fold spread cut from 16.671 to 3.686 yr.

**Until this runs, RES is *unknown*, not *unchanged*, and Stage 5 cannot be assessed.**

---

## 0. A bug fixed first, because it would have wasted the whole run

`run_multi_local.py` built its `DataConfig` **without `bulk_integrity_gate`**, so it defaulted OFF.
Running the retrain as it stood would have spent hours producing **pre-C-7 labels while every log
line said C-7** — the same inert-flag failure that cost the first C-7 build (`e6fc183`).

Now wired: a module-level `BULK_INTEGRITY_GATE` (default `False`, so no other runner changes
behaviour), set from `CELLFATE_BULK_GATE` by `run_loocv.py`, which **prints the gate state before
building**. Pinned by `tests/test_c7_reaches_the_retrain.py` (12 static checks, no build required).

---

### 0b. RUN 1 (2026-08-15) WAS INVALID — the same defect, a third time

The §0 fix was necessary and not sufficient. The flag reached `DataConfig`, `run` applied it to the
injected sources, the header printed `bulk_integrity_gate = ON` — and **the gate still did nothing.**

`run_multi_local.py` called `gill.plan()` to list donors ~30 lines *before* the `DataConfig`
carrying the flag existed. `plan()` reads and **caches** the matrix, and the gate's screen lives
inside that read (`sources.py:_load`, which early-returns on its cache). By the time
`apply_source_flags` set the flag, the unscreened 124-column matrix was already cached forever.

Measured, all six folds:

| | required | run 1 produced |
|---|---|---|
| cells | 42,600 | **42,605** |
| ΔAge labels masked | 19, all N2, `no_control_baseline` | **0** |
| Gill columns | 119 | **124** |

`n_age_labeled` came out **equal to** `n_samples` — every cell carried a label, which is the
unambiguous signature of "nothing was masked" regardless of how cells are counted.

**Fixes (both shipped):**
1. `bulk_integrity_gate` is now a **property** on `GillReprogrammingSource`; its setter drops the
   cached read when the value changes (idempotent when it does not). The flag now means the same
   thing no matter when it is set.
2. `run_multi_local.py` sets the gate on the source **before** `plan()`, so the donor list also
   comes from the gated corpus.

**And the check moved into the run.** Three C-7 failures have now been "the flag was on and nothing
happened", each caught only by inspecting artefacts afterwards. `run_multi_local.py` now **aborts**
if the gate is on and either no bulk column was rejected or no ΔAge label was masked. §3 below is
now a confirmation, not the only line of defence.

---

## 1. Choose the suffix — and do NOT reuse `_c7`

The existing `cellfate_loocv_*_c7` folds are **dataset-only** builds, and every recorded analysis of
this arc reads them: `stage3a_diagnose`, `stage3a_bis_resolvability`, `stage3a_regime_e`,
`p3_progress`, `p4_two_heads`. `run_multi_local` starts each fold with `shutil.rmtree(ROOT)`.

**Retraining into `_c7` would delete and rebuild the exact folds those results were computed on.**
The data *should* come out identical, but "should" is not a guarantee, and this arc has already
found three errors in one measurement.

**Use `_c7t`** (C-7, trained). `_c7` stays frozen as the basis of the recorded analyses.

---

## 2. Retrain — arm A, gate ON

Arm **A** is the normal arm (`AGE_MASKED` empty, no shuffle). Arms B/C/D are the label-ablation
arms from Stage 1.5.3 step 6 and are **not** what ships.

```bash
cd /d/cellfate-rx && CELLFATE_FOLD_SUFFIX=_c7t CELLFATE_BULK_GATE=1 \
  /d/.venv-cellfate/Scripts/python.exe local_runners/run_loocv.py "D:\GSE242423" "D:\Gill" --arm A
```

**Check the first ten lines of output before walking away:**

```
[C-7] bulk_integrity_gate = ON   (CELLFATE_BULK_GATE=1)
```

If it says `off`, stop — the run is worthless.

**That header is necessary, not sufficient — run 1 printed it and was still pre-C-7.** The line
that actually proves the gate bit appears at the end of fold 1's build:

```
[C-7] gate BIT: rejected 5 bulk column(s) [...]; 19 ΔAge label(s) masked of 42600 cells
```

If the gate did not bite, the run now aborts on fold 1 within minutes rather than completing.

**Cost:** six full builds plus training. The runner's own docstring says *"a few hours; run it
overnight."*

---

## 3. Verify the retrain used C-7 labels — before trusting any metric

Cheap, decisive, and it catches a gate that was on in the log but inert in the build:

- `n_samples` should be **42,600**, not 42,605 — exactly the five rejected columns.
- Donor **N2** should have **19 cells masked** for ΔAge in every fold, with reason
  `no_control_baseline`.
- The `_c7t` ΔAge labels should match the frozen `_c7` folds. **If they differ, stop and find out
  why before reading a single metric.**

---

## 4. Snapshot

```bash
cd /d/cellfate-rx && CELLFATE_FOLD_SUFFIX=_c7t \
  /d/.venv-cellfate/Scripts/python.exe scorecard.py snapshot --tag c7_A_keep_hff
```

`scorecard.py` honours the same `CELLFATE_FOLD_SUFFIX`, so it reads the trained folds.

---

## 5. Compare against the RIGHT baseline

**Not** the newest snapshot. `gc2_D_stratshuffle_hff_s0` is arm **D**, a *stratified-shuffle
ablation* — comparing against it would confound C-7 with a destroyed label pairing.

The correct comparator is the same arm, pre-C-7: **`gc2_A_keep_hff`** (2026-08-03).

```bash
cd /d/cellfate-rx && /d/.venv-cellfate/Scripts/python.exe scorecard.py compare gc2_A_keep_hff c7_A_keep_hff
```

Same arm, same k, same seeds, same regime — **the only difference is C-7.** That is what makes the
diff interpretable.

---

## 6. How to read the diff — decided in advance

| observation | reading |
|---|---|
| guards move but every paired CI includes 0 | C-7 changed labels **without** degrading the model; Stage 1's guard record survives the relabel |
| a guard **degrades** with a CI excluding 0 | C-7's cleaner labels cost real performance — record it, do **not** revert the labels to recover a number |
| a guard **improves** with a CI excluding 0 | the pre-C-7 number was inflated by the degenerate control |
| **RES** | currently **unknown**. Whatever it reads is its **first measurement** under these labels, not a change from a known value — say so |

**The guard record legitimately restarts.** C-7 moves labels, so the four-run `+0.000` history does
not carry across, and that was stated when C-7 was adopted.

---

## 7. What this does NOT do

It does not revisit Stage 3 — the forward tool is closed pending a second dense single-cell line,
and no retrain changes that. It does not touch the harmonizer. It does not settle Stage 2's
per-donor-offset premise, which M-E4 left inconclusive. And a good scorecard is **not** a licence to
ship: Stage 5 is a separate decision that this run merely makes *assessable*.
