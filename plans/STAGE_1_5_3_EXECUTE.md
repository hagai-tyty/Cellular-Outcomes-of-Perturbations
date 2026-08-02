# STAGE 1.5.3 — EXECUTE: the code changes Stage 1.5.2 forces

> ## ✅ **STEPS 1–4 EXECUTED 2026-08-01. No label moved.**
>
> *The pre-registration below is unchanged. This box records what happened.*
>
> | step | shipped | gate |
> |---|---|---|
> | **1** | C-6 `age_mask_reason`, C-3 HFF donor metadata | **IDENTICAL, max\|Δ\| = 0.00e+00** |
> | **2** | **C-1** — `age_mask` can finally address HFF | **IDENTICAL, 0.0** |
> | **3** | C-2 — the clock's `age_range` carried, flag off | **IDENTICAL, 0.0** |
> | **4** | C-4 option (a) + PART B.2's 7 annotations | `res.py` untouched; **zero deletions** in `plans/` |
>
> Gate: `plan_tests/verify_age_mask_identical.py`, baseline captured **before any `src/` edit**,
> 7 chunks / 1944 cells, and it **self-tests that it can fail** before trusting a pass.
> **676 tests pass**, ruff clean.
>
> **Demonstrated live:** HFF and Gill are now separable inside one chunk —
> `hff_sc → age_mask=False, reason=dataset_policy` beside `gill_bulk → age_mask=True`.
> **G-c step 2 is runnable; it was not before.**
>
> **One deviation:** C-6 ships the TOLERANT shard reader, not the strict one — see the note in
> C-6 and in `CHANGES.md`. Steps **5–7 remain open by design**, and both policy flags are still
> off.

**Status:** 🔵 **PRE-REGISTERED — steps 1–4 now EXECUTED (see the box above); steps 5–7 open.**
**Depends on:** `STAGE_1_5_2_LABEL_ANCHOR.md` (✅ closed) and `STAGE_1_5_1_REV_FINAL.md` §11 (✅ closed).
**Blocking for:** Stage 2's premise, Stage 3's recommendation rule, Stage 5's claims.
**Not blocking for:** the fate/safety head, which consumes no ΔAge and is untouched throughout.

**Scope:** 6 `src/` changes (**PART A**) and 7 downstream plan annotations (**PART B**, with the
paste-ready text in **B.2**). One of the six is blocking, four are mechanical, and the one product
decision — **C-4 — was taken on 2026-07-31: option (a)**.

**Decisions still outstanding, and neither blocks starting:**

| | decision | needed by |
|---|---|---|
| **C-5** | which of three fixes for the age head's label starvation | **step 5** — the only one on the critical path to the retrain |
| **C-2** | whether to enable the clock-range check (masks the two neonatal donors) | a **separate experiment after step 6**; the code ships at step 3 with the flag off and needs no decision to do so |

**Steps 1–4 need no decision at all.**

| where to look | |
|---|---|
| **§1** | the 27 verified facts everything below rests on |
| **PART A** | the `src/` changes — defect, evidence, ready-to-apply code, tests, guard |
| **PART B / B.2** | the plan annotations, with the exact anchor line and the text to paste |
| **PART C** | **G-c step 2 designed in full** — arms, metric, bar, cost, pre-registered outcomes |
| **PART D–F** | change manifest, commands per step, rollback per change |
| **§3–§7** | holes I went looking for, scope limits, order of operations, verification, falsification |

> **Why this stage exists.** Stage 1.5.2 was measurement-only and kept `git diff --stat src/` empty
> for every measurement — deliberately, so no guard could move while a verdict was being decided.
> **That was correct and it left a bill.** This document is the bill, itemised, with the line number
> and the evidence for each item.
>
> **Every claim below was verified against the file at the stated line on 2026-07-31, not quoted
> from another plan.** Where a number is cited it was recomputed from the run artefacts; §1 lists
> what was checked and how.

---

## 0. The one-paragraph version

Stage 1.5.2 established that the transcriptomic clock **is not calibratable** and that HFF's ΔAge
labels are ambiguous enough to require a retrain comparison (G-c step 2). **That retrain cannot be
run today**, because `age_mask` is a function of `source` alone and Gill and HFF share the same
`source` string. Fixing that is the blocking change. The rest follow from the same finding reaching
parts of the system that still treat ΔAge as trustworthy — most importantly the deployed
`rejuvenation_efficacy_score`, which is **multiplied** by a term derived entirely from ΔAge.

---

## 1. The evidence base — what was verified, and how

Every row was checked directly. Nothing here is inherited from a previous document's summary.

| # | fact | how it was verified |
|---|---|---|
| **E1** | `age_mask` is a function of **`source` only** | read `src/cellfate/data/aging.py:219` |
| **E2** | Gill **and** HFF both report `source = "reprogramming"` | `sources.py:423` (Gill) and `sources.py:549` (GSE242423); `grep -n 'name = "reprogramming"'` returns 301/423/549 |
| **E3** | `dataset_id` **already** distinguishes them — `gill_bulk` vs `hff_sc` — and `age_mask` ignores it | `sources.py:523` and `sources.py:716` |
| **E4** | the clock's own `age_range` is **written and never read** | written at `clock_fit.py:90`; `grep -rn "age_range" src/` returns that one line and nothing else |
| **E5** | the clock's fitted range is **[1.0, 96.0]** | read from `configs/clocks/fleischer_clock.json` → `meta.age_range` |
| **E6** | HFF's chunk builder passes **no** `extra=`, so G-b's `donor_age`/`batch` never reach HFF | `sources.py:715-716` |
| **E7** | RES is **multiplied** by a ΔAge-derived term | `res.py:38-41`: `R_eff = max(0, −(mu_age + z_conf·sigma_age))`, `g = R_eff/(R_eff+kappa)`, `res = phi · S**k · g · exp(−lam·P_loss)` |
| **E8** | `REJECTED_NO_REJUVENATION` fires **purely** on ΔAge | `res.py:44-45`: `elif R_eff == 0.0` |
| **E9** | the response carries **no** ΔAge-validity field | `schema.py:28-41` — `delta_age_mean`, `delta_age_interval`, and no flag |
| **E10** | the OOD gate **cannot** know the clock is out of domain | `ood.py:24-44` is a latent-space Mahalanobis distance; it has no notion of the clock, cell type or donor age |
| **E11** | training-split composition: **33 613 HFF / 33 688 = 99.7774 %**, non-HFF **75** | recomputed from `runs/cellfate_loocv_O1/manifest.parquet` joined to `splits/holdout.json` |
| **E12** | per-line age-valid TRAIN labels in the O1 fold: N2 14, N3 16, O2 18, Y1 13, Y2 14 | same computation. (O1 is the held-out donor, hence absent) |
| **E13** | **N2 and N3 are neonatal (donor age 0)** — also outside the clock's `[1, 96]` range | GEO `GSE165176_series_matrix` `donor age` characteristic, now parsed by G-b |
| **E14** | `batch_size = 512` | `configs/train/default.yaml` |
| **E15** | the masked age loss returns a **hard zero** for a batch with no age-valid cell | `losses.py:55-57` |
| **E16** | the conformal quantile is fitted on **cross-donor** residuals, **n = 103** | `train_model.py:254` uses `xstats.abs_residuals`; `runs/cellfate_loocv_O1/bundle/xdonor_stats.npz` → `abs_residuals.shape = (103,)` |
| **E17** | the inner-LODO **skips the HFF fold**, so those 103 residuals are **Gill-only** | `xdonor_calib.py:52` `MIN_INNER_TRAIN_FRAC = 0.5`, with the reason in the comment above it; 124 non-HFF cells − 21 held-out O1 = **exactly 103** |
| **E18** | the deployed ΔAge interval half-width is **q = 34.64 yr** | `runs/cellfate_loocv_O1/bundle/conformal.json` |
| **E19** | `sigma_scale = 12.67`, `z_conf = 1.0` | `bundle/conformal.json` and `bundle/res_params.json` |
| **E20** | the **two out-of-range donors have 3.01× the ΔAge MAE**: N2 21.79 / N3 29.70 vs O1 5.39 / O2 7.54 / Y1 7.28 / Y2 14.06 | `scorecard/baseline.json` → `folds.*.dage_mae_model` |
| **E21** | **their conformal coverage is exactly 0.000**, against 0.601 on the four adults | same file → `folds.*.conformal_coverage` |
| **E22** | it is **not a sample-size effect** — N2 and N3 carry 21 test cells each, the same as O1/O2/Y2 | same file → `folds.*.n_cells` |
| **E23** | it is **not a model artefact** — the **ridge baseline** shows the same split, 23.93 vs 9.11 (2.63×) | same file → `folds.*.dage_mae_ridge` |
| **E24** | **ranking barely degrades** on those donors: 0.910 vs 0.967 | same file → `folds.*.rank_model_dage` |
| **E25** | `DataConfig` has no flag governing clock-range enforcement | `build_dataset.py:76-97` — full field list read |
| **E26** | the training loader is `shuffle=True` with no sampler | `train.py:117` `train_dl = loader(train_ds, cfg.batch_size, shuffle=True)` |
| **E27** | the deconfounder is a 2-parameter OLS, `ΔAge ~ a·cc + b` | `proliferation.py:37-39`, `47` |

---

## 2. Ledger of changes

| # | change | class | why now |
|---|---|---|---|
| **C-1** | `age_mask` must be able to address **HFF specifically** | 🔴 **BLOCKING** | G-c step 2 is unrunnable without it |
| **C-2** | read the clock's declared `age_range` instead of discarding it | 🟠 substantive | E4/E5/E13 — the range exists, is ignored, and would change which labels are valid |
| **C-3** | stamp `donor_age` / `batch` on HFF too | 🟡 mechanical | G-b reached Gill and not HFF (E6); C-2 cannot fire for HFF without it |
| **C-4** | give the deployed response a way to say **"ΔAge is not validated here"** | ✅ **DECIDED 2026-07-31 — option (a)** | E7–E10 — the product's score is multiplied by an unvalidated term and cannot report that |
| **C-5** | make the age head **and the cell-cycle deconfounder** survive a 450×-smaller label set | 🟠 substantive | E11/E14/E15 + `build_dataset.py:449-457` — arithmetic below |
| **C-6** | record **why** a label was masked, not just that it was | 🟡 mechanical | three distinct reasons will exist after C-1/C-2; a bool cannot carry them |

**Nothing in `src/cellfate/models/`, `training/train.py`, or the fate/safety path is touched.**

---

# PART A — the `src/` changes

## C-1 🔴 BLOCKING — `age_mask` cannot address HFF

### The defect

`src/cellfate/data/aging.py:219`

```python
age_mask = np.full(age.shape[0], source not in C.CANCER_SOURCES, dtype=bool)
```

`source` is the **class attribute** `name`, and it is `"reprogramming"` for both reprogramming
sources (E2). `CANCER_SOURCES` is `frozenset({"tahoe"})` (`constants.py:73`). So today:

* every reprogramming cell gets `age_mask = True`, HFF and Gill alike;
* **there is no expressible policy that masks HFF and keeps Gill.**

Stage 1.5.2's G-c step 2 is exactly that policy. **It is not merely unrun — it is currently
impossible to run**, and no document before this one said so.

### Why the fix is small

The discriminator already exists and is already in `obs`: `dataset_id` is `"gill_bulk"`
(`sources.py:523`) and `"hff_sc"` (`sources.py:716`). `delta_age` already receives `obs`. So the
change is to consult a **policy** rather than a single frozenset.

### The change

**`src/cellfate/common/constants.py`, after line 73** — a sibling of `CANCER_SOURCES`, not a
replacement:

```python
# Datasets whose ΔAge labels are withheld from the age head by POLICY rather than by
# cell type. Stage 1.5.2 G-c: the label may be an artefact of an out-of-domain clock.
# Empty by default -- populating it is a pre-registered experiment, never a silent default.
AGE_MASKED_DATASETS: frozenset[str] = frozenset()
```

**`src/cellfate/data/aging.py`** — replace the single line 219 with a small pure helper plus one
call. The helper is separate so it is unit-testable **without constructing a clock or an expression
matrix**, which is the pattern every `diag_*` script in this project follows.

```python
def age_label_policy(
    n: int,
    source: str,
    obs: pd.DataFrame,
    *,
    masked_datasets: frozenset[str] = frozenset(),
    clock_age_range: tuple[float, float] | None = None,
) -> tuple[np.ndarray, list[str | None]]:
    """Which cells get a usable ΔAge label, and -- when they do not -- WHY. Pure.

    Three independent reasons, checked in decreasing order of certainty. The first is the
    pre-existing rule and is never weakened by the other two; a cell excluded for more than
    one reason reports the first that applied, so the string is stable under reordering of
    the later rules.

      1. ``cancer_source``   -- the clock is out of distribution on transformed lines.
                                Today's ONLY rule (``CANCER_SOURCES``), unchanged.
      2. ``dataset_policy``  -- Stage 1.5.2 gate G-c: this dataset's labels are withheld by
                                decision, not by cell type. Empty by default.
      3. ``donor_out_of_clock_range`` -- the donor's chronological age is outside the range
                                the clock was FITTED on (``configs/clocks/*.json`` ->
                                ``meta.age_range``). An UNKNOWN age never masks: absence of
                                evidence is recorded, not acted on.

    Returns ``(age_mask, reasons)`` with ``reasons[i] is None`` exactly where ``age_mask[i]``
    is True, which is the invariant ``Sample`` validation depends on.
    """
    mask = np.full(n, True, dtype=bool)
    reasons: list[str | None] = [None] * n

    def _exclude(bad: np.ndarray, why: str) -> None:
        newly = bad & mask
        mask[newly] = False
        for i in np.flatnonzero(newly):
            reasons[i] = why

    if source in C.CANCER_SOURCES:
        _exclude(np.ones(n, dtype=bool), "cancer_source")

    # FAIL LOUD, NEVER OPEN. If a withholding policy is switched on and the column it needs is
    # absent, the silent outcome is to KEEP labels that were meant to be withheld -- the unsafe
    # direction, and invisible. This is not hypothetical: C-3 records that G-b reached Gill and
    # never reached HFF, so `donor_age` is missing on HFF *today*. The step order (C-3 in step 1,
    # C-2 in step 3) protects us, but correctness must not depend on the step order.
    if masked_datasets:
        if "dataset_id" not in obs.columns:
            raise KeyError(
                "age_label_policy: masked_datasets was requested but obs has no 'dataset_id' "
                "column, so the policy cannot be applied. Refusing to silently keep labels "
                "that were meant to be withheld.")
        _exclude(obs["dataset_id"].isin(masked_datasets).to_numpy(), "dataset_policy")

    if clock_age_range is not None:
        if "donor_age" not in obs.columns:
            raise KeyError(
                "age_label_policy: clock_age_range was requested but obs has no 'donor_age' "
                "column (see C-3: G-b reached Gill, not HFF). Refusing to silently treat every "
                "cell as in-range.")
        a = pd.to_numeric(obs["donor_age"], errors="coerce").to_numpy(dtype=float)
        lo, hi = clock_age_range
        # NaN comparisons are False, so an UNKNOWN donor age cannot mask. Deliberate, and
        # distinct from the raise above: a missing COLUMN means the policy is inapplicable and
        # is an error; a missing VALUE in a present column is recorded absence, and absence of
        # evidence is not acted on.
        _exclude((a < lo) | (a > hi), "donor_out_of_clock_range")

    return mask, reasons
```

**Two tests this forces** (add to C-1's test list):

```python
def test_masked_datasets_without_the_column_raises_rather_than_keeping_labels():
    with pytest.raises(KeyError, match="dataset_id"):
        age_label_policy(3, "reprogramming", pd.DataFrame({"x": [1, 2, 3]}),
                         masked_datasets=frozenset({"hff_sc"}))


def test_clock_age_range_without_donor_age_raises():
    with pytest.raises(KeyError, match="donor_age"):
        age_label_policy(3, "reprogramming", pd.DataFrame({"x": [1, 2, 3]}),
                         clock_age_range=(1.0, 96.0))
```

and at the call site (today's line 219):

```python
age_mask, age_mask_reason = age_label_policy(
    age.shape[0], source, obs,
    masked_datasets=C.AGE_MASKED_DATASETS,
    clock_age_range=getattr(clock, "age_range", None) if enforce_clock_age_range else None,
)
return d, age_mask, age_mask_reason
```

### ⚠️ The signature change, and why it is worth it

`delta_age` currently returns a 2-tuple. Adding `age_mask_reason` makes it a 3-tuple, which touches
**two** call sites — `build_dataset.py:137` and `:145` — and the tests that call it directly
(`test_data_units.py:240`, `:250`; `test_clock_and_reprogramming.py:114`; `test_harmonize.py:104`).

**The alternative — an out-parameter like G-a's `census`** — was considered and rejected here.
G-a's census is *diagnostic*: the pipeline works identically without it. `age_mask_reason` is
**load-bearing** (C-6 persists it to the manifest, and Stage 5's claims audit reads it), so hiding
it in an optional dict would make the important thing look optional. Two call sites is a small
price for a return value that cannot be forgotten.

### The guard, and why it is credible

**With `AGE_MASKED_DATASETS` empty and `enforce_clock_age_range=False` — both defaults —
`age_mask` is bit-identical to today's**, because rules 2 and 3 do not fire and rule 1 reproduces
the existing expression exactly. Same discipline G-a shipped under, verified the same way:
`np.array_equal`, not `allclose`.

`tests/test_data_units.py:246` (`test_delta_age_masks_cancer_sources`) pins the existing behaviour.
It will need its unpacking widened to 3 values — **that is the only permitted edit**. If its
*assertions* need changing, the change is wrong and must be reverted rather than reconciled.

### The tests, written out

```python
def test_default_policy_is_bit_identical_to_the_old_expression():
    """The guard. Rules 2 and 3 are off by default, so nothing may move."""
    obs = _obs(["A"] * 6, [True] + [False] * 5, dataset_id=["hff_sc"] * 6, donor_age=[0.0] * 6)
    mask, reasons = age_label_policy(6, "reprogramming", obs)
    assert np.array_equal(mask, np.full(6, True))     # even though dataset_id and age would
    assert reasons == [None] * 6                       # qualify under rules 2 and 3

def test_a_masked_dataset_is_masked_and_its_neighbour_in_the_same_chunk_is_not():
    """The blocking capability: HFF and Gill separable inside ONE chunk."""
    obs = _obs(["HFF"] * 2 + ["O1"] * 2, [False] * 4,
               dataset_id=["hff_sc", "hff_sc", "gill_bulk", "gill_bulk"])
    mask, reasons = age_label_policy(4, "reprogramming", obs,
                                     masked_datasets=frozenset({"hff_sc"}))
    assert list(mask) == [False, False, True, True]
    assert reasons == ["dataset_policy", "dataset_policy", None, None]

def test_the_cancer_rule_wins_and_is_reported_first():
    obs = _obs(["A"] * 2, [False] * 2, dataset_id=["hff_sc"] * 2)
    mask, reasons = age_label_policy(2, "tahoe", obs, masked_datasets=frozenset({"hff_sc"}))
    assert not mask.any()
    assert reasons == ["cancer_source"] * 2      # NOT "dataset_policy" -- order is stable

def test_an_unknown_donor_age_never_masks():
    """Absence of evidence is recorded, not acted on. NaN must not silently exclude."""
    obs = _obs(["A"] * 3, [False] * 3, donor_age=[float("nan"), 53.0, 0.0])
    mask, _ = age_label_policy(3, "reprogramming", obs, clock_age_range=(1.0, 96.0))
    assert list(mask) == [True, True, False]

def test_a_missing_column_is_not_an_error():
    """A source that stamps neither dataset_id nor donor_age behaves exactly as today."""
    obs = _obs(["A"] * 3, [False] * 3)
    mask, _ = age_label_policy(3, "reprogramming", obs,
                               masked_datasets=frozenset({"hff_sc"}),
                               clock_age_range=(1.0, 96.0))
    assert mask.all()

def test_the_reason_is_none_exactly_where_the_mask_is_true():
    """The invariant Sample validation depends on (schemas.py:126-130)."""
    obs = _obs(["A"] * 4, [False] * 4, donor_age=[0.0, 53.0, 200.0, 30.0])
    mask, reasons = age_label_policy(4, "reprogramming", obs, clock_age_range=(1.0, 96.0))
    assert [r is None for r in reasons] == list(mask)
```

Plus the real-data guard, in the same shape G-a used: run all six Gill donors and one HFF chunk
through `delta_age` with defaults, and assert `np.array_equal` against the pre-change values.

> ⚠️ **Every test name in this document is lowercase, and must stay that way.** `pyproject.toml:59`
> selects the `N` rules and `:71` ignores only `N812, N818, N803, N806, N815` — **`N802` is not
> among them**, so a capital inside a test function name is a lint error, and CI runs
> `ruff check src/ tests/ scripts/` **before** `pytest`. That exact mistake kept CI red from
> 2026-07-26 until 2026-08-01 and hid the whole suite behind it. Three names in this spec were
> capitalised for emphasis and were lowercased on 2026-08-01; the emphasis belongs in the docstring.

---

## C-2 🟠 The clock's declared validity range is written and never read

### The defect

`src/cellfate/data/clock_fit.py:90` writes:

```python
"age_range": [float(ages.min()), float(ages.max())],
```

`configs/clocks/fleischer_clock.json` therefore carries `age_range: [1.0, 96.0]`. **`grep -rn
"age_range" src/` returns that single write and no read** (E4). `LinearClock.from_json`
(`aging.py:62-71`) loads `weights` and `intercept` and **drops `meta` entirely** — its last line is
`return cls(weights, intercept=float(d.get("intercept", 0.0)))`.

So the clock states the range over which it was fitted, and the pipeline discards it.

### Why it matters, with the number

HFF is a neonatal foreskin fibroblast line. GSE113957 contains **no sample below age 1**, and HFF's
day-0 baseline reads **84.5 yr** (`diag_d2_replication_results.json`, `timepoints.D0.pseudobulk_age`
= 84.459). That is the definition of extrapolation, and nothing in the code notices.

### 🔴 The hole this opens, which the naive G-c design misses

**If "outside the clock's fitted age range" is the criterion for masking HFF, it masks N2 and N3
too.** Both are declared `donor age: 0` in GEO (E13) — neonatal, below the range's lower bound of 1.

Measured on the O1 fold (E12):

| policy | age-valid TRAIN labels |
|---|---|
| today | **33 688** |
| mask HFF | **75** |
| mask HFF **and** neonatal donors | **45** |

**A stage that masks HFF for being out of range and keeps N2/N3 is internally inconsistent**, and
no previous document noticed because no previous document had `donor_age` in `src/` to check
against. G-b put it there; this is the first consequence.

### 🔵 The evidence that turns this from hygiene into a finding

**The project's own scorecard already shows that the range check would have predicted its two worst
folds.** From `scorecard/baseline.json`, split by whether the donor is inside the clock's declared
`[1, 96]`:

| metric | **N2, N3** *(age 0 — OUT of range)* | **O1, O2, Y1, Y2** *(29–53 — in range)* | ratio |
|---|---:|---:|---:|
| `dage_mae_model` | **25.74 yr** | 8.57 yr | **3.01×** |
| `dage_mae_ridge` | 23.93 yr | 9.11 yr | 2.63× |
| **`conformal_coverage`** | **0.000** | 0.601 | **0×** |
| `rank_model_dage` | 0.910 | 0.967 | 0.94× |
| `n_cells` | 21, 21 | 21, 21, 19, 21 | — |

Per donor: `dage_mae_model` = N2 **21.79**, N3 **29.70**, O1 5.39, O2 7.54, Y1 7.28, Y2 14.06.

Three things make this hard to explain away:

1. **It is not a sample-size effect** (E22) — N2 and N3 carry 21 test cells each, the same as O1,
   O2 and Y2.
2. **It is not a model artefact** (E23) — the **ridge baseline** shows the same split at 2.63×. Two
   different estimators, same two donors.
3. **Ranking barely moves** (E24) — 0.910 vs 0.967. That is the signature of a **calibration**
   failure, not a signal failure: the order is right, the level is not. Exactly what extrapolating
   past a fitted range produces.

> **`conformal_coverage = 0.000` on both out-of-range donors.** The interval that carries a
> finite-sample coverage *guarantee* covers **none** of their cells. The clock's metadata said this
> was extrapolation, in a field written at fit time, and the pipeline threw it away.

This does not prove that masking N2/N3 is right — **it proves the range field is informative**, and
that discarding it cost the project its two worst folds' worth of diagnosis.

> #### ⚠️ What the range does NOT explain — Y2
>
> *Added on review, because the aggregate hides it and a reviewer will not.* The in-range group's
> mean coverage of 0.601 is not uniform:
>
> | donor | O1 | O2 | Y1 | **Y2** |
> |---|---:|---:|---:|---:|
> | age | 53 | 53 | 29 | **35** |
> | `conformal_coverage` | 0.810 | 0.667 | 0.737 | **0.190** |
> | `dage_mae_model` | 5.39 | 7.54 | 7.28 | **14.06** |
>
> **Y2 is comfortably inside `[1, 96]` and still fails.** Three folds have broken coverage
> (N2 0.000, N3 0.000, Y2 0.190) and the range criterion identifies **two of the three**.
>
> So the claim above is correct exactly as worded — the range predicts the two worst folds *by
> MAE*, and it is informative — but it is **not a complete account of the calibration failure**,
> and it must not be presented as one. Something else is wrong with Y2. Whatever that is, it is
> not out-of-range extrapolation, and this stage does not identify it.
>
> **Consequence for C-2:** the range check is worth wiring in on the evidence above; it is *not*
> a sufficient basis for concluding that in-range donors are trustworthy. That inference would
> be refuted by Y2 on the project's own scorecard.

### The change

**`aging.py`, `LinearClock.__init__` and `from_json`** — carry the range instead of dropping it:

```python
def __init__(self, weights: dict[str, float], intercept: float = 0.0,
             age_range: tuple[float, float] | None = None) -> None:
    self.weights = weights
    self.intercept = float(intercept)
    # The range the clock was FITTED on, from its own metadata. Carried so callers can ask
    # whether a query is extrapolation; `LinearClock` itself never enforces it, because the
    # policy of what to do about extrapolation belongs to the label pipeline, not the clock.
    self.age_range = (float(age_range[0]), float(age_range[1])) if age_range else None

@classmethod
def from_json(cls, path: str | Path) -> LinearClock:
    ...
    meta = d.get("meta") or {}
    return cls(weights, intercept=float(d.get("intercept", 0.0)),
               age_range=meta.get("age_range"))
```

**`build_dataset.py`, `DataConfig`** (field list at `76-97`, E25) — one new field, defaulted off:

```python
enforce_clock_age_range: bool = False   # Stage 1.5.3 C-2: mask donors outside the clock's
                                        # fitted age_range. OFF by default because turning it
                                        # on MOVES LABELS (N2/N3 are age 0, range starts at 1).
```

Rule 3 of `age_label_policy` (C-1) then does the work; no separate code path.

### ⚠️ The decision this forces, stated rather than assumed

Enabling the range check moves labels **even with `AGE_MASKED_DATASETS` empty**, because N2 and N3
are masked by it. **So it ships off**, and turning it on is a pre-registered change with its own bar
— not a side effect of C-1.

**And there is a real argument on each side, which is why this is a decision and not a fix:**

| keep N2/N3 | mask N2/N3 |
|---|---|
| 30 of the 75 surviving labels are theirs (E12) — masking costs 40% of what is left | their `dage_mae` is 3× and their coverage is **0.000** (E20/E21) |
| the clock's `age_range` is the range of the **training donors**, not a hard validity boundary | age 0 is below the minimum by the clock's own declaration, and 84.5 yr on a neonate is the symptom |
| their **ranking** is fine (0.910) — usable for within-donor comparison, which is Stage 5's defensible claim | the *absolute* ΔAge is what Stage 2, Stage 3 and RES consume, and that is what fails |

**Recommendation: mask them for the absolute-ΔAge path and keep them for ranking** — which is not
expressible today, and is the strongest argument for C-6's `age_mask_reason`: a consumer that only
needs order can honour `donor_out_of_clock_range` differently from `cancer_source`.

---

## C-3 🟡 G-b reached Gill and never reached HFF

### The defect

`sources.py:715-716`, `GSE242423SingleCellSource.fetch`:

```python
return self.build_chunk(chunk["id"], dense, genes, self.cell_line,
                        list(pert[idx]), list(time_h[idx]), factor_as_token=True, dataset_id="hff_sc")
```

No `extra=`. Compare `GillReprogrammingSource.fetch` (`sources.py:522-524`), whose fourth argument
line is `extra={"donor_age": donor_age, "batch": batch}`.

So `donor_age` exists for six Gill donors and is **absent** for the 42 481 HFF cells — 99.7 % of the
dataset. **C-2 cannot fire for HFF without this.**

### Why it is not simply an oversight to patch silently

HFF is a **cell line**, not a donor: GSE242423's series matrix has no `donor age` characteristic to
parse, because there is no per-sample donor. The value has to be asserted from the line's identity
(human foreskin fibroblast ⇒ neonatal ⇒ 0.0) rather than read.

**That assertion belongs in the source, visibly, with its provenance** — not buried as a default.
Add a class constant so the value is greppable and has one home:

```python
class GSE242423SingleCellSource(ReprogrammingSource):
    ...
    # HFF is a neonatal human foreskin fibroblast LINE. GSE242423 declares no per-sample donor
    # age because there is no per-sample donor, so this is ASSERTED from the line's identity
    # rather than parsed -- and it is asserted HERE, visibly, because it is the value that puts
    # every HFF cell outside the clock's fitted range of [1, 96]
    # (configs/clocks/fleischer_clock.json -> meta.age_range). A silent default would hide the
    # single most consequential fact about 99.7% of the age labels.
    DONOR_AGE_YEARS: float = 0.0
    DONOR_AGE_PROVENANCE: str = "asserted from cell-line identity (neonatal foreskin); not in GEO"
```

and in `fetch` (today `sources.py:715-716`):

```python
n = len(idx)
return self.build_chunk(chunk["id"], dense, genes, self.cell_line,
                        list(pert[idx]), list(time_h[idx]), factor_as_token=True,
                        dataset_id="hff_sc",
                        extra={"donor_age": [self.DONOR_AGE_YEARS] * n, "batch": [""] * n})
```

`batch` is `""` — **not** a fabricated value — because GSE242423 has no batch structure to record.
`census_warnings` (`aging.py`) already treats a single-valued column as uninformative rather than as
a cross-batch finding, so an empty batch cannot produce a spurious D1-style warning.

### Guard

`donor_age` and `batch` are metadata and must **never** become model input. `schema.py:19`
(`model_config = ConfigDict(extra="forbid")`) forbids it at the deployment boundary and
`tests/test_baseline_census.py::test_donor_age_and_batch_are_metadata_not_model_input` pins it.
C-3 adds no new exposure — it populates an existing metadata column for a source that was skipping
it.

### Test

```python
def test_hff_stamps_its_asserted_neonatal_age_and_an_empty_batch():
    src = GSE242423SingleCellSource(samples=..., genes_file=...)
    raw = src.fetch(src.plan()[0])
    assert set(raw.obs["donor_age"]) == {0.0}
    assert set(raw.obs["batch"]) == {""}      # empty, not fabricated

def test_hff_age_is_outside_the_shipped_clocks_fitted_range():
    """The whole point of C-3: without it, C-2's rule 3 cannot fire for 99.7% of the data."""
    clock = LinearClock.from_json("configs/clocks/fleischer_clock.json")
    lo, _hi = clock.age_range
    assert GSE242423SingleCellSource.DONOR_AGE_YEARS < lo
```

---

## C-4 ✅ DECIDED — the product multiplies by a term it cannot vouch for

> ### ✅ DECISION TAKEN 2026-07-31 — **option (a), declare it.**
>
> Recorded here rather than left as three costed options. The three-way comparison below is kept
> intact, because a decision without the alternatives it was chosen over is not auditable.
>
> **Two consequences follow immediately, and the second is the useful one:**
>
> 1. **C-4 moves from step 7 to step 4** in §5's order of operations. It was scheduled last *only*
>    because the choice between (a)/(b)/(c) depended on how G-c step 2 came out. **Once the choice
>    is (a), that dependency disappears** — (a) changes no arithmetic, adds only defaulted fields,
>    and its default is the conservative answer. Nothing about it needs the retrain's result.
> 2. **It ships with the other no-decision changes**, so the deployed contract tells the truth about
>    ΔAge *before* any retrain runs rather than after.
>
> ⚠️ **What (a) does not do, stated so the record is accurate.** RES is still computed and returned
> unchanged; a caller who ignores `age_validated` gets exactly the number they got before. **(a)
> makes the documentation honest, not the output.** Closing that gap is option (c), and it remains
> **Stage 3's** question (`STAGE_3_TOOL.md:24-26`) — this decision does not pre-empt it either way.

### The defect, in three lines of shipped code

`src/cellfate/inference/res.py:38-41`

```python
R_eff = max(0.0, -(float(mu_age) + p.z_conf * float(sigma_age)))
g = R_eff / (R_eff + p.kappa)
phi = float(_sigmoid((float(S) - p.tau_safe) / p.w))
res = float(phi * (float(S) ** p.k) * g * np.exp(-p.lam * float(P_loss)))
```

`g` is a **multiplicative factor**, so `RES ∝ g(R_eff)`, and `R_eff` is a function of `mu_age` and
`sigma_age` alone. And `res.py:44-45`:

```python
elif R_eff == 0.0:
    status = REJECTED_NO_REJUVENATION
```

**So both the score and one of the three rejection statuses rest entirely on ΔAge** — the label
Stage 1.5.2 proved is produced by an instrument that is out of domain on exactly the cells this tool
is for.

`Response` (`schema.py:28-41`) returns `delta_age_mean` and `delta_age_interval` and has **no field
that can say the number is not validated** (E9). `warning` is set only for OOD (`service.py:41`), and
`OODDetector` is a latent-space Mahalanobis test that has no notion of the clock's domain (E10).

### The arithmetic that shows this is already biting

`sigma_scale = 12.67`, `z_conf = 1.0` (E19). The raw ensemble spread is ≈ 2.4 yr, so the honest
`sigma_age` served at inference is ≈ **30 yr**. `R_eff > 0` therefore requires `mu_age < −30 yr`,
against a measured effect of ≈ **−11 yr**.

> **`R_eff` is zero for essentially every query, so RES is zero and the status is
> `REJECTED_NO_REJUVENATION` — not because the treatment failed, but because the instrument cannot
> resolve it.**

This is `MASTER_PLAN.md` §5b-ter's own table ("honest, uncorrected model → `R_eff` 0.0, `g` 0.00"),
arrived at from the code. **Stage 1.5.2 supplies the missing explanation: the labels the age head
learned are artefact-laden, which is why the cross-donor residual is 34.64 yr (E18).**

### The three options, and what each costs

| option | change | cost | honest? |
|---|---|---|---|
| **(a) declare it** | add `age_validated: bool` to `Response`; set it from a bundle-level flag; leave RES arithmetic untouched | one schema field, one bundle field | ✅ says what is true, changes no number |
| **(b) refuse** | new status `REJECTED_AGE_UNVALIDATED`, `RES = 0` when the age label class is unvalidated | changes ranking behaviour | ✅ but conflates "cannot tell" with "no effect" — the exact conflation `REJECTED_NO_REJUVENATION` already makes |
| **(c) drop the term** | rank on fate/safety only when ΔAge is unvalidated; renormalise RES | changes the product | ⚠️ largest change; makes the tool honest but different |

### Option (a) in full — the recommended first move

Three files, no arithmetic touched.

**1. `src/cellfate/common/schemas.py`** — a bundle-level provenance field, defaulted so every
existing bundle stays loadable:

```python
@dataclass
class AgeProvenance:
    """What the ΔAge labels this bundle was trained on actually rest on. Stage 1.5.3 C-4.

    `validated` is False whenever the age labels come from a clock that has not been shown to
    track a ground truth on the cells it was applied to. Stage 1.5.2 measured exactly that and
    returned NOT CALIBRATABLE, so a bundle trained on RNA-clock labels must set it False and
    say so, rather than shipping a number that reads as if it were validated.
    """
    validated: bool = False
    basis: str = "rna_clock_uncalibrated"     # rna_clock_uncalibrated | methylation | mixed
    note: str = ("Stage 1.5.2: the transcriptomic clock is NOT calibratable against "
                 "methylation (rho_partial 0.267/0.516 vs a 0.50 bar). ΔAge is reported "
                 "for ranking WITHIN a donor; absolute values are not validated.")
```

**2. `src/cellfate/inference/schema.py`** — two fields on `Response` (after line 38):

```python
    age_validated: bool                # False => delta_age_* is for WITHIN-donor ranking only
    age_basis: str                     # what the label rests on; see AgeProvenance
```

**3. `src/cellfate/inference/service.py:41`** — extend the existing `warning`, which today only
fires for OOD:

```python
warnings_ = []
if not s["in_dist"]:
    warnings_.append("Out-of-distribution: prediction not trustworthy.")
if not pred.age_provenance.validated:
    warnings_.append(
        "delta_age is NOT validated in absolute terms (Stage 1.5.2: the transcriptomic clock "
        "is not calibratable). Use it to RANK conditions within this donor; do not read the "
        "number as years.")
return Response(..., age_validated=pred.age_provenance.validated,
                age_basis=pred.age_provenance.basis,
                warning=" ".join(warnings_) or None)
```

**Backward compatibility.** `Response` uses `extra="forbid"`, so adding **required** fields breaks
every existing construction site — `tests/test_inference.py:200-212` builds one field-by-field, and
`service.py:30` is the only production site. Give both fields defaults
(`age_validated: bool = False`, `age_basis: str = "unknown"`) so that:

* existing callers keep working with no edit;
* **the default is the conservative answer** — an old bundle reports "not validated", which is
  exactly what it is. A default of `True` would silently vouch for every bundle ever built.

`pred.age_provenance` loads alongside the other bundle artefacts, next to
`self.res_params = load_res_params(self.paths)` (`predictor.py:94`), with the same
missing-file-means-default behaviour.

### Why not (b)

`REJECTED_NO_REJUVENATION` already fires for the "we cannot resolve it" reason (the `sigma_scale`
arithmetic above), so adding a second status that also means "we do not know" gives the contract
**two** statuses for one condition and still no way to say which. That is worse than the status quo,
not better.

### Why (c) is the real answer, later

Option (c) — rank on fate/safety alone when ΔAge is unvalidated — is where this ends up, because it
is the only option that makes the *output* honest rather than the *documentation*. It is deferred to
Stage 3 for one reason: **Stage 3 already owns this question** (`STAGE_3_TOOL.md:24-26`, "the tool
should **decline to report ΔAge** when it cannot"), and solving it twice in two places would leave
two conventions.

**This was the one item in this stage that is a product decision rather than an engineering one.
It was costed, then decided: option (a), 2026-07-31.** The alternatives are left above so the choice
can be re-examined rather than merely trusted.

---

## C-5 🟠 The age head — and the deconfounder — starve if HFF is masked

### The arithmetic

`batch_size = 512` (E14) over **33 688** training cells. Masking HFF leaves **75** age-valid labels
(E11). Expected age labels per batch:

```
75 × 512 / 33 688 = 1.14
```

`losses.py:55-57` returns a differentiable **hard zero** when a batch contains none (E15), so with
Poisson-ish arrivals at λ = 1.14, **roughly e^(−1.14) ≈ 32 % of batches contribute nothing to the
age task at all**, and the rest contribute a Huber loss over one or two cells.

Worse, `MultiTaskLoss` (`losses.py:61+`) learns `s_age = log σ²_age` by Kendall & Gal weighting from
those same one-or-two-cell batches. **A learned task weight estimated from ~1 observation per step
is not a weighting; it is noise.**

### ✅ The part that turns out to be fine — checked rather than assumed

I expected the conformal calibration to collapse too. **It does not, and the reason is worth
recording**: `train_model.py:254` fits `q` on `xstats.abs_residuals`, which is **cross-donor**
inner-LODO, and `xdonor_calib.py:52` (`MIN_INNER_TRAIN_FRAC = 0.5`) **skips the HFF fold entirely**
because holding HFF out leaves 0.2 % of the training data.

`abs_residuals.shape = (103,)`, and 124 non-HFF cells − 21 held-out O1 = **exactly 103** (E16, E17).

> **The ΔAge conformal quantile and `sigma_scale` are already fitted on Gill cells only, with zero
> HFF contribution. Masking HFF cannot degrade them, because they never used HFF.**

That is a genuine argument *for* masking being less disruptive than it looks, and it was found by
checking a worry rather than by assuming one. It also sharpens the diagnosis: **the model learns
ΔAge from 33 613 HFF labels and is then graded against 103 Gill residuals** — which is precisely why
`q` is 34.64 yr.

### The change

No single fix is obviously right, so this ships as a **pre-registered choice with a bar**, not a
default:

The loader today is `train_dl = loader(train_ds, cfg.batch_size, shuffle=True)` (`train.py:117`,
E26) — plain shuffling, no sampler. Three options, in increasing order of intrusiveness:

**Option 1 — age-aware sampling (recommended).** Keep the batch size; change *which* cells land in
each batch so the age labels are spread rather than clumped:

```python
# Stage 1.5.3 C-5. With 75 age labels among 33,688 cells, uniform shuffling puts ~1.14 of
# them in a 512-cell batch and NONE in ~32% of batches -- and losses.py:55-57 makes such a
# batch contribute exactly zero to the age task. Weighting the sampler so age-valid cells are
# drawn at a rate that guarantees several per batch fixes the occupancy without touching the
# loss, the batch size, or the fate task's effective sample size.
w = np.where(age_mask_train, age_weight, 1.0)
sampler = WeightedRandomSampler(w, num_samples=len(train_ds), replacement=True)
train_dl = DataLoader(train_ds, batch_size=cfg.batch_size, sampler=sampler)
```

⚠️ **`replacement=True` changes the fate task's sampling too**, so this is not free: the fate head
would see the 75 age-valid cells oversampled. That must be measured, not waved through — the fate
guards (`fate_prauc`, `fate_roc`, `fate_ece`) exist precisely to catch it.

**Option 2 — gradient accumulation for the age term only.** Leave sampling alone; accumulate the
age loss across a window of *W* batches before stepping it. Preserves the fate task exactly, at the
cost of a more complicated training loop and an effective age-batch of `1.14 × W`.

**Option 3 — pin the task weight.** `MultiTaskLoss` (`losses.py:61`) learns `s_age = log σ²_age`.
With ~1 observation per step that estimate is noise. Replace it with a constant when the age-label
count falls below a threshold. Smallest change; does nothing about the 32% empty batches.

### The bar, pre-registered

**Whichever option is chosen, the metric it is judged on goes through `audit_metrics.bar_verdict`
BEFORE the run** (`REF_GROUND_RULES.md` §5b), and gets a row in `tests/test_bars_resolvable.py`.
Given **four hairline margins already on this project's record** (E1b 0.009, D2 0.014, M-2a 0.016,
G-c 0.084), this is not the place to skip it.

**The guard that must not move:** the fate metrics. If Option 1 is chosen, `fate_prauc`, `fate_roc`
and `fate_ece` are re-read before/after, and a move in any of them is a finding to explain — because
oversampling 75 cells is exactly the kind of change that helps one head by quietly reshaping the
other's training distribution.

> ### ✅ STEP 5 EXECUTED 2026-08-02 — the bar is registered, and **it overturned the recommendation above**
>
> `python plan_tests/register_c5_bar.py` → `results/register_c5_bar_results.json`. No `src/` file
> touched, no label moved, no retrain. **699 tests pass**, ruff clean.
>
> *(Correction, 2026-08-02: the committed state of step 5 is **703** passing, not 699 — the figure
> was written mid-step and four more tests went into the same commit before it landed. Left as
> written above per the annotate-never-rewrite rule; measured with
> `pytest --ignore=tests/test_c5_deeper_tests.py` at `724e359`. Nothing else in the block changes.)*
>
> #### What the bar had to be, given the gate
>
> Step 5's gate is *"bar RESOLVABLE **before any retrain**"*, and that decides what it can measure:
> `dage_mae_model` needs step 6's run, so the bar must grade the **mechanism**, not the outcome.
> Two bars, because "non-zero" is not "usable":
>
> | | | |
> |---|---|---|
> | **B1** | P(update contributes **any** age gradient) | ≥ 0.95 |
> | **B2** | P(that gradient uses **≥ 4 cells**) | ≥ 0.95 |
>
> **B1 alone would have been too easy.** C-5's diagnosis is not only the 32 % empty batches — it is
> also that the survivors carry *"a Huber loss over one or two cells"*. A mechanism can clear B1 and
> still feed the optimiser per-cell noise. `k = 4` is the smallest value that halves the per-update
> standard error relative to a single cell (SE ∝ 1/√m), and any `k ≥ 2` is already a bar only an
> oversampling or accumulating mechanism can meet.
>
> #### The result
>
> | candidate | mean cells/update | B1 | B2 | |
> |---|---:|---:|---:|---|
> | status quo (uniform shuffling) | 1.15 | 68.9 % | 2.9 % | ❌ **FAIL** |
> | **Option 3** — pin `s_age` only | 1.14 | 68.4 % | 2.8 % | ❌ **FAIL** |
> | **Option 2** — accumulate, W = 8 | **9.13** | **100 %** | **98.2 %** | ✅ **PASS** |
> | **Option 1** — sampler, w = 7.1 | 7.97 | 100 % | 96.2 % | ✅ **PASS** |
>
> **Resolvable:** the dense regime — every cell age-labelled, i.e. today, before masking — clears
> both at **100 %**. **Discriminating:** the bar separates the candidates, and the script exits
> non-zero if it ever stops doing so. A bar everything passes decides nothing.
>
> #### 🔴 The recommendation above is WITHDRAWN. Option 2, not Option 1.
>
> C-5 recommended Option 1 on intuition. Measured, **Option 2 wins on both axes that matter**:
>
> 1. **It scores higher on the harder bar** — 98.2 % vs 96.2 % on B2.
> 2. **It costs the fate task nothing.** Option 1 needs `w = 7.1`, which oversamples the 75 age
>    cells **7.0×** — from 0.223 % to 1.563 % of every batch, a **1.34 %** shift in the fate head's
>    training mix. C-5 flagged that as *"not free"*; this is the number, and Option 2's is zero
>    because it changes no sampling at all.
>
> Option 1's only advantage was simplicity, and it buys that by putting Stage 1's guard record — the
> `+0.000` bit-identical run — at risk for no measured gain.
>
> **Option 3 is dead, and now provably so:** it is `weight=1, accumulate=1`, i.e. *identical to the
> status quo by construction*. Pinning `s_age` does nothing about occupancy. That was C-5's
> criticism of it; it is now a measurement.
>
> #### What still rides on step 6
>
> The bar grades the mechanism. Whether the age head actually **learns** from 75 labels is
> `dage_mae_model` at step 6, and no simulation can answer it. **The fate guards
> (`fate_prauc`, `fate_roc`, `fate_ece`) must still read "noise"** — with Option 2 there is no
> resampling to disturb them, which is exactly why it is the safer choice.
>
> Registered in `tests/test_bars_resolvable.py` (6 rows) with 12 unit tests on the pure functions,
> including the closed-form checks: the uniform mean reproduces C-5's 1.14, and the empty-batch rate
> matches both the exact binomial `(1−p)^512` and the plan's `e^−1.14` estimate.

> ### ✅ STEP 5b EXECUTED 2026-08-02 — deeper tests before committing to the option
>
> `python plan_tests/c5_deeper_tests.py` → `results/c5_deeper_tests_results.json`. READ-ONLY:
> no `src/` file, no label moved, no training. 30 unit tests
> (`tests/test_c5_deeper_tests.py` + `tests/test_register_c5_bar.py`), ruff clean.
>
> #### Why a second script, when step 5 already chose
>
> B1/B2 grade **occupancy** — does an update get an age gradient, over how many cells. Deciding a
> design on that alone is deciding on the one axis that happened to get measured. Seven axes it
> cannot see were tested; a fourth candidate (a hybrid) was added so the comparison was not forced
> between two extremes. **Two of the seven changed the reading, and one of my own claims was
> weaker than I had stated it.**
>
> | candidate | eff cells | dup | **grad upd** | cover | donor | **reps/ep** | fate churn |
> |---|---:|---:|---:|---:|---:|---:|---:|
> | status quo (shuffle) | 1.14 | 1.00 | 2 660 | 98.8 % | 1.01 | 0.99 | 0.0 % |
> | Option 1 — sampler w = 7.1 | 7.59 | 1.05 | **3 900** | 99.9 % | 1.04 | **6.93** | **36.4 %** |
> | **Option 2 — accumulate W = 8** | **9.07** | **1.00** | 480 | 98.5 % | **1.01** | **0.98** | **0.0 %** |
> | Option 4 — hybrid w = 3, W = 3 | 9.41 | 1.07 | 1 260 | **94.8 %** | **1.06** | 2.92 | 36.1 % |
>
> #### D6 — the diagnostic that actually settled it: *information* vs *repetition*
>
> A sampler weight does not create labels. **There are 75 and there will be 75.** Weight `w` runs
> `w` age-epochs inside every fate-epoch, so across the run:
>
> | | passes over the same 75 labels |
> |---|---:|
> | status quo / **Option 2** | **59** — one pass per epoch, i.e. what "60 epochs" means |
> | Option 4 | 175 |
> | Option 1 | **416** |
>
> Option 1's extra gradient updates are bought entirely by **re-showing the same 75 labels 7× per
> epoch**. That is not more data; it is 416 effective epochs over 75 examples — a memorisation
> regime. Worse for step 6 specifically: it changes *three* things at once (delivery, exposure, and
> the fate head's sampling), so a `dage_mae_model` move could not be attributed to the fix. **Step 6
> is a diagnostic retrain whose entire purpose is attribution**, and the one-change rule applies.
>
> Option 2 changes **delivery only** — same labels, same one pass per epoch, same fate training set,
> regrouped so no update is empty. That is the minimal intervention, and it is asserted as a test
> (`test_accumulation_changes_delivery_and_not_exposure`): if it ever stops being true, this
> decision must be reopened.
>
> #### D7 — I would have overstated Option 2's cost by 47 %
>
> The status quo does **not** get 3 900 age updates. 32 % of its batches hit the hard zero at
> `losses.py:55-57`, so it gets **2 660**. Comparing Option 2's 480 against 3 900 — which is what
> the update count alone invites — inflates the apparent cost of accumulation by nearly half.
>
> #### 🟡 D5 came out WEAKER than step 5 implied, and is corrected here
>
> Step 5 quantified Option 1's fate cost as a 1.34 % shift in batch composition and I expected the
> bootstrap to be the larger, unmeasured cost. Per epoch it is — **36.4 %** of fate cells are missed.
> But a bootstrap **re-rolls its misses every epoch**: over 60 epochs `P(a cell is never seen)` is
> **1.8 × 10⁻²⁶**. Nothing is deleted. The real cost is variance — **59.3 ± 7.7 visits, CV 13 %** —
> against a permutation's exact 60. That is a genuine cost, but it is *sampling variance*, not
> *data loss*, and reporting the 36 % alone would have been the same kind of overclaim this stage
> keeps catching elsewhere.
>
> #### Option 4 (the hybrid) loses on its own merits, not by exclusion
>
> It was added to break a forced choice, and it is **dominated**: it still pays Option 1's full
> bootstrap cost (36.1 %), still repeats labels 3×, and posts the **worst label coverage (94.8 % —
> it misses 4 of the 75 labels in an average epoch)** and the **worst donor balance (1.06)** of any
> candidate, for 1 260 updates. There is no axis on which it is the best choice.
>
> #### 🔵 W = 8, and W = 7 is rejected for a measured reason
>
> W = 8 was chosen for comfort, not derived, so the whole range was swept. **W = 7 is the smallest
> that clears B2** (95.6 % vs the 95 % bar) and buys 540 updates instead of 480. It is still the
> wrong choice — 75 is not a constant, it is what survives C-1 masking *on this fold*, and other
> LODO folds hold out other donors:
>
> | n_age | W = 7 | W = 8 |
> |---:|---:|---:|
> | 75 | 95.6 % ✅ | 98.1 % ✅ |
> | 70 | **93.7 % ❌** | 97.2 % ✅ |
> | 65 | 91.3 % ❌ | 95.6 % ✅ |
> | 60 | 88.0 % ❌ | 93.4 % ❌ |
>
> W = 7 falls below its own bar as soon as the label count moves at all. **12.5 % more updates is not
> worth sitting 0.6 pp above the bar.** W = 8 holds to n_age ≥ 65; below that, C-5 needs revisiting
> regardless of W — recorded here so it is a known boundary rather than a surprise at step 6.
>
> #### ✅ Decision: **Option 2, W = 8** — confirmed, now on seven axes rather than one
>
> Step 5's choice survives the deeper tests and is strengthened: Option 2 is best or tied-best on
> every axis except update count, and D7 shows that gap is 5.5× rather than the 8× it appears.
>
> #### 🟠 The one residual risk, pre-registered
>
> **480 age updates may be too few to converge**, and no simulation can tell — it is `dage_mae_model`
> at step 6. Pre-registered contingency so it is not decided after seeing the answer: if the age
> head is still underfit at the end of the run (training age loss visibly falling at the final
> epoch), the remedy is a higher **age learning rate**, *not* a smaller W — W is pinned by the
> sensitivity table above, and lowering it trades a bar-margin for updates. Any such change is its
> own pre-registered Change with its own bar.
>
> #### ⚠️ The implementation trap step 6 must avoid — pinned now, before any code is written
>
> `huber_age_loss` (`src/cellfate/models/losses.py:48-58`) ends in
> `F.huber_loss(age_pred[m], age_true[m], delta=delta)` — **`reduction='mean'` by default**, over
> the valid cells *in that batch*. So the obvious implementation of Option 2 — average the per-batch
> age losses over the window — is **wrong**: it weights a 1-cell batch exactly as heavily as a
> 9-cell batch, which is the very defect C-5 exists to remove, moved up one level.
>
> **Requirement, not a suggestion:** the window's age loss must be `Σ(per-cell losses) / Σ(valid
> cells)` across all W batches — one mean over the window's cells, not a mean of means. That needs
> `reduction='sum'` plus a running valid-cell count.
>
> Two second-order points that follow, left to step 6 to resolve but recorded so they are not
> discovered late:
> * the normaliser is not known until the window closes, so either hold the graph across W batches
>   (memory ×W) or normalise by the expected count `W × 1.14 ≈ 9.1` and accept a small bias — the
>   second is standard gradient accumulation and almost certainly right here;
> * `MultiTaskLoss` (`losses.py:61`) sums the fate and age terms, and **the fate term must keep
>   stepping every batch**. Only the age term accumulates. If both accumulate, Option 2 has silently
>   become "train 8× less", the fate guards will move, and its whole claim to cost the fate task
>   nothing is void.
>
> A test asserting the mean-of-means error is *not* present belongs with that change.

> ### 🛑 READINESS AUDIT 2026-08-02 — step 6 is **NOT** ready. Two problems, both found by asking.
>
> Asked "are we ready for step 6", I checked rather than answered, and the answer is no.
>
> #### Problem 1 — no step actually implements C-5
>
> The step table runs 1, 2, 3, 4, 5, 5b, 6. Step 5 is *"C-5 **design** + its bar"* and 5b chose the
> option. **PART D's manifest lists `training/train.py` as a file this stage changes, but no step
> schedules that change.** `train.py:117` is still `train_dl = loader(train_ds, cfg.batch_size,
> shuffle=True)` — plain shuffling, exactly as E26 recorded it. C-5 is graded and unbuilt.
>
> This is not a bookkeeping gap. **Step 6's arm B *is* the starved regime C-5 exists to fix** — 75
> labels, 1.14 per batch, 32 % of updates a hard zero. Running step 6 as it stands measures
> "do HFF's labels help?" **confounded with** "is the age head trainable at 75 labels with the
> current loader?" — and the pre-registered reading *"A better ⇒ HFF's labels help, keep them"*
> would then be wrong for a reason the outcome table cannot express.
>
> **A new step 5c is required: implement C-5 Option 2, in both arms, before step 6 runs.**
>
> #### Problem 2 — 🔴 a fixed W = 8 biases step 6 **toward its own treatment**
>
> This one I got wrong in 5b and it is the more serious of the two. I pinned W = 8 by asking what
> the **masked** regime needs. Step 6 runs **two** arms, and arm A is not masked:
>
> | | age-valid cells | age cells per batch | age updates/epoch at fixed W = 8 | vs today |
> |---|---:|---:|---:|---|
> | **arm A** (control) | **33 688 of 33 688** | ~512 | 8 | **65 → 8, an 8× cut for no reason** |
> | **arm B** (treatment) | 75 of 33 688 | 1.14 | 8 | 44 → 8, but each is usable |
>
> Arm A has **no occupancy problem** — every batch is full. Fixed W = 8 buys it nothing and costs it
> 8× its age optimisation. **So the mechanism handicaps the control and helps the treatment**, which
> pushes `dage_mae_model` toward *"B better, CI excludes 0"* — one of the three pre-registered
> outcomes, and the one that would conclude *"99.7 % of the labels were net-negative."* A mechanism
> that tilts the result toward the treatment conclusion is a validity threat, not a detail.
>
> #### The fix — one rule, not one constant
>
> Make the trigger the **accumulated age-cell count**, not a batch count: *step the age term once the
> window holds ≥ k age cells, or after W_max batches, whichever comes first.*
>
> | | behaviour | effect |
> |---|---|---|
> | **arm A** | first batch already holds ~512 ≥ k ⇒ W = 1 | **identical to today** — the control is left alone, and `scorecard/baseline.json` stays meaningful |
> | **arm B** | ~7–8 batches to reach k ⇒ W ≈ 8 | exactly the regime 5b validated |
>
> It is **one policy applied identically to both arms**; it only *behaves* differently because the
> data differ, which is the definition of a controlled comparison rather than a confound. It also
> satisfies B2 **by construction** rather than at 98.1 % probability, and `W_max = 8` from 5b's
> sensitivity table becomes the cap for when even 8 batches cannot reach `k`.
>
> **5b's W = 8 analysis is not withdrawn** — it is still what fixes `W_max`, and the n_age ≥ 65
> boundary still holds. What changes is that W = 8 is a **ceiling**, not a constant.
>
> #### Consequences to carry into 5c
>
> * `k` needs registering through `bar_verdict` like every other bar. `k = 4` is B2's existing
>   threshold and the natural default; a larger `k` trades updates for per-update SNR.
> * The **data-dependent stop is deterministic** given the shuffle seed, so reproducibility holds —
>   but it must be asserted, not assumed.
> * Every label is still used **exactly once per epoch**, so the stopping rule selects *windows*,
>   not *labels*, and introduces no selection bias. This also needs a test.

> ### ✅ STEP 5c EXECUTED 2026-08-02 — C-5 Option 2 is built, and it ships **inert**
>
> `python plan_tests/register_c5c_bar.py` → `results/register_c5c_bar_results.json`, then the code.
> **743 tests pass** (+22), ruff clean. **No label moved, no retrain, nothing rebuilt.**
>
> #### The bar went first, and it failed — which is the point of putting it first
>
> §5b of the ground rules: bar before the change. Attempt 1 forced the window to close at each
> epoch's last batch, so every label would be consumed inside its own epoch. It scored **93.9 %**
> against A2's 95 % and **failed**.
>
> I did not lower the bar. Attributing the shortfall: the epoch-end window was **4.44 pp** of the
> 6.12 pp gap and the irreducible `W_max` limit only **1.67 pp** — so the *mechanism* was wrong, not
> the bar. Letting the window **carry across the epoch boundary** removes the artificial partial
> window entirely, and is *simpler code* (one fewer special case). Re-run: **98.2 %, PASS.**
>
> | bar | what it grades | result |
> |---|---|---|
> | **A1** | control arm closes every window at W = 1 — an **equality**, not a rate | **1.0000** ✅ |
> | **A2** | P(window holds ≥ 4 age cells) in the masked arm | **98.2 %** ✅ (bar 0.95) |
> | **A3** | masked arm gets *more* age updates than fixed W = 8 would | **980 vs 480** ✅ |
>
> #### A3 is a bonus I did not expect: the bias fix also **doubles** the age optimisation
>
> Triggering on cells rather than batches closes a window as soon as it is worth stepping on, so the
> masked arm gets **16.3 updates/epoch (980 over the run)** instead of fixed-W's 8/epoch (480) — at
> the *same* per-update quality. **That directly reduces the "480 updates may be too few to
> converge" risk that 5b had to leave open**, without touching the learning rate.
>
> #### What shipped
>
> | file | change |
> |---|---|
> | `models/losses.py` | `+ huber_age_window()` — one Huber over the window's cells, `Σloss/Σcells` |
> | `models/__init__.py` | export it |
> | `training/train.py` | `+ _AgeWindow`, and 6 lines in the batch loop |
> | `training/train_model.py` | `+ age_window_k: int = 1`, `+ age_window_max_batches: int = 8` |
>
> **`age_window_k = 1` is the default, and 1 means OFF — the pre-1.5.3 path, bit for bit.** It ships
> inert on purpose: this stage's guard is that nothing moves until step 6 turns it on deliberately,
> in **both** arms. That also makes the rollback a one-value edit rather than a revert.
>
> #### The gate, proved rather than asserted
>
> `test_arm_a_is_bit_identical_when_every_cell_is_age_valid` runs `train_member` twice — mechanism
> off, then on — and compares **every parameter tensor** with `torch.equal`. It passes, and holds
> for k ∈ {2, 4, 8, 16}.
>
> A test that only asserts invariance can pass on a no-op, so two more sit beside it: one confirming
> the mechanism **does** move a sparsely-labelled run, and — the real check — **I re-injected the
> exact bug the readiness audit found** (a fixed-W window that ignores the cell count) and confirmed
> it fails **both** arm-A identity tests plus the drift check, then restored. The guard catches the
> thing it was built to catch.
>
> 18 tests in `tests/test_c5c_age_accumulation.py` + 5 rows in `tests/test_bars_resolvable.py`,
> covering all five gates: arm-A identity, `Σloss/Σcells` (constructed so mean-of-means gives a
> visibly different answer), the fate head still stepping on a held-back batch, determinism under a
> fixed seed, and windows-not-labels. One test drives the **shipped** `_AgeWindow` against the bar
> script's `close_windows` over 30 random sequences, so the simulation the decision rests on cannot
> drift from the code that ships.
>
> #### ⚠️ Still open, and unchanged by this step
>
> Whether the age head **learns** from 75 labels is `dage_mae_model` at step 6. 5c improves the
> odds (980 updates, not 480) and removes a bias; it settles nothing about the outcome.

### 🔴 The second consequence, which is easy to miss: the cell-cycle deconfounder moves too

`build_dataset.py:449-457` fits the deconfounder on **age-valid TRAIN cells only**:

```python
for i, cell in enumerate(aux.cell_ids):
    if aux.age_mask[i] and cell in train_ids:
        d_tr.append(aux.d_age_raw[i]); cc_tr.append(aux.cc[i])
coef = fit_deconfounder(...) if len(d_tr) >= 2 else (0.0, 0.0)
```

So masking HFF moves the fitting set for `coef` from **33 613 single-cell HFF** cells to **75 bulk
Gill** samples. Two things change at once, and neither is a sample-size problem:

* **the modality changes** — `cell_cycle_score` on a **bulk** RNA-seq sample is a population average,
  not a per-cell phase score. It is a different quantity with a different variance;
* **the guard is `len(d_tr) >= 2`** — there is no minimum that would notice a 450× drop, and no
  warning is emitted.

`coef` is then applied to **every** shard (`build_dataset.py:462`), including the masked HFF cells,
whose `y_age` is subsequently NaN'd (`463-466`) — so the masked cells are unaffected in the end, but
**every surviving Gill label is deconfounded by a coefficient estimated from bulk data.**

**Add to step 6's pre-registration:** report `coef` before and after, and treat a large move as a
finding to explain rather than a number to absorb. If the deconfounder is judged inapplicable at
bulk resolution, `cfg.deconfound = False` is the honest alternative and must be chosen explicitly.

---

## C-6 🟡 `age_mask` records *that* a label was withheld, never *why*

### The defect

`schemas.py:55` — `age_mask: bool`. After C-1 and C-2 there will be **three** distinct reasons a
label is withheld:

1. cancer / transformed line (today's only reason);
2. dataset masked by G-c policy;
3. donor outside the clock's fitted age range.

A boolean cannot distinguish them, and the downstream consumers need to: `STAGE_5_PUBLICATION.md`'s
claims audit must state *which* cells are excluded and why, and `STAGE_6_NEW_DATA.md:31` already
instructs acquirers to mask labels — with no way to record the reason they did.

### Why a bool is not enough, concretely

Two consumers already need the distinction and cannot get it:

* **`STAGE_5_PUBLICATION.md`'s claims audit** must state which cells are excluded and why. "33,613
  cells were excluded" is not a claim a reviewer can assess; "33,613 excluded by G-c policy, 30 by
  the clock's fitted range, 0 as cancer" is.
* **C-2's recommendation** — mask the out-of-range donors for *absolute* ΔAge but keep them for
  *ranking* — is **not expressible** with a boolean. A consumer that only needs order can honour
  `donor_out_of_clock_range` differently from `cancer_source`; it cannot honour `False` differently
  from `False`.

### The change, all four sites

```python
# schemas.py:55 -- Sample
age_mask: bool                          # True iff y_age is a valid label
age_mask_reason: str | None = None      # why it is NOT, when age_mask is False; see
                                        # aging.age_label_policy. None iff age_mask is True.

# schemas.py -- the existing validator at 124-130 gains one line, so the invariant is
# enforced rather than documented:
if self.age_mask and self.age_mask_reason is not None:
    raise ValueError("age_mask=True requires age_mask_reason to be None")

# io.py:134 -- _SHARD_SCHEMA
("age_mask", pa.bool_()),
("age_mask_reason", pa.string()),       # nullable

# io.py:236 -- _MANIFEST_SCHEMA
("age_mask", pa.bool_()), ("age_mask_reason", pa.string()),
```

`shard_to_numpy` (`io.py:194-228`) gains one key; `assemble_samples` (`assemble.py:31-50`) passes it
straight through alongside `age_mask`, which it already computes per row.

### ⚠️ The compatibility consequence, stated up front

`io.py:194-228` reads a **fixed key list** from `table.to_pydict()`, so a shard written before this
change has no `age_mask_reason` column and **the new reader raises `KeyError` on it**. There is no
schema-evolution layer in this repo — I checked; `_SHARD_SCHEMA` is a literal `pa.schema([...])` with
no version field.

**Two honest ways to handle it, and the choice is not free:**

| | approach | consequence |
|---|---|---|
| **(i)** | tolerate absence: `d.get("age_mask_reason", [None] * n)` | old shards keep working; the cost is that "column missing" and "all reasons None" become indistinguishable |
| **(ii)** | require it, and rebuild | clean invariant; **every existing shard must be regenerated** |

**Take (ii).** C-1 and C-2 change the labels, so a rebuild is forced regardless — and (i) would buy
compatibility with artefacts that are about to be regenerated anyway, at the price of a permanent
ambiguity in the one field whose whole purpose is to be unambiguous.

> **This reasoning is valid only because a rebuild is already happening. `age_mask_reason` must not
> ship in a release that does not rebuild** — in that situation (i) would be the right call, and
> shipping (ii) would break every bundle on disk for no benefit.

---

# PART B — the downstream plan changes

Stage 1.5.2's verdict has **not propagated**: `grep -n "1\.5\.2\|NOT CALIBRATABLE\|calibratab"` over
`MASTER_PLAN.md`, `STAGE_2`, `STAGE_4`, `STAGE_5`, `STAGE_6` and `REF_GROUND_RULES.md` returns
**zero hits** in all six. Each annotation below is **purely additive** — per the standing rule, no
existing line is edited, and the correction sits beside the claim it corrects.

| # | file | line | what is there now | why it needs a note |
|---|---|---|---|---|
| **P-1** | `STAGE_3_TOOL.md` | **450** | `rec = min(ok, key=lambda o: o.delta_age)` | **the tool's entire recommendation is an argmin over ΔAge.** If ΔAge is unvalidated, so is every recommendation the tool makes. This is the single most consequential line in the downstream plans |
| **P-2** | `STAGE_3_TOOL.md` | **385-391** | `StoppingOption` carries `delta_age`, `delta_age_lo/hi` and no validity field | mirrors C-4 at the plan level: the option object cannot express "not validated" either |
| **P-3** | `MASTER_PLAN.md` | **164-176** | §5b-ter's table showing `R_eff = 0.0` for the honest model, and the 2026-07-26 note that the −11 yr effect is artefact-suppressed | §5b-ter diagnosed the *symptom* (honest uncertainty kills `R_eff`). **1.5.2 supplies the cause** — the labels come from an out-of-domain instrument, which is why the cross-donor residual is 34.64 yr |
| **P-4** | `STAGE_2_LEVEL_CORRECTION.md` | **135-145** | §0's non-code prerequisite: k ≈ 3 reference cells per donor | Stage 2 corrects a **per-donor offset in ΔAge**. 1.5.2 says the quantity being offset is not validated. **Correcting the offset of an unvalidated quantity is not a fix**, and this must be read before any wet-lab spend |
| **P-5** | `STAGE_5_PUBLICATION.md` | **160-180** | the claims audit, and §2's "what cannot be claimed" list | two rows need company: the ΔAge-derived claims now have a named upstream limit, and **"the RNA clock is calibratable against methylation" is now a *disproved* claim**, which belongs in §2's ❌ list |
| **P-6** | `STAGE_6_NEW_DATA.md` | **31** and **33-35** | already says unlabelled data must be masked, and that methylation is worth more than more RNA | **1.5.2 promotes both from advice to requirement**, and REV FINAL §11.2 adds the sizing rule: ≈16 transient-arm pairs, sized against a **±7 yr between-donor** instrument error |
| **P-7** | `REF_GROUND_RULES.md` | **91-118** (§5b) | the resolvability rule | 1.5.2 produced two lessons §5b does not yet carry: **(i)** a pre-committed *fallback* must be resolvability-checked too — 1.5.2's was not, and was itself unresolvable at 92.3 %; **(ii)** a bar near the empirical ceiling of the instrument is not a bar, and the ceiling must be measured (meth↔meth ρ_partial +0.568 against bars assuming 0.70) |

**Additional note for P-1/P-2:** `STAGE_3_TOOL.md:8-28` already contains an internal-control design
question ("the tool should **decline to report ΔAge** when it cannot"). **1.5.2 makes that question
load-bearing rather than a refinement**, and C-4 option (c) is its natural landing place. The
annotation should connect them rather than open a third thread.

---

## PART B.2 — the annotation text, ready to paste

Each block is **appended at the stated anchor**, never inserted over an existing line. They are
written out here so PART B is executable rather than a to-do list, and so the wording is reviewed
once rather than improvised six times.

### P-1 + P-2 → `STAGE_3_TOOL.md`, after line 451 (end of the recommendation rule)

```markdown
> ## 🆕 ADDED 2026-07-31 — the recommendation rule sorts on a label Stage 1.5.2 could not validate
>
> *Additive; §4's rule and the `StoppingOption` dataclass above are unmodified.*
>
> `rec = min(ok, key=lambda o: o.delta_age)` makes **ΔAge the sole ordering key** of this tool's
> output. `STAGE_1_5_2_LABEL_ANCHOR.md` returned **NOT CALIBRATABLE**: the transcriptomic clock does
> not track methylation age (ρ_partial +0.267 / +0.516 against a pre-frozen 0.50 bar), and §12-R
> confirmed the anchor itself is sound, so the failure is the RNA clock's.
>
> **What this does and does not invalidate:**
>
> | | |
> |---|---|
> | ❌ **absolute** ΔAge as a quantity to compare across donors | the number is not validated in years |
> | ✅ **within-donor ordering** of withdrawal times | `rank_model_dage` is **0.91–0.99** across all six folds (`scorecard/baseline.json`), and 1.5.2 §17 found the RNA clock reaches 91% of the meth↔meth ceiling against one of the two references |
>
> **So the rule survives as a *ranking* rule and fails as a *reporting* rule** — which is precisely
> the distinction §0.3's internal-control note above already asks this stage to make. The two should
> be resolved together, not separately.
>
> **Concretely, for 3d:** `StoppingOption` (line 385) needs a validity field beside `delta_age`, and
> `_to_day`'s output should be presented as an ordered shortlist rather than a ΔAge in years, unless
> and until the label basis changes. `STAGE_1_5_3_EXECUTE.md` C-4 costs the equivalent change at the
> inference boundary; do not solve it twice with two conventions.
```

### P-3 → `MASTER_PLAN.md`, after line 176 (the existing 2026-07-26 note in §5b-ter)

```markdown
> **🆕 ADDED 2026-07-31 — §5b-ter diagnosed the symptom; Stage 1.5.2 found the cause.**
> *Additive; the table and the 2026-07-26 note above are unmodified.*
>
> §5b-ter's row *"honest, uncorrected model → R_eff 0.0, g 0.00"* is correct and its arithmetic
> holds. What it could not say is **why** the honest uncertainty is ~39 yr in the first place.
>
> `STAGE_1_5_2_LABEL_ANCHOR.md`: the ΔAge labels are produced by a clock that is **out of domain on
> reprogramming cells and NOT calibratable** against methylation. The deployed cross-donor conformal
> half-width is **q = 34.64 yr** (`bundle/conformal.json`) — the same order as §5b-ter's ~39 yr, and
> now with a mechanism rather than a measurement.
>
> **The consequence for this section:** the fix §5b-ter implies (better uncertainty, level
> correction) cannot recover `R_eff`, because the problem is not that the uncertainty is badly
> estimated — it is that **the target it is estimating uncertainty about is artefact-laden**. That
> makes §5c's per-donor level correction (Stage 2) a correction to an unvalidated quantity; see the
> note added to `STAGE_2_LEVEL_CORRECTION.md` §0.
```

### P-4 → `STAGE_2_LEVEL_CORRECTION.md`, after line 145 (end of §0's prerequisite)

```markdown
> ## 🛑 ADDED 2026-07-31 — a SECOND non-code prerequisite, and it is upstream of the first
>
> *Additive; §0's wet-lab prerequisite above is unmodified and still applies.*
>
> §0 asks whether k ≈ 3 reference cells per donor is experimentally acceptable. **Ask this first:**
>
> > **Stage 2 corrects a per-donor OFFSET in ΔAge. `STAGE_1_5_2_LABEL_ANCHOR.md` established that
> > ΔAge itself is not validated — the clock producing it is NOT CALIBRATABLE. Correcting the offset
> > of an unvalidated quantity does not make it validated.**
>
> **This does not cancel Stage 2.** The offset may well be real, and §1's evidence for it stands.
> What changes is the order of operations and what a success would mean:
>
> | | |
> |---|---|
> | ✅ still worth doing | the offset is measurable and the k≈3 recovery result (T7.4.3, T16) is unaffected |
> | ⚠️ changed | a corrected ΔAge is a **better-centred unvalidated number**, not a validated one |
> | 🛑 **do not spend wet-lab budget** | until `STAGE_1_5_3_EXECUTE.md`'s G-c step 2 has run. If HFF's labels are masked, the population Stage 2 corrects changes from 33,613 cells to ~75, and the reference-cell arithmetic must be redone |
>
> **The decisive number:** the two donors outside the clock's fitted age range have `dage_mae` **3×**
> the adults' and **conformal coverage of exactly 0.000** (`scorecard/baseline.json`). A per-donor
> offset correction is the right instrument for a level shift; it is the wrong instrument for
> extrapolation past a clock's fitted range.
```

### P-5 → `STAGE_5_PUBLICATION.md`, after line 180 (end of §2's ❌ list)

```markdown
> ## 🆕 ADDED 2026-07-31 — one claim is now DISPROVED, and three need a stated limit
>
> *Additive; §1's audit table and §2's list above are unmodified.*
>
> **Add to §2's ❌ list — this is stronger than "cannot be claimed", it is measured and false:**
>
> - ❌ **"the transcriptomic clock can be calibrated against methylation"** — tested directly on
>   paired samples and **refuted**. `STAGE_1_5_2_LABEL_ANCHOR.md` M-2a: ρ_partial **+0.267** and
>   **+0.516** against a pre-frozen 0.50 bar, on 68 paired conditions; SPLIT ⇒ NOT CALIBRATABLE.
>   The negative verdict passed its own falsification check (§12-R).
>
> **Limits that must travel with the ΔAge rows in §1:**
>
> | claim | limit to attach |
> |---|---|
> | anything quoting ΔAge **in years** | the label is RNA-clock-derived and not validated in absolute terms; the deployed cross-donor interval is **±34.64 yr** |
> | *"deep matches but does not beat linear for clock-ΔAge"* | still true, and now **narrower**: it is a statement about predicting *this clock's output*, which is not a statement about age |
> | *"donor-level calibration is the binding constraint (±12.7 yr)"* | true, and **not the only one** — the label's own validity is upstream of its calibration |
>
> **And one finding worth publishing that this creates:** a pre-registered, paired-ground-truth
> demonstration that a widely-used transcriptomic clock **cannot** be rescued by calibration on
> reprogramming cells — with the methylation anchor verified sound in the same experiment, and the
> per-arm ceiling (§17) showing the reference's own precision varies 4× with the axis under study.
> Negative, decisive, and reusable by anyone applying an RNA clock to reprogramming data.
```

### P-6 → `STAGE_6_NEW_DATA.md`, after line 35 (end of the prioritisation note)

```markdown
> ## 🆕 ADDED 2026-07-31 — both notes above are promoted from advice to REQUIREMENT, and sized
>
> *Additive; the acquisition checklist and the prioritisation note above are unmodified.*
>
> `STAGE_1_5_2_LABEL_ANCHOR.md` closed with the RNA-clock route **shut**: five repair attempts
> failed, the last being calibration against paired methylation. So the two notes above are no longer
> preferences.
>
> 1. **Methylation is not "worth more" — for the age arm it is the only thing that works.** RNA can
>    still serve the fate/safety head, whose labels do not run through a clock.
> 2. **`STAGE_1_5_1_REV_FINAL.md` §11.2 supplies the sizing**, which this stage previously lacked:
>    **≈16 transient-arm pairs** to settle the retention question — against the 9 that exist.
> 3. **Size for the instrument, not just for n.** 1.5.2 §12-R measured **±7 yr donor-level
>    methylation-clock error** (two donors of identical age 53 read 44.0 and 58.5). The retention
>    effect is −6 to −9 yr — **the same size as the error of the instrument measuring it.** More
>    donors improve the pairing; they do not sharpen the instrument. An acquisition that adds donors
>    without adding replicates per donor will not settle it.
> 4. **HFF methylation remains the single highest-value acquisition** — it is 99.7% of the age labels
>    and no public series has it (`REV FINAL` §11.3).
```

### P-7 → `REF_GROUND_RULES.md`, after line 118 (end of §5b)

```markdown
## 5b-bis. Two things §5b does not yet say — both learned in Stage 1.5.2

*Additive; §5b above is unmodified.*

**(i) A pre-committed FALLBACK is a bar, and must be resolvability-checked too.**
`STAGE_1_5_2` §6 anticipated that its primary criterion might be unresolvable and pre-registered a
fallback — which is good practice, and it registered the fallback **without checking it**. Measured,
**the fallback was itself UNRESOLVABLE at 92.3%**. On the originally-planned geometry the stage would
have had *no valid decisive criterion at all* and would not have known.

> **Every branch of a pre-registration is a bar. Check the ones you hope not to use.**

**(ii) A bar near the instrument's empirical ceiling is not a bar.**
§5b asks whether a system meeting the intent *exactly* passes. It does not ask whether **anything**
can reach the bar on this data. Stage 1.5.2 set ρ ≥ 0.50 against a null simulated at ρ_true = 0.70 —
and then measured that **two clocks of the same modality, on the same samples, agree with each other
at only ρ_partial +0.568** (§12-R), varying **+0.233 to +0.936** by cell state (§17).

> **Measure the ceiling — the agreement achievable between two instruments of the same kind on the
> same data — before setting a bar near it. A bar the reference standard cannot itself clear tests
> the data, not the candidate.**

Cheapest form: score a second instance of the *reference* modality by the identical criterion, on
the identical samples. If it fails too, the bar is the problem.
```

---

## PART C — G-c step 2, designed in full

This is the experiment everything above exists to enable. **It is not run by this stage**, but it is
specified here so C-1's design can be checked against its actual consumer.

### C.1 The question, and the two arms

> **Do HFF's 33,613 ΔAge labels help the age head, or does it learn artefact from them while 75
> usable labels are drowned out?**

| arm | config |
|---|---|
| **A (control)** | today's build — `AGE_MASKED_DATASETS = frozenset()` |
| **B (treatment)** | `AGE_MASKED_DATASETS = frozenset({"hff_sc"})` |

**One change between them.** `enforce_clock_age_range` stays `False` in both, so N2/N3 are in both
arms and the comparison isolates the HFF policy — the C-2 decision is a **separate** experiment and
must not ride along.

### C.2 The metric, chosen before the run

Judged on `scorecard.py`'s existing metrics, on the **non-HFF held-out donors**, because those are
the only cells with labels in both arms:

| role | metric | why |
|---|---|---|
| **primary** | `dage_mae_model` | the direct question: does the age head predict Gill's ΔAge better with or without HFF's labels? |
| **secondary** | `rank_model_dage` | 1.5.2 §17 and the baseline both show ranking survives where absolute values do not; a split between these two is itself informative |
| **guard** | `fate_prauc`, `fate_roc`, `fate_ece` | must **not** move. The fate head consumes no ΔAge, so a move means the change reached something it must not — ~~especially under C-5 Option 1's resampling~~ **(2026-08-02: Option 2 was chosen, so there is no resampling and this guard is now a strictly cleaner test — under Option 1 a move would have been ambiguous between the fix and the bootstrap; under Option 2 the fate head's input is untouched, so any move is unambiguous)** |
| **context** | `conformal_coverage` | already 0.000 on N2/N3; watch whether masking changes that |

**Baseline to beat (arm A, `scorecard/baseline.json`):** `dage_mae_model` = N2 21.79, N3 29.70,
O1 5.39, O2 7.54, Y1 7.28, Y2 14.06.

### C.3 The bar, and its resolvability check

`scorecard.py`'s comparison rule is a **paired 95% CI across folds excluding 0**, and
`audit_metrics.sensitivity_multiplier(6) ≈ 1.05` — so the minimum detectable mean effect is
`1.05 × SD(per-fold change)`.

**Before the run**, simulate the paired-difference null at 6 folds and confirm the intended effect
size is detectable; register the row in `tests/test_bars_resolvable.py`. **A pre-registered
three-way outcome, so no result can be reinterpreted afterwards:**

| outcome on `dage_mae_model` | reading | action |
|---|---|---|
| **B better**, CI excludes 0 | HFF's labels were **hurting** the age head | mask them; record that 99.7% of the labels were net-negative |
| **CI includes 0** | HFF's labels are **not contributing** | mask them anyway — the simpler model wins on a tie, and it removes an unvalidated 99.7% |
| **A better**, CI excludes 0 | HFF's labels **help**, artefact notwithstanding | **keep them**, and record the tension with G-c step 1 explicitly rather than explaining it away |

⚠️ **The tie rule is a judgement and is stated in advance** because it is the likely outcome and it
is where a post-hoc rationalisation would otherwise slip in.

### C.4 What it costs

`xdonor_calib.py`'s docstring: *"One extra ensemble per training donor — 5 for a 6-donor LOOCV fold,
so roughly 6× the training time."* Both arms are full LOOCV retrains ⇒ **≈2× a full retrain**, plus
two scorecard snapshots. Arm B is *cheaper* per fold than arm A on the age task (fewer labels), but
the fate task is unchanged, so plan for symmetric cost.

### C.5 The mandatory before/after record

`y_age` moves, so **every Stage 1 guard legitimately moves and the guard record restarts** — stated
in advance, per `STAGE_1_5_2_LABEL_ANCHOR.md` §8.2. Snapshot before, exercise the rollback rather
than assume it, and report:

* the deconfounder `coef` in both arms (C-5's second consequence);
* the age-valid label count per split in both arms;
* `n_age_labeled` from `dataset_summary.json`, plus the new `baseline_census` and
  `age_mask_reason` tallies.

---

## PART D — the change manifest

Every file this stage touches, so the diff can be checked against an expectation rather than read
cold.

| file | change | ~lines | class |
|---|---|---|---|
| `src/cellfate/common/constants.py` | `+ AGE_MASKED_DATASETS` | +5 | C-1 |
| `src/cellfate/data/aging.py` | `+ age_label_policy()`, `LinearClock.age_range`, `delta_age` returns 3-tuple | +55 / ~10 | C-1, C-2 |
| `src/cellfate/data/sources.py` | `+ DONOR_AGE_YEARS`, HFF `extra=` | +8 | C-3 |
| `src/cellfate/data/build_dataset.py` | `+ enforce_clock_age_range`, 2 call sites unpack 3 values, reason plumbed to `assemble_samples` | ~12 | C-1, C-2, C-6 |
| `src/cellfate/data/assemble.py` | pass `age_mask_reason` through | +3 | C-6 |
| `src/cellfate/common/schemas.py` | `+ Sample.age_mask_reason`, `+ AgeProvenance`, one validator line | +18 | C-4, C-6 |
| `src/cellfate/common/io.py` | 2 schemas + `shard_to_numpy` gain one column | +5 | C-6 |
| `src/cellfate/inference/schema.py` | `+ age_validated`, `+ age_basis` (defaulted) | +2 | C-4 |
| `src/cellfate/inference/service.py` | warning list instead of a single string | ~10 | C-4 |
| `src/cellfate/inference/predictor.py` | load `AgeProvenance` | +3 | C-4 |
| `src/cellfate/models/losses.py` | `+ huber_age_window()` (Σloss/Σcells over a window) | +20 | C-5 |
| `src/cellfate/training/train.py` | ~~sampler (**only if** C-5 Option 1)~~ ~~→ age-loss accumulation, W = 8~~ → **`_AgeWindow`, threshold on accumulated age CELLS, `W_max = 8`** (step 5c) | +75 | C-5 |
| `src/cellfate/training/train_model.py` | `+ age_window_k = 1` (OFF), `+ age_window_max_batches = 8` | +12 | C-5 |
| `tests/` | 4 new files + widened unpacking in 4 existing | +250 | all |

**Untouched:** `models/`, `training/train_model.py`'s calibration path, `evaluation/`, every
`experiments/diag_*` script, and every existing plan line.

**`res.py` is NOT in this table.** C-4 option (a) deliberately changes no arithmetic — the RES
formula at `res.py:38-41` is left exactly as it is, and only the *reporting* around it changes.
Option (c) would touch it, and is deferred to Stage 3.

---

## PART E — the commands, per step

```bash
# STEP 1-3: the code changes. After EACH one:
python -m pytest tests/ -q
python -m ruff check src/ tests/ scripts/ plan_tests/

# the bit-identity guard on real data (steps 1-3 must all pass this)
python plan_tests/verify_age_mask_identical.py "D:\Gill" "D:\GSE242423"

# STEP 4: C-4 option (a) + plan annotations -- both additive only
python -m pytest tests/test_inference.py -q   # the two new Response fields are DEFAULTED,
                                              # so every existing construction site still works
git diff --stat plans/          # every file must show insertions and ZERO deletions

# STEP 5: register C-5's bar BEFORE any retrain
python -m pytest tests/test_bars_resolvable.py -q

# STEP 6: G-c step 2 -- arm A then arm B, snapshot each
python retrain_stage1.py && python scorecard.py snapshot --tag gc2_A_keep_hff
#   ... set AGE_MASKED_DATASETS = {"hff_sc"} ...
python retrain_stage1.py && python scorecard.py snapshot --tag gc2_B_mask_hff
python scorecard.py compare gc2_A_keep_hff gc2_B_mask_hff
```

`plan_tests/verify_age_mask_identical.py` does not exist yet — **writing it is part of step 1**, and
it is the same shape as the G-a real-data check: run every Gill donor and one HFF chunk through
`delta_age` with defaults and assert `np.array_equal` against values captured before the change.

---

## PART F — rollback, per change

| change | rollback | data consequence |
|---|---|---|
| C-1, C-2, C-3 | `git revert` | **none** — defaults make them no-ops |
| C-4 (a) | `git revert` | none; bundles built meanwhile carry an extra file that is ignored |
| C-5 | `git revert` | none until a retrain runs |
| C-6 | `git revert` **+ rebuild** | shards written with the new column are unreadable by the old reader. This is the one irreversible-without-rebuild step, which is why it rides with C-1's rebuild |
| **G-c step 2** | restore the pre-run snapshot | `y_age` moves; **the guard record restarts**. Stated in advance, not absorbed |


---

## 3. Holes I looked for, and what I found

Stated including the ones that came back clean, because a check that found nothing is still evidence.

| # | worry | finding |
|---|---|---|
| 1 | masking HFF collapses the conformal calibration | ❌ **no** — it is already Gill-only, n=103, HFF fold skipped by `MIN_INNER_TRAIN_FRAC` (E16/E17). Checked, not assumed |
| 2 | masking HFF for being neonatal is inconsistent unless N2/N3 go too | ✅ **real** — 30 of the 75 remaining labels are neonatal (E12/E13). Addressed in C-2 |
| 3 | 75 labels over 33 688 cells starves the age head | ✅ **real and quantified** — 1.14 labels/batch, ~32 % of batches contribute nothing (C-5) |
| 4 | adding `age_mask_reason` breaks old shards | ✅ **real** — `io.py:194-228` reads fixed keys. Acceptable only because a rebuild is forced anyway (C-6) |
| 5 | the fate/safety head is affected | ❌ **no** — it consumes no ΔAge; `losses.py` keeps the two tasks separate and `age_mask` gates only the age term |
| 6 | `donor_age` could leak into the model as a feature | ❌ **no** — `schema.py:19` forbids extra fields and `tests/test_baseline_census.py` pins it |
| 7 | C-1 silently changes labels on the default config | ❌ **no by construction** — `AGE_MASKED_DATASETS` is empty by default and the guard is bit-identity |
| 8 | Stage 2 becomes pointless | ⚠️ **not pointless, but re-premised** — the offset may be real; what is unvalidated is the quantity it offsets. P-4 says exactly that rather than cancelling the stage |
| 9 | the ΔAge concept itself is dead | ❌ **no** — REV FINAL §4 measured real rejuvenation on methylation with an inert control. **The concept is vindicated; the RNA labels are not** (`REV FINAL` §4.5) |
| 10 | masking a cell leaves a stale `y_age` behind, tripping the schema | ❌ **no** — `assemble.py:44` is already `y_age=(float(y_age[i]) if masked else None)`, which is exactly what `schemas.py:126-130` requires in both directions. **C-1 needs no companion change here**, and that was checked rather than assumed |
| 11 | the cell-cycle deconfounder is unaffected | ✅ **real, and easy to miss** — `build_dataset.py:449-457` fits it on age-valid TRAIN cells, so masking HFF moves it from 33 613 single-cell to 75 **bulk** samples, past a guard that only checks `>= 2`. See C-5 |
| 12 | Stage 1.5.2's own `src/`-untouched guarantee is broken by this stage | ❌ **no** — 1.5.2 is closed. This is a **separate** stage whose whole purpose is the `src/` change, which is why it is a new document rather than an appendix to that one |

---

## 4. What this stage does NOT do

* **It does not run G-c step 2.** C-1 makes it *possible*; running it is a pre-registered experiment
  with its own bar, and it needs a rebuild + retrain.
* **It does not implement C-4 option (c).** Option (a) is decided and ships at step 4; making the
  *output* honest rather than the documentation is Stage 3's question, and this stage does not
  pre-empt it.
* **It does not touch the fate/safety head, the model architecture, or Stage 1's calibration.**
* **It does not change any label on the default configuration.** Every label-moving change ships
  behind a flag that is off.
* **It does not edit a single existing line of any plan.** PART B is annotation only.

---

## 5. Order of operations

The sequence is chosen so that **nothing that can move a label runs before the guard that proves it
did not**.

| step | action | gate before proceeding |
|---|---|---|
| step | action | gate before proceeding | moves a label? | rebuild? |
|---|---|---|---|---|
| **1** | C-6 (`age_mask_reason`), C-3 (HFF metadata) | full suite green; bit-identity asserted | ❌ no | ⚠️ **yes — shard schema gains a column** |
| **2** | C-1 (`AGE_MASKED_DATASETS`, empty) | **ΔAge and `age_mask` bit-identical**; `tests/test_data_units.py:246` assertions unmodified | ❌ no | rides step 1's |
| **3** | C-2 (`age_range`, flag off) | same bit-identity guard, flag off | ❌ no | rides step 1's |
| **4** | **C-4 option (a)** + PART B annotations (text in **B.2**) | `Response` gains two **defaulted** fields; annotations additive, `git diff --stat plans/` shows **zero deletions** | ❌ no | ❌ no |
| **5** ✅ | C-5 design + its bar — **DONE 2026-08-02, and it chose Option 2 over Option 1** | bar RESOLVABLE on the dense regime, DISCRIMINATES between options, 6 rows in `tests/test_bars_resolvable.py` | ❌ no | ❌ no |
| **5b** ✅ | deeper tests D1–D7 before committing — **DONE 2026-08-02: Option 2 CONFIRMED, W = 8 pinned** | ranking stable on 7 axes, not 1; W chosen for margin against a shrinking label set, not for the minimum | ❌ no | ❌ no |
| **5c** ✅ | **C-5 Option 2 IMPLEMENTED — DONE 2026-08-02, ships inert (`age_window_k = 1`)** — `training/train.py`, threshold on accumulated age cells with `W_max = 8`. Added 2026-08-02: no step previously scheduled it, and a fixed W biases arm A. See the readiness audit in C-5 | `k` registered via `bar_verdict`; **arm-A behaviour bit-identical to today**; window loss is `Σloss/Σcells` not a mean of means; fate term still steps every batch; determinism asserted | ❌ no | ❌ no |
| **6** | **G-c step 2** (PART C) | snapshot first; **rollback exercised, not assumed**; **5c must have shipped**; **`age_window_k = 4` set in BOTH arms** (see 🆕 below); **the arm-comparison bar registered and its MDE reported** | ✅ **yes** | ✅ yes, ×2 arms |

> ### 🆕 ADDED 2026-08-02 — two things step 6 was missing, found on review
>
> #### (a) `age_window_k = 4` must be set explicitly, in BOTH arms
>
> **5c ships inert at `age_window_k = 1`, and 1 means OFF.** The gate above previously said only
> *"5c must have shipped"*, and no command in PART E sets `k`. **Run as written, both arms would use
> `k = 1`, arm B would be starved, and problem #1 from the readiness audit would return silently** —
> the exact confound 5c was built to remove, reintroduced by a default.
>
> `k = 4` is not a new choice: it is **B2's registered value** (the smallest `k` that halves the
> per-update standard error, `1/√4 = 0.50`). Arm A holds ~512 age cells per batch, so it closes at
> `W = 1` and stays bit-identical to today; arm B accumulates. **One policy, both arms** — they
> differ only because the data differ.
>
> #### (b) The comparison itself had NO registered bar — now it does
>
> Every bar registered for this stage grades a **mechanism** (B1/B2 for C-5, A1/A2/A3 for C-5c).
> The comparison step 6 actually decides on — arm A vs arm B on `dage_mae_model`, paired across 6
> donor folds — had none. Ground rule §5b: *"a bar with no such test is not considered
> pre-registered."*
>
> Registered by **`plan_tests/register_gc_step2_bar.py`** → `results/register_gc_step2_bar_results.json`,
> with 3 rows in `tests/test_bars_resolvable.py`. **Δ\* = 3.57 yr**, derived from Stage 2 §12's
> existing *"≥ 25 % drop in `dage_mae_model`"* applied to the 14.29 yr recorded baseline — an
> established threshold for this same metric, not a number invented here.
>
> | SD(per-fold difference) | MDE | P(detect Δ\*) | |
> |---|---|---|---|
> | 0.5 yr | 0.52 | 1.0000 | ✅ RESOLVABLE |
> | **1.0 yr** | **1.05** | **1.0000** | ✅ **RESOLVABLE** |
> | 2.0 yr | 2.10 | 0.9338 | ❌ UNRESOLVABLE |
> | 3.0 yr | 3.15 | 0.6476 | ❌ UNRESOLVABLE |
> | 5.0 yr | 5.25 | 0.2955 | ❌ UNRESOLVABLE |
> | 13.7 yr *(arms independent)* | 14.38 | 0.0752 | ❌ **almost pure noise** |
>
> *(false-positive rate at a true effect of 0: **0.0508** — the CI is honest, it is only weak.)*
>
> **Δ\* is detectable at ≥95 % only if the two arms track each other to within ~1 yr per fold**, on a
> metric whose baseline already ranges 5.39 → 29.69 across folds. That is a demanding requirement and
> **it is not known to hold.**
>
> #### The reading rule, pre-registered because the null is the dangerous outcome
>
> The two directions are not symmetric. *"B better"* is self-limiting. **"CI includes 0" is the
> trap** — read as *"HFF's labels contribute nothing, discard them"*, it would throw away **99.7 % of
> the age labels on a null that may simply be underpowered.**
>
> | observed | licensed conclusion |
> |---|---|
> | CI excludes 0 **and** \|effect\| > MDE | conclusive in that direction |
> | CI includes 0 **and** MDE ≤ Δ\* | genuine null — the labels really do not help |
> | CI includes 0 **and** MDE > Δ\* | 🔴 **INCONCLUSIVE. Licenses NOTHING**, and specifically does **not** license discarding HFF's labels. Report as underpowered and say so |
>
> **The run must report its own observed SD and MDE beside the effect.** Without them the outcome
> table cannot be applied, because which row you are in depends on the MDE.

> **C-4 moved from step 7 to step 4 on 2026-07-31**, when option (a) was chosen. It sat last only
> because the choice between the three options depended on G-c step 2's result; **(a) has no such
> dependency** — it changes no arithmetic, adds only defaulted fields, and defaults to the
> conservative answer. Shipping it at step 4 means the deployed contract stops overstating ΔAge
> **before** the retrain rather than after it.

**Steps 1–5 cannot move a number.** ⚠️ **But step 1 is not free**: C-6 adds a column to
`_SHARD_SCHEMA`, so existing shards must be regenerated even though no *label* changes. That
regeneration must produce **bit-identical `y_age` and `age_mask`** — which is exactly what the
step-1 gate checks, and the reason C-6 is scheduled with C-1 rather than on its own (see C-6's
compatibility note and PART F).

**Step 6 is the first step that moves a label**, and it is the one Stage 1.5.2 pre-registered.
Commands for every step are in **PART E**; per-change rollback is in **PART F**.

---

## 6. Verification

| requirement | how |
|---|---|
| no label moves in steps 1–5 | `np.array_equal` on **ΔAge, `age_mask` and the new `age_mask_reason`** (which must be all-`None` at defaults), in a unit test **and** on all six real Gill donors plus one HFF chunk — the pattern G-a already used. Written as `plan_tests/verify_age_mask_identical.py` in step 1 — it is a per-stage gate, which is what that folder is for |
| the cancer rule is untouched | `tests/test_data_units.py:246` — only its tuple unpacking may widen to 3 values; **if an assertion needs changing, revert** |
| the rebuild in step 1 is faithful | rebuilt shards compared against the pre-change ones on every column except the new one |
| every new pure function is unit-tested with no repo data | the pattern of the five `diag_*` scripts |
| every new bar is registered | one row per bar in `tests/test_bars_resolvable.py`; a bar without one is not pre-registered (§5b) |
| full suite green | currently **564 passing** |
| lint | `ruff check src/ tests/ scripts/ plan_tests/` clean — **`plan_tests/` was added to CI's scope on 2026-08-01**, so the command is wider than it was when this stage was written. *(12 pre-existing errors remain in older `experiments/` files, which are deliberately still out of scope)* |
| the record | `CHANGES.md` + `experiments/DELTAAGE_LAB_NOTEBOOK.md`, prediction before result |

**Rollback.** Steps 1–5 are revertible by `git revert` with no data consequence. Step 6 changes
`y_age` and therefore every guard: it requires a full snapshot **and an exercised rollback**, per
`STAGE_1_5_2_LABEL_ANCHOR.md` §8.2. The guard record restarts, and that is stated in advance rather
than absorbed.

---

## 7. What would falsify this stage's own reasoning

* **C-1 is wrong** if some other mechanism can already mask HFF. Checked: `grep -rn "age_mask" src/`
  returns **39** sites, and `aging.py:219` is the only one that *assigns* it from data — every other
  site reads, stores, serialises or propagates the value.
* **C-4 is overstated** if RES does not actually depend on ΔAge. Checked: `g` is a multiplicative
  factor at `res.py:41` and `R_eff` at `res.py:38` has no other input.
* **C-5 is overstated** if the age head does not need many labels. Not claimed — the arithmetic
  states the batch-occupancy problem and stops there; whether 75 labels *suffice* is exactly what
  G-c step 2 measures, and this stage does not pre-judge it.
* **The whole stage is unnecessary** if Stage 1.5.2's verdict is wrong. That verdict passed its own
  falsification check (§12-R), and §17 re-audited every number in it against the artefacts.

---
