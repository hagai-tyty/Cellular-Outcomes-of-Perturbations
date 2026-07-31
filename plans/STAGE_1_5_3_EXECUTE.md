# STAGE 1.5.3 — EXECUTE: the code changes Stage 1.5.2 forces

**Status:** 🔵 **PRE-REGISTERED — NOT EXECUTED.** No `src/` file is edited by this document.
**Depends on:** `STAGE_1_5_2_LABEL_ANCHOR.md` (✅ closed) and `STAGE_1_5_1_REV_FINAL.md` §11 (✅ closed).
**Blocking for:** Stage 2's premise, Stage 3's recommendation rule, Stage 5's claims.
**Not blocking for:** the fate/safety head, which consumes no ΔAge and is untouched throughout.

**Scope:** 6 `src/` changes (**PART A**) and 7 downstream plan annotations (**PART B**).
One of the six is blocking, one is a decision only the project owner can make, four are mechanical.

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

---

## 2. Ledger of changes

| # | change | class | why now |
|---|---|---|---|
| **C-1** | `age_mask` must be able to address **HFF specifically** | 🔴 **BLOCKING** | G-c step 2 is unrunnable without it |
| **C-2** | read the clock's declared `age_range` instead of discarding it | 🟠 substantive | E4/E5/E13 — the range exists, is ignored, and would change which labels are valid |
| **C-3** | stamp `donor_age` / `batch` on HFF too | 🟡 mechanical | G-b reached Gill and not HFF (E6); C-2 cannot fire for HFF without it |
| **C-4** | give the deployed response a way to say **"ΔAge is not validated here"** | 🔴 **DECISION REQUIRED** | E7–E10 — the product's score is multiplied by an unvalidated term and cannot report that |
| **C-5** | make the age head **and the cell-cycle deconfounder** survive a 450×-smaller label set | 🟠 substantive | E11/E14/E15 + `build_dataset.py:445-451` — arithmetic below |
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

**`src/cellfate/data/aging.py:219`** — `source` stays the first gate so the cancer rule is
untouched; `dataset_id` is consulted only when the column exists:

```python
valid = np.full(age.shape[0], source not in C.CANCER_SOURCES, dtype=bool)
if "dataset_id" in obs.columns and C.AGE_MASKED_DATASETS:
    valid &= ~obs["dataset_id"].isin(C.AGE_MASKED_DATASETS).to_numpy()
age_mask = valid
```

### The guard, and why it is credible

**With `AGE_MASKED_DATASETS` empty — the default — `age_mask` is bit-identical to today's**, because
the second line is a no-op. That is the same record-only discipline G-a shipped under, and it is
testable the same way: `np.array_equal`, not `allclose`.

`tests/test_data_units.py:246` (`test_delta_age_masks_cancer_sources`) pins the existing behaviour
and **must keep passing unmodified**. If it needs editing, the change is wrong.

### New tests

| test | asserts |
|---|---|
| default policy is a no-op | `age_mask` bit-identical with and without the new branch |
| a masked dataset is masked | `AGE_MASKED_DATASETS={"hff_sc"}` ⇒ HFF False, Gill True, **in the same chunk** |
| the cancer rule still wins | a Tahoe cell stays masked regardless of `dataset_id` |
| an absent column is not an error | a source that never stamps `dataset_id` behaves as today |

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

### The change

**`aging.py`, `LinearClock`** — carry the range instead of dropping it:

```python
def __init__(self, weights, intercept=0.0, age_range=None):
    ...
    self.age_range = tuple(age_range) if age_range else None

@classmethod
def from_json(cls, path):
    ...
    return cls(weights, intercept=float(d.get("intercept", 0.0)),
               age_range=(d.get("meta") or {}).get("age_range"))
```

**`aging.py:219`, in the same expression as C-1** — a donor outside the range is out of the clock's
own stated scope:

```python
if clock_range is not None and "donor_age" in obs.columns:
    a = pd.to_numeric(obs["donor_age"], errors="coerce").to_numpy(dtype=float)
    # NaN (unknown age) does NOT mask -- absence of evidence is recorded, not acted on
    valid &= ~((a < clock_range[0]) | (a > clock_range[1]))
```

### ⚠️ Decision this forces, stated rather than assumed

Enabling the range check changes labels **even with `AGE_MASKED_DATASETS` empty**, because N2/N3
would be masked. **Therefore it ships behind its own flag** (`DataConfig.enforce_clock_age_range`,
default `False`), and turning it on is a pre-registered change with its own bar — not a side effect
of C-1.

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

**That assertion belongs in the source, visibly, with its provenance** — not buried as a default:

```python
# HFF is a neonatal foreskin fibroblast LINE: GSE242423 declares no per-sample donor age
# because there is no per-sample donor. 0.0 is asserted from the line's identity, and is
# recorded here rather than defaulted silently, because it is the value that puts HFF
# outside the clock's fitted range of [1, 96] (configs/clocks/fleischer_clock.json).
extra={"donor_age": [0.0] * len(idx), "batch": [""] * len(idx)}
```

### Guard

`donor_age` and `batch` are metadata and must never become model input. `schema.py:19`
(`model_config = ConfigDict(extra="forbid")`) already forbids it at the deployment boundary, and
`tests/test_baseline_census.py` pins that. C-3 adds no new exposure.

---

## C-4 🔴 DECISION REQUIRED — the product multiplies by a term it cannot vouch for

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

**Recommendation: (a) now, and (c) considered at Stage 3**, because Stage 3 owns the decision rule
and already has an open design question about declining to report ΔAge (`STAGE_3_TOOL.md:24-26`).
**(b) is not recommended**: adding a second status that means "we don't know" next to one that means
"no effect" — when the existing one *already* fires for the "we don't know" reason — makes the
contract harder to read, not easier.

**This is the one item in this stage that is a product decision, not an engineering one. It is
listed, not decided.**

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

| option | mechanism |
|---|---|
| **age-aware batch sampling** | ensure ≥ k age-valid cells per batch (a `WeightedRandomSampler` or a two-stream loader) |
| **accumulate** | compute the age loss over a window of batches rather than per batch |
| **fixed task weight** | replace the learned `s_age` with a constant when the label count is below a threshold |

**Whichever is chosen must be justified before the run, and the metric it is judged on must go
through `audit_metrics.bar_verdict` first** — `REF_GROUND_RULES.md` §5b. Given four hairline margins
already on this project's record, this is not the place to skip that.

### 🔴 The second consequence, which is easy to miss: the cell-cycle deconfounder moves too

`build_dataset.py:445-451` fits the deconfounder on **age-valid TRAIN cells only**:

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

`coef` is then applied to **every** shard (`build_dataset.py:456`), including the masked HFF cells,
whose `y_age` is subsequently NaN'd (`457-460`) — so the masked cells are unaffected in the end, but
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

### The change

Add `age_mask_reason: str | None` to `Sample` (`schemas.py:55`) and to the shard schema
(`io.py:134`) and manifest schema (`io.py:236`), `None` when `age_mask` is True.

### ⚠️ The compatibility consequence, stated up front

`io.py:194-228` (`shard_to_numpy`) reads a fixed key list from `table.to_pydict()`, so **adding a
column makes previously-written shards unreadable by the new code without a migration**.

**This is acceptable here and only here**, because C-1/C-2 change the labels and therefore *force* a
rebuild anyway. **It would not be acceptable as a standalone change**, and it must not be shipped in
a release that does not already rebuild.

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
| 11 | the cell-cycle deconfounder is unaffected | ✅ **real, and easy to miss** — `build_dataset.py:445-451` fits it on age-valid TRAIN cells, so masking HFF moves it from 33 613 single-cell to 75 **bulk** samples, past a guard that only checks `>= 2`. See C-5 |
| 12 | Stage 1.5.2's own `src/`-untouched guarantee is broken by this stage | ❌ **no** — 1.5.2 is closed. This is a **separate** stage whose whole purpose is the `src/` change, which is why it is a new document rather than an appendix to that one |

---

## 4. What this stage does NOT do

* **It does not run G-c step 2.** C-1 makes it *possible*; running it is a pre-registered experiment
  with its own bar, and it needs a rebuild + retrain.
* **It does not decide C-4.** That is a product decision and is presented as three costed options.
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
| **1** | C-6 (`age_mask_reason`), C-3 (HFF metadata) | full suite green; no label moves — assert bit-identity |
| **2** | C-1 (`AGE_MASKED_DATASETS`, empty) | **ΔAge and `age_mask` bit-identical**; `tests/test_data_units.py:246` unmodified and passing |
| **3** | C-2 (`age_range`, flag off) | same bit-identity guard with the flag off |
| **4** | PART B annotations | additive only; `git diff` shows no deletions |
| **5** | C-5 design + its bar through `audit_metrics.bar_verdict` | bar RESOLVABLE **before** any retrain |
| **6** | **G-c step 2** — the retrain comparison | snapshot first; rollback path exercised, not assumed |
| **7** | C-4 | after step 6, because the answer changes which option is right |

**Steps 1–4 cannot move a number.** Step 5 is planning. **Step 6 is the first step that moves
anything**, and it is the one Stage 1.5.2 pre-registered.

---

## 6. Verification

| requirement | how |
|---|---|
| no label moves in steps 1–4 | `np.array_equal` on ΔAge and `age_mask`, in a unit test **and** on all six real Gill donors — the pattern G-a already used |
| the cancer rule is untouched | `tests/test_data_units.py:246` passes **unmodified** |
| every new pure function is unit-tested with no repo data | the pattern of the five `diag_*` scripts |
| every new bar is registered | one row per bar in `tests/test_bars_resolvable.py`; a bar without one is not pre-registered (§5b) |
| full suite green | currently **564 passing** |
| lint | `ruff check src/ tests/ scripts/` clean. *(12 pre-existing errors in four older `experiments/`+`tests/` files are not introduced here and are not in scope)* |
| the record | `CHANGES.md` + `experiments/DELTAAGE_LAB_NOTEBOOK.md`, prediction before result |

**Rollback.** Steps 1–4 are revertible by `git revert` with no data consequence. Step 6 changes
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
