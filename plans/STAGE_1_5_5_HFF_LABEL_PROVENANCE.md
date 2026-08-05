# STAGE 1.5.5 — What IS HFF's ΔAge label? The within-timepoint provenance test

**Status:** ✅ **EXECUTED 2026-08-03.** Primary verdict: **NOT_DOMINATED** — HFF's cell-level ΔAge is
**not** an identity readout. See §5.

**Implements:** the question arm C left open — it proved the labels carry *"consistent exploitable
structure of unknown provenance"* (STEP6 report). This measures the provenance.
**Scope:** 1 new script, 1 test file, **0 lines changed in `src/`**, **no label moves under any
outcome.**

---

## 1. Why this exists

HFF supplies **33,613 of the project's 33,688 age labels — 99.8 %** — and what they *mean* has never
been measured. Two established facts pointed opposite ways:

| | |
|---|---|
| **arm C** (step-6 follow-up) | permuting HFF's labels collapsed `rank_model_dage` **0.9476 → 0.5765**, **5.4× the entire A–B gap**. The labels carry structure the model exploits — provenance unknown |
| **G-c step 1** | the trajectory has the *shape* of real rejuvenation: ρ(day, ΔAge) = **−0.905**, stable under leave-one-timepoint-out, reaching **−24.0 yr** against methylation's verified **−24.1**. But the identity artefact produces that same monotone shape, so shape cannot separate them |

The artefact has a known signature — `corr(age, pluripotency) = −0.62` and
`corr(age, fibroblast identity) = +0.62` in Gill, and **+36.5 yr on cells that could not have aged.**
So the question is answerable:

> **Within a single timepoint, how much of a cell's ΔAge is explained by its identity markers?**

**Within-timepoint is the entire design.** Between timepoints ΔAge and pluripotency move together
trivially — every cell is further along — so a pooled correlation is uninformative *by construction*.
Holding day constant asks the discriminating question: among cells at the **same** point in the
protocol, does the one that looks more pluripotent read younger?

---

## 2. Why this one could resolve when nothing else did

Every dead end in this project is **donor**-limited — 3 donors, 6 folds, MDE = 1.049 × SD, M3's
95 % CI of *[9 %, 100 %]*. **HFF is the one place that is not.**

**44,473 cells across 8 timepoints, ~5,800 per timepoint.** SE(ρ) ≈ 0.013 — roughly **two orders of
magnitude** more power than anything else in the dataset. Any bar in [0.1, 0.9] is trivially
resolvable, which is why no separate resolvability simulation was needed. It is the only
well-powered question left in the data, and nobody had asked it.

---

## 3. Design and pre-registration

Registered in the script docstring before the first run.

| | |
|---|---|
| **primary** | Spearman(`y_age`, pluripotency) **within each timepoint** |
| **bar** | \|ρ\| ≥ **0.50** in at least **6 of 8** timepoints ⇒ `IDENTITY_DOMINATED` |
| **also** | R² of `y_age` on [pluripotency, somatic identity] within timepoint — the fraction of the label that *is* identity |
| **reported** | every timepoint, always; none selected on its value |

Pluripotency is the existing `OSKM_PLURIPOTENCY` signature (18 genes present) and somatic identity is
`DEFAULT_SIGNATURES["safe"]` (4 genes) — both **reused verbatim** so neither can be tuned for this
stage. iPSC (day 21) is excluded: it is a cell-type change, not a point on the trajectory.

**Honest limit, stated before the run:** a *real* rejuvenation signal would **also** correlate with
pluripotency at cell level, because a cell further along is both more pluripotent and more
rejuvenated. So a strong ρ would not by itself prove artefact. What the R² adds is the quantity that
does matter — **if identity explains essentially all the within-timepoint variance, there is nothing
left for age to be.**

---

## 4. What each outcome licenses

| | |
|---|---|
| **IDENTITY_DOMINATED** | the label is an identity readout at cell level; arm C's "exploitable structure" is pluripotency structure, and calling it age is unjustified |
| **NOT_DOMINATED** | ΔAge carries within-timepoint structure identity does not explain. That residual is the *candidate* age signal — **candidate, not confirmed** |

**Neither outcome moves a label.** This measures what the existing labels are.

---

## 5. ✅ RESULT — `NOT_DOMINATED`

**0 of 8 timepoints reach the 0.50 bar.**

| day | n | ρ(age, pluripotency) | ρ(age, somatic) | **R² identity** |
|---|---:|---:|---:|---:|
| 0 | 5,981 | +0.079 | +0.171 | 0.033 |
| 2 | 3,947 | −0.079 | +0.265 | 0.082 |
| 4 | 5,735 | +0.036 | +0.111 | 0.034 |
| 6 | 5,632 | −0.158 | +0.281 | 0.081 |
| 8 | 5,663 | −0.130 | +0.133 | 0.020 |
| 10 | 5,903 | −0.128 | +0.251 | 0.076 |
| 12 | 5,830 | −0.270 | +0.397 | **0.162** |
| 14 | 5,782 | −0.194 | +0.256 | 0.068 |

**Identity explains 2–16 % of within-timepoint ΔAge variance, typically ~7 %. So 84–98 % of it is
something else.**

### 5.1 The confound the design was built to remove, visible in one number

**Pooled** ρ(age, pluripotency) = **−0.216**, against within-timepoint values near zero. Between
timepoints they move together; hold day constant and the association largely vanishes. Reported as
**descriptive only** and never graded — pooling here would have manufactured the artefact the test
exists to detect.

### 5.2 What this establishes, and what it does not

**Establishes:** the identity artefact — the mechanism behind +36.5 yr on non-responders and
`corr(age, pluripotency) = −0.62` in Gill — **is not what HFF's cell-level ΔAge is made of.** That
hypothesis is now measured and rejected rather than argued about.

**Does not establish that the residual is age.** The honest alternative is **clock noise**: `cv_mae`
is 12.27 yr and single-cell profiles are sparse. *"Not identity"* is not *"real."* §6 tests the
next-most-likely mundane explanation.

### 5.3 A second-order observation

ρ(age, **somatic**) is consistently **positive** (+0.11 to +0.40) and larger in magnitude than the
pluripotency term at most timepoints — more fibroblast-like cells read *older*. Same direction as
1.5.1's `corr(age, fibroblast identity) = +0.62`, but far weaker once day is held constant.

---

## 6. Is the residual structured, or is it noise?

The obvious mundane explanation for 84–98 % unexplained variance in single-cell data is **technical**:
library depth, detected genes, mitochondrial fraction. So the same within-timepoint decomposition is
run against those covariates, and against identity + technical combined.

| | reading |
|---|---|
| ΔAge tracks library depth / gene count | **technical artefact** — the label is a sequencing-depth readout |
| neither identity nor technical explains it | the residual is neither, and arm C showed a model can **exploit** it — which noise, by definition, cannot be |

### 6.1 ✅ Result — it is neither

| day | R² identity | R² technical | **R² both** | ρ(age, lib) | ρ(age, n_genes) |
|---|---:|---:|---:|---:|---:|
| 0 | 0.033 | 0.086 | 0.100 | +0.211 | +0.192 |
| 2 | 0.082 | 0.075 | 0.123 | +0.201 | +0.241 |
| 4 | 0.034 | 0.030 | 0.054 | +0.147 | +0.159 |
| 6 | 0.081 | **0.008** | 0.083 | +0.009 | +0.023 |
| 8 | 0.020 | **0.007** | 0.026 | −0.002 | +0.014 |
| 10 | 0.076 | 0.030 | 0.114 | +0.168 | +0.168 |
| 12 | 0.162 | **0.013** | 0.167 | −0.006 | +0.014 |
| 14 | 0.068 | 0.010 | 0.112 | +0.066 | +0.057 |

**Identity and technical covariates together explain 2.6–16.7 % of within-timepoint ΔAge variance —
so 83–97 % is explained by neither.**

The label is **not** a sequencing-depth readout. And the two effects are largely independent rather
than the same thing twice: at days 6, 8 and 12 the technical R² is essentially **zero**
(0.008 / 0.007 / 0.013) while identity is at its strongest (ρ −0.158 / −0.130 / −0.270).

### 6.2 Where that leaves the label

Two leading mundane explanations for HFF's ΔAge are now **measured and rejected**:

| candidate | status |
|---|---|
| identity artefact — the mechanism behind +36.5 yr and `corr(age, pluri) = −0.62` | ❌ **rejected**, R² ≤ 0.16 |
| technical / sequencing depth | ❌ **rejected**, R² ≤ 0.09 |
| clock noise | ⚠️ **still live** |
| real biological signal | ⚠️ **still live** |

**This does not establish the label is age.** It establishes that the two cheapest explanations do
not account for it, which is the first positive characterisation these labels have ever had.

### 6.3 Two limitations that bound the claim

1. **Identity is measured by 22 genes** — 18 pluripotency + 4 somatic. That is a thin proxy. The
   strong form of the identity hypothesis is not "these genes" but "the cell's position on the
   reprogramming manifold", which a small marker set under-measures. A richer axis (leading
   transcriptome PCs) could explain materially more, and this result does **not** rule that out.
2. **R² is linear on ranks.** It is monotone-robust but would miss a genuinely non-monotone
   relationship between ΔAge and identity.

Both push the same way: **the true identity share is a floor, not a ceiling.** The honest reading is
"identity does not dominate", not "identity is absent".

---

## 7. Artefacts

| file | role |
|---|---|
| `experiments/diag_hff_label_identity.py` | the measurement; read-only; pure logic separated from I/O |
| `tests/test_diag_hff_label_identity.py` | 10 tests, no repo data required, including a regression pinning the observed values against the bar |
| `results/diag_hff_label_identity_results.json` | full per-timepoint output |

**One implementation note worth recording:** the source's `cells_per_run=None` puts every cell in a
single densified chunk — ~48,000 × 36,601 × 4 bytes ≈ **7 GB** — and the first attempt simply
thrashed. Its own docstring says the parameter "bounds peak RAM". Set to 4,000. Batching cannot
change a result here because every statistic is computed **per timepoint after the chunks are
concatenated**, so where the boundaries fall is irrelevant.
