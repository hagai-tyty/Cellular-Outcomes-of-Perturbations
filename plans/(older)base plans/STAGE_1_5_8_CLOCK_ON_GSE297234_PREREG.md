# CLOCK ON GSE297234 (ages 22 vs 96) — PRE-REGISTRATION. Committed BEFORE the run.

**Status:** 🔒 pre-registration. Committed before `experiments/clock_gse297234.py` was executed.
**Decides:** whether our clock can order two donors **74 years apart**, and whether the **+30 yr
absolute bias** measured on `GSE165177` is cross-study, modality-specific, or source-specific.

---

## 1. Why this dataset, and why now

M-E1 found the clock reads fresh 38–53 yr fibroblasts as 72–82 — a **+30 yr floor**. Our own data
cannot diagnose it: 3 donors spanning 15 years, against a `cv_mae` of 12.27, gives
`P(correct ordering) = 0.755`. Nothing is testable there.

`GSE297234` (Lu et al. 2025) has **two untreated adult donors 74 years apart** — `GM23815` (22) and
`GM00731` (96) — same lab, same Sendai-OSKM protocol, 10x scRNA-seq. At a 74 yr gap the same error
model gives **`P(correct ordering) ≈ 0.9997`**. If the clock cannot order *those*, the problem is
not sampling noise.

### 🔑 The discriminating fact this dataset supplies

Fleischer's training lines are **Coriell** repository stock (`AG09599`, `AG04054`, …, NIA / NIGMS).
`GM23815` and `GM00731` are **also Coriell**. `GSE165177`'s cells were purchased from **Lonza**.

> **So if the +30 yr floor is driven by cell source / supplier / handling, this dataset should show
> a SMALLER bias. If the bias is unchanged, source is not the driver.** That is a real, falsifiable
> discrimination between hypotheses, and it is only possible because of the supplier difference.

---

## 2. What will be run

`experiments/clock_gse297234.py`, **READ-ONLY** with respect to the repo — no build, no retrain,
`src/` untouched. New external data, two files only:

- `GSM8986586_GM00731_D0_filtered_feature_bc_matrix.h5` (28 MB) — age **96**, day 0, untreated
- `GSM8986590_GM23815_D0_filtered_feature_bc_matrix.h5` (34 MB) — age **22**, day 0, untreated

from `ftp.ncbi.nlm.nih.gov/geo/samples/…`. The 266 MB `RAW.tar` and the 1.2/3.5 GB Seurat `.rds`
files are **not** downloaded — day 0 is all this test needs. New dependency: `h5py`.

**Clocking.** The clock is bulk-trained, so the primary estimate is **pseudobulk**: sum raw counts
across cells within a sample → `normalize_counts(target_sum=1e4)` → `log1p` → `predict_age`. This
is the same normalisation path as every other ΔAge in the project. A **per-cell** distribution is
reported as a secondary, clearly-labelled quantity, because per-cell dropout is a domain shift the
clock was never fitted for.

---

## 3. 🔒 PRE-REGISTERED OUTCOMES

### N1 — does the clock ORDER a 74-year gap correctly?

| result | meaning | action |
|---|---|---|
| **`pred(96) > pred(22)`** | the clock carries real age signal on this data | ordering is usable; calibration is the remaining problem |
| **`pred(96) ≤ pred(22)`** | ❌ the clock fails at ~4× its own error | **escalate.** Not a noise story. Any ΔAge interpretation that leans on the clock tracking age would need re-examining |

*This is a binary observation on n = 2 donors, not a significance test, and is reported as such.*

### N2 — absolute calibration at both ends

| result | meaning |
|---|---|
| both \|pred − true\| ≤ 12.27 (`cv_mae`) | **PASS-CALIBRATION** — absolute ages usable here |
| otherwise | **FAIL-CALIBRATION** — consistent with M-E1; absolute ages remain unusable |

### N3 — 🔑 is the +30 yr floor SOURCE-specific? *(the discriminating test)*

Compare the mean bias here against `GSE165177`'s day-0 bias of **+30.0 yr**.

| result | reading | action |
|---|---|---|
| mean bias **≤ +15 yr** (less than half) | **SOURCE MATTERS.** Coriell-to-Coriell transfers much better than Coriell-to-Lonza | supplier/handling becomes a first-class variable; the acquisition spec should require it be matched |
| mean bias **within ±10 yr of +30** | **SOURCE IS NOT THE DRIVER.** The offset is a general cross-study/platform failure | supports Gill's ComBat-harmonisation route as the only fix; drop supplier from the spec |
| mean bias **≥ +45 yr** | single-cell adds its **own** penalty on top | pseudobulk-from-scRNA is its own domain shift and must be handled separately |

### N4 — the SLOPE, which two points can give and three clustered donors cannot

`b = (pred(96) − pred(22)) / 74`.

| result | meaning |
|---|---|
| `b ≈ 1` (0.7–1.3) | the clock tracks age at the right **rate**; the failure is a pure **intercept** offset, which a single additive correction could fix |
| `0 < b < 0.7` | **compressed dynamic range** — a scale error as well as an offset; an intercept fix is not enough |
| `b ≤ 0` | no age signal at all on this data |

*Two points determine a slope exactly, so this has **no error bar**. A within-donor cell bootstrap
is reported alongside it, and it measures **pseudobulk stability only** — it is NOT a donor-level
interval and must never be quoted as one.*

### N5 — gene coverage, for comparison with `GSE165177`

Report clock genes present and the fraction of total |weight| reachable, against `GSE165177`'s
**57.1 % / 89.2 %**. A large difference would mean coverage is confounded with source in N3.

---

## 4. 🔒 PRE-REGISTERED EXPECTATIONS (stated so they can fail)

- **P-N1** the clock **will** order 22 < 96 correctly. At 74 yr this should be near-certain if the
  clock works at all.
- **P-N2** calibration **will** FAIL, with a positive bias, as it did in M-E1.
- **P-N3** the bias **will be smaller than +30** because these are Coriell lines like the clock's
  own training stock. *This is the prediction most worth being wrong about — if the bias is
  unchanged, the supplier hypothesis dies and ComBat-style harmonisation becomes the only route.*
- **P-N4** the slope **will be well below 1** — a clock this badly offset is more likely compressed
  than merely shifted.

## 5. Declared limits

1. **2 donors.** No donor-level confidence interval is possible. N1 is an observation, not a test.
2. **Pseudobulk from scRNA-seq is not bulk.** Depth, dropout and cell-composition differ from the
   clock's training domain; this is a stated confound in N3, not something corrected for.
3. **Day 0 only.** Nothing here touches reprogramming, ΔAge, or the safety target.
4. **Different protocol** from `GSE165177` (Sendai vs lentiviral) — irrelevant at day 0, since no
   factors have been delivered, but recorded.
5. **This does not re-open Stage 1.5** and produces no new labels.

## 6. Recording

`results/clock_gse297234_results.json`; write-up to `CHANGES.md` and the lab notebook; unit tests
in `tests/test_clock_gse297234.py`. Every bar in §3 and every expectation in §4 graded as written,
including those that fail.

---

## 7. RESULT — 2026-08-14. The clock **orders 74 years correctly but compresses them ~3.4×**, and the Coriell hypothesis is dead.

*Graded against §3/§4 as written. Artefacts: `experiments/clock_gse297234.py`,
`results/clock_gse297234_results.json`, `tests/test_clock_gse297234.py` (13 tests).*
**Run on the corrected normalisation** — the first version carried the same double-`log1p` bug
recorded as §9 of the ΔAge pre-registration.

| line | true age | cells | predicted (pseudobulk) | cell-bootstrap 95 % | bias |
|---|---|---|---|---|---|
| GM23815 | **22** | 7,782 | **84.1** | [83.6, 84.6] | **+62.1** |
| GM00731 | **96** | 5,021 | **106.1** | [105.5, 106.6] | **+10.1** |

| # | outcome | verdict |
|---|---|---|
| **N1** | ordering at a 74-yr gap | ✅ **CORRECT** — 106.1 > 84.1. The clock does carry real age signal |
| **N2** | absolute calibration | ❌ **FAIL** — biases +62.1 and +10.1 against `cv_mae` 12.27 |
| **N3** | is the +30 floor SOURCE-specific? | ❌ **SOURCE IS NOT THE DRIVER** — mean bias **+36.1** here vs `GSE165177`'s +30, inside the ±10 band |
| **N4** | slope | ⚠️ **0.297** — **COMPRESSED**, not merely shifted |
| **N5** | coverage | 62.8 % genes / **92.7 %** weight mass, vs 57.1 % / 89.2 % — comparable, so coverage does **not** confound N3 |

### 🔑 The two findings

**1. The Coriell hypothesis is dead — and it was mine.** `P-N3` predicted a smaller bias here
because `GM23815`/`GM00731` are Coriell stock like Fleischer's own training lines, while
`GSE165177` used Lonza. **The bias is not smaller — it is slightly larger (+36.1 vs +30).**
Supplier/source is **not** what drives the offset. That kills a hypothesis cheaply and removes
"match the supplier" from the acquisition spec. **It also leaves ComBat-style harmonisation — what
Gill actually did — as the only route on the table.**

**2. The clock COMPRESSES age, it does not merely shift it.** A slope of **0.297** means a real
74-year difference is rendered as **22 years**. An additive intercept correction — the obvious fix
for a pure offset — **cannot repair this**; the clock's *sensitivity* to age is down ~3.4× out of
domain, on top of a large positive offset.

**Independent support that this is real and not another artefact:** pseudobulk gives 0.297 and the
per-cell route gives **0.307** — two different aggregations agreeing closely. Before the
double-`log1p` fix they disagreed badly (0.160 vs 0.307), which is what a normalisation bug looks
like and what its repair looks like.

### The question this raises about ΔAge, stated and NOT answered

If the clock's sensitivity to *chronological* age is ~0.3 out of domain, is its sensitivity to
*reprogramming-induced* change attenuated too? **This run cannot say** — between-donor
chronological discrimination and a within-dataset same-day contrast are different quantities. But
it is now a live question about the magnitude of every ΔAge this project reports, and it is
recorded as open rather than resolved.

### Expectations, graded: 3 of 4

| | expectation | held? |
|---|---|---|
| **P-N1** | orders 22 < 96 correctly | ✅ YES |
| **P-N2** | calibration FAILS with a positive bias | ✅ YES |
| **P-N3** | bias smaller than +30 (Coriell hypothesis) | ❌ **NO** — the one flagged as most worth being wrong about |
| **P-N4** | slope well below 1 | ✅ YES — 0.297 |

### What is NOT claimed

That 2 donors establish a slope with any uncertainty — two points fix a line exactly and the
cell-bootstrap measures **pseudobulk stability only**, never donor-level error. That the compression
factor transfers to ΔAge. Or that ordering 22 vs 96 correctly means the clock can order donors
closer together — at 15 years apart it demonstrably cannot.
