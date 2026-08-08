# CHANGE C-7 — Bulk sample integrity: reject degenerate columns before they become labels

**Status:** 🔵 **PRE-REGISTERED 2026-08-08. NOT IMPLEMENTED. `src/` untouched, no label moved.**

**Why a file and not a section:** two machines pushed a `## 5.8` into `STAGE_1_5_6_SPARSE_CLOCK.md`
concurrently and neither saw the other's. This change is not part of 1.5.6 — it is a data-integrity
gate that 1.5.6 *surfaced* — and it takes the next free **change ID** (C-1 … C-6 are Stage 1.5.3's)
rather than a section number, so it cannot collide.

**Owner:** whoever holds `src/`. **Cost:** one assertion + tests to write; **one rebuild to adopt.**
**Blocking:** nothing. **Blocked by:** nothing. It can be written before 3b/3c report.

---

## 1. The defect, verified twice on two machines

`GSE165176_Log2_RPM_Sendai_reprogramming` contains columns that are **not transcriptomes**.

| column | log2 range | linear RPM sum | role |
|---|---:|---:|---|
| `Y1_d7_CD13_Sendai_Exp1` | **0.15** | **2.148e+09** | perturbation |
| `N3_d21_SSEA4_Sendai_Exp2` | 2.47 | 2.379e+08 | perturbation |
| `O2_d9_SSEA4_Sendai_Exp1` | 2.15 | 1.628e+08 | perturbation |
| **`N2_Fib_Sendai_Exp2`** | **1.74** | **1.030e+08** | 🔴 **CONTROL** |
| `N2_d21_CD13_Sendai_Exp2` | 7.26 | 1.694e+07 | perturbation |
| *— the other 119 —* | **9.00 – 15.26** | **2.859e+05 – 3.880e+06** | |

`Y1_d7_CD13` is **entirely constant** across ~20 000 genes. `N2_Fib` sits **0.0008 log2 above its own
floor** — a mean that close to the minimum is something RNA-seq cannot produce, because the
highly-expressed minority always pulls the mean well clear.

**`apply_qc` passes every one of them.** Its gates are `min_genes` and `max_mito_frac`
(`QCConfig`) — single-cell gates that a constant bulk column clears trivially.

**Reach.** `gill_bulk` is *both* the harmonization reference *and* a training source, so these
columns are simultaneously (a) age-labelled training rows and (b) inputs to `σ_ref`. `N2_Fib` is
`is_control`, so it enters `Harmonizer.fit` in **five of six LOOCV folds** — including **O1**, the
fold behind July's −24.02 (`STAGE_1_5_6_SPARSE_CLOCK.md` §5.7, §5.11).

**Stage 1.5.2's gate G-a made `n = 1` visible. Nothing checks whether that `n = 1` is sound.**

---

## 2. The gate — two conditions, both justified by units, not by this cohort's quantiles

An earlier proposal thresholded `mean − min` at ⅕ of the cohort median. **That was rejected and the
rejection is recorded** (§5.10): it cuts a continuous distribution 8 % from its neighbour, so on a
new cohort it would flag or miss arbitrarily. Both conditions below come from what the numbers
*mean*, and this cohort only confirms that they separate.

### G1 — library size

> The matrix is **Reads Per Million**. A sound column's linear values must therefore sum to ≈ **1e6**
> *by definition*. Accept `[1e5, 1e7]` — a decade either side of the value the units mandate.

| | |
|---|---|
| the 5 degenerate | 1.694e+07 – 2.148e+09 — **all above** |
| the other 119 | 2.859e+05 – 3.880e+06 — **all inside**, and all within 3.9× of 1e6 |
| margin below the ceiling | 1e7 / 3.880e+06 = **2.58×** |
| margin above it | 1.694e+07 / 1e7 = **1.69×** |

### G2 — dynamic range

> Any real transcriptome spans several orders of magnitude between its least- and most-expressed
> gene. Require **log2(max) − log2(min) ≥ 8**, i.e. at least a **256-fold** spread.

| | |
|---|---|
| the 5 degenerate | 0.15 – 7.26 — **all below** |
| the other 119 | 9.00 – 15.26 — **all above** |
| separation | **no overlap** |

### Why both

Each condition alone flags exactly the same 5 with 0 false positives on this cohort. **They are kept
as two because they fail differently:** G1 catches a mis-scaled library, G2 catches a collapsed
distribution. A future cohort that defeats one is unlikely to defeat both. **A column must satisfy
both to be admitted.**

### 🔴 What the gate does NOT catch, stated plainly

Seven further columns look poor on `mean − min` but pass both G1 and G2: `O2_d40`, `O2_d34`,
`O1_d34`, `N2_d11_CD13_Exp2`, `Y2_d34`, `O1_d11_CD13_Exp2`, and `Y1_Fib`.

**`Y1_Fib` is a control.** Its library (1.51e+06) and range (14.43) are entirely normal and only its
`mean − min` is low, on a continuum (0.2745 → 0.2967 → 0.3069). §5.10 downgraded it from "defective"
to "not established" and that downgrade stands. **C-7 does not resolve these seven.** They are
recorded as open, and Y1's unexplained floor ratio (§5.5) stays unexplained.

---

## 3. 🔴 The consequence that makes this more than one assertion

**`N2_Fib` is N2's only control.** Rejecting it leaves donor N2 with **no control at all**, and
`aging.py:88` then does this:

```python
ref = values[ctrl] if ctrl.any() else values[in_line]     # <- silent zero-point switch
```

**N2 would fall through to self-centring — which subtracts its own mean perturbation effect and
forces its mean ΔAge toward 0, with no warning, no counter and no mask.** That is the exact
behaviour Stage 1.5's Group D pinned as a defect. **Rejecting the sample without deciding the donor
would trade a known-bad control for a silent fallback, which is worse: the first is visible.**

### The three options, and the recommendation

| | consequence | |
|---|---|---|
| **(a) reject the sample AND the donor** | a donor with no sound control cannot carry control-relative ΔAge. Drops N2's **21** columns | ✅ **recommended** — consistent, and nothing lands in the fallback |
| (b) reject the sample, let the fallback fire | N2's ΔAge forced toward 0, silently | ❌ Group D pinned this |
| (c) reject it from `σ_ref` but keep it as N2's baseline | the same column is "too broken to estimate variance, sound enough to define zero" | ❌ incoherent |

**Under (a) the corpus goes 124 → 100 Gill columns and 6 donors → 5.**

### And that reaches further than this change

* **LOOCV goes from 6 folds to 5.** Every Stage 1 guard, and step 6's MDE arithmetic, are computed
  over folds. **This change cannot be adopted without re-reporting them.**
* **It reaches C-2.** N2 is donor age 0. If N2 leaves on integrity grounds, C-2's *"masks the two
  neonatal donors"* becomes *"masks N3"* — plus HFF (`00_START_HERE.md`, C-2 section).
* **It reaches §4.7's whole finding.** The 16.67 yr fold spread is defined over six folds, one of
  which is N2's. Removing N2 does not "fix" that spread — **it removes the fold the spread was
  measured on.** 3b and 3c must therefore report **before** C-7 is adopted, or their question
  disappears rather than gets answered.

> **Sequencing consequence, and it is not negotiable: C-7 is WRITTEN now and ADOPTED after 3b and
> 3c report.** Writing the gate costs nothing and cannot be undone by their results; adopting it
> first would delete their evidence.

---

## 4. Pre-registered bars

| | bar |
|---|---|
| **B1 — separation** | on the Gill matrix, G1 ∧ G2 flag **exactly the 5** columns in §1, with **0** of the other 119 flagged |
| **B2 — no silent fallback** | after rejection, **no `cell_line` reaches `_control_baseline` with zero controls.** Asserted, not logged: a donor losing its last control must **raise**, per the `age_label_policy` fail-open precedent |
| **B3 — the gate can fail** | fed a synthetic constant column, the check must reject it; fed a synthetic sound column, it must admit it. **Both branches execute in the test suite** — a branch that never runs is not a check (the `verify_1a` lesson) |
| **B4 — bit-identical when off** | with the flag off, `ΔAge` is **bit-identical** to the current build. The gate ships **off** and is enabled by its own pre-registered run, exactly as C-2 did |
| **resolvability (§5b)** | B1 is a **deterministic** classification on a fixed matrix, not an estimate — there is no sampling null, so `bar_verdict` records it as *deterministic, resolvability N/A* rather than simulating one. **B2/B3 are likewise deterministic.** No bar here needs a power calculation, and claiming one would be theatre |

---

## 5. Decision branches, fixed in advance

| outcome | meaning | what happens |
|---|---|---|
| **B1–B4 all pass** | the gate is correct and inert until switched on | ships **off**. Adoption is a separate run **after 3b/3c report**, with a snapshot and every Stage 1 guard re-reported over **5** folds |
| **B1 fails — extra columns flagged** | a threshold is wrong, or the cohort is not what §1 measured | **do not widen the band to fit.** Re-derive from the units or record that the units do not separate this cohort |
| **B2 fails** | a donor reaches the fallback | **blocking.** The donor-level decision (§3 option a) is not optional; fix that first |
| **B4 fails** | the gate is not inert when off | **blocking, and a bug in this change** — nothing about a disabled flag may move a label |

---

## 6. What this change does NOT license

* **It is not a re-analysis.** §§4.5, 4.6, 1c and 1d were computed on a contaminated `σ_gill`
  (§5.11) and **C-7 does not re-measure them.** That costs one HFF stream and belongs to 3b/3c.
* **It does not withdraw anything.** No recorded result is retracted on the strength of this.
* **It does not decide whether the defect is GEO's deposit or our read of it.** Either way the
  pipeline consumes it. If it is our read, C-7 becomes redundant rather than wrong — **and that is
  worth ten minutes against the GEO supplementary file before the gate is implemented.**
* **It does not touch HFF.** The gate is for bulk columns; the single-cell path has `apply_qc`.

---

## 7. Artefacts

| file | role |
|---|---|
| `plans/STAGE_1_5_6_SPARSE_CLOCK.md` §5.7 | the discovery (second machine) |
| `plans/STAGE_1_5_6_SPARSE_CLOCK.md` §5.9 / §5.10 / §5.11 | the census, its correction, and the contamination reach |
| *(to be written)* `src/cellfate/data/qc.py` or the Gill source | where G1/G2 live |
| *(to be written)* `tests/test_bulk_sample_integrity.py` | B1–B4 |
