# STAGE 1.5.6 — The clock's density IS the defect. Sparsify it.

**Status:** ✅ **MEASURED 2026-08-04, validated leave-one-donor-out. Not yet applied to any label.**

**Scope of what has run:** 3 new read-only scripts, **0 lines changed in `src/`**, **no label moved.**
**Scope of what this proposes:** one config change, gated, with a pre-registered bar.

---

## 1. The finding

The Fleischer clock is a dense RidgeCV over **33,155 genes fitted from 133 samples**. Restricting it
to its **~100 largest-|weight| genes** and changing nothing else:

| | MAE vs methylation | bias | ρ | sign agreement |
|---|---:|---:|---:|---:|
| **full clock (33,155 genes)** | **16.61 yr** | **−14.10** | +0.703 | 0.62 |
| **top-100 genes** | **5.36 yr** | **−1.61** | **+0.835** | **0.94** |

*68 conditions, Gill transient arm, Horvath multi-tissue as truth, ΔAge vs ΔAge, replicates averaged.*

**MAE 5.36 yr is below the reference instrument's own donor-level error of ±7 yr** (1.5.2 §12-R: two
donors of identical age 53 read 44.0 and 58.5). The RNA clock now agrees with methylation about as
closely as methylation agrees with itself.

### 1.1 The mechanism, visible in the sweep

| k | 20 | 50 | **100** | 150 | 300 | 1000 | all 33,155 |
|---|---:|---:|---:|---:|---:|---:|---:|
| MAE | 6.19 | 5.99 | **5.36** | 5.42 | 9.99 | 12.33 | **16.61** |
| bias | +5.06 | +3.22 | **−1.61** | −1.99 | −7.69 | −10.15 | **−14.10** |

**The bias crosses zero exactly where MAE bottoms out.** That is not a tuning coincidence — it is the
mechanism. Thousands of near-zero weights each contribute a little drift; summed over 33,155 genes
they become a **−14 yr systematic offset**. Dropping them removes the offset rather than merely
shrinking noise.

### 1.2 It generalises — leave-one-donor-out

| held-out donor | k chosen on the other two | held-out MAE | full clock |
|---|---:|---:|---:|
| O1 | 50 | **6.70** | 16.59 |
| O2 | 100 | **6.84** | 16.48 |
| O3 | 100 | **5.45** | 16.75 |

**k is stable at 50–100 and the ~2.5–3× improvement survives donor-level cross-validation.** The full
clock is worse than every sparse variant for **every donor individually**, so this is not selection.

### 1.3 The Sendai arm agrees, independently

22 conditions, absolute age (GSE165178 has no untreated control, so ΔAge cannot be formed there).
MAE minimum at **k = 50** — 28.52 yr against the full clock's 65.63. Different protocol, partly
different donors, same conclusion about where the optimum lies.

---

## 2. Why nine months of correlation tests never saw this

Every earlier test in this arc scored `ρ_partial`. **Spearman is shift- and scale-free, so a uniform
−14 yr offset is completely invisible to it.** M-2a, 1.5.4 and the first variant sweep all measured
ordering and all reported "weak but present". None of them could see the bias, because the statistic
was blind to it by construction.

**MAE and bias against a real instrument are what exposed it.** The lesson generalises: a
correlation is the wrong summary for a quantity whose *units* are the claim.

---

## 3. What is NOT fixed

**Horvath skin & blood does not come right at any k.** Its MAE improves (17.84 → 6.69 at k=50) but
ordering stays poor throughout (ρ ≤ 0.43, sign agreement 0.41–0.68, sometimes below chance). So the
sb/mt asymmetry recorded in 1.5.2 and 1.5.4 **survives** — sparsification fixes the *offset*, not the
disagreement between the two reference clocks.

**This does not clear M-2a's SPLIT rule.** A variant passing on one clock and not the other is still
a SPLIT. What has changed is that we now know *what kind* of failure it is: not a weak signal, a
biased one on one axis and a mis-ordered one on the other.

---

## 4. 🔑 Why this plausibly fixes HFF too

**The −14 yr bias is a property of the clock's dense weights, not of any dataset.** It is applied to
every cell the pipeline labels — including HFF's **33,613 labels, 99.8 % of the total**. Sparsifying
the clock changes all of them, not just the 90 with paired methylation.

**A concrete, falsifiable prediction:** G-c step 1 measured HFF's trajectory reaching **−24.0 yr at
day 14** under the full clock. If the −14 yr offset is the same artefact, the sparse clock should
read HFF's day-14 at roughly **−10 yr**, with the same ρ ≈ −0.9 trajectory shape.

**That is checkable today, on data already on disk, with no retrain.** It is step 1 below.

### 4.1 ⚠️ RUN 2026-08-04 — the prediction is VOID, and the reason is a bigger finding

| k | ρ(day, ΔAge) | day-14 ΔAge |
|---|---:|---:|
| top50 | −0.905 | −8.36 |
| **top100** | −0.881 | **−10.72** |
| top150 | −0.857 | −11.22 |
| **all 33,155** | −0.857 | **−10.62** |

**The shape survives sparsification** (ρ −0.857 → −0.881) — that was the falsification condition and
it passed. **But the predicted magnitude shift did not occur, because the full clock in this run
already reads −10.62, not the −24.0 the prediction was built on.**

**My script initially reported CONFIRMED. That was wrong** — it checked the result against the
prediction without checking that the baseline the prediction was predicated on still held. Verdict
logic corrected to `BASELINE_NOT_REPRODUCED`.

### 4.2 🔑 The 13.4-year gap, which is the real finding here

| day | pipeline `y_age` (built shards) | clock applied directly to counts |
|---:|---:|---:|
| 2 | **+3.85** | −0.52 |
| 4 | +3.54 | +3.58 |
| 6 | −5.80 | −0.35 |
| 12 | −8.23 | −1.26 |
| **14** | **−24.02** | **−10.62** |

G-c step 1 read **built shards** — the pipeline's `y_age`, which adds **harmonization (Gill
Projection), cell-cycle deconfounding and control re-centring**. This script applies the clock
directly and does none of that.

> **About half of HFF's apparent ΔAge magnitude is contributed by pipeline processing, not by the
> clock.** −10.62 becomes −24.02. Nobody has audited that 13.4 yr.

That reorders the plan. A sparse clock addresses a −14 yr bias *in the clock*; it cannot address a
+13.4 yr contribution from *downstream processing*. **Both are the same size, and only one of them
has been measured.**

**Note also the direct route's day-4 = +3.58** — cells reading 3.6 years OLDER four days into
reprogramming, in both routes. That is the same class of impossibility as the +36.5 yr
non-responder reading, surviving in HFF and unexplained.

---

## 5. The plan

| # | step | gate | cost |
|---|---|---|---|
| **1** | ✅ **DONE.** Shape survives (ρ −0.881), but the baseline was not reproduced — see §4.1/§4.2 | — | done |
| **1b** | 🔴 **NEW, and now the priority: audit the 13.4 yr the pipeline adds.** Decompose `y_age` into clock → harmonize → deconfound → re-centre, and report each step's contribution on HFF | each step's ΔAge contribution reported in years | free, no retrain |
| **2** | **Pre-register the bar** for adopting a sparse clock: MAE ≤ 8 yr **and** sign agreement ≥ 0.80 vs methylation, on **both** arms, k fixed at 100 in advance | `bar_verdict` row in `tests/test_bars_resolvable.py` | free |
| **3** | **Write `configs/clocks/fleischer_clock_top100.json`** — the same coefficients, 33,055 zeroed. Provenance in `meta`, original untouched | ships as a **new file**; nothing switches automatically | free |
| **4** | **One rebuild + LOOCV under the sparse clock**, full scorecard, snapshot and rollback | every Stage 1 guard reported before/after | one retrain |
| **5** | Only then decide on the label change | — | — |

**Step 1 is the falsifier and it costs nothing. Do it first.**

---

## 6. What this does not license

* **It is not yet a label change.** Nothing in `src/` or `configs/` has moved.
* **It does not resolve the sb/mt split** (§3).
* **It does not make ΔAge "correct"** — it makes it agree with methylation to within methylation's
  own error on one arm, one clock. That is a large improvement on a 16.6 yr disagreement, and it is
  not the same as validation.
* **k = 100 is chosen from a 12-point sweep.** LODO says the choice generalises, but k should be
  **frozen at 100 before step 4**, not re-tuned on the retrain.

---

## 7. Artefacts

| file | role |
|---|---|
| `experiments/diag_dage_variants.py` | step 1 of the search — 9 variants vs the FACS outcome (all null; ceiling effect at AUC 0.9221) |
| `experiments/diag_dage_variants_meth.py` | the same 9 vs methylation, both arms |
| `experiments/diag_dage_ledger.py` | the per-condition ledger: truth / expected / actual / error |
| `experiments/diag_dage_ksweep.py` | the k-sweep that located the optimum |
| `results/DAGE_LEDGER.md` | the readable record, 68 per-condition rows |
| `results/dage_ledger.csv` | full table, 90 rows × 60 columns |
| `results/diag_dage_ksweep_results.json` | every k, both clocks, both arms |
