# ΔAge ON GSE165177 — PRE-REGISTRATION. Written and committed BEFORE the run.

**Status:** 🔒 pre-registration. Committed before `experiments/dage_gse165177.py` was executed.
**Why now:** regime E established that `p_unsafe` is **not expressible in bulk at any
replication**, and explicitly left open that ΔAge — continuous per sample — does **not** have that
problem. `GSE165177` offers three things this project has never had at once, and they map directly
onto defects D1/D2/D3 of the original zero-point fix plan:

| what `GSE165177` has | the defect it addresses |
|---|---|
| **33 contemporaneous negative controls**, 2–3 per donor **per timepoint** | **D1** — every `gill_bulk` ΔAge is measured against a **day-0, cross-batch** zero-point; ~50 % of samples are `age(Exp1) − age(Exp2 baseline)` |
| **replicated** controls (n = 2–3, not n = 1) | **D2** — `gill_bulk`'s baseline is one unreplicated sample per donor, with no error bar and nothing recording that |
| donors aged **53, 53, 38 — all inside the clock's fitted `[1, 96]`** | **D3** — HFF is neonatal (age 0), so 99.7 % of our age labels are extrapolations past the clock's declared validity |

---

## 1. Normalisation — the thing most likely to silently corrupt this

`configs/clocks/fleischer_clock.json` declares `normalization: log1p_cp10k`. `GSE165177` ships
**Log2 RPM**. Applying the clock directly to Log2 RPM is a **normalisation mismatch** and would put
every age on the wrong scale (≈ 1/ln 2 ≈ 1.44× inflated, plus a pseudocount nonlinearity).

**The pipeline already solves this for the identically-formatted `gill_bulk`**
(`sources.py:506`): `rpm = 2**log2 − 1.0` → linear RPM → the pipeline's own
`normalize_counts(…, target_sum=1e4)` → `log1p`, valid because CP10k(RPM) == CP10k(counts).
**This run reuses that exact path** rather than re-deriving it. Any deviation would make the
numbers incomparable with every recorded ΔAge in the project.

---

## 2. What will be measured

`experiments/dage_gse165177.py`, **READ-ONLY** — no build, no retrain, `src/` untouched. Clock via
`LinearClock.from_json('configs/clocks/fleischer_clock.json')`, ages via `predict_age`.

Arms are grouped by name into **control** (`negative_control*`), **failed** (`fail*`),
**transient** (`transient*`), plus the day-0 fibroblasts.

---

## 3. 🔒 PRE-REGISTERED MEASUREMENTS AND BARS

### M-E0 — the C-7 integrity gate, on this dataset

Run `integrity.screen_bulk_matrix` on the raw Log2 columns — the same gate that caught the
degenerate `N2_Fib_Sendai_Exp2`, which corrupted five of six `gill_bulk` folds.

| result | meaning | action |
|---|---|---|
| **0 columns rejected** | the dataset is clean by the project's own standard | proceed; report it |
| **1+ rejected** | the same defect class is present here | **exclude them and report which**; every downstream number is computed without them |
| **> 10 % rejected** | a systemic quality problem | **stop and report**; do not read M-E1…M-E5 |

### M-E1 — does the clock read absolute age on this data?

Clock the **untreated** samples (negative controls + day-0 fibroblasts). True donor ages 53/53/38,
**mean 48.0**.

| result | meaning | action |
|---|---|---|
| **PASS-CALIBRATION** — \|mean predicted − 48.0\| ≤ **12.27 yr** (one `cv_mae`) | the clock reads this dataset in the right absolute range | absolute ages are usable, with the error bar stated |
| **FAIL-CALIBRATION** — otherwise | the clock does not read absolute age here | **only control-relative ΔAge may be used**; every absolute age is reported as invalid |

> **The 53-vs-38 contrast is NOT gated and must not be claimed either way.** The gap is 15 yr
> against a 12.27 yr `cv_mae` — about 1.2×, badly underpowered. It is reported as *indicative
> only*, exactly as M1's middle contrast was in the Stage 1.5 audit.

### M-E2 — the contemporaneous-control ΔAge (what `gill_bulk` could never compute)

For every treated sample: `ΔAge = clock(sample) − mean(clock(negative controls of the SAME donor
at the SAME day))`. Descriptive; reported per (donor, day, arm). **No bar** — this is the
quantity, not a test.

### M-E3 — does transient reprogramming rejuvenate? (external validation against Gill 2022)

Paired across (donor, day) cells that carry both arms, on the contemporaneous-control ΔAge:

| result | meaning | action |
|---|---|---|
| **REPRODUCED** — mean ΔAge(transient) < 0 and the 95 % CI **excludes 0** | the published direction is recovered by our clock **with proper contemporaneous controls** | a genuine external validation; record it as such |
| **NOT REPRODUCED** — CI includes 0 | our pipeline does not recover it at this n | record; **do not** claim the effect is absent — state the power |
| **CONTRADICTED** — CI excludes 0 on the **positive** side | our pipeline reads transient reprogramming as **ageing** | **escalate**; this would question the ΔAge target itself |

### M-E4 — how noisy is the zero-point really? (the question D2 could never answer)

SD of the clock reading **within** each (donor, day) control group (n = 2–3), pooled.

| result | meaning | action |
|---|---|---|
| **pooled control SD ≥ 12.27 yr** (`cv_mae`) | the zero-point wobbles as much as the clock's own error | the ±12.7 yr per-donor offset is **indistinguishable from measurement noise** — Stage 2's premise stays void as stated |
| **pooled control SD < ½ × 12.27 yr** | the clock is far more reproducible on replicates than its CV implies | the per-donor offset is **more likely real biology**; Stage 2's premise is strengthened and this is the first evidence for it |
| in between | inconclusive at this n; report the number and the CI | |

### M-E5 — the `exp1` / `exp2` batch offset (D1, measured directly)

At matched (donor, day, arm) present in both batches, the mean clock difference.

| result | meaning | action |
|---|---|---|
| \|offset\| ≥ 12.27 yr | batch alone moves ΔAge by a full clock error | **D1 is confirmed as severe**; any cross-batch ΔAge is unusable |
| \|offset\| < 12.27 yr | batch is present but sub-error | record the number; treat as a stated caveat |
| **no matched cells exist** | the design does not permit the comparison | say so plainly rather than substituting a weaker contrast |

---

## 4. Declared limits — stated before any number is seen

1. **Bulk.** Every sample is a population average. This is fine for ΔAge (continuous per sample)
   and was fatal for `p_unsafe` — that asymmetry is the entire reason this run exists.
2. **3 donors, 2 distinct ages (53, 53, 38).** Any age-calibration claim is weak by construction.
3. **No harmonizer.** ΔAge here is raw control-relative, not σ-harmonised, so it is **not**
   numerically comparable to the project's harmonised ΔAge figures. Direction and magnitude only.
4. **Arms are grouped by name.** `failing_to_transiently_reprogram_intermediate` and
   `failed_to_transiently_reprogram` are pooled as *failed*; the grouping is reported so it can be
   disputed.
5. **This does not re-open Stage 1.5.** It measures a dataset that was never in a training config.

---

## 5. Recording

`results/dage_gse165177_results.json`; write-up to `CHANGES.md` and
`experiments/DELTAAGE_LAB_NOTEBOOK.md`; unit tests in `tests/test_dage_gse165177.py`. Every bar in
§3 graded as written, including any that fail.
