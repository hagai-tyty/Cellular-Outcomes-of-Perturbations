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

---

## 6. RESULT — 2026-08-12. **Gill 2022's rejuvenation claim REPRODUCED**; the clock's absolute age FAILS; and two of my own bars were badly posed.

*Graded against §3 exactly as written. Artefacts: `experiments/dage_gse165177.py`,
`results/dage_gse165177_results.json`, `tests/test_dage_gse165177.py` (18 tests).*
**READ-ONLY** — no build, no retrain, `src/` untouched. 93 samples × 35,720 genes; groups
control 33, failed 33, transient 24, day-0 3.

### M-E0 — the C-7 integrity gate: ✅ **CLEAN**

**0 of 93 columns rejected.** `GSE165177` passes the gate that rejected five `gill_bulk` columns
and whose single worst case (`N2_Fib_Sendai_Exp2`) corrupted five of six folds. By the project's
own standard this dataset is cleaner than the one every ΔAge to date was computed on.

### M-E1 — absolute age: ❌ **FAIL-CALIBRATION**

| donor | true age | n untreated | predicted | 95 % CI | error |
|---|---|---|---|---|---|
| O1 | 53 | 12 | 89.2 | [85.0, 93.5] | **+36.2** |
| O2 | 53 | 12 | 97.3 | [93.6, 100.9] | **+44.3** |
| O3 | 38 | 12 | 95.9 | [92.4, 99.4] | **+57.9** |

Pooled **94.1 against a true mean of 48.0 → |Δ| = 46.1 yr**, nearly **4× one `cv_mae`**. Every
untreated fibroblast reads 65–105 yr. **Per the pre-registration, absolute ages on this dataset are
invalid and only control-relative ΔAge is used below.**

*Gene coverage does not excuse it.* 18,928 of 33,155 clock genes are present (57.1 %) — but that
is **89.2 % of the clock's total |weight|**, so the missing genes are overwhelmingly low-weight.
Coverage is a contributing bias, not an explanation for +46 yr.

*The 53-vs-38 contrast was pre-registered as NOT gated and is reported as indicative only —* and
the indication is unfavourable: **O3, the youngest donor at 38, reads the second-highest.** At
15 yr against a 12.27 yr `cv_mae` this resolves nothing, and no claim is made either way.

### M-E2 — contemporaneous-control ΔAge: computed, for the first time in this project

Every value is `clock(sample) − mean(clock(controls of the SAME donor at the SAME day))`. `failed`
runs −6.2 to −11.0 yr; `transient` runs −3.1 to −24.3 yr. **This is the quantity `gill_bulk` is
structurally incapable of producing** — its only baseline is one day-0 sample per donor, from a
different batch for ~half the data.

### M-E3 — does transient reprogramming rejuvenate? ✅ **REPRODUCED**

| | mean | 95 % CI | n |
|---|---|---|---|
| ΔAge(transient) vs its own contemporaneous control | **−17.88 yr** | **[−21.13, −14.64]** | 12 |
| paired transient − failed | **−9.58 yr** | **[−12.77, −6.39]** | 12 |

Both CIs exclude zero, and **11 of 12 (donor, day) cells are negative** (the sole exception is
O3 day 15 at +3.29). Gill 2022's central claim is recovered by this project's own clock, on data
that was in no training config, against **replicated contemporaneous controls**.

> **Why M-E3 can succeed while M-E1 fails, and this is the load-bearing point:** ΔAge is a
> *difference*. The clock's +46 yr bias, and every missing-gene term, appear in the treated sample
> and in its control alike and **cancel exactly**. Absolute age needs the clock to be *accurate*;
> ΔAge needs it only to be *consistent*. This is the first direct evidence in the project that the
> control-relative design does the job it was chosen for.

### M-E4 — the zero-point's own noise: bar fired, ⚠️ **but the inference I attached to it does not follow**

Pooled within-(donor, day) control SD over 12 groups of n = 2–3: **5.04 yr**, against
`cv_mae` 12.27. By §3 as written this is the *"< ½ × cv_mae"* branch, whose pre-registered reading
was *"the per-donor offset is MORE LIKELY REAL BIOLOGY; first evidence for Stage 2's premise."*

**That reading is wrong, and this same run disproves it.** The two quantities are not comparable:
`cv_mae` is cross-validated error against **true age across 133 donors**; the 5.04 yr is
**replicate scatter within one condition**. A clock can be highly reproducible and badly biased —
and M-E1 shows *exactly that*: 5 yr replicate scatter alongside a **+46 yr** bias. **Reproducibility
is not accuracy**, so low replicate SD says nothing about whether the per-donor offset is biology.

**Recorded as: the bar fired as written; the conclusion attached to it is withdrawn.** Stage 2's
premise is **not** strengthened by this. The honest statement is narrower and still useful: *the
zero-point is reproducible to ~5 yr when controls are contemporaneous and replicated*, which is a
real improvement over an unreplicated single sample and is worth having on its own.

### M-E5 — the batch offset: pooled bar says "sub-error", ⚠️ **but the pooled statistic was the wrong one**

Pooled `exp1 − exp2` = **−2.91 yr**, CI [−6.35, +0.54] → by §3, *"batch present but sub-error"*.
**Stratifying by arm — not pre-registered — shows that summary is misleading:**

| arm | mean `exp1 − exp2` | 95 % CI | n |
|---|---|---|---|
| **control** | **−8.52 yr** | **[−10.13, −6.92]** | 9 |
| failed | −4.29 yr | [−7.52, −1.06] | 9 |
| transient | +7.59 yr | [−3.60, +18.77] | 6 |

The control and transient arms move in **opposite directions**, so the pooled mean cancels an
effect that is real in both. The control CI excludes zero: **the ZERO-POINT itself shifts −8.52 yr
between batches — 0.69 `cv_mae`.** Because ΔAge is measured *against* that control, an
arm-dependent offset **does not cancel**, which is precisely D1's concern. **D1 is confirmed as
material for any cross-batch comparison**, and my pre-registered pooled bar could not have seen it.

### Summary of what this establishes

| | |
|---|---|
| `GSE165177` for **ΔAge** | ✅ **valuable, and now demonstrated** — clean by C-7, contemporaneous replicated controls, and it reproduces the published effect |
| Gill 2022's rejuvenation direction | ✅ **reproduced** by our clock, −17.88 yr [−21.13, −14.64] |
| the control-relative ΔAge design | ✅ **vindicated** — it works precisely where absolute age fails |
| the clock's **absolute** age on bulk fibroblasts | ❌ **unusable** — +46 yr, ~4 `cv_mae`; **not** explained by gene coverage (89.2 % of weight mass present) |
| D1 (cross-batch zero-point) | ✅ **confirmed material** — the control arm alone shifts −8.52 yr [−10.13, −6.92] |
| D2 (unreplicated baseline) | **quantified at last** — replicate SD 5.04 yr |
| Stage 2's premise (the ±12.7 yr offset is biology) | **NOT strengthened.** My M-E4 bar conflated reproducibility with accuracy; withdrawn above |

### What is NOT claimed

That the clock is calibrated on this data (it is not). That the 53-vs-38 contrast resolves (it does
not, and it points the wrong way). That these ΔAge values are comparable to the project's
harmonised figures — **no harmonizer was applied**, by design, so they are directional only. And
nothing here re-opens Stage 1.5 or touches `p_unsafe`, which regime E settled separately.
