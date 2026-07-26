# STAGE 1.5.1 — NEW CHANGES (review findings + revised execution plan)

**Companion to** `STAGE_1_5_1_CLOCK_PRECISION.md`, which is **unmodified** and remains the
authoritative pre-registration. Everything here is additive: findings from an independent review,
the tests that settled them, and the plan changes those tests imply.

**Status:** review EXECUTED (2026-07-26, on the real GSE113957). Steps below are PLANNED.

---

## 0. ⚠️ THE PROBLEM — what goes wrong if you follow `STAGE_1_5_1_CLOCK_PRECISION.md` literally

**Read this before executing that document.** Its bars, discipline, revalidation suite and fallback
are all correct and unchanged (§1). But three things in its *execution steps* are now known to be
wrong, and one of them is not merely wasteful:

| # | Problem in the plan as written | Where | Severity |
|---|---|---|---|
| **P1** | **Step 1 budgets an audit for a leak that does not exist.** R4 claims the CV is optimistic because data is "standardised before the split". There is **no cross-sample scaler in `clock_fit.py` at all** — proven three ways (§2). Re-measuring "leak-free" reproduces 12.27 and finds nothing. | §1 R4, §2 Step 1 | ⚠️ **Wasted time.** Harmless to correctness |
| **P2** | **Step 2 lists a candidate already eliminated.** C3 (slope recalibration) was measured at **12.78 yr — 1% *worse* than the control** (§4). Compression is a *symptom* of R1, not a patchable defect. | §2 Step 2 (C3) | ⚠️ **Wasted time.** Harmless to correctness |
| **P3** | **No step pins the training sample set — and it is currently unknown.** The artefact records `n_samples = 133`; GSE113957 contains **143**. **Which ten were excluded is recorded nowhere** (§5). | missing entirely | 🔴 **CORRECTNESS.** See below |

### Why P3 is the one that actually matters

Every bar in §3 is a comparison against `cv_mae = 12.27`. If the refit trains on a different sample
set than the original did, then "12.27 → X" is **not a like-for-like comparison**, and a PASS could
come from the 7.5% extra data rather than from the new estimator. The stage would then have proven
nothing about C1/C2 while appearing to succeed — the exact failure mode this project's discipline
exists to prevent.

It also propagates: whatever clock ships inherits the same undocumented sample set, so Stage 5
cannot state its training data — a question a reviewer *will* ask.

**Fix:** run **Step 0** (§6) before anything else — identify the original 133, or explicitly define
and record the set the refit will use. Minutes of work; without it the rest is uninterpretable.

### What to do

**Follow `STAGE_1_5_1_CLOCK_PRECISION.md` for the bars, revalidation, fallback and discipline —
those stand.** Take the execution steps from **§6 of this file**, which drops P1 and P2 and adds
Step 0 for P3.

> **One caution that cuts against my own finding.** Refuting R4 does **not** mean the leakage worry
> was baseless. The original's §3 guard — *"any feature selection inside each fold"* — is
> **correct and load-bearing for C1/C2**: unlike per-row normalisation, selecting genes by
> age-correlation **is** a cross-sample statistic and **will** leak if done outside the fold. R4
> named the wrong culprit; the guard it motivated must still be honoured.

---

## 1. Verdict on the original plan: **endorsed**

Checked against the tree, not accepted on report.

| Claim in `STAGE_1_5_1_CLOCK_PRECISION.md` | Check | Result |
|---|---|---|
| dense ridge, 33,155 genes, no selection | read the artefact | ✅ all 33,155 weights non-zero |
| `cv_mae` 12.27, `cv_pearson` 0.837, range [1,96] | read the artefact | ✅ exact |
| primary bar arithmetic (`√2·cv_mae ≤ 5.5` ⇒ `cv_mae ≤ 3.9`) | recomputed | ✅ correct |
| **R1** dense ridge at p/n≈250 overfits | new evidence, §3 | ✅ **strongly corroborated** |
| **R2** shrinkage compression | measured, §3 | ✅ **confirmed and quantified** |

**The SNR≈1 diagnosis is right, the bars are sound, and this is the correct next stage.** Nothing
below changes that. The changes are: one hypothesis is refuted, two Step-1 items are already done,
and one candidate is eliminated.

---

## 2. ❌ R4 is refuted — tested three ways, not argued

> **R4 as written:** *"`cross_val_predict` is run on `Xn`, which is standardised **before** the
> split. If the scaler sees all samples, `cv_mae = 12.27` is itself an optimistic estimate."*

### Test 1 — does any cross-sample statistic exist?

```
normalize_counts(X[:5])   vs   normalize_counts(X)[:5]      ->   max |difference| = 0.0
```

A sample's normalised value is **identical** whether or not the other 138 samples are present,
because of one line in `normalize.py`:

```python
lib = counts.sum(axis=1, keepdims=True)   # axis=1 = across GENES, within one sample
```

Each row is divided by **its own** total; no column statistic is ever computed. A split cannot leak
through a function that never looks across samples. `grep StandardScaler|scaler|fit_transform` over
`clock_fit.py` → **no match**. **The scaler R4 describes does not exist.**

### Test 2 — is there a *different* leak R4 might have been sensing?

The plausible alternative is **group leakage**: if one donor contributed several samples,
`KFold(shuffle=True)` would split them across folds and the CV *would* be optimistic — same
symptom, nothing to do with scalers. Measured on the real series matrices:

```
143 samples   ->   143 unique `cell id`s   ->   0 donors appearing twice
```

Every sample is a distinct donor. **Not that either.**

### Test 3 — reproduce the CV directly (decisive)

Exact protocol (`KFold(5, shuffle, seed=0)`, `RidgeCV` refit inside each fold, per-row
normalisation), run on the real data:

| | replicated | shipped artefact |
|---|---|---|
| `cv_mae` | **12.67 yr** | 12.27 yr |
| `cv_pearson` | **0.841** | 0.837 |

It reproduces. (Gap expected: 143 samples here vs the artefact's 133, and a newer NCBI annotation.)

### Conclusion, and why the error was a reasonable one

**`cv_mae = 12.27` is an honest, leak-free number. R4 is wrong in mechanism *and* conclusion.**

Scaler-before-split is *the* classic CV leak, and the code does read `Xn = normalize_counts(X)`
directly above `cross_val_predict`. But **normalisation ≠ standardisation**:

| | computes | leaks? |
|---|---|---|
| library-size **normalisation** (what this code does) | one row's own total | **no** |
| per-gene **standardisation** | a column mean/SD *across samples* | **yes** |

They look nearly identical in code and behave completely differently. Right instinct, wrong
identification — and the original hedged it correctly ("possible", "to be audited, not assumed").

> **This does not weaken the stage — it strengthens it.** R4 was the one route by which the SNR
> problem could have been *overstated*. It is closed. The ruler really is ±12.3 yr.

---

## 3. What the same run established (Step 1, effectively complete)

Total compute: ~3 seconds.

| Quantity | Measured | Reading |
|---|---|---|
| predicted-vs-true **slope** on held-out folds | **0.717** | **R2 confirmed.** Bar S1 wants [0.85, 1.15] — real compression, not a rounding artefact |
| `alpha` chosen by `RidgeCV` | **0.272** | near the **bottom** of `logspace(-1, 4)` → penalty barely binding → **R1 confirmed** |
| in-sample MAE vs CV MAE | **0.77 vs 12.67 = 16×** | **memorisation**, not mild over-regularisation |

The 16× gap is the most direct evidence in the record for R1 and was not available when the
original was written (it came from §9's reproduction check).

### Step 1 status

| Original Step 1 item | Status |
|---|---|
| re-measure `cv_mae` leak-free (tests R4) | ✅ **done — no leak exists; 12.27 is honest** |
| predicted-vs-true slope (tests R2) | ✅ **done — 0.717, compression confirmed** |
| error profile by age decile (feeds S4) | ⏳ **still to run** |

---

## 4. ❌ C3 eliminated — tested, not predicted

The original lists **C3 (ridge + slope recalibration)** as "the cheapest possible fix for R2 alone".
Tested directly, with the recalibration fitted **out-of-fold** so it cannot cheat:

| candidate | `cv_mae` | vs §3 bar |
|---|---|---|
| C4 dense ridge (control) | 12.67 yr | FAIL |
| **C3 + slope recalibration** | **12.78 yr (−1%, i.e. slightly worse)** | **FAIL** |

Correcting the slope of a model that is *memorising* rescales its errors without removing them. The
compression is a **symptom** of R1, not an independent defect that can be patched downstream.

**C3 is answered and should not consume Step 2 time.** It remains useful only as a reported control.

---

## 5. ⚠️ Loose end found: `n_samples` 133 vs 143

The artefact records `n_samples = 133`; GSE113957 contains **143** samples (13 on GPL16791 + 130 on
GPL18573). **Ten samples were excluded when the clock was fit, and which ten is not recorded
anywhere.**

Consequences:
- The §9 reproduction (MAE 0.77 yr) covered **133 of 143 in-sample**, plus ~10 that were genuinely
  held out — so that number is *slightly* less purely in-sample than stated, without changing the
  conclusion.
- Any refit must **state its sample set explicitly**, or the new clock inherits the same ambiguity.
- If the ten were dropped for a quality reason, that reason should be re-applied deliberately; if
  they were dropped by an accident of parsing, they are 7.5% more training data for a p≫n problem.

**Action:** Step 0 below.

---

## 6. Revised execution plan

Changes only. Everything not mentioned stands as written in `STAGE_1_5_1_CLOCK_PRECISION.md`.

### Step 0 (new, minutes) — pin the sample set
Determine which 143→133 samples the current clock used, or, failing that, **define and record** the
sample set the refit will use. No refit is interpretable against an unknown training set.

### Step 1 (reduced) — only the age-decile profile remains
R4 and the slope are done (§2, §3). Run the error-by-decile profile to feed bar **S4**.

### Step 2 (narrowed) — C1 and C2 versus C4
Drop C3 to a reported control (§4). Evaluate under one leak-free harness:
- **C1** ElasticNet / Lasso — sparse selection, targets R1 + R3.
- **C2** Ridge on a pre-filtered gene set, selection **inside** each fold.
- **C4** current dense ridge — control.

> **The in-fold guard still matters — for C1/C2, not for R4.** The original's §3 guard ("any feature
> selection inside each fold") is *correct and load-bearing*: unlike per-row normalisation, gene
> selection by age-correlation **is** a cross-sample statistic and **will** leak if done outside the
> fold. R4 was the wrong target; the guard it motivated is right.

### Steps 3–8 — unchanged
Bars (§3), revalidation (§5), fallback (§6), discipline (§7), deliverables (§8) all stand.

---

## 7. Honest expectation, recorded before Step 2

`cv_mae` must fall **12.27 → ≤ 4.0** — a **3× improvement** — from n=133 samples. Sparse methods
address p≫n but cannot create information that is not in 133 samples.

| outcome | my estimate |
|---|---|
| **PASS** (≤ 4.0 yr) | ~25–35% |
| **MARGINAL** (4.0–6.0 yr) | ~40% — the most likely single outcome |
| **FAIL** (> 6.0 yr) | ~30% |

Recorded now so a MARGINAL result is read as the pre-registered outcome it is (ranking-only
framing, §6.3 of the original) rather than argued about afterwards. **The fallback is a real
finding, not a defeat** — "transcriptomic age on n=133 fibroblasts cannot support per-cell
rejuvenation claims" is publishable and directly shapes Stage 5.

---

## 8. Also withdrawn here: the OOD-detector proposal

Proposed in discussion before §10's D2 result; recorded rather than quietly dropped. The idea was to
gate ΔAge on the model's existing OOD detector, treating reprogramming intermediates as
out-of-domain. **Three independent reasons it fails:**

1. **It measures the wrong distribution.** The detector is a Gaussian over the *model's* latent `z`
   fitted on `train_ds` (`train_model.py:291`) — and the model's training set **contains** the Gill
   reprogramming intermediates. They are in-distribution by construction; it would flag nothing.
2. **It is already known to be uninformative** — measured OOD AUC ≈ **0.47** (chance), with
   `train_model.py:288-290` recording "disable the gate" as the anticipated outcome.
3. **§10 removed the premise.** D2 found the trajectory sign **flips** between datasets (Gill
   +0.205, GSE242423 −0.214, each clearing its bar by hundredths). A stable domain effect gives the
   *same* sign in both; this is noise at SNR≈1 — which is precisely this stage's thesis.

Even absent all three, gating would flag the reprogramming cells — the entire use case — making it
option C (retreat) in a domain-condition costume. **The precision fix is the correct target.**
