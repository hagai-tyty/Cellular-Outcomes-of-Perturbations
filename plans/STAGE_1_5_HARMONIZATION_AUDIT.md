# STAGE 1.5 — Harmonization & ΔAge zero-point audit (measurement only)

**Implements:** nothing new — it *validates a claim four existing plan docs already make*.
**Depends on:** Stage 1 scored and closed (so no snapshot mixes two changes).
**Blocking for:** Stage 2 — its premise (the per-donor offset is real biology) is what this audit
tests. Not blocking for the tool.
**Scope:** 1 new test file, 1 new verifier script, **0 lines changed in `src/`**.

**Status:** ✅ **EXECUTED and INDEPENDENTLY VERIFIED** (2026-07-24). The audit passed on its own
terms; it surfaced a second, still-open artefact (§5.2) and a fix plan (§5.4, tightened in §6.2).

### Status ledger — what is done, what is open

*Every ✅ below was re-checked against the tree by an independent review (§6.1), not accepted on
report. The pre-registration text in §0–§4 is left **exactly as written**; only this ledger and
§5–§6 record what happened.*

| Item | Status |
|---|---|
| §2 **Group A** — intercept cancellation, additive-batch immunity, `mu_ref` drops out | ✅ pass (with one correction: cancellation is numerical ~1e-14, **not** bit-identical as §2 specified) |
| §2 **Group B** — the true scale-gain invariant, closed form | ✅ pass |
| §2 **Group C** — fit / leak-safety, variance floor, gene space, raises, `_align`, round-trip | ✅ pass |
| §2 **Group D** — ΔAge zero-point incl. the silent fallback pinned | ✅ pass |
| §2 **Group E** — real-build replay | ✅ **PASS — 51/51 chunks carry ≥1 control; the fallback never fired.** All six LOOCV donors covered, so the PASS is not vacuous |
| §2 `verify_stage1_5.py` with pure `decide_verdict()` | ✅ built, every branch unit-tested |
| §3 Groups A–D pass with no repo data | ✅ 21/21 |
| §3 full suite green | ✅ **303 passed** |
| §3 `git diff --stat src/` empty | ✅ empty across all 5 commits |
| §3 Group E reports whether the fallback fired | ✅ reported: it did not |
| **Tests have real detection power** | ✅ mutation-tested — 4 injected defects, 4 caught (§6.1) |
| §0.1 "batch-immune by construction" is an overstatement | ✅ **confirmed false as written**; exact invariant now pinned |
| §0.3 is the ±12.7 yr offset the *fallback* artefact? | ✅ **answered: no** |
| Is the offset *biology*? | ⏳ **still open** — §5.2 found a second candidate artefact (`n=1` baseline) |
| §5.4 **Phase 1** — M1/M2/M3 measurements | ✅ **EXECUTED — M1 FAILED. ACTION: ESCALATE** (§7). M2's verdict was a stub and is being fixed; M3 indeterminate as predicted |
| §5.4 Phases 2–4 | ⛔ **BLOCKED by the M1 failure** — do not proceed until the clock's validity is settled |
| Does the clock read chronological age on this data? | ❌ **NO at the extremes** — 11.8 yr contrast vs a 53 yr true gap; two age-0 donors read 62 yr apart (§7) |

---

## 0. Why this stage exists

Four plan documents assert cross-modality harmonization is already validated:

| Claim | Where | Kind of claim |
|---|---|---|
| "parameter-free, leak-safe, **unit-tested**" | `MASTER_PLAN.md:48`, `REF_ARCHITECTURE.md:20` | capability table |
| "unit-tested; **intercept cancellation proven**" | `STAGE_5_PUBLICATION.md:127` | **a claim made to a reviewer** |
| "**intercept-cancellation unit test still passes**" | `STAGE_6_NEW_DATA.md:143` | **an acceptance gate** |

**No test imports or exercises `src/cellfate/data/harmonize.py`.** Zero hits across all test
files. So the Stage 6 gate names a test that does not exist — it can never fail — and the Stage 5
row promises a reviewer a proof that was never written. The intent of this stage is to make the
existing claim **true**, not to weaken it: the plan says "unit-tested", so we write the unit test.

Reading the module against the claim surfaced two concrete, testable facts.

### 0.1 "Batch-immune by construction" is only half-true

For two cells in the same dataset `d`, harmonization Z-scores against `d`'s control stats and then
projects back through the reference dataset's scale ([harmonize.py:118-132](src/cellfate/data/harmonize.py)),
before ΔAge subtracts the per-line control baseline ([aging.py:144](src/cellfate/data/aging.py)):

```
ΔAge = Σ_g (x_pert,g − x_ctrl,g) · sigma_ref,g / (sigma_d,g + EPS) · w_g
```

`mu_d`, `mu_ref`, and the clock intercept **all cancel** — that part of the claim is real and
worth proving. But `sigma_d` does **not** cancel: it survives as a per-dataset multiplicative
**gain** `sigma_ref / (sigma_d + EPS)`. HFF cells carry gain `sigma_gill / sigma_hff`; Gill cells
carry ≈1 because `gill_bulk` is the reference. So ΔAge is immune to **additive** batch effects but
carries a **scale** factor by design. That is fine as an architecture; the docstring's unqualified
"batch-immune by construction" overstates it and would not survive review. This stage replaces the
overstatement with the exact invariant.

### 0.2 A silent zero-point switch in ΔAge

ΔAge is control-relative: `ΔAge = age − mean(age over that line's vehicle controls)`. But the
baseline has a silent fallback ([aging.py:88](src/cellfate/data/aging.py)):

```python
ref = values[ctrl] if ctrl.any() else values[in_line]
```

ΔAge is computed **per chunk** ([build_dataset.py:306](src/cellfate/data/build_dataset.py)), and
chunks come from `plan_all` per source. If a donor lands in a chunk with **no vehicle controls**,
its zero-point silently flips from *control-relative* to *self-centred* — and self-centring
subtracts that donor's own mean perturbation effect, dragging its mean ΔAge toward 0. No warning,
no counter, no mask records that it happened.

### 0.3 Why this must run before Stage 2, and why it is 1.5 not 0

Stage 2 spends wet-lab resources (k≈3 reference cells per donor) on the premise that the
**±12.7 yr per-donor offset is real donor-response biology**. `delta_age` baselines per
`cell_line`, so an *additive* donor baseline cancels by design. An offset that survives is
therefore *either* that biology *or* the fallback in 0.2 firing — and **nothing currently
distinguishes the two.** Settling it is cheap and belongs before the money is spent.

It is **1.5, not 0**, because it should run *after* Stage 1's measurement is closed (so no snapshot
mixes two changes) but *before* Stage 2. Placement is safe for Stage 1: runs 1–3 all use the same
shards, so Stage 1's guards are a relative comparison throughout, and a finding here would qualify
how ΔAge is *interpreted* without re-running anything.

---

## 1. What it checks (and the discipline)

**Measurement only. Zero lines change in `src/`.** The tests pin *current* behaviour. If the
real-data check shows the fallback fired, that becomes its **own** pre-registered Change with its
own snapshot — silently "fixing" the ΔAge zero-point mid-audit would move every target and
invalidate the guards, which is the "be careful not to hurt any biology" line. `git diff --stat
src/` must be empty when this stage finishes.

---

## 2. Tests needed for the audit to be clear

`tests/test_harmonize.py` (new), organised so each group answers one question and a failure points
at one cause.

### Group A — the proof the plans already promise *(synthetic, no repo data, CPU-seconds)*
The intercept-cancellation proof `STAGE_5` and `STAGE_6` name. On a synthetic 2-dataset fixture:
- **intercept cancellation** — perturbing `LinearClock.intercept` leaves ΔAge **bit-identical**.
- **additive batch immunity** — add a per-gene offset to *all* of one dataset's cells, refit the
  harmonizer, and ΔAge is unchanged to tolerance.
- **`mu_ref` drops out** — the reference mean cancels in any control-relative difference through
  `project_to_clock`.

### Group B — the *true* scope (makes the paper claim precise)
- **scale is a gain, not immune** — scaling one dataset per-gene by `c` changes ΔAge by exactly
  the predicted factor `sigma_ref / (sigma_d + EPS)`. This asserts the real invariant from 0.1 in
  place of the overstated one, so the manuscript sentence can be corrected to something provable.

### Group C — fit / leak-safety (`Harmonizer.fit`)
- cells **not** passed in `controls` cannot move `mu`/`sigma` — the held-out-donor guarantee.
- variance floor holds: `sigma >= median(sigma)`.
- the common gene space is the **sorted intersection** of per-dataset admissible sets.
- `MIN_REPLICATES` violation and an unknown `dataset_id` **both raise** (not silently degrade).
- `_align` places permuted / missing genes in the correct columns.
- `to_json` / `from_json` round-trips a harmonizer unchanged.

### Group D — the ΔAge zero-point (`aging.py`)
- `_control_baseline` is genuinely **per-line**: two lines with different baselines both land at
  ~0 after subtraction.
- **the fallback, pinned explicitly** — a line with *no* controls is self-centred, so its mean
  ΔAge is forced to 0. Asserting this makes the behaviour *visible in a test* instead of silent,
  so any future change to it is a deliberate, reviewed act.
- `recenter_on_controls` restores the control-zero after `deconfound_age` re-centres the population.

### Group E — real-data diagnostic *(needs the data machine; `pytest.mark.skipif` when absent)*
- replay `plan_all(sources)` and assert that **every `(chunk, cell_line)` group containing
  non-control cells also contains ≥1 control** — i.e. the 0.2 fallback **never actually fired** on
  the real build. This is the test that converts "±12.7 yr is real biology" from an assumption
  into a checked fact.

### `verify_stage1_5.py` (new)
A runnable gate mirroring [verify_1a.py](verify_1a.py): a **pure `decide_verdict()`** separated
from I/O so every branch is unit-testable (the `verify_1a` lesson — a decision function whose only
exercised path is the one that says PASS is not a gate). Prints a PASS/FAIL table. **Reuses**
`plan_all` ([chunking.py:30](src/cellfate/data/chunking.py)) and `_control_baseline`
([aging.py:81](src/cellfate/data/aging.py)) rather than reimplementing either.

---

## 3. Acceptance

| Check | Bar |
|---|---|
| Group A–D | pass on synthetic data with **no repo data present** (CPU, seconds) |
| Full suite | stays green at its current count |
| `git diff --stat src/` | **empty** — the guarantee that no guard can move |
| Group E / `verify_stage1_5.py` | reports whether the no-control fallback fired on the real build |

A **Group E failure is a finding, not a bug to patch here**: it means some donor's ΔAge zero-point
was self-centred, so part of the ±12.7 yr offset is an artefact rather than biology — recorded in
`STAGE_1_DEVIATIONS.md`, pre-registered as its own Change, and fixed under its own snapshot. A
Group A–D failure means harmonization does not do what four plan docs claim, which is a
publication-blocking finding in its own right.

---

## 4. Documentation (appends only; existing plans untouched)

- `CHANGES.md` — per the standing changelog rule.
- `experiments/DELTAAGE_LAB_NOTEBOOK.md` — pre-register predictions **before** running.
- `STAGE_1_DEVIATIONS.md` — **only if** a test fails, recording plan-claim vs measured.

**Not touched:** `MASTER_PLAN.md`, `REF_ARCHITECTURE.md`, `STAGE_1..6_*.md`. If Group A–D pass,
the overstated "batch-immune by construction" wording in `harmonize.py` and the two reviewer-facing
rows (`STAGE_5:127`, `STAGE_6:143`) should be *corrected to the exact invariant from 0.1* — but
that is a wording change proposed to the user, not made unilaterally by this stage.

---
---

# 5. RESULTS, AND THE FIX PLAN THEY PRODUCED

> **Everything above this line is the original pre-registration, written before anything ran, and is
> left exactly as written.** Everything below was added **after** execution on 2026-07-24. It is
> appended, never substituted, so the plan and what actually happened stay auditable side by side.
> **Status: the plan in §5.4 is PLAN ONLY — none of it has been executed.**

## 5.1 What was run, and what it returned

| Part | Result |
|---|---|
| Groups A–D (synthetic) | ✅ **21/21 pass**; full suite **303**; ruff clean; `git diff --stat src/` **empty** |
| Group E (`verify_stage1_5.py`, data machine) | ✅ **PASS — 51/51 chunks carry ≥1 control.** The `aging.py:88` fallback **never fired** |

Per-chunk census: GSE242423 HFF = 45 stratified batches, **111–112 controls** of ~980 cells each
(`_batch_indices` works exactly as its docstring claims); Gill = 6 donor chunks, **exactly 1
control** of 19–21 cells each. The pre-registered pre-QC caveat did not bite — Gill cell counts are
unchanged post-QC.

**§3 acceptance: all four checks met. Stage 1.5 passes on its own terms.**

Two claims the tests corrected (both measured, not argued):
1. **"batch-immune by construction"** (`harmonize.py:9`) is **false as written**. ΔAge is immune to
   *additive* batch effects but carries a per-dataset multiplicative gain —
   `ΔAge = Σ_g δ_g · sigma_ref,g / (sigma_d,g + EPS) · w_g`, now pinned as a closed form.
2. **"intercept cancellation is bit-identical"** (§2 Group A) is **not exact**. The cancellation is
   *numerical*, not symbolic; immune to ~1e-12, and `np.array_equal` fails.

## 5.2 The finding the gate was not asking about

Every Gill donor's zero-point rests on **exactly one unreplicated control sample**. That is a PASS
by the pre-registered rule — a control exists — but `n = 1` has zero degrees of freedom, so any
error in that single day-0 measurement propagates **1:1 into every ΔAge for that donor**, landing as
a per-donor *additive offset*: structurally the same shape as the effect Stage 2 is premised on.

**Scope of the PASS.** Ruled out: the offset is an artefact of the *self-centring fallback*. **Not**
ruled out: the offset is noise in a single unreplicated baseline. Not addressed at all: whether the
offset is biology.

## 5.3 Following that into the Gill metadata — three findings

**The number that makes this urgent.** The clock carries its own cross-validated error in its
metadata (`configs/clocks/fleischer_clock.json`, Fleischer 2018 / GSE113957, 133 samples):
**`cv_mae_years = 12.27`**. The per-donor offset Stage 2 exists to correct is **±12.7 yr (ridge) /
13.12 yr (model mean |shift|)**.

> The offset attributed to donor biology is the magnitude of **one** clock measurement's error — and
> every donor's zero-point **is** exactly one clock measurement. That does not prove it is noise; it
> proves the two are **currently indistinguishable**.

**D1 — the zero-point is CROSS-BATCH (a real defect).** All six baselines are
`*_Fib_Sendai_`**`Exp2`**, while treatment samples span both experiments:

| donor | treatment n | Exp1 | Exp2 | share measured against a cross-batch baseline |
|---|---|---|---|---|
| N2 / N3 / O1 / O2 / Y2 | 20 | 10 | 10 | **50%** |
| Y1 | 18 | 10 | 8 | **56%** |

So for ~half of every donor's samples, `ΔAge = age(Exp1) − age(Exp2 baseline)`. The batch term sits
**inside the definition of `y_age`**, not in a downstream metric — everything computed since
inherits it, and nothing records which batch a baseline came from.

**D2 — baseline replication is INVISIBLE.** `_control_baseline` (`aging.py:81-90`) averages whatever
controls exist and records neither count nor composition. Stage 1.5 made `n=0` visible; **`n=1` is
still silent.** Same class as the two defects that already cost real time.

**D3 — `donor age` ground truth is UNUSED.** `grep -rn "donor age" src/ local_runners/ scripts/` →
**zero hits**; `_parse_series` reads only `days of reprogramming` and `cell type`. GEO declares:

| donor | N2 | N3 | Y1 | Y2 | O1 | O2 |
|---|---|---|---|---|---|---|
| chronological age | 0 | 0 | 29 | 35 | 53 | 53 |

It is the only external ground truth able to test what the whole ΔAge target rests on: **does this
clock read age on this data at all?** At `cv_mae ≈ 12.3 yr` the test is well powered at the extremes
(0 vs 53 ≈ 4× the error) and deliberately underpowered in the middle (29 vs 35 is half the error) —
only the extreme contrast may be claimed.

## 5.4 THE FIX PLAN — plan only, nothing executed

Sequenced so the **cheap measurements decide whether the expensive change is needed at all.** It is
entirely possible the answer is *"the clock and the baseline are the only real problems"* and Phase
3 shrinks to uncertainty-propagation plus documentation. The plan permits that answer rather than
assuming the large fix.

### Deliberately LEFT ALONE, and why

1. **The ΔAge definition (control-relative).** The design is right; the data feeding it is the
   problem. Redefining the target invalidates every prior result and all four guards.
2. **The clock — validate, never refit.** A frozen external artefact with published provenance.
   Reweighting it to improve our numbers is fitting the test.
3. **`models/`, `training/`, `evaluation/`, all Stage 1 calibration.** Stage 1 is closed at
   **PARTIAL**; re-opening it mixes changes and destroys the four-run `+0.000` guard record.
4. **Do NOT drop the Exp1 samples to "solve" D1** — a silent 50% selection is worse than the
   confound. Model or match the batch; never delete.
5. **Do NOT recruit `Failing to reprogram fibroblast` (47 samples) as baseline.** They have been
   through reprogramming. Tempting for replication, wrong biologically.
6. **Every prior record**, including §0–§4 above. Annotate, never rewrite.

### Phase 1 — measurements only; nothing rebuilt. *The decisive phase.*

New read-only `experiments/diag_zero_point.py`, in the shape of the existing diagnostics
(`dump_pool_diag.py`, `diag_calibrators.py`): pure functions, printed table, JSON dump. **Predictions
pre-registered in the lab notebook before it runs.**

- **M1 — does the clock track chronological age?** `LinearClock.predict_age` on the six day-0
  baselines vs `[0, 0, 29, 35, 53, 53]`, judged against the clock's own 12.27 yr CV MAE.
- **M2 — is there an Exp1/Exp2 batch effect?** Matched `(donor, day, marker)` comparison — measures
  the offset D1 injects.
- **M3 — bound the share** of the ±12.7 yr explainable by the single unreplicated Exp2 baseline.

| M1 | M2 | Action |
|---|---|---|
| separates 0 from 53 | no batch effect | baselines informative and unconfounded → **Phase 2 only**; Stage 2 proceeds |
| separates 0 from 53 | **batch effect present** | D1 real and quantified → **Phase 2 + Phase 3** |
| **does NOT separate 0 from 53** | either | **escalate** — the clock does not read age on this data; ΔAge's target is unvalidated. Reaches past Stage 1.5 into Stage 4, and Stage 2's premise is void as stated |

### Phase 2 — instrumentation. Value-neutral, worth doing regardless.

| File | Change |
|---|---|
| `src/cellfate/data/sources.py` (`_parse_series`) | also parse `donor age`, stamp into `obs` |
| `src/cellfate/data/aging.py` (`_control_baseline` / `delta_age`) | record per-line baseline **count + composition** (n, batch, marker). Recording only — arithmetic untouched |
| `src/cellfate/data/build_dataset.py` | persist that composition into chunk metrics |
| `verify_stage1_5.py` | flag `n_baseline < k` **and cross-batch baselines** — turn D1/D2 into a gate that can fail |

**Hard guard:** ΔAge must be **bit-identical** before/after (`max|Δ| == 0.00e+00`). It records, it
does not compute. If any ΔAge moves the change is wrong — revert, do not rationalise. This is why
Phase 2 needs no re-score and is safe even if Phase 3 never happens.

### Phase 3 — the zero-point fix. Only the option Phase 1 licenses.

A batch-matched baseline is **impossible** — no Exp1 day-0 sample exists. Realistic candidates:

- **(a)** estimate the Exp1↔Exp2 offset from matched samples and remove it before ΔAge — targets D1;
- **(b)** shrinkage baseline `λ·(donor day-0) + (1−λ)·(age-anchored grand mean)`, λ from
  between/within variance, using the newly parsed `donor age` — targets the `n=1` variance;
- **(c)** propagate baseline uncertainty into `sigma_age`, which today covers prediction spread but
  **not** the error of the zero-point it is measured against, so intervals are overconfident by
  construction. Cheapest honest option; possibly sufficient alone.

Exactly **one** ships, as its own pre-registered Change with its own bar and snapshot tag. It
changes `y_age`, so it needs a **rebuild + full re-score**, and the four guards will legitimately
move — the `+0.000` record restarts, stated in advance.

### Phase 4 — re-score, then rule on Stage 2's framing

Which this work has already shifted: from *"correct a known biological offset"* to *"replicate the
baseline so we can determine whether the offset exists."* Stage 2's k≈3 reference cells are the
right intervention either way; the justification changes, not the action.

## 5.5 Verification

| Phase | Verified by |
|---|---|
| 1 | read-only script + JSON; pure functions unit-tested (`tests/test_diag_zero_point.py`) over **every** branch — a branch that never executes is not a check |
| 2 | full `pytest -q` green; **`y_age` bit-identical** on a rebuilt fold; `verify_stage1_5.py` shows the new baseline-composition columns; `ruff check src/ tests/ scripts/` clean |
| 3 | rebuild + `scorecard.py snapshot --tag <new>` + `compare baseline <new>`; the new bar passes `audit_metrics.bar_verdict` **before** the run (ground rule §5b) and is registered in `tests/test_bars_resolvable.py` |

---

# 6. INDEPENDENT REVIEW OF §5 (2026-07-24) — verified, then tightened

> §5 was produced on the data machine. Everything in it was **re-checked against the tree rather
> than taken on trust**, including by breaking the code to confirm the new tests can fail. §5 is
> left exactly as written; this section records what verification found and the gaps it closed.

## 6.1 Verification — every checkable claim held

| Claim in §5 | How checked | Result |
|---|---|---|
| clock `cv_mae_years = 12.27` | read `configs/clocks/fleischer_clock.json` | ✅ **12.2688**, 133 samples, GSE113957 |
| `donor age` unused, 0 hits | `grep -rn` over `src/ local_runners/ scripts/` | ✅ exact; `_parse_series` reads only `days of reprogramming` + `cell type` |
| Exp1/Exp2 identity discarded | grep `sources.py` | ✅ appears only in a docstring example, never parsed into `obs` |
| `git diff --stat src/` empty | diff over the 5 commits | ✅ **`src/` untouched**; only docs, tests, verifier |
| Groups A–D 21/21, suite 303 | ran both | ✅ **21/21**, ✅ **303 passed** |
| Group E 51/51, fallback never fired | read `verify_stage1_5_results.json` | ✅ and **all six LOOCV donors present** — the PASS is not vacuous |
| every Gill donor has exactly 1 control | per-chunk census in the JSON | ✅ N2/N3/O1/O2/Y2 = 1 of 21, Y1 = 1 of 19 |

**The tests were mutation-tested — they are not decorative.** Four deliberate defects were injected
and `src/` restored after each; each was caught by the right test:

| Injected defect | Caught by |
|---|---|
| variance floor removed | `test_variance_floor_lifts_every_sigma_to_at_least_the_median` |
| control branch killed (always self-centre) | `test_control_baseline_matches_the_raw_control_mean_when_controls_exist` |
| `sigma_ref` dropped from the Gill Projection | `test_the_gain_actually_differs_between_datasets_so_it_is_not_immune` |
| `_align` made positional (ignores gene names) | `test_align_places_permuted_and_missing_genes_in_the_right_columns` |

**§5 corrected this document, and the correction is right.** §2 Group A specified intercept
cancellation as **bit-identical**; it is not. `(age+b) − mean(age_ctrl+b)` re-rounds, so the
cancellation is numerical (~1e-14), not symbolic. Independently reproduced. **§2 was wrong; the
implementation is right.**

**One concern raised and dismissed by checking.** The verifier counts controls **per chunk**, while
production `_control_baseline` groups per `cell_line` *within* a chunk — so a mixed-line chunk could
mask a fallback. Checked: every source emits one chunk per cell line by construction
(`sources.py:364`, `:459`, `:507`), so chunk↔line is 1:1 and the check is exactly equivalent. **Not
a defect** — but the invariant is nowhere asserted (see T4).

## 6.2 Gaps found in the §5.4 plan, and the tightenings that close them

**T1 — Phase 1 does not comply with the ground rule this project just adopted.** §5.5 routes only
*Phase 3* through `audit_metrics.bar_verdict`. But M1/M2/M3 each carry an implicit bar ("separates
0 from 53"), and `REF_GROUND_RULES.md §5b` requires **every** bar to be shown resolvable *before*
the run. M1 is the one that matters: with 2 samples at age 0 and 2 at 53, `SE(diff) = 12.27·√(1/2+1/2)
= 12.27 yr` against a 53 yr contrast — ~4.3σ, comfortably powered. **That is the calculation §5.4
asserts qualitatively and must instead register:** each of M1–M3 gets a pre-registered bar, a
`bar_verdict` check, and an entry in `tests/test_bars_resolvable.py` **before** `diag_zero_point.py`
runs. If a measurement has no resolvable bar, it is a description, not a test.

**T2 — M3 is measured but decides nothing.** The §5.4 decision table is M1 × M2 only; M3 ("bound
the share of the ±12.7 yr explained by the single unreplicated baseline") has no row. M3 is the
quantity that should *size* Phase 3, so it needs a decision role:

| M3 result | Consequence |
|---|---|
| baseline noise explains **most** of the offset | Stage 2's premise is reframed, not merely re-justified: there may be no donor-biology offset to correct. Phase 3 becomes **required**, and option (b) leads |
| explains **little** | the offset survives as biology-or-batch; Phase 3 is driven by M2/D1 instead, option (a) leads |
| **indeterminate** at n=6 | say so and stop — an underpowered bound is not a finding. Record it and let Stage 2's extra donors settle it |

**T3 — option (c) is partly redundant with work Stage 1 already did.** §5.4 offers "propagate
baseline uncertainty into `sigma_age`" as the cheapest option. But `sigma_scale_factor`
(`xdonor_calib.py:374`) already fits `sigma_age` to the **true out-of-donor residuals**, and those
residuals are `|pred − y_age|` where `y_age` *already contains* the baseline error. So the baseline
error is **already absorbed in magnitude**, on average, by `sigma_scale`. Option (c) therefore adds
nothing as stated. It adds value **only if made per-donor** — scaling each donor's interval by the
quality of *its own* baseline (n, batch match), which is exactly what a single global multiplier
cannot express. **Restated that way it stays on the menu; as written it should be struck.**

**T4 — two unstated preconditions.**
- **Option (a) may not be estimable.** It needs matched `(donor, day, marker)` samples spanning
  Exp1/Exp2. If no such pairs exist, the Exp1↔Exp2 offset is unidentifiable and (a) is off the
  menu regardless of M2. **M2 must report pair counts first**, and the plan must permit "(a) is
  impossible" as an outcome.
- **Phase 3 reopens Stage 1's closed verdict.** §5.4 states the four guards will move. It does not
  state that changing `y_age` also moves **both Stage 1 targets** — `conformal_coverage` (PASS) and
  `fate_ece` (MISS) are computed against `y_age`. Stage 1's PARTIAL verdict would need re-stating,
  not just its guards. That is acceptable but must be declared **before** Phase 3, not discovered
  after.

**T5 — cheap hardening for the gate (do with Phase 2).** `verify_stage1_5.py` should assert the
chunk↔line invariant it silently relies on (group by `raw.obs["cell_line"]` rather than the chunk's
metadata label), so the gate cannot weaken silently if a future source emits mixed-line chunks. One
line; no behaviour change today.

## 6.3 Standing verdict

§5 is **accepted as sound work**: the tests are real, the Group E result is meaningful and
non-vacuous, the discipline held (`src/` untouched), and the reasoning corrected this document
where it was wrong. The §5.4 plan is **directionally right and now concrete** with T1–T5 folded in.

**Phase 1 remains the correct next action** — it is read-only, cheap, and genuinely decisive: M1
can escalate past this entire stage if the clock does not read age on this data.

---

# 7. PHASE 1 EXECUTED (2026-07-24) — **M1 FAILED. ESCALATE.**

Numbers, per-donor table and full reasoning are in `experiments/DELTAAGE_LAB_NOTEBOOK.md` under
*RESULT — PHASE 1*; kept in one place so the two cannot drift. Summary:

| Measurement | Verdict |
|---|---|
| **M1** — does the clock read chronological age? | ❌ **FAIL** — extreme contrast **11.8 yr** vs bar **20.2** (true gap 53 yr). N2 (age 0) predicts **98.7**, older than both 53-year-olds; the two age-0 donors read **62 yr apart** |
| **M2** — Exp1/Exp2 batch offset | ⚠️ reported `NOT_ESTIMABLE`, but that was a **stub** (`m2_verdict([])`) and its claim "option (a) is impossible" is **false** — matched pairs exist in the series matrix. Being fixed |
| **M3** — share of offset variance from one baseline | ⏳ **INDETERMINATE** as pre-registered — 56%, 95% CI [9%, 100%] |

**The Phase 1 prediction was FALSIFIED.** It predicted `PHASE_2_AND_3` with M1 clearing.

**Consequence, per the §5.4 pre-registered branch table:** the clock does not separate the age
extremes on this data, so **ΔAge's target is unvalidated**. This reaches past Stage 1.5 into
**Stage 4 validation**, and **Stage 2's premise is void as stated**. Phases 2–4 are blocked.

**The failure is structured, not random** — and that is the lead for the escalation. O1/O2 (both 53)
agree to **0.4 yr**; across the four *adult* donors the old-vs-young separation is ≈18 yr against a
true 21 yr gap. The catastrophe is confined to the **neonatal** donors, and
`fleischer_clock.json` was fit on **adult** dermal fibroblasts (GSE113957) — so age 0 is
extrapolation outside its fitted domain. Every donor is also over-predicted (+22.7 to +98.7),
so a positive bias sits on top. **Hypothesis for the escalation to settle:** the clock may be
usable on the adults and invalid on the neonates — materially different from "the clock is broken",
and it would mean two of six LOOCV folds carry an unvalidated target.

**Also found, unrelated to the measurement:** `run_multi_local.py:53` points `CLOCK` at
`local_runners/configs/clocks/fleischer_clock.json`, **which does not exist**. `build_clock` fails
loud, so a rebuild would abort at the clock step — the "we can always harmonize again" fallback is
currently broken. The only tracked clock is `configs/clocks/fleischer_clock.json`.

## 7.1 Immediate follow-ups (before the escalation is scoped)

1. **Fix M2** to parse `(donor, day, marker, Exp)` from the series-matrix titles and actually measure
   the Exp1−Exp2 offset, then re-run and re-record. Its current verdict text must not stand.
2. **Fix the runner's clock path** to the tracked `configs/clocks/fleischer_clock.json` so a rebuild
   is possible at all.

Both are corrections to *diagnostics and wiring*, not to `src/` model or data code — `git diff
--stat src/` stays empty for this stage.
