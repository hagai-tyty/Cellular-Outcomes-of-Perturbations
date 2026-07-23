# STAGE 1.5 — Harmonization & ΔAge zero-point audit (measurement only)

**Implements:** nothing new — it *validates a claim four existing plan docs already make*.
**Depends on:** Stage 1 scored and closed (so no snapshot mixes two changes).
**Blocking for:** Stage 2 — its premise (the per-donor offset is real biology) is what this audit
tests. Not blocking for the tool.
**Scope:** 1 new test file, 1 new verifier script, **0 lines changed in `src/`**.

**Status:** PLANNED — not yet run. This document is the pre-registration.

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
