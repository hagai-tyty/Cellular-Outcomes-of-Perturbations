# Stage 10 — Is the pluripotency component of ΔAge CONTAMINATION or MEDIATION?

**Status:** PLAN. Nothing in `src/` moves until §10.5 rules. All diagnostics read-only.
**Opened because** the previous session reported "pluripotency has to come out of the ΔAge readout"
on the strength of one number (`resid_pluri` 22.69 → 13.00 MAE vs methylation). That was a
recommendation that had not been earned, and the objection to it is correct: **OSKM induces
pluripotency, and pluripotency induction may BE the mechanism of rejuvenation.** Removing it could
be deleting the signal, not the noise.

> ### 🟩 STATUS AT CLOSE — EXECUTED 2026-08-16, CLOSED
>
> Verdict **MEDIATION, not contamination** (3/3 tests). §10.5 ruled: the removal recommendation is **WITHDRAWN** -- removing the signature deletes signal (rho 0.770 -> 0.354). `src/` correctly never moved.
>
> *Commit: `3b71412`. The status line above is the pre-execution claim and is left unedited, per the append-never-rewrite rule.*

---

## 10.0 The two readings, stated so they can be told apart

| | claim | what follows |
|---|---|---|
| **CONTAMINATION** | pluripotency-correlated expression adds error to the age readout for incidental reasons | removing it purifies ΔAge; a `resid_pluri` variant is a genuine improvement |
| **MEDIATION** | OSKM → pluripotency → rejuvenation. The pluripotency component IS the effect | removing it deletes real signal; `resid_pluri` is a mutilated estimator that only *looks* better |

Both predict the observed MAE improvement is *possible*. They differ on **what else** must be true,
and that is what §10.2–10.4 test.

---

## 10.1 What is actually being removed — established, recorded here

Measured before the plan was written:

- The signature is **five genes**: `NANOG`, `POU5F1`, `LIN28A`, `SOX2`, `ZFP42`
  (`constants.py:138`, `DEFAULT_SIGNATURES["loss"]`), scored by `labels.py:74`.
- They carry **0.0005 % of the clock's total |w| mass**; ranks 17,829 / 22,836 / 24,862 / 26,826 /
  17,829 out of 33,155. **The clock barely reads them.**
- Therefore `resid_pluri` does NOT remove the clock's direct reading of pluripotency genes. It
  regresses ΔAge on a pluripotency SCORE and removes the **co-varying component across samples** —
  a far more invasive operation, and the one the objection is about.
- **Two of the five (`POU5F1` = OCT4, `SOX2`) are OSKM transgenes.** In Sendai-transduced samples
  their measured RNA includes vector transcript, so the score is partly a **vector-dose** readout
  rather than endogenous biology. Any interpretation must carry this caveat.
- `top100` **excludes all five by construction** (all rank > 17,000), so truncation removes them
  automatically — plus 33,050 other genes. The two interventions are not the same operation and
  must not be treated as interchangeable evidence.

---

## 10.2 TEST A — the falsifier for MEDIATION: does the ARM ORDERING survive?

**Logic.** The arms encode the outcome: `negative_control` (no OSKM) → `failed_to_transiently_
reprogram` → `transiently_reprogrammed` / `transient_reprogramming_intermediate`. Real rejuvenation
should order them. If pluripotency **mediates** that rejuvenation, regressing it out must **collapse
the ordering** — you have removed the causal path. If it is **contamination**, the ordering survives
or sharpens, because only noise was removed.

**Measure.** The transient-vs-control ΔAge gap, per donor, before and after `resid_pluri`.

**PRE-REGISTERED**
- **MEDIATION SUPPORTED** if the gap shrinks by **> 50 %**.
- **CONTAMINATION SUPPORTED** if it shrinks by **< 50 %**, is unchanged, or grows.
- Reported per donor (3), clustered on donor, because arms within a donor share material.

---

## 10.3 TEST B — the falsifier for CONTAMINATION: does pluripotency predict ΔAge in CONTROLS?

**Logic, and this is the sharpest of the three.** `negative_control` samples receive **no OSKM**.
There is no pluripotency induction, so there is **nothing for pluripotency to mediate**. If the
pluripotency score still predicts ΔAge among controls, that association cannot be causal
rejuvenation — it is baseline covariation between the two gene sets, i.e. contamination.

**PRE-REGISTERED**
- **CONTAMINATION SUPPORTED** if |Spearman(pluri score, ΔAge)| among control samples is
  materially non-zero (≥ 0.5) with a consistent sign.
- **MEDIATION SUPPORTED** if the association is ~0 in controls and appears only in OSKM arms —
  i.e. it is reprogramming-specific, as a causal path would be.
- n is small (controls only). An inconclusive answer is a legitimate outcome and will be reported
  as such rather than forced.

---

## 10.4 TEST C — does `resid_pluri` still track METHYLATION on the arm contrast?

**Logic.** Methylation clocks read CpGs. They **cannot see RNA pluripotency expression at all**. So
methylation's arm ordering is an outside witness to whatever rejuvenation is real. If `resid_pluri`
still reproduces that ordering, the removed component was not needed to explain what methylation
sees → contamination. If `resid_pluri` loses the methylation agreement that `raw` had, the removed
component was carrying the real signal → mediation.

**PRE-REGISTERED**
- **CONTAMINATION** if `resid_pluri`'s Spearman against methylation ΔAge on the arm contrast is
  ≥ `raw`'s.
- **MEDIATION** if it is materially lower.

---

## 10.5 Decision rule, and the ONLY thing that licenses a `src/` change

Count the three tests.

| outcome | verdict | action |
|---|---|---|
| ≥ 2 of 3 support CONTAMINATION | contamination | a `resid_pluri` ΔAge variant may be proposed as its own Change, pre-registered, with its own bar and snapshot. **Still not adopted by this stage.** |
| ≥ 2 of 3 support MEDIATION | mediation | **the previous session's recommendation is withdrawn** and recorded as withdrawn. No `src/` change. |
| split / inconclusive | undetermined | no `src/` change; record what would resolve it |

**`src/` is not touched by Stage 10 under any outcome.** The most a positive result buys is the
*right to propose* a change as a separate, separately-pre-registered Change — the project's
one-change rule. This is stated here so a favourable result cannot be read as authorisation.

---

## 10.6 What this stage cannot settle

- **Partial mediation.** If ΔAge = a·(true rejuvenation) + b·(spurious pluripotency), and the true
  rejuvenation is itself partly pluripotency-driven, all three tests can point to contamination
  while some real signal is still being removed. The tests bound the interpretation; they do not
  decompose it.
- **The transgene confound.** With OCT4 and SOX2 in the signature, "pluripotency score" and "vector
  dose" cannot be separated in this data. A cleaner signature excluding the delivered factors would
  be the follow-up, and is noted rather than done here.
- **n.** 3 donors in GSE165177. Everything here is descriptive at that scale.

---

## 10.7 Verification

| item | how |
|---|---|
| diagnostic | `experiments/diag_stage10_pluri.py`, read-only, writes `results/diag_stage10_pluri_results.json` |
| tests | `tests/test_diag_stage10_pluri.py` — every decision branch of §10.2–10.5 exercised on constructed input, since a branch that never runs is not a check |
| record | `CHANGES.md` entry stating the verdict and, if mediation, the explicit withdrawal of the earlier recommendation |
| suite | full `pytest` green, `ruff` clean, before commit |
