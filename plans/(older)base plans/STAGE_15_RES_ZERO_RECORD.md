# Stage 15 — why RES is identically zero

**Status:** 🟩 **RETROSPECTIVE RECORD, written 2026-08-20. EXECUTED 2026-08-17, commit `aec6f57`.**

> ⚠️ **This is not a plan and must not be read as one.** Stage 15 ran without a plan file: its
> pre-registration lived in the docstring of `experiments/diag_stage15_res_zero.py`, which was
> written before the numbers, and its result went straight to `CHANGES.md`. This document exists so
> `plans/` is not silently missing a stage — it records what was asked and what came back, after
> the fact, and claims no pre-registration it does not have.

---

## The question

`res_median`, `res_max` and `res_approvals` had read **0.000 in every snapshot ever taken**, and
the cause was never established. An earlier "Spearman 0.40 over RES" headline turned out to be a
correlation over floating-point residue, which is what made the zero worth attributing rather than
shrugging at.

## The method

`RES = φ(S) · S^k · g(R_eff) · exp(−λ·P_loss)` is a product of four factors, so the zero can be
**attributed** rather than inferred: check all four, on all 119 held-out cells, and see which one
is zero.

## The result

| factor | outcome |
|---|---|
| `φ(S)` = sigmoid | **> 0 for all 119** — a sigmoid has no zero |
| `S^k` | **> 0 for all 119** |
| `exp(−λ·P_loss)` | **exactly 1.000** — `lam: 0.0` ships, so this factor is inert |
| **`g(R_eff)`** | **zero for 119 of 119** |

`R_eff = max(0, −(µ_age + z_conf·σ_age))`, and **σ_age is 2.0–4.5× larger than |µ_age|** on every
fold. The closest any cell comes to the credit threshold is **+2.00 yr** (N3).

**RES is not broken. It is working exactly as designed and correctly reporting that there is no
confident rejuvenation to credit.** The σ values are honest: conformal coverage is 0.714–1.000
against a nominal 0.90.

**Re-tuning the gate would not rescue it.** The largest `z_conf` at which *any* cell would qualify
is 0.235–0.898 against a shipped 1.0 — even `z_conf = 0.9` lights up one cell in one fold.

## Two further findings

**The zero is OVER-DETERMINED.** Three independent gates are all closed: `R_eff = 0` (119/119),
`REJECTED_UNSAFE` (11–12 of ~20 in five of six folds), `REJECTED_OOD` (6–7 per fold; 16 of 18 for
Y1). Fixing any one still leaves RES at zero.

**The status field understates the failure it reports.** `compute_res_batch` applies precedence
OOD → UNSAFE → NO_REJUVENATION, so the counts show **1–3 cells** as `REJECTED_NO_REJUVENATION`
when in fact **100 %** have `R_eff = 0`. Pinned by a test so it cannot mislead later.

## What it closed and what it opened

**Closed:** why RES is zero. Fully attributed.
**Open:** RES cannot become non-zero until `σ_age < |µ_age|` for some cell. That is a statement
about the model's uncertainty, not about the RES formula — the same signal-versus-noise wall as
the ΔAge work. **No `src/` change was proposed, because there is no defect here to fix.**

## Artefacts

`experiments/diag_stage15_res_zero.py` · `results/diag_stage15_res_zero_results.json` ·
`tests/test_diag_stage15_res_zero.py` (19 tests) · `CHANGES.md` 2026-08-17.
