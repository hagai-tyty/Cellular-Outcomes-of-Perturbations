# STAGE 6 (REV) — Acquisition, rewritten against the data you actually hold

**Status:** 🔵 **PLAN.** Supersedes `STAGE_6_NEW_DATA.md` for inventory, sizing and ordering.
**The original is left byte-unmodified** — its species blocker (§0) and its `input_ablation` gate
(§5) are still correct and are carried forward here unchanged.

**Why a rewrite:** the original names **2 datasets — GSE165176 and GSE242423.** You now hold **six**,
and four of the questions it was written to answer have since been closed by measurement. Executed as
written it would buy data to answer questions the data already answers.

---

## 1. What is now closed, and what that removes from the shopping list

| question | status | what it removes |
|---|---|---|
| Can the RNA clock be calibrated against methylation? | ❌ **NO** — 1.5.2 M-2a, SPLIT | — |
| Can a model *learn* methylation age from RNA? | ❌ **NO** — 1.5.4, SPLIT on all 3 families | **removes "acquire more RNA for the age arm"** |
| Is HFF's ΔAge the identity artefact? | ❌ **NO** — 1.5.5, R² ≤ 0.16 within timepoint | **removes "discard HFF" as the obvious move** |
| Is it a sequencing-depth artefact? | ❌ **NO** — 1.5.5, R² ≤ 0.09 | as above |
| Do HFF's labels carry exploitable structure? | ✅ **YES** — arm C, ranking fell 5.4× the A–B gap when permuted | — |
| Does masking HFF change ΔAge MAE? | ⚠️ **INCONCLUSIVE** — +0.661 yr, MDE 5.045 > Δ\* 3.572 | — |

**The RNA route is closed on measurement, not inference.** That is the single biggest change from the
original document, which still treats more RNA as the default purchase.

**Still open, and narrow:** whether HFF's residual structure is real biological signal or clock noise.
Everything cheaper than that has been eliminated.

---

## 2. The real inventory

| dataset | what it is | n | age labels? |
|---|---|---|---|
| **GSE113957** | Fleischer fibroblast RNA — the clock's own training set | 133 | the clock itself |
| **GSE165176** | Gill Sendai RNA — **the age arm trains on this** | 124 | 75 (6 donors) |
| **GSE165177** | Gill transient RNA — **on disk, in no training config** | 95 | 3 donors, ages 38/53 |
| **GSE165178** | Gill Sendai methylation — joins GSE165176 **22/22** | 22 | working instrument |
| **GSE165179** | Gill transient methylation — joins GSE165177, 68 conditions | 96 | working instrument |
| **GSE242423** | HFF scRNA, D0→D14 + iPSC — **99.8 % of all age labels** | 44,473 cells | RNA clock only |

**~90 samples already carry BOTH the model's input (RNA) and a working instrument for the label
(methylation)** — 68 from 165177×165179, 22 from 165176×165178. That set has never been used to
train anything.

---

## 3. The binding constraint is DONORS, and here is the number

Every unresolved statistic in this project is donor-limited, not sample-limited. Sizing derived from
the **measured** instrument errors, not guessed:

| open question | measured spread | donors needed |
|---|---|---|
| step 6's arm comparison, σ treated as known | SD 4.808 | **10** |
| step 6's arm comparison, **σ-robust** *(the honest bar)* | SD 4.808 | **19** |
| contrast B retention, multi-tissue | SD 11.6 | **17** |
| contrast B retention, skin & blood | SD 18.2 | **38** |
| M3 — is the per-donor offset real or `n=1` noise? | σ CI factor 3.93 at n=6 | **~20** to reach a factor of 2 |

**σ-robust means the conclusion survives σ being an estimate** — the bar the step-6 rerun correctly
*refused* to claim at n=4, when the χ² interval on σ put the MDE at 6.704 against Δ\* 3.572.

---

## 4. 🔑 The finding that changes what to buy

**Restricting to donors inside the clock's fitted range collapses the requirement.** When N2/N3
(donor age 0, outside `[1, 96]`) are dropped, the per-fold difference SD falls **4.808 → 1.130**, a
factor of 4.3 measured on 2 of 6 folds:

| | σ-robust donors needed |
|---|---|
| all donors, as measured | **19** |
| **in-range donors only** | **6** |

**You already hold 4 in-range donors (O1, O2, Y1, Y2). The gap is 2 — not 15.**

> **Fixing the instrument's domain is worth more than tripling the donor count.** C-2 — the
> `enforce_clock_age_range` switch — is already built and shipped inert. Turning it on is a
> pre-registered label change, and this is the arithmetic that justifies paying for it.

That reframes the acquisition entirely: **the ask is 2 adult in-range donors with paired methylation,
not a large cohort.**

---

## 5. Carried forward unchanged from the original

**§0's species blocker still holds.** Every public dose-varying or dense-time reprogramming dataset
is **mouse** (scTF-seq, Schiebinger, OSKM-stoichiometry). Mouse is usable for **fate** and for
**method validation**; it cannot give ΔAge, because that needs a species-matched clock. Do not
orthology-map and assume transfer — ±12.7 yr of shift was measured between *human donors within one
study*.

**§5's decisive gate still holds.** For any dose dataset, `python experiments/test_suite.py
input_ablation`: if `u_only` does not separate from `x_only`, the dose axis is not learnable there
either — stop before spending compute.

---

## 6. ⛔ What to acquire — GATED

**This section is deliberately not final.** The remaining open question — is HFF's residual structure
real signal or clock noise? — is being settled by the **stratified shuffle** (permute within
donor/timepoint, preserving the artefact's between-stratum structure while destroying within-stratum
pairing). Its result changes the target:

| stratified shuffle says | acquisition target |
|---|---|
| HFF's structure is **systematic artefact** | HFF must be **replaced**, not supplemented. Buy adult in-range donors at volume |
| HFF's structure is **real signal** | HFF is an asset. Buy only the **2 in-range donors** §4 sizes, to close step 6 and M3 |

**Writing a target now would be guessing at the one number this stage exists to get right.** Filled
in when the shuffle lands.

---

## 7. Honest ordering

1. **Turn on C-2** (`enforce_clock_age_range`) as a pre-registered label change. §4 shows it is worth
   more than 13 extra donors. **Free — the code is already written and inert.**
2. **Integrate GSE165177.** 95 adult in-range samples, all methylation-paired, in no training config.
   Free.
3. **The stratified shuffle.** Settles §6's gate. Free.
4. **Then acquire**, sized by §4 and targeted by §6.

**Steps 1–3 are free and must precede step 4**, because each one changes what step 4 should buy.

---

## 8. What no acquisition fixes

* **HFF cannot be anchored by more donors.** It is one cell line; its ΔAge question is settled by
  methylation *on HFF*, not by cohort size.
* **The RNA route stays closed.** 1.5.4 measured it. More RNA does not reopen it.
* **Contrast B on skin & blood needs 38 donors** — likely out of reach, and that should be stated as
  a limit rather than quietly carried as a plan.
