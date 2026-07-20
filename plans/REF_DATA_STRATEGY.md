# REFERENCE — Data strategy

**The highest-leverage item in the entire program**, and the one no amount of code addresses.
Covers `MASTER_PLAN.md` §6, §11.

---

## 1. The two structural limits

| Limit | Consequence | Fixable by code? |
|---|---|---|
| **One cocktail (OSKM), time-varying only** | the cross-perturbation thesis is **unposed** — not failed, unasked | **No** |
| **6 donors** | per-donor shift unlearnable without reference cells; the fate edge rests on one fold | **No** |

Plus two smaller ones:

- **Proving time-redundancy is benign** — state and time are collinear here (T11.1)
- **Statistical power on the fate edge** — 4/6 folds saturate near PR-AUC 1.0

## 2. What the next dataset needs, in priority order

| # | Requirement | Unlocks |
|---|---|---|
| **1** | **Multiple distinct perturbations** (≥5 cocktails/regimes) | makes the **core claim falsifiable for the first time** |
| **2** | **≥15 donors** | per-donor generalization becomes learnable; the fate edge becomes testable |
| **3** | **Reference/control samples per donor** | makes Stage 2 automatic instead of a protocol ask |
| **4** | **Varied time at fixed state** | would let T11.1 *prove*, not merely suggest, that time-redundancy is benign |

**Note:** T16 already shows **k ≈ 3 labelled reference cells per donor** suffices for the level
fix. That is a small, concrete ask — not a redesign.

## 3. For dose–response specifically

The identifiability wall: **"best dose and time" splits in half. Time is answerable; dose is not.**

To cross it you need a **dose × time grid**:

```
4 dose levels  ×  5 timepoints  ×  6+ donors  ≈  120 conditions
```

with paired controls per donor. Routine for a funded lab, and **the bottleneck is the wet lab, not
the modelling**.

**Critical warning on shortcuts.** Pooling across studies to manufacture a dose axis **confounds
dose with lab/batch**. Study A at 1 µg/ml vs Study B at 2 µg/ml — a difference could be dose or
could be the site. We measured **±12.7 yr between donors within a single study**; across studies it
will be worse. Your harmonization method mitigates but does not eliminate this.

**What actually works:** within-study dox titration, or within-study factor-number variation
(OSKM vs OSK vs OK), with a time course at each level.

## 4. Public datasets worth evaluating

| Dataset | What it provides | Caveat |
|---|---|---|
| **scTF-seq** (Nat Comms 2025 / bioRxiv 2024.01.30.577921) | **dose measured per cell** — dox-inducible barcoded TF overexpression, 39,187 cells, 384 TFs | mouse MSCs, single TFs not OSKM combos |
| **Schiebinger 2019 (Waddington-OT)** | dense time axis **plus a withdrawal event** (dox day 0–8, profiled to day 18) | one dox level |
| **OSKM stoichiometry** (PMC10592962) | an **OSKM ratio axis** — closest thing to dose for your exact factors | early timepoints, multiome |

Accessions are in the papers.

**Before committing compute:** run `python test_suite.py input_ablation` on any candidate dataset.
On Gill, `x+u ≈ x_only` because dose never varied. On a genuine dose dataset, `u_only` should
separate. **If it doesn't, the dose axis isn't learnable there either — stop before building.**

## 5. The strategic read

CellFate-Rx is a well-built instrument that **has not yet been pointed at an experiment capable of
testing it**. Acquiring requirement #1 is worth more than any modelling change in this program.

**The model is not the moat. The dataset is.** That is not a bad position — if you obtain a
varied-perturbation dataset, the modelling is comparatively straightforward, and the parts that
matter for a go/no-go decision are already built and validated: the fate head, the calibration
layer, the uncertainty accounting, and the per-donor correction.
