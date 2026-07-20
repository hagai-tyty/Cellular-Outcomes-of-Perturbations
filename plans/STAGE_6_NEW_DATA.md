# STAGE 6 — Acquiring and integrating new data

**Why this is a STAGE and not a reference note:** "we need more data" is not a conclusion, it is
work with downloads, checks, gates and acceptance criteria. This file is that work.

**Depends on:** nothing — it can run in parallel with Stages 1–5.
**Blocking for:** the cross-perturbation thesis, and possibly the tool (see §2).

---

## 0. READ THIS FIRST — the species blocker

**Every public dose-varying or dense-time reprogramming dataset I could find is MOUSE.**

| Dataset | Species | Scale | What it has |
|---|---|---|---|
| **scTF-seq** (*Nat Genet* 2025, s41588-025-02343-7) | **mouse** MSCs (C3H10T1/2) | 39,187 cells, 384 TFs | **dose measured per cell** |
| **Schiebinger 2019** (Waddington-OT) | **mouse** embryonic fibroblasts | **165,892 cells, 39 timepoints** | dense time + **dox withdrawal at day 8** |
| **OSKM stoichiometry** (PMC10592962) | **mouse** fibroblast | — | OSKM ratio axis |

**Your clock is the Fleischer human dermal fibroblast clock. None of these yields ΔAge directly.**

Do **not** orthology-map and assume the human clock transfers. We measured **±12.7 yr** of shift
between **human donors within one study**; cross-species will be worse, and there is no reference
to calibrate against. That path produces confident nonsense.

**What is still usable, and for what:**

| Use | Mouse data OK? | Why |
|---|---|---|
| **Fate** (safe / loss / death) | ✅ **yes** | identity and pluripotency markers have clear orthologs; no clock needed |
| **Forward dynamics as METHOD validation** | ✅ **yes** | whether Δt is learnable is a modelling question, not a species question |
| **ΔAge** | ❌ **no** | needs a species-matched clock |
| **The product itself** | ❌ **no** | it must run on human cells with a human clock |

---

## 1. The opportunity hiding in the blocker

Recall the load-bearing problem from the gate: your donors have ~21 cells spread across the time
course, so per-timepoint uncertainty may exceed the effect size.

**Schiebinger has ~4,250 cells per timepoint** (165,892 / 39). That is **~2,400× more per
timepoint**, and it includes a **dox withdrawal event at day 8** — structurally *exactly* the
stopping-time question.

> **So Schiebinger can answer "does the forward stopping-time METHOD work at all?" free of the
> statistical-power confound.** If the method fails there, with thousands of cells per timepoint,
> it will certainly fail on 21. If it succeeds, you have separated *"the method is wrong"* from
> *"we lack cells"* — which the Gill data alone cannot do.

That makes Stage 6 worth running **before** committing to Stage 3 on thin data.

---

## 2. Sub-stage 6a — Method validation on Schiebinger (mouse, fate + dynamics only)

### 6a.1 Acquire

```powershell
# GEO accession is in Schiebinger et al. 2019 (Cell 176:928-943). Broad also hosts a
# processed version via the Waddington-OT documentation.
# Download to D:\schiebinger\  -- do NOT put it in D:\cellfate-rx\
```

**Confirm before proceeding:**
```python
import scanpy as sc
a = sc.read_h5ad("D:/schiebinger/reprogramming.h5ad")
print(a.shape)                                   # expect ~165k cells
print(sorted(a.obs.columns))                     # need a day/time column
print(a.obs["day"].value_counts().sort_index())  # expect ~39 timepoints
print("cells/timepoint:", a.n_obs / a.obs["day"].nunique())
```

**Gate 6a.1:** ≥1,000 cells per timepoint and ≥10 timepoints. If not, this dataset does not solve
the power problem and 6a is pointless.

### 6a.2 Build fate labels (no clock required)

The existing fate definition uses marker genes. Map to mouse orthologs:

| Class | Human markers | Mouse orthologs |
|---|---|---|
| pluripotent / loss-of-identity | `POU5F1`, `NANOG`, `LIN28A` | `Pou5f1`, `Nanog`, `Lin28a` |
| somatic identity retained | `COL1A1`, `THY1`, `S100A4` | `Col1a1`, `Thy1`, `S100a4` |
| death / stress | apoptosis panel | same panel, mouse symbols |

**Mouse gene symbols are title-case** — a case-sensitive lookup silently returns zero matches and
you get an all-one-class dataset. Assert:
```python
assert (a.var_names.str.match(r"^[A-Z][a-z]")).mean() > 0.5, "expected mouse symbol casing"
```

### 6a.3 Run the gate on Schiebinger

```powershell
python test18_forward_gate.py --data D:\schiebinger --targets fate
```

*(needs a `--data` flag adding to the gate; ΔAge parts will be skipped for lack of a clock)*

| Outcome | What it means | Action |
|---|---|---|
| **Part C passes with thousands of cells/tp** | the forward safety method **works** when power is adequate | your Gill limitation is **cells, not method** → pool timepoints (Stage 3) or get more human cells |
| **Part C fails even here** | the method does not work regardless of power | **abandon the forward tool.** Ship the calibrated readout |

> **This is the highest-information run in the whole program.** It separates a power problem from a
> method problem, and those have completely different responses.

---

## 3. Sub-stage 6b — The human dose dataset does not exist publicly

Searched and not found: a **human**, within-study, dose-titrated reprogramming scRNA-seq dataset
with a time course.

**Therefore the dose axis requires generating data.** The minimum design:

```
4 dose levels  ×  5 timepoints  ×  6+ donors  ≈  120 conditions
+ paired day-0 controls per donor     (which ALSO makes Stage 2 automatic)
```

**Do not manufacture a dose axis by pooling studies.** Study A at 1 µg/ml vs Study B at 2 µg/ml
confounds dose with site. We measured ±12.7 yr between donors *within* one study; across studies it
is worse, and the harmonization method mitigates but does not eliminate it.

**Acceptable shortcuts, in order of preference:**
1. within-study **dox titration** with a time course
2. within-study **factor-number variation** (OSKM vs OSK vs OK)
3. nothing else

---

## 4. Sub-stage 6c — Integrating any new human dataset

Only for a **human** dataset. The pipeline steps, in order:

| # | Step | Gate |
|---|---|---|
| 1 | align to the 2000-gene panel | ≥80% of panel genes present; **record which are missing** — T0 showed the panel already captures only ~47% of clock signal |
| 2 | apply the Gill Projection harmonization | intercept-cancellation unit test still passes |
| 3 | verify the clock applies | human fibroblast-like; otherwise ΔAge is invalid |
| 4 | build fate labels from markers | all three classes present, none below ~5% |
| 5 | **run `python test_suite.py input_ablation`** | **the decisive gate — see §5** |
| 6 | re-run `scorecard.py snapshot --tag <dataset>` | compare against the Gill baseline |

## 5. The decisive gate for ANY new dataset

```powershell
python test_suite.py input_ablation
```

On Gill, `x+u ≈ x_only` because **dose never varied**. On a genuine dose dataset, `u_only` should
**separate** — the perturbation channel should carry information that cell state does not.

| Result | Meaning | Action |
|---|---|---|
| `u_only` separates from `x_only` | the dose axis is **real and learnable** | proceed to re-pair and retrain |
| `u_only ≈ x_only` | the axis is **not learnable here either** | **stop.** Do not spend compute |

**One hour, on data you can download today, and it answers "weeks or blocked" definitively.**

## 6. What each acquisition unlocks

| Acquisition | Unlocks | Does NOT unlock |
|---|---|---|
| Schiebinger (mouse) | **method validation** free of the power confound; forward-fate proof of principle | ΔAge, the product |
| scTF-seq (mouse) | whether dose is learnable **in principle** | ΔAge, human transfer |
| **Human dose × time grid** (must be generated) | **the core cross-perturbation thesis**, dose–response, the full product | — |
| **≥15 human donors** | per-donor generalization; statistical power on the fate edge | dose |

## 7. Honest ordering

1. **6a Schiebinger method validation** — cheap, public, and separates *power* from *method*.
   **Do this before committing to Stage 3.**
2. **6c integration** for any human dataset that appears.
3. **6b generating a dose grid** — the highest-value item and the slowest. Start conversations
   early; the bottleneck is bench time, not modelling.

## 8. Done when

Either 6a has returned a verdict on whether the forward method works when power is adequate, or a
new human dataset has passed the §5 gate and been integrated — **and the result is recorded in the
lab notebook, including if the answer is "this data does not help."**
