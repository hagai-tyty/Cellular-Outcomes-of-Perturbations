# STAGE 3 — The stopping-time tool

**Implements:** `MASTER_PLAN.md` §15 (forward model) **and** §5b-ter (condition-level RES) — the
same object, see §0.3.
**Depends on:** Stage 1 **required**; Stage 2 optional.
**Scope:** 3 new modules (~600 lines), 4 internal sub-stages.

| Sub-stage | Produces | Blocking for |
|---|---|---|
| **3a Gate** | GO / WEAK GO / STOP | everything below |
| **3b Data** | `forward_pairs.py` | 3c |
| **3c Training** | `train_forward.py` | 3d |
| **3d Decision** | `stopping.py` | Stage 4 |

---

## 0. Framing

### 0.1 What is being built

A researcher submits a transcriptome of their culture mid-protocol and receives, for each
candidate withdrawal time: predicted ΔAge with honest error bars, a calibrated P(identity loss),
and a recommended withdrawal window.

### 0.2 Why this product survives what the data cannot support

It asks a **relative** question along the **one axis this data varies**. Comparing "day 15 vs day
21 **for the same donor**" cancels the ±12.7 yr per-donor level shift. **No dose axis, no
cross-perturbation generalization** — the two things the data cannot pose.

### 0.3 The convergence

`MASTER_PLAN.md` §5b-ter argues RES must move from per-cell to **condition-level** scoring. The
arithmetic that forces it — `R_eff = max(0, −(mu + u))`, credit requires the *upper* bound negative:

| uncertainty used | mu | u | R_eff | g |
|---|---|---|---|---|
| ensemble spread (miscalibrated) | −11.0 | 2.4 | **8.6** | 0.63 |
| honest, uncorrected model | −11.0 | ~39 | **0.0** | 0.00 |
| honest, after level correction | −11.0 | ~19 | **0.0** | 0.00 |

**With honest per-cell uncertainty, R_eff = 0** — uncertainty (~19 yr) exceeds the real effect
(N2's true median ΔAge is **−11.35 yr**). But uncertainty on a *mean* shrinks by √n:

| n cells | SE of mean (q = 17–21) | vs effect −11.35 |
|---|---|---|
| 5 | **7.6–9.4 yr** | marginal |
| 10 | 5.4–6.6 yr | detectable |
| 21 | **3.7–4.6 yr** | **comfortably detectable** |

**That redesign and this tool are the same object.** Built once, in 3d.

### 0.4 Hard constraint discovered in the code

```python
# src/cellfate/common/constants.py:66
N_DOSE_TIME: int = 2               # [log10(dose_uM), log(time_h)]
# src/cellfate/common/schemas.py:76  -- a validator RAISES on any other length
```

So the forward model gets its **own bundle**, reusing the width-2 tensor but reinterpreting
column 1:

| | column 0 | column 1 |
|---|---|---|
| scoring bundle (existing, untouched) | log10 dose | **log absolute time** |
| forward bundle (new) | log10 dose | **log Δt** |

A `mode` field records which, **asserted at load** (3d.6).

> ⚠️ **BLOCKER FOUND IN AUDIT.** `BundleMeta` is declared `model_config = ConfigDict(extra="forbid")`
> (`schemas.py:269`). **You cannot simply set `bundle_meta.mode`** — pydantic will reject the extra
> key and the bundle will fail to write. The field must be **declared** first:
>
> ```python
> # src/cellfate/common/schemas.py, inside class BundleMeta
> mode: str = "scoring"          # "scoring" | "forward"; defaulted so OLD BUNDLES STILL LOAD
> ```
>
> A declared field with a default is compatible with `extra="forbid"`, and every existing bundle
> keeps validating. **Do this before writing any forward bundle.**

---

# SUB-STAGE 3a — The gate

**Run before writing any tool code:**

```powershell
python test18_forward_gate.py
```

**What it tests.** Test 11.1 showed the current model ignores its time input entirely (ΔAge shifts
**0.035 yr** across a full sweep) because time is redundant with state along one trajectory. The
hope is that training on `(t_i → t_j)` **pairs** fixes this: with pairs the same starting state
appears with *different* Δt and *different* targets, so Δt is no longer redundant.

**That is a hypothesis.** The gate checks it with **ridge** — because this project has established
repeatedly (T3, T5, T6, T9) that at this scale ridge matches or beats flexible models. **If ridge
cannot find a Δt signal, a neural net will not either.**

| Verdict | Condition | Action |
|---|---|---|
| **GO** | Δt improves prediction beyond noise **and** the sweep moves >2 yr | build 3b–3d |
| **WEAK GO** | sweep moves but Δt does not beat state-only | build, tempered |
| **STOP** | neither | **do not write tool code.** Ship the scoring model; go to Stage 5 |

**A STOP is a real result** — this dataset cannot support forward prediction, which is worth
knowing and worth reporting.

---

# SUB-STAGE 3b — The data layer

## 3b.1 New file: `src/cellfate/data/forward_pairs.py`

```python
"""Forward (t_i -> t_j) pairs for the stopping-time model.

Sampling is destructive - the same cell is never observed twice - so a cell-to-cell pairing does
not exist. The unit is a POPULATION SNAPSHOT: mean expression at t_i paired with the mean outcome
at t_j.

From 12 timepoints this yields up to 66 ordered pairs per donor (~396 across six donors), versus
the 6 examples a day-0 -> endpoint pairing would give.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.data import gather_split
from cellfate.common.constants import LOSS_IDX, DEATH_IDX

MIN_CELLS_PER_TP = 1          # raise to 3 if the data allows
REGIME = "holdout"


def _resolve_root(name: str) -> str:
    """Fold folders may sit at the repo root or under runs/ - check both."""
    from pathlib import Path
    for base in (".", "runs", ".."):
        if (Path(base) / name).exists():
            return str(Path(base) / name)
    return name


@dataclass(frozen=True)
class ForwardPair:
    donor: str
    x_i: np.ndarray          # (G,) population-mean expression at t_i
    t_i: float               # log time_h of the source timepoint
    dt: float                # log-gap to the target timepoint  <-- THE NEW SIGNAL
    y_age_j: float           # mean TRUE ΔAge at t_j        (age-valid cells only)
    y_frac_unsafe_j: float   # fraction loss/death at t_j   (all cells)
    n_i: int
    n_j: int


def _timepoint_rows(donor: str, regime: str = REGIME):
    """Per-(donor, timepoint) population summaries."""
    te = gather_split(ArtifactPaths.of(_resolve_root(f"cellfate_loocv_{donor}")), regime, "test")
    t_all = np.asarray(te.dose_time[:, 1], float)
    X_all = np.asarray(te.X, float)
    y_all = np.asarray(te.y_age, float)
    cls   = te.y_cls.astype(int)
    am    = te.mask.astype(bool)

    rows = []
    for tp in np.unique(np.round(t_all, 6)):
        sel = np.isclose(t_all, tp)
        if sel.sum() < MIN_CELLS_PER_TP:
            continue
        sel_age = sel & am                       # NOTE: two different masks - do not conflate
        unsafe = (cls[sel] == LOSS_IDX) | (cls[sel] == DEATH_IDX)
        rows.append({
            "t": float(tp),
            "x": X_all[sel].mean(0),             # ALL cells -> the state
            "y_age": float(y_all[sel_age].mean()) if sel_age.any() else np.nan,
            "frac_unsafe": float(unsafe.mean()),
            "n": int(sel.sum()),
            "n_age": int(sel_age.sum()),
        })
    return sorted(rows, key=lambda r: r["t"])


def build_forward_pairs(donors, regime: str = REGIME) -> list[ForwardPair]:
    out = []
    for d in donors:
        rows = _timepoint_rows(d, regime)
        if len(rows) < 3:
            continue
        for i, ri in enumerate(rows):
            for rj in rows[i + 1:]:              # ordered: t_j > t_i, guaranteed by the sort
                if not np.isfinite(rj["y_age"]):
                    continue                     # target has no age-valid cells
                out.append(ForwardPair(
                    donor=d, x_i=ri["x"], t_i=ri["t"], dt=rj["t"] - ri["t"],
                    y_age_j=rj["y_age"], y_frac_unsafe_j=rj["frac_unsafe"],
                    n_i=ri["n"], n_j=rj["n"]))
    return out


def pairs_to_arrays(pairs):
    """-> (X, dose_time, y_age, y_unsafe, donor).  dose_time[:,1] = log Δt, NOT absolute time."""
    X  = np.vstack([p.x_i for p in pairs]).astype(np.float32)
    dt = np.zeros((len(pairs), 2), np.float32)
    dt[:, 1] = [p.dt for p in pairs]             # column 1 REINTERPRETED - see §0.4
    return (X, dt,
            np.array([p.y_age_j for p in pairs], np.float32),
            np.array([p.y_frac_unsafe_j for p in pairs], np.float32),
            np.array([p.donor for p in pairs]))
```

## 3b.2 The leakage rule — the single most important constraint in this stage

> **Split by DONOR. Never split by pair.**

Pairs within a donor share timepoints: `(t1→t5)` and `(t2→t5)` share the **target**; `(t1→t5)` and
`(t1→t9)` share the **source**. A random pair split therefore puts the *same measurement* on both
sides and will produce a beautiful, meaningless result.

```python
# CORRECT
train = [p for p in pairs if p.donor != held_out]
test  = [p for p in pairs if p.donor == held_out]

# CATASTROPHIC - never do this
train, test = train_test_split(pairs, test_size=0.2)     # shares timepoints across the split
```

**Stage 4 re-verifies this independently.** It is the most likely explanation for a suspiciously
good result.

## 3b.3 Sanity checks the builder must print

```python
def audit(pairs):
    import collections
    per = collections.Counter(p.donor for p in pairs)
    dts = np.array([p.dt for p in pairs])
    assert (dts > 0).all(),            "ordering bug: non-positive Δt"
    assert dts.std() > 0,              "all Δt identical - Δt cannot be learned"
    assert all(v >= 3 for v in per.values()), f"donor with <3 pairs: {per}"
    print(f"pairs={len(pairs)} donors={len(per)} "
          f"Δt min/med/max = {dts.min():.2f}/{np.median(dts):.2f}/{dts.max():.2f}")
    print("per donor:", dict(per))
```

| Check | Failure means |
|---|---|
| `dt > 0` everywhere | ordering bug |
| `dt` spans a usable range | if all gaps are similar, Δt cannot be learned — **stop** |
| ≥3 pairs per donor | that donor cannot be held out |
| no pair crosses donors | leakage |

## 3b.4 The weighting decision — make it before coding

A pair built from 1 cell at `t_i` is far noisier than one from 20. Options: (a) unweighted,
(b) weight by `min(n_i, n_j)`, (c) drop pairs below a threshold.

**Recommendation: (c) with a low threshold, then unweighted.** Weighting adds a hyperparameter this
dataset is too small to tune honestly — and this project has repeatedly found added flexibility
hurts at n≈100.

## 3b.5 Acceptance for 3b

- [ ] `build_forward_pairs(DONORS)` returns a non-empty list
- [ ] `audit()` passes every assertion
- [ ] Δt range is reported and non-degenerate
- [ ] **no training code is written until this holds**

---

# SUB-STAGE 3c — Training

## 3c.1 New file: `src/cellfate/training/train_forward.py`

Mirrors `train_model.py` with **three** differences. Everything else — `CellFateNet`, encoders,
heads, scalers, the bundle writer — is **reused unchanged**.

| # | Difference | Why |
|---|---|---|
| 1 | input is 3b pairs; `dose_time[:,1] = dt` | the forward signal |
| 2 | targets are `y_age_j` and `y_frac_unsafe_j` | population outcomes |
| 3 | calibration via Stage 1's `crossdonor_stats` | never in-distribution |

## 3c.2 The ordering rule

```
   train model on pairs
        ↓
   inner-LODO over training donors  ->  pooled cross-donor residuals / logits / features
        ↓
   fit temperature, q, sigma_scale, OOD on THOSE
        ↓
   write bundle   (bundle_meta.mode = "forward")
```

> **Calibration is fitted last, to the model that will ship.** Reusing the existing bundle's
> calibration, or fitting before the final training run, reproduces exactly the defect Test 14
> measured (coverage 0.40 vs 0.90).

## 3c.3 The head decision

The existing fate head is 3-class trained on per-cell labels. 3b targets are **population
fractions**.

| Option | Implication |
|---|---|
| (a) regress the unsafe fraction — single sigmoid | simplest; **loses the loss-vs-death distinction** |
| **(b) keep 3-class, train on soft labels** (class proportions) | preserves it; needs soft-target cross-entropy |

**Recommendation: (b).** Losing identity and dying are different failures with different
consequences, and soft-target cross-entropy is a two-line change:

```python
# hard labels:  F.cross_entropy(logits, targets.argmax(1))
# soft labels:  -(targets * F.log_softmax(logits, dim=1)).sum(1).mean()
```

## 3c.4 What must NOT change

- `N_DOSE_TIME = 2` — hard contract with a validator that raises
- the existing scoring bundle — untouched, still measured by `scorecard.py`
- `CellFateNet` — same architecture. This project's own evidence (T3, T6, T9) is that added
  capacity **hurts** at this scale; what changes is what it trains on, not how big it is

## 3c.5 Expected effects — record before running

| Metric | Now | After 3c | Basis |
|---|---|---|---|
| **`dt_response`** | **0.035 yr** | **> 2 yr** | the 3a gate threshold |
| forward coverage | — | 0.85–0.95 | Stage 1's method |
| forward ΔAge MAE | — | comparable to per-cell MAE | — |

> **If `dt_response` stays near zero after training on pairs, 3c has FAILED** even if the losses
> look fine — the model has again learned to read time off state. **Check this before 3d.**

---

# SUB-STAGE 3d — The decision layer

## 3d.1 New file: `src/cellfate/inference/stopping.py`

```python
"""Stopping-time recommendations from the forward bundle.

Answers: "given my culture as it is now, when should I withdraw?" -- a RELATIVE question along
the one axis this data varies, which is why it survives the ±12.7 yr per-donor level shift.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Sequence
import numpy as np

from cellfate.inference import Predictor


@dataclass
class StoppingOption:
    dt_hours: float
    withdraw_day: float
    delta_age: float
    delta_age_lo: float      # calibrated, POPULATION-level
    delta_age_hi: float
    p_unsafe: float          # calibrated
    in_range: bool           # False => Δt outside the trained range


@dataclass
class StoppingReport:
    options: list[StoppingOption]
    recommendation: StoppingOption | None
    reason: str
    warnings: list[str] = field(default_factory=list)


def recommend_stopping(X_now, bundle_root, risk_threshold=0.20,
                       dt_grid_hours: Sequence[float] | None = None,
                       dt_range=None, donor_offset=None) -> StoppingReport:
    pred = Predictor(bundle_root)

    # --- 1. ASSERT the right bundle. Loading the scoring bundle here silently produces
    #        nonsense, because column 1 means absolute time there and Δt here.
    #  AUDIT FIX: the attribute is `pred.meta` (predictor.py:67), NOT `pred.bundle_meta`.
    #  Using getattr(pred, "bundle_meta", {}) would silently yield {} -> mode=None -> the
    #  assert would fire on EVERY bundle, blocking the tool entirely. ---
    mode = getattr(pred.meta, "mode", "scoring")
    if mode != "forward":
        raise ValueError(f"stopping requires a FORWARD bundle, got mode={mode!r}")

    X_now = np.atleast_2d(np.asarray(X_now, np.float32))
    n_cells = X_now.shape[0]
    warnings: list[str] = []

    # --- 2. population mean: the tool scores CONDITIONS, not cells (§0.3) ---
    x_pop = X_now.mean(0, keepdims=True)

    grid = list(dt_grid_hours) if dt_grid_hours is not None else _default_grid(dt_range)
    options = []
    for dt in grid:
        dose_time = np.array([[0.0, np.log(dt)]], np.float32)   # col 1 = log Δt
        # AUDIT FIX: Predictor has no `.arch`; the architecture lives on the members.
        # NOTE: arch holds BOTH "n_fp" (fingerprint width) and "n_pert" (actual input width;
        # network.py:43 sets n_pert = len(TF_VOCAB) for tf-kind, else n_fp). Use n_pert -
        # "n_fp" is wrong for TF cocktails, which is what this dataset uses.
        fp = np.zeros((1, pred.members[0].arch["n_pert"]), np.float32)
        row = pred.predict_encoded(x_pop, fp, dose_time)[0]

        mu = float(row["mu_age"])
        if donor_offset is not None and donor_offset.applied:
            mu -= donor_offset.d                                # Stage 2, optional

        # --- 3. population interval: sigma/sqrt(n) is a LOWER BOUND (cells are correlated) ---
        half = pred.q / np.sqrt(max(n_cells, 1))
        options.append(StoppingOption(
            dt_hours=dt, withdraw_day=_to_day(dt),
            delta_age=mu, delta_age_lo=mu - half, delta_age_hi=mu + half,
            p_unsafe=float(1.0 - row["S"]),
            in_range=(dt_range is None or dt_range[0] <= np.log(dt) <= dt_range[1]),
        ))

    # --- 4. recommendation rule ---
    ok = [o for o in options if o.p_unsafe <= risk_threshold and o.in_range]
    rec = min(ok, key=lambda o: o.delta_age) if ok else None
    reason = ("most rejuvenation with P(unsafe) below threshold" if rec
              else f"no withdrawal time meets P(unsafe) <= {risk_threshold:.0%}")

    # --- 5. mandatory warnings ---
    warnings.append("intervals are population-level and assume independent cells "
                    "(LOWER BOUND on true uncertainty)")
    if donor_offset is None or not donor_offset.applied:
        warnings.append("absolute ΔAge carries a ±12.7 yr per-donor offset; comparisons WITHIN "
                        "this donor are reliable, absolute values are not")
    if n_cells < 10:
        warnings.append(f"only {n_cells} cells - uncertainty estimate unreliable below ~10")
    if any(not o.in_range for o in options):
        warnings.append("some options are outside the trained Δt range and were not recommended")
    if rec is None:
        warnings.append("NO SAFE WITHDRAWAL WINDOW FOUND")

    return StoppingReport(options, rec, reason, warnings)
```

## 3d.2 Refuse-to-extrapolate

`in_range = dt_min_observed <= log(dt) <= dt_max_observed`, taken from 3b.

Out-of-range options are **computed and shown but never recommended**, and flagged. The model has
no information beyond the observed gaps, and **silently extrapolating a safety prediction is the
single most dangerous thing this tool could do.**

## 3d.3 The √n caveat is not decoration

`sigma/√n` assumes independent cells. **Cells from one culture are correlated**, so this is a
**lower bound** on true uncertainty. It is stated in the report, printed to the user, and repeated
in the manuscript limitations.

## 3d.4 Mandatory warnings

| Condition | Warning |
|---|---|
| always | population intervals are a lower bound |
| no donor offset applied | absolute ΔAge carries ±12.7 yr; **within-donor comparisons are the reliable output** |
| `n_cells < 10` | uncertainty unreliable |
| any option out of range | not recommended, shown only |
| OOD flag fires | "this culture is unlike the training trajectories" |
| no option meets threshold | **no safe window found** |

**`recommendation = None` is a correct, useful answer** — it means this culture has no safe
window, which is exactly what a researcher needs to hear when true.

## 3d.5 Output the researcher sees

```
Culture: donor D, sampled day 10, 847 cells

  withdraw day 12   ΔAge  −6 [−10, −2]   P(unsafe)  8%
  withdraw day 15   ΔAge  −9 [−13, −5]   P(unsafe) 14%
  withdraw day 18   ΔAge −11 [−15, −7]   P(unsafe) 29%   ← above threshold
  withdraw day 21   ΔAge −12 [−16, −8]   P(unsafe) 47%   ← above threshold

  RECOMMENDATION: withdraw day 15
    most rejuvenation with P(unsafe) below your 20% threshold

  ⚠ intervals are population-level, assume independent cells (lower bound)
  ⚠ ΔAge comparisons within this donor are reliable; absolute values carry ±12.7 yr
```

## 3d.6 Training/serving skew — the silent failure

The live sample **must** be preprocessed exactly as 3b did: same gene panel, same scaler, same
population-mean step. **A mismatch here will not announce itself** — it produces confident wrong
answers.

**AUDIT FIX:** `Predictor` has no `.panel`. The panel identity lives in the metadata, and the
Predictor already cross-checks it against the scalers at construction (`predictor.py:74`). Use:

```python
# gene count: the scalers carry the panel width
assert X_now.shape[1] == len(pred.scalers.params.x_mean), \
    f"expected {len(pred.scalers.params.x_mean)} genes, got {X_now.shape[1]}"
# panel identity: meta vs scalers is already asserted inside Predictor.__init__,
# so a successful load already guarantees they agree. Re-assert only against YOUR data:
assert forward_pairs_panel_hash == pred.meta.gene_panel_hash, \
    "the pairs were built on a different gene panel than the bundle was trained on"
```

## 3d.7 Explicit non-goals

- **No per-cell recommendations** — per-cell uncertainty exceeds the effect size
- **No dose recommendations** — identifiability wall, one dose in the data
- **No cross-donor absolute claims** without reference cells

## 3d.8 Failure modes

| Symptom | Cause | Action |
|---|---|---|
| `ValueError: requires a FORWARD bundle` | wrong bundle loaded | correct — this assert exists to prevent silent nonsense |
| every option flagged out-of-range | `dt_range` not passed from 3b | pass it; do not disable the guard |
| `recommendation` always `None` | risk threshold too strict, or fate uncalibrated | check Stage 1 acceptance passed |
| ΔAge identical across all Δt | `dt_response` ≈ 0 → **3c failed** | return to 3c; do not ship |
| intervals absurdly narrow | `q` from the scoring bundle, not the forward one | check `bundle_meta.mode` |

## 3d.9 Acceptance for Stage 3

```powershell
python scorecard.py snapshot --tag C_forward
python scorecard.py compare B_percalib C_forward
```

| Role | Metric | Bar |
|---|---|---|
| **TARGET** | `dt_response` | **> 2 yr** (from 0.035) |
| **TARGET** | `forward_coverage` | 0.85–0.95 |
| **GUARD** | the existing scoring metrics | must be **noise** — the forward bundle is a sibling, not a replacement |

## 3d.10 Done when

`recommend_stopping` runs on a held-out donor's day-N sample, returns a report whose intervals are
calibrated, and every §3d.4 warning fires when it should.
