# C-7 — IMPLEMENTATION SPEC

**Status:** 🔵 **SPEC FOR REVIEW. Nothing written. `src/` untouched.**
**Implements:** `plans/CHANGE_C7_BULK_SAMPLE_INTEGRITY.md` §2 (the gate), §9 (option (c) + rule 4),
§10 (predicate host + B2′), §11 (B1 verified, the five blockers).
**Ships:** OFF. Enabling it is a separate pre-registered run, exactly as C-2 did.

> **Five decisions are marked 🟡 DECIDE. They change what gets built and I have not chosen them
> unilaterally.** Everything else is determined by code that already exists.

---

## 0. What is being built, in one paragraph

Four pieces, which must land **together** because any three of them leave the system in a state
C-7's own bars forbid: **(A)** a pure integrity gate that rejects degenerate bulk columns; **(B)** a
global per-`cell_line` control census computed in an **unconditional** pass; **(C)** a fourth rule in
`age_label_policy` that masks ΔAge for any line the census shows has no controls; **(D)** the B2′
assertion that no line both falls back *and* keeps its label. With the flag off, ΔAge must be
**bit-identical** to today.

---

## 1. The code as it actually is — what the spec has to fit

| fact | location | why it constrains us |
|---|---|---|
| `work = plan_all(sources)` | `build_dataset.py:379` | **unconditional** — the only pre-pass hook that always runs |
| `harmonizer = fit_harmonizer(cfg, work) if cfg.harmonize else None` | `build_dataset.py:383` | **conditional** — §9 wanted the census here; §10 showed it cannot live here |
| `fit_harmonizer` pools `controls.setdefault(str(ds), …)` | `build_dataset.py:360` | keyed on **`dataset_id`**, not `cell_line` — the tally does not exist even when it runs |
| `delta_age(...)` called **twice** | `build_dataset.py:189` (harmonized) and `:199` (raw) | rule 4 must reach **both** paths, not just the harmonized one |
| `_control_baseline` called **twice** | `aging.py:245` (`recenter_on_control_arrays`, S4) and `aging.py:301` (`delta_age`, S2) | **B2′ must cover both.** §9/§10 only ever discussed the `delta_age` one |
| `age_label_policy(n, source, obs, *, masked_datasets, clock_age_range)` | `aging.py:135` | pure; three rules; rule 4 must not break purity |
| `enforce_clock_age_range: bool = False` | `build_dataset.py:126` | the **ships-off precedent** to copy verbatim |

**The `recenter_on_control_arrays` site is the one nobody has mentioned.** It is S4, the
post-deconfounder re-centring, and it calls `_control_baseline` **without** a census and **without**
any mask. A B2′ that only guards `delta_age` would pass while S4 silently self-centres the same
orphaned line.

---

## 2. Component A — the integrity gate

**New file: `src/cellfate/data/integrity.py`.** Pure, no I/O, no config — so it is testable with
synthetic arrays and cannot be weakened by a data change.

```python
G1_LIBRARY_BAND = (1e5, 1e7)     # linear RPM must sum to ~1e6 BY DEFINITION of the units
G2_MIN_LOG2_RANGE = 8.0          # >= 256-fold spread between least- and most-expressed gene

def bulk_column_verdict(log2_values: np.ndarray) -> tuple[bool, str | None]:
    """Admit or reject ONE bulk sample from its raw log2-RPM column. Pure.

    Returns (admitted, reason). `reason` is None exactly when admitted, mirroring
    `age_label_policy`'s (mask, reasons) invariant so the two read the same way.
    Reasons: "library_out_of_band" | "dynamic_range_collapsed".
    """
```

**Verified before writing (§11):** on all 124 Gill columns this rejects **exactly 5**, with **0**
false positives of 119; G1 margins 2.58× / 1.69×, G2 gap 7.26 vs 9.00, no overlap. Each condition
independently rejects all five.

### 🟡 DECIDE A1 — where the gate is *called*

| option | effect | cost |
|---|---|---|
| **(i) in the bulk source's `fetch`** | the degenerate column never becomes a sample; `n_samples` drops by 5 | rejection is invisible to anything that does not read the source's log |
| **(ii) in `build_dataset`'s chunk loop** | central, one place, easy to census | the source has already materialised the row; needs the raw log2 values, which `fetch` has already transformed |

**My recommendation: (i)**, because the gate's claim is *"this column is not a transcriptome"* — it
should not become a cell in the first place. But (i) makes the rejection **per-source**, so a second
bulk source would need wiring too.

### 🟡 DECIDE A2 — does the gate apply to single-cell sources?

G1/G2 are **bulk-column** conditions. HFF is single-cell; a single cell legitimately has a small
library and a narrow range. **Recommendation: bulk only, asserted** — the gate must *refuse to run*
on a single-cell source rather than silently pass it.

---

## 3. Component B — the global control census

**Where: a new unconditional pre-pass in `build_dataset.py`, immediately after
`work = plan_all(sources)` (`:379`) and BEFORE `fit_harmonizer` (`:383`).**

```python
def control_census(cfg, work, gate_on: bool) -> dict[str, int]:
    """Per-`cell_line` count of ADMISSIBLE control samples, over the whole corpus. Global,
    not chunk-local -- Stage 1.5's Group E is the chunk-local case and must NOT trip rule 4.
    """
```

Two properties it must have, both from §10:

* **unconditional** — computed whether or not `cfg.harmonize` is set, so `harmonize=False` builds
  (arms B/C/D, any single-dataset build) get the same protection;
* **per `cell_line`** — not per `dataset_id`.

It must apply the **same gate** as component A, so a line whose only control is rejected reads
`0` here. That is the whole point.

### 🟡 DECIDE B1 — the cost of the extra pass

`fit_harmonizer` already does a full `src.fetch` sweep when `harmonize=True`. A separate census pass
means **a second full read** on those builds, and a **first** one on `harmonize=False` builds.

| option | cost |
|---|---|
| **(i) always a dedicated census pass** | one extra full read of the corpus per build |
| **(ii) census pass, and `fit_harmonizer` reuses its result** | one pass total when harmonized; needs `fit_harmonizer` refactored to accept pre-read controls |
| **(iii) census from `obs` only, if sources can yield `obs` without counts** | cheap — but **I have not verified any source can do this**, and inventing the capability is scope creep |

**Recommendation: (i) for correctness now, with (ii) as a follow-up if the build time hurts.** I
would rather ship a slow correct census than a fast one entangled with the harmonizer — that
entanglement is exactly what §10 rejected.

---

## 4. Component C — rule 4 in `age_label_policy`

```python
def age_label_policy(
    n, source, obs, *,
    masked_datasets=frozenset(),
    clock_age_range=None,
    lines_without_controls: frozenset[str] = frozenset(),   # <- NEW, rule 4
) -> tuple[np.ndarray, list[str | None]]:
```

New rule, in the docstring's own idiom:

> **4. `no_control_baseline`** — the cell's `cell_line` has **zero admissible controls in the whole
> corpus**, so it has no zero-point, so its ΔAge is **undefined**. Keyed on data integrity, not on
> identity. Empty by default, so nothing changes unless a caller opts in.

Passing a `frozenset` keeps the function **pure** — it takes the global fact as an argument rather
than reaching for it. Default-empty preserves B4.

### 🟡 DECIDE C1 — where rule 4 sits in the reason order

The docstring states reasons are checked *"in decreasing order of certainty"* and a cell reports the
**first** that applied, so the string is stable under reordering of the later rules.

| position | argument |
|---|---|
| **first, before `cancer_source`** | "no zero-point" is the most certain of all — ΔAge is not merely untrusted, it is undefined |
| **second, after `cancer_source`** | does not perturb any existing reason string |

**Recommendation: second.** The gain from being first is cosmetic; the cost is changing recorded
reasons for any cell that is both. **This must be stated in the docstring either way**, because the
docstring currently promises the first rule "is never weakened by the other two".

---

## 5. Component D — B2′, and it must guard **two** sites

**The amendment, to be recorded in C-7 §4 before the run (§10, §5b):**

> **B2 → B2′, amended 2026-08-08.** *Reason: the original conflated the **mechanism** (raise) with
> the **invariant** (no unmasked fallback label).*
> **B2′ — no `cell_line` may reach `_control_baseline`'s fallback AND retain its ΔAge label.**
> Formally: `assert not (fell_back and not masked)`.

`_control_baseline` (`aging.py:90`) must report whether it fell back — it already builds a `census`,
so the natural shape is an added per-line `"fell_back": bool`. Then:

| site | guard |
|---|---|
| `delta_age` (`aging.py:301`) | after computing `age_mask`, assert no line with `fell_back` has any unmasked cell |
| `recenter_on_control_arrays` (`aging.py:245`) | **the site nobody named.** It has no mask and no census. It needs either the mask threaded in, or an assertion that it is never called for a line with zero global controls |

### 🟡 DECIDE D1 — how to guard the S4 re-centring site

| option | |
|---|---|
| **(i) thread `age_mask` into `recenter_on_control_arrays`** | signature change; two callers to update |
| **(ii) assert at the call site in `build_dataset`** | keeps `aging.py` pure, but the guard lives away from the thing it guards |

**Recommendation: (i)** — the guard belongs with the fallback, and a signature change is cheap
compared with a silent self-centring.

---

## 6. Config surface

```python
# Stage 1.5.6 / C-7: reject bulk samples that are not transcriptomes (G1 library band,
# G2 dynamic range) and mask ΔAge for any cell_line left with no controls.
# OFF by default because turning it on MOVES LABELS -- N2's zero-point disappears and its
# 21 ΔAge labels are masked. Enabling it is a pre-registered change with its own bar.
bulk_integrity_gate: bool = False
```

One flag, not two: the gate and rule 4 **must not be independently switchable**, because the gate
alone creates the orphaned line B2′ forbids (§9's sequencing point, §11 item 3).

---

## 7. Test matrix

| bar | test | must fail when |
|---|---|---|
| **B1** | `bulk_column_verdict` over the 124 recorded Gill columns → exactly the 5, 0 of 119 | a threshold is widened to fit |
| **B3a** | synthetic **constant** column → rejected | the gate cannot reject |
| **B3b** | synthetic **sound** column → admitted | the gate rejects everything |
| **B3c** | synthetic corpus where a line loses its last control → rule 4 masks it | rule 4 is inert |
| **B3d** | synthetic corpus where controls exist but **not in this chunk** → rule 4 does **NOT** fire | rule 4 confuses Group E for Group D |
| **B2′a** | orphaned line **masked** → passes | — |
| **B2′b** | orphaned line **unmasked** → raises | the invariant is unenforced |
| **B2′c** | same, via `recenter_on_control_arrays` | the S4 site is unguarded |
| **B4** | flag off → `y_age` **bit-identical** to the current build, `max\|Δ\| == 0.00e+00` | anything leaks when disabled |

**B3d is the test that distinguishes this spec from a naive one** — it is the global-vs-chunk-local
distinction, and without it rule 4 fires on Stage 1.5's Group E and blocks C-7 for the wrong reason.

Every branch must **execute** in the suite (the `verify_1a` lesson): a branch that never runs is not
a check.

---

## 8. Order of work

1. `integrity.py` + B1/B3a/B3b — pure, no build, no data dependency
2. `control_census` + B3d — the global predicate, and the test that it is *not* chunk-local
3. rule 4 + B3c — `age_label_policy`, default-empty
4. `_control_baseline` fell-back reporting + B2′a/b/c — **both** sites
5. wire the flag; **B4** — bit-identical when off, on a rebuilt fold
6. **stop.** Adoption (flag on) is a separate pre-registered run with its own snapshot and a full
   guard re-report over **6** folds

Steps 1–4 touch `src/` but cannot move a label while the flag is off. **Step 5's B4 is the gate on
everything after it.**

---

## 9. What this spec does NOT do

* **It does not adopt C-7.** Ships off; enabling is a separate run.
* **It does not decide drop-vs-re-quantify.** Option (c) keeps the donor; SRA `SRP302546` remains an
  optional upgrade that would restore N2's ΔAge labels (§9).
* **It does not resolve the seven `mean − min` outliers**, including `Y1_Fib`, a control whose
  library and range are both normal. §5.10's *not established* stands and **Y1's floor ratio stays
  unexplained**.
* **It does not touch Stage 3a.** 3a stays blocked until C-7 is adopted and the labels are clean.

---

## 10. Open questions for review

| # | question | my recommendation |
|---|---|---|
| **A1** | gate called in the source's `fetch`, or in the build loop? | source `fetch` |
| **A2** | does the gate refuse to run on single-cell sources? | yes, assert |
| **B1** | dedicated census pass, reuse `fit_harmonizer`'s, or `obs`-only? | dedicated pass; optimise later |
| **C1** | rule 4 first or second in the reason order? | second |
| **D1** | guard S4 by threading the mask, or asserting at the call site? | thread the mask |

**None of these is a close call in my view, but each changes the code, so they are yours.**
