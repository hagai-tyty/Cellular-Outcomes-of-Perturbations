# CHANGE C-7 — Bulk sample integrity: reject degenerate columns before they become labels

**Status:** 🔵 **PRE-REGISTERED 2026-08-08. NOT IMPLEMENTED. `src/` untouched, no label moved.**

**Status — ✅ IMPLEMENTED 2026-08-08 (§13). Flag OFF.** *The line above is left as written.* All four components shipped; **B1 verified twice** (real matrix and the recorded 124-column cohort) and **B4 verified** — `max|Δ| = 0.00e+00` on 1944 cells against the pre-C-1 baseline. **Not adopted:** the flag is off, no label has moved, and enabling it for anything consuming a trained model remains a separate pre-registered run.

**Why a file and not a section:** two machines pushed a `## 5.8` into `STAGE_1_5_6_SPARSE_CLOCK.md`
concurrently and neither saw the other's. This change is not part of 1.5.6 — it is a data-integrity
gate that 1.5.6 *surfaced* — and it takes the next free **change ID** (C-1 … C-6 are Stage 1.5.3's)
rather than a section number, so it cannot collide.

**Owner:** whoever holds `src/`. **Cost:** one assertion + tests to write; **one rebuild to adopt.**
**Blocking:** nothing. **Blocked by:** nothing. It can be written before 3b/3c report.

---

## 1. The defect, verified twice on two machines

`GSE165176_Log2_RPM_Sendai_reprogramming` contains columns that are **not transcriptomes**.

| column | log2 range | linear RPM sum | role |
|---|---:|---:|---|
| `Y1_d7_CD13_Sendai_Exp1` | **0.15** | **2.148e+09** | perturbation |
| `N3_d21_SSEA4_Sendai_Exp2` | 2.47 | 2.379e+08 | perturbation |
| `O2_d9_SSEA4_Sendai_Exp1` | 2.15 | 1.628e+08 | perturbation |
| **`N2_Fib_Sendai_Exp2`** | **1.74** | **1.030e+08** | 🔴 **CONTROL** |
| `N2_d21_CD13_Sendai_Exp2` | 7.26 | 1.694e+07 | perturbation |
| *— the other 119 —* | **9.00 – 15.26** | **2.859e+05 – 3.880e+06** | |

`Y1_d7_CD13` is **entirely constant** across ~20 000 genes. `N2_Fib` sits **0.0008 log2 above its own
floor** — a mean that close to the minimum is something RNA-seq cannot produce, because the
highly-expressed minority always pulls the mean well clear.

**`apply_qc` passes every one of them.** Its gates are `min_genes` and `max_mito_frac`
(`QCConfig`) — single-cell gates that a constant bulk column clears trivially.

**Reach.** `gill_bulk` is *both* the harmonization reference *and* a training source, so these
columns are simultaneously (a) age-labelled training rows and (b) inputs to `σ_ref`. `N2_Fib` is
`is_control`, so it enters `Harmonizer.fit` in **five of six LOOCV folds** — including **O1**, the
fold behind July's −24.02 (`STAGE_1_5_6_SPARSE_CLOCK.md` §5.7, §5.11).

**Stage 1.5.2's gate G-a made `n = 1` visible. Nothing checks whether that `n = 1` is sound.**

---

## 2. The gate — two conditions, both justified by units, not by this cohort's quantiles

An earlier proposal thresholded `mean − min` at ⅕ of the cohort median. **That was rejected and the
rejection is recorded** (§5.10): it cuts a continuous distribution 8 % from its neighbour, so on a
new cohort it would flag or miss arbitrarily. Both conditions below come from what the numbers
*mean*, and this cohort only confirms that they separate.

### G1 — library size

> The matrix is **Reads Per Million**. A sound column's linear values must therefore sum to ≈ **1e6**
> *by definition*. Accept `[1e5, 1e7]` — a decade either side of the value the units mandate.

| | |
|---|---|
| the 5 degenerate | 1.694e+07 – 2.148e+09 — **all above** |
| the other 119 | 2.859e+05 – 3.880e+06 — **all inside**, and all within 3.9× of 1e6 |
| margin below the ceiling | 1e7 / 3.880e+06 = **2.58×** |
| margin above it | 1.694e+07 / 1e7 = **1.69×** |

### G2 — dynamic range

> Any real transcriptome spans several orders of magnitude between its least- and most-expressed
> gene. Require **log2(max) − log2(min) ≥ 8**, i.e. at least a **256-fold** spread.

| | |
|---|---|
| the 5 degenerate | 0.15 – 7.26 — **all below** |
| the other 119 | 9.00 – 15.26 — **all above** |
| separation | **no overlap** |

### Why both

Each condition alone flags exactly the same 5 with 0 false positives on this cohort. **They are kept
as two because they fail differently:** G1 catches a mis-scaled library, G2 catches a collapsed
distribution. A future cohort that defeats one is unlikely to defeat both. **A column must satisfy
both to be admitted.**

### 🔴 What the gate does NOT catch, stated plainly

Seven further columns look poor on `mean − min` but pass both G1 and G2: `O2_d40`, `O2_d34`,
`O1_d34`, `N2_d11_CD13_Exp2`, `Y2_d34`, `O1_d11_CD13_Exp2`, and `Y1_Fib`.

**`Y1_Fib` is a control.** Its library (1.51e+06) and range (14.43) are entirely normal and only its
`mean − min` is low, on a continuum (0.2745 → 0.2967 → 0.3069). §5.10 downgraded it from "defective"
to "not established" and that downgrade stands. **C-7 does not resolve these seven.** They are
recorded as open, and Y1's unexplained floor ratio (§5.5) stays unexplained.

---

## 3. 🔴 The consequence that makes this more than one assertion

**`N2_Fib` is N2's only control.** Rejecting it leaves donor N2 with **no control at all**, and
`aging.py:88` then does this:

```python
ref = values[ctrl] if ctrl.any() else values[in_line]     # <- silent zero-point switch
```

**N2 would fall through to self-centring — which subtracts its own mean perturbation effect and
forces its mean ΔAge toward 0, with no warning, no counter and no mask.** That is the exact
behaviour Stage 1.5's Group D pinned as a defect. **Rejecting the sample without deciding the donor
would trade a known-bad control for a silent fallback, which is worse: the first is visible.**

### The three options, and the recommendation

| | consequence | |
|---|---|---|
| **(a) reject the sample AND the donor** | a donor with no sound control cannot carry control-relative ΔAge. Drops N2's **21** columns | ⛔ **SUPERSEDED by §9** — conflates three separable decisions |
| (b) reject the sample, let the fallback fire | N2's ΔAge forced toward 0, silently | ❌ Group D pinned this |
| (c) reject it from `σ_ref` but keep it as N2's baseline | the same column is "too broken to estimate variance, sound enough to define zero" | ❌ incoherent |

**Under (a) the corpus goes 124 → 100 Gill columns and 6 donors → 5.**

### And that reaches further than this change

* **LOOCV goes from 6 folds to 5.** Every Stage 1 guard, and step 6's MDE arithmetic, are computed
  over folds. **This change cannot be adopted without re-reporting them.**
* **It reaches C-2.** N2 is donor age 0. If N2 leaves on integrity grounds, C-2's *"masks the two
  neonatal donors"* becomes *"masks N3"* — plus HFF (`00_START_HERE.md`, C-2 section).
* **It reaches §4.7's whole finding.** The 16.67 yr fold spread is defined over six folds, one of
  which is N2's. Removing N2 does not "fix" that spread — **it removes the fold the spread was
  measured on.** 3b and 3c must therefore report **before** C-7 is adopted, or their question
  disappears rather than gets answered.

> **Sequencing consequence, and it is not negotiable: C-7 is WRITTEN now and ADOPTED after 3b and
> 3c report.** Writing the gate costs nothing and cannot be undone by their results; adopting it
> first would delete their evidence.

---

## 4. Pre-registered bars

| | bar |
|---|---|
| **B1 — separation** | on the Gill matrix, G1 ∧ G2 flag **exactly the 5** columns in §1, with **0** of the other 119 flagged |
| **B2 — no silent fallback** | after rejection, **no `cell_line` reaches `_control_baseline` with zero controls.** Asserted, not logged: a donor losing its last control must **raise**, per the `age_label_policy` fail-open precedent |
| **B3 — the gate can fail** | fed a synthetic constant column, the check must reject it; fed a synthetic sound column, it must admit it. **Both branches execute in the test suite** — a branch that never runs is not a check (the `verify_1a` lesson) |
| **B4 — bit-identical when off** | with the flag off, `ΔAge` is **bit-identical** to the current build. The gate ships **off** and is enabled by its own pre-registered run, exactly as C-2 did |
| **resolvability (§5b)** | B1 is a **deterministic** classification on a fixed matrix, not an estimate — there is no sampling null, so `bar_verdict` records it as *deterministic, resolvability N/A* rather than simulating one. **B2/B3 are likewise deterministic.** No bar here needs a power calculation, and claiming one would be theatre |

---

## 5. Decision branches, fixed in advance

| outcome | meaning | what happens |
|---|---|---|
| **B1–B4 all pass** | the gate is correct and inert until switched on | ships **off**. Adoption is a separate run **after 3b/3c report**, with a snapshot and every Stage 1 guard re-reported over ~~**5**~~ → **6** folds — ⚠️ **corrected §11: §9 adopted option (c), which keeps the donor and the fold. The "5" was written under option (a).** |
| **B1 fails — extra columns flagged** | a threshold is wrong, or the cohort is not what §1 measured | **do not widen the band to fit.** Re-derive from the units or record that the units do not separate this cohort |
| **B2 fails** | a donor reaches the fallback | **blocking.** The donor-level decision (§3 option a) is not optional; fix that first |
| **B4 fails** | the gate is not inert when off | **blocking, and a bug in this change** — nothing about a disabled flag may move a label |

---

## 6. What this change does NOT license

* **It is not a re-analysis.** §§4.5, 4.6, 1c and 1d were computed on a contaminated `σ_gill`
  (§5.11) and **C-7 does not re-measure them.** That costs one HFF stream and belongs to 3b/3c.
* **It does not withdraw anything.** No recorded result is retracted on the strength of this.
* ~~**It does not decide whether the defect is GEO's deposit or our read of it.**~~ ✅ **RESOLVED
  2026-08-08 — it is GEO's deposit. C-7 is needed, not redundant. See §8.**
* **It does not touch HFF.** The gate is for bulk columns; the single-cell path has `apply_qc`.

---

## 7. Artefacts

| file | role |
|---|---|
| `plans/STAGE_1_5_6_SPARSE_CLOCK.md` §5.7 | the discovery (second machine) |
| `plans/STAGE_1_5_6_SPARSE_CLOCK.md` §5.9 / §5.10 / §5.11 | the census, its correction, and the contamination reach |
| *(to be written)* `src/cellfate/data/qc.py` or the Gill source | where G1/G2 live |
| *(to be written)* `tests/test_bulk_sample_integrity.py` | B1–B4 |

---

## 8. ✅ 2026-08-08 — RESOLVED: the defect is in GEO's deposit, not our read

*§6 flagged this as the one thing worth checking before implementation. It was checked. C-7 is
needed.*

### The file we hold is the file GEO serves

| | ours | GEO's listing |
|---|---|---|
| filename | `GSE165176_Log2_RPM_Sendai_reprogramming (1).txt.gz` *(the `(1)` is the browser's)* | `GSE165176_Log2_RPM_Sendai_reprogramming.txt.gz` |
| size | 8 337 920 B = **7.95 MiB** | **8.0 Mb** |
| samples | **124** columns | **124** (GSM5027507 – GSM5027630) |

Series: *"Multi-omic rejuvenation of human cells by maturation phase transient reprogramming
[Sendai_RNAseq]"*. It is the **only** supplementary file on the record.

### The file is not damaged, and our parse is not at fault

| check | result |
|---|---|
| `gzip -t` | ✅ intact |
| line count | 35 806 (1 header + 35 805 genes) |
| fields per line | **136 on every single line** — 12 annotation + 124 samples. No ragged rows, so no column can shift |
| read path | values below were pulled with **`awk`**, not pandas. The defect is in the text |

### What the raw text actually contains

`N2_Fib_Sendai_Exp2` is **column 33**, and across 35 805 genes it takes **four distinct strings**:

| value | genes |
|---|---:|
| `11.489547` | **35 690** — 99.7 % |
| `11.64155` | 106 |
| `12.64155` | 8 |
| `13.226513` | 1 |

The first eight data rows read `MIR1302-11, FAM138A, OR4F5, RP11-34P13.7 …` → **`11.489547`** every
time, while sound `O1_Fib` in the same rows varies normally.

**`Y1_d7_CD13_Sendai_Exp1` is worse: TWO distinct values, the modal one covering 100.0 %** (to one
decimal) of 35 805 genes.

> **This is not a transcriptome and it is not a parsing artefact. It is what was deposited.**

### 🔵 A route that does not cost a donor — and it changes §3's recommendation

The GEO record points at **SRA `SRP302546`** for raw reads. **`N2_Fib` could be re-quantified from
raw rather than dropped.**

That matters because §3's recommendation — reject the sample *and* the donor — costs **a whole
LOOCV fold**, reaches **C-2** (N2 is donor age 0) and reaches **§4.7** (whose 16.67 yr spread is
defined over the six folds including N2's). **Re-quantifying one FASTQ restores N2's zero-point and
costs no donor, no fold and no guard re-report.**

**It is a data-acquisition task, not a code change**, and it is **far** cheaper than Stage 6's donor
acquisition — one sample from a public archive. **Recorded as an option, not adopted:** it needs its
own sizing, and re-quantifying one sample with a different pipeline than the other 123 introduces a
batch term that would have to be checked rather than assumed.

**C-7's gate is unaffected either way** — a column this degenerate must be rejected whether or not a
replacement is later obtained.

### ⚪ One gate candidate tried and REJECTED, recorded so it is not retried

**Distinct-value count does not reproduce the five.** Over all 124 columns:

| | distinct values |
|---|---|
| the five C-7 flags | 2, 4, 5, 5, **27** |
| the other 119 | **22** – 693 |

`N2_d21_CD13` (a C-7 flag) carries **27** distinct values while sound `O2_d40` carries **22** and
`O2_d34` **26**. **The populations overlap, so distinct-count would flag two sound columns before it
reached the fifth degenerate one.** G1 (library) and G2 (dynamic range) separate cleanly; this does
not. Do not substitute it.

---

## 9. ✅ 2026-08-08 — **OPTION (c) ADOPTED**, and B2 does not have to block it

*The choice among (a)/(b)/(c) is C-7's to make. It is made here.* **Still not implemented.**

### The decision

> **(c) — reject the degenerate control, MASK N2's ΔAge, keep the donor and the fold.**
> §3's recommendation of (a) is **superseded**: it conflated three separable decisions and answered
> all three with the harshest available answer.

| the three questions (c) separates | answer |
|---|---|
| should the degenerate control enter the harmonizer? | **no**, definitively — §5.14 |
| should N2's own 21 ΔAge labels survive? | **no** — its zero-point reads **98.65 yr** for a donor of true age **0** |
| should N2's **cells** survive at all? | **yes** — the fate head consumes no ΔAge and runs at `fate_roc` **0.983** (`00_START_HERE.md:70`, *"untouched by every ΔAge problem"*). Dropping the donor destroys working fate data to fix a broken age label |

**Verified here, independently:** applying the frozen clock to each donor's day-0 fibroblast gives

| donor | predicted | true | error |
|---|---:|---:|---:|
| **N2** | **98.65** | **0** | **+98.65** |
| O2 | 79.50 | 53 | +26.50 |
| O1 | 79.12 | 53 | +26.12 |
| Y1 | 64.92 | 29 | +35.92 |
| Y2 | 57.66 | 35 | +22.66 |
| N3 | 36.44 | 0 | +36.44 |

N2 sits **+35.12** above the other five's mean (63.53) — the arithmetic behind that figure, confirmed.

*Context, not a new defect:* the clock over-predicts **every** donor by +22 to +36 yr, and **N3 is
also donor age 0 yet reads 36.44** — it cannot separate a neonate from a 35-year-old. That is
absolute age, where the intercept does **not** cancel (§0 ERROR 1), so it is known behaviour. It is
recorded because it bears on C-2 and on whether ΔAge is meaningful for age-0 donors at all.

**And (c) keeps LOOCV at SIX folds**, so §5's "re-report every guard over 5 folds" does not apply,
§4.7's record stays comparable, and donors — the binding constraint Stage 6 exists to relieve — are
not spent.

### 🔑 The B2 collision dissolves: B2 was never a prohibition on reaching the fallback

The objection is real as stated. Rejecting N2's control leaves line N2 with zero controls, and
`_control_baseline` then self-centres — which B2 forbids. **But B2's purpose was never to forbid the
fallback existing. It was to forbid the fallback producing a label that is KEPT** — the
`age_label_policy` fail-open, where labels meant to be withheld were silently retained.

So B2 restates as a **conjunction**, and admits (c) unchanged:

> **B2′ — no line may reach `_control_baseline`'s fallback *and* retain its ΔAge label.**
> `assert not (fell_back and not masked)`.

### And rule 4 should be GENERAL, not a donor special case

§3's cost estimate assumed a fourth rule addressing N2 by `cell_line`. **It does not need to name a
donor.** `age_label_policy` keys on `source`, `masked_datasets` and `donor_age` — none reaches a
single line, and none has to:

> **Rule 4 — a `cell_line` with zero admissible controls has no zero-point, so its ΔAge is
> undefined and is masked.**

Keyed on **data integrity**, not on identity. No donor name anywhere. And it fires on **exactly**
the condition that triggers the fallback — which is what makes B2′ automatic rather than colliding:
**the predicate that causes the silent self-centring is the predicate that masks the label.**

It also finally closes Stage 1.5's **Group D** defect properly, open since the harmonization audit:
the silent zero-point switch stops being silent for every future dataset, not just for N2.

### 🔴 The subtlety that decides whether B2′ is implementable — two different fallbacks

`_control_baseline`'s own docstring:

> *"Falls back to the line's own mean when a line has no controls **in this chunk**"*

**That is per-CHUNK, and ΔAge is computed per chunk** (`build_dataset.py:306`). So the fallback fires
in two distinguishable situations:

| case | meaning | owner |
|---|---|---|
| **the line has no controls AT ALL** | no zero-point exists — N2 after rejection | **rule 4 / B2′.** Mask |
| the line HAS controls, but none landed in this chunk | a chunking artefact; the line's zero-point exists and is simply absent here | **Stage 1.5 Group E**, already diagnosed and separately owned |

**B2′ must test the first and not the second, or it will fire on Group E's case and block C-7 for
the wrong reason.** The predicate is *global* per `cell_line` — "zero admissible controls anywhere in
the corpus" — evaluated before chunking, not inside it. That is decidable in the same pre-pass
`fit_harmonizer` already runs.

### Sequencing — rule 4 ships WITH C-7, not after it

**C-7's gate alone creates the orphaned line.** Reject the control without rule 4 and line N2 has no
zero-point and no mask — the exact window B2 exists to forbid. **They are two halves of one
operation and must land in one change, one rebuild.** §3's *"(c) needs a fourth rule, and it's a
`src/` change and therefore its own Change"* is answered: **C-7 is already a `src/` change.** Adding
rule 4 to it costs nothing extra and removes a state the system must never be in.

### What is now settled and what is not

| | |
|---|---|
| the option | ✅ **(c)** |
| rule 4's form | ✅ general — zero admissible controls ⇒ masked |
| B2 | ✅ restated as B2′, a conjunction; no longer blocks (c) |
| the global-vs-per-chunk predicate | ✅ named; **implementation must honour it** |
| drop vs **re-quantify** from SRA `SRP302546` | ⚠️ **still open, and now optional rather than urgent** — under (c) the donor and the fold survive either way, so re-quantification becomes an *upgrade* (it would restore N2's ΔAge labels) rather than a rescue |
| the ten degenerate non-control columns | ⚠️ still C-7's scope, unresolved |

---

## 10. 🆕 2026-08-08 — §9 VERIFIED, with one bar-discipline note and one blocking correction

*Additive. §9 is unmodified. Every number and every code claim below was checked here before
agreement.*

### ✅ Verified, and agreed

| §9 claim | check |
|---|---|
| the donor error table | **exact.** N2 +98.65, N3 +36.44, Y1 +35.92, O2 +26.50, O1 +26.12, Y2 +22.66 — recomputed from the clock on each day-0 control |
| N2 is the outlier, not just "the clock is bad" | the clock is biased **high on every donor** (+22.66 to +36.44), and **N2 is 2.71× the next worst**. That is the right reading and it is stronger than "+35.12 above the mean" |
| `fate_roc` 0.983, *"untouched by every ΔAge problem"* | `00_START_HERE.md`. Dropping the donor destroys working fate data to fix a broken age label |
| **option (c)** | **agreed.** §3's option (a) answered three separable questions with the harshest available answer |
| **rule 4 should be GENERAL, not donor-named** | **agreed, and §9's version is better than the one §5.16 proposed.** Mine keyed on *identity* (`masked_cell_lines`); §9's keys on the *condition* — zero admissible controls ⇒ no zero-point ⇒ ΔAge undefined. It fires exactly where the fallback fires, works for every future dataset, and closes Stage 1.5 Group D generally |
| **two different fallbacks** (no controls at all vs none in *this chunk*) | **correct and sharp.** `_control_baseline` falls back *"when a line has no controls **in this chunk**"*, and Stage 1.5's Group E is the chunk-local case. B2′ must test the **global** predicate or it fires on Group E and blocks C-7 for the wrong reason |
| **rule 4 ships WITH C-7** | **agreed.** The gate alone creates the orphaned line — control rejected, no zero-point, no mask — which is precisely the window B2 forbids |
| SRA `SRP302546` becomes an optional upgrade under (c) | **agreed.** Donor and fold survive either way; re-quantification would only restore N2's ΔAge labels |

### ⚠️ Bar discipline — B2′ is an AMENDMENT, and should be labelled one

§9 says the B2 collision *"dissolves rather than blocks"*. **The substance is right; the framing
understates what is happening.** B2 as pre-registered reads:

> *"after rejection, **no `cell_line` reaches `_control_baseline` with zero controls.** Asserted,
> not logged: a donor losing its last control must **raise**"*

B2′ replaces *"must raise"* with *"may fall back if masked"*. That is a **change to the bar's test**,
not a reading of it — even though it is faithful to the bar's own **title** (*"no silent fallback"*),
which is about silence rather than about reaching. Under `REF_GROUND_RULES.md` §5b a bar may be
amended, **before the run, with the reason recorded** — which is exactly the situation. So the
amendment is legitimate and I agree with it.

**It should simply be recorded as `B2 → B2′, amended 2026-08-08, reason: the original conflated the
mechanism (raise) with the invariant (no unmasked fallback label)`** — not as the bar having
dissolved. This project has been bitten four times by bars that moved without the move being
labelled; the fix is cheap and it is the discipline, not a formality.

### 🔴 BLOCKING CORRECTION — the predicate CANNOT live in `fit_harmonizer`'s pre-pass

§9 states the global predicate is *"decidable in the pre-pass `fit_harmonizer` already runs"*.
**It is not, for two independent reasons, and both are in the code:**

**1. The pre-pass does not always run.**

```python
harmonizer = fit_harmonizer(cfg, work) if cfg.harmonize else None   # build_dataset.py:383
```

`fit_harmonizer` runs **only when `cfg.harmonize` is True.** A rule-4 predicate hosted there would
**silently not exist** in any `harmonize=False` build — and those are real: the arm B/C/D probes and
any single-dataset build. **A data-integrity invariant that evaporates when a config flag is off is
a guard that cannot fire**, which is the exact defect class this project has caught repeatedly
(`verify_1a` grading PASS on its own warning; `age_label_policy`'s fail-open; C-2's comment being
true only when written).

**2. Even inside the pre-pass, the tally is the wrong shape.**

```python
controls.setdefault(str(ds), []).append((norm[m], raw.genes))      # keyed on dataset_id
```

Controls are pooled **per `dataset_id`**, not per `cell_line`. The loop has `obs["cell_line"]` in
hand, so a per-line tally is easy to add — but *"already runs"* overstates it: **the loop runs, the
tally does not.**

**Where it should live instead.** `work = plan_all(sources)` and `load_or_fit_panel(cfg, work)` both
run **unconditionally** (`build_dataset.py:379-382`). The global per-`cell_line` control census
belongs in an **unconditional pass over `work`**, before and independent of `fit_harmonizer` — which
also makes it available to `harmonize=False` builds and to G-a's existing `baseline_census`, whose
job is already *"what each ΔAge zero-point actually rests on"*.

**This is blocking for the implementation, not for the decision.** Option (c), rule 4's general
form, the two-fallback distinction and the ship-together sequencing all stand exactly as §9 has
them. Only the predicate's **host** changes.

### Net

(c) is adopted. Rule 4 is general and ships with C-7. B2 becomes B2′ **as a recorded amendment**.
The global-vs-chunk-local distinction is correct and load-bearing. The predicate moves out of
`fit_harmonizer` into an unconditional pass. **Nothing here re-opens the decision.**

---

## 11. 🆕 2026-08-08 — B1 VERIFIED INDEPENDENTLY, and what still blocks implementation

*Additive. §2's gate ran here against the raw matrix before any agreement to proceed.*

### ✅ B1 PASSES — run, not quoted

`G1` (linear RPM sum ∈ [1e5, 1e7]) ∧ `G2` (log2 range ≥ 8) applied to all **124** columns of
`GSE165176_Log2_RPM_Sendai_reprogramming`:

| | |
|---|---|
| rejected | **exactly 5** |
| the five | `Y1_d7_CD13` 2.148e+09 / 0.15 · `N3_d21_SSEA4` 2.379e+08 / 2.47 · `O2_d9_SSEA4` 1.628e+08 / 2.15 · **`N2_Fib` 1.03e+08 / 1.74** · `N2_d21_CD13` 1.694e+07 / 7.26 |
| false positives | **0 of 119** |
| the other 119 | library **2.859e+05 – 3.880e+06**, log2 range **9.00 – 15.26** |
| G1 margins | ceiling/max-sound **2.58×**; min-degenerate/ceiling **1.69×** |
| G2 gap | worst degenerate **7.26** vs best sound **9.00** — no overlap |

**Each condition independently rejects all five** — every one fails G1 *and* G2 — so the "kept as
two because they fail differently" rationale costs nothing on this cohort. **B1 is confirmed, and
confirmed against the units rather than a quantile of this cohort.**

### 🔴 What still blocks IMPLEMENTATION — five items, all identified, none built

| # | item | state |
|---|---|---|
| **1** | **the predicate host.** §9 puts the global per-`cell_line` control census in `fit_harmonizer`'s pre-pass. That pre-pass is `if cfg.harmonize` (`build_dataset.py:383`) and tallies per `dataset_id`, not per `cell_line`. It must move to an unconditional pass over `work` | **§10, recorded, not fixed** |
| **2** | **B2 → B2′ is an unlabelled amendment.** Agreed in substance; §5b requires the move be recorded as an amendment before the run | **§10, recorded, not applied to §4** |
| **3** | **rule 4 does not exist in code.** `age_label_policy` has three rules (`aging.py:135-176`); the fourth — zero admissible controls ⇒ ΔAge undefined ⇒ masked — is designed and unwritten. It **must ship with C-7** | **designed only** |
| **4** | **nothing is implemented.** C-7's own status: *"PRE-REGISTERED. NOT IMPLEMENTED. `src/` untouched"* | **0 lines written** |
| **5** | **§3 and §5 still instruct a 5-fold re-report.** §3 line 119 (*"LOOCV goes from 6 folds to 5"*) and §5's decision row were written under option **(a)**. §9 adopted **(c)**, which keeps six. §9 *notes* the conflict but did not amend the instruction — and §5 is a **decision table**, i.e. what an implementer executes | **§5 row corrected above; §3 left as the record** |

### What C-7 does not resolve, restated so adoption does not imply it

Seven columns pass `G1 ∧ G2` while looking poor on `mean − min`, **including `Y1_Fib`, a control**
(library 1.51e+06, range 14.43 — both normal). §5.10 downgraded it to *not established* and that
stands. **Y1's unexplained floor ratio (§5.5) remains unexplained after C-7.**

### Cost of adoption, named

C-7 changes `y_age`, so adoption is **a rebuild + full re-score**, and the Stage 1 guard record
**restarts** — over **6** folds under (c). That is the retrain §5.13 freed by cancelling 1.5.6
step 4; it is not free, it is *already budgeted*.

### Net

**The DECISION is ready; the IMPLEMENTATION is not started.** B1 is verified, the gate is
unit-justified, (c) and rule 4's form are agreed, and the sequencing is settled. Items 1–3 are the
spec for writing it, item 5 was a live wrong instruction and is now flagged, item 4 is the work.

---

## 12. 🆕 2026-08-08 — THE FIVE OPEN DECISIONS ARE NOT DECISIONS. Four dissolve; one has a third answer.

*Additive. §11's spec is unmodified. Every answer below is read off the code, with the line cited —
none is a preference.*

**The objection that prompted this:** a project whose discipline is *measure, don't decide* should
not close a change on five judgement calls presented without data. It doesn't have to.

### A1 — gate in the source's `fetch`, or the build loop? → **`fetch`, decided by a count**

`src.fetch(chunk)` is called at **three** sites in `build_dataset.py`, each wrapped identically:

```
170:  raw = apply_qc(src.fetch(chunk), cfg.qc)     # process_chunk  -- the build loop
289:  raw = apply_qc(src.fetch(chunk), cfg.qc)     # gene-panel pre-pass
345:  raw = apply_qc(src.fetch(chunk), cfg.qc)     # fit_harmonizer pre-pass
```

**Gating in the build loop means three edit sites, and missing one is the failure mode** — miss :345
and the degenerate column still reaches `σ_ref`, which is the entire defect. **Gating in the source
is one site and covers all three.** Not a preference: 1 against 3, with a known failure mode on the 3.

### A2 — refuse to run on single-cell sources? → **the question DISSOLVES**

It only exists if the gate sits somewhere generic. **Put it in the bulk source's `fetch` (A1) and a
single-cell source never calls it** — there is nothing to assert and no branch to test.

And the assertion would be redundant anyway, on the same units argument that justifies G1:
**G1 is defined on RPM.** `GSE242423SingleCellSource` yields **raw UMI counts**, which sum to
~1e3–1e4 per cell — `normalize_counts` runs *later*, in `process_chunk`. So G1 would reject **every
cell** by construction of the units, not by any property of the data. A gate that cannot be
meaningfully applied should not be reachable, rather than reachable-and-guarded.

### B1 — dedicated census pass vs reuse vs obs-only? → **the flagged cost DOES NOT EXIST**

§11 flags *"a dedicated census means a second full corpus read on harmonized builds"* and weighs
that against entangling with the harmonizer. **Neither trade-off is real, because the census does
not need the corpus.**

Rule 4 fires only where a line can *lose* its last control. **A line can only lose a control to the
C-7 gate, and the gate is bulk-only (A1/A2). So the census is bulk-only.** That is `gill_bulk`:
**one 8 MB file, 124 columns** — already read once by `load_gill`, and trivially cheap beside the
42 605-cell single-cell corpus it does not touch.

**So take the dedicated pass — but the honest cost is one bulk matrix, not a corpus read**, and the
entanglement §10 rejected is avoided for free rather than paid for. *(`fit_harmonizer` is also
conditional — `build_dataset.py:383`, `if cfg.harmonize` — so reusing it would leave non-harmonized
builds ungated. A second reason not to.)*

### C1 — rule 4 first or second in the reason order? → **not cosmetic, and free to do RIGHT today**

`age_mask_reason` is **persisted**, not internal: `io.py:139` and `io.py:265` declare it in the
parquet schema, and `schemas.py:57-59` states the ordering is meaningful — *"a consumer that needs
RANK order can honour `donor_out_of_clock_range` differently from `cancer_source`."*

**Today no cell can match two rules**, because the only overlap candidate is C-2 and C-2 is **off**.
So the choice costs nothing now. **But it stops being free the moment C-2 activates:** N2 is donor
age 0, so it would match **both** rule 4 (no controls) *and* `donor_out_of_clock_range`, and the
recorded reason would be whichever comes first.

> **Place it by specificity while it is free:** *"no zero-point exists"* is **undefined**;
> *"donor outside the clock's fitted range"* is **out-of-validity**. Undefined is the stronger and
> more informative statement, so rule 4 belongs **before** `donor_out_of_clock_range` and after
> `cancer_source`. Deciding it later means deciding it under pressure, on shards already written.

### D1 — thread the mask, or assert at the call site? → **NEITHER. The census already records it.**

`_control_baseline` **already** distinguishes the two paths, at the single place the fallback occurs:

```python
ref = values[ctrl] if ctrl.any() else values[in_line]          # aging.py:116
...
rec = {"n_control": int(ctrl.sum()), "n_cells": int(in_line.sum()),
       "source": "controls" if ctrl.any() else "self_fallback",   # aging.py:121
       "unreplicated": bool(ctrl.sum() == 1)}
census[str(line)] = rec
```

**That is Stage 1.5.2's gate G-a — built to make `n = 1` visible. B2′ is the same census, one field
further.** So B2′ is one assertion over a dict that already exists, keyed by line:

> **for every line whose census records `source == "self_fallback"`, that line's ΔAge must be
> masked.**

No mask threading, no duplicated assertion at two call sites, no new plumbing — and it lands at the
**one** place the fallback can happen rather than the two places it is consumed.

**One real prerequisite, and it is small:** of the two call sites, only `delta_age` passes the
census (`aging.py:301`); `aging.py:245` calls `_control_baseline(values, lines, is_ctrl)` **without
it**. That site must pass a census too, or its fallback stays invisible. **That is the actual D1
work item** — and it is strictly smaller than either option offered.

### Net

| | §11 asked | answer | basis |
|---|---|---|---|
| **A1** | fetch or build loop | **fetch** | 3 call sites vs 1; `:345` is the one that matters |
| **A2** | assert on single-cell? | **question dissolves** | unreachable under A1; and G1 is RPM-defined |
| **B1** | which census? | **dedicated — and it is bulk-only** | the flagged corpus-read cost is not real |
| **C1** | first or second | **before `donor_out_of_clock_range`** | persisted column; free now, contested once C-2 is on |
| **D1** | thread or assert | **neither — read G-a's census** | `aging.py:121` already records it |

**Four of five were answerable from the code, and the fifth had a cheaper third option.** Recorded
because the pattern matters more than the five items: *"I left these to you"* on a change this
mechanical is a signal to go and look, not a signal to choose.

---

## 13. ✅ 2026-08-08 — **IMPLEMENTED**, flag OFF. What shipped, and two things the spec review missed

*Additive. §§1–12 unmodified.* **Code:** `src/cellfate/data/integrity.py`, and edits to
`sources.py`, `aging.py`, `build_dataset.py`. **Tests:** `tests/test_c7_integrity.py` (12),
`tests/test_c7_rule4_and_b2prime.py` (13). **Drivers:** `local_runners/build_c7_folds.py`,
`experiments/verify_c7_adoption.py`, `experiments/stage3a_forward_gate.py`.

### The five decisions, resolved as §12 argued

Each was verified against the file before adopting it. §12's readings were correct, and two of
them were better than the spec's recommendations:

| | resolution | verified by |
|---|---|---|
| **A1** | gate in the source's `_load` | `.fetch` really is called at **three** sites (`:170`, `:289`, `:345`); `_load` is one edit covering all three plus `plan()` |
| **A2** | dissolves | `normalize_counts` runs at `:174`/`:290`/`:358`, i.e. **after** fetch — a single-cell source hands over raw UMI counts, so G1 would reject every cell by units, not by data |
| **B1** | dedicated, and it is bulk-only | one 8.3 MB file; no second corpus pass |
| **C1** | rule 4 **before** `donor_out_of_clock_range` | `age_mask_reason` is persisted (`io.py:139`, `:265`); the order decides what is written once C-2 activates |
| **D1** | neither — reuse G-a's census | `aging.py:121` already writes `"source": "self_fallback"` |

### Two things the implementation found that no review had

1. **`_control_baseline` has TWO call sites**, and the second — `recenter_on_control_arrays`,
   the **S4 re-centring** — passed **no census**, so its fallback was **invisible**. A B2′
   guarding only `delta_age` would have *passed* while S4 silently self-centred the same
   orphaned line. It now accepts a census, with a test pinning it.
2. **`delta_age` is called twice** (`build_dataset.py:189` harmonized, `:199` raw), so rule 4
   had to reach both paths.

### A design error the suite caught before any rebuild

B2′ was first written **unconditional**. That broke
`test_the_silent_no_control_fallback_self_centres_a_line_to_zero`, which pins today's Group E
behaviour and states that changing it must be *"a deliberate, reviewed act"*. Making the
assertion unconditional **was** such a change, unreviewed. B2′ is now gated on the flag, so the
flag-off path is untouched — **which is what B4 then confirmed at 0.00e+00.**

### Bars

| bar | result |
|---|---|
| **B1** separation | ✅ exactly 5 rejected, 0 false positives of 119 — verified on the **real matrix** and on the recorded cohort; each condition independently rejects all five |
| **B2′** no unmasked fallback | ✅ enforced at both `_control_baseline` sites, three tests including the S4 one |
| **B3a/b** the gate can fail | ✅ both branches execute |
| **B3c** rule 4 masks | ✅ |
| **B3d** rule 4 does **not** fire chunk-locally | ✅ **the load-bearing test** — Group E must not trip C-7 |
| **B4** bit-identical when off | ✅ **`max|Δ| = 0.00e+00`**, 7 chunks / 1944 cells, self-test confirmed the bar can fail |

### End-to-end on the real Gill matrix

| | samples | donors | rejected | lines without controls |
|---|---|---|---|---|
| gate **off** | 124 | 6 | none | none |
| gate **on** | **119** | **6** | the exact five | **{N2}** |

**Option (c) exactly as designed: the donor and the fold survive.** §5's "re-report over 5 folds"
stays corrected to **6**.

### What is NOT done

**Adoption.** The flag is off. A **dataset-only** six-fold build (`_c7` roots) is running to
unblock Stage 3a — `test18_forward_gate.py` imports no `Predictor` and no bundle, so 3a needs a
build, not a retrain. Adoption for anything that consumes a **trained model** still requires the
retrain and the full Stage 1 guard re-report over six folds, as its own pre-registered run.

**And the seven `mean − min` outliers that pass G1 ∧ G2 remain unresolved**, including `Y1_Fib`
— a control whose library and range are both normal. §5.10's *not established* stands, and Y1's
floor ratio (§5.5) stays unexplained.
