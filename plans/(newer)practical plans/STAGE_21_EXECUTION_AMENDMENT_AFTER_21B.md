# Stage 21 — EXECUTION AMENDMENT after 21B

**Status:** 🔒 PRE-REGISTERED AMENDMENT. Frozen before any public-dataset search begins.
**Additive.** `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md` and the 21A/21B results and records
are **not** rewritten. This document only splits the remaining execution.

---

## 1. Why this amendment exists

`STAGE_21_..._V2.md` §2 laid out three phases (21A local audit → 21B public qualification → 21C
download/reconstruct) on the assumption that the local geometry was still open. **21B resolved the
local questions instead**, and closed both local prospective routes:

```text
GSE242423 -> LINEAGE_ABSENT_PROVEN
GSE165176 -> ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY
```

So the phase that was going to be "public qualification" now needs its own number, and the
download/reconstruct phase moves after it. Renumbering keeps the record unambiguous — otherwise
there would be two different things called 21B.

## 2. The execution sequence, restated

```text
21A — Local Data Audit                       ✅ DONE   (f6a0056)
21B — Local Source/Design Resolution         ✅ DONE   (3b8b644)

21C — Public Prospective Dataset Qualification
21D — Acquire + Reconstruct Qualified Dataset(s)

21D PASS -> STAGE 22 — Prospective Benchmark Build
```

**Nothing in 21A/21B is reopened by this amendment.** Their verdicts stand as frozen.

---

# STAGE 21C — Public Prospective Dataset Qualification

**No modelling. No large reconstruction. No `src/` changes.** The only purpose is to decide which
public datasets genuinely satisfy the required prospective geometry.

## 3. The two roles

### ROLE A — core reprogramming prospective dataset (MANDATORY)

```text
X_before  +  reprogramming intervention U  ->  independently measured future outcome Y
```

Preferred: pre-treatment RNA · clone/lineage linkage · orthogonal future success/failure ·
enough independent biological units · reconstructable public files.

### ROLE B — multi-perturbation interaction dataset (HIGH VALUE, NOT BLOCKING)

```text
X_before  +  multiple treatments U1, U2, …  ->  independent future response Y
```

Ideally the same clone / sister population meets several treatments, so that Stage 25 can later
ask *does treatment preference depend on starting molecular state?*

**If Role B is not found inside the frozen budget, the paper is NOT held.** Scope the
treatment-ranking contribution down and proceed on a strong Role A.

## 4. Frozen search budget — exactly two passes

### PASS 1 — six named candidate families

| role | family | accession(s) |
|---|---|---|
| A | Rewind | `GSE227151` / `GSE243933` |
| A | CellTag | `GSE99915` |
| A | CellTag-multi | `GSE216518` / `GSE216521` |
| B | ReSisTrace | `GSE223003` |
| B | multi-treatment melanoma lineage | `GSE279162` |
| B | sequential-treatment / barcode | `GSE253739` |

Accessions and geometry are to be **verified against GEO / the source papers**, not inherited from
any summary in this repository — including summaries written by me.

### PASS 2 — limited expansion

Runs **only for a role still unresolved after Pass 1**. Maximum **3 additional serious candidates
per unresolved role**. A candidate counts toward Pass 2 only if its abstract/methods already make a
plausible case for pre-outcome molecular state **and** real biological linkage **and** a future
outcome.

### After Pass 2

```text
HARD STOP
```

No third pass. No open-ended "one more dataset" loop.

## 5. Cheap rejection cascade

Applied in order, so an obvious failure is found before anything large is downloaded.

```text
Q1. Is X measured before Y?                      NO -> reject
Q2. Real biological linkage early X -> later Y?  NO -> reject
Q3. Is Y independently measured, not f(X_future)? NO -> SURROGATE_ONLY
Q4. Can public files reconstruct the X->Y link?  NO -> MISSING_REQUIRED_FILE / reject as primary
Q5. Enough independent biological units?         NO -> PILOT_OR_REPLICATION_ONLY
Q6. Does treatment vary meaningfully?            required for ROLE B
```

**Do not download large raw data to discover a Q1/Q2/Q3 failure.**

## 6. Required qualification table

Per candidate:

```text
dataset · accession · paper · role_A_candidate · role_B_candidate
biological_system · species · cell_type
pre_outcome_rna_available · time_of_X
lineage_or_clone_link · linkage_method
future_outcome · outcome_measurement_method · outcome_orthogonal_to_rna
treatments · n_treatments · same_clone_or_related_population_across_treatments
n_independent_experiments · n_independent_clones_or_lineages · n_positive · n_negative
outer_split_unit · minimum_attainable_p · resolvable
processed_files_available · lineage_mapping_available · outcome_mapping_available · reconstructable
missing_files · evidence · verdict
```

Every important field keeps the Stage-21 standard — `value` / `status` / `evidence` — with the
status one of:

```text
PRESENT
ABSENT_PROVEN
UNKNOWN_REQUIRES_SOURCE_FILE
```

**"Not found" is never turned into "absent".** Carried forward unchanged from 21A/21B.

## 7. Qualification verdicts

```text
QUALIFIED_ROLE_A · QUALIFIED_ROLE_B · QUALIFIED_BOTH
QUALIFIED_SURROGATE · QUALIFIED_REPLICATION_ONLY
NOT_QUALIFIED · UNKNOWN_REQUIRES_SOURCE_FILE
```

**Role A qualifies only with all of:** pre-outcome X · valid biological linkage · orthogonal future
Y · reconstructable public data · resolvable independent units.

**Role B additionally requires** perturbation variation sufficient to later test `X × U`.

## 8. Stage 21C decision rule

```text
FULL_DATA_PATH   Role A qualified AND Role B qualified
CORE_DATA_PATH   Role A qualified, Role B not found inside the budget
                 -> proceed; treatment-ranking contribution is SCOPED DOWN
DATA_BLOCKED     Role A not found inside the budget
                 -> do NOT build the prospective architecture yet
```

If Role A qualifies and Role B does not, **do not keep searching**. The main paper continues as a
prospective reprogramming prediction paper; Stage 25 state-conditioned multi-treatment ranking
becomes optional.

## 9. Artifacts

```text
plans/(newer)practical plans/STAGE_21_EXECUTION_AMENDMENT_AFTER_21B.md   (this file)
experiments/diag_stage21c_public_qualification.py
tests/test_diag_stage21c_public_qualification.py
results/diag_stage21c_public_qualification_results.json
plans/(newer)practical plans/RECORDs/stage_21C_RECORD.md
```

The record stays compact — goal · candidates checked · source files/pages inspected · files
added/modified · bugs/corrections · final qualification table · final verdict · what this proves ·
what it does not prove · next action. **An audit record, not a diary.**

## 10. Engineering constraints

```text
src/ byte-unchanged          no predictive modelling      no target tuning
_s16 unchanged               no logistic regression       no outcome fishing
no architecture work         no neural network            no change to 21A/21B frozen verdicts
```

Full relevant tests and CI-scope ruff at the end.

## 11. After 21C

**Do not automatically download everything.** Report the qualification result first. Only once the
best Role-A (and, if available, Role-B) dataset(s) are frozen does **21D — Acquisition +
Reconstruction** open:

```text
download the selected data
reconstruct clone/lineage links
reconstruct future outcome labels
reproduce one source-study sanity result
build the frozen X_before + U -> Y_future table
run leakage guards
```

Only after 21D passes does **Stage 22 — Prospective Benchmark Build** begin.
