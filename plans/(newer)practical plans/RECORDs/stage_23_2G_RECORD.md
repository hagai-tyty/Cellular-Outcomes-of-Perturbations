# stage_23_2G_ RECORD — step 1 of 5 (source qualification)

## Goal
Execute step 1 of the confirmation protocol's frozen execution order (V4 §19.0): qualify candidate
confirmation evidence against the §11 source criteria using declared metadata, author code and file
organisation only, and **fix the qualifying replicate set** before anything is reconstructed,
downloaded or fitted.

## Inputs
- `STAGE_23_2_ROLE_A_CONFIRMATION_V4.md` (live; V1–V3 archived unchanged)
- 23.2A reserved-confirmation ledger — declared GEO metadata only
- `GSE227151_family.xml`, the author Zenodo materials, and the frozen Stage-21C search record

## Files added
- `results/stage23_2/stage23_2g_qualification.json`

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2g-qualify`
- `tests/test_stage23_2_role_a_resolution.py` — 8 further contracts (110 total)

## What changed
- Nothing outside the additive qualification artifact

## What did NOT change
- **No reserved matrix downloaded. No outcome value read for any reserved candidate.**
- `src/` unchanged · no Stage-23 or earlier-substage artifact touched · no model fitted

## Tests
- 110 Stage-23.2 contracts passed, **0 skipped** · ruff clean

## Result

**`QUALIFYING_SET_EMPTY_FROM_FROZEN_SEARCH_SPACE`** — zero qualifying non-R1 biological replicates.
Gate 18.1 (`>= 2`) fails. Stage 24 remains **BLOCKED**. Runtime 0.01 min.

### The decisive criterion was the plainest one on the list

Not power, not geometry, not the ≥140 floor — **"independently measured later outcome"**.

```text
                    pre-state hiFT    iPS OUTCOME    declared outcome tables
biological rep 1          2                6                 (gDNA table)
biological rep 2          2                0                      none
biological rep 3          3                0                      none
```

Biological replicates 2 and 3 have pretreatment transcriptomes and **no later outcome measurement
of any kind**. Every one of the 13 GSMs in GSE227151 declares only 10X RNA files — barcodes,
features, matrix. **No GSM declares an outcome table**, and all six iPS outcome-side samples belong
to replicate 1.

The single outcome table that exists is not a GEO supplementary file at all: it comes from the
author's Zenodo materials, carries **one** selection unit (`SampleNum = 3`), and belongs to the
replicate Stage 23 already consumed.

A pre-state measurement without a later outcome cannot confirm a prospective claim, however good
the transcriptomes are. All five reserved candidates therefore fail §11 on four criteria at once:
later outcome, clone-lineage linkage for grouped evaluation, sufficiency to reconstruct the
endpoint, and establishable outcome-unit structure.

`GSM7092520` and `GSM7092521` fail on a fifth, independent ground as well — sorted for cycling
status, so a different population, failing the same-scientific-claim requirement. That was
anticipated in the 23.2A ledger and needed no new information.

### The frozen Stage-21C search space, reasoned over rather than re-searched

Stage 21C ran a pre-registered two-pass public search whose Pass-1 families are committed. Rather
than open an unbounded new search, each is dispositioned against the confirmation criteria:

```text
GSE227151   QUALIFIED_ROLE_A                   CONSUMED    R1's own series
GSE279162   QUALIFIED_ROLE_B                   CONSUMED    Stage-23 Role B
GSE99915    SURROGATE_ONLY                     DISQUALIFIED
GSE216518   SURROGATE_ONLY                     DISQUALIFIED
GSE216521   SURROGATE_ONLY                     DISQUALIFIED
GSE223003   RETAINED_PROSPECTIVE_REPLICATION   DISQUALIFIED_DIFFERENT_CLAIM
GSE253739   RETAINED_SECONDARY_SEQUENTIAL      DISQUALIFIED_DIFFERENT_CLAIM
```

The two *retained* candidates deserve their own note, because at first glance `GSE223003` reads like
exactly what is needed — Stage 21C retained it as a "prospective replication candidate". It is a
ReSisTrace drug-resistance lineage dataset. Its later outcome is **drug resistance**, not
reprogramming priming. It cannot confirm this Role-A hypothesis because it measures a different
biological event; retaining it for Role-B-style replication was correct, and reusing it here would
be a claim substitution rather than a confirmation. `GSE253739` fails for the same reason.

### The external search pass is authorised but unspent

```text
status      NOT_SPENT
authority   confirmation protocol V4 §12 -- "a single external search pass"
blocker     §12 authorises the pass but does not define its budget. Stage 21C's own
            Pass 2 was recorded NOT TRIGGERED and remains unspent.
```

Spending it is a scope decision rather than a mechanical step, so it was **not** taken unilaterally.
Every Stage-21C Pass-1 family is already dispositioned above; none yields a non-R1 Role-A
confirmation replicate.

### Status

```text
qualifying non-R1 biological replicates          0
gate 18.1 (>= 2 independent, none R1)            FAILS
BIOLOGICAL_REPLICATION_LIMITATION                SUPPORTED  (unchanged)
projected exit if the external pass is not
  spent, or is spent and fails                   ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
Stage 24                                         BLOCKED
```

## Bugs found
- **Eight contracts had begun silently skipping.** When V1–V3 were archived to `arcive/`, every
  contract keyed to a protocol version resolved to a non-existent path and its `skipif` guard turned
  it off. One test failed outright, which is how it was noticed — the other seven were skipping
  quietly, and a skipped contract protects nothing. Added a resolver that finds a version in either
  its live or archived location, pointed "current" at V4 so the V3 clarification contracts now
  verify the **live** protocol rather than a frozen ancestor, and verified the module reports
  **0 skipped**

## Scientific interpretation

**Proves:** the reserved confirmation evidence identified in 23.2A cannot confirm the corrected
Role-A hypothesis, and the reason is structural rather than statistical. Biological replicates 2
and 3 were never assayed for a later reprogramming outcome, so no reconstruction, cohort size or
analysis choice could make them usable. The frozen Stage-21C search space contains no other
candidate measuring this endpoint.

The 23.2F feasibility warning is therefore superseded by something stronger. That warning was about
*event count* — whether reserved replicates could reach 140 positives. The actual obstacle is that
they have **no positives at all**, because they have no outcome.

**Does NOT prove:**
- **That no confirming dataset exists anywhere.** The §12 external pass is unspent, and this
  qualification reasons over the frozen Stage-21C Pass-1 families rather than the whole literature.
- **That Role A is unconfirmable in principle.** It establishes that the *currently identified*
  evidence cannot do it.
- **Anything about the corrected hypothesis itself.** It remains untested on independent evidence,
  exactly as 23.2F left it.

## Next action
**STOPPED for review before spending the §12 external search pass.** Steps 2–5 of the V4 §19.0
order (per-unit reconstruction, floor and geometry checks, fold materialisation, model fitting) are
unreachable while the qualifying set is empty, and were not attempted.
