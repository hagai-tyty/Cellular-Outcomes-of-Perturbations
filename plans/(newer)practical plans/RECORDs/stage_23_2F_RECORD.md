# stage_23_2F_ RECORD

## Goal
Combine 23.2A–23.2E into a multi-label failure decomposition without forcing a single-cause story,
freeze exactly one mechanically-indicated corrected Role-A hypothesis, and write the confirmatory
protocol and Stage-24 handoff **before** any untouched evidence is inspected. **23.2F fits nothing.**

## Inputs
- 23.2A protocol `78edd5d7f9900349…` and the four substage artifacts, each pinned by SHA-256
- no raw data, no matrix, no model

## Files added
- `results/stage23_2/stage23_2_diagnostic_synthesis.json`
- `results/stage23_2/stage23_2_handoff_to_stage24.json`
- `plans/(newer)practical plans/STAGE_23_2_HANDOFF_TO_STAGE_24.md`
- `plans/(newer)practical plans/STAGE_23_2_ROLE_A_CONFIRMATION_V1.md`

## Files modified
- `experiments/run_stage23_2_role_a_resolution.py` — `--stage 23.2f`; 23.2A–23.2E untouched
- `tests/test_stage23_2_role_a_resolution.py` — 15 further contracts (81 total)
- `tests/test_stage23_learnability_gate.py` — two globs scoped so they stop claiming the Stage-23.2
  subdirectory as a Stage-23 artifact (see Bugs found). No Stage-23 assertion weakened

## What changed
- Nothing outside the additive 23.2F artifacts

## What did NOT change
- `src/` unchanged · no Stage-23 or earlier-substage artifact rewritten
- the historical `ROLE_A_SIGNAL_FAIL` is restated as permanent, not revisited

## Tests
- **2021 passed, 0 failed** across the whole suite · ruff clean (CI scope + both stage modules)
- 81 Stage-23.2 contracts, 15 new in 23.2F

## Result

### FROZEN 23.2F LEDGER

```text
MODEL_SELECTION_NULL_INFLATION              UNRESOLVED
RESIDUAL_DEPTH_STRUCTURE                    SUPPORTED
OUTCOME_LABEL_LIMITATION                    UNRESOLVED
WITHIN_R1_EVENT_COUNT_LIMITATION            SUPPORTED
BIOLOGICAL_REPLICATION_LIMITATION           SUPPORTED
ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA    UNRESOLVED

same-data status        MECHANISM_DIAGNOSED_ON_EXISTING_REWIND
projected 23.2G exit    AWAITING_23_2G_CONFIRMATION
stage_24_ready          false
```

**Two mechanisms are established, two are undetermined, and one is a design fact.** The decomposition
deliberately does not collapse to a single cause: the factor-interaction term in 23.2C was
indistinguishable from zero, so depth and selection act additively rather than one masking the other.

### Why "no biological signal" is not an available conclusion

`ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA = NOT_SUPPORTED` requires the label **and** event-count
limitations to be `NOT_SUPPORTED`. Neither is:

- `OUTCOME_LABEL_LIMITATION` cannot reach `NOT_SUPPORTED` without independent outcome-assay
  replication (V2 §8.6), which the Rewind materials do not contain;
- `WITHIN_R1_EVENT_COUNT_LIMITATION` is `SUPPORTED` — 29% power at 35 positives.

So the stage ends `UNRESOLVED` on the robust-signal axis. That is the designed behaviour, recorded
in V2 §14 in advance: **Stage 23.2 cannot conclude that Rewind has no signal, and a null result at
this event count is not evidence of absence.**

### Decomposition table

```text
                          B0            Bdepth
full K x C search     μ00 +0.00350    μ01 +0.00233
no-K-selection        μ10 +0.00299    μ11 +0.00161
53.9% of the historical null centre removed by the two corrections together

observed cells        O00 +0.01050    O10 +0.00962
                      O01 +0.00694    O11 +0.00872

corrected same-data   O11 +0.00872 vs q95 +0.00835, p_diag 0.0547  ->  NEGATIVE

technical retention   ρ ≈ 0.33-0.34 on all three pre-declared pairs
power ladder          35 pos 0.290 | 70 pos 0.520 | 140 pos 0.940   (oracle AUC 0.66)
```

### The corrected hypothesis — one, and mechanically derived

```text
status        ONE_CORRECTION_INDICATED
name          depth_complete_nuisance_control
```

> Under the historical Rewind outcome and the frozen Stage-22/23 evaluation geometry, pretreatment
> transcriptional state predicts the Role-A outcome beyond a **depth-complete** nuisance baseline
> `Bdepth = [log1p(n_pretreatment_cells), n_lanes, log1p(total_raw_GE_UMI),
> log1p(n_detected_GE_features_in_raw_pseudobulk)]`.

It is the only correction the ledger indicates, and the two alternatives are explicitly rejected
rather than silently dropped:

```text
search_matched_pipeline              rejected: MODEL_SELECTION_NULL_INFLATION is UNRESOLVED
alternative_outcome_representation   rejected: OUTCOME_LABEL_LIMITATION is UNRESOLVED,
                                     and V2 §8.7 forbids adopting one on that basis
```

**The correction was not chosen by effect size.** A contract re-derives the indicated set from the
ledger and fails if the frozen hypothesis disagrees, and further contracts assert that an
unindicated correction has not been smuggled into the hypothesis — the search grid must still read
"the frozen historical K × C grid, unchanged", and the outcome "unchanged".

### The benchmark-change firewall does not trigger

The correction alters only the model's **nuisance block**. That is not one of the material
benchmark semantics V2 §10.7 protects (outcome definition, positive/negative ontology, experimental
unit, source reconstruction, leakage firewall, split rule), so no benchmark versioning or Stage-22/23
re-gating is required. Had the label limitation been `SUPPORTED`, an alternative outcome would have
been indicated and the firewall *would* have fired — which is why the two are checked together.

### Minimum design requirement for confirmation

```text
 35 positive clones -> power 0.290
 70 positive clones -> power 0.520
140 positive clones -> power 0.940      <- smallest TESTED cohort reaching 0.80

requirement: >= 140 positive clones at oracle AUC 0.66
```

Read off tested rungs only. The ladder is coarse and V2 §9.5 forbids interpolating a precise
required N between them. A smaller confirmation cohort may be run, but a null from it is not
evidence against the hypothesis and must be reported as underpowered.

### Documents frozen

`STAGE_23_2_ROLE_A_CONFIRMATION_V1.md` — 14 sections covering hypothesis, outcome, allowed `X`,
nuisance block, model family and grid, grouping unit, metric, permutation design, PASS threshold,
positive-count floor, source qualification, search budget, forbidden already-inspected data, and
the Stage-27 firewall. It names `GSM7092515`/`GSM7092516` and the whole Stage-22 Rewind benchmark
as **consumed and forbidden**, and states plainly that whether reserved replicates 2/3 carry a
usable outcome is `UNVERIFIED`.

It also declares the null's known conservatism in advance: 23.2C measured ρ ≈ 0.34 donor–recipient
technical similarity surviving the permutation, so the confirmation retains that null unchanged and
puts the correction in the nuisance block — the null is not weakened to make confirmation easier.

`STAGE_23_2_HANDOFF_TO_STAGE_24.md` + JSON — Role A, Role B and Global blocks, marked
**`stage_24_ready: false`** with the reason. Role B is carried forward untouched, including the
limitations that must stay visible: Doxorubicin worse on both endpoints, Cisplatin's C1 gain
numerically negligible, abundance still ~3.45× the state contribution.

## Bugs found
None in the synthesis itself. Four contract failures surfaced, all of them integration faults I had
introduced earlier in Stage 23.2 and none affecting a result:

1. **A 23.2A contract was too blunt.** It asserted no reserved accession appears anywhere in the
   module source, and fired once 23.2F's confirmation protocol legitimately began naming
   `GSM7092517`–`GSM7092521` as *forbidden until qualified*. Scoped to data-access code: reserved
   accessions may appear in the documentation writer, nowhere else, and the ledger must still be
   derived from `family.xml`.
2. **My module broke the repo-wide `_RESULTS` convention.** `tests/test_results_paths.py` requires
   the literal form `_RESULTS = Path(__file__).resolve().parents[N] / "results"` so the constant is
   provably `__file__`-relative rather than CWD-relative. I had written `_RESULTS = ROOT /
   "results"`, which is equivalent at runtime but unverifiable by the contract. Adopted the
   convention.
3. **`results/stage23_2/` is caught by Stage-23's `stage23_*` globs.** Two Stage-23 tests —
   the large-artifact check and the determinism-coverage check — began claiming a *different
   stage's* subdirectory as a Stage-23 artifact. Both globs are now scoped to Stage-23's own files,
   excluding directories and the `stage23_2` prefix. Stage-23 never owned that subtree; the tests
   simply predated its existence.
4. **The CI-condition suite failed downstream of 2 and 3**, and went green once they were fixed.

Worth naming as a pattern: none of these were caught by the Stage-23.2 contracts, which all passed
throughout. They were caught only by running the *whole* suite, because the faults were in how this
stage sits alongside the rest of the repository rather than in anything it computes.

## Scientific interpretation

**Proves:** Stage 23's Role-A permutation failure has an identified, measured, partly methodological
cause. Residual technical depth structure surviving the abundance-preserving permutation is
established at the outcome level and mechanistically; together with model-selection breadth it
accounts for 53.9% of the positive null centre. Independently, the experiment was underpowered:
29% detection probability for an AUC-0.66 signal at 35 positives. Exactly one correction follows
mechanically from that ledger, and it has been frozen with a confirmation protocol that forbids the
evidence already consumed.

**Does NOT prove:**
- **That Role A has a signal.** The corrected same-data analysis is `NEGATIVE` (`p_diag = 0.0547`),
  and even a positive would have been `MECHANISM_DIAGNOSED_ON_EXISTING_REWIND`, never a
  confirmation.
- **That Rewind has no signal.** Structurally unavailable as a conclusion here, and correctly so.
- **That the correction is right.** It is the mechanically indicated one, not a validated one. Its
  entire status is "a hypothesis worth testing on evidence that has not been touched".
- **That confirmation evidence exists.** Reserved replicates 2 and 3 are declared metadata only.
  Whether they carry a reconstructable Role-A outcome is unverified, and 23.2G is where that gets
  established — or fails to.

## Next action
**23.2G — independent confirmation / roadmap resolution. NOT STARTED, and stopped here for review.**
Stage 24 remains blocked; `stage_24_ready` is `false`.
