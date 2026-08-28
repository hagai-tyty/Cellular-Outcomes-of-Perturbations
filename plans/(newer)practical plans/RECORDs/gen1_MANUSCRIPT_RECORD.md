# gen1_MANUSCRIPT_RECORD — the last Generation-1 stage

## Goal
Execute §9 of the frozen ship plan — `MANUSCRIPT + REPRODUCIBILITY PACKAGE -> PREPRINT /
SUBMISSION` — under the six obligations the claim lock handed forward: re-verify the evidence lock,
carry every claim with its qualifiers, report `p_perm` as a floor, state that replication is
Generation 2, scan any new abstract-level sentence against the nine, and bind to both digests.

## What is different about this stage
Every previous stage checked something someone else wrote, or something a machine produced. This
one checks **prose I wrote myself** — the weakest position an instrument can be in, because the
author and the reviewer are the same process. The only defence is a checker that refuses
mechanically and has been shown to do so.

## Inputs
- `results/gen1_handoff_to_manuscript.json` — `GEN1_CLAIMS_LOCKED`
- evidence digest `455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9`
- claim digest `0b3c7f038a41c3b58b4c47cc5769d1e5e7be27e4a01d6c5f790a8e7da296ae5f` (was `23ea00b8...` when this stage first ran)

## Files added
- `plans/(newer)practical plans/GEN1_MANUSCRIPT_PACKAGE_V1.md`
- `experiments/run_gen1_manuscript.py`
- `tests/test_gen1_manuscript.py`
- `results/manuscript/MANUSCRIPT.md`, `REPRODUCIBILITY.md`
- `results/manuscript/manuscript_compliance.json`, `manuscript_controls.json`,
  `GEN1_MANUSCRIPT.json`, `GEN1_PACKAGE_DIGEST.json`

## What did NOT change
No locked artifact, no number, no analysis. Both locks verified clean before a sentence was
checked, and again after.

---

## Result

```text
  GEN1_MANUSCRIPT_READY

  package digest   7a467e7f02ea4cb6a669a1300d048369e39778371efe6d2bb3d9fafdde6b66c6
  evidence digest  455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9
  claim digest     0b3c7f038a41c3b58b4c47cc5769d1e5e7be27e4a01d6c5f790a8e7da296ae5f

  MS-A  both locks verify              6 checks
  MS-C  compliance                     10 checks, 0 forbidden hits
  MS-D  number traceability            17 numbers, 0 untraceable
  MS-E  the package                    7 checks, every documented command exists
  MS-F  the checker refuses            4 controls + 1 positive, all fire
```

### The checker was fired at four broken copies of the manuscript

```text
  plant a forbidden claim from the claim lock's corpus   -> caught
  drop the outcome-semantics qualifier                   -> caught
  quote p as `p = 0.000999` instead of `p < 0.001`       -> caught
  change delta_RANK by one digit                         -> caught
  the unmodified manuscript                              -> passes all four
```

All on in-memory copies; the file on disk is never touched, and a contract asserts that. The last
line matters as much as the first four: controls that fire on everything prove nothing.

### Every number traces to a locked artifact

Seventeen figures, each pinned to the words around it and matched against the JSON it came from —
`delta_RANK` and its interval, all three `R(W)` values, the null p95 and null max, the permutation
count, `delta_TOP1`, the eligible and excluded clone counts, the design width, the adversarial
refusal count, and Role A's power both as recorded (0.64) and as audited (0.45).

A manuscript is where a figure gets retyped, and a retyped figure is how a paper ends up
disagreeing with its own data.

### The manuscript states the uncomfortable things

Contracts assert these survive into the document rather than being quietly dropped:

```text
  abundance remains ~3.45x the whole state contribution
  Cisplatin is negligible on C1; Doxorubicin is negative on both endpoints
  Role A's gate 18.3 FAILED at 0.64, and our own audit put the true power lower, at 0.45
  the outcome is a detection proxy and is not death
  no independent biological replication has been performed
```

`Limitations` and `What this does not show` are separate sections on purpose: what the result
cannot support, and what may not be said, are different things, and a reader should not have to
infer either.

---

## Bugs found — one, mine, and it is the third of its kind

**A line wrap defeated a number-trace check.** The excluded-clone count is written
`... — 472 were` / `never detected under any condition`, and the pattern required `472 were never
detected` on one line. The number was correct; the checker could not see it.

This is the **third** time a hard wrap has broken a text check in this project:

```text
  Stage 23.2   the leak scan tripped on its own explanatory note
  claim lock   a newline treated as a clause boundary orphaned "not a / clinical
               recommendation" from its negation
  here         a wrap orphaned a number from the words that identify it
```

Fixed the same way each time should have been fixed: whitespace is normalised before the patterns
run, so the checker does not care where the author's editor broke the line. Three other patterns
were also too tight — a 20-space column gap against a 12-character allowance, a `\D` run that
could not cross `1,000`, and `audited` where the text says `audit`. All four were pattern defects,
not manuscript defects; no number was wrong.

## A second bug, found in the final close-out pass

**The test suite dirtied the working tree.** `test_the_controls_never_modify_the_manuscript` proved
its point by actually running the negative controls — and `negative_controls()` wrote
`manuscript_controls.json`, whose `runtime_seconds` differs on every run. So `pytest` left
`git status` non-clean, every time.

Nothing was wrong with the manuscript or the checks. But a test that modifies a committed results
artifact is a side effect, not a check, and in any setting that asserts a clean tree after tests it
would fail for no real reason. `negative_controls(write=False)` now exists for callers that only
want the result, and the contract additionally asserts that **the controls JSON is byte-identical
after the test runs**, so the same thing cannot return quietly.

Fixing it changed the executor and the contracts, both of which the package digest covers, so the
digest moved `e4df73af... -> 68a1fca2...`, then to `3a593709...` when the package document gained
the note below, then to `7a467e7f...` when a ruff `N802` fix in the claim-lock contracts moved
the claim digest this manuscript quotes. Earlier values are recorded rather than erased.

**A ruff error was also found in the close-out pass, and only one of three was fixed.** Five
lint findings sat in Gen-1 modules: `N802` on a helper in `tests/test_gen1_claim_lock.py`, and
four unused names across three executors. CI lints `src/ tests/ scripts/ plan_tests/`, so only
the first would have failed a build, and it is fixed. Three unused imports in
`run_gen1_manuscript.py` are fixed too, that file being covered by the package digest which
was moving anyway.

The remaining two -- an unused `sys` import in `run_gen1_evidence_lock.py` and a dead
`mixed_ok` variable in `run_stage26_scope_lock.py` -- are **left as they are, deliberately**.
Both files are locked evidence artifacts. Rewriting them would move the evidence digest that
the README, the manuscript, the claim-lock plan and four records all quote, and would
invalidate a lock for two cosmetic findings with no behavioural effect. `experiments/` is not
in CI's lint scope, by a decision recorded in the workflow itself. The findings are recorded
here instead: that is what a lock is for.

**And it is now enforced rather than merely fixed.** CI snapshots
`git status --porcelain --untracked-files=all` *before* the suite and diffs it after, failing the
build on any change and printing what moved. Snapshotting before, rather than asserting a clean
tree, means a line-ending or checkout quirk cannot fail the step — only something the suite itself
did. Verified in both directions locally: a real run leaves the snapshot identical, and appending
one byte to a committed artifact makes the comparison fire and name the file.

## Tests
- 21 manuscript contracts, 0 skipped
- both locks re-verify clean, before and after
- a full-suite run leaves the working tree clean, verified by hashing every file under `results/`
  before and after

---

## Scientific interpretation

**Proves:** the manuscript and package are consistent with everything locked beneath them. No
forbidden claim appears unnegated; every mandatory qualifier is present, and the abstract carries
system, vocabulary and outcome by itself because an abstract travels alone; every number traces to
a locked source; every documented command exists and accepts the flag it is given; and the checker
was shown to refuse four specific corruptions of the document it just approved.

**Does NOT prove:**
- **that the science is good.** `READY` is a consistency verdict. It says the paper agrees with its
  own data and stays inside its own ceiling. It says nothing about whether the question was worth
  asking, whether the design was the right one, or whether a reviewer will agree.
- **that the writing is clear or fair.** A document can mislead through emphasis, ordering, or what
  it chooses to make prominent, without tripping a single pattern.
- **that the checker is complete.** It catches phrasings it has patterns for. A claim nobody thought
  to forbid, or one carried by implication across two sentences, would pass.
- **that a reader can reproduce the result.** The package proves every documented command exists
  and every referenced path is real. Actually re-running the 10.7 h null on another machine is a
  separate claim nobody has tested.

## Generation 1 is complete

```text
  Stage 22   benchmark
  Stage 23   learnability and interaction gate
  Stage 23.2 Role A -- underpowered, gate FAILED, recorded as supporting only
  Stage 23.5 claim revision and the Stage-25 preregistration
  Stage 24   frozen tool, reproducible to 5e-16
  Stage 25   STAGE_25_RANKING_SUPPORTED -- the one load-bearing new result
  Stage 26   KNOWN_TREATMENT_ONLY_SCOPED_LIMIT
  evidence   GEN1_EVIDENCE_LOCKED, 54 artifacts
  claims     GEN1_CLAIMS_LOCKED, 3 allowed, 9 forbidden
  manuscript GEN1_MANUSCRIPT_READY
```

## Next action
`PREPRINT / SUBMISSION`. Nothing in Generation 1 remains open.

Generation 2 — independent new-system biological replication, unseen-condition transfer,
calibration, out-of-distribution validation — is future work and was never a Generation-1 gate. The
first of those is the one that matters: everything here rests on a single lineage-traced system,
and one system is one system however carefully it is evaluated.

Verify the whole chain at any time:

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
  python experiments/run_gen1_manuscript.py --verify
```
