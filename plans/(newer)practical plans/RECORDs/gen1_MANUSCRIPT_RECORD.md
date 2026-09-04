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
- evidence digest `e206bfd37c5a93998a773b8bd058eac5e5e144cd2a8ee5d78e9907911a956bc5`
- claim digest `a81ee43b07fae32f9bb45b4a4133de0b1f3979eeda3de7d700a7ca6897affb77` (was `23ea00b8...` when this stage first ran)

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

  package digest   1e514de4f59570d7a67ed15c4573a12bed6a2250410dc1b1ac58ab9dba732284
  evidence digest  e206bfd37c5a93998a773b8bd058eac5e5e144cd2a8ee5d78e9907911a956bc5
  claim digest     a81ee43b07fae32f9bb45b4a4133de0b1f3979eeda3de7d700a7ca6897affb77

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

## A literature check, and the citation it forced

Before anything shipped, the claims were checked against the published record rather than against
our own artifacts alone. It found the most consequential omission in the manuscript.

**The manuscript cited nothing at all — including the paper that generated its primary dataset.**
GSE279162 is the data of Schaff DL, White PE, Cote CJ, Watterson GE, Lin KZ, Fasse AJ, Zhang NR and
Shaffer SM, *Pre-existing cell states predict resistance to multiple treatments*, Cell Genomics
6(6):101191, 2026, doi:10.1016/j.xgen.2026.101191 (PMID 41916275). Their design — a barcode library
into WM989 A6-G3, 350,000 uniquely barcoded cells, split across dabrafenib, trametinib, CoCl2,
acidic media, cisplatin and doxorubicin — is exactly the six conditions reanalysed here. Reanalysing
someone's dataset without citing them is not a stylistic lapse. Both accessions are now cited and
the Data section states plainly that we generated no new data.

**The novelty framing was also wrong, and is corrected.** The Introduction read as though prospective
clone-level evaluation were the gap this work fills. It is not: Schaff et al. built the prospective
system and already showed that pre-existing state predicts which clones resist, identifying CD44 in
treatment-naive cells as a marker of resistance across multiple conditions.

**What the check did NOT do is weaken the result — it sharpened it.** Their finding is a general
propensity: some clones survive many things. That is a state *main effect*, and it is exactly what
W4 contains. Two facts already in the frozen design separate our claim from theirs:

```text
  R(W4) 0.692176 sits BELOW R(W1) 0.692654 -- the additive/general state term adds
  nothing to ordering, so the entire gain is the interaction

  within-clone AUROC compares the six scores of ONE clone, so any quantity acting on
  that clone as a whole shifts all six equally and cannot change their order. A purely
  clone-level propensity signal contributes exactly zero, by construction
```

Both are now stated in Results rather than left implicit. The result is complementary to theirs,
not a rediscovery of it: they establish that state predicts *how resistant*, this establishes that
state also carries information about *which condition* — the part a general axis cannot supply.

## A correction I nearly introduced myself

Checking the older ΔAge claims, the `diag_clock_circularity` artifact appeared to show the blanket
word "circular" overstating a mixed result — one arm reading `NOT CIRCULAR`. I began softening the
README and ARCHITECTURE on that basis.

**That reading was wrong: I had seen only the first of two arm sets.** Under **C-7 all five arms
verdict CIRCULAR** (ridge-vs-label ρ 0.965–0.995). The dissenting arm, N3, belongs to the earlier
`pre-C-7` set and is itself CIRCULAR under C-7. The original claim was better supported than my
correction to it, and the softening was reverted before it shipped. Both documents now carry the
per-arm detail so the question cannot be re-litigated from memory.

The claim itself needs no literature support: the clock is an elastic-net linear model on
log1p-CP10K expression, 1,956 of the 2,000 panel genes carry clock weights, so predicting ΔAge from
that same expression recovers a linear functional of the input. That is arithmetic.

One wording was genuinely too strong and is fixed: the **7.30 yr** figure is the disagreement
between the two reference methylation clocks *on our samples*, not a published constant. The
published clocks report ~3.6 yr (Horvath) and ~3.9 yr (Hannum) MAE against chronological age in
their own cohorts. Calling 7.30 an "instrument floor" without that qualifier invited a reader to
take it as a field constant.

## Why CI was red on every commit

Every workflow run since the manuscript stage opened had failed, and the cause was **one test of
mine that could only ever pass on the machine that wrote it**.

`test_every_evidence_lock_input_is_a_real_path` asserted that every path the handoff names exists on
disk. One of them is `stage24_w5_artifact.npz` — 44 MB, gitignored, **deliberately absent from a
fresh clone**, rebuilt by `--stage 24c`. That gap is recorded in the plan, in the evidence lock, in
the manuscript and in the reproducibility package. The test contradicted all four.

It passed locally because the artifact was sitting in my working tree. Reproduced by hiding the
file and re-running: exactly one failure across the whole suite, and it was this.

The exemption is now **derived from the evidence manifest's own `git_ignored` list** rather than
hardcoded, and the manifest separately gates that list down to exactly one entry — so a second
unbuildable path still fails. Verified in both directions: the suite is green with the artifact
present and with it absent, and the tree-check step passes in both.

**Lesson, and it is not a new one here.** A check that has only ever run in one environment has
only ever been tested in one environment. The locks are verified from a fresh clone by design;
their own contracts were not.

### The cascade this forced

The fixed file is `tests/test_stage26_scope_lock.py`, which the evidence lock covers under `code` —
so the lock refused, correctly, and the whole chain moved:

```text
  evidence  455892ff -> 2edc73c5     the fixed contract file
  claim     0b3c7f03 -> 0453a1af     the claim plan quotes the evidence digest
  package   6211abc5 -> ba2c9989     the manuscript quotes both
```

Nothing scientific changed: no number, no claim, no qualifier, no forbidden entry. Every earlier
digest is recorded above rather than erased, and a sweep confirms no document quotes a superseded
one.

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

---

## Ground-up verification — 2026-08-28

A full audit before submission, against the published reporting standards rather than against my
own checklist. The strongest check is the first one, because it does not trust this project's code
at all.

### The headline result was re-derived independently

A separate implementation, importing **no project module**, reading only the locked out-of-fold
table, recomputing eligibility, within-clone AUROC and the mean from first principles:

```text
  clones 1401   eligible 892   never detected 472   always detected 37    -- all match

  R(W1)       mine 0.692653836572   recorded 0.692653836572   |diff| 3.3e-16
  R(W4)       mine 0.692175822123   recorded 0.692175822123   |diff| 1.1e-16
  R(W5)       mine 0.743781141006   recorded 0.743781141006   |diff| 1.1e-16
  delta_RANK  mine 0.051605318884   recorded 0.051605318884   |diff| 2.2e-16

  p_perm from the exported draws: 0 of 1000 at-or-above -> 0.000999001, matches
  bootstrap CI from my own per-clone values: EXACT, 0.00e+00 on both endpoints
```

### The shipped artifact regenerates every frozen prediction

All **8,406** prediction cells over all 1,401 clones, each scored by the fold component that did not
train on it, against the frozen `pred_W5` column: max |diff| **6.66e-16**. Stage 24C recorded
4.996e-16 over the same table computed through its own path; mine is marginally larger because it
goes through the public API, whose design matrix has a different row count. Both are four orders
inside the frozen 1e-12 bound, and the difference between them is not a discrepancy in substance.

### Two misstatements found in the manuscript, both about someone else's experiment

**1. The clone count misdescribed the source experiment.** The abstract said *"1,401 barcoded clones
were split and exposed to six conditions"*. That is not what happened: 350,000 barcoded cells were
isolated and the experiment recovered many thousands of clones. **1,401 is our analysable subset** —
the clones carrying a pretreatment observation, which is what a prospective question requires.
Stating our subset as the experiment's size misrepresents work that is not ours. Corrected in the
abstract and in Data.

**2. The treatment schedule was asserted more precisely than the sources support.** We wrote that
cisplatin and doxorubicin ran *"two weeks followed by a two-week holiday"*. The GEO summary says
2 + 2 for both; the paper's methods give doxorubicin as 2.5 + 1.5. Two sources disagree and the
detail is not load-bearing here, so the manuscript no longer asserts a split: four weeks per arm,
treat-then-recover for the two chemotherapies, exact schedules cited to [1]. Figure 1 carried the
same error and was regenerated.

Neither error touched a number in the result. Both were in prose describing the source experiment,
which is exactly where a reanalysis has the least excuse to be sloppy.

### Everything else checked

```text
  figures            9 locked values present; every 6-dp number drawn traces to the verdict
  tool refusals      Vemurafenib / acid / Carboplatin / "" -> UNSUPPORTED_TREATMENT, no score
                     missing nuisance -> MISSING_REQUIRED_NUISANCE
                     ranking_status without the verdict file -> NOT_SUPPORTED
  digests            no stale 64-hex value in any live document
  bundle             384 files, BUNDLE_INTACT
  locks              evidence e206bfd37c5a93998a773b8bd058eac5e5e144cd2a8ee5d78e9907911a956bc5
                     claim    a81ee43b07fae32f9bb45b4a4133de0b1f3979eeda3de7d700a7ca6897affb77
                     package  1e514de4f59570d7a67ed15c4573a12bed6a2250410dc1b1ac58ab9dba732284
```

The evidence lock refused mid-audit when Figure 1's generator was corrected — a locked artifact
changed and the chain would not proceed until it was re-locked. That is the machinery working, not
a fault.

---

## Release-preparation corrections — 2026-08-28

Four problems raised in external review, all verified before acting on, all real. Two were
defects rather than tidying. Recorded here because none of them was written down when it was
fixed — the corrections landed while a summary of them was interrupted, and an unrecorded fix is
the failure mode this project's records exist to prevent.

### 1. The published step order would have invalidated the locks

`export_gen1_source_data.py` writes **four locked files**: the two per-draw CSVs and
`environment_lock.txt` (evidence lock), and `figure_source_data.json` (package digest). The release
checklist had it running *after* the three locks, which would have invalidated the digests just
computed. Corrected to run before. The lock demonstrated the point unprompted by refusing mid-work
for exactly this reason.

### 2. The bundle could not verify after unpacking

`environment_lock.txt` is evidence-locked and was **absent from the archive**, so `--verify` on an
unpacked bundle would have reported it missing — defeating the only thing the archive exists to
allow. `LICENSE` and `requirements.txt` were missing too.

`--check` was also weaker than it read: it compared working-tree files and ZIP **filenames**, never
the bytes inside the archive, and printed the bundle's own SHA-256 without comparing it to
anything. It now hashes every member from inside the zip, pins the zip against its recorded
checksum, and records the git commit.

A build-time guard now **refuses** if any artifact hashed by any lock is absent from the archive.
Negative control: run against the previous member set it names `environment_lock.txt` — it would
have caught the original bug.

### 3. The pseudobulk cache was optional and must not be

The documented rebuild path `--stage 24c` fails without it (`23A pseudobulk cache missing`), so an
archive lacking it can regenerate nothing and cannot re-run the Stage-25 null. Now required; the
build refuses without it.

### 4. Wording, and the novelty claim

```text
  README   "which conditions it survives"  ->  "under which conditions it remains detected"
  README   "reproducible bit-for-bit"      ->  hash-verification of locked artifacts (holds
                                               on any machine) distinguished from REFITTING
                                               (depends on BLAS, threading, library versions;
                                               not claimed)
```

And the prior-work paragraph was **narrowed**, which matters more than the wording. Schaff et al.'s
deposited analysis does report condition-associated markers and signatures. The manuscript no longer
implies the literature establishes only a general propensity. The contribution is now stated as what
it is — a frozen, clone-held-out, preregistered test of clone-specific *ordering* with abundance held
fixed — not priority on condition-specific expression analysis.

## The licence was contradictory, and unsound as written

Three files disagreed: `pyproject.toml` and `CITATION.cff` said MIT while `LICENSE` declared
"GPLv3 for academic and non-commercial use" plus a required commercial licence.

Beyond the contradiction, that clause **cannot be offered**: GPLv3 §7 forbids adding field-of-use
restrictions, so "GPLv3, non-commercial only" is not a licence GPLv3 permits anyone to grant. This
was surfaced rather than resolved unilaterally, because the choice is the copyright holder's.

Resolved to **PolyForm Noncommercial License 1.0.0** (SPDX `PolyForm-Noncommercial-1.0.0`), the
official text unmodified, with a separate `COMMERCIAL-LICENSING.md`. It expresses the intended
structure precisely, and its *Noncommercial Organizations* clause covers educational institutions
and public research organisations **regardless of source of funding** — so grant money does not make
academic work commercial.

```text
  LICENSE  pyproject.toml  CITATION.cff  .zenodo.json  README  SUBMISSION.md  MANUSCRIPT.md
  all now say PolyForm-Noncommercial-1.0.0; no "MIT" remains anywhere
```

Citation is requested strongly in the README and the commercial notice but is deliberately **not**
in the licence grant, so it never becomes a condition on running the software. This is
source-available, not OSI open source. It satisfies BMC Bioinformatics' requirement that software be
freely available to non-commercial researchers with no gate, and the submission pack tells the
editor so directly rather than letting it surface in review.

## Three gaps this record pass itself found

Writing the above surfaced what the interruption had cost:

```text
  MANUSCRIPT.md         had NO software-availability or licence statement at all, which a
                        BMC submission requires
  REPRODUCIBILITY.md    said nothing about terms, so someone unpacking the bundle had no
                        statement of what they may do with it
  SUBMISSION.md         the pre-flight checklist had no licence-coherence step
```

All three added. The lesson is the one this project keeps relearning: the work is not finished when
the code is correct, it is finished when the record says what happened.

---

## The claim lock changed its own digest on every run — 2026-08-28

The most serious defect found in the release-preparation work, and it was inside the lock system
rather than in anything the locks were guarding.

`runtime_seconds` sat inside `GEN1_CLAIM_LOCK.json`, a file the claim digest **covers**. Re-running
`--stage all` with nothing substantively changed therefore produced a different digest every time:

```text
  run 1  528849f2...      run 2  9666ca97...      run 3  fef69804...
```

**This was not cosmetic.** The claim digest is quoted in the manuscript, the README and the
submission pack as the identity of what may be said. It was valid for exactly one execution.
Anyone reproducing the stage would have re-run it, got a different number, and concluded the claim
set had moved. A lock that cannot survive its own stage being re-run is not a lock.

Found by running the stage three times and watching the number change — not by any check, because
no check asserted the property.

### The fix

The digest now hashes **content, not timing**: JSON members are normalised with volatile fields
removed and keys sorted before hashing; other text keeps the canonical-LF rule. Three consecutive
runs now yield one digest, and all three locks were verified stable across a re-run.

The **evidence manifest is deliberately not changed**. Its hashes must equal the bytes git stores,
which the clone-portability check depends on; introducing JSON normalisation there would break that
correspondence. The evidence lock covers frozen upstream outputs that are never regenerated, so it
was never exposed to this.

### What is still true, and stated rather than hidden

Re-running the three stages still rewrites **10 files on disk** — purely the timing fields. Verified
mechanically: stripping volatile keys, every one of those files is byte-identical in content to its
committed version. The digests are unaffected, so verification is unaffected, but a reproducer who
re-runs the pipeline will see a dirty working tree and should know it is benign.

The release bundle is a snapshot of one tree, and re-running a stage invalidates its checksums. That
is correct behaviour and the bundle checker catches it: during this audit the determinism test
itself produced a `BUNDLE_MISMATCH`, which is the checker doing its job on me.

### Two contracts so the class cannot return

```text
  stable_sha256 must ignore timing and key order, and must still catch a real content change
  no file any self-digest covers may carry an unstripped volatile field
```

Neither re-runs a stage, because doing so would dirty the tree that CI checks — the fix for one
problem must not manufacture another.

### Why this reached a record late

The fix was committed with a full explanation in the commit message and did not reach a record until
a later verification pass looked for it. That is the second time in this sequence that a correction
landed without being recorded, and the reason is the same both times: the work felt finished when
the code was correct. It is finished when the record says what happened.

---

## Tree stability, and a verifier that was not verifying — 2026-09-03

The section above ends by stating, honestly, that re-running the three stages still rewrote **10
files on disk**, that the churn was purely timing, and that a reproducer should know it was benign.
That statement was true when written. It is no longer true: the churn is now **zero files**, and
chasing the last of it turned up something considerably worse than churn.

### The churn had a second cause, and it was not timing

Stripping volatile timing at the three lock writers removed most of it. Sixteen files still moved.
The residue was one field: the evidence manifest recorded each artifact's size as `st_size`.

`st_size` counts a CRLF as two bytes. `.gitattributes` declares `results/** text eol=lf`, so a file
written by a stage on Windows lands CRLF while the same file in a fresh checkout is LF — identical
content, two different sizes. The manifest therefore disagreed with itself across checkouts. This is
the same class as the raw-versus-canonical hashing bug this lock already fixed once for hashes; the
size field was simply never brought along.

Measured on a clean clone: **34 of 62 artifacts** had an `st_size` different from the recorded size.
Sizes are now taken from the same canonical-LF content that is hashed. Hashes never moved.

Checked and found sound, rather than assumed: Stage 22's CSV writers emit LF, so the recorded sizes
Stage 23 compares against are portable; Stage 26 measures `len(text.encode())` after a universal-
newline read, which is already LF. Neither needed changing.

### The serious one: `--verify` never looked at the verdict

While cascading, the manuscript stage recorded `GEN1_MANUSCRIPT_REFUSED` — `MANUSCRIPT.md` still
quoted the previous claim digest, and the compliance check "both digests are quoted" was `False`.
The gate worked. `--verify` returned **clean anyway**, exit 0.

It did so because all three verifiers asked only one question: *have the covered bytes moved?* A
refused run still writes its outputs, so the bytes matched what was recorded, and the verifier
called that intact. CI runs exactly these three commands. **CI would have stayed green over a
refused package**, and so would the verification a reviewer runs on the unpacked Zenodo bundle.

Byte-intactness is not the same as a passing stage. Each verifier now reads the recorded verdict and
reports `EVIDENCE_STAGE_REFUSED` / `CLAIMS_STAGE_REFUSED` / `PACKAGE_STAGE_REFUSED`, distinct from
`*_MOVED`, with a non-zero exit. A guard plants a refusal in each of the three records and requires
each verifier to catch it, paired with a positive control so it cannot be satisfied by a verifier
that always refuses.

A refused run also printed `"next": "PREPRINT / SUBMISSION. Generation 1 is complete."` beside its
own `REFUSED` verdict. All three now drop that line when they refuse.

### A count that went stale while the check stayed correct

CL-A's label read `all 54 locked artifacts still verify`. The lock covers 62. The verification was
never wrong — it checks `live["clean"]`, not a number — but the label describing it was, and it had
been shipped that way. The label is now derived from the manifest, and its test cross-checks the
number against what was actually verified.

### What was changed in the plan, and what was left alone

`GEN1_CLAIM_LOCK_V1.md` pins the evidence digest it was written against, so a cascade moves it. The
new digest is written in, and **Amendment V1.1 records what V1 said** — the original digest and the
original count of 54 — so the change is visible rather than silent. No gate was relaxed.

### Final state

```text
  evidence  79f96dd66fe7a586df585126bc920d9433746e78ca25f7c453e2ed0ab110dfdb
  claim     5aab50db44fbbfeb655325e1d329acf1955f383b0ca7c265a6e88a2b803de737
  package   4eb0abcd7d722154df07cf92f03336242686f4edeea5abfd9b7f614ce6694343
```

Three consecutive full cascades leave **0 files modified**. A fresh clone verifies clean on all
three layers. The bundle was rebuilt, unpacked into a scratch directory, and verified there by the
three commands the Zenodo notes give a reviewer — all `INTACT`, with the 44 MB model artifact
present, so the archive verifies without a rebuild.

### Why this one is worth remembering

The previous record closed by saying the work is finished when the record says what happened. This
one adds a second edge: a check that cannot fail is not a check. The manuscript gate did fail, said
so in its own record, and the verifier built to police it read straight past the verdict and
reported clean. Everything downstream — CI, the bundle, the reviewer's three commands — inherited
that blind spot. It was found only because a single file kept dirtying the working tree.

### Follow-up: the bundle cited a commit that no longer existed

Immediately after the push above, the release bundle was checked and found to record
`git_commit: 9d6745c470942986668c1935b1886ecbf8cec78c+dirty`. Both halves were wrong for a
published artifact. It was cut from a tree with uncommitted changes, so the commit it names does not
describe the bytes in the archive; and those checkpoint commits were squashed before the push, so
the commit does not exist in the repository at all. `git cat-file -e` on it fails.

Nothing downstream caught this. The bundle checker compares the archive against its own manifest and
was perfectly happy — the archive was internally consistent, it simply carried false provenance.

This had to be fixed before an upload rather than after, because a Zenodo DOI is immutable: an
archive published with a dangling commit reference stays that way.

`make_release_bundle.py` now refuses to build from a dirty tree, naming both failure modes, with
`--allow-dirty` for throwaway builds. The correct release order is therefore: cascade, commit, then
build — the bundle is written to a gitignored directory, so building never dirties the tree it
records.

The digests in the section above were correct at commit `ff08966`. Adding this guard moved them, as
any change to a locked file does:

```text
  evidence  6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081
  claim     a5a92ad356a571fb30aa26254745bf70d445c9496a79b569388832aa924d5366
  package   28ecf262b97c6dd460011b6be4c5c43b30fe9109e114821ad7134d84eb088396
```

Separately, and left alone deliberately: seven evidence-locked JSONs from the frozen upstream stages
(23.2H, 24F, 25, 26) still carry a `runtime_minutes` or `runtime_seconds` field. Re-running those
stages would churn the tree. They are never re-run — the cascade touches only the three Gen-1 lock
stages — and rewriting a frozen artifact to strip a field is exactly what the no-rewrite rule
forbids. The count is stated here so the exposure is known rather than discovered later.

### The front page was stale, and nothing was checking it

A sweep of every tracked document for superseded digests found two that assert current state rather
than history:

`README.md` still displayed the release's **first** evidence, claim and package digests — stale
across several cascades. A reader following the front page would have been given three digests that
verify against nothing.

Worse, the same paragraph still carried the uncorrected sentence: *"1,401 barcoded clones were split
across six observed experimental conditions."* That is the misstatement about Schaff et al.'s
experiment that was corrected in the manuscript during the earlier verification pass — it attributes
our analysable subset to their experimental design. The correction reached `MANUSCRIPT.md` and never
reached the README. It now matches the manuscript: a barcoded population was split across the six
conditions, and of the clones that experiment recovered, 1,401 carry a pretreatment profile.

The README could not simply be added to the package digest, because it inlined the package digest —
a file covered by a digest cannot quote that digest. It now quotes the evidence and claim digests
and points at `results/manuscript/GEN1_PACKAGE_DIGEST.json` for the third, which is the same
no-cycle rule the two upstream records already follow. It is covered by the package digest, and a
new MS-C/MS-D check requires it to quote both current digests. Verified by planting a wrong digest:
the stage refuses, and restoring it returns the package digest to the identical value.

`GEN1_MANUSCRIPT_PACKAGE_V1.md` pinned the same two stale digests in its Entry block. Re-pinned,
with Amendment V1.1 recording what V1 said.

Both were found by sweeping for stale digests rather than by any check, which is why the check now
exists. The manuscript document had been protected this way since it was written; the front page
most people will actually read had not.

```text
  evidence  6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081
  claim     a5a92ad356a571fb30aa26254745bf70d445c9496a79b569388832aa924d5366
  package   db2c8aa7bd83dd770a6d6a60f57de03d3b651fbdcc942b6dd76be73ef8192e6a
```

### The same correction had not reached the abstract that gets submitted

Sweeping for the corrected sentence rather than for digests found a third copy. `SUBMISSION.md`
carries an "Abstract as submitted" block — the text that is actually pasted into a preprint server
or a journal's submission form — and it still read:

> in which 1,401 barcoded clones were split and exposed to six observed experimental conditions

That is the pre-correction wording. The manuscript had been fixed; its submission copy had not. Of
the three places this sentence lived, the correction had reached one.

The submitted abstract is now regenerated from the manuscript's own abstract, and a check requires
them to be identical after whitespace normalisation, so wrapping cannot hide a divergence. Proven by
restoring the old wording: the stage refuses, and putting it back returns the package digest to the
identical value.

This is the same failure as the README, one step further along: a correction applied to the primary
document and not to its copies. Both copies are now checked against the primary rather than
maintained by hand.

```text
  evidence  6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081
  claim     a5a92ad356a571fb30aa26254745bf70d445c9496a79b569388832aa924d5366
  package   9ff5419237d78c67b363f7a59f38db0b398e6c76b39b5ec049438492fdc7fac6
```
