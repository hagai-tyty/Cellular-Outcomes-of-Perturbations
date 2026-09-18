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

### The package document still described the old verification semantics

`REPRODUCIBILITY.md` told a reviewer that `EVIDENCE_INTACT` and `CLAIMS_INTACT` mean every hashed
file is byte-for-byte what was locked, and that *"anything else means something moved"*. After the
verifier fix that second half is wrong: `*_STAGE_REFUSED` means nothing moved and the stage itself
refused. A reviewer reading the old text would have gone looking for a changed file that does not
exist.

The document now states both conditions and distinguishes the two failure modes. This is the third
copy-of-a-fact problem in this pass — the fix changed behaviour, and the document describing that
behaviour was a separate copy that did not follow. The pattern is consistent enough to name: every
correction in this repository should end with a search for the other places the same fact is
written down.

```text
  package   ed4321ea1ca369e46c9421af8c44ff631b4c6925ed93975baaac6a09a3e782b5
```

### Acting on the pattern instead of the instance

Three separate copy-of-a-fact failures in one pass — the README's digests, the submitted abstract's
wording, the package document's description of verification — is a pattern, not three accidents. The
common shape: a fact is corrected in the primary document and the copies are maintained by hand.

A sweep of every headline number across all live documents found **no** further divergence, which is
luck rather than design: nothing was checking. The manuscript is held to a strong contract — it must
*state* each of the seventeen pinned numbers with the words around it — but that tracer cannot be
pointed at the README or the submission pack, which quote only a subset.

So they are held to the complementary contract: any number they quote must be a locked value, and
the `p` floor applies to them exactly as it does to the manuscript. DOIs are stripped before the
scan, because `10.1016/j.xgen.2026.101191` ends in a six-decimal fragment that is not a number this
project reports — a detail worth naming, since a checker that flagged it would have been switched
off within a day.

Both proven by planting a failure: one digit changed in the README's confidence interval refuses,
and `p = 0.000999` in the submission pack refuses. Restoring each returns the package digest to the
identical value.

The three secondary documents are now checked against the primary rather than maintained beside it:
the README against the current digests, the submitted abstract against the manuscript's abstract,
and any number either of them quotes against the locked artifacts.

```text
  package   4ae3404f987b3dee5a1b91dacfe78dc5e90a77172908d17954f0fcdb6b3eb362
```

One more instance of the same shape, found while checking the documents against each other: the
Zenodo notes gave a reviewer all three verify commands, but `REPRODUCIBILITY.md` — the document a
stranger is told to follow first — listed only two. Anyone following it never checked the package
layer at all. The third command is now documented there, with its digest read from
`GEN1_PACKAGE_DIGEST.json` rather than inlined, because this document is one of the files that
digest covers.

```text
  package   67db3051099a77b9ba822c35dde6bc55f31c4eb46ddfb463b58917c76e12db3d
```

### The environment lock could never have verified after a commit

Running the full suite left the tree dirty on exactly one file: `environment_lock.txt`. The diff was
a single line —

```text
  -e git+https://github.com/.../Cellular-Outcomes-of-Perturbations.git@<sha>#egg=cellfate_rx
```

`pip freeze` renders this project's own editable install that way, and it resolves the sha by asking
git for **HEAD at freeze time**. The file is one of the evidence-locked artifacts. So every commit
silently invalidated it: regenerate the export at any later commit and the evidence lock refuses on
a file whose dependency set never changed.

The exporter now normalises that line to a comment naming the package. The commit is already
recorded twice over — in git history and in the release bundle's `git_commit` — and what belongs in
an environment lock is the dependency set, which is stable. Verified: two consecutive exporter runs
now produce a byte-identical lock, and the normalised live `pip freeze` matches the normalised
committed body exactly, so the sha was the only volatile element. The 1,000 null draws, the 2,000
bootstrap replicates and the figure source data all regenerated byte-identically in the same run.

A note on how this was investigated, because the first attempt was wrong. A bisect script was run in
the background to find which test rewrote the file, while the exporter was being run in the
foreground to test the fix. It reported two unrelated test files as culprits — both artifacts of the
two processes writing the same file. Worse, the bisect's own `git checkout` reverted the fix
mid-work. It was stopped and the state re-established. Concurrent investigation of a file while
editing that file produces confident nonsense; the decisive check is simply whether a full suite
leaves the tree clean.

```text
  evidence  e9326ab66143337b0bc547f38bbf6599daa34f722bdcd2290280ebfd0078f491
  claim     7cd881265a94cc6d7847d2f5889e5dd9fa95e6a6b647c3b3ad739b5b9fa41a99
  package   9f6d9d009e5807652026aecfe4e2b8fb6488485b516fa621639f9248a0ca5876
```

The README needed re-pinning in this cascade and the new MS-C/MS-D check caught it, refusing the
manuscript stage until it was current. That is the check earning its place on the first cascade
after it was added.

### Correction: the environment-lock fix above was not actually in the tree

The section above states the exporter now normalises the editable-install line and reports the
resulting digests. The exporter change was real and was committed. **The regenerated
`environment_lock.txt` was not.** The bisect script described above holds a `git checkout` on that
file; it reverted the regenerated lock before the commit, and `git add -A` therefore staged nothing
for it. `git log -- environment_lock.txt` shows no commit of mine touching it at all.

The repository was left in a state where the exporter normalised and the committed lock did not —
internally consistent to every verifier, because the lock hashed the old file and the old file was
what was there, and broken the moment anyone ran the exporter.

It was caught by unpacking the published bundle and grepping its `environment_lock.txt` for the
normalised line, which was absent while the working tree appeared correct.

Two things this changes about the account above. The digests recorded there
(`e9326ab6…` / `7cd88126…` / `9f6d9d00…`) are the digests of a tree that did not contain the fix.
And the claim that a test re-runs the exporter is **not established**: with the normalising exporter
and the old committed lock, a full suite left the tree clean, which it could not have done if a test
had regenerated that file. What originally dirtied it during a suite run remains unexplained. The
honest position is that the volatile line is gone and the lock is now stable across runs, and that
the trigger was never identified.

Corrected digests, from a tree that does contain the fix:

```text
  evidence  de6429c8d7075249b03b508000ec14eafcfe048101f6923ce32541a9394e6c95
  claim     670f89a136c11cde54aeb62279ecf14e920fbcdfb7010986cae899b47a9c7aab
  package   5f5cf37d0b39aaa9006c5e08a2cf20db9e8c209b269b6de4f9aa7d825cf0ec89
```

The general lesson is the one the bisect already taught, in a second form: a background process that
writes to the repository does not merely produce misleading readings, it can silently undo work that
is then committed as done. Nothing that mutates the tree should run unattended alongside editing.

## The cascade became a script, and it found a stale pin on its first run

The release sequence was being carried out by hand: run the exporter, run three locks in a forced
order, and re-pin four documents that quote the digests those locks produce. Re-pinning by hand is
the same failure this pass spent its time removing — a fact maintained beside the primary rather
than checked against it — and it had already cost one refused cascade when the README was missed.

`experiments/cascade_gen1.py` now does it. The ordering is forced by what each layer hashes, and
the docstring says so rather than leaving it to be rediscovered: the evidence lock hashes its own
module and the exported source data; the claim lock hashes the plan that pins the evidence digest;
the manuscript hashes the documents that quote both. Out of order, a lock refuses against a digest
that has already moved.

Two things it deliberately does not do. It does not run the exporter — that is slow, it is the only
step that touches the derived source data, and regenerating it should be a decision rather than a
side effect. It does not build the bundle — the bundle records the commit it was cut from and
refuses a dirty tree, so the sequence is cascade, commit, then build.

**On its first `--check` it reported `PINS_STALE`.** `GEN1_MANUSCRIPT_PACKAGE_V1.md` pins both
digests in its Entry block, and had been quoting superseded values for several cascades. Every lock
verified clean throughout, because the package digest hashes that file — and hashing a file proves
only that it has not changed, never that what it says is true. The scratch helper used during this
pass did not know about the document at all; that is precisely why the helper had to stop being a
scratch helper.

Two contracts now hold it: every document the tool is responsible for must quote the live digests,
and any tracked non-record `.md` that quotes a lock digest must be registered with the tool. The
second one is the important half — it catches the next document that starts quoting a digest
without being added, which is how every stale pin here happened.

```text
  evidence  213593c3b71e7db7064a5bb704288d3d51bddf725d5002aeefaa5b936c0e53b8
  claim     a5b014a08c898d6f8aa2a627e2bb8f2f6c78efcc64c7d83b807d8a53b5fce2f9
  package   8f9394cfa719f25a7a274ca98ecf1bd490e8c215058e1e9608aef7778fe7f28d
```

---

## Phase 1 of the release: the manuscript in BMC Research-article form — 2026-09-12

Decision, 2026-09-12: release through Zenodo and bioRxiv, defer any journal, and hold the manuscript
to BMC Bioinformatics' Research-article structure regardless, so that it is journal-ready.

### What was done

- **Frozen plan.** Amendment V1.2 appended to `GEN1_MANUSCRIPT_PACKAGE_V1.md`. §2's V1 section list is
  left as written; the amendment records the new structure and the source it was read from.
- **Checks.** The required-section list follows V1.2. Five structure checks were added: section order,
  the eight Declarations headings, the abstract's three labels, at most 350 words, three to ten
  keywords. The references partition now accepts a top-level `## References`. Each new check has a
  negative control.
- **Manuscript.** Restructured by moving the V1 sections, not retyping them. Of 75 V1 paragraphs, 72
  are carried verbatim and 3 changed, all declared in advance: the pointer to the renamed
  Availability section, the verify-command block gaining the third command, and the "what is not"
  block gaining a line saying the archive includes the model file. 17 of 17 numbers still trace;
  the abstract is 243 words.
- **New prose**, each piece scanned clean by the claim lock's instrument before it was written: the
  author block, keywords, Use of AI assistance, the Discussion opening, Conclusions (assembled only
  from locked claims), the list of abbreviations, the Declarations, and figure legends taken from the
  figures' own on-image titles and panel labels.
- **Submission pack.** Declarations and keywords single-sourced in the manuscript. The submitted
  abstract regenerated and parity-checked. The AI note corrected: disclosure belongs in Methods under
  Springer Nature's policy, and bioRxiv's guidance asks for no separate disclosure. The fee recorded:
  £2,290 / $3,090 / €2,590 plus VAT, due only on acceptance, waiver requestable only at submission.
  Cover letter finished with the preprint and licence disclosures; pre-flight checklist rewritten.
- **Release.** README gains a DOI slot. Version 1.0.0 in `pyproject.toml` and
  `src/cellfate/__init__.py`, which had both said 0.1.0. The release bundle refuses any FILL marker in
  its six release documents; LATER markers are allowed and recorded in `BUNDLE_CONTENTS.json`. Six
  tests cover the gate.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     19 of 19 pass
  negative controls     10 of 10 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  digest pins           PINS_CURRENT
  lint, CI paths        clean
  bundle gate tests     6 passed
  full test suite       pytest's own exit code 0, no failures, 2176 tests collected;
                        the run left the working tree unchanged
```

The test count comes from `--collect-only`, not from pytest's summary line, and there is a reason.
`pyproject.toml` already sets `addopts = "-q"`; adding a second `-q` suppresses the summary line
entirely, which is why no full-suite log in this pass showed one. It also broke the first attempt to
count: at that verbosity `--collect-only` prints per-file counts rather than test IDs, so counting
test IDs gave 0 and the finish script refused before committing anything. Pass or fail was read
from pytest's own exit code throughout, and verbosity does not change that.

### A negative control that could not fire

The first cascade after the restructure was REFUSED at MS-F: "an unstructured abstract is caught"
stayed silent. The control built its broken copy by deleting `**Results.** ` with a trailing space.
The restructured manuscript puts each label on its own line, so the deletion matched nothing and an
unbroken copy was handed to the checker. The checker was right; the control never exercised it.

The step-2 smoke test had passed because its synthetic abstract put a space after the label — a test
document shaped differently from the real one. It is the lesson this project has recorded before: a
check has to read what the gate reads.

Fixed two ways. The label is removed whatever follows it, and every V1.2 control now requires its
copy to differ from the original, so a mutation that silently does not happen can never read as a
catch.

### Digests

```text
  evidence  0763f9229665c50c1cb74769e2271f0fc5664b4fb807679ee873ef87f984efa2
  claim     939a26854adf2942776bd96cbd38f1ea8973699d3769b8bc4d8b5fdc97197228
  package   3f0855ab579fed2096a4ebf3ddf76da498bb8916dbcb50d4a2a4c3f613ee92b3
```

### Deliberately not done yet

The release bundle now refuses to build: 25 FILL markers remain across `README.md`, `CITATION.cff`,
`.zenodo.json`, `MANUSCRIPT.md` and `SUBMISSION.md`. They are author details, declarations, the
AI-use description and the Zenodo DOI, all Phase 2 inputs. The restructured manuscript has not yet
been read by its author.

---

## Phase 2A of the release: corrections from an external review — 2026-09-13

An external review of the release workflow made five claims. Each was checked against the repository
and against primary sources before anything changed. All five held; two held more narrowly than
they were stated.

### What was checked, and what was found

- **Freeze wording.** The Stage-23 predictive analyses of WM989 are recorded on 2026-08-22. The
  Stage-23.5 plan that preregistered the ranking test was committed on 2026-08-26 and states that no
  Stage-25 ranking statistic existed. "Frozen before any result existed" was therefore true of the
  ranking statistics and false of the analysis as a whole. Counted with the gate's own pattern, the
  overstated form appeared 13 times: 6 in `MANUSCRIPT.md`, 4 in `SUBMISSION.md`, and one each in
  `README.md`, `.zenodo.json` and `CITATION.cff`. A first count from a line-based grep said 11 — a
  phrase broken across a line is invisible to grep — and it was corrected before anything was
  written.
- **Competing interests.** BMC counts financial interests that may be affected by publication. The
  author holds the copyright and offers commercial licences, so "none" must not be the default.
- **Licence wording.** PolyForm Noncommercial 1.0.0 permits personal use only "without any
  anticipated commercial application". `COMMERCIAL-LICENSING.md`, the README, the manuscript,
  `REPRODUCIBILITY.md` and `.zenodo.json` describe a commercial licence as needed only for
  revenue-generating use, and the notice says evaluation needs none. The licence governs, not the
  notice. **Not changed in this phase**: how to resolve it is the author's decision, in Phase 2B.
- **Figures.** The manuscript embedded none, and pandoc's default PDF route needs LaTeX, which is not
  installed.
- **Archive and checkout.** `REPRODUCIBILITY.md` said the model was not in "this package". The Zenodo
  archive contains it, and the pseudobulk cache as well.

### A verification that tested the wrong code

Run inside an unpacked archive, a plain `python -m cellfate.gen1_cli` imported
`D:/cellfate-rx/src` — the checkout's editable install — and exited 0. With `PYTHONPATH=src` it
imported the archive's own copy, and the documented example exited 0 with a score for each of the six
conditions. A verification of an archive that silently runs a different copy of the code is not a
verification. The flag is now required by the package document and by an MS-E check.

### What was done

- **Freeze wording** made precise in all 13 places: clone-held-out folds, features and exclusions
  fixed before any model was fitted; the ranking test frozen before any ranking statistic was
  computed, after earlier predictive analyses of the same data. `GEN1_CLAIMS.md` is a lock and keeps
  its wording; the manuscript says less than the lock, never more. MS-C check and MS-F control added.
- **Placeholder prompts** for competing interests and AI use sharpened. Nothing was filled on the
  author's behalf.
- **`REPRODUCIBILITY.md`** separates a checkout from the archive, and a new section verifies the
  unpacked archive with the three verifiers and the predictor under `PYTHONPATH=src`. MS-E check and
  MS-F control added. Recorded as Amendment V1.3.
- **`experiments/render_gen1_submission.py`** builds the BMC `.docx` with separate figure PDFs, and a
  bioRxiv `.docx` with the figures embedded above their legends, exported to PDF by Word. It refuses
  while any FILL marker remains unless run as a draft, and writes a manifest of output hashes. Seven
  tests, none of which needs pandoc. Not yet run: the tools are not installed.
- **`.gitattributes`** gains binary rules. `results/** text eol=lf` forced text handling on binary
  files. `example_clone_expression.npy` had been committed under it; its attribute moved from
  `text eol=lf` to `-text`, and git showed no change, so its committed bytes were unaffected. A
  rendered `.docx` or `.pdf` would have been corrupted on commit.
- **Draft renders** are gitignored, and the release bundle skips them: it walks the tree on disk, not
  git, so a draft would otherwise have been archived.
- **`SUBMISSION.md`** sections 5 and 7 now carry the render step, page inspection, the archive
  predictor, the licensing decision, the competing-interests disclosure, and the rule that the release
  tag never moves.

### Two script bugs the gates caught before anything was written

- The step-1 script re-wrapped any paragraph that grew long. `.zenodo.json` keeps its description on
  one line and `CITATION.cff` has no blank lines, so both would have collapsed into prose. The parse
  checks refused it; those two files are now edited without re-wrapping.
- A patch anchor containing a newline escape lost its backslash in the shell and matched nothing. The
  patch stopped before saving.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2189 tests collected
  FILL markers left     25, all Phase 2B inputs
```

### Digests

```text
  evidence  e467de64fd0cb0f9ca00d67db3c9ff99a2b04a3db4424a1207cf7eb53e01f71f
  claim     c6895e963baf43b024b5701eb004b723ff95ef0de422d0d426ebe32403112e9a
  package   b15ab5ddb010175b5bb99ff6a9f60e3935bfb76d0b3a8e319348b90c3dcbcff9
```

### Not done in this phase

The licensing wording, which waits on the author's choice; the FILL markers; the installs; rendering;
and the Zenodo DOI.

---

## Phase 3a of the release: the author's details, and the licence wording aligned — 2026-09-14

### What the author decided

```text
  name                  Hagai Aviv
  affiliation           Independent researcher
  correspondence        hagai.aviv.home@gmail.com
  ORCID                 0009-0004-4503-6629
  competing interests   holds the copyright and offers commercial licences; no income to date
  funding               none
  preprint licence      CC BY
  AI tools              Claude Opus 5, Claude Opus 4.8, Claude Haiku 4.5, Gemini 3.1 Pro,
                        GPT 5.6 sol, GPT 6 astra -- names confirmed by the author as written
  AI role               the Opus models wrote most of the code and documentation; the other
                        models checked that work; the logic and direction came mainly from the author
  licensing             "just do the same" -- read as option A: make the wording match the licence
```

The licensing line is an interpretation of the author's words, not a quotation, and is recorded as
one. Option A grants and removes nothing; it corrects text that described the licence inaccurately.

### Checked before it was written

- **ORCID.** The check digit is valid, and ORCID's public record for that iD gives the name
  "hagai aviv".
- **Email.** `LICENSE` and `COMMERCIAL-LICENSING.md` gave `hagay.aviv.home@gmail.com`, with a y. The
  author was shown the difference and chose `hagai.aviv.home@gmail.com`, with an i. Both files now
  carry that address; the previous one is recorded here as history.

### What changed

- **Author details** written into the manuscript's title page, `CITATION.cff` and `.zenodo.json`,
  and the cover letter's signature.
- **Declarations.** Competing interests disclose the copyright and the commercial-licence offer.
  Funding: none. The acknowledgements placeholder was removed, since the author named no one else.
- **Authors' contributions made consistent with the AI disclosure.** The single-author template said
  the author implemented the analysis and wrote the manuscript. The disclosure says the Opus models
  did most of the coding and documentation. The contribution statement now matches the disclosure.
- **Use of AI assistance** now names every tool and its role, as the author gave them.
- **The commercial boundary** now follows the licence in `COMMERCIAL-LICENSING.md`, the README,
  the manuscript, `REPRODUCIBILITY.md` and `.zenodo.json`. Commercial use needs a separate licence
  whether or not it has earned revenue yet, and every summary says the licence governs where the two
  differ. The notice no longer says that evaluating the software needs no licence. The body of
  `LICENSE` is unchanged; only its Required Notice line changed, for the email.
- **`CITATION.cff`**: an abstract line left running long by the Phase 2A edit is re-wrapped, with
  its content checked unchanged.

### Deliberately left open

Two confirmations the author can only give after reading the manuscript -- that he reviewed the work
and takes full responsibility for it, and that he approved the final manuscript -- remain FILL
markers, as do the Zenodo DOI and publication date. The release bundle refuses to build until all
of them are filled.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2189 tests collected
  FILL markers left     8
```

### Digests

```text
  evidence  ce77f340e6ae1a5ff4341952927fe49384f0a59254da09a13ad411e10ef2239b
  claim     cdd652f74900b3336c2bf05bfd80bf41e964eb7db30c0150e7c3a131d4126786
  package   67d604f26ac623793dcbdd010f37e4d9448f46a30ee2c647d414eaabf0296ecb
```

---

## Phase 3b of the release: an external manuscript review, checked before it was applied — 2026-09-14

### What was proposed

A review by another AI assistant, passed on by the author with the instruction to read it first and
change only what held up (paraphrased), proposed:

```text
  A    fill the Zenodo DOI slots, the release date and the author's two sign-offs
  B1   name the ranking test's evaluable population (892) in the abstract
  B2   say, after the strata, that condition is not among them and the interaction is not uniform
  B3   add to the abstract that pretreatment state "carries different information for different
       perturbations, rather than indexing a single resistance axis"
  C1   group the six conditions by mechanism in the Discussion: two MAPK-pathway inhibitors and two
       microenvironmental stresses carry the interaction; the two genotoxic chemotherapeutics do not
  D    an order of work, including: register C1 in the claim set first, because the claim lock
       would otherwise refuse it
```

### What was checked, and against what

- **The per-condition result (C1, B2).** Stage 23's treatment-level table (`stage_23_RECORD.md`),
  C1 log-loss gain of W5 over W4: Trametinib +0.03326, Dabrafenib +0.02592, CoCl2 +0.01894, Acid
  +0.01097, Cisplatin +0.00002, Doxorubicin -0.00332. The four that carry meaningful interaction are
  Acid, CoCl2, Dabrafenib and Trametinib, as the locked permitted phrasing in `GEN1_CLAIMS.md` has it.
- **The classes (C1).** Schaff et al. [1], full text through PubMed Central (PMC13261651), chose the
  six as three pairs: two targeted inhibitors (dabrafenib against V600E BRAF, trametinib against
  MEK), two agents simulating selective pressures (CoCl2 for hypoxia, acidic media for extracellular
  acidosis), and two chemotherapeutics (cisplatin, DNA cross-linking and breaks; doxorubicin,
  topoisomerase inhibition). ChEMBL agrees on each drug: B-raf V600E inhibitor, MEK1/2 inhibitor,
  DNA inhibitor, DNA topoisomerase II alpha inhibitor.
- **The schedule (C1).** The same paper gives dabrafenib, trametinib, CoCl2 and acidic media for four
  weeks, and cisplatin and doxorubicin for a treatment period followed by untreated recovery. The
  manuscript's own Methods says the same. **The class grouping and the schedule grouping are the same
  split.** The review did not mention the schedule.
- **The strata (B2).** Ship plan §8.9 lists the descriptive breakdowns as by outer fold and by
  pretreatment-depth bin, alongside a pairwise condition-ranking matrix; condition is not a stratum.
- **The single-axis clause (B3).** A single resistance axis whose effect differs by condition also
  changes the within-clone order, so the result does not show that the state is not a single axis.
  The manuscript already states that a clone-wide quantity contributes exactly zero to this metric by
  construction, so the metric cannot test that alternative. "Different perturbations" also reaches
  beyond the six observed conditions. The title and the first sentence of the abstract's Results
  already say the information is condition-specific.
- **The claim lock (D).** `run_gen1_claim_lock.py` does not read the manuscript. Its "FAILED STAGE"
  wording concerns a sentence in the adversarial corpus that the scanner misses, not new prose in the
  manuscript. And the lock cannot take a new claim: the ceiling may be lowered, never raised. New
  prose is held by the manuscript stage's scan, with the same instrument.

### What was decided

```text
  B1   APPLIED      one sentence added to the abstract's Background
  B2   APPLIED      reworded: "The breakdown is by fold and by depth only", pointing to Limitation 3
  C1   APPLIED      corrected: class and treatment schedule named together, as coinciding, and
                    written as a non-claim, after the Limitation 2 paragraph
  B3   NOT APPLIED
  D    NOT FOLLOWED on registering a claim; the cascade ran the locks in order, as after any edit
  A    NOT APPLIED  the author asked for the DOI slots to be left alone, and the two sign-offs wait
                    until the author has read the manuscript
```

`SUBMISSION.md`'s abstract as submitted was regenerated from the manuscript. No existing manuscript
line was altered: the edit script refused to write if one was, and it added 15 lines.

### The instrument refused the first wording

The first C1 draft called cisplatin and doxorubicin "chemotherapeutics". The scanner fired on it
(`3_clinical_recommendation`, cue word "therapeutic") and the edit script refused to write. The pair
is now described by mechanism, as "DNA-damaging agents". The refusal is recorded because it happened.

### Noticed, not changed

The source paper gives two different recovery lengths for the treat-then-recover arms of the barcoded
experiment. One passage has cisplatin for two weeks then two untreated, and doxorubicin for 2.5 then
1.5; another has both for two weeks then three untreated. The manuscript's "each arm spanning four
weeks in total" matches the first. Both passages agree on which arms were continuous, which is all C1
depends on.

### Tools installed, as the author approved

```text
  pandoc      3.11     winget, JohnMacFarlane.Pandoc
  svglib      2.2.0    pip, into the project environment
  reportlab   5.0.1    pip
  cffconvert  2.0.0    pip
```

### Draft render

The first run of the renderer, now that its tools are installed.

```text
  command      python experiments/render_gen1_submission.py --draft
  exit code    0
  tools        pandoc 3.11; the PDF exported by Word 16.0
  wrote        MANUSCRIPT_BMC.docx, MANUSCRIPT_bioRxiv.docx, MANUSCRIPT_bioRxiv.pdf (15 pages),
               COVER_LETTER.docx, figure_1_design.pdf, figure_2_primary.pdf, figure_3_robustness.pdf
  where        results/manuscript/submission/draft/, ignored by git; nothing from it entered the tree
```

Checked after rendering, not assumed from the exit code:

- Each of the three new passages appears exactly once in the text of both manuscript .docx files, and
  the manuscript's three FILL markers show in both, as a draft should.
- The bioRxiv .docx embeds its three figures. The BMC .docx embeds none, as intended; its figures are
  the separate PDFs.

**A defect the render exposed.** On the PDF's title page each digest line wraps: the last 11
characters of the 64-character digest fall onto a second line, so a code line holds 71 characters.
Code blocks are set in Consolas 11 pt on a Letter page with one-inch margins, from pandoc's default
styles; the renderer sets no smaller code style. 43 of the manuscript's 96 code-block lines are longer
than 71 characters, in 12 of its 16 blocks, among them both digest blocks, the limitations, the
tool's refusals and the references; at the observed capacity every one of them wraps. The draft
outputs are not committed. The renderer is corrected in the next step rather than here, so that this
entry records what the first render actually produced.

Visual inspection was partial: the in-app browser drew the top of page 1 and then timed out, with the
app window hidden. The count above does not depend on it.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2189 tests collected: 2188 passed, 1 skipped
  FILL markers left     8
```

### Digests

```text
  evidence  ce77f340e6ae1a5ff4341952927fe49384f0a59254da09a13ad411e10ef2239b   unchanged
  claim     cdd652f74900b3336c2bf05bfd80bf41e964eb7db30c0150e7c3a131d4126786   unchanged
  package   2c3ee19613bba0e4a1a1633f0f6144692b4226373952c7ebcc415c4c6eb237dc
```

---

## Phase 3c of the release: code blocks no longer wrap in the rendered documents — 2026-09-14

### What the first render showed

Phase 3b's draft render wrapped code blocks. Word, asked directly through COM about that render's
`MANUSCRIPT_bioRxiv.docx`, reported the page and the damage:

```text
  page                 US Letter, 612 pt wide, side margins 90 pt (1.25 in): text width 432 pt
  code style           Consolas 11 pt, from pandoc's default reference document
  line capacity        71 characters, which is where both digest lines broke
  code blocks wrapped  12 of 16, as text lines -> laid-out lines:
                         title-page digests  2 -> 4     tool refusals    4 -> 6
                         Data                2 -> 4     limitations     12 -> 21
                         Models              3 -> 4     never-claims     9 -> 10
                         the null            5 -> 6     Generation 2     8 -> 12
                         verify digests      2 -> 4     licensing        4 -> 5
                         what is not         7 -> 12    references      23 -> 36
  not wrapped          the ranking scores, the strata table, delta_TOP1, the verify commands
```

Word's twelve are exactly the twelve blocks that a character count had found longer than 71.

### What changed

- `render_gen1_submission.py` sets the code character style to 8 pt in every .docx it writes, by
  changing that one size in the document's `styles.xml`; every other part of the file is copied as it
  is. At 8 pt a line holds 98 characters, and the manuscript's longest code line is 93. 8.5 pt would
  hold 92, one short; 8 pt is the largest half-point size at which every line fits.
- Before writing anything, the renderer refuses if a code line in the manuscript is longer than a
  line holds.
- Before exporting the PDF, Word counts the code blocks it lays out on more lines than they have. The
  count goes into `RENDER_MANIFEST.json`, and the renderer exits 4 if the count is not zero, or if Word
  did not report one: an unmeasured layout is not a pass.
- Six tests: the capacity formula reproduces the 71 characters Word showed at 11 pt; every manuscript
  code line fits at 8 pt; the style patch changes only the code size, and refuses when the style is
  missing or has no single size; the renderer refuses an over-long code line before writing; and a
  Word step that reports nothing is not read as zero. None of them needs pandoc or Word.

The renderer is in the evidence inventory, so the locks were cascaded. The draft outputs remain
uncommitted.

### Draft render after the change

```text
  command      python experiments/render_gen1_submission.py --draft
  exit code    0
  manifest     code_point_size 8.0, code_blocks_wrapped_in_word 0, pdf_exported_by_word true
  wrote        the same seven files as Phase 3b's render, into the same ignored folder
```

Measured separately through Word afterwards, not taken from the renderer's own report:

```text
                             code blocks   wrapped   code size   pages
  MANUSCRIPT_BMC.docx             16           0        8 pt       13
  MANUSCRIPT_bioRxiv.docx         16           0        8 pt       14
  MANUSCRIPT_bioRxiv.pdf                                           14   (Phase 3b's render: 15)
```

The renderer's new check and the separate measurement agree. Before the cascade, the renderer's own
tests ran on their own: 13 passed, 7 existing and 6 new.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  lint, as CI runs it   ruff 0.16.0 check src/ tests/ scripts/ plan_tests/: all checks passed, exit 0
  FILL markers left     8
```

### Digests

```text
  evidence  861fb144badb1a886bb7a76645bcca7f09db6e00d694f88b4515464af1101b3c
       was  ce77f340e6ae1a5ff4341952927fe49384f0a59254da09a13ad411e10ef2239b
  claim     36014402e26bec4f6af8da61d8ccc9509d4370379b82a805e7b4c81e1cd1bc81
       was  cdd652f74900b3336c2bf05bfd80bf41e964eb7db30c0150e7c3a131d4126786
  package   28ddd3968510aeb162197d5c056119912bfca6a2085156a6fba96d16013ce6d2
       was  2c3ee19613bba0e4a1a1633f0f6144692b4226373952c7ebcc415c4c6eb237dc
```

---

## Phase 3d of the release: the reserved Zenodo DOI written in — 2026-09-15

### What the author gave

```text
  DOI   10.5281/zenodo.22769563
```

### Checked before it was written

- **Form.** Zenodo's prefix and a record number, `10.5281/zenodo.<n>`. `CITATION.cff` with the DOI
  written in validates against the CFF 1.2.0 schema (cffconvert 2.0.0), checked on a scratch copy with
  a stand-in release date, so that the still-open date could not hide another fault.
- **Whether it is published.** Zenodo's public records API and doi.org's handle service both answer
  404 for it. That is what a DOI reserved on a saved, unpublished draft returns: Zenodo registers the
  DOI only when the record is published. It also means the number could not be confirmed from outside
  the author's Zenodo account; it is written exactly as the author gave it, and resolves once the
  record is published.
- **Nothing else waited on the stub.** No test, plan test or script refers to it or to a count of
  placeholders.

### What changed

```text
  README.md       citation paragraph            the DOI, linked to doi.org
  CITATION.cff    doi, identifiers[doi].value   the DOI
  MANUSCRIPT.md   Availability of data          the DOI, linked to doi.org
  SUBMISSION.md   cover letter                  "archived at Zenodo under DOI ..."
```

The pre-flight checklist in `SUBMISSION.md` now ticks four items that were already done and recorded:
the author block and declarations, the licensing choice, and its application (Phase 3a), and the
installs (Phase 3b). The Zenodo item stays unticked: that the repository's GitHub integration is off
cannot be seen from here.

The locks were cascaded: `README.md` and `CITATION.cff` are in the archived repository, and the
manuscript changed.

### Deliberately left open

```text
  CITATION.cff date-released   the Zenodo publication date; written when the bundle is cut
  the two sign-offs            wait until the author has read the manuscript
```

The release bundle still refuses to build while these three remain.

### Checks

```text
  placeholders, four files    FILL 8 -> 3: the date and the two sign-offs; LATER 2 -> 2
  DOI occurrences             README 2 (text and link), CITATION.cff 2, manuscript 2 (text and
                              link), SUBMISSION.md 1; the stub is gone everywhere
  manuscript                  structure intact; no overstated freeze wording; the submitted
                              abstract still matches; no number drift in README or SUBMISSION.md
  lines replaced              README 1, CITATION.cff 2, MANUSCRIPT.md 1, SUBMISSION.md 5 (the DOI
                              line and four checklist ticks); nothing else moved
  cffconvert, as written      exit 1; the error names date-released, which is still FILL, and not doi
  cffconvert, stand-in date   exit 0: valid against the CFF 1.2.0 schema
  cascade                     exit 0, PINS_CURRENT
  draft render                exit 0; Word counted 0 wrapped code blocks; the PDF exported
  the rendered documents      both manuscript .docx files carry the DOI once, linked to
                              https://doi.org/10.5281/zenodo.22769563; the cover letter carries it
                              once; none carries the stub; the manuscripts show only the two
                              sign-off markers
```

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     3
```

### Digests

```text
  evidence  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
       was  861fb144badb1a886bb7a76645bcca7f09db6e00d694f88b4515464af1101b3c
  claim     712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
       was  36014402e26bec4f6af8da61d8ccc9509d4370379b82a805e7b4c81e1cd1bc81
  package   164c9c5b7947897df7c423ea1cc29b8982704b3fddb9ce7a41c028ab7452fb1a
       was  28ddd3968510aeb162197d5c056119912bfca6a2085156a6fba96d16013ce6d2
```

---

## Phase 3e of the release: rebuilding from the public data, shown and written down — 2026-09-15

### Why

Reading the manuscript, the author objected that others still could not replicate the work. Its
"What is not" section ended: "Naming a gap is not closing it. Both remain open." Checked: a GitHub
checkout can verify the locks and run the predictor, but nothing told a reader how to get from the
published data to the analysis inputs, and the Zenodo archive that carries the pseudobulk is not yet
public. The author approved closing the gap by showing the rebuild and writing it down, using the
data already on this machine and downloading nothing.

### What a rebuild needs, established from the code and the sources

- **Inputs.** Stage 22 reads, and records the size and SHA-256 of, 29 files for GSE279162 and 11 for
  GSE227151, plus the authors' code: five Schaff et al. files and two Rewind R1 scripts.
- **Both datasets, even for the primary result.** The Stage-22 builder reconstructs Role A before
  Role B and stops if Role A fails; the Stage-23A input audit checks both and writes
  `PROTOCOL_BLOCKED` if any file is missing or differs.
- **Sources, from listings and metadata only.** GEO's file lists: GSE279162_RAW.tar holds exactly the
  27 recorded sample files, at the recorded sizes; GSE227151_RAW.tar holds the six recorded GSM7092515
  and GSM7092516 files among its 39. The three Rewind barcode tables are not on GEO; Jain et al. 2024
  state that their processed data are deposited publicly and list the links in the paper's key
  resources table, which the text extraction did not carry. The Rewind code is Zenodo record 7707418
  (arjunrajlaboratory/iPSC_Rewind). Schaff et al. state that their code is on Zenodo; a search finds
  dylanschaff/Schaff_manuscript, 10.5281/zenodo.13935305. The SHA-256 checks decide whether any
  given copy is the recorded one.
- **GEO metadata files.** The series matrix and family XML of both datasets are required and hashed,
  and no stage reads their contents: searched across `experiments/`, `src/` and `scripts/`.

### The rebuild, attempted twice

**First attempt: refused.** A fresh clone of `8ccdbe1`, with `PYTHONPATH` on the clone's own `src`
and no analysis cache. Stage 22 stopped after 13 s: `GSE227151: BENCHMARK_BLOCKED_LINKAGE`. The cause
is not the code: the builder is byte-identical to the one that produced the committed outputs, whose
manifests record its SHA-256 as `785b9811cee22fe6...`. On 2026-08-24 the three Rewind barcode tables
were moved from `D:\GSE227151_Rewind\` into `r1\`, and the frozen loader reads them from the root.
`stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md` recorded that breakage then and left it unrepaired,
because the loader is frozen; it is not repaired here either.

**Before the second attempt, every input was hashed against the Stage-22 record.**

```text
  GSE227151   11 of 11 files and 2 of 2 author scripts byte-identical
  GSE279162   29 of 29 files and 5 of 5 author scripts byte-identical
```

**Second attempt: passed.** The Rewind inputs were copied, not moved, into the loader's layout at
`D:\cellfate-repro-data\GSE227151_Rewind\`; the originals were only read, and WM989's folder was used
where it is. A new fresh clone of `8ccdbe1`, again with `PYTHONPATH` on its own `src` and no cache:

```text
  step   exit   time    last line
  22     0      102 s   GSE227151 and GSE279162 BENCHMARK_READY_WITH_DECLARED_MISSINGNESS;
                        gates: all pass; OVERALL: STAGE_23_READY
  23a    0       74 s   OVERALL: PROTOCOL_FROZEN
  23b    0       57 s   OVERALL: ROLE_A_SIGNAL_PASS
  23c    0       80 s   OVERALL: ROLE_B_ADDITIVE_PASS
  23d    0       69 s   OVERALL: INTERACTION_PASS_MULTI_TREATMENT
  23f    0        3 s   next: STAGE 23R
  24b    0      141 s   C1_W5 BYTE_IDENTICAL, R1-R3 true
  24c    0       29 s
  25a    0       10 s
  total         591 s   including the copy and the clone
```

### What came back

Six tracked files differ from the commit. Every other committed output is byte-identical, including
all the Stage-22 tables and every Stage-23 manifest, result and out-of-fold file.

```text
  stage22_rewind_benchmark_manifest.json     local_source_path_used: the copy's folder
  stage22_prospective_benchmark_results.json that manifest's recorded bytes and SHA-256
  stage24/stage24b_reproduction.json         runtime_minutes only
  stage24/stage24c_serialization.json        runtime_minutes only
  stage25/stage25a_observed.json             runtime_minutes only
  stage24/stage24_w5_artifact.json           line endings only (CRLF); identical once normalised
```

```text
  WM989 clone pseudobulk   rebuilt c61b9fd1d429e071...  = the original
  model artifact           rebuilt 954cef7cff296d99...  = its locked SHA-256
  observed statistic       identical apart from runtime: delta_RANK 0.05160531888390629,
                           R(W1) 0.692654, R(W4) 0.692176, R(W5) 0.743781, 892 clones
```

### Not re-run

The 23E permutation nulls (about 5.6 h) and the 25b ranking null (10.7 h). Each is over the three-hour
limit this project sets for runs made on the author's behalf, and both refit on the inputs shown
identical above.

### What changed

- `REPRODUCIBILITY.md` gains §2.1, "Rebuilding from the public data": every input and where it comes
  from, the folder layout the loaders read, the commands, and what this rebuild showed. §2 and §5 now
  point to it; §5 keeps saying that the Zenodo archive includes the model and the pseudobulk.
- The manuscript's "What is not" section describes the raw-data row as recorded and rebuildable, and
  its closing line no longer says the gaps remain open, because they do not.

The edit script refused its first dry run: one new line of §2.1 was 101 characters, over the
100-character limit it enforces on new lines. The paragraph was rewrapped; nothing had been written.

**The full suite failed once, and the test was right.** After the first apply it ran 1 failed,
2,193 passed, 1 skipped, with pytest's own exit code 1:
`test_the_package_names_what_it_does_not_contain` requires `REPRODUCIBILITY.md` to keep saying
"Naming a gap is not closing it", and the §5 edit had rewritten the paragraph that opened with that
sentence and dropped it. The sentence is the package's commitment to naming its gaps, and it is still
true of the raw data, which stays outside every copy. It was restored in a form that says so ("Naming a gap
is not closing it, so each one has a way through."), the test was not changed, and the locks were
cascaded again before the suite was re-run.

The draft render after the change exited 0; Word counted 0 wrapped code blocks and exported the PDF.

The clone and the copied inputs on `D:\` were removed once these results were read.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     3
```

### Digests

```text
  evidence  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
       was  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
  claim     712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
       was  712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
  package   239b83389086ae4e717d38a9e55baea454c0fecd631297356f431cd7aafc95e8
       was  164c9c5b7947897df7c423ea1cc29b8982704b3fddb9ce7a41c028ab7452fb1a
```

---

## Phase 3f of the release: the Conclusions say what the result means — 2026-09-15

### Why

Reading the manuscript, the author found that the Conclusions dwelt on what the work cannot do, when
the Discussion already carries a whole section of limitations. The author asked for the Conclusions
to say what was found, what it means and what it could mean, and how others could take it further,
without referring to Generation 2. The replacement wording for the Conclusions section and for the
abstract's Conclusions was shown to the author before anything was written, and approved on
2026-09-15.

### The plan rule, and the amendment

Amendment V1.2 said the Conclusions "may be assembled only from the locked claims and their
qualifiers". Amendment V1.4, appended to `GEN1_MANUSCRIPT_PACKAGE_V1.md`, lets them state the locked
result with its qualifiers, then its significance and the directions it opens, written as implication
or possibility and never as a finding. The claim ceiling, the qualifiers, the Discussion's three
subsections and MS-A to MS-F are unchanged. V1.2 is not edited; V1.4 is added after it.

### Checked before it was written

```text
  claim scanner              clean on all five new paragraphs: the finding, what it means, the
                             directions, the release sentence, and the abstract's Conclusions
  overstated wording         none in the manuscript, the submission pack or the new amendment
  abstract                   284 -> 318 words, within the 350 limit; still carries WM989, the six
                             observed experimental conditions and "not death"
  qualifiers                 every qualifier marker is still present in the manuscript
  Generation 2               absent from both Conclusions; still in the Discussion alongside
                             "biological replication", so the replication check holds
  submitted abstract         regenerated, and equal to the manuscript's
  scope of the edit          16 manuscript lines replaced, all inside the two Conclusions blocks
  the release sentences      the script refuses unless REPRODUCIBILITY.md already has §2.1, since
                             two new sentences say the result can be reproduced from the public
                             data; Phase 3e had added it
  draft render               exit 0; Word counted 0 wrapped code blocks and exported the PDF
```

The edit script refused its first dry run on its own gate, not on the new text. It had run the
overstated-freeze check over the whole package plan, and the plan quotes "before any result existed"
on purpose: Amendment V1.3 records that wording as what it corrected. The manuscript stage holds only
release documents to that rule, so the gate was narrowed to the manuscript, the submission pack and the
new amendment text, and it now prints the plan's quotation as history. Nothing had been written.

### What changed

```text
  MANUSCRIPT.md                    the abstract's Conclusions; the Conclusions section, now four
                                   paragraphs: the finding, what it means, the directions it opens,
                                   and what is released
  SUBMISSION.md                    the abstract as submitted, regenerated from the manuscript
  GEN1_MANUSCRIPT_PACKAGE_V1.md    Amendment V1.4
```

Moved out of the abstract, not out of the manuscript: that the six conditions are the whole supported
vocabulary, that no claim is made about unseen conditions, other cell lines or patients, and that the
model emits no calibrated probability. All three stay in "What this does not show". The abstract keeps
the outcome qualifier, and both Conclusions say that the result has not yet been tested in an
independent system.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     3
```

### Digests

```text
  evidence  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
       was  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
  claim     712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
       was  712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
  package   bc6166a1242e4df99fae0796d73c4364001b6bf0139ce577402a4be9b036344d
       was  239b83389086ae4e717d38a9e55baea454c0fecd631297356f431cd7aafc95e8
```

---

## Phase 3g of the release: the README's counts, corrected — 2026-09-15

### What was stale

The author asked for the README's counts to be corrected. Three were checked against the repository
before anything was written; the fourth count in the README, 1,401 clones, is a traced headline number
and is correct.

```text
  README said                 written                    now
  2,259 tests                 1d2cb1d, 2026-08-28        2,195, as `pytest -q` collects them
  re-hash 77 files            396b11c, 2026-09-03        82
  the 54 locked artifacts     1d2cb1d, 2026-08-28        64
```

### How the new numbers were established

- **Files re-hashed by the three `--verify` commands.** The evidence lock's artifacts, the claim
  digest's `covers` and the package digest's `covers`. The three sets do not overlap, so the count of
  distinct files equals their sum: 64 + 5 + 13 = 82. At `396b11c` the same computation gives
  62 + 5 + 10 = 77, which is the number the README gave then; it was right when written and went stale
  as the locks grew.
- **Locked artifacts.** The evidence lock verifies 64 (`n_checked` in its `--verify` output). The
  README's 54 was already stale by `396b11c`, when the lock held 62.
- **Tests.** CI runs `pytest -q`, whose configured `testpaths` is `tests/`; it collects 2,195.
  `plan_tests/` holds scripts, not pytest tests, and collects none. Since `1d2cb1d` the number of
  test functions defined under `tests/` grew from 1,764 to 1,801 and no test file was deleted, so
  2,259 is not a count of this suite as CI runs it now; how it was taken is not recorded, and no
  explanation is invented here.

The edit script derives all three numbers when it runs, from the lock manifests and from pytest's
own collection, rather than taking them from this entry.

### What changed

Three numbers in `README.md`, one line each. Nothing else in the README moved.

### Registered results

```text
  README edit, dry run  every gate passed: numbers derived as 2,195 / 82 / 64; 3 lines changed;
                        no number drift; no overstated wording
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     3
```

### Digests

```text
  evidence  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
       was  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
  claim     712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
       was  712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
  package   e2d9d0dc3ce906bf49552753b5bc15720bf8aa4a4e0341ac695dfbb47e1d2eec
       was  bc6166a1242e4df99fae0796d73c4364001b6bf0139ce577402a4be9b036344d
```

---

## Phase 3h of the release: what the ordering metric is blind to, stated narrowly — 2026-09-15

### The correction, and why it is right

An external reviewer, passed on by the author before final approval, objected to the Results claim
that a general resistance-propensity axis "enters a model as an additive state term". That is too
categorical, and it conflates two different things: a clone-level *quantity*, and a contribution
*shared equally across conditions*. Only the second is what W4 carries and what a within-clone
ordering metric cannot see. A single biological programme whose influence varies by condition would
enter as an interaction, so this design cannot tell one programme from several.

The objection was checked and accepted. It is the same reading this project used on 2026-09-14 to
reject the opposite suggestion, that the result shows the state is "not a single resistance axis":
both directions claim more than a within-clone ordering test can support. The wording the reviewer
supplied is adopted:

```text
  A state contribution shared additively across all conditions cannot change their ordering within
  a clone. Allowing state effects to vary by condition improved that ordering in WM989.
```

### Where the claim was made, and what it says now

Eight passages made or implied the categorical version; all eight were changed in one pass, so the
document says one thing.

```text
  MANUSCRIPT.md  Background                 "any effect that acts on a clone as a whole"
                 Relation to prior work     "Any quantity acting on a clone as a whole"
                 Results                    "A general resistance-propensity axis ... enters a
                                            model as an additive state term"
                 abstract, Conclusions      "an analysis that averages state across conditions"
                 Conclusions, paragraph 2   "any analysis in which state acts on a clone as a whole"
                 Conclusions, paragraph 3   "programmes that separate one condition from another",
                                            which implied separate programmes
  SUBMISSION.md  section 3.2                "A clone-level propensity cannot, by construction"
                 cover letter               "a general resistance propensity cannot, by
                                            construction, satisfy"
```

The Results paragraph now also states what is not excluded: one shared programme whose influence
differs between conditions would appear as an interaction too, so the measurement locates
condition-dependent state effects and not the number of programmes behind them. The figure note
"R(W4) sits BELOW R(W1): the additive state term adds nothing to ordering" is unchanged: it was
already about the additive term, not about propensity.

Nothing here raises the claim ceiling. The locked ranking claim is an improvement over a
non-interactive additive model, which is exactly what the narrow statement says.

### Checked before it was written

```text
  claim scanner              clean on all nine new passages
  old categorical wording    none of the seven phrasings survives in either document
  the narrow statement       8 occurrences across the manuscript and the submission pack
  abstract                   318 -> 325 words, within the 350 limit; still carries WM989, the six
                             observed experimental conditions and "not death"
  qualifiers                 every qualifier marker still present in the manuscript
  Generation 2               still absent from both Conclusions
  overstated wording         none; no number drift in the submission pack
  submitted abstract         regenerated, and equal to the manuscript's
  scope of the edit          manuscript paragraphs changed = the six targeted; submission-pack
                             paragraphs changed = section 3.2, the cover letter and the regenerated
                             abstract; no other paragraph moved, word for word
  lines                      no new line over 100 characters, and none that would start markdown syntax
  draft render               exit 0; Word counted 0 wrapped code blocks and exported the PDF
```

The first dry run refused on two points, and nothing was written. The claim scanner fired on the
Background and Relation-to-prior-work sentences: "all six conditions" and "every condition" are its
cue words for forbidden claim 7, uniform benefit across all six conditions. The same words passed in
sentences carrying "cannot", which it reads as negated. Neither sentence claims uniform benefit, but
the wording was changed rather than the instrument: "across the six conditions" and "is the same for
each of them". The script's own line-length gate also fired, on two lines of the regenerated
submitted abstract; those lines are manuscript abstract lines behind a "> " prefix, which the
manuscript's own limit already holds, so the gate now exempts exact mirrors of manuscript lines.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     3
```

### Digests

```text
  evidence  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
       was  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
  claim     712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
       was  712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
  package   0b30b1444f37324d7cdc6b1ffae026aade9236f301bee13704cfea026307c91c
       was  e2d9d0dc3ce906bf49552753b5bc15720bf8aa4a4e0341ac695dfbb47e1d2eec
```

---

## Phase 3i of the release: the last placeholders filled, and the final render — 2026-09-15

### What the author decided

On 2026-09-15 the author approved the manuscript and set the release date, in these words: "date
tomorrow / manuscript approved fill every single placeholder left".

```text
  release date   2026-09-16, the Zenodo publication date
  manuscript     approved
```

### What was filled

```text
  CITATION.cff    date-released            "2026-09-16"
  MANUSCRIPT.md   Use of AI assistance     "The author reviewed the work and takes full
                                           responsibility for its content."
                  Authors' contributions   "H.A. reviewed the work and approved the final
                                           manuscript."
```

Each confirmation is the sentence its placeholder had quoted since Phase 3a; the edit script refused
to write any other. It also refused unless 2026-09-16 was the day after the day it ran, since the
author's word was "tomorrow".

The submission checklist ticks three items that are now true: the manuscript read and approved; every
FILL marker filled; `CITATION.cff` validated (cffconvert 2.0.0, exit 0, CFF schema 1.2.0).

### Not filled, and why

```text
  SUBMISSION.md   <<LATER: bioRxiv DOI>>             exists only once the preprint is posted
  SUBMISSION.md   <<LATER: request an APC waiver>>   belongs to a journal submission, which is
                                                     deferred
```

Filling either now would mean writing something untrue. Both sit in the cover letter, which is for a
later journal submission; the release bundle allows LATER markers and lists them in its contents.

### The numbers still reproduce, checked without rewriting the environment lock

The checklist's `python experiments/export_gen1_source_data.py` was not run as a script, because its
last export rewrites `environment_lock.txt` from the current interpreter. Since that lock was written,
20 packages were installed for rendering and validation on 2026-09-14 -- cffconvert, reportlab,
svglib and their dependencies -- and no existing package changed version. Rewriting the lock now
would enter those into the record of the environment that produced the results, and change a locked
artifact to do it.

Its three data exports were run instead, each of which refuses unless it reproduces the recorded
verdict, with their output sent to the gitignored `dist/export_check/`:

```text
  1,000 null draws             reproduce the verdict; identical to the committed file
  2,000 bootstrap replicates   reproduce CI95 [0.037197, 0.065571]; identical to the committed file
  figure source data           identical to the committed file
```

Compared after normalising line endings. The repository tree was unchanged by the check. That checklist
item is left unticked, because the script it names was not the thing that ran.

### A formatting slip from Phase 3h, corrected before the release commit

Reading the abstract back for the upload guide showed "**Conclusions.**" joined to the first line of its
text, while Background and Results keep their labels on their own lines. The cause was the Phase 3h
edit script: it flattened the paragraph before rewrapping it, so its check for a label on its own line
never saw the line break. The structure checks passed throughout, because they look for the labels in
order, not for their line breaks; and the rendered documents were unaffected, since a single line
break renders as a space. The label was put back on its own line, with the abstract's words proven
unchanged and the submitted abstract regenerated, and the cascade and the render below were run after
that correction.

### The final render

Rendered after the cascade, so the title page carries the final digests, with
`python experiments/render_gen1_submission.py`, not `--draft`. The files are committed in
`results/manuscript/submission/`.

```text
  exit code        0
  manifest         draft false; fill_markers_remaining 0; pdf_exported_by_word true;
                   code set at 8 pt; code_blocks_wrapped_in_word 0; pandoc 3.11
  outputs          MANUSCRIPT_bioRxiv.pdf (14 pages), MANUSCRIPT_bioRxiv.docx, MANUSCRIPT_BMC.docx,
                   COVER_LETTER.docx, figure_1_design.pdf, figure_2_primary.pdf,
                   figure_3_robustness.pdf, RENDER_MANIFEST.json
```

Checked in the text of each file, not assumed from the exit code:

```text
                            DOI   FILL   LATER   both confirmations   evidence digest   figures
  MANUSCRIPT_bioRxiv.docx   yes   0      0       yes                  yes               3 embedded
  MANUSCRIPT_BMC.docx       yes   0      0       yes                  yes               separate PDFs
  COVER_LETTER.docx         yes   0      2       --                   --                --
```

The two LATER markers in the cover letter are the bioRxiv DOI and the APC-waiver line, as intended.
The PDF's pages were not inspected visually here; the author read and approved the manuscript, and
the checks above are what was verified mechanically.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped (run after the label fix and the re-render)
  FILL markers left     0
```

### Digests

```text
  evidence  dd655cc7656b286f3a2d71ac3680264d31bfa2b4a1865b3864e300795fb824be
       was  7ed8c8b12cd287ab830e8f1ea1a6f821e33d11b8da7481279e27ae00385fe245
  claim     a748a55e540a4652e00af27ec93acc601ae0b88b630cee0bdc0f43ea953b25af
       was  712d30837e1f416fa7dc108d78ab56d47ba27f4fa87e9825cd463cc817f2d296
  package   00ce90658a1e7d725e1312c9ef8ddb59eceb6c73dc551224780b6ff34b2abff7
       was  0b30b1444f37324d7cdc6b1ffae026aade9236f301bee13704cfea026307c91c
```

---

## Phase 3j of the release: the release bundle, built and verified — 2026-09-15

### Built from the release commit

```text
  commit            e9f6184768d16ae2114eb9cf8befb9b0df1430ed   (Phase 3i, pushed; CI success)
  bundle            dist/cellfate-rx-gen1-bundle.zip
  files             401
  size              81.7 MB uncompressed, 74.5 MB compressed
  SHA-256           679a0db6c716e9f07f4a29e9b6e91d4f9e4a0e6ffd635f6ed47b1eb35e3d0b23
  lock digests      evidence dd655cc7656b286f3a2d71ac3680264d31bfa2b4a1865b3864e300795fb824be
                    claim    a748a55e540a4652e00af27ec93acc601ae0b88b630cee0bdc0f43ea953b25af
                    package  00ce90658a1e7d725e1312c9ef8ddb59eceb6c73dc551224780b6ff34b2abff7
  carries           results/stage24/stage24_w5_artifact.npz and
                    _cc_cache/stage23/GSE279162_pseudobulk.npz, which git does not
  pending LATER     the bioRxiv DOI and the APC-waiver line, both in the cover letter
```

`make_release_bundle.py` refused nothing: the tree was clean, no FILL marker remained, both gitignored
files were present, and every file a lock hashes was in the archive. `make_release_bundle.py --check`
then re-verified every member against `SHA256SUMS.txt`: BUNDLE_INTACT, 401 files checked, none changed, missing or mismatched, and the zip's hash equal to the recorded one. `BUNDLE_CONTENTS.json` records the
commit it was cut from and the zip's SHA-256; both were re-read against the commit and a fresh hash of
the zip, and matched.

### The unpacked archive verifies with its own code

As `REPRODUCIBILITY.md` §5.1 tells a downloader: unpacked into `D:\cellfate-bundle-verify\`, with no
git metadata inside, and run from its root with `PYTHONPATH=src`, so that the archive's own `cellfate`
is imported rather than the checkout's editable install.

```text
  git metadata inside the archive   absent
  cellfate imported from            the archive's own src/cellfate
  run_gen1_evidence_lock --verify   exit 0   EVIDENCE_INTACT
  run_gen1_claim_lock --verify      exit 0   CLAIMS_INTACT
  run_gen1_manuscript --verify      exit 0   PACKAGE_INTACT
  predictor, cellfate.gen1_cli      exit 0   a score for each of the six conditions, all
                                             SUPPORTED_KNOWN_CONDITION; ranking_status
                                             NOT_SUPPORTED because no verdict file was
                                             supplied, which is the tool working as specified
```

The unpacked copy was removed afterwards. `dist/` is gitignored, so the tree was unchanged.

### For the uploads

The GitHub release `gen1-v1.0.0` belongs on `e9f6184768d16ae2114eb9cf8befb9b0df1430ed`, the commit the bundle was cut from. This
entry is committed after it and is not inside the archive. The exact steps for Zenodo, the GitHub
release and bioRxiv were given to the author in `dist/UPLOAD_GUIDE.md`.

No test suite was run for this record-only commit, since nothing but this file changed. The three
`--verify` commands ran on it before the push, and CI ran after.

---

## Phase 3k of the release: Figure 2B's null labels no longer overlap — 2026-09-15

### What the author found

On page 14 of the rendered PDF, in Figure 2B, the labels "p95 0.0087" and "max of 1,000 0.0137" were
drawn over each other. The author asked for them to be staggered onto separate lines and the PDF
regenerated.

### Measured before anything was changed

The figure's own coordinates confirm it. Both labels sat on one baseline (y 410, 8.5 pt), centred at
x 250.2 and x 290.1, under markers only 0.005 apart on the axis. A deliberately narrow width estimate
(characters x size x 0.5, so that a reported overlap is a real one) puts the overlap at 21.7 px. The
same measurement over all three figures found no other overlapping pair.

### What changed

`experiments/make_gen1_figures.py`, Figure 2 only. The "max of 1,000" label moves one line down, 12 px,
still centred under its own marker; the two caption lines below it move down 12 px; the figure grows
from 476 to 488 px tall. The patched builder was run in memory before it was written:

```text
  figure_1_design.svg       byte-identical to the committed file
  figure_3_robustness.svg   byte-identical to the committed file
  figure_2_primary.svg      5 lines changed -- the height, twice; the lower label; the two caption
                            lines -- with every string and number unchanged
  overlapping label pairs   Figure 2: 1 -> 0; Figures 1 and 3: 0
```

The builder is in the evidence inventory and the figures are in the package lock, so the locks were
cascaded, and the final render was regenerated from the new figures.

After the builder was written, `python experiments/make_gen1_figures.py` regenerated the three files
on disk, and they were checked again there rather than trusted from the dry run:

```text
  builder diff              8 lines in, 5 out
  figure_1_design.svg       unchanged (canonical-LF SHA-256 equal before and after)
  figure_3_robustness.svg   unchanged
  figure_2_primary.svg      changed
  overlapping label pairs   0 in each of the three written figures
  Figure 2B labels          "p95 0.0087" at x 250.2, y 410; "max of 1,000 0.0137" at x 290.1, y 422
  cascade                   exit 0, PINS_CURRENT
```

### The Phase 3j bundle is superseded

The bundle recorded in Phase 3j, cut from `e9f6184` with SHA-256 `679a0db6c716e9f0...`, carries the
overlapping figure. It was never uploaded, and it is superseded: the bundle is rebuilt from this
phase's commit in the next entry, and the GitHub release belongs on that commit, not on `e9f6184`. The
Phase 3j entry stands as the record of what was built then.

### The final render

Regenerated after the cascade with `python experiments/render_gen1_submission.py`, not `--draft`, and
committed in `results/manuscript/submission/`.

```text
  exit code        0
  manifest         draft false; fill_markers_remaining 0; pdf_exported_by_word true;
                   code_blocks_wrapped_in_word 0
  PDF              14 pages
```

```text
                            DOI   FILL   LATER   both confirmations   evidence digest   figures
  MANUSCRIPT_bioRxiv.docx   yes   0      0       yes                  yes               3 embedded
  MANUSCRIPT_BMC.docx       yes   0      0       yes                  yes               separate PDFs
  COVER_LETTER.docx         yes   0      2       --                   --                --
```

The Figure 2 embedded in the bioRxiv document is the new one: its SVG carries the 488 px height. The
two labels now sit on baselines 12 px apart at 8.5 pt, so they cannot meet whatever text metrics Word
uses to draw them.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 of 20 pass
  negative controls     13 of 13 fire
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0, 2195 tests collected: 2194 passed, 1 skipped
  FILL markers left     0
```

### Digests

```text
  evidence  60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa
       was  dd655cc7656b286f3a2d71ac3680264d31bfa2b4a1865b3864e300795fb824be
  claim     fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61
       was  a748a55e540a4652e00af27ec93acc601ae0b88b630cee0bdc0f43ea953b25af
  package   fef252d4c2bbaf52dcb6f9755779a0f49af9b805f54dee282169c03eda343375
       was  00ce90658a1e7d725e1312c9ef8ddb59eceb6c73dc551224780b6ff34b2abff7
```

---

## Phase 3l of the release: the release bundle, rebuilt and verified — 2026-09-15

The Phase 3j bundle was superseded by the Figure 2B correction in Phase 3k before anything was
uploaded. This entry records its replacement; the Phase 3j entry stands as the record of what was built
then.

### Built from the Phase 3k commit

```text
  commit            65ac157785db4bed4547e996158eca79e7ee6be1   (Phase 3k, pushed; CI success)
  bundle            dist/cellfate-rx-gen1-bundle.zip
  files             401
  size              81.7 MB uncompressed, 74.5 MB compressed
  SHA-256           1b32277604ffefa47c582e6ad83f612c103418c8c054ea235cbf455103f3b0b6
  replaces          the Phase 3j bundle from e9f6184, SHA-256 679a0db6c716e9f0..., never uploaded
  lock digests      evidence 60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa
                    claim    fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61
                    package  fef252d4c2bbaf52dcb6f9755779a0f49af9b805f54dee282169c03eda343375
  carries           results/stage24/stage24_w5_artifact.npz and
                    _cc_cache/stage23/GSE279162_pseudobulk.npz, which git does not
  pending LATER     the bioRxiv DOI and the APC-waiver line, both in the cover letter
```

`make_release_bundle.py` refused nothing: the tree was clean, no FILL marker remained, both gitignored
files were present, and every file a lock hashes was in the archive. `make_release_bundle.py --check`
then re-verified every member against `SHA256SUMS.txt`: BUNDLE_INTACT, 401 files checked, none changed, missing or mismatched, and the zip's hash equal to the recorded one. `BUNDLE_CONTENTS.json` records the
commit it was cut from and the zip's SHA-256; both were re-read against the commit and a fresh hash of
the zip, and matched.

### The unpacked archive verifies with its own code

As `REPRODUCIBILITY.md` §5.1 tells a downloader: unpacked into `D:\cellfate-bundle-verify\`, with no
git metadata inside, and run from its root with `PYTHONPATH=src`, so that the archive's own `cellfate`
is imported rather than the checkout's editable install.

```text
  git metadata inside the archive   absent
  cellfate imported from            the archive's own src/cellfate
  run_gen1_evidence_lock --verify   exit 0   EVIDENCE_INTACT
  run_gen1_claim_lock --verify      exit 0   CLAIMS_INTACT
  run_gen1_manuscript --verify      exit 0   PACKAGE_INTACT
  predictor, cellfate.gen1_cli      exit 0   a score for each of the six conditions, all
                                             SUPPORTED_KNOWN_CONDITION, identical to the Phase 3j
                                             bundle's; ranking_status NOT_SUPPORTED because no
                                             verdict file was supplied, the tool as specified
```

The unpacked copy was removed afterwards. `dist/` is gitignored, so the tree was unchanged.

### For the uploads

The GitHub release `gen1-v1.0.0` belongs on `65ac157785db4bed4547e996158eca79e7ee6be1`, the commit this bundle was cut from, and
no longer on `e9f6184`. The upload guide in `dist/UPLOAD_GUIDE.md` was regenerated with this bundle's
commit, size and SHA-256, and sent to the author again.

No test suite was run for this record-only commit, since nothing but this file changed. The three
`--verify` commands ran on it before the push, and CI ran after.

---

## Phase 4a of the release: the Zenodo record and the GitHub release, published — 2026-09-16

Both uploads were made by the author, from `dist/` as built in Phase 3l. This entry records what was
checked around them, what went wrong once, and what is still open.

### Zenodo, before publishing

The files were compared, not assumed. Zenodo shows an MD5 for each file it received; each was
compared with the file on this machine:

```text
  file                          bytes      MD5 on Zenodo's draft               MD5 here
  cellfate-rx-gen1-bundle.zip   74546611   d6de0f183ca43a4d1d464297d26c40a7   equal
  SHA256SUMS.txt                   42689   f4186ac658dd3c49aef76bb5aebabb46   equal
  BUNDLE_CONTENTS.json             47960   8b2bfaec177942278be06285dfcf7cae   equal
```

So the draft held the Phase 3l bundle, with the Figure 2B correction, SHA-256
`1b32277604ffefa47c582e6ad83f612c103418c8c054ea235cbf455103f3b0b6`, and not the superseded Phase 3j
one.

The form was reviewed field by field from the author's screenshots. The draft was the one holding
the reserved DOI, `10.5281/zenodo.22769563`, and not a new upload. The licence was PolyForm
Noncommercial License 1.0.0, chosen from Zenodo's own list. One field was wrong: the related work
for GSE227151 was entered as "Is derived from" instead of "References". The author was asked to
change it before publishing. The published metadata has not been read back, so that change is not
confirmed here.

### Zenodo, after publishing

```text
  doi.org handle API     200, responseCode 1
  10.5281/zenodo.22769563  -> https://zenodo.org/doi/10.5281/zenodo.22769563
  Zenodo record API      one request: 504
```

DataCite registers a Zenodo DOI only when the record is published, so the resolving handle confirms
the publication. The record itself could not be read: Zenodo was timing out, and its site warned of
outages caused by automated traffic.

**A mistake, recorded as it happened.** After that 504, a background check was started that would
have retried Zenodo every 30 seconds for up to 10 minutes. The author stopped it, rightly: against a
service struggling with bots and scrapers, that is the traffic it is trying to block. It had sent
one request, also answered 504, before it was stopped. The rule since then: one request per check to
a service that is not ours, never a timed loop, and the author decides when to try again.

**Zenodo's GitHub integration.** With Zenodo down, this was checked from GitHub's side. The author
opened the repository's Settings, then Webhooks, and found no `zenodo.org` webhook. Without that
webhook, publishing a GitHub release cannot create a second Zenodo record with a different DOI. The
setting on Zenodo's account page itself was not seen.

### The GitHub release

Published by the author with the three files attached. Checked with a single request to GitHub's
public releases API and a single `git ls-remote`:

```text
  tag                   gen1-v1.0.0
  remote tag points to  65ac157785db4bed4547e996158eca79e7ee6be1   the bundle's commit
  target_commitish      65ac157785db4bed4547e996158eca79e7ee6be1
  title                 CellFate-Rx Generation 1 (gen1-v1.0.0)
  draft / prerelease    false / false
  published             2026-09-16T16:56:58Z
  notes                 carry the DOI 10.5281/zenodo.22769563

  asset                         bytes      state      SHA-256 reported by GitHub
  BUNDLE_CONTENTS.json             47960   uploaded   equal to the file in dist/
      bc0030b333f384401474608f63f1e607b76837a753e1ec467aaca08562a1ba1b
  cellfate-rx-gen1-bundle.zip   74546611   uploaded   equal to the file in dist/
      1b32277604ffefa47c582e6ad83f612c103418c8c054ea235cbf455103f3b0b6
  SHA256SUMS.txt                   42689   uploaded   equal to the file in dist/
      a1895c68aac54c02485594dc47cde18127e5351cddd070b0bcd2bde755debc41
```

GitHub computes those digests itself from what it received, so the release's copy of the bundle is
confirmed byte-identical to the verified one. The tag is on `65ac157`, not on the tip of `main`: the
one commit after it, `efcd791`, only adds the Phase 3l entry to this record.

### The pre-flight checklist is not ticked in this commit

`SUBMISSION.md` is in the manuscript package lock, so ticking a box changes the package digest and
means cascading the locks again. Several boxes describe work that is done and recorded but not
ticked: the Phase 3 steps from Phases 3i to 3l, with the note in Phase 3i that the source-data
exports were run in-process instead of through `export_gen1_source_data.py`; the Zenodo line in
Phase 2B; and the GitHub release above. They will be ticked in the one new commit that fills the
bioRxiv DOI, so that the locks are cascaded once, not twice.

### Still open

```text
  Zenodo      read the published record back: its metadata, including the GSE227151 relation, and
              its files' checksums
  Zenodo      download the ZIP from the published record; its SHA-256 must equal the one in
              BUNDLE_CONTENTS.json (checklist, Phase 4, first line)
  bioRxiv     the submission, by the author: MANUSCRIPT_bioRxiv.pdf, New Results, Bioinformatics,
              CC BY
  afterwards  the bioRxiv DOI into the Zenodo record's related works (author), and into the cover
              letter's LATER field in a new commit
  LATER       the APC-waiver line stays until a journal submission is actually made
```

No test suite was run for this record-only commit, since nothing but this file changed. The three
`--verify` commands ran on it before the push, and CI ran after.

---

## Phase 4b of the release: the bioRxiv submission — 2026-09-16

The author submitted the preprint and approved it. It is BIORXIV/2026/752106, version 1, and is
waiting for bioRxiv's screening, so it has no DOI yet. The form was reviewed from the author's
screenshots before approval. This entry records what was checked, what was corrected, and one
difference between the web abstract and the manuscript.

### The Phase 4a commit's CI

The Phase 4a commit, `1276bc35dcf3c63d35a790f468899085e8d01206`, was still running CI at the first
check. A second check, delayed by two and a half minutes, was stopped by the author; it printed
nothing. A single request afterwards found it `completed / success`.

### Zenodo: one open item closed

On the published Zenodo record, the author confirmed that the related work for GSE227151 reads
"References". The Phase 4a entry had left this unconfirmed. The published record's file checksums,
and a download of the ZIP checked against `BUNDLE_CONTENTS.json`, are still open.

### Corrected on the form before approval

```text
  author dialog         last name "aviv" -> "Aviv"; institution "none" -> "Independent researcher",
                        the manuscript's affiliation; "Mark as Corresponding Author" ticked
  title                 entered without "carries" and with "condition specific" unhyphenated;
                        replaced with the manuscript's exact title
  competing interests   "no competing interest" was selected; changed to "details given below",
                        with the manuscript's statement word for word
  external data         three links added: the Zenodo DOI, the GitHub repository, GSE279162
```

Every one of these was then read back on bioRxiv's own pages: the author list showed "Hagai Aviv
(Corresponding Author), Independent researcher" with the ORCID linked. The proofing page showed the
exact title, the three links, "Yes" with the competing-interests statement, funders "None", CC BY
4.0, and one file, `MANUSCRIPT_bioRxiv.pdf`, as the main manuscript file.

### The PDF bioRxiv built is ours, byte for byte

Before uploading, the page size was checked, because bioRxiv asks for US Letter: all 14 page objects
of `MANUSCRIPT_bioRxiv.pdf` have MediaBox 612 x 792 pt, which is 8.5 x 11 in. After uploading, the
author downloaded the PDF bioRxiv generated, `BIORXIV-2026-752106v1-Aviv.pdf`:

```text
  bioRxiv's PDF       386757 bytes  d215bd2bdce37e5a104dbd948887adf981109a93b9e87fa814c7d7a9baabc8f2
  ours, from 65ac157  386757 bytes  d215bd2bdce37e5a104dbd948887adf981109a93b9e87fa814c7d7a9baabc8f2
  cmp                 byte-identical
```

bioRxiv passed the file through unchanged. Every check recorded for that file in Phases 3k and 3l
therefore holds for the preprint: 14 pages, three embedded figures, Figure 2B with its labels
staggered, no placeholder left, and no code line wrapped.

### The web abstract uses brackets where the manuscript uses dashes

The abstract was pasted into the form as written in the upload guide. When the author viewed it,
bioRxiv had lost both em dashes in the Results paragraph: the first was shown as a line break in the
middle of a sentence, and the second as a semicolon. Neither was visible from the form itself.

A version using the HTML entity for the dash was prepared first, and discarded before it was given
to the author. The author pointed out that bioRxiv's own lists of accented and special characters do
not include the long dash, and asked why it was needed at all. It is not: there it only marks off an
aside, and brackets do the same. The sentence was rewritten:

```text
  before   ... preregistered in full — the metric, ... of the same data — a frozen ...
  after    ... preregistered in full (the metric, ... of the same data), a frozen ...
```

Proof, against the abstract in the upload guide, which is the manuscript's own:

```text
  ASCII only                       yes
  words, in the same order         324, identical
  every character-level change     "— " -> "(" and " —" -> "),", nothing else
```

The first run of that proof reported the words as different. The fault was in the check, not the
text: it compared against the guide's lines with their `>` quote markers still attached. With the
markers removed, it reported them identical, as shown. The author then viewed the abstract again on
bioRxiv: three paragraphs, the Results one unbroken, the brackets in place, and every number as in
the manuscript.

The PDF and `MANUSCRIPT.md` keep their dashes, which display correctly there. So the abstract on
bioRxiv's web page differs from the manuscript's, and from the "abstract as submitted" in
`SUBMISSION.md`, by exactly those two punctuation marks. `SUBMISSION.md` is not changed here.

### Still open

```text
  bioRxiv     screening; the preprint's DOI once it is posted
  afterwards  the bioRxiv DOI into the Zenodo record's related works (author), and into the cover
              letter's LATER field in a new commit, with the pre-flight checklist ticked
  Zenodo      the published record's file checksums, and the downloaded ZIP's SHA-256
  LATER       the APC-waiver line stays until a journal submission is actually made
```

No test suite was run for this record-only commit, since nothing but this file changed. The three
`--verify` commands ran on it before the push, and CI ran after.

---

## Phase 4c of the release: bioRxiv declined, and the preprint moves to Zenodo — 2026-09-18

### What bioRxiv decided

bioRxiv declined BIORXIV/2026/752106 during screening. The reason given is that it requires authors to
have an organisational affiliation: an organisation that provides oversight of research activities and
can adjudicate ethical issues or disputes. Nothing in the manuscript was questioned, and no scientific
objection was raised. The submission had been accepted, converted and acknowledged two days earlier,
as Phase 4b records; the decline is about affiliation only, and an independent researcher cannot meet
it without joining an organisation.

### What was checked before choosing a new route

```text
  BMC preprint policy      a posted preprint is not prior publication, is not counted against the
                           advance a study provides, and its DOI is disclosed at submission
  BMC fees                 no submission fee; the APC falls due only after editorial acceptance
  BMC waivers              discretionary waivers must be requested at submission; requests during
                           review or after acceptance cannot be considered. The submission-form tick
                           is not the request: a separate form follows within 14 days, and the
                           decision usually comes within two working days
  arXiv q-bio              no affiliation requirement, but a first submission needs an endorsement
                           from an established arXiv author in that area
  affiliation route        the Ronin Institute takes applications from independent scholars, about
                           $100 a year; IGDORE has suspended applications since 2024
```

Checked on 2026-09-18 against the servers' own support pages. Springer Nature's main pages now redirect
automated requests to a login, so the APC amount was not re-read today; it stands as recorded on
2026-09-12, £2,290 / $3,090 / €2,590 plus tax.

The author chose: a Zenodo preprint record now, then BMC Bioinformatics with the waiver requested at
submission. The preprint comes first because review takes months, a preprint costs nothing, and BMC
does not hold it against the submission.

### What changed in the repository

`SUBMISSION.md`, eight edits, and Amendment V1.5 appended to the package plan:

```text
  header block        PREPRINT SERVER is Zenodo, with the reason bioRxiv declined; RELEASE is the
                      new order, dated, and says what the 2026-09-12 plan had been
  form-entry notes    the bioRxiv funding note becomes BMC's; the bioRxiv AI-policy sentence goes
  render section      why the single-file PDF exists, and why it keeps the name MANUSCRIPT_bioRxiv.pdf
  cover letter        the preprint line now names Zenodo, with <<LATER: Zenodo preprint DOI>>
  checklist, 2B       the Zenodo line ticked: no webhook, DOI reserved on the kept draft
  checklist, 3        all eight steps ticked, with the note that the source-data exports were run
                      in-process so the environment lock was not rewritten
  checklist, 4        Zenodo published and the GitHub release ticked; the download check left open;
                      bioRxiv replaced by the Zenodo preprint, the preprint DOI, and the BMC line
```

The manuscript, the figures and every number are untouched. The published copies are untouched: the
PDF keeps its file name, so the Zenodo archive, the GitHub release and this checkout still name the
same file.

### The edit script's gates, including one that was wrong

The edit ran as a dry run first. Its first line-length gate flagged 21 lines, of which only one was
ours: it measured every line in both files rather than the lines the edit adds, and both files already
contained longer lines from before. The gate was corrected to measure only added lines, and the one
real offender, a new line of 102 characters, was rewrapped. A second gate demanded that the number of
changed hunks equal the number of edits; that was wrong too, since one edit can leave unchanged lines
inside it. It was replaced by a stronger proof: every line that leaves must belong to some edit's old
text, and every line that arrives to its new text.

```text
  lines out / in           24 / 35, all accounted for by the eight edits
  claim scanner            no hit on any added line
  line length, added only  none over 100
  the plan                 appended to, not edited: the new text starts after the old file, whole
```

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 / 20
  negative controls     13 / 13
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0; 2,195 results: 2,194 passed, 1 skipped
  FILL markers left     0
  LATER markers left    2: the Zenodo preprint DOI, and the APC-waiver line
```

### Digests

```text
  evidence  60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa   unchanged
  claim     fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61   unchanged
  package   c4e62ecd27a860d5a7bca2c5b34890036f2cd8907b180978f60aaeec7cc2046b
       was  fef252d4c2bbaf52dcb6f9755779a0f49af9b805f54dee282169c03eda343375
```

### The committed render is one file behind, deliberately

`COVER_LETTER.docx` in `results/manuscript/submission/` still carries the bioRxiv sentence, because the
cover letter will change again when the preprint DOI arrives. It is regenerated then, once, rather than
twice. The manuscript outputs are unaffected: `MANUSCRIPT.md` did not change, and `RENDER_MANIFEST.json`
records its digest `c5bba1ce...`, still current, so `MANUSCRIPT_bioRxiv.pdf` stays byte-identical to the
copy on Zenodo and on the GitHub release.

### Still open

```text
  Zenodo      the preprint record: MANUSCRIPT_bioRxiv.pdf, type Preprint, CC BY, linked to the
              archive record, by the author
  then        the preprint DOI into the cover letter's LATER field, the cover letter re-rendered,
              and the DOI added to the archive record's related works
  BMC         the submission, with the APC waiver requested at submission
  Zenodo      the published archive's ZIP downloaded back and its SHA-256 compared
```

---

## Phase 4d of the release: the AI disclosure is renamed, and stays in the Methods — 2026-09-18

### What the author asked, and what was checked

The author asked why the AI section sits so early in the PDF, whether it should move much further
down, and whether it could be renamed "AI use disclosure".

It was already the last subsection of the Methods, on page 6 of 14; it reads as early only because the
Methods themselves are early. Springer Nature's policy, which BMC follows, says that use of a large
language model should be documented in the Methods section, and in a suitable alternative part only
where there is no Methods section. Moving the disclosure to the end of the paper would therefore be the
first thing a screener could query. That was put to the author with the alternative — a pointer line in
the Methods and the full text in the Declarations — and the author chose to rename only and keep it in
place.

### The edit

```text
  MANUSCRIPT.md    "### Use of AI assistance" -> "### AI use disclosure"   one line
  SUBMISSION.md    the form-entry note now names the new heading, and says where the policy asks for
                   the disclosure: in the Methods, another part only if there is no Methods section
```

Nothing in the locks, the checks or the tests names that heading: it is not in `REQUIRED_SECTIONS`,
and the only references to it anywhere were the two prose lines in `SUBMISSION.md`. The dry run proved
the rest:

```text
  matches per edit             exactly one
  the disclosure's own text    unchanged, byte for byte
  its neighbours               still "### The preregistered ranking test" before, "## Results" after
  lines out / in               1 / 1 in MANUSCRIPT.md, 2 / 3 in SUBMISSION.md, all from the edits
  claim scanner                no hit on any added line
  old heading name left        0 mentions anywhere
```

### Re-rendered, and what that costs

The manuscript changed, so the submission files were rebuilt with
`python experiments/render_gen1_submission.py`, not `--draft`.

```text
  exit code                     0
  draft                         false
  fill_markers_remaining        0
  pdf_exported_by_word          true
  code_blocks_wrapped_in_word   0
  PDF                           14 pages, every page US Letter, 612 x 792 pt
  the heading in the PDF        "AI use disclosure", page 6; the old name appears nowhere
  MANUSCRIPT_bioRxiv.pdf        386758 bytes
                                009ea3c4aa6cb90cfd577007ec1e50c3268d13077a8caef0f139d060fc6f00be
                          was   386757 bytes, d215bd2b...c8f2, the copy on Zenodo and the release
```

So the repository's PDF is no longer the one archived at Zenodo and attached to the GitHub release.
That is expected for a correction made after a release: the archive is the frozen v1.0.0 snapshot and
does not change, while the checkout moves on. The preprint, when it goes up, carries this newer PDF.
The three figure PDFs were rewritten too, but only in their timestamps: 74 bytes differ, all inside
`/CreationDate`, `/ModDate` and the `/ID` reportlab derives from them.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 / 20
  negative controls     13 / 13
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0; 2,195 results: 2,194 passed, 1 skipped
  FILL markers left     0
  LATER markers left    2: the Zenodo preprint DOI, and the APC-waiver line
```

### Digests

```text
  evidence  60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa   unchanged
  claim     fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61   unchanged
  package   5641c0f6d136c01e1c0f03abeb59cd818c5de881702c826068ef81850fc213d6
       was  c4e62ecd27a860d5a7bca2c5b34890036f2cd8907b180978f60aaeec7cc2046b
  manuscript, canonical LF   51eba6f5ffaeabaf518b7da493599afeebcd51269fe3e5fe1b521649a6d0797a
                       was   c5bba1ce73d4671b6dc43da1067d16cc42b023402ea9aec6a61573d92bbcc972
```

The claim digest is unchanged, as it must be: a heading's name is not a claim, and no sentence of the
disclosure moved.

---

## Phase 4e of the release: the preprint is published, and the abstract states the measured half — 2026-09-18

### The Zenodo preprint record

The author published it. Read back with one request to Zenodo's public API:

```text
  version DOI     10.5281/zenodo.22829748          published, state done
  concept DOI     10.5281/zenodo.22829747          resolves to the newest version
  type / version  Preprint / v1        language eng        CC BY 4.0        open
  creator         Aviv, Hagai, Independent researcher, ORCID 0009-0004-4503-6629
  keywords        9, the same as the archive record
  related works   isSupplementedBy 10.5281/zenodo.22769563 (software)
                  isSupplementedBy the GitHub repository (software)
                  isDerivedFrom 10.1016/j.xgen.2026.101191 (journal article)
                  references GSE279162, references GSE227151 (datasets)
  file            CellFate-Rx-Gen1-preprint.pdf, 386758 bytes
                  Zenodo's MD5 fd8fbb27f7d73a5d882566188ca48aef equals this machine's
```

Before publishing, the draft had the file under its repository name, `MANUSCRIPT_bioRxiv.pdf`, on a
record that has nothing to do with bioRxiv. The content was right -- the MD5 matched -- but the name
would have been public and permanent, so it was replaced with the identical bytes under the name
above.

The same request closed an older item. The **archive** record, 10.5281/zenodo.22769563, was read back
in full: its metadata is as intended, and Zenodo's own MD5 for each of the three archived files equals
this machine's. The checklist line asking for a re-download of the ZIP is replaced by that check,
which tests the stored bytes rather than the transfer.

### What the author's reviewer proposed, and what was verified

Two additions to the abstract were proposed, plus a wording fix. Each was checked against the locked
results before any edit.

```text
  W4 below W1      TRUE.  R(W1) 0.692654, R(W4) 0.692176, R(W5) 0.743781. The Results body already
                   says the additive term contributes nothing and the whole gain is the interaction
  delta_TOP1       TRUE.  +0.115471, CI95 [+0.082960, +0.145740]; 82.8% of evaluable clones against
                   71.3%. But it was preregistered as a directional-consistency check that could
                   withhold support and never grant it
  abstract limit   324 words against the checked limit of 350: 26 words of headroom
```

The W4 sentence was added, in 21 words. The abstract already carried the argument's theoretical half
-- an additively shared state contribution cannot change within-clone ordering -- and now carries the
measured half beside it. The top-choice number was not added: honestly qualified it costs about 24
more words, which would have meant cutting qualifiers to fit, and putting a check that cannot grant
support in the abstract would misstate its status. It keeps its own Results subsection and Figure 3C.

### The edit

```text
  MANUSCRIPT.md    one sentence at the end of the abstract's Results paragraph
  SUBMISSION.md    the mirrored "abstract as submitted", kept word for word identical
  SUBMISSION.md    Zenodo is a repository, not a preprint server: the header label, and the sentence
                   that explained the single-file PDF by bioRxiv's conversion engine
  SUBMISSION.md    the cover letter's <<LATER: Zenodo preprint DOI>> filled with the concept DOI,
                   10.5281/zenodo.22829747, which stays right across versions, naming the submitted
                   version 10.5281/zenodo.22829748 beside it
  SUBMISSION.md    checklist: the preprint published, the DOI written, the archive read back
```

The new sentence:

> An additive expression term did not improve ordering over condition identity alone (R(W4) 0.692176,
> R(W1) 0.692654): the gain is the interaction.

### A gate that was wrong again

The dry run refused, reporting that the mirrored abstract no longer matched the manuscript's. The text
was fine; the check was not. Its extractor for the mirrored copy stopped at the next `---`, but the
mirrored block ends with its own quoted `> ---`, so it read on into the KEYWORDS line below and
compared 335 words against 324. The extractor was corrected to stop at the quoted rule, and then
validated on the **unedited** files first, where it must report the two copies identical: it did, 324
words each. Only then was the edit re-run. This is the second gate in two phases to be wrong in the
same way -- measuring more than the thing under test -- and the lesson is the same: prove a new check
against the unchanged file before trusting its refusal.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 / 20
  negative controls     13 / 13
  numbers traced        17
  abstract              345 words, limit 350; the two copies identical word for word
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0; 2,195 results: 2,194 passed, 1 skipped
  render                exit 0; draft false; 0 FILL; Word export true; 0 code blocks wrapped
  PDF                   14 pages; the new sentence present; 386917 bytes
                        43a99da0aa0e8f1cfe49f02774116307c138ed6299a08f8ba18b3ef4bf4dd076
  LATER markers left    1: the APC-waiver line
```

### Digests

```text
  evidence  60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa   unchanged
  claim     fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61   unchanged
  package   c48419e9eecf0d37398a39546f0710e290f43632571811c8a5518ea12f81cc1d
       was  5641c0f6d136c01e1c0f03abeb59cd818c5de881702c826068ef81850fc213d6
```

The claim digest does not move: both numbers were already reported in the Results, traced to locked
artifacts, and the sentence adds no claim above the lock.

### Still open

```text
  Zenodo      the published v1 carries the abstract without this sentence. The record can take a new
              version, which keeps v1 and gives v2 its own DOI, with the concept DOI resolving to the
              newest. The author decides whether to publish v2 with the PDF rendered here
  Zenodo      on the archive record, add the related work: Is supplement to, the preprint DOI
  BMC         the submission, with the APC waiver requested at submission and the preprint DOI
              disclosed
```

---

## Phase 4f of the release: audited against BMC's own rules, and the data sources named — 2026-09-19

Before preparing a journal submission, the manuscript was checked against BMC Bioinformatics' current
instructions, read from the journal's own pages rather than from memory. Seven things would have been
flagged. All seven are fixed here. The author also raised the one that matters most.

### What the journal requires, and where the manuscript stood

```text
  abstract        350 words max, no citations, Background/Results/Conclusions   PASSED (345, none)
  keywords        three to ten                                                  PASSED (9)
  sections        Background, Methods, Results, Discussion, Conclusions,
                  List of abbreviations, Declarations, References              PASSED
  declarations    all eight subheadings, "Not applicable" where it does not
                  apply                                                        PASSED
  figures         separate files, an accepted format, numbered in order of
                  mention; titles <= 15 words, legends <= 300 words            PASSED (PDF; 6-9
                                                                               and 22-40 words)
  LLM use         documented in the Methods                                     PASSED
  manuscript file editable, not a PDF                                           PASSED (.docx)
  abstract        minimise abbreviations                                        FAILED -- R(W4)
  references      Vancouver, numbered in order of first citation                FAILED -- both
  URLs            every web link gets a reference number and an access date     FAILED
  tables          numbered Table 1..n and cited in sequence                     FAILED -- none
  competing       use the author's initials                                     FAILED
  contributions   initials, and "read and approved the final manuscript"        FAILED
  title page      institutional address                                         no city or country
  data            where the data actually came from                             INCOMPLETE
```

### The data statement was wrong by omission, and the author caught it

The availability statement said the analysed data are in GEO. Not all of them are. Checked against the
frozen Stage-22 manifests and the GEO family XML itself:

```text
  GSE279162, Role B    all 29 analysed files are from GEO
  Schaff et al. code   five preprocessing files, from Zenodo 10.5281/zenodo.13935305, not GEO
  GSE227151, Role A    the two GSM sample sets and the series files are from GEO
  three tables         filtered10XCells.txt, stepThreeStarcodeShavedReads_BC_10X.txt and
                       stepThreeStarcodeShavedReads_BC_gDNA.txt are NOT in the GEO deposit: the
                       series' only supplementary file is GSE227151_RAW.tar, and neither name
                       appears anywhere in GSE227151_family.xml. They come from the authors' shared
                       data package, linked from the key resources table of Jain et al., taken on
                       2026-08-21
  Rewind code          two R1 scripts, from Zenodo 10.5281/zenodo.7707418, not GEO
```

A second shared package, taken on 2026-08-24, supplied the R2 and R3 folders for the re-examination
recorded in `stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md`. No file from it appears in either Stage-22
manifest, so nothing in the manuscript rests on it. Both packages are now named, with their URLs and
dates, in `REPRODUCIBILITY.md` §2.1; the manuscript's availability statement names the first, because
that is the data the findings rest on, and says plainly that a shared folder is not a persistent
identifier: what protects a reader is the recorded size and SHA-256 of every input file, which the
pipeline refuses to proceed without.

### What changed

```text
  title page      "Independent researcher, Ma'ale Adumim, Israel"
  abstract        the model labels R(W4) and R(W1) are gone: "The additive model did not itself
                  improve that ordering over condition identity alone (0.692176 against 0.692654)"
  Tables 1-5      model specifications; ranking scores; the null; the strata; the top-choice
                  diagnostic -- each captioned and cited in the text
  competing       "HA holds the copyright ... HA has received no income from it to date"
  contributions   "HA conceived ... The author read and approved the final manuscript."
  availability    rewritten: every source, in GEO or not, each with a reference
  references      6 entries -> 13, in Vancouver style, renumbered in order of first citation
```

The old numbering was not in citation order either: reference [2] was first cited in the Methods,
after [6] in the Background. Renumbering was done by replacing each citation with its key, ordering
the keys by first appearance and numbering them from there, with a gate that every key's citation
count survives the edit and that the numbers rise with position.

### The reference list stays fenced in the source, and leaves as a list

`_partition_references` needs the fenced block to exempt a cited paper's title from the claim scan --
"cancer drug resistance" in a title is not a claim this manuscript makes -- and its loophole control
depends on that shape. Rewriting the block as a Markdown list in the source would have disabled that
exemption. So the source keeps the fence and the renderer converts the block to a numbered list on
the way out, which is what a journal sees:

```text
  MANUSCRIPT_BMC.docx   5 real Word tables; no reference sits in a monospace paragraph
  the transform         every entry reappears whole, its words in order
```

Two tests were added for it. The first one failed on the real manuscript and was right to: the
transform left the closing fence behind. That is a bug the dry run could not have caught, because the
dry run never ran the transform. Two further iterations were the tests' own fault, not the code's --
a normalisation that stripped the brackets from "[software]" as if they were citation brackets, and a
comparison that did not collapse double spaces. Both were corrected before the code was trusted.

### One checker pattern was loosened, deliberately

The null table puts " | " between a label and its value, which the old pattern's `\s+` cannot cross:

```text
  before   ("null max", r"largest of [\d,]+ draws\s+@@")
  after    ("null max", r"largest of [\d,]+ draws\D{0,12}@@")
```

`\D{0,12}` is the same form the R(W1), R(W4), R(W5) and delta_TOP1 patterns already use, so the number
still has to sit next to its own label. The control that a changed number is caught still fires; all
13 negative controls pass.

### The machine ran out of disk, and the cascade refused

The first cascade after these edits failed with `OSError: [Errno 28] No space left on device`: drive
C: had 3.7 MB free, and the evidence lock's negative controls copy files into the temp directory.
Nothing was written in a half state -- the lock refused and the cascade stopped before anything
downstream ran. Re-run with `TMP`, `TEMP` and `TMPDIR` set to `D:\tmp-cellfate`, it passed. The
author's standing rule, that large temporary data belongs on D:, now applies to these runs too.

### Registered results

```text
  manuscript stage      GEN1_MANUSCRIPT_READY
  compliance checks     20 / 20
  negative controls     13 / 13
  numbers traced        17
  abstract              345 words, limit 350; the mirrored copy identical
  three --verify        EVIDENCE_INTACT, CLAIMS_INTACT, PACKAGE_INTACT
  full test suite       pytest's own exit code 0; 2,197 results: 2,196 passed, 1 skipped
                        (two more than Phase 4e: the two new renderer tests)
  render                exit 0; draft false; 0 FILL; Word export true; 0 code blocks wrapped
  PDF                   16 pages, up from 14, the tables taking the space
                        403199 bytes
                        45d0abd55f7a0ec423701ce9f8d398b1ad0e25706e8918d1b002a878c3cb012a
  LATER markers left    1: the APC-waiver line
```

### Digests

```text
  evidence  c84a4d7a2e2d254ed92e43ccf1f91d74da0ab57c8328b1be822d73e5c62ec350
       was  60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa
  claim     8c820412ba325cd053f0bc9d907d65ef1fd8e1a82aa93e0f947e85957be69796
       was  fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61
  package   4d08ef4a6d3aa61e744949acee84098ef11ecaf08a6d68ff5658f79a835b0f19
       was  c48419e9eecf0d37398a39546f0710e290f43632571811c8a5518ea12f81cc1d
```

The evidence digest moves because the renderer and the manuscript checker are in its inventory, and
the claim digest follows it. No claim changed: the claim lock's allowed set is untouched, and the
scanner is clean.

### Still open

```text
  Zenodo      the published preprint is v1, without any of this. A new version carries it: v1
              stays, v2 gets its own DOI, and the concept DOI keeps resolving to the newest
  Zenodo      on the archive record, add the related work: Is supplement to, the preprint DOI
  BMC         the submission itself, with the APC waiver requested at submission
  venue        the scope question stands: BMC Bioinformatics asks for computational methods, models
              and tools, so the cover letter has to lead with the method and the frozen tool, not
              with the reanalysis. The fallbacks named in SUBMISSION.md are unchanged
```
