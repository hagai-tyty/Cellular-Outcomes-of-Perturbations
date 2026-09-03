# gen1_EVIDENCE_LOCK_RECORD — freezing what the Generation-1 claims are made of

## Goal
Execute §9 of the frozen ship plan — *freeze benchmark, tool, OOF predictions, ranking verdict and
limitations* — plus the three things Stage 26 handed forward: hash every artifact and refuse if one
has moved, carry the nine forbidden claims into the claim lock unchanged, and record that
independent biological replication is Generation 2.

A file listing hashes is a manifest. It becomes a lock only when something **refuses** on the
strength of it. The failure this exists to prevent is unglamorous: the manuscript gets written
months from now against files that have quietly moved, and nothing notices because the only thing
checking them is a document asserting they were checked.

## Inputs
- `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md` @ canonical-LF `8da16fca...bced48`, FROZEN
- `results/stage26_handoff_to_evidence_lock.json` — `KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`
- every artifact named in §1 of the lock plan

## Files added
- `plans/(newer)practical plans/GEN1_EVIDENCE_LOCK_V1.md`
- `experiments/run_gen1_evidence_lock.py`
- `tests/test_gen1_evidence_lock.py`
- `results/evidence_lock/` — manifest, controls, numbers, claim-lock input, verdict, and
  `GEN1_EVIDENCE_LOCK.md`
- `results/gen1_handoff_to_claim_lock.json`

## What did NOT change
Nothing. The lock fitted nothing, regenerated nothing, and edited none of the 54 artifacts it
covers. Live verification after the run: 54 checked, 0 moved, 0 missing.

---

## Result

```text
  GEN1_EVIDENCE_LOCKED

  lock digest   a4f81d40d56760346b5c291a3a0fa0a84ca46a56843ca09d754bffea46e78e90
  artifacts     54
  runtime       11.8 s

  EL-A  inventory              54 artifacts across 9 classes, portable to a clone
  EL-B  chain of custody       3 chains, closed
  EL-C  verifier can refuse    3 negative controls, all fire
  EL-D  numbers                14 patterns, 0 disagreements, canary live
  EL-F  claim lock input       9 forbidden claims carried verbatim from the frozen plan
  --    live verification      54 checked, clean
```

The first lock issued was `99c35793...`, then `6c8420b5...` and `6ef9d81c...` as completeness gaps
and the self-audit below closed. Each is a real lock over a different set of bytes; the last one is
the one that ships. Earlier digests are recorded rather than erased.

### What is locked

```text
  benchmark            3    Stage-22 clone table, two frozen Stage-23 out-of-fold files
  predictions          1    the Stage-24 table Stage 25 actually consumed
  tool                10    API, CLI, artifact + metadata, model card, schema, example clone
  verdicts             6    Stage-25, Stage-26, and every inter-stage handoff
  limitations          1    GEN1_SCOPE_LIMIT.md
  supporting_role_A    5    Rewind verdict, confirmation, power, and the power AUDIT
  protocol             7    the frozen plans and the Stage-23.5 protocol JSON
  records              8    the prose the manuscript will be written from
  code                13    every executor and test that produced or constrains a locked artifact
```

`lock_digest` is a SHA-256 over `path  sha256` lines, sorted by path and LF-joined. One number that
names the whole evidence base, reproducible from the manifest by anyone.

### The chain of custody is closed

An artifact's hash must be the same number **everywhere a stage recorded it**. If Stage 24 handed
Stage 25 a table hashed X and the lock hashes Y, the analysis did not consume what is being locked,
and no amount of documentation repairs that.

```text
  out-of-fold table   identical across 4 independent records
                      on disk == 24F freeze == 24->25 handoff == 25A observed
  model artifact      identical across 4 independent records
  ship-plan digest    identical across 6 independent records
                      on disk == frozen value == 23.5 protocol == 24->25 == 25C == 26E
```

Nobody had checked this before. Each stage asserted its own inputs; none had checked that the
numbers *agreed across stages*. They do.

### The verifier can refuse, and that was tested before any lock was issued

```text
  INTACT     change nothing            -> reports clean
  MOVED      flip one bit, one byte    -> caught, and names the file
  MISSING    delete one file           -> caught, and names the file
```

All three run on **copies in a scratch directory**. Mutating a locked artifact to test the lock
would be a spectacular own goal. A verifier that has never failed is an assumption, not a check.

The verifier also **locks itself**: `run_gen1_evidence_lock.py` is in its own manifest. This is not
circular — the manifest is written after hashing, so a later edit makes `--verify` report the
verifier moved. Tested directly: appending one comment line to a copy is detected.

### The numbers the manuscript will use

```text
  eligible clones        892
  R(W1) / R(W4) / R(W5)  0.692654 / 0.692176 / 0.743781
  delta_RANK             +0.051605   CI95 [+0.037197, +0.065571]
  null p95               0.008672
  permutation            0 of 1000 draws reached the observed value
  p_perm                 0.000999    report as p < 0.001 (0 of 1,000), never as a point estimate
  delta_TOP1             +0.115471
  adversarial refusals   56 of 56
  design columns         309
```

Every one of these was checked **against the JSON it came from**, in the record where the
manuscript will read it. The records are prose typed by hand; 13 substrings, 0 disagreements. The
check also asserts `p_perm` is reported as a floor and equals the finite-sample formula exactly.

---

## Bugs found — two, both mine, both in the lock itself

**1. The forbidden-claims check counted the code-fence tag as a tenth claim.**
`ship.split("## 3.5 ...")[1].split("```")[1]` returns the fence contents *including* the `text`
language tag, so the parse produced ten entries and refused a list that was in fact identical. The
nine were correct all along. Fixed by dropping the tag line. Same class as the Stage-26 scanner
tripping on its own description — reading a document mechanically means reading its markup too.

**2. The inventory gated on "untracked" when it meant "gitignored".**
`GEN1_EVIDENCE_LOCK_V1.md` was untracked for the ordinary reason that it had not been committed
yet, and the lock refused. But a file awaiting its first commit is fine; a file git will **never**
carry is the actual hole in the reproducibility package. The check now gates on `git check-ignore`
and reports pending-commit files separately, without gating.

Both were caught by the lock refusing on its first run, which is the behaviour the stage is for.

## Two completeness gaps I had left open, and closed

The first run locked 34 artifacts. Re-reading the plan against what the manuscript will actually
say surfaced two omissions:

- **Role-A supporting evidence was absent.** §3.3 of the frozen plan permits one supporting
  sentence about Rewind, and the standing limitation *"gate 18.3 FAILED at 0.64 (audited ~0.45)"*
  is quoted in the scope document and in the Stage-25 record. A claim the manuscript makes is a
  claim whose evidence must be locked — including, and especially, the audit that lowered the power
  from 0.64 to 0.45.
- **The records were absent.** EL-D checks headline numbers against the records, so a record free to
  change after the lock makes that check meaningless.

Locking both took the manifest from 34 to 54 artifacts and changed the lock digest, which is
correct: a more complete lock is a different lock.

## What this lock does NOT contain

```text
  stage24_w5_artifact.npz   44 MB, gitignored. A fresh clone does NOT contain it. Its
                            hash IS locked, and the rebuild is one command:
                              python experiments/run_stage24_gen1_tool.py --stage 24c
                            Anyone verifying this lock from a clone must rebuild first.

  raw sequencing data       GSE279162 (WM989), GSE227151 (Rewind). Accessions locked,
                            bytes not vendored.
```

Naming a gap is not closing it. Both stay open and both are in the lock document, not only here.

## A note on the two digest forms

The ship plan appears twice with different numbers and they do not disagree. The manifest hashes
raw bytes — what a file manifest must do. The frozen protocol identity `8da16fca...` is a
*canonical-LF* digest of the same file, CRLF normalised to LF so the protocol has one identity on
every platform. The plan is `59f22e9a...` as a file and `8da16fca...` as a protocol; both are
checked, and the lock document now says so rather than leaving a reader to reconcile them.

---

## Self-audit before the claim lock — four more findings

**A. The lock was verifiable only on this machine.** This is the serious one, and it is the exact
failure the lock exists to prevent, inverted.

```text
  the repo runs core.autocrlf=true
  -> a text file's bytes in the working tree are NOT its bytes in the repository
  -> measured: 28 of 53 tracked artifacts differed between the two
  -> a clone on Linux, or on Windows with different settings, computes a different
     digest and the lock REFUSES -- for everyone except me
```

The lock's whole audience is someone who did not build it, on a machine that is not this one, and
for that audience the lock cried wolf universally. Fixed by hashing text **canonical-LF** and
binary raw — the same rule this project already uses to give a frozen protocol one identity on
every platform, and the rule `.gitattributes` had already written down for `results/**` for exactly
this reason. Measured after the fix: **28 disagreements → 0**.

This is now a **gate**, not something I check by hand. EL-A compares every committed artifact's
locked hash against the bytes git actually stores, skipping and reporting files with uncommitted
edits. Only two artifacts were skipped at lock time — this executor and its plan, both being edited
in this pass.

**B. A number could be satisfied by an accident.** EL-D used bare substrings, and `"56"` occurs
nine times in the Stage-26 record — inside `SHA-256`, inside `frozen_24F_sha256`. A record that
never stated the refusal count would have passed. Every number is now pinned to the words around it
(`eligible clones 892`, `56 / 56 adversarial`, `design columns 309 =`), and each pattern is
controlled by perturbing its value by one digit and requiring the pattern to then find **nothing**.
Matching the right value and rejecting a wrong one, together, is what makes it a live check.

**C. The canary placeholder was regex-escaped, so it never fired.** The templates were written as
`\{n\}` — an escaped brace — while the substitution looked for `{n}`. Nothing was ever replaced, so
all fourteen checks failed at once, including ones that trivially hold. Caught on the first run
precisely because the failure was total; a placeholder that half-worked would have been far worse.
The placeholder is now `@@`, which has no regex meaning of its own.

**D. The plan said seven classes; the code had nine.** `supporting_role_A` and `records` were added
to the module during the completeness pass below without going back to §1. Same class of drift as
Stage 26's finding A. The plan now lists all nine, states which artifacts are deliberately excluded,
and carries the hashing rule as §1.1 and the number-pinning rule as §4.1.

## Tests
- 27 evidence-lock contracts, 0 skipped (22 before the audit; 5 added for findings A-C)
- Stage-25, Stage-26 and the path-convention suites re-run green

---

## Scientific interpretation

**Proves:** the Generation-1 evidence base is in a known, verifiable state. 54 artifacts are
hashed under one digest; the out-of-fold table, the model artifact and the protocol digest are
identical across every stage that recorded them; and the verifier that enforces this has been shown
to catch a one-bit change, a deleted file, and an edit to itself.

**Does NOT prove:**
- **anything scientific.** Locking grants no claim. It fixes what the existing claims are made of.
  No number moved and no analysis ran.
- **that the evidence is sufficient.** A lock says *these files, unchanged*. It says nothing about
  whether they support the claim — that was Stage 25's job, under its own preregistration.
- **that the package is self-contained.** It is not: the 44 MB model artifact and the raw
  sequencing data are both outside the repository, and both are named above.
- **that a clone reproduces the analysis.** The lock proves a clone can reproduce the *hashes*.
  Re-running Stage 25's 1,000-draw null on another machine is a separate claim this does not make.
- **that the records are correct.** EL-D checks that the numbers in them match their JSON sources.
  It does not check that the prose around those numbers is right.

## Next action
`GEN-1 CLAIM LOCK`. Consume `GEN1_CLAIM_LOCK_INPUT.json` unchanged — it **may narrow, never
widen** — re-verify the lock digest before writing an abstract-level claim, and state that
independent biological replication is Generation 2 and not a Generation-1 gate. Then manuscript and
reproducibility package.

Verify this lock at any time:

```text
  python experiments/run_gen1_evidence_lock.py --verify
```

No lock outcome reopens an earlier stage.
