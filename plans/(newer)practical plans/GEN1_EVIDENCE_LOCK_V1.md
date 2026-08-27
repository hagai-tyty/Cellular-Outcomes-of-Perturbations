# GEN-1 EVIDENCE LOCK

**Status** V1.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Entry** `results/stage26_handoff_to_evidence_lock.json`, verdict
`KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`.

**Mandate** §9 of the frozen ship plan, verbatim:

```text
GEN-1 EVIDENCE LOCK
  freeze benchmark, tool, OOF predictions, ranking verdict and limitations
```

plus the three things Stage 26 handed forward:

```text
hash every artifact and refuse to proceed if one has moved
carry the nine forbidden claims into the claim lock unchanged
record that independent biological replication is Generation 2, not a Gen-1 gate
```

---

# 0. What a lock is, and what a manifest is

A file listing hashes is a manifest. It becomes a lock only when something **refuses** on the
strength of it. So the deliverable here is not a list — it is a verifier that fails loudly, plus
proof that it can fail.

The failure this stage exists to prevent is specific and unglamorous: the manuscript is written
months from now against files that have quietly moved, and nobody notices because the only thing
checking them is a document that says they were checked.

## 0.1 Authority

```text
MAY      read every artifact; hash it; write results/evidence_lock/*
MAY      record a gap honestly and proceed with it named

MAY NOT  fit, refit, re-run or regenerate any analysis
MAY NOT  edit any artifact it locks -- including to fix a typo
MAY NOT  soften a limitation, widen a claim, or drop a forbidden claim
MAY NOT  proceed past a moved artifact under any justification
```

If an artifact has moved, the lock **refuses**. It does not re-hash and carry on. The repair is to
find out why it moved, not to record the new value.

## 0.2 Compute budget

Seconds. Hashing a 44 MB file and some text.

---

# 1. EL-A — the inventory

Nine classes. Every artifact carries `path`, `sha256`, how it was hashed, `bytes`, and whether git
tracks or ignores it.

```text
BENCHMARK          the Stage-22 clone table
                   the two frozen Stage-23 out-of-fold files

PREDICTIONS        the Stage-24 out-of-fold table Stage 25 actually consumed

TOOL               the prediction API and the CLI
                   the serialized model artifact and its metadata
                   the model card, the io schema, the example clone

VERDICTS           the Stage-25 ranking verdict
                   the Stage-26 scope-lock verdict
                   every inter-stage handoff

LIMITATIONS        GEN1_SCOPE_LIMIT.md, the authoritative scope document

SUPPORTING_ROLE_A  the Rewind verdict, confirmation, power, and the power AUDIT

PROTOCOL           the frozen plans and the Stage-23.5 protocol JSON

RECORDS            the stage records the manuscript will be written from

CODE               every executor that produced a locked artifact
                   every test file that constrains one
```

`CODE` is in the lock because a result without the code that made it is not reproducible evidence,
and because a silent edit to an executor is exactly as damaging as a silent edit to a result. The
lock's own executor is included: the manifest is written after hashing, so a later edit makes
`--verify` report the verifier itself moved.

`SUPPORTING_ROLE_A` is in the lock because §3.3 of the frozen plan permits one supporting sentence
about Rewind, and the standing limitation *"gate 18.3 FAILED at 0.64 (audited ~0.45)"* is quoted in
the scope document. A claim the manuscript makes is a claim whose evidence must be locked —
including, and especially, the audit that lowered that number.

`RECORDS` is in the lock because EL-D checks headline numbers against the records. A record free to
change afterwards makes that check meaningless.

Deliberately **not** locked, being outputs of this run: the manifest, the lock document, the lock
verdict, and the evidence-lock record.

## 1.1 How an artifact is hashed

```text
BINARY  .npz .npy        raw bytes
TEXT    everything else  canonical-LF: CRLF normalised to LF before hashing
```

Raw bytes are the wrong unit for text here and this is not a stylistic choice. The repository runs
`core.autocrlf=true`, so a text file's bytes in the working tree are not its bytes in the
repository — measured: 28 of 53 tracked artifacts differ between the two. A lock built on raw text
bytes is a property of one working tree on one platform, and would refuse for **everyone who cloned
the repository**, which is precisely the audience a lock exists to serve.

Canonical-LF is the same rule this project already uses to give a frozen protocol one identity on
every platform. Binary stays raw: normalising a float array would mangle any `0D 0A` byte pair that
happens to fall inside the data.

The lock must record which rule it applied to each artifact, so nobody has to guess later.

## 1.1 The lock digest

```text
lock_digest = SHA-256 over the canonical manifest
              (path + sha256, sorted by path, LF-joined)
```

One number that names the entire Generation-1 evidence base. It goes in the manuscript.

---

# 2. EL-B — chain of custody

An artifact's hash must be the **same number** everywhere a stage recorded it. If Stage 24 handed
Stage 25 a table hashed `X` and the lock now hashes `Y`, the analysis did not consume what is being
locked, and no amount of documentation fixes that.

```text
OUT-OF-FOLD TABLE     on disk == 24F freeze == 24->25 handoff == 25A observed
MODEL ARTIFACT        on disk == 24F freeze == 24->25 handoff == its own metadata
SHIP PLAN DIGEST      on disk == 23.5 protocol == 24->25 handoff == 25C == 26E
```

All three chains gate. A disagreement is `EVIDENCE_LOCK_REFUSED`, not a note.

---

# 3. EL-C — the verifier must be able to refuse

A verifier that has never failed is an assumption. Three negative controls, on **copies** in a
scratch directory, never on the artifacts themselves:

```text
MOVED     flip one byte in a copy      -> the verifier must report it, naming the file
MISSING   delete a copy                -> the verifier must report it, naming the file
INTACT    change nothing               -> the verifier must report clean
```

If any of the three does not behave, the lock is `REFUSED` regardless of whether the real artifacts
are intact. An instrument that cannot fail cannot pass anything either.

---

# 4. EL-D — the numbers

The manuscript will be written from the records, and the records are prose typed by hand. Every
headline number in them is therefore checked against the machine-readable source it came from.

```text
delta_RANK, its CI, R(W1), R(W4), R(W5), the null p95, the permutation count,
p_perm, delta_TOP1, the eligible-clone count, the adversarial-refusal count,
the design-column count, and both verdict strings
```

Each must appear in the record **formatted exactly as the JSON holds it**. A transcription slip is
a `NUMBERS_DISAGREE` refusal.

## 4.1 Pinned to meaning, not to a substring

A bare substring check is satisfied by an accident. `"56"` occurs inside `SHA-256` and inside
`frozen_24F_sha256`, so a record that never states the refusal count would still pass. Every number
is therefore matched by a pattern that carries the words around it — `eligible clones 892`,
`56 / 56 adversarial`, `design columns 309 =`.

And the patterns are themselves controlled: each value is perturbed by one digit and the pattern
must then find **nothing**. A pattern keyed to the surrounding words alone would keep matching and
would prove nothing about the number.

---

# 5. EL-E — what is NOT in the lock

Recorded plainly rather than quietly omitted.

```text
THE MODEL ARTIFACT IS GITIGNORED
  stage24_w5_artifact.npz is 44 MB and is not in the repository. A fresh clone does
  not contain it. Its hash is locked and its exact rebuild command is recorded, but
  anyone verifying the lock from a clone must rebuild it first.

RAW DATA IS EXTERNAL
  GSE279162 (WM989) and GSE227151 (Rewind) are not vendored. Accessions are locked;
  bytes are not.
```

Naming a gap is not the same as closing it. Both stay open, and the record says so.

---

# 6. EL-F — the claim lock input

Emitted as machine-readable JSON for the claim lock, carried through **unchanged**:

```text
the nine forbidden claims, verbatim from §3.5 of the frozen ship plan
the allowed primary claim, verbatim from §3.1
the allowed ranking claim, verbatim from §3.2, in the form the Stage-25 verdict selects
the allowed supporting claim, verbatim from §3.3
independent biological replication is GENERATION 2 and is not a Generation-1 gate
```

The claim lock may narrow these. It may not widen them.

---

# 7. Verdict

```text
GEN1_EVIDENCE_LOCKED        EL-A through EL-F all pass
GEN1_EVIDENCE_LOCK_REFUSED  any artifact moved, any chain broken, any control failed,
                            any number disagreed, any input missing
```

A refusal is not a scientific result. It means the evidence base is not in the state the records
claim, and the discrepancy is investigated before anything else happens.

---

# 8. Anti-rescue firewall

```text
no lock outcome reopens Stage 25, Stage 26, or any earlier stage
no lock outcome changes a recorded number
no lock outcome authorizes new data, a new condition, or a new model
locking grants no claim; it fixes what the existing claims are made of
```

---

# 9. Handoff

```text
GEN1_EVIDENCE_LOCKED  ->  GEN-1 CLAIM LOCK  ->  MANUSCRIPT + REPRODUCIBILITY PACKAGE
```

Generation 2 — new-system biological replication, unseen-treatment transfer, broad calibration and
OOD validation — remains future work and is not a Generation-1 gate.
