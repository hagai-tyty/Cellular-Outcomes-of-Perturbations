# GEN-1 MANUSCRIPT + REPRODUCIBILITY PACKAGE

**Status** V1. The last Generation-1 stage.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Entry** `results/gen1_handoff_to_manuscript.json`, verdict `GEN1_CLAIMS_LOCKED`.

```text
  evidence digest  9245e605f6272aa809858d3f32dbe55ed53864df90e0b55750ab0a7d577da400
  claim digest     42e4e1a027022405d16035a2e1ea24f6a42209847ba5d7d39891ffbf6f4d404a
```

**Mandate** §9 of the frozen ship plan: `MANUSCRIPT + REPRODUCIBILITY PACKAGE -> PREPRINT /
SUBMISSION`, under the six obligations the claim lock handed forward:

```text
re-verify the evidence lock before submission
carry each allowed claim with its mandatory qualifiers
report p_perm as p < 0.001 (0 of 1,000), never a point estimate
state that independent biological replication is Generation 2, not a Gen-1 gate
scan any new abstract-level sentence against the same nine forbidden claims
bind to BOTH digests
```

---

# 0. What is different about this stage

Every previous stage checked something someone else wrote, or something a machine produced. This
one checks **prose I wrote myself**, which is the weakest position an instrument can be in — the
author and the reviewer are the same process.

The only defence is to make the checker refuse things mechanically, and to prove it can. So the
manuscript is not merely written and declared compliant: it is scanned by the claim lock's own
instrument, its every number is traced back to a locked artifact, and the checker is fired at
deliberately non-compliant copies of the manuscript to show it says no.

## 0.1 Authority

```text
MAY      write results/manuscript/*; write the manuscript and the package document
MAY      state fewer claims than the claim lock permits

MAY NOT  state more, or state one without its qualifiers
MAY NOT  edit any locked artifact, or either lock
MAY NOT  run an analysis, fit anything, or produce a new number
MAY NOT  proceed if either lock fails to verify
```

## 0.2 Compute budget

Seconds.

---

# 1. MS-A — preflight

Before a sentence is written or checked:

```text
the evidence lock verifies over all 54 artifacts, digest 455892ff...
the claim lock verifies, digest 0b3c7f03...
both digests equal the values in the handoff
the frozen ship-plan digest still holds
```

Any failure is `GEN1_MANUSCRIPT_REFUSED`.

---

# 2. MS-B — the manuscript

One document, `results/manuscript/MANUSCRIPT.md`. Required sections:

```text
TITLE            system-bounded; no claim beyond the lock
ABSTRACT         assembled from locked claims and their qualifiers
INTRODUCTION     the question, and why clone-level prospective evaluation is the hard part
DATA             Role B primary, Role A supporting, with accessions
METHODS          benchmark, models, endpoints, evaluation, the preregistered ranking test
RESULTS          the locked numbers, and only the locked numbers
THE TOOL         what it does, what it refuses
LIMITATIONS      every standing limitation, carried verbatim
WHAT THIS DOES NOT SHOW   the nine forbidden claims, stated as prohibitions
AVAILABILITY     both digests, accessions, the rebuild command, what is not vendored
GENERATION 2     what would actually test this, named as future work
```

`LIMITATIONS` and `WHAT THIS DOES NOT SHOW` are separate sections on purpose. The first is what the
result cannot support; the second is what may not be said. A reader should not have to infer either.

---

# 3. MS-C — compliance

The manuscript is scanned with the claim lock's **extended patterns under clause-scoped negation**
— the same instrument, not a copy of it.

```text
no forbidden claim appears unnegated, anywhere in the document
every allowed claim the manuscript makes appears with all five qualifiers present
p_perm appears as `p < 0.001` and the point estimate 0.000999 appears NOWHERE
the replication statement is present and says Generation 2
both digests are quoted, and both verify
```

## 3.1 The qualifier rule is checked, not trusted

It is not enough that the qualifiers appear somewhere in a long document. Each of the five must be
present, and the abstract — the part that travels alone — must carry the system, the vocabulary and
the outcome semantics by itself.

---

# 4. MS-D — every number traces to a locked artifact

Each numeric claim in the manuscript is matched against the locked source it came from, pinned to
the words around it, exactly as the evidence lock pins numbers in the records.

```text
delta_RANK, its CI, R(W1)/R(W4)/R(W5), the null p95, null max, the permutation count,
delta_TOP1, the eligible-clone count, the exclusions, the design columns,
the adversarial-refusal count, the Role-A power before and after audit
```

A number in the manuscript with no locked source is a `NUMBERS_UNTRACEABLE` refusal.

---

# 5. MS-E — the reproducibility package

`results/manuscript/REPRODUCIBILITY.md`. It must be executable-in-principle by a stranger:

```text
every command it names must exist -- the script file present, the --stage value accepted
every artifact it references must be in the evidence manifest
the environment is recorded: python version, key package versions, platform
runtimes are stated honestly, including the 10.7 h permutation run
what is NOT in the package is named: the 44 MB artifact and the raw sequencing data
the two verification commands are given first, before anything else
```

## 5.1 The commands are checked, not asserted

Each `--stage` value named in the package is checked against the executor's own argument parser. A
package that documents a flag the code does not accept is worse than no package.

---

# 6. MS-F — the checker must refuse

Negative controls, on **copies** in a scratch directory:

```text
PLANT A FORBIDDEN CLAIM      insert one sentence from the claim lock's adversarial corpus
                             -> must be caught, and name it
DROP A QUALIFIER             remove the outcome-semantics qualifier
                             -> must be caught
QUOTE p AS A POINT ESTIMATE  replace `p < 0.001` with `p = 0.000999`
                             -> must be caught
BREAK A NUMBER               change delta_RANK by one digit
                             -> must be caught
```

All four must fire. A checker that has never refused its own document is decoration.

---

# 7. Verdict

```text
GEN1_MANUSCRIPT_READY     MS-A through MS-F all pass
GEN1_MANUSCRIPT_REFUSED   a lock failed, a forbidden claim appeared, a qualifier was
                          missing, a number did not trace, a documented command did not
                          exist, or a negative control did not fire
```

`GEN1_MANUSCRIPT_READY` means the document is consistent with everything locked beneath it. It is
**not** a judgement that the science is good, that the writing is clear, or that a reviewer will
agree — none of which an instrument of this kind can assess.

---

# 8. Anti-rescue firewall

```text
no manuscript outcome reopens any earlier stage
no manuscript outcome changes a recorded number, a locked artifact, or either lock
the manuscript may state FEWER claims than the lock permits, never more
```

---

# 9. After this

Generation 1 is complete. Generation 2 — independent new-system biological replication,
unseen-condition transfer, calibration, and out-of-distribution validation — is future work and was
never a Generation-1 gate.
