# GEN-1 CLAIM LOCK

**Status** V1.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Entry** `results/gen1_handoff_to_claim_lock.json`, verdict `GEN1_EVIDENCE_LOCKED`,
lock digest `6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081`.

**Mandate** §9 of the frozen ship plan, verbatim:

```text
GEN-1 CLAIM LOCK
  freeze abstract-level allowed/forbidden claims
```

plus the three things the evidence lock handed forward:

```text
consume GEN1_CLAIM_LOCK_INPUT.json unchanged; it may narrow, never widen
re-verify the lock digest before writing an abstract-level claim
state that independent biological replication is Generation 2, not a Gen-1 gate
```

---

# 0. What this stage is for

Every previous stage produced a number. This one produces **sentences**, and sentences are where a
project of this kind actually fails. Not in the statistics — in the abstract, where "six observed
experimental conditions in one melanoma line" becomes "treatments in cancer", and a detection proxy
becomes "response".

The evidence lock fixed what the claims are made of. This fixes what may be said about it.

## 0.1 The failure mode this exists to stop

An abstract assembled entirely from permitted fragments can still read as a broader claim than any
of them. So it is not enough to list allowed sentences: the lock must also demonstrate where the
boundary is, by naming the tempting sentences that are **not** permitted and showing that the
instrument catches them.

## 0.2 Authority

```text
MAY      narrow a claim; add a scope qualifier; refuse a claim entirely
MAY      write results/claim_lock/*

MAY NOT  widen any claim from the evidence lock's input
MAY NOT  add a claim not entailed by locked evidence
MAY NOT  drop, soften, or reword any of the nine forbidden claims
MAY NOT  proceed if the evidence lock does not verify
MAY NOT  run any analysis, fit anything, or produce a new number
```

## 0.3 Compute budget

Seconds. This stage hashes, scans text, and writes.

---

# 1. CL-A — the evidence must verify first

The very first thing, before a single sentence is written:

```text
re-run the evidence-lock verifier over every locked artifact
the lock digest must equal 6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081
```

A moved artifact is `GEN1_CLAIM_LOCK_REFUSED`. Writing claims against evidence that has shifted is
the failure the previous stage existed to prevent, and it would be absurd to commit it here.

---

# 2. CL-B — the allowed claims

Each locked claim carries three things. A claim missing any of them is not a claim.

```text
TEXT       the exact sentence permitted at abstract level
EVIDENCE   the locked artifacts that support it, by path
NUMBERS    the locked figures it rests on, or "none"
```

The allowed set comes from `GEN1_CLAIM_LOCK_INPUT.json` — primary, ranking, supporting — and may
only be **narrowed** here. Narrowing means adding scope, not removing it.

## 2.1 Scope qualifiers are mandatory, not decorative

Every allowed claim is bound to a qualifier that must travel with it wherever it appears:

```text
system      one BRAF-V600E melanoma cell line, WM989, 1,401 lineage-traced clones
vocabulary  six observed experimental conditions, named
outcome     an observed post-treatment clone DETECTION proxy, not death or response
evaluation  clone-held-out, frozen before any result
replication NONE. Generation 2.
```

A claim quoted without its qualifier is a claim the lock did not grant.

---

# 3. CL-C — the forbidden nine, carried unchanged

The nine from §3.5 of the frozen ship plan pass through **verbatim**, parsed from the plan itself
rather than from a copy. Count must be nine. Text must match exactly.

The claim lock may add to this list. It may not subtract from it or reword an entry.

---

# 4. CL-D — the adversarial abstract corpus

The claim-lock equivalent of firing `Vemurafenib` at the tool.

A corpus of sentences a well-meaning author would plausibly write, each of which is **forbidden**.
Every one must be caught by the same scanner Stage 26 used on the shipped tool — the same nine
claim patterns, the same twelve negation tokens, the same canary. Reusing the instrument is
deliberate: the claim text and the tool are held to one standard, not two.

**Fifteen sentences, verbatim.** Not paraphrases: the executor asserts that every sentence it runs
appears in this list exactly, so the corpus cannot drift from the plan that declared it.

```text
GENERALISATION   "The model generalises to new treatments."
GENERALISATION   "Our approach works for any drug."
CROSS_SYSTEM     "We predict outcomes in cancer cells."
CROSS_SYSTEM     "Validated on other cell lines."
CROSS_SYSTEM     "The model predicts patient response."
CLINICAL         "The tool identifies the best treatment for each clone."
CLINICAL         "This supports clinical decision-making."
CLINICAL         "A therapeutic recommendation tool for melanoma."
CAUSAL           "We estimate the causal effect of each treatment."
CALIBRATION      "The model outputs a calibrated probability of death."
REPLICATION      "Independently replicated in an external cohort."
UNIFORMITY       "The interaction helps uniformly across all six conditions."
ROLE_A           "The finding was confirmed in a second system."
ROLE_A           "Rewind confirms the result."
SINGLE_CELL      "The tool scores a single cell."
```

Each is declared here before it is run, so the corpus cannot be trimmed after seeing which entries
the scanner misses. **A missed entry is a failed stage**, not a note — it means the instrument
cannot see a claim the manuscript could make.

## 4.1 The near-miss table

For each forbidden sentence, the lock records the **nearest permitted phrasing**. A writer who is
told only "no" learns nothing; a writer shown the permitted neighbour learns the boundary.

Every permitted phrasing is itself scanned. A table that swapped one refused sentence for another
would be worse than no table.

## 4.2 The prose extension — and why it was needed

Fired at the unmodified Stage-26 patterns, **three of the fifteen walked straight through**:

```text
  "We predict outcomes in cancer cells."          `cancer` is REVIEW_ONLY in Stage 26,
                                                  reported but not gating
  "The finding was confirmed in a second system." no pattern for confirmation without the
                                                  word Rewind or Role-A beside it
  "The tool scores a single cell."                the pattern reads `score a single cell`;
                                                  the sentence says `scores`
```

Stage 26's patterns were calibrated on the tool's own surfaces — code, a model card, a JSON schema.
The tool's documentation would never say "cancer cells". An abstract easily would.

So the claim lock **adds** patterns for prose. Three rules govern the addition:

```text
ADDED, NEVER SUBSTITUTED     the nine claims and the twelve negation tokens are unchanged
THE LOCKED MODULE IS NOT EDITED
                             run_stage26_scope_lock.py is a locked evidence artifact;
                             extending it in place would break the evidence lock and
                             rewrite a completed stage's basis
TURNED BACK ON THE OLD TEXT  the extended set is re-run over every surface Stage 26 already
                             passed. A stricter instrument pointed only at new text is not
                             an instrument, and if it finds something in the shipped tool
                             that is a real finding, reported, not buried
```

The extension must also still refuse to fire on a negation — checked against every sentence in the
corpus, negated.

## 4.3 Negation is scoped to the clause, not to a window

Stage 26 excuses a forbidden phrase when a negation token appears anywhere within ±160 characters.
On code and a model card that is fine. On prose it is close to toothless, because prose is full of
legitimate negations and proximity cannot tell which clause they govern:

```text
  "The model is not calibrated for abundance, and outputs a calibrated probability of death."
  "We make no claim about dosing; the tool identifies the best treatment for each clone."
  "This was not replicated internally, but was independently replicated in an external cohort."
```

Each makes a plainly forbidden claim. **All three pass the window rule.**

Prose is therefore scanned clause-scoped: a negation excuses a hit only in the same clause, split
on `. ; :` and on a comma before a coordinating conjunction. Same twelve tokens, same nine claims —
a tighter scope, not a different instrument.

```text
WHITESPACE IS NORMALISED FIRST
  a newline is not a clause boundary. Treating it as one split "...and not a /
  clinical recommendation." in the shipped predictor's docstring and reported a
  negated sentence as a forbidden claim.

A WORD MAY NEGATE ITSELF
  `uncalibrated` contains `calibrated`, and the Stage-26 pattern has no word
  boundary. A match immediately preceded by un-/non-/de- is excused. The prefix
  must be contiguous, so "run calibrated" is untouched.

STRUCTURED TEXT KEEPS THE WINDOW RULE
  clause scoping is right for prose and wrong for JSON, where it splits a key from
  the value that negates it: `"calibrated_probability": "NEVER emitted"`. The locked
  surfaces are gated with the extended patterns under Stage 26's window rule, and the
  clause-scoped reading is reported beside it. Only the negation SCOPE differs, matched
  to the kind of text being read.
```

---

# 5. CL-E — the abstract-level statement

One worked example: a complete, permitted, abstract-level paragraph assembled only from locked
claims and their mandatory qualifiers. It is scanned by the same instrument and must come back
clean.

This is a demonstration, not a mandate. The manuscript may write its own abstract; it may not write
one this instrument would refuse.

---

# 5.1 The claims must have an identity of their own

The evidence lock hashes 62 artifacts, and **none of them is this stage's output**. That is the
correct layering — the claim lock is downstream, and re-locking to include it would make CL-A
circular, since CL-A pins the evidence digest. But it leaves the document the manuscript is
actually written from with no identity at all: it could change and nothing would notice, which is
the exact failure the layer below exists to prevent.

So this stage computes its own digest, by the same canonical-LF rule, over the five files that
constitute it:

```text
  the plan, the executor, the contracts, the claims document, the verdict JSON
  -> claim_digest, carried in the verdict and in the handoff to the manuscript
```

`--verify` re-checks it. The manuscript binds to two numbers: the evidence digest for what the
claims are made of, and the claim digest for what may be said about it.

---

# 6. Verdict

```text
GEN1_CLAIMS_LOCKED          CL-A through CL-E all pass
GEN1_CLAIM_LOCK_REFUSED     evidence moved, a claim widened, a forbidden entry changed,
                            an adversarial sentence went undetected, or the worked
                            abstract failed its own scan
```

---

# 7. Anti-rescue firewall

```text
no claim-lock outcome reopens any earlier stage
no claim-lock outcome changes a recorded number or a locked artifact
locking a claim grants nothing; it fixes the ceiling of what may be said
the ceiling may be lowered later. It may not be raised.
```

---

# 8. Handoff

```text
GEN1_CLAIMS_LOCKED  ->  MANUSCRIPT + REPRODUCIBILITY PACKAGE
```

Generation 2 — new-system biological replication, unseen-treatment transfer, broad calibration and
OOD validation — remains future work and is not a Generation-1 gate.

---

# 9. Amendment V1.1 — 2026-09-03

The two numbers this plan pins are re-pinned here rather than quietly overwritten, so that what
V1 said stays visible.

**Entry lock digest.** V1 was written against evidence lock digest
`e206bfd37c5a93998a773b8bd058eac5e5e144cd2a8ee5d78e9907911a956bc5`.
It is now `6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081`.
Nothing about the evidence changed: the digest moved because the lock's own manifest recorded
each artifact's `st_size`, which counts a CRLF as two bytes. `results/**` is `text eol=lf`, so a
file written by a stage on Windows and the same file in a fresh checkout have identical content
and different sizes -- the manifest disagreed with itself across checkouts and re-dirtied the
working tree on every run. Sizes are now measured on the same canonical-LF content that is
hashed. The artifact hashes themselves were already canonical-LF and did not move.

**Artifact count.** V1 said 54. It is 62: the three figures, the source-data export, the release
bundle builder and the two upstream Gen-1 records were added to the lock after V1 was written.
The count is no longer written into the CL-A check label, which now derives it from the manifest
-- a hard-coded count in a label that describes a verification is a claim that can go stale
while the verification stays correct, which is what happened here.

Neither amendment relaxes a gate. CL-A still refuses if a single covered byte has moved.
