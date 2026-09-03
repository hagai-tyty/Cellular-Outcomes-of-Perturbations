# gen1_CLAIM_LOCK_RECORD — freezing what may be said

## Goal
Execute §9 of the frozen ship plan — *freeze abstract-level allowed/forbidden claims* — plus the
three things the evidence lock handed forward: consume `GEN1_CLAIM_LOCK_INPUT.json` unchanged (it
may narrow, never widen), re-verify the lock digest before writing a claim, and state that
independent biological replication is Generation 2.

Every previous stage produced a number. This one produces **sentences**, and sentences are where a
project of this kind actually fails. Not in the statistics — in the abstract, where "six observed
experimental conditions in one melanoma line" becomes "treatments in cancer", and a detection proxy
becomes "response".

## Inputs
- `results/gen1_handoff_to_claim_lock.json` — `GEN1_EVIDENCE_LOCKED`
- evidence lock digest `a4f81d40d56760346b5c291a3a0fa0a84ca46a56843ca09d754bffea46e78e90`
- `results/evidence_lock/GEN1_CLAIM_LOCK_INPUT.json`
- `experiments/run_stage26_scope_lock.py` — imported for its scanner, **not modified**

## Files added
- `plans/(newer)practical plans/GEN1_CLAIM_LOCK_V1.md`
- `experiments/run_gen1_claim_lock.py`
- `tests/test_gen1_claim_lock.py`
- `results/claim_lock/GEN1_CLAIM_LOCK.json`, `GEN1_CLAIMS.md`, `claim_lock_adversarial.json`
- `results/gen1_handoff_to_manuscript.json`

## What did NOT change
No locked artifact. No number. No analysis ran and nothing was fitted. The evidence lock still
verifies clean over all 54 artifacts, before and after.

---

## Result

```text
  GEN1_CLAIMS_LOCKED

  evidence lock digest   a4f81d40d56760346b5c291a3a0fa0a84ca46a56843ca09d754bffea46e78e90
  claim digest           f69bd7f682ab3738ce73171ddb148af7e786f5813f717554f4932e1137cc9817
  allowed claims          3    each bound to its evidence and its mandatory qualifiers
  forbidden claims        9    verbatim from §3.5, parsed from the plan not from a copy
  adversarial sentences  15 of 15 caught
  worked abstract         passes the same instrument, all qualifiers present
```

The manuscript binds to **two** numbers: the evidence digest for what the claims are made of, and
the claim digest for what may be said about it. Both have a `--verify`.

The claim digest was `23ea00b8...` when this stage closed. It moved to `0b3c7f03...` in the
final close-out pass, when a ruff `N802` on a helper in `tests/test_gen1_claim_lock.py` was
fixed -- the contracts are one of the five files the digest covers. No claim, qualifier or
forbidden entry changed; the verdict was re-derived and is unchanged. The earlier value is
recorded rather than erased.

### The three allowed claims

Verbatim from the evidence lock's input — the lock may narrow, never widen, and it did not narrow.

```text
  PRIMARY      pretreatment expression carries treatment-specific information about future
               clonal detection, beyond treatment identity and clone abundance
  RANKING      the interaction model improves clone-specific ordering of the six conditions
               delta_RANK +0.051605, CI95 [+0.037197, +0.065571], 0 of 1,000 draws
  SUPPORTING   Rewind gave positive but underpowered evidence -- and must always travel with
               the word SUPPORTING, because its own confirmation gate FAILED
```

### The qualifiers are mandatory, not decorative

A claim quoted without them is a claim the lock did not grant.

```text
  system       one BRAF-V600E melanoma cell line, WM989, 1,401 lineage-traced clones
  vocabulary   six named observed experimental conditions
  outcome      an observed post-treatment clone DETECTION proxy -- not death, sensitivity,
               resistance, clinical response or patient benefit
  evaluation   clone-held-out, five outer folds, frozen before any result existed
  replication  NONE. Generation 2.
```

---

## The finding: the Stage-26 scanner was not enough for prose

This is the reason the stage was worth running rather than writing.

The adversarial corpus — fifteen sentences a well-meaning author would plausibly write, all
forbidden, all declared in the plan **before** being run — was fired at the unmodified Stage-26
scanner, the same instrument that passed the shipped tool. **Three walked straight through:**

```text
  "We predict outcomes in cancer cells."           `cancer` is REVIEW_ONLY in Stage 26 --
                                                   reported, deliberately not gating
  "The finding was confirmed in a second system."  no pattern for confirmation without the
                                                   word Rewind or Role-A beside it
  "The tool scores a single cell."                 the pattern reads `score a single cell`;
                                                   the sentence says `scores`
```

Stage 26's patterns were calibrated on the tool's own surfaces: code, a model card, a JSON schema.
**The tool's documentation would never say "cancer cells". An abstract easily would.** That is not
a defect in Stage 26 — it is a different surface with different failure modes, and it would have
been invisible without an adversarial corpus written for prose.

### How the extension was made, and what it was not allowed to do

```text
  ADDED, NEVER SUBSTITUTED   the nine claims and the twelve negation tokens are unchanged;
                             Stage 26's patterns survive intact as a prefix of every list
  THE LOCKED MODULE UNTOUCHED
                             run_stage26_scope_lock.py is a locked evidence artifact.
                             Editing it in place would break the evidence lock and rewrite
                             a completed stage's basis. It is imported, never modified.
  TURNED BACK ON THE OLD TEXT
                             the extended set was re-run over every surface Stage 26 had
                             already passed
```

That last one mattered and could have gone badly. **It came back empty** — the shipped tool, model
card, schema, scope document and Stage-25 verdict contain no forbidden claim even under the
stricter prose patterns. Stage 26's verdict holds under an instrument it was never tested against.

The extension also had to keep refusing to fire on a negation, checked against all fifteen corpus
sentences negated. It does.

## The near-miss table

A writer told only "no" learns nothing. Every forbidden sentence is paired with its nearest
permitted phrasing, and **every permitted phrasing was itself scanned clean** — otherwise the table
would teach one refused sentence to be swapped for another.

```text
  FORBIDDEN   The tool identifies the best treatment for each clone.
  PERMITTED   The tool orders six observed experimental conditions for each clone; this is
              experimental-condition selection, not treatment recommendation.

  FORBIDDEN   Rewind confirms the result.
  PERMITTED   Rewind supports the result without confirming it; gate 18.3 FAILED at 0.64,
              audited ~0.45.
```

---

## Self-audit — three more findings, one of them the sharpest in the stage

**A. The negation rule was proximity, not scope — and on prose it was close to toothless.**

Stage 26 excuses a forbidden phrase when a negation token appears anywhere within ±160 characters.
That is fine on code and a model card. Prose is full of legitimate negations, and proximity cannot
tell which clause they govern. Three sentences, each making a plainly forbidden claim, each
containing an unrelated negation:

```text
  "The model is not calibrated for abundance, and outputs a calibrated probability of death."
  "We make no claim about dosing; the tool identifies the best treatment for each clone."
  "This was not replicated internally, but was independently replicated in an external cohort."

  all three PASS the window rule
```

This mattered directly: the worked abstract is full of disclaimers, so almost anything in it would
have been excused. Prose is now scanned **clause-scoped** — a negation excuses a hit only in the
same clause, split on `. ; :` and on a comma before a coordinating conjunction. Same twelve tokens,
same nine claims; a tighter scope, not a different instrument.

**The abstract was re-checked under the strict rule and came back with zero hits** — it was
genuinely clean, not merely excused. The six forbidden phrases it contains ("unseen treatments",
"other cell lines", "patients", "clinical", "calibrated", "Independent biological replication")
each sit in a clause that negates them.

Two refinements were needed, and both came from real false positives:

```text
  WHITESPACE IS NORMALISED FIRST
    a newline is not a clause boundary. Treating it as one split "...and not a /
    clinical recommendation." across a line wrap in the shipped predictor's docstring
    and reported a correctly-negated sentence as a forbidden claim.

  A WORD MAY NEGATE ITSELF
    `uncalibrated` contains `calibrated` and the Stage-26 pattern has no word boundary,
    so "outputs an uncalibrated score" read as claiming a calibrated probability. A
    match immediately preceded by un-/non-/de- is excused. The prefix must be
    contiguous, so "run calibrated" is still caught.

  STRUCTURED TEXT KEEPS THE WINDOW RULE
    clause scoping splits a JSON key from the value that negates it:
    `"calibrated_probability": "NEVER emitted in Generation 1"` reads as a claim. The
    locked surfaces are gated with the EXTENDED patterns under Stage 26's window rule,
    and the clause-scoped reading is reported beside it. Only the negation scope
    differs, matched to the kind of text being read. Both readings are now empty.
```

**B. The claims document had no identity at all.** The evidence lock hashes 54 artifacts and none
of them is this stage's output — correct layering, since re-locking to include them would make CL-A
circular. But it left the document the manuscript is actually written from free to change with
nothing noticing: the exact failure the layer below exists to prevent, one level up. The claim lock
now computes its own digest over its plan, executor, contracts, claims document and verdict, by the
same canonical-LF rule, stored **outside** the verdict it covers so it does not hash itself.
`--verify` re-checks it.

**C. Plan/code drift, for the third time in this sequence.** The plan listed the adversarial corpus
as abbreviated fragments; seven of the fifteen sentences did not appear verbatim. The old test
spot-checked three and passed. The plan now carries all fifteen verbatim with their groups, and the
test asserts every sentence appears in the plan exactly — counts and text, not samples.

## Tests
- 28 claim-lock contracts, 0 skipped (21 before the audit; 7 added for findings A-C)
- the evidence lock re-verifies clean, 54 artifacts, before and after — digest unchanged

---

## Scientific interpretation

**Proves:** the ceiling of what Generation 1 may say is fixed, written down, bound to locked
evidence, and enforced by an instrument demonstrated to catch fifteen tempting forbidden sentences
including three the previous stage's scanner missed.

**Does NOT prove:**
- **anything scientific.** Locking a claim grants nothing. It fixes what may be said about evidence
  locked separately. No analysis ran.
- **that the manuscript will be honest.** A scanner catches phrasings it has patterns for. A claim
  nobody thought to forbid is not covered by a corpus of the ones that were, and an abstract can
  mislead through emphasis, omission or figure choice without tripping a single pattern.
- **that clause scoping is the right granularity.** It is tighter than a ±160-character window and
  it caught three sentences the window let through, but it is still a lexical rule with no
  understanding of scope. A forbidden claim spread across two clauses, or carried by implication
  rather than by phrase, would pass.
- **that the allowed claims are true.** That was Stage 25's job, under its own preregistration.
  This stage checks only that nothing was widened past it.
- **that three claims are the right three.** They are the three the frozen plan permits. Whether
  the paper should make fewer is a judgment the manuscript stage may still exercise downward.

## Next action
`MANUSCRIPT + REPRODUCIBILITY PACKAGE`. The handoff binds it to five things: re-verify the evidence
lock before submission, carry every claim with its mandatory qualifiers, report `p_perm` as
`p < 0.001 (0 of 1,000)` and never as a point estimate, state that independent biological
replication is Generation 2, and scan any new abstract-level sentence against the same nine.

**The ceiling may be lowered. It may not be raised.** No claim-lock outcome reopens an earlier
stage.
