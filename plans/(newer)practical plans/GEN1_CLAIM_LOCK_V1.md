# GEN-1 CLAIM LOCK

**Status** V1.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Entry** `results/gen1_handoff_to_claim_lock.json`, verdict `GEN1_EVIDENCE_LOCKED`,
lock digest `455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9`.

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
re-run the evidence-lock verifier over all 54 artifacts
the lock digest must equal 455892ff50de483fe6e82097f0ab7b96476781d6037e56d93106643045a8b1a9
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

```text
GENERALISATION   "the model generalises to new treatments"
                 "our approach works for any drug"
CROSS-SYSTEM     "we predict outcomes in cancer cells"
                 "validated on other cell lines"
                 "predicts patient response"
CLINICAL         "identifies the best treatment for each clone"
                 "supports clinical decision-making"
                 "a therapeutic recommendation tool"
CAUSAL           "we estimate the causal effect of each treatment"
CALIBRATION      "outputs a calibrated probability of death"
REPLICATION      "independently replicated in an external cohort"
UNIFORMITY       "the interaction helps uniformly across all six conditions"
ROLE A           "confirmed in a second system"
                 "Rewind confirms the result"
SINGLE CELL      "scores a single cell"
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

---

# 5. CL-E — the abstract-level statement

One worked example: a complete, permitted, abstract-level paragraph assembled only from locked
claims and their mandatory qualifiers. It is scanned by the same instrument and must come back
clean.

This is a demonstration, not a mandate. The manuscript may write its own abstract; it may not write
one this instrument would refuse.

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
