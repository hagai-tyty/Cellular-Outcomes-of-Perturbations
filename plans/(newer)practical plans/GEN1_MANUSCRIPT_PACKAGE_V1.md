# GEN-1 MANUSCRIPT + REPRODUCIBILITY PACKAGE

**Status** V1. The last Generation-1 stage.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Entry** `results/gen1_handoff_to_manuscript.json`, verdict `GEN1_CLAIMS_LOCKED`.

```text
  evidence digest  c84a4d7a2e2d254ed92e43ccf1f91d74da0ab57c8328b1be822d73e5c62ec350
  claim digest     8c820412ba325cd053f0bc9d907d65ef1fd8e1a82aa93e0f947e85957be69796
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

---

# Amendment V1.1 — 2026-09-04

The Entry block above is re-pinned rather than quietly overwritten. V1 was written against
evidence digest `c84a4d7a2e2d254ed92e43ccf1f91d74da0ab57c8328b1be822d73e5c62ec350` and claim
digest `a81ee43b07fae32f9bb45b4a4133de0b1f3979eeda3de7d700a7ca6897affb77`.

Both moved during the release-verification pass: the evidence manifest had been recording each
artifact's `st_size`, which counts a CRLF as two bytes and therefore differed between a file a
stage had just written and the same file in a fresh checkout. Sizes now come from the same
canonical-LF content that is hashed. No artifact content changed and no gate was relaxed.

This plan is covered by the package digest, so it is re-pinned in the same pass that rebuilds
the package.

---

# Amendment V1.2 — 2026-09-12

**What changes is the shape of the manuscript, not what it says.** §2 above pins eleven required
sections, and it stays as written. From V1.2 the manuscript takes the structure of a BMC
Bioinformatics Research article, so that one text serves the preprint and any later journal
submission. On 2026-09-12 the decision was to release through Zenodo and bioRxiv only, defer any
journal, and hold the manuscript to a journal's structure regardless.

Required top-level structure from V1.2, in this order:

```text
ABSTRACT               structured under Background / Results / Conclusions; at most 350 words
KEYWORDS               three to ten
BACKGROUND             was INTRODUCTION, including "Relation to prior work"
METHODS                DATA becomes its first subsection; adds "Use of AI assistance"
RESULTS                THE TOOL becomes its last subsection
DISCUSSION             LIMITATIONS, WHAT THIS DOES NOT SHOW and GENERATION 2 as subsections, in that order
CONCLUSIONS            new
LIST OF ABBREVIATIONS  new
DECLARATIONS           the eight BMC headings; AVAILABILITY becomes "Availability of data and materials"
REFERENCES             moved from inside DATA to the end
FIGURE LEGENDS         new
```

Source for the structure: the BMC Bioinformatics Research-article submission guidelines, read on
2026-09-12 (https://link.springer.com/journal/12859/submission-guidelines/research-article). Source
for placing AI disclosure in Methods: Springer Nature's editorial policy on large language models.

**What does not change.**

```text
every V1 section's content is carried into the new structure -- moved, not rewritten
LIMITATIONS and WHAT THIS DOES NOT SHOW stay separate sections, in the same order
MS-A through MS-F are unchanged; no check is removed or loosened
the anti-rescue firewall (§8) is unchanged: no recorded number, locked artifact or lock
  changes, and the manuscript may state fewer claims than the lock permits, never more
```

**New prose, and how it is held.** The title-page author block, Keywords, a short opening to the
Discussion, Conclusions, Use of AI assistance, List of abbreviations, the Declarations text (carried
from `SUBMISSION.md` §4), and Figure legends (from the figure builder's own titles and
descriptions). All of it is scanned by the same instrument as the rest of the document. The
Conclusions may be assembled only from the locked claims and their qualifiers.

**Checks added by V1.2.** Each must be shown to refuse a broken copy:

```text
the top-level sections appear in the required order
all eight Declarations headings are present
the abstract carries the three structured labels, in order, and is at most 350 words
there are three to ten keywords
```

**Placeholders.** Fields only a human can supply are marked FILL and must block the release bundle.
Fields that cannot exist until after archiving -- the preprint DOI -- are marked LATER and do not.

---

# Amendment V1.3 — 2026-09-13

Two corrections from an external review of the release workflow. Each was checked against the
repository before it was accepted.

**1. When the protocol was frozen, said precisely.** The manuscript said the evaluation was frozen
"before any result existed". The ranking test was frozen by digest before any ranking statistic was
computed -- the Stage-23.5 plan, committed 2026-08-26, states that no Stage-25 statistic exists -- but
after the Stage-23 predictive analyses of the same data, recorded 2026-08-22. The overstated form
appeared 13 times across the release documents. Each now says what is true: clone-held-out
folds, features and exclusions fixed before any model was fitted, and the ranking test frozen before
any ranking statistic was computed, after earlier predictive analyses. An MS-C check refuses the
overstated wording in any release document, with a negative control. The locked claim file
`GEN1_CLAIMS.md` keeps the old wording, because a lock is not edited; the manuscript may say less
than the lock, never more, and here it says less.

**2. A checkout is not the archive (§5).** §5 required the package to name what it does not contain.
The Zenodo archive does contain the 44 MB model artifact, and the pseudobulk cache the null refits on,
so "not in this package" was true of a GitHub checkout and false of the archive. `REPRODUCIBILITY.md`
now separates the two and verifies the unpacked archive with the three verifiers and the predictor,
run with `PYTHONPATH=src`. That flag is required: inside an unpacked archive, a plain
`python -m cellfate.gen1_cli` imported the checkout's editable install and exited 0 -- a pass on the
wrong code. MS-E gains a check for both points, with a negative control.

No gate is loosened. MS-A to MS-F otherwise stand as written, with the structure of V1.2.

---

# Amendment V1.4 — 2026-09-15

**What changes: the Conclusions say what the result means.** V1.2 said the Conclusions "may be
assembled only from the locked claims and their qualifiers". Read that way, they had become a second
list of limitations: the Discussion already carries the limitations in full, and the author asked for
the Conclusions to say what was found, what it means and what others can do with it, and not to refer
to Generation 2. From V1.4 the Conclusions state the locked result with its qualifiers, then its
significance and the directions it opens; the abstract's Conclusions do the same in brief.

Significance and directions are written as implication or possibility -- "what this means", "could",
"if the result holds" -- and never as a finding.

**What does not change.**

```text
the ceiling: nothing above the claim lock, and every new sentence is scanned by the same instrument
the qualifiers travel with the finding -- WM989, the six observed experimental conditions, a detection
  proxy that is not death, clone-held-out evaluation, and no independent replication yet
LIMITATIONS, WHAT THIS DOES NOT SHOW and GENERATION 2 stay in the Discussion, in that order
MS-A to MS-F are unchanged; no check is removed or loosened
```


---

# Amendment V1.5 — 2026-09-18

**What changes: the preprint venue.** V1.2 and the plan's release section named bioRxiv as the
preprint server. On 2026-09-18 bioRxiv declined the submission, BIORXIV/2026/752106, because it
requires authors to hold an organisational affiliation that can adjudicate ethical disputes. The
manuscript was not read on its merits, and nothing in it was questioned. An independent researcher
cannot meet that requirement without joining an organisation.

From V1.5 the preprint is a Zenodo record of its own, separate from the archive record, carrying the
same `MANUSCRIPT_bioRxiv.pdf` under CC BY, with its own DOI and a related-works link to the archive.
The journal route stays BMC Bioinformatics, whose policy states that a posted preprint is not prior
publication and is not counted against the advance a study provides. The preprint DOI is disclosed
at submission, and the APC waiver is requested there too, since it cannot be requested later.

**What does not change.**

```text
the manuscript, the figures and the locks: nothing here touches a claim or a number
the archived PDF keeps its file name, so the Zenodo archive, the GitHub release and this checkout
  continue to name the same file
the release order: archive first, preprint second, journal third
```
