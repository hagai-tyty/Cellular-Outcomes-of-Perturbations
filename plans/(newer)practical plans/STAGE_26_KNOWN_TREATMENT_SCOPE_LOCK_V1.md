# STAGE 26 — KNOWN-TREATMENT-ONLY SCOPE LOCK

**Status** V1, the executable form of one line in the frozen Stage-23.5 ship plan.
**Parent** `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, FROZEN.
**Mandate** §9 of that plan, verbatim:

```text
STAGE 26
  record KNOWN_TREATMENT_ONLY_SCOPED_LIMIT
  no unseen-treatment claim and no rescue experiment
```

---

# 0. What this stage is, and what it is not

Stage 25 produced a positive result. The predictable failure mode of a positive result is that the
claim quietly grows: six observed conditions become "treatments", one melanoma line becomes
"cancer", and a detection proxy becomes "response". Stage 26 exists to make that growth mechanically
impossible before the evidence lock, the claim lock and the manuscript are written on top of it.

So Stage 26 is **not** a prose exercise. Asserting a scope limit in a document is worth nothing if
the shipped tool will happily score `Vemurafenib`. This stage verifies the limit **in the code that
ships**, adversarially, and only then records it.

## 0.1 Authority

```text
MAY      read every frozen artifact; run the shipped tool on adversarial inputs
MAY      write results/stage26/*, and append a delimited Stage-26 section to the model card
MAY      fix the TOOL if a scope hole is found, then re-run Stage 26 from the top

MAY NOT  fit anything, on any data, for any reason
MAY NOT  touch the model artifact, the out-of-fold table, or the Stage-25 verdict
MAY NOT  add a condition, widen a vocabulary, or relax a refusal
MAY NOT  reopen Stage 25 or any earlier stage under any outcome
```

## 0.2 Compute budget

Seconds. Stage 26 runs the shipped tool a few hundred times on one clone and scans text. If this
stage ever needs a long run, something is wrong with its design.

---

# 1. The nine forbidden claims

Taken verbatim from §3.5 of the frozen plan. These are the scan targets in §3 below.

```text
 1  unseen-treatment generalization
 2  cross-cell-line or cross-patient generalization
 3  clinical treatment recommendation
 4  causal treatment-effect estimation
 5  calibrated probability unless calibration is separately frozen and passed
 6  independent biological replication of Role B
 7  uniform benefit across all six conditions
 8  confirmed Role-A prediction
 9  single-cell input equivalence when the model was trained on clone pseudobulk
```

---

# 2. 26A — vocabulary closure, adversarially

The frozen vocabulary is exactly:

```text
Acid  Cisplatin  CoCl2  Dabrafenib  Doxorubicin  Trametinib      reference = Acid
```

## 2.1 The corpus is fixed here, before it is run

Every string below must return `support_status = UNSUPPORTED_TREATMENT` and
`future_detection_score = null`. The corpus is declared in this plan so it cannot be trimmed after
seeing which entries fail.

```text
CASE            acid  ACID  cisplatin  CisPlatin  cocl2  COCL2  dabrafenib  trametinib
WHITESPACE      " Acid"  "Acid "  "Acid\t"  "Cis platin"  "Co Cl2"
PHARMACOLOGY    Vemurafenib  PLX4720  Cobimetinib  Selumetinib  Encorafenib  Binimetinib
                Carboplatin  Oxaliplatin  Paclitaxel  Docetaxel  Etoposide  Nilotinib
                Pembrolizumab  Nivolumab  Temozolomide  5-FU
FORMAT          Cisplatin_1uM  "Cisplatin (1 uM)"  Cisplatin-high  Acid+Cisplatin
CONTROL         DMSO  Vehicle  Control  Untreated  None  NA  null
CONFUSABLE      Cyrillic-A Acid   Cyrillic-i Acid   CoCl-subscript-2   CoCI2 (capital i)
STRUCTURAL      ""  " "  "*"  "all"  "Acid,Cisplatin"  "[Acid]"
```

The pharmacology block is the one that matters. `Vemurafenib` is *the* drug for a BRAF-V600E
melanoma line, and `Carboplatin` is a platinum agent one substitution away from a condition that IS
in the vocabulary. A user who expects the tool to be helpful will try exactly these. The tool must
refuse both without hesitation and without a nearest-neighbour suggestion.

## 2.2 The reference-leak test

`Acid` is the reference level: it is encoded as five zero dummies. An implementation that builds
dummies by equality over unrestricted input would therefore produce **all zeros — the Acid row — and
return the Acid score under a different name**. This is the single most dangerous silent failure
available to this tool.

```text
PASS requires   no unknown condition ever reaches the dummy encoder
PASS requires   for a fixed clone, no refused condition returns any number at all
PASS requires   the number of scored rows equals the number of KNOWN requested conditions
PASS requires   a refused condition's row is not numerically equal to the reference row
```

## 2.3 Structural closure

```text
len(treatment_vocabulary) == 6
design columns == 309 == 50 PCs + 4 nuisance + 5 dummies + 250 interaction
```

A seventh condition cannot be added without changing 309, which cannot happen without refitting,
which Stage 26 may not do. The vocabulary is closed by the geometry of the frozen model, not only
by a list.

---

# 3. 26B — claim-surface scan

Every user-facing surface is scanned for the nine forbidden claims:

```text
src/cellfate/gen1_predictor.py            src/cellfate/gen1_cli.py
results/stage24/tool/MODEL_CARD.md        results/stage24/tool/io_schema.json
results/stage24/tool/example_clone_README.md
results/stage24/stage24_w5_artifact.json     known_limitations, support_flags
results/stage25/stage25_verdict.json         consequence, standing_limitations
results/stage26/GEN1_SCOPE_LIMIT.md
```

## 3.1 The rule

A forbidden-claim phrase may appear **only inside a negation**. For each hit, a negating token must
occur within 160 characters:

```text
not  never  no  cannot  without  forbidden  refus  unsupported  withheld  limit  outside  only
```

A hit with no nearby negation is a `CLAIM_SURFACE_VIOLATION` and fails the stage.

## 3.2 The scanner does not scan itself

The Stage-26 module, its tests and this plan are excluded from their own scan. A scanner that trips
on its own description of what it is looking for produces a false failure and teaches nothing — that
mistake was already made once in Stage 23.2 and will not be repeated here.

---

# 4. 26C — scope-limit propagation

The limit must reach the caller on **every** path, including the failing ones. A caller who gets a
refusal and no limitations has been told less than a caller who got a score.

```text
predict()          SUPPORTED_KNOWN_CONDITION      known_limitations present
predict()          UNSUPPORTED_TREATMENT          known_limitations present
predict()          MISSING_REQUIRED_NUISANCE      known_limitations present
predict()          UNSUPPORTED_FEATURE_SCHEMA     known_limitations present
rank_conditions()  with and without verdict       known_limitations present
CLI                exit 0 / 2 / 3                 support_status present on stdout or stderr
```

Additional invariants:

```text
the known-treatment-only limitation is IN the shipped known_limitations list
no output key named calibrated_probability exists on any path
CLI exit code is 2 whenever any requested condition was refused, never 0
```

---

# 5. 26D — no rescue experiment

```text
the Stage-26 module contains no model fitting
every Stage-24 deliverable hash is re-verified, before and after Stage 26 runs
the Stage-25 verdict file is byte-identical to the one Stage 25 wrote
no file under results/stage22, stage23*, stage24 or stage25 is modified
```

The model card is the one deliberate exception, handled in §6.2, under a byte-level proof that only
an appended section changed.

---

# 6. 26E — the record

## 6.1 `results/stage26/GEN1_SCOPE_LIMIT.md`

The authoritative scope document that the evidence lock, claim lock and manuscript consume. It
states what Generation 1 may claim, what it may not, and the exact vocabulary the tool serves.

## 6.2 Model-card Stage-26 section

The shipped model card was written before Stage 25 ran and speaks about the ranking verdict in the
future tense. A tool user reading it today cannot tell that the verdict exists. Stage 26 appends one
delimited section stating the actual verdict and the scope limit.

This moves a hash that 24F recorded, so it is done under proof:

```text
the pre-Stage-26 bytes must be an exact PREFIX of the post-Stage-26 bytes
the old hash and the new hash are both recorded
nothing above the delimiter is altered by a single byte
appending twice is idempotent -- a re-run must not stack a second section
```

## 6.3 Verdict

```text
KNOWN_TREATMENT_ONLY_SCOPED_LIMIT      all of 26A-26D pass
STAGE_26_SCOPE_HOLE_FOUND              any check fails
```

A `STAGE_26_SCOPE_HOLE_FOUND` is **not** a scientific result and does not touch any claim. It means
the tool leaks scope, the tool gets fixed, and Stage 26 re-runs from the top. There is no path in
which a scope hole is resolved by widening the scope.

---

# 7. Anti-rescue firewall

```text
no Stage-26 outcome reopens Stage 25, Stage 24, or any earlier stage
no Stage-26 outcome changes delta_RANK, the ranking verdict, or any recorded number
no Stage-26 outcome authorizes new data, a new condition, or a new model
a PASS grants no claim; it only records that the existing claim is enforced in code
```

---

# 8. Handoff

```text
KNOWN_TREATMENT_ONLY_SCOPED_LIMIT  ->  GEN-1 EVIDENCE LOCK
```

Evidence lock freezes benchmark, tool, out-of-fold predictions, ranking verdict and limitations.
Claim lock and manuscript follow. Independent biological replication remains Generation 2 and is
not a Generation-1 gate.
