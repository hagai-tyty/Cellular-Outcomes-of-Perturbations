# stage_26_RECORD — the known-treatment-only scope lock

## Goal
Execute the Stage-23.5 §9 mandate: `record KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`, with
`no unseen-treatment claim and no rescue experiment`.

Stage 25 came back positive. The predictable failure mode of a positive result is that the claim
quietly grows: six observed conditions become "treatments", one melanoma line becomes "cancer", a
detection proxy becomes "response". Stage 26 makes that growth mechanically impossible before the
evidence lock, the claim lock and the manuscript are written on top of it.

So Stage 26 was **not** written as prose. A scope limit asserted in a document is worth nothing if
the shipped tool will happily score `Vemurafenib`.

## Inputs
- `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md` @ canonical-LF `8da16fca...bced48`, FROZEN, digest
  re-asserted before anything was read
- `results/stage24/stage24_w5_artifact.npz` + `.json` — the shipped model, hash-checked
- `results/stage24/stage24f_tool_freeze.json` — the 24F deliverable hashes
- `results/stage25/stage25_verdict.json` — `STAGE_25_RANKING_SUPPORTED`
- `results/stage24/tool/example_clone_*` — one real benchmark clone, `L100615`

## Files added
- `plans/(newer)practical plans/STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md`
  (canonical-LF `1d5caa3296be553c628f117ee60d1ce823e8833cbcf3c57907a2bdd403f12cf8`)
- `experiments/run_stage26_scope_lock.py`
- `tests/test_stage26_scope_lock.py`
- `results/stage26/` — 26A-26E JSONs, `GEN1_SCOPE_LIMIT.md`, `stage26_verdict.json`
- `results/stage26_handoff_to_evidence_lock.json`

## Files modified
- `results/stage24/tool/MODEL_CARD.md` — one appended section, under proof. See below.

## What did NOT change
- No model, coefficient, prediction, benchmark row or recorded number
- No Stage-22/23/23.2/23.5/24/25 artifact — every 24F hash re-verified before and after
- The out-of-fold table Stage 25 consumed still hashes to `abdc2999...`
- Nothing was fitted. Stage 26 has no authority to fit and the check that proves it is gating.

---

## Result

```text
  KNOWN_TREATMENT_ONLY_SCOPED_LIMIT        all five substages pass

  26A  vocabulary closure      56 / 56 adversarial strings refused
  26B  claim surface           8 surfaces, 9 gating claims, 0 violations
  26C  propagation             10 / 10 checks, Python API and CLI
  26D  no rescue               fits nothing, every frozen hash holds
  26E  model-card append       byte-identical base, byte-idempotent re-run

  total runtime                5.5 s
```

### 26A — the vocabulary was attacked, not asserted

56 strings across seven groups, all declared in the plan **before** the run so the corpus could not
be trimmed after seeing which entries failed:

```text
  case          9    acid  ACID  cisplatin  CisPlatin  cocl2  ...
  whitespace    6    " Acid"  "Acid "  "Acid\t"  "Cis platin"  ...
  pharmacology 16    Vemurafenib  PLX4720  Cobimetinib  Carboplatin  Oxaliplatin
                     Paclitaxel  Pembrolizumab  Temozolomide  ...
  format        5    Cisplatin_1uM  "Cisplatin (1 uM)"  Acid+Cisplatin  ...
  control       8    DMSO  Vehicle  Control  Untreated  None  NA  null  0
  confusable    5    Cyrillic-А Acid, Cyrillic-і Acid, CoCl-subscript-2, CoCI2, zero-width space
  structural    7    ""  " "  "*"  "all"  "Acid,Cisplatin"  "[Acid]"  ...

  every one     ->   UNSUPPORTED_TREATMENT, future_detection_score = null
```

The pharmacology block is the one that carries the argument. `Vemurafenib` is *the* drug for a
BRAF-V600E melanoma line and `Carboplatin` is a platinum agent one substitution from a condition
that IS supported. A user who expects the tool to be helpful will type exactly these.

### The reference-leak hazard is real, and the guard is what stops it

This is the sharpest thing Stage 26 found, and it is a property of the frozen design rather than a
bug:

```text
  Acid is the REFERENCE level -> five zero dummies

  p._dummies(["Vemurafenib"])            -> [0, 0, 0, 0, 0]      sum 0.0
  component.score(x, b, those dummies)   -> 0.673109077807343
  component.score(x, b, Acid's dummies)  -> 0.673109077807343    IDENTICAL
```

The dummy encoder, on its own, maps **any** unknown string to the Acid row and would return the
Acid score under another name. Nothing in the linear algebra objects. The vocabulary filter inside
`predict()` is the only thing standing between a user and a silently wrong number, so Stage 26
demonstrates the hazard and then proves the guard blocks it, rather than asserting the guard exists.

### Structural closure

```text
  design columns  309 = 50 PCs + 4 nuisance + 5 dummies + 250 interaction
```

A seventh condition cannot be added without changing 309, which cannot happen without refitting,
which Stage 26 may not do. The vocabulary is closed by geometry, not only by a list.

### 26B — the claim surface, with a working instrument

Eight shipped surfaces scanned for all nine §3.5 forbidden claims. Zero unnegated hits.

That number means nothing on its own, so the scan carries a **canary**: each of the nine claims is
stated plainly and must be caught, and its negated twin must not be.

```text
  claims probed                    9
  caught                           9
  missed                           0
  negated-twin false alarms        0
```

Broad topic words (`cancer`, `death`, `resistance`, `sensitivity`, `response`) are **reported, not
gating** — 11 hits, 0 unnegated. A gate tuned on broad words is a gate tuned until it passes.

### 26C — the limit reaches the caller on the failing paths too

```text
  SUPPORTED_KNOWN_CONDITION     known_limitations present
  UNSUPPORTED_TREATMENT         known_limitations present
  MISSING_REQUIRED_NUISANCE     known_limitations present    (3 ways of triggering it)
  UNSUPPORTED_FEATURE_SCHEMA    known_limitations present

  CLI  all six known             exit 0
  CLI  one unknown of two        exit 2, the unknown row prints null, the known row prints 0.2466
  CLI  no nuisance               exit 2
  CLI  unreadable input          exit 3, support_status on stderr
```

A caller who gets a refusal and no limitations has been told *less* than a caller who got a score,
which is why the refusals are checked rather than only the successes.

The Stage-25 verdict was verified to unlock a **claim, not a computation**: with the verdict file
`ranking_status = SUPPORTED` and the order is exposed; without it `NOT_SUPPORTED` and the order is
withheld; **the six scores are byte-identical either way.**

```text
  validated order, clone L100615, lowest predicted detection first
    Trametinib > Dabrafenib > CoCl2 > Cisplatin > Doxorubicin > Acid
```

### 26E — the model card, appended under proof

```text
  base   8902d4999ad98b9d...   byte-identical to the card 24F froze
  after  489f826266b8edd7...   2,940 -> 5,158 bytes
  append-only proof            the frozen bytes are an exact PREFIX of the new file
  re-run                       byte-idempotent: after == before
  delimiter count in file      1
```

The card's **Ranking** section was written before Stage 25 ran and speaks in the future tense. It
is left byte-for-byte as frozen; the appended section states the actual verdict. Leaving a stale
future tense in a shipped deliverable and correcting it in an appended section is the same rule
this project applies to records.

---

## Bugs found — three, all mine, all in the instrument

**1. The mixed-request check asserted bit-identity across batch sizes.** First run of 26A failed.
The cause is real and had to be named rather than shrugged at (§7.1 R4): the design matrix passed
to the GEMM has 1, 3 or 6 rows depending on how many known conditions were requested, and BLAS
selects a different kernel by shape. Measured difference `1.11e-16` — four orders of magnitude
inside the frozen `1e-12` cell bound. The check was **replaced by a bounded one, never removed**,
and split into three: routing correctness (gating), agreement inside `1e-12` (gating), and §7.1 R3
within-clone ordering identical across batch sizes (gating, because the ranking claim is a function
of orderings and nothing else).

**2. 26D reported its own search terms as evidence of fitting.** `GroupKFold` and `cross_val`
appear in this module exactly once each — inside `FIT_TOKENS`, the list of things to look for. The
exclusion used `src.split("]\n")`, which cut at the first such bracket anywhere in the file rather
than at the end of the token list. Fixed by excising the declaration by line range. Same class of
mistake as a scanner tripping on its own description, which is why the claim scanner was written to
exclude itself from the start.

**3. The model-card append stacked a newline per run while claiming idempotence.** The strip
removed the delimiter but not the separator that preceded it, so `base` drifted by one `\n` every
run and the append-only proof compared against a moving target. "No second section" was true while
the file quietly changed. Fixed by anchoring `base` to the **24F recorded hash** instead of to
whatever the file happened to contain — `base_sha256 == frozen_24F_sha256` is now gating, and a
re-run is byte-identical.

None of the three was a scope hole. All three were checks that would have passed for the wrong
reason, which is exactly what an instrument stage exists to catch in itself.

## Tests
- 23 Stage-26 contracts, 0 skipped
- `tests/test_results_paths.py` green (335 checks) — the convention suite that has broken CI before

---

## Scientific interpretation

**Proves:** the Generation-1 scope limit is enforced in the code that ships, not only asserted in
documentation. The condition vocabulary is closed against 56 adversarial strings including sixteen
real oncology drugs; the reference-leak hazard that would silently return the `Acid` score for an
unknown condition is real and is blocked; none of the nine forbidden claims appears unnegated on any
shipped surface, verified with an instrument shown to fire on all nine; and every response path,
including every refusal, carries the limitations.

**Does NOT prove:**
- **anything scientific.** Stage 26 grants no claim. It records that the existing claim is enforced.
  It ran no analysis, fitted no model, and moved no number.
- **that the tool is safe outside WM989.** It is not applicable there at all: `B` counts a clone's
  cells in WM989's three naive libraries, so data from another lab, cell line or library design
  cannot produce a valid `B`. The scope lock records that boundary; it does not extend it.
- **that the documentation is complete.** Nine claims are gated. A claim nobody thought to forbid is
  not covered by a scan for the nine that were.
- **that `Acid`-as-reference is a good design.** It is the frozen design, and Stage 26's finding is
  that its safety rests entirely on one filter in `predict()`. That is recorded as a standing
  fragility for Generation 2, not repaired here — repairing it would mean refitting.

## Standing asymmetry, recorded rather than fixed

`known_limitations` in the frozen artifact metadata still reads *"ranking is NOT validated until
Stage 25 records RANKING_SUPPORTED"*. That sentence is conditional and still literally true, and
Stage 25 has since satisfied it. The metadata is left byte-frozen because changing it moves a hash
the reproduction path checks. A user who loads the tool without the verdict file therefore sees
`NOT_SUPPORTED` — an **underclaim**, the safe direction. Recorded, not repaired.

## Next action
`GEN-1 EVIDENCE LOCK`. Freeze benchmark, tool, out-of-fold predictions, ranking verdict and
limitations under hashes that refuse to proceed if one has moved; carry the nine forbidden claims
into the claim lock unchanged; record that independent biological replication is Generation 2 and
not a Generation-1 gate. Then claim lock, manuscript, reproducibility package.

No Stage-26 outcome reopens an earlier stage.
