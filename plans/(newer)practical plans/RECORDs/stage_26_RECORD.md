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
  (canonical-LF `f34a85055071e59be0d169cace92ba38892b49b5d9f22e4b19248b9b3535e58d` after the
  self-audit below; it was `1d5caa32...f12cf8` on the first pass)
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
  KNOWN_TREATMENT_ONLY_SCOPED_LIMIT        all eight gates pass

  26A  vocabulary closure      11 checks   56 / 56 adversarial strings refused
  26B  claim surface            4 checks   8 surfaces, 9 gating claims, 0 violations
  26C  propagation             11 checks   Python API and CLI
  26D  no rescue                7 checks   fits nothing, every frozen hash holds
  26E  model card                          append-only against the 24F hash
  26E  hashes after the run                all six re-verified at the end, not only the start
  26E  module stamps                       all four substages from one executor
  26E  evidence-lock inputs                all ten paths exist

  total runtime                5.4 s
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
  CLI  every printed row         known_limitations present, refusals included
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

## A note on editing the plan after running against it

Stage 23.2 established that a frozen protocol must not be edited while a run against it executes.
That rule protects a **preregistered inference test**, where a mid-run edit can move a threshold
toward a result the analyst has already glimpsed.

Stage 26 is not that. It computes no statistic, has no null, and has no degree of freedom that could
be tuned toward a preferred answer — every check is a deterministic pass/fail on frozen bytes. The
plan was therefore edited during the self-audit, under three rules held throughout:

```text
  every change makes the gate STRICTER or fixes a plan/code disagreement, never looser
  the effect of each change on the verdict is measured and reported (finding B)
  the whole stage is re-run from the top and the digest re-recorded
```

The plan digest moved from `1d5caa32...` to `f34a8505...` and the verdict was re-derived, not
carried over. Stage 25's frozen plan, digest `8da16fca...bced48`, was not touched and is asserted
on every run.

---

## Self-audit before the evidence lock — five more findings

The stage was re-read end to end before handing anything to the evidence lock. Five inconsistencies
between the plan and the module, none of which changed the verdict, all of which are now gated.

**A. The plan declared 50 adversarial strings; the module ran 56.** A superset, so nothing was
cherry-picked away — but the plan's whole purpose in §2.1 is that the corpus is fixed in writing,
and a module that disagrees with it defeats that. The existing test spot-checked four sentinel
strings and sailed past a six-string gap. Fixed: §2.1 now lists all 56 with a per-group count and
spells out the five unicode confusables by codepoint; the module carries `EXPECTED_GROUP_SIZES` and
fails if a group shrinks; the test compares counts parsed from the plan against the module rather
than checking sentinels.

**B. The negation list in code was 26 tokens; the plan declared 12.** This is the one that mattered,
because a longer negation list is a **looser** gate than the one written down. Measured before
changing anything:

```text
  violations under the 26-token code list     0
  violations under the 12-token plan list     0
  hits rescued by exactly one token           25, every one of them by `not` or `never`
```

None of the fourteen undeclared tokens ever rescued a hit. The verdict was unaffected — but the code
was tightened to the plan's twelve anyway, because a gate that quietly loosens is not a gate. The
plan now also spells out that `"no "` carries a trailing space, since bare `no` matches inside
*not*, *none*, *know* and *cannot*.

**C. The handoff told the evidence lock to hash `"...artifact.npz + .json"`.** A human-readable
string in a machine-readable field. The evidence lock would have tried to hash a path that does not
exist. `evidence_to_lock` is now a list of real paths per group — ten of them, including the API,
the CLI, the model card and the schema, which the first version omitted entirely — and every one is
checked to exist before the handoff is written.

**D. §5 said hashes are re-verified "before and after"; the code only verified before.** Checking
only at the start proves nothing about what the run then did. 26E now re-verifies all six at the
end and gates on it.

**E. Nothing tied the substage JSONs to the executor that produced them.** 26E merges the 26A/26C/26D
files left on disk. During this very stage the module was repaired between substage runs, so a
verdict assembled from a mixture of versions was not hypothetical — it was the actual situation for
part of the session. Every substage now stamps the SHA-256 of the executor and 26E refuses a verdict
unless all four stamps match the running module.

**Finding E caught a sixth bug on its first execution.** `write_json` stamped the file but returned
the unstamped `obj`, and 26E used the *returned* 26B dict while reading 26A/26C/26D from disk — so
26B had no stamp and the check failed immediately. `write_json` now returns exactly what it wrote.

Two things audited and found clean: `_scan_text` indexes a lowercased copy with offsets from that
same copy, which would break on any character whose `lower()` changes length — zero such characters
exist in any scanned surface; and every other path the handoff names resolves.

## Tests
- 27 Stage-26 contracts, 0 skipped (23 before the audit; 4 added for findings A-E)
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
