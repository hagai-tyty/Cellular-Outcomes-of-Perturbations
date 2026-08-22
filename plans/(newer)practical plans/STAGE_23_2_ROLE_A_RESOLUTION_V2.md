# STAGE 23.2 — Role-A Resolution / Failure Decomposition

**Status:** DRAFT, ROADMAP-REVIEWED — commit this plan alone, then complete an independent repository/data-machine pre-execution audit before running any 23.2B–23.2G decomposition. No Stage-24 model may be fitted while Stage 23.2 is unresolved.  
**Plan version:** V2  
**Roadmap dependency:** `CELLFATE_RX_CURRENT_VISION_ROADMAP_V5.md`  
**Historical dependency:** Stage 23 is closed. Its Role-A permutation failure is permanent and may not be rewritten by this stage.  
**Stage-23 closure anchor at plan drafting:** `2e04ccf` (audit must verify the committed closure record/artifacts rather than trusting this string alone).  
**Naming compatibility:** closed Stage-23 artifacts may still contain the historical next-stage label `STAGE 23R`. In this plan and the current roadmap, the canonical name is **Stage 23.2**. Treat `STAGE 23R` only as a historical alias for Stage 23.2; do not rewrite closed Stage-23 artifacts merely to rename it.

---

# V2 change log

V1 is preserved unmodified as historical at commit `c06fc98` and archived under `arcive/`. It was
never executed. Its independent pre-execution audit verified all 30 Stage-23 numeric anchors in
section 2 exactly, confirmed the closure anchor, and **proved the historical permutation mappings
are exactly recoverable** — two replayed draws reproduced the committed null values bit-for-bit
(`|diff| = 0.00e+00`). Those findings carry forward unchanged.

The audit also found one executable contradiction and several factual errors. V2 corrects them and
applies seven design changes:

```text
BLOCKER
  B1  V1 section 8.5 required per-GSM gDNA concordance. The gDNA table carries
      SampleNum = 3 on all 49,554 rows -- one pooled library -- so outcome support
      is not attributable to either control GSM. Removed in V2 (see 8.5).

FACTUAL
  F1  the benchmark DOES record replicate structure: biological_replicate = R1
      for all 3,905 cells, generalization_scope = within_R1_clone_heldout,
      and GEO titles read "biol rep 1, sample 1" / "biol rep 1, sample 2".
  F2  V1's evidence hierarchy ranked GEO title above raw-file naming, which would
      have inverted the frozen SampleNum -> GSM mapping. Tie-break added.
  F3  the per-permutation D00 array was never committed; only summary statistics
      are in git. V2 commits the 200 values as a compact artifact.
  F4  there is no separate historical bootstrap artifact; it lives inside
      stage23_rewind_results.json.
  F5  GSE227151 contains biological replicates 2 and 3 (GSM7092517-21) that
      Stage 23 never used. V1 never mentioned them.

DESIGN CHANGES
  1  three fixed-K arms (10/20/50) replace the single K=20 matched arm; their
     equal-weight paired average is the primary no-K-selection reference.
  2  cross-GSM gDNA concordance removed from 23.2D entirely.
  3  OUTCOME_LABEL_LIMITATION = NOT_SUPPORTED is no longer reachable from
     multinomial / top-N stability alone.
  4  replicate semantics corrected to one biological replicate, R1.
  5  power interpretation split into within-R1 event-count detectability versus
     biological-replication limitation.
  6  the 200 historical D00 values become a committed reproducibility artifact;
     mappings stay cache-only with a committed digest.
  7  biological reps 2/3 enter a metadata-only reserved confirmation ledger.
```

Every other V1 rule is preserved.

---

# 0. Purpose and non-negotiable interpretation

Stage 23.2 exists because Stage 23 produced a specific, scientifically important contradiction:

```text
Rewind / Role A:
  bootstrap learnability candidate = positive
  full-refit permutation gate      = FAIL

WM989 / Role B:
  additive state signal            = PASS
  explicit X×U interaction         = PASS
  structural controls              = PASS

roadmap gate:
  STAGE_24_BLOCKED_ROLE_A
```

The purpose of Stage 23.2 is **not** to keep trying Rewind models until one passes.

It is to determine which combination of the following five explanations is supported:

```text
1. NO ROBUST BIOLOGICAL SIGNAL
2. MODEL-SELECTION NULL INFLATION
3. RESIDUAL DEPTH / SAMPLING STRUCTURE
4. OUTCOME-LABEL / MEASUREMENT LIMITATION
5. POWER / EXPERIMENTAL-UNIT LIMITATION
     5a. within-R1 event-count detectability   (simulable, see 9.5.1)
     5b. biological-replication limitation     (design fact, see 9.5.2)
```

These explanations are **not mutually exclusive**.

Stage 23.2 must end with one of two useful outcomes:

```text
A. the Role-A problem is resolved well enough to define and independently
   confirm a corrected prospective Role-A claim, after which Stage 24 may
   reopen under a frozen handoff contract;

or

B. the exact reason Stage 24 must remain blocked is known, together with the
   specific new data / benchmark redesign / roadmap revision required next.
```

The historical Stage-23 Role-A verdict remains:

```text
ROLE_A_SIGNAL_FAIL
```

forever.

No Stage-23.2 result may replace that historical verdict with a retrospective PASS.

---

# 1. Lifecycle of this plan

This plan follows the same freeze-before-fit discipline used for Stage 23.

Execution order:

```text
1. commit STAGE_23_2_ROLE_A_RESOLUTION_V2.md alone
   (V1 was committed at c06fc98, audited, and is archived unmodified;
    steps 1-3 below have already run once against V1)

2. independent pre-execution audit
   - read this plan end to end
   - verify the repository contains the frozen Stage-23 closure anchored by
     commit `2e04ccf` or an additive descendant that preserves its numbers
   - compare every frozen Stage-23 anchor against committed artifacts
   - verify the legacy `STAGE 23R` text, if present in closed artifacts, is
     only a naming alias for current Stage 23.2 and not a gate contradiction
   - verify all required raw/source files exist
   - verify the proposed diagnostics are mechanically executable
   - do NOT run 23.2B–23.2G
   - do NOT fit diagnostic models
   - do NOT change Stage-23 artifacts
   - do NOT edit the roadmap merely to make the audit pass
   - report whether V1 is executable exactly as written

3. if the audit finds an executable contradiction
   - preserve the committed plan version as historical
   - create the next version
   - independently review and commit it
   - do not execute until the corrected plan is frozen

   ALREADY EXERCISED ONCE: the V1 audit found blocker B1 (section 8.5),
   factual errors F1-F5, and the K=20 confound. V1 is preserved at c06fc98
   and archived; this document is the corrected version. V2 must itself be
   reviewed before step 4 begins.

4. execute 23.2A
   - source-design / protocol / provenance freeze
   - no decomposition model fitting

5. execute 23.2B
   - model-selection null decomposition

6. execute 23.2C
   - residual depth / nuisance decomposition

7. execute 23.2D
   - outcome-label reliability
   - NO predictive model selection on alternate labels

8. execute 23.2E
   - power / identifiability / sample-size planning

9. execute 23.2F
   - diagnostic synthesis
   - freeze confirmatory protocol
   - freeze Stage-24 handoff draft

10. execute 23.2G only if suitable untouched confirmation evidence exists
    - independent confirmation
    - final roadmap resolution

11. full tests / lint / determinism / clean-tree audit
    commit + push
    STOP

12. Stage 24 may start only if the final roadmap status explicitly permits it
```

---

# 2. Frozen historical evidence from Stage 23

23.2A must verify these against the committed Stage-23 artifacts before any diagnostic run.

## 2.1 Rewind benchmark

Expected frozen Role-A benchmark:

```text
retained clones             3,147
positive clones                35
negative clones             3,112
positive prevalence       1.112...%
outer folds                    5
positive clones / fold         7
Gene Expression features  36,601 before fold-local filtering
```

The two control GSMs used by the benchmark are:

```text
GSM7092515
GSM7092516
```

The existing project treated them operationally as lanes.

The existing record also reports:

```text
306 retained clones span both GSMs
```

**Corrected in V2.** V1 said the record does not establish whether these are biological replicates
or technical lanes. It does. Two independent sources agree:

```text
benchmark    stage22_rewind_cells.csv
             biological_replicate  = R1   for all 3,905 cells
             generalization_scope  = within_R1_clone_heldout

GEO          GSM7092515  "biol rep 1, sample 1, hiFT cells, ungated"
             GSM7092516  "biol rep 1, sample 2, hiFT cells, ungated"
```

The frozen Rewind benchmark is therefore **one biological replicate, R1**, observed as two samples.
It is not two biological replicates, and `REPLICATE_STRUCTURE_BIOLOGICAL` is not an available
finding for these two GSMs.

What remains genuinely open is finer, and is the only question 23.2A resolves here:

```text
are "sample 1" and "sample 2" two 10X lanes of one cell suspension,
or two separately handled libraries / cultures within biological replicate R1?
```

23.2A answers only that within-R1 question. Section 5.2 defines the allowed findings.

## 2.2 Frozen Stage-23 Role-A models

Historical Stage-23 nuisance block:

```text
B0 =
  log1p(n_pretreatment_cells)
  n_lanes
```

Historical models:

```text
R0 = outer-training prevalence
R1 = B0
R2 = PCA(X)
R3 = PCA(X) + B0
```

Primary Stage-23 statistic:

```text
ΔAP_state = AP(R3) - AP(R1)
```

Historical logistic grid:

```text
C ∈ {0.01, 0.1, 1, 10}
```

Historical expression dimension grid:

```text
K ∈ {10, 20, 50}
```

Therefore:

```text
R1 inner candidate count = 4
R3 inner candidate count = 12
```

Historical inner CV:

```text
3-fold StratifiedKFold
shuffle=True
random_state=23023
selection metric = mean Average Precision
unit = clone
```

Historical outer folds remain immutable.

## 2.3 Stage-23 Role-A observed result

Expected Stage-23B pooled OOF result:

```text
R0 AP       0.01112
R1 AP       0.01035
R2 AP       0.01923
R3 AP       0.02085

R3 ROC-AUC  0.6628

ΔAP_state   +0.01050
95% bootstrap CI
             [+0.00397, +0.02258]

fold direction R3-R1
             positive in 5 / 5 folds
             diagnostic only
```

The bootstrap result is historical evidence about the **fitted OOF predictions**.

It is not evidence that the model-selection procedure is exceptional under the null.

## 2.4 Stage-23 Role-A permutation failure

Historical Stage-23E expression permutation:

```text
200 permutations
base seed = 23323

whole clone-level X profile permuted
training profiles remain inside outer training
test profiles remain inside outer test

Rewind strata:
  n_pretreatment_cells = 1 / 2 / 3+
  crossed with n_lanes
```

The implementation always crosses the two factors; "where possible" in V1 was imprecise. The
realized partition, verified during the V1 audit, has exactly five non-empty cells:

```text
1|1   2,584 clones
2|1     220
2|2     196
3+|2    110
3+|1     37
```

There is no `1|2` cell: a clone with one pretreatment cell cannot span two lanes. 23.2C reasons
about how much technical structure survives inside these cells, so the realized partition -- not
the nominal 3x2 grid -- is the object of that analysis.

Every permutation reran the nested model-selection pipeline.

Expected frozen null:

```text
observed ΔAP         +0.01050

null mean            +0.00350
null sd              +0.00631
null p95             +0.01455
null max             +0.05144

null >= observed     16 / 200
p_perm               17 / 201 = 0.084577...
```

Historical gate:

```text
ROLE_A_PERMUTATION_PASS = false
ROLE_A_SIGNAL_FAIL      = true
```

This failure is the object to explain.

## 2.5 What Stage 23E did NOT identify

Stage 23E did not separate at least these two mechanisms:

```text
A. R3 has a broader model-selection path than R1
   (12 candidates versus 4)

B. the coarse abundance-preserving expression permutation may retain
   residual sampling/depth information that B0 does not absorb
```

Stage 23.2B and 23.2C are specifically designed to separate them.

The V1 audit confirmed both mechanisms remain live and that Stage 23E cannot distinguish them: the
historical `R3` search spans 12 `(K, C)` candidates against `R1`'s 4, and the frozen strata are
coarse enough (2,584 of 3,147 clones in a single cell) that within-stratum depth similarity can
survive the permutation.

## 2.6 Role B remains frozen positive evidence

Stage 23.2 is about Role A.

It must not refit, reinterpret, or weaken the successful WM989 result.

The Stage-23 final ledger is expected to preserve:

```text
ROLE_B_ADDITIVE_PASS
INTERACTION_PASS_MULTI_TREATMENT
STRUCTURAL_CONTROLS_PASS
```

Stage 23.2 may read those verdicts for the Stage-24 handoff.

It may not use Role B to silently substitute for the failed Role-A gate.

---

# 3. Global Stage-23.2 rules

## 3.1 Stage 23 is read-only

The following are immutable:

```text
Stage-21 artifacts
Stage-22 benchmark artifacts
Stage-23A–23F artifacts
Stage-23 OOF predictions
Stage-23 null statistics
Stage-23 verdicts
Stage-23 outer folds
Stage-23 benchmark labels
```

Stage 23.2 writes only additive Stage-23.2 artifacts.

If a historical artifact is discovered to be factually wrong, stop and classify the issue before changing anything.

Do not silently repair history inside Stage 23.2.

## 3.2 Diagnosis is not confirmation

Existing Rewind data may be used to diagnose failure mechanisms.

Any corrected analysis designed after inspecting Stage-23 results is:

```text
EXPLORATORY ON REWIND
```

even if its p-value is excellent.

It cannot by itself reopen Stage 24.

## 3.3 No Stage-24 architecture work

Forbidden in Stage 23.2:

```text
neural-network architecture search
new CellFateNet fitting
treatment-ranking model development
Stage-24 hyperparameter search
Stage-25 ranking optimization
Stage-26 unseen-treatment evaluation
Stage-28 calibration / OOD tuning
```

Only simple diagnostic models necessary to decompose the Role-A failure are allowed.

## 3.4 Outer folds remain frozen

Every Rewind predictive diagnostic that uses outcomes must preserve the Stage-22 / Stage-23 outer fold assignment.

No new random outer splits.

No leave-one-positive-out redesign.

No pooled random CV.

## 3.5 Clone is the independent unit

All outcome inference remains at clone grain.

Cells are never treated as independent outcome replicates.

## 3.6 Feature firewall

Primary molecular state remains:

```text
pretreatment Gene Expression only
```

Forbidden predictive inputs remain:

```text
clone_id
cell_uid
GSM / sample identity unless explicitly used as a diagnostic nuisance
outer_fold
outcome-rule fields
gDNA outcome counts / ranks
future outcome fields
source paths
```

Technical scalar summaries may enter **diagnostic nuisance models** only when this plan explicitly allows them.

They are not automatically promoted into the future Stage-24 production baseline.

## 3.7 Alternate labels cannot be model-shopped

23.2D may study the reliability of alternate thresholds, soft labels, or continuous gDNA measurements.

But:

```text
NO R2/R3 predictive model may be fit across a menu of alternate outcome labels
for the purpose of selecting the label that gives the best prediction.
```

That would be direct outcome shopping.

## 3.8 Existing 200 permutation mappings are the primary paired diagnostic basis

23.2B and 23.2C must use the **same 200 Stage-23 permutation mappings** whenever they can be reconstructed or recovered.

Reason:

```text
paired null comparisons remove Monte-Carlo differences between mechanisms
```

The audit must determine whether the exact mapping was stored.

If it was not stored, reconstruct it deterministically from:

```text
base seed
outer-fold partition
frozen strata
deterministic merge rules
```

Then verify reconstruction by rerunning the historical Stage-23 null procedure on a small audited subset and, before accepting the mapping for the decomposition, reproducing the committed historical null statistics within the frozen numerical tolerance.

If the exact historical mappings cannot be recovered, stop and create a revised plan.

Do not substitute a new random 200-permutation sample and call it paired decomposition.

---

# 4. Provenance design — fix the Stage-23 builder-hash problem

Stage 23 repeatedly encountered a bookkeeping failure because the protocol artifact hashed an implementation builder that continued to grow.

Stage 23.2 must not repeat that design.

## 4.1 Immutable protocol identity

23.2A creates:

```text
results/stage23_2/stage23_2_protocol.json
```

It contains the frozen scientific protocol surface only.

Canonical protocol digest:

```text
SHA-256 over canonical JSON:
  UTF-8
  sorted object keys
  compact separators
  LF line endings
  no timestamps
  no absolute paths
  no git commit field inside the hashed payload
```

All Stage-23.2 result artifacts pin:

```text
stage23_2_protocol_sha256
```

They do **not** use the hash of the whole builder/source file as their scientific protocol identity.

## 4.2 Source-code provenance is separate

Every substage records, separately:

```text
git commit
source file hashes
dependency versions
machine/runtime information
```

These may change as later substage code is added.

A later source change must not invalidate an already-frozen earlier scientific protocol if the protocol JSON is unchanged.

## 4.3 Append-only determinism registry

The determinism manifest may grow as new Stage-23.2 artifacts are added.

Rules:

```text
existing artifact digest entries may never silently disappear
new artifact entries may be appended
prior scientific artifacts are not rewritten merely because the registry grew
```

Large matrices / 200-draw detailed caches stay outside git under a gitignored Stage-23.2 cache directory.

Committed artifacts store content hashes sufficient to verify regeneration.

**V2 addition — the historical null array is committed, the mappings are not.**

The V1 audit found that the 200 per-permutation historical `D00` values existed only in the
gitignored Stage-23 cache, while git held summary statistics alone. The whole 23.2B/C paired design
rests on that array, so it must not depend on a machine-local cache:

```text
COMMITTED   results/stage23_2/stage23_2_historical_null_d00.json
              the 200 historical D00 values in permutation_id order,
              full float64 repr, plus their SHA-256
              (~200 floats -- compact by any measure)

CACHE-ONLY  _cc_cache/stage23_2/permutation_mappings/
              the recovered recipient -> donor tables
              (5 folds x 3,147 rows x 200 draws)
              committed as a mapping-set SHA-256 only
```

The array is *replayed*, not copied on trust: 23.2A regenerates all 200 draws from the frozen seed
and pipeline and must reproduce the committed Stage-23 summary statistics exactly before the array
is written. The V1 audit already demonstrated this is bit-exact on two draws.

---

# 5. Stage 23.2A — Resolution protocol + source-design freeze

23.2A fits **no diagnostic predictor**.

Its job is to prove that the failure-decomposition experiment is executable and to settle the experimental-unit semantics as far as the source materials permit.

## 5.1 Mandatory historical-artifact preflight

Verify:

```text
Stage-23 final synthesis exists
Stage-23 closure record exists and states Stage 23 is formally closed
Stage-23 Role-A final verdict = FAIL
Stage-23 Role-B final verdicts unchanged
STRUCTURAL_CONTROLS_PASS = true
ROADMAP GATE = STAGE_24_BLOCKED_ROLE_A

3,147 Rewind clones
35 positives
3,112 negatives
7 positives in each outer fold

historical R1/R3 OOF files hash correctly
historical permutation artifact hashes correctly
historical 200-draw null count = 200
```

**Corrected in V2 (F4).** V1 listed a "historical bootstrap artifact" as if it were a separate
file. There is none: the Role-A bootstrap block lives inside `stage23_rewind_results.json`, whose
hash is already checked. Verify the block's presence and its recorded fields
(`point`, `ci95_low`, `ci95_high`, `replicates = 2000`, `seed = 23123`) rather than a filename.

Also verify, before any decomposition:

```text
results/stage23_permutation_results.json holds SUMMARY statistics only
the 200 per-permutation D00 values are NOT yet committed anywhere
```

Both are expected at plan time. 23.2A is the substage that fixes it, per 4.3.

Also verify:

```text
no Stage-24 model fitted
no Stage-23.2 result artifact already exists except this plan
```

## 5.2 Rewind source-design audit

**Scope corrected in V2.** The biological-replicate question is already settled (see 2.1): the
benchmark is one biological replicate, R1. 23.2A resolves only the finer within-R1 structure of the
two samples.

Resolve it using this evidence hierarchy:

```text
1. source-study paper Methods / supplement
2. GEO sample title + characteristics + relation metadata
3. author-provided metadata
4. author code / comments
5. raw-file naming and library construction metadata
```

**Tie-break added in V2 (F2).** These sources conflict on sample numbering, and V1's ordering would
have silently contradicted the frozen benchmark:

```text
GEO title      GSM7092515 = "sample 1"     GSM7092516 = "sample 2"
file naming    GSM7092515_1_2_control_*    GSM7092516_1_1_control_*
benchmark      SampleNum 2 -> GSM7092515   SampleNum 1 -> GSM7092516
```

The frozen benchmark followed the file naming, which is also what the author's own barcode tables
key on through `SampleNum`. Therefore:

```text
where GEO title numbering and author SampleNum / file naming disagree,
the author SampleNum convention wins for any joins or identifiers,
because the benchmark and the source barcode tables are keyed on it.

the conflict must be recorded verbatim in stage23_2_source_design.json.
the GSM <-> SampleNum mapping is NOT re-derived or changed:
section 3.1 makes the Stage-22 benchmark immutable.
```

This tie-break governs identity and joins only. It has no effect on any Stage-23 number: the two
GSM labels are interchangeable names for lane-composition diagnostics, and `n_lanes` is unaffected.

Clone overlap, barcode reuse, or count patterns may support interpretation but may **not** alone
prove independent handling.

Output exactly one:

```text
WITHIN_R1_TECHNICAL_LANES
    one suspension, two 10X lanes

WITHIN_R1_SEPARATE_LIBRARIES
    separately handled libraries / cultures inside biological replicate R1

WITHIN_R1_STRUCTURE_UNRESOLVED
    source evidence insufficient
```

`REPLICATE_STRUCTURE_BIOLOGICAL` is **removed as a possible finding** for these two GSMs. Nothing in
23.2 may report the current benchmark as containing more than one biological replicate.

Record:

```text
verbatim source evidence snippets / locations
interpretation
what is established
what remains uncertain
```

## 5.2.1 Reserved candidate confirmation ledger — metadata only

**New in V2 (F5, design change 7).** `GSE227151_family.xml` lists thirteen samples. Stage 23 used
two. The series therefore already contains material that Stage 23 has never touched:

```text
GSM7092515  biol rep 1, sample 1, hiFT ungated        USED BY STAGE 23
GSM7092516  biol rep 1, sample 2, hiFT ungated        USED BY STAGE 23
GSM7092517  biol rep 2, sample 4, hiFT ungated        RESERVED CANDIDATE
GSM7092518  biol rep 2, sample 5, hiFT ungated        RESERVED CANDIDATE
GSM7092519  biol rep 3, sample 1, hiFT ungated        RESERVED CANDIDATE
GSM7092520  biol rep 3, sample 2, sorted fast cycling RESERVED, DIFFERENT DESIGN
GSM7092521  biol rep 3, sample 2, sorted slow cycling RESERVED, DIFFERENT DESIGN
GSM7092522-27  iPS cells, DMSO / LSD1i / DOT1Li       OUTCOME-SIDE, NOT PRE-STATE
```

23.2A writes a **metadata-only** ledger:

```text
results/stage23_2/stage23_2_reserved_confirmation_candidates.json
```

recording for each candidate: accession, title, platform, library strategy, declared biological
replicate, declared gating/sorting condition, and whether a matching future-outcome measurement is
*declared* to exist.

Hard restrictions, pre-registered now:

```text
DO NOT download raw matrices for reserved candidates during 23.2A-23.2E
DO NOT compute any expression, barcode, outcome or performance quantity on them
DO NOT inspect their outcome values
DO NOT use them to choose any correction, model, label, nuisance block or threshold
```

Only declared GEO metadata may be read. Reading a title is not inspecting evidence; reading a count
matrix is. The ledger exists so 23.2F can freeze a confirmation protocol knowing what candidate
evidence plausibly exists -- not so 23.2B-23.2E can peek at it.

Whether reps 2/3 carry a reconstructable Role-A outcome is **unverified** and must stay unverified
until 23.2F freezes `STAGE_23_2_ROLE_A_CONFIRMATION_V1.md`. Sorted fast/slow-cycling samples are
flagged separately because gating changes the population and they may not satisfy the same claim.

## 5.3 Verify outcome-rule reconstruction inputs

Before 23.2D, verify access to the exact source inputs needed to reconstruct the gDNA label:

```text
gDNA barcode count source
SampleNum / lane mapping
filtered10XCells source
exact barcode key used by author code
exact top-N selection semantics
tie handling
special-barcode exclusion
retained clone mapping
```

Historical anchor expected from the source reconstruction, all three confirmed by the V1 audit:

```text
top-N = 100
source tie behavior yields 101 selected barcodes
  rank 100 and rank 101 are tied at counts = 2365
gDNA table: 49,554 rows, 1,936 distinct barcodes, total counts N = 782,826
```

**Clarified in V2 (M4).** Two naming facts that must not later be "discovered" as bugs:

```text
the gDNA support column is `counts`, not `nUMI`
  (`nUMI` is a 10X-side column; the gDNA table has cellID / counts /
   BC50StarcodeD8 / SampleNum)

the frozen rule groups by (BC50StarcodeD8, SampleNum), but the gDNA table
carries SampleNum = 3 on every row, so grouping by SampleNum is a no-op
  -- it is equivalent to grouping by barcode alone, and reproduces 101 barcodes
```

Record both verbatim in the source-design artifact.

Do not yet run alternate thresholds.

## 5.4 Verify technical-depth quantities for 23.2C

For every retained Rewind clone, recompute from raw pretreatment Gene Expression counts:

```text
n_pretreatment_cells
n_lanes
total_raw_GE_UMI
n_detected_GE_features_in_raw_pseudobulk
```

These are outcome-free scalar summaries.

Define the frozen expanded diagnostic nuisance block:

```text
Bdepth =
  log1p(n_pretreatment_cells)
  n_lanes
  log1p(total_raw_GE_UMI)
  log1p(n_detected_GE_features_in_raw_pseudobulk)
```

`Bdepth` is diagnostic.

It is not automatically a proposed Stage-24 baseline.

If 23.2A proves the two GSMs are purely technical lanes, 23.2C may additionally run the pre-declared **secondary** lane-composition sensitivity:

```text
log1p(n_cells_GSM7092515)
log1p(n_cells_GSM7092516)
```

If replicate structure is biological or unresolved, those two sample-specific counts are not allowed in the primary depth-nuisance diagnostic because they could absorb biological-unit signal.

## 5.5 Recover the historical permutation mappings

Create a mapping digest for each historical permutation:

```text
permutation_id
outer_fold
recipient_clone_id
donor_expression_clone_id
```

The mapping table itself is cache-only. Commit its digest and shape:

```text
mapping-set SHA-256
number of rows
number of permutations
fixed-clone counts by fold
```

Alongside it, commit the compact historical null array defined in 4.3:

```text
results/stage23_2/stage23_2_historical_null_d00.json
```

**Recoverability is already proven, not assumed.** The V1 audit established that a mapping is a
pure function of `(SEED_PERMUTATION + permutation_id)` replayed through the frozen
`permute_within` five times per draw, once per outer fold, in fold order, and that nothing else
consumes that generator. Two draws were replayed end to end:

```text
draw 0   replayed -0.00004568070915419   committed -0.00004568070915419   |diff| 0.00e+00
draw 1   replayed +0.00156705757555709   committed +0.00156705757555709   |diff| 0.00e+00
```

Draw-0 fixed-clone counts by fold were `{0:13, 1:12, 2:11, 3:10, 4:4}`. Regenerating a draw twice
gives an identical mapping; different draws differ; every mapping was verified train-to-train,
test-to-test, within-stratum and bijective.

23.2A must extend this from two draws to all 200 and require that the replayed array reproduce the
committed Stage-23 summary statistics exactly:

```text
mean +0.00350   sd +0.00631   p95 +0.01455   max +0.05144
count >= observed  16 / 200      p_perm  17/201 = 0.084577...
```

Any mismatch is a stop condition under section 21.

The mapping must obey:

```text
train -> train only
test  -> test only
whole-profile mapping
historical strata / merge rule
```

## 5.6 Freeze primary decomposition design

The core Stage-23.2B/C experiment is a paired 2×2 diagnostic. **Design change 1 in V2** replaces
V1's single `K=20` matched arm with three fixed-K arms whose equal-weight average is the reference.

```text
                     historical B0         expanded Bdepth

FULL R3 SEARCH       cell 00                cell 01
K x C selected
12 candidates

NO-K-SELECTION       cell 10                cell 11
equal-weight mean
of three fixed-K arms
K = 10 / 20 / 50
4 C candidates each
```

Definitions:

```text
cell 00:
  historical R1(B0)
  historical full R3:
    K ∈ {10,20,50}   selected by inner CV
    C ∈ {0.01,0.1,1,10}

cell 01:
  R1depth(Bdepth)
  R3depth full search:
    K ∈ {10,20,50}   selected by inner CV
    C ∈ {0.01,0.1,1,10}

cell 10:
  R1(B0) unchanged
  three separate arms, K fixed at 10, 20, 50 in turn:
    C ∈ {0.01,0.1,1,10}   still selected by inner CV
  D10_j = (1/3) * [ dAP(K=10)_j + dAP(K=20)_j + dAP(K=50)_j ]

cell 11:
  R1depth(Bdepth)
  the same three fixed-K arms under Bdepth
  D11_j = (1/3) * [ dAPdepth(K=10)_j + dAPdepth(K=20)_j + dAPdepth(K=50)_j ]
```

Why three arms rather than one fixed `K`:

```text
V1 fixed K at 20 "by protocol position". The V1 audit found that the historical
R3 selected K = 50, 10, 10, 10, 10 across the five outer folds -- K = 20 was
chosen in ZERO of five. A single K=20 arm would therefore have conflated two
different changes:

    (a) removing K selection            <- the mechanism under study
    (b) moving to a subspace the data never favoured   <- a confound

Averaging the three fixed-K arms with equal weight removes K selection without
privileging or penalising any particular K. No arm is chosen by performance,
and the weights are fixed at 1/3 before execution.
```

Each individual arm is reported as a diagnostic:

```text
per-arm paired means         mean(dAP(K)_j)         for K ∈ {10,20,50}, under B0 and Bdepth
per-arm paired CIs           same bootstrap as the primary
arm dispersion               max arm mean - min arm mean
```

Large arm dispersion is itself informative: it means the K choice matters materially, and it must
be reported next to the primary contrast rather than hidden inside the average. Arm-level results
are **descriptive only** and may not be substituted for the primary equal-weight reference.

Every cell uses:

```text
same outer folds
same inner StratifiedKFold
same C grid
same AP selection metric
same preprocessing
same permutation mapping per permutation_id
```

## 5.7 23.2A verdict

Pass:

```text
STAGE_23_2_PROTOCOL_FROZEN
```

only if:

```text
historical artifacts verified
source-design audit recorded, within-R1 status emitted
reserved confirmation ledger written, metadata only
gDNA reconstruction inputs available
Bdepth exactly computable for all 3,147 clones
historical permutation mappings exactly recoverable for ALL 200 draws
replayed D00 array reproduces the committed Stage-23 summary statistics
committed historical-null artifact written, mapping digest recorded
2×2 design, including all three fixed-K arms, executable in every fold
immutable protocol digest created
```

Otherwise:

```text
STAGE_23_2_INPUT_BLOCKED
```

and stop.

---

# 6. Stage 23.2B — Model-selection null decomposition

Question:

```text
How much of the positive Stage-23 permutation-null center is attributable
to R3's broader K×C model-selection path?
```

This is a diagnostic question.

It does not change the historical null.

## 6.1 Primary paired comparison

Use the historical 200 permutation mappings.

For each permutation `j`, compute:

```text
D00_j =
  ΔAP under cell 00
  = full R3(B0) - R1(B0)

for each K ∈ {10,20,50}:
  A_j(K) =
    ΔAP with K fixed, C selected from {0.01,0.1,1,10}
    = R3(K fixed, B0) - R1(B0)

D10_j =
  ΔAP under cell 10
  = (1/3) * [ A_j(10) + A_j(20) + A_j(50) ]

S_j =
  D00_j - D10_j
```

`S_j` is therefore the paired amount by which allowing the inner CV to *choose* K raises the null
statistic above a fixed-K reference that spans the same grid, with C selection retained in both.

**Corrected in V2 (F3).** V1 said "do not recompute them if the committed null array exists". No
such committed array existed: git held summary statistics only. Per 4.3 and 5.5, 23.2A replays all
200 draws, verifies they reproduce the committed Stage-23 summary statistics exactly, and commits
the array as `stage23_2_historical_null_d00.json`. 23.2B then reads that committed artifact and
does **not** recompute `D00`.

## 6.2 Primary statistic

Model-selection contribution:

```text
selection_shift = mean(S_j)
```

Also report:

```text
median(S_j)
fraction S_j > 0
95% CI for mean(S_j)
```

CI:

```text
10,000 paired bootstrap resamples of permutation_id
seed = 23421
```

The permutation draw is the paired unit.

Report additionally, as diagnostics (design change 1):

```text
per-arm paired mean          mean(A_j(K))        for K = 10, 20, 50
per-arm paired 95% CI        same bootstrap, same permutation IDs
arm dispersion               max_K mean(A_j(K)) - min_K mean(A_j(K))
per-arm S contribution       mean(D00_j - A_j(K))
```

These describe how sensitive the reference is to K. They do not replace `selection_shift`, and no
arm may be promoted to the primary reference after seeing results.

## 6.3 Diagnostic status

```text
MODEL_SELECTION_NULL_INFLATION = SUPPORTED
```

if:

```text
95% CI lower bound of selection_shift > 0
```

```text
MODEL_SELECTION_NULL_INFLATION = NOT_SUPPORTED
```

if:

```text
95% CI upper bound <= 0
```

otherwise:

```text
MODEL_SELECTION_NULL_INFLATION = UNRESOLVED
```

Report the descriptive fraction of the historical positive null mean removed:

```text
fraction_null_mean_explained_by_search =
  selection_shift / mean(D00)
```

only when `mean(D00) > 0`.

Do not gate on an arbitrary fraction threshold.

## 6.4 Search-width ladder — unconditional in V2

V1 ran a width ladder only when the primary status was `UNRESOLVED`. In V2 the three fixed-K arms
of 5.6 are computed for the primary reference anyway, so the ladder costs one extra 8-candidate
search and is run **unconditionally**:

```text
 4-candidate   K fixed          (the three arms of cell 10, reported individually)
 8-candidate   K ∈ {10,20}      selected by inner CV
12-candidate   K ∈ {10,20,50}   selected by inner CV  = historical cell 00
```

Compare the paired null means across the ladder. A monotonic increase in the null center with
search width is descriptive support for model-selection inflation; a flat ladder alongside a
positive `selection_shift` is a signal to re-read the arm dispersion before interpreting.

The ladder remains descriptive. It does not replace the primary paired CI rule, and no ladder rung
may be substituted for the equal-weight reference.

## 6.5 Observed-data sensitivity

For context only, refit the **observed** Role-A comparison under the no-K-selection reference using
the frozen outer folds: the three fixed-K arms and their equal-weight average.

Report:

```text
observed ΔAP per fixed-K arm        K = 10, 20, 50
observed no-K-selection ΔAP         equal-weight average of the three
historical full-search observed ΔAP +0.01050
difference
```

Do not bootstrap a new "rescue p-value" and do not reinterpret the historical Stage-23 failure.

This is diagnostic effect attribution only.

## 6.6 23.2B artifact

Write:

```text
results/stage23_2/stage23_2_model_selection_decomposition.json
```

Large per-permutation details may live in cache.

Committed JSON contains:

```text
protocol hash
historical null artifact hash
mapping-set hash
D00 mean / sd
D10 mean / sd                      equal-weight three-arm reference
per-arm mean / sd / CI             K = 10, 20, 50
arm dispersion
selection_shift
CI
status
search-width ladder means          4 / 8 / 12 candidate
observed no-K-selection sensitivity
```

---

# 7. Stage 23.2C — Residual depth / nuisance decomposition

Question:

```text
How much of the positive Stage-23 permutation null remains because
the permuted expression profile retains technical depth / sparsity structure
that B0 does not represent?
```

## 7.1 Primary expanded nuisance block

Use exactly:

```text
Bdepth =
  log1p(n_pretreatment_cells)
  n_lanes
  log1p(total_raw_GE_UMI)
  log1p(n_detected_GE_features_in_raw_pseudobulk)
```

All continuous nuisance features are standardized training-only.

No gDNA quantity enters Bdepth.

No outcome quantity enters Bdepth.

## 7.2 Full-search depth comparison

For each historical permutation `j`:

```text
D01_j =
  full-search R3depth(Bdepth) - R1depth(Bdepth)

depth_shift_full_j =
  D00_j - D01_j
```

Primary depth contribution:

```text
depth_shift_full = mean(depth_shift_full_j)
```

CI:

```text
10,000 paired bootstrap resamples of permutation_id
seed = 23422
```

`D01_j` keeps full `K x C` selection, so the depth contrast is measured with the selection
mechanism held at its historical setting. The no-K-selection construction enters only through
cell 11 below, exactly as under `B0`.

## 7.3 Joint corrected cell

Also compute, using the **same three fixed-K arms** as 5.6 (design change 1):

```text
for each K ∈ {10,20,50}:
  Adepth_j(K) =
    R3depth(K fixed, Bdepth) - R1depth(Bdepth)

D11_j =
  (1/3) * [ Adepth_j(10) + Adepth_j(20) + Adepth_j(50) ]
```

This is the joint "no K selection + depth expanded" null. The construction is identical under both
nuisance blocks, so the 2×2 stays balanced and the interaction term below remains interpretable.

Report the per-arm `Adepth_j(K)` means, CIs and dispersion exactly as in 6.2.

The full 2×2 null summary is:

```text
μ00 = mean(D00)  full K x C search   + B0
μ10 = mean(D10)  no-K-selection mean + B0
μ01 = mean(D01)  full K x C search   + Bdepth
μ11 = mean(D11)  no-K-selection mean + Bdepth
```

Report:

```text
selection main contrast = μ00 - μ10
depth main contrast     = μ00 - μ01

joint residual null     = μ11

factor interaction      =
  μ00 - μ10 - μ01 + μ11
```

Bootstrap paired CIs for every contrast using the same permutation IDs.

Do not force an additive explanation if the factor interaction is large.

### Observed 2×2 diagnostic

Fit the corresponding four **observed-label** comparisons once under the frozen Stage-23 outer folds:

```text
O00 = historical observed full-search    + B0     ΔAP
O10 = observed no-K-selection mean       + B0     ΔAP
O01 = observed full-search               + Bdepth ΔAP
O11 = observed no-K-selection mean       + Bdepth ΔAP
```

Each `O` cell whose definition involves the no-K-selection reference is the equal-weight average of
its three fixed-K arms, computed on observed labels under the frozen outer folds.

`O00` must equal the historical `+0.01050` within the frozen numerical tolerance.

For the jointly corrected diagnostic cell report:

```text
diagnostic_null_centered_effect = O11 - μ11

q95_11 = 95th percentile(D11)

p_diag_11 =
  (1 + number of D11_j >= O11) / 201
```

Define:

```text
CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC = POSITIVE
```

only when:

```text
O11 > 0
O11 > q95_11
p_diag_11 <= 0.05
```

otherwise:

```text
CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC = NEGATIVE
```

This is an **exploratory same-Rewind diagnostic**.

Even a POSITIVE result cannot emit `ROLE_A_CONFIRMATORY_SUPPORTED`.

## 7.4 Direct technical-structure retention diagnostic

For every historical permutation mapping, compute outcome-free technical alignment between the donor expression profile and recipient clone.

At minimum:

```text
donor pseudobulk nonzero-gene count
  vs recipient log1p(total_raw_GE_UMI)

donor pseudobulk nonzero-gene count
  vs recipient log1p(n_detected_GE_features_in_raw_pseudobulk)

donor total_raw_GE_UMI
  vs recipient total_raw_GE_UMI
```

Use Spearman correlation within each outer partition, then aggregate by permutation.

The donor quantities are properties of the permuted profile source clone.

The recipient quantities are properties of the clone receiving that profile.

This directly measures how much continuous technical similarity survives the coarse Stage-23 strata.

Report:

```text
median correlation
95% bootstrap CI across permutation_id
```

No outcome labels are used in this diagnostic.

## 7.5 Optional lane-composition sensitivity

Run only if 23.2A establishes:

```text
WITHIN_R1_TECHNICAL_LANES
```

Add diagnostic sample-composition nuisance:

```text
log1p(n_cells_GSM7092515)
log1p(n_cells_GSM7092516)
```

This is secondary.

**Updated in V2 (design change 4).** The gate is now the within-R1 finding, since
`REPLICATE_STRUCTURE_BIOLOGICAL` no longer exists as an outcome. It remains forbidden under
`WITHIN_R1_SEPARATE_LIBRARIES` or `WITHIN_R1_STRUCTURE_UNRESOLVED`, because separately handled
libraries could carry biological-unit structure that these two counts would absorb.

## 7.6 Residual-depth status

Primary outcome-level evidence:

```text
depth_shift_full = μ00 - μ01
```

Mechanistic evidence:

```text
technical-retention correlations
```

Status:

```text
RESIDUAL_DEPTH_STRUCTURE = SUPPORTED
```

if:

```text
95% CI lower bound(depth_shift_full) > 0
and
at least one predeclared technical-retention correlation has
a 95% CI excluding 0 in the expected positive direction
```

```text
RESIDUAL_DEPTH_STRUCTURE = NOT_SUPPORTED
```

if:

```text
95% CI upper bound(depth_shift_full) <= 0
and
all predeclared technical-retention correlation CIs include 0
```

otherwise:

```text
RESIDUAL_DEPTH_STRUCTURE = UNRESOLVED
```

The joint residual `μ11` is reported, not separately threshold-tuned.

## 7.7 Critical interpretation

Even if Bdepth removes the positive null:

```text
Bdepth is NOT automatically the correct future production baseline.
```

It includes technical/sparsity summaries chosen specifically to diagnose the Stage-23 null.

Promotion of any such feature into Stage 24 requires the later frozen Stage-24 plan.

---

# 8. Stage 23.2D — Outcome-label reliability

Question:

```text
Is the frozen Role-A top-gDNA-barcode outcome sufficiently stable to support
a hard 35-positive / 3,112-negative classification claim?
```

This stage studies the **measurement**, not predictive performance.

No alternate outcome is allowed to become a new prediction target inside 23.2D.

## 8.1 Exact source-rule reconstruction

Reproduce the historical label from raw/source tables exactly.

Record:

```text
gDNA grouping key
SampleNum behavior
nUMI aggregation
top-N selection
tie behavior
join to retained expression clones
special exclusions
clone-level positive aggregation
```

Require exact reproduction of:

```text
35 positive retained clones
```

before any reliability analysis.

## 8.2 Cutoff geometry

For each source selection unit, report the gDNA count/rank neighborhood:

```text
ranks 80 through 120
```

including:

```text
count
rank
tie size
gap to previous rank
gap to next rank
ratio to rank-100 cutoff
```

Report the exact rank-100 tie behavior.

No predictive outcome model is fitted.

## 8.3 Primary gDNA sampling-stability model

Use a conditional multinomial UMI resampling diagnostic.

**Scope, stated plainly (V2).** The gDNA table has exactly **one** source selection unit --
`SampleNum = 3`, all 49,554 rows, `N = 782,826` total counts over 1,936 barcodes. This model
therefore captures **sequencing-count sampling noise only**. It does not capture colony-level
biological sampling, PCR duplication, barcode collision, or assay-to-assay variation. It is a
*lower bound* on label instability, and section 8.6 is constrained accordingly.

For each independent source selection unit used by the author's top-N rule:

```text
observed barcode counts = c_i
total UMI               = N
p_i                     = c_i / N

draw:
  c* ~ Multinomial(N, p)
```

Number of resamples:

```text
5,000
seed = 23431
```

For each resample:

```text
reapply the exact source top-100 + tie rule
reapply the frozen join / exclusion rule
derive retained clone-positive set
```

For every retained clone report:

```text
P(selected positive under resampling)
```

Aggregate:

```text
mean positive-set size
distribution of positive-set size

for the frozen 35 positives:
  mean retention probability
  median retention probability
  number with P < 0.8
  number with P < 0.5

for frozen negatives:
  maximum promotion probability
  number with P > 0.2
  number with P > 0.5

Jaccard(resampled positive set, frozen positive set):
  mean
  median
  5th / 95th percentile
```

This is a measurement-sampling model, not proof of total biological label error.

## 8.4 Secondary cutoff sensitivity

Without fitting any predictor, reconstruct clone-positive sets for:

```text
top N ∈ {80, 90, 100, 110, 120}
```

using the same tie behavior.

Report:

```text
positive clone count
Jaccard vs frozen N=100 set
number of frozen positives lost
number of frozen negatives gained
```

The N=100 outcome remains historical.

No alternate N may be selected because it improves model performance.

## 8.5 Cross-GSM gDNA concordance — REMOVED IN V2

**Not executable, and not identifiable in principle.** V1 asked for
`Spearman(gDNA support in GSM7092515, gDNA support in GSM7092516)`. The V1 audit found:

```text
stepThreeStarcodeShavedReads_BC_gDNA.txt   SampleNum = 3   on all 49,554 rows
stepThreeStarcodeShavedReads_BC_10X.txt    SampleNum ∈ {1, 2}
```

The gDNA outcome is a **single pooled library**. The two control GSMs are the 10X *pre-state* RNA
samples (`SampleNum 1 → GSM7092516`, `SampleNum 2 → GSM7092515`). There is no gDNA support
attributable to either GSM, so per-GSM outcome concordance is not a hard measurement -- it is a
category error. No repair inside 23.2D is possible.

This subsection number is retained rather than renumbered so that cross-references from V1 remain
resolvable.

**Where the residual idea goes instead.** Concordance of the *lineage-barcode assay* between the
two 10X control samples is a real, computable quantity, but it measures the **input** side, not the
outcome label. It therefore moves to 23.2A / QC as a source-design descriptive, and:

```text
it carries NO outcome-label-reliability weight
it may not appear in the 8.6 status rule
it may not appear in the 10.1 ledger
it is reported only in stage23_2_source_design.json
```

23.2D's label-reliability evidence is 8.3 and 8.4 alone.

## 8.6 Outcome-label limitation status

Use both the multinomial stability and cutoff sensitivity.

```text
OUTCOME_LABEL_LIMITATION = SUPPORTED
```

if **both** are true:

```text
A. multinomial stability:
   mean retention probability of the frozen 35 positives < 0.80
   OR at least 7 / 35 frozen positives have P(selected) < 0.50

and

B. cutoff sensitivity:
   min(
     Jaccard(top90, top100),
     Jaccard(top110, top100)
   ) < 0.80
```

```text
OUTCOME_LABEL_LIMITATION = NOT_SUPPORTED
```

**Design change 3 in V2 — this status is NOT reachable from 8.3 + 8.4 alone.**

V1 allowed `NOT_SUPPORTED` whenever multinomial retention and top-N Jaccard both looked stable.
That inference does not hold. As 8.3 now states explicitly, the multinomial model captures
sequencing-count sampling noise from a single pooled library; the top-N ladder captures sensitivity
to the cutoff position. Neither observes colony-level biological sampling, assay repetition, or any
independent re-measurement of the same clones. Passing both means *"the label is stable under the
two noise sources we could actually model"* -- which is not evidence that the label is a reliable
measurement of the biological event.

Therefore:

```text
NOT_SUPPORTED requires ALL of the V1 stability criteria:

  mean frozen-positive retention >= 0.90
  at most 3 / 35 frozen positives have P(selected) < 0.80
  Jaccard(top90, top100)  >= 0.90
  Jaccard(top110, top100) >= 0.90

AND, in addition and without exception:

  genuinely independent outcome-assay replication of the SAME clones exists
  and agrees -- an independent gDNA library, a repeat colony assay, or an
  equivalent independent measurement of the same reprogramming outcome.
```

No such replication is available in the current Rewind materials. Until one exists,
`NOT_SUPPORTED` is unreachable and the correct status when the stability criteria pass is:

```text
OUTCOME_LABEL_LIMITATION = UNRESOLVED
```

Instability still yields `SUPPORTED` on the 8.3 + 8.4 criteria alone -- the asymmetry is
deliberate. These diagnostics can demonstrate that a label *is* fragile; they cannot demonstrate
that it is sound.

Otherwise:

```text
OUTCOME_LABEL_LIMITATION = UNRESOLVED
```

These thresholds are diagnostic conventions frozen before the analysis.

They do not redefine the historical target.

## 8.7 If label limitation is supported

Do **not** immediately choose a new label.

23.2D may propose candidate future outcome formulations, for example:

```text
continuous gDNA support
soft probabilistic positive membership
a margin-separated hard label
a source-author validated alternate outcome
```

but every candidate remains:

```text
EXPLORATORY_PROPOSAL_ONLY
```

until 23.2F freezes a confirmatory protocol and new evidence is obtained.

---

# 9. Stage 23.2E — Power / identifiability analysis

Question:

```text
Given the actual rare-positive design and the appropriate null,
what effect sizes are identifiable, and how much new independent evidence
would be needed to test a Role-A signal with useful power?
```

This is not a retrospective "observed power" calculation.

## 9.1 Historical test geometry

Report directly from Stage 23:

```text
observed ΔAP           +0.01050
historical null mean   +0.00350
historical null sd     +0.00631
historical null p95    +0.01455
historical p_perm      0.0846
```

Derived descriptive quantities:

```text
null-centered observed separation =
  observed - null mean

distance to historical p95 =
  observed - null p95
```

Do not convert these alone into a claim of underpowering.

## 9.2 Experimental-unit interpretation

Power analysis must use the 23.2A within-R1 status.

**Corrected in V2 (design change 4).** The benchmark is one biological replicate, R1, on both the
benchmark's own record and the GEO titles. The number of biological replicates in the current
Rewind evidence is therefore **exactly one**, and no 23.2A finding can change that number. The
within-R1 status affects only how the two samples are described:

```text
WITHIN_R1_TECHNICAL_LANES
    two 10X lanes of one suspension inside R1
    -> the 7.5 lane-composition sensitivity may run

WITHIN_R1_SEPARATE_LIBRARIES
    separately handled libraries / cultures inside R1
    -> 7.5 forbidden; report that within-R1 handling structure exists

WITHIN_R1_STRUCTURE_UNRESOLVED
    -> 7.5 forbidden; within-R1 handling described as unresolved
```

In every case the biological-replication count for power purposes is `n_biological_replicates = 1`.

## 9.2.1 Two separate power questions — design change 5

V1 collapsed detectability and replication into a single `POWER_LIMITATION`. These are different
questions with different remedies, and only one of them is answerable by simulation:

```text
Q1  WITHIN-R1 EVENT-COUNT DETECTABILITY
    Given this feature geometry, this null, and a controlled signal scale,
    how many positive clones are needed before the historical pipeline
    detects it?
    -> answerable by the 9.3 semi-synthetic study
    -> remedy: more clones / more events in the same biological context

Q2  BIOLOGICAL-REPLICATION LIMITATION
    Does the evidence rest on a single biological replicate, so that a result
    could reflect one culture's idiosyncrasy rather than a reproducible
    biological effect?
    -> NOT answerable by any simulation
    -> remedy: an independent biological replicate, i.e. new data
```

**Synthetic cohort scaling cannot create or estimate biological replicate diversity.** Resampling
clones from R1 with replacement multiplies rows inside one biological context; it does not
manufacture a second one. Every scale-2 and scale-4 cohort in 9.3 remains, biologically, replicate
R1. A power curve that reaches 0.80 at scale 4 answers Q1 only, and says nothing whatever about Q2.

Both statuses are reported separately and neither substitutes for the other.

## 9.3 Primary semi-synthetic full-pipeline power study

The purpose is planning, not confirmation.

Use the **historical Role-A pipeline** and a label-free synthetic state direction.

### Synthetic state direction

Construct one simulation-generating score `z` without using `y_primed`:

```text
1. use the 3,147 clone-level CP10K/log1p X matrix
2. regress each gene on Bdepth using all clones
3. form residualized X
4. fit one PCA on the residualized all-clone matrix
5. orient PC1 deterministically:
     find the loading with largest absolute magnitude
     ties break by stable feature ID
     force that loading to be positive
6. use oriented PC1 as synthetic state score z
7. standardize z
```

Using all clones here is allowed **only because z is a label-free simulation generator**.

It is never used to score the real Role-A outcome and never enters the fitted evaluation pipeline.

This avoids fold-specific PC sign/axis ambiguity while keeping the synthetic direction outcome-free.

### Cohort-size / event-count ladder

Preserve the historical positive prevalence approximately rather than increasing positives inside a fixed 3,147-clone cohort.

Primary synthetic cohort sizes:

```text
scale 1:
  N = 3,147
  positives = 35
  positives / fold = 7

scale 2:
  N = 6,294
  positives = 70
  positives / fold = 14

scale 4:
  N = 12,588
  positives = 140
  positives / fold = 28
```

For scales 2 and 4, construct synthetic covariate cohorts by empirical resampling **with replacement within each original outer fold**.

Each synthetic row receives a new synthetic clone ID.

The resampling preserves the empirical feature/nuisance distribution approximately and keeps five evaluation folds, but it is a planning approximation — not evidence that duplicated empirical covariates are biological replicates.

### Synthetic positive assignment

Within each synthetic outer fold, generate exactly the planned number of positive clones by weighted sampling without replacement:

```text
weight_i = exp(beta * z_i)
```

This preserves exact fold-level class counts.

Choose `beta` by deterministic numerical calibration so the **oracle z score** has the requested median AUC across a separate calibration simulation stream.

Primary target signal:

```text
oracle AUC = 0.66
```

This approximates the scale of the historical observed R3 discrimination but uses a label-free direction.

Secondary planning sensitivity:

```text
oracle AUC = 0.70
```

### Null critical values

For each cohort scale:

```text
200 synthetic-null outcome allocations
```

with positives sampled uniformly within each synthetic outer fold at the exact planned count.

For each null allocation, rerun the full historical R1/R3 nested pipeline.

This produces:

```text
q95_null(scale)
```

### Alternative simulations

For each:

```text
cohort scale
target oracle AUC
```

run:

```text
100 alternative outcome simulations
```

and rerun the full historical nested R1/R3 pipeline.

Estimated power:

```text
fraction of alternative ΔAP values > q95_null(scale)
```

Seeds:

```text
covariate-resample base seed  23440
null base seed                23441
alternative base seed         23442
beta calibration seed         23443
```

No Stage-23 observed label is used to generate simulated outcomes.

## 9.4 Runtime rule

This study may be expensive.

Process-level parallelism across independent seed streams is allowed.

Forbidden after execution starts:

```text
reducing simulation counts because results look obvious
dropping difficult m values
changing target AUC after seeing power
reusing outer-test outcomes for tuning
changing the synthetic direction to get a nicer curve
```

If runtime is infeasible, stop 23.2E and revise the plan before weakening it.

## 9.5 Power statuses — two separate labels

**Design change 5 in V2.** V1's single `POWER_LIMITATION` is replaced by two statuses. The first is
estimated from the 9.3 curve; the second is determined by the design and is not simulated at all.

### 9.5.1 Within-R1 event-count detectability

Use the primary AUC=0.66 curve.

```text
WITHIN_R1_EVENT_COUNT_LIMITATION = SUPPORTED
```

if:

```text
estimated power at scale 1 (N=3,147 / 35 positives) < 0.50
and
estimated power reaches >= 0.80 at scale 2 or scale 4
```

```text
WITHIN_R1_EVENT_COUNT_LIMITATION = NOT_SUPPORTED
```

if:

```text
estimated power at scale 1 >= 0.80
```

otherwise:

```text
WITHIN_R1_EVENT_COUNT_LIMITATION = UNRESOLVED
```

If power remains <0.80 even at scale 4, report that the tested signal scale would require more
information than the tested cohort-size range provides; do not extrapolate a precise required N
from a failed range.

The AUC=0.70 curve is secondary planning sensitivity only.

### 9.5.2 Biological-replication limitation

This status is **not** estimated from any simulation. It follows from the design:

```text
BIOLOGICAL_REPLICATION_LIMITATION = SUPPORTED
```

whenever:

```text
n_biological_replicates in the evidence used for the claim = 1
```

which is the current state of the Rewind benchmark and cannot be changed by 23.2A-23.2F. It may
become `NOT_SUPPORTED` only when a Role-A claim is supported by evidence spanning more than one
independent biological replicate -- i.e. after 23.2G succeeds on untouched confirmation evidence.

```text
BIOLOGICAL_REPLICATION_LIMITATION = NOT_SUPPORTED
    requires >= 2 independent biological replicates supporting the claim

BIOLOGICAL_REPLICATION_LIMITATION = UNRESOLVED
    is not an available value; the replicate count is known
```

A high within-R1 power curve **never** clears this status. Any projected clone-count requirement is
not a substitute for independent biological replication, and 23.2G still requires untouched
confirmation evidence.

## 9.6 What this power study does not prove

A synthetic power curve cannot prove that the true biological effect equals AUC 0.66 or 0.70.

It answers:

```text
if a state signal of this controlled scale existed in this feature geometry,
how often would the historical pipeline detect it at different rare-positive counts?
```

That is sufficient for planning new confirmation evidence, and insufficient for anything else.

Explicitly, the study **cannot**:

```text
create biological replicate diversity
    scales 2 and 4 resample clones from R1 with replacement; every synthetic
    cohort is still biological replicate R1

estimate between-replicate variance
    with one replicate there is nothing to estimate it from

establish that a corrected same-data result would replicate
    that is 23.2G's question and needs untouched evidence

be read as effective sample size
    duplicated covariate rows are not independent observations; effective N is
    below nominal N and the reported curve is therefore optimistic
```

The duplicate-row caveat must be printed next to the curve in the 23.2E artifact and record, not
only stated here.

---

# 10. Stage 23.2F — Diagnostic synthesis + confirmatory protocol freeze

23.2F fits no new diagnostic model.

It reads 23.2A–23.2E and creates the failure-decomposition ledger.

## 10.1 Multi-label diagnostic ledger

At minimum:

```text
MODEL_SELECTION_NULL_INFLATION
RESIDUAL_DEPTH_STRUCTURE
OUTCOME_LABEL_LIMITATION
WITHIN_R1_EVENT_COUNT_LIMITATION
BIOLOGICAL_REPLICATION_LIMITATION
ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA
```

Each is:

```text
SUPPORTED
NOT_SUPPORTED
UNRESOLVED
```

except `BIOLOGICAL_REPLICATION_LIMITATION`, which per 9.5.2 admits only `SUPPORTED` or
`NOT_SUPPORTED` because the replicate count is known rather than estimated.

The first five use their frozen substage rules. `OUTCOME_LABEL_LIMITATION = NOT_SUPPORTED` is
unreachable without independent outcome-assay replication (8.6), and
`BIOLOGICAL_REPLICATION_LIMITATION = SUPPORTED` holds for as long as the claim rests on one
biological replicate (9.5.2). Both asymmetries are intentional: they stop an absence of measurable
problems from being recorded as positive evidence of soundness.

`ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA` is a synthesis status based on the predeclared corrected same-data diagnostic, not a new confirmatory claim.

Use:

```text
SUPPORTED
```

only if:

```text
CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC = POSITIVE
and
OUTCOME_LABEL_LIMITATION != SUPPORTED
```

Use:

```text
NOT_SUPPORTED
```

only if:

```text
CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC = NEGATIVE
and
OUTCOME_LABEL_LIMITATION = NOT_SUPPORTED
and
WITHIN_R1_EVENT_COUNT_LIMITATION = NOT_SUPPORTED
```

Note the consequence of 8.6 and 9.5.2: `OUTCOME_LABEL_LIMITATION = NOT_SUPPORTED` requires
independent outcome-assay replication that the current materials do not contain. Until such
evidence exists, `ROBUST_STATE_SIGNAL_COMPATIBLE_WITH_DATA = NOT_SUPPORTED` is therefore also
unreachable, and a negative corrected diagnostic yields `UNRESOLVED`. This is the intended
behaviour: Stage 23.2 must not conclude "no biological signal" from the current Rewind data.

Otherwise:

```text
UNRESOLVED
```

This keeps measurement/power limitations from being misread as biological absence.

A `SUPPORTED` status means only:

```text
a robust state signal remains compatible with the already-inspected Rewind data
under the frozen diagnostic correction
```

It still requires independent 23.2G confirmation.

No hidden weighted score is allowed.

## 10.2 Required decomposition table

Report:

```text
historical observed ΔAP

historical null:
  μ00

no-K-selection null (equal-weight 3-arm mean):
  μ10

expanded-depth null:
  μ01

no-K-selection + expanded-depth null:
  μ11

per-arm diagnostics, both nuisance blocks:
  mean / CI for K = 10, 20, 50
  arm dispersion

search-width ladder:
  4 / 8 / 12 candidate paired null means

observed cells:
  O00
  O10
  O01
  O11

selection contrast
depth contrast
2×2 interaction
joint residual

corrected same-data:
  O11 - μ11
  q95_11
  p_diag_11
  CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC

outcome-label stability metrics
power curve
replicate-structure conclusion
```

## 10.3 No same-data rescue

Even if:

```text
μ11 ≈ 0
and
a corrected observed comparison looks strong
```

the result is still:

```text
MECHANISM_DIAGNOSED_ON_EXISTING_REWIND
```

not:

```text
ROLE_A_CONFIRMATORY_SUPPORTED
```

## 10.4 Freeze the corrected Role-A hypothesis

If the synthesis supports a confirmable Role-A hypothesis, define exactly one corrected confirmatory hypothesis.

Examples of what may be frozen depending on findings:

```text
historical hard label + search-matched pipeline

historical hard label + depth-complete nuisance control

a pre-specified alternative outcome representation justified by 23.2D

a combination of the above
```

The choice must be a mechanical consequence of the diagnostic ledger.

Do not choose between several corrected hypotheses based on which one gave the largest same-data Rewind effect.

If multiple corrections remain equally plausible:

```text
ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
```

and Stage 24 stays blocked.

## 10.5 Confirmatory protocol file

Before any untouched confirmation dataset/outcome is inspected for performance, create and commit:

```text
STAGE_23_2_ROLE_A_CONFIRMATION_V1.md
```

It must freeze:

```text
confirmed scientific hypothesis
outcome definition
allowed input X
nuisance block
model family / grid
preprocessing
grouping unit
primary metric
primary statistic
permutation / null design
PASS threshold
minimum positive-count / design requirement
source-qualification criteria
search budget for confirmation data
what data are forbidden because already inspected
```

## 10.6 Stage-24 handoff draft

Create:

```text
STAGE_23_2_HANDOFF_TO_STAGE_24.md
results/stage23_2/stage23_2_handoff_to_stage24.json
```

The handoff must contain:

```text
ROLE A
  historical Stage-23 failure
  Stage-23.2 diagnostic ledger
  corrected confirmatory hypothesis
  outcome version
  experimental unit
  replicate status
  nuisance variables
  primary metric / gate
  allowed claims
  forbidden claims
  evidence consumed by diagnosis
  evidence reserved for confirmation

ROLE B
  frozen Stage-23 additive verdict
  frozen Stage-23 interaction verdict
  strongest frozen simple baseline
  C1 / C2 endpoint roles
  nuisance block
  treatment coding
  treatment-level limitations / exceptions

GLOBAL
  feature universe
  benchmark versions
  split/group policy
  datasets already inspected
  datasets reserved for Stage 27
  unresolved limitations
  exact Stage-24 opening rule
```

This handoff is not marked Stage-24-ready until 23.2G confirms Role A.

## 10.7 Material benchmark-change firewall

If the proposed correction changes any material benchmark semantic:

```text
outcome definition
positive / negative ontology
experimental unit
source reconstruction
leakage firewall
split / grouping rule
```

then Stage 23.2 cannot directly confirm and hand off the changed benchmark.

Required route:

```text
version revised benchmark
rerun affected Stage-22 qualification/construction contracts
rerun affected Stage-23 learnability + structural + permutation gates
then return to confirmation
```

This prevents Stage 23.2 from bypassing the original benchmark gates by redefining the problem after failure.

---

# 11. Stage 23.2G — Independent confirmation / roadmap resolution

23.2G runs only after:

```text
23.2F complete
STAGE_23_2_ROLE_A_CONFIRMATION_V1.md committed
confirmation evidence remains untouched by diagnostic design
```

## 11.1 Confirmation evidence eligibility

Preferred evidence:

```text
1. a new independent biological Rewind replicate
2. an independent Rewind-like dataset with pre-state RNA and future lineage outcome
3. a prospectively held-out clone cohort never inspected during Stage 23.2 design
4. another independent future-outcome anchor satisfying the same scientific claim
```

**The 5.2.1 reserved ledger is the first place to look under category 1.** GSE227151 declares
biological replicates 2 and 3 (`GSM7092517`-`GSM7092521`) that Stage 23 never used. At 23.2F, and
only then, evaluate each ledger entry against 11.2's eligibility criteria.

Two cautions are pre-registered now:

```text
declared metadata is not a usable outcome
    whether reps 2/3 carry a reconstructable gDNA / colony outcome for the same
    claim is UNVERIFIED. The ledger records what GEO declares, nothing more.

gated samples are a different population
    GSM7092520 / GSM7092521 are sorted for fast / slow cycling. Sorting changes
    the population, so they may fail the "same scientific claim" requirement even
    if their data are complete.
```

If a ledger entry qualifies, it becomes confirmation evidence under 11.3 and is thereafter subject
to the 11.4 Stage-27 firewall.

The evidence must not have been used to:

```text
choose the correction
choose the model/grid
choose the label definition
choose the nuisance block
choose the metric
choose the PASS threshold
```

## 11.2 Confirmation-source qualification

If no reserved confirmation set already exists, perform a bounded search only after the confirmatory protocol is frozen.

Eligibility must require:

```text
pre-intervention molecular measurement
independently measured later outcome
clone / lineage / biological-unit linkage sufficient for grouped evaluation
raw or processed data sufficient to reconstruct the endpoint
no same-state outcome leakage
enough positives / events to meet the 23.2E planning minimum
```

If no qualifying dataset is found inside the frozen search budget:

```text
ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
```

Stage 24 remains blocked.

Do not lower inclusion criteria because search results are sparse.

## 11.3 Confirmation analysis

Run exactly the frozen confirmation protocol.

No architecture search.

No multiple corrected hypotheses.

No post-hoc threshold change.

No alternate endpoint rescue.

The confirmation result must have its own structural/leakage audit and exact null test as frozen in `STAGE_23_2_ROLE_A_CONFIRMATION_V1.md`.

## 11.4 Stage-27 firewall

Any dataset used in 23.2G to choose or reopen Stage 24 is now development/confirmation evidence.

It is **not** an untouched Stage-27 independent replication set.

Stage 27 must preserve another independent biological test of the eventual Stage-24 model.

## 11.5 Final roadmap statuses

Exactly one final roadmap status:

### `ROLE_A_CONFIRMATORY_SUPPORTED`

Require:

```text
new pre-registered independent confirmation PASS
structural controls PASS
no material benchmark semantic left un-re-gated
Stage-24 handoff complete
```

Then:

```text
STAGE 24 MAY OPEN
```

### `ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE`

Use when:

```text
failure mechanism is plausible / diagnosable
but no untouched confirmation has passed
```

Then:

```text
STAGE 24 REMAINS BLOCKED
```

### `ROLE_A_REDESIGN_REQUIRED`

Use when:

```text
current outcome / experimental unit / benchmark semantics are inadequate
for the intended Role-A claim
```

Then:

```text
version + re-gate benchmark
obtain new evidence
STAGE 24 REMAINS BLOCKED
```

### `ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL`

Use when:

```text
after decomposition there is no defensible Role-A signal worth preserving
as the mandatory prospective anchor
```

Then:

```text
do not open Stage 24 under the current roadmap
explicitly revise the roadmap / claim architecture
decide separately whether Role B becomes the primary scientific anchor
```

---

# 12. Decision matrix

Stage 23.2 must not collapse the diagnostics into a simplistic single-cause rule.

Examples:

```text
selection SUPPORTED
depth NOT_SUPPORTED
label UNRESOLVED
within-R1 event count SUPPORTED
biological replication SUPPORTED
    -> plausible model-selection + event-count explanation
    -> freeze one corrected confirmatory hypothesis
    -> require new confirmation

selection SUPPORTED
depth SUPPORTED
label UNRESOLVED
within-R1 event count SUPPORTED
biological replication SUPPORTED
    -> multiple methodological / event-count mechanisms
    -> same-data corrected result remains exploratory
    -> require new confirmation

label SUPPORTED
regardless of model-selection/depth
    -> hard-label claim is measurement-limited
    -> outcome redesign may be necessary
    -> material outcome change triggers benchmark re-gating

selection NOT_SUPPORTED
depth NOT_SUPPORTED
label NOT_SUPPORTED          <- requires independent outcome-assay replication (8.6)
within-R1 event count NOT_SUPPORTED
biological replication NOT_SUPPORTED   <- requires >= 2 biological replicates (9.5.2)
and no corrected signal remains
    -> strong evidence toward NO ROBUST ROLE-A SIGNAL

    NOTE: with the current Rewind materials this row is UNREACHABLE, because two
    of its five conditions require evidence the benchmark does not contain. That
    is deliberate. Stage 23.2 cannot conclude "no biological signal" from one
    biological replicate measured by one un-replicated outcome assay.

all major mechanisms UNRESOLVED
    -> ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE
    -> do not invent a rescue
```

---

# 13. What Stage 23.2 may establish

Stage 23.2 may establish:

```text
whether model-selection width contributes to the positive null center
whether residual technical depth contributes to the positive null center
whether the hard gDNA-derived label is unstable under declared diagnostics
whether the rare-positive design has low within-R1 detection power for a controlled signal scale
what the within-R1 sample structure of the two Rewind control GSMs is, if source evidence resolves it
what corrected Role-A hypothesis is scientifically defensible enough to confirm
what new independent evidence is required
whether Stage 24 may reopen
```

---

# 14. What Stage 23.2 cannot establish

Stage 23.2 cannot establish from the existing Rewind data alone:

```text
a new final Role-A PASS
causal biological mechanism
clinical / laboratory utility
unseen-treatment generalization
independent biological replication of the final Stage-24 model
calibrated production probabilities
OOD performance
treatment ranking
```

It also cannot claim:

```text
"no biological signal"
```

merely because the historical Stage-23 permutation failed.

That stronger statement requires the failure decomposition to exclude plausible methodological,
measurement, and power explanations -- and per 8.6 and 9.5.2, two of those exclusions are not
achievable from the current Rewind materials at all. Stage 23.2 is therefore expected to end with
`UNRESOLVED` on at least one axis, and that is a correct result, not a failed stage.

Additionally, Stage 23.2 cannot establish:

```text
that biological replicates 2 and 3 of GSE227151 contain a usable Role-A outcome
    the 5.2.1 ledger records declared metadata only; nothing is downloaded,
    inspected, or evaluated before 23.2F freezes the confirmation protocol
```

---

# 15. Engineering deliverables

Suggested committed artifacts:

```text
STAGE_23_2_ROLE_A_RESOLUTION_V2.md          (V1 archived under arcive/)

results/stage23_2/stage23_2_protocol.json
results/stage23_2/stage23_2_source_design.json
results/stage23_2/stage23_2_historical_null_d00.json          <- new in V2 (4.3, 5.5)
results/stage23_2/stage23_2_reserved_confirmation_candidates.json  <- new in V2 (5.2.1)
results/stage23_2/stage23_2_model_selection_decomposition.json
results/stage23_2/stage23_2_depth_decomposition.json
results/stage23_2/stage23_2_label_reliability.json
results/stage23_2/stage23_2_power_identifiability.json
results/stage23_2/stage23_2_diagnostic_synthesis.json
results/stage23_2/stage23_2_handoff_to_stage24.json

stage_23_2A_RECORD.md
stage_23_2B_RECORD.md
stage_23_2C_RECORD.md
stage_23_2D_RECORD.md
stage_23_2E_RECORD.md
stage_23_2F_RECORD.md
stage_23_2G_RECORD.md
stage_23_2_RECORD.md

STAGE_23_2_ROLE_A_CONFIRMATION_V1.md
STAGE_23_2_HANDOFF_TO_STAGE_24.md
```

`STAGE_23_2_ROLE_A_CONFIRMATION_V1.md` is created only at 23.2F if a confirmable hypothesis exists.

If no hypothesis is defensible:

```text
do not create a fake confirmation protocol
```

Large cache-only artifacts:

```text
_cc_cache/stage23_2/
  recovered permutation mappings          (digest committed; array itself is not)
  2×2 per-permutation diagnostic arrays   (including all three fixed-K arms)
  multinomial gDNA bootstrap draws
  semi-synthetic power simulation draws
```

The 200 historical `D00` values are the one exception: they are compact and load-bearing, so per
4.3 they are committed rather than cached.

Committed JSON stores their content hashes and compact summaries.

Do not commit giant clone×gene matrices or simulation dumps.

---

# 16. Required automated contracts

## 16.1 23.2A

Tests must assert:

```text
Stage-23 historical artifacts unchanged
Stage-23 closure record remains closed
Role-A historical final = FAIL
Role-B historical positives unchanged
ROADMAP GATE remains STAGE_24_BLOCKED_ROLE_A
legacy `STAGE 23R` is accepted only as a historical alias for Stage 23.2

the committed D00 array has exactly 200 values and its digest matches
the replayed D00 array reproduces the committed Stage-23 summary statistics
  (mean, sd, p95, max, count >= observed, p_perm) exactly
the reserved confirmation ledger contains metadata fields ONLY --
  no expression, barcode, outcome or performance quantity for any reserved GSM
the source-design status is one of the three WITHIN_R1_* values;
  REPLICATE_STRUCTURE_BIOLOGICAL must not appear anywhere in Stage-23.2 output
n_biological_replicates recorded for the current benchmark = 1
3,147 / 35 / 3,112 counts exact
7 positives per outer fold
B0 exact
Bdepth exact
all Bdepth fields finite
no gDNA field enters Bdepth
historical 200 permutation mappings recover exactly
train/test boundaries preserved in every mapping
protocol canonical JSON digest stable across LF/CRLF source text
scientific protocol hash does not depend on whole builder file hash
```

## 16.2 23.2B

Tests must assert:

```text
D00 is historical full-search null
D10 uses exactly K=20 + four C values
R1 candidate path unchanged
same permutation_id paired in every contrast
paired CI uses permutation_id resampling
status rule is mechanical
conditional width ladder runs only when primary status UNRESOLVED
```

Mutation tests should make the status fail if:

```text
K=10 or K=50 is silently substituted for K=20
a permutation mapping is regenerated with a new seed
outer test crosses into train
```

## 16.3 23.2C

Tests must assert:

```text
Bdepth has exactly four declared fields
Bdepth values are pre-outcome
μ00/μ10/μ01/μ11 formulas exact
O00/O10/O01/O11 formulas exact
O00 reproduces historical observed ΔAP
p_diag_11 uses the 200 D11 draws with +1 correction
CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC rule exact
factorial contrast formula exact
technical-retention diagnostics use no y
sample-specific lane counts are forbidden unless 23.2A says TECHNICAL_LANES
```

## 16.4 23.2D

Tests must assert:

```text
historical top-100 reconstruction yields exactly 35 positives
multinomial resampling preserves each source selection unit's total UMI
5,000 draws
seed exact
tie behavior exact
alternate N grid = 80/90/100/110/120
no predictive estimator is fitted on alternate labels
status thresholds exact
```

A grep/AST contract should fail if 23.2D routes alternate labels into the Role-A model-fitting function.

## 16.5 23.2E

Tests must assert:

```text
synthetic z construction never reads y_primed
global simulation-generator PC1 orientation deterministic
target AUC values exact
cohort scales N=3147/6294/12588 exact
positive counts 35/70/140 exact
covariate resampling stays within original outer fold
synthetic clone IDs are unique
null/alternative simulation counts exact
fold positive allocation deterministic
historical full nested pipeline rerun
power status rule exact
no "observed power" shortcut
```

## 16.6 23.2F

Tests must assert:

```text
no estimator fit
all substage statuses read from frozen artifacts
no weighted hidden synthesis score
same-data diagnostic result cannot emit ROLE_A_CONFIRMATORY_SUPPORTED
handoff includes every mandatory field
material benchmark-change flag blocks direct Stage-24 opening
```

## 16.7 23.2G

Tests must assert:

```text
confirmation protocol committed before performance read
confirmation evidence not in diagnostic-consumed-data set
only one frozen confirmatory hypothesis
PASS rule exact
Stage-27 reserved evidence not consumed
Stage 24 opens only on ROLE_A_CONFIRMATORY_SUPPORTED
```

## 16.8 Determinism

Fresh clean-tree regeneration must verify every committed Stage-23.2 result artifact.

The determinism registry may grow.

A test must assert:

```text
it can never shrink without an explicit migration
```

Scientific protocol digests must remain stable when later substage implementation code is appended but `stage23_2_protocol.json` is unchanged.

---

# 17. Suggested implementation boundary

Prefer a dedicated Stage-23.2 experiment package rather than extending the Stage-23 builder indefinitely.

For example:

```text
experiments/
  stage23_2/
    protocol.py
    source_design.py
    selection_null.py
    depth_null.py
    label_reliability.py
    power.py
    synthesis.py
    confirmation.py
```

or an equivalent repository-native structure.

The exact filenames may follow repository conventions.

The required design property is:

```text
scientific protocol identity is separate from mutable implementation identity
```

Do not route these diagnostic experiments through production `CellFateNet`.

---

# 18. Runtime / resource expectations

23.2B/C intentionally reuse the same 200 historical permutation mappings.

Likely expensive components:

```text
23.2B no-K-selection null (three fixed-K arms)
23.2C full + no-K-selection Bdepth nulls
23.2E semi-synthetic full-pipeline power study
```

**Measured basis (V1 audit, this machine: 8 cores).** One historical null draw costs **0.48 min**
uncontended. The Stage-23 figure of 1.55 min/draw reflected three-way CPU contention, so it should
not be used for planning.

```text
23.2A  replay all 200 historical draws                    ~96 min
       (mapping regeneration alone needs no model fitting)

23.2B  cell 10, three fixed-K arms x 4 C = 12 fits        ~96 min
       NOTE: design change 1 makes cell 10 cost the SAME as a full
       search, not one third. V1's single K=20 arm would have been ~43 min.
       search-width ladder, 8-candidate rung                ~64 min

23.2C  cell 01, full search + Bdepth                      ~96 min
       cell 11, three fixed-K arms + Bdepth               ~96 min

       23.2B + 23.2C subtotal                          ~350-450 min  (6-7.5 h)
       cell 00 is NOT recomputed: it is the committed D00 array

23.2E  200 nulls x 3 cohort scales
       + 100 alternatives x 2 target AUCs x 3 scales
       = 1,200 full nested pipeline runs
       scale 1 ~0.48 min, scale 2 ~1.0 min, scale 4 ~2.1 min per run
                                                       ~1,350 min (~22 h)
                                                       single-threaded,
                                                       ~8 h at 3-way parallelism
       plus beta calibration draws
```

**Memory budget**, from the measured per-fold retained-gene counts (13,589-13,627):

```text
scale 1   N =  3,147   train  2,517   dense train block   261 MB   (~523 MB with the scaled copy)
scale 2   N =  6,294   train  5,035   dense train block   523 MB   (~1.0 GB)
scale 4   N = 12,588   train 10,070   dense train block  1,045 MB  (~2.1 GB peak per inner fit)
```

Additionally, 9.3 step 2 regresses every gene on `Bdepth` over the **unfiltered** 3,147 x 36,601
matrix: roughly 921 MB dense before any filtering. Budget for it explicitly or stream the
regression gene-block-wise.

23.2E is the dominant cost and the only substage where the frozen design may prove infeasible. If
it does, 9.4 applies: stop and revise the plan rather than weakening it.

Parallel execution across independent permutation/simulation IDs is allowed.

Do not parallelize operations in a way that changes deterministic seeds or floating-point reduction order inside one fitted comparison.

Every long-running substage should support:

```text
resume from completed permutation/simulation IDs
validate cache digest before resume
refuse mixed-protocol cache
```

No early stopping based on scientific result.

---

# 19. Stage-24 readiness contract

Stage 23.2 does not need to guarantee that Role A can be rescued.

It must guarantee that **the next decision is unambiguous**.

Stage 24 may open only when all are true:

```text
final Stage-23.2 status = ROLE_A_CONFIRMATORY_SUPPORTED

confirmation evidence:
  independent
  pre-registered
  structural controls PASS
  primary null / inference gate PASS

benchmark semantics:
  unchanged from already re-gated version
  OR materially changed version has completed required Stage-22/23 re-gating

STAGE_23_2_HANDOFF_TO_STAGE_24.md complete
stage23_2_handoff_to_stage24.json complete

Stage-27 independent replication evidence remains reserved
```

Then the next required file is:

```text
STAGE_24_PROSPECTIVE_MODEL_V1.md
```

before any Stage-24 model fit.

If these conditions are not met:

```text
Stage 24 does not start.
```

---

# 20. Completion report format

The final `stage_23_2_RECORD.md` must report at least:

```text
PLAN
  plan version
  plan commit
  protocol SHA-256
  execution commits

HISTORICAL ANCHOR
  Stage-23 Role-A observed ΔAP
  bootstrap CI
  permutation null mean/sd/p95/max
  permutation p
  historical final verdict

23.2A
  replicate-structure status
  exact source evidence
  permutation mapping recovery
  Bdepth audit
  protocol verdict

23.2B
  μ00
  μ10
  selection shift + CI
  status
  width-ladder result if triggered

23.2C
  μ01
  μ11
  O01
  O11
  depth shift + CI
  2×2 interaction
  corrected same-data q95 / p
  technical-retention correlations
  status

23.2D
  exact label reconstruction
  cutoff geometry
  multinomial label stability
  top-N sensitivity
  cross-GSM concordance
  status

23.2E
  historical identifiability geometry
  replicate interpretation
  power at cohort scales 1 / 2 / 4 for AUC=0.66
  corresponding N / positive counts
  AUC=0.70 sensitivity
  status

23.2F
  multi-label diagnostic ledger
  corrected hypothesis or "none"
  benchmark-change flag
  confirmation protocol path/hash
  Stage-24 handoff draft path/hash

23.2G
  confirmation evidence source
  qualification result
  primary confirmation statistic
  null/inference result
  structural controls
  final roadmap status

ENGINEERING
  tests
  lint
  determinism
  cache hashes
  deviations
  files modified

FINAL
  exact Stage-24 gate
  exact next action
```

---

# 21. Stop conditions

Stop immediately and do not reinterpret the plan if any occurs:

```text
historical Stage-23 null cannot be reproduced / paired
permutation mappings cannot be recovered exactly
gDNA source rule cannot reproduce 35 positives
source evidence contradicts the benchmark reconstruction
a proposed label change materially changes benchmark semantics
power simulation cannot execute as frozen
confirmation evidence was already used to design the correction
Stage-27 reserved evidence would be consumed
a diagnostic requires looking at Stage-24 performance
```

The correct response to a stop condition is:

```text
record the blocker
create a new plan version if needed
do not improvise a rescue
```

---

# 22. Final Stage-23.2 principle

The successful outcome of Stage 23.2 is **not**:

> "we made Rewind significant."

The successful outcome is:

> **we know why the mandatory prospective anchor failed, we know which parts of the failure are methodological versus measurement/power limitations, and we have frozen the exact independent evidence required to justify — or reject — reopening Stage 24.**

That is the standard required before CellFate-Rx proceeds to the prospective architecture stage.
