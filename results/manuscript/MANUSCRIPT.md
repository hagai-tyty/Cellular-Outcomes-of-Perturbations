# Pretreatment transcriptional state carries condition-specific information about future clonal detection in a lineage-traced melanoma line

**CellFate-Rx, Generation 1.**

**Author:** Hagai Aviv · ORCID [0009-0004-4503-6629](https://orcid.org/0009-0004-4503-6629)

**Affiliation:** Independent researcher, Ma'ale Adumim, Israel

**Correspondence:** hagai.aviv.home@gmail.com

```text
  evidence lock   60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa
  claim lock      fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61
```

Both digests are verifiable from the repository. See **Availability of data and materials**.

---

## Abstract

**Background.**
Whether a cell's molecular state before a perturbation predicts what happens to it afterwards is
usually asked retrospectively, after outcome and state have been measured in the same cells. We ask
it prospectively at clone level, in one BRAF-V600E melanoma cell line (WM989, GSE279162), in which
a barcoded population was split across six observed experimental conditions: Acid, Cisplatin, CoCl2,
Dabrafenib, Doxorubicin and Trametinib. Of the clones that experiment recovered, 1,401 carry a
pretreatment profile and are therefore analysable prospectively; those are the clones used here.
The ranking test uses the 892 of them detected under at least one condition and undetected under at
least one other.

**Results.**
Within this system, pretreatment gene expression contains condition-specific information about
future clonal detection beyond condition identity and captured pretreatment clone abundance, under
clone-held-out evaluation whose folds, features and exclusions were fixed
before any model was fitted. Under a test preregistered in full — the metric, population, weighting,
comparator, null and verdict rule all fixed by digest before any ranking statistic was computed,
though after earlier predictive analyses of the same data — a frozen state-by-condition interaction
model improves clone-specific ordering of the six conditions over a non-interactive additive model:
+0.051605 in equal-clone-weighted within-clone AUROC, 95% CI [+0.037197, +0.065571], with 0 of 1000
full-refit permutation draws reaching the observed value (p < 0.001). The additive model did not
itself improve that ordering over condition identity alone (0.692176 against 0.692654): the gain
is the interaction.

**Conclusions.**
A state contribution shared additively across all conditions cannot change their ordering within a
clone. Allowing state effects to vary by condition improved that ordering in WM989: part of what a
clone's pretreatment state carries bears on which of the six conditions it is still detected after,
not only on how detectable it is overall. The outcome is a detection proxy and is not death,
sensitivity, resistance or clinical response. The design and the artifacts behind it are released,
so that the result can be reproduced from the public data and tested in other lineage-traced
systems, where it has not yet been tested.

---

## Keywords

clonal barcoding; lineage tracing; drug resistance; preregistration; reproducibility;
within-clone ranking; permutation test; melanoma; reanalysis

---

## Background

The question is easy to state and hard to evaluate honestly: before a perturbation is applied, does
a cell's transcriptional state say anything about which perturbation that cell is still detected
after?

Barcoded lineage tracing solves the hardest part of asking it. A clone is split, part of it is
profiled before anything is done to it, and the remaining parts are exposed to different
conditions, so the pretreatment profile is genuinely prior to the outcome rather than a consequence
of it.

**That design, and the data used here, are not ours.** Schaff et al. built exactly such a system
across six conditions in parallel and showed that pre-existing state predicts which clones go on to
resist, identifying high CD44 expression in treatment-naive cells as a marker of resistance across
multiple conditions [1]. This work is a reanalysis of their data and takes no credit for the
experiment.

What that establishes is a **general propensity**: some clones resist more of the conditions than
others do. It leaves a different question open — whether pretreatment state says *which* condition a
given clone is still detected after, rather than how resistant that clone is overall. The two are
separable, and separating them is what this work does.

Two things make the ordering question hard to evaluate rather than merely to observe. The first is
that clone abundance dominates: a clone that was large before treatment is more likely to be
detected after it for reasons that have nothing to do with state, so any comparison that does not
hold abundance fixed will find a signal that is really a headcount. The second is that
clone-specific ordering is not something a model can produce by being generally right about a clone
— it requires an explicit interaction between state and condition, because a state contribution
shared additively across the six conditions shifts all six of a clone's scores together and leaves
their order untouched.

This work evaluates exactly that, once, under a protocol frozen in advance.

### Relation to prior work

The finding that pre-existing, non-genetic single-cell state predicts which cells resist therapy is
established in this system and is not claimed here. Shaffer et al. showed that rare transcriptional
states in WM989 predict which cells resist vemurafenib and are stabilised by drug exposure [2];
Emert et al. resolved substructure within those rare states and linked it to distinct resistant
outcomes [3]; Goyal et al. showed that clonal fates after drug are largely predetermined by
pre-treatment molecular differences and are diverse rather than binary [4]; and Schaff et al.
extended clonal tracing to six conditions in parallel, reporting cross-condition resistance
correlation and CD44 as a marker of resistance across several of them [1].

**We do not claim to be the first to look at condition-specific expression in this system.** Schaff
et al.'s own deposited analysis reports condition-associated markers and signatures alongside the
cross-condition results, and earlier work resolves substructure that distinguishes resistant
outcomes [4,5]. What has not been done is the specific thing tested here: a **frozen,
clone-held-out, preregistered test of whether an explicit state-by-condition interaction improves
clone-specific ORDERING** over a non-interactive additive model, with captured pretreatment
abundance held fixed in every model including the null.

The distinction the design turns on is not rhetorical. A state contribution shared additively across
the six conditions — including whatever part of a general resistance propensity, or of a per-clone
marker level, is the same for each of them — shifts all six of that clone's predicted scores
together and therefore contributes **exactly zero** to a within-clone ordering metric. The
comparator was chosen to enforce that separation before any ranking statistic was computed, and the
outcome is that the additive state term contributes nothing to ordering while the interaction
contributes all of the gain. The contribution is the evaluation and its preregistration, not the
observation that state carries condition-relevant information.

The methodological posture is borrowed rather than invented. Kapoor and Narayanan catalogue eight
kinds of leakage across 294 papers in seventeen fields, and observe that complex models frequently
fail to beat logistic regression once the leakage is corrected [5]. That is the failure mode this
design is built against: the comparator is a simpler model of the same family, every preprocessing
step is refitted inside the training fold, the permutation null refits the whole pipeline rather
than shuffling labels, and the metric, population, comparator and verdict rule were fixed in a
digest-frozen protocol before any ranking statistic was computed. That protocol was written after,
and in the light of, earlier predictive analyses of the same data, which the stage records keep in
full.

---

## Methods

### Data

```text
ROLE B, primary       GSE279162 (WM989)      benchmark, model, tool and ranking analysis
ROLE A, supporting    GSE227151 (Rewind)     historical supporting evidence only
```

Role B is one BRAF-V600E melanoma cell line. The source experiment recovered many thousands of
barcoded clones; **1,401 of them carry a pretreatment observation**, which is what a prospective
clone-level question requires, and those 1,401 are the analysis population here. It carries the
whole primary claim.

**Role B is a reanalysis. We generated no new data.** The experiment was performed by Schaff et
al. [1], who introduced a barcode library into WM989 A6-G3, isolated 350,000 uniquely barcoded
cells, expanded them for approximately six doublings, and divided the population across twelve
treatment arms — two replicates of each of the six conditions reanalysed here. Dabrafenib,
trametinib, CoCl2 and acidic media were applied continuously; cisplatin and doxorubicin were applied
for a treatment period followed by a recovery period, each arm spanning four weeks in total. Exact
concentrations and schedules are in [1].

Role A is a separately reconstructed reprogramming system [6]. It contributes one supporting
sentence and nothing else; its own confirmation gate failed. It is not a replication of Role B, and
it does not provide the same multi-condition task or the same outcome.

No additional dataset was searched, downloaded, qualified or used. Raw sequencing data is not
vendored; accessions are given above. **Figure 1** summarises the design and the evaluable
population.

### Benchmark construction

Clone assignments, outcome construction, the five outer folds, feature rules, condition aliases and
exclusions were all fixed before any model was fitted, and were not revised afterwards.

Two endpoint families were built. **C1** is post-treatment clone detection: an observed zero means
no assigned post-treatment cell was seen for that clone-condition row. **C2** is a clone-balanced
abundance endpoint. Everything reported below is C1.

Expression is clone-level pseudobulk: raw pretreatment counts are summed over a clone's cells, then
CP10K-normalised and `log1p`-transformed exactly once. Applying the transform twice, or summing
already-normalised cells, produces a different feature space and a model the benchmark never
evaluated.

### Models

**Table 1** The three preregistered model specifications.

| Model | Terms | What it adds |
|---|---|---|
| W1 | B + U | nuisance and condition identity only |
| W4 | X + B + U | plus an additive expression term |
| W5 | X + B + U + X*U | plus an explicit state-by-condition interaction |

`X` is the clone expression profile reduced to 50 principal components on a train-only basis. `B`
is the captured-abundance nuisance block: `log1p` cell counts, total and per pretreatment library.
`U` is condition identity, encoded as five non-reference indicators with Acid as reference. The
design has 309 columns: 50 components, 4 nuisance terms, 5 indicators and 250 interaction terms.

`B` is not optional. Abundance is the confounder that makes a naive version of this question
trivial, so it is in every model including the null, and the tool refuses to score without it.

### Evaluation

Five outer folds, held out **by clone**, so no model ever scores a clone it trained on. Gene
filtering, the PCA basis and every scaler are refitted inside each training fold. Hyperparameters
are selected by an inner grouped split within the training folds only.

### The preregistered ranking test

The primary question is whether W5 orders the six conditions *within a clone* better than W4. The
metric is equal-clone-weighted within-clone AUROC: for each clone, the mean over positive/zero
condition pairs with ties scoring exactly 0.5; then a plain mean over clones, so a clone
contributing many pairs cannot outweigh one contributing few.

The evaluable population is clones with at least one detected and one undetected condition; a clone
with no contrast has an undefined within-clone AUROC. This yields 892 of 1,401 clones — 472 were
never detected under any condition and 37 were always detected.

The comparator is W4, not W1. An additive expression term cannot by construction create
clone-specific ordering, so W4 isolates the interaction as the only thing that could.

The null is a **full refit**: profiles are permuted within stratum on each side of the outer-fold
boundary and never across it, and the entire pipeline — filtering, PCA, scaler fitting and
hyperparameter selection — is re-run inside every draw. Observed-data hyperparameters are never
reused. 1000 draws, no early stopping, with a completeness assertion that refuses an incomplete
null rather than silently reporting a smaller one.

Uncertainty on the observed statistic is a 2,000-replicate clone bootstrap, conditional on the
fitted models.

### AI use disclosure

Generative AI assistants were used throughout this work. Claude Opus 5 and Claude Opus 4.8
(Anthropic) wrote most of the code and documentation. Gemini 3.1 Pro (Google), Claude Haiku 4.5
(Anthropic), GPT 5.6 sol and GPT 6 astra (OpenAI) were used to check that work. The study's logic
and direction came mainly from the author. Each protocol was fixed by cryptographic digest before
the statistics it governs were computed, and every number reported here is produced by the frozen
code and traced mechanically to a locked artifact. The author reviewed the work and takes full
responsibility for its content.

---

## Results

### The interaction improves clone-specific ordering

**Table 2** Ranking score by model, and the preregistered difference; **Table 1** defines the
models.

| Model | Ranking score R | Terms |
|---|---|---|
| W1 | 0.692654 | nuisance + condition |
| W4 | 0.692176 | + additive X |
| W5 | 0.743781 | + explicit X x U |
| delta_RANK | +0.051605 | W5 minus W4; CI95 [+0.037197, +0.065571] |

**Figure 2** shows the three models and the observed statistic against its null.

`R(W4)` sits *below* `R(W1)` by 0.0005. The additive expression term contributes nothing to
ordering, which is precisely why W4 was preregistered as the comparator. The entire ordering gain is
the interaction.

**This also says what the gain is not.** A state contribution shared additively across all six
conditions cannot change their ordering within a clone: it moves all six of that clone's scores
together, and within-clone AUROC compares those six scores, so such a contribution counts **exactly
zero** here by construction. W4 carries exactly that contribution, and it adds nothing to the
ordering. Allowing state effects to vary by condition, as W5 does, is what improved it. That does
not rule out one shared biological programme behind those effects: a multi-condition signal of the
kind CD44 marks in this system [1] would be such a programme, and one whose influence differs
between conditions would appear as an interaction too. What this measurement locates is
condition-dependent state effects, not the number of programmes behind them.

### The separation, not the p-value, is the result

**Table 3** The observed statistic against its permutation null.

| Quantity | Value |
|---|---|
| null p95 | 0.008672 |
| observed / null p95 | 6.0x |
| observed vs null mean | 11.8 null standard deviations above it |
| largest of 1,000 draws | 0.013722 -- the observed value exceeds EVERY null draw |
| draws reaching observed | 0 of 1000 full-refit permutations |

`p < 0.001` is the floor of a 1,000-draw permutation test and is reported as such, never as a point
estimate. The number that carries weight is the separation: nothing the null produced came close.

### It holds in every stratum it was broken down by

**Table 4** delta_RANK within each stratum, by outer fold and by pretreatment depth.

| Outer fold | delta_RANK | Pretreatment depth | delta_RANK |
|---|---|---|---|
| fold 0 | +0.0435 | 1 cell | +0.0535 |
| fold 1 | +0.0548 | 2 cells | +0.0528 |
| fold 2 | +0.0658 | 3-4 | +0.0314 |
| fold 3 | +0.0506 | 5-9 | +0.0462 |
| fold 4 | +0.0435 | 10+ | +0.0779 |

Positive in all five folds and all five depth strata (**Figure 3A, 3B**). These were preregistered
as descriptive and could not have rescued a failed primary gate; they were not asked to.
The breakdown is by fold and by depth only; across the six conditions the interaction is not uniform
(Limitation 3).

### Choosing the lowest-scoring condition

**Table 5** The top-choice diagnostic, preregistered as a directional-consistency check.

| Diagnostic | Value | 95% CI |
|---|---|---|
| delta_TOP1 | +0.115471 | [+0.082960, +0.145740] |

Selecting each clone's lowest predicted detection score finds a genuine zero for 82.8% of evaluable
clones under W5 against 71.3% under W4 (**Figure 3C**). This was preregistered as a directional-consistency check,
not a significance test: it could withhold support, never grant it. It did not withhold.

### The tool

A frozen predictor ships with this work. For one starting clone it returns a
`future_detection_score` for each of the six observed conditions, reproducing the frozen
out-of-fold predictions to within 5e-16.

What it refuses is as much of the specification as what it returns:

```text
  an unknown condition        -> UNSUPPORTED_TREATMENT, and no score
  a missing nuisance block    -> MISSING_REQUIRED_NUISANCE, never imputed
  a wrong feature schema      -> UNSUPPORTED_FEATURE_SCHEMA
  a validated ordering        -> withheld unless the preregistered verdict file is supplied
```

The refusal was tested adversarially rather than asserted. 56 of 56 hostile condition strings were
refused: case and whitespace variants, dose formats, unicode look-alikes, controls, and sixteen real
oncology drugs including Vemurafenib — the drug for this exact mutation — and Carboplatin, one
substitution from a condition that *is* supported.

That test found something worth stating. Acid is the reference level and is encoded as five zero
indicators, so the indicator encoder alone maps *any* unrecognised string to the Acid row and would
return the Acid score under another name. The vocabulary filter is the only thing preventing this,
and it was verified to hold rather than assumed to.

---

## Discussion

The result is narrow by design, and this section keeps apart three things that are easy to run
together: what the result cannot support, what may not be said about it in any form, and what would
actually test it. The first two are separate on purpose; neither should have to be inferred from the
other.

### Limitations

Carried verbatim from the preregistered verdict.

```text
  1  No independent biological replication of the Role-B finding. Clone-held-out folds and
     two endpoint families are not replication.
  2  Captured pretreatment clone abundance remains ~3.45x the whole state contribution.
     The ordering is abundance first, then condition-specific state.
  3  Four of six conditions carry meaningful interaction. Cisplatin is negligible on C1
     and Doxorubicin is negative on both endpoint families.
  4  C1 is an observed detection proxy -- not death, sensitivity or clinical response.
  5  The bootstrap interval is conditional on the fitted models; only the null refits.
  6  Role A remains positive-but-underpowered supporting evidence. Its confirmation gate
     18.3 FAILED at 0.64 against a 0.80 threshold, and a later audit of the instrument
     put the true power at 0.45 -- lower still. Its effect size must not be quoted as an
     estimate, because an underpowered design that reaches significance inflates it.
```

Limitation 2 is the one most likely to be misread. This work does not show that state dominates
outcome. It shows that state adds something specific, on top of an abundance term that is several
times larger.

Limitation 3 has a pattern that this experiment cannot interpret. Schaff et al. chose the six
conditions as three pairs [1]: two targeted inhibitors, dabrafenib against BRAF and trametinib
against MEK; two non-clinical stresses, CoCl2 mimicking hypoxia and acidic media mimicking acidosis;
and two DNA-damaging agents, cisplatin cross-linking DNA and doxorubicin inhibiting topoisomerases.
The four conditions that carry meaningful interaction are the first two pairs, and the two that do
not are the DNA-damaging pair. The same four were also the ones applied continuously, and the same
two were given for a treatment period followed by recovery (Methods), so class and treatment
schedule divide the six conditions identically and cannot be told apart here. With two conditions in
each class, and a grouping noticed only after the results, neither reading is a finding and neither
is claimed. Separating them needs conditions in which class and schedule vary independently.

Limitation 6 is the honest position on Role A. Its own gate failed, we audited our own power
calculation and found it had been too generous, and we report the worse number.

### What this does not show

Separate from the limitations above, these are claims this work may not make, in any form.

```text
  1  NEVER  unseen-condition generalization
  2  NEVER  cross-cell-line or cross-patient generalization
  3  NEVER  clinical treatment recommendation
  4  NEVER  causal treatment-effect estimation
  5  NEVER  a calibrated probability
  6  NEVER  independent biological replication of Role B
  7  NEVER  uniform benefit across the six conditions
  8  NEVER  a confirmed Role-A result
  9  NEVER  single-cell input equivalence, the model being trained on clone pseudobulk
```

The result is bounded to one cell line, six observed experimental conditions, and one observed
detection proxy. It is not evidence about therapy, and the six conditions include non-clinical
stress contexts that no one would administer to anything.

### Generation 2

What would actually test this, none of which was a gate on the present result:

```text
  independent biological replication in a different lineage-traced system
  transfer to conditions the model has never seen, which requires a condition
    representation the present design does not have
  a dataset-independent nuisance definition, since the current abundance block counts
    cells in this experiment's three specific pretreatment libraries and therefore
    cannot be computed anywhere else
  calibration, frozen and tested separately
  out-of-distribution behaviour
```

The first is the one that matters. Everything here rests on a single lineage-traced system, and one
system is one system however carefully it is evaluated.

---

## Conclusions

In the WM989 lineage-traced melanoma system, a clone's gene expression before treatment carries
information about which of six observed experimental conditions it is still detected after:
information specific to each condition, beyond condition identity and captured pretreatment clone
abundance, under clone-held-out evaluation. Under a test preregistered in full, all of the gain in
ordering a clone's six conditions came from an explicit state-by-condition interaction; an additive
state term contributed nothing.

What this means is that pre-existing state in this system is not only a matter of how detectable a
clone is overall. Part of it bears on which condition the clone is still detected after. A state
contribution shared additively across all conditions cannot change their ordering within a clone, so
an analysis that admits state only that way cannot see this part, as the additive model here could
not; allowing state effects to vary by condition improved the ordering. That does not mean separate
biology behind each condition: one programme whose effect differs between conditions would show the
same pattern.

This opens directions others can take. The interaction terms of the frozen model are a starting
point for finding the pretreatment programmes behind these condition-dependent effects, whether one
or several. Any lineage-traced experiment that splits clones across several perturbations can ask
the same question with the same design: ordering within a clone, abundance held fixed in every
model, a null that refits the whole pipeline, and a protocol fixed before the statistic is computed.
The most informative of these is a test in an independent lineage-traced system, which has not yet
been done. If the result holds there, pretreatment profiling could help choose which perturbations
to test on which subpopulations in experimental work. The outcome measured here is a detection
proxy, not death or clinical response; that is where the result points, not what it shows.

Everything needed to reproduce the result from the public data, and to apply the design to new data,
is released with this work.

---

## List of abbreviations

- **AUROC** — area under the receiver operating characteristic curve
- **B** — the captured-abundance nuisance block
- **C1, C2** — the two endpoint families: post-treatment clone detection, and clone-balanced abundance
- **CI** — confidence interval
- **CoCl2** — cobalt(II) chloride
- **CP10K** — counts per 10,000
- **delta_RANK** — W5 minus W4 in equal-clone-weighted within-clone AUROC
- **delta_TOP1** — W5 minus W4 in the lowest-predicted-score top-choice diagnostic
- **GEO** — Gene Expression Omnibus
- **PCA** — principal component analysis
- **U** — the condition-identity indicators
- **W1, W4, W5** — the three preregistered model specifications
- **X** — the clone expression profile, as 50 train-only principal components

---

## Declarations

### Ethics approval and consent to participate

Not applicable. This study is a computational reanalysis of published, publicly available data from
an immortalised cell line. No human participants, human material or animals were involved.

### Consent for publication

Not applicable.

### Availability of data and materials

The datasets analysed during the current study are all public, but they are not all in one
archive. **No new data were generated for this study.**

Role B, which carries the primary result, is in the Gene Expression Omnibus under accession
GSE279162 [7], generated and deposited by Schaff et al. [1]. Five files of
their own analysis code, read to reconstruct their preprocessing rules rather than guess at them,
are archived at Zenodo [8].

Role A is split between two places. The sequencing data of the two samples used here are in GEO
under accession GSE227151 [6]. The three processed barcode tables the reconstruction
needs -- `filtered10XCells.txt`, `stepThreeStarcodeShavedReads_BC_10X.txt` and
`stepThreeStarcodeShavedReads_BC_gDNA.txt` -- are **not part of that GEO deposit**, whose only
supplementary file is `GSE227151_RAW.tar`. The original authors published them in a shared data
package [9], linked from the key resources table of Jain et al. [10], and they
were taken from there on 21 August 2026. The Rewind authors' two R1 scripts are archived at Zenodo
[11].

A shared folder is not a persistent identifier and may move. The size and SHA-256 of every input
file are recorded in the two Stage-22 manifests inside this study's archive, and the pipeline
refuses any file whose bytes differ, so a reader can check whatever copy they obtain against what
was analysed here. Those files belong to the original authors and are not redistributed here.

All analysis code, frozen protocols, stage records, out-of-fold predictions, the serialized
predictor and the verification tooling are archived at Zenodo [12] and developed openly on
GitHub [13].

#### Verify before reading anything else

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
  python experiments/run_gen1_manuscript.py --verify
```

The first re-hashes every locked artifact and refuses if one has moved. The second does the same
for the claim set. Both were shown to refuse a one-bit change before either was issued. The
third does the same for this manuscript and the package around it. All three also refuse if the
stage that produced the files did not pass.

```text
  evidence lock digest   60602531449079b8f86debea9ffd69753f8bfad033be4c476751d39da08c78aa
  claim lock digest      fc6dc2f221acf1c7661e36071b9e377571c8f903d38398b55bedc164db7efa61
```

#### Licensing

```text
  software and frozen model   PolyForm Noncommercial License 1.0.0
                              SPDX: PolyForm-Noncommercial-1.0.0
  manuscript text and figures CC BY 4.0
  GSE279162, GSE227151        original depositors' terms; NOT relicensed here
```

Academic, educational, nonprofit and personal research use requires **no permission request,
registration, payment or signed agreement**. Any use outside the licence's permitted noncommercial
purposes — including developing or evaluating a commercial product or service, before it has
earned any revenue — requires a separate commercial license (`COMMERCIAL-LICENSING.md`); where that
notice and the licence differ, the licence governs. This is source-available rather than
OSI-approved open source, because commercial use is restricted; it is stated here rather than left
to be discovered.

#### What is in the repository

Benchmark tables, frozen out-of-fold predictions, the serialized model metadata, the prediction API
and CLI, the model card and schema, every stage protocol, every stage record, and every executor
and contract file. Full inventory and per-file hashes: `results/evidence_lock/`.

#### What is not

```text
  stage24_w5_artifact.npz   44 MB, gitignored. A fresh clone does NOT contain it. Its hash
                            is locked and it rebuilds in about half a minute:
                              python experiments/run_stage24_gen1_tool.py --stage 24c
                            The Zenodo archive does include it, so an unpacked archive
                            verifies without a rebuild.
  raw sequencing data       GSE279162, GSE227151. Not vendored. The size and SHA-256 of every
                            input file are recorded, and the pipeline rebuilds from the public
                            files, refusing on any byte that differs.
```

Neither is a gap any longer. The artifact is one command away, and a fresh clone rebuilt the
clone pseudobulk and the model byte for byte from the public files and recovered the observed
statistic exactly (`REPRODUCIBILITY.md` §2.1).

Full reproduction instructions, environment and runtimes: `results/manuscript/REPRODUCIBILITY.md`.

### Competing interests

HA holds the copyright in CellFate-Rx and offers commercial licences for it. HA has received no
income from it to date.

### Funding

This work received no external funding.

### Authors' contributions

HA conceived the study, set its logic and direction, and directed the AI-assisted implementation
and documentation described in Methods. The author read and approved the final manuscript.

### Acknowledgements

We thank Schaff et al. for generating and openly depositing GSE279162, without which this
reanalysis would not be possible.

### Authors' information

Not applicable.

---

## References

1. Schaff DL, White PE, Cote CJ, Watterson GE, Lin KZ, Fasse AJ, et al.
    Pre-existing cell states predict resistance to multiple treatments.
    Cell Genomics. 2026;6(6):101191. doi:10.1016/j.xgen.2026.101191

2. Shaffer SM, Dunagin MC, Torborg SR, Torre EA, Emert B, et al.
    Rare cell variability and drug-induced reprogramming as a mode of cancer drug resistance.
    Nature. 2017;546(7658):431-5. doi:10.1038/nature22794

3. Emert BL, Cote CJ, Torre EA, Dardani IP, Jiang CL, Jain N, et al.
    Variability within rare cell states enables multiple paths toward drug resistance.
    Nat Biotechnol. 2021;39(7):865-76. doi:10.1038/s41587-021-00837-3

4. Goyal Y, Busch GT, Pillai M, Li J, Boe RH, et al.
    Diverse clonal fates emerge upon drug treatment of homogeneous cancer cells.
    Nature. 2023;620(7974):651-9. doi:10.1038/s41586-023-06342-8

5. Kapoor S, Narayanan A.
    Leakage and the reproducibility crisis in machine-learning-based science.
    Patterns. 2023;4(9):100804. doi:10.1016/j.patter.2023.100804

6. Jain N, Goyal Y, Dunagin MC, Cote CJ, Mellis IA, Emert B, et al.
    Retrospective identification of cell-intrinsic factors that mark pluripotency potential
    in rare somatic cells [dataset]. Gene Expression Omnibus, GSE227151. 2023.
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227151  Accessed 21 Aug 2026.

7. Schaff DL, White PE, Cote CJ, Watterson GE, Lin KZ, Fasse AJ, et al.
    Pre-existing cell states predict resistance to multiple treatments [dataset].
    Gene Expression Omnibus, GSE279162. 2024.
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279162  Accessed 12 Jul 2026.

8. Schaff DL (dylanschaff).
    dylanschaff/Schaff_manuscript: Schaff_manuscript_first_submission [software].
    Zenodo. 2024. doi:10.5281/zenodo.13935305
    https://doi.org/10.5281/zenodo.13935305  Accessed 21 Aug 2026.

9. Jain N, Goyal Y, Dunagin MC, Cote CJ, Mellis IA, Emert B, et al.
    Processed barcode data for iPSC Rewind, GSE227151 [data package].
    filtered10XCells.txt, stepThreeStarcodeShavedReads_BC_10X.txt and
    stepThreeStarcodeShavedReads_BC_gDNA.txt for GSE227151. Dropbox. 2024.
    https://www.dropbox.com/sh/ulu6728tcp49dv2/AAAPwLYQiVLloH_JL38lvTj6a?dl=0
    Accessed 21 Aug 2026. Linked from the key resources table of Jain et al.

10. Jain N, Goyal Y, Dunagin MC, Cote CJ, Mellis IA, Emert B, et al.
    Retrospective identification of cell-intrinsic factors that mark pluripotency potential
    in rare somatic cells. Cell Syst. 2024;15(2):109-133.e10.
    doi:10.1016/j.cels.2024.01.001

11. Jain N, et al. (goldengopherforlife).
    arjunrajlaboratory/iPSC_Rewind: Final Release [software].
    Zenodo. 2024. doi:10.5281/zenodo.7707418
    https://doi.org/10.5281/zenodo.7707418  Accessed 21 Aug 2026.

12. Aviv H.
    CellFate-Rx Gen-1: frozen model and reproducibility artifacts, version 1.0.0 [software].
    Zenodo. 2026. doi:10.5281/zenodo.22769563
    https://doi.org/10.5281/zenodo.22769563  Accessed 16 Sep 2026.

13. Aviv H.
    CellFate-Rx: cellular outcomes of perturbations [software repository].
    GitHub. 2026. Archived at doi:10.5281/zenodo.22769563
    https://github.com/hagai-tyty/Cellular-Outcomes-of-Perturbations  Accessed 18 Sep 2026.

---

## Figure legends

Figures are generated from the locked result files by
`python experiments/make_gen1_figures.py`; no number in them is typed by hand.

**Figure 1. Clone-level prospective design and evaluable population.** The design in WM989, and
the funnel from all clones to the evaluable subset.

**Figure 2. Preregistered clone-specific ranking result.** **(A)** Ranking score by model: W1, W4
and W5. **(B)** The observed ΔRANK against 1,000 full-refit permutations.

**Figure 3. Robustness across strata, and the top-choice diagnostic.** **(A)** ΔRANK by held-out
fold. **(B)** ΔRANK by pretreatment clone depth. **(C)** Choosing each clone's lowest predicted
detection score under W5 and W4 — a consistency check, not a significance test.
