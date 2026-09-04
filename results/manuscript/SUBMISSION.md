# Submission pack — Generation 1

Everything a submission form asks for, in one place, plus the full related-work section. Nothing
here is uploaded by this repository; this is the copy-paste source.

**Fields marked `<<FILL>>` need a human.** They are things I cannot know or must not invent.

---

## 1. Identity

```text
TITLE
  Pretreatment transcriptional state carries condition-specific information about
  future clonal detection in a lineage-traced melanoma line

SHORT TITLE
  Condition-specific state information in lineage-traced melanoma clones

AUTHORS            <<FILL: full name, ORCID>>
AFFILIATION        <<FILL>>
CORRESPONDING      <<FILL: email>>

ARTICLE TYPE       Research Article
PREPRINT SERVER    bioRxiv
bioRxiv CATEGORY   Bioinformatics   (secondary: Cancer Biology)
TARGET JOURNAL     BMC Bioinformatics
BMC SECTION        Methods / benchmarking and reanalysis
LICENCE            manuscript + figures  CC BY 4.0
                   software + model      PolyForm-Noncommercial-1.0.0
```

### Honest note on venue fit

BMC Bioinformatics publishes reanalyses and benchmark/methodology papers, which this is. Two things
a handling editor will notice immediately, and neither should be hidden in a cover letter:

- **No new data.** The primary dataset is Schaff et al.'s. The contribution is the evaluation
  design, the preregistered test, and the reproducibility apparatus.
- **One biological system.** No replication. The manuscript says so in the abstract.

If the editor considers the increment too small for a Research Article, the natural fallbacks are a
BMC Bioinformatics *Software* article (the frozen predictor plus its refusal semantics), or *GigaScience*
/ *Bioinformatics Advances*, both of which take reanalysis-plus-resource submissions. Deciding that
is not a Generation-1 gate.

---

## 2. Abstract as submitted

> Whether a cell's molecular state before a perturbation predicts what happens to it afterwards is
> usually asked retrospectively, after outcome and state have been measured in the same cells. We
> ask it prospectively at clone level, in one BRAF-V600E melanoma cell line (WM989, GSE279162), in
> which a barcoded population was split across six observed experimental conditions: Acid,
> Cisplatin, CoCl2, Dabrafenib, Doxorubicin and Trametinib. Of the clones that experiment
> recovered, 1,401 carry a pretreatment profile and are therefore analysable prospectively; those
> are the clones used here.
>
> Within this system, pretreatment gene expression contains condition-specific information about
> future clonal detection beyond condition identity and captured pretreatment clone abundance,
> under clone-held-out evaluation frozen before any result existed. Under a test preregistered in
> full — the metric, population, weighting, comparator, null and verdict rule all fixed before the
> numbers existed — a frozen state-by-condition interaction model improves clone-specific ordering
> of the six conditions over a non-interactive additive model: +0.051605 in equal-clone-weighted
> within-clone AUROC, 95% CI [+0.037197, +0.065571], with 0 of 1000 full-refit permutation draws
> reaching the observed value (p < 0.001).
>
> The outcome is an observed post-treatment clone-detection proxy and is **not death**,
> sensitivity, resistance or clinical response. The six conditions are the entire supported
> vocabulary. We make no claim about unseen conditions, other cell lines, or patients, and the
> model emits no calibrated probability. Independent biological replication has not been performed
> and is Generation 2 work, not a gate on this result.
>
> ---

**KEYWORDS**
clonal barcoding; lineage tracing; drug resistance; preregistration; reproducibility;
within-clone ranking; permutation test; melanoma; reanalysis

---

## 3. Related work

Four strands, and the paper sits at the intersection of the first and the last.

### 3.1 Non-genetic single-cell state predicts which cells resist

The founding result in this exact system: Shaffer et al. showed that rare, transient transcriptional
states in WM989 predict which cells resist vemurafenib, and that drug exposure epigenetically
stabilises those states into durable resistance [3]. Emert et al. developed Rewind to capture the
rare precursors directly, resolving substructure within them that predicts distinct downstream
resistant behaviours [4]. Goyal et al., with FateMap, showed that resistant fates are diverse rather
than binary and are largely **predetermined** by pre-treatment molecular differences rather than by
extrinsic factors [5]. Schaff et al. extended clonal tracing to six conditions in parallel and
reported cross-condition resistance correlation plus CD44 as a marker of resistance across several
conditions [1] — the dataset reanalysed here.

**We do not claim priority on condition-specific expression analysis.** Schaff et al.'s deposited
analysis reports condition-associated markers and signatures, and earlier work resolves substructure
distinguishing resistant outcomes. What that literature does not contain is a frozen, clone-held-out,
preregistered test of clone-specific *ordering* across conditions with abundance held fixed.

### 3.2 The gap this addresses

A clone-level propensity cannot, by construction, order conditions *within* a clone: any quantity
acting on the clone as a whole shifts all of its predicted scores together and leaves their order
unchanged. So "state predicts resistance" and "state predicts which condition" are separate claims
requiring separate tests, and the second needs an explicit state-by-condition interaction. This work
tests the second, with the first entering the model as an additive term so that it cannot supply the
answer. Empirically it does not: the additive model scores *below* the no-state baseline on ordering.

### 3.3 Confounding by capture depth

Detection-based clone outcomes are dominated by how many cells a clone contributed before treatment.
Any comparison that does not hold that fixed measures a headcount. Here abundance is a mandatory
model term, present in every model including the permutation null, and the shipped tool refuses to
score without it. It remains roughly 3.45× the whole state contribution.

### 3.4 Leakage and preregistration in ML-based science

Kapoor and Narayanan document eight kinds of leakage across 294 papers in seventeen fields and show
that complex models frequently fail to beat logistic regression once leakage is corrected [6]. The
design here is a direct response: the comparator is a simpler model of the same family; gene
filtering, PCA and scalers are refitted inside each training fold; hyperparameters are selected in
an inner split of the training folds only; the null refits the entire pipeline rather than shuffling
labels; and the metric, population, weighting, comparator, null and verdict rule were fixed in a
digest-frozen protocol before any number existed. Preregistration of a computational analysis is
still uncommon, and the protocol digests make the claim checkable rather than assertable.

### 3.5 References

```text
[1] Schaff DL, White PE, Cote CJ, Watterson GE, Lin KZ, Fasse AJ, Zhang NR, Shaffer SM.
    Pre-existing cell states predict resistance to multiple treatments.
    Cell Genomics 6(6):101191, 2026. doi:10.1016/j.xgen.2026.101191  PMID 41916275

[2] GEO GSE227151. Retrospective identification of cell-intrinsic factors that mark
    pluripotency potential in rare somatic cells (scRNA-seq). Human hiF-T fibroblasts.

[3] Shaffer SM, Dunagin MC, Torborg SR, Torre EA, Emert B, et al.
    Rare cell variability and drug-induced reprogramming as a mode of cancer drug resistance.
    Nature 546(7658):431-435, 2017. doi:10.1038/nature22794  PMID 28607484

[4] Emert BL, Cote CJ, Torre EA, Dardani IP, Jiang CL, Jain N, Shaffer SM, Raj A.
    Variability within rare cell states enables multiple paths toward drug resistance.
    Nature Biotechnology 39(7):865-876, 2021. doi:10.1038/s41587-021-00837-3  PMID 33619394

[5] Goyal Y, Busch GT, Pillai M, Li J, Boe RH, et al.
    Diverse clonal fates emerge upon drug treatment of homogeneous cancer cells.
    Nature 620(7974):651-659, 2023. doi:10.1038/s41586-023-06342-8  PMID 37468627

[6] Kapoor S, Narayanan A.
    Leakage and the reproducibility crisis in machine-learning-based science.
    Patterns 4(9):100804, 2023. doi:10.1016/j.patter.2023.100804  PMID 37720327
```

---

## 4. Declarations

### Availability of data and materials

> The dataset analysed here is publicly available from the Gene Expression Omnibus under accession
> GSE279162, generated and deposited by Schaff et al. [1]. Supporting Role-A evidence uses
> GSE227151. **No new data were generated for this study.**
>
> All analysis code, frozen protocols, stage records, out-of-fold predictions, the serialized
> predictor and the verification tooling are archived at Zenodo, DOI `<<FILL: 10.5281/zenodo.XXXXXXX>>`,
> and developed openly at `https://github.com/hagai-tyty/Cellular-Outcomes-of-Perturbations`.
> The archive contains a manifest of SHA-256 digests and three verification commands; the evidence,
> claim and package digests recorded in the manuscript can be re-derived from it.

### Software availability and licensing

> CellFate-Rx is source-available under the PolyForm Noncommercial License 1.0.0 (SPDX:
> `PolyForm-Noncommercial-1.0.0`). Academic, educational, nonprofit and other noncommercial use
> requires no permission request, registration, payment or signed agreement. Commercial use as part
> of a revenue-generating product or paid service requires a separate license; see
> `COMMERCIAL-LICENSING.md`. Manuscript text and figures are licensed CC BY 4.0.

**Note for the editor.** BMC Bioinformatics' software policy requires that software be freely
available to non-commercial researchers, which this license satisfies without any gate. It is,
however, **source-available rather than OSI-approved open source**, since commercial use is
restricted. State this plainly in the cover letter rather than letting it surface in review.

### Competing interests

> `<<FILL — if none: "The author declares no competing interests.">>`

### Funding

> `<<FILL — if none: "This research received no specific grant from any funding agency.">>`

### Ethics approval and consent to participate

> Not applicable. This study is a computational reanalysis of published, publicly available data
> from an immortalised cell line. No human participants, human material or animals were involved.

### Consent for publication

> Not applicable.

### Authors' contributions

> `<<FILL>>` — single-author template: "H.A. designed the evaluation, implemented the analysis and
> verification tooling, and wrote the manuscript."

### Acknowledgements

> We thank Schaff et al. for generating and openly depositing GSE279162, without which this
> reanalysis would not be possible. `<<FILL: any AI-assistance disclosure your venue requires>>`

**Note on AI assistance.** BMC and bioRxiv both require disclosure where generative AI contributed
to the work. This project was developed with substantial AI assistance for implementation, checking
and drafting. Say so plainly in the Acknowledgements; do not list a model as an author, which every
major publisher forbids.

---

## 5. Figures

```text
Figure 1  Clone-level prospective design and evaluable population
Figure 2  Preregistered clone-specific ranking result (models; observed vs null)
Figure 3  Robustness across strata, and the top-choice diagnostic
```

Vector SVG at `results/manuscript/figures/`. Regenerate with
`python experiments/make_gen1_figures.py`; every number is read from a locked result file, none is
typed into the script. Convert to PDF/EPS/TIFF at submission time if the venue demands it.

---

## 5.1 DOME self-assessment

DOME (Data, Optimization, Model, Evaluation) is the community reporting standard for supervised ML
in biology [7,8]. **Caveat on this assessment:** the verbatim questionnaire sits behind a paywall
and in the DOME registry wizard, so this is scored against the four categories and the sub-areas the
open sources name, not against a numbered item list. Anyone submitting should re-score in the
registry itself.

```text
DATA
  provenance          GSE279162, generated by Schaff et al., cited; NO new data generated
  splits              5 outer folds, held out BY CLONE, fixed before any model was fitted
  test independence   no component scores a clone it trained on; fold isolation verified
                      per component, and every training clone set recorded
  preprocessing       gene filter, PCA basis and all scalers refitted INSIDE each training
                      fold; hyperparameters chosen in an inner split of training folds only
  population          892 of 1,401 evaluable, with both exclusions counted and reported
  GAP                 raw sequencing not redistributed; accessions given, derived clone
                      pseudobulk included in the archive

OPTIMIZATION
  search              inner GroupKFold over the frozen grid; selection rule fixed in Stage 23
  seeds               recorded and reported (bootstrap 23501, permutation 23523)
  reuse               observed-data hyperparameters NEVER reused inside a null draw
  GAP                 the environment lock is post-hoc, not captured per stage; bit-identical
                      reproduction on a different stack is not claimed

MODEL
  specification       W5 = X + B + U + X*U, logistic, 309 design columns, fully enumerated
  availability        serialized and shipped; regenerates every frozen prediction to 6.7e-16
  interpretability    linear in 50 PCs, 4 nuisance terms, 5 indicators and 250 interactions;
                      coefficients are in the artifact
  refusal semantics   unknown condition, missing nuisance and wrong schema each documented,
                      each tested adversarially

EVALUATION
  metric              equal-clone-weighted within-clone AUROC, preregistered
  comparator          W4, preregistered, chosen so an additive term cannot supply the answer
  baseline            W1 reported alongside
  uncertainty         2,000-replicate clone bootstrap; stated as CONDITIONAL on the fitted
                      models, which the permutation null is not
  significance        1,000-draw full-refit permutation; p reported as a floor, never a point
                      estimate; all 1,000 per-draw values published
  calibration         none performed, and calibrated-probability claims are forbidden
  GAP                 NO external or independent validation. One cell line. Stated in the
                      abstract, not buried in the discussion.
```

Two DOME gaps are real and neither is hidden: no independent validation, and a post-hoc
environment lock. Both are in the manuscript.

```text
[7] Walsh I, Fishman D, Garcia-Gasulla D, Titma T, Pollastri G, et al.
    DOME: recommendations for supervised machine learning validation in biology.
    Nature Methods 18:1122-1127, 2021.  doi:10.1038/s41592-021-01205-4

[8] Ghiandoni GM, et al. DOME Registry: implementing community-wide recommendations for
    reporting supervised machine learning in biology.
    GigaScience, 2024.  doi:10.1093/gigascience/giae094   PMID 39661723
```

---

## 6. Cover-letter skeleton

> We submit *Pretreatment transcriptional state carries condition-specific information about future
> clonal detection in a lineage-traced melanoma line* for consideration as a Research Article.
>
> Prior work in this system has established that pre-existing single-cell state predicts *whether* a
> clone resists treatment. We ask the adjacent question of *which* condition a clone is still
> detected after — a clone-specific ordering claim that a general resistance propensity cannot, by
> construction, satisfy. Using the publicly deposited six-condition clonal-tracing dataset of Schaff
> et al. (GSE279162), and a protocol frozen by cryptographic digest before any result existed, we
> find that an explicit state-by-condition interaction improves within-clone ordering over a
> non-interactive additive model, exceeding all 1,000 full-refit permutation draws.
>
> The work generates no new data and makes a deliberately bounded claim: one cell line, six observed
> conditions, an observed detection proxy, and no independent biological replication. Those limits
> are stated in the abstract, not only the discussion.
>
> The complete analysis, frozen protocols, stage-by-stage records including negative and failed
> results, and a verification tool that refuses on any modified artifact are archived at Zenodo
> `<<FILL DOI>>`.

---

## 7. Pre-flight checklist

```text
[ ] licence coherent: LICENSE, pyproject.toml, CITATION.cff and .zenodo.json all say
    PolyForm-Noncommercial-1.0.0, and COMMERCIAL-LICENSING.md is present
[ ] GitHub CI green on the commit being archived
[ ] fresh clone verifies: evidence, claim and package digests
[ ] Zenodo DOI reserved
[ ] DOI written into MANUSCRIPT.md, README.md, CITATION.cff, this file
[ ] export_gen1_source_data.py run BEFORE the locks -- it writes four LOCKED files
    (the two per-draw CSVs, figure_source_data.json, environment_lock.txt), so running
    it afterwards invalidates the very digests just computed
[ ] locks re-run in order: evidence -> claim -> manuscript, digests propagated
[ ] bundle rebuilt, then --check: BUNDLE_INTACT (hashes bytes inside the zip, pins the
    zip against its recorded sha256, and refuses if any locked artifact is absent)
[ ] Zenodo record published (immutable)
[ ] GitHub release tagged, release notes link the Zenodo DOI
[ ] bioRxiv submission (PDF + figures + declarations above)
[ ] BMC Bioinformatics submission
```

**Order matters at steps 3–6.** Reserving the DOI first, then writing it into the documents, then
re-running the locks, is the only sequence in which the archived bundle contains its own DOI *and*
its digests are correct. Writing the DOI after locking silently invalidates all three digests.
