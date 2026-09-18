# Submission pack — Generation 1

Everything a submission form asks for, in one place, plus the full related-work section. Nothing
here is uploaded by this repository; this is the copy-paste source.

**Fields marked FILL, in double angle brackets, need a human** — things that cannot be known here
or must not be invented. The release bundle refuses to build while any remains. Fields marked
LATER can only be filled after archiving, such as the preprint DOI, and do not block it.

---

## 1. Identity

```text
TITLE
  Pretreatment transcriptional state carries condition-specific information about
  future clonal detection in a lineage-traced melanoma line

SHORT TITLE
  Condition-specific state information in lineage-traced melanoma clones

AUTHORS            as on the MANUSCRIPT.md title page -- one source, not a copy
AFFILIATION        as on the MANUSCRIPT.md title page
CORRESPONDING      as on the MANUSCRIPT.md title page

ARTICLE TYPE       Research Article
PREPRINT           Zenodo, a repository rather than a preprint server: a record of its own,
                   separate from the archive record
                   bioRxiv declined the submission on 2026-09-18: it requires an organisational
                   affiliation that can adjudicate ethical disputes, which an independent
                   researcher does not have
RELEASE            Zenodo archive, then a Zenodo preprint, then BMC Bioinformatics with the APC
                   waiver requested at submission (2026-09-18; the 2026-09-12 plan was Zenodo and
                   bioRxiv, with the journal deferred)
FORMATTED FOR      BMC Bioinformatics Research article, ready if a journal is chosen later
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

### Fee, if a journal is chosen

BMC Bioinformatics charges an article-processing charge of £2,290 / $3,090 / €2,590, plus VAT or
local taxes, with country-tiered pricing. There is no submission fee: the charge falls due only after
editorial acceptance. Discretionary waivers are considered case by case and can be requested **only
at submission**, never later. Checked on 2026-09-12 against the journal's own "How to publish with
us" page and BMC's article-processing-charge page.

---

## 2. Abstract as submitted

> **Background.**
> Whether a cell's molecular state before a perturbation predicts what happens to it afterwards is
> usually asked retrospectively, after outcome and state have been measured in the same cells. We ask
> it prospectively at clone level, in one BRAF-V600E melanoma cell line (WM989, GSE279162), in which
> a barcoded population was split across six observed experimental conditions: Acid, Cisplatin, CoCl2,
> Dabrafenib, Doxorubicin and Trametinib. Of the clones that experiment recovered, 1,401 carry a
> pretreatment profile and are therefore analysable prospectively; those are the clones used here.
> The ranking test uses the 892 of them detected under at least one condition and undetected under at
> least one other.
>
> **Results.**
> Within this system, pretreatment gene expression contains condition-specific information about
> future clonal detection beyond condition identity and captured pretreatment clone abundance, under
> clone-held-out evaluation whose folds, features and exclusions were fixed
> before any model was fitted. Under a test preregistered in full — the metric, population, weighting,
> comparator, null and verdict rule all fixed by digest before any ranking statistic was computed,
> though after earlier predictive analyses of the same data — a frozen state-by-condition interaction
> model improves clone-specific ordering of the six conditions over a non-interactive additive model:
> +0.051605 in equal-clone-weighted within-clone AUROC, 95% CI [+0.037197, +0.065571], with 0 of 1000
> full-refit permutation draws reaching the observed value (p < 0.001). The additive model did
> not itself improve that ordering over condition identity alone (0.692176 against 0.692654):
> the gain is the interaction.
>
> **Conclusions.**
> A state contribution shared additively across all conditions cannot change their ordering within a
> clone. Allowing state effects to vary by condition improved that ordering in WM989: part of what a
> clone's pretreatment state carries bears on which of the six conditions it is still detected after,
> not only on how detectable it is overall. The outcome is a detection proxy and is not death,
> sensitivity, resistance or clinical response. The design and the artifacts behind it are released,
> so that the result can be reproduced from the public data and tested in other lineage-traced
> systems, where it has not yet been tested.
>
> ---

**KEYWORDS** — one source: the `## Keywords` section of `MANUSCRIPT.md`.

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

A state contribution shared additively across all conditions cannot, by construction, order
conditions *within* a clone: it shifts all of that clone's predicted scores together and leaves
their order unchanged. A single programme acting unequally across conditions is not excluded; it
would appear as an interaction. So "state predicts resistance" and "state predicts which condition"
are separate claims requiring separate tests, and the second needs an explicit state-by-condition
interaction. This work tests the second, with the shared-across-conditions part of state entering
the model as an additive term so that it cannot supply the answer. Empirically it does not: the
additive model scores *below* the no-state baseline on ordering.

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
digest-frozen protocol before any ranking statistic was computed, and after earlier predictive
analyses of the same data. Preregistration of a computational analysis is still uncommon, and the
protocol digests make the claim checkable rather than assertable.

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

The declarations live in one place: the `## Declarations` section of `MANUSCRIPT.md`, under the eight
headings BMC requires (Amendment V1.2). They are not repeated here, so there is no second copy to fall
out of date.

Form-entry notes:

- **BMC** asks during submission for the organisations that funded the work, with any grant
  numbers. Copy them from the manuscript's Funding declaration. Zenodo has a funding field too, and
  it stays empty for the same reason.
- **AI use** is documented in the manuscript's Methods, under *AI use disclosure* — where Springer
  Nature's editorial policy asks for it: in the Methods, and in another part only if there is no
  Methods section. No AI tool is listed as an author.
- **Competing interests.** Do not default to "none". The author holds the copyright in
  CellFate-Rx and offers commercial licences for it, which is the kind of financial interest BMC
  asks authors to declare. Disclose it, and say whether any income has been received from it.
- **Software licence.** BMC Bioinformatics requires software to be freely available for non-commercial
  use, without restrictions such as a material transfer agreement, and recommends but does not require
  an open-source licence. PolyForm Noncommercial 1.0.0 meets the requirement. It is source-available
  rather than OSI open source, and the cover letter says so.

---

## 5. Figures

```text
Figure 1  Clone-level prospective design and evaluable population
Figure 2  Preregistered clone-specific ranking result (models; observed vs null)
Figure 3  Robustness across strata, and the top-choice diagnostic
```

Vector SVG at `results/manuscript/figures/`, regenerated with `python experiments/make_gen1_figures.py`;
every number is read from a locked result file and none is typed into the script. The legends are in
the manuscript's `## Figure legends` section.

For submission, `python experiments/render_gen1_submission.py` builds two manuscript files from the
one Markdown source. `MANUSCRIPT_bioRxiv.pdf` has the three figures embedded above their legends,
because a preprint is distributed as one PDF, and the servers that take one do not accept SVG.
Zenodo converts nothing: it stores the file as uploaded. The name is kept
that name although bioRxiv is no longer the venue: the same file is archived at Zenodo and attached
to the GitHub release, and renaming it here would make those copies disagree with this one.
`MANUSCRIPT_BMC.docx` keeps the legends in the text and ships the figures as separate PDF files,
which Springer Nature accepts. The PDF is exported through Word, since pandoc's own PDF route needs
LaTeX. Every page is inspected before it is approved.

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

## 6. Cover letter

For a journal submission, if one is made. Paste it unchanged apart from the marked fields.

> Dear Editor,
>
> We submit *Pretreatment transcriptional state carries condition-specific information about future
> clonal detection in a lineage-traced melanoma line* for consideration as a Research Article.
>
> Prior work in this system has established that pre-existing single-cell state predicts *whether* a
> clone resists treatment. We ask the adjacent question of *which* condition a clone is still
> detected after — an ordering within a clone that a state contribution shared additively across
> conditions cannot, by construction, produce. Using the publicly deposited six-condition
> clonal-tracing dataset of Schaff et al. (GSE279162), and a ranking protocol frozen by
> cryptographic digest before any ranking statistic was computed, we find that an explicit
> state-by-condition interaction improves within-clone ordering over a non-interactive additive
> model, exceeding all 1,000 full-refit permutation draws.
>
> The work generates no new data and makes a deliberately bounded claim: one cell line, six observed
> conditions, an observed detection proxy, and no independent biological replication. Those limits
> are stated in the abstract, not only the discussion.
>
> The complete analysis, frozen protocols, stage-by-stage records including negative and failed
> results, and a verification tool that refuses on any modified artifact are archived at Zenodo
> under DOI `10.5281/zenodo.22769563`.
>
> **Preprint.** This manuscript is posted as a preprint at Zenodo, DOI
> `10.5281/zenodo.22829747`, under a CC BY licence. That DOI resolves to the newest version
> of the record; the version submitted here is `10.5281/zenodo.22829748`.
>
> **Software licence.** The software and frozen model are released under the PolyForm Noncommercial
> License 1.0.0. They are free for any non-commercial use without registration or agreement, which
> meets the journal's software-availability policy, but the licence is source-available rather
> than OSI-approved open source, because commercial use is restricted. We say so here rather than
> leave it to be discovered in review.
>
> **Use of AI.** Generative AI assistance is disclosed in the Methods section.
>
> `<<LATER: request an APC waiver here if one is needed -- it can only be requested at submission>>`
>
> Sincerely,
> Hagai Aviv

---

## 7. Pre-flight checklist

```text
PHASE 1 -- structure, in the repository
[x] manuscript in BMC Research-article structure (Amendment V1.2); structure checks added, each
    shown to refuse a broken copy
[x] declarations single-sourced in MANUSCRIPT.md; the submitted abstract checked against it
[x] the release bundle refuses any FILL marker; LATER markers are allowed and listed

PHASE 2A -- review corrections, in the repository (Amendment V1.3)
[x] freeze wording precise: folds and features fixed before any model was fitted; the ranking test
    frozen before any ranking statistic, after earlier predictive analyses; checked
[x] REPRODUCIBILITY.md separates a checkout from the Zenodo archive, and verifies the archive with
    its own code (PYTHONPATH=src)
[x] render script for the preprint PDF, the BMC .docx and the figure files; binary outputs protected
    from line-ending conversion; draft renders kept out of git and out of the archive

PHASE 2B -- human inputs and decisions
[x] the manuscript read and approved, including the new prose
[x] author block, declarations and the AI-use description filled -- competing interests disclose the
    copyright and the commercial-licence offer rather than "none"
[x] licensing: A (wording matches the licence) or B (an explicit additional permission, reviewed)
[x] installs approved: pandoc; svglib, reportlab, cffconvert
[x] Zenodo: GitHub integration OFF for this repository; DOI reserved on a saved draft, draft kept

PHASE 3 -- lock, render, build
[x] the licensing decision applied wherever the commercial boundary is described
[x] every FILL marker filled; the DOI written into MANUSCRIPT.md, README.md, CITATION.cff, this file
[x] python experiments/export_gen1_source_data.py    -- the numbers still reproduce; its three
    exports were run in-process instead, leaving the environment lock unrewritten (record, Phase 3i)
[x] python experiments/cascade_gen1.py               -- locks in order, digests re-pinned
[x] the three --verify commands; the full test suite read by pytest's own exit code
[x] CITATION.cff validated
[x] python experiments/render_gen1_submission.py; every page of the PDF inspected: three figures,
    their legends, the DOI, no placeholder
[x] commit, push, green CI on that exact commit
[x] python experiments/make_release_bundle.py, then --check
[x] the unpacked archive passes the three --verify commands and the PYTHONPATH=src predictor

PHASE 4 -- publish, in this order
[x] Zenodo record published; its DOI resolves at doi.org, and the three files' MD5s matched this
    machine before publishing
[x] the published archive record read back: Zenodo's own MD5 for each of the three files
    equals this machine's, so the stored files are the verified ones
[x] GitHub release gen1-v1.0.0 on the archived commit, linking the DOI -- the tag never moves
[x] a Zenodo preprint record: 10.5281/zenodo.22829748, type Preprint, CC BY, linked to the
    archive record -- bioRxiv declined on 2026-09-18 for want of an organisational affiliation
[x] the preprint DOI written into the cover letter
[ ] the preprint DOI added to the archive record's related works: Is supplement to
[ ] BMC Bioinformatics: the APC waiver requested at submission, and the preprint DOI disclosed
```

**Order matters in Phase 3.** The DOI and every FILL field must be written before the locks are
re-run, or the archived bundle does not contain its own DOI. The locks must be re-run before the
commit, and the bundle built after it: the bundle records the commit it was cut from and refuses a
dirty tree. A Zenodo record's files cannot be changed once it is published, and the release tag is
never moved afterwards -- later corrections are new commits.
