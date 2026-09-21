# Pre-flight checklist — Generation 1

Split out of `SUBMISSION.md` on 2026-09-19. It lives in its own file because it is the one document
here that is MEANT to change: every tick records progress through the release and the submission.

`SUBMISSION.md` is hashed into the package digest, so while the checklist sat inside it, ticking a
box moved a digest that is supposed to be frozen, and the next stage run or `--verify` would report
`PACKAGE_MOVED`. That happened: the archive published on 2026-09-19 carries a package digest one
cover-letter edit out of date. This file is deliberately NOT in `PACKAGE_FILES`, so ticking a box
costs nothing. It is deliberately NOT added to `RELEASE_DOCUMENTS` in `make_release_bundle.py` either: that
module is itself inside the evidence lock, so editing it to scan this file would move the two
digests printed on the manuscript's title page -- which is the very churn this split exists to
stop. This file is an internal progress tracker, not a document a reader of the archive meets
first, so it stays outside both.

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
[x] archive v1.0.1 published: 10.5281/zenodo.22849702, concept 10.5281/zenodo.22769562, read back
    with all three MD5s equal to this machine's -- superseded: its package digest was stale
[x] archive v1.0.2 published: 10.5281/zenodo.22849793 (2026-09-20), concept unchanged, read back
    with all three MD5s equal to this machine's and PACKAGE_INTACT in the tree it was cut from
[x] GitHub release gen1-v1.0.0 on the archived commit, linking the DOI -- the tag never moves
[x] GitHub release gen1-v1.0.2 on 94ca278, the commit the 1.0.2 bundle was cut from --
    read back: tag on that commit, three assets whose SHA-256 equal dist/ and Zenodo
[x] a Zenodo preprint record: concept 10.5281/zenodo.22829747, newest version v3
    10.5281/zenodo.22850050 (2026-09-20); type Preprint, CC BY, linked to the
    archive record -- bioRxiv declined on 2026-09-18 for want of an organisational affiliation
[x] the preprint DOI written into the cover letter -- by concept DOI only, so republishing the
    preprint no longer requires editing the letter
[x] the preprint DOI added to the archive record's related works: Is supplement to
    10.5281/zenodo.22829747, the preprint's concept DOI, so it follows the newest version
[x] Springer Nature asked, before submitting, whether the waiver can be decided first and whether
    withdrawal is free if it is refused -- enquiry sent 2026-09-21, reply awaited
[ ] BMC Bioinformatics: submit ONLY once that reply is in; waiver checkbox ticked; preprint disclosed
```

**Order matters in Phase 3.** The DOI and every FILL field must be written before the locks are
re-run, or the archived bundle does not contain its own DOI. The locks must be re-run before the
commit, and the bundle built after it: the bundle records the commit it was cut from and refuses a
dirty tree. A Zenodo record's files cannot be changed once it is published, and the release tag is
never moved afterwards -- later corrections are new commits.
