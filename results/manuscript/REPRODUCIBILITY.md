# CellFate-Rx Generation 1 — reproducibility package

Everything needed to check this work, in the order a stranger should do it.

---

## 1. Verify the locks first

Before reproducing anything, confirm the artifacts are the ones the manuscript was written against.

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
  python experiments/run_gen1_manuscript.py --verify
```

```text
  evidence lock digest   338bff1073130f9a237a409a86920357749bc12b6c6405f03633c57bb71185c5
  claim lock digest      489d0f20f57b86d98e85e80667253d62ba16437629cd01f33f7fec4db2f1aeef
  package digest         results/manuscript/GEN1_PACKAGE_DIGEST.json
```

The third command checks the manuscript and this document itself. Its digest is not printed here
because this file is one of the files it covers — a document cannot quote a digest taken over its
own bytes without changing it. Read it from the file named above.

`EVIDENCE_INTACT`, `CLAIMS_INTACT` and `PACKAGE_INTACT` mean two things at once: every hashed file
is byte-for-byte what was locked, **and** the stage that produced it recorded a passing verdict.
The two are checked separately because they came apart once — a refused run still writes its
outputs, so its bytes match what was recorded, and a verifier that only re-hashed them reported
clean over a refused package.

A failure therefore names which of the two it is. `*_MOVED` means a hashed file changed;
`*_STAGE_REFUSED` means nothing moved but the stage itself refused, so the recorded state is
not one to build on.
Either way the right response is to find out why — not to re-lock.

Text is hashed canonical-LF and binary raw, so the digests are the same on Windows, macOS and
Linux regardless of line-ending settings.

**In a GitHub checkout, one step first.** The model artifact is not in the repository (see §5), so
rebuild it before verifying. The Zenodo archive already contains it; skip this step there (§5.1).

```text
  python experiments/run_stage24_gen1_tool.py --stage 24c
```

---

## 2. Environment

`environment_lock.txt` at the repository root is the authoritative record — a full freeze of the
interpreter that produced the locked results.

```text
  Python 3.11.0    numpy 2.4.6    pandas 3.0.3    scipy 1.17.1    scikit-learn 1.9.0
```

**It supersedes `requirements.txt` for Generation 1.** That file claimed scikit-learn 1.8.0 while
the models were fitted under 1.9.0; a pinned file disagreeing with the machine that produced the
numbers is worse than none. **Honest caveat:** the lock captures the environment as it stands, not
retroactively at each stage's execution. Bit-identical reproduction of the Stage-25 null on a
different stack is not claimed.

Set `PYTHONUTF8=1` on Windows if the console codepage is not UTF-8, or the stage scripts will fail
on non-ASCII output rather than on anything meaningful.

No GPU. No network access is required by any stage; the raw data is not downloaded by these
scripts. §2.1 says where each input file comes from and how to rebuild from it.

## 2.1 Rebuilding from the public data

§3 starts from the clone pseudobulk. This section builds that pseudobulk, and the benchmark tables
before it, from the data the original authors published, so that no analysis input has to be taken
from this repository on trust. The scripts download nothing; the files below are fetched by hand,
once.

**What is needed.** Both datasets, even for the primary result: the frozen Stage-22 builder
reconstructs Role A before Role B, and the Stage-23 input audit checks both. The size and SHA-256 of
every file below are recorded in the two Stage-22 manifests, and Stage 23A refuses unless each one
matches, so a wrong or partial download is caught rather than analysed.

```text
  WM989, GSE279162 -- Role B, the primary result
    GEO     GSE279162_RAW.tar (762,593,280 bytes): for each of the nine samples, its
            filtered barcodes, features and matrix -- 27 files
    GEO     GSE279162_series_matrix.txt.gz, and GSE279162_family.xml from the MINiML download
    Zenodo  Schaff et al.'s analysis code, repository dylanschaff/Schaff_manuscript
            (a Zenodo search finds 10.5281/zenodo.13935305): the five files listed under
            author_code_files in results/stage22_wm989_benchmark_manifest.json

  Rewind, GSE227151 -- Role A, supporting
    GEO     from GSE227151_RAW.tar, only the six files of GSM7092515 and GSM7092516
    GEO     GSE227151-GPL18573_series_matrix.txt.gz, and GSE227151_family.xml
    authors filtered10XCells.txt, stepThreeStarcodeShavedReads_BC_10X.txt and
            stepThreeStarcodeShavedReads_BC_gDNA.txt. These three are NOT in the GEO deposit,
            whose only supplementary file is GSE227151_RAW.tar. They are in the authors' shared
            data package, linked from the key resources table of Jain et al. (Cell Systems 2024,
            doi:10.1016/j.cels.2024.01.001):
              https://www.dropbox.com/sh/ulu6728tcp49dv2/AAAPwLYQiVLloH_JL38lvTj6a?dl=0
            taken from there on 2026-08-21
    Zenodo  the Rewind authors' code, record 7707418 (arjunrajlaboratory/iPSC_Rewind): the two
            R1 scripts listed under author_code_files in
            results/stage22_rewind_benchmark_manifest.json
```

**A second package exists, and the frozen results do not use it.** The outcome materials for
replicates 2 and 3 were never deposited in GEO either; they are in a second shared package,

```text
  https://www.dropbox.com/sh/zz958910t4fkj9w/AAAgTVwO5yAKZ1TpSQVfV6Qga?dl=0   taken 2026-08-24
```

which supplied the R2 and R3 folders used in the re-examination recorded in
`stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md`. No file from it appears in either Stage-22 manifest,
so nothing reported in the manuscript rests on it. A rebuild does not need it.

**Where they go.** The loaders read fixed paths under two folders:

```text
  <wm989-root>/                           <rewind-root>/
    GSE279162_family.xml                    filtered10XCells.txt
    GSE279162_series_matrix.txt.gz          stepThreeStarcodeShavedReads_BC_10X.txt
    GSM8562999_Naive1_filtered_*.gz         stepThreeStarcodeShavedReads_BC_gDNA.txt
      ... all 27 files side by side         GSE227151-GPL18573_series_matrix.txt.gz
    author_code_Schaff_manuscript/          GSE227151_family.xml
      the five files                        GSM7092515/GSM7092515_1_2_control_*.gz
                                            GSM7092516/GSM7092516_1_1_control_*.gz
                                            author_code_zenodo7707418/plotScripts/rewind10X/R1/
                                              the two R1 scripts
```

Both roots default to the folders on the machine that produced the results, `D:\GSE279162` and
`D:\GSE227151_Rewind`; anywhere else, pass them. The layout is not a formality. On that machine the
three Rewind barcode tables were later moved into a subfolder, and the first rebuild attempt stopped
at Stage 22 with `BENCHMARK_BLOCKED_LINKAGE` until it was given a copy with them at the root.

**The commands**, with `R` and `W` the two roots:

```text
  python experiments/build_stage22_prospective_benchmarks.py --rewind-root R --wm989-root W
  python experiments/run_stage23_learnability_gate.py --stage 23a --rewind-root R --wm989-root W
  python experiments/run_stage23_learnability_gate.py --stage 23b --rewind-root R --wm989-root W
  python experiments/run_stage23_learnability_gate.py --stage 23c --rewind-root R --wm989-root W
  python experiments/run_stage23_learnability_gate.py --stage 23d --rewind-root R --wm989-root W
  python experiments/run_stage23_learnability_gate.py --stage 23f --rewind-root R --wm989-root W
  python experiments/run_stage24_gen1_tool.py --stage 24b --wm989-root W
```

Then §3, from `--stage 24c`.

**What it was shown to do.** On 2026-09-15, in a fresh clone of commit `8ccdbe1` with no analysis
cache, and with input files whose sizes and SHA-256 all match the record, this chain and §3 through
`--stage 25a` ran in under ten minutes on the machine that produced the results:

```text
  22    102 s     23a    74 s     23b    57 s     23c    80 s     23d    69 s
  23f     3 s     24b   141 s     24c    29 s     25a    10 s
```

Every committed output came back byte-identical except in three ways, none of them a number the
analysis reports: the Rewind manifest records the folder it was read from, and the Stage-22 summary
records that manifest's size and hash; three result files record their own runtime; and the model's
metadata file was written with Windows line endings, identical once those are normalised. The clone
pseudobulk came back byte-identical, the rebuilt model matches its locked SHA-256, and the observed
statistic is the recorded one, ΔRANK +0.051605. The two permutation nulls, 23E and 25b, were not
re-run: together they take about 16 hours, and they refit on the inputs shown identical here.

The GEO series matrix and family XML are checked for presence, size and SHA-256 only; no stage
reads their contents. GEO regenerates those two files when a record's metadata changes, so a fresh
download may no longer match the recorded bytes. If so, Stage 23A refuses and names them, and the
refusal concerns the record of those two files, not the data the analysis uses.

---

## 3. Reproducing the result

Stages run in order. Each writes its own JSON and refuses if its inputs have moved.

```text
  python experiments/run_stage24_gen1_tool.py --stage 24c     rebuild the model artifact
  python experiments/run_stage25_ranking.py --stage 25a       the observed statistic
  python experiments/run_stage25_ranking.py --stage smoke     sharding equivalence proof
  python experiments/run_stage25_ranking.py --stage 25b       one shard of the null
  python experiments/run_stage25_ranking.py --stage 25c       merge, verdict
  python experiments/run_stage26_scope_lock.py --stage all    the scope lock
  python experiments/run_gen1_evidence_lock.py --stage all    the evidence lock
  python experiments/run_gen1_claim_lock.py --stage all       the claim lock
  python experiments/run_gen1_manuscript.py --stage all       this package
```

### Runtimes, stated honestly

```text
  24c rebuild            ~0.5 min
  25a observed           seconds
  25b the null           10.7 h wall across three shards, measured; 115 s per draw per
                         shard. Single-process it is 66 s per draw, so three shards buy
                         about 1.7x, not 3x -- the rest is memory-bandwidth contention.
  25c merge and verdict  seconds
  26, evidence, claim    seconds each
```

The 10.7 h figure is the real measurement, not an estimate. It came in under the 19-20 h budget
that was accepted before the run started.

The null writes per-shard cache files and asserts completeness before reading a number. A missing
draw is an integrity stop, not a smaller null — an earlier design lost one draw of 300 to a race
between shards appending to a shared file, which is why every later design writes one file per
shard.

---

## 4. Tests

```text
  python -m pytest -q
```

Every stage ships contracts, and the sharpest of them are negative: that an incomplete null is
refused rather than silently shrinking, that a scope hole is not resolved by widening the scope,
that a lock refuses a one-bit change, and that the manuscript checker refuses a manuscript with a
forbidden claim planted in it.

**The suite must leave the working tree unchanged.** Run it and confirm:

```text
  python -m pytest -q
  git status --porcelain          # must print nothing
```

A test that writes a committed artifact is a side effect, not a check — one of them did, and the
close-out pass caught it. CI now snapshots the tree before the suite and fails the build if the
suite moved anything.

---

## 4.1 Per-draw source data

The distributions behind the headline statistic, not just their summaries:

```text
  results/stage25/stage25_null_draws.csv             1,000 full-refit permutation draws
  results/stage25/stage25_bootstrap_replicates.csv   2,000 clone-bootstrap replicates
  results/manuscript/figures/figure_source_data.json every value each panel draws
```

Regenerate with `python experiments/export_gen1_source_data.py`, which refuses to write any file
that does not reproduce the recorded verdict statistic exactly. The null draws previously existed
only in a gitignored shard cache — 10.7 h of compute on a single machine, summarised to six numbers.

## 5. What a GitHub checkout lacks, and what the Zenodo archive adds

```text
  results/stage24/stage24_w5_artifact.npz
      44 MB, gitignored. A fresh clone does NOT contain it. Its hash is locked and the
      rebuild is one command, above. The Zenodo archive includes it.

  _cc_cache/stage23/GSE279162_pseudobulk.npz
      the clone pseudobulk the permutation null refits on. Gitignored, so a fresh clone
      lacks it too; §2.1 rebuilds it from the public data. The Zenodo archive includes it.

  raw sequencing data
      GSE279162 (WM989, Role B primary) and GSE227151 (Rewind, Role A supporting) are in
      neither. Accessions are locked; bytes are not vendored.
```

Naming a gap is not closing it, so each one has a way through. In a checkout the model is one
command away, and the pseudobulk is the §2.1 rebuild away. In the archive only the raw sequencing
data remains outside, and §2.1 names where each file comes from and checks every byte of it.

## 5.1 Verifying the Zenodo archive

Check the download before unpacking it: its SHA-256 must equal `bundle_sha256` in
`BUNDLE_CONTENTS.json`, published beside it. Then, from inside the unpacked `cellfate-rx-gen1/`:

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
  python experiments/run_gen1_manuscript.py --verify
  PYTHONPATH=src python -m cellfate.gen1_cli \
      --artifact results/stage24/stage24_w5_artifact.npz \
      --meta results/stage24/stage24_w5_artifact.json \
      --expression results/stage24/tool/example_clone_expression.npy \
      --nuisance results/stage24/tool/example_clone_nuisance.txt
```

`PYTHONPATH=src` is not decoration. If CellFate-Rx is also installed from a checkout, a plain
`python -m cellfate.gen1_cli` run inside the archive imports that installed copy instead of the
archive's own code, and exits cleanly while testing the wrong thing -- this was observed, not
supposed. In Windows PowerShell set it first with `$env:PYTHONPATH = "src"`. The predictor should exit
`0` and print a score for each of the six conditions; its ordering is reported as not validated unless
the preregistered verdict file is supplied, which is the tool working as specified.

---

## 6. Using the tool

```text
  python -m cellfate.gen1_cli \
      --artifact results/stage24/stage24_w5_artifact.npz \
      --meta results/stage24/stage24_w5_artifact.json \
      --expression results/stage24/tool/example_clone_expression.npy \
      --nuisance results/stage24/tool/example_clone_nuisance.txt
```

A real clone from the benchmark ships with it, so the tool can be run without any setup beyond the
rebuild. Exit codes distinguish outcomes: `0` every condition scored, `2` at least one refused,
`3` input unreadable. A refusal that exited `0` would let a caller treat a missing score as a real
one.

Inputs, outputs and refusal semantics: `results/stage24/tool/io_schema.json` and
`results/stage24/tool/MODEL_CARD.md`.

---

## 6.1 Terms

```text
  software + frozen model     PolyForm Noncommercial License 1.0.0  (see LICENSE)
  manuscript + figures        CC BY 4.0
  GSE279162 / GSE227151       original depositors' terms; not relicensed
```

Academic, educational, nonprofit and personal research need no permission request, registration or
agreement. Any commercial use needs a separate license — including developing or evaluating a
commercial product or service before it has earned any revenue: `COMMERCIAL-LICENSING.md`. Where
that notice and the licence differ, the licence governs.

## 7. Reading the record

Each stage has a record in `plans/(newer)practical plans/RECORDs/` stating what it did, what it
proved, what it did **not** prove, and what went wrong. Errors are recorded where they happened and
are not edited out later; where a conclusion changed, the original stands and the correction is
added beside it.

The scientific result is `stage_25_RECORD.md`. The rest is how it was made trustworthy.
