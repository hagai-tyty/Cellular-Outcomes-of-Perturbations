# CellFate-Rx Generation 1 — reproducibility package

Everything needed to check this work, in the order a stranger should do it.

---

## 1. Verify the locks first

Before reproducing anything, confirm the artifacts are the ones the manuscript was written against.

```text
  python experiments/run_gen1_evidence_lock.py --verify
  python experiments/run_gen1_claim_lock.py --verify
```

```text
  evidence lock digest   2b29d8c198e56ccb4823d0eed7cd15d5dd67c51694108fb9920a16c2b3e11cae
  claim lock digest      7db3c35948267d32552ad80aab374160ecb620f170faff870a7ee9e4066ca0fa
```

`EVIDENCE_INTACT` and `CLAIMS_INTACT` mean every hashed file is byte-for-byte what was locked.
Anything else means something moved, and the right response is to find out why — not to re-lock.

Text is hashed canonical-LF and binary raw, so the digests are the same on Windows, macOS and
Linux regardless of line-ending settings.

**One step first.** The model artifact is not in the repository (see §5), so rebuild it before
verifying:

```text
  python experiments/run_stage24_gen1_tool.py --stage 24c
```

---

## 2. Environment

```text
  Python        3.11 (developed on 3.11, Windows 10)
  numpy         2.x
  pandas        2.x
  scikit-learn  1.x
  pytest        8.x
```

Set `PYTHONUTF8=1` on Windows if the console codepage is not UTF-8, or the stage scripts will fail
on non-ASCII output rather than on anything meaningful.

No GPU. No network access is required by any stage; the raw data is not downloaded by these
scripts.

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

## 5. What is NOT in this package

```text
  results/stage24/stage24_w5_artifact.npz
      44 MB, gitignored. A fresh clone does NOT contain it. Its hash is locked and the
      rebuild is one command, above.

  raw sequencing data
      GSE279162 (WM989, Role B primary) and GSE227151 (Rewind, Role A supporting) are
      not vendored. Accessions are locked; bytes are not.
```

Naming a gap is not closing it. Both stay open.

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

## 7. Reading the record

Each stage has a record in `plans/(newer)practical plans/RECORDs/` stating what it did, what it
proved, what it did **not** prove, and what went wrong. Errors are recorded where they happened and
are not edited out later; where a conclusion changed, the original stands and the correction is
added beside it.

The scientific result is `stage_25_RECORD.md`. The rest is how it was made trustworthy.
