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
  evidence lock digest   6e0c805592d515214fe8795d852b01c7778680762c9deae15d601faa0189e081
  claim lock digest      a5a92ad356a571fb30aa26254745bf70d445c9496a79b569388832aa924d5366
```

`EVIDENCE_INTACT` and `CLAIMS_INTACT` mean two things at once: every hashed file is byte-for-byte
what was locked, **and** the stage that produced it recorded a passing verdict. The two are checked
separately because they came apart once — a refused run still writes its outputs, so its bytes match
what was recorded, and a verifier that only re-hashed them reported clean over a refused package.

A failure therefore names which of the two it is. `*_MOVED` means a hashed file changed;
`*_STAGE_REFUSED` means nothing moved but the stage itself refused, so the recorded state is
not one to build on.
Either way the right response is to find out why — not to re-lock.

Text is hashed canonical-LF and binary raw, so the digests are the same on Windows, macOS and
Linux regardless of line-ending settings.

**One step first.** The model artifact is not in the repository (see §5), so rebuild it before
verifying:

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

## 6.1 Terms

```text
  software + frozen model     PolyForm Noncommercial License 1.0.0  (see LICENSE)
  manuscript + figures        CC BY 4.0
  GSE279162 / GSE227151       original depositors' terms; not relicensed
```

Academic, educational, nonprofit and personal research need no permission request, registration or
agreement. Commercial deployment needs a separate license: `COMMERCIAL-LICENSING.md`.

## 7. Reading the record

Each stage has a record in `plans/(newer)practical plans/RECORDs/` stating what it did, what it
proved, what it did **not** prove, and what went wrong. Errors are recorded where they happened and
are not edited out later; where a conclusion changed, the original stands and the correction is
added beside it.

The scientific result is `stage_25_RECORD.md`. The rest is how it was made trustworthy.
