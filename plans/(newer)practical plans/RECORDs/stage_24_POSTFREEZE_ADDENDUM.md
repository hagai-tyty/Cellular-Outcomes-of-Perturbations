# stage_24_POSTFREEZE_ADDENDUM — two tool changes made after 24F froze the deliverables

`stage_24_RECORD.md` is unchanged and stays the record of Stage 24 as executed. This addendum
covers two commits made to the shipped tool *after* the 24F freeze, so that nothing in the tool
exists without a record behind it.

## Why an addendum rather than an edit
24F hashed the §6.5 deliverables. Both commits below changed hashed files, so both re-hashed and
re-recorded. Editing `stage_24_RECORD.md` would have hidden that the freeze moved; this file makes
the movement explicit.

---

## Commit `a4e0a7a` — a runnable example clone

```text
  results/stage24/tool/example_clone_expression.npy   36,601 float64, one real benchmark clone
  results/stage24/tool/example_clone_nuisance.txt     its four-value nuisance block
  results/stage24/tool/example_clone_README.md        the exact CLI invocation
```

`example_clones.csv` already existed and was hashed by 24F, but it is a *table*, not something the
CLI can be pointed at. Without a `.npy` in the shipped shape, nobody could actually run the tool
they had been handed. No hashed 24F file was modified by this commit.

## Commit `a17f1c6` — §6.2 form A, `clone_input_from_cells`

The plan's §6.2 defines two input forms: form A (raw pretreatment cells) and form B (an
already-aggregated clone vector plus the nuisance block). Stage 24 shipped only form B, which
forced the caller to hand-compute `B` — four `log1p` cell counts the tool is perfectly capable of
counting itself. Form A now does it.

```python
clone_input_from_cells(counts, samples) -> (X, B)
```

Verified against the raw GEO matrices: the `(X, B)` this returns for a benchmark clone matches the
frozen Stage-22 row exactly, so form A and form B reach the same model input.

**This did not widen scope, and the code says so in three places.** `B` counts cells per WM989
naive library. Those three libraries are the structure of one experiment. Sample labels that merely
*look* like `Naive1/2/3` — another lab, another depth, another library design — produce a `B` the
model never saw. Unknown labels raise rather than being counted into a plausible-looking vector.

Files modified: `src/cellfate/gen1_predictor.py`, `tests/test_gen1_predictor.py`,
`results/stage24/tool/MODEL_CARD.md`, `results/stage24/tool/io_schema.json`,
`results/stage24/stage24_w5_artifact.json`, and the four Stage-24 substage JSONs that carry the
re-hashed values.

## Freeze state after both commits

```text
  MODEL_CARD.md                8902d4999ad98b9d...   re-hashed, recorded in 24F
  io_schema.json               ceb798affe8d1a7e...   re-hashed, recorded in 24F
  example_clones.csv           330cc7c8a68d77c0...   unchanged
  stage24_oof_for_stage25.csv  abdc2999c0ee388f...   unchanged -- Stage 25 consumed this exact file
  stage24_w5_artifact.json     bda2c2de80187d0f...   re-hashed
  stage24_w5_artifact.npz      954cef7cff296d99...   unchanged
```

**No model, no coefficient, no prediction and no benchmark table changed.** The out-of-fold table
Stage 25 ran against still hashes to the value Stage 24 handed it. Both commits are input-plumbing
and documentation.

## What this does NOT change
Nothing in Stage 24's verdict, Stage 25's inputs, or any recorded number. `STAGE_24_GEN1_TOOL_READY`
stands as recorded.
