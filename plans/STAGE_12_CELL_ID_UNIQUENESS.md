# Stage 12 — `cell_id` is not unique, and the split map is keyed on it

**Status:** PLAN, then EXECUTE. This is the first stage in this arc that **changes `src/`**.
**Scope of that change is deliberately narrow** — see §12.5 for what it explicitly does NOT do.

---

## 12.1 The defect, measured

`cell_id` is built as `f"{source}:{cell_line}:{index_within_chunk}"`
(`sources.py:208` and `sources.py:363`). The chunk is **not** in the key. HFF is planned as **45
chunks**, so `reprogramming:HFF:0` exists 45 times.

Measured on a built fold:

- 42,481 HFF cells carry **981 distinct `cell_id`s**.
- `splits/holdout.json` therefore holds **1,100 entries for 42,600 cells**.
- The same index→split decision is applied to all 45 shards.

`chunking.py:33` already states the invariant that makes the fix obvious: *"Chunk ids must be
globally unique across sources; a collision is a configuration error and is raised eagerly."*
`CellChunk.id` exists and is unique. The cell_id construction simply does not use it.

## 12.2 The harm, measured — this is not cosmetic

Each shard contains **all 9 timepoints**, and within a shard **D0 occupies indices 0–111
exclusively** (the other 8 days are interleaved from index 112 on).

So D0's split assignment is decided by **~112 index-slots**, not by ~4,700 D0 cells. The observed
consequence:

| split | share of D0 cells |
|---|---|
| calib | **9.0 %** |
| train | 11.9 % |
| val | **13.3 %** |

`calib` is depleted of the control timepoint by roughly 23 % relative to `val`. **Calib is what
conformal intervals are computed on**, and D0 is the control anchor. The effective sample size of
the split assignment is **981, not 42,481** — and for D0 specifically, **112**.

The ±4 % spread is exactly the sampling noise of n=112, which is what identifies the mechanism
rather than merely correlating with it.

## 12.3 The fix

Include the chunk id in the cell key at both construction sites. `chunk["id"]` is in scope at
`sources.py:208`; `chunk_id` is already in scope at `sources.py:363` (it is used in that function's
error messages).

    before  f"{source}:{cell_line}:{i}"          ->  reprogramming:HFF:0        (x45)
    after   f"{chunk_id}:{i}"                    ->  reprogramming:HFF:b0:0     (unique)

## 12.4 The guard

A fix that can silently regress is half a fix. Add a **build-time** uniqueness assertion: if the
assembled cell_ids are not distinct, raise. This defect survived because nothing ever checked.

**The guard must be BUILD-time, not READ-time.** Existing folds on disk carry colliding ids; a
read-time assertion would make every recorded artefact unloadable and destroy the ability to
re-read past results.

## 12.5 What this stage explicitly does NOT do

- **No rebuild. No re-score. No retrain.** Those are hours of compute and a separate decision.
- **Existing folds are unaffected and remain readable.** `gather_split` reads cell_ids from the
  shard and splits from the splits file — both written together, both old-format, so they stay
  mutually consistent. The fix is **forward-only**: it changes what NEW builds write.
- **Every recorded result stands exactly as measured.** This stage does not revise a single
  number. It makes the next build correct.
- **The size of the effect on model metrics is UNKNOWN** and is not claimed. Quantifying it needs
  a rebuild under the fix and a paired comparison — its own Change, separately pre-registered.

## 12.6 Why this is worth changing `src/` for when the earlier stages were not

Stages 10 and 11 were **inference** questions — what does a number mean. Their answers change
interpretation, and interpretation should not be hard-coded before it is settled.

This is a **correctness** defect with a known mechanism, a measured harm, and a fix that cannot
make anything worse: a unique key is unambiguously more correct than a colliding one, independent
of what any downstream result turns out to be.

## 12.7 Verification

| item | how |
|---|---|
| uniqueness | new tests asserting cell_ids are distinct across chunks for both source paths |
| the guard fires | a test that constructs a colliding assembly and requires it to raise |
| the guard is build-time only | a test that existing read paths carry no such assertion |
| no behaviour change beyond the key | full `pytest` green, `ruff` clean |
| record | `CHANGES.md`, stating the measured harm and that no result is revised |
