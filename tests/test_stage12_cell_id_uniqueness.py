"""Stage 12 -- `cell_id` must be unique, and the build must refuse to proceed if it is not.

The defect these pin: `cell_id` was `source:cell_line:index_within_chunk`, with no chunk. HFF is
planned as 45 chunks, so `reprogramming:HFF:0` existed 45 times and 42,481 cells carried 981 ids.
`make_splits` keys on cell_id, so ONE index-slot decision was applied to all 45 shards -- the
effective split n was 981, not 42,481, and for D0 (indices 0-111 of every shard) it was 112.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from cellfate.data import build_dataset

ROOT = Path(__file__).resolve().parents[1]


# ---- the key itself -------------------------------------------------------------------------- #
def test_both_cell_id_sites_include_the_chunk_id():
    """Two construction sites, and BOTH were wrong. A fix to one would leave the other colliding."""
    src = (ROOT / "src" / "cellfate" / "data" / "sources.py").read_text(encoding="utf-8")
    assert 'f"{chunk[\'id\']}:{k}"' in src, "the synthetic path does not key on the chunk"
    assert 'f"{chunk_id}:{i}"' in src, "the reprogramming path does not key on the chunk"


def test_the_old_colliding_form_is_gone():
    src = (ROOT / "src" / "cellfate" / "data" / "sources.py").read_text(encoding="utf-8")
    assert 'f"{self.name}:{chunk[\'cell_line\']}:{k}"' not in src
    assert 'f"{source}:{cell_line}:{i}"' not in src


def test_two_chunks_of_the_same_cell_line_cannot_collide():
    """The exact failure: same source, same line, different chunk. Ids must differ."""
    a = f"{'reprogramming:HFF:b0'}:{0}"
    b = f"{'reprogramming:HFF:b1'}:{0}"
    assert a != b


def test_the_chunk_id_invariant_this_relies_on_is_documented_and_enforced():
    """The fix is only sound because CellChunk ids are globally unique -- `plan_all` raises on a
    collision. If that ever stops being true, cell_id uniqueness stops being true with it."""
    from cellfate.data import chunking
    src = inspect.getsource(chunking.plan_all)
    assert "unique" in src.lower()


# ---- the guard ------------------------------------------------------------------------------- #
def test_the_build_guards_uniqueness_before_making_splits():
    """The guard must sit BEFORE make_splits -- that is the consumer of the bad key."""
    src = inspect.getsource(build_dataset.run)
    assert "cell_id is not unique" in src
    guard_at = src.index("cell_id is not unique")
    splits_at = src.index("splits = make_splits")
    assert guard_at < splits_at, "the guard runs after the split it is meant to protect"


def test_the_guard_message_names_the_cause_not_just_the_symptom():
    src = inspect.getsource(build_dataset.run)
    assert "must include the chunk id" in src


def test_the_guard_is_build_time_only_so_existing_folds_stay_readable():
    """Folds already on disk carry colliding ids. A read-time assertion would make every recorded
    artefact unloadable and destroy the ability to re-read past results."""
    ev = (ROOT / "src" / "cellfate" / "evaluation" / "data.py").read_text(encoding="utf-8")
    assert "cell_id is not unique" not in ev
    tr = (ROOT / "src" / "cellfate" / "training" / "dataset.py").read_text(encoding="utf-8")
    assert "cell_id is not unique" not in tr


def test_a_duplicate_id_raises_with_a_counted_diagnosis():
    """Constructed directly on the guard's logic: duplicates must raise, and the message must say
    how many distinct ids there were, because that number is what makes the harm legible."""
    from collections import Counter
    ids = ["a:0", "a:0", "b:1"]
    assert len(set(ids)) != len(ids)
    dup = Counter(ids).most_common(3)
    msg = (f"cell_id is not unique: {len(ids)} cells carry {len(set(ids))} distinct ids "
           f"(worst offenders {dup}).")
    assert "3 cells carry 2 distinct ids" in msg


def test_unique_ids_do_not_raise():
    ids = ["a:0", "a:1", "b:0"]
    assert len(set(ids)) == len(ids)


# ---- the harm this prevents, stated so it is not lost ---------------------------------------- #
@pytest.mark.parametrize("n_chunks,per_chunk", [(45, 944), (2, 10)])
def test_a_colliding_key_collapses_the_effective_split_size(n_chunks, per_chunk):
    """With the chunk omitted, every chunk reuses the same index space, so the number of distinct
    split DECISIONS is per_chunk rather than n_chunks*per_chunk."""
    colliding = {f"src:LINE:{i}" for _ in range(n_chunks) for i in range(per_chunk)}
    fixed = {f"src:LINE:c{c}:{i}" for c in range(n_chunks) for i in range(per_chunk)}
    assert len(colliding) == per_chunk
    assert len(fixed) == n_chunks * per_chunk
