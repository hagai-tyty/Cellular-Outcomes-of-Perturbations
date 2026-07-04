"""Chunk planning + fault tolerance (Document 2, S4).

``ProgressTracker`` is implemented in ``cellfate.common`` (shared infrastructure)
and re-exported here so it is importable as ``cellfate.data.chunking.ProgressTracker``.
A *chunk* is the unit of resumable work: each is fetched, processed, and written
to its own shard independently, so a crash only loses the in-flight chunk.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypedDict

from cellfate.common.progress import ProgressTracker  # re-export

if TYPE_CHECKING:  # avoid an import cycle (sources imports CellChunk)
    from .sources import DataSource

__all__ = ["ProgressTracker", "CellChunk", "plan_all"]


class CellChunk(TypedDict):
    """One independently-processable unit of source data."""
    id: str           # stable, e.g. "tahoe:HEPG2:rapamycin"
    uri: str          # how to fetch it (HF path / census query / GEO id)
    source: str       # "tahoe" | "sciplex" | "synthetic" | ...
    cell_line: str
    pert_ids: list[str]


def plan_all(sources: list[DataSource]) -> list[tuple[DataSource, CellChunk]]:
    """Flatten every source's plan into a single ordered work list.

    Chunk ids must be globally unique across sources; a collision is a
    configuration error and is raised eagerly.
    """
    work: list[tuple[DataSource, CellChunk]] = []
    seen: set[str] = set()
    for src in sources:
        for chunk in src.plan():
            cid = chunk["id"]
            if cid in seen:
                raise ValueError(f"duplicate chunk id across sources: {cid!r}")
            seen.add(cid)
            work.append((src, chunk))
    return work
