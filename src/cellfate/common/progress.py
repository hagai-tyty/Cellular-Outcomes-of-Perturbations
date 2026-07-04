"""ProgressTracker: a persistent chunk-completion ledger that makes the ETL
resumable. Lives in ``common`` because fault tolerance is shared infrastructure;
``cellfate.data.chunking`` re-exports it so Document 2 finds it where expected.
"""

from __future__ import annotations

from pathlib import Path

from .io import read_json, write_json
from .schemas import ChunkDone, ChunkFailed, ProgressState


class ProgressTracker:
    """Tracks which chunks are done/failed, checkpointing after every update."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        if self.path.exists():
            self.state = ProgressState.model_validate(read_json(self.path))
        else:
            self.state = ProgressState()
            self._flush()

    def _flush(self) -> None:
        write_json(self.path, self.state.model_dump())

    def is_done(self, chunk_id: str) -> bool:
        return chunk_id in self.state.done

    def mark_done(self, chunk_id: str, n_samples: int) -> None:
        self.state.done[chunk_id] = ChunkDone(n=n_samples)
        self.state.failed.pop(chunk_id, None)
        self._flush()

    def mark_failed(self, chunk_id: str, err: str) -> None:
        self.state.failed[chunk_id] = ChunkFailed(err=err[:500])
        self._flush()

    def pending(self, chunk_ids: list[str]) -> list[str]:
        return [c for c in chunk_ids if not self.is_done(c)]

    @property
    def n_done(self) -> int:
        return len(self.state.done)

    @property
    def n_samples(self) -> int:
        return sum(c.n for c in self.state.done.values())
