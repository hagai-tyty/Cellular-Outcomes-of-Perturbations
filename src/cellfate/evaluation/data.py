"""Materialise one split of one regime as plain arrays (Document 5).

Reads shards + the regime's split assignment directly; keeps ``X`` in the training
``Sample.X`` space (unscaled) so the model scales it internally and baselines can
standardise it themselves. The test split is touched only by the evaluator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellfate.common import io
from cellfate.common.io import ArtifactPaths, load_splits

_ARRAY_KEYS = ("X", "u_chem_fp", "dose_time", "y_cls", "y_age", "age_mask",
               "scaffold_id", "cell_line", "cell_id")


@dataclass
class SplitData:
    X: np.ndarray            # (N, G) log-normalised panel expression (unscaled)
    fp: np.ndarray           # (N, 2048) fingerprint bits
    dose_time: np.ndarray    # (N, 2) [log10 dose, log time]
    y_cls: np.ndarray        # (N,) class in {0=safe, 1=loss, 2=death}
    y_age: np.ndarray        # (N,) ΔAge vs control (meaningful where mask)
    mask: np.ndarray         # (N,) bool: age is valid
    scaffold_id: np.ndarray
    cell_line: np.ndarray
    cell_id: np.ndarray

    @property
    def n(self) -> int:
        return len(self.X)

    @property
    def y1h(self) -> np.ndarray:
        oh = np.zeros((self.n, 3), dtype=np.float64)
        if self.n:
            oh[np.arange(self.n), self.y_cls.astype(int)] = 1.0
        return oh


def gather_split(paths: ArtifactPaths, regime: str, split: str) -> SplitData:
    assign = load_splits(paths, regime)
    wanted = {cid for cid, sp in assign.items() if sp == split}
    acc: dict[str, list] = {k: [] for k in _ARRAY_KEYS}
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        arr = io.shard_to_numpy(io.read_shard(shard))
        if arr["u_chem_fp"] is None:
            continue
        ids = arr["cell_id"]
        keep = np.fromiter((c in wanted for c in ids), bool, len(ids))
        if not keep.any():
            continue
        for k in _ARRAY_KEYS:
            acc[k].append(np.asarray(arr[k])[keep])

    def cat(k, dtype=None):
        if not acc[k]:
            return np.array([], dtype=dtype)
        out = np.concatenate(acc[k])
        return out.astype(dtype) if dtype else out

    # y_cls is stored as a soft (N,3) distribution; the hard label is its argmax
    y_cls_soft = cat("y_cls", np.float64)
    y_cls = (np.argmax(y_cls_soft, axis=1).astype(np.int64)
             if y_cls_soft.ndim == 2 and len(y_cls_soft) else np.array([], dtype=np.int64))

    return SplitData(
        X=cat("X", np.float32),
        fp=cat("u_chem_fp", np.float32),
        dose_time=cat("dose_time", np.float32),
        y_cls=y_cls,
        y_age=cat("y_age", np.float64),
        mask=cat("age_mask", bool),
        scaffold_id=cat("scaffold_id"),
        cell_line=cat("cell_line"),
        cell_id=cat("cell_id"),
    )
