"""Artefact input/output: the only module that reads/writes the on-disk
contracts. Everything is crash-safe (atomic write-then-rename) and goes through
the canonical paths in :class:`ArtifactPaths`, so no path string is duplicated.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import constants as C
from .errors import BundleError, ContractViolation, GenePanelMismatch, ShardIOError
from .panel import GenePanel
from .schemas import (
    BundleMeta,
    ConformalParams,
    ManifestRow,
    ResParams,
    Sample,
    TemperatureParams,
)


# --------------------------------------------------------------------------- #
# Atomic primitives                                                           #
# --------------------------------------------------------------------------- #
def atomic_write_bytes(path: str | Path, data: bytes) -> None:
    """Write bytes atomically: temp file in the same dir, then ``os.replace``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def write_json(path: str | Path, obj: object) -> None:
    atomic_write_bytes(path, json.dumps(obj, indent=2, default=str).encode("utf-8"))


def read_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sanitize_id(raw: str) -> str:
    """Turn a chunk id like ``tahoe:HEPG2:rapamycin`` into a filesystem-safe stem."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", raw).strip("_")


# --------------------------------------------------------------------------- #
# Canonical paths                                                             #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ArtifactPaths:
    """All standard paths derived from a single run directory."""

    root: Path

    @classmethod
    def of(cls, root: str | Path) -> ArtifactPaths:
        return cls(Path(root))

    # dataset artefacts
    @property
    def shards_dir(self) -> Path: return self.root / C.SHARDS_DIRNAME
    @property
    def manifest_parts_dir(self) -> Path: return self.root / C.MANIFEST_PARTS_DIRNAME
    @property
    def manifest_file(self) -> Path: return self.root / C.MANIFEST_FILENAME
    @property
    def splits_dir(self) -> Path: return self.root / C.SPLITS_DIRNAME
    @property
    def scalers_file(self) -> Path: return self.root / C.SCALERS_FILENAME
    @property
    def progress_file(self) -> Path: return self.root / C.PROGRESS_FILENAME
    @property
    def reports_dir(self) -> Path: return self.root / C.REPORTS_DIRNAME

    def shard_file(self, shard_id: str) -> Path:
        return self.shards_dir / f"{shard_id}.parquet"

    def split_file(self, regime: str) -> Path:
        return self.splits_dir / f"{regime}.json"

    # bundle artefacts
    @property
    def bundle_dir(self) -> Path: return self.root / C.BUNDLE_DIRNAME
    @property
    def bundle_members_dir(self) -> Path: return self.bundle_dir / C.BUNDLE_MEMBERS_DIRNAME
    @property
    def bundle_ood_dir(self) -> Path: return self.bundle_dir / C.BUNDLE_OOD_DIRNAME
    @property
    def bundle_meta_file(self) -> Path: return self.bundle_dir / C.BUNDLE_META_FILENAME
    @property
    def bundle_temperature_file(self) -> Path: return self.bundle_dir / C.BUNDLE_TEMPERATURE_FILENAME
    @property
    def bundle_conformal_file(self) -> Path: return self.bundle_dir / C.BUNDLE_CONFORMAL_FILENAME
    @property
    def bundle_res_file(self) -> Path: return self.bundle_dir / C.BUNDLE_RES_FILENAME
    @property
    def bundle_scalers_file(self) -> Path: return self.bundle_dir / C.SCALERS_FILENAME


# --------------------------------------------------------------------------- #
# Shard schema + IO                                                           #
# --------------------------------------------------------------------------- #
_F32 = pa.list_(pa.float32())
_U8 = pa.list_(pa.uint8())

SHARD_SCHEMA = pa.schema([
    ("cell_id", pa.string()),
    ("X", _F32),
    ("u_modality", pa.string()),
    ("u_chem_fp", _U8),
    ("u_gene_emb", _F32),
    ("u_tf_emb", _F32),
    ("dose_time", _F32),
    ("y_cls", _F32),
    ("y_age", pa.float32()),
    ("age_mask", pa.bool_()),
    ("sig_scores", _F32),
    ("cell_line", pa.string()),
    ("pert_id", pa.string()),
    ("scaffold_id", pa.string()),
    ("source", pa.string()),
])


def write_shard(path: str | Path, samples: list[Sample]) -> None:
    """Write a list of validated :class:`Sample` to a Parquet shard atomically."""
    if not samples:
        raise ShardIOError("refusing to write an empty shard")
    cols: dict[str, list] = {f.name: [] for f in SHARD_SCHEMA}
    for s in samples:
        cols["cell_id"].append(s.cell_id)
        cols["X"].append(s.X)
        cols["u_modality"].append(s.u_modality.value)
        cols["u_chem_fp"].append(s.u_chem_fp)
        cols["u_gene_emb"].append(s.u_gene_emb)
        cols["u_tf_emb"].append(s.u_tf_emb)
        cols["dose_time"].append(s.dose_time)
        cols["y_cls"].append(s.y_cls)
        cols["y_age"].append(s.y_age)
        cols["age_mask"].append(s.age_mask)
        cols["sig_scores"].append(s.sig_scores)
        cols["cell_line"].append(s.cell_line)
        cols["pert_id"].append(s.pert_id)
        cols["scaffold_id"].append(s.scaffold_id)
        cols["source"].append(s.source)
    table = pa.table(cols, schema=SHARD_SCHEMA)
    buf = pa.BufferOutputStream()
    pq.write_table(table, buf, compression="zstd")
    atomic_write_bytes(path, buf.getvalue().to_pybytes())


def read_shard(path: str | Path) -> pa.Table:
    return pq.read_table(path)


def rewrite_shard_yage(path: str | Path, y_age: np.ndarray) -> None:
    """Overwrite only the ``y_age`` column of an existing shard, in row order.

    All other columns are preserved byte-for-byte. Non-finite entries (NaN) are
    written as null so the round-trip stays consistent with masked rows. Used by
    the ETL to re-apply a train-fit cell-cycle deconfounder to every shard.
    """
    table = pq.read_table(path)
    d = table.to_pydict()
    if len(y_age) != len(d["cell_id"]):
        raise ShardIOError(
            f"rewrite_shard_yage length mismatch: {len(y_age)} != {len(d['cell_id'])}"
        )
    d["y_age"] = [None if not np.isfinite(v) else float(v) for v in y_age]
    new = pa.table(d, schema=SHARD_SCHEMA)
    buf = pa.BufferOutputStream()
    pq.write_table(new, buf, compression="zstd")
    atomic_write_bytes(path, buf.getvalue().to_pybytes())


def shard_to_numpy(table: pa.Table) -> dict[str, object]:
    """Materialise a shard into dense arrays for fast training/eval access.

    Returns X (N,G), dose_time (N,2), y_cls (N,3), sig_scores (N,3),
    y_age (N,) float32 with NaN where masked, age_mask (N,) bool,
    u_chem_fp (N,2048) uint8 (or None if not chem), plus id/meta lists.
    """
    d = table.to_pydict()
    n = len(d["cell_id"])

    def stack(key: str) -> np.ndarray:
        return np.asarray(d[key], dtype=np.float32).reshape(n, -1)

    fp = d["u_chem_fp"]
    chem = np.asarray(fp, dtype=np.uint8) if all(v is not None for v in fp) else None
    tf = d["u_tf_emb"]
    tf_emb = (np.asarray(tf, dtype=np.float32).reshape(n, -1)
              if all(v is not None for v in tf) else None)
    y_age = np.asarray([np.nan if v is None else v for v in d["y_age"]], dtype=np.float32)
    return {
        "cell_id": d["cell_id"],
        "X": stack("X"),
        "dose_time": stack("dose_time"),
        "y_cls": stack("y_cls"),
        "sig_scores": stack("sig_scores"),
        "y_age": y_age,
        "age_mask": np.asarray(d["age_mask"], dtype=bool),
        "u_chem_fp": chem,
        "u_tf_emb": tf_emb,
        "u_modality": d["u_modality"],
        "cell_line": d["cell_line"],
        "pert_id": d["pert_id"],
        "scaffold_id": d["scaffold_id"],
        "source": d["source"],
    }


# --------------------------------------------------------------------------- #
# Manifest (fragments -> consolidated)                                        #
# --------------------------------------------------------------------------- #
_MANIFEST_SCHEMA = pa.schema([
    ("cell_id", pa.string()), ("cell_line", pa.string()), ("pert_id", pa.string()),
    ("scaffold_id", pa.string()), ("source", pa.string()), ("age_mask", pa.bool_()),
    ("shard_id", pa.string()), ("row_idx", pa.int64()),
])


def write_manifest_part(paths: ArtifactPaths, shard_id: str, rows: list[ManifestRow]) -> None:
    """Write one shard's manifest fragment (crash-safe; consolidated later)."""
    cols = {f.name: [] for f in _MANIFEST_SCHEMA}
    for r in rows:
        for f in _MANIFEST_SCHEMA:
            cols[f.name].append(getattr(r, f.name))
    table = pa.table(cols, schema=_MANIFEST_SCHEMA)
    buf = pa.BufferOutputStream()
    pq.write_table(table, buf)
    atomic_write_bytes(paths.manifest_parts_dir / f"{shard_id}.parquet", buf.getvalue().to_pybytes())


def consolidate_manifest(paths: ArtifactPaths) -> pa.Table:
    """Merge all manifest fragments into a single ``manifest.parquet`` and return it."""
    parts = sorted(paths.manifest_parts_dir.glob("*.parquet"))
    if not parts:
        raise ContractViolation("no manifest fragments found; did ETL run?")
    table = pa.concat_tables([pq.read_table(p) for p in parts])
    buf = pa.BufferOutputStream()
    pq.write_table(table, buf)
    atomic_write_bytes(paths.manifest_file, buf.getvalue().to_pybytes())
    return table


def load_manifest(paths: ArtifactPaths) -> pa.Table:
    if paths.manifest_file.exists():
        return pq.read_table(paths.manifest_file)
    return consolidate_manifest(paths)


def manifest_rows(table: pa.Table) -> list[ManifestRow]:
    return [ManifestRow.model_validate(r) for r in table.to_pylist()]


# --------------------------------------------------------------------------- #
# Splits                                                                       #
# --------------------------------------------------------------------------- #
def write_splits(paths: ArtifactPaths, regime: str, mapping: dict[str, str]) -> None:
    """Persist a ``cell_id -> split`` mapping for one regime."""
    write_json(paths.split_file(regime), {"regime": regime, "map": mapping})


def load_splits(paths: ArtifactPaths, regime: str) -> dict[str, str]:
    return read_json(paths.split_file(regime))["map"]


# --------------------------------------------------------------------------- #
# Gene-panel + config + deps fingerprints                                     #
# --------------------------------------------------------------------------- #
def assert_gene_panel(expected_hash: str, panel: GenePanel, detail: str = "") -> None:
    if panel.hash() != expected_hash:
        raise GenePanelMismatch(expected_hash, panel.hash(), detail)


def to_container(cfg: object) -> object:
    """Best-effort convert an OmegaConf/pydantic object to a plain container."""
    try:
        from omegaconf import OmegaConf

        if OmegaConf.is_config(cfg):
            return OmegaConf.to_container(cfg, resolve=True)
    except ImportError:
        pass
    if hasattr(cfg, "model_dump"):
        return cfg.model_dump()
    return cfg


def hash_config(cfg: object) -> str:
    """Stable 12-char hash of a resolved config (order-independent)."""
    import hashlib

    blob = json.dumps(to_container(cfg), sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def deps_hash() -> str:
    """Fingerprint of the installed dependency set (name==version)."""
    import hashlib
    from importlib import metadata

    items = sorted(f"{d.metadata['Name']}=={d.version}" for d in metadata.distributions())
    return hashlib.sha256("\n".join(items).encode("utf-8")).hexdigest()[:12]


# --------------------------------------------------------------------------- #
# Bundle helpers (the sole input to inference)                                #
# --------------------------------------------------------------------------- #
def save_bundle_meta(paths: ArtifactPaths, meta: BundleMeta) -> None:
    write_json(paths.bundle_meta_file, meta.model_dump())


def load_bundle_meta(paths: ArtifactPaths) -> BundleMeta:
    return BundleMeta.model_validate(read_json(paths.bundle_meta_file))


def save_temperature(paths: ArtifactPaths, t: TemperatureParams) -> None:
    write_json(paths.bundle_temperature_file, t.model_dump())


def load_temperature(paths: ArtifactPaths) -> TemperatureParams:
    return TemperatureParams.model_validate(read_json(paths.bundle_temperature_file))


def save_conformal(paths: ArtifactPaths, c: ConformalParams) -> None:
    write_json(paths.bundle_conformal_file, c.model_dump())


def load_conformal(paths: ArtifactPaths) -> ConformalParams:
    return ConformalParams.model_validate(read_json(paths.bundle_conformal_file))


def save_res_params(paths: ArtifactPaths, r: ResParams) -> None:
    write_json(paths.bundle_res_file, r.model_dump())


def load_res_params(paths: ArtifactPaths) -> ResParams:
    return ResParams.model_validate(read_json(paths.bundle_res_file))


def assert_bundle_complete(paths: ArtifactPaths) -> None:
    """Verify a bundle has all required files before inference trusts it."""
    required = [
        paths.bundle_meta_file, paths.bundle_temperature_file,
        paths.bundle_conformal_file, paths.bundle_res_file, paths.bundle_scalers_file,
        paths.bundle_members_dir, paths.bundle_ood_dir,
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise BundleError(f"incomplete bundle at {paths.bundle_dir}: missing {missing}")
    if not any(paths.bundle_members_dir.glob("*.pt")):
        raise BundleError(f"no ensemble members (*.pt) in {paths.bundle_members_dir}")
