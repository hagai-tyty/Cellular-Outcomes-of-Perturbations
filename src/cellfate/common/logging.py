"""Structured logging and run-directory management.

Every stage writes JSON-lines logs (machine-parseable) and runs inside a unique
directory ``artifacts/<timestamp>-<git_sha>-<cfg_hash>/`` so runs are traceable
and reproducible (Document 1, S8).
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path


class _JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "ts": round(record.created, 3),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        extra = getattr(record, "extra_fields", None)
        if extra:
            payload.update(extra)
        if record.exc_info:
            payload["exc"] = self.formatException(record.exc_info)
        return json.dumps(payload, default=str)


def get_logger(name: str = "cellfate", level: int = logging.INFO) -> logging.Logger:
    """Return a process-wide singleton logger emitting JSON lines to stdout."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(_JsonFormatter())
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.propagate = False
    return logger


def log_event(logger: logging.Logger, msg: str, **fields: object) -> None:
    """Log a message with structured extra fields."""
    logger.info(msg, extra={"extra_fields": fields})


def git_sha(default: str = "nogit") -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5, check=True,
        )
        return out.stdout.strip() or default
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return default


def make_run_dir(base: str | Path, stage: str, config_hash: str) -> Path:
    """Create and return ``base/<UTC-timestamp>-<git_sha>-<cfg_hash>-<stage>``."""
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    run = Path(base) / f"{stamp}-{git_sha()}-{config_hash}-{stage}"
    run.mkdir(parents=True, exist_ok=True)
    return run


def write_run_metadata(run_dir: str | Path, *, stage: str, config_hash: str,
                       deps_hash: str, extra: dict | None = None) -> None:
    """Persist run provenance to ``run_dir/run_metadata.json``."""
    meta = {
        "stage": stage,
        "config_hash": config_hash,
        "deps_hash": deps_hash,
        "git_sha": git_sha(),
        "python": sys.version.split()[0],
        "created_at": time.time(),
        **(extra or {}),
    }
    Path(run_dir, "run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
