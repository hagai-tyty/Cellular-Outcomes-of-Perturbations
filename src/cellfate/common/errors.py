"""Typed exceptions. Every package raises these (never bare ``Exception``) so
callers can handle failures precisely and messages are actionable.
"""

from __future__ import annotations


class CellFateError(Exception):
    """Base class for all CellFate-Rx errors."""


class SchemaError(CellFateError):
    """A data object violates its schema/contract (wrong shape, bad value)."""


class ContractViolation(CellFateError):
    """An on-disk artefact does not match the expected contract / schema version."""


class GenePanelMismatch(CellFateError):
    """The gene panel of an artefact does not match what a consumer expects.

    Raised, e.g., when a deployment bundle was trained on a different feature
    order than the input ``X`` provides. Carries both hashes for debugging.
    """

    def __init__(self, expected: str, got: str, detail: str = "") -> None:
        msg = (
            f"Gene-panel hash mismatch: expected {expected!r}, got {got!r}. "
            f"The model's input feature order differs from the data's. {detail}"
        ).strip()
        super().__init__(msg)
        self.expected = expected
        self.got = got


class ShardIOError(CellFateError):
    """Failure reading or writing a Parquet shard."""


class BundleError(CellFateError):
    """A deployment bundle is missing files or is internally inconsistent."""


class ConfigError(CellFateError):
    """Invalid or missing configuration."""


class DataSourceError(CellFateError):
    """A data-source connector failed to plan or stream a chunk."""


class NotImplementedInFoundation(CellFateError):
    """Raised by skeleton stubs whose implementation lives in Documents 2-5."""

    def __init__(self, what: str, document: str) -> None:
        super().__init__(
            f"{what} is not part of the foundation (Document 1). "
            f"Implement it per {document}."
        )
