"""The frozen gene panel: the ordered list of genes that defines the model
input ``X``. It is a *contract* shared by data (which produces ``X``) and
inference (which must validate that a bundle was trained on the same order).

The heavy ``fit_gene_panel`` (which needs scanpy) lives in ``cellfate.data``
and returns a :class:`GenePanel`; this module keeps the panel dependency-light.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


class GenePanel:
    """An ordered, immutable list of gene symbols with a stable content hash."""

    def __init__(self, genes: list[str]) -> None:
        if not genes:
            raise ValueError("GenePanel cannot be empty.")
        if len(set(genes)) != len(genes):
            dupes = sorted({g for g in genes if genes.count(g) > 1})
            raise ValueError(f"GenePanel has duplicate genes: {dupes[:5]}...")
        self._genes: tuple[str, ...] = tuple(genes)

    @property
    def genes(self) -> tuple[str, ...]:
        return self._genes

    def __len__(self) -> int:
        return len(self._genes)

    def __iter__(self):
        return iter(self._genes)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, GenePanel) and other._genes == self._genes

    def hash(self) -> str:
        """16-char content hash; ties an artefact's ``X`` order to the model."""
        joined = "\n".join(self._genes).encode("utf-8")
        return hashlib.sha256(joined).hexdigest()[:16]

    # -- persistence -------------------------------------------------------- #
    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f"# CellFate-Rx gene panel | n={len(self)} | hash={self.hash()}\n"
        path.write_text(header + "\n".join(self._genes) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> GenePanel:
        """Load a panel file. Lines starting with ``#`` and blank lines are ignored."""
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        genes = [ln.strip() for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
        return cls(genes)
