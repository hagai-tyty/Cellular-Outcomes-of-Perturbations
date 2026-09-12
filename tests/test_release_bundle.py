"""Contracts for the release bundle's placeholder gate.

A Zenodo record's files cannot be changed after it is published, so a FILL marker left in a document
the archive carries would be permanent. The builder must refuse one. LATER markers -- fields that can
only exist after archiving, such as the preprint DOI -- must be allowed and reported separately.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import make_release_bundle as B  # noqa: E402


def _tree(root: Path, files: dict[str, str]) -> Path:
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return root


def test_a_fill_marker_is_found_in_every_release_document(tmp_path):
    for i, rel in enumerate(B.RELEASE_DOCUMENTS):
        root = _tree(tmp_path / str(i), {rel: "before <<FILL: author name>> after\n"})
        fill, later = B.placeholders(root)
        assert fill == [f"{rel}: <<FILL: author name>>"], rel
        assert later == [], rel


def test_a_marker_wrapped_across_lines_is_still_found(tmp_path):
    root = _tree(tmp_path, {"README.md": "see <<FILL: describe what the tools\ndid and did not do>> here\n"})
    fill, _ = B.placeholders(root)
    assert fill == ["README.md: <<FILL: describe what the tools did and did not do>>"]


def test_later_markers_are_allowed_and_reported_separately(tmp_path):
    root = _tree(tmp_path, {"results/manuscript/SUBMISSION.md": "DOI `<<LATER: bioRxiv DOI>>`\n"})
    fill, later = B.placeholders(root)
    assert fill == []
    assert later == ["results/manuscript/SUBMISSION.md: <<LATER: bioRxiv DOI>>"]


def test_prose_that_explains_the_convention_is_not_a_marker(tmp_path):
    """The submission pack explains its markers; the explanation must not trip the gate itself."""
    root = _tree(tmp_path, {"results/manuscript/SUBMISSION.md":
                            "Fields marked FILL, in double angle brackets, need a human.\n"})
    assert B.placeholders(root) == ([], [])


def test_the_gate_sits_in_build_before_the_archive_is_written():
    """A helper nothing calls protects nothing."""
    src = inspect.getsource(B.build)
    assert "placeholders()" in src
    assert src.index("placeholders()") < src.index("zipfile.ZipFile")


def test_the_gate_is_not_bypassed_by_allow_dirty():
    """--allow-dirty exists for throwaway builds; it must not also let a FILL marker through."""
    src = inspect.getsource(B.build)
    gate = src[src.index("fill, later = placeholders()"):src.index("zipfile.ZipFile")]
    assert "allow_dirty" not in gate
