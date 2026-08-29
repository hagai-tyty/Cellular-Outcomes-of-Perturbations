"""Contracts on the test suite itself — the things that are green locally and red in CI.

CI is the only place this project ever runs from a fresh checkout, so any assertion that depends on
the *working tree* rather than on committed content is invisible on a developer machine and fails
on every build. That is not hypothetical here. An `st_mtime` comparison in the Stage-23.2H contracts
kept CI red continuously, across commits that never touched Stage 23.2, because git does not
preserve modification times: a checkout stamps every file with the same instant in arbitrary order,
so the assertion was a coin flip on whichever filesystem it happened to run on.

Two classes are guarded here, both of which have actually bitten this repository:

    filesystem metadata git does not carry   -> mtime/ctime assertions
    filename case                            -> Windows finds `RECORDS/`, Linux does not

Neither guard is clever. They exist because the expensive part was not fixing either bug, it was
the weeks of red builds before anyone traced them.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"

# `st_size` is deliberately NOT here: git preserves content, so size assertions are portable.
NON_PORTABLE = {
    "st_mtime": "git does not preserve modification times",
    "st_ctime": "git does not preserve creation/change times",
    "getmtime": "git does not preserve modification times",
    "getctime": "git does not preserve creation/change times",
    "st_ino": "inode numbers are not stable across checkouts",
    "st_uid": "ownership is not preserved by git",
}


def _test_sources() -> list[Path]:
    return [p for p in sorted(TESTS.glob("test_*.py")) if p.name != Path(__file__).name]


def _code_lines(text: str) -> list[str]:
    """Comments and docstring prose mention these tokens legitimately; code must not use them."""
    out, in_doc = [], False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.count('"""') == 1:
            in_doc = not in_doc
            continue
        if in_doc or stripped.startswith("#") or stripped.startswith('"""'):
            continue
        out.append(line.split("#", 1)[0])
    return out


def test_no_test_asserts_on_filesystem_metadata_git_discards():
    """The exact bug that kept CI red across dozens of commits."""
    offenders = []
    for src in _test_sources():
        for i, line in enumerate(_code_lines(src.read_text(encoding="utf-8")), 1):
            for token, why in NON_PORTABLE.items():
                if token in line:
                    offenders.append(f"{src.name}:{i}: {token} — {why}")
    assert not offenders, "non-portable assertions:\n  " + "\n  ".join(offenders)


def test_the_metadata_guard_actually_fires():
    """A guard that has never caught anything is decoration."""
    planted = 'def t():\n    assert A.stat().st_mtime <= B.stat().st_mtime\n'
    lines = _code_lines(planted)
    assert any("st_mtime" in ln for ln in lines), "the scanner must see real code"
    # ...and must NOT fire on prose that merely discusses the token
    prose = '# we used to compare st_mtime here, which git does not preserve\nx = 1\n'
    assert not any("st_mtime" in ln for ln in _code_lines(prose)), \
        "a comment explaining the rule must not trip it"


def _case_exact(rel: str) -> bool:
    """True only if every path component matches the on-disk name byte for byte.

    `Path.exists()` is case-insensitive on Windows and case-sensitive on the Linux runner, so a
    mis-cased literal passes locally and fails in CI. This walks the real directory listings.
    """
    cur = ROOT
    for part in Path(rel).parts:
        try:
            names = {p.name for p in cur.iterdir()}
        except (NotADirectoryError, FileNotFoundError, PermissionError):
            return False
        if part not in names:
            return False
        cur = cur / part
    return True


@pytest.mark.parametrize("module_name", ["run_gen1_evidence_lock", "run_gen1_claim_lock",
                                         "run_gen1_manuscript"])
def test_every_path_a_lock_names_matches_the_filesystem_case(module_name):
    """Windows finds `RECORDS/`; the Linux runner does not."""
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    mod = __import__(module_name)

    paths: set[str] = set()
    for attr in ("INVENTORY",):
        if hasattr(mod, attr):
            for group in getattr(mod, attr).values():
                paths.update(group)
    for attr in ("CLAIM_LOCK_FILES", "PACKAGE_FILES"):
        if hasattr(mod, attr):
            paths.update(getattr(mod, attr))

    assert paths, f"{module_name} names no paths — the guard would be vacuous"
    for rel in sorted(paths):
        if not (ROOT / rel).exists():
            continue          # genuinely absent (the gitignored artifact); checked elsewhere
        assert _case_exact(rel), \
            f"{rel} differs in case from the file on disk; the Linux runner will not find it"


def test_documented_commands_use_forward_slashes():
    """A backslash path in a documented command is a Windows-only instruction."""
    doc = ROOT / "results" / "manuscript" / "REPRODUCIBILITY.md"
    if not doc.exists():
        pytest.skip("the package document has not been written")
    for m in re.finditer(r"python (experiments\S+)", doc.read_text(encoding="utf-8")):
        assert "\\" not in m.group(1), f"{m.group(1)} is not runnable on Linux"
