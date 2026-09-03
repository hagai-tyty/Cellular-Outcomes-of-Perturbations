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

# `st_size` is deliberately NOT banned outright -- a loose bound ("under 2 MB") is portable and
# several stages use one. But it is NOT a property of the content: `results/**` is `text eol=lf`,
# so a checkout lands LF while a Python `write_text` on Windows lands CRLF, and the same file
# then has two different sizes. Recording `st_size` as if it described the file is what made the
# evidence manifest churn on every re-run. A size COMPARED against a recorded value must be
# measured on the same normalised content that was hashed.
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

# ============================================================================================== #
# Re-running a stage must not dirty the working tree
# ============================================================================================== #
VOLATILE = {"runtime_seconds", "runtime_minutes", "total_runtime_seconds"}

GEN1_STAGE_OUTPUTS = [
    "results/evidence_lock", "results/claim_lock", "results/manuscript",
]


def _volatile_paths(obj, path=""):
    hits = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in VOLATILE:
                hits.append(f"{path}.{k}".lstrip("."))
            hits += _volatile_paths(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for v in obj:
            hits += _volatile_paths(v, path)
    return hits


def test_no_committed_gen1_stage_output_carries_a_timestamp():
    """Timing written into a committed output makes every re-run dirty the tree, and in the claim
    lock it moved the digest itself. Stage 23.2 established this convention first -- its contract
    forbids timing inside a hashed payload and its docstring records that the class bit three
    times. The Gen-1 locks reintroduced it; this stops a third occurrence."""
    import json
    offenders = {}
    for d in GEN1_STAGE_OUTPUTS:
        base = ROOT / d
        if not base.is_dir():
            continue
        for p in sorted(base.rglob("*.json")):
            hits = _volatile_paths(json.loads(p.read_text(encoding="utf-8")))
            if hits:
                offenders[p.relative_to(ROOT).as_posix()] = hits
    assert not offenders, f"volatile timing in committed stage outputs: {offenders}"


def test_each_gen1_lock_writer_strips_timing():
    """Directly, at the writer, so the property holds by construction rather than by discipline."""
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    for name in ("run_gen1_evidence_lock", "run_gen1_claim_lock", "run_gen1_manuscript"):
        mod = __import__(name)
        assert hasattr(mod, "_strip_volatile"), f"{name} has no volatile stripper"
        stripped = mod._strip_volatile({"a": 1, "runtime_seconds": 9.9,
                                        "n": [{"runtime_minutes": 1, "b": 2}]})
        assert stripped == {"a": 1, "n": [{"b": 2}]}, f"{name} strips incompletely"


def test_manifest_sizes_are_measured_on_hashed_content():
    """The manifest's `bytes` must describe the content, not this filesystem.

    `st_size` counts CRLF as two bytes. Because `results/**` is `text eol=lf`, a file written by a
    stage on Windows (CRLF) and the same file in a fresh checkout (LF) have identical content and
    different `st_size` -- so a manifest built from `st_size` disagrees with itself across checkouts
    and re-dirties the tree on every run. Sizes must come from the bytes that were hashed.
    """
    import json
    import sys
    manifest = ROOT / "results" / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json"
    if not manifest.is_file():
        pytest.skip("the evidence lock has not been built")
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_evidence_lock as EL

    wrong = []
    for rel, meta in sorted(json.loads(manifest.read_text(encoding="utf-8"))["artifacts"].items()):
        if not (ROOT / rel).is_file():
            continue                      # gitignored artifact; its absence is checked elsewhere
        measured = EL.artifact_bytes(rel)
        if meta["bytes"] != measured:
            wrong.append(f"{rel}: manifest {meta['bytes']} != hashed-content {measured}")
    assert not wrong, "manifest sizes do not describe the hashed content: " + "; ".join(wrong)


def test_the_size_measurement_is_line_ending_independent(tmp_path):
    """A guard that would also pass under the old buggy implementation is worthless."""
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_evidence_lock as EL

    crlf, lf = tmp_path / "a.json", tmp_path / "b.json"
    crlf.write_bytes(b'{"a":1}\r\n{"b":2}\r\n')
    lf.write_bytes(b'{"a":1}\n{"b":2}\n')
    assert crlf.stat().st_size != lf.stat().st_size, "the fixture must differ on disk"

    original = EL.ROOT
    try:
        EL.ROOT = tmp_path
        assert EL.artifact_bytes("a.json") == EL.artifact_bytes("b.json"), \
            "identical content measured to different sizes -- the CRLF bug is back"
    finally:
        EL.ROOT = original
