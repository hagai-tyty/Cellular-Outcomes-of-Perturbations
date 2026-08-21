"""The suite must pass on CI, and that must be checkable BEFORE pushing.

CI runs on ubuntu-latest. Two things are true there and false on the development machine:

    1. `D:\\` does not exist, so every dataset root is absent.
    2. paths are POSIX, so a Windows-separator string stored in a committed artifact
       (`results\\table.tsv`) is one filename containing a backslash, not a path.

Both have now produced a red X more than once — Stage 21A (tests reading gitignored fold
directories), and Stage 22 (two Stage-21B tests calling the audit against `D:\\Gill`, plus
`clone_table` written by `str(Path.relative_to(...))` on Windows). The pattern is always the same:
a test that passes locally because the machine happens to have something CI does not.

This module closes that loop. `test_the_whole_suite_survives_the_ci_condition` re-runs the
data-dependent test modules in a subprocess with `CELLFATE_NO_LOCAL_DATA=1`, which points every
local dataset root at a directory that does not exist. If a new test quietly depends on `D:\\`, it
fails HERE, on the development machine, instead of on GitHub.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

# Modules whose tests exercise code that reaches for a local dataset root. A module belongs here
# the moment it imports an `experiments/` script that defines a `D:\` constant.
DATA_DEPENDENT_MODULES = [
    "tests/test_diag_stage21_data_audit.py",
    "tests/test_diag_stage21b_source_design.py",
    "tests/test_diag_stage21d_public_reconstruction.py",
    "tests/test_stage22_prospective_benchmarks.py",
]

# Experiment modules that own a local root and must therefore honour the switch.
SWITCHED_MODULES = [
    "experiments/diag_stage21b_source_design.py",
    "experiments/diag_stage21d_public_reconstruction.py",
]


def test_the_whole_suite_survives_the_ci_condition():
    """Re-run the data-dependent modules with every local dataset root made absent.

    Passing or skipping is fine; failing is not. This is the exact condition GitHub runs under, so
    a green result here means the push will be green too.
    """
    env = dict(os.environ, CELLFATE_NO_LOCAL_DATA="1", PYTHONUTF8="1",
               PYTEST_DISABLE_PLUGIN_AUTOLOAD="")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", *DATA_DEPENDENT_MODULES,
         "-q", "-p", "no:cacheprovider", "--tb=line", "-o", "addopts="],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=900)
    tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-25:])
    assert proc.returncode == 0, (
        "a test depends on local data that CI does not have.\n"
        "Reproduce with:  CELLFATE_NO_LOCAL_DATA=1 pytest " + " ".join(DATA_DEPENDENT_MODULES)
        + "\n\n" + tail)


@pytest.mark.parametrize("rel", SWITCHED_MODULES)
def test_every_local_root_honours_the_switch(rel):
    """A `D:\\` literal is fine; a module-level root that ignores the switch is not.

    Checked on the syntax tree rather than line by line: a root may be written as a multi-line
    conditional, and a line-wise scan flags its `else` branch as if it were unguarded.
    """
    import ast

    src = (ROOT / rel).read_text(encoding="utf-8")
    assert "_NO_LOCAL_DATA" in src, f"{rel} defines a local root but cannot simulate CI"
    for node in ast.parse(src).body:
        if not isinstance(node, ast.Assign):
            continue
        literals = [n.value for n in ast.walk(node)
                    if isinstance(n, ast.Constant) and isinstance(n.value, str)
                    and n.value.startswith("D:")]
        if not literals:
            continue
        names = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
        target = ", ".join(t.id for t in node.targets if isinstance(t, ast.Name))
        assert "_NO_LOCAL_DATA" in names or "_ABSENT_ROOT" in names, (
            f"{rel}: `{target}` points at {literals} unconditionally, so CELLFATE_NO_LOCAL_DATA "
            f"cannot simulate CI for it. Wrap it: `X = _ABSENT_ROOT if _NO_LOCAL_DATA else ...`")


@pytest.mark.parametrize("name", sorted(p.name for p in RESULTS.glob("*.json")))
def test_no_committed_result_stores_a_windows_path_separator(name):
    """`str(Path.relative_to(ROOT))` yields `results\\x.tsv` on Windows, which is not a path on
    Linux. Repo-relative references must be POSIX. This is what broke the Stage-22 CI run."""
    raw = (RESULTS / name).read_text(encoding="utf-8")
    try:
        doc = json.loads(raw)
    except json.JSONDecodeError:
        pytest.skip(f"{name} is not JSON")

    bad = []

    def walk(node, path="$"):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{path}.{k}")
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")
        elif isinstance(node, str) and "\\" in node:
            # Two things are legitimate and must not be flagged:
            #   * an absolute local path -- provenance, "the file lived here when this ran"
            #   * prose that happens to name a directory, e.g. "checked 11 paths under D:\..."
            # What is NOT legitimate is a bare REPO-RELATIVE reference, because something will
            # later try to open it. Those are single tokens: no whitespace, no drive letter.
            looks_absolute = len(node) > 2 and node[1] == ":"
            is_prose = any(c.isspace() for c in node)
            if not looks_absolute and not is_prose:
                bad.append((path, node))

    walk(doc)
    assert not bad, (
        f"{name} stores repo-relative path(s) with a Windows separator; use `.as_posix()`:\n  "
        + "\n  ".join(f"{p} = {v}" for p, v in bad))


def test_the_committed_clone_tables_are_reachable_by_their_recorded_path():
    """The concrete failure: Stage 22 read `clone_table` out of the Stage-21D results and opened
    `ROOT / that`. On Linux it was a single filename containing a backslash."""
    f = RESULTS / "diag_stage21d_public_reconstruction_results.json"
    if not f.exists():
        pytest.skip("Stage 21D has not been run")
    rec = json.loads(f.read_text(encoding="utf-8"))
    for gse in ("GSE227151", "GSE279162"):
        rel = rec[gse]["clone_table"]
        assert "\\" not in rel, f"{gse}: {rel!r} is not portable"
        assert (ROOT / rel).exists(), f"{gse}: {rel} does not resolve"


def test_the_ci_workflow_still_runs_the_full_suite():
    """A green X obtained by narrowing what CI runs would be worse than the red one."""
    wf = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert "pytest -q" in wf, "CI must run the whole suite, unfiltered"
    assert "ruff check src/ tests/ scripts/ plan_tests/" in wf
    for narrowing in ("--ignore", "-k ", "--deselect", "continue-on-error"):
        assert narrowing not in wf, f"CI must not be narrowed with {narrowing!r}"
