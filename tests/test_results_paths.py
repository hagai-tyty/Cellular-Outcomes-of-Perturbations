"""Every script that writes a `*_results.json` must write it into `results/`, correctly.

WHY THIS FILE EXISTS
--------------------
On 2026-08-01 the 19 `*_results.json` files were moved out of the repo root into `results/`, and
the writers were repointed by a regex: `Path("x.json")` -> `_RESULTS / "x.json"`. That rewrite was
wrong in **20 places** in a way no existing test could catch:

    _RESULTS / "x.json".write_text(...)      # `.write_text` binds to the STRING, not the Path
    (_RESULTS / "x.json").write_text(...)    # correct

Python's `.` binds tighter than `/`, so the broken form raises
`AttributeError: 'str' object has no attribute 'write_text'` -- but **only when `main()` actually
runs**. The unit tests exercise the pure functions and read the recorded JSON; none of them calls
`main()`, so all 567 passed against code that could not write its own output.

These tests are static: they read the source text rather than executing the scripts, so they need
no data, no GPU and no network, and they run in CI in milliseconds.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Every place a *_results.json is produced. Kept as a glob rather than a list so a new script is
# covered the moment it lands.
SCRIPT_DIRS = [ROOT / "experiments", ROOT / "plan_tests", ROOT]

# Any call that could put a file on disk. Used to tell a genuine WRITER from a script that
# merely reads a `*_results.json` as input -- the distinction the `_RESULTS` check turns on.
_WRITE_CALL = re.compile(r"\.write_text\s*\(|\.write_bytes\s*\(|json\.dump\s*\(")


def _scripts() -> list[Path]:
    out: list[Path] = []
    for d in SCRIPT_DIRS:
        out.extend(p for p in d.glob("*.py") if p.is_file())
    return sorted(set(out))


def _writers() -> list[Path]:
    return [p for p in _scripts() if "_results.json" in p.read_text(encoding="utf-8")]


def test_there_are_writers_to_check():
    """If this fails the glob is wrong and every other test here is vacuous."""
    assert len(_writers()) >= 15


@pytest.mark.parametrize("path", _writers(), ids=lambda p: p.name)
def test_no_unparenthesised_path_division_before_a_method_call(path: Path):
    """The exact 2026-08-01 defect: `_RESULTS / "x.json".write_text(...)`.

    Matches any `<name> / "<something>.json".<method>` -- the bug is the missing parentheses,
    not the particular constant name, so a future `_OUT / "x.json".unlink()` is caught too.
    """
    bad = re.findall(r'\b\w+ / "[^"]+\.json"\.\w+', path.read_text(encoding="utf-8"))
    assert not bad, (
        f"{path.name}: `.` binds tighter than `/`, so these call the method on the STRING:\n  "
        + "\n  ".join(bad) + "\nWrap the path: (_RESULTS / \"x.json\").method(...)")


@pytest.mark.parametrize("path", _writers(), ids=lambda p: p.name)
def test_results_json_is_never_written_to_a_bare_relative_path(path: Path):
    """`Path("x_results.json")` resolves against the CWD, so the file lands wherever the user
    happened to be standing. That is what put 19 JSONs in the repo root."""
    bad = re.findall(r'Path\("[a-z0-9_]+_results\.json"\)', path.read_text(encoding="utf-8"))
    assert not bad, (f"{path.name}: writes {bad} relative to the working directory; "
                     "use the module's `_RESULTS` constant so the location is __file__-relative")


@pytest.mark.parametrize("path", _writers(), ids=lambda p: p.name)
def test_every_writer_defines_a_results_constant_pointing_at_the_repo_root(path: Path):
    """`_RESULTS` must resolve to `<repo>/results` from the script's own location, not the CWD."""
    t = path.read_text(encoding="utf-8")
    if "_RESULTS" not in t:
        # The skip below used to be unconditional, with the message "reads results but does not
        # write any" -- an ASSUMPTION, not a check. `verify_rev_final_4_4.py` mentioned a
        # `*_results.json`, defined no `_RESULTS`, and wrote to `root / "..._results.json"`; it was
        # skipped under a message asserting it did not write. The tidy-up had moved its output to
        # `results/`, so the next run would have dropped a stray JSON back into the repo root and
        # turned `test_no_results_json_is_left_in_the_repo_root` red -- a latent failure that this
        # test existed to prevent and silently waved through.
        assert not _WRITE_CALL.search(t), (
            f"{path.name}: mentions a *_results.json AND writes files, but defines no _RESULTS "
            "constant, so this test cannot check where it writes. Adopt the convention:\n"
            '  _RESULTS = Path(__file__).resolve().parents[N] / "results"')
        pytest.skip("mentions a results file but writes nothing")
    m = re.search(r'_RESULTS = Path\(__file__\)\.resolve\(\)\.(parent|parents\[(\d+)\]) / "results"',
                  t)
    assert m, f"{path.name}: _RESULTS must be __file__-relative, not CWD-relative"
    depth = int(m.group(2)) if m.group(2) else 0
    expected = len(path.relative_to(ROOT).parts) - 1      # 0 at root, 1 in a subdirectory
    assert depth == expected, (
        f"{path.name} lives {expected} level(s) below the repo root but resolves _RESULTS with "
        f"parents[{depth}] -- it would write to the wrong directory")


@pytest.mark.parametrize("path", _writers(), ids=lambda p: p.name)
def test_every_writer_still_parses(path: Path):
    """A syntax check over the scripts CI never imports, so a bad edit cannot sit unnoticed."""
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_no_results_json_is_left_in_the_repo_root():
    """The tidy-up's invariant. A new one appearing in root means a writer regressed."""
    stray = sorted(p.name for p in ROOT.glob("*_results.json"))
    assert not stray, f"these belong in results/: {stray}"


def test_the_results_directory_exists_and_holds_the_recorded_runs():
    assert (ROOT / "results").is_dir()
    assert len(list((ROOT / "results").glob("*_results.json"))) >= 15
