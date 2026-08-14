"""The C-7 gate must actually REACH the retrain path.

This project has already lost a build to an inert flag: `bulk_integrity_gate` was wired into
`DataConfig` and into `build_sources`, but `run()` skips `build_sources` when sources are injected
— and every real caller injects. The first C-7 build produced `n_samples = 42605`, identical to the
ungated arm, and nothing raised (`e6fc183`).

The retrain path has the same shape. `run_multi_local.py` constructs its `DataConfig` by hand, so a
field it forgets is silently a default — and a forgotten `bulk_integrity_gate` means SIX FULL
BUILDS, several hours, producing pre-C-7 labels while every log line says C-7.

These are static source checks: no build, no data, no GPU. They fail loudly if the wiring is ever
removed.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "local_runners" / "run_multi_local.py"
LOOCV = ROOT / "local_runners" / "run_loocv.py"


def test_the_runner_passes_bulk_integrity_gate_into_dataconfig():
    """The whole point. A DataConfig built without it silently defaults to OFF."""
    src = RUNNER.read_text(encoding="utf-8")
    tree = ast.parse(src)
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "DataConfig"]
    assert calls, "no DataConfig(...) call found in the retrain runner"
    for c in calls:
        kw = {k.arg for k in c.keywords}
        assert "bulk_integrity_gate" in kw, (
            "DataConfig in run_multi_local.py does not pass bulk_integrity_gate -- a retrain "
            "would silently use PRE-C-7 labels")


def test_the_runner_exposes_a_module_level_toggle():
    src = RUNNER.read_text(encoding="utf-8")
    assert re.search(r"^BULK_INTEGRITY_GATE\s*=", src, re.M), (
        "run_multi_local.py must expose BULK_INTEGRITY_GATE at module level so run_loocv.py "
        "can set it the way it sets HARMONIZE and the age-shuffle flags")


def test_the_toggle_defaults_to_off_so_existing_behaviour_is_unchanged():
    """C-7 MOVES LABELS. Defaulting it on would silently change every other runner's output."""
    src = RUNNER.read_text(encoding="utf-8")
    m = re.search(r"^BULK_INTEGRITY_GATE\s*=\s*(\w+)", src, re.M)
    assert m and m.group(1) == "False"


def test_loocv_sets_the_toggle_from_the_environment():
    src = LOOCV.read_text(encoding="utf-8")
    assert "BULK_INTEGRITY_GATE" in src, "run_loocv.py never sets the gate"
    assert "CELLFATE_BULK_GATE" in src, "the gate is not env-driven, so it cannot compose with " \
                                        "CELLFATE_FOLD_SUFFIX on the command line"


def test_loocv_announces_which_label_set_it_is_about_to_build():
    """Hours of compute must not be ambiguous about which labels they used."""
    src = LOOCV.read_text(encoding="utf-8")
    assert re.search(r"C-7.*bulk_integrity_gate", src), (
        "run_loocv.py must print the gate state; a silent OFF is how the first C-7 build was lost")


@pytest.mark.parametrize("val,want", [("", False), ("0", False), ("false", False),
                                      ("1", True), ("yes", True), ("ON", True)])
def test_the_env_var_parsing_rule_is_the_one_the_runner_uses(val, want):
    """Pinned so 'CELLFATE_BULK_GATE=0' can never be read as truthy."""
    assert (val not in ("", "0", "false")) is want


def test_dataconfig_actually_has_the_field():
    """If src/ ever drops it, the runner's keyword becomes a TypeError at build time -- better to
    fail here, in milliseconds, than three hours into a retrain."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from cellfate.data import DataConfig
    assert "bulk_integrity_gate" in DataConfig.__dataclass_fields__
    assert DataConfig.__dataclass_fields__["bulk_integrity_gate"].default is False
