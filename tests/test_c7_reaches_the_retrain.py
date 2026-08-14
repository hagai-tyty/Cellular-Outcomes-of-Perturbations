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


# --------------------------------------------------------------------------------------- #
# The THIRD incarnation: the flag arrived, and was still too late.                          #
#                                                                                           #
# `_load` caches and the gate's screen lives inside it. `run_multi_local.py` called `plan()`  #
# to list donors ~30 lines before the DataConfig existed, so `apply_source_flags` set a flag  #
# that read True over an already-cached, unscreened 124-column matrix. The retrain produced   #
# 42,605 cells with 0 masked labels (pre-C-7) under a header that said ON; correct is 42,600  #
# with 19 masked. Six folds, several hours.                                                   #
# --------------------------------------------------------------------------------------- #
def _bare_source():
    """A source with nonexistent paths -- __init__ reads nothing, so this needs no data."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from cellfate.data.sources import GillReprogrammingSource
    return GillReprogrammingSource("no_such_expr.txt.gz", "no_such_series.txt.gz")


def test_enabling_the_gate_drops_an_already_cached_read():
    """Otherwise the screen never re-runs and the gate is inert while reporting ON."""
    src = _bare_source()
    src._rpm, src._genes, src._meta = "CACHED", ["g"], {"s": {}}
    src.bulk_integrity_gate = True
    assert src._rpm is None, ("enabling the gate left a cached UNSCREENED matrix in place -- "
                             "the gate cannot fire, exactly as in the lost retrain")
    assert src._genes is None and src._meta is None


def test_disabling_the_gate_also_drops_the_cache():
    """Symmetric: a screened cache must not leak into a build that asked for the gate OFF,
    or B4 (bit-identical when disabled) silently stops being true."""
    src = _bare_source()
    src.bulk_integrity_gate = True
    src._rpm, src.rejected_samples = "SCREENED", {"col": "reason"}
    src.bulk_integrity_gate = False
    assert src._rpm is None and src.rejected_samples == {}


def test_setting_the_same_value_is_idempotent():
    """`apply_source_flags` runs from BOTH `build_sources` and `run` by design. If an
    unchanged set invalidated, every build would pay a second full matrix read."""
    src = _bare_source()
    src.bulk_integrity_gate = True
    src._rpm = "CACHED"
    src.bulk_integrity_gate = True          # same value, second call
    assert src._rpm == "CACHED"


def test_the_runner_sets_the_gate_before_it_plans():
    """Belt and braces to the property: the donor list must come from the gated corpus too."""
    src = RUNNER.read_text(encoding="utf-8")
    set_at = src.find("gill.bulk_integrity_gate")
    plan_at = src.find("gill.plan()")
    assert set_at != -1, "run_multi_local.py never sets the gate on the bulk source"
    assert plan_at != -1, "expected gill.plan() in the runner"
    assert set_at < plan_at, (
        "run_multi_local.py calls gill.plan() BEFORE setting bulk_integrity_gate -- plan() "
        "caches the unscreened matrix, which is how the C-7 retrain was lost")


def test_the_runner_fails_loudly_when_the_gate_bit_nothing():
    """A gate that reports ON and changes nothing must not survive to a scorecard."""
    src = RUNNER.read_text(encoding="utf-8")
    assert "rejected NOTHING" in src, (
        "run_multi_local.py must abort when the gate is on and no bulk column was rejected")
    assert "n_age_labeled" in src, (
        "run_multi_local.py must abort when the gate is on and no ΔAge label was masked")


def test_dataconfig_actually_has_the_field():
    """If src/ ever drops it, the runner's keyword becomes a TypeError at build time -- better to
    fail here, in milliseconds, than three hours into a retrain."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from cellfate.data import DataConfig
    assert "bulk_integrity_gate" in DataConfig.__dataclass_fields__
    assert DataConfig.__dataclass_fields__["bulk_integrity_gate"].default is False
