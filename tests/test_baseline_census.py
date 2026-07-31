"""STAGE 1.5.2 gates G-a and G-b — the baseline census and donor-age wiring.

Both gates are **record-only**. The single most important test in this file is
`test_delta_age_is_bit_identical_with_and_without_the_census`: the plan's hard guard is that
"ΔAge values must come out **bit-identical** before/after. It records, it does not compute. If any
ΔAge moves, the change is wrong — revert, do not rationalise." That is asserted here, not assumed
from reading the diff.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from verify_stage1_5 import ChunkControlStat, decide_verdict  # noqa: E402

from cellfate.data.aging import (  # noqa: E402
    LinearClock,
    _control_baseline,
    census_warnings,
    delta_age,
)
from cellfate.data.sources import ReprogrammingSource, _maybe_float  # noqa: E402


def _obs(lines, ctrl, **extra):
    d = {"cell_line": list(lines), "is_control": list(ctrl)}
    d.update({k: list(v) for k, v in extra.items()})
    return pd.DataFrame(d)


# ------------------------------------------------------------------ THE HARD GUARD ---- #
def test_delta_age_is_bit_identical_with_and_without_the_census():
    """G-a records; it must not compute. Any drift here means the change is wrong."""
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(8)]
    clock = LinearClock({g: float(w) for g, w in zip(genes, rng.normal(size=8), strict=True)},
                        intercept=41.5)
    expr = rng.normal(3.0, 1.0, size=(12, 8))
    obs = _obs(["A"] * 6 + ["B"] * 6, [True, False, False, True, False, False] * 2,
               batch=["E1", "E2"] * 6, donor_age=[53.0] * 12)
    d_plain, m_plain = delta_age(clock, expr, genes, obs, source="reprogramming")
    census: dict = {}
    d_census, m_census = delta_age(clock, expr, genes, obs, source="reprogramming", census=census)
    assert np.array_equal(d_plain, d_census)          # bit-identical, not "close"
    assert np.array_equal(m_plain, m_census)
    assert census, "the census must actually have been filled, or the test proves nothing"


def test_control_baseline_values_do_not_depend_on_the_composition_argument():
    values = np.array([10.0, 20.0, 30.0, 40.0])
    lines = np.array(["A", "A", "B", "B"])
    ctrl = np.array([True, False, True, False])
    a = _control_baseline(values, lines, ctrl)
    b = _control_baseline(values, lines, ctrl, census={},
                          composition={"batch": np.array(["x", "y", "x", "y"])})
    assert np.array_equal(a, b)


# ---------------------------------------------------------------------- the census ---- #
def test_census_records_the_unreplicated_case_that_was_silent_before():
    values = np.array([10.0, 20.0, 30.0])
    census: dict = {}
    _control_baseline(values, np.array(["A"] * 3), np.array([True, False, False]), census=census)
    assert census["A"]["n_control"] == 1
    assert census["A"]["unreplicated"] is True
    assert census["A"]["source"] == "controls"


def test_census_records_the_no_control_fallback():
    census: dict = {}
    _control_baseline(np.array([1.0, 2.0]), np.array(["A", "A"]), np.array([False, False]),
                      census=census)
    assert census["A"]["source"] == "self_fallback"
    assert census["A"]["n_control"] == 0


def test_census_records_only_the_batches_the_baseline_itself_came_from():
    """Finding D1 in one assertion: controls all Exp2 while the line spans Exp1 and Exp2."""
    census: dict = {}
    _control_baseline(
        np.array([1.0, 2.0, 3.0, 4.0]), np.array(["A"] * 4), np.array([False, True, False, True]),
        census=census, composition={"batch": np.array(["Exp1", "Exp2", "Exp1", "Exp2"])})
    assert census["A"]["batch"] == ["Exp2"]
    assert census["A"]["batch_in_line"] == ["Exp1", "Exp2"]
    assert census["A"]["n_cells"] == 4 and census["A"]["n_control"] == 2


def test_census_warnings_name_each_problem_exactly_once():
    census = {
        "ok": {"n_control": 5, "n_cells": 20, "source": "controls", "unreplicated": False},
        "solo": {"n_control": 1, "n_cells": 20, "source": "controls", "unreplicated": True},
        "none": {"n_control": 0, "n_cells": 20, "source": "self_fallback", "unreplicated": False},
    }
    w = " | ".join(census_warnings(census))
    assert "solo" in w and "n=1" in w
    assert "none" in w and "NO controls" in w
    assert "ok:" not in w


def test_census_warning_flags_a_single_batch_baseline_under_a_multi_batch_line():
    census = {"A": {"n_control": 2, "n_cells": 20, "source": "controls", "unreplicated": False,
                    "batch": ["Exp2"], "batch_in_line": ["Exp1", "Exp2"]}}
    assert any("cross-batch" in x for x in census_warnings(census))


def test_a_column_that_is_constant_within_a_line_never_warns():
    """`donor_age` is a per-donor constant, so a single baseline value is the ONLY possible
    answer. An earlier version warned on every donor -- noise that trains the reader to
    ignore the warnings that matter."""
    census = {"A": {"n_control": 2, "n_cells": 20, "source": "controls", "unreplicated": False,
                    "donor_age": ["53.0"], "donor_age_in_line": ["53.0"],
                    "batch": ["Exp2"], "batch_in_line": ["Exp1", "Exp2"]}}
    w = census_warnings(census)
    assert any("cross-batch" in x for x in w)
    assert not any("donor_age" in x for x in w)


# --------------------------------------------------- verify_stage1_5's extended census ---- #
def test_the_stage_1_5_pass_rule_is_unchanged_by_the_new_flags():
    """Four runs are recorded against what a Stage 1.5 PASS means. It must not move."""
    stats = [ChunkControlStat("c1", "A", 20, 1, control_batches=("E2",),
                              cell_batches=("E1", "E2")),
             ChunkControlStat("c2", "B", 20, 5)]
    v = decide_verdict(stats)
    assert v["status"] == "PASS"                      # still PASS: no fallback fired
    assert v["baseline_warnings"]                     # but the problems are now visible
    assert v["unreplicated_chunks"][0]["cell_line"] == "A"
    assert v["cross_batch_chunks"][0]["control_batches"] == ["E2"]


def test_a_fallback_chunk_still_fails():
    v = decide_verdict([ChunkControlStat("c1", "A", 20, 0)])
    assert v["status"] == "FAIL"


def test_no_batch_information_is_not_reported_as_verified_single_batch():
    """A source that does not stamp `batch` must not be silently cleared of D1."""
    s = ChunkControlStat("c1", "A", 20, 5)            # no batch info at all
    assert s.cross_batch_baseline is False
    assert decide_verdict([s])["cross_batch_chunks"] == []


# ------------------------------------------------------------------------ G-b wiring ---- #
@pytest.mark.parametrize("raw,expect", [("53", 53.0), (" 0 ", 0.0), ("35.0", 35.0),
                                        ("", None), ("N/A", None), (None, None)])
def test_maybe_float_never_defaults_a_missing_age_to_zero(raw, expect):
    """0.0 is a REAL age here (N2/N3 are neonatal), so a missing value must read as None."""
    assert _maybe_float(raw) == expect


def test_build_chunk_carries_extra_metadata_into_obs():
    raw = ReprogrammingSource.build_chunk(
        "t:A", np.zeros((3, 2)), ["G0", "G1"], "A",
        ["control", "OSKM", "OSKM"], [0.0, 24.0, 48.0],
        extra={"donor_age": [53.0] * 3, "batch": ["Exp2", "Exp1", "Exp2"]})
    assert list(raw.obs["donor_age"]) == [53.0] * 3
    assert list(raw.obs["batch"]) == ["Exp2", "Exp1", "Exp2"]


def test_build_chunk_rejects_a_mis_sized_extra_column():
    """Silently broadcasting or truncating would put the wrong age on a cell."""
    with pytest.raises(Exception, match="extra column"):
        ReprogrammingSource.build_chunk(
            "t:A", np.zeros((3, 2)), ["G0", "G1"], "A",
            ["control", "OSKM", "OSKM"], [0.0, 24.0, 48.0], extra={"donor_age": [53.0]})


def test_donor_age_and_batch_are_metadata_not_model_input():
    """They must never reach the model: the deployed request schema forbids extra fields."""
    from pydantic import ValidationError

    from cellfate.inference.schema import Request
    with pytest.raises(ValidationError, match="donor_age"):
        Request(X_raw=[0.0], u_modality="tf", u_descriptor="OSKM", dose_uM=1.0, time_h=24.0,
                donor_age=53.0)
