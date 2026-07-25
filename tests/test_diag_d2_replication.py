"""STAGE 1.5 §10 D2 — every branch of the independent-replication diagnostic.

Per the `verify_1a` lesson: a branch that never executes is not a check. Nothing here touches data;
`diag_d2_replication` keeps all repo-data imports inside `timepoint_ages`.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "diag_d2_replication", _ROOT / "experiments" / "diag_d2_replication.py")
d2 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = d2
_SPEC.loader.exec_module(d2)


# ------------------------------- day_of_label ------------------------------ #
def test_day_of_label_parses_the_day_timepoints():
    assert d2.day_of_label("D0") == 0.0
    assert d2.day_of_label("D14") == 14.0
    assert d2.day_of_label("d8") == 8.0                     # case-insensitive


def test_day_of_label_returns_none_for_ipsc():
    """iPSC is an endpoint, not a day -- it must be excluded from the primary window."""
    assert d2.day_of_label("iPSC") is None
    assert d2.day_of_label("IPSC") is None


def test_day_of_label_returns_none_for_unrecognised_labels():
    assert d2.day_of_label("Fib") is None
    assert d2.day_of_label("") is None


# --------------------------------- spearman -------------------------------- #
def test_spearman_is_minus_one_for_a_clean_decrease():
    assert d2.spearman([0, 2, 4, 6], [40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_is_plus_one_for_a_clean_increase():
    assert d2.spearman([0, 2, 4, 6], [10, 20, 30, 40]) == pytest.approx(1.0)


def test_spearman_is_nan_without_enough_points_or_spread():
    assert np.isnan(d2.spearman([0, 2], [1, 2]))            # too few
    assert np.isnan(d2.spearman([0, 2, 4], [5, 5, 5]))      # no spread in y


def test_spearman_drops_non_finite_pairs():
    assert d2.spearman([0, 2, np.nan, 6], [10, 20, 99, 40]) == pytest.approx(1.0)


# -------------------------------- d2_verdict ------------------------------- #
def test_d2_replicates_when_age_rises_matching_gill():
    v = d2.d2_verdict(0.31, 8)
    assert v["status"] == "REPLICATES" and "not a Gill design artefact" in v["reason"]


def test_d2_replicates_even_for_a_weak_positive_because_the_bar_is_direction():
    """n=8 cannot resolve significance, so the pre-registered bar is direction alone."""
    assert d2.d2_verdict(0.05, 8)["status"] == "REPLICATES"


def test_d2_contradicts_when_age_clearly_falls():
    v = d2.d2_verdict(-0.55, 8)
    assert v["status"] == "CONTRADICTS" and "dissolves" in v["reason"]


def test_d2_ambiguous_between_the_contradiction_boundary_and_zero():
    assert d2.d2_verdict(-0.1, 8)["status"] == "AMBIGUOUS"
    assert d2.d2_verdict(0.0, 8)["status"] == "AMBIGUOUS"   # exactly zero is not a replication


def test_d2_contradiction_boundary_is_inclusive_at_the_pre_committed_value():
    assert d2.d2_verdict(d2.CONTRADICTS_AT, 8)["status"] == "CONTRADICTS"
    assert d2.d2_verdict(d2.CONTRADICTS_AT + 0.01, 8)["status"] == "AMBIGUOUS"


def test_d2_cannot_verify_on_too_few_points_or_undefined_rho():
    assert d2.d2_verdict(0.5, 2)["status"] == "CANNOT_VERIFY"
    assert d2.d2_verdict(float("nan"), 8)["status"] == "CANNOT_VERIFY"


# ----------------------------------- bars ---------------------------------- #
def test_bars_registers_d2_and_states_the_power_limitation():
    b = {x["id"]: x for x in d2.bars()}
    assert set(b) == {"D2"}
    assert "UNDERPOWERED" in b["D2"]["power_note"]


def test_pre_committed_constants_are_what_the_notebook_registered():
    assert d2.MAX_CELLS_PER_TIMEPOINT == 2000
    assert d2.N_PSEUDO_REPLICATES == 5
    assert d2.CONTRADICTS_AT == -0.2
    assert d2.GILL_E1B_RHO == pytest.approx(0.205)
