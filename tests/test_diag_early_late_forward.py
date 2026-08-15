"""Unit tests for the early->late forward diagnostic.

The transition rule is tested BEFORE it is trusted on real data, on trajectories where the answer
is known by construction -- especially the 'sustained' clause, which is the only thing separating
a real transition from the 40-50 yr swings the early single samples show.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "delf", ROOT / "experiments" / "diag_early_late_forward.py")
delf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(delf)

DAYS = [7.0, 9.0, 11.0, 13.0, 15.0, 21.0, 29.0, 34.0, 40.0, 47.0, 54.0]


# ---- transition_day (the pre-registered rule) --------------------------------------------- #
def test_a_clean_step_returns_the_first_low_day():
    ages = [100, 100, 100, 100, 100, 100, 100, 40, 40, 40, 40]
    assert delf.transition_day(DAYS, ages) == 34.0


def test_a_late_step_is_found_at_the_right_day():
    ages = [100, 100, 100, 100, 100, 100, 100, 100, 100, 40, 40]
    assert delf.transition_day(DAYS, ages) == 47.0


def test_a_transient_dip_that_recovers_is_not_the_transition():
    """The load-bearing clause. O1 swings d9=121 -> d11=77 with no transition; a first-crossing
    rule without 'sustained' would report day 11."""
    ages = [100, 130, 20, 110, 100, 100, 100, 40, 40, 40, 40]
    assert delf.transition_day(DAYS, ages) == 34.0


def test_a_trajectory_that_never_settles_low_is_undefined():
    ages = [100, 100, 100, 100, 100, 100, 100, 40, 40, 40, 200]
    assert np.isnan(delf.transition_day(DAYS, ages))


def test_a_flat_trajectory_is_undefined():
    """mid equals every value, so nothing is strictly below it -- must be nan, not day 7."""
    assert np.isnan(delf.transition_day(DAYS, [50.0] * len(DAYS)))


def test_missing_early_or_late_window_is_undefined():
    assert np.isnan(delf.transition_day([7.0, 9.0], [100.0, 90.0]))          # no late
    assert np.isnan(delf.transition_day([40.0, 47.0], [40.0, 38.0]))         # no early


def test_days_need_not_be_supplied_in_order():
    """Real per-donor series come from a groupby and must not depend on row order."""
    ages = [40, 100, 100, 40, 100, 100, 100, 40, 100, 40, 100]
    d = [34.0, 7.0, 9.0, 40.0, 11.0, 13.0, 15.0, 47.0, 21.0, 54.0, 29.0]
    assert delf.transition_day(d, ages) == 34.0


def test_the_rule_uses_the_donors_own_midpoint_not_a_fixed_threshold():
    """A donor whose whole trajectory sits low must still yield a transition, and one sitting
    high must too -- the rule is scale-free per donor."""
    low = [50, 50, 50, 50, 50, 50, 50, 20, 20, 20, 20]
    high = [500, 500, 500, 500, 500, 500, 500, 200, 200, 200, 200]
    assert delf.transition_day(DAYS, low) == 34.0
    assert delf.transition_day(DAYS, high) == 34.0


# ---- partial_corr -------------------------------------------------------------------------- #
def test_partial_equals_full_when_the_covariate_is_unrelated():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    y = x + rng.normal(0, 0.1, 200)
    z = rng.normal(size=200)
    assert delf.partial_corr(x, y, z) == pytest.approx(delf.pearson(x, y), abs=0.05)


def test_a_fully_mediating_covariate_drives_the_partial_to_zero():
    """The exact scenario under test: x and y correlate only because both track z."""
    rng = np.random.default_rng(1)
    z = rng.normal(size=400)
    x = z + rng.normal(0, 0.05, 400)
    y = z + rng.normal(0, 0.05, 400)
    assert abs(delf.partial_corr(x, y, z)) < 0.2
    assert delf.pearson(x, y) > 0.9          # the raw correlation is large and misleading


def test_partial_is_undefined_when_the_covariate_explains_a_variable_exactly():
    z = np.linspace(0, 1, 50)
    assert np.isnan(delf.partial_corr(z.copy(), np.sin(z), z))


def test_partial_is_undefined_on_degenerate_input():
    z = np.linspace(0, 1, 50)
    assert np.isnan(delf.partial_corr(np.ones(50), np.sin(z), z))


def test_pearson_guards_short_and_constant_input():
    assert np.isnan(delf.pearson(np.array([1.0, 2.0]), np.array([1.0, 2.0])))
    assert np.isnan(delf.pearson(np.ones(10), np.arange(10.0)))


# ---- the pre-registered constants ---------------------------------------------------------- #
def test_thresholds_are_stated_constants_not_fitted_after_the_fact():
    assert delf.SPEARMAN_CRIT_N6 == 0.886
    assert delf.PEARSON_CRIT_N6 == 0.811
    assert delf.PARTIAL_CRIT_DF3 == 0.878
    assert (delf.EARLY_LO, delf.EARLY_HI, delf.LATE_LO) == (7.0, 29.0, 34.0)
