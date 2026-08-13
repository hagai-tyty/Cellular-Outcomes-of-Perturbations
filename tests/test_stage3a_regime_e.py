"""Unit tests for regime E — pure functions only, no repo data and no GEO files.

The one that matters most is `parse_sample_name`. The first implementation matched only the
treated samples' `{donor}_{arm}_{N}days_{exp}` pattern, and GSE165177 names its DAY-0 fibroblasts
`O1 Fib` — space-separated, no `days` token. All three were dropped **silently**, which cut every
trajectory's first timepoint and changed the pair count from 10 per donor to 6. Nothing failed;
the run simply answered a smaller question. That is the same class of invisible filter that has
cost this project real time before, so the day-0 branch is pinned here.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


E = _load("stage3a_regime_e", "experiments/stage3a_regime_e.py")


# ------------------------------------------------------------------ parse_sample_name ---- #
def test_parses_a_treated_sample():
    got = E.parse_sample_name("O1_negative_control_intermediate_13days_exp2")
    assert got == {"sample": "O1_negative_control_intermediate_13days_exp2", "donor": "O1",
                   "arm": "negative_control_intermediate", "day": 13.0, "exp": "exp2"}


def test_parses_the_day_zero_fibroblast_that_the_first_implementation_dropped():
    """`O1 Fib` — the regression this file exists for."""
    got = E.parse_sample_name("O1 Fib")
    assert got is not None, "day-0 fibroblasts must not be dropped silently"
    assert got["donor"] == "O1"
    assert got["day"] == 0.0
    assert got["arm"] == "day0_fibroblast"


@pytest.mark.parametrize("col", ["O1 Fib", "O2 Fib", "O3 Fib"])
def test_every_donor_day_zero_sample_is_recovered(col):
    got = E.parse_sample_name(col)
    assert got is not None and got["day"] == 0.0


def test_day_zero_fibroblast_is_not_marked_as_a_control_arm():
    """It is the trajectory's starting state. Marking it control deletes the only pre-transition
    timepoint, and the 33 contemporaneous `negative_control` samples are the better reference."""
    assert "negative_control" not in E.parse_sample_name("O1 Fib")["arm"]


@pytest.mark.parametrize("col", ["iPSC 13", "iPSC 21"])
def test_ipsc_lines_are_excluded(col):
    """No donor attribution — the series files them under donor `iPSC`, day 51. Not a point on
    any donor's trajectory."""
    assert E.parse_sample_name(col) is None


@pytest.mark.parametrize("col", ["Probe", "Chromosome", "Feature Strand", "Distance", "ID", ""])
def test_annotation_columns_are_rejected(col):
    assert E.parse_sample_name(col) is None


def test_multi_word_arm_names_survive_intact():
    got = E.parse_sample_name("O3_failing_to_transiently_reprogram_intermediate_10days_exp2")
    assert got["arm"] == "failing_to_transiently_reprogram_intermediate"
    assert got["day"] == 10.0 and got["donor"] == "O3"


# -------------------------------------------------------------------------- pairs_of ---- #
def _rows(days, ns=None, us=None):
    ns = ns or [4] * len(days)
    us = us or [0.5] * len(days)
    return [{"day": float(d), "x": np.array([float(d)]), "n": int(n), "u": float(u)}
            for d, n, u in zip(days, ns, us, strict=True)]


def test_pairs_are_strictly_forward_in_time():
    for p in E.pairs_of(_rows([0, 10, 13, 15, 17])):
        assert p["dt"] > 0


def test_the_five_gse165177_timepoints_give_ten_ordered_pairs():
    """The count the pre-registration claimed; it holds only once day 0 is parsed."""
    assert len(E.pairs_of(_rows([0, 10, 13, 15, 17]))) == 10


def test_dropping_day_zero_silently_costs_four_pairs():
    """Pins the size of the bug: 10 -> 6, with nothing raising."""
    assert len(E.pairs_of(_rows([10, 13, 15, 17]))) == 6


def test_pairs_carry_the_endpoint_sample_count_for_the_binomial_draw():
    ps = E.pairs_of(_rows([0, 10], ns=[1, 6]))
    assert len(ps) == 1 and ps[0]["n_j"] == 6


def test_day_zero_is_never_an_endpoint():
    """It is the minimum, so its n=1 never enters a target."""
    assert all(p["day_j"] > 0 for p in E.pairs_of(_rows([0, 10, 13])))


# ------------------------------------------------------------------------------ p_of ---- #
def test_p_of_reproduces_the_hff_curve_at_alpha_one_up_to_the_clip():
    """alpha=1 is HFF's own curve, except that it is clipped into (0.01, 0.99)."""
    assert E.p_of(E.HFF_DAYS, 1.0) == pytest.approx(np.clip(E.HFF_CURVE, 0.01, 0.99), abs=1e-9)


def test_p_of_clips_the_terminal_point_so_the_logit_stays_finite():
    """HFF's day-21 value is 0.9996; an unclipped p would make logit(p) blow up and the
    binomial draw degenerate. The clip is deliberate, so it is pinned rather than tolerated."""
    assert E.HFF_CURVE[-1] > 0.99
    assert E.p_of([21.0], 1.0)[0] == pytest.approx(0.99)


def test_p_of_is_flat_at_alpha_zero():
    v = E.p_of([0, 10, 13, 15, 17], 0.0)
    assert v.min() == pytest.approx(v.max())


def test_p_of_scales_the_amplitude_about_the_curve_mean():
    mid = E.HFF_CURVE.mean()
    day = 14.0                      # inside the curve, so the clip does not confound the check
    full = E.p_of([day], 1.0)[0]
    half = E.p_of([day], 0.5)[0]
    assert half - mid == pytest.approx((full - mid) / 2, rel=1e-9)


def test_p_of_stays_inside_the_unit_interval():
    v = E.p_of(np.linspace(0, 21, 50), 1.0)
    assert v.min() >= 0.0 and v.max() <= 1.0


def test_p_of_interpolates_between_hff_timepoints():
    """GSE165177's days 13/15/17 fall between HFF's 12, 14 and 21 and must interpolate."""
    assert E.HFF_CURVE[6] <= E.p_of([13.0], 1.0)[0] <= E.HFF_CURVE[7]
