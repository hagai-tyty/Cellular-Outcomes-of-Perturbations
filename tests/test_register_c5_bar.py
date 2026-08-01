"""Unit tests for STAGE 1.5.3 step 5's C-5 bar — pure functions only, no repo data.

The occupancy model is what the whole C-5 decision rests on, so it is tested against closed-form
values rather than only against itself: a binomial mean is `n*p`, and the two knobs must reproduce
the three candidate mechanisms exactly.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


C5 = _load("register_c5_bar", "plan_tests/register_c5_bar.py")

N_TRAIN, N_AGE, BATCH = 33_688, 75, 512


# ------------------------------------------------------------ age_cells_per_update ---- #
def test_uniform_sampling_reproduces_the_plans_1_14_figure():
    """C-5's headline arithmetic: 75 * 512 / 33688 = 1.14 age cells in an average batch."""
    sim = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, n_sim=40_000,
                                  rng=np.random.default_rng(0))
    assert sim.mean() == pytest.approx(N_AGE * BATCH / N_TRAIN, rel=0.05)
    assert sim.mean() == pytest.approx(1.14, abs=0.08)


def test_the_empty_batch_rate_matches_the_poisson_estimate():
    """The plan estimates ~e^-1.14 = 32% of batches contribute nothing. The exact binomial
    answer is (1-p)^512; both should agree to within a point or so."""
    sim = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, n_sim=40_000,
                                  rng=np.random.default_rng(1))
    p = N_AGE / N_TRAIN
    assert float((sim == 0).mean()) == pytest.approx((1 - p) ** BATCH, abs=0.02)
    assert float((sim == 0).mean()) == pytest.approx(np.exp(-1.14), abs=0.03)


def test_a_sampler_weight_raises_occupancy_monotonically():
    rng = np.random.default_rng(2)
    means = [C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, weight=w, n_sim=8000,
                                     rng=rng).mean() for w in (1.0, 2.0, 4.0, 8.0)]
    assert means == sorted(means)
    assert means[-1] > 4 * means[0]


def test_accumulation_scales_the_cells_per_update_linearly():
    """Option 2's whole mechanism: one update over W batches sees W times the draws."""
    rng = np.random.default_rng(3)
    one = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, accumulate=1, n_sim=20_000, rng=rng)
    eight = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, accumulate=8, n_sim=20_000, rng=rng)
    assert eight.mean() == pytest.approx(8 * one.mean(), rel=0.08)


def test_option_3_is_exactly_the_status_quo_by_construction():
    """Pinning `s_age` touches neither sampling nor accumulation, so its occupancy is
    identical -- which is precisely C-5's criticism, encoded rather than argued."""
    rng = np.random.default_rng(4)
    a = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, weight=1.0, accumulate=1,
                                n_sim=5000, rng=np.random.default_rng(9))
    b = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, n_sim=5000, rng=np.random.default_rng(9))
    assert np.array_equal(a, b)
    assert rng is not None      # (the fixed seeds above are the point of the assertion)


# ---------------------------------------------------------------- required_weight ---- #
def test_required_weight_hits_its_target_when_simulated():
    """Solve for the weight, then MEASURE that it lands -- closed form checked against the
    simulation rather than trusted."""
    w = C5.required_weight(N_AGE, N_TRAIN, BATCH, target_cells=8.0)
    sim = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, weight=w, n_sim=40_000,
                                  rng=np.random.default_rng(5))
    assert sim.mean() == pytest.approx(8.0, rel=0.05)


def test_a_weight_of_one_is_the_status_quo_target():
    """Sanity floor on the algebra: w=1 must correspond to the uniform mean."""
    target = N_AGE * BATCH / N_TRAIN
    assert C5.required_weight(N_AGE, N_TRAIN, BATCH, target_cells=target) == pytest.approx(1.0,
                                                                                           rel=1e-6)


def test_an_unreachable_target_is_infinite_rather_than_silently_clipped():
    """Asking for more age cells than the batch holds has no answer; returning a finite
    number would invent one."""
    assert C5.required_weight(N_AGE, N_TRAIN, BATCH, target_cells=BATCH) == float("inf")


# --------------------------------------------------------------- oversampling_cost ---- #
def test_the_fate_cost_is_reported_as_a_fold_change_and_a_share():
    """Option 1's price. The plan flags it as 'not free'; this is the number."""
    w = C5.required_weight(N_AGE, N_TRAIN, BATCH, target_cells=8.0)
    cost = C5.oversampling_cost(N_AGE, N_TRAIN, w)
    assert cost["fold_oversampled"] == pytest.approx(7.0, abs=0.5)
    assert 0.0 < cost["extra_share_of_batch"] < 0.05
    assert cost["age_share_after"] > cost["age_share_before"]


def test_weight_one_costs_the_fate_task_nothing():
    cost = C5.oversampling_cost(N_AGE, N_TRAIN, 1.0)
    assert cost["fold_oversampled"] == pytest.approx(1.0)
    assert cost["extra_share_of_batch"] == pytest.approx(0.0, abs=1e-12)


# ------------------------------------------------------------------- discriminates ---- #
def test_a_bar_everything_passes_does_not_discriminate():
    same = {"a": {"B1": {"verdict": "RESOLVABLE"}, "B2": {"verdict": "RESOLVABLE"}},
            "b": {"B1": {"verdict": "RESOLVABLE"}, "B2": {"verdict": "RESOLVABLE"}}}
    assert C5.discriminates(same) is False


def test_a_bar_that_separates_candidates_does():
    mixed = {"a": {"B1": {"verdict": "RESOLVABLE"}, "B2": {"verdict": "RESOLVABLE"}},
             "b": {"B1": {"verdict": "UNRESOLVABLE"}, "B2": {"verdict": "UNRESOLVABLE"}}}
    assert C5.discriminates(mixed) is True
