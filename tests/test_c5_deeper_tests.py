"""Unit tests for STAGE 1.5.3 step 5b — the D1..D7 diagnostics that chose C-5's option.

These decided a design choice, so they are tested against **closed forms and constructions with a
known answer**, not against themselves. Two of them (D6, D7) overturned readings taken from the
occupancy bar alone, which is exactly the class of number this project has been burned by before.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter
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


D = _load("c5_deeper_tests", "plan_tests/c5_deeper_tests.py")

N_TRAIN, N_AGE, BATCH, STEPS = 33_688, 75, 512, 65


def _epoch(**kw):
    kw.setdefault("weight", 1.0)
    kw.setdefault("accumulate", 1)
    kw.setdefault("replacement", False)
    kw.setdefault("rng", np.random.default_rng(0))
    return D.simulate_epoch(N_AGE, N_TRAIN, BATCH, STEPS, **kw)


# ------------------------------------------------------------ draw_probabilities ---- #
def test_draw_probabilities_is_a_distribution_with_the_requested_ratio():
    p = D.draw_probabilities(75, 1000, 7.0)
    assert p.sum() == pytest.approx(1.0)
    assert p[0] / p[-1] == pytest.approx(7.0)          # an age cell is 7x likelier than a non-age
    assert (p[:75] == p[0]).all() and (p[75:] == p[-1]).all()


def test_weight_one_is_the_uniform_distribution():
    p = D.draw_probabilities(75, 1000, 1.0)
    assert p == pytest.approx(np.full(1000, 1 / 1000))


# ------------------------------------------------------------------ D1: duplicates ---- #
def test_a_permutation_can_never_duplicate_but_a_bootstrap_can():
    """D1's premise. Without it, 'effective cells' and 'raw cells' would be the same column."""
    perm = _epoch()
    assert perm["raw_cells_per_update"] == pytest.approx(perm["effective_cells_per_update"])
    boot = _epoch(weight=7.1, replacement=True, accumulate=8,
                  rng=np.random.default_rng(1))
    assert boot["raw_cells_per_update"] > boot["effective_cells_per_update"]


# -------------------------------------------------------------- D2/D7: update count ---- #
def test_accumulation_divides_the_update_count_exactly():
    assert _epoch(accumulate=8)["age_updates_per_epoch"] == STEPS // 8
    assert _epoch(accumulate=1)["age_updates_per_epoch"] == STEPS


def test_the_status_quo_has_fewer_nonzero_updates_than_updates():
    """D7. `losses.py:55-57` returns a hard zero for an age-free batch, so counting updates
    overstates the status quo's age optimisation by the empty-batch rate (~32%)."""
    r = _epoch(rng=np.random.default_rng(2))
    assert r["nonzero_age_updates_per_epoch"] < r["age_updates_per_epoch"]
    assert r["nonzero_age_updates_per_epoch"] / STEPS == pytest.approx(1 - np.exp(-1.14), abs=0.06)


def test_accumulation_at_w8_makes_every_update_nonzero():
    """Option 2's actual claim: not 'more cells' but 'no wasted updates'."""
    r = _epoch(accumulate=8, rng=np.random.default_rng(3))
    assert r["nonzero_age_updates_per_epoch"] == r["age_updates_per_epoch"]


# --------------------------------------------------------------------- D3: coverage ---- #
def test_a_permutation_visits_each_label_at_most_once():
    """The definition of an epoch, and the baseline D6 is measured against."""
    r = _epoch(rng=np.random.default_rng(4))
    assert max(r["visits_per_label"].values()) == 1
    assert r["total_age_draws"] <= N_AGE


def test_oversampling_costs_coverage_when_it_is_also_a_bootstrap():
    """A weighted bootstrap can miss a label entirely; a permutation only misses the tail."""
    weak = _epoch(weight=3.0, replacement=True, rng=np.random.default_rng(5))
    assert weak["label_coverage"] < 1.0


# ----------------------------------------------------------------- D4: donor balance ---- #
def test_donor_balance_is_flat_when_every_label_is_visited_equally():
    donors = {"a": 10, "b": 20, "c": 5}
    visits = Counter({i: 4 for i in range(35)})
    bal = D.donor_balance(visits, donors)
    assert bal["max_over_min"] == pytest.approx(1.0)
    assert bal["worst_deviation"] == pytest.approx(0.0, abs=1e-12)


def test_donor_balance_detects_an_over_visited_donor():
    """Constructed so the answer is known: donor 'a' gets twice the per-label visits."""
    donors = {"a": 10, "b": 10}
    visits = Counter({**{i: 2 for i in range(10)}, **{i: 1 for i in range(10, 20)}})
    bal = D.donor_balance(visits, donors)
    assert bal["max_over_min"] == pytest.approx(2.0)
    assert bal["per_donor"]["a"]["ratio"] == pytest.approx(4 / 3)


def test_donor_balance_survives_a_donor_with_no_visits_at_all():
    """max/min must not raise or silently return inf garbage on a zero denominator."""
    bal = D.donor_balance(Counter({0: 5}), {"a": 1, "b": 1})
    assert np.isfinite(bal["max_over_min"])


# -------------------------------------------------------------------- D5: bootstrap ---- #
def test_bootstrap_loss_is_one_minus_the_seen_fraction():
    assert D.bootstrap_loss(100, 63) == pytest.approx(0.37)
    assert D.bootstrap_loss(100, 100) == pytest.approx(0.0)


def test_bootstrap_spread_matches_a_direct_simulation():
    """D5's honest framing is a closed form, so it is checked against the thing it models
    rather than trusted -- the per-epoch miss rate is re-rolled, and the run-level spread is
    what actually survives."""
    n_train, draws, epochs = 200, 200, 10
    got = D.bootstrap_visit_spread(n_train, draws, epochs)
    rng = np.random.default_rng(7)
    counts = [int((rng.integers(0, n_train, size=draws * epochs) == 0).sum()) for _ in range(3000)]
    assert float(np.mean(counts)) == pytest.approx(got["expected_visits_over_run"], rel=0.05)
    assert float(np.std(counts)) == pytest.approx(got["sd_visits_over_run"], rel=0.08)


def test_a_bootstrap_loses_nobody_over_a_long_run():
    """The claim that overturned D5's headline: 36% per epoch, ~0 over 60 epochs."""
    got = D.bootstrap_visit_spread(N_TRAIN, STEPS * BATCH, 60)
    assert got["p_never_seen_in_whole_run"] < 1e-20
    assert 0.05 < got["cv"] < 0.25


def test_one_epoch_of_bootstrap_does_lose_about_a_third():
    """...but the per-epoch cost is real, and must not be rounded away either."""
    got = D.bootstrap_visit_spread(N_TRAIN, STEPS * BATCH, 1)
    assert got["p_never_seen_in_whole_run"] == pytest.approx(np.exp(-1.0), abs=0.02)


# ------------------------------------------------------ D6: information vs repetition ---- #
def test_a_permutation_is_exactly_one_pass_per_epoch():
    """D6's zero point. 1.0 means 'an epoch', and is what the status quo and Option 2 both do."""
    r = _epoch(rng=np.random.default_rng(8))
    v = D.visits_per_label_per_epoch(r["total_age_draws"], N_AGE)
    assert v == pytest.approx(1.0, abs=0.05)
    assert v <= 1.0                       # a permutation cannot exceed one pass, ever


def test_a_sampler_weight_of_w_runs_about_w_passes_per_epoch():
    """The measurement that decided the option: weight `w` buys repetition, not information.

    Graded against the EXACT expectation `draws*w / (n_age*w + rest)` rather than against `w`
    itself -- the two agree only because `n_age*w << rest`, and asserting the approximation
    would be asserting the rounding. One epoch is too noisy to resolve either (sd ~ 0.2 on a
    mean of 2.95), so the mean over 20 epochs is used.
    """
    rng = np.random.default_rng(9)
    rest = N_TRAIN - N_AGE
    for w in (3.0, 7.1):
        vs = [D.visits_per_label_per_epoch(
            _epoch(weight=w, replacement=True, rng=rng)["total_age_draws"], N_AGE)
            for _ in range(20)]
        exact = (STEPS * BATCH) * w / (N_AGE * w + rest)
        assert float(np.mean(vs)) == pytest.approx(exact, rel=0.05)
        # and the interpretive claim the plan actually leans on: ~w passes, i.e. repetition
        assert float(np.mean(vs)) == pytest.approx(w, rel=0.15)


def test_accumulation_changes_delivery_and_not_exposure():
    """Option 2's defining property, stated as a test: same labels, same number of passes,
    only regrouped. If this ever fails, Option 2 is no longer the minimal intervention and the
    step 5b decision must be revisited."""
    plain = _epoch(rng=np.random.default_rng(10))
    accum = _epoch(accumulate=8, rng=np.random.default_rng(10))
    assert plain["total_age_draws"] == accum["total_age_draws"]
    assert (D.visits_per_label_per_epoch(accum["total_age_draws"], N_AGE)
            == pytest.approx(D.visits_per_label_per_epoch(plain["total_age_draws"], N_AGE)))
    assert accum["effective_cells_per_update"] > 7 * plain["effective_cells_per_update"]
