"""`diag_calibrators`: the ECE estimator, the LODO harness, and the effective-n decomposition.

These decide which calibrator gets proposed next, so the properties that make their output
trustworthy are asserted rather than assumed -- above all that `lodo_scores` is genuinely
out-of-sample (a leak would make every family look good and the comparison worthless).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from diag_calibrators import (  # noqa: E402
    FAMILIES,
    _fit_identity,
    _fit_isotonic,
    _fit_logit_platt,
    ece,
    effective_n,
    lodo_scores,
)


# ---------------------------------------------------------------------- ece ----- #
def test_ece_matches_scorecard_on_a_worked_example():
    """Same binning as scorecard.py:112 -- if this drifts, we stop measuring the graded quantity."""
    p = np.array([0.95, 0.95, 0.95, 0.95])
    y = np.array([1.0, 1.0, 1.0, 0.0])           # one bin, conf 0.95, acc 0.75
    assert ece(p, y) == pytest.approx(0.20)


def test_ece_is_zero_for_a_perfectly_calibrated_bin():
    p = np.full(10, 0.5)
    y = np.array([1.0] * 5 + [0.0] * 5)
    assert ece(p, y) == pytest.approx(0.0)


def test_ece_weights_bins_by_occupancy():
    # 9 cells perfectly calibrated at 0.05, 1 cell wrong at 0.95 -> 0.1 weight on a 0.95 gap
    p = np.array([0.05] * 9 + [0.95])
    y = np.array([0.0] * 9 + [0.0])
    assert ece(p, y) == pytest.approx(0.1 * 0.95 + 0.9 * 0.05, abs=1e-9)


# --------------------------------------------------------------- lodo_scores ---- #
def test_lodo_is_actually_out_of_sample():
    """The load-bearing property. A memoriser must score badly out-of-fold and 0 in-sample.

    Isotonic on a perfectly separable-but-inconsistent pool: each donor's label is the OPPOSITE
    of what the other donors imply, so anything that leaked the held-out donor would look good.
    """
    donor = np.repeat([1, 2], 20)
    p = np.tile(np.linspace(0.01, 0.99, 20), 2)
    y = np.concatenate([(np.linspace(0.01, 0.99, 20) > 0.5).astype(float),
                        (np.linspace(0.01, 0.99, 20) < 0.5).astype(float)])
    oos, ins = lodo_scores(p, y, donor, _fit_isotonic)
    assert ins < 0.05, "isotonic should interpolate its own fitting data"
    assert oos > 0.3, "held-out donors contradict each other; a leak would hide that"
    assert oos > ins


def test_identity_scores_the_same_in_and_out_of_sample():
    """It has nothing to fit, so there is no optimism -- a sanity check on the harness itself."""
    rng = np.random.default_rng(0)
    donor = np.repeat([1, 2, 3], 15)
    p = rng.random(45)
    y = (rng.random(45) < p).astype(float)
    oos, ins = lodo_scores(p, y, donor, _fit_identity)
    assert oos == pytest.approx(ins)


def test_every_row_gets_an_out_of_fold_prediction():
    donor = np.repeat([1, 2, 3], 10)
    rng = np.random.default_rng(1)
    p = rng.uniform(0.1, 0.9, 30)
    y = (rng.random(30) < p).astype(float)
    for fit in FAMILIES.values():
        oos, _ = lodo_scores(p, y, donor, fit)
        assert np.isfinite(oos)


def test_single_class_donor_passes_through_instead_of_crashing():
    """A donor whose training remainder has one class cannot fit a calibrator; must not raise."""
    donor = np.repeat([1, 2], 10)
    p = np.linspace(0.1, 0.9, 20)
    y = np.ones(20)                               # no variation anywhere
    oos, ins = lodo_scores(p, y, donor, _fit_logit_platt)
    assert np.isfinite(oos) and np.isfinite(ins)


def test_a_miscalibrated_but_well_ranked_pool_is_fixable_by_calibration():
    """The regime run 3 is in: near-perfect ranking, badly offset probabilities.

    Calibration should beat identity out-of-sample here -- if it did not, the harness would be
    unable to detect the very thing it exists to measure.
    """
    rng = np.random.default_rng(3)
    donor = np.repeat([1, 2, 3, 4, 5], 30)
    z = rng.normal(0.8, 1.2, 150)
    y = (rng.random(150) < 1 / (1 + np.exp(-z))).astype(float)
    p = 1 / (1 + np.exp(-3.0 * z))                # same ranking, over-sharpened
    base, _ = lodo_scores(p, y, donor, _fit_identity)
    fixed, _ = lodo_scores(p, y, donor, _fit_logit_platt)
    assert fixed < base


# -------------------------------------------------------------- effective_n ----- #
def test_icc_near_one_when_donor_offsets_dominate():
    donor = np.repeat([1, 2, 3, 4, 5], 20)
    p = np.full(100, 0.5)
    y = np.repeat([0.0, 0.25, 0.5, 0.75, 1.0], 20)   # all variance is BETWEEN donors
    e = effective_n(p, y, donor)
    assert e["icc"] > 0.95
    assert e["n_eff"] < 10, "100 rows over 5 donors with no within-variance is ~5 points"


def test_icc_near_zero_when_error_is_independent_of_donor():
    rng = np.random.default_rng(5)
    donor = np.repeat([1, 2, 3, 4, 5], 40)
    p = np.full(200, 0.5)
    y = (rng.random(200) < 0.5).astype(float)        # no donor structure
    e = effective_n(p, y, donor)
    assert e["icc"] < 0.15
    assert e["n_eff"] > 100


def test_effective_n_reports_shape():
    donor = np.repeat([1, 2], 10)
    e = effective_n(np.full(20, 0.5), np.zeros(20), donor)
    assert e["donors"] == 2 and e["n"] == 20
