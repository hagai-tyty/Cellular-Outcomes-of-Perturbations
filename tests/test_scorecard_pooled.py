"""`scorecard.pooled_fate_ece` — the repaired form of the calibration target.

The per-fold `fate_ece` it sits beside is measured on ~21 cells in 10 bins, where a PERFECTLY
calibrated model scores 0.183 and clears the 0.169 bar only 26.9% of the time. Pooling raises
that to 99.6%. These tests pin the properties the repair depends on, above all that `excess`
does not reward a calibrator for merely sharpening.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import scorecard as sc  # noqa: E402


def _folds(pairs: dict):
    return {d: {"_fate_S": list(map(float, s)), "_fate_y": list(map(int, y))}
            for d, (s, y) in pairs.items()}


def test_old_snapshots_without_the_field_return_none_not_an_error():
    """baseline.json predates pooled scoring and must still load."""
    assert sc.pooled_fate_ece({"N2": {"fate_ece": 0.28}, "N3": {"fate_ece": 0.3}}) is None


def test_error_folds_are_skipped():
    f = _folds({"N2": (np.full(60, 0.5), np.ones(60))})
    f["N3"] = {"_error": "boom", "_fate_S": [0.9] * 60, "_fate_y": [0] * 60}
    r = sc.pooled_fate_ece(f, trials=200)
    assert r["n"] == 60


def test_too_few_cells_returns_none():
    assert sc.pooled_fate_ece(_folds({"N2": (np.full(5, 0.5), np.ones(5))}), trials=50) is None


def test_it_actually_pools_rather_than_averaging():
    """Two folds with OPPOSITE errors. Averaging per-fold ECEs adds their magnitudes; pooling
    lets them offset, which is the point -- the pooled number describes the METHOD across the
    donor population rather than the worst donor twice."""
    p = np.full(40, 0.3)
    f = _folds({"N2": (p, np.zeros(40)), "N3": (p, np.ones(40))})
    # per-fold: |0.3-0| = 0.3 and |0.3-1| = 0.7  ->  mean 0.500
    assert np.mean([sc._ece(p, np.zeros(40)), sc._ece(p, np.ones(40))]) == pytest.approx(0.5)
    r = sc.pooled_fate_ece(f, trials=200)
    # pooled: 80 cells at p=0.3 with accuracy 0.5  ->  |0.3-0.5| = 0.2
    assert r["n"] == 80
    assert r["ece"] == pytest.approx(0.2, abs=1e-9)


def test_a_well_calibrated_pool_scores_near_zero_and_a_biased_one_does_not():
    rng = np.random.default_rng(0)
    p = rng.uniform(0.05, 0.95, 4000)
    good = _folds({"N2": (p, (rng.random(len(p)) < p).astype(int))})
    bad = _folds({"N2": (p, (rng.random(len(p)) < np.clip(p + 0.25, 0, 1)).astype(int))})
    assert sc.pooled_fate_ece(good, trials=100)["ece"] < 0.02
    assert sc.pooled_fate_ece(bad, trials=100)["ece"] > 0.15


def test_floor_is_positive_and_shrinks_as_the_pool_grows():
    """The whole reason pooling repairs the criterion."""
    rng = np.random.default_rng(1)
    small = rng.uniform(0.1, 0.9, 21)
    big = rng.uniform(0.1, 0.9, 500)
    fs = sc.pooled_fate_ece(_folds({"N2": (small, (rng.random(21) < small).astype(int))}),
                            trials=800)["floor"]
    fb = sc.pooled_fate_ece(_folds({"N2": (big, (rng.random(500) < big).astype(int))}),
                            trials=800)["floor"]
    assert fs > fb > 0


def test_excess_does_not_reward_pure_sharpening():
    """THE property the repair exists for.

    Sharpening a perfectly calibrated model makes it strictly worse. Raw ECE can be misleading
    because the floor drops underneath it; `excess` must report the degradation.
    """
    rng = np.random.default_rng(2)
    p = rng.uniform(0.15, 0.85, 400)
    y = (rng.random(len(p)) < p).astype(int)
    z = np.log(p / (1 - p))
    sharp = 1 / (1 + np.exp(-3.0 * z))

    base = sc.pooled_fate_ece(_folds({"N2": (p, y)}), trials=600)
    worse = sc.pooled_fate_ece(_folds({"N2": (sharp, y)}), trials=600)

    assert worse["floor"] < base["floor"], "sharpening must lower the floor -- the confounder"
    assert worse["excess"] > base["excess"], "excess must still call the sharpened model worse"


def test_percentile_flags_a_genuinely_miscalibrated_pool():
    rng = np.random.default_rng(3)
    p = rng.uniform(0.1, 0.9, 300)
    perfect = _folds({"N2": (p, (rng.random(len(p)) < p).astype(int))})
    broken = _folds({"N2": (p, (rng.random(len(p)) < np.clip(p + 0.3, 0, 1)).astype(int))})
    assert sc.pooled_fate_ece(perfect, trials=600)["pctile"] < 0.95
    assert sc.pooled_fate_ece(broken, trials=600)["pctile"] > 0.99


def test_result_is_deterministic_for_a_given_seed():
    f = _folds({"N2": (np.linspace(0.05, 0.95, 100), (np.arange(100) % 2))})
    assert sc.pooled_fate_ece(f, trials=200) == sc.pooled_fate_ece(f, trials=200)
