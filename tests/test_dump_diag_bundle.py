"""`donor_ids_from_counts` and the pool reader.

The reconstruction is load-bearing: leave-one-donor-out WITHIN the pool is the only honest way to
compare calibrators without touching the graded folds, and it needs per-row donor labels the pool
does not store. Getting it silently wrong would mis-assign rows to donors and quietly invalidate
every comparison built on it -- so the guard that refuses to guess is tested as hard as the
happy path.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dump_diag_bundle import donor_ids_from_counts  # noqa: E402


def test_labels_follow_insertion_order_of_the_counts():
    ids = donor_ids_from_counts({1: 3, 2: 2, 5: 1}, 6)
    assert ids.tolist() == [1, 1, 1, 2, 2, 5]


def test_matches_run_3_geometry():
    counts = {1: 21, 2: 21, 3: 21, 4: 19, 5: 21}      # N2's actual pool
    ids = donor_ids_from_counts(counts, 103)
    assert len(ids) == 103
    assert [int((ids == d).sum()) for d in counts] == [21, 21, 21, 19, 21]
    assert len(np.unique(ids)) == 5


def test_refuses_when_totals_disagree():
    """The condition that rules out an am-mask mismatch. Must return None, not guess."""
    assert donor_ids_from_counts({1: 21, 2: 21}, 50) is None     # fate rows > residual rows
    assert donor_ids_from_counts({1: 21, 2: 21}, 21) is None     # and the other direction


def test_refuses_on_empty_counts():
    assert donor_ids_from_counts({}, 0) is None
    assert donor_ids_from_counts({}, 10) is None


def test_dtype_is_integral_so_it_can_index():
    ids = donor_ids_from_counts({7: 2, 9: 2}, 4)
    assert np.issubdtype(ids.dtype, np.integer)


def test_a_single_donor_pool_still_reconstructs():
    # degenerate but not invalid -- the caller decides whether 1 donor is usable, not this fn
    assert donor_ids_from_counts({3: 4}, 4).tolist() == [3, 3, 3, 3]


# ------------------------------------------------------------------ round trip ---- #
def test_pool_reader_round_trips_a_real_xstats(tmp_path):
    """Write a pool with the real save_xstats, read it back through _pool."""
    import dump_diag_bundle

    from cellfate.training.xdonor_calib import XDonorStats, save_xstats

    rng = np.random.default_rng(0)
    n = 12
    xs = XDonorStats(
        abs_residuals=rng.random(n),
        logits=rng.random((n, 3)),
        targets=np.eye(3)[rng.integers(0, 3, n)],
        probs_mean=rng.dirichlet(np.ones(3), n),
        sigma_pred=rng.random(n),
        sigma_pred_mc=rng.random(n),
        n_donors=3,
        feats=rng.random((n, 4)),
        residuals_per_donor={1: 4, 2: 4, 3: 4},
    )
    save_xstats(tmp_path / "bundle", xs)

    arrays, meta = dump_diag_bundle._pool(tmp_path)
    assert meta["pool_n"] == n
    assert meta["pool_n_donors"] == 3
    assert meta["pool_donor_ids_reconstructed"] is True
    assert arrays["pool_donor_id"].tolist() == [1] * 4 + [2] * 4 + [3] * 4
    np.testing.assert_allclose(arrays["pool_probs_mean"], xs.probs_mean)
    np.testing.assert_allclose(arrays["pool_abs_residuals"], xs.abs_residuals)


def test_pool_reader_flags_unreconstructible_labels_instead_of_crashing(tmp_path):
    """Counts describing fewer rows than the fate arrays: labels must be withheld, not invented."""
    import dump_diag_bundle

    from cellfate.training.xdonor_calib import XDonorStats, save_xstats

    n = 10
    xs = XDonorStats(
        abs_residuals=np.zeros(6), logits=np.zeros((n, 3)), targets=np.zeros((n, 3)),
        probs_mean=np.zeros((n, 3)), sigma_pred=np.zeros(6), sigma_pred_mc=np.zeros(6),
        n_donors=2, feats=np.zeros((n, 2)),
        residuals_per_donor={1: 3, 2: 3},          # 6 residual rows vs 10 fate rows
    )
    save_xstats(tmp_path / "bundle", xs)

    arrays, meta = dump_diag_bundle._pool(tmp_path)
    assert meta["pool_donor_ids_reconstructed"] is False
    assert "pool_donor_id" not in arrays
    assert "NOT reconstructible" in meta["pool_donor_ids_note"]
