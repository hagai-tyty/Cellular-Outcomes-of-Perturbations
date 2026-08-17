"""Stage 13 -- the scorecard's verdicts must be computed on the quantity they claim to judge.

Two defects, both of which INVERTED decisions rather than merely degrading a number:

A. `"abs"` metrics were judged signed. `abs()` was applied once, to the aggregate, at display
   time -- so the column measured how far the donor panel CANCELS (0.230) instead of how large
   the error is (12.72 yr), the paired CI was built on signed differences, and `_verdict`'s
   better_is_down then read `-28 -> -22` as an increase. Both level-shift rows of the C-7
   comparison printed REGRESSION when the correct verdict is noise.

B. The two columns could average different fold sets, because each snapshot was aggregated over
   whatever folds were valid in ITSELF while the paired test used the intersection. 13 of 18
   rows printed a 6-fold mean beside a 5-fold one, and a verdict computed from neither.

The real snapshots are committed, so the defects are pinned against actual data, not only
synthetic constructions.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import scorecard as sc  # noqa: E402

GC2 = ROOT / "scorecard" / "gc2_A_keep_hff.json"
C7 = ROOT / "scorecard" / "c7_A_keep_hff.json"


def _real():
    return (json.loads(GC2.read_text(encoding="utf-8"))["folds"],
            json.loads(C7.read_text(encoding="utf-8"))["folds"])


def _folds(vals: dict, key="level_shift_model"):
    return {d: {key: v, "n_cells": 20} for d, v in vals.items()}


# ---- A1: the aggregate must accumulate magnitude, not cancel ------------------------------- #
def test_opposite_shifts_do_not_cancel_to_zero():
    """THE defect. Two donors shifted +10 and -10 are both badly wrong; the panel is not fine."""
    f = _folds({"N2": 10.0, "N3": -10.0})
    assert sc._agg(f, "level_shift_model") == pytest.approx(0.0)          # the old, signed path
    assert sc._agg(f, "level_shift_model", magnitude=True) == pytest.approx(10.0)


def test_magnitude_and_signed_agree_when_every_fold_has_the_same_sign():
    f = _folds({"N2": 4.0, "N3": 6.0})
    assert sc._agg(f, "level_shift_model", magnitude=True) == pytest.approx(5.0)
    assert sc._agg(f, "level_shift_model") == pytest.approx(5.0)


def test_on_real_data_the_ridge_level_shift_was_understated_55x():
    """`gc2_A` printed 0.230. The mean magnitude is 12.72 yr -- and `MASTER_PLAN.md` records the
    per-donor level shift, derived independently in Test 7.4.3, as +-12.7 yr. The corrected
    statistic reproduces the project's own number; the printed one erased it."""
    A, _ = _real()
    signed = sc._agg(A, "level_shift_ridge")
    mag = sc._agg(A, "level_shift_ridge", magnitude=True)
    assert abs(signed) == pytest.approx(0.230, abs=0.001)
    assert mag == pytest.approx(12.723, abs=0.001)
    assert mag / abs(signed) > 50


def test_on_real_data_the_model_level_shift_was_understated_too():
    A, B = _real()
    assert sc._agg(A, "level_shift_model", magnitude=True) == pytest.approx(13.120, abs=0.001)
    assert sc._agg(B, "level_shift_model", magnitude=True) == pytest.approx(9.621, abs=0.001)


# ---- A2/A3: the paired statistic and the verdict ------------------------------------------- #
def test_a_shift_moving_toward_zero_reads_as_improvement_not_increase():
    """-28 -> -22 is a 6 yr improvement in magnitude. Signed, it is +6, and better_is_down then
    calls it a regression. This is exactly what happened to the C-7 comparison."""
    A = _folds({"N2": -28.0, "N3": -26.0, "O1": -30.0, "O2": -27.0})
    B = _folds({"N2": -22.0, "N3": -20.0, "O1": -24.0, "O2": -21.0})
    md_signed, _, _ = sc._paired(A, B, "level_shift_model")
    md_mag, (lo, hi), n = sc._paired(A, B, "level_shift_model", magnitude=True)
    assert md_signed == pytest.approx(+6.0)
    assert md_mag == pytest.approx(-6.0)
    assert n == 4
    assert sc._verdict("abs", md_mag, lo, hi) == "ACCEPT (better)"
    assert sc._verdict("abs", md_signed, *sc._paired(A, B, "level_shift_model")[1]) == "REGRESSION"


def test_a_shift_moving_away_from_zero_is_still_caught_as_a_regression():
    """The fix must not simply flip signs -- a genuine worsening must still read REGRESSION."""
    A = _folds({"N2": -5.0, "N3": -4.0, "O1": -6.0, "O2": -5.0})
    B = _folds({"N2": -25.0, "N3": -24.0, "O1": -26.0, "O2": -25.0})
    md, (lo, hi), _ = sc._paired(A, B, "level_shift_model", magnitude=True)
    assert md == pytest.approx(+20.0)
    assert sc._verdict("abs", md, lo, hi) == "REGRESSION"


def test_the_two_false_regressions_on_the_real_c7_comparison_are_gone():
    """The regression test for this stage: the exact rows the user saw read REGRESSION."""
    A, B = _real()
    for key in ("level_shift_model", "level_shift_ridge"):
        md_old, (lo_old, hi_old), _ = sc._paired(A, B, key)
        assert sc._verdict("abs", md_old, lo_old, hi_old) == "REGRESSION", "the defect, reproduced"
        md, (lo, hi), _ = sc._paired(A, B, key, magnitude=True)
        assert sc._verdict("abs", md, lo, hi) == "noise (CI incl. 0)"
        assert md < 0, "and the point estimate is in the IMPROVING direction"


def test_the_corrected_magnitudes_match_the_hand_computation():
    A, B = _real()
    md, (lo, hi), n = sc._paired(A, B, "level_shift_model", magnitude=True)
    assert n == 5
    assert md == pytest.approx(-3.118, abs=0.001)
    assert lo == pytest.approx(-9.100, abs=0.001) and hi == pytest.approx(+2.865, abs=0.001)


# ---- B: fold-set alignment ------------------------------------------------------------------ #
def test_common_folds_excludes_errored_and_missing_values():
    A = {"N2": {"k": 1.0}, "N3": {"k": 2.0}, "O1": {"_error": "boom"}, "O2": {"k": None}}
    B = {"N2": {"k": 5.0}, "N3": {"_error": "boom"}, "O1": {"k": 3.0}, "O2": {"k": 4.0}}
    assert sc._common_folds(A, B, "k") == ["N2"]


def test_columns_are_averaged_over_the_common_folds_not_each_snapshots_own():
    """A has an extra fold with an extreme value. It must not enter A's column, because there is
    nothing in B to compare it against."""
    A = _folds({"N2": 1.0, "N3": 1.0, "O1": 100.0})
    B = _folds({"N2": 3.0, "N3": 3.0})
    common = sc._common_folds(A, B, "level_shift_model")
    assert common == ["N2", "N3"]
    assert sc._agg(A, "level_shift_model") == pytest.approx(34.0)          # the old, unaligned path
    assert sc._agg(A, "level_shift_model", only=common) == pytest.approx(1.0)
    assert sc._agg(B, "level_shift_model", only=common) == pytest.approx(3.0)


def test_on_real_data_the_unaligned_columns_disagreed_with_the_verdict_statistic():
    """`dage_mae_model` displayed 14.291 -> 15.713 (a visible +1.42) while the verdict was driven
    by +2.92. Aligned, the column difference IS the verdict statistic."""
    A, B = _real()
    assert sc._agg(A, "dage_mae_model") == pytest.approx(14.291, abs=0.001)        # 6 folds
    assert sc._agg(B, "dage_mae_model") == pytest.approx(15.713, abs=0.001)        # 5 folds
    common = sc._common_folds(A, B, "dage_mae_model")
    assert len(common) == 5 and "N2" not in common
    va = sc._agg(A, "dage_mae_model", only=common)
    assert va == pytest.approx(12.791, abs=0.001)
    md, _, _ = sc._paired(A, B, "dage_mae_model")
    assert sc._agg(B, "dage_mae_model", only=common) - va == pytest.approx(md, abs=1e-9)


# ---- the invariant that can only hold if BOTH fixes are right ------------------------------- #
def test_column_difference_equals_the_mean_diff_on_every_real_metric():
    """`col_B - col_A == mean diff`, for every metric and every direction. False for 13 of 18
    rows before this stage. For an "abs" metric it holds only if the columns and the paired
    statistic take the magnitude at the SAME point in the computation."""
    A, B = _real()
    checked = 0
    for key, (direction, _) in sc.METRICS.items():
        mag = direction == "abs"
        common = sc._common_folds(A, B, key)
        va = sc._agg(A, key, only=common, magnitude=mag)
        vb = sc._agg(B, key, only=common, magnitude=mag)
        md, _, _ = sc._paired(A, B, key, magnitude=mag)
        if va is None or vb is None or md is None:
            continue
        assert vb - va == pytest.approx(md, abs=1e-9), f"{key} violates the invariant"
        checked += 1
    assert checked == 18, "every metric must actually have been exercised"


@pytest.mark.parametrize("direction", ["lower", "higher", "abs", "neutral"])
def test_the_invariant_holds_synthetically_for_every_direction(direction):
    rng = np.random.default_rng(13)
    a = dict(zip(sc.DONORS, rng.normal(0, 15, len(sc.DONORS)), strict=True))
    b = dict(zip(sc.DONORS, rng.normal(0, 15, len(sc.DONORS)), strict=True))
    A, B = _folds(a), _folds(b)
    B["Y2"] = {"_error": "boom"}                      # force an unequal fold set
    mag = direction == "abs"
    common = sc._common_folds(A, B, "level_shift_model")
    va = sc._agg(A, "level_shift_model", only=common, magnitude=mag)
    vb = sc._agg(B, "level_shift_model", only=common, magnitude=mag)
    md, _, _ = sc._paired(A, B, "level_shift_model", magnitude=mag)
    assert len(common) == 5
    assert vb - va == pytest.approx(md, abs=1e-12)


# ---- what must NOT have changed ------------------------------------------------------------- #
def test_measure_fold_still_stores_a_signed_level_shift():
    """The per-donor sign is real information and is read by diag_zero_point. Taking abs() at
    measurement time would destroy it permanently and make every snapshot unreadable in its own
    terms -- the fix belongs at aggregation, not measurement."""
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert 'out["level_shift_model"] = float(np.median(mu[m]) - np.median(y))' in src
    assert 'out["level_shift_ridge"] = float(np.median(r_te[m]) - np.median(y))' in src


def test_no_snapshot_on_disk_was_modified_by_this_stage():
    """Nine snapshots are committed. This stage re-judges them; it does not rewrite them."""
    A, _ = _real()
    assert A["N2"]["level_shift_model"] == pytest.approx(15.0246, abs=1e-4)
    assert A["N3"]["level_shift_model"] == pytest.approx(-28.3459, abs=1e-4)


def test_the_compare_path_no_longer_takes_abs_of_an_aggregate():
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "va, vb = abs(va), abs(vb)" not in src


def test_the_signed_mean_is_still_reported_for_abs_metrics():
    """Replacing one number with another would lose the cancellation question. Both are kept."""
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "signed mean" in src
    assert "(context, never judged)" in src


def test_non_abs_verdicts_are_untouched_by_this_stage():
    """Only the two "abs" rows change verdict. If a "lower"/"higher" verdict moved, the fix
    reached further than it should have."""
    A, B = _real()
    for key, (direction, _) in sc.METRICS.items():
        if direction == "abs":
            continue
        old = sc._paired(A, B, key)
        new = sc._paired(A, B, key, magnitude=False)
        assert old == new


def test_agg_without_only_still_spans_every_fold():
    """Backwards compatibility: the default must remain the historical behaviour."""
    f = _folds({"N2": 1.0, "N3": 2.0, "O1": 3.0})
    assert sc._agg(f, "level_shift_model") == pytest.approx(2.0)


def test_agg_returns_none_when_nothing_is_usable():
    assert sc._agg({"N2": {"_error": "boom"}}, "level_shift_model") is None
    assert sc._agg({}, "level_shift_model", only=["N2"]) is None


def test_paired_needs_at_least_two_folds():
    A, B = _folds({"N2": 1.0}), _folds({"N2": 2.0})
    md, (lo, hi), n = sc._paired(A, B, "level_shift_model", magnitude=True)
    assert md is None and lo is None and hi is None and n == 1
