"""The Stage 1a gate's decision table.

`verify_1a.py` printed the right warning on run 1's data and then graded it PASS, because only
the PASS branch was ever exercised. The operator followed the PASS and lost 3.5 h of GPU time to
a void experiment. These tests drive EVERY branch, so a verdict that should stop the run cannot
silently become one that doesn't.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "verify_1a", Path(__file__).resolve().parents[1] / "verify_1a.py")
verify_1a = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(verify_1a)

EXPECTED = verify_1a.EXPECTED_TRAIN_DONORS       # 5 = six donors minus the held-out one


def _fold(counts: dict[int, int], n_cols: int = 7) -> dict:
    """One fold's check_fold() output, with `counts` the per-donor cell counts in train."""
    return {
        "splits": {
            "train": {"n_cols": n_cols, "n_rows": sum(counts.values()),
                      "dtype_ok": True, "len_ok": True,
                      "n_donors": len(counts), "codes": sorted(counts), "counts": counts},
        },
        "empty_n_cols": n_cols,
        "empty_n_rows": 0,
    }


def _decide(results: dict, all_cols_ok: bool = True):
    ok = list(results)
    bulk, usable = verify_1a.bulk_and_usable(results, ok, lambda c: f"D{c}")
    return verify_1a.decide_verdict(ok, all_cols_ok, usable, bulk)


def test_pass_on_five_clean_donors():
    results = {f"F{i}": _fold({c: 20 for c in range(EXPECTED)}) for i in range(2)}
    status, reason = _decide(results)
    assert status == "PASS", reason


def test_pass_when_a_bulk_corpus_is_present_but_skipped():
    """THE RUN-1 GEOMETRY. Six labels, one of them a corpus -> five usable -> proceed.

    A STOP here would block a configuration the calibration code handles correctly; a PASS
    without the skip would repeat run 1. The distinction is the whole point of the gate.
    """
    counts = {0: 33613}                      # the corpus, ~99.8%
    counts.update({c: 15 for c in range(1, EXPECTED + 1)})
    status, reason = _decide({"N2": _fold(counts)})
    assert status == "PASS", reason
    assert "skipped as bulk corpora" in reason
    assert "D0" in reason                     # and it must NAME the corpus


def test_stop_when_cell_line_is_finer_grained_than_donor():
    """The dangerous direction: more usable donors than folds. Holding out one such group is
    not holding out a donor, so residuals understate cross-donor error and `q` comes out too
    small -- a failure that LOOKS like success."""
    status, reason = _decide({"N2": _fold({c: 20 for c in range(EXPECTED + 3)})})
    assert status == "STOP", reason
    assert "1:1" in reason


def test_stop_when_too_few_donors_survive_the_skip():
    counts = {0: 10_000, 1: 20}               # corpus + a single real donor
    status, reason = _decide({"N2": _fold(counts)})
    assert status == "STOP", reason
    assert "inner-LODO cannot run" in reason


def test_stop_when_folds_disagree_on_donor_count():
    """Folds must be homogeneous; a mixture means the split is not what it claims."""
    results = {
        "N2": _fold({c: 20 for c in range(EXPECTED)}),
        "N3": _fold({c: 20 for c in range(EXPECTED - 1)}),
    }
    assert _decide(results)[0] == "STOP"


def test_fail_on_wrong_column_count():
    results = {"N2": _fold({c: 20 for c in range(EXPECTED)}, n_cols=6)}
    assert _decide(results, all_cols_ok=False)[0] == "FAIL"


def test_cannot_verify_with_no_folds():
    assert _decide({})[0] == "CANNOT_VERIFY"


@pytest.mark.parametrize("corpus_frac", [0.99, 0.90, 0.75, 0.51])
def test_corpus_is_skipped_across_the_dominant_range(corpus_frac):
    """Anything holding more than the floor is skipped, not just the extreme case."""
    others = 100
    corpus = int(others * corpus_frac / (1 - corpus_frac))
    counts = {0: corpus}
    counts.update({c: others // EXPECTED for c in range(1, EXPECTED + 1)})
    bulk, usable = verify_1a.bulk_and_usable({"N2": _fold(counts)}, ["N2"], lambda c: f"D{c}")
    assert usable["N2"] == EXPECTED, f"corpus at {corpus_frac:.0%} was not skipped"
    assert "D0" in bulk["N2"]


def test_a_donor_just_under_the_floor_is_not_skipped():
    """The documented gap, pinned so a change to the threshold is a deliberate act.

    A donor at 49% of the split is KEPT -- holding it out still leaves 51%, above the floor --
    yet it goes on to supply ~49% of the pooled residuals, so `q` is shaped largely by that one
    donor while tripping neither the skip nor the >50% pool warning. Whether 50% is the right
    floor is a threshold decision, not an accident; this pins the current behaviour so any
    change to it is deliberate.
    """
    counts = {0: 49, 1: 26, 2: 25}           # nobody exceeds the floor
    bulk, usable = verify_1a.bulk_and_usable({"N2": _fold(counts)}, ["N2"], lambda c: f"D{c}")
    assert bulk == {}, "no donor exceeds the floor, so none may be skipped"
    assert usable["N2"] == 3

    # ...and the boundary itself: 51% IS skipped, because holding it out leaves 49% < 50%
    bulk2, usable2 = verify_1a.bulk_and_usable(
        {"N2": _fold({0: 51, 1: 49})}, ["N2"], lambda c: f"D{c}")
    assert bulk2 == {"N2": {"D0": 51}}
    assert usable2["N2"] == 1
