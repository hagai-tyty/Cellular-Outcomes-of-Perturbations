"""Change C-7, components C and D — rule 4 (`no_control_baseline`) and bar B2'.

Covers **B3c** (rule 4 masks a line with no zero-point), **B3d** (it does NOT fire on the
chunk-local case), **B2'a/b/c** (the invariant holds, fails, and covers the S4 site), the **C1**
reason-ordering decision, and **B4** (nothing moves while the flag is off).

**B3d is the test that matters most.** `_control_baseline` falls back in two different
situations and they have different owners:

  * **no controls ANYWHERE** for that line -- Stage 1.5 Group D, what rule 4 is for;
  * **controls exist but none landed in THIS chunk** -- Stage 1.5 Group E, separately owned.

A rule 4 that cannot tell them apart fires on Group E and blocks C-7 for the wrong reason. The
census is global by construction (`GillReprogrammingSource.lines_without_controls` reads every
column of the matrix), so a Group E line never enters the set -- and B3d pins that.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cellfate.data.aging import (
    LinearClock,
    age_label_policy,
    assert_no_unmasked_fallback,
    delta_age,
    recenter_on_control_arrays,
)


def _obs(lines, ctrl, donor_age=None):
    d = {"cell_line": list(lines), "is_control": list(ctrl)}
    if donor_age is not None:
        d["donor_age"] = list(donor_age)
    return pd.DataFrame(d)


# --------------------------------------------------------------------------- B3c
def test_b3c_rule4_masks_a_line_with_no_zero_point():
    obs = _obs(["N2", "N2", "O1"], [False, False, True])
    mask, reasons = age_label_policy(3, "reprogramming", obs,
                                     lines_without_controls=frozenset({"N2"}))
    assert list(mask) == [False, False, True]
    assert reasons[0] == reasons[1] == "no_control_baseline"
    assert reasons[2] is None


def test_b3c_rule4_is_keyed_on_the_condition_not_on_identity():
    """Any line in the set is masked; the rule never names a donor."""
    obs = _obs(["A", "B", "C"], [False, False, False])
    mask, reasons = age_label_policy(3, "s", obs, lines_without_controls=frozenset({"B"}))
    assert list(mask) == [True, False, True]
    assert reasons[1] == "no_control_baseline"


# --------------------------------------------------------------------------- B3d
def test_b3d_rule4_does_not_fire_on_the_chunk_local_fallback():
    """Group E: this chunk has NO controls, but the line has them elsewhere in the corpus.

    The global census therefore does not list the line, and the label must survive. If this
    test ever fails, rule 4 has started firing on Group E and C-7 is blocking for the wrong
    reason.
    """
    obs = _obs(["X", "X", "X"], [False, False, False])   # no controls IN THIS CHUNK
    mask, reasons = age_label_policy(3, "s", obs, lines_without_controls=frozenset())
    assert list(mask) == [True, True, True]
    assert reasons == [None, None, None]


def test_b3d_the_two_fallbacks_are_distinguishable_end_to_end():
    """The same chunk, the same absence of controls -- opposite verdicts, decided globally."""
    clock = LinearClock({"G0": 1.0}, intercept=40.0)
    expr = np.array([[5.0], [7.0], [9.0]])
    obs = _obs(["X", "X", "X"], [False, False, False])

    _, mask_e, _ = delta_age(clock, expr, ["G0"], obs, source="s",
                             lines_without_controls=frozenset())          # Group E
    _, mask_d, _ = delta_age(clock, expr, ["G0"], obs, source="s",
                             lines_without_controls=frozenset({"X"}))     # Group D
    assert mask_e.all()
    assert not mask_d.any()


# --------------------------------------------------------------------------- C1 ordering
def test_c1_rule4_precedes_donor_out_of_clock_range():
    """N2 is donor_age 0 and the clock starts at 1.0, so it matches BOTH rules.

    `age_mask_reason` is persisted (`io.py:139`, `:265`), so the order decides what is written
    to the shard. "No zero-point exists" is undefined; "outside the fitted range" is
    out-of-validity. Undefined is the stronger claim and is the one recorded.
    """
    obs = _obs(["N2"], [False], donor_age=[0.0])
    _, reasons = age_label_policy(1, "s", obs,
                                  clock_age_range=(1.0, 96.0),
                                  lines_without_controls=frozenset({"N2"}))
    assert reasons[0] == "no_control_baseline"


def test_c1_cancer_source_still_wins_over_rule4():
    """Rule 1 is documented as never weakened by the later rules. That must stay true."""
    import cellfate.common.constants as C
    src = next(iter(C.CANCER_SOURCES))
    obs = _obs(["L"], [False])
    _, reasons = age_label_policy(1, src, obs, lines_without_controls=frozenset({"L"}))
    assert reasons[0] == "cancer_source"


# --------------------------------------------------------------------------- B2'
def test_b2prime_passes_when_the_fallback_line_is_masked():
    census = {"X": {"n_control": 0, "n_cells": 3, "source": "self_fallback"}}
    lines = np.array(["X", "X", "X"])
    assert_no_unmasked_fallback(census, lines, np.zeros(3, dtype=bool))   # must not raise


def test_b2prime_raises_when_the_fallback_line_keeps_its_label():
    census = {"X": {"n_control": 0, "n_cells": 3, "source": "self_fallback"}}
    lines = np.array(["X", "X", "X"])
    with pytest.raises(AssertionError, match="B2' violated"):
        assert_no_unmasked_fallback(census, lines, np.ones(3, dtype=bool))


def test_b2prime_ignores_lines_that_used_real_controls():
    census = {"X": {"n_control": 2, "n_cells": 3, "source": "controls"}}
    lines = np.array(["X", "X", "X"])
    assert_no_unmasked_fallback(census, lines, np.ones(3, dtype=bool))    # must not raise


def test_b2prime_covers_the_s4_recentring_site():
    """`recenter_on_control_arrays` is `_control_baseline`'s SECOND call site.

    It passed no census, so its fallback was invisible -- a B2' guarding only `delta_age` would
    pass while this site silently self-centred the same orphaned line. It now accepts a census,
    and this pins that the fallback is actually recorded there.
    """
    census: dict = {}
    vals = np.array([1.0, 2.0, 3.0])
    lines = np.array(["X", "X", "X"])
    out = recenter_on_control_arrays(vals, lines, np.zeros(3, dtype=bool), census=census)
    assert census["X"]["source"] == "self_fallback"
    assert out.mean() == pytest.approx(0.0)
    with pytest.raises(AssertionError, match="B2' violated"):
        assert_no_unmasked_fallback(census, lines, np.ones(3, dtype=bool))


# --------------------------------------------------------------------------- B4
def test_b4_defaults_change_nothing_in_the_policy():
    """Rule 4 is empty by default, so every existing caller is unaffected."""
    obs = _obs(["A", "B"], [True, False])
    mask, reasons = age_label_policy(2, "s", obs)
    assert list(mask) == [True, True]
    assert reasons == [None, None]


def test_b4_delta_age_does_not_enforce_b2prime_by_default():
    """With the C-7 flag off, `delta_age` must not raise on the Group E fallback.

    `tests/test_harmonize.py::test_the_silent_no_control_fallback_self_centres_a_line_to_zero`
    pins that behaviour deliberately -- "any future change to this behaviour a deliberate,
    reviewed act". Enforcing B2' unconditionally WAS such a change and is gated for this reason.
    """
    clock = LinearClock({"G0": 1.0}, intercept=40.0)
    expr = np.array([[5.0], [7.0], [9.0]])
    obs = _obs(["X", "X", "X"], [False, False, False])
    d, mask, _ = delta_age(clock, expr, ["G0"], obs, source="s")   # no enforcement
    assert d.mean() == pytest.approx(0.0)
    assert mask.all()


def test_b4_recenter_without_a_census_is_unchanged():
    vals = np.array([1.0, 2.0, 3.0])
    lines = np.array(["X", "X", "X"])
    a = recenter_on_control_arrays(vals, lines, np.zeros(3, dtype=bool))
    b = recenter_on_control_arrays(vals, lines, np.zeros(3, dtype=bool), census={})
    assert np.array_equal(a, b), "adding a census must not move a single value"
