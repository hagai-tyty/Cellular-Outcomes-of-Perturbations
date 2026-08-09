"""Change C-7, component A — the bulk integrity gate.

Covers pre-registered bars **B1** (separation on the recorded cohort) and **B3a/B3b** (the gate
can fail in both directions). Both B3 branches execute here: a check whose rejection branch
never runs is not a check — the `verify_1a` lesson.

B1 runs against `results/diag_gill_control_integrity_results.json`, the committed 124-column
census of `GSE165176_Log2_RPM_Sendai_reprogramming`. Reading the recorded artefact rather than
the raw matrix keeps these tests runnable in CI with **no data, no GPU and no network**, which
is the standing convention for `tests/` — and it makes B1 a check on the *recorded* cohort, so
a future change to the gate is graded against the same 124 columns every time.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cellfate.data.integrity import (
    G1_LIBRARY_BAND,
    G2_MIN_LOG2_RANGE,
    REASON_LIBRARY,
    REASON_RANGE,
    bulk_column_verdict,
    linear_library_size,
    log2_dynamic_range,
    screen_bulk_matrix,
)

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "results" / "diag_gill_control_integrity_results.json"

# The five degenerate columns, established in STAGE_1_5_6_SPARSE_CLOCK §5.7 and re-verified in
# §11. `N2_Fib_Sendai_Exp2` is the one that matters: it is a CONTROL, and therefore both donor
# N2's ΔAge zero-point and one of the samples `sigma_gill` is fitted on.
EXPECTED_REJECTS = {
    "Y1_d7_CD13_Sendai_Exp1",
    "N3_d21_SSEA4_Sendai_Exp2",
    "O2_d9_SSEA4_Sendai_Exp1",
    "N2_Fib_Sendai_Exp2",
    "N2_d21_CD13_Sendai_Exp2",
}


def _cohort() -> dict[str, dict]:
    if not COHORT.exists():
        pytest.skip(f"{COHORT.name} not present")
    return json.loads(COHORT.read_text(encoding="utf-8"))["per_column"]


def _verdict_from_summary(rec: dict) -> tuple[bool, str | None]:
    """Apply G1/G2 to a recorded (library, range) pair.

    The gate itself consumes a column; the census records only the two summary statistics the
    gate depends on. Applying them here tests the same thresholds on the same cohort without
    shipping a 8 MB matrix into the test suite.
    """
    lo, hi = G1_LIBRARY_BAND
    if not (lo <= rec["implied_library_size"] <= hi):
        return False, REASON_LIBRARY
    if rec["log2_range"] < G2_MIN_LOG2_RANGE:
        return False, REASON_RANGE
    return True, None


# --------------------------------------------------------------------------- B1
def test_b1_separation_exactly_five_rejects_and_no_false_positives():
    per_col = _cohort()
    assert len(per_col) == 124, f"cohort is {len(per_col)} columns, expected 124"
    rejected = {n for n, rec in per_col.items() if not _verdict_from_summary(rec)[0]}
    assert rejected == EXPECTED_REJECTS, (
        f"missing={sorted(EXPECTED_REJECTS - rejected)} "
        f"false_positives={sorted(rejected - EXPECTED_REJECTS)}")


def test_b1_the_rejected_control_is_the_one_the_defect_is_about():
    per_col = _cohort()
    rec = per_col["N2_Fib_Sendai_Exp2"]
    assert rec["is_control"] is True, "N2_Fib must be a control -- that is why it reaches HFF"
    assert not _verdict_from_summary(rec)[0]


def test_b1_each_condition_independently_rejects_all_five():
    """G1 and G2 are kept as two because they fail differently, not because either is weak."""
    per_col = _cohort()
    lo, hi = G1_LIBRARY_BAND
    for name in EXPECTED_REJECTS:
        rec = per_col[name]
        assert not (lo <= rec["implied_library_size"] <= hi), f"{name} passes G1"
        assert rec["log2_range"] < G2_MIN_LOG2_RANGE, f"{name} fails to fail G2"


def test_b1_margins_are_not_hairline():
    """A threshold that only just separates would flag or miss arbitrarily on a new cohort."""
    per_col = _cohort()
    sound = [r for n, r in per_col.items() if n not in EXPECTED_REJECTS]
    bad = [per_col[n] for n in EXPECTED_REJECTS]
    max_sound_lib = max(r["implied_library_size"] for r in sound)
    min_bad_lib = min(r["implied_library_size"] for r in bad)
    assert G1_LIBRARY_BAND[1] / max_sound_lib > 2.0, "G1 ceiling is too close to a sound column"
    assert min_bad_lib / G1_LIBRARY_BAND[1] > 1.5, "G1 ceiling is too close to a reject"
    assert min(r["log2_range"] for r in sound) - max(r["log2_range"] for r in bad) > 1.0


# --------------------------------------------------------------------------- B3
def test_b3a_a_constant_column_is_rejected():
    """The rejection branch. Y1_d7_CD13 is literally this: 20k genes, one value."""
    col = np.full(20_000, 11.49)
    ok, why = bulk_column_verdict(col)
    assert ok is False
    assert why in (REASON_LIBRARY, REASON_RANGE)


def test_b3b_a_sound_column_is_admitted():
    """The admission branch. Both must execute or the gate is untested in one direction."""
    rng = np.random.default_rng(0)
    # log-normal-ish: most genes near the floor, a minority carrying the signal -- and scaled
    # so the linear sum lands at ~1e6, which is what the units mandate.
    col = np.clip(rng.normal(2.0, 2.2, size=20_000), 0.0, None)
    col[rng.integers(0, 20_000, 40)] = 14.0
    ok, why = bulk_column_verdict(col)
    assert ok is True, f"sound synthetic column rejected: {why}"
    assert why is None


def test_b3_reason_is_none_exactly_when_admitted():
    """The (admitted, reason) invariant `age_label_policy` also holds for (mask, reasons)."""
    for col in (np.full(5_000, 3.0), np.linspace(0.0, 15.0, 5_000), np.zeros(5_000)):
        ok, why = bulk_column_verdict(col)
        assert (why is None) == ok


def test_g2_catches_a_collapsed_column_whose_library_is_fine():
    """G2 must earn its place: a column can sit inside the RPM band and still be collapsed."""
    # ~1e6 linear over 20k genes with a range of only 1 log2 unit.
    col = np.full(20_000, np.log2(1e6 / 20_000 + 1.0))
    col[:10] += 1.0
    assert G1_LIBRARY_BAND[0] <= linear_library_size(col) <= G1_LIBRARY_BAND[1]
    assert log2_dynamic_range(col) < G2_MIN_LOG2_RANGE
    ok, why = bulk_column_verdict(col)
    assert ok is False and why == REASON_RANGE


# --------------------------------------------------------------------------- screen
def test_screen_returns_empty_for_a_clean_matrix():
    rng = np.random.default_rng(1)
    m = np.clip(rng.normal(2.0, 2.2, size=(20_000, 3)), 0.0, None)
    m[rng.integers(0, 20_000, 40), :] = 14.0
    assert screen_bulk_matrix(m, ["a", "b", "c"]) == {}


def test_screen_names_only_the_bad_columns():
    """Column b is constant with an IN-BAND library, so G2 is what catches it.

    Constant at a high value would fail G1 first (that is what `N2_Fib` does: constant at
    11.49 over ~36k genes gives a library of 1.03e+08). Pinning the library in-band here makes
    this test grade G2 specifically rather than accidentally re-testing G1.
    """
    rng = np.random.default_rng(2)
    m = np.clip(rng.normal(2.0, 2.2, size=(20_000, 3)), 0.0, None)
    m[rng.integers(0, 20_000, 40), :] = 14.0
    m[:, 1] = np.log2(1e6 / 20_000 + 1.0)
    assert screen_bulk_matrix(m, ["a", "b", "c"]) == {"b": REASON_RANGE}


def test_a_constant_column_at_a_high_value_fails_the_library_condition_first():
    """`N2_Fib`'s actual signature: constant AND mis-scaled. G1 is the first to fire."""
    col = np.full(36_000, 11.489547)
    assert linear_library_size(col) > G1_LIBRARY_BAND[1]
    ok, why = bulk_column_verdict(col)
    assert ok is False and why == REASON_LIBRARY


def test_screen_refuses_a_mismatched_name_list():
    """A mis-paired verdict would reject the WRONG sample, so this must raise, not guess."""
    m = np.zeros((10, 3))
    with pytest.raises(ValueError, match="refusing to guess"):
        screen_bulk_matrix(m, ["a", "b"])


# --------------------------------------------------------------------------- wiring
def test_the_gate_flag_reaches_injected_sources_not_only_constructed_ones():
    """The bug this pins actually happened, and it made the gate a no-op in production.

    `run` takes `sources: list | None`. When a caller INJECTS sources -- which
    `local_runners/run_multi_local.py` does, and so does every test with a synthetic source --
    `build_sources` is never called. Setting the flag only there meant the first C-7 build came
    out byte-identical to arm A: 42605 samples, N2's degenerate control still present, the gate
    silently inert. A guard that cannot fire is the exact defect C-7 exists to remove.
    """
    from cellfate.data.build_dataset import DataConfig, apply_source_flags

    class _FakeSource:
        def __init__(self):
            self.bulk_integrity_gate = False

    cfg_on = DataConfig(out="x", gene_panel="p", bulk_integrity_gate=True)
    cfg_off = DataConfig(out="x", gene_panel="p", bulk_integrity_gate=False)

    injected = [_FakeSource(), _FakeSource()]
    apply_source_flags(cfg_on, injected)
    assert all(s.bulk_integrity_gate for s in injected), "injected sources did not get the flag"

    apply_source_flags(cfg_off, injected)
    assert not any(s.bulk_integrity_gate for s in injected), "flag-off must also propagate"


def test_apply_source_flags_ignores_sources_that_do_not_have_the_attribute():
    """Single-cell sources never gate, so they simply do not carry the attribute."""
    from cellfate.data.build_dataset import DataConfig, apply_source_flags

    class _NoAttr:
        pass

    apply_source_flags(DataConfig(out="x", gene_panel="p", bulk_integrity_gate=True), [_NoAttr()])
