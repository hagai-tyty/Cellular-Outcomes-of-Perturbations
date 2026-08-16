"""Unit tests for the out-of-cohort transfer test.

The two failures this script already hit were both SILENT -- a cohort vanishing to an empty frame
rather than raising -- so those are what is pinned.
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_age_transfer.py"
spec = importlib.util.spec_from_file_location("dat", SRC)
dat = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dat)


def test_the_fib_pattern_matches_both_naming_conventions():
    """Gill writes `N3_Fib_Sendai_Exp2`; GSE165177 writes `O1 Fib`. Both are day-0 fibroblasts."""
    pat = re.compile(r"(^|_| )Fib($|_| )")
    for c in ("N3_Fib_Sendai_Exp2", "O1 Fib", "O2 Fib"):
        assert pat.search(c), c
    for c in ("O1_negative_control_13days_exp1", "N3_d21_SSEA4_Sendai_Exp2", "Fibroblast_x"):
        assert not pat.search(c), c


def test_title_normalisation_bridges_space_and_underscore():
    """GSE165177's series matrix says `O1_Fib` while its expression header says `O1 Fib`. Keying
    on the raw title made the age lookup miss and the ENTIRE cohort disappear with no error."""
    assert dat._norm("O1 Fib") == dat._norm("O1_Fib") == "o1_fib"
    assert dat._norm("  N3_Fib_Sendai_Exp2 ") == "n3_fib_sendai_exp2"


def test_normalisation_does_not_merge_distinct_samples():
    assert dat._norm("O1 Fib") != dat._norm("O2 Fib")


def test_the_part2_matrix_is_the_one_that_carries_the_fibroblasts():
    """GSE165177 splits its samples across two files; the main matrix has no Fib column at all, so
    pointing at it yields an empty cohort rather than an error."""
    src = SRC.read_text(encoding="utf-8")
    assert "part2" in src


def test_zscore_centres_and_scales_per_feature():
    rng = np.random.default_rng(0)
    a = rng.normal(5.0, 3.0, size=(20, 4))
    z = dat.zscore(a)
    assert np.allclose(z.mean(0), 0.0, atol=1e-12)
    assert np.allclose(z.std(0), 1.0, atol=1e-12)


def test_zscore_survives_a_constant_feature():
    a = np.column_stack([np.ones(10), np.arange(10.0)])
    assert np.all(np.isfinite(dat.zscore(a)))


def test_the_bars_are_stated_constants():
    assert dat.TRANSFER_MAE_BAR == 20.0
    assert dat.TRANSFER_RHO_BAR == 0.6
    assert dat.ALPHAS == (1.0, 10.0, 100.0, 1000.0, 10000.0)


def test_the_c7_gate_is_applied_to_the_held_out_cohorts():
    """N2_Fib_Sendai_Exp2 is not a transcriptome; it must be rejected by the gate, not by hand."""
    src = SRC.read_text(encoding="utf-8")
    assert "bulk_column_verdict" in src
