"""Unit tests for the robustness sweep.

The sweep's value rests on two structural properties rather than on any computation of its own:
that every variant runs the IDENTICAL procedure as the headline result, and that each variant
changes exactly ONE axis. Both are checkable without touching the data.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_residual_robustness.py"
spec = importlib.util.spec_from_file_location("drr", SRC)
drr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(drr)

BASE_MARKERS = ("CD13", "SSEA4", "Fib")
BASE_WINDOW = (7.0, 29.0)
BASE_FEATURES = "panel"


def test_the_procedure_is_imported_not_reimplemented():
    """If the sweep had its own copy of the LOO / permutation code it could drift from the result
    it is meant to be testing, and a 'robust' verdict would mean nothing."""
    src = SRC.read_text(encoding="utf-8")
    assert "diag_residual_expression.py" in src
    for fn in ("loo_spearman", "permutation_null"):
        assert f"dre.{fn}" in src, f"{fn} must come from the headline module, not a local copy"
    assert f"def {'loo_spearman'}" not in src, "the sweep must not define its own loo_spearman"


def test_the_baseline_is_the_published_configuration():
    name, mk, win, fs = drr.RUNS[0]
    assert name.startswith("BASELINE")
    assert (mk, win, fs) == (BASE_MARKERS, BASE_WINDOW, BASE_FEATURES)


def test_every_variant_changes_exactly_one_axis():
    """A grid would be multiple-testing soup; one-axis-at-a-time is what makes each row readable."""
    for name, mk, win, fs in drr.RUNS[1:]:
        changed = sum([mk != BASE_MARKERS, win != BASE_WINDOW, fs != BASE_FEATURES])
        assert changed == 1, f"{name!r} changes {changed} axes, not 1"


def test_all_three_named_axes_are_actually_exercised():
    """The three checks that were promised: early window, marker, feature set."""
    axes = set()
    for _, mk, win, _fs in drr.RUNS[1:]:
        axes.add("marker" if mk != BASE_MARKERS else
                 "window" if win != BASE_WINDOW else "features")
    assert axes == {"window", "marker", "features"}


def test_the_feature_sets_include_a_gene_specificity_probe():
    assert any(fs == "random" for _, _, _, fs in drr.RUNS), (
        "a random gene set asks whether the effect is specific to a gene set or just to donors "
        "resembling one another")


def test_the_robustness_bar_is_a_stated_constant_and_a_majority():
    assert drr.MIN_ROBUST == 6
    assert drr.MIN_ROBUST > len(drr.RUNS) / 2, "the bar must be a majority of runs"
    assert drr.MIN_ROBUST < len(drr.RUNS), (
        "unanimity is the wrong bar: SSEA4-only and the shortest window are expected to be weak")


def test_the_random_gene_seed_is_fixed():
    assert isinstance(drr.RANDOM_GENES_SEED, int)
    assert drr.N_RANDOM == 500 and drr.N_HVG == 500


def test_the_sweep_reuses_the_headline_thresholds():
    """A variant scored against different alphas or a different percentile would not be
    comparable to the result it is testing."""
    src = SRC.read_text(encoding="utf-8")
    assert "dre.ALPHAS" in src and "dre.PERM_PCTILE" in src
    assert "dre.MIN_ALPHAS_PASSING" in src
