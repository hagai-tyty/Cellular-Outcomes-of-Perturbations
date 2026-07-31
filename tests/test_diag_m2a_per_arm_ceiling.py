"""Unit tests for STAGE 1.5.2 §17's per-arm ceiling re-audit — pure functions, no repo data.

The point of the re-audit is that a correlation's numerator is meaningless without its denominator,
so the tests focus on the `interpretable` flag: an arm where the two references disagree with each
other must never be readable as evidence about a third instrument, in either direction.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

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


C = _load("diag_m2a_per_arm_ceiling", "experiments/diag_m2a_per_arm_ceiling.py")


def _rows(arm, meth_a, meth_b, rna):
    return ([{"donor": "O1", "arm": arm, "day": i, "age_meth": v, "age_rna": r}
             for i, (v, r) in enumerate(zip(meth_a, rna, strict=True))],
            [{"donor": "O1", "arm": arm, "day": i, "age_meth": v, "age_rna": r}
             for i, (v, r) in enumerate(zip(meth_b, rna, strict=True))])


# ------------------------------------------------------------------ per_arm_ceiling ---- #
def test_a_sharp_reference_with_an_inverted_rna_clock_is_interpretable():
    """The load-bearing case: the two references agree, so the RNA clock's failure is real."""
    a, b = _rows("x", [1, 2, 3, 4, 5], [1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
    out = C.per_arm_ceiling(a, b)["x"]
    assert out["meth_vs_meth"] == pytest.approx(1.0)
    assert out["rna_mean"] == pytest.approx(-1.0)
    assert out["interpretable"] is True


def test_a_blunt_reference_is_flagged_uninterpretable_even_with_a_high_rna_correlation():
    """The trap §11 fell into: a good-looking RNA number in an arm where the references
    do not agree with each other says nothing."""
    a, b = _rows("x", [1, 2, 3, 4, 5], [3, 1, 5, 2, 4], [1, 2, 3, 4, 5])
    out = C.per_arm_ceiling(a, b)["x"]
    assert out["meth_vs_meth"] < C.SHARP_CEILING
    assert out["rna_mean"] > 0.4               # looks fine on its own...
    assert out["interpretable"] is False       # ...and is not readable


def test_arms_below_the_minimum_n_are_dropped_not_reported_with_a_nan():
    a, b = _rows("tiny", [1, 2], [1, 2], [1, 2])
    assert C.per_arm_ceiling(a, b) == {}


def test_only_conditions_present_in_both_clocks_are_used():
    a, b = _rows("x", [1, 2, 3, 4, 5], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
    b = b[:4]                                   # one condition missing from the second clock
    out = C.per_arm_ceiling(a, b)["x"]
    assert out["n"] == 4


def test_reprogramming_arms_are_labelled_from_the_vocabulary_not_from_the_numbers():
    a, b = _rows("transiently_reprogrammed", [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4])
    assert C.per_arm_ceiling(a, b)["transiently_reprogrammed"]["reprogramming"] is True
    a, b = _rows("negative_control", [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4])
    assert C.per_arm_ceiling(a, b)["negative_control"]["reprogramming"] is False


# ---------------------------------------------------------------- reading_correction ---- #
def test_no_interpretable_arm_blocks_every_reading():
    arms = {"a": {"meth_vs_meth": 0.3, "interpretable": False, "rna_mean": 0.9}}
    d = C.reading_correction(arms)
    assert d["status"] == "NO_ARM_IS_INTERPRETABLE"
    assert "either direction" in d["reason"]


def test_the_sharpest_and_bluntest_arms_are_both_named():
    arms = {
        "sharp": {"meth_vs_meth": 0.94, "interpretable": True, "rna_mean": -0.16},
        "blunt": {"meth_vs_meth": 0.23, "interpretable": False, "rna_mean": 0.15},
        "mid": {"meth_vs_meth": 0.86, "interpretable": True, "rna_mean": 0.40},
    }
    d = C.reading_correction(arms)
    assert d["sharpest_arm"] == "sharp" and d["bluntest_arm"] == "blunt"
    assert d["n_interpretable"] == 2 and d["n_arms"] == 3
    assert "sharp" in d["rna_fails_in_interpretable_arms"]
    assert "mid" not in d["rna_fails_in_interpretable_arms"]


# ------------------------------------------------------------------ the recorded run ---- #
def test_the_recorded_result_is_what_the_rules_produce():
    import json
    p = ROOT / "diag_m2a_per_arm_ceiling_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    r = json.loads(p.read_text(encoding="utf-8"))
    assert C.reading_correction(r["arms"])["status"] == r["reading_correction"]["status"]


def test_the_finding_that_corrects_section_11():
    """§11 read the per-arm table as a clean in-domain/out-of-domain boundary. Pin the two
    facts that make that reading wrong: the SHARPEST reference is in a REPROGRAMMING arm, and
    at least one NON-reprogramming arm with a sharp reference also fails."""
    import json
    p = ROOT / "diag_m2a_per_arm_ceiling_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    arms = json.loads(p.read_text(encoding="utf-8"))["arms"]
    sharpest = max(arms.items(), key=lambda kv: kv[1]["meth_vs_meth"])
    assert sharpest[1]["reprogramming"] is True
    assert sharpest[1]["rna_mean"] < 0
    interp_non_reprog = [v for v in arms.values()
                         if v["interpretable"] and not v["reprogramming"]]
    assert any(v["rna_mean"] < 0.2 for v in interp_non_reprog)
