"""Unit tests for STAGE 1.5.2 M-2b — pure functions only, no repo data.

The join logic gets the most attention here. §9-R5 is "I am wrong about the join", and §10 step 1
answers it with an ABORT rather than a warning — so the abort has to actually fire, on every way
the join can be wrong. A branch that never executes is not a check.
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


M = _load("diag_m2b_contrast_agreement", "experiments/diag_m2b_contrast_agreement.py")

METH = ["Y2_d11_SSEA4", "Y2_d11_CD13"]
RNA = ["Y2_d11_SSEA4_Sendai_Exp1", "Y2_d11_SSEA4_Sendai_Exp2", "Y2_d11_CD13_Sendai_Exp1"]


# ------------------------------------------------------------------------ join_key ---- #
@pytest.mark.parametrize("title,expect", [
    ("Y2_d11_SSEA4", "Y2_d11_SSEA4"),                       # methylation title IS the key
    ("Y2_d11_SSEA4_Sendai_Exp1", "Y2_d11_SSEA4"),           # RNA title = key + batch suffix
    ("N2_d9_CD13_Sendai_Exp2", "N2_d9_CD13"),
    ("O1_transiently_reprogrammed_17days_exp1", None),      # the TRANSIENT vocabulary
    ("iPSC_control", None),
    ("Y2_d11_SSEA5", None),                                 # unknown marker
    ("", None),
])
def test_join_key_accepts_only_the_shape_it_documents(title, expect):
    assert M.join_key(title) == expect


def test_a_foreign_series_produces_no_join_at_all_rather_than_a_wrong_one():
    """GSE165177's titles must not accidentally key against GSE165178's."""
    j = M.verify_join(METH, ["O1_transiently_reprogrammed_17days_exp1"], expected=2)
    assert j["n_matched"] == 0 and j["verdict"] == "ABORT"


# ---------------------------------------------------------------------- verify_join ---- #
def test_a_complete_join_passes():
    j = M.verify_join(METH, RNA, expected=2)
    assert j["verdict"] == "OK"
    assert j["n_matched"] == 2 and j["arm_counts"] == {"CD13": 1, "SSEA4": 1}
    assert j["rna_replicates_per_key"]["Y2_d11_SSEA4"] == 2      # replicates counted, not dropped


def test_a_short_join_aborts():
    j = M.verify_join(METH, RNA[:1], expected=2)
    assert j["verdict"] == "ABORT" and j["n_matched"] == 1


def test_an_unmatched_methylation_sample_aborts_even_at_the_expected_count():
    """The subtle failure: the right NUMBER of matches, but not the right ones."""
    meth = [*METH, "O9_d99_CD13"]
    j = M.verify_join(meth, RNA, expected=2)
    assert j["n_matched"] == 2                # count is right...
    assert j["unmatched_meth"] == ["O9_d99_CD13"]
    assert j["verdict"] == "ABORT"            # ...and it still aborts


# -------------------------------------------------------------------- pair_contrasts ---- #
def test_replicates_are_averaged_not_treated_as_independent():
    """1.5.1's unit-of-analysis rule. Averaging two SSEA4 reads must give their mean."""
    keys = ["Y2_d11_SSEA4", "Y2_d11_CD13"]
    out = M.pair_contrasts({"Y2_d11_SSEA4": 10.0, "Y2_d11_CD13": 4.0}, keys)
    assert len(out) == 1
    assert out[0]["delta"] == pytest.approx(6.0)


def test_an_unpaired_condition_is_dropped_not_half_counted():
    out = M.pair_contrasts({"Y2_d11_SSEA4": 10.0}, ["Y2_d11_SSEA4", "Y2_d11_CD13"])
    assert out == []


def test_pairs_never_cross_donors_or_days():
    vals = {"Y2_d11_SSEA4": 10.0, "O1_d15_CD13": 1.0}
    assert M.pair_contrasts(vals, list(vals)) == []


# ------------------------------------------------------------------- sign_agreement ---- #
def test_sign_agreement_counts_matched_directions():
    s = M.sign_agreement([-1.0, -2.0, +3.0], [-5.0, +6.0, +7.0])
    assert s["n_pairs"] == 3 and s["n_agree"] == 2


def test_sign_agreement_ignores_non_finite_pairs():
    s = M.sign_agreement([-1.0, float("nan")], [-5.0, -6.0])
    assert s["n_pairs"] == 1


# ---------------------------------------------------------------------- m2b_verdict ---- #
def test_exactly_on_the_bar_is_reported_as_fragile():
    """Three hairline margins have already misled this project (0.009, 0.014, 0.016)."""
    v = M.m2b_verdict({"n_pairs": 11, "n_agree": 7}, 0.5, bar=7)
    assert v["status"] == "AGREE_FRAGILE" and v["fragile"] is True
    assert "EXACTLY on the moved bar" in v["caveat"]


def test_clearing_the_bar_outright_is_not_fragile():
    v = M.m2b_verdict({"n_pairs": 11, "n_agree": 10}, 0.5, bar=7)
    assert v["status"] == "AGREE" and v["fragile"] is False


def test_below_the_bar_disagrees_and_carries_no_weakened_bar_caveat():
    v = M.m2b_verdict({"n_pairs": 11, "n_agree": 6}, 0.1, bar=7)
    assert v["status"] == "DISAGREE"
    assert v["caveat"] == ""          # the caveat is about a PASS on a loosened bar


def test_too_few_pairs_cannot_verify():
    assert M.m2b_verdict({"n_pairs": 1, "n_agree": 1}, float("nan"))["status"] == "CANNOT_VERIFY"


# ------------------------------------------------------------------ agreement_by_day ---- #
def test_by_day_splits_a_pooled_number_that_would_otherwise_hide_the_structure():
    rows = [{"day": 9, "same_sign": False, "delta_rna_years": 40.0, "delta_meth_years": -2.0},
            {"day": 9, "same_sign": False, "delta_rna_years": 60.0, "delta_meth_years": -1.0},
            {"day": 15, "same_sign": True, "delta_rna_years": -50.0, "delta_meth_years": -70.0}]
    b = M.agreement_by_day(rows)
    assert b["9"]["n_agree"] == 0 and b["9"]["n"] == 2
    assert b["9"]["mean_rna"] == pytest.approx(50.0)
    assert b["15"]["n_agree"] == 1


def test_the_recorded_result_still_matches_the_rules():
    """Replay 2026-07-31's measured pairs through the pure functions."""
    import json
    p = ROOT / "results" / "diag_m2b_contrast_agreement_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    res = json.loads(p.read_text(encoding="utf-8"))
    assert res["join"]["n_matched"] == 22 and res["join"]["verdict"] == "OK"
    for blk in res["clocks"].values():
        rows = blk["pairs"]
        s = M.sign_agreement([r["delta_rna_years"] for r in rows],
                             [r["delta_meth_years"] for r in rows])
        assert M.m2b_verdict(s, blk["rho"])["status"] == blk["status"]
        # the finding: day 9 -- where methylation says nothing has happened -- agrees 0 of 3
        assert blk["by_day"]["9"]["n_agree"] == 0
        assert blk["by_day"]["9"]["mean_rna"] > 30.0
        assert blk["by_day"]["15"]["n_agree"] == blk["by_day"]["15"]["n"]
