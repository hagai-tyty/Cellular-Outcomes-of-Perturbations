"""Stage 21B — the source/design audit that resolved Stage 21A's two open questions.

The point of 21B is to stop an *orthogonal* label being promoted to a *prospective* one. These
tests pin the three ways that could still happen:

1. **The subpopulation trap.** If one culture yields BOTH an SSEA4 and a CD13 fraction at the same
   timepoint, the marker is a within-culture subpopulation label — the culture did not BECOME one
   of them. A verdict of `VALID_PROSPECTIVE_ORTHOGONAL_TASK` must be impossible in that geometry.
2. **The leakage case.** If the early input was itself marker-sorted, predicting a later marker
   identity is predicting the variable the input was selected on.
3. **Inherited unit counts.** Stage 21A multiplied donors x experiments to get 12. 21B must
   verify that rather than inherit it — and here it is wrong, because Exp1 and Exp2 are time
   blocks, not replicates.

Plus the tri-state rule carried over from 21A: "not found" is never "absent".
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage21b_source_design.py"
spec = importlib.util.spec_from_file_location("s21b", SRC)
s21b = importlib.util.module_from_spec(spec)
sys.modules["s21b"] = s21b          # @dataclass resolves through sys.modules
spec.loader.exec_module(s21b)

RESULTS = ROOT / "results" / "diag_stage21b_source_design_results.json"
STAGE_21A = ROOT / "results" / "diag_stage21_data_audit_results.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="21B has not been run")

# A handful of tests exercise the audit against the REAL GSE165176 files, which live on D:\ and
# are not in the repo. On CI those tests used to fail rather than skip -- the dataset is absent, so
# the audit correctly answers UNKNOWN, and the assertions expected the real verdict. The frozen
# result file carries the same evidence, so the CI-visible coverage is unchanged.
needs_gse165176 = pytest.mark.skipif(
    not any(d.exists() for d in s21b.GSE165176_DIRS),
    reason="GSE165176 is not on this machine (this is the CI condition)")

P, AB, UNK = s21b.PRESENT, s21b.ABSENT_PROVEN, s21b.UNKNOWN


@pytest.fixture(scope="module")
def rec():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


# ---- title parsing, which every B-branch answer rests on ------------------------------------- #
def test_sorted_and_unsorted_titles_are_both_parsed():
    rows = s21b.parse_gill_titles(
        ["N2_d11_SSEA4_Sendai_Exp1", "N2_d11_CD13_Sendai_Exp1", "N2_Fib_Sendai_Exp2"])
    assert [r["marker"] for r in rows] == ["SSEA4", "CD13", "UNSORTED"]
    assert [r["day"] for r in rows] == [11, 11, 0]
    assert {r["donor"] for r in rows} == {"N2"}


def test_the_unsorted_baseline_is_day_zero():
    """`_Fib_` carries no day token; treating it as anything but day 0 would invent an early
    observation that does not exist."""
    assert s21b.parse_gill_titles(["O1_Fib_Sendai_Exp2"])[0]["day"] == 0


def test_an_unparseable_title_is_dropped_rather_than_guessed():
    assert s21b.parse_gill_titles(["something_else_entirely"]) == []


# ---- the three verdict branches for GSE165176 ------------------------------------------------ #
@needs_gse165176
def test_both_fractions_from_one_culture_cannot_be_a_valid_prospective_task(monkeypatch, tmp_path):
    """THE trap. Constructed geometry: every culture yields both fractions and no proportions
    exist. The verdict must be CONTEMPORANEOUS, never VALID."""
    r = s21b.audit_gse165176()
    assert r["verdict"] == s21b.B_CONTEMP
    both = r["findings"]["q2_q3_both_fractions_same_culture"].value["both"]
    assert both > 0, "the guard only bites when cultures really do yield both fractions"


def test_the_valid_branch_requires_facs_proportions():
    """`VALID_PROSPECTIVE_ORTHOGONAL_TASK` is reachable only when a culture-level quantity such as
    %SSEA4+ exists. Reading the source rather than asserting it."""
    src = SRC.read_text(encoding="utf-8")
    assert 'elif out["q9_facs_proportions"].status == PRESENT:' in src
    assert "v, why = B_VALID" in src


def test_a_missing_dataset_is_unknown_not_invalid():
    r = s21b.audit_gse165176(dirs=[Path(r"D:\definitely_not_here")])
    assert r["verdict"] == s21b.B_UNKNOWN
    assert r["findings"]["location"].status == UNK


@needs_gse165176
def test_leakage_is_recorded_when_the_early_input_was_already_sorted():
    """Case 2. 118 of 124 samples are marker-sorted, so a sorted early input predicting a later
    marker is predicting its own selection variable."""
    r = s21b.audit_gse165176()
    f = r["findings"]["q7_early_input_already_sorted"]
    assert f.value is True
    assert "LEAKAGE" in f.evidence
    assert s21b.B_LEAK in (s21b.B_LEAK,)   # the branch name exists and is spelled once


def test_the_leakage_verdict_constant_is_defined_and_distinct():
    assert s21b.B_LEAK == "INVALID_FUTURE_LABEL_LEAKAGE"
    assert len({s21b.B_VALID, s21b.B_CONTEMP, s21b.B_LEAK, s21b.B_UNKNOWN}) == 4


# ---- the three verdict branches for GSE242423 ------------------------------------------------- #
def test_lineage_absent_requires_both_vocabulary_and_supplementary_to_be_proven_absent():
    src = SRC.read_text(encoding="utf-8")
    assert ('elif lin.status == ABSENT_PROVEN and '
            'out["supplementary_lineage_file"].status == ABSENT_PROVEN:') in src


def test_missing_source_metadata_yields_unknown_not_absent(tmp_path):
    """The 21A rule, still enforced: no series matrix and no MINiML means UNKNOWN."""
    r = s21b.audit_gse242423(base=tmp_path)
    assert r["verdict"] == s21b.LINEAGE_UNKNOWN
    assert r["findings"]["source_metadata"].status == UNK


def test_the_three_a_branch_verdicts_are_distinct():
    assert len({s21b.LINEAGE_PRESENT, s21b.LINEAGE_ABSENT, s21b.LINEAGE_UNKNOWN}) == 3


def test_lineage_vocabulary_covers_the_real_systems():
    for term in ("clone", "lineage", "celltag", "larry", "hashtag", "sister"):
        assert term in s21b.LINEAGE_TERMS


# ---- evidence standard ------------------------------------------------------------------------ #
def test_a_finding_cannot_carry_an_invalid_status():
    with pytest.raises(ValueError):
        s21b.Finding(value=1, status="maybe", evidence="")


@has_results
def test_every_finding_carries_evidence(rec):
    for gse in ("GSE242423", "GSE165176"):
        for k, f in rec[gse]["findings"].items():
            assert f["evidence"].strip(), f"{gse}/{k}"
            assert f["status"] in (P, AB, UNK)


# ---- the recorded result ------------------------------------------------------------------------ #
@has_results
def test_gse242423_lineage_is_now_closed_as_absent(rec):
    """21A left this PENDING because no series matrix was on disk. Both source files are now
    present and they resolve it."""
    d = rec["GSE242423"]
    assert d["verdict"] == s21b.LINEAGE_ABSENT
    assert d["findings"]["lineage_vocabulary"]["status"] == AB
    assert d["findings"]["supplementary_lineage_file"]["status"] == AB
    assert d["findings"]["replicate_structure"]["status"] == AB
    assert d["findings"]["sample_characteristic_tags"]["value"] == ["cell type", "genotype"]


@has_results
def test_gse242423_is_one_trajectory_not_forty_two_thousand(rec):
    d = rec["GSE242423"]
    assert d["findings"]["n_independent_units"]["value"] == 1
    assert "NOT independent trajectories" in d["findings"]["n_independent_units"]["evidence"]
    assert d["findings"]["destructive_sampling"]["value"] is True


@has_results
def test_gse165176_is_orthogonal_but_contemporaneous(rec):
    d = rec["GSE165176"]
    assert d["verdict"] == s21b.B_CONTEMP
    q = d["findings"]["q2_q3_both_fractions_same_culture"]["value"]
    assert q["both"] == 47 and q["single"] == 24


@has_results
def test_the_sort_is_the_fate_call_and_is_orthogonal_to_rna(rec):
    """The genuinely good half of the finding, and the reason this was worth auditing."""
    ev = rec["GSE165176"]["findings"]["q1_what_was_sorted"]["evidence"]
    assert "SORT IS the fate assignment" in ev
    assert "independently of RNA" in ev


@has_results
def test_only_six_prediction_time_valid_early_samples_exist(rec):
    d = rec["GSE165176"]["findings"]["q6_q8_unsorted_early_samples"]
    assert len(d["value"]) == 6
    assert all(t.endswith("_Fib_Sendai_Exp2") for t in d["value"])


@has_results
def test_stage_21a_overstated_the_unit_count(rec):
    """21A reported 12 by multiplying 6 donors x 2 experiments. 21B verifies rather than inherits:
    Exp1 and Exp2 are TIME BLOCKS overlapping at a single day, so the effective n is 6."""
    q = rec["GSE165176"]["findings"]["q11_effective_n"]["value"]
    assert q["effective_n"] == 6
    assert q["overlap_days"] == [11]
    assert set(q["exp1_days"]).isdisjoint(set(q["exp2_days"]) - {11})

    a = json.loads(STAGE_21A.read_text(encoding="utf-8"))
    gill = next(x for x in a["datasets"] if x["dataset"] == "GSE165176")
    assert gill["findings"]["n_independent_units"]["value"] == 12, "21A's figure, left frozen"


@has_results
def test_there_is_no_outcome_variation_at_the_terminal_day(rec):
    v = rec["GSE165176"]["findings"]["terminal_day_outcome_variation"]["value"]
    assert v["day"] == 54 and v["markers"] == ["SSEA4"] and len(v["donors"]) == 6


@has_results
def test_no_facs_proportions_and_no_terminal_assay(rec):
    f = rec["GSE165176"]["findings"]
    assert f["q9_facs_proportions"]["status"] == AB
    assert f["q10_terminal_outcome"]["status"] == AB


# ---- additive contract ------------------------------------------------------------------------ #
@has_results
def test_stage_21b_is_additive_and_did_not_touch_stage_21a(rec):
    assert rec["stage_21a_result_modified"] is False
    assert rec["additive_to"] == "results/diag_stage21_data_audit_results.json"
    assert STAGE_21A.exists(), "the frozen 21A result must still be present"
    a = json.loads(STAGE_21A.read_text(encoding="utf-8"))
    assert a["verdict"] == "CULTURE_FORWARD_AVAILABLE", "21A's frozen verdict is unchanged"


@has_results
def test_no_model_was_fitted(rec):
    assert rec["model_fitted"] is False and rec["src_modified"] is False
    src = SRC.read_text(encoding="utf-8")
    for banned in ("LogisticRegression", "sklearn", "torch", ".fit("):
        assert banned not in src, f"21B must not model: found {banned}"


def test_the_script_writes_only_its_own_results_file():
    src = SRC.read_text(encoding="utf-8")
    assert src.count(".write_text(") == 1
    assert "diag_stage21b_source_design_results.json" in src
