"""Stage 21A — the local dataset geometry audit.

Per plan §20, every verdict branch must be hit by a constructed test. Beyond that, these pin the
two rules the audit exists to enforce, because both are easy to violate silently:

* **every classification carries its evidence** — a bare `TRAJECTORY_FORWARD` is unauditable;
* **"not found" is never "proven absent"** — a dataset must not be discarded because a
  supplementary file was not downloaded.

They also pin four real bugs found while building it, each of which produced a WRONG verdict on
real data:

1. day tokens appear as `d11` AND `13days`; missing the second made GSE165177 look like it had no
   timecourse at all;
2. `GM00731_D0` + `GM23815_D0` is two donors at ONE timepoint, not two timepoints;
3. an UNKNOWN `n_timepoints` was silently driving `INVALID` instead of `PENDING`;
4. asserting "no independent outcome" everywhere hid the SSEA4/CD13 antibody sorts, which are the
   only non-RNA phenotype in the local corpora.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage21_data_audit.py"
spec = importlib.util.spec_from_file_location("s21a", SRC)
s21 = importlib.util.module_from_spec(spec)
# `@dataclass` resolves annotations through `sys.modules[cls.__module__]`, so the module must be
# registered BEFORE exec_module or the decorator raises on a module loaded purely by path.
sys.modules["s21a"] = s21
spec.loader.exec_module(s21)

RESULTS = ROOT / "results" / "diag_stage21_data_audit_results.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="21A has not been run")


def _audit(**findings):
    a = s21.DatasetAudit(dataset="X", path="/tmp/x", present=True)
    for k, (v, st, ev) in findings.items():
        a.add(k, v, st, ev)
    return a


P, AB, UNK = s21.PRESENT, s21.ABSENT_PROVEN, s21.UNKNOWN


# ---- every verdict branch (plan §20) --------------------------------------------------------- #
def test_lineage_plus_independent_outcome_gives_strict_lineage():
    a = s21.classify(_audit(lineage_link=(True, P, "clone_id column"),
                            independent_outcome=(True, P, "imaging endpoint"),
                            n_timepoints=(4, P, "4 days")))
    assert a.level == s21.LEVEL_3


def test_culture_links_only_gives_culture_forward():
    a = s21.classify(_audit(lineage_link=(False, AB, "no clone column"),
                            independent_outcome=(True, P, "sorted phenotype"),
                            n_independent_units=(12, P, "6 donors x 2 experiments"),
                            n_timepoints=(11, P, "11 days")))
    assert a.level == s21.LEVEL_2


def test_timecourse_with_no_link_gives_trajectory_forward():
    a = s21.classify(_audit(lineage_link=(False, AB, "none"),
                            independent_outcome=(False, AB, "RNA only"),
                            n_independent_units=(1, P, "one culture"),
                            n_timepoints=(9, P, "9 days")))
    assert a.level == s21.LEVEL_1


def test_cross_sectional_dataset_is_invalid():
    a = s21.classify(_audit(lineage_link=(False, AB, "none"),
                            independent_outcome=(False, AB, "RNA only"),
                            n_independent_units=(143, P, "143 donors"),
                            n_timepoints=(1, P, "one collection")))
    assert a.level.startswith(s21.LEVEL_0)


def test_d0_only_dataset_is_invalid():
    """Two donors sampled once is not a timecourse."""
    a = s21.classify(_audit(lineage_link=(False, AB, "none"),
                            independent_outcome=(False, AB, "RNA only"),
                            n_independent_units=(2, P, "2 donors"),
                            n_timepoints=(1, P, "D0 only")))
    assert a.level.startswith(s21.LEVEL_0)


def test_a_future_derived_predictor_is_a_hard_failure():
    """Plan §19.4 — an outcome-derived column among the predictors is an error, not a warning."""
    assert s21.check_forbidden_predictors(["time_h", "donor", "final_fate"]) == ["final_fate"]
    assert s21.check_forbidden_predictors(["X_0", "survivor_label"]) == ["survivor_label"]
    assert s21.check_forbidden_predictors(["time_h", "dose_uM"]) == []


# ---- RULE 2: not found is not proven absent -------------------------------------------------- #
def test_unknown_lineage_does_not_rule_out_strict_lineage():
    """THE rule. A missing series matrix must leave LEVEL 3 open, flagged for one more download."""
    a = s21.classify(_audit(lineage_link=(None, UNK, "series matrix not on disk"),
                            independent_outcome=(False, AB, "RNA only"),
                            n_independent_units=(1, P, "one culture"),
                            n_timepoints=(9, P, "9 days")))
    assert a.level == s21.LEVEL_1 + s21.PENDING
    assert "NOT RULED OUT" in a.ruled_out[s21.LEVEL_3]
    assert "lineage_link" in a.level_reason


def test_proven_absent_does_rule_it_out():
    """The contrast: when the file that WOULD carry it exists and does not, the level closes."""
    a = s21.classify(_audit(lineage_link=(False, AB, "372 characteristics, none mention clone"),
                            independent_outcome=(False, AB, "RNA only"),
                            n_independent_units=(6, P, "6 donors"),
                            n_timepoints=(11, P, "11 days")))
    assert a.level == s21.LEVEL_1
    assert s21.PENDING not in a.level
    assert "ruled out" in a.ruled_out[s21.LEVEL_3]


def test_unknown_timepoints_forces_pending_rather_than_invalid():
    """Bug 3: a parsing failure was driving a REJECTION. It must ask for evidence instead."""
    a = s21.classify(_audit(lineage_link=(False, AB, "none"),
                            independent_outcome=(False, AB, "none"),
                            n_independent_units=(4, P, "4 donors"),
                            n_timepoints=(None, UNK, "no day token matched")))
    assert a.level == s21.LEVEL_0 + s21.PENDING
    assert "n_timepoints" in a.level_reason


def test_a_missing_dataset_is_unknown_not_absent():
    a = s21.audit_dataset("NOT_A_REAL_GSE", candidates=[r"D:\definitely_not_here"])
    assert a.present is False
    assert a.findings["location"].status == UNK
    assert "NOT the same as proven absent" in a.findings["location"].evidence


def test_a_dataset_can_be_found_under_an_alias_path():
    """GSE165176 lives at D:\\Gill. Checking only the canonical path would report a present
    dataset as missing -- the same error class as calling UNKNOWN 'absent'."""
    assert r"D:\Gill" in s21.DEFAULT_LOCATIONS["GSE165176"]


# ---- RULE 1: evidence on every finding -------------------------------------------------------- #
def test_a_finding_cannot_be_created_without_a_valid_status():
    with pytest.raises(ValueError):
        s21.Finding(value=1, status="probably", evidence="vibes")


@has_results
def test_every_recorded_finding_carries_non_empty_evidence():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    for d in r["datasets"]:
        for k, f in d["findings"].items():
            assert f["evidence"].strip(), f"{d['dataset']}/{k} has no evidence"
            assert f["status"] in (P, AB, UNK)


@has_results
def test_every_dataset_records_why_each_level_was_ruled_out():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    for d in r["datasets"]:
        if d["level"].startswith(s21.LEVEL_3):
            continue
        assert d["level_reason"].strip()


# ---- the four parsing bugs -------------------------------------------------------------------- #
def test_day_tokens_are_parsed_in_both_formats():
    """Bug 1: GSE165177 writes `13days`; GSE165178 writes `d11`. Missing either format reports a
    real timecourse as having none."""
    assert s21.parse_days(["Y2_d11_SSEA4", "O1_d15_CD13"]) == [11, 15]
    assert s21.parse_days(["O1_failed_to_transiently_reprogram_13days_exp1"]) == [13]
    assert s21.parse_days(["no_day_here"]) == []


def test_donor_prefix_is_stripped_before_counting_timepoints():
    """Bug 2: `GM00731_D0` + `GM23815_D0` is 2 donors x 1 timepoint. Counting labels makes a
    D0-only corpus look like a 2-point timecourse and would license a forward task."""
    src = SRC.read_text(encoding="utf-8")
    assert "TWO\n    # DONORS AT ONE TIMEPOINT" in src or "TWO DONORS AT ONE TIMEPOINT" in src
    assert 'donor_re = re.compile(r"^(GM\\d+|[NOY]\\d+)_")' in src


def test_surface_marker_sorts_are_recognised_as_an_orthogonal_phenotype():
    """Bug 4: SSEA4/CD13 are ANTIBODY phenotypes, not RNA. Asserting 'no independent outcome'
    everywhere hid the only non-RNA readout in the local corpora."""
    assert s21.find_orthogonal_phenotype(["Y2_d11_SSEA4"], []) == ["SSEA4"]
    assert s21.find_orthogonal_phenotype(["N2_d21_CD13_Sendai_Exp2"], []) == ["CD13"]
    assert s21.find_orthogonal_phenotype(["plain_sample"], []) == []


def test_plain_10x_barcodes_prove_only_that_the_barcode_file_carries_no_tag(tmp_path):
    """The distinction that keeps GSE242423's LEVEL 3 open: encoding is provable, study design
    is not."""
    import gzip
    p = tmp_path / "b.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("AAACCCAAGAAACACT-1\nAAACCCAAGAAACCCA-1\n")
    enc, why = s21.barcode_encoding(p)
    assert enc == "PLAIN_10X" and "matches" in why
    doc = SRC.read_text(encoding="utf-8")
    assert "not that the STUDY lacked lineage tracing" in doc


# ---- the recorded run -------------------------------------------------------------------------- #
@has_results
def test_all_seven_datasets_were_located():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert len(r["datasets"]) == 7
    assert all(d["present"] for d in r["datasets"]), "every local corpus should be found"


@has_results
def test_gse242423_keeps_level_3_open_because_its_metadata_is_missing():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    d = next(x for x in r["datasets"] if x["dataset"] == "GSE242423")
    assert d["findings"]["series_matrix"]["status"] == UNK
    assert d["findings"]["clone_id_in_barcode"]["status"] == AB   # encoding IS provable
    assert d["findings"]["lineage_link"]["status"] == UNK          # study design is NOT
    assert s21.PENDING in d["level"]


@has_results
def test_gse297234_is_d0_only():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    d = next(x for x in r["datasets"] if x["dataset"] == "GSE297234")
    assert d["findings"]["n_timepoints"]["value"] == 1
    assert d["level"].startswith(s21.LEVEL_0)


@has_results
def test_the_surface_marker_finding_is_recorded_with_its_contemporaneous_caveat():
    """The lead this audit surfaced — and the reason it is not yet a forward task."""
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    d = next(x for x in r["datasets"] if x["dataset"] == "GSE165176")
    assert d["findings"]["independent_outcome"]["value"] is True
    assert d["findings"]["outcome_is_rna_surrogate"]["value"] is False
    assert d["findings"]["outcome_is_contemporaneous"]["value"] is True
    assert "AT COLLECTION" in d["findings"]["outcome_is_contemporaneous"]["evidence"]


@has_results
def test_the_audit_fitted_no_model_and_left_src_alone():
    r = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert r["src_modified"] is False
    src = SRC.read_text(encoding="utf-8")
    for banned in ("LogisticRegression", "fit(", "sklearn"):
        assert banned not in src, f"Stage 21A must not model: found {banned}"
