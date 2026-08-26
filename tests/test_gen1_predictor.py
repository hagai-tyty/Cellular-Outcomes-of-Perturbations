"""Contracts for the Gen-1 predictor — the Stage-23.5 §6 tool contract.

Three of these are about what the tool REFUSES to do, and those matter more than the scoring:
a tool that quietly imputes a missing nuisance block, or maps an unknown drug onto a known one, or
presents an unvalidated sort order as condition selection, would report frozen-benchmark numbers for
something the benchmark never evaluated.

The artifact itself is 44 MB and gitignored (it rebuilds in ~0.5 min via `--stage 24c`), so the
scoring tests skip when it is absent. The contract tests that need no artifact always run.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
ART_NPZ = ROOT / "results" / "stage24" / "stage24_w5_artifact.npz"
ART_META = ROOT / "results" / "stage24" / "stage24_w5_artifact.json"
SRC = ROOT / "src" / "cellfate" / "gen1_predictor.py"

has_artifact = pytest.mark.skipif(
    not (ART_NPZ.exists() and ART_META.exists()),
    reason="the 24C artifact is gitignored; rebuild with --stage 24c")
has_meta = pytest.mark.skipif(not ART_META.exists(), reason="24C has not been run")


@pytest.fixture(scope="module")
def predictor():
    from cellfate.gen1_predictor import Gen1Predictor
    return Gen1Predictor.load(ART_NPZ, ART_META)


@pytest.fixture(scope="module")
def one_clone():
    """A real benchmark clone, its nuisance vector, and its frozen out-of-fold row."""
    import sys

    import pandas as pd
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_stage23_learnability_gate as S23
    X, clones = S23._load_wm989_x()
    ck = pd.read_csv(ROOT / "results" / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis = np.column_stack([np.log1p(ck[c].to_numpy(float)) for c in
                            ("n_naive_cells", "n_naive1_cells", "n_naive2_cells",
                             "n_naive3_cells")])
    frozen = pd.read_csv(ROOT / "results" / "stage23_wm989_interaction_oof.csv")
    cid = frozen["clone_id"].iloc[0]
    i = clones.index(cid)
    return {"clone_id": cid, "x": np.asarray(X[i].todense()).ravel(), "b": nuis[i],
            "rows": frozen[frozen.clone_id == cid]}


# ============================================================================================== #
# What the tool refuses
# ============================================================================================== #
@has_artifact
def test_a_missing_nuisance_block_fails_closed_and_is_never_imputed(predictor, one_clone):
    """§6.2. Expression alone is not equivalent to the evaluated model."""
    rows = predictor.predict(one_clone["x"], None)
    assert rows, "a refusal must still return one row per requested condition"
    for r in rows:
        assert r["support_status"] == "MISSING_REQUIRED_NUISANCE"
        assert r["future_detection_score"] is None, "a refused row must carry no score"
    assert "never imputed" in rows[0]["detail"]


@has_artifact
def test_an_incomplete_or_nonfinite_nuisance_block_also_fails_closed(predictor, one_clone):
    for bad in (one_clone["b"][:2], np.array([np.nan, 1.0, 2.0, 3.0])):
        rows = predictor.predict(one_clone["x"], bad)
        assert all(r["support_status"] == "MISSING_REQUIRED_NUISANCE" for r in rows)
        assert all(r["future_detection_score"] is None for r in rows)


@has_artifact
def test_an_unknown_condition_is_refused_not_mapped(predictor, one_clone):
    """§6.3. Never embedded, nearest-neighboured, or silently mapped to a known condition."""
    rows = predictor.predict(one_clone["x"], one_clone["b"],
                             treatments=["Dabrafenib", "Aspirin", "Pembrolizumab"])
    by = {r["condition"]: r for r in rows}
    assert by["Dabrafenib"]["support_status"] == "SUPPORTED_KNOWN_CONDITION"
    for unknown in ("Aspirin", "Pembrolizumab"):
        assert by[unknown]["support_status"] == "UNSUPPORTED_TREATMENT"
        assert by[unknown]["future_detection_score"] is None
    # and the known condition's score is unaffected by the unknown ones being present
    alone = predictor.predict(one_clone["x"], one_clone["b"], treatments=["Dabrafenib"])
    assert by["Dabrafenib"]["future_detection_score"] == alone[0]["future_detection_score"]


@has_artifact
def test_a_wrong_feature_schema_is_refused(predictor, one_clone):
    rows = predictor.predict(np.zeros(10), one_clone["b"])
    assert all(r["support_status"] == "UNSUPPORTED_FEATURE_SCHEMA" for r in rows)


@has_artifact
def test_the_condition_order_is_withheld_until_stage_25_validates_it(predictor, one_clone):
    """§6.1/§6.4. Sorting six scores is not condition selection until Stage 25 says it is."""
    r = predictor.rank_conditions(one_clone["x"], one_clone["b"])
    assert r["ranking_status"] == "NOT_SUPPORTED"
    assert r["validated_condition_order"] is None
    assert len(r["scores"]) == 6, "the six scores are always returned"
    assert "not a validated condition ranking" in r["detail"]


@has_artifact
def test_the_order_appears_only_on_a_real_stage_25_supported_verdict(predictor, one_clone, tmp_path):
    """The gate is the recorded verdict, not a flag anyone can set."""
    from cellfate.gen1_predictor import Gen1Predictor
    for verdict, expect_order in (("STAGE_25_RANKING_NOT_SUPPORTED", False),
                                  ("STAGE_25_RANKING_SUPPORTED", True)):
        v = tmp_path / f"{verdict}.json"
        v.write_text(json.dumps({"verdict": verdict}), encoding="utf-8")
        p = Gen1Predictor.load(ART_NPZ, ART_META, stage25_verdict=v)
        r = p.rank_conditions(one_clone["x"], one_clone["b"])
        got = r["validated_condition_order"]
        assert (got is not None) is expect_order, verdict
        if expect_order:
            assert sorted(got) == sorted(p.conditions)


# ============================================================================================== #
# What the tool computes
# ============================================================================================== #
@has_artifact
def test_the_shipped_api_reproduces_the_frozen_out_of_fold_predictions(predictor, one_clone):
    """The whole artifact is worthless if it does not regenerate the frozen column."""
    g = one_clone["rows"]
    f = int(g["outer_fold"].iloc[0])
    rows = predictor.predict(one_clone["x"], one_clone["b"],
                             treatments=list(g["treatment"]), component=f"fold{f}")
    got = np.array([r["future_detection_score"] for r in rows])
    assert np.abs(got - g["pred_W5"].to_numpy()).max() < 1e-12


@has_artifact
def test_every_response_carries_its_provenance_and_limitations(predictor, one_clone):
    for r in predictor.predict(one_clone["x"], one_clone["b"]):
        assert r["model_version"] and r["feature_contract_version"]
        assert r["ranking_status"] in ("SUPPORTED", "NOT_SUPPORTED")
        assert r["known_limitations"], "limitations travel with every response"


@has_artifact
def test_both_component_families_are_shipped(predictor):
    """Fold components reproduce the benchmark; deployment scores a clone in no fold."""
    assert "deployment" in predictor.components
    assert {f"fold{i}" for i in range(5)} <= set(predictor.components)


# ============================================================================================== #
# The artifact's own declarations
# ============================================================================================== #
@has_meta
def test_the_deployment_component_is_declared_as_packaging_not_validation():
    m = json.loads(ART_META.read_text(encoding="utf-8"))["deployment_component"]
    assert "PACKAGING, not a new model" in m["status"]
    assert "NOT validated on held-out data" in m["validation"]
    assert "ESTIMATED by the frozen out-of-fold result" in m["validation"]


@has_meta
def test_the_limitations_name_the_things_that_would_otherwise_be_overclaimed():
    lim = " ".join(json.loads(ART_META.read_text(encoding="utf-8"))["known_limitations"])
    assert "known conditions only" in lim
    assert "may not be imputed" in lim
    assert "single cell is not an equivalent input" in lim
    assert "no calibrated probability" in lim
    assert "no independent biological replication" in lim
    assert "3.45x the state contribution" in lim


def test_the_module_documents_its_three_refusals():
    doc = SRC.read_text(encoding="utf-8")
    assert "MISSING_REQUIRED_NUISANCE" in doc
    assert "UNSUPPORTED_TREATMENT" in doc
    assert "never embedded, nearest-neighboured" in doc
    assert "not a calibrated probability" in doc


def test_pseudobulk_normalisation_is_applied_exactly_once():
    from cellfate.gen1_predictor import PredictionError, clone_pseudobulk_from_counts
    counts = np.array([[1.0, 0.0, 3.0], [2.0, 2.0, 2.0]])
    got = clone_pseudobulk_from_counts(counts)
    want = np.log1p(counts.sum(axis=0) * (10_000.0 / counts.sum()))
    assert np.allclose(got, want)
    # summing already-normalised cells would double-normalise; a zero clone is a blocking condition
    with pytest.raises(PredictionError):
        clone_pseudobulk_from_counts(np.zeros((2, 3)))
