"""Unit tests for Stage 10 -- pluripotency: contamination or mediation?

Every decision branch of the plan's §10.2-10.5 is exercised on constructed input, because a
branch that never runs is not a check, and this stage's verdict withdraws a previous
recommendation.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage10_pluri.py"
spec = importlib.util.spec_from_file_location("s10", SRC)
s10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s10)


# ---- residualise: the operation under scrutiny ---------------------------------------------- #
def test_residualising_removes_the_linear_component_entirely():
    x = np.linspace(0, 10, 50)
    y = 3.0 * x + 7.0
    assert np.allclose(s10.residualise(y, x), 0.0, atol=1e-9)


def test_residualising_leaves_an_orthogonal_component_untouched():
    """If y carries signal NOT aligned with x, residualising must preserve it -- otherwise the
    operation would destroy signal even under the contamination reading."""
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    extra = rng.normal(size=200)
    r = s10.residualise(2.0 * x + extra, x)
    assert abs(np.corrcoef(r, extra)[0, 1]) > 0.9


def test_residualising_against_a_constant_only_centres():
    y = np.array([1.0, 2.0, 3.0])
    assert np.allclose(s10.residualise(y, np.ones(3)), y - y.mean())


# ---- pluri_score --------------------------------------------------------------------------- #
def test_the_score_is_the_mean_of_z_scored_signature_genes():
    expr = np.array([[0.0, 10.0], [2.0, 20.0], [4.0, 30.0]])
    genes = [s10.PLURI[0], s10.PLURI[1]]
    out = s10.pluri_score(expr, genes)
    assert out.shape == (3,) and out[0] < out[1] < out[2]
    assert abs(float(out.mean())) < 1e-12          # z-scored, so centred


def test_a_gene_that_is_off_contributes_no_variance():
    """The finding that decided Test B: in untreated fibroblasts the signature genes are OFF, so
    the score is exactly constant and cannot covary with anything."""
    expr = np.zeros((5, 2))
    out = s10.pluri_score(expr, [s10.PLURI[0], s10.PLURI[1]])
    assert float(np.std(out)) == 0.0


def test_missing_signature_genes_raise_rather_than_score_silently():
    with pytest.raises(ValueError):
        s10.pluri_score(np.zeros((3, 2)), ["NOT_A_GENE", "ALSO_NOT"])


# ---- verdict_from: every branch of §10.5 ---------------------------------------------------- #
def test_two_contamination_votes_carry_it():
    assert s10.verdict_from("CONTAMINATION", "CONTAMINATION", "MEDIATION") == "CONTAMINATION"


def test_two_mediation_votes_carry_it():
    assert s10.verdict_from("MEDIATION", "MEDIATION", "CONTAMINATION") == "MEDIATION"


def test_the_observed_result_is_unanimous_mediation():
    assert s10.verdict_from("MEDIATION", "MEDIATION", "MEDIATION") == "MEDIATION"


def test_one_vote_plus_two_undetermined_is_undetermined():
    """A single test must not decide a stage that withdraws a recommendation."""
    assert s10.verdict_from("CONTAMINATION", "UNDETERMINED", "UNDETERMINED") == "UNDETERMINED"
    assert s10.verdict_from("MEDIATION", "UNDETERMINED", "UNDETERMINED") == "UNDETERMINED"


def test_a_one_all_split_is_undetermined():
    assert s10.verdict_from("CONTAMINATION", "MEDIATION", "UNDETERMINED") == "UNDETERMINED"


def test_all_undetermined_is_undetermined():
    assert s10.verdict_from("UNDETERMINED", "UNDETERMINED", "UNDETERMINED") == "UNDETERMINED"


# ---- the bars, and the stage's own restraint ------------------------------------------------ #
def test_the_bars_are_stated_constants():
    assert s10.GAP_COLLAPSE_BAR == 0.50
    assert s10.CONTROL_RHO_BAR == 0.50
    assert s10.PLURI == ("NANOG", "POU5F1", "LIN28A", "SOX2", "ZFP42")


def test_the_stage_does_not_touch_src():
    """Plan §10.5: no outcome authorises a src change. The most a positive result buys is the
    right to PROPOSE one as a separate, separately pre-registered Change."""
    src = SRC.read_text(encoding="utf-8")
    assert "src/ is NOT changed" in src
    # exactly one write, and it is the results JSON
    assert src.count(".write_text(") == 1
