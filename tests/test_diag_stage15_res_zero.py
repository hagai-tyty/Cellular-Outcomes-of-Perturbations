"""Unit tests for the Stage 15 RES-zero diagnosis.

The claim is an ATTRIBUTION -- that `g(R_eff)` and nothing else produces the zero -- so the tests
have to establish that the other three factors genuinely cannot be zero, not merely that they
happened not to be on these folds.

They also pin the status-precedence artefact, because it is the part most likely to mislead a
future reader: the status field reports 1-3 cells failing on rejuvenation when in fact every cell
does, and a reader trusting those counts would conclude rejuvenation is a minor issue.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage15_res_zero.py"
spec = importlib.util.spec_from_file_location("s15", SRC)
s15 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s15)

RESULTS = ROOT / "results" / "diag_stage15_res_zero_results.json"
has_results = pytest.mark.skipif(not RESULTS.exists(), reason="diagnostic has not been run")


@pytest.fixture(scope="module")
def recorded():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


class P:
    tau_safe, w, k, kappa, z_conf, lam = 0.85, 0.03, 2.0, 5.0, 1.0, 0.0


# ---- the attribution: only g can be zero ---------------------------------------------------- #
def test_phi_is_never_exactly_zero_even_for_a_maximally_unsafe_cell():
    """A sigmoid has no zero. S=0 is the worst possible input and phi is still positive."""
    d = s15.factor_decomposition(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                                 np.array([1.0]), np.array([True]), P)
    assert d["phi_min"] > 0.0
    assert d["n_phi_zero"] == 0


def test_s_to_the_k_is_zero_only_if_s_is_exactly_zero():
    d = s15.factor_decomposition(np.array([0.0, 0.5]), np.zeros(2), np.zeros(2),
                                 np.ones(2), np.ones(2, bool), P)
    assert d["n_s_k_zero"] == 1


def test_the_loss_term_is_identically_one_because_lam_ships_at_zero():
    """One of the four factors of the RES formula is switched off in the shipped config."""
    d = s15.factor_decomposition(np.array([0.5]), np.array([0.99]), np.array([0.0]),
                                 np.array([1.0]), np.array([True]), P)
    assert d["loss_term_min"] == d["loss_term_max"] == 1.0
    cfg = (ROOT / "configs" / "infer" / "default.yaml").read_text(encoding="utf-8")
    assert "lam: 0.0" in cfg


def test_g_is_zero_exactly_when_the_upper_age_bound_is_non_negative():
    mu = np.array([-100.0, -1.0, 5.0])
    sig = np.array([1.0, 10.0, 1.0])          # UB = -99, +9, +6
    d = s15.factor_decomposition(np.full(3, 0.9), np.zeros(3), mu, sig, np.ones(3, bool), P)
    assert d["n_g_zero"] == 2
    assert d["g_max"] > 0.0


# ---- the headroom arithmetic ---------------------------------------------------------------- #
def test_z_required_is_the_largest_z_at_which_some_cell_still_qualifies():
    mu = np.array([-10.0, -1.0])
    sig = np.array([20.0, 2.0])               # -mu/sig = 0.5 and 0.5
    h = s15.headroom(mu, sig, z_conf=1.0)
    assert h["z_required_for_any_cell"] == pytest.approx(0.5)
    assert h["n_upper_bound_negative"] == 0


def test_a_cell_that_qualifies_shows_a_negative_upper_bound_and_z_required_above_the_ship():
    mu = np.array([-30.0])
    sig = np.array([10.0])                    # UB = -20, -mu/sig = 3.0
    h = s15.headroom(mu, sig, z_conf=1.0)
    assert h["min_upper_bound"] == pytest.approx(-20.0)
    assert h["n_upper_bound_negative"] == 1
    assert h["z_required_for_any_cell"] > h["z_conf_shipped"]


def test_z_required_is_none_when_no_cell_is_predicted_to_rejuvenate():
    h = s15.headroom(np.array([5.0, 8.0]), np.array([1.0, 1.0]), z_conf=1.0)
    assert h["z_required_for_any_cell"] is None
    assert h["n_mu_negative"] == 0


# ---- the recorded result -------------------------------------------------------------------- #
@has_results
def test_res_is_zero_on_every_fold(recorded):
    assert recorded["folds"]
    for d, f in recorded["folds"].items():
        assert f["res_all_zero"] is True, f"{d}: RES is not all zero"
        assert f["res_max"] == 0.0


@has_results
def test_the_zero_is_attributed_to_g_alone_on_every_cell_of_every_fold(recorded):
    """THE finding. Not 'g was probably it' -- the other three factors are non-zero for every
    single cell, so the attribution is exhaustive."""
    for d, f in recorded["folds"].items():
        x = f["factors"]
        assert x["n_phi_zero"] == 0, f"{d}: phi contributed a zero"
        assert x["n_s_k_zero"] == 0, f"{d}: S^k contributed a zero"
        assert x["loss_term_min"] == x["loss_term_max"] == 1.0, f"{d}: the lam term is not inert"
        assert x["n_g_zero"] == x["n_cells"], f"{d}: g is not zero for every cell"


@has_results
def test_sigma_dwarfs_mu_on_every_fold(recorded):
    """The mechanism: the model's uncertainty is 2-4.5x its signal, so the confidence-gated
    credit can never open. RES is working exactly as designed; the model has no confidence."""
    ratios = [f["headroom"]["sigma_over_abs_mu_median"] for f in recorded["folds"].values()]
    assert all(r > 1.9 for r in ratios)
    assert max(ratios) > 4.0


@has_results
def test_no_cell_anywhere_comes_within_two_years_of_qualifying(recorded):
    """The closest miss across all six folds is N3 at +1.998 yr. Everything else is further --
    the next nearest is Y1 at +5.40 and the worst is O1 at +27.32."""
    ubs = [f["headroom"]["min_upper_bound"] for f in recorded["folds"].values()]
    assert min(ubs) > 1.99
    assert min(ubs) == pytest.approx(1.998, abs=0.001)
    assert sorted(ubs)[1] == pytest.approx(5.40, abs=0.01)


@has_results
def test_lowering_z_conf_would_not_rescue_this(recorded):
    """`z_required` is 0.235-0.898 against a shipped 1.0. Even z_conf = 0.9 would light up ONE
    cell in ONE fold -- so this is not a near miss for the system, and re-tuning the gate is not
    the fix."""
    zs = [f["headroom"]["z_required_for_any_cell"] for f in recorded["folds"].values()]
    assert all(z is not None for z in zs)
    assert max(zs) < 1.0
    assert max(zs) == pytest.approx(0.898, abs=0.005)
    assert sum(1 for z in zs if z > 0.5) == 1


@has_results
def test_the_status_field_understates_the_rejuvenation_failure(recorded):
    """The precedence artefact, pinned so it cannot mislead later. `compute_res_batch` reports
    OOD first, then UNSAFE, then NO_REJUVENATION -- so only the cells that pass the first two
    gates get labelled with the third. R_eff is zero for 100% of cells; the status field says
    1-3 per fold."""
    for d, f in recorded["folds"].items():
        labelled = f["status_counts"].get("REJECTED_NO_REJUVENATION", 0)
        actual = f["factors"]["n_g_zero"]
        assert labelled < actual, f"{d}: expected the status field to undercount"
        assert labelled <= 3
        assert actual == f["factors"]["n_cells"]


@has_results
def test_res_zero_is_over_determined_by_three_independent_gates(recorded):
    """Fixing any ONE would still leave RES at zero, which is why this is a structural finding
    rather than a bug report."""
    for f in recorded["folds"].values():
        sc = f["status_counts"]
        n = f["factors"]["n_cells"]
        assert f["factors"]["n_g_zero"] == n                       # gate 1: every cell
        assert sc.get("REJECTED_OOD", 0) > 0                       # gate 2
        assert sc.get("REJECTED_UNSAFE", 0) > 0                    # gate 3
        assert sum(sc.values()) == n
        assert "APPROVED" not in sc


@has_results
def test_most_cells_are_rejected_unsafe_in_five_of_six_folds(recorded):
    """S < tau_safe - 3w = 0.76. Recorded as an observation with a hypothesis attached, not as a
    claim: fate PR-AUC is 0.965-0.992, so the head RANKS well while its probabilities sit low --
    consistent with miscalibration rather than genuinely unsafe cells. Untested here."""
    unsafe = {d: f["status_counts"].get("REJECTED_UNSAFE", 0) / f["factors"]["n_cells"]
              for d, f in recorded["folds"].items()}
    assert sum(1 for v in unsafe.values() if v > 0.5) == 5
    assert unsafe["Y1"] < 0.1, "Y1 is the exception -- it is 89% out-of-distribution instead"


@has_results
def test_the_shipped_params_are_the_ones_measured_against(recorded):
    assert recorded["params"] == {"tau_safe": 0.85, "w": 0.03, "k": 2.0,
                                  "kappa": 5.0, "z_conf": 1.0, "lam": 0.0}


@has_results
def test_all_six_folds_were_measured_including_the_one_without_age_labels(recorded):
    """N2 has no ΔAge under C-7 but still produces predictions, so it can and must be included --
    RES does not depend on having a ground-truth label."""
    assert set(recorded["folds"]) == {"N2", "N3", "O1", "O2", "Y1", "Y2"}
    assert recorded["errors"] == {}


# ---- contract -------------------------------------------------------------------------------- #
def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_it_uses_the_shipped_res_implementation_not_a_reimplementation():
    """The decomposition must import the real `_sigmoid`, or it could agree with a formula the
    product does not actually use."""
    assert "from cellfate.inference.res import _sigmoid" in SRC.read_text(encoding="utf-8")
    assert "compute_res_batch" in SRC.read_text(encoding="utf-8")
