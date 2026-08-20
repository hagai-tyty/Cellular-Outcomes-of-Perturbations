"""Unit tests for Stage 16 — the safety floor rejecting demonstrably safe cells.

Three ways this stage could report a false success, each pinned:

* the ORACLE calibrator being read as a result. It is fitted on the cells it scores and is an
  upper bound only; tests assert it is named and reported as such, and that the deployable arm
  never touches the test labels.
* a calibrator "fixing" false rejections by shoving every probability upward, trading them for
  false approvals. Both directions are asserted to be reported.
* the verdict being decided by anything other than the pre-registered rule in plan §16.4.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage16_safety_floor.py"
spec = importlib.util.spec_from_file_location("s16", SRC)
s16 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s16)

RESULTS = ROOT / "results" / "diag_stage16_safety_floor_results.json"
def _find_plan(name: str) -> Path:
    """Locate a plan by FILENAME, wherever it currently sits under plans/.

    Plans were reorganised into `(older)base plans/` and `(newer)practical plans/` on 2026-08-20,
    which broke every test that read one by a fixed path. Searching by name survives the next
    reshuffle too.
    """
    hits = sorted(ROOT.joinpath("plans").rglob(name))
    if not hits:
        raise AssertionError(f"plan {name!r} not found anywhere under plans/")
    return hits[0]


has_results = pytest.mark.skipif(not RESULTS.exists(), reason="diagnostic has not been run")


@pytest.fixture(scope="module")
def recorded():
    return json.loads(RESULTS.read_text(encoding="utf-8"))


# ---- the confusion, in BOTH directions ------------------------------------------------------ #
def test_a_truly_safe_cell_below_the_bar_is_a_false_rejection():
    c = s16.confusion(np.array([0.5]), np.array([True]), 0.76)
    assert c["false_rejections"] == 1 and c["false_approvals"] == 0
    assert c["false_rejection_rate"] == 1.0


def test_a_truly_unsafe_cell_above_the_bar_is_a_false_approval():
    c = s16.confusion(np.array([0.9]), np.array([False]), 0.76)
    assert c["false_approvals"] == 1 and c["false_rejections"] == 0


def test_both_error_directions_are_reported_so_the_trade_is_visible():
    """THE trap: a calibrator that shoves every probability up 'fixes' false rejections by
    creating false approvals. Reporting only one direction would hide that."""
    S = np.array([0.5, 0.5, 0.9, 0.9])
    truth = np.array([True, False, True, False])
    c = s16.confusion(S, truth, 0.76)
    assert c["false_rejections"] == 1 and c["false_approvals"] == 1
    assert set(c) >= {"false_rejections", "false_approvals", "sensitivity", "specificity",
                      "balanced_accuracy"}


def test_pushing_everything_above_the_bar_is_caught_by_balanced_accuracy():
    S = np.array([0.5, 0.5, 0.6, 0.6])
    truth = np.array([True, True, False, False])
    before = s16.confusion(S, truth, 0.76)
    after = s16.confusion(np.full(4, 0.99), truth, 0.76)
    assert after["false_rejections"] < before["false_rejections"]     # "improved"
    assert after["false_approvals"] > before["false_approvals"]       # by trading
    assert after["balanced_accuracy"] == before["balanced_accuracy"] == 0.5


def test_a_fold_with_no_unsafe_cells_reduces_balanced_accuracy_to_sensitivity():
    """N2 is 19/19 truly safe. Specificity is undefined there and must not poison the number."""
    c = s16.confusion(np.array([0.9, 0.5]), np.array([True, True]), 0.76)
    assert c["n_unsafe"] == 0
    assert np.isnan(c["specificity"])
    assert c["balanced_accuracy"] == pytest.approx(c["sensitivity"]) == pytest.approx(0.5)


def test_a_perfect_separation_scores_one():
    S = np.array([0.9, 0.95, 0.1, 0.2])
    truth = np.array([True, True, False, False])
    assert s16.confusion(S, truth, 0.76)["balanced_accuracy"] == pytest.approx(1.0)


# ---- Platt: no leak in the deployable arm --------------------------------------------------- #
def test_platt_fitted_on_one_set_and_applied_to_another_is_the_deployable_form():
    rng = np.random.default_rng(16)
    p_fit = rng.uniform(0.2, 0.8, 200)
    y_fit = (rng.random(200) < p_fit).astype(int)
    out, fitted = s16.platt(p_fit, y_fit, np.array([0.3, 0.7]))
    assert fitted is True and out.shape == (2,)


def test_platt_passes_through_when_the_fit_set_has_a_single_class():
    """Unidentifiable boundary. It must pass through, not silently emit a constant."""
    out, fitted = s16.platt(np.array([0.4, 0.6]), np.array([1, 1]), np.array([0.4, 0.6]))
    assert fitted is False
    assert np.allclose(out, [0.4, 0.6])


def test_platt_is_monotone_so_it_cannot_change_the_ranking():
    """Calibration must move probabilities, never the order -- otherwise a 'calibration' gain
    would really be a ranking change and PR-AUC would move with it."""
    rng = np.random.default_rng(17)
    p_fit = rng.uniform(0.1, 0.9, 300)
    y_fit = (rng.random(300) < p_fit).astype(int)
    score = np.linspace(0.05, 0.95, 50)
    out, _ = s16.platt(p_fit, y_fit, score)
    assert np.all(np.diff(out) >= -1e-12) or np.all(np.diff(out) <= 1e-12)


# ---- the oracle threshold sweep -------------------------------------------------------------- #
def test_best_threshold_recovers_a_known_optimum():
    S = np.array([0.1, 0.2, 0.8, 0.9])
    truth = np.array([False, False, True, True])
    r = s16.best_threshold(S, truth)
    assert 0.2 < r["oracle_threshold"] <= 0.8
    assert r["oracle_balanced_accuracy"] == pytest.approx(1.0)


def test_best_threshold_is_undefined_with_only_one_class():
    r = s16.best_threshold(np.array([0.1, 0.9]), np.array([True, True]))
    assert r["oracle_threshold"] is None and r["oracle_balanced_accuracy"] is None


# ---- the pre-registered verdict rule --------------------------------------------------------- #
def test_the_bar_matches_the_plan():
    plan = _find_plan("STAGE_16_SAFETY_FLOOR_MISCALIBRATION.md").read_text("utf-8")
    assert "50 % reduction" in plan or "50 %" in plan
    assert s16.FALSE_REJECTION_DROP_BAR == 0.50


@pytest.mark.parametrize("dep,orc,expect", [
    (True, True, "H1"),
    (False, True, "H4"),
    (False, False, "H3"),
    (True, False, "H3"),          # deployable beating the oracle is incoherent -> not H1
])
def test_every_verdict_branch_is_reachable_and_matches_the_plan(dep, orc, expect):
    v = s16.verdict_from({"deployable_clears_bar": dep, "oracle_clears_bar": orc})
    assert v.startswith(expect)


# ---- the post-hoc prior-gap check ------------------------------------------------------------ #
def _fold(test_frac, fr_rate, n_safe=10):
    return {"raw": {"n_safe": n_safe, "false_rejection_rate": fr_rate},
            "priors": {"train_safe_frac": 0.532, "calib_safe_frac": 0.515,
                       "test_safe_frac": test_frac}}


def test_prior_gap_correlation_is_positive_when_a_bigger_gap_costs_more():
    folds = {"a": _fold(0.45, 0.10), "b": _fold(0.60, 0.30),
             "c": _fold(0.80, 0.60), "d": _fold(1.00, 0.90)}
    r = s16.prior_gap_correlation(folds)
    assert r["spearman"] == pytest.approx(1.0)
    assert r["n"] == 4


def test_prior_gap_correlation_reports_its_own_power_rather_than_a_p_value():
    """n=6 cannot establish this. The critical value is carried in the result so nobody has to
    look it up, and `significant` is computed against it rather than asserted."""
    folds = {c: _fold(0.5 + 0.08 * i, 0.1 * i) for i, c in enumerate("abcdef")}
    r = s16.prior_gap_correlation(folds)
    assert r["critical_rho_n6_alpha05"] == 0.886
    assert r["significant"] is (abs(r["spearman"]) >= 0.886)
    assert "POST-HOC" in r["note"]


def test_prior_gap_correlation_declines_to_report_on_too_few_folds():
    assert s16.prior_gap_correlation({"a": _fold(0.5, 0.1)})["spearman"] is None


def test_a_fold_with_no_safe_cells_is_skipped_rather_than_producing_a_nan():
    folds = {"a": _fold(0.45, 0.10), "b": _fold(0.60, 0.30),
             "c": _fold(0.80, 0.60), "dead": _fold(0.9, float("nan"), n_safe=0)}
    r = s16.prior_gap_correlation(folds)
    assert r["n"] == 3 and "dead" not in r["folds"]


# ---- the recorded run ------------------------------------------------------------------------ #
@has_results
def test_all_six_folds_were_measured(recorded):
    assert set(recorded["folds"]) == {"N2", "N3", "O1", "O2", "Y1", "Y2"}
    assert recorded["errors"] == {}


@has_results
def test_the_gate_holds_most_held_out_cells_really_are_safe(recorded):
    """If this ever fails, the whole stage is void -- the rejections would be correct."""
    p = recorded["pooled"]["raw"]
    assert p["n_safe"] / p["n"] > 0.7
    assert recorded["folds"]["N2"]["raw"]["n_unsafe"] == 0
    assert recorded["folds"]["N2"]["raw"]["n_safe"] == 19


@has_results
def test_the_head_separates_the_classes_it_just_scores_the_safe_ones_low(recorded):
    """The finding rests on this: median S for truly-unsafe cells is far below truly-safe, so
    the head is not confused -- its safe class merely sits under the bar."""
    pl = recorded["pooled"]
    assert pl["median_S_true_unsafe"] < 0.4 < pl["median_S_true_safe"]
    assert pl["median_S_true_safe"] < pl["threshold"], "the safe class sits BELOW the bar"


@has_results
def test_a_majority_of_truly_safe_cells_fall_below_the_shipped_bar(recorded):
    assert recorded["pooled"]["frac_true_safe_below_threshold"] > 0.5


@has_results
def test_the_verdict_is_the_mechanical_result_of_the_pre_registered_rule(recorded):
    assert recorded["verdict"] == s16.verdict_from(recorded["pooled"])


@has_results
def test_the_oracle_is_undefined_on_n2_which_is_why_it_is_not_a_per_fold_upper_bound(recorded):
    """N2 is 19/19 safe, so a Platt fitted on it has an unidentifiable boundary and passes
    through unchanged -- its oracle column equals its raw column. The DEPLOYABLE calibrator is
    fitted on all six folds, because calib always carries both classes. That is a practical
    advantage of the deployable arm, not a defect of it."""
    assert recorded["folds"]["N2"]["platt_fitted_oracle"] is False
    assert recorded["folds"]["N2"]["raw"]["n_unsafe"] == 0
    assert (recorded["folds"]["N2"]["oracle"]["false_rejections"]
            == recorded["folds"]["N2"]["raw"]["false_rejections"])
    assert all(f["platt_fitted_deployable"] for f in recorded["folds"].values())


@has_results
def test_the_deployable_arm_matches_the_oracle_on_the_metric_that_matters(recorded):
    """H4 refuted: if calibration failed to transfer from calib to a held-out donor, the
    deployable arm would lag the oracle badly. It does not -- 26 vs 27 false rejections."""
    pl = recorded["pooled"]
    assert abs(pl["deployable"]["false_rejections"]
               - pl["oracle"]["false_rejections"]) <= 2
    assert pl["deployable_clears_bar"] is True


@has_results
def test_the_repair_is_not_bought_by_waving_unsafe_cells_through(recorded):
    """The trap. False approvals may rise, but only slightly, and balanced accuracy must
    improve -- otherwise the calibrator merely shifted the operating point."""
    pl = recorded["pooled"]
    assert pl["deployable"]["false_approvals"] <= pl["raw"]["false_approvals"] + 5
    assert pl["deployable"]["balanced_accuracy"] > pl["raw"]["balanced_accuracy"]
    assert pl["deployable"]["sensitivity"] > 2 * pl["raw"]["sensitivity"]


@has_results
def test_the_shipped_sensitivity_is_the_headline_failure(recorded):
    """PR-AUC 0.965-0.992, and the deployed verdict approves under a third of genuinely safe
    cells. The head's quality is not reaching the decision."""
    assert recorded["pooled"]["raw"]["sensitivity"] < 0.31


@has_results
def test_the_false_rejection_drop_is_computed_against_the_raw_arm(recorded):
    pl = recorded["pooled"]
    fr0 = pl["raw"]["false_rejections"]
    for arm in ("deployable", "oracle"):
        expect = (fr0 - pl[arm]["false_rejections"]) / fr0
        assert pl[f"{arm}_fr_drop"] == pytest.approx(expect)
        assert pl[f"{arm}_clears_bar"] == (expect >= s16.FALSE_REJECTION_DROP_BAR)


# ---- contract -------------------------------------------------------------------------------- #
def test_the_oracle_arm_is_labelled_everywhere_it_appears():
    """It is fitted on the cells it scores. If it is ever reported without that word attached,
    someone will read it as a deployable result."""
    src = SRC.read_text(encoding="utf-8")
    assert "ORACLE" in src
    assert "not deployable" in src.lower()


@has_results
def test_the_results_file_names_the_oracle_arm_explicitly(recorded):
    assert "oracle" in recorded["pooled"]
    assert "deployable" in recorded["pooled"]


def test_the_deployable_arm_never_sees_the_test_labels():
    """The no-leak assertion, read off the source: `platt(S_cal, truth_cal, S)` fits on calib."""
    src = SRC.read_text(encoding="utf-8")
    assert "platt(S_cal, truth_cal, S)" in src


def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_the_train_prior_is_read_without_loading_the_expression_matrix():
    """gather_split reads every shard's full 2000-dim X; the class prior needs only y_cls."""
    src = SRC.read_text(encoding="utf-8")
    assert 'columns=["cell_id", "y_cls"]' in src
    assert 'gather_split(paths, "holdout", "train")' not in src
