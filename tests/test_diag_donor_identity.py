"""Unit tests for REV FINAL §6.3's donor-identity check — pure functions only, no repo data.

The verdict logic gets the most attention: it is the difference between "the donor labels denote
the same people" and "every cross-series statement keyed on them is void". Both conditions of the
pre-registered bar are exercised independently, because the second one (margin separation) exists
precisely to stop the first from passing by chance at 1-in-9.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
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


D = _load("diag_donor_identity", "experiments/diag_donor_identity.py")


# ------------------------------------------------------------------ between_donor_f ---- #
def test_f_is_high_for_a_probe_that_separates_donors_and_is_stable_within_one():
    beta = np.array([[0.0, 0.0, 1.0, 1.0],        # perfect donor separator
                     [0.0, 1.0, 0.0, 1.0]])       # varies within donor, separates nothing
    f = D.between_donor_f(beta, ["A", "A", "B", "B"])
    assert f[0] > f[1] * 100


def test_f_is_zero_for_a_constant_probe_rather_than_dividing_by_zero():
    f = D.between_donor_f(np.full((1, 4), 0.5), ["A", "A", "B", "B"])
    assert np.isfinite(f[0]) and f[0] == pytest.approx(0.0)


# ------------------------------------------------------------------- trimodal_score ---- #
def test_a_genotype_shaped_probe_qualifies():
    """Betas at the three genotype modes, at least two of them occupied."""
    tri, n = D.trimodal_score(np.array([[0.02, 0.51, 0.98, 0.49]]))
    assert tri[0] and n[0] == 3


def test_an_invariant_probe_is_rejected_even_though_every_sample_sits_on_a_mode():
    """It satisfies the mode test trivially and fingerprints nothing — the trap the
    two-modes-occupied condition exists to catch."""
    tri, n = D.trimodal_score(np.array([[0.01, 0.02, 0.0, 0.01]]))
    assert not tri[0] and n[0] == 1


def test_a_continuously_varying_probe_is_rejected():
    """Real methylation variation sits between the modes; genotype variation does not."""
    tri, _ = D.trimodal_score(np.array([[0.30, 0.35, 0.65, 0.70]]))
    assert not tri[0]


def test_the_tolerance_is_honoured_in_both_directions():
    b = np.array([[0.0, 0.5 + 0.14, 1.0, 0.0]])
    assert D.trimodal_score(b, tol=0.15)[0][0]
    assert not D.trimodal_score(b, tol=0.10)[0][0]


# --------------------------------------------------------------------------- assign ---- #
def _profiles(seed=0, n=200):
    rng = np.random.default_rng(seed)
    return {d: rng.normal(size=n) for d in ("O1", "O2", "O3")}


def test_a_donor_matches_itself_when_the_panel_carries_identity():
    tgt = _profiles()
    qry = {d: v + np.random.default_rng(1).normal(0, 0.1, v.size) for d, v in tgt.items()}
    a = D.assign(qry, tgt)
    assert all(a[d]["correct"] for d in tgt)
    assert all(a[d]["margin"] > 0.3 for d in tgt)


def test_a_query_with_no_counterpart_is_flagged_and_never_counted_correct():
    tgt = {k: v for k, v in _profiles().items() if k != "O3"}
    qry = {"Y1": _profiles(seed=7)["O1"]}
    a = D.assign(qry, tgt)
    assert a["Y1"]["has_true_counterpart"] is False
    assert a["Y1"]["correct"] is False
    assert a["Y1"]["best"] in tgt          # it still reports a best match — that is the point


def test_margin_is_best_minus_runner_up():
    tgt = _profiles()
    a = D.assign({"O1": tgt["O1"]}, tgt)["O1"]
    ranked = sorted(a["correlations"].values(), reverse=True)
    assert a["margin"] == pytest.approx(ranked[0] - ranked[1])


# ------------------------------------------------------------------ identity_verdict ---- #
def _asg(correct, shared_margin, ctrl_margin):
    out = {}
    for q in ("O1", "O2"):
        out[q] = {"correct": correct, "has_true_counterpart": True, "margin": shared_margin}
    for q in ("Y1", "Y2"):
        out[q] = {"correct": False, "has_true_counterpart": False, "margin": ctrl_margin}
    return out


def test_correct_and_separated_is_the_only_full_pass():
    v = D.identity_verdict(_asg(True, 0.12, 0.06))
    assert v["status"] == "SAME_DONORS"


def test_correct_but_not_separated_from_the_controls_is_only_weak():
    """1-in-9 by chance is why condition 2 exists; without separation the pass is not evidence."""
    v = D.identity_verdict(_asg(True, 0.05, 0.09))
    assert v["status"] == "SAME_DONORS_WEAK"


def test_no_shared_label_matching_itself_voids_the_cross_series_keying():
    v = D.identity_verdict(_asg(False, 0.2, 0.01))
    assert v["status"] == "DIFFERENT_DONORS"
    assert "void" in v["reason"]


def test_a_partial_match_is_reported_as_inconsistent_not_rounded_either_way():
    a = _asg(True, 0.2, 0.01)
    a["O2"]["correct"] = False
    assert D.identity_verdict(a)["status"] == "INCONSISTENT"


def test_no_shared_labels_cannot_verify():
    assert D.identity_verdict({}, shared=())["status"] == "CANNOT_VERIFY"


# ---------------------------------------------------------------------- the null sims ---- #
def test_chance_of_getting_both_right_is_one_in_nine():
    assert float(D.sim_random_assignment(2, 3).mean()) == pytest.approx(1 / 9, abs=0.02)


def test_a_perfect_panel_loses_correlation_as_the_sample_count_drops():
    """The §5b check's premise: correlation between two noisy means of the SAME vector is
    bounded below 1, and the bound tightens as n falls."""
    hi = D.sim_stability(400, 0.32, 0.10, 7, 7, n_sim=200).mean()
    lo = D.sim_stability(400, 0.32, 0.10, 7, 2, n_sim=200).mean()
    assert 1.0 > hi > lo


# ------------------------------------------------------------------ the recorded run ---- #
def test_the_recorded_result_is_what_the_rules_produce():
    import json
    p = ROOT / "results" / "diag_donor_identity_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    r = json.loads(p.read_text(encoding="utf-8"))
    if "assignment" not in r:
        pytest.skip("run did not reach the assignment")
    assert D.identity_verdict(r["assignment"])["status"] == r["verdict"]["status"]
    # the finding: the shared donors separate from the no-counterpart controls
    assert min(r["verdict"]["shared_margins"]) > max(r["verdict"]["control_margins"])


def test_the_stability_gate_actually_fired_and_is_recorded():
    """The panel failing its own gate on reprogrammed cells is a result, not a mishap —
    it is why the test is restricted to non-reprogramming cells. Pin it."""
    import json
    p = ROOT / "results" / "diag_donor_identity_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    r = json.loads(p.read_text(encoding="utf-8"))
    if "panel_stability" not in r:
        pytest.skip("pre-registration only")
    assert r["panel_stability"]["verdict"] == "MOVES"
    assert r["stability_bar_check"]["verdict"] == "RESOLVABLE"   # the bar was NOT loosened
    assert min(r["diagnosis"]["stability_vs_failed_arm"].values()) >= 0.95
