"""Unit tests for STAGE 1.5.2 §12's anchor-reliability check — pure functions only, no repo data.

`STAGE_1_5_2_LABEL_ANCHOR.md` §10 requires "every pure function unit-tested with no repo data
present, per the pattern of the four existing `diag_*` scripts". The decision rules matter more
than the arithmetic here, because they are what turns four numbers into "M-2a's negative verdict
stands" or "it is withdrawn" — and a branch that never executes is not a check (the `verify_1a`
lesson).
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


R1 = _load("diag_r1_anchor_reliability", "experiments/diag_r1_anchor_reliability.py")
DMA = _load("diag_methylation_anchor", "experiments/diag_methylation_anchor.py")


# ----------------------------------------------------------------- lodo_age_errors ---- #
def test_lodo_recovers_a_perfect_clock_exactly():
    """If the linear predictors are exactly consistent with the ages, LODO error is 0."""
    ages = {"A": 30.0, "B": 50.0, "C": 70.0}
    k = -0.4471
    lp = {d: float(DMA.trafo(a)) - k for d, a in ages.items()}
    out = R1.lodo_age_errors(lp, ages, DMA.trafo, DMA.anti_trafo)
    assert out["n_folds"] == 3
    assert out["mae_years"] == pytest.approx(0.0, abs=1e-9)
    assert all(f["intercept"] == pytest.approx(k) for f in out["folds"])


def test_lodo_never_uses_the_held_out_donor_for_its_own_intercept():
    """The non-circularity claim, asserted rather than trusted."""
    ages = {"A": 30.0, "B": 50.0, "C": 70.0}
    lp = {"A": 0.1, "B": 0.5, "C": 0.9}
    out = R1.lodo_age_errors(lp, ages, DMA.trafo, DMA.anti_trafo)
    for f in out["folds"]:
        assert f["held_out"] not in f["intercept_from"]
        assert len(f["intercept_from"]) == 2


def test_lodo_error_grows_with_a_donor_level_offset():
    """A clock biased on one donor must show it on exactly that fold, not smear it away."""
    ages = {"A": 40.0, "B": 40.0, "C": 40.0}
    lp = {d: float(DMA.trafo(a)) for d, a in ages.items()}
    lp["C"] += 10.0 / 21.0                      # +10 yr on one donor
    out = R1.lodo_age_errors(lp, ages, DMA.trafo, DMA.anti_trafo)
    errs = {f["held_out"]: f["abs_err"] for f in out["folds"]}
    assert errs["C"] > errs["A"] and errs["C"] > errs["B"]


# -------------------------------------------------------------- intercept_free_gap ---- #
def test_gap_is_exactly_intercept_free():
    """Adding ANY constant to every linear predictor must not move the recovered gap."""
    old, young = [0.9, 1.1], [0.3]
    base = R1.intercept_free_gap(old, young)
    for k in (-5.0, -0.4471, 0.0, 0.696, 12.3):
        shifted = R1.intercept_free_gap([x + k for x in old], [y + k for y in young])
        assert shifted == pytest.approx(base, abs=1e-12)


def test_gap_uses_the_published_adult_slope():
    assert R1.intercept_free_gap([1.0], [0.0]) == pytest.approx(21.0)
    assert R1.ADULT_SLOPE == 21.0


# ------------------------------------------------------------------ the null models ---- #
def test_sim_lodo_mae_is_centred_where_the_analytic_value_says():
    """MAE of |N(0, sigma*sqrt(1.5))| is sigma*sqrt(1.5)*sqrt(2/pi) for 3 donors."""
    sim = R1.sim_lodo_mae(3, sigma=4.0, n_sim=20000)
    expected = 4.0 * np.sqrt(1.5) * np.sqrt(2.0 / np.pi)
    assert float(sim.mean()) == pytest.approx(expected, rel=0.03)


def test_a_zero_error_clock_passes_every_bar():
    """Sanity floor: the null must not be so wide that a perfect clock fails."""
    assert float(R1.sim_lodo_mae(3, sigma=1e-9, n_sim=200).max()) < 1e-6
    assert float(R1.sim_gap_abs_err(2, 1, sigma=1e-9, n_sim=200).max()) < 1e-6


# ------------------------------------------------------------------ the decision rules ---- #
def _r1a(ok):
    return {"verdict": "PASS" if ok else "FAIL", "detail": "d"}


def _r1b(ok):
    return {"verdict": "PASS" if ok else "FAIL", "detail": "d"}


def _r1d(ok):
    return {"verdict": "PASS" if ok else "FAIL", "detail": "d"}


@pytest.mark.parametrize("a,b,cov,expect", [
    (True, True, True, "ANCHOR_READS_AGE"),
    (False, True, True, "ANCHOR_QUESTIONABLE"),
    (True, False, True, "ANCHOR_QUESTIONABLE"),
    (True, True, False, "ANCHOR_QUESTIONABLE"),
])
def test_per_clock_decision_covers_every_branch(a, b, cov, expect):
    assert R1.r1_decide(_r1a(a), _r1b(b), cov)["action"] == expect


def test_r1d_failure_withdraws_the_verdict_even_when_everything_else_passes():
    """The load-bearing branch: an anchor that cannot arbitrate itself cannot arbitrate RNA."""
    per_clock = {"c1": {"decision": R1.r1_decide(_r1a(True), _r1b(True), True)}}
    out = R1.overall_decision(per_clock, _r1d(False))
    assert out["action"] == "M2A_NEGATIVE_VERDICT_WITHDRAWN"


def test_all_pass_accepts_the_verdict():
    per_clock = {"c1": {"decision": R1.r1_decide(_r1a(True), _r1b(True), True)},
                 "c2": {"decision": R1.r1_decide(_r1a(True), _r1b(True), True)}}
    assert R1.overall_decision(per_clock, _r1d(True))["action"] == "M2A_NEGATIVE_VERDICT_STANDS"


def test_a_single_questionable_clock_downgrades_to_the_caveat_branch():
    """One clock failing an under-powered check must not be silently absorbed."""
    per_clock = {"c1": {"decision": R1.r1_decide(_r1a(True), _r1b(True), True)},
                 "c2": {"decision": R1.r1_decide(_r1a(False), _r1b(True), True)}}
    out = R1.overall_decision(per_clock, _r1d(True))
    assert out["action"] == "M2A_NEGATIVE_VERDICT_STANDS_WITH_CAVEAT"


def test_r1d_reuses_m2a_bar_verbatim():
    """R1d's whole value is that it is the SAME criterion M-2a was graded on. Pin it."""
    assert R1.R1D_RHO_BAR == 0.50


def test_the_recorded_result_is_what_the_rules_produce():
    """Replay 2026-07-31's measured numbers through the pure rules, so a later edit to either
    the numbers or the rules that breaks their correspondence fails here."""
    import json
    p = ROOT / "results" / "diag_r1_anchor_reliability_results.json"
    if not p.exists():
        pytest.skip("results file not present")
    m = json.loads(p.read_text(encoding="utf-8"))["measurement"]
    assert R1.overall_decision(m["clocks"], m["R1d_interclock"])["action"] == m["overall"]
    # and the ceiling finding itself: meth<->meth barely clears the bar RNA was graded on
    assert 0.50 <= m["R1d_interclock"]["rho_partial"] < 0.70
