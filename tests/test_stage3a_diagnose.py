"""Unit tests for the Stage 3a diagnosis — pure functions only, no repo data.

Three things here are load-bearing and are pinned as PROPERTIES rather than as observations,
because each one is a claim the recorded 3a diagnosis rests on:

  * `FrozenRidge` must equal `sklearn.linear_model.Ridge`. D4 refits the same folds 2000 times
    per cell with only the target changing, which is only legitimate if the fast path IS the
    estimator it replaces. If it drifts, every pass-rate in D4 is measuring the wrong model.
  * Part B's Δt swing must be independent of `x0` for ANY linear model. The diagnosis says the
    five identical −269.13 rows are arithmetic, not agreement; that is a theorem about the
    estimator, so it is tested as one rather than by re-observing the number.
  * The bounded predictors must actually stay in [0,1], and the logit link must survive a target
    saturated at exactly 0 and exactly 1 — which is what ~90 % of the real target is.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


S = _load("stage3a_diagnose", "experiments/stage3a_diagnose.py")


def _design(n_tr=40, n_te=12, p=25, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n_tr, p)), rng.normal(size=(n_te, p))


# ------------------------------------------------------------------------ FrozenRidge ---- #
def test_frozen_ridge_equals_sklearn_ridge():
    """The fast path IS the estimator, not an approximation of it."""
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    Xtr, Xte = _design()
    y = np.random.default_rng(1).normal(size=len(Xtr))
    fr = S.FrozenRidge(Xtr, Xte, alpha=S.ALPHA)
    sc = StandardScaler().fit(Xtr)
    ref = Ridge(alpha=S.ALPHA).fit(sc.transform(Xtr), y).predict(sc.transform(Xte))
    assert np.abs(fr.predict(y) - ref).max() < 1e-8


def test_frozen_ridge_equals_sklearn_when_features_outnumber_samples():
    """The real geometry is p >> n (2000 genes, ~264 pairs); the dual form must hold there."""
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    Xtr, Xte = _design(n_tr=20, n_te=8, p=300, seed=5)
    y = np.random.default_rng(6).normal(size=len(Xtr))
    fr = S.FrozenRidge(Xtr, Xte, alpha=S.ALPHA)
    sc = StandardScaler().fit(Xtr)
    ref = Ridge(alpha=S.ALPHA).fit(sc.transform(Xtr), y).predict(sc.transform(Xte))
    assert np.abs(fr.predict(y) - ref).max() < 1e-7


def test_frozen_ridge_prediction_is_affine_in_the_target():
    """D4 reuses one H across thousands of targets; that is only valid if the map is affine."""
    Xtr, Xte = _design()
    rng = np.random.default_rng(2)
    a, b = rng.normal(size=len(Xtr)), rng.normal(size=len(Xtr))
    fr = S.FrozenRidge(Xtr, Xte)
    lhs = fr.predict(3.0 * a - 2.0 * b)
    rhs = 3.0 * fr.predict(a) - 2.0 * fr.predict(b)
    assert np.abs(lhs - rhs).max() < 1e-9


def test_frozen_ridge_reproduces_a_constant_target_exactly():
    """A constant target must come back as that constant -- the intercept path."""
    Xtr, Xte = _design()
    out = S.FrozenRidge(Xtr, Xte).predict(np.full(len(Xtr), 0.37))
    assert np.abs(out - 0.37).max() < 1e-9


def test_frozen_ridge_handles_a_constant_feature_column():
    """StandardScaler maps a zero-variance column to scale 1; it must not divide by zero."""
    Xtr, Xte = _design()
    Xtr[:, 3] = 2.0
    Xte[:, 3] = 2.0
    out = S.FrozenRidge(Xtr, Xte).predict(np.random.default_rng(3).normal(size=len(Xtr)))
    assert np.isfinite(out).all()


# ---------------------------------------------------------------- the Part B arithmetic ---- #
def test_the_dt_swing_is_independent_of_the_start_state_for_any_linear_model():
    """D3's claim as a theorem: `pred(x0, hi) - pred(x0, lo)` has no x0 term.

    This is why 3a's five folds printed −269.13 five times. It is not a property of that data;
    it holds for every x0 and every fitted linear model, so a per-fold Part B column produced
    this way can never disagree across folds.
    """
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    rng = np.random.default_rng(7)
    G = 15
    X = rng.normal(size=(60, G))
    dt = rng.uniform(0.2, 10.0, size=60)
    Z = np.hstack([X, dt.reshape(-1, 1), (dt ** 2).reshape(-1, 1)])
    y = 3.0 * dt + X[:, 0] + rng.normal(scale=0.1, size=60)
    sc = StandardScaler().fit(Z)
    m = Ridge(alpha=1.0).fit(sc.transform(Z), y)

    lo, hi = 0.5, 9.0
    swings = []
    for _ in range(6):
        x0 = rng.normal(size=G) * 5.0                     # wildly different start states
        p = [float(m.predict(sc.transform(np.hstack([x0, [q, q ** 2]]).reshape(1, -1)))[0])
             for q in (lo, hi)]
        swings.append(p[1] - p[0])
    assert np.std(swings) < 1e-9

    zl = (np.array([lo, lo ** 2]) - sc.mean_[G:]) / sc.scale_[G:]
    zh = (np.array([hi, hi ** 2]) - sc.mean_[G:]) / sc.scale_[G:]
    assert abs(swings[0] - float(m.coef_[G:] @ (zh - zl))) < 1e-9


# -------------------------------------------------------------------- bounded predictors ---- #
def test_logit_round_trips_values_strictly_inside_the_unit_interval():
    p = np.array([0.2, 0.5, 0.8])
    assert np.abs(S.expit(S.logit(p)) - p).max() < 1e-12


@pytest.mark.parametrize("v", [0.0, 1.0])
def test_logit_survives_a_target_saturated_at_the_bounds(v):
    """~90 % of the real Part C target is exactly 0 or exactly 1; logit(0) must not be -inf."""
    z = S.logit(np.array([v]))
    assert np.isfinite(z).all()
    assert abs(float(S.expit(z)[0]) - v) <= S.LOGIT_EPS + 1e-12


def test_expit_output_is_always_inside_the_unit_interval():
    z = np.array([-1e6, -50.0, 0.0, 50.0, 1e6])
    out = S.expit(z)
    assert out.min() >= 0.0 and out.max() <= 1.0


def test_expit_does_not_overflow_on_extreme_input():
    with np.errstate(over="raise"):
        assert np.isfinite(S.expit(np.array([-1e308, 1e308]))).all()


def test_mae_matches_the_definition():
    assert S.mae([1.0, 2.0, 3.0], [1.0, 4.0, 7.0]) == pytest.approx(2.0)


# --------------------------------------------------------------------------- lodo_folds ---- #
def _per(pairs_by_donor: dict[str, int]) -> dict:
    return {d: {"rows": [], "pairs": [{"dt": 1.0}] * n} for d, n in pairs_by_donor.items()}


def test_lodo_folds_keeps_every_fold_when_all_are_large_enough():
    assert S.lodo_folds(_per({"A": 10, "B": 10, "C": 10})) == ["A", "B", "C"]


def test_lodo_folds_drops_a_fold_with_too_few_held_out_pairs():
    """`>= 3 held-out pairs` is test18's own rule; the branch must actually fire."""
    assert S.lodo_folds(_per({"A": 10, "B": 10, "C": 2})) == ["A", "B"]


def test_lodo_folds_drops_a_fold_whose_training_side_is_too_small():
    """`>= 8 training pairs`: holding out A leaves 4, which is not enough."""
    assert S.lodo_folds(_per({"A": 20, "B": 4})) == ["B"]


def test_lodo_folds_returns_nothing_when_no_fold_qualifies():
    assert S.lodo_folds(_per({"A": 3, "B": 3})) == []


# ------------------------------------------------------------------------ oracle_on_tj ---- #
def _rows(ts, vals):
    return {"rows": [{"t": float(t), "u": float(v), "y": float(v), "n": 2}
                     for t, v in zip(ts, vals, strict=True)]}


def _ci(diffs):
    """test18's paired_ci contract, reduced to what the oracle needs."""
    d = [x for x in diffs if np.isfinite(x)]
    if len(d) < 2:
        return float("nan"), (float("nan"), float("nan")), len(d)
    md = float(np.mean(d))
    se = float(np.std(d, ddof=1)) / np.sqrt(len(d))
    return md, (md - 2.776 * se, md + 2.776 * se), len(d)


def test_oracle_is_perfect_when_every_donor_follows_the_same_step():
    """If the target is a shared function of t, predicting from the other donors is exact."""
    ts, step = [0.0, 1.0, 2.0, 3.0], [0.0, 0.0, 1.0, 1.0]
    per = {d: _rows(ts, step) for d in ("A", "B", "C")}
    r = S.oracle_on_tj(per, ["A", "B", "C"], "u", _ci)
    assert r["oracle_mean"] == pytest.approx(0.0)
    assert r["pooled_mean"] > 0.0
    assert r["verdict"] == "t_j HELPS"


def test_oracle_gains_nothing_when_the_target_is_pure_donor_identity():
    """A target constant within a donor and unrelated to t: t_j carries no information."""
    per = {"A": _rows([0.0, 1.0, 2.0, 3.0], [0.0] * 4),
           "B": _rows([0.0, 1.0, 2.0, 3.0], [1.0] * 4),
           "C": _rows([0.0, 1.0, 2.0, 3.0], [0.0] * 4)}
    r = S.oracle_on_tj(per, ["A", "B", "C"], "u", _ci)
    assert r["oracle_mean"] == pytest.approx(r["pooled_mean"])
    assert r["verdict"] == "tied"


def test_oracle_never_uses_the_held_out_donors_own_values():
    """Leakage check: changing ONLY the held-out donor's targets must not change its prediction.

    The oracle's whole standing as a ceiling depends on this -- a ceiling computed with the
    answer in hand is not a ceiling.
    """
    ts = [0.0, 1.0, 2.0, 3.0]
    base = {"A": _rows(ts, [0.0, 0.0, 1.0, 1.0]), "B": _rows(ts, [0.0, 0.0, 1.0, 1.0]),
            "C": _rows(ts, [0.0, 1.0, 0.0, 1.0])}
    alt = {k: _rows(ts, [r["u"] for r in v["rows"]]) for k, v in base.items()}
    alt["C"] = _rows(ts, [1.0, 1.0, 1.0, 1.0])
    a = S.oracle_on_tj(base, ["C"], "u", _ci)
    b = S.oracle_on_tj(alt, ["C"], "u", _ci)
    assert a["preds"]["C"] == b["preds"]["C"]
    # and the prediction really is the other donors' value at t_j -- A and B both step at t=2
    assert a["preds"]["C"] == [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def test_oracle_skips_a_fold_whose_target_is_not_finite():
    """`per_all` carries y = nan by construction; that fold must be skipped, not averaged in."""
    ts = [0.0, 1.0, 2.0, 3.0]
    per = {"A": _rows(ts, [0.0, 0.0, 1.0, 1.0]), "B": _rows(ts, [0.0, 0.0, 1.0, 1.0]),
           "C": _rows(ts, [np.nan] * 4)}
    r = S.oracle_on_tj(per, ["A", "B", "C"], "u", _ci)
    assert [row[0] for row in r["rows"]] == ["A", "B"]
