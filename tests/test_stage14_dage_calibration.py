"""Stage 14 -- the ΔAge scale correction, at the reporting boundary only.

The ways this could be wrong, each pinned:

* it could silently REPLACE the raw number instead of sitting beside it, hiding an untested
  cross-cohort assumption behind a single figure;
* it could leak upstream into `res.py`, where `kappa` is a half-saturation IN YEARS, silently
  reinterpreting the RES formula -- the same class of defect as the Huber knee that ruled out
  rescaling the training target;
* it could reach the scorecard, where a 1/k drop in ΔAge MAE would read as a large ACCEPT for a
  change that improved nothing;
* the flattering factor (`k_LS`, which beats the instrument floor by under-reporting magnitude
  40 %) could be shipped instead of the honest one.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from cellfate.inference import dage_calibration as dc

ROOT = Path(__file__).resolve().parents[1]


# ---- the factor, and its provenance --------------------------------------------------------- #
def test_the_shipped_factor_is_the_variance_matched_one():
    assert dc.K_VAR == pytest.approx(0.5991, abs=1e-4)
    assert dc.K_LS < dc.K_VAR, "least squares is strictly smaller -- it shrinks"


def test_the_factors_match_stage_11s_recorded_lodo_fits():
    """Neither number may drift from the measurement it came from."""
    import numpy as np
    r = json.loads((ROOT / "results" / "diag_stage11_scale_results.json").read_text("utf-8"))
    raw = r["variants"]["raw"]
    assert np.mean(list(raw["k_var_per_donor"].values())) == pytest.approx(dc.K_VAR, abs=1e-3)
    assert np.mean(list(raw["k_per_donor"].values())) == pytest.approx(dc.K_LS, abs=1e-3)


def test_the_rejected_alternative_is_the_one_that_shrinks_the_spread():
    """k_LS reaches MAE 6.78 against a 7.30 floor -- the quotable headline -- but SD ratio 0.597.
    For a REPORTING transform the objective is unbiased magnitude, so the worse MAE ships."""
    r = json.loads((ROOT / "results" / "diag_stage11_scale_results.json").read_text("utf-8"))
    raw = r["variants"]["raw"]
    assert raw["sd_ratio_scaled"] < 0.7, "k_LS under-reports magnitude"
    assert 0.95 < raw["sd_ratio_var_matched"] < 1.05, "k_var preserves it"
    assert raw["mae_scaled"] < raw["mae_var_matched"], "and k_LS does win on MAE"


# ---- the arithmetic ------------------------------------------------------------------------- #
def test_calibrate_scales_and_preserves_sign():
    assert dc.calibrate(10.0) == pytest.approx(5.991)
    assert dc.calibrate(-10.0) == pytest.approx(-5.991)
    assert dc.calibrate(0.0) == 0.0


def test_calibrating_an_interval_scales_both_ends_so_it_stays_coherent():
    lo, hi = dc.calibrate_interval(-30.0, 10.0)
    assert lo == pytest.approx(-30.0 * dc.K_VAR)
    assert hi == pytest.approx(10.0 * dc.K_VAR)
    assert lo < hi, "the interval must not invert"


def test_a_positive_rescale_cannot_reorder_anything():
    """The property that makes this safe to ship: every ranking result stands unchanged."""
    vals = [-30.0, -12.0, 0.0, 4.0, 19.0]
    cal = [dc.calibrate(v) for v in vals]
    assert sorted(range(len(vals)), key=lambda i: vals[i]) == \
        sorted(range(len(cal)), key=lambda i: cal[i])


def test_the_caveat_states_the_untested_transfer_and_names_the_factor():
    assert "UNTESTED" in dc.CAVEAT
    assert "DISJOINT" in dc.CAVEAT
    assert "0.5991" in dc.CAVEAT
    assert "alongside" in dc.CAVEAT


# ---- containment: it must not leak anywhere it would change a decision ---------------------- #
def test_res_is_computed_on_raw_values_not_calibrated_ones():
    """`kappa = 5.0` is a rejuvenation half-saturation IN YEARS. Calibrating upstream of RES
    would silently reinterpret it -- exactly the defect that ruled out rescaling the target."""
    svc = (ROOT / "src" / "cellfate" / "inference" / "service.py").read_text(encoding="utf-8")
    call = svc[svc.index("compute_res("):svc.index("lo, hi = interval")]
    assert "calibrate" not in call, "RES must see raw mu_age/sigma_age"
    assert 's["mu_age"], s["sigma_age"]' in call


def test_the_predictor_emits_raw_values_only():
    """If calibration moved into the Predictor, RES, the scorecard and every evaluation path
    would inherit it silently."""
    src = (ROOT / "src" / "cellfate" / "inference" / "predictor.py").read_text(encoding="utf-8")
    assert "dage_calibration" not in src
    assert "K_VAR" not in src


def test_the_scorecard_is_untouched_by_calibration():
    """A 1/k drop in ΔAge MAE would print as a large ACCEPT for a change that improved nothing."""
    src = (ROOT / "scorecard.py").read_text(encoding="utf-8")
    assert "dage_calibration" not in src and "K_VAR" not in src


def test_calibration_lives_in_exactly_one_place():
    """One factor, one site -- so it stays trivially reversible."""
    hits = [p for p in (ROOT / "src").rglob("*.py")
            if "dage_calibration" in p.read_text(encoding="utf-8") and p.name != "dage_calibration.py"]
    assert [p.name for p in hits] == ["service.py"]


# ---- the response contract ------------------------------------------------------------------ #
def test_the_response_carries_both_numbers_never_one():
    from cellfate.inference.schema import Response
    f = Response.model_fields
    for name in ("delta_age_mean", "delta_age_interval", "epistemic_std",
                 "delta_age_calibrated", "delta_age_interval_calibrated",
                 "epistemic_std_calibrated", "delta_age_calibration_k"):
        assert name in f, name


def test_the_calibrated_fields_default_to_none_so_old_responses_stay_valid():
    from cellfate.inference.schema import Response
    r = Response(status="APPROVED", rejuvenation_efficacy_score=0.0, p_identity_preserved=0.9,
                 p_identity_loss=0.05, p_apoptosis=0.05, delta_age_mean=-3.0,
                 delta_age_interval=[-10.0, 4.0], in_distribution=True, epistemic_std=2.0,
                 predictive_entropy=0.4)
    assert r.delta_age_calibrated is None
    assert r.delta_age_calibration_k is None


def test_the_raw_field_is_still_the_raw_number():
    """The load-bearing guarantee: `delta_age_mean` must NOT quietly become calibrated."""
    svc = (ROOT / "src" / "cellfate" / "inference" / "service.py").read_text(encoding="utf-8")
    assert 'delta_age_mean=round(s["mu_age"], 2),' in svc
    assert 'delta_age_calibrated=round(calibrate(s["mu_age"]), 2),' in svc


def test_the_existing_unvalidated_age_warning_is_not_weakened():
    """Calibration makes the number better scaled, not validated. The pre-existing caveat about
    absolute interpretation must survive intact alongside the new one."""
    svc = (ROOT / "src" / "cellfate" / "inference" / "service.py").read_text(encoding="utf-8")
    assert "do not read the number as years" in svc
    assert "DAGE_CAVEAT" in svc


def test_a_built_response_carries_the_calibrated_values_and_the_caveat():
    """End to end on a synthetic summary, so the wiring is exercised rather than asserted."""
    from cellfate.inference.service import build_response

    class _P:
        q = 20.0
        res_params = __import__("cellfate.common.schemas", fromlist=["ResParams"]).ResParams()
        age_provenance = None

    s = {"S": 0.9, "P_loss": 0.05, "P_death": 0.05, "mu_age": -10.0, "sigma_age": 4.0,
         "in_dist": True, "entropy": 0.3}
    r = build_response(_P(), s)
    assert r.delta_age_mean == pytest.approx(-10.0)
    assert r.delta_age_calibrated == pytest.approx(round(-10.0 * dc.K_VAR, 2))
    assert r.delta_age_interval == [-30.0, 10.0]
    assert r.delta_age_interval_calibrated == [round(-30.0 * dc.K_VAR, 2),
                                               round(10.0 * dc.K_VAR, 2)]
    assert r.delta_age_calibration_k == pytest.approx(dc.K_VAR)
    assert "UNTESTED" in r.warning


def test_the_calibrated_interval_still_brackets_the_calibrated_point():
    from cellfate.inference.service import build_response

    class _P:
        q = 20.0
        res_params = __import__("cellfate.common.schemas", fromlist=["ResParams"]).ResParams()
        age_provenance = None

    for mu in (-25.0, 0.0, 18.0):
        r = build_response(_P(), {"S": 0.9, "P_loss": 0.05, "P_death": 0.05, "mu_age": mu,
                                  "sigma_age": 3.0, "in_dist": True, "entropy": 0.3})
        lo, hi = r.delta_age_interval_calibrated
        assert lo <= r.delta_age_calibrated <= hi
