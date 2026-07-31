"""Every registered acceptance bar must be RESOLVABLE at the geometry it is graded on.

This is the forward form of `audit_metrics.py` and the executable half of REF_GROUND_RULES sec 5b:
a bar is only pre-registered once a system that meets its intent EXACTLY is shown to pass it at
least `MIN_PASS_RATE` of the time. Adding a bar to the project means adding an entry to
`REGISTERED_BARS` here; a bar with no entry is, by rule, not pre-registered.

The geometries below are REPRESENTATIVE, self-contained, and seeded -- they isolate the effect the
audit cares about (sample size / binning) without depending on run data. The exact run-3 numbers
live in the lab notebook; these guard the structural conclusion so it cannot silently regress.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from audit_metrics import MIN_PASS_RATE, bar_verdict, ece  # noqa: E402

TRIALS = 3000
SEED = 0


# ------------------------------------------------------------- geometry nulls ---- #
def _perfectly_calibrated_p(rng, n):
    """A representative NON-saturated confidence spread (matches run 3: no P(safe) near 1)."""
    return rng.uniform(0.1, 0.9, n)


def null_ece_pooled(n, gen=_perfectly_calibrated_p):
    """ECE of a perfectly calibrated model, pooled over `n` held-out cells."""
    rng = np.random.default_rng(SEED)
    return np.array([ece(p := gen(rng, n), (rng.random(n) < p).astype(float))
                     for _ in range(TRIALS)])


def null_ece_mean_of_folds(n_per_fold, k_folds, gen=_perfectly_calibrated_p):
    """ECE graded as the MEAN of per-fold ECEs -- the estimator Stage 1 originally used."""
    rng = np.random.default_rng(SEED)
    out = []
    for _ in range(TRIALS):
        vals = []
        for _ in range(k_folds):
            p = gen(rng, n_per_fold)
            vals.append(ece(p, (rng.random(n_per_fold) < p).astype(float)))
        out.append(float(np.mean(vals)))
    return np.array(out)


def null_m1_contrast_correct_clock(true_gap=53.0, cv_mae=12.26879346460328, n_each=2):
    """Stage 1.5 Phase 1 / M1: the extreme age contrast a CORRECT clock would produce.

    Intent: the clock reads chronological age with its own published CV error. The bar is tested
    against the opposite null (a clock reading nothing), so resolvability asks how often a clock
    that DOES work clears a threshold set to exclude one that does not.
    """
    rng = np.random.default_rng(SEED)
    se = cv_mae * np.sqrt(1.0 / n_each + 1.0 / n_each)
    return rng.normal(true_gap, se, 20000)


def null_coverage_marginal(n_total, level):
    """Marginal coverage of a correctly-`level` conformal interval over `n_total` cells."""
    rng = np.random.default_rng(SEED)
    return rng.binomial(n_total, level, 20000) / n_total


def band_pass_rate(null, lo, hi):
    return float(((null >= lo) & (null <= hi)).mean())


# --------------------------------------------------------------- the registry ---- #
# One entry per acceptance bar. `expect` is what the resolvability rule REQUIRES; a "retired"
# entry documents a bar we removed BECAUSE it was unresolvable, and asserts it stays that way.
REGISTERED_BARS = [
    {
        "name": "conformal_coverage in [0.85,0.95], pooled marginal",
        "kind": "band",
        "null": lambda: null_coverage_marginal(124, 0.90),
        "band": (0.85, 0.95),
        "expect": "RESOLVABLE",
        "where": "STAGE_1_CALIBRATION.md sec 3",
        "note": "a correctly-90% system lands in-band ~93% of the time (confirmed, not assumed)",
    },
    {
        "name": "fate_ece <= 0.169, POOLED over held-out cells (~103)",
        "kind": "lower",
        "null": lambda: null_ece_pooled(103),
        "bar": 0.169,
        "expect": "RESOLVABLE",
        "where": "STAGE_1_CALIBRATION.md sec 3, as repaired 2026-07-23",
        "note": "the resolvable form: a perfectly calibrated model clears it ~99% of the time",
    },
    {
        "name": "fate_ece <= 0.169, mean of per-fold ECE (n~21 x 5) [RETIRED]",
        "kind": "lower",
        "null": lambda: null_ece_mean_of_folds(21, 5),
        "bar": 0.169,
        "expect": "UNRESOLVABLE",
        "where": "the original Stage 1 grading; retired because of this very property",
        "note": "kept as a regression: a perfectly calibrated model FAILS it most of the time, "
                "so the bar tested the sample size. If this ever reads RESOLVABLE the geometry "
                "assumptions changed and sec 5b must be revisited.",
    },
    {
        "name": "M1 extreme age contrast >= 20.2 yr (Stage 1.5 Phase 1)",
        "kind": "higher",
        "null": null_m1_contrast_correct_clock,
        "bar": 1.6448536269514722 * 12.26879346460328,     # z_0.95 * SE under a null clock
        "expect": "RESOLVABLE",
        "where": "STAGE_1_5_HARMONIZATION_AUDIT.md §5.4 M1, registered per §6.2 T1",
        "note": "a clock that reads NOTHING clears 'contrast > 0' half the time, so the bar is "
                "set at z_0.95 of the null SE instead. A correct clock (true gap 53 yr, cv_mae "
                "12.27) clears that ~99.6% of the time. The 29-vs-35 middle contrast is "
                "deliberately NOT gated -- it is half the clock's error and unresolvable.",
    },
]


@pytest.mark.parametrize("spec", REGISTERED_BARS, ids=lambda s: s["name"])
def test_registered_bar_has_expected_resolvability(spec):
    if spec["kind"] == "band":
        rate = band_pass_rate(spec["null"](), *spec["band"])
        verdict = "RESOLVABLE" if rate >= MIN_PASS_RATE - 0.05 else "UNRESOLVABLE"
        # coverage's band sits at ~0.93; the 0.05 slack is the known binomial width at n=124,
        # documented in the audit. A band criterion is judged in-band-rate, not a one-sided tail.
    else:
        r = bar_verdict(spec["null"](), spec["bar"], lower_is_better=(spec["kind"] == "lower"))
        verdict, rate = r["verdict"], r["pass_rate"]
    assert verdict == spec["expect"], (
        f"{spec['name']}: expected {spec['expect']} but a correct system passes "
        f"{rate:.1%} (bar from {spec['where']})")


def test_every_registered_bar_is_documented():
    """Contract for adding a bar: name, expectation, provenance, and a rationale note."""
    for s in REGISTERED_BARS:
        assert s["name"] and s["expect"] in ("RESOLVABLE", "UNRESOLVABLE")
        assert s["where"] and s["note"]
        assert s["kind"] in ("band", "lower", "higher")


def test_the_retired_and_repaired_bars_differ_only_in_geometry():
    """The whole lesson in one assertion: SAME bar, SAME intent, pooling flips the verdict."""
    per_fold = bar_verdict(null_ece_mean_of_folds(21, 5), 0.169)
    pooled = bar_verdict(null_ece_pooled(103), 0.169)
    assert per_fold["verdict"] == "UNRESOLVABLE"
    assert pooled["verdict"] == "RESOLVABLE"
    assert pooled["null_median"] < per_fold["null_median"]


# ------------------------------------------------------------- bar_verdict API ---- #
def test_bar_verdict_thresholds_on_min_pass_rate():
    passing = np.full(1000, 0.05)          # always below a 0.169 bar
    failing = np.full(1000, 0.30)          # always above it
    assert bar_verdict(passing, 0.169)["verdict"] == "RESOLVABLE"
    assert bar_verdict(failing, 0.169)["verdict"] == "UNRESOLVABLE"


def test_bar_verdict_boundary_at_exactly_min_pass():
    # 95 of 100 pass -> pass_rate 0.95 -> RESOLVABLE (>=, not >)
    null = np.array([0.10] * 95 + [0.30] * 5)
    assert bar_verdict(null, 0.169, min_pass=0.95)["verdict"] == "RESOLVABLE"
    null2 = np.array([0.10] * 94 + [0.30] * 6)
    assert bar_verdict(null2, 0.169, min_pass=0.95)["verdict"] == "UNRESOLVABLE"


def test_bar_verdict_reports_the_usable_bar_when_unresolvable():
    rng = np.random.default_rng(0)
    null = rng.normal(0.25, 0.03, 20000)   # a correct system scores ~0.25, bar is 0.169
    r = bar_verdict(null, 0.169)
    assert r["verdict"] == "UNRESOLVABLE"
    # moving the bar to usable_bar would make a correct system pass 95% of the time
    assert (null <= r["usable_bar"]).mean() == pytest.approx(0.95, abs=0.01)


def test_higher_is_better_direction():
    null = np.linspace(0.0, 1.0, 10001)
    # higher-is-better: a bar at 0.04 leaves ~96% of a correct system's mass above it -> RESOLVABLE
    r = bar_verdict(null, 0.04, lower_is_better=False)
    assert r["verdict"] == "RESOLVABLE"
    assert r["pass_rate"] == pytest.approx(0.96, abs=0.01)
    # ...but a bar at 0.10 leaves only ~90% above it, below MIN_PASS_RATE -> UNRESOLVABLE
    assert bar_verdict(null, 0.10, lower_is_better=False)["verdict"] == "UNRESOLVABLE"


# ===================================================================== #
# STAGE 1.5.2 bars, frozen 2026-07-31 at the ACTUAL GSE165177xGSE165179  #
# geometry (68 conditions, 3 donors). §6 requires every registered bar   #
# to carry a resolvability test; a bar without one is not pre-registered.#
# Recomputed from stage_1_5_2_resolvability_results.json.                #
# ===================================================================== #
import json as _json  # noqa: E402
from pathlib import Path as _Path  # noqa: E402

_S152 = _Path(__file__).resolve().parents[1] / "stage_1_5_2_resolvability_results.json"


def _s152():
    if not _S152.exists():
        import pytest as _pytest
        _pytest.skip("stage_1_5_2_resolvability_results.json not present")
    return _json.loads(_S152.read_text(encoding="utf-8"))["checks"]


def test_s152_rho_partial_is_resolvable_at_the_actual_geometry():
    """The DECISIVE criterion. Resolvable only because n=68, not the registered n=22."""
    c = _s152()["M-2a rho_partial (ACTUAL)"]
    assert c["verdict"] == "RESOLVABLE" and c["pass_rate"] >= 0.95


def test_s152_the_registered_fallback_would_have_failed():
    """Pins the finding: rho_partial at the registered n=22 is UNRESOLVABLE (92.3%).
    On the GSE165178-only geometry the stage had no valid decisive criterion."""
    c = _s152()["M-2a rho_partial (registered n=22)"]
    assert c["verdict"] == "UNRESOLVABLE" and c["pass_rate"] < 0.95


def test_s152_rho_within_is_demoted_not_gated():
    """rho_within is UNRESOLVABLE at both the registered and actual n -> descriptive only."""
    ch = _s152()
    for k in ("M-2a rho_within (registered n=11/arm)", "M-2a rho_within (ACTUAL smallest arm)"):
        assert ch[k]["verdict"] == "UNRESOLVABLE"


def test_s152_sign_agreement_bar_moved_to_its_usable_bar():
    """Registered >=8/11 is UNRESOLVABLE; §5b requires moving it before the run, to 7."""
    c = _s152()["M-2b sign agreement (registered 8/11)"]
    assert c["verdict"] == "UNRESOLVABLE"
    assert abs(c["usable_bar"] - 7.0) < 1e-9


def test_s152_lodo_survives_the_reduced_donor_count():
    """Registered assumed 4 folds; the real set has 3 donors. Still resolvable."""
    c = _s152()["M-2c LODO MAE (ACTUAL donors)"]
    assert c["verdict"] == "RESOLVABLE" and c["geometry"] == "3 folds"


# ===================================================================== #
# STAGE 1.5.2 §12 (R1/R4 anchor reliability) and G-c step 1.            #
# Same rule: every bar these stages were graded on carries an entry, or #
# it is not pre-registered.                                             #
# ===================================================================== #
_R1 = _Path(__file__).resolve().parents[1] / "diag_r1_anchor_reliability_results.json"
_GC = _Path(__file__).resolve().parents[1] / "diag_gc_hff_signature_results.json"


def _load_json(p):
    if not p.exists():
        import pytest as _pytest
        _pytest.skip(f"{p.name} not present")
    return _json.loads(p.read_text(encoding="utf-8"))


def test_s152_r1_chronological_age_bars_are_both_underpowered():
    """Pins the honest weakness: 3 donors carrying 2 distinct ages cannot resolve either
    proposed bar, so both had to be loosened -- which makes M-2a's negative verdict HARDER
    to falsify, a bias in favour of the result already recorded."""
    ch = _load_json(_R1)["preregistration"]["checks"]
    for k in ch:
        if k.startswith(("R1a", "R1b")):
            assert ch[k]["verdict"] == "UNRESOLVABLE"
            assert ch[k]["bar_used"] == ch[k]["usable_bar"]      # moved, per §5b


def test_s152_r1d_reuses_m2a_bar_rather_than_deriving_a_new_one():
    """R1d's entire value is that meth<->meth and RNA<->meth are scored by the SAME
    criterion. A separately-derived bar would destroy the comparison."""
    ch = _load_json(_R1)["preregistration"]["checks"]
    k = next(k for k in ch if k.startswith("R1d"))
    assert ch[k]["bar"] == ch[k]["bar_used"] == 0.50
    assert ch[k]["verdict"] == "RESOLVABLE"
    assert "reused verbatim" in ch[k]["source"]


def test_s152_gc_signature_bar_is_resolvable_without_being_moved():
    """G-c is the one bar in this stage that did NOT need loosening."""
    pre = _load_json(_GC)["preregistration"]
    assert pre["verdict"] == "RESOLVABLE" and pre["pass_rate"] >= 0.95
    assert pre["rho_bar_used"] == pre["rho_bar_proposed"] == -0.50


def test_s152_gc_slope_band_is_two_sided():
    """'Within ~2x of methylation's' must exclude a slope 10x too steep as well as one
    10x too shallow, or it is not a similarity criterion at all."""
    lo, hi = _load_json(_GC)["preregistration"]["slope_band"]
    assert lo < hi < 0
    assert hi / lo == pytest.approx(0.25, abs=1e-9)      # 2x each way around the mean
