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
# --------------------------------------------------------------------------- #
# STAGE 1.5.3 step 6 / G-c step 2 -- the ARM COMPARISON.                       #
# --------------------------------------------------------------------------- #
def null_gc_step2_paired(sd: float, delta: float = 3.5723, n_folds: int = 6,
                         n_sim: int = 20_000, seed: int = 0) -> np.ndarray:
    """|t| for a paired `n_folds`-fold comparison whose TRUE effect is `delta`.

    `delta` defaults to DELTA* -- Stage 2 sec 12's registered ">=25% drop in dage_mae_model"
    applied to the 14.29 yr recorded baseline. Graded against the two-sided t critical value,
    which is exactly the outcome table's "CI excludes 0" rule.

    Registered because step 6 decides whether 99.7% of the age labels are discarded and its
    criterion had NO bar. See plan_tests/register_gc_step2_bar.py.
    """
    rng = np.random.default_rng(seed)
    d = rng.normal(delta, sd, size=(n_sim, n_folds))
    se = d.std(axis=1, ddof=1) / np.sqrt(n_folds)
    return np.abs(d.mean(axis=1)) / np.maximum(se, 1e-12)


_T_CRIT_5 = 2.5705818366147395      # t(.975, df=5)

# HFF's measured unsafe fraction by timepoint (day 0 -> day 21), GSE242423 pooled over 42,481
# cells at ~4,700 per timepoint. Recorded in STAGE_3_TOOL.md "3a-bis". Used as the shape of a
# REAL forward safety curve so the 3a bar is simulated against something the data actually shows
# rather than an invented effect size.
_HFF_UNSAFE_CURVE = np.array([0.0835, 0.3509, 0.3508, 0.4226, 0.4708, 0.3964, 0.4659, 0.6791,
                              0.9996])


def null_forward_gate_paired(cells_per_tp: int, n_folds: int, alpha: float = 1.0):
    """3a's paired forward comparison when the Δt signal is REAL, at a given held-out precision.

    Both arms are given the best predictor their information set allows -- the Δt arm knows
    `p(t_j)` exactly, the state arm the best value available from `t_i` alone -- so this is a
    system that meets 3a's intent exactly. The only thing that varies is the binomial noise on
    the held-out target at `cells_per_tp` cells. Returns the paired 95 % CI upper end, which 3a
    requires to be < 0.

    This isolates the PRECISION axis only. `STAGE_3_TOOL.md` "3a-bis" measured the other axis --
    the held-out line/modality shift -- on the real data, and THAT is what makes the Gill-donor
    geometry unresolvable at any cell count. It is not synthesisable here and is not claimed to be.
    """
    from scipy.stats import t as _t

    rng = np.random.default_rng(SEED)
    p = np.clip(_HFF_UNSAFE_CURVE.mean() + alpha * (_HFF_UNSAFE_CURVE - _HFF_UNSAFE_CURVE.mean()),
                0.01, 0.99)
    T = len(p)
    ii = np.array([i for i in range(T) for j in range(T) if j > i])
    jj = np.array([j for i in range(T) for j in range(T) if j > i])
    blind = np.array([p[jj[ii == i]].mean() if (ii == i).any() else p.mean() for i in range(T)])
    tc = float(_t.ppf(0.975, n_folds - 1))
    out = []
    for _ in range(TRIALS):
        d = []
        for _f in range(n_folds):
            y = rng.binomial(cells_per_tp, p)[jj] / cells_per_tp
            d.append(float(np.abs(p[jj] - y).mean()) - float(np.abs(blind[ii] - y).mean()))
        d = np.asarray(d)
        out.append(d.mean() + tc * d.std(ddof=1) / np.sqrt(n_folds))
    return np.asarray(out)


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
    {
        "name": "G-c step 2 arm comparison, DELTA*=3.57 yr at SD(diff)=1.0 (step 6)",
        "kind": "higher",
        "null": lambda: null_gc_step2_paired(sd=1.0),
        "bar": _T_CRIT_5,
        "expect": "RESOLVABLE",
        "where": "plan_tests/register_gc_step2_bar.py; STAGE_1_5_3_EXECUTE.md step 6",
        "note": "the ONLY regime in which step 6 can detect an effect worth acting on. The arms "
                "must track each other to within ~1 yr per fold for the paired CI to resolve "
                "DELTA*.",
    },
    {
        "name": "G-c step 2 arm comparison at SD(diff)=3.0 -- UNDERPOWERED (step 6)",
        "kind": "higher",
        "null": lambda: null_gc_step2_paired(sd=3.0),
        "bar": _T_CRIT_5,
        "expect": "UNRESOLVABLE",
        "where": "plan_tests/register_gc_step2_bar.py",
        "note": "a REAL effect of DELTA* is detected only ~65% of the time here, so a null is not "
                "evidence of absence. Registered so the underpowered regime is a recorded fact "
                "rather than something rediscovered after the run -- the fate_ece lesson.",
    },
    {
        "name": "G-c step 2 arm comparison at SD(diff)=13.7 (arms independent) (step 6)",
        "kind": "higher",
        "null": lambda: null_gc_step2_paired(sd=13.7),
        "bar": _T_CRIT_5,
        "expect": "UNRESOLVABLE",
        "where": "plan_tests/register_gc_step2_bar.py",
        "note": "sqrt(2) x the 9.67 yr baseline fold SD -- the pessimistic bound if pairing "
                "cancels nothing. Detection collapses to ~7.5%, barely above the 5% false-positive "
                "rate: the test would be almost pure noise.",
    },
    {
        "name": "3a forward gate, held-out target on ~470 cells/tp, full HFF curve (10 folds)",
        "kind": "lower",
        "null": lambda: null_forward_gate_paired(472, 10, alpha=1.0),
        "bar": 0.0,
        "expect": "RESOLVABLE",
        "where": "STAGE_3_TOOL.md '3a-bis'; experiments/stage3a_bis_resolvability.py regime B",
        "note": "the ONLY held-out geometry in which 3a's rule can register its own data's "
                "forward curve. Measured on the real folds at 1.000. Registered so that any "
                "future 3a run states which side of this line it sits on.",
    },
    {
        "name": "3a forward gate, held-out target on 2 cells/tp, HALF-amplitude curve (10 folds)",
        "kind": "lower",
        "null": lambda: null_forward_gate_paired(2, 10, alpha=0.5),
        "bar": 0.0,
        "expect": "UNRESOLVABLE",
        "where": "STAGE_3_TOOL.md '3a-bis'; regime C, the precision counterfactual",
        "note": "the SENSITIVITY FLOOR. At 2 cells/timepoint a real effect of half HFF's measured "
                "amplitude is missed almost always, while the same effect at ~470 cells/tp is "
                "caught every time. NOTE what this does NOT show: at FULL amplitude 2 cells/tp "
                "still resolves (regime C, 0.965), so held-out precision is NOT what makes the "
                "Gill-donor geometry unresolvable -- the held-out line/modality shift is, and "
                "that was measured on the real folds (regime D: 0.000 even at 472 cells/tp). "
                "Kept as the regression for the precision axis only.",
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

_S152 = _Path(__file__).resolve().parents[1] / "results" / "stage_1_5_2_resolvability_results.json"


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
_R1 = _Path(__file__).resolve().parents[1] / "results" / "diag_r1_anchor_reliability_results.json"
_GC = _Path(__file__).resolve().parents[1] / "results" / "diag_gc_hff_signature_results.json"


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


# ===================================================================== #
# STAGE 1.5.3 C-5 — the age head's per-update occupancy bars.           #
# Registered BEFORE any retrain, which is what step 5 requires: the bar #
# grades the MECHANISM, because the outcome metric needs step 6's run.  #
# ===================================================================== #
_C5 = _Path(__file__).resolve().parents[1] / "results" / "register_c5_bar_results.json"


def test_s153_c5_bars_are_resolvable_on_the_dense_regime():
    """The system that meets the intent exactly is today's DENSE regime -- before masking,
    every cell carries an age label. If THAT failed, the bars would be measuring the sample
    size rather than the mechanism."""
    r = _load_json(_C5)["reference_dense"]
    assert r["B1"]["verdict"] == "RESOLVABLE" and r["B1"]["pass_rate"] == 1.0
    assert r["B2"]["verdict"] == "RESOLVABLE" and r["B2"]["pass_rate"] == 1.0


def test_s153_c5_bars_actually_discriminate_between_the_options():
    """A bar every candidate passes decides nothing. This is the check that makes it a bar
    rather than a formality, and the registration script exits non-zero without it."""
    assert _load_json(_C5)["discriminates"] is True


def test_s153_the_status_quo_fails_both_bars():
    """Pins C-5's premise: 75 labels among 33,688 cells at batch 512 is 1.14 per update, so
    ~31% of updates contribute NOTHING and almost none clears 4 cells."""
    c = _load_json(_C5)["candidates"]["status quo (uniform shuffling)"]
    assert c["B1"]["verdict"] == "UNRESOLVABLE" and c["B1"]["pass_rate"] < 0.75
    assert c["B2"]["verdict"] == "UNRESOLVABLE" and c["B2"]["pass_rate"] < 0.10
    assert c["mean_cells"] < 1.3


def test_s153_option_3_is_indistinguishable_from_the_status_quo():
    """Pinning `s_age` does nothing about occupancy -- exactly C-5's criticism of it, now
    measured rather than asserted."""
    cs = _load_json(_C5)["candidates"]
    quo = cs["status quo (uniform shuffling)"]
    opt3 = cs["Option 3 (pin s_age only)"]
    assert opt3["B1"]["verdict"] == quo["B1"]["verdict"] == "UNRESOLVABLE"
    assert abs(opt3["mean_cells"] - quo["mean_cells"]) < 0.2


def test_s153_both_surviving_options_clear_both_bars():
    cs = _load_json(_C5)["candidates"]
    survivors = [k for k, v in cs.items()
                 if v["B1"]["verdict"] == v["B2"]["verdict"] == "RESOLVABLE"]
    assert len(survivors) == 2
    assert any("Option 1" in s for s in survivors) and any("Option 2" in s for s in survivors)


def test_s153_option_2_scores_higher_and_costs_the_fate_task_nothing():
    """The finding that overturned the plan's original recommendation of Option 1.

    Option 2 (accumulate) clears B2 by more than Option 1 (sampler) AND changes no sampling,
    so the fate head's training distribution is untouched. Option 1 oversamples the 75 age
    cells ~7x, which the plan itself flagged as 'not free'."""
    d = _load_json(_C5)
    cs = d["candidates"]
    opt1 = next(v for k, v in cs.items() if k.startswith("Option 1"))
    opt2 = next(v for k, v in cs.items() if k.startswith("Option 2"))
    assert opt2["B2"]["pass_rate"] > opt1["B2"]["pass_rate"]
    assert d["option1_fate_cost"]["fold_oversampled"] > 5.0     # the cost Option 2 avoids


# ===================================================================== #
# Stage 1.5.3 STEP 5c -- C-5's THRESHOLD, registered before the code.   #
# 5b pinned a fixed W=8; the readiness audit found that would have      #
# handicapped step 6's CONTROL arm and tilted the result toward the     #
# treatment's own conclusion. These grade the rule that replaced it.    #
# ===================================================================== #
_C5C = _Path(__file__).resolve().parents[1] / "results" / "register_c5c_bar_results.json"


def test_s153_c5c_control_arm_never_accumulates():
    """Bar A1, and the most important row in this file for step 6.

    The unmasked arm has every cell age-valid, so the first batch already clears `k` and the
    window must close at W=1 every single time -- i.e. training is bit-identical to today. If
    this is ever below 1.0, `scorecard/baseline.json` has stopped being a valid reference and
    the step-6 comparison is confounded by the mechanism meant to de-confound it.

    An EQUALITY, deliberately, not a >= 0.95 rate: "almost never accumulates" is not identity.
    """
    a1 = _load_json(_C5C)["bars"]["A1"]
    assert a1["value"] == 1.0 and a1["pass"] is True
    assert _load_json(_C5C)["arm_a"]["mean_batches_per_window"] == 1.0


def test_s153_c5c_masked_arm_clears_the_cell_count_bar():
    """Bar A2 -- B2's >= 4 cells, restated for the adaptive rule, where it now holds by
    construction rather than by luck. The residual shortfall is only the W_max forced close."""
    a2 = _load_json(_C5C)["bars"]["A2"]
    assert a2["value"] >= a2["bar"] and a2["pass"] is True


def test_s153_c5c_beats_the_fixed_w_design_it_replaced():
    """Bar A3. The 5c redesign exists to remove a bias, but it must not cost arm B its age
    optimisation in the process -- otherwise it trades one defect for another."""
    a3 = _load_json(_C5C)["bars"]["A3"]
    assert a3["value"] > a3["bar"] and a3["pass"] is True


def test_s153_c5c_attempt_1_is_kept_as_a_regression():
    """Forcing a close at each epoch's end manufactured one partial window per epoch and failed
    A2 on its own (4.44 pp of a 6.12 pp shortfall). Kept measurable so the record shows why it
    was dropped. If this ever passes, the geometry changed and the choice must be revisited."""
    d = _load_json(_C5C)
    assert d["attempt1_close_at_epoch_end"]["frac_windows_ge_k"] < d["bars"]["A2"]["bar"]


def test_s153_c5c_all_bars_pass():
    assert _load_json(_C5C)["all_pass"] is True


_GC2 = _Path(__file__).resolve().parents[1] / "results" / "register_gc_step2_bar_results.json"


def test_s153_step6_crossover_is_solved_not_read_off_the_grid():
    """The sweep grid jumps 1.0 -> 2.0, so its largest PASSING gridpoint is 1.0 while the true
    crossover is ~1.91. Reading the gridpoint understated the usable SD by about 2x and would
    have declared any run with SD in (1.0, 1.91] INCONCLUSIVE while it was in fact >=95% powered
    -- discarding a real result on a reporting artefact. Corrected 2026-08-02.

    Pinned here so the grid can be re-tuned for display without silently re-introducing it.
    """
    d = _load_json(_GC2)
    solved, gridpoint = d["max_resolvable_sd_years"], d["max_resolvable_sd_gridpoint"]
    assert solved > gridpoint, "the crossover was read off the grid again"
    assert 1.7 < solved < 2.1, f"crossover moved to {solved}; the power model changed"


def test_s153_step6_crossover_actually_delivers_the_power_it_claims():
    """Independent of the script: simulate at the reported crossover and just above it."""
    import numpy as _np
    from scipy.stats import t as _t
    sd = _load_json(_GC2)["max_resolvable_sd_years"]
    delta, n, sims = _load_json(_GC2)["delta_star_years"], 6, 60_000
    tcrit = float(_t.ppf(0.975, n - 1))

    def _power(s):
        rng = _np.random.default_rng(4)
        d = rng.normal(delta, s, size=(sims, n))
        m, sdv = d.mean(axis=1), d.std(axis=1, ddof=1)
        return float(((_np.abs(m) > tcrit * sdv / _np.sqrt(n)) & (m > 0)).mean())

    assert _power(sd) >= MIN_PASS_RATE - 0.01     # powered at the crossover
    assert _power(sd * 1.6) < MIN_PASS_RATE       # and genuinely not, well above it
