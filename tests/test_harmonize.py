"""STAGE 1.5 — the harmonization + ΔAge zero-point audit as unit tests.

Four plan documents assert cross-modality harmonization is "unit-tested" and its
"intercept cancellation proven" (`MASTER_PLAN.md:48`, `STAGE_5_PUBLICATION.md:127`,
`STAGE_6_NEW_DATA.md:143`) — yet **no test imported `harmonize.py`**. This file makes
that claim true rather than weaker: it proves the part that cancels, states the exact
scale-gain that does *not* cancel (correcting the docstring's "batch-immune by
construction"), and makes the silent ΔAge zero-point fallback in `aging.py:88`
**visible in a test** so any future change to it is deliberate.

Every invariant below is derived from the real code, not assumed:

    transform(x)       = (align(x) - mu_d) / (sigma_d + EPS)         # harmonize.py:118
    project_to_clock(z)= z * sigma_ref + mu_ref                       # harmonize.py:127
    age                = x_proj @ w + intercept                       # aging.py:47
    ΔAge               = age - per-line-control-mean(age)             # aging.py:144

Composing them, for a perturbed cell whose raw profile is a control cell plus δ, in
the SAME line, ΔAge collapses to the closed form pinned in Group B:

    ΔAge = Σ_g  δ_g · sigma_ref,g / (sigma_d,g + EPS) · w_g

mu_d, mu_ref and the clock intercept all cancel (control-relative difference); sigma_d
does NOT — it is a per-dataset multiplicative gain. Groups A–D run on synthetic data
with no repo artefacts. Group E replays the real build and is skipped when its data is
absent; its logic (`decide_verdict`) is unit-tested here regardless.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellfate.data.aging import (
    LinearClock,
    _control_baseline,
    delta_age,
    recenter_on_control_arrays,
    recenter_on_controls,
)
from cellfate.data.harmonize import (
    EPS,
    MIN_REPLICATES,
    Harmonizer,
)

# verify_stage1_5.py lives at the repo root; tests/ has no __init__.py, so pytest does not
# put the repo root on sys.path. Load it by path — exactly as tests/test_verify_1a.py loads
# verify_1a.py — and reuse its PURE decision logic (importing it runs no data code, no I/O).
_SPEC = importlib.util.spec_from_file_location(
    "verify_stage1_5", Path(__file__).resolve().parents[1] / "verify_stage1_5.py")
verify_stage1_5 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = verify_stage1_5   # register BEFORE exec so @dataclass can resolve
_SPEC.loader.exec_module(verify_stage1_5)
ChunkControlStat = verify_stage1_5.ChunkControlStat
decide_verdict = verify_stage1_5.decide_verdict

# --------------------------------------------------------------------------- #
# Synthetic 2-dataset fixture builder                                         #
# --------------------------------------------------------------------------- #
# `d_ds` is the non-reference dataset (plays HFF single-cell); `ref_ds` is the
# reference (plays gill_bulk). Control expression is centred well above the
# admissibility floor so the common gene space is the full gene list, which keeps
# the closed-form indices aligned.
D_DS = "hff_sc"
REF_DS = "gill_bulk"
LINE = "DONOR_A"


def _controls(rng: np.random.Generator, n: int, g: int, *, loc: float, scale: float) -> np.ndarray:
    """A block of `n` control cells over `g` genes, mean ~`loc` (>> expr_floor)."""
    return rng.normal(loc, scale, size=(n, g))


def _fit_fixture(seed: int = 0, g: int = 6, *, d_scale: float = 0.7, ref_scale: float = 1.3):
    """Return (harmonizer, genes, clock). Two datasets with DIFFERENT control
    spreads so `sigma_d != sigma_ref` and the projection gain is non-trivial."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(g)]
    ctrl_d = _controls(rng, 60, g, loc=2.0, scale=d_scale)
    ctrl_ref = _controls(rng, 60, g, loc=2.0, scale=ref_scale)
    h = Harmonizer.fit(
        {REF_DS: [(ctrl_ref, genes)], D_DS: [(ctrl_d, genes)]},
        ref_dataset=REF_DS,
    )
    # A deterministic clock over exactly the harmonizer's gene space.
    w = rng.normal(0.0, 1.0, size=len(h.genes))
    clock = LinearClock({gene: float(wi) for gene, wi in zip(h.genes, w, strict=True)},
                        intercept=40.0)
    return h, genes, clock


def _dage(h: Harmonizer, clock: LinearClock, X: np.ndarray, ds: str,
          is_ctrl: list[bool], *, source: str = "synth") -> np.ndarray:
    """Run the exact production ΔAge path for a chunk of cells `X` (on h.genes order)
    belonging to one line: transform -> project_to_clock -> delta_age."""
    z = h.transform(X, h.genes, ds)
    x_clock = h.project_to_clock(z)
    obs = pd.DataFrame({"cell_line": [LINE] * len(X), "is_control": is_ctrl})
    d, _mask = delta_age(clock, x_clock, h.genes, obs, source)
    return d


# ===================================================================== #
# GROUP A — the proof the plans already promise                         #
# ===================================================================== #
def test_intercept_cancels():
    """Perturbing the clock intercept leaves ΔAge unchanged: ΔAge is a control-relative
    difference, so the additive intercept cancels.

    It cancels NUMERICALLY, not symbolically — `age + b`, then subtracting a control mean,
    re-rounds, so a large `b` shifts the low bits. It is immune to ~1e-12, not bit-exact;
    `np.array_equal` here would be exactly the kind of overclaim STAGE 1.5 exists to correct
    ('unit-tested', 'batch-immune by construction'). This asserts what is true."""
    h, _genes, clock = _fit_fixture()
    rng = np.random.default_rng(1)
    X = rng.normal(2.0, 0.6, size=(8, len(h.genes)))
    is_ctrl = [True, True, True, False, False, False, False, False]

    clock_hi = LinearClock(dict(clock.weights), intercept=clock.intercept + 137.0)
    d0 = _dage(h, clock, X, D_DS, is_ctrl)
    d1 = _dage(h, clock_hi, X, D_DS, is_ctrl)
    assert np.allclose(d0, d1, rtol=0, atol=1e-9)


def test_additive_batch_offset_is_immune():
    """Adding a per-gene offset to ALL of one dataset's cells and refitting leaves
    ΔAge unchanged: mu_d absorbs the offset, so the Z-score is identical."""
    h, genes, clock = _fit_fixture()
    rng = np.random.default_rng(2)
    ctrl_d = _controls(rng, 60, len(genes), loc=2.0, scale=0.7)
    ctrl_ref = _controls(rng, 60, len(genes), loc=2.0, scale=1.3)
    h0 = Harmonizer.fit({REF_DS: [(ctrl_ref, genes)], D_DS: [(ctrl_d, genes)]}, ref_dataset=REF_DS)
    # rebuild the clock on h0's gene space (same genes, so identical mapping)
    clock0 = LinearClock({g: clock.weights[g] for g in h0.genes}, intercept=40.0)

    X = rng.normal(2.0, 0.6, size=(6, len(genes)))
    is_ctrl = [True, True, False, False, False, False]
    d_before = _dage(h0, clock0, X, D_DS, is_ctrl)

    offset = rng.normal(0.0, 1.0, size=len(genes))          # per-gene additive batch effect
    h1 = Harmonizer.fit(
        {REF_DS: [(ctrl_ref, genes)], D_DS: [(ctrl_d + offset, genes)]}, ref_dataset=REF_DS)
    clock1 = LinearClock({g: clock.weights[g] for g in h1.genes}, intercept=40.0)
    d_after = _dage(h1, clock1, X + offset, D_DS, is_ctrl)

    assert np.allclose(d_before, d_after, atol=1e-9)


def test_reference_mean_drops_out():
    """mu_ref enters project_to_clock as a constant added to every cell, so it
    cancels in the control-relative ΔAge. Shifting it changes nothing."""
    h, _genes, clock = _fit_fixture()
    rng = np.random.default_rng(3)
    X = rng.normal(2.0, 0.6, size=(6, len(h.genes)))
    is_ctrl = [True, True, False, False, False, False]

    d_before = _dage(h, clock, X, D_DS, is_ctrl)
    h._stats[REF_DS].mu[:] = h._stats[REF_DS].mu + rng.normal(0.0, 5.0, size=len(h.genes))
    d_after = _dage(h, clock, X, D_DS, is_ctrl)
    assert np.allclose(d_before, d_after, rtol=0, atol=1e-9)


# ===================================================================== #
# GROUP B — the TRUE scope: sigma_d is a gain, not immune                #
# ===================================================================== #
def test_scale_is_a_gain_matches_the_closed_form():
    """The exact invariant the manuscript should state instead of "batch-immune by
    construction": for a perturbed cell = control + δ in one line,

        ΔAge = Σ_g δ_g · sigma_ref,g / (sigma_d,g + EPS) · w_g
    """
    h, _genes, clock = _fit_fixture()
    rng = np.random.default_rng(4)
    x_ctrl = rng.normal(2.0, 0.6, size=len(h.genes))
    delta = rng.normal(0.0, 0.5, size=len(h.genes))
    x_pert = x_ctrl + delta
    X = np.vstack([x_ctrl, x_pert])

    d = _dage(h, clock, X, D_DS, [True, False])

    sigma_d = h._stats[D_DS].sigma
    sigma_ref = h._stats[REF_DS].sigma
    w = np.array([clock.weights[g] for g in h.genes])
    expected = float(np.sum(delta * sigma_ref / (sigma_d + EPS) * w))

    assert d[0] == pytest.approx(0.0, abs=1e-9)           # the control sits at its own baseline
    assert d[1] == pytest.approx(expected, rel=1e-9, abs=1e-9)


def test_the_gain_actually_differs_between_datasets_so_it_is_not_immune():
    """The same raw δ produces a DIFFERENT ΔAge depending on the dataset's sigma —
    which is exactly what "carries a scale factor by design" means. A truly
    batch-immune transform would give the same ΔAge; this asserts it does not."""
    h, _genes, clock = _fit_fixture(d_scale=0.7, ref_scale=1.3)
    rng = np.random.default_rng(5)
    x_ctrl = rng.normal(2.0, 0.6, size=len(h.genes))
    delta = rng.normal(0.0, 0.5, size=len(h.genes))
    X = np.vstack([x_ctrl, x_ctrl + delta])

    d_in_d = _dage(h, clock, X, D_DS, [True, False])[1]     # gain sigma_ref/sigma_d
    d_in_ref = _dage(h, clock, X, REF_DS, [True, False])[1]  # gain sigma_ref/sigma_ref ~ 1

    # ratio of ΔAge equals the ratio of gains — the proof it is a scale factor, not immune
    sigma_d = h._stats[D_DS].sigma
    sigma_ref = h._stats[REF_DS].sigma
    w = np.array([clock.weights[g] for g in h.genes])
    num = float(np.sum(delta * sigma_ref / (sigma_d + EPS) * w))
    den = float(np.sum(delta * sigma_ref / (sigma_ref + EPS) * w))
    assert not np.isclose(d_in_d, d_in_ref)                 # NOT immune
    assert d_in_d / d_in_ref == pytest.approx(num / den, rel=1e-9)


# ===================================================================== #
# GROUP C — fit / leak-safety (Harmonizer.fit)                          #
# ===================================================================== #
def test_fit_uses_only_the_controls_it_is_handed():
    """The held-out-donor guarantee: mu/sigma are exactly the moments of the passed
    control cells (variance-floored). A cell not in `controls` cannot move them."""
    rng = np.random.default_rng(6)
    genes = [f"G{i}" for i in range(5)]
    ctrl_d = _controls(rng, 40, 5, loc=2.0, scale=0.8)
    ctrl_ref = _controls(rng, 40, 5, loc=2.0, scale=1.1)
    h = Harmonizer.fit({REF_DS: [(ctrl_ref, genes)], D_DS: [(ctrl_d, genes)]}, ref_dataset=REF_DS)

    raw_sigma = ctrl_d.std(axis=0)
    expected_sigma = np.maximum(raw_sigma, np.median(raw_sigma))
    assert np.allclose(h._stats[D_DS].mu, ctrl_d.mean(axis=0))
    assert np.allclose(h._stats[D_DS].sigma, expected_sigma)


def test_variance_floor_lifts_every_sigma_to_at_least_the_median():
    """sigma is floored at median(sigma), so a near-constant gene cannot produce an
    exploding Z-score."""
    rng = np.random.default_rng(7)
    genes = [f"G{i}" for i in range(6)]
    ctrl = _controls(rng, 50, 6, loc=2.0, scale=1.0)
    ctrl[:, 0] = 2.0                                        # gene 0: zero variance
    ref = _controls(rng, 50, 6, loc=2.0, scale=1.0)
    h = Harmonizer.fit({REF_DS: [(ref, genes)], D_DS: [(ctrl, genes)]}, ref_dataset=REF_DS)

    sigma = h._stats[D_DS].sigma
    assert np.all(sigma >= np.median(ctrl.std(axis=0)) - 1e-12)
    assert sigma[0] == pytest.approx(np.median(ctrl.std(axis=0)))   # floored, not ~0


def test_common_gene_space_is_the_sorted_intersection_of_admissible_sets():
    """G = genes admissible (mean control expr >= floor) in EVERY dataset, sorted."""
    genes = ["B", "A", "C", "D"]                # deliberately unsorted input order
    # columns are [B, A, C, D]. dataset d: B,A,C admissible, D below floor;
    #                            ref:       B,C,D admissible, A below floor.
    d = np.array([[3.0, 3.0, 3.0, 0.0]] * 5)    # B=3 A=3 C=3 D=0  -> {B,A,C}
    r = np.array([[3.0, 0.0, 3.0, 3.0]] * 5)    # B=3 A=0 C=3 D=3  -> {B,C,D}
    h = Harmonizer.fit({REF_DS: [(r, genes)], D_DS: [(d, genes)]}, ref_dataset=REF_DS)
    assert h.genes == ["B", "C"]                # sorted( {B,A,C} ∩ {B,C,D} ) == ["B","C"]


def test_fit_raises_below_min_replicates():
    genes = [f"G{i}" for i in range(4)]
    too_few = np.full((MIN_REPLICATES - 1, 4), 2.0)
    ref = np.full((10, 4), 2.0) + np.random.default_rng(8).normal(0, 0.1, size=(10, 4))
    with pytest.raises(ValueError, match="control observations"):
        Harmonizer.fit({REF_DS: [(ref, genes)], D_DS: [(too_few, genes)]}, ref_dataset=REF_DS)


def test_transform_raises_on_unknown_dataset():
    h, genes, _clock = _fit_fixture()
    with pytest.raises(KeyError, match="no harmonization stats"):
        h.transform(np.zeros((2, len(h.genes))), h.genes, "not_a_dataset")


def test_fit_raises_when_reference_not_present():
    genes = [f"G{i}" for i in range(4)]
    block = np.full((10, 4), 2.0) + np.random.default_rng(9).normal(0, 0.1, size=(10, 4))
    with pytest.raises(ValueError, match="reference dataset"):
        Harmonizer.fit({D_DS: [(block, genes)]}, ref_dataset="gill_bulk")


def test_align_places_permuted_and_missing_genes_in_the_right_columns():
    """transform aligns an input in ANY gene order onto h.genes; a gene absent from
    the input becomes a zero column, never a mis-slotted one."""
    h, _genes, _clock = _fit_fixture(g=4)
    permuted = list(reversed(h.genes))                     # reverse order
    row = np.arange(1.0, len(h.genes) + 1.0)               # distinct per-gene values
    x_perm = row[::-1][None, :]                            # values matching `permuted`
    z_perm = h.transform(x_perm, permuted, D_DS)
    z_ref = h.transform(row[None, :], h.genes, D_DS)
    assert np.allclose(z_perm, z_ref)                      # same cell, order-independent

    drop = h.genes[1]
    kept = [g for g in h.genes if g != drop]
    x_missing = row[[h.genes.index(g) for g in kept]][None, :]
    z_missing = h.transform(x_missing, kept, D_DS)
    # the dropped gene's aligned input is 0 -> its Z-score is (0 - mu)/(sigma+EPS)
    st = h._stats[D_DS]
    j = h.genes.index(drop)
    assert z_missing[0, j] == pytest.approx((0.0 - st.mu[j]) / (st.sigma[j] + EPS))


def test_to_json_from_json_round_trips(tmp_path):
    h, _genes, _clock = _fit_fixture()
    p = tmp_path / "harmonization.json"
    h.to_json(p)
    h2 = Harmonizer.from_json(p)
    assert h2.genes == h.genes and h2.ref_dataset == h.ref_dataset
    x = np.random.default_rng(10).normal(2.0, 0.6, size=(4, len(h.genes)))
    assert np.allclose(h.transform(x, h.genes, D_DS), h2.transform(x, h2.genes, D_DS))


# ===================================================================== #
# GROUP D — the ΔAge zero-point (aging.py), incl. the silent fallback   #
# ===================================================================== #
def test_control_baseline_is_per_line_both_lines_land_at_zero():
    """Two lines with different absolute ages both centre on ~0 after subtracting
    their OWN control baseline — the per-line control-relative zero-point."""
    ages = np.array([10.0, 12.0, 14.0, 100.0, 104.0, 108.0])   # line A ~12, line B ~104
    lines = np.array(["A", "A", "A", "B", "B", "B"])
    is_ctrl = np.array([True, True, False, True, True, False])
    centred = recenter_on_control_arrays(ages, lines, is_ctrl)
    assert centred[:2].mean() == pytest.approx(0.0)            # A controls
    assert centred[3:5].mean() == pytest.approx(0.0)          # B controls
    # the perturbed cell in each line is measured against ITS line's controls
    assert centred[2] == pytest.approx(14.0 - 11.0)
    assert centred[5] == pytest.approx(108.0 - 102.0)


def test_the_silent_no_control_fallback_self_centres_a_line_to_zero():
    """THE audit's core check, made visible: a line with NO controls in its chunk
    falls back to self-centring (aging.py:88), which forces that line's MEAN ΔAge to
    0 by construction — dragging a real per-donor offset toward zero with no warning.
    Pinning it makes any future change to this behaviour a deliberate, reviewed act."""
    clock = LinearClock({"G0": 1.0}, intercept=40.0)
    expr = np.array([[5.0], [7.0], [9.0]])                    # ages 45, 47, 49
    obs = pd.DataFrame({"cell_line": ["X", "X", "X"], "is_control": [False, False, False]})
    d, _mask = delta_age(clock, expr, ["G0"], obs, source="synth")
    assert d.mean() == pytest.approx(0.0)                     # self-centred, not control-relative
    # and it is genuinely the within-line mean that was subtracted, not a control mean
    assert np.allclose(d, np.array([45.0, 47.0, 49.0]) - 47.0)


def test_control_baseline_matches_the_raw_control_mean_when_controls_exist():
    """When controls ARE present the baseline is their mean, not the line mean —
    the two branches of aging.py:88 must be distinguishable, which is the whole point."""
    values = np.array([1.0, 3.0, 100.0])                     # 2 controls, 1 far perturbed
    lines = np.array(["L", "L", "L"])
    is_ctrl = np.array([True, True, False])
    base = _control_baseline(values, lines, is_ctrl)
    assert np.allclose(base, 2.0)                            # mean(1,3), NOT mean(1,3,100)


def test_recenter_on_controls_dataframe_form_matches_the_array_form():
    """`recenter_on_controls` is the obs-DataFrame wrapper called after cell-cycle
    deconfounding re-centres the population; it must reproduce the array form exactly,
    restoring the control-relative zero-point the rejuvenation score depends on."""
    values = np.array([10.0, 12.0, 14.0, 100.0, 104.0])
    obs = pd.DataFrame({
        "cell_line": ["A", "A", "A", "B", "B"],
        "is_control": [True, True, False, True, False],
    })
    got = recenter_on_controls(values, obs)
    want = recenter_on_control_arrays(
        values, obs["cell_line"].to_numpy(), obs["is_control"].to_numpy().astype(bool))
    assert np.array_equal(got, want)


# ===================================================================== #
# GROUP E — real-build fallback scan: LOGIC unit-tested here (decide_verdict was  #
# loaded at the top); the real replay lives in verify_stage1_5.py and runs on the #
# data machine.                                                                   #
# ===================================================================== #
def test_decide_verdict_passes_when_every_chunk_has_a_control():
    stats = [
        ChunkControlStat("gse:HFF:d0", "HFF", n_cells=100, n_control=20),
        ChunkControlStat("gill:N2", "N2", n_cells=10, n_control=2),
    ]
    v = decide_verdict(stats)
    assert v["status"] == "PASS" and v["fallback_chunks"] == []


def test_decide_verdict_fails_and_names_a_control_less_chunk():
    """A chunk with perturbed cells but zero controls fired the silent fallback."""
    stats = [
        ChunkControlStat("gill:N2", "N2", n_cells=10, n_control=2),
        ChunkControlStat("gill:O1", "O1", n_cells=8, n_control=0),   # fallback fired
    ]
    v = decide_verdict(stats)
    assert v["status"] == "FAIL"
    assert [c["chunk_id"] for c in v["fallback_chunks"]] == ["gill:O1"]


def test_decide_verdict_ignores_all_control_chunks():
    """A chunk that is entirely controls (no perturbed cells) cannot fire the
    fallback — there is nothing whose zero-point could flip."""
    stats = [ChunkControlStat("gill:N2:ctrl_only", "N2", n_cells=5, n_control=5)]
    assert decide_verdict(stats)["status"] == "PASS"


def test_decide_verdict_reports_cannot_verify_on_no_chunks():
    assert decide_verdict([])["status"] == "CANNOT_VERIFY"
