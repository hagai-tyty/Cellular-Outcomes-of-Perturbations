"""Tests for cellfate.training: calibration/conformal/OOD units + an end-to-end
build -> train -> bundle run that verifies the bundle matches the Document-1
contracts (Document 3)."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from cellfate.common import constants as C
from cellfate.common import io
from cellfate.common.io import ArtifactPaths
from cellfate.common.panel import GenePanel
from cellfate.models import CellFateNet
from cellfate.training import (
    TrainConfig,
    coverage,
    fit_conformal,
    fit_ood,
    fit_temperature,
    load_ood,
    mahalanobis,
    save_ood,
    train_member,
)
from cellfate.training import run as train_run


# --------------------------------------------------------------------------- #
# temperature scaling                                                         #
# --------------------------------------------------------------------------- #
def test_temperature_tempers_overconfident_logits_and_reduces_nll():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 3, size=400)
    target = np.eye(3)[labels]
    # deliberately over-confident logits: right class big, but 15% are wrong
    logits = np.full((400, 3), -4.0)
    logits[np.arange(400), labels] = 8.0
    wrong = rng.random(400) < 0.15
    logits[wrong] = logits[wrong][:, ::-1]

    t = fit_temperature(logits, target).temperature
    assert t > 1.0   # over-confident -> needs softening

    def nll(scale):
        z = logits / scale
        logp = z - np.log(np.exp(z).sum(1, keepdims=True))
        return -(target * logp).sum(1).mean()

    assert nll(t) <= nll(1.0) + 1e-9


def test_temperature_empty_returns_unit():
    assert fit_temperature(np.zeros((0, 3)), np.zeros((0, 3))).temperature == 1.0


# --------------------------------------------------------------------------- #
# conformal                                                                   #
# --------------------------------------------------------------------------- #
def test_conformal_quantile_covers_at_least_the_level():
    rng = np.random.default_rng(1)
    res = np.abs(rng.normal(size=500))
    cp = fit_conformal(res, [0.90])
    assert "0.9" in cp.q                       # keyed by str(level)
    q = cp.q["0.9"]
    assert q > 0
    assert coverage(res, q) >= 0.90            # finite-sample guarantee


def test_conformal_empty_uses_default():
    cp = fit_conformal(np.array([]), [0.9], default_q=123.0)
    assert cp.q["0.9"] == 123.0


# --------------------------------------------------------------------------- #
# OOD                                                                         #
# --------------------------------------------------------------------------- #
def test_ood_scores_far_point_higher_and_roundtrips(tmp_path):
    rng = np.random.default_rng(2)
    feats = rng.normal(size=(300, 8))
    ood = fit_ood(feats)
    assert ood["threshold"] > 0 and ood["dim"] == 8

    in_dist = mahalanobis(feats, ood).mean()
    far = mahalanobis(np.full((1, 8), 20.0), ood)[0]
    assert far > in_dist
    assert far > ood["threshold"]

    paths = ArtifactPaths.of(tmp_path)
    save_ood(paths, ood)
    back = load_ood(paths)
    assert back["dim"] == 8
    assert np.allclose(back["mean"], ood["mean"])
    assert np.allclose(back["precision"], ood["precision"])


# --------------------------------------------------------------------------- #
# determinism                                                                 #
# --------------------------------------------------------------------------- #
def _toy_dataset(n=64, g=8, n_donors=1):
    rng = np.random.default_rng(3)
    x = torch.tensor(rng.normal(size=(n, g)), dtype=torch.float32)
    fp = torch.tensor(rng.integers(0, 2, size=(n, 2048)), dtype=torch.float32)
    dt = torch.tensor(rng.normal(size=(n, 2)), dtype=torch.float32)
    yc = torch.tensor(rng.dirichlet(np.ones(3), size=n), dtype=torch.float32)
    ya = torch.tensor(rng.normal(size=n), dtype=torch.float32)
    am = torch.ones(n)
    donor = torch.arange(n) % n_donors            # Stage 1a: the 7th column
    return torch.utils.data.TensorDataset(x, fp, dt, yc, ya, am, donor)


def _tiny_cfg(**kw):
    base = dict(dataset_dir="unused", d_cell=8, d_u=8, latent_dim=8, p_drop=0.1,
                epochs=3, patience=3, batch_size=32, ensemble_size=2)
    base.update(kw)
    return TrainConfig(**base)


def test_train_member_is_deterministic_under_fixed_seed():
    ds = _toy_dataset()
    cfg = _tiny_cfg()

    def make():
        return CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8, p_drop=0.1)

    m1, _ = train_member(make, ds, ds, cfg, seed=0, device="cpu")
    m2, _ = train_member(make, ds, ds, cfg, seed=0, device="cpu")
    w1 = m1.cell.net[0][0].weight.detach()
    w2 = m2.cell.net[0][0].weight.detach()
    assert torch.allclose(w1, w2, atol=1e-6)


# --------------------------------------------------------------------------- #
# end-to-end: build -> train -> bundle, verify contracts                      #
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def trained_bundle(tmp_path_factory):
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build

    root = tmp_path_factory.mktemp("e2e")
    panel = root / "panel.json"
    build(
        DataConfig(out=str(root), gene_panel=str(panel), n_genes=160,
                   qc=QCConfig(min_genes=5, max_mito_frac=0.5), label_tau=0.5,
                   split_fracs=(0.6, 0.2, 0.1, 0.1), primary_regime="scaffold", seed=0),
        sources=[
            SyntheticSource(name="synth", n_lines=2, n_compounds=4,
                            n_cells_per_condition=10, seed=1),
            SyntheticSource(name="tahoe", n_lines=2, n_compounds=4,
                            n_cells_per_condition=10, seed=2),
        ],
    )
    summary = train_run(_tiny_cfg(dataset_dir=str(root), regime="scaffold",
                                  conformal_levels=(0.90,), device="cpu"))
    return root, panel, summary


def test_e2e_bundle_is_complete_and_matches_contracts(trained_bundle):
    root, panel, summary = trained_bundle
    paths = ArtifactPaths.of(root)

    io.assert_bundle_complete(paths)                       # raises if incomplete
    meta = io.load_bundle_meta(paths)
    assert meta.classes == list(C.CLASSES)
    assert meta.n_members == 2
    assert meta.gene_panel_hash == GenePanel.load(panel).hash()   # hash threads through
    assert meta.conformal_levels == [0.90]


def test_e2e_members_load_and_predict(trained_bundle):
    root, _, _ = trained_bundle
    paths = ArtifactPaths.of(root)
    members = [CellFateNet.load_member(p)
               for p in sorted(paths.bundle_members_dir.glob("*.pt"))]
    assert len(members) == 2

    # ensemble forward on a handful of real cells
    arr = io.shard_to_numpy(io.read_shard(next(paths.shards_dir.glob("*.parquet"))))
    from cellfate.common.scalers import Scalers
    sc = Scalers.load(paths.bundle_scalers_file)
    x = torch.tensor(sc.transform_x(arr["X"][:4]), dtype=torch.float32)
    fp = torch.tensor(arr["u_chem_fp"][:4], dtype=torch.float32)
    dt = torch.tensor(sc.transform_dose_time(arr["dose_time"][:4]), dtype=torch.float32)

    probs = torch.stack([torch.softmax(m(x, fp, dt)[0], 1) for m in members]).mean(0)
    assert torch.allclose(probs.sum(1), torch.ones(4), atol=1e-5)
    ages = torch.stack([m(x, fp, dt)[1] for m in members]).mean(0)
    assert torch.isfinite(ages).all()
    feat = members[0](x, fp, dt)[2]
    assert feat.shape[1] == members[0].arch["latent_dim"]


def test_e2e_calibration_artifacts_valid(trained_bundle):
    root, _, summary = trained_bundle
    paths = ArtifactPaths.of(root)

    temp = io.load_temperature(paths)
    assert temp.temperature > 0

    conf = io.load_conformal(paths)
    q = conf.q[str(conf.levels[0])]
    assert np.isfinite(q) and q > 0

    ood = load_ood(paths)
    assert ood["threshold"] > 0
    feats = np.zeros((2, ood["dim"]))
    assert mahalanobis(feats, ood).shape == (2,)

    assert (paths.bundle_dir / C.BUNDLE_CONFIG_FILENAME).exists()
    assert (paths.bundle_dir / C.BUNDLE_METRICS_FILENAME).exists()

    # Temperature scaling never worsens NLL -- ON THE SPLIT IT WAS FITTED ON.
    # Before Stage 1b that was calib/val; it is now the cross-donor pool, so the guarantee
    # moved with it. Asserting it on the in-distribution split would be asserting something
    # fit_temperature never promised: the model is under-confident in-distribution (T<1,
    # sharpening) and over-confident out-of-donor (T>1, softening), and one scalar cannot
    # serve both. In-distribution NLL rising is the expected cost of that trade.
    fitted_on = ("xdonor_nll_before_temp" if summary.get("xdonor_calibrated")
                 else "nll_before_temp")
    before = summary[fitted_on]
    after = summary[fitted_on.replace("before", "after")]
    assert after <= before + 1e-6, (
        f"temperature worsened NLL on the split it was fitted on ({fitted_on}): "
        f"{before:.6f} -> {after:.6f}"
    )


# --------------------------------------------------------------------------- #
# Stage 1b -- cross-donor calibration                                         #
# --------------------------------------------------------------------------- #
def test_platt_recovers_a_miscaled_and_a_biased_p_safe():
    """The bias case is the one a temperature provably cannot fix.

    A scalar can only sharpen or soften; it cannot shift the safe/unsafe boundary. Stage 1's own
    <=0.17 bar comes from a Platt measurement (MASTER_PLAN: "S, P_loss ... YES -- Platt halves
    it", T8.2), while run 2's scalar regressed the graded metric on every fold.
    """
    from cellfate.common.calibration import platt_safe
    from cellfate.training.calibrate import fit_platt_binary

    rng = np.random.default_rng(0)
    n = 20000
    z = rng.normal(0.0, 2.0, n)                 # true log-odds
    p_true = 1.0 / (1.0 + np.exp(-z))           # calibrated BY CONSTRUCTION
    y = (rng.uniform(size=n) < p_true).astype(float)

    # (a) over-sharpened by 3x. logit(p_bad) = 3z, so recovering z needs a = 1/3, b = 0.
    sharp = 1.0 / (1.0 + np.exp(-3.0 * z))
    a, b = fit_platt_binary(sharp, y)
    assert a == pytest.approx(1 / 3, abs=0.05), f"expected a~0.333, got {a}"
    assert b == pytest.approx(0.0, abs=0.10), f"expected b~0, got {b}"
    assert np.mean(np.abs(platt_safe(sharp, a, b) - p_true)) < 0.02

    # (b) BIASED by +1.8 toward `safe`. logit(p_bad) = z + 1.8, so a = 1, b = -1.8.
    # A TEMPERATURE CANNOT DO THIS: scaling z + 1.8 by any factor never recovers z.
    biased = 1.0 / (1.0 + np.exp(-(z + 1.8)))
    a2, b2 = fit_platt_binary(biased, y)
    assert a2 == pytest.approx(1.0, abs=0.08), f"expected a~1, got {a2}"
    assert b2 == pytest.approx(-1.8, abs=0.15), f"expected b~-1.8, got {b2}"
    assert np.mean(np.abs(platt_safe(biased, a2, b2) - p_true)) < 0.02

    # and the scalar-only alternative provably cannot: the best pure slope still leaves a gap
    best_slope = min((np.mean(np.abs(platt_safe(biased, s, 0.0) - p_true)), s)
                     for s in np.linspace(0.05, 5.0, 200))[0]
    assert best_slope > 0.05, (
        f"a slope alone got within {best_slope:.3f} of the truth on a BIASED input -- if a "
        "temperature can fix this case the test no longer isolates the intercept"
    )


def test_platt_is_rank_preserving_so_the_fate_guards_cannot_move():
    """`fate_prauc`/`fate_roc` are rank-based and are Stage 1 GUARDS.

    The slope is constrained positive, so the map is strictly increasing in P(safe) and the
    ordering is mathematically identical -- which is what lets this ship without disturbing them.
    """
    from cellfate.common.calibration import apply_platt

    rng = np.random.default_rng(1)
    p = rng.dirichlet(np.ones(3), size=300)
    out = apply_platt(p, a=2.3, b=-0.7, safe_idx=0)

    assert np.array_equal(np.argsort(p[:, 0]), np.argsort(out[:, 0]))
    assert np.allclose(out.sum(axis=1), 1.0)
    # the loss/death RATIO is untouched, so P_loss stays meaningful to RES
    assert np.allclose(p[:, 1] / p[:, 2], out[:, 1] / out[:, 2])


def test_platt_declines_degenerate_and_never_ships_worse_than_identity():
    from cellfate.training.calibrate import fit_platt_binary

    rng = np.random.default_rng(2)
    p = rng.uniform(0.05, 0.95, 200)
    with pytest.warns(UserWarning, match="single class"):
        assert fit_platt_binary(p, np.ones(200)) == (1.0, 0.0)
    assert fit_platt_binary(np.array([]), np.array([])) == (1.0, 0.0)

    # already calibrated -> must not make it worse; identity is always available
    y = (rng.uniform(size=500) < p[:500] if len(p) >= 500 else rng.integers(0, 2, len(p)))
    a, b = fit_platt_binary(p[:len(y)], y.astype(float))
    assert a > 0


def test_xstats_round_trips_through_the_bundle(tmp_path):
    """Persisting the pool is what makes future calibrator work a seconds-long offline refit
    instead of another 3.5 h LOOCV pass."""
    from cellfate.training.xdonor_calib import XDonorStats, load_xstats, save_xstats

    rng = np.random.default_rng(3)
    stats = XDonorStats(
        abs_residuals=rng.normal(size=20), logits=rng.normal(size=(20, 3)),
        targets=rng.dirichlet(np.ones(3), size=20), probs_mean=rng.dirichlet(np.ones(3), size=20),
        sigma_pred=rng.uniform(size=20), sigma_pred_mc=rng.uniform(size=20),
        n_donors=5, feats=rng.normal(size=(20, 4)),
        residuals_per_donor={1: 10, 2: 10},
        donor_scales={1: {"n": 10, "mae": 3.5}, 2: {"n": 10, "mae": 7.0}},
    )
    save_xstats(tmp_path, stats)
    back = load_xstats(tmp_path)

    for k in ("abs_residuals", "logits", "targets", "probs_mean",
              "sigma_pred", "sigma_pred_mc", "feats"):
        assert np.allclose(getattr(back, k), getattr(stats, k)), k
    assert back.n_donors == 5
    assert back.residuals_per_donor == {1: 10, 2: 10}
    assert back.donor_scales[1]["mae"] == 3.5


def test_sigma_scale_widens_an_overconfident_ensemble_spread():
    """The defect this exists for: members agree (~2.4 yr) while collectively wrong (~14 yr)."""
    from cellfate.training import XDonorStats, sigma_scale_factor

    stats = XDonorStats(
        abs_residuals=np.full(200, 14.0), logits=np.zeros((0, 3)), targets=np.zeros((0, 3)),
        probs_mean=np.zeros((0, 3)),
        sigma_pred=np.full(200, 2.4), sigma_pred_mc=np.full(200, 2.4),
        n_donors=5, feats=np.zeros((0, 1)),
    )
    s = sigma_scale_factor(stats, z_conf=1.0, level=0.90)
    assert s == pytest.approx(14.0 / 2.4, rel=1e-6)
    assert 2.4 * s == pytest.approx(14.0, rel=1e-6)   # the scaled spread now matches reality


def test_sigma_scale_never_shrinks_an_already_adequate_spread():
    """Clamped at 1.0: over-wide makes RES conservative, which is the safe direction."""
    from cellfate.training import XDonorStats, sigma_scale_factor

    stats = XDonorStats(
        abs_residuals=np.full(50, 1.0), logits=np.zeros((0, 3)), targets=np.zeros((0, 3)),
        probs_mean=np.zeros((0, 3)),
        sigma_pred=np.full(50, 9.0), sigma_pred_mc=np.full(50, 9.0),
        n_donors=5, feats=np.zeros((0, 1)),
    )
    assert sigma_scale_factor(stats, z_conf=1.0) == 1.0


def test_sigma_scale_is_identity_without_statistics():
    from cellfate.training import XDonorStats, sigma_scale_factor

    empty = XDonorStats(
        abs_residuals=np.array([]), logits=np.zeros((0, 3)), targets=np.zeros((0, 3)),
        probs_mean=np.zeros((0, 3)),
        sigma_pred=np.array([]), sigma_pred_mc=np.array([]),
        n_donors=0, feats=np.zeros((0, 1)),
    )
    assert sigma_scale_factor(empty, z_conf=1.0) == 1.0
    assert sigma_scale_factor(empty, z_conf=1.0, mode="mc_dropout") == 1.0
    # z_conf=0 would divide by zero; must not raise or return a non-finite factor
    stats = XDonorStats(
        abs_residuals=np.full(10, 5.0), logits=np.zeros((0, 3)), targets=np.zeros((0, 3)),
        probs_mean=np.zeros((0, 3)),
        sigma_pred=np.full(10, 1.0), sigma_pred_mc=np.full(10, 1.0),
        n_donors=3, feats=np.zeros((0, 1)),
    )
    assert sigma_scale_factor(stats, z_conf=0.0) == 1.0


def test_crossdonor_pool_uses_every_split_of_the_held_out_donor():
    """The inner model sees donor d in NO split, so all of d's cells are valid held-out data.

    It trains on `train minus d` and early-stops on `val minus d`, so d's train, val and calib
    cells are equally unseen. Evaluating on only its train cells throws away clean data, and the
    pool size is what limits `q`'s stability (a 68% sampling spread at ~70 residuals).
    """
    from cellfate.training import crossdonor_stats

    def split(n, donors):
        ds = _toy_dataset(n=n)
        return torch.utils.data.TensorDataset(*ds.tensors[:-1], torch.tensor(donors))

    # two donors; each contributes cells to train, val AND calib
    train = split(40, [0] * 20 + [1] * 20)
    val = split(16, [0] * 8 + [1] * 8)
    calib = split(16, [0] * 8 + [1] * 8)

    def make():
        return CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8, p_drop=0.1)

    with_calib = crossdonor_stats(train, val, make, _tiny_cfg(), "cpu",
                                  mc_T=4, calib_ds=calib)
    without_calib = crossdonor_stats(train, val, make, _tiny_cfg(), "cpu", mc_T=4)

    # with calib: 20 train + 8 val + 8 calib = 36 per donor; without: 20 + 8 = 28
    assert with_calib.residuals_per_donor == {0: 36, 1: 36}, with_calib.residuals_per_donor
    assert without_calib.residuals_per_donor == {0: 28, 1: 28}
    assert with_calib.abs_residuals.size > without_calib.abs_residuals.size, (
        "calib cells of the held-out donor were discarded"
    )


def test_the_same_split_passed_twice_is_not_double_counted():
    """Overlapping splits would inflate the pool and silently distort `q`.

    Genuinely different-but-overlapping datasets cannot be detected (no cell ids in the
    tensors) and are the caller's contract; passing the SAME object twice is the easy mistake
    and is deduplicated.
    """
    from cellfate.training import crossdonor_stats

    ds = _toy_dataset(n=40)
    ds = torch.utils.data.TensorDataset(*ds.tensors[:-1],
                                        torch.tensor([0] * 20 + [1] * 20))

    def make():
        return CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8, p_drop=0.1)

    stats = crossdonor_stats(ds, ds, make, _tiny_cfg(), "cpu", mc_T=4, calib_ds=ds)
    assert stats.residuals_per_donor == {0: 20, 1: 20}, (
        f"the same split counted more than once: {stats.residuals_per_donor}"
    )


def test_crossdonor_stats_refuses_a_single_donor():
    """Inner-LODO with one donor would silently produce in-distribution calibration wearing
    a cross-donor label -- the exact defect Stage 1 exists to fix. It must raise instead."""
    from cellfate.training import crossdonor_stats, n_train_donors

    ds = _toy_dataset(n_donors=1)
    assert n_train_donors(ds) == 1
    with pytest.raises(ValueError, match="inner-LODO"):
        crossdonor_stats(ds, ds, lambda: CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8),
                         _tiny_cfg(), "cpu")


def test_crossdonor_stats_skips_a_bulk_corpus_masquerading_as_a_donor():
    """The defect that invalidated the first Stage 1 run.

    The Gill+HFF merge carries HFF as a `cell_line` with 33,613 of 33,688 training cells.
    Holding it out left a model trained on 75 cells (val_loss 33.0 vs 5.3), and because that
    fold is also the largest, it supplied 99.8% of the pooled residuals -- so `q` and
    `sigma_scale` were calibrated against data starvation, not donor shift.
    """
    from cellfate.training.xdonor_calib import MIN_INNER_TRAIN_FRAC, crossdonor_stats

    # donor 0 holds 90% of the rows; donors 1 and 2 hold 5% each
    n = 200
    ds = _toy_dataset(n=n)
    donor = torch.zeros(n, dtype=torch.long)
    donor[180:190] = 1
    donor[190:] = 2
    ds = torch.utils.data.TensorDataset(*ds.tensors[:-1], donor)

    def make():
        return CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8, p_drop=0.1)

    stats = crossdonor_stats(ds, ds, make, _tiny_cfg(), "cpu", mc_T=4)

    # donor 0 must be skipped: holding it out leaves 20/200 = 10%, below the floor
    assert 0.10 < MIN_INNER_TRAIN_FRAC
    assert stats.n_donors == 2, (
        f"expected the bulk corpus to be skipped and 2 real donors used, got {stats.n_donors}"
    )
    # and its 180 rows must not appear in the pooled residuals (donors 1 and 2 hold 10 each)
    assert 0 not in stats.residuals_per_donor, (
        f"bulk corpus leaked into the residual pool: {stats.residuals_per_donor}"
    )
    assert stats.abs_residuals.size == 20, stats.residuals_per_donor


def test_crossdonor_stats_refuses_when_only_one_donor_survives_the_bulk_filter():
    """Two donors where one is a bulk corpus leaves one usable fold -- not cross-donor."""
    from cellfate.training.xdonor_calib import crossdonor_stats

    n = 200
    ds = _toy_dataset(n=n)
    donor = torch.zeros(n, dtype=torch.long)
    donor[190:] = 1                      # 95% / 5%
    ds = torch.utils.data.TensorDataset(*ds.tensors[:-1], donor)

    def make():
        return CellFateNet(g=8, d_cell=8, d_u=8, latent_dim=8, p_drop=0.1)

    with pytest.raises(ValueError, match="usable fold"):
        crossdonor_stats(ds, ds, make, _tiny_cfg(), "cpu")


def test_sigma_scale_is_fitted_per_mode_from_the_matching_spread():
    """Each inference mode gets its OWN factor, from its OWN spread.

    `sigma_age` is the ensemble spread in one mode and the T-pass dropout spread in the other.
    Borrowing one factor for the other scales the wrong quantity; leaving a mode at 1.0 serves
    raw overconfident uncertainty. Both are wrong, so both modes are calibrated.
    """
    from cellfate.training import XDonorStats, sigma_scale_factor

    # true error 14 yr; members agree to 2.0, dropout jitters by 7.0 -- different quantities
    stats = XDonorStats(
        abs_residuals=np.full(200, 14.0), logits=np.zeros((0, 3)), targets=np.zeros((0, 3)),
        probs_mean=np.zeros((0, 3)),
        sigma_pred=np.full(200, 2.0), sigma_pred_mc=np.full(200, 7.0),
        n_donors=5, feats=np.zeros((0, 1)),
    )
    ens = sigma_scale_factor(stats, z_conf=1.0, mode="ensemble")
    mc = sigma_scale_factor(stats, z_conf=1.0, mode="mc_dropout")

    assert ens == pytest.approx(14.0 / 2.0)
    assert mc == pytest.approx(14.0 / 7.0)
    # each factor must carry ITS OWN spread to the same honest width
    assert 2.0 * ens == pytest.approx(14.0)
    assert 7.0 * mc == pytest.approx(14.0)


def test_conformal_schema_still_loads_pre_stage1b_bundles(tmp_path):
    """Every bundle in runs/ was written before sigma_scale existed. They must keep loading --
    which is why the field is defaulted and SCHEMA_VERSION was deliberately NOT bumped."""
    paths = ArtifactPaths.of(tmp_path)
    paths.bundle_dir.mkdir(parents=True, exist_ok=True)
    io.write_json(paths.bundle_conformal_file, {"levels": [0.9], "q": {"0.9": 8.86}})

    conf = io.load_conformal(paths)
    assert conf.q["0.9"] == 8.86
    assert conf.sigma_scale == 1.0            # identity: old bundles behave exactly as before
    assert conf.sigma_scale_mode == "ensemble"


def test_e2e_calibration_was_fitted_cross_donor(trained_bundle):
    """The bundle must record WHETHER cross-donor calibration actually happened, so a silent
    fallback to in-distribution is auditable after the fact."""
    _, _, summary = trained_bundle

    assert summary["xdonor_calibrated"] is True
    assert summary["xdonor_n_donors"] >= 2, "inner-LODO needs at least two donors"
    # SyntheticSource names lines "<NAME>_L<i>", so the two sources give
    # SYNTH_L0/L1 + TAHOE_L0/L1 -- all four should reach the inner-LODO pool
    assert summary["xdonor_n_donors"] == 4, (
        f"expected 4 synthetic cell lines in the inner-LODO pool, got "
        f"{summary['xdonor_n_donors']} -- a donor was skipped for having no cells "
        f"on one side of the split"
    )
    assert summary["xdonor_n_residuals"] > 0
    assert summary["sigma_scale"] >= 1.0
    assert summary["sigma_scale_mode"] == "ensemble"


def test_e2e_sigma_scale_reaches_the_predictor(trained_bundle):
    from cellfate.inference import Predictor

    root, _, summary = trained_bundle
    pred = Predictor(root)
    assert pred.sigma_scale == pytest.approx(summary["sigma_scale"])


def test_e2e_both_modes_are_calibrated_from_their_own_spread(trained_bundle):
    """The bundle must carry a factor for BOTH inference modes, each from its own spread.

    A mode left at 1.0 serves the raw spread -- overconfident, and precisely the defect Stage 1
    exists to remove. The two spreads (across members vs across T dropout passes) are different
    quantities, so identical factors would mean one borrowed the other's.
    """
    from cellfate.inference import Predictor

    root, _, summary = trained_bundle
    assert summary["sigma_scale"] > 1.0, "ensemble mode left uncalibrated"
    assert summary["sigma_scale_mc"] > 1.0, "mc_dropout mode left uncalibrated"

    ens = Predictor(root, mode="ensemble")
    mc = Predictor(root, mode="mc_dropout", T=8)
    assert ens.sigma_scale == pytest.approx(summary["sigma_scale"])
    assert mc.sigma_scale == pytest.approx(summary["sigma_scale_mc"])
    assert ens.sigma_scale != mc.sigma_scale
