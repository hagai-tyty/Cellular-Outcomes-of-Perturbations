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
def _toy_dataset(n=64, g=8):
    rng = np.random.default_rng(3)
    x = torch.tensor(rng.normal(size=(n, g)), dtype=torch.float32)
    fp = torch.tensor(rng.integers(0, 2, size=(n, 2048)), dtype=torch.float32)
    dt = torch.tensor(rng.normal(size=(n, 2)), dtype=torch.float32)
    yc = torch.tensor(rng.dirichlet(np.ones(3), size=n), dtype=torch.float32)
    ya = torch.tensor(rng.normal(size=n), dtype=torch.float32)
    am = torch.ones(n)
    return torch.utils.data.TensorDataset(x, fp, dt, yc, ya, am)


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
    # temperature scaling never worsens NLL
    assert summary["nll_after_temp"] <= summary["nll_before_temp"] + 1e-6
