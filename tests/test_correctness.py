"""Rigorous / adversarial correctness tests for Document 3.

These go beyond "it runs": they assert the network actually *learns*, the
calibration statistics hold *out-of-sample*, OOD genuinely *separates*
distributions, the age mask cannot *leak*, the Kendall weighting *responds* to a
noisy task, the scalers are *really applied*, temperature *recovers* a known
mis-scaling, and ``run`` is *reproducible*.
"""

from __future__ import annotations

import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from torch.utils.data import TensorDataset

from cellfate.models import (
    CellFateNet,
    MultiTaskLoss,
    class_balanced_weights,
    focal_loss,
    huber_age_loss,
    mc_dropout_predict,
)
from cellfate.training import (
    TrainConfig,
    coverage,
    fit_conformal,
    fit_ood,
    fit_temperature,
    mahalanobis,
    member_outputs,
    train_member,
)
from cellfate.training.dataset import AM_I, YA_I, YC_I

G = 24
# Planted task weights are FIXED across splits — train and test must share the
# same labelling function, otherwise generalisation is impossible by construction.
_W_CLS = np.random.default_rng(100).normal(size=(G, 3))
_W_AGE = np.random.default_rng(101).normal(size=G)


def _planted_dataset(n, seed):
    """A linearly-separable 3-class task + a linear age signal, both functions of X.

    Fingerprint and dose/time are held constant so the only signal is in the cell
    pathway -- if the network learns, the whole encoder->trunk->head stack works.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, G)).astype(np.float32)
    cls = (x @ _W_CLS).argmax(1)
    yc = np.eye(3, dtype=np.float32)[cls]
    age = (x @ _W_AGE).astype(np.float32)
    age = (age - age.mean()) / (age.std() + 1e-8)
    fp = np.zeros((n, 2048), np.float32)
    dt = np.zeros((n, 2), np.float32)
    ds = TensorDataset(
        torch.from_numpy(x), torch.from_numpy(fp), torch.from_numpy(dt),
        torch.from_numpy(yc), torch.from_numpy(age), torch.ones(n),
    )
    return ds, cls, age


def _cfg(**kw):
    base = dict(dataset_dir="unused", d_cell=32, d_u=8, latent_dim=32, p_drop=0.05,
                lr=3e-3, epochs=40, patience=40, batch_size=64, ensemble_size=1,
                focal_gamma=2.0, huber_delta=2.0)
    base.update(kw)
    return TrainConfig(**base)


def _make(g=G, **kw):
    p = dict(g=g, d_cell=32, d_u=8, latent_dim=32, p_drop=0.05)
    p.update(kw)
    return lambda: CellFateNet(**p)


# --------------------------------------------------------------------------- #
# 1. the network actually LEARNS a planted signal                             #
# --------------------------------------------------------------------------- #
def test_network_learns_classification_and_age():
    train_ds, _, _ = _planted_dataset(700, seed=0)
    test_ds, test_cls, test_age = _planted_dataset(300, seed=1)

    model, _ = train_member(_make(), train_ds, test_ds, _cfg(), seed=0, device="cpu")
    logits, age, _ = member_outputs(model, test_ds, "cpu")

    acc = (logits.argmax(1).numpy() == test_cls).mean()
    assert acc > 0.65, f"3-class separable task only reached acc={acc:.2f} (chance=0.33)"

    r = float(np.corrcoef(age.numpy(), test_age)[0, 1])
    assert r > 0.5, f"age head failed to track the linear signal (r={r:.2f})"


def test_training_reduces_loss_vs_untrained():
    train_ds, _, _ = _planted_dataset(400, seed=2)
    cw = torch.tensor(class_balanced_weights(train_ds.tensors[YC_I].numpy().sum(0)))

    def ce(model):
        lg, ag, _ = member_outputs(model, train_ds, "cpu")
        c = focal_loss(lg, train_ds.tensors[YC_I], cw, 2.0)
        a = huber_age_loss(ag, train_ds.tensors[YA_I], train_ds.tensors[AM_I], 2.0)
        return float(c + a)

    torch.manual_seed(0)
    untrained = CellFateNet(**{"g": G, "d_cell": 32, "d_u": 8, "latent_dim": 32, "p_drop": 0.05})
    trained, _ = train_member(_make(), train_ds, train_ds, _cfg(epochs=30), seed=0, device="cpu")
    assert ce(trained) < ce(untrained)


# --------------------------------------------------------------------------- #
# 2. Kendall multi-task weighting responds to a noisy task                     #
# --------------------------------------------------------------------------- #
def test_kendall_downweights_the_higher_loss_task():
    # constant losses: age task is consistently harder than the cls task
    mtl = MultiTaskLoss()
    opt = torch.optim.Adam(mtl.parameters(), lr=0.05)
    l_cls, l_age = torch.tensor(0.1), torch.tensor(2.0)
    for _ in range(2000):
        opt.zero_grad()
        mtl(l_cls, l_age).backward()
        opt.step()
    # Optima follow the asymmetric weighting: the classification term
    # exp(-s)*L + 0.5*s is minimised at s = ln(2L), while the regression term
    # 0.5*exp(-s)*L + 0.5*s is minimised at s = ln(L). Either way, the larger
    # loss -> larger log-variance -> the task is down-weighted.
    assert mtl.log_var_age.item() > mtl.log_var_cls.item()
    assert abs(mtl.log_var_age.item() - np.log(2.0)) < 0.15        # age: ln(L_age)
    assert abs(mtl.log_var_cls.item() - np.log(2 * 0.1)) < 0.15    # cls: ln(2*L_cls)


# --------------------------------------------------------------------------- #
# 3. the age mask is leak-proof                                                #
# --------------------------------------------------------------------------- #
def test_masked_ages_cannot_leak_into_loss_or_gradient():
    pred = torch.randn(8, requires_grad=True)
    true = torch.randn(8)
    mask = torch.tensor([1., 1., 1., 0., 0., 0., 0., 0.])

    clean = huber_age_loss(pred, true, mask, 2.0)
    poisoned = true.clone()
    poisoned[mask == 0] = 1e9                      # garbage in the masked cells
    after = huber_age_loss(pred, poisoned, mask, 2.0)
    assert torch.allclose(clean, after)            # value unaffected

    after.backward()
    assert torch.allclose(pred.grad[mask == 0], torch.zeros(5))  # no gradient leak


# --------------------------------------------------------------------------- #
# 4. conformal coverage holds OUT-OF-SAMPLE (the real guarantee)              #
# --------------------------------------------------------------------------- #
def test_conformal_coverage_holds_on_fresh_data():
    rng = np.random.default_rng(7)
    calib = np.abs(rng.normal(size=2000))
    fresh = np.abs(rng.normal(size=8000))          # exchangeable, never seen at fit time
    for level in (0.80, 0.90, 0.95):
        q = fit_conformal(calib, [level]).q[str(level)]
        cov = coverage(fresh, q)
        assert level - 0.04 <= cov <= level + 0.04, f"level {level}: out-of-sample cov {cov:.3f}"


def test_conformal_interval_widens_with_level():
    rng = np.random.default_rng(8)
    res = np.abs(rng.normal(size=1000))
    cp = fit_conformal(res, [0.80, 0.90, 0.95])
    assert cp.q["0.8"] <= cp.q["0.9"] <= cp.q["0.95"]


# --------------------------------------------------------------------------- #
# 5. OOD genuinely separates distributions (AUROC)                            #
# --------------------------------------------------------------------------- #
def test_ood_separates_shifted_and_scaled_distributions():
    rng = np.random.default_rng(9)
    d = 16
    ood = fit_ood(rng.normal(size=(1000, d)))

    in_test = rng.normal(size=(1000, d))
    shifted = rng.normal(size=(1000, d)) + 3.0           # mean shift
    scaled = rng.normal(size=(1000, d)) * 3.0            # variance inflation

    for name, out in (("shift", shifted), ("scale", scaled)):
        scores = np.concatenate([mahalanobis(in_test, ood), mahalanobis(out, ood)])
        labels = np.concatenate([np.zeros(1000), np.ones(1000)])
        auroc = roc_auc_score(labels, scores)
        assert auroc > 0.9, f"OOD ({name}) AUROC only {auroc:.3f}"


# --------------------------------------------------------------------------- #
# 6. temperature RECOVERS a known mis-scaling                                 #
# --------------------------------------------------------------------------- #
def _calibrated_logits(n, seed):
    rng = np.random.default_rng(seed)
    p = rng.dirichlet(np.ones(3) * 2.0, size=n)
    labels = np.array([rng.choice(3, p=row) for row in p])
    return np.log(p + 1e-9), np.eye(3)[labels]


def test_temperature_recovers_overconfidence_and_underconfidence():
    logits, target = _calibrated_logits(3000, seed=11)
    # over-confident by 3x -> T should come back near 3; under-confident by 3x -> near 1/3
    t_over = fit_temperature(3.0 * logits, target).temperature
    t_under = fit_temperature(logits / 3.0, target).temperature
    assert 2.0 < t_over < 4.5, f"expected ~3, got {t_over:.2f}"
    assert 0.2 < t_under < 0.55, f"expected ~0.33, got {t_under:.2f}"


# --------------------------------------------------------------------------- #
# 7. scalers are REALLY applied by the loader (via the data pipeline)          #
# --------------------------------------------------------------------------- #
def test_loader_actually_standardizes_train_features(tmp_path):
    from cellfate.common.io import ArtifactPaths
    from cellfate.common.scalers import Scalers
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build
    from cellfate.training import load_split_tensors

    build(
        DataConfig(out=str(tmp_path), gene_panel=str(tmp_path / "panel.json"), n_genes=120,
                   qc=QCConfig(min_genes=5, max_mito_frac=0.5), label_tau=0.5,
                   split_fracs=(0.7, 0.15, 0.1, 0.05), primary_regime="scaffold", seed=0),
        sources=[SyntheticSource(name="synth", n_lines=2, n_compounds=4,
                                 n_cells_per_condition=12, seed=1)],
    )
    paths = ArtifactPaths.of(tmp_path)
    sc = Scalers.load(paths.scalers_file)
    ds = load_split_tensors(paths, sc, "scaffold", "train")
    x = ds.tensors[0].numpy()
    # scalers were fit on the train split, so standardized train features are ~N(0,1)
    assert abs(x.mean()) < 0.15
    assert 0.7 < x.std() < 1.3
    # dose/time column standardized too
    dt = ds.tensors[1 + 1].numpy()  # DT column index is 2
    assert abs(dt.mean()) < 0.4

    # -- Stage 1a: the donor column, without which inner-LODO calibration is impossible.
    # Sourced from `cell_line`; SyntheticSource above built n_lines=2, so we expect 2 codes.
    from cellfate.training.dataset import DONOR_I

    donor = ds.tensors[DONOR_I]
    assert len(ds.tensors) == 7, "donor column missing from the training tensors"
    assert donor.dtype == torch.long, "donor codes must be integer"
    assert len(donor) == len(ds.tensors[0]), "donor column length mismatch"
    assert len(set(donor.tolist())) == 2, "expected one code per synthetic cell line"

    # the empty-split branch builds its own tensors -- it must grow the column too
    empty = load_split_tensors(paths, sc, "scaffold", "__no_such_split__")
    assert len(empty.tensors) == 7, "empty-split branch still returns 6 columns"
    assert len(empty) == 0


# --------------------------------------------------------------------------- #
# 8. run() is reproducible end-to-end                                         #
# --------------------------------------------------------------------------- #
def test_run_is_reproducible(tmp_path):
    from cellfate.common.io import ArtifactPaths
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build
    from cellfate.training import run as train_run

    ds_dir = tmp_path / "ds"
    build(
        DataConfig(out=str(ds_dir), gene_panel=str(ds_dir / "panel.json"), n_genes=120,
                   qc=QCConfig(min_genes=5, max_mito_frac=0.5), label_tau=0.5,
                   split_fracs=(0.6, 0.2, 0.1, 0.1), primary_regime="scaffold", seed=0),
        sources=[SyntheticSource(name="synth", n_lines=2, n_compounds=4,
                                 n_cells_per_condition=10, seed=1),
                 SyntheticSource(name="tahoe", n_lines=2, n_compounds=4,
                                 n_cells_per_condition=10, seed=2)],
    )
    common = dict(dataset_dir=str(ds_dir), regime="scaffold", d_cell=16, d_u=8,
                  latent_dim=16, p_drop=0.1, epochs=4, patience=4, batch_size=64,
                  ensemble_size=2, device="cpu")
    out_a, out_b = tmp_path / "a", tmp_path / "b"
    s1 = train_run(TrainConfig(out=str(out_a), **common))
    s2 = train_run(TrainConfig(out=str(out_b), **common))

    assert s1["temperature"] == s2["temperature"]
    assert s1["conformal_q"] == s2["conformal_q"]
    m1 = CellFateNet.load_member(ArtifactPaths.of(out_a).bundle_members_dir / "member_0.pt")
    m2 = CellFateNet.load_member(ArtifactPaths.of(out_b).bundle_members_dir / "member_0.pt")
    w1 = m1.cell.net[0][0].weight.detach()
    w2 = m2.cell.net[0][0].weight.detach()
    assert torch.allclose(w1, w2, atol=1e-6)


# --------------------------------------------------------------------------- #
# 9. numerical stability + MC-dropout edge behaviour                          #
# --------------------------------------------------------------------------- #
def test_focal_loss_stable_under_extreme_logits():
    logits = torch.tensor([[1e3, -1e3, -1e3], [-1e3, -1e3, 1e3]])  # confident, 2nd is wrong
    target = torch.tensor([[1., 0., 0.], [1., 0., 0.]])
    loss = focal_loss(logits, target, torch.ones(3), 2.0)
    assert torch.isfinite(loss)


def test_eval_forward_deterministic_but_mc_varies():
    net = CellFateNet(g=G, d_cell=16, d_u=8, latent_dim=16, p_drop=0.3).eval()
    x = (torch.randn(5, G), torch.randint(0, 2, (5, 2048)).float(), torch.randn(5, 2))
    a = net(*x)[0].detach()
    b = net(*x)[0].detach()
    assert torch.allclose(a, b)                                  # eval = deterministic
    probs, _ = mc_dropout_predict(net, *x, n_samples=1)
    assert probs.shape == (1, 5, 3)
    probs_many, _ = mc_dropout_predict(net, *x, n_samples=30)
    assert probs_many.std(0).sum() > 0                           # dropout = stochastic


# --------------------------------------------------------------------------- #
# 10. robustness on degenerate inputs (regression guards)                     #
# --------------------------------------------------------------------------- #
def test_member_outputs_handles_empty_dataset():
    from torch.utils.data import TensorDataset
    net = CellFateNet(g=G, d_cell=16, d_u=8, latent_dim=16, p_drop=0.0)
    empty = TensorDataset(torch.empty(0, G), torch.empty(0, 2048), torch.empty(0, 2),
                          torch.empty(0, 3), torch.empty(0), torch.empty(0))
    logits, ages, feats = member_outputs(net, empty, "cpu")
    assert logits.shape == (0, 3)
    assert ages.shape == (0,)
    assert feats.shape == (0, 16)


def test_conformal_warns_when_level_unattainable_for_n():
    import warnings as _w
    # n=5 cannot guarantee 0.99 coverage (max ~5/6); must warn, not silently under-cover
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        cp = fit_conformal(np.abs(np.random.default_rng(0).normal(size=5)), [0.99])
    assert any("calibration points" in str(c.message) for c in caught)
    assert np.isfinite(cp.q["0.99"])  # still returns the widest finite interval


def test_focal_soft_label_reduces_to_calibrated_ce_at_gamma_zero():
    # the documented escape hatch: gamma=0 gives exact soft CE (min at p=target)
    tgt = torch.tensor([[0.6, 0.3, 0.1]])
    lg = torch.log(tgt.clamp_min(1e-6)).clone().requires_grad_(True)
    focal_loss(lg, tgt, torch.ones(3), gamma=0.0).backward()
    assert lg.grad.norm() < 1e-5   # gradient ~0 at p=target when gamma=0
