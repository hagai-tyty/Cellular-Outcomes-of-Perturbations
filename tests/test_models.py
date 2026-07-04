"""Unit tests for cellfate.models (network + losses, Document 3)."""

from __future__ import annotations

import numpy as np
import torch
import torch.nn.functional as F

from cellfate.models import (
    CellFateNet,
    MultiTaskLoss,
    class_balanced_weights,
    focal_loss,
    huber_age_loss,
    mc_dropout_predict,
)


# --------------------------------------------------------------------------- #
# losses                                                                      #
# --------------------------------------------------------------------------- #
def test_focal_reduces_to_soft_cross_entropy_at_gamma_zero():
    torch.manual_seed(0)
    logits = torch.randn(5, 3)
    target = F.softmax(torch.randn(5, 3), dim=1)
    w = torch.ones(3)
    fl0 = focal_loss(logits, target, w, gamma=0.0)
    ce = -(target * F.log_softmax(logits, dim=1)).sum(1).mean()
    assert torch.allclose(fl0, ce, atol=1e-6)


def test_focal_is_lower_when_prediction_matches_target():
    target = torch.tensor([[1.0, 0.0, 0.0]])
    good = torch.tensor([[10.0, 0.0, 0.0]])
    bad = torch.tensor([[0.0, 0.0, 10.0]])
    w = torch.ones(3)
    assert focal_loss(good, target, w, 2.0) < focal_loss(bad, target, w, 2.0)


def test_class_balanced_weights_upweight_rare_classes():
    w = class_balanced_weights(np.array([1000.0, 1000.0, 50.0]), beta=0.999)
    assert w[2] > w[0] and w[2] > w[1]            # the rare class gets more weight
    assert np.isclose(w.mean(), 1.0, atol=1e-6)   # normalised to mean 1


def test_huber_age_loss_masks_invalid_cells():
    pred = torch.tensor([1.0, 2.0, 3.0, 4.0], requires_grad=True)
    true = torch.tensor([1.5, 2.5, 0.0, 0.0])
    none_valid = huber_age_loss(pred, true, torch.zeros(4), delta=2.0)
    assert float(none_valid.detach()) == 0.0
    # only the first two (valid) cells contribute
    half = huber_age_loss(pred, true, torch.tensor([1.0, 1.0, 0.0, 0.0]), delta=2.0)
    ref = F.huber_loss(pred[:2], true[:2], delta=2.0)
    assert torch.allclose(half, ref)


def test_multitask_loss_is_scalar_and_learns_weights():
    mtl = MultiTaskLoss()
    out = mtl(torch.tensor(1.0), torch.tensor(2.0))
    assert out.ndim == 0
    out.backward()
    assert mtl.log_var_cls.grad is not None and mtl.log_var_age.grad is not None


# --------------------------------------------------------------------------- #
# network                                                                     #
# --------------------------------------------------------------------------- #
def _net():
    return CellFateNet(g=16, n_fp=2048, n_dt=2, d_cell=16, d_u=16, latent_dim=16, p_drop=0.3)


def _inputs(n=4, g=16):
    torch.manual_seed(1)
    return (torch.randn(n, g), torch.randint(0, 2, (n, 2048)).float(), torch.randn(n, 2))


def test_forward_shapes():
    net = _net().eval()
    logits, age, feat = net(*_inputs())
    assert logits.shape == (4, 3)
    assert age.shape == (4,)
    assert feat.shape == (4, 16)


def test_mc_dropout_varies_and_restores_eval_mode():
    net = _net().eval()
    probs, ages = mc_dropout_predict(net, *_inputs(), n_samples=8)
    assert probs.shape == (8, 4, 3)
    assert torch.allclose(probs.sum(-1), torch.ones(8, 4), atol=1e-5)
    assert probs.std(0).sum() > 1e-6        # dropout produced variation
    assert not net.training                 # eval mode restored


def test_member_roundtrip_reproduces_outputs(tmp_path):
    net = _net().eval()
    x = _inputs()
    out1 = net(*x)[0].detach()
    path = tmp_path / "m.pt"
    net.save_member(path)
    net2 = CellFateNet.load_member(path)
    out2 = net2(*x)[0].detach()
    assert torch.allclose(out1, out2, atol=1e-6)
    assert net2.arch["g"] == 16
