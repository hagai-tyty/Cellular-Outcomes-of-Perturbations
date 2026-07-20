"""Training loop + ensembling (Document 3, S4).

Each ensemble member is an independently-seeded network trained with the
class-balanced focal + masked-Huber objective, balanced by the Kendall
multi-task weighting, with gradient clipping and early stopping on the val split.
Deep ensembling (training K members) is the backbone of the epistemic-uncertainty
estimate; MC-dropout (in ``cellfate.models``) layers on top at inference.
"""

from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

from cellfate.common.logging import get_logger, log_event
from cellfate.common.seeding import set_global_seed
from cellfate.models import (
    CellFateNet,
    MultiTaskLoss,
    class_balanced_weights,
    focal_loss,
    huber_age_loss,
)

from .dataset import AM_I, DT_I, FP_I, X_I, YA_I, YC_I, loader

log = get_logger("cellfate.training")


def class_mass(ds: TensorDataset) -> np.ndarray:
    """Per-class soft-label mass over a dataset (drives class-balanced weights)."""
    if len(ds) == 0:
        return np.ones(3, dtype=np.float64)
    return ds.tensors[YC_I].numpy().sum(axis=0)


@torch.no_grad()
def member_outputs(model: CellFateNet, ds: TensorDataset, device: str,
                   batch_size: int = 2048):
    """Eval-mode (dropout OFF) outputs over a dataset: (logits, age, trunk_feature)."""
    model.eval()
    if len(ds) == 0:  # empty split -> correctly-shaped empty tensors
        return (torch.empty(0, model.arch["n_classes"]), torch.empty(0),
                torch.empty(0, model.arch["latent_dim"]))
    logits, ages, feats = [], [], []
    for x, fp, dt, *_ in DataLoader(ds, batch_size=batch_size):
        lg, ag, ft = model(x.to(device), fp.to(device), dt.to(device))
        logits.append(lg.cpu())
        ages.append(ag.cpu())
        feats.append(ft.cpu())
    return torch.cat(logits), torch.cat(ages), torch.cat(feats)


def ensemble_logits(members, ds, device) -> torch.Tensor:
    acc = None
    for m in members:
        lg = member_outputs(m, ds, device)[0]
        acc = lg if acc is None else acc + lg
    return acc / len(members)


def ensemble_age(members, ds, device) -> torch.Tensor:
    acc = None
    for m in members:
        ag = member_outputs(m, ds, device)[1]
        acc = ag if acc is None else acc + ag
    return acc / len(members)


def _eval_loss(model, dl, class_w, gamma, huber_delta, device) -> float:
    """Fixed-weight validation objective for model selection.

    Deliberately Kendall-free: the multitask loss's learned log-variances can
    lower the loss *number* by inflating uncertainty, which would let early
    stopping / best-checkpoint selection be gamed independently of predictive
    quality. A fixed-weight focal + Huber objective tracks quality instead.
    """
    model.eval()
    tot, n = 0.0, 0
    with torch.no_grad():
        for batch in dl:   # indexed, not unpacked: the schema grows (donor column)
            x, fp, dt = batch[X_I], batch[FP_I], batch[DT_I]
            yc, ya, am = batch[YC_I], batch[YA_I], batch[AM_I]
            lg, ag, _ = model(x.to(device), fp.to(device), dt.to(device))
            l_cls = focal_loss(lg, yc.to(device), class_w, gamma)
            l_age = huber_age_loss(ag, ya.to(device), am.to(device), huber_delta)
            tot += (l_cls + l_age).item() * x.size(0)
            n += x.size(0)
    return tot / max(n, 1)


def train_member(make_model, train_ds, val_ds, cfg, seed: int, device: str):
    """Train one member; return (model in eval mode, best monitored loss)."""
    set_global_seed(seed)
    model = make_model().to(device)
    mtl = MultiTaskLoss().to(device)
    params = list(model.parameters()) + list(mtl.parameters())
    opt = torch.optim.Adam(params, lr=cfg.lr, weight_decay=cfg.wd)
    class_w = torch.tensor(class_balanced_weights(class_mass(train_ds), cfg.class_weight_beta),
                           device=device)

    train_dl = loader(train_ds, cfg.batch_size, shuffle=True)
    monitor_dl = loader(val_ds, cfg.batch_size, shuffle=False) if len(val_ds) else train_dl

    best, best_state, bad = float("inf"), None, 0
    for _epoch in range(cfg.epochs):
        model.train()
        for batch in train_dl:   # indexed, not unpacked: the schema grows (donor column)
            x, fp, dt, yc, ya, am = (batch[i].to(device)
                                     for i in (X_I, FP_I, DT_I, YC_I, YA_I, AM_I))
            lg, ag, _ = model(x, fp, dt)
            l_cls = focal_loss(lg, yc, class_w, cfg.focal_gamma)
            l_age = huber_age_loss(ag, ya, am, cfg.huber_delta)
            loss = mtl(l_cls, l_age)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, cfg.clip)
            opt.step()
        cur = _eval_loss(model, monitor_dl, class_w, cfg.focal_gamma, cfg.huber_delta, device)
        if cur < best - cfg.min_delta:
            best, bad = cur, 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= cfg.patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model, best


def train_ensemble(make_model, train_ds, val_ds, cfg, device):
    """Train ``cfg.ensemble_size`` independently-seeded members."""
    members, val_losses = [], []
    for i in range(cfg.ensemble_size):
        model, vloss = train_member(make_model, train_ds, val_ds, cfg, cfg.base_seed + i, device)
        members.append(model)
        val_losses.append(vloss)
        log_event(log, "member.trained", idx=i, seed=cfg.base_seed + i, val_loss=round(vloss, 5))
    return members, val_losses
