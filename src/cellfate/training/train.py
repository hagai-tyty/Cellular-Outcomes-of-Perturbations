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
    huber_age_window,
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


def ensemble_probs(members, ds, device) -> torch.Tensor:
    """Mean over members of ``softmax(member_logits)`` -- exactly ``Predictor``'s ``pbar`` at T=1.

    NOT ``softmax(ensemble_logits(...))``: averaging logits then softmaxing is a different
    quantity by Jensen, and calibrating one while inference produces the other is the fit/apply
    mismatch that Stage 1 run 2 shipped. Every calibrator fitted on P(safe) must use this.
    """
    acc = None
    for m in members:
        p = torch.softmax(member_outputs(m, ds, device)[0], dim=-1)
        acc = p if acc is None else acc + p
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


class _AgeWindow:
    """Stage 1.5.3 C-5 Option 2 — accumulate age cells until a window is worth stepping on.

    THE DEFECT. With HFF's ΔAge labels masked (C-1), 75 age-valid cells sit among 33 688 training
    cells, so uniform shuffling puts ~1.14 of them in a 512-cell batch: ~32% of batches hit the
    hard zero in ``huber_age_loss`` and the survivors carry a Huber over one or two cells, from
    which ``MultiTaskLoss`` then learns ``s_age``. That is noise, not a gradient.

    THE RULE, registered at step 5c before this code existed
    (``plan_tests/register_c5c_bar.py`` -> ``results/register_c5c_bar_results.json``):

        hold age cells back until the window carries >= ``k`` of them, or ``w_max`` batches have
        passed, then take ONE Huber over the whole window.

    WHY A CELL COUNT AND NOT A BATCH COUNT. Step 6 compares two arms and only one is masked. The
    control arm has every cell age-valid (~512 per batch), so a fixed window of W batches would cut
    its age updates 65 -> 8 for no benefit while helping the treatment arm -- the mechanism would
    handicap the control and tilt the result toward the treatment's own conclusion. Triggering on
    ``k`` cells makes the control close at W = 1 on the first batch every time, i.e. **exactly
    today's behaviour**, while the masked arm accumulates. One policy, both arms; it only behaves
    differently because the data differ, which is what a controlled comparison is. Bar A1 asserts
    that identity, and ``tests/test_c5c_age_accumulation.py`` asserts it against the real loop.

    WHY BUFFERED CELLS ARE RE-RUN, NOT RE-USED FROM THEIR ORIGINAL GRAPH. The optimiser steps on
    every batch (the fate task must not be slowed down), so a gradient held from an earlier batch
    would be stale -- computed against parameters that have since moved. The buffer therefore
    stores detached INPUTS and the window re-runs them through the current model, costing one extra
    forward over ~3 cells. Correctness over cleverness; and in the control arm the buffer is always
    empty, so there is no extra forward at all.

    The window CARRIES across the epoch boundary. Forcing a close at each epoch's end was tried
    first (step 5c attempt 1) and manufactured one deliberately-partial window per epoch, which
    was 4.44 pp of a 6.12 pp shortfall and failed bar A2 on its own.
    """

    __slots__ = ("k", "w_max", "x", "fp", "dt", "ya", "batches", "n_closed", "n_short")

    def __init__(self, k: int, w_max: int) -> None:
        self.k, self.w_max = int(k), int(w_max)
        self.x: list[torch.Tensor] = []
        self.fp: list[torch.Tensor] = []
        self.dt: list[torch.Tensor] = []
        self.ya: list[torch.Tensor] = []
        self.batches = 0
        self.n_closed = self.n_short = 0

    @property
    def n_cells(self) -> int:
        return sum(int(t.shape[0]) for t in self.ya)

    def offer(self, x, fp, dt, ya, am) -> bool:
        """Add a batch's age-valid cells; return True when the window should close now.

        The current batch's cells are counted but NOT buffered -- on a close they are consumed
        straight from the live graph, and only on a hold are they detached and kept.
        """
        self.batches += 1
        m = am.bool()
        n_here = int(m.sum())
        if self.n_cells + n_here >= self.k or self.batches >= self.w_max:
            return True
        if n_here:
            self.x.append(x[m].detach())
            self.fp.append(fp[m].detach())
            self.dt.append(dt[m].detach())
            self.ya.append(ya[m].detach())
        return False

    def close(self, model, ag, ya, am, delta: float):
        """One Huber over every age cell in the window: buffered ones re-run, current ones live."""
        m = am.bool()
        preds, trues = [], []
        if self.ya:                       # empty in the control arm -> no extra forward
            _, buf_ag, _ = model(torch.cat(self.x), torch.cat(self.fp), torch.cat(self.dt))
            preds.append(buf_ag)
            trues.append(torch.cat(self.ya))
        if bool(m.any()):
            preds.append(ag[m])
            trues.append(ya[m])
        n = sum(int(t.shape[0]) for t in trues)
        self.n_closed += 1
        self.n_short += int(n < self.k)
        self.reset()
        return huber_age_window(preds, trues, ag, delta)

    def reset(self) -> None:
        self.x.clear()
        self.fp.clear()
        self.dt.clear()
        self.ya.clear()
        self.batches = 0


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

    # C-5 Option 2. `age_window_k <= 1` restores the pre-1.5.3 path exactly -- the documented
    # rollback, and what every test written before this change continues to exercise.
    win = (_AgeWindow(cfg.age_window_k, cfg.age_window_max_batches)
           if getattr(cfg, "age_window_k", 1) > 1 else None)

    best, best_state, bad = float("inf"), None, 0
    for _epoch in range(cfg.epochs):
        model.train()
        for batch in train_dl:   # indexed, not unpacked: the schema grows (donor column)
            x, fp, dt, yc, ya, am = (batch[i].to(device)
                                     for i in (X_I, FP_I, DT_I, YC_I, YA_I, AM_I))
            lg, ag, _ = model(x, fp, dt)
            l_cls = focal_loss(lg, yc, class_w, cfg.focal_gamma)
            if win is None:
                l_age = huber_age_loss(ag, ya, am, cfg.huber_delta)
            elif win.offer(x, fp, dt, ya, am):
                l_age = win.close(model, ag, ya, am, cfg.huber_delta)
            else:
                # held back: the age task contributes nothing to THIS step, exactly as an
                # age-free batch already does today. The fate task steps regardless.
                l_age = ag.sum() * 0.0
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
