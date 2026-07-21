"""Cross-donor calibration statistics via inner leave-one-donor-out (Stage 1b).

The bundle's calibrators are fitted on val/calib/train -- all from donors the model trained
alongside. Out-of-donor they fail: conformal coverage 0.40 vs a nominal 0.90 (0.00 on N2 and N3),
fate ECE 0.28, OOD AUC 0.47. One architectural mistake -- calibrating against donors already
seen -- with four manifestations.

This module produces statistics from the regime deployment actually faces: for each donor in the
TRAINING set, train on the others and predict on it, then pool the results. `temperature`, `q`
and `sigma_scale` are fitted on those.

THREE of the four, not four. The OOD reference is deliberately NOT fitted here -- see the note on
`XDonorStats.feats` and the call site in `train_model.run`. Its refit is deferred to the Stage 3d
decision about whether to keep the gate at all.

COST. One extra ensemble per training donor -- 5 for a 6-donor LOOCV fold, so roughly 6x the
training time. The inner ensembles must use the SAME `ensemble_size` as the deployed model,
because `sigma_scale` calibrates the ensemble spread and a spread over k members is a different
quantity for a different k. There is deliberately no knob to shrink them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch.utils.data import TensorDataset

from cellfate.common.logging import get_logger, log_event

from .dataset import AM_I, DONOR_I, YA_I, YC_I
from .train import ensemble_age, ensemble_logits, member_outputs, train_ensemble

log = get_logger("cellfate.training")

MIN_DONORS = 2

# An inner fold is only a proxy for the DEPLOYED model if it is trained on a comparable amount
# of data. Holding out a donor that is most of the training set produces a data-starved model
# whose residuals measure "what happens with almost no training data", not "what happens on an
# unseen donor" -- and because those residuals are pooled, one such fold can swamp every honest
# one. Measured on the Gill+HFF merge: holding out HFF left 75 of 33,688 cells (0.2%), and that
# single fold supplied 33,613 of 33,688 pooled residuals (99.8%).
MIN_INNER_TRAIN_FRAC = 0.5

# The inference mode `sigma_scale` is valid for. NOT a configurable value: `sigma_pred` is
# collected as the spread across ENSEMBLE MEMBERS, so an ensemble-derived factor is the only
# thing this module can produce. It is a constant so the label written into the bundle always
# describes what was actually computed, rather than whatever the caller declared.
SIGMA_SCALE_MODE = "ensemble"


@dataclass
class XDonorStats:
    abs_residuals: np.ndarray   # (M,)  |ΔAge error| pooled over held-out donors -> fits q
    logits: np.ndarray          # (M,3) fate logits, out-of-donor -> fits temperature
    targets: np.ndarray         # (M,3) matching soft labels
    sigma_pred: np.ndarray      # (M,)  ENSEMBLE spread     -> fits sigma_scale
    sigma_pred_mc: np.ndarray   # (M,)  MC-DROPOUT spread   -> fits sigma_scale_mc
    n_donors: int               # inner-LODO donors actually used
    # DIAGNOSTIC ONLY -- deliberately not used to fit the OOD reference. These come from
    # independently-seeded INNER models, whose latent bases differ by arbitrary rotation, while
    # OODDetector compares the DEPLOYED member[0]'s z against the stored Gaussian. Pooling them
    # would make the Mahalanobis distance meaningless. See train_model.run and the Stage 1
    # deviation log.
    feats: np.ndarray           # (M,D) trunk features, out-of-donor

    @property
    def has_residuals(self) -> bool:
        return self.abs_residuals.size > 0

    @property
    def has_logits(self) -> bool:
        return self.logits.shape[0] > 0


def _subset(ds: TensorDataset, mask: torch.Tensor) -> TensorDataset:
    return TensorDataset(*[t[mask] for t in ds.tensors])


MC_ROW_BUDGET = 8192          # tiled rows per forward; batch = budget // T


@torch.no_grad()
def mc_dropout_spread(model, ds: TensorDataset, device: str, T: int,
                      row_budget: int = MC_ROW_BUDGET) -> np.ndarray:
    """Std of ΔAge across T dropout passes -- what `sigma_age` IS in mode="mc_dropout".

    Mirrors ``Predictor._raw_batch``'s mc_dropout branch exactly: only Dropout modules go to
    train mode, the T passes are ONE tiled forward, and the spread is ``std(0, unbiased=False)``.
    A mismatch here would calibrate a quantity inference never produces.

    ⚠ T COUPLING. ``std(0, unbiased=False)`` over T samples is biased LOW, and the bias grows as
    T shrinks (~4% at T=8, <1% at T=50). The factor is therefore only exact for the T it was
    fitted at -- ``TrainConfig.mc_dropout_T``, recorded in the bundle metrics. Running
    ``Predictor(T=...)`` far from that value shifts sigma by a few percent: second-order against
    the ~5x miscalibration this exists to remove, but real. Keep them equal when it matters.

    The tiled input is T x batch rows, so the batch is sized from a ROW budget rather than
    fixed -- a fixed batch of 256 at T=50 is 12,800 rows, which can exhaust a small GPU on a
    large held-out fold.
    """
    from torch.utils.data import DataLoader

    from .dataset import DT_I, FP_I, X_I

    model.eval()
    for m in model.modules():                      # dropout ON, everything else stays eval
        if isinstance(m, torch.nn.Dropout):
            m.train()
    try:
        out = []
        for batch in DataLoader(ds, batch_size=max(1, row_budget // max(T, 1))):
            x, fp, dt = (batch[i].to(device) for i in (X_I, FP_I, DT_I))
            n = x.shape[0]
            _, ag, _ = model(x.repeat(T, 1), fp.repeat(T, 1), dt.repeat(T, 1))
            out.append(ag.view(T, n).std(0, unbiased=False).cpu())
    finally:
        model.eval()
    return torch.cat(out).numpy() if out else np.array([])


def n_train_donors(train_ds: TensorDataset) -> int:
    """How many distinct donors the training split carries (the inner-LODO precondition)."""
    if len(train_ds) == 0:
        return 0
    return len(set(train_ds.tensors[DONOR_I].tolist()))


def crossdonor_stats(train_ds: TensorDataset, val_ds: TensorDataset,
                     make_model, cfg, device: str, mc_T: int = 50) -> XDonorStats:
    """Inner leave-one-donor-out over the training donors; pool the out-of-donor statistics.

    For each donor d: train an ensemble on the other donors, predict on d, keep the ΔAge
    residuals, fate logits, trunk features and ensemble spread. Pool across d.

    The monitoring split passed to each inner ensemble is the outer val split with donor d
    REMOVED -- see the leakage note in the module tests and STAGE_1 deviation log.
    """
    donors = train_ds.tensors[DONOR_I]
    uniq = sorted(set(donors.tolist()))
    if len(uniq) < MIN_DONORS:
        raise ValueError(
            f"inner-LODO needs >={MIN_DONORS} training donors, found {len(uniq)}. "
            "Check the donor column (Stage 1a) rather than bypassing this."
        )

    n_total = len(train_ds)
    res, log_, tgt, fts, sig, sig_mc = [], [], [], [], [], []
    used, skipped = 0, []
    for d in uniq:
        hold = donors == d
        inner_tr, inner_te = _subset(train_ds, ~hold), _subset(train_ds, hold)
        if len(inner_te) == 0 or len(inner_tr) == 0:
            continue

        # A donor that IS the training set cannot be held out and still leave a model worth
        # calibrating against. Skipping is the only honest option: including it pools residuals
        # from a data-starved model with residuals from real ones, and since it is also the
        # largest fold it dominates the pooled quantile.
        if len(inner_tr) < MIN_INNER_TRAIN_FRAC * n_total:
            skipped.append(int(d))
            log.warning(
                "inner-LODO: SKIPPING donor %d -- holding it out leaves %d of %d training "
                "cells (%.1f%%, below the %.0f%% floor). It is a bulk corpus, not a donor; "
                "calibrating against it would measure data starvation, not donor shift.",
                int(d), len(inner_tr), n_total, 100.0 * len(inner_tr) / max(n_total, 1),
                100.0 * MIN_INNER_TRAIN_FRAC,
            )
            continue

        # Early stopping must NOT see the held-out donor, or the residuals we collect from it
        # are best-case rather than honest out-of-donor. Deployment gets no such privilege.
        inner_val = (_subset(val_ds, val_ds.tensors[DONOR_I] != d) if len(val_ds) else val_ds)

        members, _ = train_ensemble(make_model, inner_tr, inner_val, cfg, device)
        used += 1

        # -- ΔAge residuals on the held-out donor (age-valid rows only) --
        age = ensemble_age(members, inner_te, device).numpy()
        ya = inner_te.tensors[YA_I].numpy()
        am = inner_te.tensors[AM_I].numpy().astype(bool)
        if am.any():
            res.append(np.abs(age[am] - ya[am]))
            per = np.stack([member_outputs(m, inner_te, device)[1].numpy() for m in members])
            # ddof=0, matching Predictor's age.std(0, unbiased=False)
            sig.append(per.std(axis=0)[am])
            # ...and the SAME rows under mc_dropout, so that mode gets its own honest factor
            # instead of borrowing one calibrated for a different spread. Dropout consumes RNG,
            # but train_member re-seeds at the top of every member, so the next fold is
            # unaffected and `run()` stays reproducible.
            sig_mc.append(mc_dropout_spread(members[0], inner_te, device, mc_T)[am])

        log_.append(ensemble_logits(members, inner_te, device).numpy())
        tgt.append(inner_te.tensors[YC_I].numpy())
        fts.append(member_outputs(members[0], inner_te, device)[2].numpy())

        log_event(log, "xdonor.fold", donor=int(d), n_train=len(inner_tr),
                  n_held=len(inner_te), n_age=int(am.sum()))

    if used < MIN_DONORS:
        raise ValueError(
            f"inner-LODO ran only {used} usable fold(s) of {len(uniq)} donors "
            f"(skipped as bulk corpora: {skipped}). Calibration from fewer than {MIN_DONORS} "
            "held-out donors is not cross-donor. Check the donor column (Stage 1a)."
        )

    stats = XDonorStats(
        abs_residuals=np.concatenate(res) if res else np.array([]),
        logits=np.vstack(log_) if log_ else np.zeros((0, 3)),
        targets=np.vstack(tgt) if tgt else np.zeros((0, 3)),
        feats=np.vstack(fts) if fts else np.zeros((0, 1)),
        sigma_pred=np.concatenate(sig) if sig else np.array([]),
        sigma_pred_mc=np.concatenate(sig_mc) if sig_mc else np.array([]),
        n_donors=used,
    )
    log_event(log, "xdonor.done", n_donors=used, n_skipped=len(skipped), skipped=skipped,
              n_residuals=int(stats.abs_residuals.size), n_logits=int(stats.logits.shape[0]))
    return stats


def sigma_scale_factor(stats: XDonorStats, z_conf: float, level: float = 0.90,
                       mode: str = "ensemble") -> float:
    """Multiplier s such that ``mu +- z_conf*(s*sigma_age)`` attains nominal coverage.

    WHY THIS EXISTS, and why refitting the conformal `q` alone does not replace it
    (MASTER_PLAN S5b-bis): `sigma_age` is the ENSEMBLE SPREAD (~2.4 yr) while the true
    out-of-donor error is ~14 yr. RES consumes sigma via ``R_eff = max(0, -(mu + z*sigma))``,
    NOT the conformal q -- so a refit of q fixes the reported interval and leaves RES exactly
    as broken. The ensemble agrees with itself while being collectively wrong, which is the
    classic failure of ensemble uncertainty under distribution shift.

    MODE DEPENDENCY. ``Predictor._raw_batch`` produces sigma_age two different ways:
        mode="ensemble"   (DEFAULT) -> spread across ensemble members
        mode="mc_dropout"           -> spread across T dropout passes of member[0]
    These are different quantities of different magnitude, so EACH GETS ITS OWN FACTOR,
    measured on the same held-out rows. Borrowing one for the other calibrates the wrong
    spread; serving an uncalibrated one is the very defect Stage 1 exists to remove. Both are
    written to the bundle and `Predictor` selects by mode.

    Clamped at 1.0: this may widen sigma, never shrink it. Over-wide makes RES conservative,
    which is the safe direction.
    """
    spread = stats.sigma_pred if mode == "ensemble" else stats.sigma_pred_mc
    if stats.abs_residuals.size == 0 or spread.size == 0:
        return 1.0
    need = float(np.quantile(stats.abs_residuals, level))    # half-width for `level` coverage
    have = float(np.median(spread)) * float(z_conf)
    if not np.isfinite(need) or not np.isfinite(have) or have <= 0:
        return 1.0
    return max(1.0, need / have)
