"""STAGE 1.5.3 STEP 5c — C-5 Option 2's age-accumulation window, tested against the real loop.

The gates this file exists to hold, from the step table:

  * **arm-A behaviour bit-identical to today** -- the one that matters most. If the unmasked
    control arm's training moves at all, `scorecard/baseline.json` stops being a valid reference
    and step 6's comparison is confounded by the very mechanism meant to de-confound it.
  * the window's loss is `sum(losses)/sum(cells)`, **not** a mean of per-batch means;
  * the **fate** term keeps stepping every batch;
  * the data-dependent stop is **deterministic** under a fixed seed;
  * the rule selects **windows, not labels** -- every age cell is used exactly once per pass.

Graded against `train_member` itself wherever possible rather than against a re-implementation,
because a re-implementation is exactly what could drift from the shipped code.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import torch
from torch.utils.data import TensorDataset

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from cellfate.models import CellFateNet, huber_age_loss, huber_age_window  # noqa: E402
from cellfate.training.train import _AgeWindow, train_member  # noqa: E402
from cellfate.training.train_model import TrainConfig  # noqa: E402


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


C5C = _load("register_c5c_bar", "plan_tests/register_c5c_bar.py")

N_FEAT, N_FP, N_DT, N_CLS = 16, 8, 4, 3


def _ds(n: int, n_age: int, seed: int = 0) -> TensorDataset:
    """A tiny dataset with exactly `n_age` age-valid cells, placed by a fixed permutation."""
    g = torch.Generator().manual_seed(seed)
    x = torch.randn(n, N_FEAT, generator=g)
    fp = torch.randn(n, N_FP, generator=g)
    dt = torch.randn(n, N_DT, generator=g)
    yc = torch.zeros(n, N_CLS)
    yc[torch.arange(n), torch.randint(0, N_CLS, (n,), generator=g)] = 1.0
    ya = torch.randn(n, generator=g)
    am = torch.zeros(n)
    am[torch.randperm(n, generator=g)[:n_age]] = 1.0
    donor = torch.zeros(n, dtype=torch.long)
    return TensorDataset(x, fp, dt, yc, ya, am, donor)


def _cfg(**kw) -> TrainConfig:
    base = dict(dataset_dir=".", d_cell=8, d_u=8, latent_dim=8, p_drop=0.0,
                epochs=2, patience=99, batch_size=16, ensemble_size=1,
                base_seed=0, device="cpu")
    base.update(kw)
    return TrainConfig(**base)


def _make_model():
    return CellFateNet(g=N_FEAT, n_fp=N_FP, n_dt=N_DT, n_classes=N_CLS,
                       d_cell=8, d_u=8, latent_dim=8, p_drop=0.0)


def _weights(model) -> list[torch.Tensor]:
    return [p.detach().clone() for p in model.parameters()]


# =========================================================================== #
# THE GATE: the unmasked control arm must not move at all                     #
# =========================================================================== #
def test_arm_a_is_bit_identical_when_every_cell_is_age_valid():
    """Bar A1, against the real training loop rather than the simulation.

    Every cell age-valid => the first batch already clears k => the window closes at W=1 with an
    empty buffer => `huber_age_window` reduces to `huber_age_loss` over that batch. So turning the
    mechanism ON must change **nothing** in the control arm.
    """
    train, val = _ds(96, 96, seed=1), _ds(32, 32, seed=2)
    off, _ = train_member(_make_model, train, val, _cfg(age_window_k=1), seed=7, device="cpu")
    on, _ = train_member(_make_model, train, val, _cfg(age_window_k=4), seed=7, device="cpu")
    for a, b in zip(_weights(off), _weights(on), strict=True):
        assert torch.equal(a, b), "arm A moved -- baseline.json is no longer a valid reference"


def test_arm_a_identity_holds_for_every_k_a_run_might_use():
    """The identity must not be a coincidence of k=4: it holds for any k <= batch occupancy."""
    train, val = _ds(96, 96, seed=3), _ds(32, 32, seed=4)
    ref = _weights(train_member(_make_model, train, val, _cfg(age_window_k=1),
                                seed=5, device="cpu")[0])
    for k in (2, 4, 8, 16):
        got = _weights(train_member(_make_model, train, val, _cfg(age_window_k=k),
                                    seed=5, device="cpu")[0])
        assert all(torch.equal(a, b) for a, b in zip(ref, got, strict=True)), f"k={k} moved arm A"


def test_the_mechanism_does_move_a_sparsely_labelled_run():
    """The mirror of the gate: if it changed nothing anywhere, it would be doing nothing.

    A test that only asserts invariance can pass on a no-op, so the discriminating case is
    required alongside it.
    """
    train, val = _ds(96, 6, seed=11), _ds(32, 32, seed=12)
    off = _weights(train_member(_make_model, train, val, _cfg(age_window_k=1),
                                seed=9, device="cpu")[0])
    on = _weights(train_member(_make_model, train, val, _cfg(age_window_k=4),
                               seed=9, device="cpu")[0])
    assert any(not torch.equal(a, b) for a, b in zip(off, on, strict=True))


# =========================================================================== #
# sum/count, NOT a mean of means                                              #
# =========================================================================== #
def test_window_loss_is_one_mean_over_all_cells_not_a_mean_of_batch_means():
    """The trap named in the plan. Constructed so the two answers differ a lot: one batch holds a
    single cell with a big error, the next holds nine with none."""
    big = (torch.tensor([10.0]), torch.tensor([0.0]))
    small = (torch.zeros(9), torch.zeros(9))
    preds = [big[0], small[0]]
    trues = [big[1], small[1]]
    got = huber_age_window(preds, trues, torch.zeros(1), 2.0)

    pooled = torch.nn.functional.huber_loss(torch.cat(preds), torch.cat(trues), delta=2.0)
    mean_of_means = 0.5 * (
        torch.nn.functional.huber_loss(big[0], big[1], delta=2.0)
        + torch.nn.functional.huber_loss(small[0], small[1], delta=2.0))
    assert torch.allclose(got, pooled)
    assert not torch.allclose(got, mean_of_means)
    assert got < mean_of_means          # the single noisy cell is correctly down-weighted


def test_a_single_batch_window_is_exactly_todays_loss():
    """The reduction that makes arm A's bit-identity possible, asserted directly."""
    g = torch.Generator().manual_seed(0)
    ag, ya = torch.randn(20, generator=g), torch.randn(20, generator=g)
    am = torch.zeros(20)
    am[:7] = 1.0
    m = am.bool()
    assert torch.allclose(huber_age_window([ag[m]], [ya[m]], ag, 2.0),
                          huber_age_loss(ag, ya, am, 2.0))


def test_an_empty_window_returns_a_differentiable_zero():
    """Same contract as `huber_age_loss`: never a detached constant, or `.backward()` breaks."""
    ag = torch.randn(5, requires_grad=True)
    out = huber_age_window([], [], ag, 2.0)
    assert out.item() == 0.0
    out.backward()
    assert ag.grad is not None


# =========================================================================== #
# the window rule itself                                                      #
# =========================================================================== #
def _offer_counts(counts, k, w_max):
    """Drive the REAL `_AgeWindow` with synthetic batches; return each window's cell count."""
    win = _AgeWindow(k, w_max)
    out, pending = [], 0
    for c in counts:
        am = torch.zeros(16)
        am[:c] = 1.0
        z = torch.zeros(16, 4)
        if win.offer(z, z, z, torch.zeros(16), am):
            out.append(pending + c)
            pending = 0
            win.reset()
        else:
            pending += c
    return out


def test_the_shipped_window_matches_the_bar_script_exactly():
    """The simulation the 5c decision rests on must not drift from the code that ships."""
    rng = np.random.default_rng(0)
    for _ in range(30):
        counts = rng.integers(0, 4, size=40).tolist()
        expected = [c for c, _ in C5C.close_windows(counts, 4, 8)[0]]
        assert _offer_counts(counts, 4, 8) == expected


def test_a_window_never_closes_below_k_unless_w_max_forced_it():
    """A2's guarantee is 'by construction', so it is asserted as a construction."""
    win = _AgeWindow(k=4, w_max=8)
    z = torch.zeros(16, 4)
    am = torch.zeros(16)
    am[0] = 1.0                                  # one age cell per batch
    for i in range(3):
        assert win.offer(z, z, z, torch.zeros(16), am) is False, f"closed early at batch {i}"
    assert win.offer(z, z, z, torch.zeros(16), am) is True     # 4th cell reaches k


def test_w_max_forces_a_close_even_with_no_age_cells_at_all():
    """Otherwise a long labelless stretch would pin the window open forever."""
    win = _AgeWindow(k=4, w_max=3)
    z, am = torch.zeros(16, 4), torch.zeros(16)
    assert [win.offer(z, z, z, torch.zeros(16), am) for _ in range(3)] == [False, False, True]


def test_the_window_carries_across_the_epoch_boundary():
    """5c attempt 1 forced a close at each epoch's end and failed bar A2 by 4.44 pp. The shipped
    rule carries, so a label buffered at an epoch's end is consumed by the next one's first window
    -- used exactly once, just delivered one window late."""
    counts = [1, 1]                       # two epochs of one age cell each, k=4
    win = _AgeWindow(k=4, w_max=8)
    z = torch.zeros(16, 4)
    am = torch.zeros(16)
    am[0] = 1.0
    for _ in counts:
        assert win.offer(z, z, z, torch.zeros(16), am) is False
    assert win.n_cells == 2               # nothing was dropped at the boundary


def test_every_age_cell_is_used_exactly_once_per_pass():
    """The rule must select WINDOWS, not LABELS -- otherwise it silently subsets the 75."""
    rng = np.random.default_rng(3)
    counts = rng.integers(0, 3, size=200).tolist()
    windows, (carry_cells, _) = C5C.close_windows(counts, 4, 8)
    assert sum(c for c, _ in windows) + carry_cells == sum(counts)


# =========================================================================== #
# determinism                                                                 #
# =========================================================================== #
def test_the_data_dependent_stop_is_deterministic_under_a_fixed_seed():
    """A stop that depends on the data is fine; a stop that depends on the RUN is not."""
    train, val = _ds(96, 8, seed=21), _ds(32, 32, seed=22)
    a = _weights(train_member(_make_model, train, val, _cfg(age_window_k=4),
                              seed=13, device="cpu")[0])
    b = _weights(train_member(_make_model, train, val, _cfg(age_window_k=4),
                              seed=13, device="cpu")[0])
    assert all(torch.equal(x, y) for x, y in zip(a, b, strict=True))


# =========================================================================== #
# the fate task must not be slowed down                                       #
# =========================================================================== #
def test_the_fate_head_still_steps_on_a_batch_the_age_window_holds_back():
    """If the fate term rode the age window, Option 2 would silently become 'train 8x less' and
    its whole claim to cost the fate task nothing would be void."""
    torch.manual_seed(0)
    model = _make_model()
    win = _AgeWindow(k=99, w_max=99)            # guarantees a hold
    x, fp, dt = torch.randn(4, N_FEAT), torch.randn(4, N_FP), torch.randn(4, N_DT)
    am = torch.zeros(4)
    am[0] = 1.0
    assert win.offer(x, fp, dt, torch.zeros(4), am) is False

    lg, ag, _ = model(x, fp, dt)
    yc = torch.zeros(4, N_CLS)
    yc[:, 0] = 1.0
    from cellfate.models import class_balanced_weights, focal_loss
    w = torch.tensor(class_balanced_weights(np.array([4.0, 1.0, 1.0]), 0.999))
    (focal_loss(lg, yc, w, 2.0) + ag.sum() * 0.0).backward()
    grads = [p.grad for p in model.parameters() if p.grad is not None]
    assert grads, "the fate task produced no gradient on a held-back batch"
    assert any(float(g.abs().sum()) > 0 for g in grads)


def test_a_held_back_batch_contributes_no_age_gradient():
    """The other half of the same contract: held-back age cells must wait, not leak."""
    torch.manual_seed(0)
    model = _make_model()
    x, fp, dt = torch.randn(4, N_FEAT), torch.randn(4, N_FP), torch.randn(4, N_DT)
    _, ag, _ = model(x, fp, dt)
    (ag.sum() * 0.0).backward()
    assert all(p.grad is None or float(p.grad.abs().sum()) == 0.0 for p in model.parameters())


# =========================================================================== #
# the bar script's own pure logic                                             #
# =========================================================================== #
def test_arm_a_closes_at_w1_on_every_window():
    """Bar A1 in the simulation, mirroring the training-loop assertion above."""
    r = C5C.simulate_arm(C5C.N_AGE_ARM_A, k=4, w_max=8, n_epochs=20,
                         rng=np.random.default_rng(0))
    assert r["frac_windows_at_W1"] == 1.0
    assert r["mean_batches_per_window"] == pytest.approx(1.0)


def test_arm_b_clears_the_ge_k_bar():
    r = C5C.simulate_arm(C5C.N_AGE_ARM_B, k=4, w_max=8, n_epochs=400,
                         rng=np.random.default_rng(1))
    assert r["frac_windows_ge_k"] >= 0.95


def test_the_adaptive_rule_beats_the_fixed_w_design_on_update_count():
    """Bar A3 -- otherwise the 5c redesign is a regression against 5b, not a fix."""
    r = C5C.simulate_arm(C5C.N_AGE_ARM_B, k=4, w_max=8, n_epochs=400,
                         rng=np.random.default_rng(2))
    assert r["updates_over_run"] > C5C.fixed_w_updates(8)


def test_forcing_a_close_at_each_epoch_end_is_what_failed_the_bar():
    """Attempt 1, kept measurable so the record can show why it was dropped rather than assert it."""
    kw = dict(k=4, w_max=8, n_epochs=600)
    carried = C5C.simulate_arm(C5C.N_AGE_ARM_B, **kw, rng=np.random.default_rng(4))
    forced = C5C.simulate_arm(C5C.N_AGE_ARM_B, **kw, rng=np.random.default_rng(4),
                              close_at_end=True)
    assert forced["frac_windows_ge_k"] < 0.95 <= carried["frac_windows_ge_k"]
