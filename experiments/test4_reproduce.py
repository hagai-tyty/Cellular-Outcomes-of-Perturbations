"""
Test 4 (ΔAge lab notebook) — "train on X->Y, feed X back, does it reproduce Y?" —
done as a DIAGNOSTIC by varying how much of the signal the model can SEE.

Idea (from the user): train the real model on (X -> Y), then feed the same X back and
check it reproduces Y. We sharpen it: the signal always fully determines Y, but we hide
a growing fraction of the signal-carrying genes from the MODEL'S INPUT. This tests
whether imperfect reproduction is a BUG or just the known missing-genes effect (the real
panel sees only ~47% of the clock's genes).

  coverage 100% -> model sees all signal genes -> should reproduce Y almost perfectly
  coverage  75/47/25% -> signal genes hidden from input -> reproduction MUST degrade,
                         because the hidden genes carry signal the model can't access.

READ:
  - reproduces at 100%, degrades as coverage drops -> imperfect reproduction on real data
    is the MISSING-GENES effect, not a bug. Diagnosis confirmed.
  - can't reproduce even at 100% -> something is genuinely BROKEN (a real find).

Signal is SMOOTH (single-index) — the kind Test 3 proved the architecture CAN learn — so
any failure to reproduce is about visibility, not an unlearnable function.

GUARDRAIL: hyperparameters fixed a priori (Adam 1e-3, 60 epochs, batch 256). Run once.
"""
from __future__ import annotations

import numpy as np
import torch

from cellfate.common.console import install_pretty_console, render_table
from cellfate.models.network import CellFateNet

SEED = 0
N = 8_000
G = 800
K = 40              # signal-carrying genes
EPOCHS = 60
BATCH = 256
LR = 1e-3
AGE_STD = 15.0
COVERAGES = [1.00, 0.75, 0.47, 0.25]   # fraction of signal genes VISIBLE to the model


def make_full_signal(seed: int):
    """X (N,G) and Y determined by a SMOOTH function of K signal genes.
    Returns X, Y, and the indices of the signal genes (so we can hide some)."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((N, G)).astype(np.float32)
    idx = rng.choice(G, K, replace=False)
    z = X[:, idx] @ rng.standard_normal(K)
    z = (z - z.mean()) / (z.std() + 1e-9)
    y = np.sin(1.5 * z) + 0.4 * z ** 2                     # smooth single-index (MLP-friendly)
    y = (y - y.mean()) / (y.std() + 1e-9) * AGE_STD
    return X, y.astype(np.float32), idx


def reproduce_mae(Xin: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Train real CellFateNet on (Xin -> y), feed the SAME Xin back, return
    (reproduction MAE on training rows, predict-mean floor)."""
    torch.manual_seed(SEED)
    net = CellFateNet(g=Xin.shape[1], pert_kind="tf")
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    n_pert = net.arch["n_pert"]
    u = torch.zeros(len(Xin), n_pert)
    dt = torch.zeros(len(Xin), 2)
    Xt, yt = torch.tensor(Xin), torch.tensor(y)
    net.train()
    for _ in range(EPOCHS):
        perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), BATCH):
            b = perm[i:i + BATCH]
            opt.zero_grad()
            _, age, _ = net(Xt[b], u[b], dt[b])
            torch.mean((age - yt[b]) ** 2).backward()
            opt.step()
    net.eval()
    with torch.no_grad():
        _, age, _ = net(Xt, u, dt)                          # feed the SAME X back
    repro = float(np.abs(age.numpy() - y).mean())
    floor = float(np.abs(y - y.mean()).mean())
    return repro, floor


def main() -> None:
    install_pretty_console()
    X, y, idx = make_full_signal(SEED)
    rng = np.random.default_rng(SEED + 1)

    print("\nTEST 4 — feed training X back, does the model reproduce Y? (vs gene visibility)")
    print(f"(N={N:,}, {G} genes, {K} signal genes, SMOOTH signal, ΔAge std={AGE_STD})")
    print("we hide signal genes from the INPUT only; Y always fully determined by them.")

    rows = []
    for cov in COVERAGES:
        keep = int(round(K * cov))
        vis = set(rng.choice(idx, keep, replace=False).tolist()) if keep else set()
        # zero out the hidden signal genes in the model's input (keep dims fixed)
        Xin = X.copy()
        hidden = [g for g in idx if g not in vis]
        Xin[:, hidden] = 0.0
        repro, floor = reproduce_mae(Xin, y)
        read = ("reproduces" if repro < 0.2 * floor
                else "degraded" if repro < 0.8 * floor else "can't (~floor)")
        rows.append([f"{int(cov*100)}%", f"{keep}/{K}", f"{repro:.2f}", f"{floor:.2f}", read])

    print("\n" + render_table(
        ["gene coverage", "signal genes seen", "reproduce MAE", "predict-mean", "read"],
        rows, aligns=["r", "r", "r", "r", "l"]))
    print("\n   EXPECTED / READ:")
    print("     100% -> 'reproduces' (MAE near 0): the model CAN reproduce Y it was trained")
    print("             on when it sees the signal. Training path works — NOT broken.")
    print("     lower coverage -> reproduction degrades toward the predict-mean floor,")
    print("             because hidden genes carry signal the model can't access.")
    print("     -> This is exactly the real situation (panel sees ~47% of clock genes):")
    print("        imperfect reproduction on real data is MISSING GENES, not a bug.")
    print("     If 100% did NOT reproduce (MAE ~ floor) -> something IS broken (a real find).")


if __name__ == "__main__":
    main()
