"""
Test 3 (ΔAge lab notebook) — Is the model tying ridge because the signal is LINEAR
(nothing to fix) or because it CAN'T CAPTURE NONLINEARITY (a fixable model problem)?

Method: synthetic data with a KNOWN ΔAge signal in a gene-expression-like input, at
LARGE n (so sample size is not a confound — that's a later test). We train the REAL
CellFateNet regressor and compare its held-out ΔAge MAE to sklearn Ridge (alpha=1.0),
exactly the baseline the real pipeline uses. Three sub-runs, each ONCE:

  3a LINEAR    : ΔAge = linear combo of genes  -> model should TIE ridge (ridge optimal)
  3b NONLINEAR : ΔAge = gene products/thresholds -> model should BEAT ridge (ridge can't)
  3c NOISE     : ΔAge = pure noise             -> BOTH should fail (MAE ~ std)  [harness check]

Key metric = ridge_MAE - model_MAE (positive => model better). Positive controls:
3a ties, 3b model wins, 3c both fail. If 3b comes back a TIE, the architecture can't
capture nonlinearity -> a real fixable model problem.

GUARDRAIL: hyperparameters fixed a priori (Adam 1e-3, 60 epochs, batch 256). Run once,
do not tune against the result.
"""
from __future__ import annotations

import numpy as np
import torch
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.console import install_pretty_console, render_table
from cellfate.models.network import CellFateNet

SEED = 0
N = 20_000          # LARGE n (isolates architecture capability from sample size)
G = 800             # genes (input width)
K = 40              # genes carrying signal
EPOCHS = 60
BATCH = 256
LR = 1e-3
AGE_STD = 15.0      # target ΔAge spread in "years" (realistic scale)


def make_data(kind: str, seed: int):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((N, G)).astype(np.float32)          # expression-like input
    idx = rng.choice(G, K, replace=False)
    if kind == "linear":
        w = rng.standard_normal(K)
        y = X[:, idx] @ w
    elif kind == "smooth":
        # single-index SMOOTH nonlinearity: y = sin(1.5 z) + 0.4 z^2, z = X.w.
        # This is the kind of nonlinearity MLPs are GOOD at; ridge can't fit sin/square.
        z = X[:, idx] @ rng.standard_normal(K)
        z = (z - z.mean()) / (z.std() + 1e-9)
        y = np.sin(1.5 * z) + 0.4 * z ** 2
    elif kind == "products":
        # pairwise products + thresholds: mean-zero in each single gene, so a LINEAR
        # model sees ~no signal. MLPs are KNOWN to struggle with multiplication.
        a = np.zeros(N)
        for j in range(0, K - 1, 2):
            a += X[:, idx[j]] * X[:, idx[j + 1]]                # interactions
        for j in range(K):
            a += np.maximum(X[:, idx[j]] - 0.5, 0.0)            # thresholds (ReLU)
        y = a
    elif kind == "noise":
        y = rng.standard_normal(N)                             # no dependence on X
    else:
        raise ValueError(kind)
    y = (y - y.mean()) / (y.std() + 1e-9) * AGE_STD            # scale to AGE_STD
    return X, y.astype(np.float32)


def train_model(Xtr, ytr, Xte) -> tuple[np.ndarray, float]:
    """Train the REAL CellFateNet on (X -> ΔAge); return (test preds, TRAIN mae).
    Train MAE tells us if a tie is undertraining (train also high) vs genuine
    inability (train low, test high) vs can't-fit-at-all (both high)."""
    torch.manual_seed(SEED)
    net = CellFateNet(g=G, pert_kind="tf")
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    n_pert = net.arch["n_pert"]
    u_tr = torch.zeros(len(Xtr), n_pert)
    dt_tr = torch.zeros(len(Xtr), 2)
    Xtr_t, ytr_t = torch.tensor(Xtr), torch.tensor(ytr)
    net.train()
    for _ep in range(EPOCHS):
        perm = torch.randperm(len(Xtr_t))
        for i in range(0, len(Xtr_t), BATCH):
            b = perm[i:i + BATCH]
            opt.zero_grad()
            _, age, _ = net(Xtr_t[b], u_tr[b], dt_tr[b])
            loss = torch.mean((age - ytr_t[b]) ** 2)           # MSE on ΔAge only
            loss.backward()
            opt.step()
    net.eval()
    with torch.no_grad():
        _, age_tr, _ = net(Xtr_t, u_tr, dt_tr)
        train_mae = float(np.abs(age_tr.numpy() - ytr).mean())
        u_te = torch.zeros(len(Xte), n_pert)
        dt_te = torch.zeros(len(Xte), 2)
        _, age, _ = net(torch.tensor(Xte), u_te, dt_te)
    return age.numpy(), train_mae


def run(kind: str) -> list:
    X, y = make_data(kind, SEED)
    ntr = int(0.7 * N)
    Xtr, Xte = X[:ntr], X[ntr:]
    ytr, yte = y[:ntr], y[ntr:]

    # ridge (same as the real baseline: standardize + Ridge(alpha=1.0))
    sx = StandardScaler().fit(Xtr)
    ridge = Ridge(alpha=1.0).fit(sx.transform(Xtr), ytr)
    ridge_pred = ridge.predict(sx.transform(Xte))
    ridge_mae = float(np.abs(ridge_pred - yte).mean())

    model_pred, train_mae = train_model(Xtr, ytr, Xte)
    model_mae = float(np.abs(model_pred - yte).mean())

    baseline_mae = float(np.abs(yte - ytr.mean()).mean())      # predict-the-mean floor
    gap = ridge_mae - model_mae
    read = ("model WINS" if gap > 0.15 * ridge_mae
            else "TIE" if abs(gap) <= 0.15 * ridge_mae else "model LOSES")
    return [kind, f"{baseline_mae:.2f}", f"{ridge_mae:.2f}", f"{model_mae:.2f}",
            f"{train_mae:.2f}", f"{gap:+.2f}", read]


def main() -> None:
    install_pretty_console()
    print("\nTEST 3 — model vs ridge on KNOWN signal (large n, no noise)")
    print(f"(N={N:,}, {G} genes, {K} signal genes, ΔAge std={AGE_STD})")
    print("predict-mean = floor (no signal).  model-train = model's error on TRAIN "
          "(low+test-high => can't generalize; both high => can't fit at all)")
    rows = [run("linear"), run("smooth"), run("products"), run("noise")]
    print("\n" + render_table(
        ["signal", "predict-mean", "ridge MAE", "model MAE", "model-train", "gap", "read"],
        rows, aligns=["l", "r", "r", "r", "r", "r", "l"]))
    print("\n   EXPECTED (positive controls):")
    print("     linear   -> ridge optimal (~0); model close        (training works)")
    print("     smooth   -> model WINS  (MLPs handle smooth single-index nonlinearity)")
    print("     products -> model may TIE (MLPs are known to struggle with multiplication)")
    print("     noise    -> both ~= predict-mean                   (harness invents no signal)")
    print("\n   READ THE RESULT:")
    print("     - model WINS on SMOOTH -> architecture CAN capture (the right kind of)")
    print("       nonlinearity. Real-data tie then means real ΔAge is ~linear or product-like.")
    print("     - model TIES on SMOOTH too (esp. if model-train is also high) -> a real,")
    print("       fixable training/architecture problem: it can't fit nonlinearity it should.")


if __name__ == "__main__":
    main()
