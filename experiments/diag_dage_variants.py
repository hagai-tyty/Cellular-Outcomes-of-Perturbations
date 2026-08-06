"""STAGE 1.5.6 — which ΔAge definition predicts an INDEPENDENT outcome? (Gill Sendai, GSE165176)

    python experiments/diag_dage_variants.py "D:\\Gill"

READ-ONLY. Writes `results/diag_dage_variants_results.json`. `src/` untouched, no labels move.

WHY THIS EXISTS
---------------
Every measurement in this arc so far correlates one RNA-derived quantity against another, which is
why it keeps ending in "consistent with either". This does not: it scores each candidate ΔAge
definition on its ability to predict the **FACS sort marker** -- CD13 (failed to reprogram, 47) vs
SSEA4 (reprogramming, 71). That is a **protein-level flow-cytometry outcome**. The clock reads RNA
and cannot have leaked into it.

THE TEST IS INCREMENTAL, AND THAT IS THE WHOLE POINT
----------------------------------------------------
SSEA4 is a pluripotency surface marker, so pluripotency predicts it well on its own. The question is
never "does ΔAge predict outcome" -- it is:

    does ΔAge add predictive power OVER pluripotency and day?

  adds nothing  -> the label is redundant with identity; that is what justifies replacing it
  adds power    -> it carries information identity does not, on an outcome outside the RNA

Baseline model  : outcome ~ day + pluripotency
Test model      : outcome ~ day + pluripotency + dAge_variant
Score           : leave-one-DONOR-out AUC increment (6 donors). Log-loss increment reported beside it.

Grouped LODO, not random CV: samples repeat within donor (exp1/exp2, many days), so a random split
would put the same donor on both sides and score memorisation.

VARIANTS SWEPT
--------------
  raw          the shipped Fleischer clock, unchanged -- the incumbent
  covnorm      renormalise by the COVERED weight mass. Gill carries 89% of the clock's |weight|, so
               ~11% of it silently reads as zero; this rescales instead of under-reading
  top100/500/2000   only the largest-|weight| genes. A dense ridge over 33,155 genes from 133
               samples is the regime that produced the clock's problems; sparsity is the obvious
               counter-hypothesis
  ranknorm     rank-transform each sample's expression before applying weights -- immune to library
               size and to the log/linear scale choice
  resid_pluri  dAge with pluripotency regressed out
  resid_cc     with the cell-cycle signature regressed out
  resid_both   both

**No variant is selected on the answer**: every one is scored and every one is reported, ranked.

WHAT THIS CANNOT DO
-------------------
It cannot prove a variant measures AGE. It measures which definition carries the most information
about a real reprogramming outcome that the RNA did not generate. A variant that wins here is the
best candidate; confirming it needs the methylation agreement check (step 2) and a second dataset.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

ANNOT_COLS = 12
SEED = 0


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# Pure logic — unit-tested with no repo data present                           #
# --------------------------------------------------------------------------- #
def rank_normalise(X: np.ndarray) -> np.ndarray:
    """Per-sample rank transform, scaled to [0,1]. Kills library size and scale choice."""
    out = np.empty_like(X, dtype=np.float64)
    for i in range(X.shape[0]):
        order = np.argsort(X[i], kind="mergesort")
        r = np.empty(X.shape[1], float)
        r[order] = np.arange(X.shape[1], dtype=float)
        out[i] = r / max(X.shape[1] - 1, 1)
    return out


def covered_weight_fraction(genes: list[str], weights: dict[str, float]) -> float:
    """Fraction of the clock's total |weight| that this gene space actually carries."""
    tot = sum(abs(v) for v in weights.values())
    have = sum(abs(weights[g]) for g in genes if g in weights)
    return float(have / tot) if tot else 0.0


def residualise(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """y with every column of Z (plus an intercept) least-squares removed."""
    Z = np.column_stack([np.ones(len(y)), np.asarray(Z, float).reshape(len(y), -1)])
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ beta


def control_relative(age: np.ndarray, donor: np.ndarray, is_ctrl: np.ndarray) -> np.ndarray:
    """ΔAge = age − that donor's own control mean. The project's definition, held constant across
    every variant so the sweep tests the CLOCK and not the zero-point."""
    d = np.empty_like(age, dtype=np.float64)
    for k in np.unique(donor):
        m = donor == k
        ref = age[m & is_ctrl]
        d[m] = age[m] - (ref.mean() if ref.size else age[m].mean())
    return d


def auc(y: np.ndarray, s: np.ndarray) -> float:
    """Rank-based AUC. Ties get 0.5 credit."""
    y = np.asarray(y).astype(bool)
    if y.all() or (~y).all():
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    r = np.empty(len(s), float)
    r[order] = np.arange(1, len(s) + 1)
    for v in np.unique(s):
        m = s == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    n1, n0 = int(y.sum()), int((~y).sum())
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def fit_logistic(X: np.ndarray, y: np.ndarray, iters: int = 200, l2: float = 1.0) -> np.ndarray:
    """Newton-IRLS logistic with an L2 penalty. Small n and few columns, so this is exact enough
    and avoids a sklearn dependency inside a diagnostic."""
    X = np.column_stack([np.ones(len(y)), X])
    w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))
        W = np.maximum(p * (1 - p), 1e-6)
        pen = l2 * np.eye(X.shape[1])
        pen[0, 0] = 0.0
        H = X.T @ (X * W[:, None]) + pen
        g = X.T @ (y - p) - pen @ w
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        w += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return w


def predict(w: np.ndarray, X: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones(X.shape[0]), X])
    return 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))


def lodo_auc(X: np.ndarray, y: np.ndarray, donor: np.ndarray) -> tuple[float, float]:
    """Out-of-fold AUC and mean log-loss, leaving one DONOR out at a time."""
    pred = np.full(len(y), np.nan)
    for k in np.unique(donor):
        te = donor == k
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        w = fit_logistic((X[tr] - mu) / sd, y[tr].astype(float))
        pred[te] = predict(w, (X[te] - mu) / sd)
    ok = np.isfinite(pred)
    p = np.clip(pred[ok], 1e-6, 1 - 1e-6)
    ll = float(-np.mean(y[ok] * np.log(p) + (1 - y[ok]) * np.log(1 - p)))
    return auc(y[ok], pred[ok]), ll


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def load_gill(gill_dir: Path):
    import pandas as pd
    expr = next(gill_dir.glob("*Log2_RPM*.txt.gz"))
    df = pd.read_csv(expr, sep="\t", low_memory=False)
    cols = list(df.columns)
    m = df.set_index(cols[0])[cols[ANNOT_COLS:]]
    lin = np.power(2.0, m.to_numpy(dtype=np.float64)) - 1.0
    lin[lin < 0] = 0.0
    sym = [str(s) for s in m.index]
    order = np.argsort(-lin.sum(axis=1))
    seen, keep = set(), []
    for i in order:
        if sym[i] not in seen:
            seen.add(sym[i])
            keep.append(i)
    keep.sort()
    return list(m.columns), [sym[i] for i in keep], lin[keep, :].T


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    from cellfate.data.normalize import normalize_counts
    dcv = _load("dcv", ROOT / "experiments" / "diag_clock_validity.py")

    gill = Path(sys.argv[1])
    samples, genes, lin = load_gill(gill)
    norm = normalize_counts(lin)                       # log1p CP10k — the pipeline's own transform
    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}

    donor = np.array([s.split("_")[0] for s in samples])
    mark = np.array(["SSEA4" if "SSEA4" in s else ("CD13" if "CD13" in s else "") for s in samples])
    day = np.array([float(m.group(1)) if (m := re.search(r"_d(\d+)_", s)) else 0.0 for s in samples])
    is_ctrl = day == 0.0

    gi = {g: i for i, g in enumerate(genes)}
    plu = norm[:, [gi[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi]].mean(axis=1)
    cc = norm[:, [gi[g] for g in dcv.CELL_CYCLE if g in gi]].mean(axis=1)

    cov = covered_weight_fraction(genes, W)
    print(f"\n[shape before statistic] {len(samples)} samples, {len(genes)} genes, "
          f"{len(np.unique(donor))} donors")
    print(f"  clock |weight| covered by this gene space: {cov:.3f}")
    print(f"  outcome: SSEA4 {int((mark == 'SSEA4').sum())} vs CD13 {int((mark == 'CD13').sum())}"
          f"  ({int((mark == '').sum())} unmarked day-0 baselines excluded from scoring)\n")

    wv = np.array([W.get(g, 0.0) for g in genes])
    absw = np.abs(wv)
    ranked = np.argsort(-absw)

    def age_from(mat: np.ndarray, wvec: np.ndarray, scale: float = 1.0) -> np.ndarray:
        return (mat @ wvec) * scale + float(clock.get("intercept", 0.0))

    variants: dict[str, np.ndarray] = {}
    variants["raw"] = age_from(norm, wv)
    variants["covnorm"] = age_from(norm, wv, scale=1.0 / cov if cov else 1.0)
    for k in (100, 500, 2000):
        sub = np.zeros_like(wv)
        sub[ranked[:k]] = wv[ranked[:k]]
        variants[f"top{k}"] = age_from(norm, sub)
    variants["ranknorm"] = age_from(rank_normalise(norm), wv)

    base_d = control_relative(variants["raw"], donor, is_ctrl)
    variants_d = {k: control_relative(v, donor, is_ctrl) for k, v in variants.items()}
    variants_d["resid_pluri"] = residualise(base_d, plu)
    variants_d["resid_cc"] = residualise(base_d, cc)
    variants_d["resid_both"] = residualise(base_d, np.column_stack([plu, cc]))

    scored = mark != ""
    y = (mark[scored] == "SSEA4").astype(float)
    Xb = np.column_stack([day[scored], plu[scored]])
    a0, l0 = lodo_auc(Xb, y, donor[scored])
    print(f"  BASELINE  outcome ~ day + pluripotency : AUC {a0:.4f}   logloss {l0:.4f}\n")
    print(f"  {'variant':>12} {'AUC':>8} {'dAUC':>8} {'logloss':>9} {'dlogloss':>9}")

    rows = {}
    for name, dv in variants_d.items():
        X = np.column_stack([Xb, dv[scored]])
        a, ll = lodo_auc(X, y, donor[scored])
        rows[name] = {"auc": a, "d_auc": a - a0, "logloss": ll, "d_logloss": ll - l0}
        print(f"  {name:>12} {a:8.4f} {a - a0:+8.4f} {ll:9.4f} {ll - l0:+9.4f}")

    rank = sorted(rows, key=lambda k: -rows[k]["d_auc"])
    print(f"\n  RANKED by AUC increment: {', '.join(rank)}")
    print(f"  best: {rank[0]}  (dAUC {rows[rank[0]]['d_auc']:+.4f})")
    out = {"script": "diag_dage_variants", "utc": datetime.now(UTC).isoformat(),
           "n_samples": len(samples), "n_scored": int(scored.sum()),
           "covered_weight_fraction": cov,
           "baseline": {"auc": a0, "logloss": l0, "model": "day + pluripotency"},
           "variants": rows, "ranked_by_d_auc": rank}
    (_RESULTS / "diag_dage_variants_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print("  wrote results/diag_dage_variants_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
