"""Is the ΔAge regression CIRCULAR?  (read-only; ridge fit + a dot product)

    python experiments/diag_clock_circularity.py

THE QUESTION
------------
`diag_target_path` showed ridge and `x_only` are identical to two decimals on every fold, so the
perturbation contributes nothing and the "prediction" is a pure function of the cell's own
expression. But ΔAge is DEFINED as `clock(cell) - clock(control)`, and the Fleischer clock is
LINEAR in log-normalised expression. So a linear model on expression may simply be re-deriving a
linear functional of its own input -- which would make the good ridge MAE a statement about the
gene panel, not about the model, the biology, or any prediction.

It is NOT circular by construction, and that is why this has to be measured:
  * the label is computed from the **full normalised profile** (~33k genes), not the model input
    (`build_dataset.py:342-346`);
  * the model sees a 2,000-gene panel that retains only ~21% of the clock's |w| mass;
  * `X` is additionally HARMONIZED, and ΔAge is deconfounded and re-anchored afterwards.
Any of those could break the correspondence. The question is whether they do.

WHAT IS MEASURED, on the HELD-OUT donor only
  clock_panel(X) = X @ w   with w[i] = clock.weights.get(panel_gene[i], 0.0)
  -- the clock's own weights applied to the model's own input. The intercept and the ΔAge control
  offset are both CONSTANTS within a single held-out donor, so they shift the values without
  touching any correlation; they are omitted rather than guessed at.

  T1  rho(clock_panel, y_age)      Is the LABEL recoverable from the model's input by the clock?
  T2  rho(ridge_pred, clock_panel) Is RIDGE that readout?
  T3  rho(ridge_pred, y_age)       Reference: the plain predictive correlation.

PRE-REGISTERED READING (constant below, fixed before running)
  CIRCULAR              T1 >= RHO_CIRCULAR and T2 >= RHO_CIRCULAR
                        -> the target is a linear readout of the input and ridge re-derives it.
                           The ΔAge MAE numbers describe the panel, not a prediction.
  LABEL-RECOVERABLE     T1 >= RHO_CIRCULAR, T2 < RHO_CIRCULAR
                        -> the task is trivially solvable but ridge is solving it another way.
  NOT CIRCULAR          T1 < RHO_CIRCULAR
                        -> the panel does NOT reconstruct the clock; ridge's skill is not a
                           tautology and needs a different explanation.

NOT CLAIMED: that a high T1 makes ΔAge biologically meaningless. It would mean the REGRESSION is
uninformative -- predicting a quantity already determined by the input -- which is a statement
about this experiment, not about the clock or about aging.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))

from cellfate.common.io import ArtifactPaths  # noqa: E402
from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.evaluation.baselines import RidgeLinear, XOnly  # noqa: E402
from cellfate.evaluation.data import gather_split  # noqa: E402

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
ARMS = {"pre-C-7": "_armA", "C-7": "_c7t"}
REGIME = "holdout"
CLOCK_PATH = ROOT / "configs" / "clocks" / "fleischer_clock.json"

RHO_CIRCULAR = 0.95


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def clock_weights_for_panel(panel_genes: list[str], weights: dict) -> tuple[np.ndarray, dict]:
    """The clock's own weights, aligned to the panel's gene order. Genes the clock never saw get
    0.0 -- the same rule `LinearClock.predict_age` uses (`aging.py:55`), so this is the clock's
    behaviour on a reduced gene set, not a re-derivation of it."""
    w = np.array([float(weights.get(g, 0.0)) for g in panel_genes], dtype=np.float64)
    total = sum(abs(v) for v in weights.values())
    return w, {"n_panel": len(panel_genes),
               "n_with_weight": int((w != 0).sum()),
               "abs_mass_retained": float(np.abs(w).sum() / total) if total else float("nan")}


def verdict(t1: float, t2: float) -> str:
    if not np.isfinite(t1):
        return "UNDETERMINED"
    if t1 < RHO_CIRCULAR:
        return "NOT CIRCULAR"
    if np.isfinite(t2) and t2 >= RHO_CIRCULAR:
        return "CIRCULAR"
    return "LABEL-RECOVERABLE"


def audit(root: Path, w: np.ndarray) -> dict:
    paths = ArtifactPaths.of(str(root))
    tr, te = gather_split(paths, REGIME, "train"), gather_split(paths, REGIME, "test")
    m = te.mask.astype(bool)
    if te.n == 0 or m.sum() < 3:
        return {"n": int(m.sum()) if te.n else 0}
    ridge = RidgeLinear().fit(tr)
    xonly = XOnly().fit(tr)
    _, r_age = ridge.predict(te.X, te.fp, te.dose_time)
    _, x_age = xonly.predict(te.X, te.fp, te.dose_time)
    cp = np.asarray(te.X, dtype=np.float64) @ w          # clock weights on the MODEL's input
    y = np.asarray(te.y_age, dtype=np.float64)
    cp, y, r_age, x_age = cp[m], y[m], np.asarray(r_age)[m], np.asarray(x_age)[m]
    t1, t2 = pearson(cp, y), pearson(r_age, cp)
    return {"n": int(m.sum()),
            "T1_clockpanel_vs_label": t1,
            "T2_ridge_vs_clockpanel": t2,
            "T3_ridge_vs_label": pearson(r_age, y),
            "xonly_vs_clockpanel": pearson(x_age, cp),
            "sd_clockpanel": float(np.std(cp, ddof=1)), "sd_label": float(np.std(y, ddof=1)),
            "verdict": verdict(t1, t2)}


def main() -> None:
    clock = json.loads(CLOCK_PATH.read_text(encoding="utf-8"))
    res: dict = {"clock_meta": clock.get("meta", {}), "arms": {}}
    print("=" * 104)
    print("CLOCK-CIRCULARITY TEST  --  is the ΔAge label a linear readout of the model's own input?")
    print(f"  pre-registered: CIRCULAR iff T1 >= {RHO_CIRCULAR} AND T2 >= {RHO_CIRCULAR}")
    print(f"  clock: {clock['meta']['source']}  ({len(clock['weights'])} genes)")
    print("=" * 104)
    for arm, sfx in ARMS.items():
        print(f"\n[{arm}]")
        print(f"  {'fold':<6}{'n':>4}{'T1 clk~lbl':>12}{'T2 ridge~clk':>14}{'T3 ridge~lbl':>14}"
              f"{'xonly~clk':>11}{'sd_clk':>9}{'sd_lbl':>9}   verdict")
        res["arms"][arm] = {}
        for d in DONORS:
            root = ROOT / f"cellfate_loocv_{d}{sfx}"
            if not (root / "shards").is_dir():
                continue
            panel = GenePanel.load(str(root / "panel.json"))
            w, cov = clock_weights_for_panel(list(panel.genes), clock["weights"])
            res.setdefault("coverage", cov)
            try:
                r = audit(root, w)
            except Exception as exc:                                     # noqa: BLE001
                print(f"  {d:<6}  FAILED: {type(exc).__name__}: {exc}")
                continue
            res["arms"][arm][d] = r
            if r["n"] < 3:
                print(f"  {d:<6}{r['n']:>4}   (no age-valid test cells -- C-7 masks N2)")
                continue
            print(f"  {d:<6}{r['n']:>4}{r['T1_clockpanel_vs_label']:>12.4f}"
                  f"{r['T2_ridge_vs_clockpanel']:>14.4f}{r['T3_ridge_vs_label']:>14.4f}"
                  f"{r['xonly_vs_clockpanel']:>11.4f}{r['sd_clockpanel']:>9.2f}"
                  f"{r['sd_label']:>9.2f}   {r['verdict']}")

    c = res.get("coverage", {})
    print(f"\n  panel coverage of the clock: {c.get('n_with_weight')} of {c.get('n_panel')} genes "
          f"carry weights, holding {c.get('abs_mass_retained', float('nan')):.4f} of total |w| mass")
    for arm in ARMS:
        vs = [v["verdict"] for v in res["arms"].get(arm, {}).values() if v.get("n", 0) >= 3]
        if vs:
            print(f"  {arm:<9} verdicts: {vs}")

    _RESULTS.mkdir(exist_ok=True)
    out = _RESULTS / "diag_clock_circularity_results.json"
    out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
