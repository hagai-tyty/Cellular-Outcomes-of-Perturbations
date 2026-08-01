"""STAGE 1.5.2 §6 — freeze the bars at the ACTUAL geometry, BEFORE any statistic is computed.

    python experiments/stage_1_5_2_resolvability.py

Writes `stage_1_5_2_resolvability_results.json`. Reads only sample TITLES (geometry), never
expression or beta values, so no bar can be tuned to the answer. `src/` untouched.

WHY THIS RUNS FIRST. §6 fixes the procedure: simulate each metric at the geometry it will be graded
on under a system that meets the intent EXACTLY, read `bar_verdict`, and freeze or move the bar
*before* the data is opened. §6's registered geometry was GSE165178's (n=11 per arm, n=22 total,
4 donor folds). The set actually being used is GSE165177 x GSE165179, whose geometry is different in
both directions — larger overall, but FEWER donors — so the bars must be re-checked rather than
assumed to carry over.

Registered bars (§6), restated:
  M-2a rho_within   n=11 per arm,  intent rho_true 0.70,  bar rho >= 0.50
  M-2a rho_partial  n=22, 1 covar, intent rho_true 0.70,  bar rho >= 0.50
  M-2b sign agree   up to 11 pairs, Bernoulli(0.85),      bar >= 8/11
  M-2c LODO MAE     4 donor folds, residual sd 5 yr,      bar <= 8.0 yr   [gated on M-2a]
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from audit_metrics import bar_verdict  # noqa: E402

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


N_SIM = 20000
RNG = np.random.default_rng(0)


def titles(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                return [x.strip('"') for x in line.rstrip().split("\t")[1:]]
    return []


def geometry(rna_series: Path, met_series: Path) -> dict:
    """Paired (donor, arm, day) conditions — TITLES ONLY, no measurements."""
    pairs = sorted(set(titles(rna_series)) & set(titles(met_series)))
    rows = []
    for t in pairs:
        m = re.match(r"^(O\d)_(.+?)_(\d+)days_(exp\d)$", t)
        if m:
            rows.append({"donor": m.group(1), "arm": m.group(2), "day": int(m.group(3))})
    cond = {(r["donor"], r["arm"], r["day"]) for r in rows}
    return {"n_joined": len(pairs), "n_conditions": len(cond),
            "donors": sorted({r["donor"] for r in rows}),
            "days": sorted({r["day"] for r in rows}),
            "per_arm_conditions": dict(Counter(a for _, a, _ in cond))}


def sim_spearman(n: int, rho_true: float, n_sim: int = N_SIM) -> np.ndarray:
    """Sampling distribution of Spearman rho for a system whose true rho IS rho_true."""
    from scipy.stats import spearmanr
    cov = np.array([[1.0, rho_true], [rho_true, 1.0]])
    out = np.empty(n_sim)
    for i in range(n_sim):
        xy = RNG.multivariate_normal([0, 0], cov, size=n)
        out[i] = spearmanr(xy[:, 0], xy[:, 1]).correlation
    return out


def sim_partial(n: int, rho_true: float, n_sim: int = 4000) -> np.ndarray:
    """Partial correlation after removing one covariate, at the same true rho."""
    out = np.empty(n_sim)
    cov = np.array([[1.0, rho_true], [rho_true, 1.0]])
    for i in range(n_sim):
        xy = RNG.multivariate_normal([0, 0], cov, size=n)
        z = RNG.normal(size=n)
        rx = xy[:, 0] - np.polyval(np.polyfit(z, xy[:, 0], 1), z)
        ry = xy[:, 1] - np.polyval(np.polyfit(z, xy[:, 1], 1), z)
        out[i] = np.corrcoef(rx, ry)[0, 1]
    return out


def sim_sign(n_pairs: int, p: float = 0.85, n_sim: int = N_SIM) -> np.ndarray:
    return RNG.binomial(n_pairs, p, size=n_sim).astype(float)


def sim_lodo_mae(n_folds: int, resid_sd: float = 5.0, per_fold: int = 20,
                 n_sim: int = N_SIM) -> np.ndarray:
    out = np.empty(n_sim)
    for i in range(n_sim):
        e = RNG.normal(0.0, resid_sd, size=n_folds * per_fold)
        out[i] = np.mean(np.abs(e))
    return out


def main() -> int:
    g = geometry(Path(r"D:\GSE165177\GSE165177_series_matrix.txt.gz"),
                 Path(r"D:\GSE165179\GSE165179_series_matrix.txt.gz"))
    print("STAGE 1.5.2 §6 — bar resolvability at the ACTUAL geometry (titles only)\n")
    print(f"  joined pairs {g['n_joined']} -> {g['n_conditions']} (donor, arm, day) conditions")
    print(f"  donors {g['donors']}  days {g['days']}")
    for a, c in sorted(g["per_arm_conditions"].items(), key=lambda x: -x[1]):
        print(f"     {a:48s} {c}")

    n_cond = g["n_conditions"]
    n_arm_min = min(g["per_arm_conditions"].values())
    n_donors = len(g["donors"])
    out: dict = {"geometry": g, "checks": {}}

    print("\n  BAR CHECKS  (RESOLVABLE = a system meeting the intent exactly passes >= 95%)\n")
    print(f"  {'bar':<44}{'geometry':>16}{'pass rate':>11}  verdict")
    print("  " + "-" * 88)

    def report(name: str, sim, bar, lower, geom_note):
        v = bar_verdict(np.asarray(sim, float), bar, lower_is_better=lower)
        out["checks"][name] = {**v, "bar": bar, "geometry": geom_note}
        flag = "" if v["verdict"] == "RESOLVABLE" else f"   -> usable_bar {v['usable_bar']:.3f}"
        print(f"  {name:<44}{geom_note:>16}{v['pass_rate']*100:>10.1f}%  {v['verdict']}{flag}")
        return v

    # M-2a rho_within: registered at n=11/arm; check BOTH the registered and the actual smallest arm
    report("M-2a rho_within (registered n=11/arm)",
           sim_spearman(11, 0.70), 0.50, False, "n=11")
    report("M-2a rho_within (ACTUAL smallest arm)",
           sim_spearman(n_arm_min, 0.70), 0.50, False, f"n={n_arm_min}")
    # M-2a rho_partial: registered n=22; actual n = all conditions
    report("M-2a rho_partial (registered n=22)",
           sim_partial(22, 0.70), 0.50, False, "n=22")
    report("M-2a rho_partial (ACTUAL)",
           sim_partial(n_cond, 0.70), 0.50, False, f"n={n_cond}")
    # M-2b sign agreement
    report("M-2b sign agreement (registered 8/11)",
           sim_sign(11), 8.0, False, "11 pairs")
    # M-2c LODO: registered 4 folds; ACTUAL is 3 donors
    report("M-2c LODO MAE (registered 4 folds)",
           sim_lodo_mae(4), 8.0, True, "4 folds")
    report("M-2c LODO MAE (ACTUAL donors)",
           sim_lodo_mae(n_donors), 8.0, True, f"{n_donors} folds")

    (_RESULTS / "stage_1_5_2_resolvability_results.json").write_text(
        json.dumps({"script": "stage_1_5_2_resolvability",
                    "utc": datetime.now(UTC).isoformat(timespec="seconds"), **out},
                   indent=2, default=str), encoding="utf-8")
    print("\n  wrote stage_1_5_2_resolvability_results.json")
    print("  NOTE: no expression or beta value was read. Bars can be frozen without seeing the answer.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
