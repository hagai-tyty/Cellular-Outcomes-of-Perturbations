"""PHASE 3 -- the forward question on GSE165177, where donor age is constant by construction.

    python experiments/diag_phase3_within_donor_forward.py

Pre-registered in `plans/THREE_TESTS_PREREG.md` Phase 3.

WHY THIS DESIGN
---------------
Every earlier forward attempt used the 6 Sendai donors, where donor chronological age explains the
outcome (r 0.931) and cannot be separated at n=6. GSE165177 has a design that was never used for
this question: each (donor, arm) is its own trajectory over days 10/13/15/17, and WITHIN a donor
the arms differ in OUTCOME while donor age is constant by construction. Up to 3 donors x 6 arms =
18 trajectories with the confound held fixed.

Clock: top100, because these are PERTURBED samples (Phase 0 Part C licensed dense only for resting
ones). Truncated weights have no valid intercept, so scores are not years; every statistic used
here is invariant to linear rescaling.

THE PRECONDITION THAT CAN KILL THE PHASE
-----------------------------------------
The window is d10 -> d17: SEVEN DAYS. If the score does not move materially across it there is
nothing to predict, and that is the result -- not something to work around. Checked and reported
FIRST, before any model is fitted.

PRE-REGISTERED READING
  Unit is the (donor, arm) trajectory; inference clusters on DONOR (3), not trajectory (18),
  because arms within a donor share material.
  SIGNAL  leave-one-DONOR-out prediction of the late score from early expression beats a
          permutation null at the 95th percentile for a majority of the alpha grid.
  NULL    otherwise.
With 3 donor clusters this is severely underpowered. Stated before the run: a NULL bounds what the
design can show, it does not establish absence.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dre = _load("dre", "experiments/diag_residual_expression.py")
p1 = _load("p1", "experiments/diag_phase1_top100.py")

from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.normalize import normalize_counts  # noqa: E402

DDIR = Path(r"D:\GSE165177")
MATRICES = ["GSE165177_Log2_RPM_Transient_reprogramming.txt.gz",
            "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz"]
CLOCK = ROOT / "configs" / "clocks" / "fleischer_clock.json"
EARLY_DAY, LATE_DAY = 10.0, 17.0
N_ANNOT = 12
MIN_MOVE = 1.0      # score units; below this the trajectory is flat and the phase ends


def load_all() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Both matrices, merged on gene symbol, CP10k+log1p, with donor/arm/day parsed."""
    frames = []
    for m in MATRICES:
        df = pd.read_csv(DDIR / m, sep="\t", low_memory=False)
        cols = list(df.columns[N_ANNOT:])
        lin = np.power(2.0, df[cols].to_numpy(dtype=np.float64)) - 1.0
        lin[lin < 0] = 0.0
        f = pd.DataFrame(lin, columns=cols)
        f["__sym__"] = df[df.columns[0]].astype(str).to_numpy()
        frames.append(f)
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.merge(f, on="__sym__", how="inner")
    sample_cols = [c for c in merged.columns if c != "__sym__"]
    merged["__tot__"] = merged[sample_cols].sum(1)
    merged = merged.sort_values("__tot__", ascending=False).drop_duplicates("__sym__", keep="first")
    genes = merged.__sym__.tolist()
    expr = normalize_counts(merged[sample_cols].to_numpy(dtype=np.float64).T, target_sum=1e4)
    rows = []
    for j, c in enumerate(sample_cols):
        mm = re.match(r"^([A-Z]\d+)_(.+?)_(\d+)days_(exp\d+)$", c)
        if not mm:
            continue
        rows.append({"j": j, "donor": mm.group(1), "arm": mm.group(2),
                     "day": float(mm.group(3)), "exp": mm.group(4), "col": c})
    return pd.DataFrame(rows), expr, genes


def main() -> None:
    meta, expr, genes = load_all()
    clock = json.loads(CLOCK.read_text(encoding="utf-8"))
    score = expr @ p1.top_n_weights(clock["weights"], genes, p1.TOP_N)
    meta["score"] = score[meta.j.to_numpy()]

    print("=" * 100)
    print("PHASE 3 -- within-donor forward test, GSE165177 (donor age constant by construction)")
    print("=" * 100)
    print(f"  samples parsed {len(meta)}  donors {sorted(meta.donor.unique())}  "
          f"arms {meta.arm.nunique()}  days {sorted(meta.day.unique())}")

    # ---- PRECONDITION: does the score move across d10 -> d17 within a trajectory? ----------- #
    print(f"\n[PRECONDITION] per (donor, arm) movement d{EARLY_DAY:.0f} -> d{LATE_DAY:.0f}")
    print(f"  {'donor':<7}{'arm':<48}{'early':>9}{'late':>9}{'move':>9}")
    moves, traj = [], []
    for (d, a), g in meta.groupby(["donor", "arm"]):
        e = g[g.day == EARLY_DAY]
        late = g[g.day == LATE_DAY]
        if e.empty or late.empty:
            continue
        ev, lv = float(e.score.mean()), float(late.score.mean())
        moves.append(abs(lv - ev))
        traj.append({"donor": d, "arm": a, "early": ev, "late": lv,
                     "early_j": e.j.to_numpy(), "late_j": late.j.to_numpy()})
        print(f"  {d:<7}{a:<48}{ev:>9.2f}{lv:>9.2f}{lv - ev:>+9.2f}")
    med_move = float(np.median(moves)) if moves else 0.0
    print(f"\n  trajectories with BOTH endpoints: {len(traj)}   median |move| {med_move:.2f}")
    res: dict = {"n_trajectories": len(traj), "median_abs_move": med_move,
                 "trajectories": [{k: v for k, v in t.items() if not k.endswith("_j")}
                                  for t in traj]}
    if len(traj) < 6 or med_move < MIN_MOVE:
        res["verdict"] = "PRECONDITION FAILED"
        print(f"\n  -> PRECONDITION FAILED (need >=6 trajectories and median move >= {MIN_MOVE}). "
              "There is\n     nothing to predict across this window. Reported as the result.")
    else:
        print("  -> precondition met; proceeding")
        panel = [g for g in GenePanel.load(str(dre.PANEL_PATH)).genes if g in set(genes)]
        gi = {g: i for i, g in enumerate(genes)}
        cols = np.array([gi[g] for g in panel])
        X = np.vstack([expr[t["early_j"], :][:, cols].mean(0) for t in traj])
        y = np.array([t["late"] for t in traj], float)
        donor = np.array([t["donor"] for t in traj])
        # cluster on DONOR: hold out every trajectory of one donor at once
        print("\n[FORWARD] predict late score from early expression, leave-one-DONOR-out")
        print(f"  {len(traj)} trajectories, {len(set(donor))} donor clusters, "
              f"{X.shape[1]} features")
        print(f"  {'alpha':>10}{'LODO spearman':>16}{'null p95':>11}   pass")
        n_pass, per = 0, {}
        rng = np.random.default_rng(0)

        def lodo(Xa, ya, perm=None):
            yy = ya[perm] if perm is not None else ya
            pred = np.empty(len(yy))
            for d in set(donor):
                te = donor == d
                tr = ~te
                if tr.sum() < 2 or np.std(yy[tr]) == 0:
                    pred[te] = yy[tr].mean() if tr.sum() else 0.0
                    continue
                pred[te] = dre.ridge_fit_predict(Xa[tr], yy[tr], Xa[te], alpha)
            return dre.spearman(yy, pred)

        for alpha in dre.ALPHAS:
            obs = lodo(X, y)
            null = np.array([lodo(X, y, rng.permutation(len(y))) for _ in range(500)])
            null = null[np.isfinite(null)]
            p95 = float(np.percentile(null, 95)) if len(null) else float("nan")
            ok = bool(np.isfinite(obs) and np.isfinite(p95) and obs > p95)
            n_pass += ok
            per[str(alpha)] = {"observed": obs, "null_p95": p95, "pass": ok}
            print(f"  {alpha:>10.0f}{obs:>16.3f}{p95:>11.3f}   {'YES' if ok else 'no'}")
        res["forward"] = {"alphas": per, "n_pass": n_pass}
        raw_verdict = "SIGNAL" if n_pass > len(dre.ALPHAS) / 2 else "NULL"
        print(f"  -> raw verdict {raw_verdict} ({n_pass}/{len(dre.ALPHAS)} alphas)")

        # ---- CONTROL 1: is this just PERSISTENCE? ------------------------------------------ #
        # Median |move| is 3.60 against a between-trajectory spread of ~40, so "late ~ early" is
        # nearly true by default. If a ONE-FEATURE model using only the early SCORE does as well as
        # 1903 genes, the "forward signal" is a trajectory staying put -- not a prediction.
        early_score = np.array([t["early"] for t in traj], float)[:, None]
        base_rho = {}
        for alpha in dre.ALPHAS:
            base_rho[str(alpha)] = lodo(early_score, y)
        best_base = max(v for v in base_rho.values() if np.isfinite(v))
        best_full = max(v["observed"] for v in per.values() if np.isfinite(v["observed"]))
        print(f"\n[CONTROL 1] persistence -- early SCORE alone (1 feature) vs {X.shape[1]} genes")
        print(f"  early-score-only best LODO spearman : {best_base:+.3f}")
        print(f"  full-expression   best LODO spearman: {best_full:+.3f}")
        print(f"  spearman(early, late) directly       : {dre.spearman(early_score[:, 0], y):+.3f}")
        adds = best_full > best_base + 0.05
        print(f"  -> expression {'ADDS over' if adds else 'does NOT add over'} persistence")

        # ---- CONTROL 2: a null that preserves donor structure ------------------------------ #
        # The null above permutes across all 17 trajectories, destroying the donor-level offsets
        # visible in the table (O1 low, O2/O3 high). Beating that null can mean "the model knows
        # donors differ", which LODO should not reward. Permuting WITHIN donor holds those offsets
        # fixed and asks the sharper question: does it order ARMS inside a donor?
        def within_donor_perm(rng_):
            p = np.arange(len(y))
            for d in set(donor):
                idx = np.where(donor == d)[0]
                p[idx] = rng_.permutation(idx)
            return p

        print("\n[CONTROL 2] permutation null that PRESERVES donor structure (shuffle within donor)")
        print(f"  {'alpha':>10}{'LODO spearman':>16}{'null p95':>11}   pass")
        n_pass2, per2 = 0, {}
        rng2 = np.random.default_rng(1)
        for alpha in dre.ALPHAS:
            obs = lodo(X, y)
            nl = np.array([lodo(X, y, within_donor_perm(rng2)) for _ in range(500)])
            nl = nl[np.isfinite(nl)]
            p95 = float(np.percentile(nl, 95)) if len(nl) else float("nan")
            ok = bool(np.isfinite(obs) and np.isfinite(p95) and obs > p95)
            n_pass2 += ok
            per2[str(alpha)] = {"observed": obs, "null_p95": p95, "pass": ok}
            print(f"  {alpha:>10.0f}{obs:>16.3f}{p95:>11.3f}   {'YES' if ok else 'no'}")
        print(f"  -> {n_pass2}/{len(dre.ALPHAS)} alphas beat the structure-preserving null")

        res["control_persistence"] = {"early_only_best": best_base, "full_best": best_full,
                                      "expression_adds": bool(adds)}
        res["control_within_donor_null"] = {"alphas": per2, "n_pass": n_pass2}
        res["verdict"] = ("SIGNAL" if (raw_verdict == "SIGNAL" and adds
                                       and n_pass2 > len(dre.ALPHAS) / 2) else "NULL")
        print(f"\n  FINAL (all three must hold): raw {raw_verdict}, adds-over-persistence {adds}, "
              f"structure-preserving null {n_pass2}/{len(dre.ALPHAS)}")
        print(f"  -> {res['verdict']}")
        print("  NOTE: 3 donor clusters. Severely underpowered, as stated before the run.")

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_phase3_within_donor_forward_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
