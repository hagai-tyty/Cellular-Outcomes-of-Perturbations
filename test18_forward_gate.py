"""
Test 18 — THE GATE for the stopping-time tool: does a forward Δt signal exist at all?

    python test18_forward_gate.py

WHY THIS RUNS BEFORE ANY TRAINING CODE.
The product we want is: "given my culture as it is now, when should I withdraw?" That requires the
model to answer a FORWARD question — state at t_i plus a gap Δt maps to the outcome at t_j.

Test 11.1 showed the CURRENT model ignores its time input entirely (ΔAge shifts 0.035 yr across a
full time sweep) because time is redundant with state along a single trajectory. The hope is that
retraining on (t_i → t_j) PAIRS fixes this: with pairs, the same starting state appears with
DIFFERENT Δt and different targets, so Δt is no longer redundant.

That is a hypothesis, not a fact. This test checks it cheaply, with ridge, before we build
anything. Rationale: this project has established repeatedly (Tests 3, 5, 6, 9) that at this data
scale ridge matches or beats flexible models. **If ridge cannot find a Δt signal here, a neural
net will not either.**

METHOD (population-level, because sampling is destructive — the same cell is never seen twice):
  1. Group each donor's held-out cells by timepoint; compute per-(donor, timepoint) mean
     expression and mean true ΔAge.
  2. Build every ordered pair (t_i < t_j) within a donor.
  3. Leave-one-donor-out ridge on two feature sets:
        STATE      = [mean expression at t_i]
        STATE+DT   = [mean expression at t_i, Δt, Δt²]
     Target = mean true ΔAge at t_j.
  4. Part B (decisive): hold the starting state FIXED and sweep Δt across its range. Measure how
     far the prediction moves. This is the forward-response analogue of Test 11.1 Part 3.

READ:
  - STATE+DT beats STATE, AND the Δt sweep moves the prediction materially
        -> a forward signal EXISTS. Build the stopping-time model. GO.
  - STATE+DT ties STATE, and/or the sweep is flat
        -> no learnable forward signal at this data scale. The stopping-time tool is NOT
           supported by this dataset. STOP before writing training code.

HONEST POWER NOTE: ~12 timepoints × 6 donors. Pairs are not independent (they share timepoints),
so treat this as a screen, not a proof. A negative result is decisive; a positive one is
permission to proceed, not a validation.

USAGE (repo root, venv active; needs cellfate_loocv_* bundles).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from cellfate.common.constants import DEATH_IDX, LOSS_IDX
from cellfate.common.io import ArtifactPaths
from cellfate.evaluation.data import gather_split

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}
MIN_CELLS_PER_TP = 1


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def timepoint_table(donor: str):
    """Per-(donor, timepoint) population means: expression, true ΔAge, time."""
    try:
        te = gather_split(ArtifactPaths.of(resolve_root(f"cellfate_loocv_{donor}")),
                          REGIME, "test")
    except Exception:  # noqa: BLE001
        return None
    m = te.mask
    if m.sum() < 4:
        return None
    t = np.asarray(te.dose_time[:, 1], float)[m]      # log time_h
    X = np.asarray(te.X, float)[m]
    y = np.asarray(te.y_age, float)[m]
    cls = te.y_cls[m].astype(int)
    unsafe = ((cls == LOSS_IDX) | (cls == DEATH_IDX)).astype(float)   # the SAFETY target
    tps = np.unique(np.round(t, 6))
    if len(tps) < 3:
        return None
    rows = []
    for tp in tps:
        sel = np.isclose(t, tp)
        if sel.sum() < MIN_CELLS_PER_TP:
            continue
        rows.append({"t": float(tp), "x": X[sel].mean(0), "y": float(y[sel].mean()),
                     "u": float(unsafe[sel].mean()), "n": int(sel.sum())})
    return rows if len(rows) >= 3 else None


def build_pairs(rows):
    """Every ordered (t_i -> t_j) pair within a donor."""
    out = []
    for i in range(len(rows)):
        for j in range(len(rows)):
            if rows[j]["t"] <= rows[i]["t"]:
                continue
            out.append({"x_i": rows[i]["x"], "dt": rows[j]["t"] - rows[i]["t"],
                        "y_j": rows[j]["y"], "u_j": rows[j]["u"]})
    return out


def feats(pairs, with_dt: bool):
    X = np.vstack([p["x_i"] for p in pairs])
    if not with_dt:
        return X
    dt = np.array([p["dt"] for p in pairs], float).reshape(-1, 1)
    return np.hstack([X, dt, dt ** 2])


def paired_ci(diffs):
    diffs = [d for d in diffs if np.isfinite(d)]
    n = len(diffs)
    if n < 2:
        return float("nan"), (float("nan"), float("nan")), n
    md = float(np.mean(diffs))
    sd = float(np.std(diffs, ddof=1))
    se = sd / math.sqrt(n)
    t = T_CRIT.get(n - 1, 2.571)
    return md, (md - t * se, md + t * se), n


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\nTEST 18 — THE GATE: does a forward Δt signal exist in this data?")
    print("population-level (t_i -> t_j) pairs, leave-one-donor-out ridge.")
    print("If this fails, the stopping-time tool is not supported by this dataset.")

    per = {}
    for d in DONORS:
        rows = timepoint_table(d)
        if rows is None:
            continue
        per[d] = {"rows": rows, "pairs": build_pairs(rows)}
    if len(per) < 2:
        print("\n   Not enough donors with >=3 timepoints.")
        return

    print("\n  DATA AVAILABLE")
    print(render_table(["fold", "timepoints", "cells", "forward pairs"],
                       [[d, str(len(per[d]["rows"])),
                         str(sum(r["n"] for r in per[d]["rows"])),
                         str(len(per[d]["pairs"]))] for d in per],
                       aligns=["l", "r", "r", "r"]))
    total = sum(len(per[d]["pairs"]) for d in per)
    print(f"   total forward pairs across donors: {total}")

    # ---- LOAD-BEARING CHECK: cells per timepoint ----------------------------------------
    # The condition-level design (MASTER_PLAN §5b-ter) claims population uncertainty of
    # ±3.7-4.6 yr, derived as q/sqrt(n) with n=21. That assumed 21 cells AT ONE TIMEPOINT.
    # If the cells are spread across the time course, the per-timepoint n is far smaller and
    # the uncertainty can EXCEED the ~11.35 yr effect -- in which case the tool cannot resolve
    # anything, regardless of how well Parts A-C do.
    print("\n  CELLS PER TIMEPOINT — the load-bearing assumption of the whole design")
    rows_cpt, worst = [], 0.0
    for d in per:
        n_tp = len(per[d]["rows"])
        n_cells = sum(r["n"] for r in per[d]["rows"])
        cpt = n_cells / max(n_tp, 1)
        for q in (17.0, 21.0):
            se = q / np.sqrt(max(cpt, 1e-9))
            worst = max(worst, se)
        se17, se21 = 17.0 / np.sqrt(cpt), 21.0 / np.sqrt(cpt)
        rows_cpt.append([d, str(n_tp), str(n_cells), f"{cpt:.1f}",
                         f"{se17:.1f}-{se21:.1f}",
                         "OK" if se21 < 11.35 else "** EXCEEDS EFFECT **"])
    print(render_table(["fold", "timepoints", "cells", "cells/tp", "SE at one tp (yr)",
                        "vs 11.35 effect"], rows_cpt,
                       aligns=["l", "r", "r", "r", "r", "l"]))
    if worst >= 11.35:
        print("   !! WARNING: at some donors the per-timepoint uncertainty EXCEEDS the effect")
        print("      size. The condition-level argument (±3.7-4.6 yr) assumed all cells at one")
        print("      timepoint; they are spread across the course. Mitigations, in order:")
        print("        (a) POOL adjacent timepoints into wider bins (trades time resolution")
        print("            for statistical resolution) -- usually the right call")
        print("        (b) report only COARSE recommendations (early / mid / late)")
        print("        (c) if neither is acceptable, the tool cannot resolve withdrawal days")
    else:
        print("   OK: per-timepoint uncertainty stays below the effect size.")

    # ---- PART A: does adding Δt improve forward prediction? ----
    rows, d_state, d_both = [], [], []
    for hd in per:
        tr_pairs = [p for d in per if d != hd for p in per[d]["pairs"]]
        te_pairs = per[hd]["pairs"]
        if len(tr_pairs) < 8 or len(te_pairs) < 3:
            continue
        ytr = np.array([p["y_j"] for p in tr_pairs])
        yte = np.array([p["y_j"] for p in te_pairs])
        res = {}
        for lbl, wdt in (("state", False), ("state+dt", True)):
            Xtr, Xte = feats(tr_pairs, wdt), feats(te_pairs, wdt)
            sc = StandardScaler().fit(Xtr)
            pred = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr).predict(sc.transform(Xte))
            res[lbl] = float(np.abs(pred - yte).mean())
        d_state.append(res["state"])
        d_both.append(res["state+dt"])
        rows.append([hd, f"{len(te_pairs)}", f"{res['state']:.2f}", f"{res['state+dt']:.2f}",
                     f"{res['state'] - res['state+dt']:+.2f}"])
    if not rows:
        print("\n   Too few pairs for leave-one-donor-out.")
        return
    print("\n  PART A — forward ΔAge MAE (lower better): does Δt add anything?")
    print(render_table(["held-out", "pairs", "state only", "state + Δt", "gain"],
                       rows, aligns=["l", "r", "r", "r", "r"]))
    diffs = [b - a for a, b in zip(d_state, d_both, strict=True)]   # negative = Δt helps
    md, (lo, hi), n = paired_ci(diffs)
    verdict_a = ("Δt HELPS" if hi < 0 else "Δt HURTS" if lo > 0 else "tied (no Δt signal)")
    print(f"   aggregate: state={np.mean(d_state):.2f}  state+Δt={np.mean(d_both):.2f}")
    print(f"   paired (state+Δt − state): mean={md:+.2f} 95% CI=[{lo:+.2f},{hi:+.2f}] (n={n})"
          f"  -> {verdict_a}")

    # ---- PART B (decisive): sweep Δt with the starting state held FIXED ----
    all_pairs = [p for d in per for p in per[d]["pairs"]]
    ytr = np.array([p["y_j"] for p in all_pairs])
    Xtr = feats(all_pairs, True)
    sc = StandardScaler().fit(Xtr)
    mdl = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr)
    dts = np.array([p["dt"] for p in all_pairs])
    lo_dt, hi_dt = float(np.quantile(dts, 0.10)), float(np.quantile(dts, 0.90))
    rows = []
    for d in per:
        x0 = per[d]["rows"][0]["x"]                     # this donor's earliest state
        preds = []
        for dt in (lo_dt, hi_dt):
            f = np.hstack([x0, [dt, dt ** 2]]).reshape(1, -1)
            preds.append(float(mdl.predict(sc.transform(f))[0]))
        rows.append([d, f"{preds[0]:+.2f}", f"{preds[1]:+.2f}", f"{preds[1] - preds[0]:+.2f}"])
    print(f"\n  PART B — sweep Δt from {lo_dt:.2f} to {hi_dt:.2f} (log-h), START STATE HELD FIXED")
    print(render_table(["fold", "ΔAge @ short Δt", "ΔAge @ long Δt", "swing"],
                       rows, aligns=["l", "r", "r", "r"]))
    swing = float(np.mean([abs(float(r[3])) for r in rows]))
    print(f"   mean |swing| = {swing:.2f} yr   (Test 11.1's current model managed 0.035 yr)")

    # ---- PART C (THE ONE THAT MATTERS): can Δt predict the UNSAFE FRACTION forward? ----
    # The recommendation is gated by `p_unsafe <= risk_threshold`, NOT by ΔAge -- and the ΔAge
    # intervals overlap so heavily between adjacent days that they cannot separate options.
    # So forward SAFETY prediction is the capability the tool actually runs on. If Δt cannot
    # predict it, the recommender has nothing to recommend with, even if Part A passes.
    rows_u, u_state, u_both = [], [], []
    for hd in per:
        tr_pairs = [p for d in per if d != hd for p in per[d]["pairs"]]
        te_pairs = per[hd]["pairs"]
        if len(tr_pairs) < 8 or len(te_pairs) < 3:
            continue
        utr = np.array([p["u_j"] for p in tr_pairs])
        ute = np.array([p["u_j"] for p in te_pairs])
        if np.std(utr) == 0 or np.std(ute) == 0:
            rows_u.append([hd, f"{len(te_pairs)}", "n/a", "n/a", "no variation"])
            continue
        res = {}
        for lbl, wdt in (("state", False), ("state+dt", True)):
            Xtr, Xte = feats(tr_pairs, wdt), feats(te_pairs, wdt)
            sc = StandardScaler().fit(Xtr)
            pred = Ridge(alpha=1.0).fit(sc.transform(Xtr), utr).predict(sc.transform(Xte))
            res[lbl] = float(np.abs(pred - ute).mean())
        u_state.append(res["state"])
        u_both.append(res["state+dt"])
        rows_u.append([hd, f"{len(te_pairs)}", f"{res['state']:.3f}",
                       f"{res['state+dt']:.3f}", f"{res['state'] - res['state+dt']:+.3f}"])
    print("\n  PART C — forward UNSAFE-FRACTION MAE (lower better)  ***THE DECISIVE ONE***")
    print("  the recommendation is gated by p_unsafe, not ΔAge — this is what it runs on.")
    print(render_table(["held-out", "pairs", "state only", "state + Δt", "gain"],
                       rows_u, aligns=["l", "r", "r", "r", "r"]))
    if len(u_state) >= 2:
        du = [b - a for a, b in zip(u_state, u_both, strict=True)]
        mdu, (lou, hiu), nu = paired_ci(du)
        verdict_c = ("Δt HELPS safety" if hiu < 0 else "Δt HURTS" if lou > 0
                     else "tied (NO forward safety signal)")
        print(f"   aggregate: state={np.mean(u_state):.3f}  state+Δt={np.mean(u_both):.3f}")
        print(f"   paired: mean={mdu:+.3f} 95% CI=[{lou:+.3f},{hiu:+.3f}] (n={nu}) -> {verdict_c}")
        safety_ok = hiu < 0
    else:
        print("   too few folds with unsafe-fraction variation to test")
        safety_ok = False

    print("\n   VERDICT:")
    helps = hi < 0
    moves = swing > 2.0
    if not safety_ok:
        print("     => STOP (or REDESIGN). Δt cannot predict the UNSAFE FRACTION forward.")
        print("        The recommendation is gated by p_unsafe, and ΔAge intervals overlap too")
        print("        heavily to separate adjacent days on their own. Without forward safety")
        print("        prediction the recommender has nothing to recommend with -- even if")
        print("        Parts A/B look fine. Options: (a) ship a ΔAge-trajectory readout with NO")
        print("        safety recommendation, (b) get data with more unsafe-cell variation.")
    elif helps and moves:
        print("     => GO. Forward signal exists for BOTH ΔAge and safety: Δt improves prediction")
        print("        beyond noise, the prediction responds materially to Δt, and the unsafe")
        print("        fraction is forward-predictable. Build the stopping-time model.")
    elif moves and not helps:
        print("     => WEAK GO. The prediction moves with Δt, but adding Δt does not beat")
        print("        state-only beyond noise. The signal may be real but is small relative to")
        print("        state information. Proceed only with tempered expectations.")
    else:
        print("     => STOP. No usable forward Δt signal at this data scale. The stopping-time")
        print("        tool is NOT supported by this dataset. Do not write training code.")
        print("        Options: (a) a dataset with more timepoints per donor, (b) more donors,")
        print("        (c) abandon forward prediction and ship retrospective ranking only.")
    print("\n   CAVEAT: pairs share timepoints so they are not independent — a screen, not a")
    print("   proof. A NEGATIVE result is decisive; a POSITIVE one is permission to proceed.")


if __name__ == "__main__":
    main()
