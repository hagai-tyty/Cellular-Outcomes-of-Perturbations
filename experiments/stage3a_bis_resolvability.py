"""STAGE 3a-bis — is the forward gate ANSWERABLE once the 42,481 HFF cells are used?

    python experiments/stage3a_bis_resolvability.py            # operative arm (_c7)
    python experiments/stage3a_bis_resolvability.py _armA

READ-ONLY. Writes `results/stage3a_bis_resolvability_results.json`. No retrain, no rebuild,
`src/` untouched. **This does NOT grade 3a.** It runs the `REF_GROUND_RULES.md` §5b check that
must precede any bar, at the geometry a repaired 3a would actually be graded on.

WHY THIS EXISTS
---------------
`stage3a_diagnose.py` (2026-08-12) established that 3a's estimator diverges, that its bar is
UNRESOLVABLE at the geometry it was run on, and — §D0 — that a forward time signal IS present on
the five folds it graded. What it did not ask is why the geometry was that small.

`test18_forward_gate.py:74` builds every row from `gather_split(..., REGIME, "test")`. Under the
`holdout` regime the test split is, by construction, **the held-out Gill donor alone**. So 3a ran
on 18–21 bulk samples per fold at **1.7 cells per timepoint** — while the same bundle's train
split holds **33,613 HFF single cells across 9 timepoints at ~3,700 cells per timepoint**, whose
unsafe fraction runs 0.0835 → 0.9996 with a per-timepoint SE of 0.006. The gate was decided on
~0.3 % of the available cells, and the 99.7 % it ignored are the precise ones.

THE QUESTION THIS ANSWERS
-------------------------
Not "is there a signal" — §D0 already measured one. It is: **which forward question can this
corpus actually resolve, and is Gill's cells-per-timepoint the binding constraint?**

  Regime A — CROSS-LINE (what the product needs). Train on HFF's dense trajectory plus the other
             Gill donors; hold out one Gill donor. 6 folds. The held-out target is measured on
             1-2 cells per timepoint.
  Regime B — WITHIN-HFF (what the corpus measures precisely). Partition HFF's cells into K
             disjoint pseudo-replicate trajectories; hold one out. The held-out target is
             measured on ~470 cells per timepoint.

Pseudo-replicates are a TRAINING-SIDE and WITHIN-LINE device only. They are one culture split up,
not independent cultures, so regime B answers "can the forward curve be recovered at all", never
"does it transfer to a new line". Regime A is the only one that speaks to the product, and it is
reported as such.

THE §5b SIMULATION, AND WHY IT IS GROUNDED RATHER THAN ASSUMED
--------------------------------------------------------------
A system that "meets the intent exactly" here is one where the unsafe probability really is a
function of elapsed time. So the simulated truth is **HFF's own measured curve** — no invented
effect size:

    p(t_j) = mean_g + alpha * (g(t_j) - mean_g),   g = HFF's measured unsafe fraction

and the observation is drawn at the **real cell counts**, `u_j ~ Binomial(n_j, p(t_j)) / n_j`, so
Gill's 1-2 cells and HFF's ~470 enter as the measurement noise they actually are. `alpha` scales
the amplitude: alpha = 1 is exactly the curve HFF measures, alpha = 0.25 a quarter of it,
alpha = 0 no time effect at all (the false-positive check).

Under this truth `state + dt` can reach `t_j` and `state` alone cannot, so a working test MUST
detect it. 3a's rule is graded verbatim: PASS iff the paired 95 % CI upper end < 0.

PRE-REGISTERED EXPECTATIONS (written before the run; graded in the output)
-------------------------------------------------------------------------
  Q1  Regime B is RESOLVABLE at alpha = 1 for the bounded (logit) estimator. If the corpus cannot
      recover its own measured curve when the target is precise, nothing else here matters.
  Q2  Regime A is NOT resolvable at alpha = 1 for the RAW estimator — the divergence §D1 located
      is a property of the estimator, not of the fold size, so more training data should not fix
      it on its own.
  Q3  The binding constraint is the held-out cells-per-timepoint, not the training-set size:
      regime A's pass rate stays below regime B's at every alpha, for the bounded estimator.
  Q4  alpha = 0 gives a pass rate at or below 0.05 in every cell (no false-positive problem),
      as it did in the diagnosis.

Stated so they can fail. Two of five failed last time and were recorded as failures.
"""
from __future__ import annotations

import contextlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "stage3a_bis_resolvability_results.json"

K_PSEUDO = 10          # HFF pseudo-replicate trajectories
SEED = 0
SIM_TRIALS = 2000
ALPHAS = (0.0, 0.25, 0.5, 1.0)
MIN_TP = 3             # a trajectory needs >= 3 timepoints to make forward pairs


def _t18(suffix: str):
    if "test18_forward_gate" in sys.modules:
        del sys.modules["test18_forward_gate"]
    m = importlib.import_module("test18_forward_gate")
    m.resolve_root = lambda name, _s=suffix: str(REPO / f"{name}{_s}")
    return m


def _diag():
    if "stage3a_diagnose" in sys.modules:
        return sys.modules["stage3a_diagnose"]
    return importlib.import_module("stage3a_diagnose")


def trajectory_rows(t, X, cls, y, mask, sel, unsafe_idx) -> list[dict]:
    """Per-timepoint population means for one trajectory: expression, unsafe fraction, n."""
    rows = []
    for tp in np.unique(np.round(t[sel], 6)):
        s = sel & np.isclose(t, tp)
        if not s.any():
            continue
        n = int(s.sum())
        am = s & mask
        rows.append({"t": float(tp), "x": np.asarray(X[s].mean(0), float), "n": n,
                     "u": float(unsafe_idx[s].mean()),
                     "y": float(y[am].mean()) if am.any() else float("nan"),
                     "n_age": int(am.sum())})
    return rows


def partition_within_timepoints(t, sel, k: int, seed: int = SEED) -> np.ndarray:
    """Split the selected cells into `k` pseudo-replicates, balanced WITHIN each timepoint.

    Every replicate must carry the WHOLE time course, not a slice of it -- a partition that split
    cells globally would hand some replicates only early timepoints and make them incomparable as
    trajectories. Returns a replicate id per cell, `-1` where `sel` is False.
    """
    rng = np.random.default_rng(seed)
    rep = np.full(len(t), -1, int)
    for tp in np.unique(np.round(np.asarray(t)[sel], 6)):
        idx = np.flatnonzero(sel & np.isclose(t, tp))
        rep[rng.permutation(idx)] = np.arange(len(idx)) % k
    return rep


def assemble_fold(t18, donor: str) -> dict | None:
    """One fold: TRAIN trajectories (HFF pseudo-reps + the other Gill donors) and the TEST donor.

    The bundle's own split assignment supplies the leakage barrier -- `holdout/test` is the held
    out donor and nothing else, so every training trajectory here is already disjoint from it.
    """
    from cellfate.common.constants import DEATH_IDX, LOSS_IDX
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split

    paths = ArtifactPaths.of(t18.resolve_root(f"cellfate_loocv_{donor}"))
    try:
        tr_parts = [gather_split(paths, t18.REGIME, s) for s in ("train", "val", "calib")]
        te = gather_split(paths, t18.REGIME, "test")
    except Exception:                                                   # noqa: BLE001
        return None

    def stack(parts, attr, cast=float):
        return np.concatenate([np.asarray(getattr(p, attr), cast) for p in parts])

    t = np.concatenate([np.asarray(p.dose_time[:, 1], float) for p in tr_parts])
    # kept float32 (the shard dtype): 42k cells x 2000 genes is 340 MB here and 680 MB as
    # float64, and only the per-timepoint MEANS need double precision.
    X = np.concatenate([np.asarray(p.X, np.float32) for p in tr_parts])
    cls = np.concatenate([p.y_cls.astype(int) for p in tr_parts])
    y = stack(tr_parts, "y_age")
    mask = np.concatenate([np.asarray(p.mask, bool) for p in tr_parts])
    line = np.concatenate([np.asarray([str(v) for v in p.cell_line]) for p in tr_parts])
    unsafe = ((cls == DEATH_IDX) | (cls == LOSS_IDX)).astype(float)

    train: dict[str, list[dict]] = {}

    # HFF -> K pseudo-replicate trajectories, partitioned WITHIN each timepoint so every
    # replicate carries the whole time course rather than a slice of it.
    h = line == "HFF"
    if h.any():
        rep = partition_within_timepoints(t, h, K_PSEUDO)
        for k in range(K_PSEUDO):
            rows = trajectory_rows(t, X, cls, y, mask, rep == k, unsafe)
            if len(rows) >= MIN_TP:
                train[f"HFF_r{k}"] = rows

    for ln in sorted(set(line[~h])):
        rows = trajectory_rows(t, X, cls, y, mask, line == ln, unsafe)
        if len(rows) >= MIN_TP:
            train[ln] = rows

    tt = np.asarray(te.dose_time[:, 1], float)
    te_rows = trajectory_rows(tt, np.asarray(te.X, float), te.y_cls.astype(int),
                              np.asarray(te.y_age, float), np.asarray(te.mask, bool),
                              np.ones(te.n, bool),
                              ((te.y_cls.astype(int) == DEATH_IDX)
                               | (te.y_cls.astype(int) == LOSS_IDX)).astype(float))
    if len(te_rows) < MIN_TP:
        return None
    return {"train": train, "test": te_rows}


def pairs_of(rows: list[dict]) -> list[dict]:
    out = []
    for i in range(len(rows)):
        for j in range(len(rows)):
            if rows[j]["t"] <= rows[i]["t"]:
                continue
            out.append({"x_i": rows[i]["x"], "dt": rows[j]["t"] - rows[i]["t"],
                        "t_j": rows[j]["t"], "n_j": rows[j]["n"], "u_j": rows[j]["u"]})
    return out


def hff_curve(fold: dict):
    """g(t): HFF's own measured unsafe fraction, pooled over the pseudo-replicates."""
    acc: dict[float, list[tuple[float, int]]] = {}
    for name, rows in fold["train"].items():
        if not name.startswith("HFF_r"):
            continue
        for r in rows:
            acc.setdefault(round(r["t"], 6), []).append((r["u"], r["n"]))
    ts = np.array(sorted(acc), float)
    g = np.array([sum(u * n for u, n in acc[k]) / sum(n for _, n in acc[k])
                  for k in sorted(acc)], float)
    return ts, g


def main() -> int:
    from audit_metrics import MIN_PASS_RATE, bar_verdict

    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    D = _diag()

    suffix = sys.argv[1] if len(sys.argv) > 1 else "_c7"
    t18 = _t18(suffix)
    print("\n" + "=" * 92)
    print(f"STAGE 3a-bis — §5b RESOLVABILITY at the geometry a repaired 3a would use   arm:"
          f" {suffix}")
    print("=" * 92)
    print("READ-ONLY. This does NOT grade 3a. It asks which forward question this corpus can")
    print("resolve, with the 42,481 HFF cells included and the real cell counts as the noise.")

    folds = {}
    for d in t18.DONORS:
        if not (REPO / f"cellfate_loocv_{d}{suffix}").exists():
            continue
        f = assemble_fold(t18, d)
        if f is not None:
            folds[d] = f
            print(f"   assembled fold {d}: {len(f['train'])} training trajectories, "
                  f"{len(f['test'])} held-out timepoints")
    if len(folds) < 2:
        print("   not enough folds assembled")
        return 1

    out: dict = {"script": "stage3a_bis_resolvability", "arm": suffix,
                 "k_pseudo": K_PSEUDO, "trials": SIM_TRIALS, "alphas": list(ALPHAS)}

    # ---- what changed: cells per timepoint on each side -------------------------------------
    any_fold = folds[next(iter(folds))]
    geo_rows = []
    for name, rows in sorted(any_fold["train"].items()):
        ns = [r["n"] for r in rows]
        geo_rows.append([name, "train", str(len(rows)), f"{np.mean(ns):.1f}",
                         f"{len(pairs_of(rows))}"])
    for d, f in folds.items():
        ns = [r["n"] for r in f["test"]]
        geo_rows.append([d, "HELD OUT", str(len(f['test'])), f"{np.mean(ns):.1f}",
                         f"{len(pairs_of(f['test']))}"])
    print("\n  GEOMETRY — cells per timepoint is the quantity 3a never had")
    print(render_table(["trajectory", "role", "timepoints", "cells/tp", "pairs"], geo_rows,
                       aligns=["l", "l", "r", "r", "r"]))
    out["geometry"] = {"rows": geo_rows}

    ts_g, g = hff_curve(any_fold)
    print("\n  THE SIMULATED TRUTH — HFF's OWN measured unsafe-fraction curve (alpha scales it)")
    print(render_table(["log t", "day", "g(t)"],
                       [[f"{x:.3f}", f"{np.exp(x) / 24.0:.1f}", f"{v:.4f}"]
                        for x, v in zip(ts_g, g, strict=True)], aligns=["r", "r", "r"]))
    out["hff_curve"] = {"log_t": ts_g.tolist(), "g": g.tolist()}
    g_mean = float(np.average(g))

    def p_of(t_arr, alpha):
        return np.clip(g_mean + alpha * (np.interp(np.asarray(t_arr, float), ts_g, g) - g_mean),
                       0.01, 0.99)

    # ---- the two regimes ---------------------------------------------------------------------
    def build_regime(regime: str):
        """-> list of {'s','d','n_te','t_te','n_tr','t_tr'} per fold, features frozen."""
        prep = []
        for d, f in folds.items():
            if regime == "A":
                tr_rows = list(f["train"].values())
                te_rows = f["test"]
                tag = d
            else:                                    # B: hold out one HFF pseudo-replicate
                names = [n for n in f["train"] if n.startswith("HFF_r")]
                if not names:
                    continue
                # one fold per pseudo-replicate, but only from the FIRST bundle -- the HFF cells
                # and their split are the same object in every bundle, so repeating them per
                # donor would be the same fold counted six times.
                if d != next(iter(folds)):
                    continue
                for hold in names:
                    tr_rows = [f["train"][n] for n in f["train"] if n != hold]
                    prep.append(_prep(tr_rows, f["train"][hold], hold))
                continue
            prep.append(_prep(tr_rows, te_rows, tag))
        return prep

    def _prep(tr_rows, te_rows, tag):
        tr = [p for rows in tr_rows for p in pairs_of(rows)]
        te = pairs_of(te_rows)
        return {
            "tag": tag, "n_tr_pairs": len(tr), "n_te_pairs": len(te),
            "s": D.FrozenRidge(t18.feats(tr, False), t18.feats(te, False)),
            "d": D.FrozenRidge(t18.feats(tr, True), t18.feats(te, True)),
            "t_tr": np.array([p["t_j"] for p in tr], float),
            "t_te": np.array([p["t_j"] for p in te], float),
            "n_tr": np.array([max(p["n_j"], 1) for p in tr], int),
            "n_te": np.array([max(p["n_j"], 1) for p in te], int),
        }

    rng = np.random.default_rng(SEED)
    for regime, label in (("B", "WITHIN-HFF — held-out target on ~470 cells/timepoint"),
                          ("A", "CROSS-LINE — held-out target on 1-2 cells/timepoint")):
        prep = build_regime(regime)
        if len(prep) < 2:
            continue
        print(f"\n  REGIME {regime} — {label}   [{len(prep)} folds, "
              f"{prep[0]['n_tr_pairs']} training pairs]")
        rows = []
        for alpha in ALPHAS:
            hi_raw, hi_log = [], []
            lat = [{"tr": p_of(p["t_tr"], alpha), "te": p_of(p["t_te"], alpha)} for p in prep]
            for _ in range(SIM_TRIALS):
                ds_r, dd_r, ds_l, dd_l = [], [], [], []
                for p, L in zip(prep, lat, strict=True):
                    ytr = rng.binomial(p["n_tr"], L["tr"]) / p["n_tr"]
                    yte = rng.binomial(p["n_te"], L["te"]) / p["n_te"]
                    ds_r.append(D.mae(p["s"].predict(ytr), yte))
                    dd_r.append(D.mae(p["d"].predict(ytr), yte))
                    ds_l.append(D.mae(D.expit(p["s"].predict(D.logit(ytr))), yte))
                    dd_l.append(D.mae(D.expit(p["d"].predict(D.logit(ytr))), yte))
                hi_raw.append(t18.paired_ci([b - a for a, b in zip(ds_r, dd_r, strict=True)])[1][1])
                hi_log.append(t18.paired_ci([b - a for a, b in zip(ds_l, dd_l, strict=True)])[1][1])
            vr = bar_verdict(np.array(hi_raw, float), 0.0, lower_is_better=True)
            vl = bar_verdict(np.array(hi_log, float), 0.0, lower_is_better=True)
            rows.append([f"{alpha:.2f}", f"{vr['pass_rate']:.3f}", vr["verdict"],
                         f"{vl['pass_rate']:.3f}", vl["verdict"]])
            out.setdefault("regimes", {}).setdefault(regime, {})[f"alpha={alpha}"] = {
                "raw": vr, "logit": vl}
        print(render_table(["alpha (curve amplitude)", "pass (raw)", "verdict (raw)",
                            "pass (logit)", "verdict (logit)"], rows,
                           aligns=["r", "r", "l", "r", "l"]))
        out["regimes"][regime]["folds"] = [p["tag"] for p in prep]
        out["regimes"][regime]["n_train_pairs"] = prep[0]["n_tr_pairs"]
        out["regimes"][regime]["n_test_pairs"] = [p["n_te_pairs"] for p in prep]

    # ---- pre-registered questions, graded ----------------------------------------------------
    R = out.get("regimes", {})

    def rate(reg, alpha, est):
        return R.get(reg, {}).get(f"alpha={alpha}", {}).get(est, {}).get("pass_rate", float("nan"))

    q = {
        "Q1_regimeB_resolvable_at_alpha1_logit": {
            "held": bool(rate("B", 1.0, "logit") >= MIN_PASS_RATE),
            "pass_rate": rate("B", 1.0, "logit")},
        "Q2_regimeA_raw_not_resolvable_at_alpha1": {
            "held": bool(rate("A", 1.0, "raw") < MIN_PASS_RATE),
            "pass_rate": rate("A", 1.0, "raw")},
        "Q3_heldout_cells_bind_not_training_size": {
            "held": bool(all(rate("A", a, "logit") <= rate("B", a, "logit")
                             for a in ALPHAS if a > 0)),
            "A": [rate("A", a, "logit") for a in ALPHAS],
            "B": [rate("B", a, "logit") for a in ALPHAS]},
        "Q4_no_false_positive_at_alpha0": {
            "held": bool(max([rate(r, 0.0, e) for r in ("A", "B") for e in ("raw", "logit")]
                             or [1.0]) <= 0.05),
            "rates": [rate(r, 0.0, e) for r in ("A", "B") for e in ("raw", "logit")]},
    }
    print("\n  PRE-REGISTERED QUESTIONS, GRADED")
    print(render_table(["question", "held?", "evidence"],
                       [[k, "YES" if v["held"] else "NO",
                         ", ".join(f"{a}={b}" for a, b in v.items() if a != "held")]
                        for k, v in q.items()], aligns=["l", "l", "l"]))
    out["pre_registered"] = q

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   NO 3a VERDICT IS TAKEN HERE. This is the §5b precondition; grading needs a bar")
    print("   registered in tests/test_bars_resolvable.py first.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
