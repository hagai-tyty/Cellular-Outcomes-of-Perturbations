"""ΔAge ON GSE165177 — the contemporaneous, replicated, in-range controls this project never had.

    python experiments/dage_gse165177.py

READ-ONLY. Writes `results/dage_gse165177_results.json`. No build, no retrain, `src/` untouched.
Graded against `plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md`, committed BEFORE this file existed.

WHY
---
Regime E established that `p_unsafe` is not expressible in bulk at any replication, and left ΔAge
explicitly open -- it is continuous per sample, so the collapse does not apply. `GSE165177` offers
three things at once that map onto D1/D2/D3 of the zero-point fix plan:

  D1  33 CONTEMPORANEOUS negative controls, 2-3 per donor PER TIMEPOINT. Every `gill_bulk` ΔAge is
      measured against a day-0 baseline from a DIFFERENT batch for ~50% of samples.
  D2  Those controls are REPLICATED (n=2-3). `gill_bulk`'s baseline is one unreplicated sample per
      donor, with no error bar and nothing recording that.
  D3  Donors aged 53, 53, 38 -- all inside the clock's fitted [1, 96]. HFF is neonatal (age 0), so
      99.7% of the project's age labels extrapolate past the clock's declared validity.

NORMALISATION -- the thing most likely to corrupt this silently
---------------------------------------------------------------
The clock declares `log1p_cp10k`; GSE165177 ships **Log2 RPM**. Feeding Log2 RPM straight in is a
mismatch worth ~1/ln2 = 1.44x on every age. The pipeline already solves it for the identically
formatted `gill_bulk` at `sources.py:506`:

    rpm = 2**log2 - 1.0   ->   normalize_counts(rpm, target_sum=1e4)   ->   log1p

valid because CP10k(RPM) == CP10k(counts). **That exact path is reused here**, not re-derived.
"""
from __future__ import annotations

import contextlib
import gzip
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "dage_gse165177_results.json"

CLOCK = REPO / "configs" / "clocks" / "fleischer_clock.json"
TRUE_AGES = {"O1": 53.0, "O2": 53.0, "O3": 38.0}
T_CRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
          9: 2.262, 10: 2.228}


def arm_group(arm: str) -> str:
    a = arm.lower()
    if "negative_control" in a:
        return "control"
    if a.startswith("day0"):
        return "day0"
    if "fail" in a:
        return "failed"
    if "transient" in a:
        return "transient"
    return "other"


def ci(vals) -> tuple[float, float, float, int]:
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    n = len(v)
    if n < 2:
        return (float(v[0]) if n else float("nan"), float("nan"), float("nan"), n)
    m = float(v.mean())
    se = float(v.std(ddof=1)) / np.sqrt(n)
    t = T_CRIT.get(n - 1, 1.96)
    return m, m - t * se, m + t * se, n


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.data.aging import LinearClock
    from cellfate.data.integrity import screen_bulk_matrix
    from cellfate.data.normalize import normalize_counts
    install_pretty_console()

    E = importlib.import_module("stage3a_regime_e")     # reuse its verified loader

    print("\n" + "=" * 92)
    print("ΔAge ON GSE165177 — contemporaneous, replicated, in-clock-range controls")
    print("=" * 92)
    print("READ-ONLY. Graded against plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md,")
    print("committed BEFORE this script existed.")

    # ---- load the RAW log2 matrix (for M-E0) and the normalised one (for the clock) ----------
    frames = []
    for fn in E.MATRICES:
        with gzip.open(E.GSE_DIR / fn, "rt", encoding="utf-8", errors="replace") as f:
            head = f.readline().rstrip("\n").split("\t")
        use = [head[0], *head[E.ANNOT_COLS:]]
        df = pd.read_csv(E.GSE_DIR / fn, sep="\t", usecols=use, compression="gzip",
                         low_memory=False).rename(columns={head[0]: "gene"}).set_index("gene")
        frames.append(df.apply(pd.to_numeric, errors="coerce"))
    j = frames[0].join(frames[1], how="inner", lsuffix="_a", rsuffix="_b")
    j = j.assign(_m=j.mean(axis=1)).sort_values("_m", ascending=False)
    j = j[~j.index.duplicated(keep="first")].drop(columns="_m").sort_index().dropna(how="any")
    meta = [E.parse_sample_name(c) for c in j.columns]
    keep = [i for i, m in enumerate(meta) if m is not None]
    obs = pd.DataFrame([meta[i] for i in keep])
    log2 = j.to_numpy(dtype=np.float64)[:, keep]                 # genes x samples
    genes = list(j.index)
    obs["group"] = obs["arm"].map(arm_group)
    print(f"\n   {log2.shape[1]} samples x {log2.shape[0]} genes; "
          f"groups={dict(obs.group.value_counts())}")

    out: dict = {"script": "dage_gse165177",
                 "prereg": "plans/STAGE_1_5_7_DAGE_ON_GSE165177_PREREG.md",
                 "n_samples": int(log2.shape[1]), "n_genes": int(log2.shape[0]),
                 "true_ages": TRUE_AGES}

    # ---- M-E0: the C-7 integrity gate --------------------------------------------------------
    rejected = screen_bulk_matrix(log2, list(obs["sample"]))
    frac = len(rejected) / log2.shape[1]
    print("\n" + "-" * 92)
    print("M-E0 — C-7 INTEGRITY GATE (the gate that caught N2_Fib_Sendai_Exp2)")
    print("-" * 92)
    print(f"   rejected {len(rejected)} of {log2.shape[1]} columns ({100 * frac:.1f}%)"
          + (f": {rejected}" if rejected else " — CLEAN by the project's own standard"))
    out["M_E0"] = {"n_rejected": len(rejected), "rejected": rejected, "frac": frac}
    if frac > 0.10:
        print("   >10% rejected -> systemic quality problem. Per the pre-registration, STOP.")
        OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
        return 0
    if rejected:
        ok = [i for i, s in enumerate(obs["sample"]) if s not in rejected]
        log2, obs = log2[:, ok], obs.iloc[ok].reset_index(drop=True)
        print(f"   excluded; {log2.shape[1]} samples remain")

    # ---- the pipeline's own normalisation, then the clock ------------------------------------
    rpm = np.power(2.0, log2) - 1.0                      # sources.py:506, verbatim
    expr = np.log1p(normalize_counts(np.clip(rpm, 0.0, None).T, target_sum=1e4))
    clock = LinearClock.from_json(CLOCK)
    age = clock.predict_age(expr, genes)
    obs["age"] = age
    present = sum(1 for g in genes if g in clock.weights)
    # NOT pre-registered. Gene COUNT understates the problem: what matters to a linear clock is
    # how much of its total |weight| is reachable, because every missing gene silently drops a
    # term from `expr @ w` while the intercept stays fixed -- a pure BIAS on absolute age. It
    # cancels in any control-relative difference, which is why M-E1 and M-E3 can disagree.
    gset = set(genes)
    w_tot = sum(abs(v) for v in clock.weights.values())
    w_present = sum(abs(v) for g, v in clock.weights.items() if g in gset)
    print("\n   normalised via 2**log2-1 -> normalize_counts(1e4) -> log1p (sources.py:506)")
    print(f"   clock genes present: {present} of {len(clock.weights)} "
          f"({100 * present / len(clock.weights):.1f}%); "
          f"of the clock's total |weight|, {100 * w_present / w_tot:.1f}% is reachable")
    print(f"   predicted age range {age.min():.1f} .. {age.max():.1f} yr")
    out["clock_genes_present"] = present
    out["clock_weight_mass_present"] = float(w_present / w_tot)

    # ---- M-E1: absolute calibration on UNTREATED samples --------------------------------------
    cv_mae = 12.26879346460328
    untreated = obs[obs.group.isin(("control", "day0"))]
    m, lo, hi, n = ci(untreated.age)
    true_mean = float(np.mean([TRUE_AGES[d] for d in untreated.donor]))
    delta = abs(m - true_mean)
    ok1 = delta <= cv_mae
    print("\n" + "-" * 92)
    print("M-E1 — DOES THE CLOCK READ ABSOLUTE AGE HERE?")
    print("-" * 92)
    rows = []
    for d in sorted(TRUE_AGES):
        sub = untreated[untreated.donor == d]
        mm, ll, hh, nn = ci(sub.age)
        rows.append([d, f"{TRUE_AGES[d]:.0f}", str(nn), f"{mm:.1f}",
                     f"[{ll:.1f},{hh:.1f}]", f"{mm - TRUE_AGES[d]:+.1f}"])
    print(render_table(["donor", "true age", "n untreated", "predicted", "95% CI", "error"],
                       rows, aligns=["l", "r", "r", "r", "r", "r"]))
    print(f"   pooled predicted {m:.1f} vs true mean {true_mean:.1f}  ->  |Δ| = {delta:.1f} yr "
          f"against one cv_mae = {cv_mae:.2f}")
    print(f"   -> {'PASS-CALIBRATION' if ok1 else 'FAIL-CALIBRATION'}")
    if not ok1:
        print("      per the pre-registration: absolute ages are INVALID here; only")
        print("      control-relative ΔAge may be used below.")
    print("   the 53-vs-38 contrast is NOT gated (15 yr against a 12.27 yr cv_mae) and is")
    print("   reported as indicative only.")

    # NOT pre-registered. `untreated` pools two very different things -- day-0 fibroblasts that
    # never saw reprogramming media, and negative controls cultured ALONGSIDE the experiment for
    # 10-17 days. Splitting them decomposes the bias into a cross-study part and a culture part,
    # which are different questions with different owners.
    sp = {}
    for g in ("day0", "control"):
        sub = obs[obs.group == g]
        tm = float(np.mean([TRUE_AGES[d] for d in sub.donor]))
        sp[g] = {"n": int(len(sub)), "pred": float(sub.age.mean()), "true": tm,
                 "bias": float(sub.age.mean() - tm)}
    print("\n   BIAS DECOMPOSED (not pre-registered):")
    print(f"      day-0 fibroblasts, never in reprogramming media: n={sp['day0']['n']}, "
          f"predicted {sp['day0']['pred']:.1f} vs true {sp['day0']['true']:.1f} "
          f"-> bias {sp['day0']['bias']:+.1f} yr")
    print(f"      negative controls, cultured alongside 10-17 d:   n={sp['control']['n']}, "
          f"predicted {sp['control']['pred']:.1f} vs true {sp['control']['true']:.1f} "
          f"-> bias {sp['control']['bias']:+.1f} yr")
    print(f"      => a ~{sp['day0']['bias']:+.0f} yr CROSS-STUDY floor even on fresh cells, plus a"
          f" further {sp['control']['bias'] - sp['day0']['bias']:+.1f} yr that tracks TIME IN")
    print("         CULTURE. The second is why a contemporaneous control matters: it drifts with")
    print("         the experiment, and measuring ΔAge against the SAME day cancels that drift.")
    out["M_E1_bias_decomposition"] = sp
    out["M_E1"] = {"pooled_pred": m, "ci": [lo, hi], "n": n, "true_mean": true_mean,
                   "abs_error": delta, "cv_mae": cv_mae,
                   "verdict": "PASS-CALIBRATION" if ok1 else "FAIL-CALIBRATION",
                   "per_donor": rows}

    # ---- M-E2: contemporaneous-control ΔAge ---------------------------------------------------
    ctrl_mean, ctrl_sd, ctrl_n = {}, {}, {}
    for (d, day), sub in obs[obs.group == "control"].groupby(["donor", "day"]):
        ctrl_mean[(d, day)] = float(sub.age.mean())
        ctrl_sd[(d, day)] = float(sub.age.std(ddof=1)) if len(sub) > 1 else float("nan")
        ctrl_n[(d, day)] = int(len(sub))
    obs["dage"] = [ctrl_mean.get((r.donor, r.day), np.nan) for r in obs.itertuples()]
    obs["dage"] = obs["age"] - obs["dage"]

    print("\n" + "-" * 92)
    print("M-E2 — CONTEMPORANEOUS-CONTROL ΔAge  (what gill_bulk could never compute)")
    print("-" * 92)
    rows = []
    for (d, day, g), sub in obs[obs.group.isin(("failed", "transient"))].groupby(
            ["donor", "day", "group"]):
        if not np.isfinite(sub.dage).any():
            continue
        rows.append([d, f"{day:.0f}", g, str(len(sub)), str(ctrl_n.get((d, day), 0)),
                     f"{sub.dage.mean():+.2f}",
                     f"{ctrl_sd.get((d, day), float('nan')):.2f}"])
    print(render_table(["donor", "day", "arm", "n", "n ctrl", "mean ΔAge (yr)", "ctrl SD (yr)"],
                       rows, aligns=["l", "r", "l", "r", "r", "r", "r"]))
    out["M_E2"] = rows

    # ---- M-E3: transient vs failed, paired on (donor, day) ------------------------------------
    print("\n" + "-" * 92)
    print("M-E3 — DOES TRANSIENT REPROGRAMMING REJUVENATE? (vs Gill 2022)")
    print("-" * 92)
    diffs, trs, drows = [], [], []
    for (d, day), sub in obs.groupby(["donor", "day"]):
        t = sub[sub.group == "transient"]["dage"]
        f = sub[sub.group == "failed"]["dage"]
        if len(t) and np.isfinite(t).any():
            trs.append(float(t.mean()))
        if len(t) and len(f) and np.isfinite(t).any() and np.isfinite(f).any():
            diffs.append(float(t.mean() - f.mean()))
            drows.append([d, f"{day:.0f}", str(len(t)), str(len(f)), f"{t.mean():+.2f}",
                          f"{f.mean():+.2f}", f"{t.mean() - f.mean():+.2f}"])
    print(render_table(["donor", "day", "n trans", "n fail", "ΔAge trans", "ΔAge fail",
                        "trans − fail"], drows, aligns=["l", "r", "r", "r", "r", "r", "r"]))
    mt, lt, ht, nt = ci(trs)
    md, ld, hd, nd = ci(diffs)
    print(f"\n   ΔAge(transient) vs its own control: mean={mt:+.2f} CI=[{lt:+.2f},{ht:+.2f}] "
          f"(n={nt} cells)")
    print(f"   paired transient − failed:            mean={md:+.2f} CI=[{ld:+.2f},{hd:+.2f}] "
          f"(n={nd} cells)")
    if np.isfinite(ht) and ht < 0:
        v3 = "REPRODUCED — ΔAge(transient) < 0, CI excludes 0"
    elif np.isfinite(lt) and lt > 0:
        v3 = "CONTRADICTED — reads as AGEING; per the pre-registration this ESCALATES"
    else:
        v3 = "NOT REPRODUCED — CI includes 0 (this is NOT evidence of absence at this n)"
    print(f"   -> {v3}")
    out["M_E3"] = {"transient_mean": mt, "transient_ci": [lt, ht], "n_cells": nt,
                   "paired_vs_failed_mean": md, "paired_ci": [ld, hd], "n_paired": nd,
                   "verdict": v3, "rows": drows}

    # ---- M-E4: how noisy is the zero-point? ---------------------------------------------------
    sds = [v for v in ctrl_sd.values() if np.isfinite(v)]
    pooled = float(np.sqrt(np.mean(np.square(sds)))) if sds else float("nan")
    print("\n" + "-" * 92)
    print("M-E4 — HOW NOISY IS THE ZERO-POINT?  (the question D2 could never answer)")
    print("-" * 92)
    print(f"   within-(donor, day) control SD, {len(sds)} groups of n=2-3: "
          f"pooled = {pooled:.2f} yr   vs cv_mae {cv_mae:.2f}")
    if pooled >= cv_mae:
        v4 = ("zero-point wobbles as much as the clock's own error -> the ±12.7 yr per-donor "
              "offset is INDISTINGUISHABLE from measurement noise; Stage 2's premise stays void")
    elif pooled < cv_mae / 2:
        v4 = ("clock is far MORE reproducible on replicates than its CV implies -> the per-donor "
              "offset is MORE LIKELY REAL BIOLOGY; first evidence for Stage 2's premise")
    else:
        v4 = "inconclusive at this n"
    print(f"   -> {v4}")
    out["M_E4"] = {"pooled_control_sd": pooled, "n_groups": len(sds), "cv_mae": cv_mae,
                   "per_group": {f"{d}_d{day:.0f}": v for (d, day), v in ctrl_sd.items()},
                   "verdict": v4}

    # ---- M-E5: the exp1/exp2 batch offset -----------------------------------------------------
    print("\n" + "-" * 92)
    print("M-E5 — exp1 / exp2 BATCH OFFSET  (D1, measured directly)")
    print("-" * 92)
    offs, orows = [], []
    for (d, day, g), sub in obs.groupby(["donor", "day", "group"]):
        e1, e2 = sub[sub.exp == "exp1"]["age"], sub[sub.exp == "exp2"]["age"]
        if len(e1) and len(e2):
            offs.append(float(e1.mean() - e2.mean()))
            orows.append([d, f"{day:.0f}", g, str(len(e1)), str(len(e2)),
                          f"{e1.mean() - e2.mean():+.2f}"])
    if orows:
        print(render_table(["donor", "day", "arm", "n exp1", "n exp2", "exp1 − exp2 (yr)"],
                           orows, aligns=["l", "r", "l", "r", "r", "r"]))
        mo, lo5, ho5, no = ci(offs)
        v5 = ("D1 CONFIRMED SEVERE — batch alone moves age by >= one clock error"
              if abs(mo) >= cv_mae else "batch present but sub-error; record as a caveat")
        print(f"   mean offset {mo:+.2f} CI=[{lo5:+.2f},{ho5:+.2f}] (n={no}) vs cv_mae "
              f"{cv_mae:.2f}  -> {v5}")
        # NOT pre-registered. The POOLED mean is the wrong summary here and hides the thing that
        # matters: an arm-dependent batch offset does NOT cancel in ΔAge = treated - control.
        by_arm = {}
        for r, o in zip(orows, offs, strict=True):
            by_arm.setdefault(r[2], []).append(o)
        print("\n   STRATIFIED BY ARM (not pre-registered — the pooled mean above hides this):")
        strat = {}
        for a, vs in sorted(by_arm.items()):
            mm, ll, hh, nn = ci(vs)
            strat[a] = {"mean": mm, "ci": [ll, hh], "n": nn}
            print(f"      {a:<10} mean {mm:+.2f} CI=[{ll:+.2f},{hh:+.2f}] (n={nn})")
        cm = strat.get("control", {}).get("mean", float("nan"))
        print(f"   ** the CONTROL arm alone shifts {cm:+.2f} yr between batches. Because ΔAge is")
        print("      measured AGAINST that control, an arm-dependent offset does NOT cancel --")
        print("      so the pooled figure understates D1 for any cross-batch comparison. **")
        out["M_E5"] = {"mean_offset": mo, "ci": [lo5, ho5], "n": no, "verdict": v5,
                       "rows": orows, "by_arm": strat}
    else:
        v5 = ("NO MATCHED CELLS — no (donor, day, arm) appears in both batches, so the design "
              "does not permit this comparison")
        print(f"   {v5}")
        out["M_E5"] = {"verdict": v5, "rows": []}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   LIMITS: bulk; 3 donors / 2 distinct ages; NO HARMONIZER so these numbers are")
    print("   directional and NOT comparable to the project's harmonised ΔAge figures.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
