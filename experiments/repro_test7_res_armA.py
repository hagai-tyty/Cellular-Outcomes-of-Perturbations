"""
REPRODUCTION TEST — do the July RES/ranking results (Test 7, 7.1, 7.2) reproduce on arm A?

WHY THIS EXISTS (user's design)
-------------------------------
Tests 7, 7.1 and 7.2 were run in July 2026 on the other machine, on the pre-Stage-1.5.6
LOOCV builds (`runs/cellfate_loocv_<donor>`). Since then the pipeline changed a lot
(C-1..C-6, the age_mask/deconfound_mask split, the C-5 age window, control re-centring)
and those builds were deleted from disk. Arm A is the current true-label build. If the
July numbers reproduce on arm A, the ranking pipeline still works as recorded; if they
do not, something between July and now moved a load-bearing result.

WHY THIS COMPARISON IS LIKE-FOR-LIKE (and the HFF trajectory one was NOT)
------------------------------------------------------------------------
The user's correction stands: `diag_pipeline_decompose` / `diag_gc_hff_signature` ran the
clock on a SINGLE raw dataset, so comparing them to a two-dataset harmonized build is
invalid. The Test 7 family is different — it runs the LOOCV donor holdout on the BUILT
artefacts, and both sides are the same two-dataset build. Verified from metadata before
running (July values read out of git at f353526^):

    field                July runs/cellfate_loocv_N2   arm A cellfate_loocv_N2_armA
    n_samples            42605                          42605
    n_shards             51                             51
    n_age_labeled        42605                          42605
    split_sizes          852/115/117/21                 852/115/117/21
    gene_panel_hash      783f269a214aa972               783f269a214aa972
    label_distribution   22635/1095/18875               22635/1095/18875

Same cells, same panel, same splits, both harmonized over HFF + the six Gill donors. The
only metadata difference is `baseline_census`, which did not exist in July (added by G-a).
So any difference in the numbers below is the MODEL/PIPELINE, not the data.

WHAT IS COMPARED
----------------
The three July scripts are imported unmodified and their `resolve_root` is redirected to
the `_armA` fold roots. Nothing in `experiments/test7*.py` is edited.

PRE-REGISTERED BEFORE LOOKING AT ANY ARM-A NUMBER
-------------------------------------------------
Prediction: the qualitative findings reproduce. Stated basis, not a peek at this test's
output: arm A's committed scorecard already carries `rank_model_dage = 0.948`, which is
the same metric family as Test 7's `model_dAge` (0.948 in July). I therefore expect P1-a
to pass comfortably. I am less sure about model_RES, whose July fold spread was 3x larger
(+/-0.091 vs +/-0.030), and about Test 7.1's precision@5, which is a 5-item statistic on a
21-cell test split and is close to pure noise.

PRIMARY bar — the load-bearing claims. ALL must hold to call it a reproduction:
  P1-a  Test 7:   model_dAge >= 0.85 AND ridge_dAge >= 0.85          (July 0.948 / 0.955)
  P1-b  Test 7:   paired (model_RES - ridge_dAge) 95% CI entirely below 0
                                                                     (July -0.268 [-0.381,-0.155])
  P1-c  Test 7.2: paired (B_true - A_true) 95% CI entirely below 0   (July -0.300 [-0.473,-0.128])

SECONDARY bar — magnitude drift. Reported, does NOT decide the verdict:
  P2-a  |armA - July| <= 0.10 for model_dAge and ridge_dAge  (~3x their July fold SD)
  P2-b  |armA - July| <= 0.15 for model_RES                  (~1.6x its July fold SD)
  P2-c  Test 7.1: RES worse than ridge_dAge on BOTH gated and penalized targets (sign only)

OUTCOMES, registered in advance:
  O1  all PRIMARY pass, all SECONDARY pass -> arm A REPRODUCES July. The ranking pipeline
      survived the 1.5.3-1.5.6 changes intact; the July RES conclusions carry forward to
      the current build unchanged.
  O2  all PRIMARY pass, some SECONDARY fail -> REPRODUCES QUALITATIVELY. Same conclusions,
      drifted magnitudes; record the drift and which metric moved. Conclusions stand.
  O3  any PRIMARY fails -> DOES NOT REPRODUCE. Something between July and now moved a
      load-bearing result. Escalate: do not treat either the July record or arm A's
      scorecard as describing the same system until it is explained.
  O4  the run cannot complete (missing artefact, incompatible bundle) -> INCONCLUSIVE,
      reported as such. NOT counted as a pass or a failure.

Note on what a pass does and does not license: this is a REPRODUCTION check, not a
validation. Reproducing Test 7 means arm A still ranks by dAge about as well as July did
and RES still degrades that ranking. It says nothing about whether the dAge labels are
CORRECT age -- that is the open question arms C and D narrowed but did not close.

READ-ONLY. Touches no build, writes only results/repro_test7_res_armA_results.json.

USAGE (repo root, venv active):
    python experiments/repro_test7_res_armA.py
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "experiments"))

ARM = "armA"
DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
OUT = REPO / "results" / "repro_test7_res_armA_results.json"

# ---------------------------------------------------------------------------
# July reference, transcribed from experiments/DELTAAGE_LAB_NOTEBOOK.md.
# Test 7    -> "### Test 7 RESULT (user ran it)"
# Test 7.1  -> "### Test 7.1 RESULT (user ran it)"
# Test 7.2  -> "### Test 7.2 RESULT (user ran it)"
# Scripts authored 2026-07-11/12 by hagai-tyty (b1c97b6, ae0dc11, 7aa4152).
# ---------------------------------------------------------------------------
JULY = {
    "test7": {
        "per_fold": {
            "N2": {"model_RES": 0.742, "model_dAge": 0.910, "ridge_dAge": 0.957},
            "N3": {"model_RES": 0.804, "model_dAge": 0.909, "ridge_dAge": 0.925},
            "O1": {"model_RES": 0.684, "model_dAge": 0.990, "ridge_dAge": 0.960},
            "O2": {"model_RES": 0.507, "model_dAge": 0.970, "ridge_dAge": 0.952},
            "Y1": {"model_RES": 0.706, "model_dAge": 0.960, "ridge_dAge": 0.951},
            "Y2": {"model_RES": 0.674, "model_dAge": 0.947, "ridge_dAge": 0.983},
        },
        "agg": {"model_RES": 0.686, "model_dAge": 0.948, "ridge_dAge": 0.955},
        "std": {"model_RES": 0.091, "model_dAge": 0.030, "ridge_dAge": 0.017},
        "paired_res_minus_ridge": {"mean": -0.268, "lo": -0.381, "hi": -0.155},
        "res_wins": 0,
    },
    "test7_1": {
        "n_unsafe": {"N2": 0, "N3": 3, "O1": 4, "O2": 5, "Y1": 8, "Y2": 5},
        "n_cells": {"N2": 21, "N3": 21, "O1": 21, "O2": 21, "Y1": 19, "Y2": 21},
        "gated": {"model_RES": -0.005, "model_dAge": 0.292, "ridge_dAge": 0.295},
        "pen": {"model_RES": 0.137, "model_dAge": 0.414, "ridge_dAge": 0.414},
        "p@5": {"model_RES": 0.20, "model_dAge": 0.27, "ridge_dAge": 0.30},
        "paired_gated": {"mean": -0.300, "lo": -0.567, "hi": -0.034},
        "paired_pen": {"mean": -0.277, "lo": -0.525, "hi": -0.030},
        # per-fold RES gated, quoted in the notebook verdict (Y2 not quoted there)
        "res_gated_per_fold": {"N2": 0.742, "N3": 0.676, "O1": -0.257, "O2": -0.431,
                               "Y1": -0.534},
    },
    "test7_2": {
        "true": {"A": 0.955, "B": 0.654, "mean": -0.300, "lo": -0.473, "hi": -0.128},
        "safe": {"A": 0.295, "B": 0.116, "mean": -0.179, "lo": -0.420, "hi": 0.061},
    },
}

BAR = {
    "P1a_min_dage_spearman": 0.85,
    "P2a_tol_dage": 0.10,
    "P2b_tol_res": 0.15,
}


def load(mod_name: str):
    """Import a July test module unmodified and redirect it at the arm-A fold roots."""
    mod = importlib.import_module(mod_name)
    mod.resolve_root = lambda name, _m=mod: str(REPO / f"{name}_{ARM}")
    return mod


def run_folds(mod) -> dict:
    out = {}
    for d in DONORS:
        r = mod.one_fold(d)
        if r is not None:
            out[d] = r
    return out


def agg(per_fold: dict, key: str) -> float:
    v = [per_fold[d][key] for d in per_fold]
    return float(np.nanmean(v)) if v else float("nan")


def std(per_fold: dict, key: str) -> float:
    v = [per_fold[d][key] for d in per_fold]
    return float(np.nanstd(v)) if v else float("nan")


def fmt_ci(m, lo, hi) -> str:
    return f"{m:+.3f} [{lo:+.3f},{hi:+.3f}]"


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\n" + "=" * 78)
    print("REPRODUCTION — July Test 7 / 7.1 / 7.2  vs  arm A")
    print("=" * 78)
    print("Both sides are the SAME two-dataset harmonized build (42605 cells, 51 shards,")
    print("panel 783f269a214aa972, splits 852/115/117/21). Differences below are the")
    print("MODEL and PIPELINE, not the data.")

    missing = [d for d in DONORS if not (REPO / f"cellfate_loocv_{d}_{ARM}").exists()]
    if missing:
        print(f"\n   INCONCLUSIVE (O4): missing arm-A fold roots: {missing}")
        return 2

    t7 = load("test7_ranking")
    t71 = load("test7_1_safe_ranking")
    t72 = load("test7_2_res_isolation")

    checks: list[tuple[str, bool, str]] = []
    report: dict = {"arm": ARM, "bar": BAR, "july_source": "DELTAAGE_LAB_NOTEBOOK.md"}

    # ---------------------------------------------------------------- Test 7
    print("\n\n--- TEST 7 — RES ranking vs ranking-by-ridge-ΔAge (Spearman vs true ΔAge) ---")
    names = ["model_RES", "model_dAge", "ridge_dAge"]
    pf7 = run_folds(t7)
    rows = []
    for d in DONORS:
        if d in pf7:
            rows.append([d] + [f"{pf7[d][n]:+.3f}" for n in names]
                        + [f"{JULY['test7']['per_fold'][d][n]:+.3f}" for n in names])
        else:
            rows.append([d] + ["n/a"] * 6)
    print("\n" + render_table(
        ["fold", "A:RES", "A:mdAge", "A:ridge", "J:RES", "J:mdAge", "J:ridge"],
        rows, aligns=["l"] + ["r"] * 6))

    a7 = {n: agg(pf7, n) for n in names}
    s7 = {n: std(pf7, n) for n in names}
    print("\n   arm A aggregate:  " + "   ".join(f"{n}={a7[n]:.3f}±{s7[n]:.3f}" for n in names))
    print("   July  aggregate:  " + "   ".join(
        f"{n}={JULY['test7']['agg'][n]:.3f}±{JULY['test7']['std'][n]:.3f}" for n in names))
    print("   delta (armA−July):" + "   ".join(
        f"  {n}={a7[n] - JULY['test7']['agg'][n]:+.3f}" for n in names))

    diffs7 = [pf7[d]["model_RES"] - pf7[d]["ridge_dAge"] for d in pf7]
    m7, (lo7, hi7) = t7.paired_ci(diffs7)
    wins7 = sum(1 for x in diffs7 if x > 0)
    j7 = JULY["test7"]["paired_res_minus_ridge"]
    print(f"\n   paired (model_RES − ridge_dAge)  arm A {fmt_ci(m7, lo7, hi7)}"
          f"   July {fmt_ci(j7['mean'], j7['lo'], j7['hi'])}")
    print(f"   model_RES ranks better on {wins7}/{len(diffs7)} folds"
          f"   (July {JULY['test7']['res_wins']}/6)")

    checks.append(("P1-a  model_dAge & ridge_dAge >= 0.85",
                   a7["model_dAge"] >= BAR["P1a_min_dage_spearman"]
                   and a7["ridge_dAge"] >= BAR["P1a_min_dage_spearman"],
                   f"model_dAge={a7['model_dAge']:.3f}  ridge_dAge={a7['ridge_dAge']:.3f}"))
    checks.append(("P1-b  (RES − ridge) 95% CI entirely below 0",
                   bool(hi7 < 0), fmt_ci(m7, lo7, hi7)))
    for n, tol in (("model_dAge", BAR["P2a_tol_dage"]), ("ridge_dAge", BAR["P2a_tol_dage"]),
                   ("model_RES", BAR["P2b_tol_res"])):
        dd = a7[n] - JULY["test7"]["agg"][n]
        tag = "P2-a" if n != "model_RES" else "P2-b"
        checks.append((f"{tag}  |Δ{n}| <= {tol:.2f}", bool(abs(dd) <= tol),
                       f"Δ={dd:+.3f}  (armA {a7[n]:.3f} vs July {JULY['test7']['agg'][n]:.3f})"))

    report["test7"] = {"per_fold": pf7, "agg": a7, "std": s7,
                       "paired_res_minus_ridge": {"mean": m7, "lo": lo7, "hi": hi7},
                       "res_wins": wins7}

    # -------------------------------------------------------------- Test 7.1
    print("\n\n--- TEST 7.1 — ranked against SAFE REJUVENATION (RES's own objective) ---")
    pf71 = run_folds(t71)
    comp = [[d, str(pf71[d]["n"]), str(pf71[d]["n_unsafe"]),
             str(JULY["test7_1"]["n_cells"].get(d, "?")),
             str(JULY["test7_1"]["n_unsafe"].get(d, "?"))] for d in pf71]
    print("\n" + render_table(["fold", "A:cells", "A:unsafe", "J:cells", "J:unsafe"],
                              comp, aligns=["l", "r", "r", "r", "r"]))

    t71_out = {}
    for tgt, label in [("gated", "GATED (unsafe ranked last)"),
                       ("pen", "PENALIZED (soft unsafe penalty)")]:
        rows = [[d] + [f"{pf71[d][f'{n}|{tgt}']:+.3f}" for n in names] for d in pf71]
        print(f"\n  target = {label}")
        print(render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
        a = {n: agg(pf71, f"{n}|{tgt}") for n in names}
        jkey = "gated" if tgt == "gated" else "pen"
        print("   arm A aggregate: " + "   ".join(f"{n}={a[n]:+.3f}" for n in names))
        print("   July  aggregate: " + "   ".join(
            f"{n}={JULY['test7_1'][jkey][n]:+.3f}" for n in names))
        diffs = [pf71[d][f"model_RES|{tgt}"] - pf71[d][f"ridge_dAge|{tgt}"] for d in pf71]
        m, (lo, hi), nn = t71.paired_ci(diffs)
        jp = JULY["test7_1"][f"paired_{jkey}"]
        print(f"   paired (RES − ridge)  arm A {fmt_ci(m, lo, hi)} (n={nn})"
              f"   July {fmt_ci(jp['mean'], jp['lo'], jp['hi'])}")
        t71_out[tgt] = {"agg": a, "paired": {"mean": m, "lo": lo, "hi": hi, "n": nn}}
        checks.append((f"P2-c  {tgt}: RES worse than ridge_dAge (sign)",
                       bool(a["model_RES"] < a["ridge_dAge"]),
                       f"RES={a['model_RES']:+.3f} < ridge={a['ridge_dAge']:+.3f}"))

    K = t71.K_TOP
    rows = [[d] + [f"{pf71[d][f'{n}|p@{K}']:.2f}" for n in names] for d in pf71]
    print(f"\n  precision@{K} (top-{K} truly SAFE & rejuvenating)")
    print(render_table(["fold"] + names, rows, aligns=["l", "r", "r", "r"]))
    ap = {n: agg(pf71, f"{n}|p@{K}") for n in names}
    print("   arm A aggregate: " + "   ".join(f"{n}={ap[n]:.2f}" for n in names))
    print("   July  aggregate: " + "   ".join(
        f"{n}={JULY['test7_1']['p@5'][n]:.2f}" for n in names))
    t71_out["p_at_k"] = {"k": K, "agg": ap}

    print("\n   RES(gated) vs unsafe-cell count — July's pinpointing pattern"
          " (RES degrades as unsafe count rises):")
    pat = [[d, str(pf71[d]["n_unsafe"]), f"{pf71[d]['model_RES|gated']:+.3f}",
            (f"{JULY['test7_1']['res_gated_per_fold'][d]:+.3f}"
             if d in JULY["test7_1"]["res_gated_per_fold"] else "not quoted")]
           for d in sorted(pf71, key=lambda x: pf71[x]["n_unsafe"])]
    print(render_table(["fold", "unsafe", "armA RES|gated", "July RES|gated"],
                       pat, aligns=["l", "r", "r", "r"]))
    xs = [pf71[d]["n_unsafe"] for d in pf71]
    ys = [pf71[d]["model_RES|gated"] for d in pf71]
    rho_pat = float(np.corrcoef(xs, ys)[0, 1]) if len(set(xs)) > 1 else float("nan")
    print(f"   arm A corr(unsafe count, RES|gated) = {rho_pat:+.3f}"
          "   (July: strongly negative, the basis of the 'fate predictions are wrong' claim)")
    t71_out["unsafe_vs_res_corr"] = rho_pat
    t71_out["per_fold"] = pf71
    report["test7_1"] = t71_out

    # -------------------------------------------------------------- Test 7.2
    print("\n\n--- TEST 7.2 — RES FORMULA isolated (same ridge ΔAge fed to A and B) ---")
    pf72 = run_folds(t72)
    t72_out = {}
    for tgt, label in [("true", "vs TRUE ΔAge"), ("safe", "vs SAFE-REJUVENATION (gated)")]:
        rows = [[d, f"{pf72[d][f'A_{tgt}']:+.3f}", f"{pf72[d][f'B_{tgt}']:+.3f}",
                 f"{pf72[d][f'B_{tgt}'] - pf72[d][f'A_{tgt}']:+.3f}"] for d in pf72]
        print(f"\n  {label}")
        print(render_table(["fold", "A = ridge ΔAge", "B = RES(ridge ΔAge)", "B − A"],
                           rows, aligns=["l", "r", "r", "r"]))
        aA = agg(pf72, f"A_{tgt}")
        aB = agg(pf72, f"B_{tgt}")
        diffs = [pf72[d][f"B_{tgt}"] - pf72[d][f"A_{tgt}"] for d in pf72]
        m, (lo, hi), nn = t72.paired_ci(diffs)
        j = JULY["test7_2"][tgt]
        print(f"   arm A: A={aA:.3f}  B={aB:.3f}   paired(B−A) {fmt_ci(m, lo, hi)} (n={nn})")
        print(f"   July : A={j['A']:.3f}  B={j['B']:.3f}   paired(B−A) "
              f"{fmt_ci(j['mean'], j['lo'], j['hi'])}")
        t72_out[tgt] = {"A": aA, "B": aB, "paired": {"mean": m, "lo": lo, "hi": hi, "n": nn}}
        if tgt == "true":
            checks.append(("P1-c  Test 7.2 (B−A) 95% CI entirely below 0",
                           bool(hi < 0), fmt_ci(m, lo, hi)))
    t72_out["per_fold"] = pf72
    report["test7_2"] = t72_out

    # ------------------------------------------------------------- the verdict
    print("\n\n" + "=" * 78)
    print("PRE-REGISTERED BAR")
    print("=" * 78)
    rows = [[nm, "PASS" if ok else "FAIL", why] for nm, ok, why in checks]
    print(render_table(["check", "result", "observed"], rows, aligns=["l", "l", "l"]))

    primary = [c for c in checks if c[0].startswith("P1")]
    secondary = [c for c in checks if c[0].startswith("P2")]
    p_ok = all(c[1] for c in primary)
    s_ok = all(c[1] for c in secondary)

    if p_ok and s_ok:
        outcome, verdict = "O1", "REPRODUCES"
    elif p_ok:
        outcome, verdict = "O2", "REPRODUCES QUALITATIVELY"
    else:
        outcome, verdict = "O3", "DOES NOT REPRODUCE"

    print(f"\n   PRIMARY {sum(c[1] for c in primary)}/{len(primary)} pass"
          f"   SECONDARY {sum(c[1] for c in secondary)}/{len(secondary)} pass")
    print(f"\n   OUTCOME {outcome} — {verdict}")
    if outcome == "O1":
        print("   -> arm A reproduces the July RES/ranking results on identical data.")
        print("      The ranking pipeline survived the 1.5.3–1.5.6 changes intact and the")
        print("      July conclusions carry forward to the current build.")
    elif outcome == "O2":
        print("   -> same conclusions, drifted magnitudes. The failing SECONDARY rows above")
        print("      say which metric moved and by how much. Conclusions stand; record drift.")
    else:
        print("   -> ESCALATE. A load-bearing result moved between July and arm A. Do not")
        print("      treat the July record and arm A's scorecard as the same system until")
        print("      this is explained.")
    print("\n   SCOPE: this is REPRODUCTION, not validation. It shows arm A still ranks by")
    print("   ΔAge as well as July did and that RES still degrades it. It says NOTHING about")
    print("   whether the ΔAge labels are correct age — still open after arms C and D.")

    report["checks"] = [{"name": n, "pass": bool(o), "observed": w} for n, o, w in checks]
    report["primary_pass"] = bool(p_ok)
    report["secondary_pass"] = bool(s_ok)
    report["outcome"] = outcome
    report["verdict"] = verdict
    report["july"] = JULY
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0 if p_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
