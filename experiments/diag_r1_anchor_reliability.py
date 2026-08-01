"""STAGE 1.5.2 §11 / §9-R1 / §9-R4 — is the METHYLATION ANCHOR itself reliable on this data?

    python experiments/diag_r1_anchor_reliability.py                # pre-register only (no data read)
    python experiments/diag_r1_anchor_reliability.py --run          # + the measurement

READ-ONLY. Writes `diag_r1_anchor_reliability_results.json`. `src/` untouched.

WHY THIS EXISTS, AND WHY IT RUNS *AFTER* M-2a BUT *BEFORE* M-2a IS BELIEVED
--------------------------------------------------------------------------
STAGE_1_5_2 §11 states the falsification condition for its own negative result:

    "A negative verdict is falsified if the methylation ages themselves are unreliable here — so
     the clocks are checked against donor chronological age on the CD13 arm (R1) before any
     negative verdict is accepted."

M-2a returned SPLIT ⇒ NOT CALIBRATABLE on 2026-07-31 **without that check having been run.** The
verdict is therefore recorded but not yet *accepted*. This script is the missing precondition. It
can only do one of two things:

  * confirm the anchor reads age here  -> M-2a's negative verdict STANDS;
  * show the anchor does not read age  -> M-2a's negative verdict is WITHDRAWN, because a
    disagreement between two instruments says nothing about which one is wrong.

WHAT IS MEASURED, AND WHY EACH FORM WAS CHOSEN
----------------------------------------------
R1a  LODO chronological-age recovery.  `diag_methylation_anchor.py`'s G2 derives the intercept from
     the SAME day-0 samples it then grades, which its own docstring flags as partly self-fulfilling.
     Here the intercept is derived from the OTHER donors and used to predict the held-out one, so
     nothing the fold is graded on contributed to its own intercept. 3 donors -> 3 folds.

R1b  The age gap, INTERCEPT-FREE.  Donors are O1 = 53, O2 = 53, O3 = 38, so the data contains one
     real 15-yr contrast. `anti_trafo` is linear above adult_age, so `age_i − age_j = 21·(lp_i − lp_j)`
     and the intercept cancels exactly (REV FINAL §4.3, applied rather than assumed). Measured on
     the **untreated `Negative control fibroblast` arm** (n=7 per donor, 21 samples) — cells that
     never received OSKM, which is the closest thing this series has to a clean age readout.

R1c  Drift of the treated NON-RESPONDER arms.  §9-R1's own words: "CD13's absolute methylation ages
     are compared to donor chronological age — a direct check of drift that costs nothing." GSE165179
     has no CD13; its equivalents are `Failed to transiently reprogram fibroblast` and its
     intermediate. Reported as a WITHIN-DONOR difference against that donor's own day-0 fibroblast,
     which is again intercept-free.

R1d  Inter-clock agreement — THE ONLY WELL-POWERED CHECK HERE, and the decisive one.
     R1a and R1b are what §9-R1 and §11 literally ask for, but this series has **3 donors carrying
     only 2 distinct ages**, so both come back UNRESOLVABLE at their proposed bars and have to be
     loosened. Loosening a bar that gates a negative verdict makes that verdict HARDER to falsify —
     a bias in its own favour, so neither can carry the decision alone.

     R1d has no such problem. It asks the same question — "are the methylation readings reliable on
     these samples?" — as: do the **two independent Horvath clocks agree with each other**, over the
     same 68 conditions, under the same pluripotency partialling, against **the same ρ_partial ≥ 0.50
     bar M-2a was graded on** (frozen 2026-07-31, RESOLVABLE at 99.4%). No new bar is introduced.

     That symmetry is what makes it decisive: RNA↔methylation and methylation↔methylation are scored
     by an identical criterion on identical samples, so if one passes and the other fails, the
     failure is localised to an instrument rather than to the data. Reported with the CpG overlap
     between the two clocks, because agreement driven by shared probes would be partly trivial.

R4   CpG coverage per clock (§9-R4): below 90% the clock is reported as degraded, matching 1.5.1.

THE NULL, AND WHY IT IS DELIBERATELY PESSIMISTIC
------------------------------------------------
"A system that meets the intent exactly" = a methylation clock with the published Horvath accuracy,
MAE ≈ 3.0 yr. Error is modelled as a **donor-level offset with no averaging benefit**: with 7
samples per donor an independent-error model would shrink the SE by √7 and make every bar look
easily resolvable. Donor-level clock error does not shrink that way, so the conservative model is
used and the bars are checked against it. If a bar is RESOLVABLE under this null it is resolvable
in fact.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from audit_metrics import bar_verdict  # noqa: E402

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


N_SIM = 20000
RNG = np.random.default_rng(0)

# Published Horvath accuracy. MAE 3.0 yr for a half-normal => sigma = MAE / sqrt(2/pi).
PUBLISHED_MAE = 3.0
SIGMA = PUBLISHED_MAE / np.sqrt(2.0 / np.pi)

# R1a's bar is NOT invented here: it is G2_MAE_TOL = 5.0 yr, already registered in
# STAGE_1_5_1_REV_FINAL §3 and used by diag_methylation_anchor.py. Reusing it means this check
# cannot be tuned to its own answer.
R1A_MAE_BAR = 5.0
# R1b: two clock readings differ with ~sqrt(2) the single-reading error, so the tolerance on a
# DIFFERENCE is sqrt(2) x the tolerance on a level. Proposed, then checked below.
R1B_GAP_BAR = R1A_MAE_BAR * np.sqrt(2.0)
TRUE_GAP_YEARS = 15.0            # O1/O2 = 53, O3 = 38
COVERAGE_BAR = 0.90              # §9-R4
# R1d reuses M-2a's OWN frozen bar verbatim — not a new one. §6, 2026-07-31: rho_partial >= 0.50 at
# n=68, RESOLVABLE at 99.4%. Scoring meth<->meth by the identical criterion that scored RNA<->meth is
# the whole point; a bar invented for this check would destroy the comparison.
R1D_RHO_BAR = 0.50

ADULT_SLOPE = 21.0               # anti_trafo slope above adult_age; 1 lp unit = 21 yr

NC = "Negative control fibroblast"
FIB = "Fibroblast"
NON_RESPONDER_ARMS = ["Failed to transiently reprogram fibroblast",
                      "Failing to transiently reprogram intermediate"]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def lodo_age_errors(lp_by_donor: dict[str, float], age_by_donor: dict[str, float],
                    trafo, anti_trafo) -> dict:
    """Leave-one-donor-out chronological-age recovery. Pure.

    For each held-out donor: `k = mean(trafo(age_j) - lp_j)` over the OTHER donors, then
    `pred = anti_trafo(lp_i + k)`. The held-out donor contributes nothing to its own intercept,
    which is exactly what `diag_methylation_anchor.py`'s G2 could not claim.
    """
    donors = sorted(lp_by_donor)
    folds = []
    for d in donors:
        others = [o for o in donors if o != d]
        if len(others) < 1:
            continue
        k = float(np.mean([trafo(age_by_donor[o]) - lp_by_donor[o] for o in others]))
        pred = float(anti_trafo(lp_by_donor[d] + k))
        folds.append({"held_out": d, "true_age": age_by_donor[d], "pred_age": pred,
                      "abs_err": abs(pred - age_by_donor[d]), "intercept_from": others,
                      "intercept": k})
    mae = float(np.mean([f["abs_err"] for f in folds])) if folds else float("nan")
    return {"n_folds": len(folds), "mae_years": mae, "folds": folds}


def intercept_free_gap(lp_old: list[float], lp_young: list[float],
                       slope: float = ADULT_SLOPE) -> float:
    """Recovered age gap in years from linear predictors alone. No intercept anywhere. Pure.

    `age_i − age_j = slope·(lp_i − lp_j)` wherever `anti_trafo` is linear, so this is algebra, not
    an approximation.
    """
    return float((np.mean(lp_old) - np.mean(lp_young)) * slope)


def sim_lodo_mae(n_donors: int, sigma: float = SIGMA, n_sim: int = N_SIM) -> np.ndarray:
    """LODO MAE for a clock with donor-level error `sigma` and no averaging benefit."""
    e = RNG.normal(0.0, sigma, size=(n_sim, n_donors))
    # held-out error minus the mean error of the other donors (which set the intercept)
    tot = e.sum(axis=1, keepdims=True)
    others_mean = (tot - e) / (n_donors - 1)
    return np.abs(e - others_mean).mean(axis=1)


def sim_gap_abs_err(n_old: int, n_young: int, sigma: float = SIGMA,
                    n_sim: int = N_SIM) -> np.ndarray:
    """|recovered gap − true gap| for a correct clock, donor-level error, no averaging benefit."""
    old = RNG.normal(0.0, sigma, size=(n_sim, n_old)).mean(axis=1)
    young = RNG.normal(0.0, sigma, size=(n_sim, n_young)).mean(axis=1)
    return np.abs(old - young)


def r1_decide(r1a: dict, r1b: dict, coverage_ok: bool) -> dict:
    """Does M-2a's negative verdict survive on the per-clock checks? Pure.

    R1d is deliberately NOT folded in here: it is a property of the clock PAIR, not of either clock
    alone, so it is decided once in `overall_decision` rather than double-counted per clock.
    """
    reasons = []
    ok = True
    if not coverage_ok:
        ok = False
        reasons.append("CpG coverage below 90% on at least one clock (§9-R4: report as degraded)")
    if r1a["verdict"] != "PASS":
        ok = False
        reasons.append(f"LODO chronological-age recovery FAILED ({r1a['detail']})")
    if r1b["verdict"] != "PASS":
        ok = False
        reasons.append(f"intercept-free 15-yr gap NOT recovered ({r1b['detail']})")
    if ok:
        return {"action": "ANCHOR_READS_AGE",
                "reason": "this clock recovers chronological age on these samples within the "
                          "(loosened) bars, so it is not disqualified as an arbiter."}
    return {"action": "ANCHOR_QUESTIONABLE", "reason": "; ".join(reasons)}


def overall_decision(per_clock: dict, r1d: dict) -> dict:
    """Does M-2a's negative verdict survive? Pure.

    Weighting is decided here rather than left implicit, because the checks differ enormously in
    power. R1d is the only one that is RESOLVABLE at its bar; R1a/R1b run on 3 donors carrying 2
    distinct ages and had to be loosened to run at all. So:

      * R1d FAILS  -> the methylation readings do not even agree with each other. The anchor is not
        usable and M-2a's negative verdict is WITHDRAWN, whatever R1a/R1b say.
      * R1d PASSES and R1a/R1b pass -> the anchor is corroborated from two directions. STANDS.
      * R1d PASSES but R1a or R1b fails -> STANDS, but flagged: the well-powered check supports the
        anchor while the under-powered ones do not, and that tension is reported rather than
        resolved by picking the convenient one.
    """
    weak_ok = all(v["decision"]["action"] == "ANCHOR_READS_AGE" for v in per_clock.values())
    if r1d["verdict"] != "PASS":
        return {"action": "M2A_NEGATIVE_VERDICT_WITHDRAWN",
                "reason": "the two methylation clocks do not agree with each other on these samples "
                          f"({r1d['detail']}), scored by the same criterion M-2a used. The anchor "
                          "cannot arbitrate the RNA clock if it cannot arbitrate itself, so the "
                          "NOT CALIBRATABLE verdict is not evidence about the RNA clock."}
    if weak_ok:
        return {"action": "M2A_NEGATIVE_VERDICT_STANDS",
                "reason": "the methylation clocks agree with each other under the identical "
                          f"criterion that RNA failed ({r1d['detail']}), and each also recovers "
                          "chronological age. §11's falsification condition is NOT met: M-2a's "
                          "NOT CALIBRATABLE verdict is ACCEPTED."}
    return {"action": "M2A_NEGATIVE_VERDICT_STANDS_WITH_CAVEAT",
            "reason": "the well-powered check (R1d) supports the anchor "
                      f"({r1d['detail']}), but at least one under-powered chronological-age check "
                      "did not clear its loosened bar. The verdict stands on R1d; the tension is "
                      "recorded, not resolved."}


# --------------------------------------------------------------------------- #
# Phase 1 — pre-registration. No measurement is read.                          #
# --------------------------------------------------------------------------- #
def preregister() -> dict:
    print("STAGE 1.5.2 §11 — R1/R4 anchor-reliability check\n")
    print("  PHASE 1: freeze the bars. No beta value is read in this phase.\n")
    print(f"  null: a clock with the published Horvath accuracy, MAE {PUBLISHED_MAE:.1f} yr "
          f"(sigma {SIGMA:.2f}),")
    print("        error modelled as a DONOR-LEVEL offset with no averaging benefit (pessimistic).\n")
    print(f"  {'bar':<52}{'pass rate':>11}  verdict")
    print("  " + "-" * 86)
    out: dict = {"sigma": float(SIGMA), "published_mae": PUBLISHED_MAE, "checks": {}}

    def report(name, sim, bar, lower=True):
        v = bar_verdict(np.asarray(sim, float), bar, lower_is_better=lower)
        out["checks"][name] = {**v, "bar": bar}
        tail = "" if v["verdict"] == "RESOLVABLE" else f"   -> usable_bar {v['usable_bar']:.2f}"
        print(f"  {name:<52}{v['pass_rate']*100:>10.1f}%  {v['verdict']}{tail}")
        return v

    report(f"R1a LODO age MAE <= {R1A_MAE_BAR:.1f} yr (3 donors)",
           sim_lodo_mae(3), R1A_MAE_BAR)
    report(f"R1b |gap - 15| <= {R1B_GAP_BAR:.2f} yr (2 old vs 1 young)",
           sim_gap_abs_err(2, 1), float(R1B_GAP_BAR))
    print(f"\n  R1d rho_partial >= {R1D_RHO_BAR:.2f} at n=68 is NOT re-simulated here: it is M-2a's")
    print("      own bar, frozen 2026-07-31 at 99.4% RESOLVABLE by stage_1_5_2_resolvability.py.")
    print("      Re-deriving it would risk drift; reusing it verbatim is what makes the")
    print("      RNA<->meth / meth<->meth comparison a like-for-like one.")
    out["checks"][f"R1d rho_partial >= {R1D_RHO_BAR:.2f} (n=68)"] = {
        "bar": R1D_RHO_BAR, "bar_used": R1D_RHO_BAR, "verdict": "RESOLVABLE",
        "pass_rate": 0.994, "source": "stage_1_5_2_resolvability_results.json (M-2a rho_partial, "
                                      "ACTUAL n=68) — reused verbatim, not re-derived"}
    return out


# --------------------------------------------------------------------------- #
# Phase 2 — the measurement                                                    #
# --------------------------------------------------------------------------- #
def measure(pre: dict, meth_dir: Path, rna_dir: Path) -> dict:
    dma = _load("diag_methylation_anchor", ROOT / "experiments" / "diag_methylation_anchor.py")
    meta = dma.load_series(meth_dir / "GSE165179_series_matrix.txt.gz")
    bpath = meth_dir / "GSE165179_Matrix_processed_transient.txt.gz"

    print("\n  PHASE 2: measurement.\n")
    results: dict = {"clocks": {}}
    lp_by_clock: dict[str, dict[str, float]] = {}
    cpgs_by_clock: dict[str, set[str]] = {}
    for cfile, cname in dma.CLOCKS:
        clock = json.loads((ROOT / "configs" / "clocks" / f"{cfile}.json").read_text(
            encoding="utf-8"))
        W = {k: float(v) for k, v in clock["weights"].items()}
        samples, betas = dma.load_betas(bpath, set(W))
        present = [s for s in samples if s in meta]
        lp, cov = {}, []
        for s in present:
            v, n = dma.linear_predictor(betas[s], W)
            lp[s] = v
            cov.append(n)
        cov_frac = float(np.mean(cov)) / max(len(W), 1)
        lp_by_clock[cname] = dict(lp)
        cpgs_by_clock[cname] = set(W)

        by_arm_donor: dict[tuple[str, str], list[float]] = {}
        age_of: dict[str, float] = {}
        for s in present:
            m = meta[s]
            by_arm_donor.setdefault((m["ctype"], m["donor"]), []).append(lp[s])
            age_of[m["donor"]] = m["age"]

        # ---- R1a: LODO on the day-0 fibroblasts -------------------------------------- #
        d0_lp = {d: float(np.mean(v)) for (ct, d), v in by_arm_donor.items() if ct == FIB}
        lodo = lodo_age_errors(d0_lp, {d: age_of[d] for d in d0_lp}, dma.trafo, dma.anti_trafo)
        # Grade against the bar FROZEN IN PHASE 1 (7.17 after the §5b move), not the proposed 5.0.
        bar_a = float(pre["checks"][f"R1a LODO age MAE <= {R1A_MAE_BAR:.1f} yr (3 donors)"]
                      ["bar_used"])
        r1a = {"verdict": "PASS" if lodo["mae_years"] <= bar_a else "FAIL", "bar": bar_a,
               "bar_proposed": R1A_MAE_BAR,
               "detail": f"MAE {lodo['mae_years']:.2f} yr vs bar {bar_a:.2f} "
                         f"(proposed {R1A_MAE_BAR:.1f}, moved by §5b before the run)", **lodo}

        # ---- R1b: intercept-free gap on the UNTREATED negative-control arm ------------ #
        nc = {d: v for (ct, d), v in by_arm_donor.items() if ct == NC}
        old = [x for d, v in nc.items() if age_of[d] == 53.0 for x in v]
        young = [x for d, v in nc.items() if age_of[d] == 38.0 for x in v]
        gap = intercept_free_gap(old, young)
        gap_err = abs(gap - TRUE_GAP_YEARS)
        bar_b = float(pre["checks"][f"R1b |gap - 15| <= {R1B_GAP_BAR:.2f} yr (2 old vs 1 young)"]
                      ["bar_used"])
        r1b = {"verdict": "PASS" if gap_err <= bar_b else "FAIL",
               "recovered_gap_years": gap, "true_gap_years": TRUE_GAP_YEARS,
               "abs_err_years": gap_err, "bar": bar_b, "n_old": len(old), "n_young": len(young),
               "detail": f"recovered {gap:+.2f} yr vs true {TRUE_GAP_YEARS:+.1f}, "
                         f"|err| {gap_err:.2f} vs bar {bar_b:.2f}"}

        # ---- R1c: drift of the treated non-responders, within donor (intercept-free) --- #
        drift = {}
        for arm in NON_RESPONDER_ARMS:
            per_donor = {}
            for d, base in d0_lp.items():
                vals = by_arm_donor.get((arm, d))
                if vals:
                    per_donor[d] = float((np.mean(vals) - base) * ADULT_SLOPE)
            if per_donor:
                drift[arm] = {"per_donor_years": per_donor,
                              "mean_years": float(np.mean(list(per_donor.values())))}

        cov_ok = cov_frac >= COVERAGE_BAR
        results["clocks"][cname] = {
            "n_cpg": len(W), "coverage_frac": cov_frac,
            "coverage_verdict": "OK" if cov_ok else "DEGRADED",
            "R1a_lodo": r1a, "R1b_gap": r1b, "R1c_drift_vs_own_day0": drift,
            "decision": r1_decide(r1a, r1b, cov_ok)}

        print(f"  === {cname} ===  {len(W)} CpGs, coverage {cov_frac:.1%} "
              f"[{'OK' if cov_ok else 'DEGRADED'}]")
        print(f"     R1a  LODO age recovery   MAE {lodo['mae_years']:.2f} yr "
              f"(bar <= {bar_a:.2f}, proposed {R1A_MAE_BAR:.1f})  -> {r1a['verdict']}")
        for f in lodo["folds"]:
            print(f"            hold out {f['held_out']}: true {f['true_age']:.0f}  "
                  f"pred {f['pred_age']:.1f}  err {f['pred_age']-f['true_age']:+.1f}")
        print(f"     R1b  intercept-free gap  {gap:+.2f} yr vs true {TRUE_GAP_YEARS:+.1f} "
              f"(|err| {gap_err:.2f} <= {bar_b:.2f})  -> {r1b['verdict']}")
        for arm, v in drift.items():
            print(f"     R1c  drift {arm:<46} {v['mean_years']:+6.2f} yr vs own day-0")
        print(f"     ==> {results['clocks'][cname]['decision']['action']}\n")

    # ---- R1d: do the two methylation clocks agree, under M-2a's OWN criterion? --------- #
    m2a = _load("diag_m2a_calibratability", ROOT / "experiments" / "diag_m2a_calibratability.py")
    _ages, plu, n_plu = m2a.rna_ages_and_pluripotency(rna_dir)
    names = [c[1] for c in dma.CLOCKS]
    cond: dict[tuple, dict[str, list[float]]] = {}
    for t in lp_by_clock[names[0]]:
        if t in plu and t in lp_by_clock[names[1]] and (p := m2a.parse_title(t)):
            k = (p["donor"], p["arm"], p["day"])
            c = cond.setdefault(k, {"a": [], "b": [], "plu": []})
            c["a"].append(lp_by_clock[names[0]][t])
            c["b"].append(lp_by_clock[names[1]][t])
            c["plu"].append(plu[t])
    A = [float(np.mean(v["a"])) for v in cond.values()]
    B = [float(np.mean(v["b"])) for v in cond.values()]
    P = [float(np.mean(v["plu"])) for v in cond.values()]
    rho_all = m2a.spearman(A, B)
    rho_par = m2a.partial_spearman(A, B, P)
    shared = cpgs_by_clock[names[0]] & cpgs_by_clock[names[1]]
    bar_d = float(pre["checks"][f"R1d rho_partial >= {R1D_RHO_BAR:.2f} (n=68)"]["bar_used"])
    r1d = {"verdict": "PASS" if rho_par >= bar_d else "FAIL",
           "rho_all": rho_all, "rho_partial": rho_par, "bar": bar_d,
           "n_conditions": len(cond), "n_plu_genes": n_plu,
           "shared_cpgs": len(shared),
           "shared_frac_of_smaller": len(shared) / min(len(cpgs_by_clock[names[0]]),
                                                       len(cpgs_by_clock[names[1]])),
           "detail": f"rho_partial {rho_par:+.3f} vs bar {bar_d:.2f} at n={len(cond)}"}
    results["R1d_interclock"] = r1d
    print(f"  === R1d  inter-clock agreement ===   n = {len(cond)} conditions "
          f"(same set M-2a used)")
    print(f"     shared CpGs between the two clocks: {len(shared)} "
          f"({r1d['shared_frac_of_smaller']:.1%} of the smaller panel)")
    print(f"     rho_all      {rho_all:+.3f}   [descriptive]")
    print(f"     rho_partial  {rho_par:+.3f}   [same bar M-2a used, >= {bar_d:.2f}]  "
          f"-> {r1d['verdict']}")

    # ---- The empirical CEILING. Descriptive, UNREGISTERED, and reported anyway. -------- #
    # M-2a's bars were simulated against a null with rho_true = 0.70. R1d measures what two clocks
    # of the SAME modality actually achieve on these samples under the same partialling. If that
    # ceiling sits near the bar, then M-2a's criterion was near the limit of what ANY instrument
    # could reach here, and the SPLIT verdict says less about RNA than it appears to. This does NOT
    # change the pre-committed decision -- it qualifies what the decision means.
    m2a_path = ROOT / "results" / "diag_m2a_calibratability_results.json"
    ceiling = {"note": "descriptive, not pre-registered; does not alter any verdict"}
    if m2a_path.exists():
        prev = json.loads(m2a_path.read_text(encoding="utf-8"))
        ceiling["meth_vs_meth_rho_partial"] = rho_par
        ceiling["rna_vs_meth_rho_partial"] = {
            k: v["rho_partial"] for k, v in prev.get("clocks", {}).items()}
        ceiling["frac_of_ceiling"] = {
            k: (v / rho_par if rho_par else float("nan"))
            for k, v in ceiling["rna_vs_meth_rho_partial"].items()}
        print("\n  --- empirical ceiling (descriptive, NOT a registered criterion) ---")
        print(f"     two methylation clocks agree with each other at rho_partial {rho_par:+.3f}")
        for k, v in ceiling["rna_vs_meth_rho_partial"].items():
            print(f"     RNA vs {k:<28} {v:+.3f}   = {v/rho_par:5.1%} of that ceiling")
        print("     M-2a's bars were simulated against a null with rho_true = 0.70; nothing on")
        print("     this data reaches 0.70, including methylation against methylation.")
    results["ceiling"] = ceiling

    dec = overall_decision(results["clocks"], r1d)
    results["overall"] = dec["action"]
    results["overall_reason"] = dec["reason"]
    print(f"\n  ==> OVERALL: {dec['action']}\n      {dec['reason']}")
    return results


def main() -> int:
    pre = preregister()
    # record the frozen bar next to each check so phase 2 cannot silently use a different one
    for name, c in pre["checks"].items():
        c["bar_used"] = c["bar"] if c["verdict"] == "RESOLVABLE" else c["usable_bar"]
        if c["verdict"] != "RESOLVABLE":
            print(f"\n  [!] {name}: UNRESOLVABLE at the proposed bar. Per ground rule §5b the bar "
                  f"moves NOW, before the data is opened, to its usable_bar {c['usable_bar']:.2f}.")

    out = {"script": "diag_r1_anchor_reliability",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"), "preregistration": pre}

    if "--run" in sys.argv:
        pos = [a for a in sys.argv[1:] if not a.startswith("--")]
        meth_dir = Path(pos[0] if pos else r"D:\GSE165179")
        rna_dir = Path(pos[1] if len(pos) > 1 else r"D:\GSE165177")
        out["measurement"] = measure(pre, meth_dir, rna_dir)
    else:
        print("\n  Pre-registration only. Re-run with --run to measure.")

    _RESULTS / "diag_r1_anchor_reliability_results.json".write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_r1_anchor_reliability_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
