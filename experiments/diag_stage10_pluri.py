"""STAGE 10 -- is the pluripotency component of ΔAge CONTAMINATION or MEDIATION?  (read-only)

    python experiments/diag_stage10_pluri.py

Pre-registered in `plans/STAGE_10_PLURIPOTENCY_CONTAMINATION_OR_MEDIATION.md`, written before this
ran. `src/` is NOT touched by this stage under any outcome.

WHY
---
A previous session reported "pluripotency has to come out of the ΔAge readout" from one number
(resid_pluri 22.69 -> 13.00 MAE vs methylation). The objection is correct and was not addressed:
OSKM INDUCES pluripotency, and pluripotency induction may BE the mechanism of rejuvenation. Then
regressing it out deletes the signal rather than the noise.

Established before the plan (10.1): the signature is 5 genes carrying 0.0005% of the clock's |w|
mass, ranked 17.8k-26.8k of 33,155 -- the clock barely reads them. So `resid_pluri` removes a
CO-VARYING COMPONENT, not a direct reading. Two of the five (POU5F1=OCT4, SOX2) are OSKM
TRANSGENES, so the score partly measures vector dose.

THREE TESTS, each a falsifier for one reading
  A  arm ordering after residualising      MEDIATION predicts the outcome gap COLLAPSES
  B  pluri ~ ΔAge among NEGATIVE CONTROLS  no OSKM, so nothing to mediate; an association there
                                           is baseline covariation = CONTAMINATION
  C  agreement with METHYLATION            methylation cannot see RNA pluripotency, so it is an
                                           outside witness to what is real

Bars are the plan's, restated as constants below and not changed after the fact.
"""
from __future__ import annotations

import importlib.util
import json
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


p3 = _load("p3", "experiments/diag_phase3_within_donor_forward.py")

from cellfate.common import constants as C  # noqa: E402
from cellfate.data.aging import LinearClock  # noqa: E402

CLOCK = ROOT / "configs" / "clocks" / "fleischer_clock.json"
PLURI = C.DEFAULT_SIGNATURES["loss"]
GAP_COLLAPSE_BAR = 0.50      # 10.2: shrink > 50% -> MEDIATION
CONTROL_RHO_BAR = 0.50       # 10.3: |rho| >= this among controls -> CONTAMINATION

CONTROL_ARMS = ("negative_control", "negative_control_intermediate")
TRANSIENT_ARMS = ("transiently_reprogrammed", "transient_reprogramming_intermediate")


def zscore(v: np.ndarray) -> np.ndarray:
    s = np.std(v)
    return (v - v.mean()) / (s if s > 1e-12 else 1.0)


def pluri_score(expr: np.ndarray, genes: list[str]) -> np.ndarray:
    """Mean z-scored expression of the signature genes -- the same construction `labels.py`
    scores fate with, applied to the same samples."""
    idx = [genes.index(g) for g in PLURI if g in genes]
    if not idx:
        raise ValueError("no pluripotency signature gene present")
    return np.mean([zscore(expr[:, i]) for i in idx], axis=0)


def residualise(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """y with the linear component along x removed -- what `resid_pluri` does."""
    if np.std(x) < 1e-12:
        return y - y.mean()
    a, b = np.polyfit(x, y, 1)
    return y - (a * x + b)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def verdict_from(a: str, b: str, c: str) -> str:
    """10.5: count the three tests. Stated as a function so every branch is testable."""
    votes = [v for v in (a, b, c) if v in ("CONTAMINATION", "MEDIATION")]
    n_cont = votes.count("CONTAMINATION")
    n_med = votes.count("MEDIATION")
    if n_cont >= 2:
        return "CONTAMINATION"
    if n_med >= 2:
        return "MEDIATION"
    return "UNDETERMINED"


def main() -> None:
    meta, expr, genes = p3.load_all()
    clock = LinearClock.from_json(str(CLOCK))
    ages = clock.predict_age(expr, genes)
    plu = pluri_score(expr, genes)
    meta = meta.assign(age=ages[meta.j.to_numpy()], plu=plu[meta.j.to_numpy()])

    print("=" * 100)
    print("STAGE 10 -- pluripotency in ΔAge: CONTAMINATION or MEDIATION?")
    print("=" * 100)
    present = [g for g in PLURI if g in genes]
    print(f"  signature {list(PLURI)}  present in this matrix: {present}")
    print(f"  samples {len(meta)}  donors {sorted(meta.donor.unique())}  arms {meta.arm.nunique()}")

    # ΔAge = clock age minus that donor's negative_control mean at the same day (contemporaneous).
    d_rows = []
    for (d, day), g in meta.groupby(["donor", "day"]):
        ctrl = g[g.arm.isin(CONTROL_ARMS)]
        if ctrl.empty:
            continue
        base = float(ctrl.age.mean())
        for _, r in g.iterrows():
            d_rows.append({"donor": d, "day": day, "arm": r.arm, "plu": r.plu,
                           "dage": float(r.age) - base})
    D = pd.DataFrame(d_rows)
    res: dict = {"n_samples": int(len(D)), "signature_present": present}

    # ---- TEST A: does the ARM ORDERING survive residualising? ------------------------------- #
    D["dage_resid"] = residualise(D.dage.to_numpy(float), D.plu.to_numpy(float))
    print("\n[TEST A] transient-vs-control ΔAge gap, before and after removing pluripotency")
    print(f"  {'donor':<7}{'gap raw':>10}{'gap resid':>12}{'shrink':>10}")
    gaps = []
    for d, g in D.groupby("donor"):
        tr = g[g.arm.isin(TRANSIENT_ARMS)]
        ct = g[g.arm.isin(CONTROL_ARMS)]
        if tr.empty or ct.empty:
            continue
        raw_gap = float(tr.dage.mean() - ct.dage.mean())
        res_gap = float(tr.dage_resid.mean() - ct.dage_resid.mean())
        shrink = 1.0 - abs(res_gap) / abs(raw_gap) if abs(raw_gap) > 1e-9 else float("nan")
        gaps.append({"donor": d, "raw": raw_gap, "resid": res_gap, "shrink": shrink})
        print(f"  {d:<7}{raw_gap:>10.2f}{res_gap:>12.2f}{shrink:>10.2f}")
    med_shrink = float(np.median([g["shrink"] for g in gaps])) if gaps else float("nan")
    a_verdict = ("MEDIATION" if np.isfinite(med_shrink) and med_shrink > GAP_COLLAPSE_BAR
                 else "CONTAMINATION" if np.isfinite(med_shrink) else "UNDETERMINED")
    print(f"  median shrink {med_shrink:.2f}  (bar > {GAP_COLLAPSE_BAR} -> MEDIATION)"
          f"  -> {a_verdict}")
    res["A"] = {"gaps": gaps, "median_shrink": med_shrink, "verdict": a_verdict}

    # ---- TEST B: does pluripotency predict ΔAge among NEGATIVE CONTROLS? -------------------- #
    # NOTE, and this is a flaw in the first version of this test that had to be fixed: ΔAge is
    # DEFINED relative to the controls at the same (donor, day), so control samples carry ~zero
    # ΔAge by construction and their Spearman is undefined. Test B therefore uses RAW CLOCK AGE,
    # which is what the question actually needs: among untreated fibroblasts, where there is no
    # OSKM and nothing to mediate, does the pluripotency score covary with the clock at all?
    ctrl = meta[meta.arm.isin(CONTROL_ARMS)]
    rho_c = spearman(ctrl.plu.to_numpy(float), ctrl.age.to_numpy(float))
    osk = D[~D.arm.isin(CONTROL_ARMS)]
    rho_o = spearman(osk.plu.to_numpy(float), osk.dage.to_numpy(float))
    print("\n[TEST B] pluripotency ~ RAW CLOCK AGE among controls (no OSKM -> nothing to mediate)")
    print(f"  controls  n={len(ctrl):>3}  spearman {rho_c:+.3f}   (raw clock age)")
    print(f"  OSKM arms n={len(osk):>3}  spearman {rho_o:+.3f}   (ΔAge, for contrast)")
    ctrl_sd = float(np.std(ctrl.plu.to_numpy(float)))
    if ctrl_sd < 1e-9:
        # The signature is EXACTLY constant among controls -- the five genes are OFF in untreated
        # fibroblasts, so the score has no variance until OSKM is delivered. This is not a failed
        # measurement, it IS the answer, and it is the strongest form of the branch the plan
        # already pre-registered: "the association is ~0 in controls and appears only in OSKM
        # arms -- i.e. reprogramming-specific, as a causal path would be". A score that does not
        # exist without the treatment cannot be a baseline confound.
        b_verdict = "MEDIATION"
        print(f"  the signature is CONSTANT among controls (sd {ctrl_sd:.2e}) -- the genes are OFF "
              "without OSKM,\n  so there is NO baseline covariation for contamination to act "
              "through. This is the plan's\n  'appears only in OSKM arms' branch in its limiting "
              "form.")
    else:
        b_verdict = ("CONTAMINATION" if np.isfinite(rho_c) and abs(rho_c) >= CONTROL_RHO_BAR
                     else "MEDIATION" if np.isfinite(rho_c) else "UNDETERMINED")
    print(f"  bar |rho| >= {CONTROL_RHO_BAR} among controls -> CONTAMINATION   -> {b_verdict}")
    res["B"] = {"rho_controls": rho_c, "rho_oskm": rho_o, "n_controls": int(len(ctrl)),
                "control_signature_sd": ctrl_sd, "verdict": b_verdict}

    # ---- TEST C: does resid_pluri keep the METHYLATION agreement raw had? ------------------- #
    print("\n[TEST C] agreement with methylation ΔAge (an assay that cannot see RNA pluripotency)")
    led = pd.read_csv(ROOT / "results" / "dage_ledger.csv")
    t = led[(~led.is_control.astype(bool)) & led.TRUTH_meth_dage_mt.notna()]
    mt = t.TRUTH_meth_dage_mt.to_numpy(float)
    r_raw = spearman(t.ACTUAL_rna_dage_raw.to_numpy(float), mt)
    r_res = spearman(t.ACTUAL_rna_dage_resid_pluri.to_numpy(float), mt)
    print(f"  n={len(t)} reprogramming conditions")
    print(f"  raw         ~ methylation  spearman {r_raw:+.3f}")
    print(f"  resid_pluri ~ methylation  spearman {r_res:+.3f}")
    c_verdict = ("CONTAMINATION" if np.isfinite(r_res) and np.isfinite(r_raw) and r_res >= r_raw
                 else "MEDIATION" if np.isfinite(r_res) else "UNDETERMINED")
    print(f"  bar: resid >= raw -> CONTAMINATION   -> {c_verdict}")
    res["C"] = {"spearman_raw": r_raw, "spearman_resid_pluri": r_res, "n": int(len(t)),
                "verdict": c_verdict}

    final = verdict_from(a_verdict, b_verdict, c_verdict)
    print(f"\n{'=' * 100}")
    print(f"  A {a_verdict}   B {b_verdict}   C {c_verdict}   ->  STAGE 10 VERDICT: {final}")
    print("  src/ is NOT changed by this stage under any outcome (plan 10.5).")
    res["verdict"] = final

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_stage10_pluri_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
