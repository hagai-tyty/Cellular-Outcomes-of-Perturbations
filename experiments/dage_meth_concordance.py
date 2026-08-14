"""P1 — does the transcriptomic ΔAge survive comparison with an independent ageing modality?

    python experiments/dage_meth_concordance.py

READ-ONLY. Writes `results/dage_meth_concordance_results.json`. `src/` untouched.
Graded against `plans/WORK_ORDER_2026_08_14.md` P1, revision 2.

**NOT** "is −42 real". A second imperfect clock cannot answer that. The question is whether the
transcriptomic effect survives comparison with methylation measured on the SAME samples.

`GSE165179` is the methylation twin of `GSE165177`: same donors, arms, days and sample names.
Both clocks are applied by reusing `diag_methylation_anchor`'s verified functions -- its loader,
its linear predictor, its Horvath anti-transform and its implied-intercept derivation -- rather
than re-deriving any of them. Re-deriving a transform is how this project acquired a double-log1p.

STRATIFICATION IS LOAD-BEARING
------------------------------
Every treatment exists in two cell states with its OWN negative control: cells RETURNED to
fibroblast identity (Gill's MPTR claim) and cells still IN the reprogramming phase. Pooling them
cost the RNA run its verdict once already. Both modalities are stratified identically here, and
each treated stratum is measured against its own matched contemporaneous control.

THE UNIT IS THE DONOR
---------------------
n = 3, `t(0.975, df=2) = 4.303`. Sample- and pair-level figures are shown for shape only. No R² is
quoted with n equal to the number of samples.
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
OUT = _RESULTS / "dage_meth_concordance_results.json"

MEDIR = Path(r"D:\GSE165179")
CLOCKS = [("horvath_skin_blood_2018", "Horvath skin & blood 2018"),
          ("horvath_multitissue_2013", "Horvath multi-tissue 2013")]
T2 = 4.302652729911275                      # t(0.975, df=2) -- three donors

# methylation `cell type` -> the stratum names the RNA side uses
CTYPE = {
    "Transiently reprogrammed fibroblast": "transient_fib",
    "Transient reprogramming intermediate": "transient_int",
    "Failed to transiently reprogram fibroblast": "failed_fib",
    "Failing to transiently reprogram intermediate": "failed_int",
    "Negative control fibroblast": "control_fib",
    "Negative control intermediate": "control_int",
    "Fibroblast": "day0",
}


def control_for(stratum: str) -> str:
    return "control_int" if stratum.endswith("_int") else "control_fib"


def ci3(vals):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if len(v) < 2:
        return (float(v[0]) if len(v) else float("nan"), float("nan"), float("nan"), len(v))
    m = float(v.mean())
    se = float(v.std(ddof=1)) / np.sqrt(len(v))
    t = T2 if len(v) == 3 else 1.96
    return m, m - t * se, m + t * se, len(v)


def methylation_ages(M, clock_file: str, meta: dict) -> dict[str, float]:
    """Per-sample methylation age, via the recorded module's own verified path."""
    clock = json.loads((REPO / "configs" / "clocks" / f"{clock_file}.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    samples, betas = M.load_betas(MEDIR / "GSE165179_Matrix_processed_transient.txt.gz", set(W))
    present = [s for s in samples if s in meta]
    lp = {s: M.linear_predictor(betas[s], W)[0] for s in present}
    d0 = [s for s in present if meta[s]["ctype"] == M.FIB]
    imp = M.implied_intercept([lp[s] for s in d0], [meta[s]["age"] for s in d0])
    k = float(clock.get("intercept", 0.0)) or (imp.get("mean", 0.0)
                                               if imp.get("status") == "OK" else 0.0)
    return {s: float(M.anti_trafo(lp[s] + k)) for s in present}


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()
    M = importlib.import_module("diag_methylation_anchor")

    print("\n" + "=" * 92)
    print("P1 — TRANSCRIPTOMIC ΔAge vs AN INDEPENDENT MODALITY (methylation), same samples")
    print("=" * 92)
    print("Donor is the unit: n = 3, t(0.975, df=2) = 4.303. Strata matched on both sides.")

    meta = M.load_series(MEDIR / "GSE165179_series_matrix.txt.gz")
    for m in meta.values():
        m["stratum"] = CTYPE.get(m["ctype"], "other")

    out: dict = {"script": "dage_meth_concordance", "prereg": "plans/WORK_ORDER_2026_08_14.md P1"}
    per_clock: dict[str, dict] = {}

    for cf, cname in CLOCKS:
        age = methylation_ages(M, cf, meta)

        # ---- P1.1 THE ADJUDICATION FLOOR --------------------------------------------------
        ctrl_vals: dict[tuple, list[float]] = {}
        for s, a in age.items():
            st = meta[s]["stratum"]
            if st.startswith("control_"):
                ctrl_vals.setdefault((meta[s]["donor"], meta[s]["day"], st), []).append(a)
        sds = [float(np.std(v, ddof=1)) for v in ctrl_vals.values() if len(v) > 1]
        sd_meth = float(np.sqrt(np.mean(np.square(sds)))) if sds else float("nan")

        # ---- ΔAge per (donor, day, stratum) against its OWN matched control ----------------
        cmean = {k: float(np.mean(v)) for k, v in ctrl_vals.items()}
        cell: dict[tuple, list[float]] = {}
        for s, a in age.items():
            st = meta[s]["stratum"]
            if not st.startswith(("transient_", "failed_")):
                continue
            ck = (meta[s]["donor"], meta[s]["day"], control_for(st))
            if ck in cmean:
                cell.setdefault((meta[s]["donor"], meta[s]["day"], st), []).append(a - cmean[ck])
        dage = {k: float(np.mean(v)) for k, v in cell.items()}

        by_donor: dict[str, dict[str, list[float]]] = {}
        for (d, _day, st), v in dage.items():
            by_donor.setdefault(st, {}).setdefault(d, []).append(v)
        donor_mean = {st: {d: float(np.mean(v)) for d, v in dd.items()}
                      for st, dd in by_donor.items()}

        per_clock[cf] = {"name": cname, "sd_meth": sd_meth, "n_ctrl_groups": len(sds),
                         "dage_cell": {f"{d}|{day:.0f}|{st}": v for (d, day, st), v in dage.items()},
                         "donor_mean": donor_mean}

        print(f"\n  {cname}")
        rows = []
        for st in ("transient_fib", "transient_int", "failed_fib", "failed_int"):
            if st not in donor_mean:
                continue
            m, lo, hi, n = ci3(list(donor_mean[st].values()))
            rows.append([st, str(n), f"{m:+.2f}", f"[{lo:+.2f},{hi:+.2f}]",
                         "YES" if (np.isfinite(hi) and hi < 0) else "no"])
        print(render_table(["stratum", "donors", "ΔAge_meth", "95% CI (donor)", "excludes 0?"],
                           rows, aligns=["l", "r", "r", "r", "l"]))
        mean_abs = float(np.mean([abs(v) for v in dage.values()])) if dage else float("nan")
        blocked = sd_meth >= 0.5 * mean_abs
        print(f"   P1.1 FLOOR: control-replicate SD = {sd_meth:.2f} yr over {len(sds)} groups; "
              f"mean |ΔAge_meth| = {mean_abs:.2f}")
        print(f"      -> {'CANNOT ADJUDICATE' if blocked else 'floor passed'} "
              f"(blocked if SD >= {0.5 * mean_abs:.2f})")
        per_clock[cf]["floor"] = {"sd_meth": sd_meth, "mean_abs_dage": mean_abs,
                                  "blocked": bool(blocked)}

    # ---- inter-clock disagreement: the rho=0.568 ceiling, made concrete --------------------
    a, b = per_clock[CLOCKS[0][0]], per_clock[CLOCKS[1][0]]
    shared = sorted(set(a["dage_cell"]) & set(b["dage_cell"]))
    dif = [a["dage_cell"][k] - b["dage_cell"][k] for k in shared]
    inter = float(np.sqrt(np.mean(np.square(dif)))) if dif else float("nan")
    print(f"\n  INTER-CLOCK DISAGREEMENT on the same {len(shared)} cells: RMS "
          f"{inter:.2f} yr  (the rho=0.568 ceiling, in years)")
    out["inter_clock_rms"] = inter

    # ---- join to the RNA side -------------------------------------------------------------
    rna = json.loads((_RESULTS / "dage_gse165177_results.json").read_text("utf-8"))
    rna_cell = {f"{r[0]}|{int(float(r[1]))}|{r[2]}": float(r[5]) for r in rna["M_E2"]}
    print(f"\n  RNA cells available: {len(rna_cell)}")

    print("\n" + "-" * 92)
    print("P1.2-P1.5 — DIRECTION, MAGNITUDE, TRAJECTORY, SCALE")
    print("-" * 92)
    for cf, cname in CLOCKS:
        pc = per_clock[cf]
        if pc["floor"]["blocked"]:
            print(f"\n  {cname}: floor NOT passed -- P1.2-P1.5 not read, per the work order.")
            continue
        keys = sorted(set(pc["dage_cell"]) & set(rna_cell))
        if not keys:
            print(f"\n  {cname}: no matched cells")
            continue
        # donor-level, per stratum
        print(f"\n  {cname} — matched cells: {len(keys)}")
        drows = []
        for st in ("transient_fib", "transient_int"):
            ks = [k for k in keys if k.endswith(st)]
            if not ks:
                continue
            dm_r: dict[str, list[float]] = {}
            dm_m: dict[str, list[float]] = {}
            for k in ks:
                d = k.split("|")[0]
                dm_r.setdefault(d, []).append(rna_cell[k])
                dm_m.setdefault(d, []).append(pc["dage_cell"][k])
            R = {d: float(np.mean(v)) for d, v in dm_r.items()}
            Me = {d: float(np.mean(v)) for d, v in dm_m.items()}
            donors = sorted(R)
            mr, lr, hr, nr = ci3([R[d] for d in donors])
            mm, lm, hm, _ = ci3([Me[d] for d in donors])
            same = sum(1 for d in donors if R[d] < 0 and Me[d] < 0)
            # SCALE: donor-level slope through the origin, the honest n=3 estimator
            sl = [R[d] / Me[d] for d in donors if abs(Me[d]) > 1e-9]
            beta = float(np.mean(sl)) if sl else float("nan")
            drows.append([st, str(nr), f"{mr:+.2f}", f"{mm:+.2f}", f"{same}/{len(donors)}",
                          f"{beta:.2f}", f"[{lr:+.1f},{hr:+.1f}]", f"[{lm:+.1f},{hm:+.1f}]"])
            per_clock[cf].setdefault("joined", {})[st] = {
                "rna_donor": R, "meth_donor": Me, "rna_ci": [mr, lr, hr],
                "meth_ci": [mm, lm, hm], "both_negative": f"{same}/{len(donors)}",
                "per_donor_ratio": sl, "mean_ratio": beta}
        print(render_table(["stratum", "donors", "ΔAge RNA", "ΔAge meth", "both<0",
                            "RNA/meth", "RNA 95% CI", "meth 95% CI"], drows,
                           aligns=["l", "r", "r", "r", "r", "r", "r", "r"]))
        print("   RNA/meth is the donor-mean RATIO and is DESCRIPTIVE ONLY -- it is unstable when")
        print("   the denominator approaches zero. Both denominators are printed beside it.")

    out["clocks"] = {k: {kk: vv for kk, vv in v.items() if kk != "dage_cell"}
                     for k, v in per_clock.items()}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
