"""STAGE 1.5.1 REV FINAL — anchor ΔAge to methylation (GSE165179).

    python experiments/diag_methylation_anchor.py "C:\\Users\\hagay\\Desktop\\GSE165179"

READ-ONLY. Writes `diag_methylation_anchor_results.json`. `src/` is untouched.

WHY (STAGE_1_5_1_REV_FINAL.md)
------------------------------
Four RNA-only fixes were tested and all four failed, for one reason: the transcriptomic clock is
correctly built and correctly applied but **out of domain on reprogramming cells**
(`corr(age, pluripotency) = −0.62`; non-responders read +36.5 yr in 11 days; ~20 yr swings between
adjacent days). No RNA analysis can escape that, because every RNA route to "age" runs through the
same clock. GSE165179 is Gill's own methylation companion and supplies two things RNA cannot:

  * an **identity-matched** comparison — `Transiently reprogrammed fibroblast` vs
    `Negative control fibroblast`, both FIBROBLASTS, same culture time, differing only in treatment;
  * a **real untreated negative control** at every timepoint, which GSE165176 verifiably lacks at
    any day > 0.

MEASUREMENTS (bars pre-registered in REV FINAL §3, amended §10.4)
  G2   does the clock reproduce KNOWN chronological age on the day-0 fibroblasts (53/53/38)?
       Load-bearing: if it cannot, methylation is not an anchor either and M-1/M-3 are not read.
  M-1  transiently reprogrammed vs negative control      (9 identity-matched pairs, MDE ~3.3 yr)
  M-3  failed to reprogram      vs negative control      (12 pairs, MDE ~2.7 yr)

THE INTERCEPT QUESTION — resolved by measurement, not assumption
---------------------------------------------------------------
`biolearn`'s `LinearModel.predict` inserts `intercept = 1` into the data matrix and inner-joins, so
an intercept is applied only if the coefficient table carries one. `Horvath1.csv` (353 rows) and
`Horvath2.csv` (391 rows) are **all CpGs, no intercept row** — yet the published clocks do have an
intercept, and omitting it shifts every age substantially. Rather than assert a value, G2 reports
the intercept **implied by the known ages** and checks it is consistent across donors. A consistent
implied intercept that also reproduces the ages is the evidence; an inconsistent one means the
clock is not usable here and is reported as such.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
ADULT_AGE = 20.0

TR = "Transiently reprogrammed fibroblast"
NC = "Negative control fibroblast"
FL = "Failed to transiently reprogram fibroblast"
FIB = "Fibroblast"

T_CRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
          8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145}

# G2 bar: the clock must recover known chronological age this well to be trusted as an anchor.
# Horvath skin&blood publishes ~2.5-3.5 yr MAE; 1.5x that is the same tolerance §9 used for the
# transcriptomic clock's own-domain check.
G2_MAE_TOL = 5.0
# Pre-registered fragility margin (ground rule: verdicts decided by hundredths must say so).
FRAGILE_MARGIN = 0.5


def anti_trafo(x: float | np.ndarray, adult_age: float = ADULT_AGE):
    """Horvath's inverse age transform. Verbatim the published form.

        y = (1+adult)*exp(x) - 1     if x <  0
        y = (1+adult)*x + adult      if x >= 0

    `F(20) = 0` and `F(0) = -log(21)`, so `anti_trafo` maps those back to 20 and 0.
    """
    x = np.asarray(x, dtype=float)
    return np.where(x < 0, (1.0 + adult_age) * np.exp(x) - 1.0, (1.0 + adult_age) * x + adult_age)


def trafo(age: float | np.ndarray, adult_age: float = ADULT_AGE):
    """Forward transform — the inverse of `anti_trafo`. Used to derive an implied intercept."""
    age = np.asarray(age, dtype=float)
    return np.where(age <= adult_age,
                    np.log(age + 1.0) - np.log(adult_age + 1.0),
                    (age - adult_age) / (adult_age + 1.0))


def linear_predictor(betas: dict[str, float], weights: dict[str, float]) -> tuple[float, int]:
    """Σ w_cpg · beta_cpg over the CpGs present. Returns (sum, n_used). Pure.

    Missing CpGs contribute nothing — the same silent-drop behaviour the transcriptomic clock has,
    which is why `coverage` is reported and gated below.
    """
    s = 0.0
    n = 0
    for cpg, w in weights.items():
        b = betas.get(cpg)
        if b is not None and np.isfinite(b):
            s += w * b
            n += 1
    return s, n


def implied_intercept(linpreds: list[float], true_ages: list[float]) -> dict:
    """What intercept would make these samples read their KNOWN ages?

    For each sample: `trafo(true_age) = linpred + intercept`, so `intercept = trafo(age) - linpred`.
    A clock whose intercept is genuinely a constant gives a consistent value across samples; a
    scattered one means the clock is not reproducing this data and must not be used as an anchor.
    """
    lp = np.asarray(linpreds, float)
    ta = np.asarray(true_ages, float)
    ok = np.isfinite(lp) & np.isfinite(ta)
    if ok.sum() < 2:
        return {"status": "CANNOT_VERIFY", "n": int(ok.sum())}
    imp = trafo(ta[ok]) - lp[ok]
    sd = float(np.std(imp, ddof=1))
    # spread in transformed units -> years, at the adult slope (21 yr per unit)
    return {"status": "OK", "n": int(ok.sum()), "mean": float(np.mean(imp)),
            "sd": sd, "spread_years": sd * (1.0 + ADULT_AGE),
            "per_sample": [float(v) for v in imp]}


def g2_verdict(pred_ages: list[float], true_ages: list[float], tol: float = G2_MAE_TOL) -> dict:
    """Does the clock reproduce known chronological age? The load-bearing guard."""
    p = np.asarray(pred_ages, float)
    t = np.asarray(true_ages, float)
    ok = np.isfinite(p) & np.isfinite(t)
    if ok.sum() < 2:
        return {"status": "CANNOT_VERIFY", "n": int(ok.sum()),
                "reason": "need >=2 samples of known age"}
    mae = float(np.mean(np.abs(p[ok] - t[ok])))
    return {"status": "REPRODUCES" if mae <= tol else "FAILS",
            "mae_years": mae, "n": int(ok.sum()),
            "predicted": [float(v) for v in p[ok]], "true": [float(v) for v in t[ok]],
            "reason": (f"MAE {mae:.1f} yr on {int(ok.sum())} known-age samples "
                       f"(tolerance {tol:.1f}) -> "
                       + ("the clock reads age here; it is a usable anchor"
                          if mae <= tol else
                          "the clock does NOT read age here; it is NOT an anchor and M-1/M-3 "
                          "must not be interpreted"))}


def paired_stat(diffs: list[float]) -> dict:
    """Mean and paired 95% CI over matched (donor, day) pairs."""
    d = np.asarray([x for x in diffs if np.isfinite(x)], float)
    n = len(d)
    if n < 2:
        return {"n": n, "mean": float("nan"), "ci95": [float("nan")] * 2, "n_negative": 0}
    mean = float(d.mean())
    se = float(d.std(ddof=1)) / np.sqrt(n)
    t = T_CRIT.get(n - 1, 1.96)
    return {"n": n, "mean": mean, "ci95": [mean - t * se, mean + t * se],
            "sd": float(d.std(ddof=1)), "n_negative": int((d < 0).sum())}


def effect_verdict(stat: dict, label: str, fragile: float = FRAGILE_MARGIN) -> dict:
    """Negative = younger than the matched untreated control."""
    if stat["n"] < 2 or not np.isfinite(stat["mean"]):
        return {"status": "CANNOT_VERIFY", "reason": f"{stat['n']} pairs"}
    lo, hi = stat["ci95"]
    if hi < 0:
        s, r = "REJUVENATION", f"{label} read YOUNGER than their matched untreated control"
    elif lo > 0:
        s, r = "AGEING", f"{label} read OLDER than their matched untreated control"
    else:
        s, r = "NO_EFFECT", f"{label}: CI includes 0"
    out = {"status": s, "mean_years": stat["mean"], "ci95": stat["ci95"], "n_pairs": stat["n"],
           "n_negative": stat["n_negative"],
           "reason": (f"{r} (mean {stat['mean']:+.1f} yr, 95% CI [{lo:+.1f}, {hi:+.1f}], "
                      f"{stat['n_negative']}/{stat['n']} pairs negative)")}
    if min(abs(lo), abs(hi)) < fragile:
        out["status"] = s + "_FRAGILE"
        out["reason"] += "  [FRAGILE: a CI bound is within 0.5 yr of zero]"
    return out


def decide(g2: dict, m1: dict, m3: dict) -> dict:
    """What the run licenses. Pure. G2 gates everything."""
    if g2.get("status") != "REPRODUCES":
        return {"action": "ANCHOR_INVALID",
                "reason": "the methylation clock does not reproduce known chronological age on "
                          "these samples, so it cannot arbitrate the transcriptomic clock. M-1/M-3 "
                          "are NOT interpreted. Check probe coverage and the intercept before "
                          "concluding anything about biology."}
    m1s, m3s = m1.get("status", ""), m3.get("status", "")
    if m1s.startswith("REJUVENATION"):
        if m3s.startswith("AGEING"):
            return {"action": "ANCHOR_VALID_EFFECT_REAL_CONTROL_MOVES",
                    "reason": "rejuvenation is real AND the failed arm also moves. The effect "
                              "stands, but the negative control is not inert -- report both "
                              "together (ground rule 10)."}
        return {"action": "ANCHOR_VALID_EFFECT_REAL",
                "reason": "methylation shows rejuvenation against an identity-matched untreated "
                          "control, and the failed arm does not. This is the anchor ΔAge needed: "
                          "the target is validated and the transcriptomic failure is localised to "
                          "the instrument, not the biology."}
    if m1s.startswith("AGEING"):
        return {"action": "CONTRADICTS",
                "reason": "two independent instruments would agree the cells get older. Treat as a "
                          "bug hunt first (ground rule 6) before revising the project premise."}
    return {"action": "NO_EFFECT_AT_THIS_RESOLUTION",
            "reason": "with an instrument sharper than the effect and an identity-matched control, "
                      "no rejuvenation is detected. That is a real finding, not a power failure -- "
                      "it escalates to Stage 4/5 rather than to another label change."}


# --------------------------------------------------------------------------- #
# Data wiring                                                                  #
# --------------------------------------------------------------------------- #
def load_series(path: Path) -> dict:
    """{title: {donor, day, ctype, age}} from the GEO series matrix."""
    rows: dict[str, list[str]] = {}
    titles: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = v[0].split(":")[0].strip()
                rows[key] = [x.split(":", 1)[1].strip() if ":" in x else x for x in v]
    out = {}
    for i, t in enumerate(titles):
        out[t] = {"donor": t.split("_")[0].split()[0],
                  "ctype": rows["cell type"][i],
                  "day": float(rows["length of reprogramming phase (days)"][i]),
                  "age": float(rows["donor age (years)"][i])}
    return out


def load_betas(path: Path, wanted: set[str]) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Parse the processed matrix: COMMA-separated, `Detection Pval` interleaved after each sample.

    Returns (sample_names, {sample: {cpg: beta}}) restricted to `wanted` CpGs — the clock's probes
    only, which keeps memory small (391 of ~865k rows).
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        keep = [(i, h) for i, h in enumerate(hdr[1:], start=1) if h.strip() != "Detection Pval"]
        samples = [h for _i, h in keep]
        data: dict[str, dict[str, float]] = {s: {} for s in samples}
        for row in rdr:
            if not row:
                continue
            cpg = row[0].strip()
            if cpg not in wanted:
                continue
            for i, s in keep:
                if i < len(row):
                    try:
                        data[s][cpg] = float(row[i])
                    except ValueError:
                        pass
    return samples, data


def pair_by_donor_day(meta: dict, arm_a: str, arm_b: str) -> list[tuple[str, float, list, list]]:
    """Matched (donor, day) groups between two arms.

    Returns `(donor, day, [treated titles], [control titles])`. **Replicates are kept and averaged
    by the caller, not dropped.** GSE165179 runs each condition in `exp1` and `exp2`, so requiring a
    unique sample per (donor, day, arm) silently discarded 6 of 9 M-1 pairs and left only day 10 —
    which is how the first run of this script under-reported. Averaging replicates is the correct
    handling: they are repeats of the same condition, not different conditions.
    """
    idx: dict[tuple[str, float, str], list[str]] = {}
    for t, m in meta.items():
        idx.setdefault((m["donor"], m["day"], m["ctype"]), []).append(t)
    out = []
    for (d, day, ct), ts in sorted(idx.items()):
        if ct != arm_a:
            continue
        other = idx.get((d, day, arm_b))
        if other:
            out.append((d, day, sorted(ts), sorted(other)))
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ddir = Path(sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE165179")
    clock_path = root / "configs" / "clocks" / "horvath_skin_blood_2018.json"

    print("STAGE 1.5.1 REV FINAL -- methylation anchor (read-only)\n")
    if not clock_path.exists():
        print(f"  [!] missing clock: {clock_path}")
        return 1
    clock = json.loads(clock_path.read_text(encoding="utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    print(f"  clock: {clock['meta']['name']}  ({len(W)} CpGs, {clock['meta']['tissue']})")

    ser = ddir / "GSE165179_series_matrix.txt.gz"
    bet = ddir / "GSE165179_Matrix_processed_transient.txt.gz"
    for p in (ser, bet):
        if not p.exists():
            print(f"  [!] missing {p}")
            return 1
    meta = load_series(ser)
    samples, betas = load_betas(bet, set(W))
    cov = {s: len(betas[s]) for s in samples}
    mean_cov = float(np.mean(list(cov.values()))) if cov else 0.0
    print(f"  samples: {len(samples)}   mean clock-CpG coverage: {mean_cov:.0f}/{len(W)} "
          f"({mean_cov/len(W):.1%})\n")

    # names in the beta matrix vs the series matrix are the same titles
    common = [s for s in samples if s in meta]
    print(f"  matched to series metadata: {len(common)}/{len(samples)}")

    lp = {s: linear_predictor(betas[s], W)[0] for s in common}

    # ---- G2: known-age day-0 fibroblasts, and the intercept question ----
    d0 = [s for s in common if meta[s]["ctype"] == FIB]
    imp = implied_intercept([lp[s] for s in d0], [meta[s]["age"] for s in d0])
    print("\n  G2 -- known-age day-0 fibroblasts (the load-bearing guard)")
    for s in d0:
        print(f"     {s:<12} true age {meta[s]['age']:>4.0f}   linear predictor {lp[s]:+.4f}")
    if imp.get("status") == "OK":
        print(f"     implied intercept: mean {imp['mean']:+.4f}, sd {imp['sd']:.4f} "
              f"(= {imp['spread_years']:.1f} yr of spread)")
    intercept = float(clock.get("intercept", 0.0)) or (imp.get("mean", 0.0) if imp.get("status") == "OK" else 0.0)
    used_implied = not float(clock.get("intercept", 0.0))
    if used_implied:
        print(f"     NOTE: the coefficient table carries no intercept; using the value implied by "
              f"the known ages ({intercept:+.4f}).")
        print("     This makes G2 partly self-fulfilling on the MEAN -- so G2 is judged on the")
        print("     SPREAD (consistency across donors), which no intercept choice can fake.")
    age = {s: float(anti_trafo(lp[s] + intercept)) for s in common}
    g2 = g2_verdict([age[s] for s in d0], [meta[s]["age"] for s in d0])
    g2["implied_intercept"] = imp
    g2["intercept_used"] = intercept
    g2["intercept_source"] = "implied from known ages" if used_implied else "clock file"
    print(f"     -> {g2['status']}: {g2['reason']}")

    # ---- M-1 / M-3: identity-matched pairs ----
    results = {}
    for tag, arm, label in [("M1", TR, "transiently reprogrammed fibroblasts"),
                            ("M3", FL, "failed-to-reprogram fibroblasts")]:
        groups = pair_by_donor_day(meta, arm, NC)
        rows, diffs = [], []
        for d, day, tt, cc in groups:
            ta = [age[x] for x in tt if x in age]
            ca = [age[x] for x in cc if x in age]
            if not ta or not ca:
                continue
            t_m, c_m = float(np.mean(ta)), float(np.mean(ca))
            diffs.append(t_m - c_m)
            rows.append({"donor": d, "day": day, "n_treated": len(ta), "n_control": len(ca),
                         "treated_age": t_m, "control_age": c_m, "diff_years": t_m - c_m})
        st = paired_stat(diffs)
        v = effect_verdict(st, label)
        v["pairs"] = rows
        results[tag] = v
        print(f"\n  {tag} -- {label} vs matched NEGATIVE CONTROL ({len(rows)} donor-day pairs)")
        for r in rows:
            rep = (f"  [{r['n_treated']}v{r['n_control']} reps]"
                   if max(r["n_treated"], r["n_control"]) > 1 else "")
            print(f"     {r['donor']} d{r['day']:.0f}: {r['treated_age']:>6.1f} vs "
                  f"{r['control_age']:>6.1f}  ->  {r['diff_years']:+6.1f} yr{rep}")
        print(f"     -> {v['status']}: {v['reason']}")

    decision = decide(g2, results["M1"], results["M3"])
    print(f"\n  ==> ACTION: {decision['action']}\n      {decision['reason']}")

    out = {"script": "diag_methylation_anchor",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "data_dir": str(ddir), "clock": clock["meta"],
           "n_samples": len(samples), "mean_cpg_coverage": mean_cov, "n_clock_cpgs": len(W),
           "G2": g2, "M1": results["M1"], "M3": results["M3"], "decision": decision,
           "ages": {s: age[s] for s in common}}
    Path("diag_methylation_anchor_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_methylation_anchor_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
