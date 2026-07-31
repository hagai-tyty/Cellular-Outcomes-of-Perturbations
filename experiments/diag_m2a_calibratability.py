"""STAGE 1.5.2 M-2a — does the RNA clock track METHYLATION age? (paired GSE165177 x GSE165179)

    python experiments/diag_m2a_calibratability.py "D:\\GSE165177" "D:\\GSE165179"

READ-ONLY. Writes `diag_m2a_calibratability_results.json`. `src/` untouched.

WHAT THIS IS. Stage 1.5.2's control-free core (§5): "Absolute ages, no control needed." Because it
needs no control baseline, it is NOT gated by G-a (baseline-count visibility) or G-b (donor-age
wiring) — those gate M-2b, the ΔAge-shaped question. M-2b is deliberately NOT run here.

THE CONFOUND, AND WHY THE HEADLINE NUMBER IS NOT THE ANSWER (§4). Both modalities move with
reprogramming progress: 1.5.1 measured corr(age_rna, pluripotency) = -0.62, and methylation age
falls -24 to -27 yr over the same axis. So a clock carrying NO age information -- one that only
detects "is this cell reprogramming?" -- would still produce a high age_rna <-> age_meth
correlation. Three readings are therefore computed and only the confound-free ones count:

  rho_all      across all conditions        DESCRIPTIVE ONLY, never a pass criterion
  rho_within   within each arm separately   descriptive here (demoted 2026-07-31: UNRESOLVABLE
                                            at the real per-arm n; see §6's frozen bars)
  rho_partial  partialling out pluripotency DECISIVE -- bar rho >= 0.50, RESOLVABLE at n=68 (99.4%)

The pluripotency score is the existing `OSKM_PLURIPOTENCY` signature from `diag_clock_validity.py`,
reused verbatim so it cannot be tuned for this stage.

GEOMETRY (verified on download, titles only): 90 joined pairs -> 68 unique (donor, arm, day)
conditions; donors O1/O2/O3; days 10/13/15/17. `exp1`/`exp2` replicates are AVERAGED, not treated as
independent (1.5.1's unit-of-analysis rule).

BOTH Horvath clocks are run and BOTH reported, whichever way they fall. A criterion met on one and
not the other is SPLIT, which §7 counts as a failure, not a pass.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "local_runners", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

ANNOT_COLS = 12          # Probe..Distance precede the sample columns (same layout as GSE165176)
RHO_BAR = 0.50           # frozen 2026-07-31; RESOLVABLE at n=68 (99.4%)


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def parse_title(t: str):
    m = re.match(r"^(O\d)_(.+?)_(\d+)days_(exp\d)$", t)
    return {"donor": m.group(1), "arm": m.group(2), "day": int(m.group(3))} if m else None


def load_rna(rna_dir: Path):
    """Merge the two Log2-RPM parts, dedup gene symbols, return (samples, genes, linear RPM)."""
    import pandas as pd
    frames = []
    for f in sorted(rna_dir.glob("*Log2_RPM*.txt.gz")):
        df = pd.read_csv(f, sep="\t", low_memory=False)
        cols = list(df.columns)
        frames.append(df.set_index(cols[0])[cols[ANNOT_COLS:]])
    merged = frames[0].join(frames[1:], how="inner") if len(frames) > 1 else frames[0]
    lin = np.power(2.0, merged.to_numpy(dtype=np.float64)) - 1.0     # Log2 RPM -> linear RPM
    lin[lin < 0] = 0.0
    sym = [str(s) for s in merged.index]
    # dedup symbols, keeping the highest-expressed row (GillReprogrammingSource's rule)
    order = np.argsort(-lin.sum(axis=1))
    seen, keep = set(), []
    for i in order:
        if sym[i] not in seen:
            seen.add(sym[i])
            keep.append(i)
    keep.sort()
    return list(merged.columns), [sym[i] for i in keep], lin[keep, :].T   # samples x genes


def rna_ages_and_pluripotency(rna_dir: Path):
    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    dcv = _load("diag_clock_validity", ROOT / "experiments" / "diag_clock_validity.py")

    samples, genes, lin = load_rna(rna_dir)
    norm = normalize_counts(lin)                                   # log1p CP10k, per sample
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    ages = clock.predict_age(norm, genes)
    idx = [i for i, g in enumerate(genes) if g in dcv.OSKM_PLURIPOTENCY]
    plu = norm[:, idx].mean(axis=1) if idx else np.zeros(len(samples))
    return (dict(zip(samples, map(float, ages), strict=True)),
            dict(zip(samples, map(float, plu), strict=True)), len(idx))


def meth_ages(meth_dir: Path) -> dict:
    """Methylation age RANK per sample, for both Horvath clocks (reuses 1.5.1's verified code).

    Returns the LINEAR PREDICTOR, not the transformed age, and that is deliberate: age is
    `anti_trafo(lp + k)`, which is strictly monotone in `lp`, so **Spearman on lp is identical to
    Spearman on age** for any intercept k. Correlating lp therefore removes the implied-intercept
    dependency entirely rather than merely sweeping it (`REV FINAL` §4.3's argument, applied here).
    """
    dma = _load("diag_methylation_anchor", ROOT / "experiments" / "diag_methylation_anchor.py")
    meta = dma.load_series(meth_dir / "GSE165179_series_matrix.txt.gz")
    bpath = meth_dir / "GSE165179_Matrix_processed_transient.txt.gz"
    out = {}
    for cfile, cname in dma.CLOCKS:
        clock = json.loads((ROOT / "configs" / "clocks" / f"{cfile}.json").read_text(
            encoding="utf-8"))
        W = {k: float(v) for k, v in clock["weights"].items()}
        samples, betas = dma.load_betas(bpath, set(W))
        out[cname] = {s: float(dma.linear_predictor(betas[s], W)[0])
                      for s in samples if s in meta}
    return out


def spearman(x, y):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return float("nan")
    return float(spearmanr(x[ok], y[ok]).correlation)


def partial_spearman(x, y, z):
    """Spearman of x,y after removing z from both (rank-residualised)."""
    from scipy.stats import rankdata
    x, y, z = (rankdata(np.asarray(v, float)) for v in (x, y, z))
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    return float(np.corrcoef(rx, ry)[0, 1])


def main() -> int:
    rna_dir = Path(sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE165177")
    met_dir = Path(sys.argv[2] if len(sys.argv) > 2 else r"D:\GSE165179")
    print("STAGE 1.5.2 M-2a — is the RNA clock calibratable against methylation? (read-only)\n")
    print(f"  DECISIVE bar: rho_partial >= {RHO_BAR:.2f} on BOTH clocks (frozen 2026-07-31, "
          "RESOLVABLE at n=68, 99.4%)")
    print("  rho_all and rho_within are DESCRIPTIVE ONLY and are never a pass criterion.\n")

    age_rna, plu, n_plu = rna_ages_and_pluripotency(rna_dir)
    print(f"  RNA: {len(age_rna)} samples; pluripotency score over {n_plu} OSKM genes")
    ages_m = meth_ages(met_dir)
    for c, d in ages_m.items():
        print(f"  METH ({c}): {len(d)} samples  [linear predictor; rank-identical to age]")

    # pair on exact title, then average exp replicates into (donor, arm, day) conditions
    results = {"n_plu_genes": n_plu, "bar": RHO_BAR, "clocks": {}}
    for cname, am in ages_m.items():
        cond = defaultdict(lambda: {"rna": [], "met": [], "plu": []})
        for t, ar in age_rna.items():
            if t in am and (p := parse_title(t)):
                k = (p["donor"], p["arm"], p["day"])
                cond[k]["rna"].append(ar)
                cond[k]["met"].append(am[t])
                cond[k]["plu"].append(plu[t])
        rows = [{"donor": k[0], "arm": k[1], "day": k[2],
                 "age_rna": float(np.mean(v["rna"])), "age_meth": float(np.mean(v["met"])),
                 "plu": float(np.mean(v["plu"]))} for k, v in sorted(cond.items())]
        R = [r["age_rna"] for r in rows]
        M = [r["age_meth"] for r in rows]
        P = [r["plu"] for r in rows]
        r_all = spearman(R, M)
        r_par = partial_spearman(R, M, P)
        per_arm = {}
        for arm in sorted({r["arm"] for r in rows}):
            s = [r for r in rows if r["arm"] == arm]
            per_arm[arm] = {"n": len(s),
                            "rho": spearman([r["age_rna"] for r in s], [r["age_meth"] for r in s])}
        print(f"\n  === {cname} ===   n = {len(rows)} conditions")
        print(f"     rho_all      {r_all:+.3f}   [descriptive only]")
        print(f"     rho_partial  {r_par:+.3f}   [DECISIVE, bar >= {RHO_BAR:.2f}]  "
              f"-> {'PASS' if r_par >= RHO_BAR else 'FAIL'}")
        print("     rho_within (descriptive):")
        for a, v in sorted(per_arm.items(), key=lambda x: -x[1]["n"]):
            print(f"        {a:48s} n={v['n']:3d}  rho {v['rho']:+.3f}")
        results["clocks"][cname] = {"n_conditions": len(rows), "rho_all": r_all,
                                    "rho_partial": r_par, "rho_within": per_arm, "rows": rows}

    passes = [v["rho_partial"] >= RHO_BAR for v in results["clocks"].values()]
    verdict = "PASS" if all(passes) else ("SPLIT" if any(passes) else "FAIL")
    results["verdict_m2a"] = verdict
    print(f"\n  ==> M-2a VERDICT: {verdict}"
          f"{'  (both clocks agree)' if verdict != 'SPLIT' else '  -- §7 counts SPLIT as a failure'}")
    print("  M-2b NOT RUN: it is the ΔAge-shaped question and is gated on G-a (baseline visibility).")

    Path("diag_m2a_calibratability_results.json").write_text(
        json.dumps({"script": "diag_m2a_calibratability",
                    "utc": datetime.now(UTC).isoformat(timespec="seconds"), **results},
                   indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_m2a_calibratability_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
