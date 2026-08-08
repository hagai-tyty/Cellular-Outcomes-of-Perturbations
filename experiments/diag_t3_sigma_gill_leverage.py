"""STAGE 1.5.6 step 3b — apply the MAX-LEVERAGE argument to T3, using Gill alone.

    python experiments/diag_t3_sigma_gill_leverage.py "C:/Users/hagay/Desktop/Gill"

READ-ONLY. Writes `results/diag_t3_sigma_gill_leverage_results.json`. `src/` untouched, no labels
move. Needs **no HFF stream** — Gill's six control samples and the shipped harmonizer are enough.

WHY
---
§5.5's floor precheck eliminated T2 with a max-leverage bound: `d` is affine in the clamped block's
ratio `R_f`, so `d_O1 · R_f/R_O1` is T2's ceiling, and for N2 that ceiling misses by 17.00 yr on a
16.671 yr spread. **The same argument can be aimed at T3, and it costs one Gill matrix.**

`σ_gill` is the ONLY remaining harmonization channel once T1 and T2 are out. It is estimated from
**five single control samples**, and which five changes per fold. This script computes, per fold, the
T3-only counterfactual — vary `σ_raw^(gill)`, hold the mask and the floor at O1, exactly §5.3's
ladder rung — and reports the per-gene ratio-of-ratios

    rho_g^(f) = ratio_g^(f) / ratio_g^(O1)     with   ratio_g = sigma_gill,g / sigma_hff,g

`d̂_f / d̂_O1` is a `δ_g·w_g`-weighted average of `rho_g`, so **it is bounded below by min(rho) and
above by max(rho) over the genes the clock reads.** If that interval excludes the observed
`d_N2/d_O1 = 0.306`, **T3 cannot explain N2 either**, by the same reasoning that eliminated T2 —
and the residue is outside harmonization entirely.

THE HONEST LIMIT, STATED UP FRONT
----------------------------------
The bound is a *containment* interval, not a prediction: `δ_g·w_g` carries both signs, so the
weighted average is not guaranteed to lie between min and max of `rho` when the weights are mixed.
**Reported here is the interval over the clock's genes AND the |w|-weighted mean**, with the sign
composition of the weights stated, so the write-up cannot present this as tighter than it is. It
NARROWS T3; it does not close it. Closing it needs `δ`, which needs the HFF stream.

WHAT IS VALIDATED FIRST (G0-style)
-----------------------------------
Before any counterfactual, the O1 fold's recomputed `sigma_gill` is checked **gene by gene** against
`runs/cellfate_multi/harmonization.json` — the O1 fold's own shipped harmonizer. If the
recomputation does not reproduce the artifact, nothing downstream is read.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

SHIPPED = ROOT / "runs" / "cellfate_multi" / "harmonization.json"
EPS = 1e-6
OBSERVED = {"N2": -7.352, "N3": -22.12, "O1": -24.023, "O2": -22.89, "Y1": -22.049, "Y2": -23.87}
REF_FOLD = "O1"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    if not SHIPPED.exists():
        print(f"FATAL: {SHIPPED} missing -- this script gates on the shipped O1 harmonizer.")
        return 1

    from cellfate.data.normalize import normalize_counts
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    ship = json.loads(SHIPPED.read_text("utf-8"))
    G = list(ship["genes"])
    sig_gill_ship = np.asarray(ship["stats"]["gill_bulk"]["sigma"], dtype=np.float64)
    sig_hff_ship = np.asarray(ship["stats"]["hff_sc"]["sigma"], dtype=np.float64)
    print(f"\n[shipped] {len(G)} genes | floors: gill {sig_gill_ship.min():.5f} "
          f"hff {sig_hff_ship.min():.5f}")

    # ---- Gill's six control samples, one per donor ---------------------------------------- #
    g_samples, g_genes, g_lin = dv.load_gill(Path(sys.argv[1]))
    g_norm = np.asarray(normalize_counts(g_lin), dtype=np.float64)
    ctrl = [(s.split("_")[0], i) for i, s in enumerate(g_samples)
            if re.search(r"_d\d+_", s) is None and "_Fib_" in s]
    donors = [d for d, _ in ctrl]
    if len(set(donors)) != len(donors) or len(donors) < 4:
        print(f"FATAL: expected one control per donor, got {donors}")
        return 1
    print(f"[gill]    {len(donors)} control samples, one per donor: {donors}")

    gi = {g: i for i, g in enumerate(g_genes)}
    missing = [g for g in G if g not in gi]
    idx = np.array([gi[g] for g in G if g in gi])
    keep = np.array([g in gi for g in G])
    C = g_norm[[i for _, i in ctrl]][:, idx]          # (6 donors, |G ∩ gill|)
    print(f"[align]   {keep.sum()}/{len(G)} shipped genes found in the Gill matrix "
          f"({len(missing)} missing)")

    wG = np.array([W.get(g, 0.0) for g in G])[keep]
    clock_mask = wG != 0.0
    sg_ship, sh_ship = sig_gill_ship[keep], sig_hff_ship[keep]
    floor_gill_O1 = float(sig_gill_ship.min())

    # ---- G0: does the recomputed O1 sigma_gill reproduce the shipped artifact? ------------- #
    def sigma_for(heldout: str) -> np.ndarray:
        rows = [k for k, d in enumerate(donors) if d != heldout]
        return C[rows].std(axis=0)

    raw_ref = sigma_for(REF_FOLD)
    rec_ref = np.maximum(raw_ref, floor_gill_O1)
    unclamped = sg_ship > floor_gill_O1 * (1 + 1e-9)
    if unclamped.sum() == 0:
        print("FATAL: no unclamped genes to compare -- G0 cannot run.")
        return 1
    rel = np.abs(rec_ref[unclamped] - sg_ship[unclamped]) / sg_ship[unclamped]
    med_rel, p90_rel = float(np.median(rel)), float(np.percentile(rel, 90))
    g0 = med_rel <= 0.05
    print(f"\n[G0] recomputed vs shipped sigma_gill on {int(unclamped.sum())} unclamped genes: "
          f"median rel.err {med_rel:.4f}, p90 {p90_rel:.4f}  -> {'PASS' if g0 else 'FAIL'}")
    if not g0:
        print("  G0 FAILED. The recomputation does not reproduce the pipeline's own harmonizer, so\n"
              "  nothing downstream is read. Likely causes: a different control definition, a\n"
              "  different normalization, or the shipped artifact is not the O1 fold.")

    # ---- T3-only counterfactual: vary sigma_raw(gill), hold mask + floor at O1 ------------- #
    ratio_ref = np.maximum(raw_ref, floor_gill_O1) / (sh_ship + EPS)
    d_ref = OBSERVED[REF_FOLD]
    out_folds = {}
    print(f"\n  {'fold':>5} {'observed':>9} {'d/d_O1':>7} | {'rho range on clock genes':>26} "
          f"{'|w|-mean':>9} {'in range?':>10}")
    for f in donors:
        raw_f = sigma_for(f)
        ratio_f = np.maximum(raw_f, floor_gill_O1) / (sh_ship + EPS)
        rho = ratio_f / ratio_ref
        rc = rho[clock_mask]
        lo, hi = float(rc.min()), float(rc.max())
        wm = float(np.average(rc, weights=np.abs(wG[clock_mask])))
        obs_ratio = OBSERVED[f] / d_ref
        inside = lo <= obs_ratio <= hi
        out_folds[f] = {"observed_day14": OBSERVED[f], "observed_over_ref": obs_ratio,
                        "rho_min_clock": lo, "rho_max_clock": hi, "rho_absw_mean_clock": wm,
                        "rho_median_all": float(np.median(rho)),
                        "observed_inside_rho_range": bool(inside)}
        print(f"  {f:>5} {OBSERVED[f]:9.2f} {obs_ratio:7.3f} | "
              f"[{lo:11.3f}, {hi:11.3f}] {wm:9.3f} {'YES' if inside else 'NO':>10}")

    # Does T3's leverage ORDER the folds the way the observed labels are ordered? The |w|-mean is
    # the wrong magnitude estimator (§5.2 A1 -- it drops delta), but a monotone relationship in the
    # ordering is still evidence about the CHANNEL, and it is scale-free so A1's objection does not
    # reach it.
    def _rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = float(pos)
        return r

    fs = list(out_folds)
    xr = _rank([out_folds[f]["rho_absw_mean_clock"] for f in fs])
    yr = _rank([out_folds[f]["observed_over_ref"] for f in fs])
    mx, my = sum(xr) / len(xr), sum(yr) / len(yr)
    num = sum((a - mx) * (b - my) for a, b in zip(xr, yr, strict=True))
    den = (sum((a - mx) ** 2 for a in xr) * sum((b - my) ** 2 for b in yr)) ** 0.5
    spearman = num / den if den else float("nan")
    print(f"\n  Spearman(|w|-weighted rho, observed d_f/d_O1) over {len(fs)} folds = "
          f"{spearman:+.3f}")

    n_pos = int((wG[clock_mask] > 0).sum())
    n_neg = int((wG[clock_mask] < 0).sum())
    print(f"\n  clock genes in the shipped space: {int(clock_mask.sum())} "
          f"({n_pos} positive weights, {n_neg} negative) -- MIXED SIGNS, so the |w|-mean is a "
          f"summary,\n  not a bound. The rho range is the containment interval.")

    n2 = out_folds.get("N2", {})
    verdict = "T3_CANNOT_EXPLAIN_N2" if (g0 and n2 and not n2["observed_inside_rho_range"]) \
        else ("T3_STILL_LIVE" if g0 else "G0_FAILED_NOT_READ")
    print(f"\n  => {verdict}")

    out = {"script": "diag_t3_sigma_gill_leverage", "utc": datetime.now(UTC).isoformat(),
           "ref_fold": REF_FOLD, "n_shipped_genes": len(G),
           "n_genes_aligned": int(keep.sum()), "n_clock_genes": int(clock_mask.sum()),
           "clock_weight_signs": {"positive": n_pos, "negative": n_neg},
           "g0": {"median_rel_err": med_rel, "p90_rel_err": p90_rel,
                  "n_unclamped": int(unclamped.sum()), "pass": bool(g0)},
           "folds": out_folds, "spearman_absw_rho_vs_observed": spearman, "verdict": verdict,
           "scope_limit": "Holds the mask and the floor at O1 (T3-only counterfactual). The rho "
                          "range is a containment interval, not a prediction: clock weights carry "
                          "both signs. NARROWS T3, does not close it -- closing needs delta."}
    (_RESULTS / "diag_t3_sigma_gill_leverage_results.json").write_text(
        json.dumps(out, indent=2), "utf-8")
    print("  wrote results/diag_t3_sigma_gill_leverage_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
