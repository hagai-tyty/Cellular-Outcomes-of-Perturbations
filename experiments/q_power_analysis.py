"""
POWER ANALYSIS for Stage 1's conformal bar — run BEFORE run 2, to know what it can achieve.

    python experiments/q_power_analysis.py

WHY. Stage 1 pre-registers `conformal_coverage` in **0.85-0.95**, and the user ruled that
anything above 0.95 counts as a FAIL. Before spending 3.5 h of GPU on run 2, it is worth asking a
cheap question: *given this data, is that bar reachable at all?*

Two things make it doubtful, and neither is a bug:

1. **The pool is tiny.** After the bulk corpus is skipped, the inner-LODO pools ~14 age-valid
   cells from each of 5 training donors -- about 70 residuals. `q` is their 90th percentile, and
   a P90 from 70 samples carries real sampling error.
2. **Donor error scales differ 5.5x** (O1 MAE 5.39 vs N3 29.69, measured, `scorecard/baseline.json`).
   A single global `q` must serve all of them.

This simulates both. It needs no bundles, no GPU, and no shards -- only the per-fold MAEs already
in the baseline. Residuals are modelled as half-normal, which reproduces the plan's own
`P90/mean = 2.07` for a normal error; the real mixture is heavier-tailed (2.67 per MASTER_PLAN
S5d), so the true `q` is LARGER than modelled here and saturation is MORE likely, not less. The
conclusion is therefore conservative.

READ: if the per-fold coverages come out bimodal -- low-error donors pinned at 1.000, high-error
donors well under 0.85 -- then no single scalar `q` can put every fold inside 0.85-0.95, and a
run-2 "failure" on that metric is a property of donor heterogeneity rather than of the
calibration code.
"""
from __future__ import annotations

import numpy as np

# model dAge MAE per held-out fold, from scorecard/baseline.json
MAE = {"N2": 21.79, "N3": 29.69, "O1": 5.39, "O2": 7.54, "Y1": 7.28, "Y2": 14.06}
N_PER_DONOR = 14        # age-valid training cells per donor, after the corpus is skipped
LEVEL = 0.90
TRIALS = 4000
SEED = 0


def half_normal(rng, mae: float, n: int) -> np.ndarray:
    """|error| with the given mean. E|N(0,s)| = s*sqrt(2/pi)."""
    return np.abs(rng.normal(0.0, mae * np.sqrt(np.pi / 2), n))


def main() -> None:
    rng = np.random.default_rng(SEED)
    donors = list(MAE)

    print("\nSTAGE 1 CONFORMAL POWER ANALYSIS")
    print(f"  pool: {len(donors) - 1} training donors x {N_PER_DONOR} cells "
          f"= ~{(len(donors) - 1) * N_PER_DONOR} residuals")
    print("  bar : coverage in 0.85-0.95 (>0.95 pre-registered as FAIL)")

    # ---- 1. how much does q move, purely from sampling? ----
    qs = []
    for _ in range(TRIALS):
        held = rng.choice(donors)
        pool = np.concatenate([half_normal(rng, MAE[d], N_PER_DONOR)
                               for d in donors if d != held])
        k = int(np.ceil((pool.size + 1) * LEVEL))
        qs.append(np.sort(pool)[min(k, pool.size) - 1])
    qs = np.array(qs)
    lo, med, hi = (float(np.quantile(qs, x)) for x in (0.05, 0.50, 0.95))
    print(f"\n  [1] q from {TRIALS} resamples: median {med:.1f} yr, 90% range [{lo:.1f}, {hi:.1f}]")
    print(f"      spread / median = {(hi - lo) / med:.0%}  <- sampling noise alone")

    # ---- 2. what coverage does a single global q deliver per donor? ----
    print(f"\n  [2] coverage each held-out donor gets from q = {med:.1f} yr")
    covs = {}
    for d in donors:
        e = half_normal(rng, MAE[d], 20000)
        covs[d] = float(np.mean(e <= med))
        flag = ("SATURATED" if covs[d] > 0.95 else
                "under-covers" if covs[d] < 0.85 else "in range")
        print(f"      {d}: MAE {MAE[d]:5.2f} -> coverage {covs[d]:.3f}   {flag}")

    in_range = [d for d, c in covs.items() if 0.85 <= c <= 0.95]
    mean_cov = float(np.mean(list(covs.values())))
    print(f"\n      aggregate mean coverage: {mean_cov:.3f}")
    print(f"      folds INSIDE 0.85-0.95:  {len(in_range)}/{len(donors)}  {in_range}")

    print("\n   VERDICT:")
    if len(in_range) <= len(donors) // 2:
        print("     => The bar is NOT reachable per-fold with a single global q.")
        print("        Donor error scales differ ~5.5x, so one scalar cannot cover them all at")
        print("        90%: low-error donors saturate at 1.000 while high-error donors")
        print("        under-cover. Note the aggregate can still land inside the window while")
        print("        NO individual fold does -- run 1 showed exactly that (mean 0.873 from")
        print("        folds 0.381/0.857/1.000/1.000/1.000/1.000).")
        print("        A run-2 miss on this metric is then a finding about donor heterogeneity,")
        print("        NOT evidence that the calibration code is wrong. The response is a")
        print("        per-donor or conditional interval, pre-registered separately -- never")
        print("        shrinking q until the number fits, which is fitting the test.")
    else:
        print("     => The bar looks reachable; a miss in run 2 would point at the code.")
    print()


if __name__ == "__main__":
    main()
