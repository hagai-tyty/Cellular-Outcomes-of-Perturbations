"""
CAN the 0.85-0.95 coverage bar be MET? -- feasibility of an adaptive conformal interval.

    python experiments/q_adaptive_feasibility.py

THE PROBLEM (measured, `q_power_analysis.py`). A single global `q` cannot hold every fold at 90%:
donor error scales differ 5.5x (O1 MAE 5.39 vs N3 29.69), so one scalar saturates the easy donors
at 1.000 and under-covers the hard ones at 0.67. 0 of 6 folds land in the window.

THE STANDARD FIX. Normalized (locally-adaptive) split conformal: instead of scoring
``|y - mu| <= q`` , score the NORMALIZED residual

    r_i = |y_i - mu_i| / s_i      ->      interval = mu_i +- q_norm * s_i

where ``s_i`` is any per-cell difficulty estimate computable at inference. The interval then
widens for hard cells and tightens for easy ones, so coverage can hold across donors WITHOUT
needing the donor's label. Requires no new data -- only a signal we already produce.

THE QUESTION THIS ANSWERS. Not "is the idea sound" (it is, and it is standard) but:
**how well must `s` track the donor's true error scale before every fold lands in 0.85-0.95?**
If the answer is "nearly perfectly", the approach is not viable here and we should say so. If a
partial signal suffices, it tells us how good a normalizer we need to find.

CANDIDATE NORMALIZERS already available at inference, none needing new data:
  * `sigma_age`      -- the ensemble spread. Stage 1 shows its MAGNITUDE is wrong (2.4 vs ~14 yr),
                        but normalization only needs its RANKING across cells to be informative.
  * Mahalanobis distance to the training manifold (already computed for the OOD gate).
  * |mu_age| itself, if larger predicted changes carry larger errors.

`alpha` below is the exponent on the donor's true scale: alpha=0 is a signal that ignores donor
difficulty entirely (equivalent to today's global q), alpha=1 tracks it perfectly. `noise` is
per-cell log-normal scatter on top, i.e. how noisy the signal is at the individual-cell level.
"""
from __future__ import annotations

import numpy as np

MAE = {"N2": 21.79, "N3": 29.69, "O1": 5.39, "O2": 7.54, "Y1": 7.28, "Y2": 14.06}
N_CAL_PER_DONOR = 14        # inner-LODO cells per training donor, after the corpus is skipped
N_TEST = 20000              # large, so the reported coverage is the donor's TRUE coverage
LEVEL = 0.90
LO, HI = 0.85, 0.95
SEED = 0


def errors(rng, scale: float, n: int) -> np.ndarray:
    """|error| with mean == scale."""
    return np.abs(rng.normal(0.0, scale * np.sqrt(np.pi / 2), n))


def signal(rng, scale: float, n: int, alpha: float, noise: float) -> np.ndarray:
    """A difficulty estimate that tracks the donor's scale to degree ``alpha``."""
    return (scale ** alpha) * np.exp(rng.normal(0.0, noise, n))


def evaluate(rng, alpha: float, noise: float) -> dict:
    """Leave-one-donor-out: fit q_norm on 5 donors, measure coverage on the 6th."""
    donors = list(MAE)
    covs = {}
    for held in donors:
        cal_r = []
        for d in donors:
            if d == held:
                continue
            e = errors(rng, MAE[d], N_CAL_PER_DONOR)
            s = signal(rng, MAE[d], N_CAL_PER_DONOR, alpha, noise)
            cal_r.append(e / s)
        r = np.concatenate(cal_r)
        k = int(np.ceil((r.size + 1) * LEVEL))
        q_norm = np.sort(r)[min(k, r.size) - 1]

        e = errors(rng, MAE[held], N_TEST)
        s = signal(rng, MAE[held], N_TEST, alpha, noise)
        covs[held] = float(np.mean(e <= q_norm * s))
    return covs


def main() -> None:
    rng = np.random.default_rng(SEED)
    print("\nADAPTIVE CONFORMAL FEASIBILITY -- can every fold reach 0.85-0.95?")
    print(f"  {len(MAE)} donors, {N_CAL_PER_DONOR} calibration cells each, level {LEVEL}")
    print("  alpha = how well the normalizer tracks the donor's true error scale "
          "(0 = not at all, 1 = perfectly)")

    for noise in (0.0, 0.3, 0.6):
        print(f"\n  per-cell signal noise (log-normal sigma) = {noise}")
        print(f"    {'alpha':>6} | " + " ".join(f"{d:>6}" for d in MAE) + " | in-range  worst")
        print(f"    {'-' * 6}-+-" + "-" * (7 * len(MAE)) + "-+---------------")
        for alpha in (0.0, 0.25, 0.5, 0.75, 1.0):
            covs = evaluate(rng, alpha, noise)
            n_ok = sum(LO <= c <= HI for c in covs.values())
            worst = min(covs.values(), key=lambda c: min(abs(c - LO), abs(c - HI))
                        if not (LO <= c <= HI) else 1e9)
            mark = "  <== ALL PASS" if n_ok == len(MAE) else ""
            print(f"    {alpha:>6.2f} | " + " ".join(f"{covs[d]:>6.3f}" for d in MAE)
                  + f" | {n_ok}/{len(MAE)}      {worst:.3f}{mark}")

    print("\n   READ:")
    print("     alpha=0 reproduces today's global q -- the saturate/under-cover split.")
    print("     The alpha at which the in-range column first reaches 6/6 is the bar the")
    print("     normalizer must clear. If that is ~1.0, only a near-perfect difficulty signal")
    print("     works and the approach is fragile here. If 0.5-0.75 suffices, a partial signal")
    print("     -- e.g. the ensemble spread's RANKING, even with its magnitude wrong -- is")
    print("     enough, and the next step is to MEASURE which available signal reaches it.")
    print("\n   NOTE: this says nothing about whether such a signal exists in our data. It says")
    print("   how good one would have to be. Measuring the candidates is the follow-up, and it")
    print("   needs the real cross-donor residuals from run 2.\n")


if __name__ == "__main__":
    main()
