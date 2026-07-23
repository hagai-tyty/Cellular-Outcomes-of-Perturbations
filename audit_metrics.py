"""Is each scorecard criterion RESOLVABLE at the geometry it is actually measured on?

A bar is not a bar unless a system that meets the underlying scientific intent passes it
reliably. Run 3 raised the question for `fate_ece`: it is graded as the mean of per-fold ECEs
over ~21 held-out cells in 10 bins, and that estimator is biased upward by construction. If a
PERFECTLY calibrated model fails the bar most of the time, the bar tests the sample size rather
than the model, and both passing and failing it are uninformative.

This measures that directly instead of asserting it, for every criterion the dump can reach.

    python audit_metrics.py diag_dump/

TWO KINDS OF CRITERION, TWO KINDS OF QUESTION
---------------------------------------------
ABSOLUTE (the TARGETs -- `fate_ece`, `conformal_coverage`): compared against a fixed threshold.
  The question is FALSE-NEGATIVE RATE: how often does a system that genuinely satisfies the
  intent get reported as failing? Simulated from the null where the intent holds exactly
  (`y ~ Bernoulli(p)` for calibration; hits ~ Bernoulli(level) for coverage), so any departure
  measured is pure estimator behaviour.

COMPARATIVE (the GUARDs -- `dage_mae_model`, `rank_model_dage`, `fate_prauc`, `fate_roc`):
  compared between snapshots, and required to read "noise (CI incl. 0)". The question is POWER:
  how large must a real regression be before the paired 6-fold CI stops calling it noise? A guard
  that cannot detect a regression is not protecting anything, however reassuring its verdict.

Nothing here reads a calibrator choice off the graded folds; it characterises the instrument.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from cellfate.common.constants import SAFE_IDX  # noqa: E402

TRIALS = 20000
FATE_BAR = 0.169          # STAGE_1_CALIBRATION.md section 3
COV_LO, COV_HI = 0.85, 0.95
COV_LEVEL = 0.90


def ece(p, y, bins: int = 10) -> float:
    """Verbatim scorecard.py:112."""
    p, y = np.asarray(p, float), np.asarray(y, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (p >= edges[i]) & (p < hi)
        if m.sum():
            e += (m.sum() / len(p)) * abs(p[m].mean() - y[m].mean())
    return float(e)


def resolvability(sim_null: np.ndarray, bar: float, lower_is_better: bool = True) -> dict:
    """How the criterion behaves when the intent is satisfied EXACTLY.

    `pass_rate` is the probability a perfect system is reported as passing. Anything well below
    1.0 means the criterion mostly measures sampling noise. `usable_bar` is where the threshold
    would have to sit for a perfect system to pass 95% of the time.
    """
    ok = sim_null <= bar if lower_is_better else sim_null >= bar
    return {
        "null_median": float(np.median(sim_null)),
        "null_p90": float(np.percentile(sim_null, 90)),
        "pass_rate": float(ok.mean()),
        "usable_bar": float(np.percentile(sim_null, 95 if lower_is_better else 5)),
    }


def sensitivity_multiplier(n_folds: int) -> float:
    """What a PAIRED test's sensitivity actually depends on.

    scorecard's rule is "real only if the paired 95% CI across folds excludes 0". The CI is built
    on the paired DIFFERENCES, so the relevant spread is the fold-to-fold heterogeneity of the
    CHANGE -- not the fold-to-fold spread of the metric itself, which cancels in the pairing.

    A change that shifts every fold by the same amount has difference-SD 0 and is therefore
    detected at ANY magnitude; that is exactly why Stage 1's guards reading +0.000 with CI
    [+0.000, +0.000] is stronger evidence than any variance argument. A change that helps some
    folds and hurts others can be large in the mean and still read as noise.

    Returns k such that the minimum detectable MEAN effect is k * SD(effect across folds).
    """
    from scipy.stats import t

    return float(t.ppf(0.975, n_folds - 1) / np.sqrt(n_folds))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dump", nargs="?", default="diag_dump")
    ap.add_argument("--baseline", default="scorecard/baseline.json",
                    help="a snapshot, used only for fold-to-fold SD of the guards")
    args = ap.parse_args()
    dd = Path(args.dump)
    files = sorted(dd.glob("*.npz"))
    if not files:
        print(f"No .npz in {dd}/. Run dump_diag_bundle.py first.")
        return 1

    rng = np.random.default_rng(0)
    man = json.loads((dd / "manifest.json").read_text(encoding="utf-8"))
    data = {f.stem: np.load(f, allow_pickle=False) for f in files}

    P, Y, scored = {}, {}, []
    for n, d in data.items():
        if "test_S" not in d.files:
            continue
        P[n] = np.asarray(d["test_S"], float)
        Y[n] = (d["test_y_cls"].astype(int) == SAFE_IDX).astype(float)
        # scorecard reports n/a for a fold with no class variation -- match that exactly
        if 0 < Y[n].sum() < len(Y[n]):
            scored.append(n)
    scored.sort()

    print("=" * 84)
    print("  ABSOLUTE CRITERIA (TARGETs) -- false-negative rate when the intent holds EXACTLY")
    print("=" * 84)

    # ---- fate_ece, as graded: mean of per-fold ECEs -------------------------------- #
    sims = np.stack([[ece(P[n], (rng.random(len(P[n])) < P[n]).astype(float))
                      for _ in range(TRIALS)] for n in scored]).mean(axis=0)
    obs = float(np.mean([ece(P[n], Y[n]) for n in scored]))
    r = resolvability(sims, FATE_BAR)
    print(f"\n  fate_ece <= {FATE_BAR}   [as graded: mean of per-fold ECE, n~21 x {len(scored)} folds]")
    print(f"    perfectly calibrated model scores : median {r['null_median']:.3f}, "
          f"p90 {r['null_p90']:.3f}")
    print(f"    ** it PASSES the bar {r['pass_rate']:.1%} of the time **")
    print(f"    bar would need to be >= {r['usable_bar']:.3f} for a perfect model to pass 95%")
    print(f"    observed = {obs:.3f}  (percentile {float((sims < obs).mean()):.1%} of the null)")

    # ---- fate_ece, pooled across folds -------------------------------------------- #
    p_all = np.concatenate([P[n] for n in scored])
    y_all = np.concatenate([Y[n] for n in scored])
    sims_p = np.array([ece(p_all, (rng.random(len(p_all)) < p_all).astype(float))
                       for _ in range(TRIALS)])
    obs_p = ece(p_all, y_all)
    rp = resolvability(sims_p, FATE_BAR)
    print(f"\n  fate_ece <= {FATE_BAR}   [POOLED over all {len(p_all)} held-out cells]")
    print(f"    perfectly calibrated model scores : median {rp['null_median']:.3f}, "
          f"p90 {rp['null_p90']:.3f}")
    print(f"    ** it PASSES the bar {rp['pass_rate']:.1%} of the time **")
    print(f"    observed = {obs_p:.3f}  (percentile {float((sims_p < obs_p).mean()):.1%} of the null)")

    # ---- conformal coverage -------------------------------------------------------- #
    hits, tot, per_fold_cov = 0, 0, []
    for n, d in data.items():
        if "test_mu_age" not in d.files:
            continue
        m = np.asarray(d["test_mask"], bool)
        q = float(man[n]["conformal_q"])
        h = np.abs(np.asarray(d["test_y_age"], float)[m]
                   - np.asarray(d["test_mu_age"], float)[m]) <= q
        hits += int(h.sum()); tot += int(m.sum()); per_fold_cov.append(float(h.mean()))
    sims_c = rng.binomial(tot, COV_LEVEL, TRIALS) / tot
    in_band = ((sims_c >= COV_LO) & (sims_c <= COV_HI)).mean()
    print(f"\n  conformal_coverage in [{COV_LO}, {COV_HI}]   [POOLED marginal over {tot} cells]")
    print(f"    a correctly-{COV_LEVEL:.0%} system scores : median {np.median(sims_c):.3f}, "
          f"90% range [{np.percentile(sims_c,5):.3f}, {np.percentile(sims_c,95):.3f}]")
    print(f"    ** it PASSES the band {in_band:.1%} of the time **")
    print(f"    observed = {hits/tot:.3f}  -> conformal guarantees MARGINAL coverage, and this")
    print(f"       is the marginal rate, so the criterion matches what the method promises")

    # ---- guards -------------------------------------------------------------------- #
    print("\n" + "=" * 84)
    print("  COMPARATIVE CRITERIA (GUARDs) -- how big must a REAL regression be to be seen?")
    print("=" * 84)
    print("\n  The paired CI is built on DIFFERENCES, so a guard's sensitivity is set by how")
    print("  CONSISTENT a change is across folds -- not by the metric's own fold-to-fold spread,")
    print("  which cancels in the pairing. Minimum detectable MEAN effect:\n")
    for nf in (6, 5):
        k = sensitivity_multiplier(nf)
        print(f"    {nf} folds:  mean effect must exceed  {k:.2f} x SD(effect across folds)")
    print("\n    -> a UNIFORM change (SD 0) is caught at ANY magnitude. This is why Stage 1's")
    print("       guards reading +0.000 with CI [+0.000,+0.000] is strong evidence and not luck.")
    print("    -> a change that helps some folds and hurts others can be large in the mean")
    print("       and still be reported as noise. That is the blind spot to watch from Stage 2 on,")
    print("       where changes DO touch the model and will not be uniform.")

    print("\n  Calibrated against the one real, non-zero change measured so far")
    print("  (A_xdonor -> B_fatecal on fate_ece, 5 folds):")
    k5 = sensitivity_multiplier(5)
    obs_mean, obs_half = 0.115, (0.142 - 0.087) / 2          # from the run-3 compare table
    sd_eff = obs_half / k5
    print(f"    mean effect {obs_mean:.3f}, CI half-width {obs_half:.4f}"
          f"  -> SD(effect) = {sd_eff:.4f}")
    print(f"    so that change was {obs_mean / (k5 * sd_eff):.1f}x the detection threshold -- "
          f"comfortably real.")
    print(f"    a change with the same heterogeneity would need a mean above "
          f"{k5 * sd_eff:.3f} to register.")

    print("\n" + "=" * 84)
    print("  NOT REACHABLE FROM THIS DUMP (would need ridge predictions and RES recomputed):")
    print("    dage_mae_ridge, level_shift_ridge, rank_ridge_dage, rank_res, res_* ")
    print("  All are context metrics -- none is a TARGET or a GUARD, so no criterion depends")
    print("  on them. Auditing them needs a dump extension, not a rerun.")
    print("=" * 84)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
