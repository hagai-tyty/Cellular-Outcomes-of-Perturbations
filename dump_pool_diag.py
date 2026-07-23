"""Read-only: print the cross-donor-pool calibrator diagnostics run 3 computed but never showed.

`train_model.py:225-246` fits a pool-only Platt alongside the shipped one and scores BOTH on the
same held-out pool, writing the result to each fold's `bundle/metrics.json`. The console only ever
printed the slopes. This reads those files back -- no training, no inference, no writes.

    python dump_pool_diag.py

The question it answers: would the calibrator fitted on the deployment regime alone have beaten
the union fit that shipped? Run 3 missed `fate_ece` (0.249 vs the <=0.169 bar) with a union fit
that was 97.7% in-distribution, so this is the number that decides whether reverting to the
cross-donor principle is worth a snapshot.

READ THE CAVEAT IN THE OUTPUT. `xdonor_only_safe_ece_insample` is scored on the data it was fitted
to and `shipped_safe_ece_on_pool` is not, so the comparison is biased toward the pool-only fit.
It bounds the gain from above; it does not measure it.
"""

from __future__ import annotations

import json
from pathlib import Path

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
RUNS = Path("runs")

# Smallest in-sample ECE advantage that justifies spending a snapshot on the pool-only fit.
# Set against the estimator floor measured for run 3: a perfectly calibrated model scores ~0.078
# at this geometry with a 90% range of [0.057, 0.105], so a gain under ~0.02 is inside the noise
# of the thing being measured -- and it is an in-sample gain, which shrinks out-of-sample.
MIN_GAIN = 0.02

KEYS = (
    "platt_a", "platt_b",
    "xdonor_only_platt_a", "xdonor_only_platt_b",
    "shipped_safe_ece_on_pool", "xdonor_only_safe_ece_insample",
    "xdonor_safe_ece_before", "xdonor_safe_ece_after",
    "fate_calib_n",   # nested: {"total", "in_dist", "xdonor"} -- see train_model.py:205
    "xdonor_only_n", "xdonor_only_n_donors",
)


def read_fold(donor: str) -> dict:
    """Load one fold's metrics; missing file/keys are reported, never raised."""
    p = RUNS / f"cellfate_loocv_{donor}" / "bundle" / "metrics.json"
    if not p.exists():
        return {"_error": f"no metrics.json at {p}"}
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return {"_error": f"unreadable: {exc!r}"[:100]}
    return {k: d.get(k) for k in KEYS}


def verdict(shipped: list[float], pool: list[float]) -> tuple[str, str]:
    """Decide what the pool-only diagnostic implies. Pure, so every branch is testable.

    Returns (verdict, reason). The in-sample bias runs one way only -- it flatters the pool-only
    fit -- so a pool-only number that is WORSE is conclusive, while one that is better is not.
    """
    if not shipped or not pool:
        return "NO DATA", "one or both calibrators produced no pooled ECE"
    ms, mp = sum(shipped) / len(shipped), sum(pool) / len(pool)
    if mp >= ms:
        return ("CONCLUSIVE — do not revert",
                f"pool-only ({mp:.3f}) is no better than shipped ({ms:.3f}) even WITH the "
                f"in-sample advantage; the union fit is not what cost the target")
    gain = ms - mp
    # "at least MIN_GAIN" -- the tolerance is load-bearing, not decoration. Differences of
    # decimal literals do not land on the decimal value: 0.250 - 0.230 == 0.019999999999999990,
    # which without the slack reports a gain of exactly 0.02 as though it were below 0.02.
    if gain < MIN_GAIN - 1e-9:
        return ("WEAK — not worth a snapshot",
                f"pool-only better by only {gain:.3f} in-sample, which is an upper bound; "
                f"out-of-sample it would be smaller still")
    return ("WORTH TESTING — but in-sample",
            f"pool-only better by {gain:.3f} in-sample. That is an UPPER BOUND on the real gain. "
            f"Needs an honest leave-one-donor-out-within-pool refit before anything is shipped")


def main() -> int:
    rows = {d: read_fold(d) for d in DONORS}
    ok = {d: r for d, r in rows.items() if "_error" not in r}
    for d, r in rows.items():
        if "_error" in r:
            print(f"  [!] {d}: {r['_error']}")
    if not ok:
        print("\nNo fold metrics found. Run this from the repo root on the machine that trained "
              "(expects runs/cellfate_loocv_<DONOR>/bundle/metrics.json).")
        return 1

    def col(key):
        return [r[key] for r in ok.values() if isinstance(r.get(key), (int, float))]

    # ASCII only: this prints on a Windows console that is cp1252 unless PYTHONUTF8=1 is set.
    print(f"\n  POOL CALIBRATOR DIAGNOSTICS -- {len(ok)} folds  (read from bundle/metrics.json)\n")
    hdr = (f"{'fold':<6}{'ship a':>9}{'ship b':>9}{'pool a':>9}{'pool b':>9}"
           f"{'ECE ship':>11}{'ECE pool':>11}{'safe pre':>11}{'safe post':>11}")
    print(hdr)
    print("-" * len(hdr))
    for d, r in ok.items():
        def f(k, w=9, p=3):
            v = r.get(k)
            return f"{v:>{w}.{p}f}" if isinstance(v, (int, float)) else f"{'n/a':>{w}}"
        print(f"{d:<6}{f('platt_a')}{f('platt_b')}{f('xdonor_only_platt_a')}"
              f"{f('xdonor_only_platt_b')}{f('shipped_safe_ece_on_pool', 11)}"
              f"{f('xdonor_only_safe_ece_insample', 11)}{f('xdonor_safe_ece_before', 11)}"
              f"{f('xdonor_safe_ece_after', 11)}")

    ship, pool = col("shipped_safe_ece_on_pool"), col("xdonor_only_safe_ece_insample")
    pre, post = col("xdonor_safe_ece_before"), col("xdonor_safe_ece_after")
    print("-" * len(hdr))
    if ship and pool:
        print(f"{'mean':<6}{'':>36}{sum(ship)/len(ship):>11.3f}{sum(pool)/len(pool):>11.3f}", end="")
        print(f"{sum(pre)/len(pre):>11.3f}" if pre else f"{'n/a':>11}", end="")
        print(f"{sum(post)/len(post):>11.3f}" if post else f"{'n/a':>11}")

    n = next(iter(ok.values()))
    fcn = n.get("fate_calib_n") or {}
    tot, xd = fcn.get("total"), fcn.get("xdonor")
    frac = f"  ({xd / tot:.2%} of the fit)" if isinstance(tot, int) and tot and isinstance(xd, int) else ""
    print(f"\n  fit composition: total={tot} in_dist={fcn.get('in_dist')} xdonor={xd}{frac}"
          f"\n  pool: n={n.get('xdonor_only_n')} donors={n.get('xdonor_only_n_donors')}")

    v, why = verdict(ship, pool)
    print(f"\n  VERDICT: {v}\n    {why}")
    print("\n  CAVEAT (load-bearing): 'ECE pool' is scored IN-SAMPLE on the pool it was fitted to;")
    print("  'ECE ship' is not. The comparison is biased toward the pool-only fit, so treat any")
    print("  advantage it shows as an upper bound. 'safe pre/post' are the shipped calibrator's")
    print("  binary P(safe) ECE on the pool, before and after -- these ARE comparable to each other.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
