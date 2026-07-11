"""
Test 5.0 (ΔAge lab notebook) — is the model's ΔAge deficit vs ridge REAL, or noise?

Aggregate LOOCV showed model 14.29 vs ridge 14.05 — the model is WORSE by 0.24. But the
per-fold spread is large (std ~8-9). This script does the honest PAIRED test across the 6
leave-one-donor-out folds: per-fold (model - ridge) difference, its mean, a paired t-CI,
and the win/loss count. If the difference is indistinguishable from 0, "0.24 worse" is
noise and the honest statement is "statistically tied".

Reads each fold's reports/holdout.json (already on disk). Run once.

USAGE (repo root, venv active):
    python test5_ridge_gap.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
# 95% two-sided t critical values by df (df = n-1)
T_CRIT = {5: 2.571, 4: 2.776, 3: 3.182, 2: 4.303, 1: 12.706}


def find_report(donor: str) -> Path | None:
    for base in (".", "runs", ".."):
        p = Path(base) / f"cellfate_loocv_{donor}" / "reports" / "holdout.json"
        if p.exists():
            return p
    return None


def main() -> None:
    try:
        from cellfate.common.console import install_pretty_console, render_table
        install_pretty_console()
    except Exception:  # noqa: BLE001
        def render_table(h, r, **k):  # minimal fallback
            return "\n".join(["  ".join(map(str, h))] + ["  ".join(map(str, x)) for x in r])

    print("\nTEST 5.0 — is the model's ΔAge deficit vs ridge REAL or noise? (paired, 6 folds)")

    rows, diffs = [], []
    for d in DONORS:
        rp = find_report(d)
        if rp is None:
            rows.append([d, "  n/a", "  n/a", "  n/a", "report not found"])
            continue
        R = json.loads(rp.read_text())
        m = R.get("model", {}).get("reg_mae", float("nan"))
        r = R.get("ridge", {}).get("reg_mae", float("nan"))
        if not (math.isfinite(m) and math.isfinite(r)):
            rows.append([d, f"{m:.2f}" if math.isfinite(m) else "n/a",
                         f"{r:.2f}" if math.isfinite(r) else "n/a", "  n/a", "missing reg_mae"])
            continue
        diff = m - r
        diffs.append(diff)
        rows.append([d, f"{m:.2f}", f"{r:.2f}", f"{diff:+.2f}",
                     "model better" if diff < 0 else "ridge better"])

    print("\n" + render_table(
        ["fold (held-out)", "model MAE", "ridge MAE", "model−ridge", "who wins"],
        rows, aligns=["l", "r", "r", "r", "l"]))

    if len(diffs) < 2:
        print("\n   Not enough folds with both values to run the paired test.")
        print("   (Make sure the cellfate_loocv_* folders — or runs/cellfate_loocv_* — are present.)")
        return

    n = len(diffs)
    mean_d = sum(diffs) / n
    var_d = sum((x - mean_d) ** 2 for x in diffs) / (n - 1)
    sd = math.sqrt(var_d)
    se = sd / math.sqrt(n)
    tcrit = T_CRIT.get(n - 1, 2.571)
    ci_lo, ci_hi = mean_d - tcrit * se, mean_d + tcrit * se
    t_stat = mean_d / se if se > 0 else float("inf")
    wins = sum(1 for x in diffs if x < 0)      # model better
    losses = sum(1 for x in diffs if x > 0)

    print(f"\n   mean(model − ridge) = {mean_d:+.2f}   (negative = model better)")
    print(f"   std of per-fold diffs = {sd:.2f}   95% CI = [{ci_lo:+.2f}, {ci_hi:+.2f}]   "
          f"paired t = {t_stat:+.2f}")
    print(f"   model wins {wins}/{n} folds, ridge wins {losses}/{n}")

    print("\n   WHAT THIS MEANS:")
    significant = (ci_lo > 0) or (ci_hi < 0)
    if not significant:
        print("     -> The 95% CI INCLUDES 0: the model–ridge difference is NOT statistically")
        print("        distinguishable from zero. The '0.24 worse' is NOISE. Honest statement:")
        print("        the model is STATISTICALLY TIED with ridge on ΔAge MAE (as expected for")
        print("        a linear target). No systematic deficit to explain — Test 5.1 not needed.")
    elif ci_hi < 0:
        print("     -> CI is entirely BELOW 0: the model is SYSTEMATICALLY BETTER than ridge.")
        print("        (Would contradict 'linear target' — worth understanding, but not a deficit.)")
    else:
        print("     -> CI is entirely ABOVE 0: the model is SYSTEMATICALLY WORSE than ridge.")
        print("        This is a REAL deficit. A more-powerful model should not lose to linear on")
        print("        a linear target -> proceed to Test 5.1 (multi-task tradeoff?).")


if __name__ == "__main__":
    main()
