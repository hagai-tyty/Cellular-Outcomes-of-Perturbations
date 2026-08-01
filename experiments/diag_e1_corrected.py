"""E1/E1b RE-RUN with the A1 and A3 errors fixed — using the CURRENT day-0 labels.

    python experiments/diag_e1_corrected.py "D:\\Gill"

READ-ONLY. Writes `diag_e1_corrected_results.json`. No labels change, no rebuild, no GPU.
`src/` is untouched.

WHY (STAGE_1_5_1_REVISED §2 A1/A3, and its review §4)
-----------------------------------------------------
`STAGE_1_5_1_REVISED.md` identified three errors. Two are **verified and undisputed**:

  A1  `diag_e1_trajectory.py` excluded only iPSC, so the **47 "Failing to reprogram" samples were
      pooled into the treatment arm** — 42% of samples that by definition cannot rejuvenate.
  A3  E1/E1b used a **monotonic Spearman** on an effect that **dips and recovers** (~day 13). A
      fall-then-rise rank-correlates to ~0.

Its third proposal — redefining `is_control` to the non-responder arm — is **disputed** (review §2:
the arms are not identity-matched; 61% of that effect is identity; it contradicts Gill's own
non-responder direction). **This script deliberately avoids that dispute.** It fixes ONLY A1 and A3
and keeps the existing day-0 control, so it answers one clean question:

    Once we stop pooling non-responders and stop using the wrong statistic,
    do the CURRENT ΔAge labels show rejuvenation in responders?

That question is decisive either way and costs nothing:
  * **YES** -> the Stage 1.5 escalation was purely a method artefact; no label change is needed.
  * **NO**  -> the current label definition really is the problem, and the identity-anchor question
               must be faced directly rather than settled by choosing a control.

PRE-REGISTRATION HONESTY
------------------------
The peak-window responder mean (+8.2 yr) was already computed during the review, so it is **not**
blind. Everything else here — the per-window profile, the shape test, the LOO check and the
non-responder comparison — is new. This is stated so the record is not overclaimed.
"""

from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested; no repo-data imports below this block.  #
# --------------------------------------------------------------------------- #
RESPONDER = "Reprogramming fibroblast"
NON_RESPONDER = "Failing to reprogram fibroblast"
BASELINE = "Dermal fibroblast"
IPSC = "iPSC"

# Gill's reported optimum. Pre-committed BEFORE looking at other windows, and taken from the
# protocol (OSKM withdrawn ~day 13), not tuned on our numbers.
PEAK_WINDOW = (10.0, 13.0)

# Windows reported as sensitivity. Fixed here so the set cannot grow until something passes.
SENSITIVITY_WINDOWS = [(7.0, 9.0), (10.0, 13.0), (13.0, 15.0), (15.0, 21.0), (21.0, 29.0)]

T_CRIT = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365}


def paired_ci(values: list[float], conf: float = 0.95) -> dict:
    """Mean and paired 95% CI over donors. Pure. The project's standard treatment."""
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    n = len(v)
    if n < 2:
        return {"n": n, "mean": float(v.mean()) if n else float("nan"),
                "ci95": [float("nan"), float("nan")], "n_negative": int((v < 0).sum())}
    mean = float(v.mean())
    se = float(v.std(ddof=1)) / np.sqrt(n)
    t = T_CRIT.get(n - 1, 1.96)
    return {"n": n, "mean": mean, "ci95": [mean - t * se, mean + t * se],
            "n_negative": int((v < 0).sum()), "sd": float(v.std(ddof=1))}


def window_verdict(stat: dict, fragile_margin: float = 0.5) -> dict:
    """Verdict for a window contrast on ΔAge (negative = rejuvenation).

    FRAGILE is reported whenever a CI bound sits within `fragile_margin` years of zero — the §10
    lesson, where E1b (+0.009) and D2 (−0.014) were both decided by hundredths.
    """
    if not np.isfinite(stat.get("mean", float("nan"))) or stat["n"] < 2:
        return {"status": "CANNOT_VERIFY", "reason": f"n={stat.get('n', 0)} donors"}
    lo, hi = stat["ci95"]
    fragile = min(abs(lo), abs(hi)) < fragile_margin
    if hi < 0:
        s, r = "REJUVENATION", "dAge falls significantly in this window"
    elif lo > 0:
        s, r = "AGEING", "dAge RISES significantly in this window"
    else:
        s, r = "NO_EFFECT", "CI includes 0"
    out = {"status": s, "mean_years": stat["mean"], "ci95": stat["ci95"],
           "n_donors": stat["n"], "n_negative": stat["n_negative"],
           "reason": f"{r} (mean {stat['mean']:+.1f} yr, 95% CI [{lo:+.1f}, {hi:+.1f}], "
                     f"{stat['n_negative']}/{stat['n']} donors negative)"}
    if fragile:
        out["status"] = s + "_FRAGILE"
        out["reason"] += "  [FRAGILE]: a CI bound is within 0.5 yr of zero"
    return out


def leave_one_donor_out(per_donor: dict[str, float]) -> dict:
    """Does the sign survive dropping any single donor? Guards against one donor carrying it."""
    donors = [d for d, v in per_donor.items() if np.isfinite(v)]
    if len(donors) < 3:
        return {"status": "CANNOT_VERIFY", "reason": "need >=3 donors"}
    means = {d: float(np.mean([per_donor[o] for o in donors if o != d])) for d in donors}
    signs = {np.sign(m) for m in means.values()}
    return {"status": "STABLE" if len(signs) == 1 else "SIGN_FLIPS",
            "loo_means": means,
            "reason": ("sign is the same for all leave-one-out subsets"
                       if len(signs) == 1 else
                       "the sign FLIPS when some donor is dropped — not robust at this n")}


def decide(peak: dict, windows: dict, loo: dict) -> dict:
    """What the corrected re-run licenses. Pure."""
    st = peak.get("status", "")
    if st.startswith("REJUVENATION"):
        return {"action": "LABELS_ARE_FINE",
                "reason": "with A1/A3 fixed, the CURRENT day-0 labels show rejuvenation in "
                          "responders. The Stage 1.5 escalation was a method artefact. No label "
                          "change is needed — do NOT redefine is_control."}
    any_rejuv = [w for w, v in windows.items() if v.get("status", "").startswith("REJUVENATION")]
    if st.startswith("AGEING"):
        return {"action": "LABELS_INADEQUATE",
                "reason": "with A1/A3 fixed, responders still read OLDER against their own day-0 "
                          "baseline. Fixing the pooling and the statistic does NOT rescue the "
                          "current labels, so the escalation is not purely a method artefact. The "
                          "identity-anchor question must be faced directly — and note this is ALSO "
                          "evidence against the day-0 control, independent of the non-responder "
                          "proposal."
                          + (f" (Some windows do show rejuvenation: {any_rejuv} — report, do not "
                             "select on them.)" if any_rejuv else "")}
    return {"action": "INCONCLUSIVE",
            "reason": "no significant movement in the pre-registered window. "
                      + (f"Other windows moved: {any_rejuv}. " if any_rejuv else "")
                      + ("Leave-one-out is unstable, so n=6 is the binding constraint."
                         if loo.get("status") == "SIGN_FLIPS" else
                         "Consistent across leave-one-out, so this is a real null, not one donor.")}


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def load_ages(gill_dir: str):
    """(ages, cell type, day, donor) per sample via the production path + frozen clock."""
    import gzip
    root = Path(__file__).resolve().parents[1]
    for p in (root, root / "local_runners", root / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    from run_multi_local import discover_gill  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GillReprogrammingSource

    expr, series = discover_gill(gill_dir)
    src = GillReprogrammingSource(expr_tsv=expr, series_matrix=series)
    src._load()
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")

    rows: dict[str, list[str]] = {}
    titles: list[str] = []
    with gzip.open(series, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                rows[v[0].split(":")[0].strip()] = v
    ct = {titles[i]: rows["cell type"][i].split(":", 1)[1].strip() for i in range(len(titles))}
    day = {titles[i]: float(rows["days of reprogramming"][i].split(":", 1)[1].strip())
           for i in range(len(titles))}
    donor = {t: t.split("_")[0] for t in titles}

    cols = [c for c in src._rpm.columns if c in ct]
    mat = src._rpm[cols].to_numpy(dtype=np.float64).T
    ages = dict(zip(cols, clock.predict_age(normalize_counts(mat), list(src._genes)), strict=True))
    return ages, ct, day, donor, cols


def delta_by_donor(ages, ct, day, donor, cols, arm: str, window: tuple[float, float]):
    """Per-donor ΔAge = mean(arm samples in window) − that donor's day-0 baseline.

    This is EXACTLY the current label definition (`sources.py:471`: control = day 0), restricted to
    one arm. No control redefinition.
    """
    out: dict[str, float] = {}
    for d in sorted({donor[c] for c in cols}):
        base = [ages[c] for c in cols if donor[c] == d and ct[c] == BASELINE]
        sel = [ages[c] for c in cols
               if donor[c] == d and ct[c] == arm and window[0] <= day[c] <= window[1]]
        if base and sel:
            out[d] = float(np.mean(sel) - np.mean(base))
    return out


def main() -> int:
    gill = sys.argv[1] if len(sys.argv) > 1 else r"D:\Gill"
    print("E1/E1b RE-RUN -- A1 (pooling) and A3 (statistic) fixed, CURRENT day-0 labels\n")
    print(f"  pre-registered window: days {PEAK_WINDOW[0]:.0f}-{PEAK_WINDOW[1]:.0f} "
          "(Gill's optimum / OSKM withdrawal)")
    print("  arm: RESPONDERS only ('Reprogramming fibroblast'); non-responders EXCLUDED (A1)")
    print("  statistic: window contrast vs the donor's own day-0 baseline, NOT a monotonic "
          "Spearman (A3)\n")
    try:
        ages, ct, day, donor, cols = load_ages(gill)
    except Exception as exc:  # noqa: BLE001
        print(f"  [!] could not load Gill data: {exc!r}")
        return 1

    resp = delta_by_donor(ages, ct, day, donor, cols, RESPONDER, PEAK_WINDOW)
    peak_stat = paired_ci(list(resp.values()))
    peak = window_verdict(peak_stat)
    loo = leave_one_donor_out(resp)

    print("  PER-DONOR dAge (responders, pre-registered window):")
    for d, v in sorted(resp.items()):
        print(f"     {d:<4} {v:+8.1f} yr")
    print(f"\n  PRE-REGISTERED WINDOW: {peak['status']}\n     {peak['reason']}")
    print(f"  leave-one-donor-out: {loo['status']} -- {loo['reason']}")

    print("\n  SENSITIVITY -- every window, responders only (reported, not selected on):")
    print(f"     {'window':<12}{'n':>3}{'mean dAge':>12}{'95% CI':>22}  verdict")
    windows = {}
    for w in SENSITIVITY_WINDOWS:
        dd = delta_by_donor(ages, ct, day, donor, cols, RESPONDER, w)
        v = window_verdict(paired_ci(list(dd.values())))
        windows[f"{w[0]:.0f}-{w[1]:.0f}"] = v
        if v["status"] == "CANNOT_VERIFY":
            print(f"     {w[0]:.0f}-{w[1]:.0f}d       n/a")
            continue
        lo, hi = v["ci95"]
        print(f"     {w[0]:.0f}-{w[1]:.0f}d{'':<6}{v['n_donors']:>3}{v['mean_years']:>12.1f}"
              f"   [{lo:>+7.1f},{hi:>+7.1f}]  {v['status']}")

    print("\n  COMPARISON -- the same windows for NON-responders (context; not used as a control):")
    for w in SENSITIVITY_WINDOWS:
        dd = delta_by_donor(ages, ct, day, donor, cols, NON_RESPONDER, w)
        v = window_verdict(paired_ci(list(dd.values())))
        if v["status"] == "CANNOT_VERIFY":
            continue
        print(f"     {w[0]:.0f}-{w[1]:.0f}d{'':<6}{v['n_donors']:>3}{v['mean_years']:>12.1f}"
              f"  {v['status']}")

    decision = decide(peak, windows, loo)
    print(f"\n  ==> ACTION: {decision['action']}\n      {decision['reason']}")

    out = {"script": "diag_e1_corrected", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "gill_dir": gill, "peak_window": list(PEAK_WINDOW),
           "responder_delta_by_donor": resp, "peak_verdict": peak,
           "leave_one_donor_out": loo, "windows_responders": windows, "decision": decision}
    (_RESULTS / "diag_e1_corrected_results.json").write_text(json.dumps(out, indent=2, default=str),
                                                      encoding="utf-8")
    print("\n  wrote diag_e1_corrected_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
