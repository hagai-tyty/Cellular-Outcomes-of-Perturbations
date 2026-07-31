"""STAGE 1.5.2 gate G-c step 1 — do HFF's ΔAge labels carry the rejuvenation signature?

    python experiments/diag_gc_hff_signature.py                     # pre-register only
    python experiments/diag_gc_hff_signature.py --run [RUN_DIR]     # + the measurement

READ-ONLY. Writes `diag_gc_hff_signature_results.json`. `src/` untouched.

THE QUESTION G-c ASKS
---------------------
"Are HFF's ΔAge labels informative, or is the ΔAge head learning artefact from 33,613 contaminated
labels while 75 usable ones are drowned out?"

Step 1 is the cheap half, and it is deliberately cheap: "Spearman(ΔAge, reprogramming day) over the
HFF trajectory, with the iPSC endpoint excluded (a cell-type change, per the standing rule), and the
dose-response slope reported beside the methylation figures."

WHAT IT IS MEASURED AGAINST
---------------------------
`REV FINAL` §4.4(b) established what real rejuvenation looks like in this system, on methylation:

    slope  -3.30 / -3.15 yr per day,  rho  -0.885 / -0.842,  p <= 0.0006,  MONOTONE

G-c's registered rule then reads:

    HFF shows the signature (monotone, slope within ~2x of methylation's)  -> keep the labels
    ambiguous                                                              -> run step 2 (retrain)
    no signature                                                           -> mask HFF's ΔAge

WHERE THE LABELS COME FROM, AND WHY NOT FROM A PSEUDOBULK
---------------------------------------------------------
§0's G-c block cites `diag_d2_replication_results.json` (-0.36 yr/day, rho -0.214) — but that is a
**pseudobulk of 2000 sampled cells per timepoint**, computed on absolute predicted age. The labels
the model actually trains on are per-cell `y_age` in the built shards, control-relative, after
cell-cycle deconfounding. Those are what G-c is about, so those are what is read here. If the two
disagree, that is itself worth knowing.

UNIT OF ANALYSIS
----------------
`day` is a TIMEPOINT-level variable, so the decisive rho is computed over the **8 timepoints**, which
is also how the methylation figures it is compared against were computed. The per-cell rho over all
~33k cells is reported alongside as descriptive: at that n almost anything is "significant", and a
criterion graded there would be measuring the cell count.
"""
from __future__ import annotations

import glob
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from audit_metrics import bar_verdict  # noqa: E402

N_SIM = 20000
RNG = np.random.default_rng(0)

# REV FINAL §4.4(b) — the signature to match. Not re-derived here; quoted.
METH_SLOPES = (-3.30, -3.15)
METH_RHOS = (-0.885, -0.842)
METH_SLOPE = float(np.mean(METH_SLOPES))
METH_RHO = float(np.mean(METH_RHOS))

# "within ~2x of methylation's" made operational, both directions.
SLOPE_LO, SLOPE_HI = METH_SLOPE * 2.0, METH_SLOPE / 2.0        # -6.45 .. -1.61
# "monotone" made operational. Deliberately WEAKER than methylation achieved (-0.86): a bar set at
# methylation's own level would be graded on a different modality's precision. Failing an easier
# bar is the decisive direction here.
RHO_BAR = -0.50

IPSC_DAY = 21.0          # GSE242423SingleCellSource.ipsc_day — a cell-type change, excluded
TIME_FLOOR = 0.01        # perturbation._TIME_FLOOR, so day-0 controls have a finite log


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def day_from_dose_time(dose_time: np.ndarray) -> np.ndarray:
    """Invert `encode_dose_time`: column 1 is log(time_h), and time_h = day * 24. Pure.

    Day-0 controls were floored to `time_h = TIME_FLOOR` before the log, so they come back as
    ~0.0004 rather than exactly 0. Rounding to the nearest half-day recovers the real grid without
    inventing a special case for the floor.
    """
    time_h = np.exp(np.asarray(dose_time, float)[:, 1])
    return np.round(time_h / 24.0 * 2.0) / 2.0


def trajectory_stats(days: np.ndarray, y_age: np.ndarray, exclude_day: float = IPSC_DAY) -> dict:
    """Timepoint-level and per-cell trend of ΔAge against reprogramming day. Pure."""
    from scipy.stats import spearmanr
    d = np.asarray(days, float)
    y = np.asarray(y_age, float)
    keep = np.isfinite(d) & np.isfinite(y) & (d != exclude_day)
    d, y = d[keep], y[keep]
    if d.size < 4 or np.ptp(d) == 0:
        return {"status": "CANNOT_VERIFY", "n_cells": int(d.size)}
    uniq = np.unique(d)
    means = np.array([y[d == u].mean() for u in uniq])
    sems = np.array([y[d == u].std(ddof=1) / np.sqrt(max((d == u).sum(), 1)) for u in uniq])
    slope = float(np.polyfit(uniq, means, 1)[0])
    return {
        "status": "OK",
        "n_cells": int(d.size), "n_timepoints": int(uniq.size),
        "days": [float(x) for x in uniq],
        "mean_dage": [float(x) for x in means], "sem_dage": [float(x) for x in sems],
        "n_per_day": [int((d == u).sum()) for u in uniq],
        "rho_timepoint": float(spearmanr(uniq, means).correlation),
        "slope_yr_per_day": slope,
        "rho_percell": float(spearmanr(d, y).correlation),
        "slope_percell": float(np.polyfit(d, y, 1)[0]),
        # a monotone series has |rho| = 1; the count of descending steps says how it fails
        "n_descending_steps": int((np.diff(means) < 0).sum()),
        "n_steps": int(uniq.size - 1),
    }


def leave_one_timepoint_out(days: list[float], means: list[float]) -> dict:
    """Recompute rho and slope with each timepoint dropped in turn. Pure. DESCRIPTIVE.

    Added after the first run. Day 14 is the last point before the iPSC endpoint that the
    standing rule already excludes as a cell-type change, so a trend that lives entirely in
    day 14 would be the identity confound arriving one timepoint early. This asks that
    directly instead of arguing about it.
    """
    from scipy.stats import spearmanr
    d, m = np.asarray(days, float), np.asarray(means, float)
    out = {}
    for i in range(d.size):
        k = np.ones(d.size, bool)
        k[i] = False
        if k.sum() >= 3 and np.ptp(d[k]) > 0:
            out[f"drop_day_{d[i]:g}"] = {
                "rho": float(spearmanr(d[k], m[k]).correlation),
                "slope": float(np.polyfit(d[k], m[k], 1)[0])}
    rhos = [v["rho"] for v in out.values()]
    slopes = [v["slope"] for v in out.values()]
    return {"folds": out,
            "rho_range": [float(min(rhos)), float(max(rhos))],
            "slope_range": [float(min(slopes)), float(max(slopes))]}


def gc_verdict(stats: dict, rho_bar: float, slope_lo: float, slope_hi: float) -> dict:
    """G-c step 1's pre-registered three-way decision. Pure.

    Both criteria must hold to KEEP. Exactly one holding is AMBIGUOUS, which is what step 2's
    retrain comparison exists for. Neither holding is NO_SIGNATURE.
    """
    if stats.get("status") != "OK":
        return {"action": "CANNOT_VERIFY", "reason": "not enough usable cells"}
    rho, slope = stats["rho_timepoint"], stats["slope_yr_per_day"]
    rho_ok = rho <= rho_bar
    slope_ok = slope_lo <= slope <= slope_hi
    n_ok = int(rho_ok) + int(slope_ok)
    detail = (f"rho_timepoint {rho:+.3f} vs bar <= {rho_bar:+.2f} -> "
              f"{'PASS' if rho_ok else 'FAIL'}; slope {slope:+.3f} yr/day vs band "
              f"[{slope_lo:+.2f}, {slope_hi:+.2f}] -> {'PASS' if slope_ok else 'FAIL'}")
    if n_ok == 2:
        return {"action": "KEEP_HFF_LABELS", "detail": detail, "rho_ok": rho_ok,
                "slope_ok": slope_ok,
                "reason": "HFF's ΔAge labels carry the rejuvenation signature methylation "
                          "established. Keep them and record the check."}
    if n_ok == 1:
        return {"action": "RUN_STEP_2", "detail": detail, "rho_ok": rho_ok, "slope_ok": slope_ok,
                "reason": "one criterion holds and one does not — ambiguous, which G-c routes to "
                          "step 2's age_mask=True/False retrain comparison on the scorecard."}
    return {"action": "MASK_HFF_IN_PHASE_2", "detail": detail, "rho_ok": rho_ok,
            "slope_ok": slope_ok,
            "reason": "HFF's ΔAge labels show NO trace of the signature. G-c: mask them, and state "
                      "the consequence plainly — the age head then trains on ~75 labels, which may "
                      "be too few. That is a finding, not a failure."}


def sim_signature(n_timepoints: int, rho_true: float = METH_RHO, n_sim: int = N_SIM):
    """Timepoint-level Spearman for a label set that DOES carry the methylation signature."""
    from scipy.stats import spearmanr
    days = np.arange(n_timepoints, dtype=float)
    z = (days - days.mean()) / days.std()
    out = np.empty(n_sim)
    sd_noise = np.sqrt(1.0 / rho_true**2 - 1.0)
    for i in range(n_sim):
        y = rho_true * z + RNG.normal(0.0, sd_noise, size=n_timepoints)
        out[i] = spearmanr(days, y).correlation
    return out


# --------------------------------------------------------------------------- #
# Data wiring                                                                  #
# --------------------------------------------------------------------------- #
def load_hff(run_dir: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Per-cell (day, y_age) for the HFF line, plus a label-volume census. Reads built shards."""
    import pandas as pd
    files = sorted(glob.glob(str(run_dir / "shards" / "*.parquet")))
    if not files:
        raise SystemExit(f"no shards under {run_dir / 'shards'}")
    days, ages, volume = [], [], {"by_line": {}, "age_masked_by_line": {}}
    for f in files:
        d = pd.read_parquet(f, columns=["cell_line", "dose_time", "y_age", "age_mask"])
        lines = d["cell_line"].to_numpy()
        for ln in np.unique(lines):
            m = lines == ln
            volume["by_line"][str(ln)] = volume["by_line"].get(str(ln), 0) + int(m.sum())
            volume["age_masked_by_line"][str(ln)] = (
                volume["age_masked_by_line"].get(str(ln), 0)
                + int(d["age_mask"].to_numpy()[m].sum()))
        m = lines == "HFF"
        if not m.any():
            continue
        dt = np.stack(d["dose_time"].to_numpy()[m])
        days.append(day_from_dose_time(dt))
        ages.append(d["y_age"].to_numpy()[m].astype(float))
    if not days:
        raise SystemExit("no HFF cells found in these shards")
    return np.concatenate(days), np.concatenate(ages), volume


def main() -> int:
    print("STAGE 1.5.2 G-c step 1 — do HFF's ΔAge labels carry the rejuvenation signature?\n")
    print("  PHASE 1: freeze the bars. No label is read in this phase.\n")
    print(f"  the signature to match (REV FINAL §4.4(b)): slope {METH_SLOPES} yr/day, "
          f"rho {METH_RHOS}")
    print(f"  operationalised: rho_timepoint <= {RHO_BAR:+.2f} AND slope in "
          f"[{SLOPE_LO:+.2f}, {SLOPE_HI:+.2f}] yr/day\n")

    sim = sim_signature(8)
    v = bar_verdict(sim, RHO_BAR, lower_is_better=True)
    print(f"  rho bar <= {RHO_BAR:+.2f} at 8 timepoints, for a system that DOES carry the "
          f"signature (rho_true {METH_RHO:+.3f}):")
    print(f"     pass rate {v['pass_rate']:.1%}  -> {v['verdict']}")
    rho_bar = RHO_BAR
    if v["verdict"] != "RESOLVABLE":
        rho_bar = float(v["usable_bar"])
        print(f"     [!] §5b: the bar moves NOW, before the labels are opened, to "
              f"{rho_bar:+.3f}")
    out: dict = {"script": "diag_gc_hff_signature",
                 "utc": datetime.now(UTC).isoformat(timespec="seconds"),
                 "methylation_reference": {"slopes": METH_SLOPES, "rhos": METH_RHOS},
                 "preregistration": {"rho_bar_proposed": RHO_BAR, "rho_bar_used": rho_bar,
                                     "slope_band": [SLOPE_LO, SLOPE_HI], **v}}

    if "--run" not in sys.argv:
        print("\n  Pre-registration only. Re-run with --run to measure.")
        Path("diag_gc_hff_signature_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 0

    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    run_dir = Path(pos[0] if pos else "runs/cellfate_loocv_O1")
    print(f"\n  PHASE 2: measurement, from the built labels in {run_dir}\n")
    days, ages, volume = load_hff(run_dir)
    st = trajectory_stats(days, ages)
    out["label_volume"] = volume
    out["hff_trajectory"] = st

    total = sum(volume["age_masked_by_line"].values())
    hff = volume["age_masked_by_line"].get("HFF", 0)
    print(f"  age-valid labels: {total} total, HFF {hff} ({hff/max(total,1):.2%}), "
          f"non-HFF {total - hff}")
    print(f"  HFF cells used: {st['n_cells']} across {st['n_timepoints']} timepoints "
          f"(iPSC day {IPSC_DAY:.0f} excluded)\n")
    print(f"  {'day':>5}{'n cells':>10}{'mean ΔAge':>12}{'SEM':>8}")
    for d, n, m, s in zip(st["days"], st["n_per_day"], st["mean_dage"], st["sem_dage"],
                          strict=True):
        print(f"  {d:>5.1f}{n:>10d}{m:>12.3f}{s:>8.3f}")

    print(f"\n  rho_timepoint  {st['rho_timepoint']:+.3f}   (bar <= {rho_bar:+.3f})   "
          f"methylation: {METH_RHOS}")
    print(f"  slope          {st['slope_yr_per_day']:+.3f} yr/day   "
          f"(band [{SLOPE_LO:+.2f}, {SLOPE_HI:+.2f}])   methylation: {METH_SLOPES}")
    print(f"  descending steps {st['n_descending_steps']}/{st['n_steps']}  "
          f"(monotone rejuvenation would be {st['n_steps']}/{st['n_steps']})")
    print(f"  [descriptive] per-cell rho {st['rho_percell']:+.3f}, slope "
          f"{st['slope_percell']:+.3f} yr/day over {st['n_cells']} cells")

    loo = leave_one_timepoint_out(st["days"], st["mean_dage"])
    out["leave_one_timepoint_out"] = loo
    print("\n  [descriptive] leave-one-timepoint-out — is the trend carried by one point?")
    for k, v in loo["folds"].items():
        print(f"     {k:<16} rho {v['rho']:+.3f}   slope {v['slope']:+.3f}")
    print(f"     rho spans [{loo['rho_range'][0]:+.3f}, {loo['rho_range'][1]:+.3f}]; "
          f"slope spans [{loo['slope_range'][0]:+.3f}, {loo['slope_range'][1]:+.3f}]")

    dec = gc_verdict(st, rho_bar, SLOPE_LO, SLOPE_HI)
    out["verdict"] = dec
    print(f"\n  ==> G-c STEP 1: {dec['action']}\n      {dec.get('detail', '')}\n"
          f"      {dec['reason']}")
    if dec["action"] == "MASK_HFF_IN_PHASE_2":
        print(f"\n  CONSEQUENCE, stated plainly: masking HFF leaves {total - hff} age labels.")

    Path("diag_gc_hff_signature_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_gc_hff_signature_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
