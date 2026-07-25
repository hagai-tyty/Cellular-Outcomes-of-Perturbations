"""STAGE 1.5 §8.3 E1 — does the clock track age CHANGE within a donor's reprogramming trajectory?

    python experiments/diag_e1_trajectory.py                 # defaults to D:\\Gill
    python experiments/diag_e1_trajectory.py "D:\\Gill"

READ-ONLY. Writes `diag_e1_trajectory_results.json`. Nothing is rebuilt, refitted or written back;
`src/` is not touched.

WHY (STAGE_1_5_HARMONIZATION_AUDIT.md §8).
Phase 1 M1 showed the frozen clock does not read ABSOLUTE chronological age on this data. But the
model's target is ΔAge = w·(x_pert − x_base), CONTROL-RELATIVE — so the clock's intercept, any
additive per-donor baseline offset, and every gene Gill is missing all CANCEL. M1's absolute-age
failure therefore does NOT establish that ΔAge is invalid. E1 tests the quantity that actually
bears on ΔAge: does predicted age move DOWN as OSKM reprogramming proceeds (day 0 → 54), the
rejuvenation direction, WITHIN each donor?

  primary statistic : per-donor Spearman(predicted_age, day) over that donor's fibroblast samples
  aggregate         : mean of the six per-donor rho, paired 95% CI (t, n=6)
  PASS              : the CI excludes 0 AND is negative -> the clock tracks the rejuvenation axis
  WRONG_DIRECTION   : CI excludes 0 and is positive -> age RISES with reprogramming (reads backwards)
  NO_TREND          : CI includes 0 -> not distinguishable from a clock that reads nothing

CONFOUNDS, pre-registered:
  * iPSC endpoints are a CELL-TYPE change, not aging; a fibroblast clock is out of domain there.
    They are EXCLUDED from the primary and reported as a sensitivity (with-iPSC) alongside.
  * The neonatal donors N2/N3 failed M1 (age 0 is below the clock's ~1–94 yr fitted range). Their
    per-donor rho is reported separately; the aggregate is also reported adults-only.
  * Direction assumes reprogramming rejuvenates (the project premise). A clock that reads a real
    aging axis should show it; a positive or null trend is evidence it does not, on this data.

BAR RESOLVABILITY (ground rule §5b): the aggregate is a paired test over 6 donors, so like M2 its
power is CONDITIONAL on how CONSISTENT the trend is across donors. `e1_power()` reports P(pass) for
an assumed true rho; a moderate, consistent trend (rho ≈ −0.6, donor SD ≈ 0.25) passes ~99%. A real
but heterogeneous trend can be large in the mean and still read NO_TREND at n=6 — reported, never
hidden.
"""
from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Pure logic — data-free, fully unit-tested; nothing below imports repo data.  #
# --------------------------------------------------------------------------- #
DONOR_AGE: dict[str, float] = {"N2": 0.0, "N3": 0.0, "Y1": 29.0, "Y2": 35.0, "O1": 53.0, "O2": 53.0}
NEONATAL = ("N2", "N3")
MIN_POINTS_PER_DONOR = 4          # below this a per-donor rank correlation is not worth stating

# E1b (§8.5 follow-up): the REPROGRAMMING-PHASE-only cutoff. Gill 2022 is maturation-phase transient
# reprogramming (MPTR) -- OSKM is withdrawn ~day 13, after which cells re-differentiate, so the age
# trajectory is non-monotonic and E1's monotonic 0->54 Spearman conflates the dip with the recovery.
# The cutoff is set at day 15 to cover the withdrawal window, and is chosen from the PROTOCOL and the
# dense-sampling break between day 15 and day 21 -- NOT tuned on the phase-restricted ages. One
# pre-committed cutoff, one verdict; alternative cutoffs are deliberately not scanned (that is the
# fishing this guards against).
REPROG_PHASE_DAY_MAX = 15.0


def donor_trend(days: list[float], ages: list[float]) -> float:
    """Spearman(age, day) for one donor. Pure. NaN when undefined (too few points, or no spread)."""
    d = np.asarray(days, float)
    a = np.asarray(ages, float)
    ok = np.isfinite(d) & np.isfinite(a)
    d, a = d[ok], a[ok]
    if len(d) < MIN_POINTS_PER_DONOR or np.ptp(d) == 0 or np.ptp(a) == 0:
        return float("nan")
    from scipy.stats import spearmanr
    return float(spearmanr(d, a).correlation)


def restrict_to_phase(days: list[float], ages: list[float], day_max: float) -> tuple[list, list]:
    """Keep only samples with `day <= day_max` (the E1b reprogramming-phase window). Pure."""
    d = np.asarray(days, float)
    a = np.asarray(ages, float)
    m = d <= day_max
    return d[m].tolist(), a[m].tolist()


def e1_power(true_rho: float, donor_sd: float = 0.25, n: int = 6, alpha: float = 0.05) -> float:
    """P(pass) for a clock whose per-donor rho is truly `true_rho` with donor-to-donor SD `donor_sd`
    — the §5b resolvability figure. Two-sided paired-t power against mean = 0."""
    from math import erf, sqrt
    se = donor_sd / sqrt(n)
    t_crit = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447}.get(n - 1, 1.96)
    z = (abs(true_rho) - t_crit * se) / se
    return float(0.5 * (1.0 + erf(z / sqrt(2.0))))


def e1_verdict(rhos_by_donor: dict[str, float]) -> dict:
    """Aggregate the six per-donor rho into the one verdict. Pure; no I/O."""
    finite = {d: r for d, r in rhos_by_donor.items() if np.isfinite(r)}
    n = len(finite)
    if n < 3:
        return {"status": "CANNOT_VERIFY", "n_donors": n,
                "reason": f"only {n} donors with a defined trend; a paired CI needs >=3"}
    vals = np.array(list(finite.values()), float)
    mean = float(vals.mean())
    sd = float(vals.std(ddof=1))
    se = sd / np.sqrt(n)
    t_crit = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447}.get(n - 1, 1.96)
    lo, hi = mean - t_crit * se, mean + t_crit * se
    n_neg = int((vals < 0).sum())
    if hi < 0:
        status, reason = "PASS", "predicted age falls with reprogramming day; the clock tracks the rejuvenation axis"
    elif lo > 0:
        status, reason = "WRONG_DIRECTION", "predicted age RISES with reprogramming day; the clock reads the axis backwards"
    else:
        status, reason = "NO_TREND", "the age-vs-day trend is not distinguishable from a clock that reads nothing"
    return {
        "status": status, "n_donors": n,
        "mean_rho": mean, "ci95": [float(lo), float(hi)],
        "n_donors_negative": n_neg,
        "per_donor": {d: float(r) for d, r in finite.items()},
        "reason": f"{reason} (mean rho {mean:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}], {n_neg}/{n} negative)",
    }


def bars() -> list[dict]:
    """Pre-registered, resolvability-annotated (§5b)."""
    common_res = ("CONDITIONAL on trend consistency across 6 donors (T4-style): a moderate, "
                  "consistent trend passes ~99%; a heterogeneous one can read NO_TREND")
    return [
        {
            "id": "E1",
            "bar": "paired 95% CI on the mean per-donor Spearman(age, day) excludes 0 and is negative",
            "null": "a clock that reads nothing about the reprogramming axis (per-donor rho ~ 0)",
            "pass_rate_if_intent_holds": e1_power(-0.6, 0.25, 6),
            "resolvability": common_res,
        },
        {
            "id": "E1b",
            "bar": f"same test, reprogramming phase only (day <= {REPROG_PHASE_DAY_MAX:.0f})",
            "null": "no age drop during the OSKM phase (per-donor rho ~ 0)",
            "pass_rate_if_intent_holds": e1_power(-0.6, 0.25, 6),
            "resolvability": common_res + "; fewer points per donor than E1, so noisier per-donor rho",
        },
    ]


# --------------------------------------------------------------------------- #
# Real-data wiring (imports repo machinery only when actually run)            #
# --------------------------------------------------------------------------- #
def donor_trajectories(gill_dir: str, *, exclude_ipsc: bool = True) -> tuple[dict, dict]:
    """Per-donor (day, predicted_age) trajectory via the production normalisation + frozen clock.

    Works from the source's loaded matrix (like the M2 fix) so cell-type is available: the pipeline
    `obs` discards it, but the series-matrix `cell type` field carries it. iPSC samples are a
    cell-type change (clock out of domain) and are excluded from the primary.
    """
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
    src._load()                                     # populates src._rpm (genes x samples) + src._meta
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")

    cols = list(src._rpm.columns)
    mat = src._rpm[cols].to_numpy(dtype=np.float64).T             # samples x genes (linear RPM)
    ages = clock.predict_age(normalize_counts(mat), list(src._genes))

    traj: dict[str, dict] = {}
    n_ipsc = 0
    for col, age in zip(cols, ages, strict=True):
        m = src._meta.get(col)
        if m is None:
            continue
        is_ipsc = str(m.get("ctype", "")).strip().lower() == "ipsc"
        if is_ipsc:
            n_ipsc += 1
            if exclude_ipsc:
                continue
        d = traj.setdefault(m["donor"], {"days": [], "ages": [], "ctypes": []})
        d["days"].append(float(m["day"]))
        d["ages"].append(float(age))
        d["ctypes"].append(m.get("ctype", ""))
    meta = {"n_samples": len(cols), "n_ipsc": n_ipsc, "excluded_ipsc": exclude_ipsc,
            "n_donors": len(traj)}
    return traj, meta


def main() -> int:
    gill_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\Gill"
    print("STAGE 1.5 §8.3 E1 — within-donor reprogramming age trajectory (read-only)\n")
    print("  PRE-REGISTERED BAR (ground rule §5b):")
    for b in bars():
        print(f"    {b['id']}: {b['bar']}")
        print(f"        vs null: {b['null']}")
        print(f"        a correct system (rho=-0.6, sd=0.25, n=6) passes: {b['pass_rate_if_intent_holds']*100:.1f}%")
        print(f"        -> {b['resolvability']}")

    try:
        traj, tmeta = donor_trajectories(gill_dir, exclude_ipsc=True)
    except Exception as exc:  # noqa: BLE001
        print(f"\n   !! could not load Gill data ({exc!r}). Pass the dir: "
              'python experiments/diag_e1_trajectory.py "D:\\Gill"')
        return 1

    rhos = {d: donor_trend(v["days"], v["ages"]) for d, v in traj.items()}
    verdict = e1_verdict(rhos)
    # sensitivity: adults only, and with iPSC included
    adults = {d: r for d, r in rhos.items() if d not in NEONATAL}
    verdict_adults = e1_verdict(adults)
    traj_all, _ = donor_trajectories(gill_dir, exclude_ipsc=False)
    verdict_with_ipsc = e1_verdict({d: donor_trend(v["days"], v["ages"]) for d, v in traj_all.items()})

    # E1b: the reprogramming phase only (day <= REPROG_PHASE_DAY_MAX), where an MPTR dip should live.
    phase = {d: restrict_to_phase(v["days"], v["ages"], REPROG_PHASE_DAY_MAX) for d, v in traj.items()}
    rhos_b = {d: donor_trend(dd, aa) for d, (dd, aa) in phase.items()}
    verdict_e1b = e1_verdict(rhos_b)

    print(f"\n  {'donor':<7}{'chrono':>8}{'n(all)':>7}{'rho E1':>9}{'n(<=15)':>9}{'rho E1b':>9}")
    print("  " + "-" * 50)
    for d in sorted(traj, key=lambda x: DONOR_AGE.get(x, -1)):
        r, rb = rhos[d], rhos_b[d]
        print(f"  {d:<7}{str(DONOR_AGE.get(d)):>8}{len(traj[d]['days']):>7}{r:>9.3f}"
              f"{len(phase[d][0]):>9}{rb:>9.3f}")

    print(f"\n  E1  (full trajectory, iPSC excluded) : {verdict['status']}\n      {verdict['reason']}")
    print(f"  E1b (reprogramming phase, day <= {REPROG_PHASE_DAY_MAX:.0f})  : {verdict_e1b['status']}\n"
          f"      {verdict_e1b['reason']}")
    print(f"  E1 sensitivity — adults : {verdict_adults['status']} "
          f"(mean rho {verdict_adults.get('mean_rho', float('nan')):+.3f})")
    print(f"  E1 sensitivity — +iPSC  : {verdict_with_ipsc['status']} "
          f"(mean rho {verdict_with_ipsc.get('mean_rho', float('nan')):+.3f})")

    out = {"script": "diag_e1_trajectory", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "bars": bars(), "trajectory_meta": tmeta,
           "verdict_primary_ipsc_excluded": verdict,
           "verdict_e1b_reprogramming_phase": verdict_e1b,
           "e1b_day_max": REPROG_PHASE_DAY_MAX,
           "verdict_adults_only": verdict_adults,
           "verdict_with_ipsc": verdict_with_ipsc}
    Path("diag_e1_trajectory_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("\n  wrote diag_e1_trajectory_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
