"""
MATCHED-DATA REPRODUCTION — does July's HFF ΔAge trajectory reproduce on arm A?

    python experiments/repro_hff_signature_armA.py

READ-ONLY. Writes `results/repro_hff_signature_armA_results.json`. `src/` untouched, no
build touched, and the July reference `results/diag_gc_hff_signature_results.json` is
NEVER written by this script (the July script writes that path; this one does not).

WHY A SECOND ATTEMPT WAS NEEDED
-------------------------------
The first attempt compared arm A's HFF trajectory against the July reference and started
to flag a ~2.4x mismatch (arm A ≈ −9.9 against the recorded −24.02). The user stopped it:
*"arm A is from 2 data sets while the tests run on 1 data set — of course it would be
different results."* The comparison was indeed invalid, so it was withdrawn.

**But the dataset count was not the operative cause, and this script establishes that.**
There are two HFF references and they are NOT the same kind of measurement:

  results/diag_pipeline_decompose_results.json   RAW GSE242423, ONE dataset, harmonization
                                                 off, clock applied directly. day-14 ≈ −8.5
                                                 to −10.6. Comparing this to a harmonized
                                                 two-dataset build IS invalid — the user's
                                                 point, and it stands for this reference.
  results/diag_gc_hff_signature_results.json     BUILT SHARDS of a two-dataset harmonized
                                                 LOOCV build. day-14 = −24.02. This is the
                                                 number the first attempt compared against,
                                                 and it was ALREADY a two-dataset build.

So the −24.02 side was never single-dataset. The real error was mine and simpler: I read
arm A's **N2** fold and compared it to a reference produced from the **O1** fold
(`diag_gc_hff_signature.py` defaults to `runs/cellfate_loocv_O1`). That is the matched
comparison this script runs, and it also measures how much the fold choice matters —
which turns out to be the finding.

WHAT IS COMPARED
----------------
`load_hff` and `trajectory_stats` are imported from the July script UNMODIFIED, so the
arithmetic is identical on both sides; only the run directory changes.

  PART 1  MATCHED FOLD — arm A's O1 fold against the July reference, field by field.
  PART 2  FOLD STABILITY — the same trajectory on all six arm-A folds. HFF is never the
          held-out line in any of them (only a Gill donor is), and HFF is 42481 of 42605
          cells, so a stable pipeline should give near-identical HFF labels in all six.

PRE-REGISTERED, before any arm-A number was read for PART 1:
  R1  arm A's O1 fold reproduces the July reference EXACTLY (float equality) on
      rho_timepoint, slope_yr_per_day, n_cells, days, mean_dage, sem_dage, n_per_day.
      Basis: the RES reproduction already showed this build is bit-stable against July,
      and `y_age` is read straight off the shards with no model involved.
  R2  the six folds agree with each other to within a tolerance of |Δ day-14| <= 2.0 yr
      and |Δ slope| <= 0.30 yr/day. Basis: holding out one Gill donor removes ~21 of
      42605 cells (0.05%), which should not move HFF's labels materially.

  R1 fails -> the build is not the same system as July; escalate.
  R2 fails -> HFF's labels depend on which donor is held out. Since HFF supplies 99.7%
      of the age-labelled corpus, that is a defect in the label construction, not a
      curiosity, and it becomes a candidate source of the between-fold variance that made
      step 6 inconclusive (observed SD 4.808, MDE 5.045).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
JULY_FOLD = "O1"                     # diag_gc_hff_signature.py defaults to runs/cellfate_loocv_O1
JULY_REF = REPO / "results" / "diag_gc_hff_signature_results.json"
OUT = REPO / "results" / "repro_hff_signature_armA_results.json"

EXACT_FIELDS = ["n_cells", "n_timepoints", "rho_timepoint", "slope_yr_per_day",
                "rho_percell", "slope_percell", "n_descending_steps", "n_steps"]
ARRAY_FIELDS = ["days", "mean_dage", "sem_dage", "n_per_day"]

TOL_DAY14 = 2.0      # R2, years
TOL_SLOPE = 0.30     # R2, years per day


def _july_module():
    """Import the July script unmodified, for its load_hff / trajectory_stats."""
    path = REPO / "experiments" / "diag_gc_hff_signature.py"
    spec = importlib.util.spec_from_file_location("_gc_hff", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_gc_hff"] = mod
    spec.loader.exec_module(mod)
    return mod


def day14(traj: dict) -> float:
    return float(traj["mean_dage"][traj["days"].index(14.0)])


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    gc = _july_module()
    july = json.loads(JULY_REF.read_text(encoding="utf-8"))["hff_trajectory"]

    print("\n" + "=" * 78)
    print("MATCHED-DATA REPRODUCTION — HFF ΔAge trajectory, July vs arm A")
    print("=" * 78)
    print("Both sides read BUILT SHARDS of a two-dataset harmonized LOOCV build.")
    print(f"July's reference came from the {JULY_FOLD} fold (the July script's default run dir).")

    per_fold: dict[str, dict] = {}
    for d in DONORS:
        run = REPO / f"cellfate_loocv_{d}_armA"
        if not run.exists():
            print(f"   missing {run}")
            return 2
        days, ages, volume = gc.load_hff(run)
        st = gc.trajectory_stats(days, ages)
        hz = json.loads((run / "harmonization.json").read_text(encoding="utf-8"))
        st["_n_harmonized_genes"] = len(hz.get("genes", []))
        st["_n_hff_cells_total"] = int(volume["by_line"].get("HFF", 0))
        per_fold[d] = st

    # ------------------------------------------------------------ PART 1
    print(f"\n\n--- PART 1: MATCHED FOLD ({JULY_FOLD}) vs the July reference ---\n")
    arm = per_fold[JULY_FOLD]
    rows, all_exact = [], True
    for f in EXACT_FIELDS:
        a, j = arm.get(f), july.get(f)
        ok = a == j
        all_exact &= ok
        rows.append([f, f"{a}"[:24], f"{j}"[:24], "EXACT" if ok else "DIFFERS"])
    for f in ARRAY_FIELDS:
        a, j = list(arm.get(f, [])), list(july.get(f, []))
        ok = a == j
        all_exact &= ok
        rows.append([f, f"len={len(a)}", f"len={len(j)}",
                     "EXACT (elementwise)" if ok else "DIFFERS"])
    print(render_table(["field", "arm A", "July", "match"], rows,
                       aligns=["l", "r", "r", "l"]))
    print(f"\n   R1  arm A {JULY_FOLD} reproduces July EXACTLY: "
          f"{'PASS' if all_exact else 'FAIL'}")

    # ------------------------------------------------------------ PART 2
    print("\n\n--- PART 2: FOLD STABILITY — the same trajectory on all six arm-A folds ---")
    print("HFF is never the held-out line; only a Gill donor is. Holding one out removes")
    print("~21 of 42605 cells (0.05%), so these six should be near-identical.\n")
    rows = []
    for d in DONORS:
        t = per_fold[d]
        rows.append([d, f"{day14(t):+.3f}", f"{t['rho_timepoint']:+.4f}",
                     f"{t['slope_yr_per_day']:+.4f}",
                     f"{t['n_descending_steps']}/{t['n_steps']}",
                     str(t["_n_harmonized_genes"]), str(t["_n_hff_cells_total"])])
    rows.append(["JULY", f"{day14(july):+.3f}", f"{july['rho_timepoint']:+.4f}",
                 f"{july['slope_yr_per_day']:+.4f}",
                 f"{july['n_descending_steps']}/{july['n_steps']}", "-", "-"])
    print(render_table(
        ["fold", "day-14 ΔAge", "rho_timepoint", "slope yr/day", "desc", "genes", "HFF cells"],
        rows, aligns=["l", "r", "r", "r", "r", "r", "r"]))

    d14 = {d: day14(per_fold[d]) for d in DONORS}
    sl = {d: float(per_fold[d]["slope_yr_per_day"]) for d in DONORS}
    spread14 = max(d14.values()) - min(d14.values())
    spread_sl = max(sl.values()) - min(sl.values())
    worst14 = max(DONORS, key=lambda d: abs(d14[d] - np.median(list(d14.values()))))
    r2 = bool(spread14 <= TOL_DAY14 and spread_sl <= TOL_SLOPE)

    print(f"\n   day-14 spread  {spread14:.3f} yr    (tolerance {TOL_DAY14:.2f})")
    print(f"   slope  spread  {spread_sl:.3f} yr/day (tolerance {TOL_SLOPE:.2f})")
    print(f"   furthest fold from the median: {worst14} "
          f"({d14[worst14]:+.3f} vs median {np.median(list(d14.values())):+.3f})")
    print(f"\n   R2  the six folds agree within tolerance: {'PASS' if r2 else 'FAIL'}")

    print("\n   WHAT THIS MEANS:")
    if all_exact and r2:
        print("     -> HFF labels reproduce July exactly AND are stable across folds. Closed.")
    elif all_exact and not r2:
        print("     -> The MATCHED comparison reproduces July exactly, so the pipeline is the")
        print("        same system it was in July. But HFF's ΔAge labels are NOT stable across")
        print("        LOOCV folds. HFF supplies 42481 of 42605 age-labelled cells (99.7%), so")
        print("        the training signal itself changes with which Gill donor is held out.")
        print("        Candidate mechanism, NOT established here: harmonization is refit per")
        print("        fold against the gill_bulk reference, and that reference is small enough")
        print("        that dropping one donor perturbs it. Also a candidate source of the")
        print("        between-fold variance that made step 6 inconclusive (SD 4.808).")
    else:
        print("     -> R1 FAILED: the matched fold does not reproduce July. Escalate before")
        print("        reading anything else here.")

    payload = {
        "script": "repro_hff_signature_armA",
        "july_fold": JULY_FOLD,
        "tolerances": {"day14_yr": TOL_DAY14, "slope_yr_per_day": TOL_SLOPE},
        "R1_matched_fold_exact": bool(all_exact),
        "R2_folds_agree": r2,
        "matched_fold_comparison": {
            f: {"armA": arm.get(f), "july": july.get(f), "exact": arm.get(f) == july.get(f)}
            for f in EXACT_FIELDS
        },
        "per_fold": per_fold,
        "day14_by_fold": d14,
        "slope_by_fold": sl,
        "day14_spread": spread14,
        "slope_spread": spread_sl,
        "furthest_fold": worst14,
        "july_reference": july,
    }
    OUT.write_text(json.dumps(payload, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print(f"   (July reference {JULY_REF.relative_to(REPO)} was NOT written by this script)")
    return 0 if all_exact else 1


if __name__ == "__main__":
    raise SystemExit(main())
