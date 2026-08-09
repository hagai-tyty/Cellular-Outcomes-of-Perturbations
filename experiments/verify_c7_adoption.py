"""C-7 ADOPTION CHECK — what did the gate actually change, fold by fold?

    python experiments/verify_c7_adoption.py

READ-ONLY. Writes `results/verify_c7_adoption_results.json`. Compares the `_c7` folds against
`_armA`, which differ by exactly one config field (`bulk_integrity_gate`).

WHAT IT ASSERTS, AND WHY EACH ONE
----------------------------------
  A1  N2's ΔAge labels are MASKED in every `_c7` fold, with reason `no_control_baseline`.
      Rule 4 firing is the whole point of option (c); if it did not fire, the gate stripped a
      donor's zero-point and left the labels, which is the state B2' forbids.
  A2  the donor and the fold SURVIVE -- six folds, and N2's CELLS are still present. Option (c)
      differs from option (a) exactly here, and it is why the guard re-report stays at 6.
  A3  HFF's day-14 ΔAge MOVES, and moves toward the uncontaminated value. §5.14's reconstruction
      predicted -26.755 -> -8.196 when the degenerate control leaves `sigma_gill`; the built
      shards are the end-to-end test of that prediction.
  A4  HFF's day-14 ΔAge is now STABLE across folds. §4.7 measured a 16.67 yr spread; if the
      degenerate control was the carrier (§5.14, ATTRIBUTED), removing it should collapse that.

A3 and A4 are the ones that matter: they are the first end-to-end confirmation that the
reconstruction's number was right, measured on real built labels rather than predicted.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "verify_c7_adoption_results.json"

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
PREDICTED_CLEAN_DAY14 = -8.196      # §5.14's reconstruction, O1 fold minus N2's control
ARMA_DAY14 = {"N2": -7.352, "N3": -22.121, "O1": -24.023,
              "O2": -22.891, "Y1": -22.049, "Y2": -23.869}


def day_from_dose_time(dt: np.ndarray) -> np.ndarray:
    t = np.exp(np.asarray(dt, float)[:, 1])
    return np.round(t / 24.0 * 2.0) / 2.0


def fold_facts(root: Path) -> dict | None:
    import pandas as pd
    files = sorted(glob.glob(str(root / "shards" / "*.parquet")))
    if not files:
        return None
    n_cells = 0
    n2_total = n2_masked = 0
    reasons: dict[str, int] = {}
    hff_days, hff_age = [], []
    for f in files:
        d = pd.read_parquet(f, columns=["cell_line", "dose_time", "y_age",
                                        "age_mask", "age_mask_reason"])
        lines = d["cell_line"].to_numpy()
        n_cells += len(d)
        m2 = lines == "N2"
        n2_total += int(m2.sum())
        n2_masked += int((~d["age_mask"].to_numpy()[m2]).sum())
        for r in d["age_mask_reason"].to_numpy():
            if r is not None and r == r:
                reasons[str(r)] = reasons.get(str(r), 0) + 1
        mh = lines == "HFF"
        if mh.any():
            hff_days.append(day_from_dose_time(np.stack(d["dose_time"].to_numpy()[mh])))
            hff_age.append(d["y_age"].to_numpy()[mh].astype(float))
    day14 = None
    if hff_days:
        dd = np.concatenate(hff_days)
        aa = np.concatenate(hff_age)
        sel = np.isclose(dd, 14.0)
        if sel.any():
            day14 = float(aa[sel].mean())
    return {"n_cells": n_cells, "n2_cells": n2_total, "n2_masked": n2_masked,
            "reasons": reasons, "hff_day14": day14}


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    print("\n" + "=" * 78)
    print("C-7 ADOPTION CHECK — `_c7` vs `_armA`, one config field apart")
    print("=" * 78)

    c7, arma = {}, {}
    for d in DONORS:
        a = fold_facts(ROOT / f"cellfate_loocv_{d}_c7")
        b = fold_facts(ROOT / f"cellfate_loocv_{d}_armA")
        if a:
            c7[d] = a
        if b:
            arma[d] = b
    if not c7:
        print("   no `_c7` folds found -- run local_runners/build_c7_folds.py first")
        return 2

    print(render_table(
        ["fold", "cells", "N2 cells", "N2 masked", "HFF day-14 (c7)", "HFF day-14 (armA)"],
        [[d, str(c7[d]["n_cells"]), str(c7[d]["n2_cells"]), str(c7[d]["n2_masked"]),
          f"{c7[d]['hff_day14']:+.3f}" if c7[d]["hff_day14"] is not None else "n/a",
          f"{ARMA_DAY14[d]:+.3f}"] for d in DONORS if d in c7],
        aligns=["l", "r", "r", "r", "r", "r"]))

    checks: list[tuple[str, bool, str]] = []

    a1 = all(c7[d]["n2_masked"] == c7[d]["n2_cells"] and c7[d]["n2_cells"] > 0 for d in c7)
    reasons_seen = sorted({r for d in c7 for r in c7[d]["reasons"]})
    checks.append(("A1  N2's ΔAge fully masked in every fold", a1,
                   f"reasons seen: {reasons_seen or 'none'}"))
    checks.append(("A1b reason is `no_control_baseline`",
                   "no_control_baseline" in reasons_seen, str(reasons_seen)))

    a2 = len(c7) == 6 and all(c7[d]["n2_cells"] > 0 for d in c7)
    checks.append(("A2  donor and fold survive (option (c))", a2,
                   f"{len(c7)} folds, N2 cells present in all"))

    d14 = {d: c7[d]["hff_day14"] for d in c7 if c7[d]["hff_day14"] is not None}
    if d14:
        vals = np.array(list(d14.values()))
        spread = float(vals.max() - vals.min())
        near = float(np.abs(vals - PREDICTED_CLEAN_DAY14).mean())
        checks.append(("A3  HFF day-14 moved toward the reconstruction's value",
                       bool(near < 6.0),
                       f"mean |d14 - {PREDICTED_CLEAN_DAY14:+.3f}| = {near:.3f} yr"))
        checks.append(("A4  HFF day-14 spread collapsed (was 16.671 yr)",
                       bool(spread < 4.0), f"spread now {spread:.3f} yr"))

    print("\n" + render_table(["check", "result", "observed"],
                              [[n, "PASS" if ok else "FAIL", w] for n, ok, w in checks],
                              aligns=["l", "l", "l"]))

    payload = {"script": "verify_c7_adoption", "c7": c7, "armA_day14_recorded": ARMA_DAY14,
               "predicted_clean_day14": PREDICTED_CLEAN_DAY14,
               "checks": [{"name": n, "pass": bool(o), "observed": w} for n, o, w in checks],
               "all_pass": all(o for _, o, _ in checks)}
    OUT.write_text(json.dumps(payload, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(ROOT)}")
    return 0 if payload["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
