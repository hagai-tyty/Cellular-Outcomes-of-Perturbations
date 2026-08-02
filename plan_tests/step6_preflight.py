"""STAGE 1.5.3 step 6 PRE-FLIGHT — everything that must hold before spending ~10 h.

    python plan_tests/step6_preflight.py

Builds BOTH arms at reduced scale into scratch roots and checks, end to end, the things that
silently wasted the first run. Writes `results/step6_preflight_results.json`. Exits non-zero on any
failure so it can gate the real run.

WHAT THIS EXISTS TO CATCH
-------------------------
The first step-6 run cost ~10 h and returned a result that could not be interpreted, for reasons
that were all visible before it started:

  P1  The arm switch never reached the data (wrong driver reused shards). Both arms trained
      identically and the null was an artefact. -> assert the two arms' label counts DIFFER.

  P2  `age_window_k` was not plumbed, so C-5 Option 2 would have run at k = 1 = OFF whatever the
      plan said. -> assert k reaches TrainConfig in both arms.

  P3  🔴 THE ONE THAT ACTUALLY INVALIDATED IT. `y_age` depended on the training-label policy,
      because the deconfounder fit and the control re-centring both used `age_mask`. Masking HFF
      redefined the TARGET VARIABLE, so the arms differed in two ways at once (C-I in
      results/STEP6_REPORT.md). -> assert `y_age` is BIT-IDENTICAL between the arms, per cell.

  P4  The arms shared a fold root, so arm B overwrote arm A and its deconfounder coefficient had
      to be reported from a proxy. -> assert the suffixed roots are distinct and both survive.

P3 is the whole point. A unit test on synthetic `ChunkAux` already covers the logic
(`tests/test_ci_deconfounder_arm_invariance.py`); this checks the REAL pipeline, because the bug
was never in the arithmetic -- it was in which mask three different call sites happened to read.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

SMOKE = {"epochs": 2, "ensemble": 1, "max_cells": 800, "holdout": "O1"}


def _rml():
    spec = importlib.util.spec_from_file_location(
        "rml_preflight", ROOT / "local_runners" / "run_multi_local.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def build_arm(arm: str, root: str) -> dict:
    """Build one arm at smoke scale into `root`. Returns its census."""
    rml = _rml()
    rml.GSE_DIR, rml.GILL_DIR = r"D:\GSE242423", r"D:\Gill"
    rml.HARMONIZE = True
    rml.HOLDOUT_DONOR = SMOKE["holdout"]
    rml.ROOT = root
    rml.AGE_MASKED = frozenset({"hff_sc"}) if arm == "B" else frozenset()
    rml.AGE_WINDOW_K = 4
    rml.EPOCHS, rml.ENSEMBLE, rml.MAX_CELLS = (
        SMOKE["epochs"], SMOKE["ensemble"], SMOKE["max_cells"])
    rml.main()
    return json.loads(Path(root, "step6_arm_census.json").read_text(encoding="utf-8"))


def read_y_age(root: str) -> dict[str, tuple[float, bool]]:
    """cell_id -> (y_age, age_mask), read straight off the shards. Pure I/O."""
    from cellfate.common import io
    out: dict[str, tuple[float, bool]] = {}
    for sh in sorted(Path(root, "shards").glob("*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        for i, cid in enumerate(a["cell_id"]):
            out[str(cid)] = (float(a["y_age"][i]), bool(a["age_mask"][i]))
    return out


def compare_y_age(A: dict, B: dict) -> dict:
    """P3. `y_age` must match cell-for-cell; only `age_mask` may differ. Pure."""
    common = sorted(set(A) & set(B))
    both_nan = mism = 0
    worst = 0.0
    worst_cell = None
    for c in common:
        ya, yb = A[c][0], B[c][0]
        if np.isnan(ya) and np.isnan(yb):
            both_nan += 1
            continue
        if np.isnan(ya) != np.isnan(yb):
            mism += 1
            continue
        d = abs(ya - yb)
        if d > worst:
            worst, worst_cell = d, c
    mask_diff = sum(1 for c in common if A[c][1] != B[c][1])
    return {"n_common": len(common), "nan_pattern_mismatches": mism, "both_nan": both_nan,
            "max_abs_y_age_delta": worst, "worst_cell": worst_cell,
            "cells_whose_age_mask_differs": mask_diff}


def main() -> int:
    print("STAGE 1.5.3 STEP 6 PRE-FLIGHT — gating a ~10 h run\n")
    rootA, rootB = "cellfate_pf_armA", "cellfate_pf_armB"
    for r in (rootA, rootB):
        shutil.rmtree(r, ignore_errors=True)

    print(f"  building arm A -> {rootA} (smoke scale {SMOKE}) ...")
    cenA = build_arm("A", rootA)
    print(f"  building arm B -> {rootB} ...")
    cenB = build_arm("B", rootB)

    yA, yB = read_y_age(rootA), read_y_age(rootB)
    cmp_ = compare_y_age(yA, yB)

    checks: list[tuple[str, bool, str]] = []

    # P1 -- the arms must differ in WHICH labels are trainable
    p1 = cenA["n_age_valid_train"] != cenB["n_age_valid_train"] and cenB["n_age_valid_train"] > 0
    checks.append(("P1 arms differ in trainable labels", p1,
                   f"A {cenA['n_age_valid_train']:,} vs B {cenB['n_age_valid_train']:,} "
                   f"of {cenA['n_train_cells']:,}"))

    # P2 -- k reached both arms
    p2 = cenA["age_window_k"] == cenB["age_window_k"] == 4
    checks.append(("P2 age_window_k = 4 in both arms", p2,
                   f"A k={cenA['age_window_k']}, B k={cenB['age_window_k']}"))

    # P3 -- THE ONE. y_age must not move.
    p3 = (cmp_["nan_pattern_mismatches"] == 0 and cmp_["max_abs_y_age_delta"] == 0.0
          and cmp_["n_common"] > 0)
    checks.append(("P3 y_age BIT-IDENTICAL across arms", p3,
                   f"max|delta| {cmp_['max_abs_y_age_delta']:.3e} over {cmp_['n_common']:,} cells, "
                   f"{cmp_['nan_pattern_mismatches']} NaN-pattern mismatches"))

    # P3b -- and it must not be vacuous: the arms DO differ, in age_mask only
    p3b = cmp_["cells_whose_age_mask_differs"] > 0
    checks.append(("P3b the difference is in age_mask ONLY (not vacuous)", p3b,
                   f"{cmp_['cells_whose_age_mask_differs']:,} cells differ in age_mask"))

    # P4 -- separate roots, both alive
    p4 = Path(rootA).is_dir() and Path(rootB).is_dir() and rootA != rootB
    checks.append(("P4 arms kept SEPARATE builds", p4, f"{rootA} and {rootB} both present"))

    print(f"\n  {'check':<48}{'':<6}detail")
    print("  " + "-" * 100)
    for name, ok, detail in checks:
        print(f"  {name:<48}{'PASS' if ok else 'FAIL':<6}{detail}")

    all_ok = all(ok for _, ok, _ in checks)
    print(f"\n  ==> PRE-FLIGHT {'PASS -- clear to run step 6' if all_ok else 'FAIL -- DO NOT RUN'}")
    if not all_ok:
        print("      Fix the failing check before spending the compute.")

    out = {"script": "step6_preflight", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "smoke": SMOKE, "census_arm_a": cenA, "census_arm_b": cenB,
           "y_age_comparison": cmp_,
           "checks": [{"name": n, "pass": bool(o), "detail": d} for n, o, d in checks],
           "all_pass": all_ok}
    (_RESULTS / "step6_preflight_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("  wrote step6_preflight_results.json")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
