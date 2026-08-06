"""ARM D pre-flight — INTRINSIC per-stratum validation (single build, no cross-build compare).

    python plan_tests/armd_intrinsic_preflight.py

Builds ONE arm-D dataset at smoke scale and asserts, inside that build, that the stratified shuffle
is a within-stratum permutation: for every `(cell_line, time_h)` stratum, the pre-shuffle ΔAge
multiset equals the post-shuffle multiset, so each timepoint's MEAN ΔAge is preserved exactly and
the between-timepoint trajectory (rho(day, ΔAge) = -0.905) survives. Writes
`results/armd_intrinsic_preflight_results.json`; exits non-zero on any failure.

WHY INTRINSIC, NOT A vs D
-------------------------
The first arm-D pre-flight compared an arm-A build against an arm-D build and failed. The failure was
NOT arm D: comparing HFF *training-label values* across two separate smoke builds is unreliable, and
it is also the wrong thing to check -- arm A and arm D differ only in HFF training labels, while the
metric (`rank_model_dage`) is scored on the never-shuffled held-out donor. The property that makes
arm D valid is intrinsic to a single build: the shuffle must permute WITHIN each timepoint. This
captures the ΔAge pool immediately before and after the shuffle and checks exactly that -- no second
build, no confound. (The pure-permutation logic itself is unit-tested in
`tests/test_arm_d_stratified_shuffle.py`; this confirms the real stratum keys are built correctly
from `raw.obs` end to end.)
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_RESULTS = ROOT / "results"
_RESULTS.mkdir(exist_ok=True)

from cellfate.data import build_dataset as bd  # noqa: E402

_cap: dict = {}
_orig = bd._shuffle_age_labels


def _spy(cfg, aux_by_sid, ys):
    def pool(source):
        g = defaultdict(list)
        for sid in sorted(aux_by_sid):
            aux = aux_by_sid[sid]
            arr = source[sid]
            for i in np.flatnonzero(aux.shuffle_mask):
                v = arr[i]
                if not np.isnan(v):
                    g[str(aux.stratum[i])].append(float(v))
        return g

    pre = pool(ys)
    out = _orig(cfg, aux_by_sid, ys)
    _cap["pre"], _cap["post"] = pre, pool(out)
    return out


def main() -> int:
    bd._shuffle_age_labels = _spy
    root = "cellfate_armd_pf"
    shutil.rmtree(ROOT / root, ignore_errors=True)

    spec = importlib.util.spec_from_file_location(
        "rml_armd_pf", ROOT / "local_runners" / "run_multi_local.py")
    rml = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rml)
    rml.GSE_DIR, rml.GILL_DIR = r"D:\GSE242423", r"D:\Gill"
    rml.HARMONIZE = True
    rml.HOLDOUT_DONOR = "O1"
    rml.ROOT = root
    rml.AGE_MASKED = frozenset()
    rml.AGE_SHUFFLE = frozenset({"hff_sc"})
    rml.AGE_SHUFFLE_STRATA = True
    rml.AGE_SHUFFLE_SEED = 0
    rml.AGE_WINDOW_K = 4
    rml.EPOCHS, rml.ENSEMBLE, rml.MAX_CELLS = 1, 1, 800
    rml.main()

    pre, post = _cap["pre"], _cap["post"]
    rows = []
    all_eq = True
    worst = 0.0
    for k in sorted(pre):
        a, b = sorted(pre[k]), sorted(post[k])
        eq = len(a) == len(b) and bool(np.allclose(a, b))
        dmean = float(abs(np.mean(pre[k]) - np.mean(post[k])))
        all_eq &= eq
        worst = max(worst, dmean)
        rows.append({"stratum": k, "n": len(a), "multiset_eq": eq,
                     "mean": float(np.mean(pre[k])), "dmean": dmean})

    total = sum(r["n"] for r in rows)
    n_strata = len(rows)
    checks = {
        "per_stratum_multiset_preserved": all_eq,
        "trajectory_means_preserved": worst < 1e-6,
        "stratified_into_multiple_timepoints": n_strata > 1,
        "labels_present": total > 100,
    }
    ok = all(checks.values())

    print("\nARM D INTRINSIC PRE-FLIGHT (single build)\n")
    print(f"  {n_strata} strata, {total} target labels")
    for r in rows:
        print(f"    {r['stratum']:<12} n={r['n']:>4}  mean {r['mean']:+7.3f}  "
              f"multiset_eq={r['multiset_eq']}  |dmean|={r['dmean']:.1e}")
    print(f"\n  worst per-stratum |mean delta|: {worst:.2e}")
    for name, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {name}")
    print(f"\n  ==> {'PASS -- clear to run arm D' if ok else 'FAIL -- DO NOT RUN'}")

    (_RESULTS / "armd_intrinsic_preflight_results.json").write_text(
        json.dumps({"n_strata": n_strata, "total_labels": total,
                    "worst_stratum_mean_delta": worst, "strata": rows,
                    "checks": checks, "all_pass": ok}, indent=2), encoding="utf-8")
    shutil.rmtree(ROOT / root, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
