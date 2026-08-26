"""Stage 24 — Gen-1 Role-B predictor engineering, under the frozen Stage-23.5 contract.

Opened by `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`, recorded in
`results/stage23_5_protocol.json`.

This is a **bounded predictor-engineering stage**, not an architecture search. Its whole job is to
reproduce the frozen W5 result, package it, and hand Stage 25 a set of out-of-fold predictions it
can rank. Per the handoff it may NOT:

    inspect the Stage-25 ranking metric
    replace W5 because another architecture scores better on the same folds
    add a dataset
    change any Stage-22 or Stage-23 frozen quantity

The reproduction never overwrites a frozen artifact. `run_stage23_learnability_gate` writes its
WM989 outputs to module-level path constants; 24B rebinds those constants to `results/stage24/repro/`
for the duration of the call, so the frozen Stage-23 files stay exactly as committed and become the
comparison target rather than the destination. The models themselves are the frozen implementation,
called rather than re-typed (plan §5.1).
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage23_learnability_gate as S23  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
# The repo convention (tests/test_results_paths.py): __file__-relative, never CWD-relative, in
# exactly this literal form so the invariant can check where every writer writes.
_RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS = _RESULTS
OUT = _RESULTS / "stage24"
REPRO = OUT / "repro"
for _d in (OUT, REPRO):
    _d.mkdir(parents=True, exist_ok=True)

PLAN = ROOT / "plans" / "(newer)practical plans" / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"
PROTOCOL_JSON = RESULTS / "stage23_5_protocol.json"
HANDOFF_JSON = RESULTS / "stage23_5_handoff_to_stage24.json"

A_JSON = OUT / "stage24a_engineering_plan.json"
B_JSON = OUT / "stage24b_reproduction.json"

# ---- the frozen artifacts 24B reproduces AGAINST, never INTO -------------------------------- #
FROZEN = {
    "C1_W0toW4": RESULTS / "stage23_wm989_detection_oof.csv",
    "C2_W0toW4": RESULTS / "stage23_wm989_abundance_oof.csv",
    "C1_W5": RESULTS / "stage23_wm989_interaction_oof.csv",
    "C2_W5": RESULTS / "stage23_wm989_interaction_abundance_oof.csv",
    "results_W0toW4": RESULTS / "stage23_wm989_results.json",
    "results_W5": RESULTS / "stage23_wm989_interaction_results.json",
}
KEY_COLS = ["clone_id", "treatment", "outer_fold", "y"]
ROW_TOLERANCE = 1e-12          # plan §7.1 R2

# Row counts are per ENDPOINT, not global. C1 (detection) scores every clone x condition row;
# C2 (abundance) is defined only where a clone was detected, so it carries the 2,256 nonzero rows.
# An earlier version of this gate hardcoded 8,406 everywhere and reported R1=False on a
# byte-identical C2 file -- harmless there because byte-identity carried the verdict, but it would
# have forced a spurious INPUT_INTEGRITY_STOP had the reproduction been merely tolerance-clean.
EXPECTED_ROWS = {"C1": 8406, "C2": 2256}
EXPECTED_CLONES = 1401
TREATMENTS = ("Dabrafenib", "Trametinib", "CoCl2", "Acid", "Cisplatin", "Doxorubicin")


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_json(p: Path, obj: dict) -> None:
    p.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


# =============================================================================================== #
# 24A — consume the handoff and freeze the engineering plan
# =============================================================================================== #
def run_24a() -> dict:
    """Assert the Stage-23.5 freeze is intact before any engineering begins.

    A handoff that does not match the plan it claims to come from is not a handoff.
    """
    t0 = time.perf_counter()
    proto = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    handoff = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    plan_digest = S23.canonical_text_sha256(PLAN)

    checks = {
        "plan_digest_matches_protocol":
            plan_digest == proto["plan_canonical_lf_sha256"],
        "plan_digest_matches_handoff":
            plan_digest == handoff["plan_canonical_lf_sha256"],
        "plan_status_frozen": proto["plan_status"] == "FROZEN",
        "stage_24_open": handoff["stage_24_opening_status"] == "STAGE_24_OPEN_ROLE_B_PRIMARY_GEN1",
        "ranking_metric_not_inspected": handoff["ranking_metric_inspected"] is False,
        "ranking_protocol_frozen": handoff["ranking_protocol_hash_frozen"] is True,
        "no_new_datasets": handoff["new_datasets_authorized"] is False,
        "stage_27_not_a_gen1_gate": handoff["stage_27_gen1_gate"] is False,
        "audit_fully_passed":
            proto["audit"]["checklist_passed"] == proto["audit"]["checklist_items"],
        "compute_budget_accepted":
            proto["permutation_contract"]["compute_budget"]["accepted_by_decision_owner"] is True,
    }
    # every source artifact the protocol pinned must still hash the same
    drift = {p: {"pinned": h, "now": sha256_file(ROOT / p)}
             for p, h in proto["source_artifacts"].items()
             if sha256_file(ROOT / p) != h}
    checks["no_source_artifact_drift"] = not drift

    # no Stage-25 statistic may exist
    ranking_leak = sorted(p.name for p in RESULTS.rglob("*")
                          if p.is_file() and "rank" in p.name.lower())
    checks["no_ranking_artifact_exists"] = not ranking_leak

    out = {
        "stage": "24A",
        "plan": {"file": PLAN.name, "canonical_lf_sha256": plan_digest},
        "checks": checks,
        "source_artifact_drift": drift,
        "ranking_artifacts_found": ranking_leak,
        "all_checks_pass": all(checks.values()),
        "engineering_plan": {
            "model": "W5 = X + B + U + X*U, the frozen Stage-23 implementation",
            "substages": handoff["stage_24_substages"],
            "may_not": handoff["stage_24_may_not"],
            "reproduction_gate": "plan §7.1 R1-R4, implemented in run_24b()",
            "output_root": "results/stage24/",
            "frozen_artifacts_are_read_only": True,
            "how_reproduction_avoids_overwriting": (
                "run_stage23_learnability_gate writes WM989 outputs to module-level path "
                "constants. 24B rebinds those constants to results/stage24/repro/ for the "
                "duration of the call, so the frozen files are the comparison target and are "
                "never written to."),
        },
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 4),
        "source_provenance": S23.source_provenance() if hasattr(S23, "source_provenance") else {},
    }
    write_json(A_JSON, out)
    if not out["all_checks_pass"]:
        raise RuntimeError(f"24A integrity failure: {[k for k, v in checks.items() if not v]}")
    return out


# =============================================================================================== #
# 24B — reproduce W1/W4/W5 and apply the §7.1 row-level gate
# =============================================================================================== #
@contextlib.contextmanager
def _redirected_wm989_outputs(dest: Path):
    """Point the frozen builder's WM989 write targets at `dest`, then restore them.

    The builder is called, not re-implemented (plan §5.1). Only where it WRITES changes; every
    input it reads -- the Stage-22 tables, the pseudobulk cache, the protocol -- is untouched.
    """
    names = {
        "WM989_C1_OOF": dest / "stage23_wm989_detection_oof.csv",
        "WM989_C2_OOF": dest / "stage23_wm989_abundance_oof.csv",
        "WM989_RESULTS": dest / "stage23_wm989_results.json",
        "WM989_W5_OOF": dest / "stage23_wm989_interaction_oof.csv",
        "WM989_INTERACTION_RESULTS": dest / "stage23_wm989_interaction_results.json",
    }
    saved = {k: getattr(S23, k) for k in names}
    try:
        for k, v in names.items():
            setattr(S23, k, v)
        yield names
    finally:
        for k, v in saved.items():
            setattr(S23, k, v)


def _gate_r1(frozen: pd.DataFrame, repro: pd.DataFrame, endpoint: str) -> dict:
    """§7.1 R1 — shape and key, against the row count documented for THIS endpoint.

    The three sub-checks are computed independently rather than short-circuited, so a row-count
    mismatch reports the key result too instead of masking it.
    """
    expected = EXPECTED_ROWS[endpoint]
    counts_agree = len(frozen) == len(repro)
    matches_documented = len(frozen) == expected
    have_keys = set(KEY_COLS) <= set(frozen.columns) and set(KEY_COLS) <= set(repro.columns)
    key_ok = bool(have_keys and counts_agree
                  and frozen[KEY_COLS].reset_index(drop=True).equals(
                      repro[KEY_COLS].reset_index(drop=True)))
    return {"endpoint": endpoint,
            "rows_frozen": int(len(frozen)), "rows_repro": int(len(repro)),
            "expected_rows_for_endpoint": expected,
            "repro_matches_frozen_row_count": bool(counts_agree),
            "frozen_matches_documented_count": bool(matches_documented),
            "key_columns_present": bool(have_keys),
            "ordered_key_identical": key_ok,
            "pass": bool(counts_agree and matches_documented and key_ok)}


def _gate_r2(frozen: pd.DataFrame, repro: pd.DataFrame) -> dict:
    """§7.1 R2 — every prediction cell, checked cell by cell."""
    cols = [c for c in frozen.columns if c.startswith("pred_")]
    per_col, worst = {}, 0.0
    for c in cols:
        if c not in repro.columns:
            per_col[c] = {"present_in_repro": False, "pass": False}
            continue
        d = np.abs(frozen[c].to_numpy(dtype=float) - repro[c].to_numpy(dtype=float))
        n_bad = int((d > ROW_TOLERANCE).sum())
        worst = max(worst, float(d.max()))
        per_col[c] = {"present_in_repro": True, "max_abs_diff": float(d.max()),
                      "cells_over_tolerance": n_bad, "cells": int(len(d)),
                      "pass": n_bad == 0}
    return {"tolerance": ROW_TOLERANCE, "columns_checked": cols,
            "worst_abs_diff_any_column": worst, "per_column": per_col,
            "pass": bool(cols) and all(v["pass"] for v in per_col.values())}


def _within_clone_sign_matrix(df: pd.DataFrame, col: str) -> dict:
    """For each clone, the sign of every ordered condition pair difference."""
    out = {}
    for cid, grp in df.groupby("clone_id", sort=True):
        g = grp.sort_values("treatment")
        s = g[col].to_numpy(dtype=float)
        out[cid] = tuple(int(np.sign(s[i] - s[j]))
                         for i in range(len(s)) for j in range(len(s)) if i < j)
    return out


def _gate_r3(frozen: pd.DataFrame, repro: pd.DataFrame, cols: list[str]) -> dict:
    """§7.1 R3 — within-clone ordering unchanged, including tie structure.

    This is the load-bearing check: the Stage-25 ranking test is a function of within-clone
    orderings and nothing else, so a reproduction that preserved every pooled metric while
    flipping one near-tie would not have reproduced the input the test consumes.
    """
    per_col = {}
    for c in cols:
        if c not in repro.columns:
            per_col[c] = {"present_in_repro": False, "pass": False}
            continue
        a = _within_clone_sign_matrix(frozen, c)
        b = _within_clone_sign_matrix(repro, c)
        changed = sorted(k for k in a if a[k] != b.get(k))
        per_col[c] = {"present_in_repro": True, "clones_checked": len(a),
                      "clones_with_changed_ordering": len(changed),
                      "examples": changed[:5], "pass": not changed}
    return {"models_checked": cols, "per_model": per_col,
            "pass": bool(cols) and all(v["pass"] for v in per_col.values())}


def run_24b(wm989_root: Path) -> dict:
    """Reproduce the frozen W0-W5 chain into results/stage24/repro/ and gate it."""
    if not A_JSON.exists():
        raise RuntimeError("24A must run before 24B")
    a = json.loads(A_JSON.read_text(encoding="utf-8"))
    if not a["all_checks_pass"]:
        raise RuntimeError("24A did not pass; 24B may not run")

    t0 = time.perf_counter()
    print("  24B: reproducing 23C (W0-W4) then 23D (W5) into results/stage24/repro/ ...",
          flush=True)
    with _redirected_wm989_outputs(REPRO) as dest:
        t_c = time.perf_counter()
        S23.run_23c(wm989_root)
        c_min = (time.perf_counter() - t_c) / 60
        print(f"  24B: 23C reproduced in {c_min:.1f} min", flush=True)
        t_d = time.perf_counter()
        S23.run_23d(wm989_root)
        d_min = (time.perf_counter() - t_d) / 60
        print(f"  24B: 23D reproduced in {d_min:.1f} min", flush=True)
        produced = {k: v for k, v in dest.items()}

    gates: dict = {}
    for label, endpoint, frozen_path, repro_path, r3_models in (
        ("C1_W0toW4", "C1", FROZEN["C1_W0toW4"], produced["WM989_C1_OOF"], ["pred_W4"]),
        ("C2_W0toW4", "C2", FROZEN["C2_W0toW4"], produced["WM989_C2_OOF"], []),
        ("C1_W5", "C1", FROZEN["C1_W5"], produced["WM989_W5_OOF"], ["pred_W4", "pred_W5"]),
    ):
        fz = pd.read_csv(frozen_path)
        rp = pd.read_csv(repro_path)
        byte_identical = sha256_file(frozen_path) == sha256_file(repro_path)
        r1 = _gate_r1(fz, rp, endpoint)
        r2 = _gate_r2(fz, rp)
        r3 = _gate_r3(fz, rp, r3_models) if r3_models else {"models_checked": [], "pass": True,
                                                            "note": "no ranking-relevant model"}
        all_r = r1["pass"] and r2["pass"] and r3["pass"]
        # A byte-identical file that fails any sub-gate means the SUB-GATE is broken, not the
        # reproduction. Surface that rather than letting byte-identity paper over it.
        gate_self_consistent = (not byte_identical) or all_r
        gates[label] = {
            "frozen": frozen_path.name, "reproduced": str(repro_path.relative_to(ROOT)),
            "byte_identical": byte_identical,
            "R1_shape_and_key": r1, "R2_every_score": r2, "R3_within_clone_ordering": r3,
            "gate_self_consistent": gate_self_consistent,
            "verdict": ("BYTE_IDENTICAL" if byte_identical
                        else ("TOLERANCE_DECLARED" if all_r else "INPUT_INTEGRITY_STOP")),
        }
        if not gate_self_consistent:
            raise RuntimeError(
                f"{label}: files are byte-identical but a sub-gate reports failure -- "
                f"R1={r1['pass']} R2={r2['pass']} R3={r3['pass']}. The gate is defective; "
                f"fix it before trusting any reproduction verdict.")

    all_byte = all(g["byte_identical"] for g in gates.values())
    all_pass = all(g["verdict"] in ("BYTE_IDENTICAL", "TOLERANCE_DECLARED")
                   for g in gates.values())

    out = {
        "stage": "24B",
        "plan_canonical_lf_sha256": a["plan"]["canonical_lf_sha256"],
        "reproduction_root": str(REPRO.relative_to(ROOT)),
        "frozen_artifacts_untouched": {
            k: {"sha256": sha256_file(v), "path": str(v.relative_to(ROOT))}
            for k, v in FROZEN.items() if v.exists()},
        "gates": gates,
        "all_byte_identical": all_byte,
        "reproduction_verdict": (
            "BYTE_IDENTICAL" if all_byte
            else ("TOLERANCE_DECLARED" if all_pass else "INPUT_INTEGRITY_STOP")),
        "R4_cause_named": (
            None if all_byte
            else "REQUIRED: a specific environment difference must be named here before this "
                 "reproduction may be accepted as TOLERANCE_DECLARED. 'Floating point' alone is "
                 "not a cause (plan §7.1 R4)."),
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
        "substage_runtime_minutes": {"23C": round(c_min, 3), "23D": round(d_min, 3)},
    }
    write_json(B_JSON, out)
    return out


# =============================================================================================== #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 24 — Gen-1 Role-B predictor engineering")
    ap.add_argument("--stage", required=True, choices=["24a", "24b"])
    ap.add_argument("--wm989-root", type=Path, default=Path("D:/GSE279162"))
    a = ap.parse_args(argv)

    if a.stage == "24a":
        r = run_24a()
        print(json.dumps({k: r[k] for k in ("stage", "checks", "all_checks_pass")}, indent=2))
    else:
        r = run_24b(a.wm989_root)
        print(json.dumps({k: r[k] for k in
                          ("stage", "all_byte_identical", "reproduction_verdict",
                           "substage_runtime_minutes")}, indent=2))
        for label, g in r["gates"].items():
            print(f"  {label:<12} {g['verdict']:<20} "
                  f"R1={g['R1_shape_and_key']['pass']} R2={g['R2_every_score']['pass']} "
                  f"R3={g['R3_within_clone_ordering']['pass']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
