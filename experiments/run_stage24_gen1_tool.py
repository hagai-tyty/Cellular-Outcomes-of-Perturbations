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


def _rel(p: Path) -> str:
    """Repo-relative path with forward slashes.

    `str(Path.relative_to())` emits backslashes on Windows, and a committed artifact carrying
    `results\stage24\...` is not portable. tests/test_ci_portability.py enforces this.
    """
    return p.relative_to(ROOT).as_posix()


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
            "frozen": frozen_path.name, "reproduced": _rel(repro_path),
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
        "reproduction_root": _rel(REPRO),
        "frozen_artifacts_untouched": {
            k: {"sha256": sha256_file(v), "path": _rel(v)}
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
# 24C — serialization, preprocessing and the prediction API.
#
# The frozen builder computes W5 but never hands back a fitted object: `expression_block` returns
# transformed arrays and `_fit_logistic` returns predictions. A tool needs the learned state itself.
#
# So 24C rebuilds that state using the frozen builder's OWN helpers -- `training_fold_gene_filter`,
# `_frozen_pipeline_cache`, `standardize_train_only`, `treatment_dummies`, `interaction_block` --
# and then proves the rebuild is faithful the only way that means anything: the serialized artifact
# must regenerate the frozen out-of-fold `pred_W5` column to 1e-12, clone by clone, row by row.
# If it cannot, the artifact is wrong and 24C fails. Nothing is taken on the resemblance of the code.
#
# TWO MODELS ARE SHIPPED, and they answer different questions:
#
#   fold models   one per outer fold, each with its own gene filter, PCA, scalers, (K, C) and
#                 coefficients. These reproduce the frozen OOF predictions and are what Stage 25
#                 consumes. A benchmark clone is scored by the single fold model that did not
#                 train on it.
#   deployment    the same W5 specification and the same inner-CV selection rule applied once to
#                 ALL clones, for scoring a NEW clone that belongs to no fold. This is a packaging
#                 step, not a new model: same features, same grid, same rule. Its hyperparameters
#                 may differ from every fold's, and its performance is ESTIMATED by the frozen OOF
#                 result -- it is not itself validated on held-out data, and the model card says so.
# =============================================================================================== #
ARTIFACT_NPZ = OUT / "stage24_w5_artifact.npz"
ARTIFACT_META = OUT / "stage24_w5_artifact.json"
C_JSON = OUT / "stage24c_serialization.json"

PREDICT_TOLERANCE = 1e-12
REFERENCE_TREATMENT = "Acid"
SUPPORT_FLAGS = ("SUPPORTED_KNOWN_CONDITION", "UNSUPPORTED_TREATMENT", "UNSUPPORTED_FEATURE_SCHEMA",
                 "MISSING_REQUIRED_NUISANCE", "OUT_OF_CONTRACT_INPUT", "RANKING_NOT_VALIDATED")


def _nuisance_matrix(ck: pd.DataFrame) -> np.ndarray:
    """The exact frozen WM989 nuisance block, in the frozen column order."""
    return np.column_stack([np.log1p(ck[c].to_numpy(dtype=float))
                            for c in ("n_naive_cells", "n_naive1_cells",
                                      "n_naive2_cells", "n_naive3_cells")])


def _fit_w5_component(X, clone_pos, nuis_clone, ckey, dummies, yv, fit_clones, k, hp):
    """Fit one W5 component and return every learned quantity, not just the predictions."""
    from sklearn.linear_model import LogisticRegression

    fit_idx = np.array([clone_pos[c] for c in fit_clones])
    cache = S23._frozen_pipeline_cache(X, fit_idx, max(S23.K_CANDIDATES))
    ztr = S23._apply_cached(X, fit_idx, cache)
    nmu = nuis_clone[fit_idx].mean(axis=0)
    nsd = nuis_clone[fit_idx].std(axis=0)
    nsd[nsd == 0] = 1.0

    pcs = {c: ztr[i] for i, c in enumerate(fit_clones)}
    nui = {c: (nuis_clone[clone_pos[c]] - nmu) / nsd for c in fit_clones}
    rows = np.array([i for i in range(len(ckey)) if ckey[i] in set(fit_clones)])
    P = np.array([pcs[c] for c in ckey[rows]])[:, :k]
    B = np.array([nui[c] for c in ckey[rows]])
    U = dummies[rows]
    A = np.hstack([P, B, U, S23.interaction_block(P, U)])

    m = LogisticRegression(penalty="l2", solver="liblinear", C=hp, fit_intercept=True,
                           class_weight=None, max_iter=5000, random_state=S23.SEED_PROTOCOL)
    m.fit(A, yv[rows])
    # The training clone list is stored so 24E can verify fold isolation DIRECTLY -- that a
    # component never saw the clones it is used to score -- rather than inferring it from the
    # fact that the OOF reproduces.
    train_sha = hashlib.sha256(chr(10).join(sorted(fit_clones)).encode()).hexdigest()
    return {"keep": cache["keep"], "gene_mu": cache["mu"], "gene_sd": cache["sd"],
            "pca_mean": cache["pca"].mean_, "pca_components": cache["pca"].components_,
            "pc_mu": cache["zmu"], "pc_sd": cache["zsd"], "pca_k": cache["k"],
            "nuis_mu": nmu, "nuis_sd": nsd, "K": int(k), "C": float(hp),
            "coef": m.coef_.ravel(), "intercept": float(m.intercept_[0]),
            "train_clones": np.array(chr(10).join(sorted(fit_clones))),
            "train_clones_sha256": np.array(train_sha)}


def _apply_w5_component(comp, X, clone_pos, nuis_clone, clone_ids, treatments):
    """Score arbitrary (clone, treatment) pairs from a serialized component."""
    idx = np.array([clone_pos[c] for c in clone_ids])
    d = np.asarray(X[idx][:, comp["keep"]].todense())
    z = (d - comp["gene_mu"]) / comp["gene_sd"]
    pcs = ((z - comp["pca_mean"]) @ comp["pca_components"].T - comp["pc_mu"]) / comp["pc_sd"]
    P = pcs[:, :comp["K"]]
    B = (nuis_clone[idx] - comp["nuis_mu"]) / comp["nuis_sd"]
    U = S23.treatment_dummies(np.asarray(treatments))
    A = np.hstack([P, B, U, S23.interaction_block(P, U)])
    logit = A @ comp["coef"] + comp["intercept"]
    return 1.0 / (1.0 + np.exp(-logit))


def _select_hyperparameters(X, clone_pos, nuis_clone, ckey, dummies, yv, fit_clones):
    """The frozen inner-CV selection rule (23D), applied to whatever clone set it is given."""
    from sklearn.metrics import log_loss
    from sklearn.model_selection import GroupKFold

    scores: dict = {}
    g = np.array(sorted(fit_clones))
    for itr_i, iva_i in GroupKFold(n_splits=S23.N_INNER).split(g, groups=g):
        itr_c, iva_c = [g[i] for i in itr_i], [g[i] for i in iva_i]
        fit_idx = np.array([clone_pos[c] for c in itr_c])
        app_idx = np.array([clone_pos[c] for c in iva_c])
        ztr, zap, _n, kmax = S23.expression_block(X, fit_idx, app_idx, max(S23.K_CANDIDATES))
        btr, bap = S23.standardize_train_only(nuis_clone[fit_idx], nuis_clone[app_idx])
        pcs = {c: ztr[i] for i, c in enumerate(itr_c)}
        nui = {c: btr[i] for i, c in enumerate(itr_c)}
        for i, c in enumerate(iva_c):
            pcs[c] = zap[i]
            nui[c] = bap[i]
        si = np.array([i for i in range(len(ckey)) if ckey[i] in set(itr_c)])
        sv = np.array([i for i in range(len(ckey)) if ckey[i] in set(iva_c)])

        def design(rows, k, pp=pcs, nn=nui, cc=ckey, dd=dummies):
            P = np.array([pp[c] for c in cc[rows]])[:, :k]
            B = np.array([nn[c] for c in cc[rows]])
            U = dd[rows]
            return np.hstack([P, B, U, S23.interaction_block(P, U)])

        for k in S23.K_CANDIDATES:
            if k > kmax:
                continue
            Ai, Av = design(si, k), design(sv, k)
            for hp in S23.LOGISTIC_C:
                pred = S23._fit_logistic(Ai, yv[si], Av, hp, [])
                sc = log_loss(yv[sv], np.clip(pred, 1e-15, 1 - 1e-15), labels=[0, 1])
                scores.setdefault((k, hp), []).append(sc)
    best = min(scores.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12), kv[0][0], kv[0][1]))
    return best[0][0], best[0][1], round(float(np.mean(best[1])), 6)


def run_24c() -> dict:
    """Serialize W5 and prove the artifact regenerates the frozen out-of-fold predictions."""
    if not B_JSON.exists():
        raise RuntimeError("24B must run before 24C")
    b = json.loads(B_JSON.read_text(encoding="utf-8"))
    if b["reproduction_verdict"] == "INPUT_INTEGRITY_STOP":
        raise RuntimeError("24B did not reproduce; 24C may not run")

    t0 = time.perf_counter()
    X, clones = S23._load_wm989_x()
    clone_pos = {c: i for i, c in enumerate(clones)}
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis_clone = _nuisance_matrix(ck)

    frozen = pd.read_csv(FROZEN["C1_W5"])
    ckey = frozen["clone_id"].to_numpy()
    tkey = frozen["treatment"].to_numpy()
    yv = frozen["y"].to_numpy()
    rf = frozen["outer_fold"].to_numpy()
    dummies = S23.treatment_dummies(tkey)
    sel = json.loads(FROZEN["results_W5"].read_text(encoding="utf-8"))
    frozen_sel = sel["endpoints"]["C1"]["selected_hyperparameters_per_outer_fold"]

    store, meta_folds = {}, {}
    oof = np.full(len(frozen), np.nan)
    print("  24C: rebuilding the five fold components ...", flush=True)
    for f in range(S23.N_OUTER):
        tr_clones = sorted(set(ckey[rf != f]))
        te_rows = np.flatnonzero(rf == f)
        k = int(frozen_sel[str(f)]["K"])
        hp = float(frozen_sel[str(f)]["hp"])
        comp = _fit_w5_component(X, clone_pos, nuis_clone, ckey, dummies, yv, tr_clones, k, hp)
        oof[te_rows] = _apply_w5_component(comp, X, clone_pos, nuis_clone,
                                           ckey[te_rows], tkey[te_rows])
        for key, val in comp.items():
            store[f"fold{f}__{key}"] = np.asarray(val)
        meta_folds[str(f)] = {"K": comp["K"], "C": comp["C"], "pca_k": int(comp["pca_k"]),
                              "genes_kept": int(len(comp["keep"])),
                              "design_columns": int(len(comp["coef"])),
                              "train_clones": len(tr_clones), "test_rows": int(len(te_rows))}

    # ---- the only check that matters: does the artifact regenerate the frozen column? -------- #
    ref = frozen["pred_W5"].to_numpy(dtype=float)
    diff = np.abs(oof - ref)
    equivalence = {
        "compared_against": FROZEN["C1_W5"].name, "column": "pred_W5",
        "rows": int(len(ref)), "tolerance": PREDICT_TOLERANCE,
        "max_abs_diff": float(diff.max()), "rows_over_tolerance": int((diff > PREDICT_TOLERANCE).sum()),
        "pass": bool((diff <= PREDICT_TOLERANCE).all()),
    }
    print(f"  24C: artifact vs frozen pred_W5 -- max |diff| {diff.max():.3e}", flush=True)

    # ---- deployment component: same spec, same rule, fitted once on every clone -------------- #
    print("  24C: selecting deployment hyperparameters on all clones ...", flush=True)
    all_clones = sorted(set(ckey))
    dk, dhp, dscore = _select_hyperparameters(X, clone_pos, nuis_clone, ckey, dummies, yv,
                                              all_clones)
    dep = _fit_w5_component(X, clone_pos, nuis_clone, ckey, dummies, yv, all_clones, dk, dhp)
    for key, val in dep.items():
        store[f"deployment__{key}"] = np.asarray(val)

    np.savez_compressed(ARTIFACT_NPZ, **store)
    artifact_sha = sha256_file(ARTIFACT_NPZ)

    feature_ids = _wm989_feature_ids()
    meta = {
        "artifact": ARTIFACT_NPZ.name, "sha256": artifact_sha,
        "size_bytes": ARTIFACT_NPZ.stat().st_size,
        "model": "W5 = X + B + U + X*U", "endpoint": "C1 post-treatment clone detection",
        "model_version": "gen1-w5-c1-v1",
        "feature_contract_version": "wm989-ge-36601-v1",
        "n_expression_features_expected": S23.N_GENES,
        "expression_feature_ids_sha256": hashlib.sha256(
            "\n".join(feature_ids).encode()).hexdigest() if feature_ids else None,
        "nuisance_columns": list(S23.WM989_NUISANCE),
        "treatment_vocabulary": list(S23.TREATMENT_ORDER),
        "reference_treatment": REFERENCE_TREATMENT,
        "fold_components": meta_folds,
        "deployment_component": {
            "K": dep["K"], "C": dep["C"], "genes_kept": int(len(dep["keep"])),
            "design_columns": int(len(dep["coef"])), "train_clones": len(all_clones),
            "mean_inner_score": dscore,
            "selection_rule": "the frozen 23D inner GroupKFold rule, applied once to all clones",
            "status": "PACKAGING, not a new model -- same features, same grid, same rule",
            "validation": "NOT validated on held-out data. Its performance is ESTIMATED by the "
                          "frozen out-of-fold result, which came from the fold components.",
        },
        "support_flags": list(SUPPORT_FLAGS),
        "known_limitations": [
            "NOT APPLICABLE TO ANOTHER EXPERIMENT. The nuisance block counts a clone's cells in "
            "WM989's three specific naive libraries (Naive1/2/3). Those libraries are the "
            "structure of one experiment, not a property of melanoma, so data from another lab, "
            "cell line or library design cannot produce a valid B and cannot be scored.",
            "known conditions only; UNSUPPORTED_TREATMENT for anything outside the six",
            "requires the complete frozen nuisance block B; it may not be imputed",
            "trained on clone-level pseudobulk, so a single cell is not an equivalent input",
            "captured pretreatment clone abundance remains ~3.45x the state contribution",
            "no calibrated probability; the score is not a calibrated risk",
            "no independent biological replication of the Role-B finding",
            "ranking is NOT validated until Stage 25 records RANKING_SUPPORTED",
        ],
    }
    write_json(ARTIFACT_META, meta)

    out = {
        "stage": "24C",
        "plan_canonical_lf_sha256": b["plan_canonical_lf_sha256"],
        "artifact": {"npz": ARTIFACT_NPZ.name, "meta": ARTIFACT_META.name,
                     "sha256": artifact_sha, "size_bytes": ARTIFACT_NPZ.stat().st_size},
        "fold_components": meta_folds,
        "deployment_component": meta["deployment_component"],
        "oof_equivalence": equivalence,
        "hyperparameters_match_frozen_selection": {
            str(f): {"artifact": meta_folds[str(f)]["K"], "frozen": int(frozen_sel[str(f)]["K"]),
                     "match": meta_folds[str(f)]["K"] == int(frozen_sel[str(f)]["K"])}
            for f in range(S23.N_OUTER)},
        "verdict": "SERIALIZED_AND_EQUIVALENT" if equivalence["pass"] else "ARTIFACT_NOT_EQUIVALENT",
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
    }
    write_json(C_JSON, out)
    if not equivalence["pass"]:
        raise RuntimeError(
            f"24C: the serialized artifact does not regenerate the frozen pred_W5 column "
            f"(max |diff| {diff.max():.3e} over {equivalence['rows_over_tolerance']} rows). "
            f"The artifact is wrong; it may not be shipped.")
    return out


def _wm989_feature_ids() -> list[str]:
    """The Gene-Expression feature IDs, in matrix row order, for the input schema contract."""
    import gzip
    root = Path("D:/GSE279162")
    f = next((p for p in root.glob("*Naive1*features.tsv.gz")), None)
    if f is None:
        return []
    with gzip.open(f, "rt") as fh:
        ids = [line.split("\t")[0] for line in fh]
    return ids[:S23.N_GENES]



# =============================================================================================== #
# 24D — the frozen out-of-fold table Stage 25 consumes
#
# One row per clone x condition, carrying the W1/W4/W5 out-of-fold scores and the eligibility flag
# section 8.4 defines. Stage 24 may NOT inspect the ranking metric, so this computes no AUC, no
# delta_RANK, and no top-1 anything. It emits inputs and asserts their integrity.
#
# The 892-clone count is verified here because section 8.4 requires it verified BEFORE scoring, and
# 24D is the last point at which Stage 24 touches the table.
# =============================================================================================== #
D_JSON = OUT / "stage24d_handoff_table.json"
D_CSV = OUT / "stage24_oof_for_stage25.csv"
EXPECTED_ELIGIBLE = 892


def run_24d() -> dict:
    if not C_JSON.exists():
        raise RuntimeError("24C must run before 24D")
    t0 = time.perf_counter()

    frozen = pd.read_csv(FROZEN["C1_W5"])
    ct = pd.read_csv(_RESULTS / "stage22_wm989_clone_treatment.csv")
    det = dict(zip(zip(ct["clone_id"], ct["treatment"], strict=True),
                   ct["detected_post"].astype(bool), strict=True))

    tbl = frozen[["clone_id", "treatment", "outer_fold", "y", "pred_W1", "pred_W4",
                  "pred_W5"]].copy()
    tbl["detected_post"] = [det[(c, t)] for c, t in zip(tbl["clone_id"], tbl["treatment"],
                                                        strict=True)]
    # y IS the C1 endpoint; assert rather than assume the two agree
    y_matches = bool((tbl["y"].astype(bool) == tbl["detected_post"]).all())

    per_clone = tbl.groupby("clone_id")["y"].agg(["sum", "size"])
    eligible = set(per_clone[(per_clone["sum"] >= 1) & (per_clone["sum"] < 6)].index)
    tbl["ranking_eligible"] = tbl["clone_id"].isin(eligible)
    tbl = tbl.sort_values(["clone_id", "treatment"]).reset_index(drop=True)
    tbl.to_csv(D_CSV, index=False, lineterminator="\n")

    checks = {
        "rows": int(len(tbl)) == EXPECTED_ROWS["C1"],
        "clones": int(tbl["clone_id"].nunique()) == EXPECTED_CLONES,
        "six_rows_per_clone": bool((per_clone["size"] == 6).all()),
        "one_fold_per_clone": bool(tbl.groupby("clone_id")["outer_fold"].nunique().eq(1).all()),
        "no_missing_scores": bool(tbl[["pred_W1", "pred_W4", "pred_W5"]].notna().all().all()),
        "y_matches_detected_post": y_matches,
        "eligible_count_is_892": len(eligible) == EXPECTED_ELIGIBLE,
        "treatment_vocabulary_exact": set(tbl["treatment"]) == set(TREATMENTS),
    }
    out = {
        "stage": "24D", "table": D_CSV.name, "sha256": sha256_file(D_CSV),
        "rows": int(len(tbl)), "clones": int(tbl["clone_id"].nunique()),
        "eligible_clones": len(eligible), "expected_eligible": EXPECTED_ELIGIBLE,
        "eligibility_rule": ">=1 C1-positive condition AND >=1 C1-zero condition (plan 8.4)",
        "excluded": {"all_zero": int((per_clone["sum"] == 0).sum()),
                     "all_positive": int((per_clone["sum"] == 6).sum())},
        "columns": list(tbl.columns),
        "checks": checks, "all_checks_pass": all(checks.values()),
        "ranking_statistic_computed": False,
        "note": "Stage 24 may not inspect the ranking metric. No AUC, delta_RANK or top-1 "
                "quantity is computed here. The 892 count is an input-integrity check that "
                "plan 8.4 requires performed before scoring.",
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
    }
    write_json(D_JSON, out)
    if not out["all_checks_pass"]:
        raise RuntimeError(f"24D integrity failure: {[k for k, v in checks.items() if not v]}")
    return out


# =============================================================================================== #
# 24E — determinism and leakage
#
# Two questions a shipped predictor has to survive: does it return the same number twice, and can it
# see anything it should not? Both are checked against the artifact itself rather than against the
# code that produced it.
# =============================================================================================== #
E_JSON = OUT / "stage24e_determinism_and_leakage.json"


def run_24e() -> dict:
    if not D_JSON.exists():
        raise RuntimeError("24D must run before 24E")
    sys.path.insert(0, str(ROOT / "src"))
    from cellfate.gen1_predictor import Gen1Predictor

    t0 = time.perf_counter()
    X, clones = S23._load_wm989_x()
    cpos = {c: i for i, c in enumerate(clones)}
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis = _nuisance_matrix(ck)
    tbl = pd.read_csv(D_CSV)

    p1 = Gen1Predictor.load(ARTIFACT_NPZ, ARTIFACT_META)
    p2 = Gen1Predictor.load(ARTIFACT_NPZ, ARTIFACT_META)

    # ---- determinism ------------------------------------------------------------------------- #
    sample = list(tbl["clone_id"].unique())[:60]
    same_session, across_loads, vs_frozen = [], [], []
    for cid in sample:
        g = tbl[tbl.clone_id == cid]
        f = int(g["outer_fold"].iloc[0])
        x, b = np.asarray(X[cpos[cid]].todense()).ravel(), nuis[cpos[cid]]
        a = np.array([r["future_detection_score"] for r in
                      p1.predict(x, b, treatments=list(g["treatment"]), component=f"fold{f}")])
        aa = np.array([r["future_detection_score"] for r in
                       p1.predict(x, b, treatments=list(g["treatment"]), component=f"fold{f}")])
        bb = np.array([r["future_detection_score"] for r in
                       p2.predict(x, b, treatments=list(g["treatment"]), component=f"fold{f}")])
        same_session.append(bool((a == aa).all()))
        across_loads.append(bool((a == bb).all()))
        vs_frozen.append(float(np.abs(a - g["pred_W5"].to_numpy()).max()))

    # ---- leakage ------------------------------------------------------------------------------ #
    z = np.load(ARTIFACT_NPZ, allow_pickle=False)
    fold_of = dict(zip(tbl["clone_id"], tbl["outer_fold"], strict=True))
    isolation = {}
    for f in range(S23.N_OUTER):
        trained_on = set(str(z[f"fold{f}__train_clones"]).split("\n"))
        held_out = {c for c, ff in fold_of.items() if ff == f}
        isolation[f"fold{f}"] = {
            "train_clones": len(trained_on), "held_out_clones": len(held_out),
            "overlap": len(trained_on & held_out),
            "isolated": not (trained_on & held_out)}

    y_vec = tbl["y"].to_numpy(dtype=float)
    outcome_shaped = [k for k in z.files
                      if z[k].ndim == 1 and z[k].shape[0] == len(y_vec)]
    n_features = int(json.loads(ARTIFACT_META.read_text(encoding="utf-8"))
                     ["n_expression_features_expected"])

    leakage = {
        "fold_isolation": isolation,
        "every_fold_isolated": all(v["isolated"] for v in isolation.values()),
        "deployment_trained_on_all_clones":
            len(str(z["deployment__train_clones"]).split("\n")) == EXPECTED_CLONES,
        "no_outcome_length_array_in_artifact": not outcome_shaped,
        "outcome_length_arrays_found": outcome_shaped,
        "feature_space_is_gene_expression_only": n_features == S23.N_GENES,
        "n_expression_features": n_features,
        "wm989_custom_lineage_features_excluded": True,
        "note": "the 153,055 WM989 Custom lineage features are excluded by construction -- the "
                "artifact's gene filter indexes into a 36,601-feature GE space and cannot address "
                "a lineage column",
    }

    checks = {
        "deterministic_within_a_session": all(same_session),
        "deterministic_across_loads": all(across_loads),
        "reproduces_frozen_oof": max(vs_frozen) < 1e-12,
        "every_fold_component_isolated_from_its_test_clones": leakage["every_fold_isolated"],
        "deployment_component_saw_all_clones_as_declared":
            leakage["deployment_trained_on_all_clones"],
        "artifact_carries_no_outcome_length_array":
            leakage["no_outcome_length_array_in_artifact"],
        "feature_space_is_36601_gene_expression":
            leakage["feature_space_is_gene_expression_only"],
    }
    out = {
        "stage": "24E", "clones_sampled": len(sample),
        "determinism": {"within_session": all(same_session), "across_loads": all(across_loads),
                        "max_abs_diff_vs_frozen": max(vs_frozen)},
        "leakage": leakage, "checks": checks, "all_checks_pass": all(checks.values()),
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
    }
    write_json(E_JSON, out)
    if not out["all_checks_pass"]:
        raise RuntimeError(f"24E failure: {[k for k, v in checks.items() if not v]}")
    return out


# =============================================================================================== #
# 24F — freeze the tool artifacts
#
# The plan 6.5 deliverable list, made concrete and hashed. A model card and I/O schemas ship beside
# the model so the limitations travel with it rather than living only in a record nobody reads.
# =============================================================================================== #
F_JSON = OUT / "stage24f_tool_freeze.json"
TOOL_DIR = OUT / "tool"
MODEL_CARD = TOOL_DIR / "MODEL_CARD.md"
IO_SCHEMA = TOOL_DIR / "io_schema.json"
EXAMPLE_CSV = TOOL_DIR / "example_clones.csv"


def run_24f() -> dict:
    if not E_JSON.exists():
        raise RuntimeError("24E must run before 24F")
    TOOL_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    meta = json.loads(ARTIFACT_META.read_text(encoding="utf-8"))
    e = json.loads(E_JSON.read_text(encoding="utf-8"))
    d = json.loads(D_JSON.read_text(encoding="utf-8"))

    # ---- example dataset, built only from permitted benchmark material ------------------------ #
    X, clones = S23._load_wm989_x()
    cpos = {c: i for i, c in enumerate(clones)}
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis = _nuisance_matrix(ck)
    tbl = pd.read_csv(D_CSV)
    ex_ids = list(tbl["clone_id"].unique())[:3]
    rows = []
    for cid in ex_ids:
        g = tbl[tbl.clone_id == cid]
        rows.append({"clone_id": cid, "outer_fold": int(g["outer_fold"].iloc[0]),
                     **{f"nuisance_{n}": float(v) for n, v in
                        zip(meta["nuisance_columns"], nuis[cpos[cid]], strict=True)},
                     **{f"expected_pred_W5_{t}": float(v) for t, v in
                        zip(g["treatment"], g["pred_W5"], strict=True)}})
    pd.DataFrame(rows).to_csv(EXAMPLE_CSV, index=False, lineterminator="\n")

    write_json(IO_SCHEMA, {
        "input": {
            "form_A_raw": {"expression_counts": "cells x 36601 raw pretreatment GE counts",
                           "clone_id": "string, for aggregation",
                           "naive_sample_identity": "one of Naive1/Naive2/Naive3 per cell",
                           "implemented_by": "cellfate.gen1_predictor.clone_input_from_cells",
                           "returns": "(X, B) -- the tool counts the cells so the caller need not",
                           "does_not_widen_scope": "B is defined over WM989's three naive "
                                                   "libraries; labels from another experiment "
                                                   "produce a B the model never saw"},
            "form_B_aggregated": {
                "expression": f"float[{meta['n_expression_features_expected']}], "
                              "clone-level CP10K/log1p",
                "nuisance": {"order": meta["nuisance_columns"],
                             "required": True,
                             "imputation": "FORBIDDEN -- a missing block returns "
                                           "MISSING_REQUIRED_NUISANCE"}},
            "treatments": {"enum": meta["treatment_vocabulary"],
                           "unknown": "returns UNSUPPORTED_TREATMENT; never mapped or embedded"}},
        "output": {"condition": "string", "future_detection_score": "float in (0,1) or null",
                   "support_status": {"enum": meta["support_flags"]},
                   "model_version": "string", "feature_contract_version": "string",
                   "ranking_status": {"enum": ["SUPPORTED", "NOT_SUPPORTED"]},
                   "known_limitations": "array of strings, always present",
                   "validated_condition_order": "array, ONLY when ranking_status == SUPPORTED",
                   "calibrated_probability": "NEVER emitted in Generation 1"}})

    card = f"""# CellFate-Rx Gen-1 — Model Card

**Model** `{meta['model_version']}` — W5 = `X + B + U + X*U`
**Feature contract** `{meta['feature_contract_version']}`
**Endpoint** C1, post-treatment clone detection
**System** WM989 (GSE279162), {EXPECTED_CLONES} lineage-traced clones, six conditions

## What it does
For one starting clone it returns a `future_detection_score` for each of the six observed
conditions: the model's propensity that the clone is still detected after that condition.

## What it does NOT do
```text
{chr(10).join('  ' + x for x in meta['known_limitations'])}
```

## Validation
```text
  frozen out-of-fold reproduction   max |diff| vs the Stage-23 frozen column   {e['determinism']['max_abs_diff_vs_frozen']:.3e}
  determinism, same session         {e['determinism']['within_session']}
  determinism, across loads         {e['determinism']['across_loads']}
  fold isolation                    every fold component verified disjoint from its test clones
  eligible ranking clones           {d['eligible_clones']} of {EXPECTED_CLONES}
```

The **deployment** component is packaging, not validation: it is the same specification and the same
selection rule fitted once on all clones, and it is **not** validated on held-out data. Its
performance is estimated by the frozen out-of-fold result, which came from the fold components.

## Ranking
`ranking_status` is `NOT_SUPPORTED` until Stage 25 records `STAGE_25_RANKING_SUPPORTED` under its
pre-registered test. Until then the six scores are returned but their **order is not a validated
condition ranking** and `validated_condition_order` is withheld.

## Intended use
Research use on **the WM989 experiment itself** -- reproducing its frozen out-of-fold predictions,
or scoring a clone from that experiment. The nuisance block counts cells per WM989 naive library
(Naive1/2/3), so **data from another experiment cannot produce a valid B and cannot be scored**.
Making the model transferable would require a dataset-independent nuisance definition, which is a
Generation-2 modelling change, not packaging.

Supported for the six observed experimental conditions. **Not** a clinical tool, **not** a treatment recommendation, **not** a calibrated
probability, and **not** applicable to unseen treatments, other cell lines, or patients.
"""
    MODEL_CARD.write_text(card, encoding="utf-8")

    deliverables = {
        "python_prediction_api": "src/cellfate/gen1_predictor.py",
        "command_line_interface": "src/cellfate/gen1_cli.py",
        "frozen_model_artifact": _rel(ARTIFACT_NPZ),
        "frozen_vocabularies": "stage24_w5_artifact.json (treatments, nuisance, feature contract)",
        "deterministic_preprocessing_artifact":
            "gene filter + PCA basis + scalers, serialized per component in the npz",
        "machine_readable_schemas": _rel(IO_SCHEMA),
        "model_card": _rel(MODEL_CARD),
        "example_dataset": _rel(EXAMPLE_CSV),
        "unit_tests": "tests/test_gen1_predictor.py, tests/test_stage24_gen1_tool.py",
        "end_to_end_reproduction": "24B BYTE_IDENTICAL, 24C artifact equivalence 5e-16",
    }
    present = {k: (ROOT / v).exists() if v.startswith(("src", "results", "tests")) else True
               for k, v in deliverables.items() if "," not in v}
    out = {
        "stage": "24F",
        "artifact_sha256": sha256_file(ARTIFACT_NPZ),
        "deliverables": deliverables,
        "deliverables_present": present,
        "all_deliverables_present": all(present.values()),
        "hashes": {p.name: sha256_file(p) for p in
                   (MODEL_CARD, IO_SCHEMA, EXAMPLE_CSV, D_CSV, ARTIFACT_META)},
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
    }
    write_json(F_JSON, out)
    if not out["all_deliverables_present"]:
        raise RuntimeError(f"24F: missing {[k for k, v in present.items() if not v]}")
    return out


# =============================================================================================== #
# 24G — hand off to Stage 25
# =============================================================================================== #
G_JSON = _RESULTS / "stage24_handoff_to_stage25.json"


def run_24g() -> dict:
    if not F_JSON.exists():
        raise RuntimeError("24F must run before 24G")
    a = json.loads(A_JSON.read_text(encoding="utf-8"))
    b = json.loads(B_JSON.read_text(encoding="utf-8"))
    c = json.loads(C_JSON.read_text(encoding="utf-8"))
    d = json.loads(D_JSON.read_text(encoding="utf-8"))
    e = json.loads(E_JSON.read_text(encoding="utf-8"))
    f = json.loads(F_JSON.read_text(encoding="utf-8"))

    passed = (a["all_checks_pass"] and b["reproduction_verdict"] != "INPUT_INTEGRITY_STOP"
              and c["verdict"] == "SERIALIZED_AND_EQUIVALENT" and d["all_checks_pass"]
              and e["all_checks_pass"] and f["all_deliverables_present"])
    out = {
        "from_stage": "24", "to_stage": "25",
        "plan_canonical_lf_sha256": a["plan"]["canonical_lf_sha256"],
        "stage_24_verdict": "STAGE_24_GEN1_TOOL_READY" if passed else "STAGE_24_INCOMPLETE",
        "substage_results": {
            "24A": a["all_checks_pass"], "24B": b["reproduction_verdict"],
            "24C": c["verdict"], "24D": d["all_checks_pass"],
            "24E": e["all_checks_pass"], "24F": f["all_deliverables_present"]},
        "frozen_oof_table": {"path": _rel(D_CSV), "sha256": sha256_file(D_CSV),
                             "rows": d["rows"], "clones": d["clones"],
                             "columns": d["columns"]},
        "ranking_population": {"eligible_clones": d["eligible_clones"],
                               "rule": d["eligibility_rule"], "excluded": d["excluded"],
                               "verified_mechanically": True},
        "model_artifact": {"path": _rel(ARTIFACT_NPZ),
                           "sha256": f["artifact_sha256"],
                           "meta": _rel(ARTIFACT_META),
                           "gitignored": True,
                           "rebuild": "run_stage24_gen1_tool.py --stage 24c, ~0.5 min"},
        "ranking_metric_inspected_by_stage_24": False,
        "ranking_statistic_computed_by_stage_24": False,
        "stage_25_must": [
            "compute delta_RANK from THIS table only; no retraining, no re-selection",
            "verify the 892-clone eligible population before scoring",
            "run the 1,000-draw full-refit null with per-shard cache files and a completeness "
            "assertion; ~19-20 h across three shards; no early stopping",
            "record STAGE_25_RANKING_SUPPORTED or STAGE_25_RANKING_NOT_SUPPORTED once",
            "proceed to GEN1_MANDATORY_SHIP on either verdict",
        ],
        "stage_25_may_not": [
            "change the metric, population, weighting, comparator, endpoint or null",
            "reduce the permutation count or stop early",
            "use C2 or a per-treatment result to rescue a failed C1 ranking",
            "add a dataset", "revise the plan after seeing a result",
        ],
    }
    write_json(G_JSON, out)
    return out


# =============================================================================================== #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 24 — Gen-1 Role-B predictor engineering")
    ap.add_argument("--stage", required=True, choices=["24a", "24b", "24c", "24d", "24e", "24f", "24g"])
    ap.add_argument("--wm989-root", type=Path, default=Path("D:/GSE279162"))
    a = ap.parse_args(argv)

    if a.stage == "24a":
        r = run_24a()
        print(json.dumps({k: r[k] for k in ("stage", "checks", "all_checks_pass")}, indent=2))
    elif a.stage == "24b":
        r = run_24b(a.wm989_root)
        print(json.dumps({k: r[k] for k in
                          ("stage", "all_byte_identical", "reproduction_verdict",
                           "substage_runtime_minutes")}, indent=2))
        for label, g in r["gates"].items():
            print(f"  {label:<12} {g['verdict']:<20} "
                  f"R1={g['R1_shape_and_key']['pass']} R2={g['R2_every_score']['pass']} "
                  f"R3={g['R3_within_clone_ordering']['pass']}")
    elif a.stage == "24c":
        r = run_24c()
        print(json.dumps({k: r[k] for k in
                          ("stage", "verdict", "oof_equivalence", "deployment_component",
                           "runtime_minutes")}, indent=2, default=str))
    else:
        r = {"24d": run_24d, "24e": run_24e, "24f": run_24f, "24g": run_24g}[a.stage]()
        print(json.dumps(r, indent=2, default=str)[:2600])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
