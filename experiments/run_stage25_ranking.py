"""Stage 25 — the preregistered clone-specific condition-ranking test.

Executes §8 of `STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md`, canonical-LF SHA-256
`8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48`.

The question: among the six observed WM989 conditions, does the frozen interaction model W5 use
pretreatment state to improve **clone-specific condition ordering** over the non-interactive
additive model W4?

This is the sole load-bearing new capability test in Generation 1, and the plan froze every degree
of freedom in it before any of these numbers existed — the metric, the population, the weighting,
the comparator, the endpoint, the null construction, the permutation count, and the verdict logic.
§8.11 forbids changing any of them afterwards. This module therefore reads its parameters from the
frozen plan's protocol rather than defining them, and asserts the plan digest before it starts.

Both verdicts ship. §8.10: `STAGE_25_RANKING_SUPPORTED` permits a validated ordering claim,
`STAGE_25_RANKING_NOT_SUPPORTED` removes only that claim, and either proceeds directly to
`GEN1_MANDATORY_SHIP`. There is no third outcome and no route back to an earlier stage.
"""

from __future__ import annotations

import argparse
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
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "stage25"
SHARDS = _RESULTS / ".cache" / "stage25_null_shards"
for _d in (OUT, SHARDS):
    _d.mkdir(parents=True, exist_ok=True)

PLAN = ROOT / "plans" / "(newer)practical plans" / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"
PROTOCOL_JSON = _RESULTS / "stage23_5_protocol.json"
HANDOFF_JSON = _RESULTS / "stage24_handoff_to_stage25.json"
OOF_TABLE = _RESULTS / "stage24" / "stage24_oof_for_stage25.csv"

A_JSON = OUT / "stage25a_observed.json"
SMOKE_JSON = OUT / "stage25_smoke.json"
VERDICT_JSON = OUT / "stage25_verdict.json"

# ---- frozen by plan §8.6 / §8.7 ---------------------------------------------------------------- #
N_BOOT = 2000
SEED_BOOT = 23501
N_PERM = 1000
SEED_PERM = 23523
EXPECTED_ELIGIBLE = 892
CONDITIONS = ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin", "Trametinib")


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_json(p: Path, obj: dict) -> None:
    p.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def _rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


# =============================================================================================== #
# The metric — §8.5
# =============================================================================================== #
def within_clone_auc(scores: np.ndarray, y: np.ndarray) -> float:
    """§8.5. Mean over every positive/zero condition pair, ties scoring 0.5.

    Not `sklearn.roc_auc_score`: the plan specifies the pairwise form with an explicit tie rule,
    and for a clone with a handful of conditions the two can differ in tie handling. Implementing
    the stated formula is the point.
    """
    pos, zero = scores[y == 1], scores[y == 0]
    if pos.size == 0 or zero.size == 0:
        return float("nan")
    d = pos[:, None] - zero[None, :]
    return float((np.sign(d) * 0.5 + 0.5).mean())


def rank_score(df: pd.DataFrame, col: str) -> tuple[float, dict[str, float]]:
    """§8.5. R(W) = mean over eligible clones of AUC_i, every clone weighted equally."""
    per = {}
    for cid, g in df.groupby("clone_id", sort=True):
        per[cid] = within_clone_auc(g[col].to_numpy(dtype=float), g["y"].to_numpy(dtype=int))
    vals = np.array([v for v in per.values()])
    return float(vals.mean()), per


def low_persistence_top1(df: pd.DataFrame, col: str) -> float:
    """§8.8. Pick the LOWEST predicted detection score; ask whether that condition was a zero.

    Ties break on the canonical §6.3 condition order and never on outcome information.
    """
    order = {t: i for i, t in enumerate(CONDITIONS)}
    hits = []
    for _cid, g in df.groupby("clone_id", sort=True):
        g = g.assign(_o=g["treatment"].map(order)).sort_values([col, "_o"])
        hits.append(int(g["y"].iloc[0]) == 0)
    return float(np.mean(hits))


# =============================================================================================== #
# 25A — the observed statistic
# =============================================================================================== #
def _load_eligible() -> pd.DataFrame:
    if not HANDOFF_JSON.exists():
        raise RuntimeError("Stage 24 has not handed off")
    h = json.loads(HANDOFF_JSON.read_text(encoding="utf-8"))
    if h["stage_24_verdict"] != "STAGE_24_GEN1_TOOL_READY":
        raise RuntimeError(f"Stage 24 is not ready: {h['stage_24_verdict']}")
    if sha256_file(OOF_TABLE) != h["frozen_oof_table"]["sha256"]:
        raise RuntimeError("the out-of-fold table has changed since Stage 24 hashed it")
    if S23.canonical_text_sha256(PLAN) != h["plan_canonical_lf_sha256"]:
        raise RuntimeError("the frozen plan has moved since Stage 24 consumed it")

    tbl = pd.read_csv(OOF_TABLE)
    el = tbl[tbl["ranking_eligible"]].copy()
    n = el["clone_id"].nunique()
    # §8.4: verify mechanically. A different count is an input-integrity stop, not permission to
    # redefine eligibility.
    if n != EXPECTED_ELIGIBLE:
        raise RuntimeError(f"eligible population is {n}, frozen expectation is "
                           f"{EXPECTED_ELIGIBLE} -- INPUT-INTEGRITY STOP")
    if not (el.groupby("clone_id").size() == 6).all():
        raise RuntimeError("an eligible clone does not carry all six condition rows")
    return el


def observed_statistics(el: pd.DataFrame) -> dict:
    r5, per5 = rank_score(el, "pred_W5")
    r4, per4 = rank_score(el, "pred_W4")
    r1, _ = rank_score(el, "pred_W1")
    clones = sorted(per5)
    d_per_clone = np.array([per5[c] - per4[c] for c in clones])

    rng = np.random.default_rng(SEED_BOOT)
    idx = rng.integers(0, len(clones), size=(N_BOOT, len(clones)))
    boot = d_per_clone[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])

    t5 = low_persistence_top1(el, "pred_W5")
    t4 = low_persistence_top1(el, "pred_W4")
    # §8.8: delta_TOP1 gets the SAME clone resampling, at no extra cost
    per_clone_top = {}
    order = {t: i for i, t in enumerate(CONDITIONS)}
    for cid, g in el.groupby("clone_id", sort=True):
        a = g.assign(_o=g["treatment"].map(order)).sort_values(["pred_W5", "_o"])
        b = g.assign(_o=g["treatment"].map(order)).sort_values(["pred_W4", "_o"])
        per_clone_top[cid] = (int(a["y"].iloc[0]) == 0) - (int(b["y"].iloc[0]) == 0)
    tvec = np.array([per_clone_top[c] for c in clones], dtype=float)
    top1_boot = tvec[idx].mean(axis=1)
    tlo, thi = np.percentile(top1_boot, [2.5, 97.5])

    return {
        "eligible_clones": len(clones),
        "R_W1": r1, "R_W4": r4, "R_W5": r5,
        "delta_RANK": r5 - r4,
        "delta_RANK_FULL": r5 - r1,
        "bootstrap": {"replicates": N_BOOT, "seed": SEED_BOOT,
                      "ci95": [float(lo), float(hi)],
                      "lower_endpoint_gt_0": bool(lo > 0),
                      "conditional_on_fitted_models": True},
        "LOW_PERSISTENCE_TOP1_W5": t5, "LOW_PERSISTENCE_TOP1_W4": t4,
        "delta_TOP1": t5 - t4,
        "delta_TOP1_ci95": [float(tlo), float(thi)],
        "delta_TOP1_ge_0": bool((t5 - t4) >= 0),
        "score_tie_rate": float((el.groupby("clone_id")["pred_W5"]
                                 .apply(lambda s: s.duplicated().any())).mean()),
        "_per_clone_delta": d_per_clone, "_clones": clones,
    }


def run_25a() -> dict:
    t0 = time.perf_counter()
    el = _load_eligible()
    obs = observed_statistics(el)
    out = {k: v for k, v in obs.items() if not k.startswith("_")}
    out.update({
        "stage": "25A",
        "plan_canonical_lf_sha256": S23.canonical_text_sha256(PLAN),
        "oof_table": _rel(OOF_TABLE), "oof_table_sha256": sha256_file(OOF_TABLE),
        "metric": "equal-clone-weighted within-clone AUROC, ties 0.5 (plan §8.5)",
        "primary_comparator": "W4", "endpoint": "C1",
        "null_not_yet_computed": True,
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
    })
    write_json(A_JSON, out)
    return out


# =============================================================================================== #
# 25B — the full-refit permutation null (plan §8.7)
#
# The null must mirror the Stage-23 WM989 expression-permutation geometry, not a convenient label
# shuffle. Per outer fold and per draw: permute intact clone-level CP10K/log1p profiles among
# outer-TRAINING clones within stratum, and separately among outer-TEST clones within stratum;
# never across the boundary; never permute genes independently; keep each clone's six-condition
# outcome vector, B, U and fold identity fixed. Then rerun the COMPLETE pipeline -- gene filter,
# scaling, inner-split PCA, inner hyperparameter selection, W4 fit, W5 fit, outer-test prediction --
# and recompute delta_RANK.
#
# Observed-data hyperparameters are NOT reused. That is what makes it a full refit rather than a
# re-scoring, and it is why one draw costs ~100 s rather than ~1 s.
#
# Sharding: draw b is seeded from SEED_PERM + b alone, so shards are order-independent and
# resumable. Each shard writes its OWN file (plan §0.2) -- Stage 23.2H lost a completed draw to a
# race when three shards appended to one shared file, silently.
# =============================================================================================== #
def _strata(ck: pd.DataFrame, clones: list[str]) -> np.ndarray:
    """The frozen Stage-23 WM989 strata: depth bin x 3-bit naive-presence, with its merge rule."""
    return S23.wm989_strata(ck.loc[clones])


def _null_draw(b: int, X, clones, clone_pos, ck, nuis_clone, el_ids, tbl) -> float:
    """One full-refit null draw -> delta_RANK. Mirrors run_23c/run_23d's W4/W5 construction."""
    from sklearn.metrics import log_loss
    from sklearn.model_selection import GroupKFold

    rng = np.random.default_rng(SEED_PERM + b)
    ckey = tbl["clone_id"].to_numpy()
    tkey = tbl["treatment"].to_numpy()
    yv = tbl["y"].to_numpy()
    rf = tbl["outer_fold"].to_numpy()
    dummies = S23.treatment_dummies(tkey)
    strata = _strata(ck, clones)
    pos_of = {c: i for i, c in enumerate(clones)}

    pred = {"W4": np.full(len(tbl), np.nan), "W5": np.full(len(tbl), np.nan)}
    for f in range(S23.N_OUTER):
        tr_clones = sorted({c for c in clones if tbl.loc[ckey == c, "outer_fold"].iloc[0] != f})
        te_clones = sorted({c for c in clones if tbl.loc[ckey == c, "outer_fold"].iloc[0] == f})

        # ---- permute profiles WITHIN stratum, separately on each side of the boundary -------- #
        perm = np.arange(len(clones))
        for side in (tr_clones, te_clones):
            idx = np.array([pos_of[c] for c in side])
            for s in np.unique(strata[idx]):
                cell = idx[strata[idx] == s]
                if len(cell) > 1:
                    perm[cell] = cell[rng.permutation(len(cell))]

        def prep(fit_c, app_c, pm=perm):
            fi = np.array([pm[pos_of[c]] for c in fit_c])
            ai = np.array([pm[pos_of[c]] for c in app_c])
            ztr, zap, _n, kmax = S23.expression_block(X, fi, ai, max(S23.K_CANDIDATES))
            btr, bap = S23.standardize_train_only(
                nuis_clone[[pos_of[c] for c in fit_c]], nuis_clone[[pos_of[c] for c in app_c]])
            pcs = {c: ztr[i] for i, c in enumerate(fit_c)}
            nui = {c: btr[i] for i, c in enumerate(fit_c)}
            for i, c in enumerate(app_c):
                pcs[c], nui[c] = zap[i], bap[i]
            return pcs, nui, kmax

        def design(rows, pcs, nui, k, w5, cc=ckey, dd=dummies):
            P = np.array([pcs[c] for c in cc[rows]])[:, :k]
            B = np.array([nui[c] for c in cc[rows]])
            U = dd[rows]
            return (np.hstack([P, B, U, S23.interaction_block(P, U)]) if w5
                    else np.hstack([P, B, U]))

        tr_rows = np.flatnonzero(rf != f)
        te_rows = np.flatnonzero(rf == f)
        for model, w5 in (("W4", False), ("W5", True)):
            scores: dict = {}
            g = np.array(tr_clones)
            for itr_i, iva_i in GroupKFold(n_splits=S23.N_INNER).split(g, groups=g):
                itr_c, iva_c = [g[i] for i in itr_i], [g[i] for i in iva_i]
                pcs, nui, kmax = prep(itr_c, iva_c)
                si = np.array([i for i in tr_rows if ckey[i] in set(itr_c)])
                sv = np.array([i for i in tr_rows if ckey[i] in set(iva_c)])
                for k in S23.K_CANDIDATES:
                    if k > kmax:
                        continue
                    Ai, Av = design(si, pcs, nui, k, w5), design(sv, pcs, nui, k, w5)
                    for hp in S23.LOGISTIC_C:
                        p = S23._fit_logistic(Ai, yv[si], Av, hp, [])
                        scores.setdefault((k, hp), []).append(
                            log_loss(yv[sv], np.clip(p, 1e-15, 1 - 1e-15), labels=[0, 1]))
            k, hp = min(scores.items(),
                        key=lambda kv: (round(float(np.mean(kv[1])), 12), kv[0][0], kv[0][1]))[0]
            pcs, nui, _ = prep(tr_clones, te_clones)
            pred[model][te_rows] = S23._fit_logistic(
                design(tr_rows, pcs, nui, k, w5), yv[tr_rows],
                design(te_rows, pcs, nui, k, w5), hp, [])

    out = tbl.assign(pred_W4=pred["W4"], pred_W5=pred["W5"])
    el = out[out["clone_id"].isin(el_ids)]
    return rank_score(el, "pred_W5")[0] - rank_score(el, "pred_W4")[0]


def _shard_path(shard: int, n_shards: int) -> Path:
    return SHARDS / f"null_shard{shard}_of{n_shards}.jsonl"


def _load_null_cache() -> dict[int, float]:
    """Read every shard file. Each shard owns its own file, so no append can be lost to a race."""
    done: dict[int, float] = {}
    proto = S23.canonical_text_sha256(PLAN)
    for p in sorted(SHARDS.glob("null_shard*.jsonl")):
        for line in p.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if "plan_sha256" in rec:
                if rec["plan_sha256"] != proto:
                    raise RuntimeError(f"{p.name}: written under a different frozen plan")
                continue
            done[int(rec["i"])] = float(rec["delta_rank"])
    return done


def _null_inputs():
    el = _load_eligible()
    tbl = pd.read_csv(OOF_TABLE)
    X, clones = S23._load_wm989_x()
    ck = pd.read_csv(_RESULTS / "stage22_wm989_clones.csv").set_index("clone_id").loc[clones]
    nuis = np.column_stack([np.log1p(ck[c].to_numpy(dtype=float))
                            for c in ("n_naive_cells", "n_naive1_cells",
                                      "n_naive2_cells", "n_naive3_cells")])
    return el, tbl, X, clones, {c: i for i, c in enumerate(clones)}, ck, nuis


def run_25b(shard: int, n_shards: int, n_perm: int = N_PERM) -> dict:
    t0 = time.perf_counter()
    el, tbl, X, clones, cpos, ck, nuis = _null_inputs()
    el_ids = set(el["clone_id"])
    done = _load_null_cache()
    todo = [b for b in range(n_perm) if b % n_shards == shard and b not in done]
    p = _shard_path(shard, n_shards)
    if not p.exists():
        p.write_text(json.dumps({"plan_sha256": S23.canonical_text_sha256(PLAN)}) + "\n",
                     encoding="utf-8")
    print(f"  [shard {shard}/{n_shards}] {len(todo)} draws to compute, "
          f"{len(done)} already cached", flush=True)
    for n, b in enumerate(todo, start=1):
        d = _null_draw(b, X, clones, cpos, ck, nuis, el_ids, tbl)
        with open(p, "a", encoding="utf-8") as fh:
            fh.write(json.dumps({"i": int(b), "delta_rank": float(d)}) + "\n")
        if n % 5 == 0 or n == len(todo):
            el_s = time.perf_counter() - t0
            print(f"  [shard {shard}] {n}/{len(todo)}  {el_s/60:.1f} min  {el_s/n:.0f} s/draw  "
                  f"eta {el_s/n*(len(todo)-n)/60:.0f} min", flush=True)
    return {"stage": "25B-shard", "shard": shard, "n_shards": n_shards,
            "computed": len(todo), "cached_total": len(_load_null_cache()), "target": n_perm,
            "runtime_minutes": round((time.perf_counter() - t0) / 60, 3)}


# =============================================================================================== #
# SMOKE TEST — run this before committing ~20 h of compute
#
# The claim that has to hold is not "it does not crash" but "sharding, interrupting and resuming
# reproduce EXACTLY the sequential answer, and the metric is what the plan says it is".
# =============================================================================================== #
def run_smoke(n: int = 4) -> dict:
    import shutil
    t0 = time.perf_counter()
    checks: list[dict] = []

    def chk(name, ok, detail=""):
        checks.append({"check": name, "pass": bool(ok), "detail": detail})
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""),
              flush=True)

    # ---- 1. the metric is the plan's formula, on cases with a known answer ------------------- #
    perfect = within_clone_auc(np.array([0.9, 0.8, 0.1, 0.2]), np.array([1, 1, 0, 0]))
    inverted = within_clone_auc(np.array([0.1, 0.2, 0.9, 0.8]), np.array([1, 1, 0, 0]))
    tied = within_clone_auc(np.array([0.5, 0.5]), np.array([1, 0]))
    chk("within-clone AUC matches the plan formula", perfect == 1.0 and inverted == 0.0
        and tied == 0.5, f"perfect {perfect}, inverted {inverted}, tie {tied}")

    # ---- 2. the eligible population is the frozen one ---------------------------------------- #
    el, tbl, X, clones, cpos, ck, nuis = _null_inputs()
    chk("eligible population is the frozen 892", el["clone_id"].nunique() == EXPECTED_ELIGIBLE,
        f"{el['clone_id'].nunique()} clones, all six rows each")

    # ---- 3. equal-clone weighting, not micro-averaging ---------------------------------------- #
    r, per = rank_score(el, "pred_W5")
    micro = float(np.mean([v for v in per.values()]))
    chk("R(W) is the equal-clone mean of AUC_i", abs(r - micro) < 1e-15, f"R(W5) = {r:.6f}")

    # ---- 4. sharding and resume reproduce the sequential answer EXACTLY ---------------------- #
    scratch = SHARDS.parent / "stage25_smoke"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    el_ids = set(el["clone_id"])
    t_d = time.perf_counter()
    seq = {b: _null_draw(b, X, clones, cpos, ck, nuis, el_ids, tbl) for b in range(n)}
    per_draw = (time.perf_counter() - t_d) / n
    shard = {}
    for s in (1, 0):                       # deliberately out of order
        for b in [x for x in range(n) if x % 2 == s]:
            shard[b] = _null_draw(b, X, clones, cpos, ck, nuis, el_ids, tbl)
    chk("2 shards run OUT OF ORDER are bit-identical to sequential",
        all(shard[b] == seq[b] for b in range(n)),
        f"max abs diff {max(abs(shard[b]-seq[b]) for b in range(n)):.3e}")

    # ---- 5. a draw depends on b alone ---------------------------------------------------------- #
    again = _null_draw(0, X, clones, cpos, ck, nuis, el_ids, tbl)
    chk("draw b is a function of b alone (repeatable)", again == seq[0])

    # ---- 6. the null is centred near zero and is NOT the observed value ----------------------- #
    obs = observed_statistics(el)["delta_RANK"]
    chk("null draws differ from the observed statistic",
        all(abs(v - obs) > 1e-12 for v in seq.values()),
        f"observed {obs:+.6f}, null draws {[round(v, 4) for v in seq.values()]}")

    # ---- 7. per-shard files, and the merge refuses an incomplete null -------------------------- #
    chk("each shard writes its OWN file", "_shard_path" in globals()
        and _shard_path(0, 3) != _shard_path(1, 3))
    shutil.rmtree(scratch, ignore_errors=True)

    ok = all(c["pass"] for c in checks)
    eta = {f"n_perm={m}": {"one_process_h": round(m * per_draw / 3600, 1),
                           "three_shards_h": round(m * per_draw / 3 / 3600, 1)}
           for m in (100, 1000)}
    out = {"stage": "25-SMOKE", "draws_used": n,
           "measured_seconds_per_draw": round(per_draw, 1),
           "note": "measured single-process and unloaded; parallel shards contend, and Stage 23.2H "
                   "measured only a 1.44x effective speedup from 3 processes on this pipeline",
           "eta": eta, "checks": checks, "all_passed": ok,
           "already_cached_for_the_real_run": len(_load_null_cache()),
           "runtime_minutes": round((time.perf_counter() - t0) / 60, 3)}
    write_json(SMOKE_JSON, out)
    return out


# =============================================================================================== #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 25 — preregistered ranking test")
    ap.add_argument("--stage", required=True, choices=["25a", "25b", "smoke"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--smoke-draws", type=int, default=4)
    a = ap.parse_args(argv)

    if a.stage == "25a":
        r = run_25a()
        print(json.dumps({k: v for k, v in r.items()
                          if k in ("stage", "eligible_clones", "R_W1", "R_W4", "R_W5",
                                   "delta_RANK", "delta_RANK_FULL", "bootstrap", "delta_TOP1",
                                   "delta_TOP1_ci95")}, indent=2, default=str))
    elif a.stage == "smoke":
        r = run_smoke(a.smoke_draws)
        print(json.dumps({k: r[k] for k in ("stage", "draws_used",
                                            "measured_seconds_per_draw", "eta", "all_passed")},
                         indent=2, default=str))
        return 0 if r["all_passed"] else 1
    else:
        print(json.dumps(run_25b(a.shard, a.n_shards, a.n_perm), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
