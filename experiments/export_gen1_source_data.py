"""Export the per-draw source data behind the Generation-1 result, and pin the environment.

Three things existed only as summary statistics, or only on one machine:

  THE 1,000 NULL DRAWS      lived in a gitignored shard cache. The verdict kept the mean, SD,
                            p95, min, max and the count at-or-above the observed value -- enough
                            to read the result, not enough to re-plot it, recompute the p-value
                            independently, or apply a different test. Regenerating costs 10.7 h.
                            Ten point seven hours of compute, one `rm` from gone.

  THE 2,000 BOOTSTRAP       only the CI endpoints were kept. The replicates are cheap to
  REPLICATES                reproduce (they are conditional on the fitted models and read only
                            the frozen out-of-fold table), but "cheap" is not "recorded".

  THE ENVIRONMENT           `requirements.txt` claimed scikit-learn 1.8.0 while the results were
                            produced under 1.9.0. A pinned file that disagrees with the machine
                            that produced the numbers is worse than no pinned file.

Every export is checked against the recorded statistic it should reproduce. An export that does not
reproduce the verdict is refused rather than written, because a source-data file that disagrees with
the paper is a liability, not an asset.

    python experiments/export_gen1_source_data.py
"""

from __future__ import annotations

import glob
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import run_stage25_ranking as S25  # noqa: E402

RESULTS = ROOT / "results"
SHARDS = RESULTS / ".cache" / "stage25_null_shards"

NULL_CSV = RESULTS / "stage25" / "stage25_null_draws.csv"
BOOT_CSV = RESULTS / "stage25" / "stage25_bootstrap_replicates.csv"
FIGDATA = RESULTS / "manuscript" / "figures" / "figure_source_data.json"
ENVLOCK = ROOT / "environment_lock.txt"

TOL = 1e-12


def _j(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# =============================================================================================== #
def export_null() -> dict:
    """Lift the 1,000 permutation draws out of the shard cache into a committed file."""
    draws: dict[int, float] = {}
    for f in sorted(glob.glob(str(SHARDS / "*.jsonl"))):
        for line in Path(f).read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if "i" in rec:
                draws[int(rec["i"])] = float(rec["delta_rank"])

    v = _j(RESULTS / "stage25" / "stage25_verdict.json")
    p, obs = v["permutation"], v["primary"]["delta_RANK"]

    if len(draws) != p["n_perm"] or sorted(draws) != list(range(p["n_perm"])):
        raise SystemExit(f"REFUSED: recovered {len(draws)} draws, expected a complete "
                         f"0..{p['n_perm'] - 1}. An incomplete null is an integrity stop.")

    vals = np.array([draws[i] for i in range(p["n_perm"])], dtype=float)
    recomputed = {
        "n_null_ge_observed": int((vals >= obs).sum()),
        "null_mean": float(vals.mean()),
        "null_sd": float(vals.std()),
        "null_p95": float(np.percentile(vals, 95)),
        "null_min": float(vals.min()),
        "null_max": float(vals.max()),
    }
    bad = {k: (recomputed[k], p[k]) for k in recomputed
           if abs(recomputed[k] - p[k]) > TOL}
    if bad:
        raise SystemExit(f"REFUSED: the exported draws do not reproduce the verdict: {bad}")

    pd.DataFrame({"draw_index": range(p["n_perm"]), "delta_rank": vals}).to_csv(
        NULL_CSV, index=False, lineterminator="\n")
    return {"file": NULL_CSV.relative_to(ROOT).as_posix(), "n": int(len(vals)),
            "reproduces_verdict": True, **recomputed}


def export_bootstrap() -> dict:
    """Recompute the 2,000 clone-bootstrap replicates with the recorded seed."""
    el = S25._load_eligible()
    _r5, per5 = S25.rank_score(el, "pred_W5")
    _r4, per4 = S25.rank_score(el, "pred_W4")
    clones = sorted(per5)
    d_per_clone = np.array([per5[c] - per4[c] for c in clones], dtype=float)

    rng = np.random.default_rng(S25.SEED_BOOT)
    idx = rng.integers(0, len(clones), size=(S25.N_BOOT, len(clones)))
    boot = d_per_clone[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])

    rec = _j(RESULTS / "stage25" / "stage25_verdict.json")["primary"]["bootstrap_ci95"]
    if abs(lo - rec[0]) > TOL or abs(hi - rec[1]) > TOL:
        raise SystemExit(f"REFUSED: recomputed CI [{lo}, {hi}] != recorded {rec}")

    pd.DataFrame({"replicate": range(S25.N_BOOT), "delta_rank": boot}).to_csv(
        BOOT_CSV, index=False, lineterminator="\n")
    return {"file": BOOT_CSV.relative_to(ROOT).as_posix(), "n": int(S25.N_BOOT),
            "seed": S25.SEED_BOOT, "ci95": [float(lo), float(hi)], "reproduces_verdict": True}


def export_figure_source_data() -> dict:
    """One JSON carrying exactly the numbers each panel draws. Journals ask; it costs nothing."""
    v = _j(RESULTS / "stage25" / "stage25_verdict.json")
    a = _j(RESULTS / "stage26" / "stage26a_vocabulary_closure.json")
    d = v["descriptives"]
    payload = {
        "note": "Every value a figure panel draws. Generated by export_gen1_source_data.py from "
                "the locked verdicts; per-draw distributions are in the two CSVs beside this file "
                "in results/stage25/.",
        "figure_1": {
            "clones_total": v["eligible_clones"] + d["excluded_all_zero_clones"]
            + d["excluded_all_positive_clones"],
            "evaluable": v["eligible_clones"],
            "excluded_never_detected": d["excluded_all_zero_clones"],
            "excluded_always_detected": d["excluded_all_positive_clones"],
            "conditions": list(S25.CONDITIONS),
            "adversarial_refused": a["n_refused"],
            "adversarial_total": a["n_adversarial_strings"],
            "design_columns": a["structural_closure"]["design_columns"],
        },
        "figure_2": {
            "R_W1": v["secondary"]["R_W1"], "R_W4": v["primary"]["R_W4"],
            "R_W5": v["primary"]["R_W5"], "delta_RANK": v["primary"]["delta_RANK"],
            "bootstrap_ci95": v["primary"]["bootstrap_ci95"],
            "null": {k: v["permutation"][k] for k in
                     ("n_perm", "null_mean", "null_sd", "null_p95", "null_min", "null_max",
                      "n_null_ge_observed", "p_perm")},
            "per_draw_null": NULL_CSV.relative_to(ROOT).as_posix(),
            "per_replicate_bootstrap": BOOT_CSV.relative_to(ROOT).as_posix(),
        },
        "figure_3": {
            "by_outer_fold": d["by_outer_fold"],
            "by_pretreatment_depth_bin": d["by_pretreatment_depth_bin"],
            "delta_TOP1": v["delta_TOP1"],
        },
    }
    FIGDATA.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return {"file": FIGDATA.relative_to(ROOT).as_posix(), "panels": 3}


def export_environment() -> dict:
    """The environment that produced these results, not the one a stale file claims."""
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            capture_output=True, text=True).stdout
    header = (
        "# Environment lock for CellFate-Rx Generation 1.\n"
        "#\n"
        "# Captured from the interpreter that produced the locked results. This supersedes\n"
        "# requirements.txt for Generation 1: that file claimed scikit-learn 1.8.0 while the\n"
        "# models were fitted under 1.9.0, and a pinned file that disagrees with the machine\n"
        "# that produced the numbers is worse than no pinned file at all.\n"
        "#\n"
        "# HONEST CAVEAT. This is the environment as it stands now. The project ran over a long\n"
        "# period and nothing recorded the interpreter state at the moment each stage executed,\n"
        "# so this is the best available record, not a retroactive one. Bit-identical\n"
        "# reproduction of the Stage-25 null on a different stack is NOT claimed.\n"
        f"#\n# python   {platform.python_version()}\n"
        f"# platform {platform.platform()}\n"
        f"# machine  {platform.machine()}\n#\n")
    ENVLOCK.write_text(header + freeze, encoding="utf-8")
    key = {ln.split("==")[0].lower(): ln.strip() for ln in freeze.splitlines() if "==" in ln}
    return {"file": ENVLOCK.name, "packages": len(key),
            "python": platform.python_version(),
            "key_versions": {k: key[k] for k in
                             ("numpy", "pandas", "scipy", "scikit-learn") if k in key}}


def main() -> int:
    out = {
        "null_draws": export_null(),
        "bootstrap_replicates": export_bootstrap(),
        "figure_source_data": export_figure_source_data(),
        "environment": export_environment(),
    }
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
