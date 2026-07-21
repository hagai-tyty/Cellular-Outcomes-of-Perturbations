"""
STAGE 1a VERIFICATION — did the donor label column land, and is inner-LODO possible?

    python verify_1a.py

WHY THIS EXISTS. Stage 1b fits calibration on cross-donor statistics, produced by an inner
leave-one-donor-out over the TRAINING donors. That loop is impossible unless the training
tensors carry donor identity -- they did not, so 1a adds a 7th column sourced from the shard's
`cell_line` field.

THE LOAD-BEARING QUESTION (STAGE_1_CALIBRATION.md S1a.2). `cell_line` is a required shard column,
so it is always present -- but presence is not the same as USEFULNESS. If every row carries the
same value (e.g. one string for the whole study), there is nothing to group by and inner-LODO
cannot run.

    >= 2 distinct donors in a training split  ->  proceed to 1b
    exactly 1                                 ->  STOP. Report it. Do NOT work around it.

A wrong grouping here is the worst possible outcome: it silently produces in-distribution
calibration wearing a cross-donor label -- the exact defect Stage 1 exists to fix, now invisible.

WHAT ELSE IT CHECKS (S1a.5 acceptance):
  - all four splits return 7 tensors, not 6
  - the empty-split branch also returns 7 (easy to miss -- it builds its own tensors)
  - donor codes are integer dtype and length-matched to X

1a MUST CHANGE NO METRIC. It only adds a column. After this passes, run:
    python scorecard.py snapshot --tag 1a_donorlabels
    python scorecard.py compare baseline 1a_donorlabels
and expect EVERY metric to read 'noise'. Anything else means something other than the column moved.

USAGE (repo root, venv active; needs the cellfate_loocv_* dataset artefacts).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

from cellfate.common.io import ArtifactPaths
from cellfate.common.scalers import Scalers
from cellfate.training.dataset import (
    DONOR_I,
    DONOR_VOCAB,
    X_I,
    load_split_tensors,
)

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REGIME = "holdout"
SPLITS = ["train", "val", "calib", "test"]
N_COLS = 7
EXPECTED_TRAIN_DONORS = len(DONORS) - 1   # LOOCV: six donors minus the held-out one


def resolve_root(name: str) -> str:
    for base in (".", "runs", ".."):
        p = Path(base) / name
        if p.exists():
            return str(p)
    return name


def peek_shard(root: str) -> dict | None:
    """S1a.2 pre-check: what does `cell_line` actually contain in the raw shards?"""
    from cellfate.common import io

    paths = ArtifactPaths.of(root)
    shards = sorted(paths.shards_dir.glob("*.parquet"))
    if not shards:
        return None
    arr = io.shard_to_numpy(io.read_shard(shards[0]))
    lines = np.asarray(arr["cell_line"])
    return {
        "shard": shards[0].name,
        "n_rows": len(lines),
        "keys": sorted(arr.keys()),
        "cell_id_sample": list(arr["cell_id"][:3]),
        "cell_line_uniq": sorted(set(lines.tolist())),
    }


def check_fold(donor: str) -> dict:
    root = resolve_root(f"cellfate_loocv_{donor}")
    paths = ArtifactPaths.of(root)
    try:
        sc = Scalers.load(paths.scalers_file)
    except Exception as exc:  # noqa: BLE001
        return {"_error": repr(exc)[:120]}

    out: dict = {"splits": {}}
    for split in SPLITS:
        try:
            ds = load_split_tensors(paths, sc, REGIME, split)
        except Exception as exc:  # noqa: BLE001
            out["splits"][split] = {"_error": repr(exc)[:120]}
            continue
        d = ds.tensors[DONOR_I]
        codes = sorted(set(d.tolist()))
        # per-donor cell counts: a donor with only a handful of cells makes its inner-LODO
        # fold nearly useless, and the pooled statistics quietly inherit that
        counts = {int(c): int((d == c).sum()) for c in codes}
        out["splits"][split] = {
            "n_cols": len(ds.tensors),
            "n_rows": len(ds),
            "dtype_ok": d.dtype == torch.long,
            "len_ok": len(d) == len(ds.tensors[X_I]),
            "n_donors": len(codes),
            "codes": codes,
            "counts": counts,
        }

    # the empty-split branch builds its own tensors and is easy to miss
    try:
        empty = load_split_tensors(paths, sc, REGIME, "__no_such_split__")
        out["empty_n_cols"] = len(empty.tensors)
        out["empty_n_rows"] = len(empty)
    except Exception as exc:  # noqa: BLE001
        out["empty_n_cols"] = None
        out["empty_error"] = repr(exc)[:120]
    return out


def main() -> None:
    from cellfate.common.console import install_pretty_console, render_table

    install_pretty_console()
    # This machine's console codepage (cp1255, Hebrew) cannot encode the box-drawing
    # tables printed below; emit UTF-8 so a print can never abort the run. The JSON
    # report is written regardless of what the console can display.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    print("\nSTAGE 1a VERIFICATION — donor labels in the training tensors")
    print("the question that gates 1b: does `cell_line` distinguish donors, or is it constant?")

    # ---- S1a.2: what is actually in the shards? ----
    peek = None
    for d in DONORS:
        peek = peek_shard(resolve_root(f"cellfate_loocv_{d}"))
        if peek:
            print(f"\n  RAW SHARD PEEK  (fold {d}, {peek['shard']}, {peek['n_rows']} rows)")
            print(f"   keys           : {peek['keys']}")
            print(f"   sample cell_id : {peek['cell_id_sample']}")
            print(f"   cell_line values in this shard: {peek['cell_line_uniq']}")
            break
    if peek is None:
        print("\n   !! No shards found. Cannot verify 1a -- the dataset artefacts are missing.")
        return

    # ---- S1a.5: the acceptance assertions, per fold ----
    results, rows = {}, []
    for d in DONORS:
        r = check_fold(d)
        results[d] = r
        if "_error" in r:
            rows.append([d, "ERROR", "", "", "", r["_error"][:40]])
            continue
        tr = r["splits"].get("train", {})
        cols_ok = all(s.get("n_cols") == N_COLS
                      for s in r["splits"].values() if "_error" not in s)
        rows.append([
            d,
            str(tr.get("n_rows", "n/a")),
            str(tr.get("n_donors", "n/a")),
            "yes" if cols_ok and r.get("empty_n_cols") == N_COLS else "NO",
            "yes" if tr.get("dtype_ok") and tr.get("len_ok") else "NO",
            "OK" if tr.get("n_donors", 0) >= 2 else "** INNER-LODO IMPOSSIBLE **",
        ])

    def inv_code(code: int) -> str:
        return {v: k for k, v in DONOR_VOCAB.items()}.get(code, str(code))

    # ---- SAVE the verdict as JSON BEFORE any console table ----
    # The box-drawing tables below can crash on a non-UTF-8 console; the JSON must
    # survive that, so it is written here while `results` is fully populated.
    ok_folds = [d for d in DONORS if "_error" not in results.get(d, {"_error": 1})]
    n_donors = [results[d]["splits"].get("train", {}).get("n_donors", 0) for d in ok_folds]
    all_cols_ok = all(
        results[d].get("empty_n_cols") == N_COLS
        and all(s.get("n_cols") == N_COLS
                for s in results[d]["splits"].values() if "_error" not in s)
        for d in ok_folds
    )
    # A donor holding most of the training split is a bulk corpus, not a donor: holding it out
    # leaves a data-starved model, and its residuals then dominate the pooled calibration.
    # xdonor_calib skips these, but they must be VISIBLE here -- the first Stage 1 run was
    # invalidated by exactly this (HFF, 33,613 of 33,688 training cells).
    dominant = {}
    for d in ok_folds:
        tr = results[d]["splits"].get("train", {})
        counts, n_rows = tr.get("counts", {}), tr.get("n_rows", 0) or 1
        big = {inv_code(c): n for c, n in counts.items() if n > 0.5 * n_rows}
        if big:
            dominant[d] = big

    if not ok_folds:
        status, reason = "CANNOT_VERIFY", "No fold loaded."
    elif not all_cols_ok:
        status, reason = "FAIL", "A split returns != 7 columns (check the empty-split branch)."
    elif min(n_donors) < 2:
        bad = [d for d, n in zip(ok_folds, n_donors, strict=True) if n < 2]
        status, reason = "STOP", f"Folds {bad} have <2 training donors; inner-LODO cannot run."
    elif dominant:
        status = "STOP"
        reason = (f"`cell_line` mixes donors with a BULK CORPUS {sorted({k for v in "
                  f"dominant.values() for k in v})}: it holds >50% of the training split, so "
                  "holding it out leaves a data-starved model whose residuals swamp the pool. "
                  "Rotating over it is not cross-donor calibration.")
    elif max(n_donors) != EXPECTED_TRAIN_DONORS:
        status = "STOP"
        reason = (f"expected exactly {EXPECTED_TRAIN_DONORS} training donors "
                  f"(six minus the held-out one), saw {min(n_donors)}-{max(n_donors)}. "
                  "`cell_line` does not correspond 1:1 to donor; inspect it before running 1b.")
    else:
        status = "PASS"
        reason = (f"7 columns everywhere; exactly {EXPECTED_TRAIN_DONORS} training donors per "
                  "fold, none dominant. Inner-LODO is possible.")

    report = {
        "script": "verify_1a",
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "verdict": {
            "status": status,
            "reason": reason,
            "min_train_donors": min(n_donors) if n_donors else None,
            "max_train_donors": max(n_donors) if n_donors else None,
            "expected_train_donors": EXPECTED_TRAIN_DONORS,
            "all_splits_7_cols": all_cols_ok,
            "ok_folds": ok_folds,
            "dominant_donors": dominant,
        },
        "raw_shard_peek": peek,
        "folds": results,
    }
    out_path = Path("verify_1a_results.json")
    out_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\n  saved -> {out_path}   |   VERDICT: {status} - {reason}")

    print("\n  PER-FOLD CHECK  (train split unless noted)")
    print(render_table(
        ["fold", "train rows", "donors", "7 cols", "dtype/len", "inner-LODO"],
        rows, aligns=["l", "r", "r", "l", "l", "l"]))

    print(f"\n  DONOR VOCAB (cell_line -> code): {DONOR_VOCAB}")

    print("\n  CELLS PER DONOR in each fold's TRAIN split — a donor with very few cells makes")
    print("  its inner-LODO fold nearly useless, and the pooled calibration inherits that")
    inv = {v: k for k, v in DONOR_VOCAB.items()}
    rows = []
    for d in DONORS:
        r = results.get(d, {})
        if "_error" in r:
            continue
        counts = r["splits"].get("train", {}).get("counts", {})
        named = ", ".join(f"{inv.get(c, c)}={n}" for c, n in sorted(counts.items()))
        low = min(counts.values()) if counts else 0
        rows.append([d, str(len(counts)), named or "n/a",
                     "OK" if low >= 20 else f"** thin: {low} cells **"])
    print(render_table(["fold", "donors", "cells per donor (train)", "smallest"],
                       rows, aligns=["l", "r", "l", "l"]))

    print("\n  PER-SPLIT COLUMN COUNT — all four splits plus the empty branch must read 7")
    rows = []
    for d in DONORS:
        r = results.get(d, {})
        if "_error" in r:
            continue
        cells = []
        for split in SPLITS:
            s = r["splits"].get(split, {})
            cells.append("err" if "_error" in s else f"{s.get('n_cols')}/{s.get('n_rows')}")
        rows.append([d] + cells + [f"{r.get('empty_n_cols')}/{r.get('empty_n_rows')}"])
    print(render_table(["fold"] + [f"{s} (cols/rows)" for s in SPLITS] + ["empty branch"],
                       rows, aligns=["l"] + ["r"] * (len(SPLITS) + 1)))

    # ---- verdict ----
    ok_folds = [d for d in DONORS if "_error" not in results.get(d, {"_error": 1})]
    n_donors = [results[d]["splits"].get("train", {}).get("n_donors", 0) for d in ok_folds]
    cols_ok = all(
        results[d].get("empty_n_cols") == N_COLS
        and all(s.get("n_cols") == N_COLS
                for s in results[d]["splits"].values() if "_error" not in s)
        for d in ok_folds
    )

    print("\n   VERDICT:")
    if not ok_folds:
        print("     => CANNOT VERIFY. No fold loaded.")
    elif not cols_ok:
        print("     => FAIL. Some split still returns 6 columns. Check the empty-split branch.")
    elif min(n_donors) < 2:
        bad = [d for d, n in zip(ok_folds, n_donors, strict=True) if n < 2]
        print(f"     => STOP. Folds {bad} have <2 training donors -- inner-LODO cannot run.")
        print("        `cell_line` does not distinguish donors in this dataset. Do NOT substitute")
        print("        a guessed grouping: a wrong one produces in-distribution calibration")
        print("        wearing a cross-donor label, which is the defect Stage 1 exists to fix.")
        print("        Report this and stop -- 1b is not implementable as designed.")
    else:
        print(f"     => PASS. 7 columns everywhere; {min(n_donors)}-{max(n_donors)} training")
        print("        donors per fold. Inner-LODO is possible. Proceed to 1b.")
        if min(n_donors) < EXPECTED_TRAIN_DONORS:
            print(f"     !! FEWER than the expected {EXPECTED_TRAIN_DONORS} (6 donors minus the")
            print(f"        held-out one); saw {min(n_donors)}. A smaller inner-LODO pool means a")
            print("        noisier calibration fit -- STAGE_1 S1b.4 names this as the most likely")
            print("        cause of 1b failing.")
        if max(n_donors) > EXPECTED_TRAIN_DONORS:
            print(f"     !! MORE than the expected {EXPECTED_TRAIN_DONORS}; saw {max(n_donors)}.")
            print("        THIS IS THE DANGEROUS DIRECTION. It means `cell_line` is finer-grained")
            print("        than donor -- e.g. donor x timepoint, or donor x passage. Holding out")
            print("        one such group is NOT holding out a donor: cells from the same donor")
            print("        stay in training, so the residuals understate true cross-donor error")
            print("        and `q` comes out too small. That looks like success and is not.")
            print("        Check the cell_line values printed above before proceeding to 1b.")

    print("\n   NEXT: python scorecard.py snapshot --tag 1a_donorlabels")
    print("         python scorecard.py compare baseline 1a_donorlabels")
    print("   1a adds a column and nothing else, so EVERY metric must read 'noise'.")
    print("   Any ACCEPT or REGRESSION means something other than the column moved -- that is")
    print("   a bug, not a bonus.")


if __name__ == "__main__":
    main()
