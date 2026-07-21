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
from cellfate.training.xdonor_calib import MIN_INNER_TRAIN_FRAC

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
    results = {}
    for d in DONORS:
        results[d] = check_fold(d)

    def inv_code(code: int) -> str:
        return {v: k for k, v in DONOR_VOCAB.items()}.get(code, str(code))

    ok_folds = [d for d in DONORS if "_error" not in results.get(d, {"_error": 1})]
    all_cols_ok = all(
        results[d].get("empty_n_cols") == N_COLS
        and all(s.get("n_cols") == N_COLS
                for s in results[d]["splits"].values() if "_error" not in s)
        for d in ok_folds
    )

    # Report what `crossdonor_stats` will ACTUALLY DO, not what a bare donor count suggests.
    # A donor holding most of a training split is a bulk corpus (here: HFF, 33,613 of 33,688
    # cells), and xdonor_calib skips it -- holding it out would leave a data-starved model whose
    # residuals then swamp the pool. That is what invalidated Stage 1 run 1. So the number that
    # matters is the count of donors that SURVIVE the skip, and the threshold is imported rather
    # than restated, so the two files cannot drift apart.
    bulk, usable = {}, {}
    for d in ok_folds:
        tr = results[d]["splits"].get("train", {})
        counts, n_rows = tr.get("counts", {}), tr.get("n_rows", 0) or 1
        skipped = {inv_code(c): n for c, n in counts.items()
                   if (n_rows - n) < MIN_INNER_TRAIN_FRAC * n_rows}
        if skipped:
            bulk[d] = skipped
        usable[d] = len(counts) - len(skipped)

    n_usable = list(usable.values())

    # per-fold summary rows, judged on USABLE donors -- the same number the verdict uses
    rows = []
    for d in DONORS:
        r = results[d]
        if "_error" in r:
            rows.append([d, "ERROR", "", "", "", r["_error"][:40]])
            continue
        tr = r["splits"].get("train", {})
        cols_ok = all(s.get("n_cols") == N_COLS
                      for s in r["splits"].values() if "_error" not in s)
        n_use = usable.get(d, 0)
        rows.append([
            d,
            str(tr.get("n_rows", "n/a")),
            f"{tr.get('n_donors', 0)} -> {n_use}",
            "yes" if cols_ok and r.get("empty_n_cols") == N_COLS else "NO",
            "yes" if tr.get("dtype_ok") and tr.get("len_ok") else "NO",
            "OK" if n_use == EXPECTED_TRAIN_DONORS
            else f"** {n_use} usable, expected {EXPECTED_TRAIN_DONORS} **",
        ])

    # ---- SAVE the verdict as JSON BEFORE any console table ----
    # The box-drawing tables below can crash on a non-UTF-8 console; the JSON must
    # survive that, so it is written here while `results` is fully populated.
    n_donors = [results[d]["splits"].get("train", {}).get("n_donors", 0) for d in ok_folds]
    if not ok_folds:
        status, reason = "CANNOT_VERIFY", "No fold loaded."
    elif not all_cols_ok:
        status, reason = "FAIL", "A split returns != 7 columns (check the empty-split branch)."
    elif min(n_usable, default=0) < 2:
        bad = [d for d, n in usable.items() if n < 2]
        status = "STOP"
        reason = (f"Folds {bad} have <2 donors left after bulk corpora are skipped; "
                  "inner-LODO cannot run.")
    elif min(n_usable) != EXPECTED_TRAIN_DONORS or max(n_usable) != EXPECTED_TRAIN_DONORS:
        status = "STOP"
        reason = (f"expected exactly {EXPECTED_TRAIN_DONORS} usable donors per fold (six minus "
                  f"the held-out one), saw {min(n_usable)}-{max(n_usable)} after skipping bulk "
                  "corpora. `cell_line` does not correspond 1:1 to donor; inspect it.")
    else:
        status = "PASS"
        skip_note = (f" ({sorted({k for v in bulk.values() for k in v})} skipped as bulk "
                     f"corpora, as crossdonor_stats will)" if bulk else "")
        reason = (f"7 columns everywhere; exactly {EXPECTED_TRAIN_DONORS} usable training "
                  f"donors per fold{skip_note}. Inner-LODO is possible.")

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
            "bulk_corpora_skipped": bulk,
            "usable_donors_per_fold": usable,
        },
        "raw_shard_peek": peek,
        "folds": results,
    }
    out_path = Path("verify_1a_results.json")
    out_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\n  saved -> {out_path}   |   VERDICT: {status} - {reason}")

    print("\n  PER-FOLD CHECK  (train split unless noted)")
    print(render_table(
        ["fold", "train rows", "donors (labels -> usable)", "7 cols", "dtype/len", "inner-LODO"],
        rows, aligns=["l", "r", "r", "l", "l", "l"]))

    print(f"\n  DONOR VOCAB (cell_line -> code): {DONOR_VOCAB}")

    print("\n  CELLS PER DONOR in each fold's TRAIN split, and what crossdonor_stats will do")
    print(f"  with each: a donor holding >{100 * (1 - MIN_INNER_TRAIN_FRAC):.0f}% of the split "
          "is a BULK CORPUS and gets SKIPPED -- holding")
    print("  it out would leave a data-starved model whose residuals swamp the pooled quantile")
    inv = {v: k for k, v in DONOR_VOCAB.items()}
    rows = []
    for d in DONORS:
        r = results.get(d, {})
        if "_error" in r:
            continue
        tr = r["splits"].get("train", {})
        counts, n_rows = tr.get("counts", {}), tr.get("n_rows", 0) or 1
        named = ", ".join(
            f"{inv.get(c, c)}={n}" + ("(SKIP)" if (n_rows - n) < MIN_INNER_TRAIN_FRAC * n_rows
                                      else "")
            for c, n in sorted(counts.items()))
        n_use = usable.get(d, 0)
        rows.append([d, str(len(counts)), str(n_use), named or "n/a",
                     "OK" if n_use == EXPECTED_TRAIN_DONORS
                     else f"** expected {EXPECTED_TRAIN_DONORS} **"])
    print(render_table(["fold", "labels", "usable", "cells per donor (train)", "inner-LODO"],
                       rows, aligns=["l", "r", "r", "l", "l"]))

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
    # ONE source of truth: `status`/`reason` were decided above and written to the JSON. This
    # only restates them. Recomputing here is how a console verdict and a saved verdict drift
    # apart and start disagreeing.
    print(f"\n   VERDICT: {status}")
    print(f"     => {reason}")
    if bulk:
        print(f"     - bulk corpora skipped by crossdonor_stats: "
              f"{sorted({k for v in bulk.values() for k in v})}")
        print("       (they are training signal, not donors -- rotating over them would")
        print("        calibrate against data starvation, which voided Stage 1 run 1)")
    if status == "STOP":
        print("     - Do NOT substitute a guessed grouping. A wrong one produces")
        print("       in-distribution calibration wearing a cross-donor label -- the exact")
        print("       defect Stage 1 exists to fix, now invisible. Report it instead.")

    print("\n   NEXT: python scorecard.py snapshot --tag 1a_donorlabels")
    print("         python scorecard.py compare baseline 1a_donorlabels")
    print("   1a adds a column and nothing else, so EVERY metric must read 'noise'.")
    print("   Any ACCEPT or REGRESSION means something other than the column moved -- that is")
    print("   a bug, not a bonus.")


if __name__ == "__main__":
    main()
