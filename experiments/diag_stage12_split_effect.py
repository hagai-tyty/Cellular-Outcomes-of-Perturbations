"""STAGE 12 EFFECT — what the colliding `cell_id` actually did to the split, measured exactly.

Stage 12 fixed the key and deliberately claimed NOTHING about the size of the effect. This
quantifies the half that can be measured without spending a rebuild:

  MEASURABLE HERE (this script)   the split map itself -- how many decisions it really contained,
                                  and how the timepoint composition of train/val/calib differs
                                  between the colliding key and the fixed one.
  NOT MEASURABLE HERE             the effect on model metrics. That needs a rebuild + retrain +
                                  re-score under the fixed key. Pre-registered separately; it is
                                  hours of compute and is the user's to run.

No rebuild is needed for the first half, because a built fold already stores everything the split
depends on: `manifest.parquet` carries `cell_id` (the OLD, colliding key) alongside `shard_id` and
`row_idx` (which together ARE the fixed key), in build order.

THE CANARY. Re-deriving a split from stored rows is only worth anything if the derivation is
faithful, so the script first reproduces the OLD map and requires it to equal `splits/holdout.json`
**exactly**. If it does not, the reconstruction is wrong and every number below it is fiction --
the run aborts rather than reporting. This is the `diag_target_shift` lesson: a plausible-looking
number from a broken join is worse than no number.

Read-only. Reads a built fold, writes one results file.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage12_split_effect_results.json"

sys.path.insert(0, str(ROOT / "src"))

# The build parameters these folds were made with (local_runners/build_c7_folds.py:78-86).
SEED = 0
FRACS = (0.8, 0.1, 0.1, 0.0)
FOLD = "cellfate_loocv_N2_c7t"
HOLDOUT_LINE = "N2"

# `dose_time[:, 1]` is log1p(day); day 0 is the control anchor and lands on log1p(0.01)-ish.
# Taken as "the smallest distinct value" rather than a hard-coded float, so a change in the
# encoding surfaces as a different D0 count instead of silently matching nothing.
D0 = "D0"


class Row:
    """The two fields `holdout_split` reads. Avoids constructing a full validated ManifestRow
    42,600 times when only these matter."""
    __slots__ = ("cell_id", "cell_line")

    def __init__(self, cell_id: str, cell_line: str):
        self.cell_id = cell_id
        self.cell_line = cell_line


def load_fold(fold: str = FOLD) -> pd.DataFrame:
    """Manifest rows in build order, with the day of each cell joined on from its shard."""
    root = ROOT / fold
    man = pd.read_parquet(root / "manifest.parquet",
                          columns=["cell_id", "cell_line", "shard_id", "row_idx"])
    days = []
    for shard_id in man["shard_id"].drop_duplicates():
        dt = pd.read_parquet(root / "shards" / f"{shard_id}.parquet", columns=["dose_time"])
        arr = np.stack(dt["dose_time"].to_numpy())[:, 1]
        days.append(pd.DataFrame({"shard_id": shard_id, "row_idx": np.arange(len(arr)),
                                  "day_code": arr}))
    man = man.merge(pd.concat(days, ignore_index=True), on=["shard_id", "row_idx"], how="left")
    if man["day_code"].isna().any():
        raise ValueError("a manifest row has no matching shard row; the join is wrong")
    # Label the smallest day code D0 rather than trusting a literal.
    lo = man["day_code"].min()
    man["day"] = np.where(np.isclose(man["day_code"], lo), D0, "later")
    return man


def old_key(man: pd.DataFrame) -> list[str]:
    """As stored: `source:cell_line:index_within_chunk`, with no chunk."""
    return man["cell_id"].tolist()


def new_key(man: pd.DataFrame) -> list[str]:
    """Stage 12's fix: the chunk id is in the key. `shard_id` is the chunk's sanitised stem, so
    (shard_id, row_idx) carries exactly the information the fixed cell_id encodes."""
    return [f"{s}:{i}" for s, i in zip(man["shard_id"], man["row_idx"], strict=True)]


def build_split(keys: list[str], lines: list[str]) -> dict[str, str]:
    from cellfate.data.splits import holdout_split
    rows = [Row(k, ln) for k, ln in zip(keys, lines, strict=True)]
    return holdout_split(rows, {HOLDOUT_LINE}, FRACS, SEED)


def assign(man: pd.DataFrame, keys: list[str], smap: dict[str, str]) -> list[str]:
    """The split each CELL ends up in -- which is the lookup `gather_split` performs, and the
    only thing that matters. With a colliding key many cells share one entry."""
    return [smap[k] for k in keys]


def composition(man: pd.DataFrame, assigned: list[str]) -> dict:
    df = man.assign(split=assigned)
    out = {}
    for sp, grp in df.groupby("split"):
        n = len(grp)
        out[sp] = {"n": int(n),
                   "d0_n": int((grp["day"] == D0).sum()),
                   "d0_share": float((grp["day"] == D0).mean())}
    return out


def run(fold: str = FOLD) -> dict:
    man = load_fold(fold)
    lines = man["cell_line"].tolist()
    ko, kn = old_key(man), new_key(man)

    old_map = build_split(ko, lines)
    new_map = build_split(kn, lines)

    # ---- THE CANARY: the reconstruction must equal what the build actually wrote ---- #
    stored = json.loads((ROOT / fold / "splits" / "holdout.json").read_text(encoding="utf-8"))
    stored_map = stored["map"]
    canary = {"stored_entries": len(stored_map), "rebuilt_entries": len(old_map),
              "identical": stored_map == old_map}
    if not canary["identical"]:
        diff = [k for k in list(stored_map)[:2000] if old_map.get(k) != stored_map.get(k)]
        canary["first_mismatches"] = diff[:5]
        return {"fold": fold, "canary": canary, "ABORTED": True,
                "reason": "reconstruction does not reproduce the stored split map"}

    a_old, a_new = assign(man, ko, old_map), assign(man, kn, new_map)
    hff = man["cell_line"] != HOLDOUT_LINE

    return {
        "fold": fold, "canary": canary, "ABORTED": False,
        "n_cells": int(len(man)),
        "n_distinct_ids_old": int(len(set(ko))),
        "n_distinct_ids_new": int(len(set(kn))),
        "split_map_entries_old": len(old_map),
        "split_map_entries_new": len(new_map),
        "d0_cells": int((man["day"] == D0).sum()),
        "d0_distinct_ids_old": int(len({k for k, d in zip(ko, man["day"], strict=True)
                                        if d == D0})),
        "composition_old": composition(man[hff], list(np.array(a_old)[hff.to_numpy()])),
        "composition_new": composition(man[hff], list(np.array(a_new)[hff.to_numpy()])),
        "cells_that_change_split": int(sum(1 for x, y in zip(a_old, a_new, strict=True)
                                           if x != y)),
        "split_sizes_old": dict(Counter(a_old)),
        "split_sizes_new": dict(Counter(a_new)),
    }


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    r = run()

    print(f"\nSTAGE 12 EFFECT — split composition under the colliding vs fixed key ({r['fold']})")
    c = r["canary"]
    print(f"\n  CANARY  stored split map {c['stored_entries']} entries; rebuilt "
          f"{c['rebuilt_entries']}; identical: {c['identical']}")
    if r["ABORTED"]:
        print(f"\n  ABORTED — {r['reason']}")
        print(f"  first mismatches: {c.get('first_mismatches')}")
        return 1
    print("  (the reconstruction reproduces the real build exactly, so what follows is measured,")
    print("   not simulated)")

    print(f"\n  {r['n_cells']} cells carry {r['n_distinct_ids_old']} distinct ids under the old "
          f"key, {r['n_distinct_ids_new']} under the fixed one.")
    print(f"  D0: {r['d0_cells']} cells sharing {r['d0_distinct_ids_old']} ids -- so the control")
    print(f"  timepoint's split was decided {r['d0_distinct_ids_old']} times, not "
          f"{r['d0_cells']}.")

    print("\n  TIMEPOINT COMPOSITION (HFF only; the held-out donor is all test)")
    print(f"     {'split':<8}{'n (old)':>10}{'D0% (old)':>12}{'n (new)':>10}{'D0% (new)':>12}"
          f"{'shift':>9}")
    for sp in ("train", "val", "calib"):
        o = r["composition_old"].get(sp)
        n = r["composition_new"].get(sp)
        if not o or not n:
            continue
        print(f"     {sp:<8}{o['n']:>10}{o['d0_share']:>11.1%}{n['n']:>10}"
              f"{n['d0_share']:>11.1%}{n['d0_share'] - o['d0_share']:>+9.1%}")

    print(f"\n  {r['cells_that_change_split']} of {r['n_cells']} cells "
          f"({r['cells_that_change_split'] / r['n_cells']:.1%}) land in a different split.")
    print("\n  NOT CLAIMED: the effect on any model metric. That needs a rebuild + retrain +")
    print("  re-score under the fixed key, pre-registered as its own Change.")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\n  saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
