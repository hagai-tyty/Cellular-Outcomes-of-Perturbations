"""Paired target audit: what did C-7 do to the ΔAge LABELS?  (read-only, no rebuild)

    python experiments/diag_target_shift.py

THE QUESTION
------------
C-7 does not merely drop 19 N2 labels. Rejecting five bulk columns -- one of them a control --
refits the harmonizer on a different reference, and `sigma_gill / sigma_hff` is the gain applied
to HFF's labels, which are 99.7% of the age-labelled corpus. So the network is being trained
against a different target landscape, and the ΔAge MAE comparison across C-7 is a comparison of
errors against TWO DIFFERENT TARGETS. This measures the target change directly.

Pairs every cell present in BOTH label sets, by `cell_id`, WITHIN each fold -- each fold fits its
own harmonizer on its own training donors, so `_armA` fold O1 is comparable only to `_c7t` fold
O1, never across folds.

PRE-REGISTERED READINGS (fixed before running; see the constants below)
----------------------------------------------------------------------
  A  OFFSET       slope ~ 1 and the linear fit leaves little residual
                  -> target re-centred. A centering/calibration problem, not harder biology.
  B  SCALE        slope departs from 1 but the linear fit still explains the target
                  -> a different age SCALE. A linear model re-fits a scale for free; a network
                     with a trained output head does not.
  C  NONLINEAR    the linear fit does NOT explain the change, or the shift varies with timepoint
                  -> C-7 changed the biological target along the trajectory; a genuine mismatch.

These are reported PER STRATUM, never pooled-only: HFF is 99.7% of the cells, so a pooled slope
is the HFF slope wearing a cohort's name.

NOT CLAIMED: which target is scientifically correct. C-7's justification is the degenerate
control (`integrity.py`), and it is unaffected by whatever this measures. A worse fit to cleaner
labels is not an argument for dirtier labels.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
OLD_SUFFIX = "_armA"       # pre-C-7 arm A
NEW_SUFFIX = "_c7t"        # C-7, trained
CHECK_SUFFIX = "_c7"       # dataset-only C-7; labels must match NEW_SUFFIX (runbook §3, check 3)

# ---- pre-registered thresholds ------------------------------------------------------------ #
SLOPE_TOL = 0.05      # |slope - 1| <= this  -> "slope ~ 1"
R2_FLOOR = 0.90       # linear fit explains >= this share of variance -> the change IS linear
TIME_TOL = 2.0        # yr; max spread of per-timepoint mean shift before calling it time-varying

COLS = ["cell_id", "y_age", "age_mask", "cell_line", "dose_time"]


def load_labels(root: Path) -> pd.DataFrame:
    """(shard, row) -> (y_age, age_mask, cell_line, time) for one built fold. Reads only the
    label columns, never `X` -- the feature matrix is orders of magnitude larger and unused.

    `cell_id` is NOT a key. It is an index WITHIN a chunk (`reprogramming:HFF:0` occurs once per
    shard, 45 times over), so joining on it produces a 43x many-to-many explosion that pairs
    unrelated cells and yields a plausible-looking table of pure noise. The first version of this
    script did exactly that; the tell was the free `_c7`-vs-`_c7t` check reading 73.77 instead of
    the 0 that identical configs must give.
    """
    frames = []
    for shard in sorted((root / "shards").glob("*.parquet")):
        f = pd.read_parquet(shard, columns=COLS)
        f["shard"] = shard.stem
        f["row"] = np.arange(len(f))
        frames.append(f)
    df = pd.concat(frames, ignore_index=True)
    # dose_time is [dose, time], log-transformed by the builder; the second entry separates
    # timepoints. Kept raw for grouping, with a readable day for the report.
    df["time"] = df["dose_time"].map(lambda v: float(np.asarray(v)[1]))
    df["day"] = df["time"].map(lambda t: 0.0 if t < 0 else float(np.exp(t) / 24.0))
    return df.drop(columns=["dose_time"])


def pair_labels(old: pd.DataFrame, new: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Pair cells across two builds. Returns (paired, n_unpairable).

    Positional within each (shard, timepoint) group, and ONLY where the group has the same size
    in both builds. That restriction is the whole safety argument: C-7 removes bulk columns, so a
    donor that lost one has every subsequent index shifted, and pairing it positionally would
    silently compare different samples. HFF is untouched by the bulk gate (45 shards, 42,481 rows
    in both), so its groups always match; a Gill (donor, day) group that lost a column never does
    and is dropped rather than guessed at.

    Every pair is then VERIFIED on `cell_line` and `time` -- source metadata the harmonizer never
    touches -- so a mispairing raises instead of being reported.
    """
    keys = ["shard", "time"]
    o = {k: g.sort_values("row") for k, g in old.groupby(keys)}
    n = {k: g.sort_values("row") for k, g in new.groupby(keys)}
    pairs, unpairable = [], 0
    for k, og in o.items():
        ng = n.get(k)
        if ng is None or len(ng) != len(og):
            unpairable += len(og)
            continue
        p = pd.DataFrame({
            "cell_line": og["cell_line"].to_numpy(), "time": og["time"].to_numpy(),
            "day": og["day"].to_numpy(),
            "y_old": og["y_age"].to_numpy(float), "y_c7": ng["y_age"].to_numpy(float),
            "m_old": og["age_mask"].to_numpy(bool), "m_c7": ng["age_mask"].to_numpy(bool),
            "_cl_chk": ng["cell_line"].to_numpy(), "_t_chk": ng["time"].to_numpy()})
        pairs.append(p)
    if not pairs:
        return pd.DataFrame(), unpairable
    out = pd.concat(pairs, ignore_index=True)
    if not (out["cell_line"] == out["_cl_chk"]).all() or not np.allclose(out["time"], out["_t_chk"]):
        raise AssertionError("pairing mismatch: cell_line/time disagree between builds -- the "
                             "positional assumption is wrong, refusing to report")
    return out.drop(columns=["_cl_chk", "_t_chk"]), unpairable


def ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """slope, intercept, R^2 of y ~ x. R^2 is of the FIT, so it answers 'is the change linear',
    which is not the same question as corr(x, y)."""
    if len(x) < 3 or np.std(x) == 0:
        return float("nan"), float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), r2


def verdict(slope: float, r2: float, time_spread: float) -> str:
    """The pre-registered branch. Order matters: a time-varying shift is C even if it is linear
    in the pooled fit, because 'linear overall' and 'the same change at every timepoint' are
    different claims and only the second licenses a simple re-centring."""
    if not np.isfinite(slope) or not np.isfinite(r2):
        return "UNDETERMINED"
    if np.isfinite(time_spread) and time_spread > TIME_TOL:
        return "C: NONLINEAR/TIME-VARYING"
    if r2 < R2_FLOOR:
        return "C: NONLINEAR"
    if abs(slope - 1.0) <= SLOPE_TOL:
        return "A: OFFSET"
    return "B: SCALE"


def stratum_stats(old: np.ndarray, new: np.ndarray) -> dict:
    d = new - old
    slope, intercept, r2 = ols(old, new)
    corr = float(np.corrcoef(old, new)[0, 1]) if len(old) > 2 and np.std(old) > 0 else float("nan")
    return {"n": int(len(old)),
            "mean_old": float(old.mean()), "mean_c7": float(new.mean()),
            "sd_old": float(old.std(ddof=1)) if len(old) > 1 else float("nan"),
            "sd_c7": float(new.std(ddof=1)) if len(new) > 1 else float("nan"),
            "corr": corr, "slope": slope, "intercept": intercept, "r2": r2,
            "mean_shift": float(d.mean()),
            "sd_shift": float(d.std(ddof=1)) if len(d) > 1 else float("nan")}


def audit_fold(donor: str) -> dict | None:
    old_root, new_root = ROOT / f"cellfate_loocv_{donor}{OLD_SUFFIX}", ROOT / f"cellfate_loocv_{donor}{NEW_SUFFIX}"
    if not (old_root / "shards").is_dir() or not (new_root / "shards").is_dir():
        return None
    old, new = load_labels(old_root), load_labels(new_root)

    # Runbook §3 check 3, free here: the dataset-only and trained C-7 builds are the same config
    # on the same data, so their labels must agree EXACTLY. This is also the canary that caught
    # the cell_id join explosion -- a pairing bug shows up here as a large max|Δ|.
    chk_root = ROOT / f"cellfate_loocv_{donor}{CHECK_SUFFIX}"
    c7_agrees = None
    if (chk_root / "shards").is_dir():
        cp, _ = pair_labels(load_labels(chk_root), new)
        b = cp["m_old"] & cp["m_c7"]
        c7_agrees = float(np.abs(cp.loc[b, "y_old"] - cp.loc[b, "y_c7"]).max()) if b.any() else float("nan")

    j, unpairable = pair_labels(old, new)
    both = j["m_old"] & j["m_c7"]
    paired = j[both].copy()
    o, n = paired["y_old"].to_numpy(float), paired["y_c7"].to_numpy(float)

    out: dict = {"donor": donor, "n_old": int(len(old)), "n_c7": int(len(new)),
                 "n_paired": int(len(j)), "n_unpairable": int(unpairable),
                 "n_common_labeled": int(both.sum()),
                 "n_masked_by_c7": int((j["m_old"] & ~j["m_c7"]).sum()),
                 "c7_vs_c7t_max_abs_diff": c7_agrees,
                 "overall": stratum_stats(o, n)}

    is_hff = (paired["cell_line"] == "HFF").to_numpy()
    out["by_dataset"] = {}
    for name, m in (("HFF", is_hff), ("Gill", ~is_hff)):
        if m.sum() > 2:
            out["by_dataset"][name] = stratum_stats(o[m], n[m])

    # HFF timepoints -- the trajectory C-7 is accused of bending
    hff = paired[is_hff]
    tp = {}
    for d, grp in hff.groupby("day"):
        if len(grp) > 2:
            tp[f"D{d:.0f}"] = stratum_stats(grp["y_old"].to_numpy(float),
                                            grp["y_c7"].to_numpy(float))
    out["hff_by_timepoint"] = tp
    shifts = [v["mean_shift"] for v in tp.values()]
    out["hff_timepoint_shift_spread"] = float(max(shifts) - min(shifts)) if shifts else float("nan")

    out["by_gill_donor"] = {}
    for line, grp in paired[~is_hff].groupby("cell_line"):
        if len(grp) > 2:
            out["by_gill_donor"][str(line)] = stratum_stats(grp["y_old"].to_numpy(float),
                                                            grp["y_c7"].to_numpy(float))

    h = out["by_dataset"].get("HFF")
    out["verdict_hff"] = verdict(h["slope"], h["r2"], out["hff_timepoint_shift_spread"]) if h else "UNDETERMINED"
    return out


def _row(label: str, s: dict) -> str:
    return (f"  {label:<14}{s['n']:>7}{s['mean_old']:>10.3f}{s['mean_c7']:>10.3f}"
            f"{s['sd_old']:>9.3f}{s['sd_c7']:>9.3f}{s['corr']:>8.3f}{s['slope']:>8.3f}"
            f"{s['intercept']:>10.3f}{s['r2']:>7.3f}{s['mean_shift']:>10.3f}{s['sd_shift']:>9.3f}")


HEAD = (f"  {'stratum':<14}{'n':>7}{'mean_old':>10}{'mean_c7':>10}{'sd_old':>9}{'sd_c7':>9}"
        f"{'corr':>8}{'slope':>8}{'intercept':>10}{'r2':>7}{'d_mean':>10}{'d_sd':>9}")


def main() -> None:
    print("=" * 118)
    print("PAIRED TARGET AUDIT  --  ΔAge labels, pre-C-7 (_armA)  ->  C-7 (_c7t)")
    print(f"pre-registered: |slope-1|<={SLOPE_TOL} -> A(offset) | r2<{R2_FLOOR} -> C | "
          f"timepoint shift spread >{TIME_TOL}yr -> C")
    print("=" * 118)
    results = []
    for d in DONORS:
        r = audit_fold(d)
        if r is None:
            print(f"\n[fold {d}] SKIPPED -- a build is missing")
            continue
        results.append(r)
        print(f"\n[fold {d}]  cells {r['n_old']} -> {r['n_c7']}  "
              f"(paired: {r['n_paired']}, unpairable: {r['n_unpairable']}, "
              f"newly masked: {r['n_masked_by_c7']}, paired+labelled: {r['n_common_labeled']})")
        if r["c7_vs_c7t_max_abs_diff"] is not None:
            chk = r["c7_vs_c7t_max_abs_diff"]
            flag = "OK" if chk == 0 else "!! MUST BE 0 -- PAIRING IS WRONG"
            print(f"           _c7 vs _c7t labels: max|Δ| = {chk:.3e}   [{flag}]")
        print(HEAD)
        print(_row("OVERALL", r["overall"]))
        for k, v in r["by_dataset"].items():
            print(_row(k, v))
        print(f"  HFF timepoint mean-shift spread: {r['hff_timepoint_shift_spread']:.3f} yr"
              f"   ->  VERDICT (HFF): {r['verdict_hff']}")

    if results:
        print("\n" + "=" * 118)
        print("HFF PER-TIMEPOINT MEAN SHIFT (yr), by fold  --  is the change uniform along the trajectory?")
        tps = sorted({t for r in results for t in r["hff_by_timepoint"]}, key=lambda s: float(s[1:]))
        print(f"  {'fold':<6}" + "".join(f"{t:>10}" for t in tps) + f"{'spread':>10}")
        for r in results:
            print(f"  {r['donor']:<6}"
                  + "".join(f"{r['hff_by_timepoint'][t]['mean_shift']:>10.2f}"
                            if t in r["hff_by_timepoint"] else f"{'-':>10}" for t in tps)
                  + f"{r['hff_timepoint_shift_spread']:>10.2f}")
        print("\n  VERDICTS (HFF): " + ", ".join(f"{r['donor']}={r['verdict_hff']}" for r in results))

    out = _RESULTS / "diag_target_shift_results.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
