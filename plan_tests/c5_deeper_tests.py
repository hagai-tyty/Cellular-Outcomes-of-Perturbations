"""STAGE 1.5.3 STEP 5b — the seven things the occupancy bar does NOT see.

    python plan_tests/c5_deeper_tests.py

READ-ONLY. Writes `results/c5_deeper_tests_results.json`. Touches no `src/` file, moves no label,
runs no training. Pure simulation on the real geometry.

WHY THIS EXISTS
---------------
`register_c5_bar.py` graded B1/B2 -- *does an update get an age gradient, over how many cells* --
and both Option 1 (weighted sampler) and Option 2 (accumulation) cleared them, with Option 2 ahead.
Occupancy is **necessary and not sufficient**, and choosing on it alone would be choosing on the one
axis that happens to have been measured. Seven things it cannot see, each able to reverse the rank:

  D1  EFFECTIVE cells, not raw cells. `WeightedRandomSampler` needs `replacement=True` for weights
      to mean anything, so Option 1 can draw the SAME age cell twice in one batch. A duplicate adds
      no information to the Huber mean -- it only inflates the count B2 measured.

  D2  HOW OFTEN the age head is updated at all. Option 2 pays for its ~9 cells by stepping the age
      loss once per W batches. At W = 8 that is 8 updates per epoch instead of 65 -- an 8x cut in
      age optimisation steps that the occupancy bar is blind to, because it grades per-update
      quality and never counts the updates.

  D3  COVERAGE of the 75 labels. With only 75, a mechanism that never visits some of them in an
      epoch is training on fewer than 75, whatever the per-batch count says.

  D4  DONOR BALANCE. Those 75 labels are 5 donors (N2 14, N3 16, O2 18, Y1 13, Y2 14). Stage 1's
      whole subject is CROSS-DONOR generalisation, so a mechanism that systematically over-visits
      one donor is doing damage the occupancy bar cannot report.

  D5  Option 1's SECOND cost, which `register_c5_bar` did not measure. `replacement=True` with
      `num_samples=len(dataset)` makes every epoch a BOOTSTRAP of the training set. **This one came
      out weaker than expected and is reported as such:** the per-epoch 36% miss rate is re-rolled
      every epoch, so across 60 epochs no cell is lost (P ~ 1e-26) and the real cost is a 13% CV in
      how often each fate cell is trained on.

  D6  INFORMATION vs REPETITION -- the diagnostic that decided this. A sampler weight `w` does not
      create labels; it runs `w` age-epochs inside every fate-epoch. There are 75 labels and no
      more, so Option 1's 416 passes over the run are the same 75 labels 416 times.

  D7  NON-ZERO updates, not updates. Today's 3900 updates are not 3900 age steps: 32% carry a hard
      zero from `losses.py:55-57`, so the honest status-quo figure is 2660. Comparing Option 2's 480
      against 3900 would have overstated its cost by 47%.

A FOURTH CANDIDATE
------------------
Because D1/D2 pull in opposite directions, a hybrid is tested too: mild oversampling AND mild
accumulation, which should reach the same cells-per-update with less of each drawback. It is
included so the comparison is not a forced choice between two extremes. (It loses on its own
merits -- worst coverage and worst donor balance while still paying the full bootstrap cost.)
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


def _load_c5():
    """Reuse step 5's occupancy model rather than re-deriving it, so the two scripts cannot
    silently disagree about the geometry they are both arguing from."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "register_c5_bar", ROOT / "plan_tests" / "register_c5_bar.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["register_c5_bar"] = m
    spec.loader.exec_module(m)
    return m


C5 = _load_c5()

# real geometry -- every number from the plan's evidence table or configs/train/default.yaml
N_TRAIN = 33_688
BATCH = 512
# FULL batches only: 33688 / 512 = 65 full + one 408-cell remainder. The real DataLoader keeps
# the short batch (drop_last=False); it is dropped here because a permutation cannot fill it and
# because 1 batch in 66 changes no comparison between mechanisms.
STEPS_PER_EPOCH = N_TRAIN // BATCH            # 65
EPOCHS = 60
# the 75 surviving age labels, by donor (E12; O1 is the held-out donor in this fold)
DONORS = {"N2": 14, "N3": 16, "O2": 18, "Y1": 13, "Y2": 14}
N_AGE = sum(DONORS.values())                  # 75


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def draw_probabilities(n_age: int, n_train: int, weight: float) -> np.ndarray:
    """Per-cell draw probability under a sampler weight `w` on the age-valid cells. Pure.

    Age cells occupy indices ``[0, n_age)`` by construction here; the rest are index >= n_age.
    """
    w = np.ones(n_train, dtype=np.float64)
    w[:n_age] = weight
    return w / w.sum()


def simulate_epoch(n_age: int, n_train: int, batch: int, steps: int, *, weight: float,
                   accumulate: int, replacement: bool,
                   rng: np.random.Generator) -> dict:
    """One training epoch. Returns the five diagnostics the occupancy bar cannot see. Pure.

    ``replacement=False`` models today's plain shuffling (a permutation, every cell exactly
    once). ``replacement=True`` models `WeightedRandomSampler`, which is what makes weights
    work -- and what makes an epoch a bootstrap.
    """
    if replacement:
        p = draw_probabilities(n_age, n_train, weight)
        drawn = rng.choice(n_train, size=steps * batch, replace=True, p=p)
    else:
        drawn = rng.permutation(n_train)[: steps * batch]

    per_batch = drawn.reshape(steps, batch)
    raw, eff = [], []
    for b in per_batch:
        age = b[b < n_age]
        raw.append(age.size)
        eff.append(np.unique(age).size)          # D1: duplicates carry no extra information

    # D2: an accumulated update spans `accumulate` batches
    n_updates = steps // max(accumulate, 1)
    raw_u, eff_u = [], []
    for i in range(n_updates):
        b = per_batch[i * accumulate:(i + 1) * accumulate].ravel()
        age = b[b < n_age]
        raw_u.append(age.size)
        eff_u.append(np.unique(age).size)

    seen = np.unique(drawn[drawn < n_age])
    counts = Counter(int(x) for x in drawn[drawn < n_age])
    # D7: an update whose age loss is a hard zero moves no age parameter. What matters for the
    # age head is not how many updates ran but how many carried ANY age gradient.
    nonzero = int(sum(1 for x in raw_u if x > 0))
    return {
        "raw_cells_per_batch": float(np.mean(raw)),
        "effective_cells_per_update": float(np.mean(eff_u)) if eff_u else 0.0,
        "raw_cells_per_update": float(np.mean(raw_u)) if raw_u else 0.0,
        "age_updates_per_epoch": int(n_updates),
        "nonzero_age_updates_per_epoch": nonzero,             # D7
        "label_coverage": float(seen.size / n_age),           # D3
        "visits_per_label": counts,                           # feeds D4
        "total_age_draws": int(drawn[drawn < n_age].size),    # feeds D6
        "unique_cells_seen": int(np.unique(drawn).size),
    }


def visits_per_label_per_epoch(total_age_draws: float, n_age: int) -> float:
    """D6 -- how many times each of the 75 labels is revisited in ONE epoch. Pure.

    This is the number that decides whether a mechanism adds *information* or only *repetition*.
    A permutation gives exactly 1.0: one pass over the age set per epoch, which is what "an epoch"
    means. A sampler weight of `w` gives ~`w`, i.e. it runs `w` age-epochs inside every fate-epoch.
    There are only 75 labels and no new ones, so everything above 1.0 is the SAME labels again.
    """
    return total_age_draws / max(n_age, 1)


def bootstrap_visit_spread(n_train: int, draws_per_epoch: int, epochs: int) -> dict:
    """D5, honestly. Pure, closed form -- checked against the simulation in the unit tests.

    The per-epoch 'never drawn' figure overstates the damage, because the misses are re-rolled
    every epoch. Over the whole run a bootstrap loses essentially nobody; what it actually costs
    is VARIANCE in how often each cell is trained on. Reporting the per-epoch number alone would
    be the same overclaim this project keeps catching elsewhere, so both are reported.
    """
    p = 1.0 / n_train
    n = draws_per_epoch * epochs
    mean, sd = n * p, float(np.sqrt(n * p * (1 - p)))
    return {"expected_visits_over_run": mean,
            "sd_visits_over_run": sd,
            "cv": sd / max(mean, 1e-9),
            "p_never_seen_in_whole_run": float((1.0 - p) ** n)}


def donor_balance(visits: Counter, donors: dict[str, int]) -> dict:
    """Visit share per donor against its share of the 75 labels. Pure. D4.

    A mechanism that is fair leaves every ratio at ~1.0; a ratio of 2 means that donor's cells
    were drawn twice as often as their share of the labels warrants.
    """
    idx, out, total = 0, {}, sum(visits.values()) or 1
    for name, k in donors.items():
        v = sum(visits.get(i, 0) for i in range(idx, idx + k))
        out[name] = {"labels": k, "visits": v,
                     "share_of_visits": v / total,
                     "share_of_labels": k / sum(donors.values()),
                     "ratio": (v / total) / (k / sum(donors.values()))}
        idx += k
    ratios = [v["ratio"] for v in out.values()]
    return {"per_donor": out, "max_over_min": max(ratios) / max(min(ratios), 1e-9),
            "worst_deviation": max(abs(r - 1.0) for r in ratios)}


def bootstrap_loss(n_train: int, unique_seen: int) -> float:
    """Fraction of the training set a bootstrap epoch never shows the FATE head. Pure. D5."""
    return 1.0 - unique_seen / n_train


def main() -> int:
    rng = np.random.default_rng(0)
    print("STAGE 1.5.3 STEP 5b — what the occupancy bar could not see\n")
    print(f"  {N_AGE} age labels among {N_TRAIN} cells | batch {BATCH} | "
          f"{STEPS_PER_EPOCH} steps/epoch | {EPOCHS} epochs")
    print(f"  donors: {DONORS}\n")

    w1 = 7.1          # from register_c5_bar: the weight that puts ~8 age cells in a batch
    candidates = {
        "status quo (shuffle)":        dict(weight=1.0, accumulate=1, replacement=False),
        "Option 1 (sampler w=7.1)":    dict(weight=w1,  accumulate=1, replacement=True),
        "Option 2 (accumulate W=8)":   dict(weight=1.0, accumulate=8, replacement=False),
        "Option 4 (hybrid w=3, W=3)":  dict(weight=3.0, accumulate=3, replacement=True),
    }

    rows = {}
    for name, kw in candidates.items():
        reps = [simulate_epoch(N_AGE, N_TRAIN, BATCH, STEPS_PER_EPOCH, rng=rng, **kw)
                for _ in range(40)]
        agg = {k: float(np.mean([r[k] for r in reps]))
               for k in ("raw_cells_per_batch", "effective_cells_per_update",
                         "raw_cells_per_update", "label_coverage", "unique_cells_seen",
                         "nonzero_age_updates_per_epoch", "total_age_draws")}
        agg["age_updates_per_epoch"] = reps[0]["age_updates_per_epoch"]
        merged = Counter()
        for r in reps:
            merged.update(r["visits_per_label"])
        agg["donor_balance"] = donor_balance(merged, DONORS)
        agg["fate_cells_unseen_per_epoch"] = bootstrap_loss(N_TRAIN, agg["unique_cells_seen"])
        # the permutation candidates lose exactly the dropped tail batch and nothing else; naming
        # it keeps the bootstrap column honest instead of charging Option 2 for my truncation
        agg["fate_unseen_above_tail"] = max(
            agg["fate_cells_unseen_per_epoch"] - (N_TRAIN - STEPS_PER_EPOCH * BATCH) / N_TRAIN, 0.0)
        agg["total_age_updates"] = agg["age_updates_per_epoch"] * EPOCHS
        agg["nonzero_age_updates_total"] = agg["nonzero_age_updates_per_epoch"] * EPOCHS
        agg["visits_per_label_per_epoch"] = visits_per_label_per_epoch(
            agg["total_age_draws"], N_AGE)                                       # D6
        agg["effective_age_epochs_over_run"] = agg["visits_per_label_per_epoch"] * EPOCHS
        agg["bootstrap_spread"] = bootstrap_visit_spread(
            N_TRAIN, STEPS_PER_EPOCH * BATCH, EPOCHS) if kw["replacement"] else None
        agg["duplicate_inflation"] = (agg["raw_cells_per_update"] /
                                      max(agg["effective_cells_per_update"], 1e-9))
        rows[name] = {"config": kw, **agg}

    print(f"  {'candidate':<28}{'eff cells':>10}{'dup':>7}{'grad upd':>10}{'cover':>8}"
          f"{'donor':>7}{'reps/ep':>9}{'fate churn':>12}")
    print("  " + "-" * 91)
    for n, r in rows.items():
        print(f"  {n:<28}{r['effective_cells_per_update']:>10.2f}"
              f"{r['duplicate_inflation']:>7.2f}{r['nonzero_age_updates_total']:>10.0f}"
              f"{r['label_coverage']:>8.1%}{r['donor_balance']['max_over_min']:>7.2f}"
              f"{r['visits_per_label_per_epoch']:>9.2f}"
              f"{r['fate_unseen_above_tail']:>12.1%}")
    print("\n  eff cells  = UNIQUE age cells per optimiser update (duplicates removed)   [D1]")
    print("  dup        = raw/effective; 1.00 means no wasted duplicates               [D1]")
    print("  grad upd   = updates carrying a NON-ZERO age gradient, whole 60-ep run    [D2/D7]")
    print("  cover      = fraction of the 75 labels seen at least once per epoch       [D3]")
    print("  donor      = max/min visit ratio across the 5 donors; 1.00 is fair        [D4]")
    print("  reps/ep    = times each label is revisited per epoch; 1.00 = one pass     [D6]")
    print("  fate churn = fate cells missed per epoch ABOVE the dropped tail batch     [D5]")
    print(f"               (the tail alone is {(N_TRAIN - STEPS_PER_EPOCH * BATCH) / N_TRAIN:.1%} "
          "and is an artefact of this simulation, not of any candidate)")

    print("\n  D5 over the WHOLE run — a bootstrap re-rolls its misses every epoch:")
    for n, r in rows.items():
        s = r["bootstrap_spread"]
        if s is None:
            print(f"     {n:<28} permutation: every cell exactly {EPOCHS} times, no spread")
        else:
            print(f"     {n:<28} {s['expected_visits_over_run']:.1f} +/- "
                  f"{s['sd_visits_over_run']:.1f} visits (CV {s['cv']:.0%}), "
                  f"P(never seen at all) = {s['p_never_seen_in_whole_run']:.1e}")
    print("     ==> the bootstrap does NOT delete data; it adds VARIANCE to how often")
    print("         each fate cell is trained on. Weaker than the per-epoch number looks.")

    print("\n  D6 — does the mechanism add INFORMATION or only REPETITION?")
    for n, r in rows.items():
        print(f"     {n:<28} {r['effective_age_epochs_over_run']:>7.0f} passes over the same "
              f"75 labels across the run")
    print("     ==> there are 75 labels and no more. Anything above 60 is the same labels again.")

    # D2 is the only axis on which Option 2 loses, and W is the knob that sets it. W=8 was picked
    # to be comfortable, not because it was the smallest that works -- so measure the smallest.
    print("\n  Option 2's W, swept: the SMALLEST W clearing B2 keeps the most gradient updates")
    print(f"     {'W':>3}{'cells/upd':>11}{'B2 pass':>10}{'grad upd (run)':>16}")
    best_w = None
    for w in range(1, 13):
        sim = C5.age_cells_per_update(N_AGE, N_TRAIN, BATCH, accumulate=w, n_sim=40_000,
                                      rng=np.random.default_rng(100 + w))
        rate = float((sim >= C5.MIN_CELLS).mean())
        ok = rate >= C5.MIN_PASS_RATE
        if ok and best_w is None:
            best_w = w
        print(f"     {w:>3}{sim.mean():>11.2f}{rate:>10.1%}{(STEPS_PER_EPOCH // w) * EPOCHS:>16.0f}"
              f"   {'PASS' if ok else 'fail'}{'   <-- smallest' if w == best_w else ''}")
    if best_w is None:
        print("     ==> NO W clears B2. That would invalidate Option 2; investigate before step 6.")
        return 1

    # 75 is not a constant of nature -- it is what survives C-1 masking on THIS fold, and other
    # folds hold out other donors. A W that only just clears the bar at 75 is not a choice, it is
    # a coincidence. Measure how each candidate W survives a smaller label set.
    print("\n  Sensitivity: does the chosen W still clear B2 if the label count shrinks?")
    print(f"     {'n_age':>7}" + "".join(f"{'W=' + str(w):>10}" for w in (best_w, best_w + 1)))
    sens = {}
    for n_age in (75, 70, 65, 60, 55):
        cells = []
        for w in (best_w, best_w + 1):
            sim = C5.age_cells_per_update(n_age, N_TRAIN, BATCH, accumulate=w, n_sim=40_000,
                                          rng=np.random.default_rng(200 + w))
            cells.append(float((sim >= C5.MIN_CELLS).mean()))
        sens[n_age] = dict(zip((best_w, best_w + 1), cells, strict=True))
        marks = "".join(f"{c:>9.1%}{'*' if c >= C5.MIN_PASS_RATE else ' '}" for c in cells)
        print(f"     {n_age:>7}{marks}")
    print("     (* clears B2.) A fold that loses a donor's labels must not silently drop the")
    print("     mechanism below its own bar, so W is chosen for MARGIN, not for the minimum.")

    out = {"script": "c5_deeper_tests",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "geometry": {"n_train": N_TRAIN, "n_age": N_AGE, "batch": BATCH,
                        "steps_per_epoch": STEPS_PER_EPOCH, "epochs": EPOCHS,
                        "donors": DONORS},
           "candidates": rows,
           "option2_smallest_W_clearing_B2": best_w,
           "b2_sensitivity_to_label_count": sens,
           "b2_min_cells": C5.MIN_CELLS, "b2_bar": C5.MIN_PASS_RATE}
    (_RESULTS / "c5_deeper_tests_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote c5_deeper_tests_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
