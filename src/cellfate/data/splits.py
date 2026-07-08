"""Group-aware splits (Document 2, S12).

Random cell-level splits leak: the same drug (or cell line) ends up in both train
and test, so a model can memorise it. Each regime instead splits by *group* so
whole drugs / cell lines are held out:

* ``scaffold``  -- leave-scaffold-out (a Bemis-Murcko scaffold is train OR test)
* ``cell_line`` -- leave-cell-line-out (a cell line is train OR test)
* ``both``      -- test = unseen scaffold AND unseen line; train = seen-and-seen;
                   ambiguous rows (seen one axis, unseen the other) are dropped.

Splits are deterministic given the seed. Vehicle controls (scaffold ``CONTROL``)
are pinned to train -- they are the reference population, not a generalisation
target.
"""

from __future__ import annotations

import numpy as np

from cellfate.common.constants import Split
from cellfate.common.schemas import ManifestRow

_FOUR = (Split.TRAIN.value, Split.VAL.value, Split.CALIB.value, Split.TEST.value)
CONTROL_SCAFFOLD = "CONTROL"


def _partition_groups(groups: list[str], fracs: tuple[float, ...], seed: int) -> dict[str, str]:
    """Assign each unique group to one of train/val/calib/test by ``fracs``.

    Uses per-split integer counts (not cumulative floors, whose float noise can
    collapse adjacent cuts and starve a split). When there are at least as many
    groups as splits, every split is guaranteed >=1 group -- required so
    leave-cell-line-out always has a held-out test donor.
    """
    if len(fracs) != 4 or abs(sum(fracs) - 1.0) > 1e-6:
        raise ValueError("fracs must be four numbers summing to 1")
    uniq = sorted(set(groups))
    rng = np.random.default_rng(seed)
    ordered = [uniq[i] for i in rng.permutation(len(uniq))]
    n = len(ordered)

    raw = np.array(fracs, dtype=float) * n
    counts = np.floor(raw).astype(int)
    rem = int(n - counts.sum())                       # groups left by flooring
    for i in np.argsort(-(raw - counts))[:rem]:       # give them to largest remainders
        counts[i] += 1
    if n >= len(counts):                              # guarantee no split is empty
        for i in range(len(counts)):
            if counts[i] == 0:
                counts[int(counts.argmax())] -= 1
                counts[i] += 1

    assign: dict[str, str] = {}
    start = 0
    for split_name, cnt in zip(_FOUR, counts, strict=True):
        for g in ordered[start:start + int(cnt)]:
            assign[g] = split_name
        start += int(cnt)
    for g in ordered[start:]:                         # any leftover -> train
        assign[g] = Split.TRAIN.value
    return assign



def _held_out(groups: list[str], frac: float, seed: int) -> set[str]:
    """Return the subset of unique groups held out (size ~= frac)."""
    uniq = sorted(set(groups))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(uniq))
    k = int(np.floor(frac * len(uniq)))
    return {uniq[i] for i in perm[:k]}


def scaffold_split(rows: list[ManifestRow], fracs: tuple[float, ...], seed: int) -> dict[str, str]:
    scaffolds = [r.scaffold_id or r.pert_id for r in rows]
    assign = _partition_groups([s for s in scaffolds if s != CONTROL_SCAFFOLD], fracs, seed)
    assign[CONTROL_SCAFFOLD] = Split.TRAIN.value  # controls are the reference
    return {r.cell_id: assign[scaffolds[i]] for i, r in enumerate(rows)}


def cell_line_split(rows: list[ManifestRow], fracs: tuple[float, ...], seed: int) -> dict[str, str]:
    assign = _partition_groups([r.cell_line for r in rows], fracs, seed)
    return {r.cell_id: assign[r.cell_line] for r in rows}


def both_split(rows: list[ManifestRow], fracs: tuple[float, ...], seed: int) -> dict[str, str]:
    test_frac = fracs[3]
    held_scaffold = _held_out([r.scaffold_id or r.pert_id for r in rows], test_frac, seed + 1)
    held_line = _held_out([r.cell_line for r in rows], test_frac, seed + 2)
    # Sub-split the seen-and-seen rows into train/val/calib (renormalised fracs).
    seen_fracs = np.array(fracs[:3], dtype=float)
    seen_fracs = tuple(seen_fracs / seen_fracs.sum()) + (0.0,)
    seen_assign = _partition_groups(
        [r.scaffold_id or r.pert_id for r in rows
         if (r.scaffold_id or r.pert_id) not in held_scaffold and r.cell_line not in held_line],
        seen_fracs, seed + 3,
    )
    out: dict[str, str] = {}
    for r in rows:
        scaf = r.scaffold_id or r.pert_id
        s_unseen, l_unseen = scaf in held_scaffold, r.cell_line in held_line
        if s_unseen and l_unseen:
            out[r.cell_id] = Split.TEST.value
        elif not s_unseen and not l_unseen:
            out[r.cell_id] = seen_assign.get(scaf, Split.TRAIN.value)
        else:
            out[r.cell_id] = Split.DROP.value
    return out


def random_split(rows: list[ManifestRow], fracs: tuple[float, ...], seed: int) -> dict[str, str]:
    """Assign each CELL at random to train/val/calib/test by ``fracs``.

    Unlike the group-based regimes, this does NOT test generalization -- the same
    cell line / perturbation can appear in train and test. It exists for datasets
    with a single donor + single perturbation (e.g. a one-line reprogramming time
    course) where the group regimes are degenerate: it validates that the model
    *fits* the per-cell fate manifold. Generalization claims still require the
    group regimes across multiple lines/datasets.
    """
    rng = np.random.default_rng(seed)
    splits = [Split.TRAIN.value, Split.VAL.value, Split.CALIB.value, Split.TEST.value]
    probs = np.array(fracs, dtype=float)
    probs = probs / probs.sum()
    draws = rng.choice(len(splits), size=len(rows), p=probs)
    return {r.cell_id: splits[draws[i]] for i, r in enumerate(rows)}


_REGIMES = {
    "scaffold": scaffold_split,
    "cell_line": cell_line_split,
    "both": both_split,
    "random": random_split,
}


def holdout_split(
    rows: list[ManifestRow],
    holdout_cell_lines: set[str],
    fracs: tuple[float, ...],
    seed: int,
) -> dict[str, str]:
    """Leave-one-(or-more)-cell-line-out for the TEST set only.

    The named cell lines go entirely to ``test`` (a leak-free generalization probe:
    the model never trains on them). Every *other* cell is split at the cell level
    into train/val/calib, so a large single-cell line (e.g. HFF) both trains the
    model AND provides a real calibration set -- instead of one giant group landing
    wholesale in a single split and starving the others.
    """
    rng = np.random.default_rng(seed)
    tvc = np.array(fracs[:3], dtype=float)
    tvc = tvc / tvc.sum()                       # renormalise train/val/calib (drop test frac)
    names = [Split.TRAIN.value, Split.VAL.value, Split.CALIB.value]
    out: dict[str, str] = {}
    for r in rows:
        if r.cell_line in holdout_cell_lines:
            out[r.cell_id] = Split.TEST.value
        else:
            out[r.cell_id] = names[int(rng.choice(3, p=tvc))]
    return out


def make_splits(
    rows: list[ManifestRow],
    fracs: tuple[float, ...],
    regimes: tuple[str, ...],
    seed: int,
    holdout_cell_lines: tuple[str, ...] = (),
) -> dict[str, dict[str, str]]:
    """Build the cell_id -> split mapping for each requested regime.

    The pseudo-regime ``"holdout"`` uses :func:`holdout_split` with
    ``holdout_cell_lines`` (those cell lines -> test, the rest cell-level split).
    """
    known = set(_REGIMES) | {"holdout"}
    unknown = [r for r in regimes if r not in known]
    if unknown:
        raise ValueError(f"unknown regimes {unknown}; choose from {sorted(known)}")
    out: dict[str, dict[str, str]] = {}
    for regime in regimes:
        if regime == "holdout":
            out[regime] = holdout_split(rows, set(holdout_cell_lines), fracs, seed)
        else:
            out[regime] = _REGIMES[regime](rows, fracs, seed)
    return out
